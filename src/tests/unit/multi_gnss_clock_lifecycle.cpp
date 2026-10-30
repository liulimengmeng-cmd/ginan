// Controlled production SPP/clock-offset lifecycle, not a real-data PPP run.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/ephemeris.hpp"
#include "common/navigation.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/zhangSystemClockDatum.hpp"
#include "common/zhangEffectiveClockDatum.hpp"
#include "common/zhangCheckpoint.hpp"
#include "orbprop/coordinates.hpp"
#include "pea/ppp.hpp"
#include "pea/zhangReceiverCheckpoint.hpp"
#include "common/sinex.hpp"
#include "trop/tropModels.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

void updateRecClocks(Trace&, ReceiverMap&, KFState&);
void updateSppPppClkOffsets(ReceiverMap&, KFState&);

namespace
{
unsigned checks = 0;
void require(bool pass, const std::string& message)
{ ++checks; if (!pass) throw std::runtime_error(message); }
struct System { E_Sys sys; double origin; E_ObsCode first, second; int firstPrn; };
const std::vector<System> allSystems{
    {E_Sys::GPS, 100, E_ObsCode::L1C, E_ObsCode::L2W, 2},
    {E_Sys::GAL, 650, E_ObsCode::L1C, E_ObsCode::L5Q, 11},
    {E_Sys::BDS, 102722, E_ObsCode::L1P, E_ObsCode::L5P, 19}};
const std::vector<Vector3d> positions{
    {26500000,0,0}, {24500000,12000000,0}, {24500000,-12000000,5000000},
    {22500000,0,-15000000}, {20500000,10000000,15000000}, {20500000,-10000000,-12000000}};
constexpr long firstEpoch = 1405216800;
constexpr double receiverRelativeClock = 12;

void configure(const std::vector<System>& systems)
{
    acsConfig.zhangFullRank.enable = true;
    acsConfig.zhangPppAr.user_adapter = false;
    acsConfig.zhangPppAr.output_products = false;
    acsConfig.process_sys.clear();
    acsConfig.zhangFullRank.sysOpts.clear();
    acsConfig.code_priorities.clear();
    for (const auto& system : systems)
    {
        acsConfig.process_sys[system.sys] = true;
        auto& opts = acsConfig.zhangFullRank.sysOpts[system.sys];
        opts.baseline_observables = {system.first, system.second};
        opts.reference_receiver = "PIVOT";
        acsConfig.code_priorities[system.sys] = {system.first, system.second};
    }
    acsConfig.adjust_rec_clocks_by_spp = true;
    acsConfig.adjust_clocks_for_jumps_only = false;
    acsConfig.preprocOpts.preprocess_all_data = false;
    receiverMap.clear();
}

void observations(Receiver& receiver, const std::vector<System>& systems,
                  long epoch, double physicalJump, bool prefilledKfChart)
{
    receiver.obsList.clear();
    receiver.aprioriPos << 6378137, 0, 0;
    auto& recOpts = acsConfig.getRecOpts(receiver.id);
    recOpts.receiver_reference_system = systems.front().sys;
    recOpts.clk.estimate = {true}; recOpts.clk.sigma = {1000}; recOpts.clk.process_noise = {0};
    recOpts.clk_rate.estimate = {false};
    VectorPos pos = ecef2pos(receiver.aprioriPos);
    std::ostringstream discardedTrace;
    for (const auto& system : systems)
        for (int satellite = 0; satellite < static_cast<int>(positions.size()); ++satellite)
        {
            GObs obs;
            obs.Sat = SatSys(system.sys, system.firstPrn + satellite);
            obs.time.bigTime = epoch; obs.mount = receiver.id; obs.rec_ptr = &receiver;
            obs.satNav_ptr = &nav.satNavMap[obs.Sat];
            obs.satStat_ptr = &receiver.satStatMap[obs.Sat];
            obs.rSatApc = positions[satellite]; obs.rSatCom = positions[satellite];
            obs.satVel.setZero(); obs.ephPosValid = obs.ephClkValid = true;
            obs.posSource = E_Source::PRECISE;
            obs.clkSource = prefilledKfChart ? E_Source::KALMAN : E_Source::PRECISE;
            const double externalClock = 10000 + static_cast<int>(system.sys) * 100 + 10 * satellite;
            obs.satClk = (externalClock - (prefilledKfChart ? system.origin : 0)) / CLIGHT;
            obs.satClkVel = 0; obs.satClkVar = 0;
            Vector3d lineOfSight;
            const double range = geodist(obs.rSatApc, receiver.aprioriPos, lineOfSight);
            satazel(pos, lineOfSight, *obs.satStat_ptr);
            TropStates trop; TropMapping mapping; double tropVariance = 0;
            const double delay = tropModel(discardedTrace, acsConfig.sppOpts.trop_models,
                obs.time, pos, *obs.satStat_ptr, trop, mapping, tropVariance);
            const double rawRange = range + system.origin + receiverRelativeClock + physicalJump - externalClock + delay;
            for (const auto code : {system.first, system.second})
            {
                const auto frequency = code2Freq.at(system.sys).at(code);
                const double lambda = genericWavelength.at(frequency);
                obs.satNav_ptr->lamMap[frequency] = lambda;
                auto& sig = obs.sigs[frequency];
                sig.code = code; sig.P = rawRange; sig.L = rawRange / lambda;
                sig.snr = 50;
            }
            auto& satOpts = acsConfig.getSatOpts(obs.Sat);
            satOpts.clockModel.sources = {E_Source::KALMAN, E_Source::PRECISE};
            satOpts.posModel.sources = {E_Source::PRECISE};
            // External clocks actually exist at both epochs. A policy which
            // stabilizes SPP must handle already-valid preprocessor KF clocks.
            for (long offset : {-30, 0, 30, 60, 90})
            {
                GTime time; time.bigTime = firstEpoch + offset;
                Pclk clock; clock.clk = externalClock / CLIGHT; clock.clkStd = 0; clock.clkIndex = 0;
                nav.pclkMap[obs.Sat.id()][time] = clock;
            }
            // SPP is allowed to recompute the complete physical satellite
            // cache. Supply actual precise ephemeris samples, not merely a
            // prefilled position which a source-stabilizing policy must trust.
            for (int sample = -10; sample <= 10; ++sample)
            {
                Peph ephemeris; ephemeris.Sat = obs.Sat;
                ephemeris.time.bigTime = firstEpoch + 30*sample;
                ephemeris.index = ephemeris.clkIndex = 0;
                ephemeris.pos = positions[satellite]; ephemeris.vel.setZero();
                ephemeris.clk = externalClock / CLIGHT;
                nav.pephMap[obs.Sat.id()][ephemeris.time] = ephemeris;
            }
            receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
        }
}

KFState clockState(Receiver& receiver, const std::vector<System>& systems, Trace& trace)
{
    KFState state;
    state.time.bigTime = firstEpoch;
    auto add = [&](KFKey key, double value)
    { InitialState initial; initial.estimate = true; initial.x = value; initial.P = 1; initial.Q = 0; state.addKFState(key, initial); };
    for (const auto& system : systems)
    {
        add(zhangReceiverClockKey(acsConfig, receiver.id, system.sys), receiverRelativeClock);
        for (int satellite = 0; satellite < static_cast<int>(positions.size()); ++satellite)
        {
            KFKey key; key.type = KF::SAT_CLOCK;
            key.Sat = SatSys(system.sys, system.firstPrn + satellite);
            add(key, 10000 + static_cast<int>(system.sys) * 100 + 10 * satellite - system.origin);
        }
    }
    for (int axis = 0; axis < 3; ++axis)
    { KFKey key; key.type = KF::REC_POS; key.str = receiver.id; key.num = axis; add(key, receiver.aprioriPos[axis]); }
    require(state.stateTransition(trace, state.time), "fixture clock state initialization failed");
    MatrixXd factor = MatrixXd::Zero(state.x.size(), state.x.size());
    const int one = state.kfIndexMap.at(KFState::oneKey);
    for (int row = 0; row < factor.rows(); ++row)
        for (int column = 0; column < factor.cols(); ++column)
            if (row != one && column != one)
                factor(row,column) = (row == column ? 2.0 : 0) + .003*(row+1)*(column+1);
    state.P = factor * factor.transpose();
    return state;
}

bool lifecycle(Trace& trace, const std::vector<System>& systems, double jump, bool observeOnly)
{
    configure(systems);
    const std::string id = "CLOCK_LIFECYCLE_" + enum_to_string(systems.front().sys) +
        (systems.size() > 1 ? "_JOINT" : "_SINGLE") + (jump ? "_JUMP" : "_NO_JUMP");
    auto& receiver = receiverMap[id]; receiver.id = id; receiver.ready = true;
    tsync.bigTime = firstEpoch;
    observations(receiver, systems, firstEpoch, 0, false);
    KFState state = clockState(receiver, systems, trace);
    // With external clocks on the first cold SPP, no KF source is offered.
    spp(trace, receiver, nullptr, nullptr);
    require(receiver.sol.status == E_Solution::SINGLE, "controlled cold SPP failed: " + id);
    const KFState coldSpp = receiver.sol.sppState;
    updateSppPppClkOffsets(receiverMap, state);
    require(receiver.sol.clkAdjustReady, "postfilter did not bind SPP offset for " + id);
    for (const auto& system : systems)
    {
        double offset = 0;
        require(zhangSavedSppClockOffset(state,receiver,system.sys,offset) && std::abs(offset+system.origin) < 1e-4,
            "postfilter did not save each system's distinct current external-clock offset");
        KFState wrongEpoch = state;
        wrongEpoch.time.bigTime = firstEpoch-30;
        require(!zhangSavedSppClockOffset(wrongEpoch,receiver,system.sys,offset),
            "saved offset was read against a different committed KF epoch");
        KFState malformed = state;
        malformed.metaDataMap[zhangClockLifecycleMetadataKey("SPP_OFFSET",system.sys,id)] = "not-json";
        require(!zhangSavedSppClockOffset(malformed,receiver,system.sys,offset),
            "malformed source-chart offset was accepted");
    }
    KFState staleSave = state;
    staleSave.time.bigTime = firstEpoch-30;
    for (const auto& system : systems)
        staleSave.metaDataMap.erase(zhangClockLifecycleMetadataKey("SPP_OFFSET",system.sys,id));
    ReceiverMap staleReceivers = receiverMap;
    updateSppPppClkOffsets(staleReceivers,staleSave);
    for (const auto& system : systems)
    {
        double offset = 0;
        require(!zhangSavedSppClockOffset(staleSave,staleReceivers.at(id),system.sys,offset),
            "current SPP was bound to a different-epoch committed FLOAT posterior");
    }
    tsync.bigTime = firstEpoch+30;
    for (const auto& system : systems)
        require(zhangSppPhysicalClockChart(receiver.sol,system.sys).empty(),
            "preceding-epoch SPP chart was reused as current evidence");
    tsync.bigTime = firstEpoch;
    const double savedOffset = receiver.sol.sppPppClkOffset;
    const auto selected = receiver.sol.clkRefSys;
    const auto selectedKey = zhangReceiverClockKey(acsConfig, id, selected);
    const KFState before = state;
    tsync.bigTime = firstEpoch + 30;
    observations(receiver, systems, firstEpoch + 30, jump, true);
    spp(trace, receiver, &state, nullptr);
    require(receiver.sol.status == E_Solution::SINGLE, "controlled second SPP failed: " + id);
    for (const auto& observation : only<GObs>(receiver.obsList))
        require(observation.clkSource == E_Source::KALMAN && observation.ephClkValid && observation.sppValid,
            "independent SPP did not preserve the caller's cache or actual SPP-valid result");
    const auto priorCovariance = coldSpp.P;
    require(coldSpp.kfIndexMap == receiver.sol.sppState.kfIndexMap, "SPP state order changed on a clock chart-only event");
    require((priorCovariance - receiver.sol.sppState.P).norm() < 1e-7 * std::max(1.0, priorCovariance.norm()),
        "constant clock chart/jump changed full SPP geometry covariance");
    updateRecClocks(trace, receiverMap, state);
    double applied = 0;
    if (state.stateTransitionMap.contains(selectedKey) && state.stateTransitionMap.at(selectedKey).contains(KFState::oneKey) &&
        state.stateTransitionMap.at(selectedKey).at(KFState::oneKey).contains(0))
        applied = state.stateTransitionMap.at(selectedKey).at(KFState::oneKey).at(0);
    bool correct = std::abs(applied-jump) < 1e-4;
    // This injection is a physical receiver-clock jump common to the signals,
    // not a one-constellation ISB injection. Check every effective clock rather
    // than silently accepting a GPS-only adjustment in a joint solution.
    for (const auto& system : systems)
    {
        const auto key = zhangReceiverClockKey(acsConfig, id, system.sys);
        double systemApplied = 0;
        if (state.stateTransitionMap.contains(key) && state.stateTransitionMap.at(key).contains(KFState::oneKey) &&
            state.stateTransitionMap.at(key).at(KFState::oneKey).contains(0))
            systemApplied = state.stateTransitionMap.at(key).at(KFState::oneKey).at(0);
        correct = correct && std::abs(systemApplied-jump) < 1e-4;
    }
    std::cout << std::setprecision(17) << "CLOCK_LIFECYCLE system=" << enum_to_string(selected)
        << " systems=" << systems.size() << " source_chart_switch=EXTERNAL_TO_PREFILLED_VALID_KF"
        << " physical_jump_m=" << jump << " saved_offset_m=" << savedOffset
        << " second_spp_clock_m=" << receiver.sol.sppClk << " applied_adjustment_m=" << applied
        << " expected_adjustment_m=" << jump << " prediction_equivalent=" << correct << '\n';
    if (!correct && observeOnly) return false;
    require(correct, "source chart motion was reapplied or a real receiver jump absorbed: " + id);
    require(state.stateTransition(trace, tsync), "production clock transition failed");
    require(state.kfIndexMap == before.kfIndexMap, "clock adjustment changed global state order");
    require((state.P-before.P).norm() < 1e-10 * std::max(1.0,before.P.norm()),
        "deterministic clock-origin handling changed complete KF covariance");
    for (const auto& [key,index] : state.kfIndexMap)
    {
        const bool receiverClock = key.type == KF::REC_CLOCK && key.str == id && key.num == 0;
        const double expected = before.x[index] + (receiverClock ? jump : 0);
        require(std::abs(state.x[index]-expected) < 1e-4, "clock handling modified a foreign system/nuisance mean");
    }
    MatrixXd design = MatrixXd::Zero(systems.size()*positions.size(), state.x.size());
    int row = 0;
    for (const auto& system : systems)
        for (int satellite = 0; satellite < static_cast<int>(positions.size()); ++satellite, ++row)
        {
            design(row,state.kfIndexMap.at(zhangReceiverClockKey(acsConfig,id,system.sys))) = 1;
            KFKey satelliteKey; satelliteKey.type = KF::SAT_CLOCK;
            satelliteKey.Sat = SatSys(system.sys,system.firstPrn+satellite);
            design(row,state.kfIndexMap.at(satelliteKey)) = -1;
            Vector3d los = (positions[satellite]-receiver.aprioriPos).normalized();
            for (int axis = 0; axis < 3; ++axis)
            {
                KFKey positionKey; positionKey.type = KF::REC_POS; positionKey.str = id; positionKey.num = axis;
                design(row,state.kfIndexMap.at(positionKey)) = -los[axis];
            }
        }
    require(((design*state.x-design*before.x).array()-jump).abs().maxCoeff() < 1e-4,
        "full multi-system range prediction lost a common physical receiver jump");
    require((design*state.P*design.transpose()-design*before.P*design.transpose()).norm() <
        1e-10*std::max(1.0,(design*before.P*design.transpose()).norm()),
        "clock lifecycle changed cross-system predicted observation covariance");
    return true;
}

bool legacyCheckpointBoundary(Trace& trace, bool observeOnly)
{
    configure({allSystems.front()});
    const std::string id = "LEGACY_GPS_CLOCK_CHECKPOINT";
    auto& receiver = receiverMap[id]; receiver.id = id; receiver.ready = true;
    tsync.bigTime = firstEpoch;
    observations(receiver,{allSystems.front()},firstEpoch,0,false);
    KFState legacy = clockState(receiver,{allSystems.front()},trace);
    spp(trace,receiver,nullptr,nullptr);
    require(receiver.sol.status == E_Solution::SINGLE,"legacy fixture SPP failed");
    updateSppPppClkOffsets(receiverMap,legacy);
    // Model an actual pre-V2 checkpoint: neither the FLOAT core nor its SPP
    // runtime had source-chart provenance. Retaining freshly generated V2
    // metadata here would not test a legacy import at all.
    const std::string lifecyclePrefix = "ZHANG_CLOCK_LIFECYCLE_V2|";
    for (auto* state : {&legacy,&receiver.sol.sppState})
        for (auto entry = state->metaDataMap.begin(); entry != state->metaDataMap.end();)
            entry = entry->first.rfind(lifecyclePrefix,0) == 0 ? state->metaDataMap.erase(entry) : std::next(entry);
    auto& site = theSinex.mapsiteids[id]; site.sitecode = id;
    auto& equipment = theSinex.mapreceivers[id][GTime{}]; equipment.sitecode = id;
    auto& antenna = theSinex.mapantennas[id][GTime{}]; antenna.sitecode = id;
    auto& eccentricity = theSinex.mapeccentricities[id][GTime{}]; eccentricity.sitecode = id;
    eccentricity.ecc = Vector3d::Zero();
    receiver.snx.id_ptr = &site; receiver.snx.rec_ptr = &equipment;
    receiver.snx.ant_ptr = &antenna; receiver.snx.ecc_ptr = &eccentricity;
    std::string payload;
    auto exported = exportZhangReceiverRuntimeCheckpointSection(receiverMap,"CLOCK-LIFECYCLE-LEGACY",payload);
    require(exported.valid,"legacy checkpoint export: "+exported.failureReason);
    ReceiverMap restored = receiverMap;
    ZhangReceiverRuntimeCheckpointRestorePlan plan;
    auto preflight = preflightZhangReceiverRuntimeCheckpointSection(restored,"CLOCK-LIFECYCLE-LEGACY",payload,plan);
    require(preflight.valid,"legacy checkpoint preflight: "+preflight.failureReason);
    auto imported = importZhangReceiverRuntimeCheckpointSection(restored,"CLOCK-LIFECYCLE-LEGACY",plan);
    require(imported.valid,"legacy checkpoint import: "+imported.failureReason);
    auto& resumed = restored.at(id);
    // Import correctly invalidates derived option caches. Reinstall the same
    // fixture's hashed model options; they are not receiver-runtime payload.
    auto& resumedOptions = acsConfig.getRecOpts(id);
    resumedOptions.clk.estimate = {true}; resumedOptions.clk.sigma = {1000};
    resumedOptions.clk.process_noise = {0}; resumedOptions.clk_rate.estimate = {false};
    require(resumed.sol.sppPppClkOffsetSys == E_Sys::NONE,
        "legacy wire restore invented constellation provenance");
    const auto legacyKey = zhangReceiverClockKey(acsConfig,id,E_Sys::GPS);
    require(legacyKey.Sat.sys == E_Sys::NONE && legacy.kfIndexMap.contains(legacyKey),
        "single GPS checkpoint no longer uses original NONE receiver key");
    tsync.bigTime = firstEpoch+30;
    resumed.sol.sppClk += 17.25;
    const KFState before = legacy;
    updateRecClocks(trace,restored,legacy);
    double unknownApplied = 0;
    if (legacy.stateTransitionMap.contains(legacyKey) && legacy.stateTransitionMap.at(legacyKey).contains(KFState::oneKey) &&
        legacy.stateTransitionMap.at(legacyKey).at(KFState::oneKey).contains(0))
        unknownApplied = legacy.stateTransitionMap.at(legacyKey).at(KFState::oneKey).at(0);
    const bool unknownRejected = std::abs(unknownApplied) < 1e-12;
    std::cout << "CLOCK_LIFECYCLE legacy_unknown_source_chart_first_adjustment_m=" << unknownApplied
        << " correctly_skipped=" << unknownRejected << '\n';
    if (!observeOnly) require(unknownRejected,"legacy scalar offset lacking source-chart evidence was applied");
    ReceiverMap legacyImported = restored;
    legacy = before;
    // After one fresh external SPP/postfilter binding, subsequent physical
    // jumps must again work. Skipping forever is not a valid legacy policy.
    observations(resumed,{allSystems.front()},firstEpoch+30,0,false);
    spp(trace,resumed,nullptr,nullptr);
    require(resumed.sol.status == E_Solution::SINGLE,"legacy rebind SPP failed");
    require(legacy.stateTransition(trace,tsync),"legacy initial skipped transition failed");
    updateSppPppClkOffsets(restored,legacy);
    tsync.bigTime = firstEpoch+60;
    observations(resumed,{allSystems.front()},firstEpoch+60,17.25,false);
    spp(trace,resumed,nullptr,nullptr);
    require(resumed.sol.status == E_Solution::SINGLE,"legacy rebound jump SPP failed");
    updateRecClocks(trace,restored,legacy);
    require(legacy.stateTransition(trace,tsync),"legacy GPS rebound transition failed");
    require(std::abs(legacy.x[legacy.kfIndexMap.at(legacyKey)]-
        before.x[before.kfIndexMap.at(legacyKey)]-17.25) < 1e-4,
        "freshly rebound legacy single-GPS offset discarded a real receiver jump");
    require((legacy.P-before.P).norm() < 1e-10*std::max(1.0,before.P.norm()),
        "legacy deterministic resume altered full covariance");

    // The same legacy NONE state and provenance-less scalar are not three
    // system clocks when a checkpoint is presented to a joint configuration.
    configure(allSystems);
    for (const auto& system : allSystems)
        require(!zhangReceiverClockMatchesSystem(acsConfig,legacyKey,system.sys),
            "joint configuration silently imported a legacy NONE clock");
    KFState jointBoundary = before;
    const auto originalMap = jointBoundary.kfIndexMap;
    const auto originalTransition = jointBoundary.stateTransitionMap;
    updateRecClocks(trace,legacyImported,jointBoundary);
    require(jointBoundary.kfIndexMap == originalMap && jointBoundary.stateTransitionMap == originalTransition,
        "provenance-less legacy SPP offset manufactured a joint-system transition");
    require((jointBoundary.x-before.x).norm() == 0 && (jointBoundary.P-before.P).norm() == 0,
        "joint legacy boundary mutated means or covariance");
    std::cout << "CLOCK_LIFECYCLE legacy_receiver_checkpoint_actual_export_preflight_import=1"
                 " single_GPS_NONE_rebind_resume=PASS joint_legacy_import=REJECTED\n";
    return unknownRejected;
}

void configureClockOnlyReceiverModel(ReceiverOptions& opts)
{
    opts.clk.estimate = {true}; opts.clk.sigma = {1000}; opts.clk.process_noise = {0};
    opts.clk_rate.estimate = {false}; opts.pos.estimate = {false,false,false};
    opts.clockModel.enable = true; opts.eccentricityModel.enable = false;
    opts.pcoModel.enable = opts.pcvModel.enable = false;
    opts.tideModels.enable = false; opts.tropModel.enable = false;
    opts.relativity = opts.relativity2 = false;
    opts.ionospheric_component = opts.ionospheric_component2 = opts.ionospheric_component3 = false;
    opts.ionospheric_model = false; opts.eop = false;
    opts.phaseBiasModel.enable = opts.codeBiasModel.enable = false;
}

void configureClockOnlySatelliteModel(SatelliteOptions& opts)
{
    opts.clk.estimate = {true}; opts.clk.sigma = {1000}; opts.clk.process_noise = {0};
    opts.clk_rate.estimate = {false}; opts.clk_rate.apriori_value = {0};
    opts.clk_rate.process_noise = {0}; opts.clockModel.enable = true;
    opts.posModel.sources = {E_Source::PRECISE};
    opts.clockModel.sources = {E_Source::KALMAN,E_Source::PRECISE};
    opts.pcoModel.enable = opts.pcvModel.enable = false;
    opts.phaseBiasModel.enable = opts.codeBiasModel.enable = false;
}

void configureClockOnlyRows(Receiver& receiver)
{
    configureClockOnlyReceiverModel(acsConfig.getRecOpts(receiver.id));
    for (auto& observation : only<GObs>(receiver.obsList))
    {
        configureClockOnlySatelliteModel(acsConfig.getSatOpts(observation.Sat));
        for (auto& [frequency,signal] : observation.sigs)
        {
            signal.codeVar = 1;
            observation.sigsLists[frequency] = {signal};
            const std::string name = enum_to_string(signal.code);
            configureClockOnlyReceiverModel(acsConfig.getRecOpts(receiver.id,{observation.Sat.sysName(),name}));
            configureClockOnlySatelliteModel(acsConfig.getSatOpts(observation.Sat,{name}));
        }
    }
}

bool movingOrbitObservationModel(Trace& trace, const std::vector<System>& systems, bool observeOnly)
{
    configure(systems);
    const std::string id = "CLOCK_PHYSICAL_EFFECTIVE_" + enum_to_string(systems.front().sys) +
        (systems.size() > 1 ? "_JOINT" : "_SINGLE");
    auto& receiver = receiverMap[id]; receiver.id = id; receiver.ready = true;
    tsync.bigTime = firstEpoch;
    observations(receiver,systems,firstEpoch,0,false);
    const Vector3d velocity(1100,2300,-1400);
    configureClockOnlyRows(receiver);
    constexpr double configuredSatelliteRate = 0.37;
    auto enableRate = [configuredSatelliteRate](SatelliteOptions& opts)
    { opts.clk_rate.estimate = {true}; opts.clk_rate.apriori_value = {configuredSatelliteRate};
      opts.clk_rate.sigma = {2}; opts.clk_rate.process_noise = {0}; };
    for (auto& observation : only<GObs>(receiver.obsList))
    {
        enableRate(acsConfig.getSatOpts(observation.Sat));
        for (const auto& [frequency,signal] : observation.sigs)
            enableRate(acsConfig.getSatOpts(observation.Sat,{enum_to_string(signal.code)}));
        const int satellite = observation.Sat.prn - [&]()
        { for (const auto& system : systems) if (system.sys == observation.Sat.sys) return system.firstPrn; return 0; }();
        for (int sample = -10; sample <= 10; ++sample)
        {
            Peph ephemeris; ephemeris.Sat = observation.Sat; ephemeris.time.bigTime = firstEpoch+30*sample;
            ephemeris.index = 0; ephemeris.clkIndex = 0;
            ephemeris.pos = positions[satellite] + velocity*(30*sample);
            ephemeris.vel = velocity; ephemeris.clk = observation.satClk;
            nav.pephMap[observation.Sat.id()][ephemeris.time] = ephemeris;
        }
    }
    const auto originalObservations = receiver.obsList;
    KFState effective = clockState(receiver,systems,trace);
    KFState physical = effective;
    for (const auto& system : systems)
    {
        physical.x[physical.kfIndexMap.at(zhangReceiverClockKey(acsConfig,id,system.sys))] += system.origin;
        for (int satellite = 0; satellite < static_cast<int>(positions.size()); ++satellite)
        {
            KFKey key; key.type = KF::SAT_CLOCK; key.Sat = SatSys(system.sys,system.firstPrn+satellite);
            physical.x[physical.kfIndexMap.at(key)] += system.origin;
        }
    }
    const bool processCode = acsConfig.process_meas[CODE], processPhase = acsConfig.process_meas[PHAS];
    acsConfig.process_meas[PHAS] = false; acsConfig.process_meas[CODE] = true;
    KFState remote;
    KFMeasEntryList rawEntries, effectiveEntries;
    receiverUducGnss(trace,receiver,physical,rawEntries,remote);
    std::map<SatSys,Vector3d> physicalPositions;
    std::map<SatSys,double> physicalClocks;
    for (const auto& observation : only<GObs>(receiver.obsList))
    { physicalPositions[observation.Sat] = observation.rSatCom; physicalClocks[observation.Sat] = observation.satClk; }
    // Fresh copies also cover the receiverUducGnss satPosClk call rather than
    // relying on a valid preprocessor cache to happen to keep the same orbit.
    receiver.obsList.clear();
    for (const auto& observation : originalObservations)
        receiver.obsList.push_back(std::make_shared<GObs>(dynamic_cast<GObs&>(*observation)));
    receiverUducGnss(trace,receiver,effective,effectiveEntries,remote);
    require(physical.stateTransition(trace,tsync) && effective.stateTransition(trace,tsync),
        "new configured satellite rate states were not allocated");
    acsConfig.process_meas[CODE] = processCode; acsConfig.process_meas[PHAS] = processPhase;
    require(rawEntries.size() == systems.size()*positions.size()*2 && rawEntries.size() == effectiveEntries.size(),
        "moving orbit fixture did not execute all production code observations");
    double maxPosition = 0, maxClockCache = 0, maxPrediction = 0, maxDesign = 0, maxRateModel = 0;
    std::map<SatSys,double> clockEpochOffset;
    for (const auto& observation : only<GObs>(receiver.obsList))
    {
        maxPosition = std::max(maxPosition,(Vector3d(observation.rSatCom)-physicalPositions.at(observation.Sat)).norm());
        // obs.satClk is used to locate the transmission event. Its coordinate
        // is physical/external, not a freely shifted S-basis product clock.
        maxClockCache = std::max(maxClockCache,std::abs(observation.satClk-physicalClocks.at(observation.Sat))*CLIGHT);
        clockEpochOffset[observation.Sat] = (observation.time-tsync).to_double()-observation.tof;
    }
    MatrixXd rawH = MatrixXd::Zero(rawEntries.size(),physical.x.size());
    MatrixXd effectiveH = MatrixXd::Zero(effectiveEntries.size(),effective.x.size());
    for (int row = 0; row < static_cast<int>(rawEntries.size()); ++row)
    {
        require(rawEntries[row].obsKey == effectiveEntries[row].obsKey,"gauge change reordered code observations");
        auto prediction = [](const KFMeasEntry& entry)
        { double sum = 0; for (const auto& [component,details] : entry.componentsMap)
            if (component != E_Component::OBSERVED && component != E_Component::NET_RESIDUAL) sum += details.value;
          return sum; };
        maxPrediction = std::max(maxPrediction,std::abs(prediction(rawEntries[row])-prediction(effectiveEntries[row])));
        KFKey satellite; satellite.type = KF::SAT_CLOCK; satellite.Sat = rawEntries[row].obsKey.Sat;
        maxDesign = std::max(maxDesign,std::abs(effectiveEntries[row].designEntryMap.at(satellite)+1));
        KFKey rateKey = satellite; rateKey.type = KF::SAT_CLOCK_RATE;
        const double dt = clockEpochOffset.at(satellite.Sat);
        maxDesign = std::max(maxDesign,std::abs(effectiveEntries[row].designEntryMap.at(rateKey)+dt));
        const double expectedClock = effective.x[effective.kfIndexMap.at(satellite)] + configuredSatelliteRate*dt;
        maxRateModel = std::max(maxRateModel,std::abs(effectiveEntries[row].componentsMap.at(E_Component::SAT_CLOCK).value+expectedClock));
        for (const auto& [key,value] : rawEntries[row].designEntryMap) rawH(row,physical.kfIndexMap.at(key)) = value;
        for (const auto& [key,value] : effectiveEntries[row].designEntryMap) effectiveH(row,effective.kfIndexMap.at(key)) = value;
    }
    const double covarianceDifference = (rawH*physical.P*rawH.transpose()-effectiveH*effective.P*effectiveH.transpose()).norm();
    const bool correct = maxPosition < 1e-6 && maxClockCache < 1e-6 && maxPrediction < 1e-6 && maxDesign < 1e-12 &&
        maxRateModel < 1e-9 && covarianceDifference < 1e-7;
    std::cout << std::setprecision(17) << "CLOCK_LIFECYCLE moving_orbit_systems=" << systems.size()
        << " first_system=" << enum_to_string(systems.front().sys) << " physical_cache_position_difference_m=" << maxPosition
        << " physical_clock_cache_difference_m=" << maxClockCache
        << " full_prediction_difference_m=" << maxPrediction << " effective_clock_design_error=" << maxDesign
        << " new_configured_rate_model_error_m=" << maxRateModel
        << " full_predicted_covariance_difference=" << covarianceDifference << " equivalent=" << correct << '\n';
    if (!correct && observeOnly) return false;
    require(correct,"effective clock gauge changed physical transmit time/range/design/full observation covariance");
    return true;
}

bool missingExternalSource(Trace& trace, bool observeOnly)
{
    configure({allSystems.front()});
    const std::string id = "CLOCK_MISSING_EXTERNAL_SOURCE";
    auto& receiver = receiverMap[id]; receiver.id = id; receiver.ready = true;
    tsync.bigTime = firstEpoch;
    observations(receiver,{allSystems.front()},firstEpoch,0,true);
    KFState state = clockState(receiver,{allSystems.front()},trace);
    const KFState before = state;
    for (const auto& observation : only<GObs>(receiver.obsList))
    { nav.pclkMap.erase(observation.Sat.id()); nav.pephMap.erase(observation.Sat.id()); }
    spp(trace,receiver,&state,nullptr);
    const bool rejected = receiver.sol.status != E_Solution::SINGLE;
    require(state.kfIndexMap == before.kfIndexMap && (state.x-before.x).norm() == 0 && (state.P-before.P).norm() == 0,
        "failed independent SPP mutated authoritative FLOAT state");
    for (const auto& observation : only<GObs>(receiver.obsList))
        require(observation.clkSource == E_Source::KALMAN && observation.ephClkValid,
            "failed independent SPP destroyed the caller's saved ephemeris cache");
    std::cout << "CLOCK_LIFECYCLE physical_external_missing_with_valid_KF_cache rejected=" << rejected << '\n';
    if (!rejected && observeOnly) return false;
    require(rejected,"independent SPP reused valid effective KF clocks after external-source failure");
    // A valid old PPP cache must not bypass the physical-source contract even
    // when ordinary user QC satellite-health screening is explicitly off.
    observations(receiver,{allSystems.front()},firstEpoch,0,true);
    configureClockOnlyRows(receiver);
    for (const auto& observation : only<GObs>(receiver.obsList))
    { nav.pclkMap.erase(observation.Sat.id()); nav.pephMap.erase(observation.Sat.id()); }
    const bool oldSvh = acsConfig.exclude.svh;
    const bool oldCode = acsConfig.process_meas[CODE], oldPhase = acsConfig.process_meas[PHAS];
    acsConfig.exclude.svh = false;
    acsConfig.process_meas[CODE] = true; acsConfig.process_meas[PHAS] = false;
    KFState remote; KFMeasEntryList entries;
    receiverUducGnss(trace,receiver,state,entries,remote);
    acsConfig.exclude.svh = oldSvh;
    acsConfig.process_meas[CODE] = oldCode; acsConfig.process_meas[PHAS] = oldPhase;
    require(std::none_of(entries.begin(),entries.end(),[](const auto& entry) { return entry.valid; }),
        "missing physical ephemeris-time clock bypassed PPP through exclude.svh=false");
    std::cout << "CLOCK_LIFECYCLE physical_external_missing_PPP_svh_disabled valid_rows=0\n";
    return true;
}

void coldSeedLifecycle(Trace& trace, const std::vector<System>& systems)
{
    configure(systems);
    auto& pivot = receiverMap["PIVOT"]; pivot.id = "PIVOT"; pivot.ready = true;
    auto& user = receiverMap["COLD_CLOCK_USER"]; user.id = "COLD_CLOCK_USER"; user.ready = true;
    KFState state, remote;
    state.time.bigTime = firstEpoch;
    const bool processCode = acsConfig.process_meas[CODE], processPhase = acsConfig.process_meas[PHAS];
    acsConfig.process_meas[CODE] = true; acsConfig.process_meas[PHAS] = false;
    for (int epoch = 0; epoch < 2; ++epoch)
    {
        tsync.bigTime = firstEpoch+30*epoch;
        if (epoch) require(state.stateTransition(trace,tsync),"cold lifecycle prior transition failed");
        for (auto* receiver : {&pivot,&user})
        {
            observations(*receiver,systems,firstEpoch+30*epoch,0,false);
            if (receiver == &pivot)
                for (auto& observation : only<GObs>(receiver->obsList))
                    for (auto& [frequency,signal] : observation.sigs)
                    { signal.P -= receiverRelativeClock; signal.L = signal.P / observation.satNav_ptr->lamMap.at(frequency); }
            spp(trace,*receiver,&state,nullptr);
            require(receiver->sol.status == E_Solution::SINGLE,"actual cold seed reference/user SPP failed");
            receiver->aprioriClk = receiver->sol.sppClk;
            configureClockOnlyRows(*receiver);
        }
        const KFState before = state;
        prepareZhangEffectiveClockSeeds(trace,receiverMap,state,remote);
        require(state.kfIndexMap == before.kfIndexMap && (state.x-before.x).norm() == 0 && (state.P-before.P).norm() == 0,
            "seed preparation modified existing effective means/full covariance");
        require(state.procNoiseMap == before.procNoiseMap && state.gaussMarkovTauMap == before.gaussMarkovTauMap &&
                state.gaussMarkovMuMap == before.gaussMarkovMuMap,"seed preparation changed existing process model");
        for (const auto& system : systems)
        {
            double translation = 0;
            require(zhangEffectiveClockSeed(state,system.sys,translation),"current cold/warm seed evidence absent");
            require(std::abs(translation+system.origin) < 1e-4,"cold or existing-anchor chart translation incorrect");
        }
        KFMeasEntryList entries;
        receiverUducGnss(trace,pivot,state,entries,remote);
        receiverUducGnss(trace,user,state,entries,remote);
        require(entries.size() == 2*systems.size()*positions.size()*2,"cold model omitted configured code rows");
        for (const auto& entry : entries)
        {
            require(entry.valid,"cold clock model rejected a valid current seed");
            const auto found = std::find_if(systems.begin(),systems.end(),[&](const System& system)
                { return system.sys == entry.obsKey.Sat.sys; });
            const int satellite = entry.obsKey.Sat.prn-found->firstPrn;
            const double external = 10000+static_cast<int>(found->sys)*100+10*satellite;
            const double expectedReceiver = entry.obsKey.str == "PIVOT" ? 0 : receiverRelativeClock;
            require(std::abs(entry.componentsMap.at(E_Component::REC_CLOCK).value-expectedReceiver) < 1e-4,
                "cold reference/nonreference effective receiver clock is not in a common datum");
            require(std::abs(entry.componentsMap.at(E_Component::SAT_CLOCK).value+external-found->origin) < 1e-4,
                "cold satellite component did not use the seeded effective clock");
        }
        require(state.stateTransition(trace,tsync),"actual cold clock allocation/transition failed");
        for (const auto& system : systems)
        {
            require(!state.kfIndexMap.contains(zhangReceiverClockKey(acsConfig,"PIVOT",system.sys)),
                "cold reference zero was replaced by an invented random clock prior");
            const auto key = zhangReceiverClockKey(acsConfig,user.id,system.sys);
            require(std::abs(state.x[state.kfIndexMap.at(key)]-receiverRelativeClock) < 1e-4,
                "allocated cold receiver clock has wrong chart");
            require(std::abs(state.P(state.kfIndexMap.at(key),state.kfIndexMap.at(key))-1e6) < 1e-6,
                "cold deterministic translation changed configured effective-coordinate prior variance");
        }
        if (epoch)
        {
            require(state.kfIndexMap == before.kfIndexMap && (state.x-before.x).norm() < 1e-6 && (state.P-before.P).norm() < 1e-7,
                "second cold lifecycle epoch reseeded existing clocks or changed complete effective prior");
        }
        const auto core = captureZhangCheckpointKfCore(state);
        KFState restored; std::string reason;
        require(restoreZhangCheckpointKfCoreWithReceiverResolver(core,restored,
            [](const std::string& id) -> Receiver*
            { const auto found = receiverMap.find(id); return found == receiverMap.end() ? nullptr : &found->second; },
            &reason),"clock metadata KF-core checkpoint restore failed: "+reason);
        require(restored.metaDataMap == state.metaDataMap && (restored.x-state.x).norm() == 0 && (restored.P-state.P).norm() == 0,
            "existing KF-core schema lost clock metadata/means/full covariance");
    }
    // SEED evidence is for one epoch only, not a permanently cached datum.
    tsync.bigTime = firstEpoch+60;
    for (const auto& system : systems)
    {
        double translation = 0;
        require(!zhangEffectiveClockSeed(state,system.sys,translation),"stale seed from the preceding epoch was reused");
        require(zhangSppPhysicalClockChart(pivot.sol,system.sys).empty(),"stale reference SPP chart was reused");
    }
    // No current physical anchor: the helper may not invent a new SPP gauge
    // or copy a surviving system's translation into this missing system.
    for (const auto& system : systems)
        for (int satellite = 0; satellite < static_cast<int>(positions.size()); ++satellite)
        {
            const auto id = SatSys(system.sys,system.firstPrn+satellite).id();
            nav.pclkMap.erase(id); nav.pephMap.erase(id);
        }
    const KFState beforeMissing = state;
    prepareZhangEffectiveClockSeeds(trace,receiverMap,state,remote);
    for (const auto& system : systems)
    {
        double translation = 0;
        require(!zhangEffectiveClockSeed(state,system.sys,translation),"missing current anchor did not fail closed");
    }
    require((state.x-beforeMissing.x).norm() == 0 && (state.P-beforeMissing.P).norm() == 0,
        "missing-anchor handling damaged an existing posterior");
    acsConfig.process_meas[CODE] = processCode; acsConfig.process_meas[PHAS] = processPhase;
    std::cout << "CLOCK_LIFECYCLE actual_cold_seed_systems=" << systems.size()
        << " first_system=" << enum_to_string(systems.front().sys)
        << " epochs=2 complete_effective_prior_preserved=1 checkpoint_metadata_preserved=1 stale_and_missing_anchor_rejected=1\n";
}
}

int main(int argc, char** argv)
{
    std::ostringstream trace;
    const bool observeOnly = argc > 2 && std::string(argv[2]) == "--observe-only";
    try
    {
        int accepted = 0;
        for (const auto& system : allSystems)
            for (double jump : {0.0, 17.25}) accepted += lifecycle(trace, {system}, jump, observeOnly);
        for (double jump : {0.0, 17.25}) accepted += lifecycle(trace, allSystems, jump, observeOnly);
        const bool legacyAccepted = legacyCheckpointBoundary(trace,observeOnly);
        int movingAccepted = 0;
        for (const auto& system : allSystems) movingAccepted += movingOrbitObservationModel(trace,{system},observeOnly);
        movingAccepted += movingOrbitObservationModel(trace,allSystems,observeOnly);
        const bool missingRejected = missingExternalSource(trace,observeOnly);
        for (const auto& system : allSystems) coldSeedLifecycle(trace,{system});
        coldSeedLifecycle(trace,allSystems);
        if (argc > 1) { std::ofstream out(argv[1]); out << trace.str(); }
        std::cout << (accepted == 8 && movingAccepted == 4 && legacyAccepted && missingRejected ? "PASS" : "OBSERVED_FAILURE") << " multi_gnss_clock_lifecycle checks=" << checks
            << " cold_SPP_epochs=16 cases=" << accepted << "/8 actual_spp_and_offset_and_transition=1"
            << " moving_orbit_cases=" << movingAccepted << "/4"
            << " legacy_boundary=" << legacyAccepted
            << " missing_external_rejected=" << missingRejected
            << " synthetic_only=1 real_PPP_acceptance=0 observe_only=" << observeOnly << '\n';
        return (accepted == 8 && movingAccepted == 4 && legacyAccepted && missingRejected) || observeOnly ? 0 : 1;
    }
    catch (const std::exception& error)
    {
        if (argc > 1) { std::ofstream out(argv[1]); out << trace.str(); }
        std::cerr << "FAIL multi_gnss_clock_lifecycle " << error.what() << '\n';
        return 1;
    }
}
