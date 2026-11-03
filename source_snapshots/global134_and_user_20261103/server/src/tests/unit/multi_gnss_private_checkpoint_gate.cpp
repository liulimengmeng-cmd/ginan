// Production-linked adapter contract. No producer source is included and no
// private map is exposed. Public operations create the controlled histories.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangMultiGnssFloatCheckpoint.hpp"
#include "common/zhangProductRelationAdmission.hpp"
#include "common/zhangProductDatumTransition.hpp"
#include "pea/zhangPppAr.hpp"
#include "pea/ppp.hpp"
#include "pea/zhangReference.hpp"
#include "ambres/GNSSambres.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/json.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "private_checkpoint_wire_schema.hpp"
#include "private_ambres_wire_schema.hpp"
#include "private_graph_wire_schema.hpp"

namespace {
unsigned checks = 0;
void require(bool pass, const std::string& reason)
{ ++checks; if (!pass) throw std::runtime_error(reason); }
const std::string runtime = "MULTI_FLOAT_PRIVATE_GATE_CONTRACT";

void mode(bool multi)
{
    ::setenv("ZHANG_R51_ENABLE", "1", 1);
    ::setenv("ZHANG_R51_RATIO_ONLY", "1", 1);
    ::setenv("ZHANG_R51_EPOCH_ONLY_AR", "1", 1);
    ::setenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", multi ? "1" : "0", 1);
    acsConfig.process_sys.clear();
    acsConfig.process_sys[E_Sys::GPS] = true;
    if (multi) {
        acsConfig.process_sys[E_Sys::GAL] = true;
        acsConfig.process_sys[E_Sys::BDS] = true;
    }
    acsConfig.zhangFullRank.enable = true;
    acsConfig.zhangPppAr.user_adapter = false;
    acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow = false;
    acsConfig.zhangPppAr.targeted_besd_capture_shadow = false;
    acsConfig.zhangPppAr.baseline_observables = {
        {E_Sys::GPS, {E_ObsCode::L1C, E_ObsCode::L2W}},
        {E_Sys::GAL, {E_ObsCode::L1C, E_ObsCode::L5Q}},
        {E_Sys::BDS, {E_ObsCode::L1P, E_ObsCode::L5P}}};
    require(zhangMultiGnssFloatCheckpointActive(acsConfig) == multi, "MODE_SELECTION_FAILED");
}

// Independent bounded prefix reader, used ONLY to locate exact scalar wire
// offsets in an authenticated production-exported EMPTY schema-6 payload.
// It neither constructs nor validates the production envelope or gate.
struct EmptyPending {
    template<class Archive> void serialize(Archive&, unsigned)
    { throw std::runtime_error("PREFIX_READER_REQUIRES_EMPTY_PENDING_VECTOR"); }
};
struct PinsPrefix {
    ZhangCheckpointTime time;
    std::set<std::string> identities;
    std::istream* stream = nullptr;
    std::streamoff timeEnd = -1;
    template<class Archive> void serialize(Archive& ar, unsigned) {
        ar & time;
        timeEnd = stream->tellg();
        ar & identities;
    }
};
struct PppPrefix {
    std::uint32_t schema = 0;
    std::string runtimeId;
    std::vector<EmptyPending> pending;
    bool hasPins = false;
    PinsPrefix pins;
    bool configured = false, hasCapture = false;
    std::istream* stream = nullptr;
    std::streamoff hasPinsOffset = -1, configuredOffset = -1, hasCaptureOffset = -1;
    template<class Archive> void serialize(Archive& ar, unsigned) {
        ar & schema & runtimeId & pending;
        hasPinsOffset = stream->tellg(); ar & hasPins;
        pins.stream = stream; ar & pins;
        configuredOffset = stream->tellg(); ar & configured;
        hasCaptureOffset = stream->tellg(); ar & hasCapture;
    }
};
PppPrefix prefix(const std::string& bytes)
{
    std::istringstream input(bytes, std::ios::binary | std::ios::in);
    boost::archive::binary_iarchive archive(input, boost::archive::no_header);
    PppPrefix value; value.stream = &input; archive >> value;
    require(value.schema == ZHANG_PPP_AR_CHECKPOINT_SCHEMA_VERSION &&
        value.runtimeId == runtime && value.pending.empty() && !value.hasPins &&
        value.pins.identities.empty() && !value.hasCapture, "PREFIX_NOT_EMPTY_SCHEMA6");
    return value;
}
std::string scalarByte(const std::string& input, std::streamoff offset)
{
    require(offset >= 0 && static_cast<std::size_t>(offset) < input.size(), "PREFIX_OFFSET_OUT_OF_RANGE");
    require(input.at(offset) == 0, "EXPECTED_ZERO_WIRE_SCALAR");
    std::string changed = input; changed.at(offset) = 1;
    return changed;
}
void rejected(const ZhangPppArCheckpointResult& result, const std::string& field)
{
    require(!result.valid, "PRIVATE_PPP_AR_IMPORT_ACCEPTED:" + field);
    require(result.failureReason == "PPP_AR_MULTI_GNSS_FLOAT_PRIVATE_STATE_NONEMPTY:" + field,
        "WRONG_PPP_AR_REJECTION:" + result.failureReason);
}
void rejected(const ZhangAmbresCheckpointResult& result, const std::string& field)
{
    require(!result.valid, "PRIVATE_AMBRES_ACCEPTED:" + field);
    require(result.failureReason == "AMBRES_MULTI_GNSS_FLOAT_PRIVATE_STATE_NONEMPTY:" + field,
        "WRONG_AMBRES_REJECTION:" + result.failureReason);
}

using FloatWire = checkpoint_test_wire::ZhangPppArCheckpointEnvelope;
FloatWire decode(const std::string& bytes)
{
    FloatWire value;
    std::istringstream input(bytes, std::ios::binary | std::ios::in);
    boost::archive::binary_iarchive archive(input, boost::archive::no_header);
    archive >> value;
    require(input.peek() == std::char_traits<char>::eof(), "TEST_DTO_TRAILING_BYTES");
    return value;
}
template<class Wire> std::string encode(const Wire& value)
{
    std::ostringstream output(std::ios::binary | std::ios::out);
    boost::archive::binary_oarchive archive(output, boost::archive::no_header);
    archive << value;
    return output.str();
}

using AmbresWire = checkpoint_test_ambres_wire::AmbresCheckpointEnvelope;
AmbresWire decodeAmbres(const std::string& bytes)
{
    AmbresWire value;
    std::istringstream input(bytes, std::ios::binary | std::ios::in);
    boost::archive::binary_iarchive archive(input, boost::archive::no_header);
    archive >> value;
    require(input.peek() == std::char_traits<char>::eof(), "AMBRES_TEST_DTO_TRAILING_BYTES");
    require(encode(value) == bytes, "AMBRES_SCHEMA_MIRROR_NOT_BYTE_EXACT");
    return value;
}

void disabledPriorCapture(KFState& owner, const std::string& emptyAmbres)
{
    mode(true);
    acsConfig.zhangFullRank.enable = true;
    acsConfig.zhangPppAr.multi_epoch_relink_shadow = false;
    acsConfig.zhangPppAr.whitened_wl_fixed_lag_shadow = false;
    std::string payload;
    for (int epoch = 0; epoch < 2; ++epoch)
    {
        captureZhangPppArFloatPrior(owner);
        const auto result = exportZhangAmbresCheckpointSection(owner, runtime, payload);
        require(result.valid && payload == emptyAmbres,
            "DISABLED_PRIOR_CAPTURE_REGISTERED_EMPTY_RUNTIME:" + result.failureReason);
        require(importZhangAmbresCheckpointSection(owner, runtime, payload, true).valid,
            "DISABLED_PRIOR_CAPTURE_PREFLIGHT_FAILED");
    }
    acsConfig.zhangFullRank.enable = false;
    acsConfig.zhangPppAr.multi_epoch_relink_shadow = true;
    captureZhangPppArFloatPrior(owner);
    require(exportZhangAmbresCheckpointSection(owner, runtime, payload).valid && payload == emptyAmbres,
        "DISABLED_FULL_RANK_CAPTURE_REGISTERED_RUNTIME");
    acsConfig.zhangFullRank.enable = true;

    // Enabled mode with no usable graph still registers its runtime, exactly
    // as before. Turning it off must clear the inner cache, not erase the
    // outer object merely to make the new FLOAT checkpoint appear empty.
    captureZhangPppArFloatPrior(owner);
    rejected(exportZhangAmbresCheckpointSection(owner, runtime, payload), "live_zhangRelinkPriorMoments");
    acsConfig.zhangPppAr.multi_epoch_relink_shadow = false;
    captureZhangPppArFloatPrior(owner);
    captureZhangPppArFloatPrior(owner);
    rejected(exportZhangAmbresCheckpointSection(owner, runtime, payload), "live_zhangRelinkPriorMoments");
    mode(false);
    require(importZhangAmbresCheckpointSection(owner, runtime, emptyAmbres, false).valid,
        "EMPTY_PRIOR_REGISTRATION_CONTROL_CLEANUP_FAILED");

    // A controlled public graph import supplies two receivers, two satellites
    // and different valid S/product spanning trees. The one current chord has
    // L1=3.25, L2=1.5, P11=4, P22=9, P12=.5. An independently specified
    // marginal oracle is |L1|=3.25, |WL|=1.75, Var(L1)=4, Var(WL)=12,
    // Cov(L1,WL)=3.5, in each ordered satellite-pair direction.
    namespace gw = checkpoint_test_graph_wire;
    const auto savedOwner = owner;
    const auto savedProductBasisMode = acsConfig.zhangFullRank.product_basis_mode;
    acsConfig.zhangFullRank.product_basis_mode = "LEGACY_PRODUCT_TREE";
    std::string emptyGraph, reason;
    require(exportZhangGraphCheckpointSection(owner, runtime, emptyGraph, reason), reason);
    gw::ZhangGraphCheckpointPayloadV4 graph;
    graph.runtimeId = runtime;
    gw::ZhangGraphCheckpointRuntimeV4 state;
    state.common.system = static_cast<int>(E_Sys::GPS);
    auto edge = [](const char* rec, int prn) {
        return gw::ZhangGraphCheckpointEdge{rec, {static_cast<int>(E_Sys::GPS), prn}};
    };
    auto& basis = state.common.basis;
    basis.rootReceiver = "R0";
    basis.edges = {edge("R0", 13), edge("R0", 14), edge("R1", 13), edge("R1", 14)};
    basis.treeEdges = {edge("R0", 13), edge("R0", 14), edge("R1", 13)};
    basis.receivers = {"R0", "R1"};
    basis.satellites = {{static_cast<int>(E_Sys::GPS), 13}, {static_cast<int>(E_Sys::GPS), 14}};
    basis.componentCount = 1; basis.connected = true;
    state.common.activeBasis = basis;
    state.common.productBasis = basis;
    state.common.productBasis.treeEdges = {edge("R0", 13), edge("R1", 13), edge("R1", 14)};
    state.common.initialized = true;
    state.common.productInitialized = true;
    state.common.productAuthorityComponentVersion = 0;
    state.common.stateEdges = basis.edges;
    state.signalIdentityAuthoritative = true;
    for (const auto& key : basis.edges)
    {
        state.common.edgeHistory[key] = {};
        state.common.phaseArcVersions[key] = {
            {static_cast<int>(E_ObsCode::L1C), 0}, {static_cast<int>(E_ObsCode::L2W), 0}};
    }
    for (const auto& key : state.common.productBasis.treeEdges)
    {
        state.common.productArcVersions[key] = 0;
        state.common.arcVersionObservationHistory[{key, 0}] = {};
        state.productSignalArcVersions[key] = state.common.phaseArcVersions.at(key);
    }
    graph.graphStates.push_back(state);
    const auto graphBytes = std::string("GINAN_ZHANG_GRAPH_V4\n") + encode(graph);
    // This fixture deliberately has independent S/P trees. Native SINGLE_S
    // must still reject its V4 wire; only this legacy-control scope permits it.
    acsConfig.zhangFullRank.product_basis_mode = "SINGLE_S_BASIS_NATIVE";
    const bool nativeLegacyGraphValid = validateZhangGraphCheckpointSection(runtime, graphBytes, reason);
    require(!nativeLegacyGraphValid && reason ==
        "ZHANG_SINGLE_S_CHECKPOINT_REJECTED_LEGACY_V3_V4_REQUIRES_EXPLICIT_MIGRATION",
        "NATIVE_ACCEPTED_LEGACY_INDEPENDENT_PRODUCT_TREE:" + reason);
    acsConfig.zhangFullRank.product_basis_mode = "LEGACY_PRODUCT_TREE";
    const bool controlledGraphValid = validateZhangGraphCheckpointSection(runtime, graphBytes, reason);
    require(controlledGraphValid, "CONTROLLED_PRIOR_GRAPH_VALIDATION:" + reason);
    require(importZhangGraphCheckpointSection(owner, runtime, graphBytes, reason), reason);
    ZhangGraphIntegerContext context;
    require(zhangGraphIntegerContext(owner, E_Sys::GPS, context) && context.productBasis.connected,
        "CONTROLLED_PRIOR_GRAPH_LOST_AUTHORITATIVE_PRODUCT_CHART");
    const auto target = ZhangProductTargetBuilder::build(context.basis, context.productBasis);
    require(target.valid && target.currentChords.size() == 1 &&
        target.targetSatellites.size() == 1 && target.matrix.size() == 1 &&
        target.matrix.front().size() == 1 && target.matrix.front().front() != 0,
        "CONTROLLED_PRIOR_GRAPH_REQUIRES_NONZERO_ONE_CYCLE_PRODUCT_TARGET");
    owner.kfIndexMap.clear();
    owner.kfIndexMap[KFKey{KF::AMBIGUITY, SatSys(E_Sys::GPS, 14), "R1", static_cast<int>(E_ObsCode::L1C)}] = 0;
    owner.kfIndexMap[KFKey{KF::AMBIGUITY, SatSys(E_Sys::GPS, 14), "R1", static_cast<int>(E_ObsCode::L2W)}] = 1;
    owner.x << 3.25, 1.5;
    owner.P << 4, .5, .5, 9;
    acsConfig.zhangFullRank.sysOpts.clear();
    acsConfig.zhangFullRank.sysOpts[E_Sys::GPS].baseline_observables = {E_ObsCode::L1C, E_ObsCode::L2W};
    acsConfig.zhangFullRank.sysOpts[E_Sys::GPS].use_spanning_tree = false;
    acsConfig.zhangPppAr.multi_epoch_relink_shadow = true;
    const auto captureAndCheck = [&]() {
        captureZhangPppArFloatPrior(owner);
        const auto result = exportZhangAmbresCheckpointSection(owner, runtime, payload);
        require(result.valid && result.relinkPriorMoments == 2,
            "ENABLED_PRIOR_CAPTURE_LOST_MOMENTS:count=" + std::to_string(result.relinkPriorMoments) +
                ":" + result.failureReason);
        const auto schema = decodeAmbres(payload);
        require(schema.relinkPriors.size() == 2, "ENABLED_CAPTURE_DTO_MOMENTS_MISSING");
        for (const auto& moment : schema.relinkPriors)
            require(std::abs(std::abs(moment.firstMean) - 3.25) < 1e-12 &&
                std::abs(std::abs(moment.wideLaneMean) - 1.75) < 1e-12 &&
                std::abs(moment.firstVariance - 4) < 1e-12 &&
                std::abs(moment.wideLaneVariance - 12) < 1e-12 &&
                std::abs(moment.firstWideLaneCovariance - 3.5) < 1e-12,
                "ENABLED_CAPTURE_SCALAR_MARGINAL_ORACLE_CHANGED");
        return payload;
    };
    const auto enabledPayload = captureAndCheck();
    require(captureAndCheck() == enabledPayload, "ENABLED_CAPTURE_NOT_DETERMINISTIC");
    mode(true);
    rejected(exportZhangAmbresCheckpointSection(owner, runtime, payload), "live_zhangRelinkPriorMoments");
    acsConfig.zhangPppAr.multi_epoch_relink_shadow = false;
    acsConfig.zhangPppAr.whitened_wl_fixed_lag_shadow = false;
    captureZhangPppArFloatPrior(owner);
    rejected(exportZhangAmbresCheckpointSection(owner, runtime, payload), "live_zhangRelinkPriorMoments");
    mode(false);
    auto cleared = exportZhangAmbresCheckpointSection(owner, runtime, payload);
    require(cleared.valid && cleared.relinkPriorMoments == 0 && decodeAmbres(payload).relinkPriors.empty(),
        "DISABLED_EXISTING_CAPTURE_DID_NOT_CLEAR_INNER");
    require(importZhangAmbresCheckpointSection(owner, runtime, emptyAmbres, false).valid,
        "PRIOR_FIXTURE_CLEANUP_FAILED");
    mode(true);
    rejected(importZhangAmbresCheckpointSection(owner, runtime, enabledPayload, true), "relinkPriors");
    require(importZhangGraphCheckpointSection(owner, runtime, emptyGraph, reason), reason);
    acsConfig.zhangFullRank.product_basis_mode = savedProductBasisMode;
    owner = savedOwner;
    require(exportZhangAmbresCheckpointSection(owner, runtime, payload).valid && payload == emptyAmbres,
        "PRIOR_CAPTURE_TEST_LEFT_HISTORY");
}

// This tests the real production adapters and actual graph-callback public
// entry points. The DTO mirror only encodes independent negative fixtures;
// it does not contain the whitelist or call a copied reference classifier.
void floatBookkeeping(KFState& owner, const std::string& empty)
{
    mode(true);
    acsConfig.zhangPppAr.output_products = true;
    acsConfig.zhangPppAr.product_mode = "HOU_OSB_LIKE";
    acsConfig.zhangPppAr.hou_product_coordinate = "HYBRID_STABLE";
    ::setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF", "1", 1);
    const auto systems = {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS};
    for (const auto sys : systems)
    {
        const SatSys first(sys, 13), second(sys, 14);
        const auto& codes = acsConfig.zhangPppAr.baseline_observables.at(sys);
        for (const auto code : codes)
        {
            (void)zhangSatelliteDatumStatus(sys, code, first);
            (void)zhangSatelliteDatumStatus(sys, code, second);
            recordZhangExactPhaseTransforms(owner.time, sys, code, {{first, .125}, {second, -.625}});
        }
        recordZhangSatellitePhaseDiscontinuity(owner.time, sys, codes, first,
            "CONTROLLED_FLOAT_PHYSICAL_DISCONTINUITY");
        recordZhangPhaseReinitialisation(owner.time, sys, codes,
            "CONTROLLED_FLOAT_GRAPH_REINITIALISATION", {second});
    }
    std::string natural;
    auto exported = exportZhangPppArCheckpointSection(owner, runtime, natural);
    require(exported.valid, "FLOAT_GRAPH_CALLBACK_EXPORT:" + exported.failureReason);
    auto mirror = decode(natural);
    require(encode(mirror) == natural, "FULL_DTO_MIRROR_NOT_BYTE_EXACT");
    require(mirror.continuity.size() == 12 && mirror.globalContinuity.size() == 6 &&
        mirror.satelliteDatumManagers.size() == 6, "FLOAT_CALLBACK_TABLES_MISSING");
    for (const auto& manager : mirror.satelliteDatumManagers)
    {
        require(manager.currentSegments.at(SatSys(manager.system, 13)) == 1 &&
            manager.discontinuityCounters.at(SatSys(manager.system, 13)) == 1 &&
            manager.datumVersions.at(SatSys(manager.system, 13)) == 1 &&
            manager.topologyVersion == 1 && manager.alignmentGeneration >= 1,
            "PHYSICAL_EVENT_NOT_PRESERVED");
        require(!manager.alignmentKnown.empty() && manager.precisionValid.empty() &&
            manager.relations.empty() && manager.eventCounts.empty(),
            "SINGLETON_BIRTH_IS_NOT_AN_INTEGER_CERTIFICATE");
    }
    require(importZhangPppArCheckpointSection(owner, runtime, empty, false, 2).valid,
        "FLOAT_TEST_RESET_FAILED");
    require(importZhangPppArCheckpointSection(owner, runtime, natural, true, 2).valid,
        "FLOAT_NATURAL_PREFLIGHT_FAILED");
    require(importZhangPppArCheckpointSection(owner, runtime, natural, false, 2).valid,
        "FLOAT_NATURAL_COMMIT_FAILED");
    std::string after;
    require(exportZhangPppArCheckpointSection(owner, runtime, after).valid && after == natural,
        "FLOAT_NATURAL_ROUNDTRIP_DROPPED_METADATA");
    for (const auto sys : systems)
    for (const auto code : acsConfig.zhangPppAr.baseline_observables.at(sys))
    {
        const auto state = zhangSatelliteDatumStatus(sys, code, SatSys(sys, 13));
        require(state.phaseSegment == 1 && state.discontinuityCounter == 1 &&
            state.datumVersion == 1 && state.componentSize == 1 && state.componentRank == 0 &&
            !state.integerDatumContinuous && !state.integerPrecisionValid && !state.integerValid,
            "FLOAT_RESTORE_MANUFACTURED_INTEGER_AUTHORITY");
        long long difference = 0;
        require(!queryZhangSatelliteProductRelation(sys, code, SatSys(sys, 13), SatSys(sys, 14), difference),
            "FLOAT_RESTORE_CONNECTED_UNRELATED_SATELLITES");
    }
    // All allowed scalar metadata, including nonzero affine chart offsets,
    // survives. Such offsets alone are not integer precision evidence.
    auto shifted = mirror;
    auto& continuity = shifted.continuity.begin()->second;
    continuity.integerShiftCycles = 17;
    continuity.fractionalShiftCycles = -.375;
    continuity.counter = 3;
    continuity.datumVersion = 4;
    continuity.iod = 5;
    continuity.stabilizationRemaining = 2;
    continuity.validFrom = captureZhangCheckpointTime(owner.time);
    continuity.lastEpoch = captureZhangCheckpointTime(owner.time);
    auto& global = shifted.globalContinuity.begin()->second;
    global.counter = 2; global.datumVersion = 3; global.iod = 4;
    global.stabilizationRemaining = 1; global.validFrom = captureZhangCheckpointTime(owner.time);
    shifted.satelliteDatumManagers.front().alignmentCycles.begin()->second = -11;
    const auto shiftedWire = encode(shifted);
    require(importZhangPppArCheckpointSection(owner, runtime, shiftedWire, false, 2).valid,
        "FLOAT_AFFINE_METADATA_REJECTED");
    require(exportZhangPppArCheckpointSection(owner, runtime, after).valid && after == shiftedWire,
        "FLOAT_AFFINE_METADATA_DROPPED");

    auto negative = [&](const std::string& field, auto mutate)
    {
        auto polluted = shifted;
        mutate(polluted);
        const auto bytes = encode(polluted);
        rejected(importZhangPppArCheckpointSection(owner, runtime, bytes, true, 2), field);
        rejected(importZhangPppArCheckpointSection(owner, runtime, bytes, false, 2), field);
        require(exportZhangPppArCheckpointSection(owner, runtime, after).valid && after == shiftedWire,
            "REJECTED_FLOAT_FIXTURE_MUTATED_RUNTIME:" + field);
    };
    negative("continuity_hasFixedDatum", [](auto& value) { value.continuity.begin()->second.hasFixedDatum = true; });
    negative("continuity_signal", [](auto& value) {
        auto record = *value.continuity.begin();
        value.continuity.erase(value.continuity.begin());
        value.continuity[{record.first.satellite, E_ObsCode::L1X}] = record.second;
    });
    negative("continuity_invalid_metadata", [](auto& value) {
        GTime bad; bad.bigTime = std::numeric_limits<long double>::quiet_NaN();
        value.continuity.begin()->second.validFrom = captureZhangCheckpointTime(bad);
    });
    negative("globalContinuity_invalid_metadata", [](auto& value) {
        GTime bad; bad.bigTime = std::numeric_limits<long double>::infinity();
        value.globalContinuity.begin()->second.validFrom = captureZhangCheckpointTime(bad);
    });
    negative("satelliteDatumManagers_precisionValid", [](auto& value) {
        auto& manager = value.satelliteDatumManagers.front();
        manager.precisionValid.insert(manager.alignmentCycles.begin()->first);
    });
    negative("satelliteDatumManagers_relations", [](auto& value) {
        auto& manager = value.satelliteDatumManagers.front();
        manager.relations.push_back({{SatSys(manager.system, 13), 1}, {SatSys(manager.system, 14), 0}, 7, true, "POLLUTED_INTEGER_FACT"});
    });
    negative("satelliteDatumManagers_redundantRelations", [](auto& value) {
        auto& manager = value.satelliteDatumManagers.front();
        manager.redundantRelations.push_back({{SatSys(manager.system, 13), 1}, {SatSys(manager.system, 14), 0}, 7, true, "POLLUTED_REDUNDANT_FACT"});
    });
    negative("satelliteDatumManagers_eventCounts", [](auto& value) {
        value.satelliteDatumManagers.front().eventCounts[ZhangProductRelationEventType::NEW_COMPONENT_EDGE] = 0;
    });
    negative("satelliteDatumManagers_conflictCount", [](auto& value) {
        value.satelliteDatumManagers.front().conflictCount = 1;
    });
    negative("satelliteDatumManagers_invalid_metadata", [](auto& value) {
        value.satelliteDatumManagers.front().currentSegments.begin()->second = -1;
    });
    negative("productHistory", [](auto& value) {
        const auto key = value.continuity.begin()->first;
        value.productHistory[{"FLOAT", key.satellite, key.observable}] = {};
    });
    negative("productDatumVersionTrackers", [](auto& value) {
        value.productDatumVersionTrackers[E_Sys::GPS] = {};
    });
    // Publicly establish real integer authority in a legacy adapter control,
    // then ensure both live export and serialized multi import reject it.
    require(importZhangPppArCheckpointSection(owner, runtime, empty, false, 2).valid,
        "FLOAT_NEGATIVE_RESET_FAILED");
    mode(false);
    acsConfig.zhangPppAr.promotion_confirmation_epochs = 1;
    const auto promoted = promoteZhangSatelliteProductRelationDetailed(owner.time,
        E_Sys::GPS, E_ObsCode::L1C, SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 14), 9,
        "LEGACY_CONTROLLED_INTEGER_RELATION");
    require(promoted.accepted, "LEGACY_RELATION_CONTROL_NOT_ACCEPTED");
    std::string integerAuthority;
    exported = exportZhangPppArCheckpointSection(owner, runtime, integerAuthority);
    require(exported.valid, "LEGACY_RELATION_EXPORT_FAILED:" + exported.failureReason);
    mode(true);
    rejected(exportZhangPppArCheckpointSection(owner, runtime, after), "satelliteDatumManagers_precisionValid");
    mode(false);
    require(importZhangPppArCheckpointSection(owner, runtime, empty, false, 2).valid,
        "LEGACY_RELATION_TEST_CLEANUP_FAILED");
    mode(true);
    rejected(importZhangPppArCheckpointSection(owner, runtime, integerAuthority, true, 2),
        "satelliteDatumManagers_precisionValid");
}
}

static void numericalCheckpointRoundtrip(KFState& owner,const std::string& empty)
{
    mode(false);
    auto snapshot=decode(empty);
    const checkpoint_test_wire::ProductHistoryKey key{"PRODUCT_FIXED",SatSys("G13"),E_ObsCode::L1C};
    checkpoint_test_wire::ProductHistoryCheckpoint history;
    history.time=captureZhangCheckpointTime(owner.time);history.correction=100;
    snapshot.productHistory[key]=history;
    snapshot.productPublicationIdentity[key]={"PUBLISHED_SEGMENT","VERIFIED_PHYSICAL_CHART",1,0.01};
    const std::string domain="VERIFIED_PHYSICAL_CHART|G1|C0|D0";
    auto dynamics=zhangAdvanceProductNumericalDynamics({},99.75,30,120,domain);
    dynamics=zhangAdvanceProductNumericalDynamics(dynamics,100,30,120,domain);
    ZhangProductDatumTransition attestation;
    auto& endpoint=attestation.target;
    endpoint.epoch=double(owner.time.bigTime);endpoint.counter=endpoint.datum=0;
    endpoint.componentVersion=endpoint.alignment=1;endpoint.satellite="G13";endpoint.observable="L1C";
    endpoint.reference="G01";endpoint.segment="PUBLISHED_SEGMENT";endpoint.physicalFunctional="VERIFIED_PHYSICAL_CHART";
    endpoint.component="GPS_COMPONENT";endpoint.posterior="FROZEN_ROOT";endpoint.bundle="ACTUAL_BUNDLE";
    snapshot.productNumericalHistory[key]={dynamics,zhangDatumWireProof(attestation)};
    const auto bytes=encode(snapshot);
    auto result=importZhangPppArCheckpointSection(owner,runtime,bytes,true,2);
    require(result.valid,"PREDICTOR_NEW_SCHEMA_PREFLIGHT:"+result.failureReason);
    require(importZhangPppArCheckpointSection(owner,runtime,bytes,false,2).valid,"PREDICTOR_NEW_SCHEMA_COMMIT");
    std::string after;
    require(exportZhangPppArCheckpointSection(owner,runtime,after).valid && after==bytes,
        "PREDICTOR_OR_PUBLIC_ATTESTATION_LOST_NEW_SCHEMA_ROUNDTRIP");
    std::string legacyValidation;
    require(!exportZhangPppArCheckpointSection(owner,runtime,legacyValidation,6).valid && legacyValidation.empty(),
        "LEGACY_VALIDATION_EXPORT_DROPPED_LIVE_PREDICTOR_OR_ATTESTATION");
    double beforePrediction=0,afterPrediction=0;
    const auto restored=decode(after);
    require(zhangPredictProductNumericalChange(dynamics,100,90,120,domain,beforePrediction) &&
        zhangPredictProductNumericalChange(restored.productNumericalHistory.at(key).dynamics,100,90,120,domain,afterPrediction) &&
        beforePrediction==afterPrediction,"RESUMED_PREDICTION_DIFFERS_FROM_UNINTERRUPTED_PREDICTION");
    for(int bad=0;bad<5;++bad)
    {
        auto corrupted=snapshot;
        if(bad==0)corrupted.productNumericalHistory[key].dynamics.rateMetresPerSecond=std::numeric_limits<double>::infinity();
        if(bad==1)corrupted.productNumericalHistory[key].dynamics.domain="UNPROVED_DIFFERENT_CHART";
        if(bad==2)corrupted.productNumericalHistory[key].dynamics.correction=99;
        if(bad==3){auto wrong=attestation;++wrong.target.epoch;corrupted.productNumericalHistory[key].publicationProof=zhangDatumWireProof(wrong);}
        if(bad==4)corrupted.productNumericalHistory.clear();
        require(!importZhangPppArCheckpointSection(owner,runtime,encode(corrupted),true,2).valid &&
            !importZhangPppArCheckpointSection(owner,runtime,encode(corrupted),false,2).valid,
            "PREDICTOR_OR_ATTESTATION_CONTAMINATION_ACCEPTED:"+std::to_string(bad));
        require(exportZhangPppArCheckpointSection(owner,runtime,after).valid && after==bytes,
            "REJECTED_PREDICTOR_IMPORT_MUTATED_RUNTIME");
    }
    auto legacy=snapshot;legacy.schemaVersion=6;legacy.productNumericalHistory.clear();
    require(importZhangPppArCheckpointSection(owner,runtime,encode(legacy),false,2).valid,
        "LEGACY_SCHEMA6_PRODUCT_HISTORY_NO_LONGER_READABLE");
    require(exportZhangPppArCheckpointSection(owner,runtime,after).valid,"LEGACY_COLD_EXPORT_FAILED");
    const auto cold=decode(after).productNumericalHistory.at(key);
    require(!cold.dynamics.hasSample && !cold.dynamics.hasRate && cold.publicationProof.empty(),
        "LEGACY_SCHEMA6_INVENTED_PREDICTOR_OR_PUBLICATION_ATTESTATION");
    require(exportZhangPppArCheckpointSection(owner,runtime,legacyValidation,6).valid && legacyValidation==encode(legacy),
        "LEGACY_SOURCE_WIRE_FORMAT_NOT_EXACT_FOR_AFTER_RESTORE_VALIDATION");
    require(importZhangPppArCheckpointSection(owner,runtime,empty,false,2).valid,"PREDICTOR_TEST_RUNTIME_CLEANUP_FAILED");
    mode(true);
}

int main()
try {
    mode(true);
    KFState owner;
    owner.time.bigTime = 1405216800;
    owner.x = VectorXd(2); owner.x << 12, 34;
    owner.P = MatrixXd(2, 2); owner.P << 2, .25, .25, 3;
    owner.dx = VectorXd::Zero(2);
    owner.kfIndexMap.clear();
    owner.kfIndexMap[KFKey{KF::REC_POS, {}, "R0", 0}] = 0;
    owner.kfIndexMap[KFKey{KF::REC_POS, {}, "R0", 1}] = 1;
    std::string reason;
    require(bindZhangCheckpointRuntimeId(owner, runtime, &reason), reason);
    const auto originalX = owner.x;
    const auto originalP = owner.P;
    const auto originalMetadata = owner.metaDataMap;
    std::string emptyPpp, emptyAmbres;
    auto ppp = exportZhangPppArCheckpointSection(owner, runtime, emptyPpp);
    require(ppp.valid && !emptyPpp.empty(), "EMPTY_PPP_EXPORT:" + ppp.failureReason);
    auto ambres = exportZhangAmbresCheckpointSection(owner, runtime, emptyAmbres);
    require(ambres.valid && !emptyAmbres.empty(), "EMPTY_AMBRES_EXPORT:" + ambres.failureReason);
    require(importZhangPppArCheckpointSection(owner, runtime, emptyPpp, true, 2).valid,
        "EMPTY_PPP_PREFLIGHT_FAILED");
    require(importZhangAmbresCheckpointSection(owner, runtime, emptyAmbres, true).valid,
        "EMPTY_AMBRES_PREFLIGHT_FAILED");
    require(importZhangPppArCheckpointSection(owner, runtime, emptyPpp, false, 2).valid,
        "EMPTY_PPP_COMMIT_FAILED");
    require(importZhangAmbresCheckpointSection(owner, runtime, emptyAmbres, false).valid,
        "EMPTY_AMBRES_COMMIT_FAILED");

    const auto offsets = prefix(emptyPpp);
    const auto hasPins = scalarByte(emptyPpp, offsets.hasPinsOffset);
    rejected(importZhangPppArCheckpointSection(owner, runtime, hasPins, true, 2), "hasPendingSnapshotPins");
    const auto hasCapture = scalarByte(emptyPpp, offsets.hasCaptureOffset);
    rejected(importZhangPppArCheckpointSection(owner, runtime, hasCapture, true, 2), "hasFactorCapture");
    const auto hiddenPinTime = scalarByte(emptyPpp, offsets.pins.timeEnd - sizeof(long double));
    rejected(importZhangPppArCheckpointSection(owner, runtime, hiddenPinTime, true, 2),
        "inactive_payload_or_history");
    const auto callbackOnly = scalarByte(emptyPpp, offsets.configuredOffset);
    const auto callbackCheck = importZhangPppArCheckpointSection(owner, runtime, callbackOnly, true, 2);
    require(callbackCheck.valid, "CALLBACK_CONFIGURATION_ONLY_REJECTED:" + callbackCheck.failureReason);
    mode(false);
    require(importZhangPppArCheckpointSection(owner, runtime, hiddenPinTime, true, 2).valid,
        "LEGACY_FALSE_PRESENCE_BEHAVIOR_CHANGED");
    mode(true);

    // Public production operation creates a registered datum object with no
    // relation or datum rows. This must not pass by virtue of zero row count.
    (void)selectZhangE18CanonicalProductRelations(owner, E_Sys::GPS, {}, {}, 0);
    std::string ignored;
    rejected(exportZhangPppArCheckpointSection(owner, runtime, ignored),
        "live_e18PersistentProductDatumRegistries");
    require(ignored.empty(), "REJECTED_EXPORT_LEFT_PPP_PAYLOAD");
    mode(false);
    std::string registeredEmptyPpp;
    ppp = exportZhangPppArCheckpointSection(owner, runtime, registeredEmptyPpp);
    require(ppp.valid, "LEGACY_EMPTY_REGISTRY_EXPORT_FAILED:" + ppp.failureReason);
    require(importZhangPppArCheckpointSection(owner, runtime, emptyPpp, false, 2).valid,
        "LEGACY_EMPTY_REGISTRY_CLEANUP_FAILED");
    mode(true);
    rejected(importZhangPppArCheckpointSection(owner, runtime, registeredEmptyPpp, true, 2),
        "hasPersistentProductDatumRegistry");

    mode(false);
    const SatSys first(E_Sys::GPS, 13), second(E_Sys::GPS, 14);
    const auto selection = selectZhangE18CanonicalProductRelations(owner, E_Sys::GPS,
        {ZhangCanonicalSatelliteRelation::ordered(first, second)}, {first, second}, 1);
    std::string nonemptyPpp;
    ppp = exportZhangPppArCheckpointSection(owner, runtime, nonemptyPpp);
    require(ppp.valid && nonemptyPpp != registeredEmptyPpp,
        "LEGACY_NONEMPTY_DATUM_EXPORT_FAILED:" + ppp.failureReason);
    mode(true);
    rejected(exportZhangPppArCheckpointSection(owner, runtime, ignored),
        "live_e18PersistentProductDatumRegistries");
    mode(false);
    require(importZhangPppArCheckpointSection(owner, runtime, emptyPpp, false, 2).valid,
        "LEGACY_NONEMPTY_REGISTRY_CLEANUP_FAILED");
    mode(true);
    rejected(importZhangPppArCheckpointSection(owner, runtime, nonemptyPpp, true, 2),
        "hasPersistentProductDatumRegistry");

    auto& registry = zhangProductRelationAdmissionStateRegistry();
    require(registry.empty(), "ADMISSION_FIXTURE_NOT_ISOLATED");
    registry[{runtime, E_Sys::GPS}] = {};
    rejected(exportZhangAmbresCheckpointSection(owner, runtime, ignored),
        "live_zhangProductRelationAdmissionStateRegistry()");
    require(ignored.empty(), "REJECTED_EXPORT_LEFT_AMBRES_PAYLOAD");
    mode(false);
    std::string registeredEmptyAmbres;
    ambres = exportZhangAmbresCheckpointSection(owner, runtime, registeredEmptyAmbres);
    require(ambres.valid, "LEGACY_EMPTY_ADMISSION_EXPORT_FAILED:" + ambres.failureReason);
    registry.clear();
    mode(true);
    rejected(importZhangAmbresCheckpointSection(owner, runtime, registeredEmptyAmbres, true),
        "admissionStates");

    auto& state = registry[{runtime, E_Sys::GPS}];
    state.certifiedRows = {{{"PRIVATE_GATE_CONTROLLED_PHYSICAL_ARC", 1}}};
    state.certifiedValues = {7};
    state.certifiedSatellites = {"G13"};
    state.certifiedRelationIds = {"PRIVATE_GATE_CONTROLLED_RELATION"};
    require(ProductRelationAdmission::validateState(state, &reason), reason);
    rejected(exportZhangAmbresCheckpointSection(owner, runtime, ignored),
        "live_zhangProductRelationAdmissionStateRegistry()");
    mode(false);
    std::string nonemptyAmbres;
    ambres = exportZhangAmbresCheckpointSection(owner, runtime, nonemptyAmbres);
    require(ambres.valid && nonemptyAmbres != registeredEmptyAmbres,
        "LEGACY_NONEMPTY_ADMISSION_EXPORT_FAILED:" + ambres.failureReason);
    require(importZhangAmbresCheckpointSection(owner, runtime, nonemptyAmbres, true).valid,
        "LEGACY_NONEMPTY_ADMISSION_PREFLIGHT_CHANGED");
    registry.clear();
    mode(true);
    rejected(importZhangAmbresCheckpointSection(owner, runtime, nonemptyAmbres, true),
        "admissionStates");

    floatBookkeeping(owner, emptyPpp);
    disabledPriorCapture(owner, emptyAmbres);

    std::string afterPpp, afterAmbres;
    require(exportZhangPppArCheckpointSection(owner, runtime, afterPpp).valid && afterPpp == emptyPpp,
        "FAILED_NEGATIVES_MUTATED_PPP_RUNTIME");
    require(exportZhangAmbresCheckpointSection(owner, runtime, afterAmbres).valid && afterAmbres == emptyAmbres,
        "FAILED_NEGATIVES_MUTATED_AMBRES_RUNTIME");
    require((owner.x.array() == originalX.array()).all() &&
        (owner.P.array() == originalP.array()).all() && owner.metaDataMap == originalMetadata,
        "ADAPTER_CONTRACT_MUTATED_AUTHORITATIVE_X_P_METADATA");
    numericalCheckpointRoundtrip(owner,emptyPpp);
    boost::json::object result{
        {"status", "PASS"}, {"assertions", checks},
        {"actual_production_export_import", true}, {"producer_recompiled_in_driver", false},
        {"empty_carrier_export_preflight_commit", true}, {"callback_only_configuration_allowed", true},
        {"false_presence_hidden_time_rejected", true}, {"empty_live_registry_rejected", true},
        {"gec_float_graph_bookkeeping_roundtrip", true},
        {"all_allowed_affine_and_physical_metadata_preserved", true},
        {"singleton_birth_not_integer_authority", true},
        {"field_level_integer_authority_contamination_rejected", true},
        {"disabled_prior_capture_no_registration", true},
        {"disabled_existing_prior_clears_inner_preserves_registration", true},
        {"enabled_prior_capture_scalar_oracle_preserved", true},
        {"nonempty_private_history_rejected", true}, {"legacy_gps_public_adapter_controls", true},
        {"failed_import_side_effect_free", true},
        {"schema7_predictor_and_attestation_exact_roundtrip",true},
        {"schema6_predictor_cold_reset_explicit",true},
        {"scope", "CONTROLLED_PRIVATE_ADAPTER_CONTRACT_NOT_FULL_BUNDLE_OR_REAL_DATA"},
        {"uncovered", "Not every private container type independently injected; full bundle and real-data resume are separate gates"}};
    std::cout << "PRIVATE_CHECKPOINT_GATE_RESULT " << boost::json::serialize(result) << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "PRIVATE_CHECKPOINT_GATE_FAILED assertions=" << checks << " reason=" << error.what() << '\n';
    return 1;
}
