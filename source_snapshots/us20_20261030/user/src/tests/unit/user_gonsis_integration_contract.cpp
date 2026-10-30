#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/biases.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserGonsisEpoch.hpp"
#include "pea/zhangPppAr.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
bool sameMoments(const KFState& first, const KFState& second)
{
    return first.x.size() == second.x.size() && first.P.rows() == second.P.rows() &&
        first.P.cols() == second.P.cols() && first.dx.size() == second.dx.size() &&
        std::memcmp(first.x.data(), second.x.data(), first.x.size() * sizeof(double)) == 0 &&
        std::memcmp(first.P.data(), second.P.data(), first.P.size() * sizeof(double)) == 0 &&
        std::memcmp(first.dx.data(), second.dx.data(), first.dx.size() * sizeof(double)) == 0;
}

// These are synthetic committed FLOAT moments, not another implementation of
// the AR selector, WL/IF construction, ILS or output update. The actual public
// production entry point is linked below, without including its .inc source.
struct Fixture
{
    Receiver receiver;
    KFState state;
    KFMeas measurement;
    unsigned callbackCalls = 0;
    std::vector<KFKey> ambiguityKeys;

    explicit Fixture(unsigned id)
    {
        receiver.id = "GONSIS_UNIT_" + std::to_string(id);
        state.time.bigTime = 100000 + 30 * id;
        if (!bindZhangCheckpointRuntimeId(state, "gonsis-integration-" + std::to_string(id)))
            throw std::runtime_error("fixture runtime binding failed");
        const int dimension = 1 + 3 + 2 * 6;
        state.x = VectorXd::Zero(dimension);
        state.x(0) = 1;
        state.dx = VectorXd::Zero(dimension);
        state.P = MatrixXd::Identity(dimension, dimension) * 0.0001;
        state.P(0, 0) = 0;
        state.kfIndexMap.clear();
        state.kfIndexMap[state.oneKey] = 0;
        state.lastPppTransactionCommitted = true;
        state.lastPppQcConverged = true;
        state.factorCommitSequence = 77 + id;
        for (int axis = 0; axis < 3; ++axis)
        {
            KFKey key;
            key.type = KF::REC_POS;
            key.str = receiver.id;
            key.num = axis;
            key.rec_ptr = &receiver;
            state.kfIndexMap[key] = 1 + axis;
            state.stateTransitionMap[key][key][0] = 1;
            state.x(1 + axis) = axis == 0 ? 6378137 : 0;
        }
        KFMeasEntryList entries;
        for (int satelliteIndex = 0; satelliteIndex < 6; ++satelliteIndex)
        {
            const SatSys satellite(E_Sys::GPS, satelliteIndex + 1);
            auto& status = receiver.satStatMap[satellite];
            status.el = (70 - satelliteIndex * 5) * D2R;
            for (const auto frequency : {F1, F2}) status.sigStatMap[ft2string(frequency)] = {};
            for (const auto code : {E_ObsCode::L1C, E_ObsCode::L2W})
                status.sigStatMap[enum_to_string(code)] = {};
            GObs observation;
            observation.Sat = satellite;
            observation.time = state.time;
            observation.satStat_ptr = &status;
            observation.sigs[F1].code = E_ObsCode::L1C;
            observation.sigs[F2].code = E_ObsCode::L2W;
            observation.sigs[F1].L = 1e8;
            observation.sigs[F2].L = 8e7;
            observation.sigs[F1].P = observation.sigs[F2].P = 2e7;
            receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(observation));
            for (int band = 0; band < 2; ++band)
            {
                const auto code = band ? E_ObsCode::L2W : E_ObsCode::L1C;
                KFKey ambiguity;
                ambiguity.type = KF::AMBIGUITY;
                ambiguity.str = receiver.id;
                ambiguity.Sat = satellite;
                ambiguity.num = static_cast<int>(code);
                ambiguity.rec_ptr = &receiver;
                const int column = 4 + 2 * satelliteIndex + band;
                state.kfIndexMap[ambiguity] = column;
                state.stateTransitionMap[ambiguity][ambiguity][0] = 1;
                state.x(column) = (band ? 85 + satelliteIndex : 100 + 2 * satelliteIndex) +
                    (band ? .0005 : .001) * satelliteIndex;
                ambiguityKeys.push_back(ambiguity);
                BiasEntry bias;
                bias.Sat = satellite;
                bias.cod1 = code;
                bias.cod2 = E_ObsCode::NONE;
                bias.measType = PHAS;
                bias.source = "bsx";
                bias.tini.bigTime = 99900;
                bias.tfin.bigTime = 200000;
                bias.refTime.bigTime = 150000;
                bias.bias = .123;
                bias.var = .004;
                biasMaps[PHAS][satellite.id() + ":G"][code][E_ObsCode::NONE][bias.tini] = bias;
                for (const auto type : {KF::CODE_MEAS, KF::PHAS_MEAS})
                {
                    KFMeasEntry entry(&state);
                    entry.obsKey.type = type;
                    entry.obsKey.str = receiver.id;
                    entry.obsKey.Sat = satellite;
                    entry.obsKey.num = static_cast<int>(code);
                    entry.obsKey.rec_ptr = &receiver;
                    for (const auto& [key, index] : state.kfIndexMap)
                        if (key.type == KF::REC_POS)
                            entry.designEntryMap[key] = .2 + .1 * key.num;
                    if (type == KF::PHAS_MEAS)
                        entry.designEntryMap[ambiguity] = genericWavelength.at(band ? F2 : F1);
                    entry.addNoiseEntry(entry.obsKey, 1, type == KF::PHAS_MEAS ? .003 * .003 : .3 * .3);
                    entry.setValue(0);
                    entries.push_back(entry);
                }
            }
        }
        measurement = KFMeas(state, entries, state.time);
        measurement.V.setZero();
        measurement.VV.setZero();
        state.acceptedMeasurementFactorCallback = [&](auto&&...) { ++callbackCalls; return true; };
        state.stateTransitionFactorCallback = [&](auto&&...) { ++callbackCalls; return true; };
        state.exactStateTransformCallback = [&](auto&&...) { ++callbackCalls; return true; };
        state.authoritativeFactorTransactionCallback = [&](auto&&...) { ++callbackCalls; return true; };
    }
};
}

extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool value, const char* label)
    {
        ++checks;
        if (!value) throw std::runtime_error(std::string("FAIL ") + label);
    };
    try
    {
        for (const char* name : {"ZHANG_USER_AR_MIN_ARC_SECONDS", "ZHANG_USER_AR_MIN_ARC_SAMPLES",
            "ZHANG_USER_AR_CONFIRM_EPOCHS", "ZHANG_USER_AR_MAX_CONDITIONAL_SIGMA_CYCLES",
            "ZHANG_USER_AR_PAR_MIN_RANK", "ZHANG_USER_AR_PAR_MAX_TRIALS",
            "ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH"}) unsetenv(name);
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "GONSIS_UDUC", 1);
        setenv("ZHANG_USER_AR_SEARCH", "GONSIS_UDUC", 1);
        setenv("ZHANG_USER_QC", "CAUSAL_TRACKING", 1);
        setenv("ZHANG_USER_AR_FEEDBACK", "EPOCH_ONLY", 1);
        setenv("ZHANG_USER_ENFORCE_AR_ELEVATION", "1", 1);
        setenv("ZHANG_R51_RATIO_ONLY", "1", 1);
        setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS", "0", 1);
        acsConfig.zhangPppAr.user_adapter = false;
        acsConfig.zhangPppAr.output_products = false;
        acsConfig.zhangFullRank.enable = false;
        acsConfig.phaseClockOsb.enable = false;
        acsConfig.ssrInOpts.one_freq_phase_bias = false;
        acsConfig.solve_amb_for[E_Sys::GPS] = true;
        acsConfig.code_priorities[E_Sys::GPS] = {E_ObsCode::L1C, E_ObsCode::L2W};
        acsConfig.ambrOpts.mode = E_ARmode::LAMBDA_ALT;
        acsConfig.ambrOpts.elevation_mask_deg = 15;
        acsConfig.ambrOpts.ratioThres = 2;
        acsConfig.ambrOpts.succsThres = .999;
        acsConfig.ambrOpts.lambda_set = 2;
        acsConfig.ambrOpts.AR_max_itr = 100;
        acsConfig.assign_closest_epoch = false;
        acsConfig.epoch_tolerance = .5;
        auto verifyIsolation = [&](const Fixture& fixture, const KFState& before, const KFState& output)
        {
            check(sameMoments(fixture.state, before), "authoritative FLOAT moments bitwise unchanged");
            check(fixture.state.kfIndexMap == before.kfIndexMap &&
                fixture.state.factorCommitSequence == before.factorCommitSequence &&
                fixture.state.metaDataMap == before.metaDataMap,
                "FLOAT layout sequence metadata unchanged");
            check(fixture.callbackCalls == 0, "private AR invoked no authoritative factor callbacks");
            check(output.factorCommitSequence == before.factorCommitSequence &&
                output.metaDataMap == before.metaDataMap,
                "output publishes no integer history or runtime bookkeeping changes");
            check(bool(output.acceptedMeasurementFactorCallback) &&
                bool(output.stateTransitionFactorCallback) && bool(output.exactStateTransformCallback) &&
                bool(output.authoritativeFactorTransactionCallback) && output.factorCallbackOwnerOverride == nullptr,
                "private bookkeeping callbacks restored on return");
        };
        Fixture success(1);
        const KFState before(success.state);
        captureZhangUserGonsisFloatObservations(success.state, success.measurement);
        KFState output(success.state);
        std::ostringstream successTrace;
        fixEpochOnlyUserAmbiguities(successTrace, output, success.state);
        check((output.x - before.x).norm() > 1e-5 && output.P.trace() < before.P.trace(),
            "production GONSIS changes only private ambiguity moments with soft covariance reduction");
        verifyIsolation(success, before, output);
        check(successTrace.str().find("GONSIS") != std::string::npos &&
            successTrace.str().find("nis_evaluations=0") != std::string::npos,
            "production strategy and zero NIS explicit in receipt");
        check(successTrace.str().find("final_rank=10 wl_rank=5 l1_rank=5") != std::string::npos &&
            successTrace.str().find("status=FULL_SUBSET_ACCEPTED") != std::string::npos &&
            successTrace.str().find("committed_rank=10 reason=PASS") != std::string::npos,
            "complete five-pair WL and L1 search passes the actual production output transaction");

        // Each runtime is new: previous captures cannot turn a missing packet
        // into current measurement authority. Corrupt only one required gate.
        for (unsigned fault = 0; fault < 11; ++fault)
        {
            Fixture fixture(10 + fault);
            if (fault == 1) fixture.measurement.time.bigTime -= 30;
            if (fault == 2) fixture.state.lastPppTransactionCommitted = false;
            if (fault == 3) fixture.state.lastPppQcConverged = false;
            if (fault == 6)
                for (const auto& [key, index] : fixture.state.kfIndexMap)
                    if (key.type == KF::REC_POS) fixture.state.P(index, index) = 1;
            if (fault == 7) fixture.measurement.VV(1) = 1;
            if (fault == 8)
            {
                // A retained G06/L1 PHASE row passes FLOAT's 4-sigma bound,
                // but the fixed displacement worsens it beyond that bound.
                // This demonstrates rollback AFTER a complete search, rather
                // than rejection merely because the original VV was invalid.
                fixture.measurement.VV(21) = .012 - 1e-6;
                check(std::abs(fixture.measurement.VV(21)) /
                    std::sqrt(fixture.measurement.R(21,21)) < 4,
                    "residual-boundary fixture initially passes FLOAT row threshold");
            }
            if (fault == 9)
            {
                for (const auto& [key, index] : fixture.state.kfIndexMap)
                    if (key.type == KF::REC_POS) fixture.state.P(index, index) = .01;
                check(std::sqrt(.01) <= .5 && std::sqrt(3 * .01) >= .15,
                    "position fixture independently passes FLOAT but fails fixed sigma thresholds");
            }
            if (fault != 0) captureZhangUserGonsisFloatObservations(fixture.state, fixture.measurement);
            if (fault == 4) fixture.state.factorCommitSequence++;
            if (fault == 5) fixture.state.x(4) += .001;
            if (fault == 10) fixture.state.P(4,4) += .00005;
            const KFState frozen(fixture.state);
            KFState rejected(fixture.state);
            std::ostringstream trace;
            fixEpochOnlyUserAmbiguities(trace, rejected, fixture.state);
            check(sameMoments(rejected, frozen), "missing or stale packet/precision/residual gate remains FLOAT");
            verifyIsolation(fixture, frozen, rejected);
            if (fault == 7 || fault == 8)
                check(trace.str().find("GONSIS_FIXED_RESIDUAL_4SIGMA_FAILED") != std::string::npos &&
                    trace.str().find("status=ROLLED_BACK committed_rank=0") != std::string::npos,
                    "fixed CODE/PHASE residual failure rolls back complete candidate transaction");
            if (fault == 9)
                check(trace.str().find("GONSIS_FIXED_POSITION_SIGMA_FAILED") != std::string::npos &&
                    trace.str().find("status=ROLLED_BACK committed_rank=0") != std::string::npos,
                    "fixed position gate is stricter than passing FLOAT position gate");
            if (fault == 10)
                check(trace.str().find("status=SAME_EPOCH_OBSERVATION_EVIDENCE_UNAVAILABLE") != std::string::npos &&
                    trace.str().find("ZHANG_USER_GONSIS_L1_SEARCH") == std::string::npos,
                    "P-only mutation cannot reuse same mean time sequence layout observation evidence");
        }
        // Capture columns are semantic, not merely an equal-sized vector.
        Fixture layout(30);
        captureZhangUserGonsisFloatObservations(layout.state, layout.measurement);
        std::swap(layout.state.kfIndexMap[layout.ambiguityKeys[0]],
                  layout.state.kfIndexMap[layout.ambiguityKeys[1]]);
        const KFState frozenLayout(layout.state);
        KFState layoutOutput(layout.state);
        std::ostringstream layoutTrace;
        fixEpochOnlyUserAmbiguities(layoutTrace, layoutOutput, layout.state);
        check(sameMoments(layoutOutput, frozenLayout), "same-sized changed KF layout cannot reuse old H");
        verifyIsolation(layout, frozenLayout, layoutOutput);
        std::cout << "USER_GONSIS_INTEGRATION_CONTRACT checks=" << checks
                  << " status=PASS real_production_entry=1 private_only=1 nis_evaluations=0"
                  << " synthetic_contract_not_real_data_quality=1\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
