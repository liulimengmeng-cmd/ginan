#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/biases.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/trace.hpp"
#include "common/zhangUserArExperiment.hpp"
#include "common/zhangUserOutputValidation.hpp"
#include "pea/zhangPppAr.hpp"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>

// Compile the exact production adapter, using real KF/observation/BiasEntry
// types and linked production query methods. The local helpers only pin a
// deterministic candidate policy/runtime; they do not mimic the adapter.
static double zhangE27Wavelength(E_Sys system, E_ObsCode code)
{
    const auto frequency = code2Freq.at(system).at(code);
    return genericWavelength.at(frequency);
}
static bool zhangUserArElevationEnabled() { return true; }
static string zhangAmbresRuntimeId(const KFState&) { return "adapter-contract"; }
static ZhangUserArExperimentArcTracker& zhangUserExperimentArcs()
{
    static ZhangUserArExperimentPolicy policy = []
    {
        ZhangUserArExperimentPolicy value;
        value.minimumArcSeconds = 0;
        value.minimumArcSamples = 1;
        return value;
    }();
    static ZhangUserArExperimentArcTracker tracker(policy);
    return tracker;
}
#include "pea/zhangUserSharedAdapter.inc"

extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool accepted, const char* label)
    {
        ++checks;
        if (!accepted) throw std::runtime_error(string("FAIL ") + label);
    };
    try
    {
        acsConfig.zhangPppAr.user_adapter = false;
        acsConfig.zhangPppAr.output_products = false;
        acsConfig.zhangFullRank.enable = false;
        acsConfig.phaseClockOsb.enable = false;
        acsConfig.ssrInOpts.one_freq_phase_bias = false;
        acsConfig.code_priorities[E_Sys::GPS] = {
            E_ObsCode::L1C, E_ObsCode::L1W, E_ObsCode::L2W};
        acsConfig.ambrOpts.elevation_mask_deg = 15;
        acsConfig.assign_closest_epoch = false;
        acsConfig.epoch_tolerance = 0.5;
        biasMaps[PHAS].clear();
        Receiver receiver;
        receiver.id = "TEST";
        KFState state;
        state.time.bigTime = 100000;
        state.x = VectorXd::Zero(3);
        state.x(1) = 100.25;
        state.x(2) = 95.75;
        state.P = MatrixXd::Identity(3, 3);
        state.dx = VectorXd::Zero(3);
        state.kfIndexMap.clear();
        state.kfIndexMap[state.oneKey] = 0;
        const SatSys satellite("G01");
        KFKey first;
        first.type = KF::AMBIGUITY;
        first.str = receiver.id;
        first.Sat = satellite;
        first.num = static_cast<int>(E_ObsCode::L1C);
        first.rec_ptr = &receiver;
        KFKey second = first;
        second.num = static_cast<int>(E_ObsCode::L2W);
        state.kfIndexMap[first] = 1;
        state.kfIndexMap[second] = 2;
        state.stateTransitionMap[first][first][0] = 1;
        state.stateTransitionMap[second][second][0] = 1;
        auto& status = receiver.satStatMap[satellite];
        status.el = 30 * D2R;
        status.sigStatMap[ft2string(F1)] = {};
        status.sigStatMap[ft2string(F2)] = {};
        status.sigStatMap[enum_to_string(E_ObsCode::L1C)] = {};
        status.sigStatMap[enum_to_string(E_ObsCode::L2W)] = {};
        GObs observed;
        observed.Sat = satellite;
        observed.time = state.time;
        observed.satStat_ptr = &status;
        observed.sigs[F1].code = E_ObsCode::L1C;
        observed.sigs[F1].L = 1e8;
        observed.sigs[F1].P = 2e7;
        observed.sigs[F2].code = E_ObsCode::L2W;
        observed.sigs[F2].L = 8e7;
        observed.sigs[F2].P = 2e7;
        const auto observation = static_cast<std::shared_ptr<GObs>>(observed);
        receiver.obsList.push_back(observation);
        for (const auto code : {E_ObsCode::L1C, E_ObsCode::L2W})
        {
            BiasEntry entry;
            entry.Sat = satellite;
            entry.cod1 = code;
            entry.cod2 = E_ObsCode::NONE;
            entry.measType = PHAS;
            entry.source = "bsx";
            entry.tini.bigTime = 99900;
            entry.tfin.bigTime = 100900;
            entry.refTime.bigTime = 100400;
            entry.bias = 0.123;
            entry.var = 0.004;
            biasMaps[PHAS]["G01:G"][code][E_ObsCode::NONE][entry.tini] = entry;
        }
        std::ostringstream trace;
        auto endpoints = [&]
        {
            return buildZhangUserSharedEndpoints(trace, state, &receiver, E_Sys::GPS, state.time);
        };
        const auto beforeX = state.x;
        const auto beforeP = state.P;
        const auto beforeObs = observation->sigs;
        auto initial = endpoints();
        check(initial.size() == 1 && initial.front().eligible, "legal dual-band endpoint");
        check(initial.front().firstCode == E_ObsCode::L1C &&
              initial.front().secondCode == E_ObsCode::L2W, "skip same-frequency alias priority");
        check(initial.front().firstRow.size() == 3 && initial.front().firstRow(1) == 1 &&
              initial.front().secondRow(2) == 1, "full-state cycle coordinate rows");
        check(!initial.front().firstReferenceZero && !initial.front().secondReferenceZero,
              "native external ordinary states are not certified chart zeros");
        check(initial.front().firstRow.dot(state.x) == 100.25 &&
              initial.front().secondRow.dot(state.x) == 95.75, "no meter/cycle remapping");
        const string initialIdentity = initial.front().identity;
        status.el = 14 * D2R;
        auto low = endpoints();
        check(low.size() == 1 && !low.front().eligible &&
              low.front().identity == initialIdentity, "candidate-only elevation retains physical endpoint");
        status.el = 30 * D2R;
        observation->sigs[F1].LLI = 2;
        check(endpoints().empty(), "half-cycle blocks integer endpoint");
        observation->sigs[F1].LLI = 0;
        check(endpoints().size() == 1, "whole-cycle recovery");
        status.sigStatMap[ft2string(F2)].slip.GF = true;
        check(endpoints().empty(), "current physical slip invalidates endpoint");
        status.sigStatMap[ft2string(F2)].slip.any = 0;
        status.sigStatMap[enum_to_string(E_ObsCode::L1C)].phaseRejectCount = 1;
        check(endpoints().empty(), "rejected current phase not carried");
        status.sigStatMap[enum_to_string(E_ObsCode::L1C)].phaseRejectCount = 0;
        observation->time.bigTime -= 1;
        check(endpoints().empty(), "stale epoch not an endpoint");
        observation->time = state.time;
        auto& firstEntry = biasMaps[PHAS]["G01:G"][E_ObsCode::L1C][E_ObsCode::NONE].begin()->second;
        const auto actualFinish = firstEntry.tfin;
        firstEntry.tfin = state.time;
        check(endpoints().empty(), "OSB end interval is exclusive no extrapolated AR");
        firstEntry.tfin = actualFinish;
        firstEntry.source = "KALMAN";
        check(endpoints().empty(), "non-BSX source is not external phase OSB evidence");
        firstEntry.source = "bsx";
        const double actualVariance = firstEntry.var;
        firstEntry.var = -1;
        check(endpoints().empty(), "negative bias variance rejects AR only");
        firstEntry.var = actualVariance;
        state.stateTransitionMap.erase(second);
        check(endpoints().empty(), "pending removed ordinary state not materialized as zero");
        state.stateTransitionMap[second][second][0] = 1;
        check(endpoints().size() == 1, "ordinary state restored");
        const auto originalIdentity = endpoints().front().identity;
        const auto generationKey = zhangUserArcGenerationKey(false, "TEST", "G01", first.num);
        check(zhangUserArcGeneration(state.metaDataMap, generationKey) == 0, "initial generation");
        state.removeState(first);
        check(zhangUserArcGeneration(state.metaDataMap, generationKey) == 1, "actual native reset hook");
        state.removeState(first);
        check(zhangUserArcGeneration(state.metaDataMap, generationKey) == 1, "duplicate remove not another physical reset");
        state.stateTransitionMap[first][first][0] = 1;
        const auto restarted = endpoints();
        check(restarted.size() == 1 && restarted.front().identity != originalIdentity,
              "same-key remove/re-add cannot inherit old HOLD identity");
        state.metaDataMap["test.reference-chart-only"] = "G02";
        check(endpoints().front().identity == restarted.front().identity,
              "reference-neutral endpoint identity");
        check((state.x - beforeX).norm() == 0 && (state.P - beforeP).norm() == 0,
              "adapter and removal scheduling preserve exported FLOAT moments");
        check(observation->sigs[F1].L == beforeObs.at(F1).L &&
              observation->sigs[F2].P == beforeObs.at(F2).P,
              "raw observations unchanged");
        // The same generic KF removal path is also used by internal QC resets.
        // In contrast, exact chart transport must not invent a physical break.
        acsConfig.zhangPppAr.user_adapter = true;
        KFState internal = state;
        const auto physicalKey = zhangUserArcGenerationKey(true, "TEST", "G01", first.num);
        check(zhangUserArcGeneration(internal.metaDataMap, physicalKey) == 0,
              "internal physical generation starts separately");
        internal.removeState(first);
        check(zhangUserArcGeneration(internal.metaDataMap, physicalKey) == 1,
              "actual internal generic QC reset hook");
        KFState chart = state;
        map<KFKey, map<KFKey, double>> chartTransform;
        chartTransform[state.oneKey][state.oneKey] = 1;
        chartTransform[second][second] = 1;
        KFKey replacement = first;
        replacement.Sat = SatSys("G02");
        chartTransform[replacement][first] = -1;
        check(chart.applyStateTransform(trace, chartTransform, "unit reference chart exchange"),
              "actual exact reference chart transport");
        check(zhangUserArcGeneration(chart.metaDataMap, physicalKey) == 0 &&
              zhangUserArcGeneration(chart.metaDataMap,
                  zhangUserArcGenerationKey(true, "TEST", "G02", first.num)) == 0,
              "chart exchange is not a physical ambiguity reset");
        acsConfig.zhangPppAr.user_adapter = false;
        state.metaDataMap[generationKey] = "not-an-integer";
        bool invalidCounter = false;
        try { (void)endpoints(); }
        catch (const std::runtime_error&) { invalidCounter = true; }
        check(invalidCounter, "generation metadata fails closed");

        // Independent hand calculation, not the production conditioner as an
        // oracle: prior mean=(0,0), P=[[2,1],[1,2]], A=(1,-1), z=3.
        // P A'=(1,-1), A P A'=2, so conditioned mean=(1.5,-1.5)
        // and P_fixed=[[1.5,1.5],[1.5,1.5]], with eigenvalues 0 and 3.
        VectorXd fixedMean(2); fixedMean << 1.5,-1.5;
        MatrixXd fixedCovariance(2,2); fixedCovariance << 1.5,1.5,1.5,1.5;
        MatrixXd fixedRows(2,2); fixedRows << 1,-1,2,-2;
        VectorXd fixedIntegers(2); fixedIntegers << 3,6;
        const auto fixed = zhangValidateUserIntegerOutput(
            fixedMean,fixedCovariance,fixedRows,fixedIntegers);
        check(fixed.valid && fixed.maximumMeanResidual==0 && fixed.maximumCovarianceResidual==0,
              "hand-derived posterior satisfies all redundant integer rows");
        check(std::abs(fixed.minimumEigenvalue)<=fixed.psdTolerance,
              "hand-derived singular posterior zero eigenvalue is legal");
        const auto fixedMeanBefore = fixedMean;
        const auto fixedCovarianceBefore = fixedCovariance;
        auto wrongMean = fixedMean;
        wrongMean(0)+=0.25;
        const auto wrongClosure = zhangValidateUserIntegerOutput(
            wrongMean,fixedCovariance,fixedRows,fixedIntegers);
        check(!wrongClosure.valid && wrongClosure.reason=="OUTPUT_POSTERIOR_MEAN_CLOSURE_FAILED",
              "finite but unclosed accepted integer posterior rejects");
        const auto floatCovariance = MatrixXd::Identity(2,2).eval();
        const auto wrongNullspace = zhangValidateUserIntegerOutput(
            fixedMean,floatCovariance,fixedRows,fixedIntegers);
        check(!wrongNullspace.valid && wrongNullspace.reason=="OUTPUT_POSTERIOR_INTEGER_NULLSPACE_FAILED",
              "correct integer mean with unconditioned covariance rejects");
        auto indefinite = fixedCovariance;
        indefinite(0,0)=1; indefinite(1,1)=1;
        indefinite(0,1)=2; indefinite(1,0)=2; // independently: eigenvalues -1,3.
        const auto wrongPsd = zhangValidateUserIntegerOutput(
            fixedMean,indefinite,fixedRows,fixedIntegers);
        check(!wrongPsd.valid && wrongPsd.reason=="OUTPUT_POSTERIOR_NOT_PSD" &&
              std::abs(wrongPsd.minimumEigenvalue+1)<1e-12,
              "finite covariance with hand-known negative eigenvalue rejects");
        auto asymmetric = fixedCovariance;
        asymmetric(0,1)+=0.25;
        const auto wrongSymmetry = zhangValidateUserIntegerOutput(
            fixedMean,asymmetric,fixedRows,fixedIntegers);
        check(!wrongSymmetry.valid && wrongSymmetry.reason=="OUTPUT_POSTERIOR_ASYMMETRIC",
              "posterior asymmetry is rejected not silently repaired");
        auto wrongShape = fixedRows;
        wrongShape.conservativeResize(2,3);
        check(!zhangValidateUserIntegerOutput(fixedMean,fixedCovariance,wrongShape,fixedIntegers).valid,
              "posterior illegal constraint shape rejects before multiplication");
        check(!zhangValidateUserIntegerOutput(fixedMean,MatrixXd::Identity(3,3),fixedRows,fixedIntegers).valid,
              "posterior illegal covariance shape rejects");
        check(!zhangValidateUserIntegerOutput(fixedMean,fixedCovariance,fixedRows,VectorXd::Zero(1)).valid,
              "posterior RHS dimension mismatch rejects");
        auto nonfinite = fixedMean;
        nonfinite(0)=std::numeric_limits<double>::quiet_NaN();
        check(!zhangValidateUserIntegerOutput(nonfinite,fixedCovariance,fixedRows,fixedIntegers).valid,
              "nonfinite posterior mean rejects");
        auto nonfiniteP = fixedCovariance;
        nonfiniteP(1,1)=std::numeric_limits<double>::infinity();
        check(!zhangValidateUserIntegerOutput(fixedMean,nonfiniteP,fixedRows,fixedIntegers).valid,
              "nonfinite posterior covariance rejects");
        auto nonfiniteA = fixedRows;
        nonfiniteA(0,0)=std::numeric_limits<double>::quiet_NaN();
        check(!zhangValidateUserIntegerOutput(fixedMean,fixedCovariance,nonfiniteA,fixedIntegers).valid,
              "nonfinite integer functional rejects");
        auto nonfiniteZ = fixedIntegers;
        nonfiniteZ(0)=std::numeric_limits<double>::infinity();
        check(!zhangValidateUserIntegerOutput(fixedMean,fixedCovariance,fixedRows,nonfiniteZ).valid,
              "nonfinite integer RHS rejects");
        const auto emptyRows = MatrixXd::Zero(0,2).eval();
        check(zhangValidateUserIntegerOutput(fixedMean,fixedCovariance,emptyRows,VectorXd::Zero(0)).valid,
              "empty integer space validates covariance without zero-size BLAS");
        check(!zhangValidateUserIntegerOutput(VectorXd::Zero(0),MatrixXd::Zero(0,0),
                  MatrixXd::Zero(0,0),VectorXd::Zero(0)).valid,
              "empty stochastic state has no valid output posterior");
        auto roundingMean = fixedMean;
        roundingMean(0)+=std::numeric_limits<double>::epsilon();
        check(zhangValidateUserIntegerOutput(roundingMean,fixedCovariance,fixedRows,fixedIntegers).valid,
              "machine-scale affine rounding has an explicit finite tolerance");
        const auto scaledCovariance = (1e12*fixedCovariance).eval();
        check(zhangValidateUserIntegerOutput(fixedMean,scaledCovariance,fixedRows,fixedIntegers).valid,
              "large but closed PSD covariance uses scale-aware numerical bounds");
        check((fixedMean.array()==fixedMeanBefore.array()).all() &&
              (fixedCovariance.array()==fixedCovarianceBefore.array()).all(),
              "output numerical gate never mutates posterior inputs");
        std::cout << "PASS user_shared_adapter_contract checks=" << checks << '\n';
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
