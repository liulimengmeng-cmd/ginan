#pragma once

#include "common/eigenIncluder.hpp"
#include "common/zhangUserGonsisReplay.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Independently implemented from the UDUC decisions in the supplied
// GONSIS-master/src/ppp_ar.cpp (sdidx_by_ele, resamb_wl, calculate_nl,
// fix_sol_wl, resamb_nl and fix_sol_uc). Its file attributes the
// RTKLIB ancestor to T. Takasu, 2012-2015. No GONSIS source is copied here.
// This core consumes one frozen ambiguity posterior; it consumes no raw
// observations, performs no triple filtering, and stores no integer history.

struct ZhangUserGonsisSatellite
{
    std::string satellite;
    // Caller joins the physical arc, KF-state and product-datum identities.
    std::string identity;
    double elevationRadians = 0;
    // Legacy name retained for this small core API: the caller supplies the
    // source-style published-FIX lock counter, NOT ordinary tracking samples.
    int effectiveSamples = 0;
    // GONSIS UDUC ranks the sum of the TWO marginal metre variances.
    // This is not Var(N1+N2), nor the unscaled sum of cycle variances.
    double metricVarianceMetres2 = std::numeric_limits<double>::quiet_NaN();
    bool eligible = true;
};

struct ZhangUserGonsisReferenceHint
{
    std::string satellite;
    std::string identity;
    bool previousPublishedFix = false;
};

struct ZhangUserGonsisTrial;
struct ZhangUserGonsisReplayView;

struct ZhangUserGonsisOptions
{
    double wideLaneResidualCycles = 0.25;
    int minimumWideLanePairs = 4;
    double ratioThreshold = 2;
    // Zero means the finite number of available one-satellite removals.
    int maximumTrials = 0;
    double wideLaneConstraintVarianceCycles2 = 1e-8;
    double ifConstraintVarianceMetres2 = 0.003 * 0.003;
    double wavelength1Metres = 299792458.0 / 1575.42e6;
    double wavelength2Metres = 299792458.0 / 1227.60e6;
    // Safe local reuse requires complete bit-identical search moments AND
    // actual targets. No cross-epoch or satellite-name cache is permitted.
    bool reuseIdenticalIls = true;
    bool targetFailedL1SubsetForPar = true;
    // Synchronous borrowed views only: no matrix copies when not requested.
    std::function<void(const ZhangUserGonsisTrial&, const ZhangUserGonsisReplayView&)> replayTrial;
};

struct ZhangUserGonsisIlsSolution
{
    // Must certify a complete two-candidate search of EVERY requested row.
    bool complete = false;
    VectorXd integers;
    double bestSquaredCost = std::numeric_limits<double>::quiet_NaN();
    double secondSquaredCost = std::numeric_limits<double>::quiet_NaN();
    std::string reason;
    ZhangUserGonsisLambdaTopTwo replayTopTwo;
};

using ZhangUserGonsisIlsSolver = std::function<ZhangUserGonsisIlsSolution(
    const VectorXd&, const MatrixXd&)>;

struct ZhangUserGonsisTrial
{
    std::vector<int> satelliteIndices;
    int referenceIndex = -1;
    std::string referenceSatellite;
    std::string referenceIdentity;
    std::vector<int> wideLaneCandidateSatelliteIndices;
    std::vector<double> wideLaneCandidateResidualsCycles;
    // Actual dot/round values, recorded only when the replay callback is on.
    std::vector<double> wideLaneCandidateMeansCycles;
    std::vector<double> wideLaneCandidateIntegers;
    std::vector<int> wideLaneSatelliteIndices;
    int requestedL1Rank = 0;
    bool ilsComplete = false;
    double ratio = std::numeric_limits<double>::quiet_NaN();
    std::string solverReason;
    std::string reason;
    // The actual next removal, not a repeated priority-queue top.
    int removedSatelliteIndex = -1;
    std::string removalPool = "NONE";
    bool ilsCacheReused = false;
};

struct ZhangUserGonsisReplayView
{
    const VectorXd& rawFrozenMean;
    const MatrixXd& rawFrozenCovariance;
    const MatrixXd& usedFrozenCovariance;
    const MatrixXd& wideLaneRows;
    const VectorXd& wideLaneIntegers;
    const VectorXd& wlConditionedMean;
    const MatrixXd& wlConditionedCovariance;
    const MatrixXd& l1Rows;
    const MatrixXd& ifRows;
    const VectorXd& l1Mean;
    const MatrixXd& l1Covariance;
    const ZhangUserGonsisIlsSolution& solution;
};

struct ZhangUserGonsisResult
{
    bool accepted = false;
    std::string reason = "NOT_EVALUATED";
    int referenceIndex = -1;
    std::string referenceSatellite;
    std::string referenceIdentity;
    // Reference followed by the rovers in the complete accepted L1 search.
    std::vector<int> acceptedSatelliteIndices;
    std::vector<int> removedSatelliteIndices;
    std::vector<ZhangUserGonsisTrial> trials;
    // All rows address the ORIGINAL interleaved cycle vector [N1,N2,...].
    MatrixXd wideLaneRows;
    VectorXd wideLaneIntegers;
    MatrixXd l1Rows;
    VectorXd l1Integers;
    // Source-compatible affine IF-minus-fixed-WL search moments. Because WL
    // is soft, these need not exactly equal l1Rows * wlConditionedMean.
    VectorXd l1FloatMean;
    MatrixXd l1FloatCovariance;
    VectorXd wlConditionedMean;
    MatrixXd wlConditionedCovariance;
    // Private output conditioning: WL rows in cycles, then IF rows in metres.
    MatrixXd constraintRows;
    VectorXd constraintValues;
    VectorXd constraintVariances;
    int ilsCalls = 0;
    int ilsCacheHits = 0;
};

inline double zhangUserGonsisRound(double value)
{
    // Deliberately differs from std::round for negative half integers.
    return std::floor(value + 0.5);
}

namespace zhangUserGonsisDetail
{
inline bool positiveSemidefinite(const MatrixXd& matrix)
{
    if (matrix.rows() == 0 || matrix.rows() != matrix.cols() || !matrix.allFinite())
        return false;
    Eigen::SelfAdjointEigenSolver<MatrixXd> spectrum(matrix);
    if (spectrum.info() != Eigen::Success || !spectrum.eigenvalues().allFinite())
        return false;
    const double scale = spectrum.eigenvalues().cwiseAbs().maxCoeff();
    const double tolerance = std::max(1e-14,
        scale * std::numeric_limits<double>::epsilon() * matrix.rows() * 64);
    // A certified omitted datum reference can give the ambient posterior
    // exact zero rows/columns. Check PSD without clipping eigenvalues or
    // adding jitter; the actual searched L1 covariance still must be PD.
    return spectrum.eigenvalues().minCoeff() >= -tolerance;
}

inline bool positiveDefinite(const MatrixXd& matrix)
{
    if (matrix.rows() == 0 || matrix.rows() != matrix.cols() || !matrix.allFinite())
        return false;
    Eigen::LLT<MatrixXd> solve(matrix);
    return solve.info() == Eigen::Success &&
        solve.matrixL().toDenseMatrix().diagonal().allFinite();
}

inline int reference(
    const std::vector<int>& active,
    const std::vector<ZhangUserGonsisSatellite>& satellites,
    const ZhangUserGonsisReferenceHint& hint)
{
    int highest = active.front();
    int previous = -1;
    for (int index : active)
    {
        if (satellites[index].elevationRadians > satellites[highest].elevationRadians)
            highest = index;
        if (hint.previousPublishedFix && !hint.identity.empty() &&
            satellites[index].satellite == hint.satellite &&
            satellites[index].identity == hint.identity &&
            satellites[index].elevationRadians > std::acos(-1.0) / 6)
            previous = index;
    }
    if (previous < 0) return highest;
    int selected = previous;
    for (int index : active)
        if (satellites[index].effectiveSamples > 200 &&
            satellites[index].elevationRadians > satellites[selected].elevationRadians)
            selected = index;
    return selected;
}

inline bool conditionWideLane(
    const VectorXd& frozenMean, const MatrixXd& frozenCovariance,
    const MatrixXd& rows, const VectorXd& integers, double variance,
    VectorXd& mean, MatrixXd& covariance)
{
    const MatrixXd cross = (frozenCovariance * rows.transpose()).eval();
    MatrixXd innovation = (rows * cross).eval();
    innovation.diagonal().array() += variance;
    innovation = (0.5 * (innovation + innovation.transpose())).eval();
    Eigen::LLT<MatrixXd> solve(innovation);
    if (solve.info() != Eigen::Success) return false;
    // Solve, never invert or jitter. This includes correlations with every
    // original satellite, even those omitted by a PAR trial.
    const MatrixXd gain = solve.solve(cross.transpose()).transpose();
    if (solve.info() != Eigen::Success || !gain.allFinite()) return false;
    mean = frozenMean + gain * (integers - rows * frozenMean);
    // Joseph form avoids cancellation with the 1e-8-cycle WL constraint.
    const MatrixXd residual = MatrixXd::Identity(frozenMean.size(), frozenMean.size()) - gain * rows;
    covariance = residual * frozenCovariance * residual.transpose() + variance * gain * gain.transpose();
    covariance = (0.5 * (covariance + covariance.transpose())).eval();
    return mean.allFinite() && covariance.allFinite();
}
}

inline ZhangUserGonsisResult zhangResolveUserGonsis(
    const VectorXd& frozenCycleMean,
    const MatrixXd& frozenCycleCovariance,
    const std::vector<ZhangUserGonsisSatellite>& satellites,
    const ZhangUserGonsisIlsSolver& solver,
    const ZhangUserGonsisOptions& options = {},
    const ZhangUserGonsisReferenceHint& referenceHint = {})
{
    ZhangUserGonsisResult result;
    const int dimension = frozenCycleMean.size();
    auto reject = [&](const std::string& reason)
    {
        result.reason = reason;
        // A rejected WL-only or incomplete L1 decision exposes no publishable
        // rows, integers or conditioned posterior.
        result.wideLaneRows.resize(0, dimension);
        result.l1Rows.resize(0, dimension);
        result.constraintRows.resize(0, dimension);
        return result;
    };
    if (dimension != 2 * static_cast<int>(satellites.size()) || dimension == 0 ||
        frozenCycleCovariance.rows() != dimension || frozenCycleCovariance.cols() != dimension ||
        !frozenCycleMean.allFinite() || !frozenCycleCovariance.allFinite() || !solver)
        return reject("INPUT_INVALID");
    if (!std::isfinite(options.wideLaneResidualCycles) || options.wideLaneResidualCycles < 0 ||
        options.wideLaneResidualCycles > 0.5 || options.minimumWideLanePairs < 4 ||
        !std::isfinite(options.ratioThreshold) || options.ratioThreshold <= 1 ||
        options.maximumTrials < 0 ||
        !std::isfinite(options.wideLaneConstraintVarianceCycles2) || options.wideLaneConstraintVarianceCycles2 <= 0 ||
        !std::isfinite(options.ifConstraintVarianceMetres2) || options.ifConstraintVarianceMetres2 <= 0 ||
        !std::isfinite(options.wavelength1Metres) || !std::isfinite(options.wavelength2Metres) ||
        options.wavelength1Metres <= 0 || options.wavelength2Metres <= options.wavelength1Metres)
        return reject("OPTIONS_INVALID");
    const double covarianceScale = std::max(1.0, frozenCycleCovariance.cwiseAbs().maxCoeff());
    if ((frozenCycleCovariance - frozenCycleCovariance.transpose()).cwiseAbs().maxCoeff() >
        covarianceScale * 1e-10)
        return reject("COVARIANCE_NOT_SYMMETRIC");
    const MatrixXd frozenCovariance =
        (0.5 * (frozenCycleCovariance + frozenCycleCovariance.transpose())).eval();
    if (!zhangUserGonsisDetail::positiveSemidefinite(frozenCovariance))
        return reject("COVARIANCE_NOT_POSITIVE_SEMIDEFINITE");
    std::set<std::string> names;
    std::vector<int> active;
    for (int index = 0; index < static_cast<int>(satellites.size()); ++index)
    {
        const auto& satellite = satellites[index];
        if (satellite.satellite.empty() || satellite.identity.empty() ||
            !names.insert(satellite.satellite).second ||
            !std::isfinite(satellite.elevationRadians) || satellite.effectiveSamples < 0 ||
            !std::isfinite(satellite.metricVarianceMetres2) || satellite.metricVarianceMetres2 < 0)
            return reject("SATELLITE_METADATA_INVALID");
        if (satellite.eligible) active.push_back(index);
    }
    if (static_cast<int>(active.size()) <= options.minimumWideLanePairs)
        return reject("INSUFFICIENT_ELIGIBLE_SATELLITES");

    const double lambda1 = options.wavelength1Metres;
    const double lambda2 = options.wavelength2Metres;
    const double lambdaNl = lambda1 * lambda2 / (lambda1 + lambda2);
    const double lambdaWideCoefficient = lambda1 * lambda1 * lambda2 / (lambda2 * lambda2 - lambda1 * lambda1);
    ZhangUserGonsisExactIlsCache<ZhangUserGonsisIlsSolution> ilsCache;
    // IF_m = lambdaNl*N1 + lambdaWideCoefficient*(N1-N2).
    while (static_cast<int>(active.size()) > options.minimumWideLanePairs)
    {
        if (options.maximumTrials > 0 &&
            static_cast<int>(result.trials.size()) >= options.maximumTrials)
            return reject("TRIAL_BUDGET_EXHAUSTED");
        ZhangUserGonsisTrial trial;
        trial.satelliteIndices = active;
        trial.referenceIndex = zhangUserGonsisDetail::reference(active, satellites, referenceHint);
        trial.referenceSatellite = satellites[trial.referenceIndex].satellite;
        trial.referenceIdentity = satellites[trial.referenceIndex].identity;
        std::vector<VectorXd> passedRows;
        std::vector<double> passedIntegers;
        for (int index : active)
        {
            if (index == trial.referenceIndex) continue;
            VectorXd row = VectorXd::Zero(dimension);
            row(2 * index) = 1;
            row(2 * index + 1) = -1;
            row(2 * trial.referenceIndex) = -1;
            row(2 * trial.referenceIndex + 1) = 1;
            const double floating = row.dot(frozenCycleMean);
            const double integer = zhangUserGonsisRound(floating);
            if (!std::isfinite(floating) || !std::isfinite(integer) || std::abs(integer) > 9007199254740991.0)
                continue;
            trial.wideLaneCandidateSatelliteIndices.push_back(index);
            trial.wideLaneCandidateResidualsCycles.push_back(floating - integer);
            if (options.replayTrial)
            {
                trial.wideLaneCandidateMeansCycles.push_back(floating);
                trial.wideLaneCandidateIntegers.push_back(integer);
            }
            // Inclusive source threshold. No sigma or bootstrap-success gate.
            if (std::abs(floating - integer) <= options.wideLaneResidualCycles)
            {
                passedRows.push_back(std::move(row));
                passedIntegers.push_back(integer);
                trial.wideLaneSatelliteIndices.push_back(index);
            }
        }
        const int rank = passedRows.size();
        trial.requestedL1Rank = rank;
        MatrixXd wideLaneRows(rank, dimension);
        VectorXd wideLaneIntegers(rank);
        MatrixXd l1Rows = MatrixXd::Zero(rank, dimension);
        for (int row = 0; row < rank; ++row)
        {
            wideLaneRows.row(row) = passedRows[row].transpose();
            wideLaneIntegers(row) = passedIntegers[row];
            l1Rows(row, 2 * trial.wideLaneSatelliteIndices[row]) = 1;
            l1Rows(row, 2 * trial.referenceIndex) = -1;
        }
        VectorXd conditionedMean;
        MatrixXd conditionedCovariance;
        VectorXd l1Float;
        MatrixXd l1Covariance;
        const MatrixXd ifRows = lambdaNl * l1Rows + lambdaWideCoefficient * wideLaneRows;
        ZhangUserGonsisIlsSolution solution;
        if (rank < options.minimumWideLanePairs)
            trial.reason = "INSUFFICIENT_PASSED_WIDE_LANE_PAIRS";
        else if (!zhangUserGonsisDetail::conditionWideLane(
            frozenCycleMean, frozenCovariance, wideLaneRows, wideLaneIntegers,
            options.wideLaneConstraintVarianceCycles2, conditionedMean, conditionedCovariance))
            trial.reason = "WIDE_LANE_CONDITIONING_FAILED";
        else
        {
            const MatrixXd searchRows = ifRows / lambdaNl;
            l1Float = searchRows * conditionedMean - (lambdaWideCoefficient / lambdaNl) * wideLaneIntegers;
            l1Covariance = searchRows * conditionedCovariance * searchRows.transpose();
            l1Covariance = (0.5 * (l1Covariance + l1Covariance.transpose())).eval();
            if (!l1Float.allFinite() || !zhangUserGonsisDetail::positiveDefinite(l1Covariance))
                trial.reason = "L1_MOMENTS_INVALID";
            else
            {
                const ZhangUserGonsisIlsKeyView key{l1Float, l1Covariance, wideLaneRows,
                    wideLaneIntegers, l1Rows, ifRows};
                const auto* reused = options.reuseIdenticalIls ? ilsCache.find(key) : nullptr;
                if (reused)
                {
                    solution = *reused;
                    trial.ilsCacheReused = true;
                    ++result.ilsCacheHits;
                }
                else
                {
                    solution = solver(l1Float, l1Covariance);
                    ++result.ilsCalls;
                    if (options.reuseIdenticalIls) ilsCache.remember(key, solution);
                }
                trial.ilsComplete = solution.complete;
                trial.solverReason = solution.reason;
                if (!solution.complete)
                    trial.reason = "L1_ILS_INCOMPLETE";
                else if (solution.integers.size() != rank || !solution.integers.allFinite() ||
                    (solution.integers.array().abs() > 9007199254740991.0).any() ||
                    (solution.integers.array() != solution.integers.array().floor()).any())
                    trial.reason = "L1_INTEGER_CONTRACT_INVALID";
                else if (!std::isfinite(solution.bestSquaredCost) || !std::isfinite(solution.secondSquaredCost) ||
                    solution.bestSquaredCost < 0 || solution.secondSquaredCost < solution.bestSquaredCost)
                    trial.reason = "L1_TWO_CANDIDATE_COSTS_INVALID";
                else
                {
                    // The source defines zero-best-cost ratio as zero.
                    trial.ratio = solution.bestSquaredCost > 0
                        ? solution.secondSquaredCost / solution.bestSquaredCost : 0;
                    trial.reason = trial.ratio > options.ratioThreshold
                        ? "ACCEPTED" : "L1_RATIO_REJECTED";
                }
            }
        }
        auto replay = [&]()
        {
            if (options.replayTrial)
                options.replayTrial(trial, {frozenCycleMean, frozenCycleCovariance, frozenCovariance,
                    wideLaneRows, wideLaneIntegers, conditionedMean, conditionedCovariance,
                    l1Rows, ifRows, l1Float, l1Covariance, solution});
        };
        if (trial.reason == "ACCEPTED")
        {
            replay();
            result.accepted = true;
            result.reason = result.trials.empty() ? "FULL_SUBSET_ACCEPTED" : "PARTIAL_SUBSET_ACCEPTED";
            result.referenceIndex = trial.referenceIndex;
            result.referenceSatellite = trial.referenceSatellite;
            result.referenceIdentity = trial.referenceIdentity;
            result.acceptedSatelliteIndices.push_back(trial.referenceIndex);
            result.acceptedSatelliteIndices.insert(result.acceptedSatelliteIndices.end(),
                trial.wideLaneSatelliteIndices.begin(), trial.wideLaneSatelliteIndices.end());
            result.wideLaneRows = std::move(wideLaneRows);
            result.wideLaneIntegers = std::move(wideLaneIntegers);
            result.l1Rows = std::move(l1Rows);
            result.l1Integers = std::move(solution.integers);
            result.l1FloatMean = std::move(l1Float);
            result.l1FloatCovariance = std::move(l1Covariance);
            result.wlConditionedMean = std::move(conditionedMean);
            result.wlConditionedCovariance = std::move(conditionedCovariance);
            result.constraintRows.resize(2 * rank, dimension);
            result.constraintRows.topRows(rank) = result.wideLaneRows;
            result.constraintRows.bottomRows(rank) = ifRows;
            result.constraintValues.resize(2 * rank);
            result.constraintValues.head(rank) = result.wideLaneIntegers;
            result.constraintValues.tail(rank) = lambdaNl * result.l1Integers + lambdaWideCoefficient * result.wideLaneIntegers;
            result.constraintVariances.resize(2 * rank);
            result.constraintVariances.head(rank).setConstant(options.wideLaneConstraintVarianceCycles2);
            result.constraintVariances.tail(rank).setConstant(options.ifConstraintVarianceMetres2);
            result.trials.push_back(std::move(trial));
            return result;
        }
        if (static_cast<int>(active.size()) <= options.minimumWideLanePairs + 1)
        {
            replay();
            result.trials.push_back(std::move(trial));
            return reject("ALL_SUBSETS_REJECTED");
        }
        // Remove once, then rebuild the reference and both combinations from
        // frozenCycleMean/frozenCovariance in the next iteration.
        // After a real L1 target exists, deleting a satellite which failed WL
        // merely repeats the same effective ILS. Select from the ACTUAL failed
        // target endpoints, reference included. If WL itself has insufficient
        // pairs, no L1 target exists and the legacy all-active pool is explicit.
        std::vector<int> removalPool = active;
        if (options.targetFailedL1SubsetForPar && rank >= options.minimumWideLanePairs)
        {
            removalPool.clear();
            for (const int index : active)
                if (index == trial.referenceIndex || std::find(trial.wideLaneSatelliteIndices.begin(),
                    trial.wideLaneSatelliteIndices.end(), index) != trial.wideLaneSatelliteIndices.end())
                    removalPool.push_back(index);
            trial.removalPool = "ACTUAL_FAILED_WL_L1_ENDPOINTS";
        }
        else trial.removalPool = rank >= options.minimumWideLanePairs
            ? "LEGACY_ALL_ACTIVE_L1_TARGET_IGNORED" : "ALL_ACTIVE_NO_FAILED_L1_TARGET";
        const auto worst = std::max_element(removalPool.begin(), removalPool.end(), [&](int left, int right)
        {
            return satellites[left].metricVarianceMetres2 < satellites[right].metricVarianceMetres2;
        });
        trial.removedSatelliteIndex = *worst;
        result.removedSatelliteIndices.push_back(*worst);
        active.erase(std::find(active.begin(), active.end(), *worst));
        replay();
        result.trials.push_back(std::move(trial));
    }
    return reject("ALL_SUBSETS_REJECTED");
}
