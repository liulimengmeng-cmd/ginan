#include "common/zhangUserGonsis.hpp"

#include <iostream>
#include <stdexcept>

namespace
{
void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

void close(const MatrixXd& actual, const MatrixXd& expected, double tolerance, const std::string& message)
{
    require(actual.rows() == expected.rows() && actual.cols() == expected.cols(), message + " dimensions");
    require(actual.allFinite() && expected.allFinite() &&
        (actual - expected).norm() <= tolerance * std::max(1.0, expected.norm()), message);
}

struct Fixture
{
    VectorXd mean;
    MatrixXd covariance;
    std::vector<ZhangUserGonsisSatellite> satellites;
};

Fixture fixture(int count)
{
    Fixture input;
    input.mean.resize(2 * count);
    input.satellites.resize(count);
    // Independent latent [L1_i, WL_i] variables, including L1-WL covariance.
    // Original N2_i = L1_i-WL_i; this makes cross-covariance nontrivial.
    MatrixXd latent = MatrixXd::Zero(2 * count, 2 * count);
    MatrixXd transform = MatrixXd::Zero(2 * count, 2 * count);
    for (int i = 0; i < count; ++i)
    {
        input.mean(2 * i) = 10 + 0.03 * i;
        input.mean(2 * i + 1) = input.mean(2 * i) - (2 + i % 3 + 0.02 * i);
        latent(i, i) = 0.02 + 0.001 * i;
        latent(count + i, count + i) = 0.1 + 0.002 * i;
        latent(i, count + i) = latent(count + i, i) = 0.01;
        transform(2 * i, i) = 1;
        transform(2 * i + 1, i) = 1;
        transform(2 * i + 1, count + i) = -1;
        auto& satellite = input.satellites[i];
        satellite.satellite = "G" + std::to_string(i + 1);
        satellite.identity = satellite.satellite + ":arc1:state1:product1";
        satellite.elevationRadians = 1.2 - 0.1 * i;
        satellite.effectiveSamples = 100;
        const double lambda1 = ZhangUserGonsisOptions{}.wavelength1Metres;
        const double lambda2 = ZhangUserGonsisOptions{}.wavelength2Metres;
        satellite.metricVarianceMetres2 = lambda1 * lambda1 * latent(i, i) +
            lambda2 * lambda2 * (latent(i, i) + latent(count + i, count + i) - 2 * latent(i, count + i));
    }
    input.covariance = transform * latent * transform.transpose();
    return input;
}

ZhangUserGonsisIlsSolution acceptedStub(const VectorXd& mean, const MatrixXd&)
{
    ZhangUserGonsisIlsSolution solution;
    solution.complete = true;
    solution.integers = mean.unaryExpr([](double value) { return std::floor(value + 0.5); });
    solution.bestSquaredCost = 1;
    solution.secondSquaredCost = 3;
    solution.reason = "TEST_COMPLETE_TWO_CANDIDATES";
    return solution;
}

// Independent bounded integer enumeration. Completeness is certified by a
// lower bound on every point OUTSIDE the box: r'Q^-1*r >= r_i^2/Q_ii.
// This test oracle contains no LAMBDA reduction or candidate reuse.
ZhangUserGonsisIlsSolution enumerate(const VectorXd& mean, const MatrixXd& covariance)
{
    const int radius = 3;
    const int rank = mean.size();
    Eigen::LDLT<MatrixXd> solve(covariance);
    require(solve.info() == Eigen::Success && solve.isPositive(), "enumeration covariance invalid");
    const VectorXd centre = mean.unaryExpr([](double value) { return std::floor(value + 0.5); });
    VectorXd candidate(rank);
    VectorXd best;
    double first = std::numeric_limits<double>::infinity();
    double second = std::numeric_limits<double>::infinity();
    auto visit = [&](auto&& self, int row) -> void
    {
        if (row < rank)
        {
            for (int offset = -radius; offset <= radius; ++offset)
            {
                candidate(row) = centre(row) + offset;
                self(self, row + 1);
            }
            return;
        }
        const VectorXd residual = mean - candidate;
        const double cost = residual.dot(solve.solve(residual));
        if (cost < first)
        {
            second = first;
            first = cost;
            best = candidate;
        }
        else if (cost < second) second = cost;
    };
    visit(visit, 0);
    double outsideLowerBound = std::numeric_limits<double>::infinity();
    for (int row = 0; row < rank; ++row)
    {
        const double low = centre(row) - radius - 1 - mean(row);
        const double high = centre(row) + radius + 1 - mean(row);
        outsideLowerBound = std::min(outsideLowerBound,
            std::min(low * low, high * high) / covariance(row, row));
    }
    ZhangUserGonsisIlsSolution solution;
    solution.complete = outsideLowerBound > second;
    solution.integers = best;
    solution.bestSquaredCost = first;
    solution.secondSquaredCost = second;
    solution.reason = solution.complete ? "EXHAUSTIVE_BOX_WITH_OUTSIDE_BOUND" : "ENUMERATION_BOUND_INSUFFICIENT";
    return solution;
}

struct OracleMoments
{
    MatrixXd wideLaneRows;
    VectorXd wideLaneIntegers;
    MatrixXd l1Rows;
    VectorXd wlMean;
    MatrixXd wlCovariance;
    VectorXd l1Mean;
    MatrixXd l1Covariance;
};

OracleMoments moments(
    const Fixture& input, const std::vector<int>& rovers, int reference,
    const ZhangUserGonsisOptions& options)
{
    OracleMoments oracle;
    const int rank = rovers.size();
    const int dimension = input.mean.size();
    oracle.wideLaneRows = MatrixXd::Zero(rank, dimension);
    oracle.wideLaneIntegers.resize(rank);
    oracle.l1Rows = MatrixXd::Zero(rank, dimension);
    for (int row = 0; row < rank; ++row)
    {
        oracle.wideLaneRows(row, 2 * rovers[row]) = 1;
        oracle.wideLaneRows(row, 2 * rovers[row] + 1) = -1;
        oracle.wideLaneRows(row, 2 * reference) = -1;
        oracle.wideLaneRows(row, 2 * reference + 1) = 1;
        oracle.l1Rows(row, 2 * rovers[row]) = 1;
        oracle.l1Rows(row, 2 * reference) = -1;
        oracle.wideLaneIntegers(row) = std::floor(oracle.wideLaneRows.row(row).dot(input.mean) + 0.5);
    }
    const MatrixXd cross = input.covariance * oracle.wideLaneRows.transpose();
    MatrixXd innovation = oracle.wideLaneRows * cross;
    innovation.diagonal().array() += options.wideLaneConstraintVarianceCycles2;
    // Deliberately use the direct Gaussian block formula as an independent
    // oracle for the production LLT/Joseph implementation.
    const MatrixXd gain = cross * innovation.inverse();
    oracle.wlMean = input.mean + gain * (oracle.wideLaneIntegers - oracle.wideLaneRows * input.mean);
    oracle.wlCovariance = input.covariance - gain * cross.transpose();
    const double f1 = 1 / options.wavelength1Metres;
    const double f2 = 1 / options.wavelength2Metres;
    const double lambdaNl = 1 / (f1 + f2);
    const double lambdaWide = f2 / (f1 * f1 - f2 * f2);
    MatrixXd ifRows = MatrixXd::Zero(rank, dimension);
    for (int row = 0; row < rank; ++row)
    {
        ifRows(row, 2 * rovers[row]) = f1 / (f1 * f1 - f2 * f2);
        ifRows(row, 2 * rovers[row] + 1) = -f2 / (f1 * f1 - f2 * f2);
        ifRows(row, 2 * reference) = -f1 / (f1 * f1 - f2 * f2);
        ifRows(row, 2 * reference + 1) = f2 / (f1 * f1 - f2 * f2);
    }
    oracle.l1Mean = (ifRows * oracle.wlMean - lambdaWide * oracle.wideLaneIntegers) / lambdaNl;
    oracle.l1Covariance = ifRows * oracle.wlCovariance * ifRows.transpose() / (lambdaNl * lambdaNl);
    return oracle;
}

void rejectedWithoutPublication(const ZhangUserGonsisResult& result)
{
    require(!result.accepted && result.acceptedSatelliteIndices.empty(), "rejection exposed accepted satellites");
    require(result.wideLaneRows.rows() == 0 && result.wideLaneIntegers.size() == 0 &&
        result.l1Rows.rows() == 0 && result.l1Integers.size() == 0 &&
        result.constraintRows.rows() == 0 && result.constraintValues.size() == 0 &&
        result.wlConditionedMean.size() == 0, "rejection exposed a WL-only or private posterior result");
}

void testBoundariesAndNoSigmaGate()
{
    require(zhangUserGonsisRound(-1.5) == -1 && zhangUserGonsisRound(-0.5) == 0 &&
        zhangUserGonsisRound(0.5) == 1, "GONSIS negative rounding changed");
    Fixture input = fixture(5);
    input.mean.setZero();
    const double values[] = {0, 1.25, -1.25, 0.25, -0.25};
    for (int i = 0; i < 5; ++i) input.mean(2 * i) = values[i];
    // Deliberately enormous sigmas: UDUC uses no WL sigma/bootstrap gate.
    input.covariance = MatrixXd::Identity(10, 10) * 10000;
    int calls = 0;
    auto solver = [&](const VectorXd& mean, const MatrixXd& covariance)
    {
        ++calls;
        return acceptedStub(mean, covariance);
    };
    const auto accepted = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, solver);
    require(accepted.accepted && calls == 1 && accepted.wideLaneIntegers.size() == 4,
        "inclusive WL boundary or sigma-free UDUC gate failed");
    VectorXd expected(4);
    expected << 1, -1, 0, 0;
    close(accepted.wideLaneIntegers, expected, 0, "negative WL boundary integers");
    input.mean(6) = std::nextafter(0.25, 1.0);
    calls = 0;
    const auto rejected = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, solver);
    rejectedWithoutPublication(rejected);
    require(calls == 0 && rejected.trials[0].requestedL1Rank == 3 &&
        rejected.trials[0].reason == "INSUFFICIENT_PASSED_WIDE_LANE_PAIRS", "WL boundary overrun accepted");
}

void testCorrelatedEnumerationAndFullCovariance()
{
    Fixture input = fixture(5);
    MatrixXd latent = MatrixXd::Zero(10, 10);
    MatrixXd transform = MatrixXd::Zero(10, 10);
    const double l1[] = {10, 10.45, 10.52, 10.13, 9.93};
    for (int i = 0; i < 5; ++i)
    {
        latent(i, i) = i == 0 ? 0.0001 : 0.02;
        latent(5 + i, 5 + i) = 0.1;
        latent(i, 5 + i) = latent(5 + i, i) = 0.0001;
        transform(2 * i, i) = 1;
        transform(2 * i + 1, i) = 1;
        transform(2 * i + 1, 5 + i) = -1;
        input.mean(2 * i) = l1[i];
        input.mean(2 * i + 1) = l1[i] - (2 + i + 0.01 * i);
    }
    latent(1, 2) = latent(2, 1) = 0.0194;
    input.covariance = transform * latent * transform.transpose();
    const VectorXd frozenMean = input.mean;
    const MatrixXd frozenCovariance = input.covariance;
    ZhangUserGonsisOptions options;
    options.ratioThreshold = 1.01; // Numerical test of the complete correlated search.
    const auto oracle = moments(input, {1, 2, 3, 4}, 0, options);
    int calls = 0;
    auto solver = [&](const VectorXd& mean, const MatrixXd& covariance)
    {
        ++calls;
        close(mean, oracle.l1Mean, 1e-10, "L1 affine IF/WL mean differs from independent oracle");
        close(covariance, oracle.l1Covariance, 1e-8, "L1 covariance lost full cross correlations");
        return enumerate(mean, covariance);
    };
    const auto result = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, solver, options);
    require(result.accepted && calls == 1, "correlated complete enumeration not accepted");
    const auto enumerated = enumerate(oracle.l1Mean, oracle.l1Covariance);
    require(enumerated.complete, "bounded oracle failed its outside-box certificate");
    close(result.l1Integers, enumerated.integers, 0, "ILS candidate differs from bounded enumeration");
    require((result.l1Integers.array() != oracle.l1Mean.array().round()).any(),
        "correlated fixture did not exercise an integer different from coordinatewise rounding");
    close(result.wlConditionedMean, oracle.wlMean, 1e-10, "WL full posterior mean mismatch");
    close(result.wlConditionedCovariance, oracle.wlCovariance, 1e-8, "WL full posterior covariance mismatch");
    close(input.mean, frozenMean, 0, "frozen input mean mutated");
    close(input.covariance, frozenCovariance, 0, "frozen input covariance mutated");
    close(result.constraintRows.topRows(4), result.wideLaneRows, 0, "output WL rows mismatch");
    require((result.constraintVariances.head(4).array() == 1e-8).all() &&
        (result.constraintVariances.tail(4).array() == 0.003 * 0.003).all(), "source soft constraint variances changed");
    const double ln = options.wavelength1Metres * options.wavelength2Metres /
        (options.wavelength1Metres + options.wavelength2Metres);
    const double lw = options.wavelength1Metres * options.wavelength1Metres * options.wavelength2Metres /
        (options.wavelength2Metres * options.wavelength2Metres - options.wavelength1Metres * options.wavelength1Metres);
    close(result.constraintRows.bottomRows(4), ln * result.l1Rows + lw * result.wideLaneRows,
        1e-14, "physical IF/WL/L1 identity mismatch");
    close(result.constraintValues.tail(4), ln * result.l1Integers + lw * result.wideLaneIntegers,
        1e-14, "soft IF output target mismatch");
}

void testParRemovalsReferenceChangeAndFreshSnapshot()
{
    Fixture input = fixture(7);
    input.satellites[0].metricVarianceMetres2 = 90;
    input.satellites[1].metricVarianceMetres2 = 80;
    input.satellites[2].metricVarianceMetres2 = 70;
    int calls = 0;
    auto solver = [&](const VectorXd& mean, const MatrixXd& covariance)
    {
        const int reference = calls;
        std::vector<int> rovers;
        for (int i = reference + 1; i < 7; ++i) rovers.push_back(i);
        const auto oracle = moments(input, rovers, reference, {});
        close(mean, oracle.l1Mean, 1e-10, "PAR reused a failed conditioned mean");
        close(covariance, oracle.l1Covariance, 1e-8, "PAR reused a failed conditioned covariance");
        auto solution = acceptedStub(mean, covariance);
        ++calls;
        if (calls < 3) solution.secondSquaredCost = 1.5;
        return solution;
    };
    const auto retry = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, solver);
    require(retry.accepted && calls == 3 && retry.reason == "PARTIAL_SUBSET_ACCEPTED", "PAR retry failed");
    require(retry.removedSatelliteIndices == std::vector<int>({0, 1}), "PAR failed to pop distinct worst satellites");
    require(retry.trials.size() == 3 && retry.trials[0].referenceIndex == 0 &&
        retry.trials[1].referenceIndex == 1 && retry.trials[2].referenceIndex == 2,
        "PAR did not rebuild references after removals");
    require(retry.trials[0].satelliteIndices.size() == 7 && retry.trials[1].satelliteIndices.size() == 6 &&
        retry.trials[2].satelliteIndices.size() == 5, "PAR subset did not actually shrink");
    auto freshMetadata = input.satellites;
    freshMetadata[0].eligible = false;
    freshMetadata[1].eligible = false;
    const auto fresh = zhangResolveUserGonsis(input.mean, input.covariance, freshMetadata, acceptedStub);
    require(fresh.accepted && fresh.referenceIndex == retry.referenceIndex, "fresh surviving subset rejected");
    close(retry.l1FloatMean, fresh.l1FloatMean, 0, "full-to-retry differs from fresh frozen subset mean");
    close(retry.l1FloatCovariance, fresh.l1FloatCovariance, 0, "full-to-retry differs from fresh frozen subset covariance");
    close(retry.constraintRows, fresh.constraintRows, 0, "retry rows differ from fresh physical subset rows");
    require(retry.wideLaneRows.cols() == input.mean.size() &&
        retry.wideLaneRows.leftCols(4).isZero(0), "PAR changed original-column physical coordinate mapping");
}

void testDeterministicOmittedReferenceAndSearchedPositiveDefiniteness()
{
    const Fixture full = fixture(5);
    Fixture datum = full;
    MatrixXd transform = MatrixXd::Identity(10, 10);
    transform.topRows(2).setZero();
    for (int satellite = 1; satellite < 5; ++satellite)
    {
        transform(2 * satellite, 0) = -1;
        transform(2 * satellite + 1, 1) = -1;
    }
    datum.mean = transform * full.mean;
    datum.covariance = transform * full.covariance * transform.transpose();
    require(datum.mean.head(2).isZero(0) && datum.covariance.topRows(2).isZero(0) &&
        datum.covariance.leftCols(2).isZero(0), "datum fixture omitted reference is not deterministic zero");
    const auto fullResult = zhangResolveUserGonsis(full.mean, full.covariance, full.satellites, acceptedStub);
    const auto datumResult = zhangResolveUserGonsis(datum.mean, datum.covariance, datum.satellites, acceptedStub);
    require(fullResult.accepted && datumResult.accepted && datumResult.referenceIndex == 0 &&
        datumResult.acceptedSatelliteIndices.size() == 5, "legal deterministic reference PSD posterior rejected");
    close(datumResult.l1FloatMean, fullResult.l1FloatMean, 1e-11, "datum transformation changed L1 search mean");
    close(datumResult.l1FloatCovariance, fullResult.l1FloatCovariance, 1e-11,
        "datum transformation changed complete searched L1 covariance");
    close(datumResult.wideLaneIntegers, fullResult.wideLaneIntegers, 0, "datum transformation changed WL integers");
    close(datumResult.l1Integers, fullResult.l1Integers, 0, "datum transformation changed L1 integers");
    close(datumResult.wlConditionedMean, transform * fullResult.wlConditionedMean, 1e-11,
        "datum transformation changed soft WL posterior mean");
    close(datumResult.wlConditionedCovariance,
        transform * fullResult.wlConditionedCovariance * transform.transpose(), 1e-11,
        "datum transformation changed soft WL posterior covariance");
    const auto oracle = moments(datum, {1, 2, 3, 4}, 0, {});
    close(datumResult.l1FloatMean, oracle.l1Mean, 1e-11, "PSD reference mean differs from independent Gaussian oracle");
    close(datumResult.l1FloatCovariance, oracle.l1Covariance, 1e-8,
        "PSD reference covariance differs from independent Gaussian oracle");

    int calls = 0;
    auto counted = [&](const VectorXd& mean, const MatrixXd& covariance)
    {
        ++calls;
        return acceptedStub(mean, covariance);
    };
    datum.covariance.setZero();
    const auto unobservable = zhangResolveUserGonsis(datum.mean, datum.covariance, datum.satellites, counted);
    rejectedWithoutPublication(unobservable);
    require(calls == 0 && unobservable.trials.size() == 1 &&
        unobservable.trials[0].reason == "L1_MOMENTS_INVALID", "singular searched L1 covariance accepted");
    datum.covariance(0, 0) = -1e-4;
    const auto negative = zhangResolveUserGonsis(datum.mean, datum.covariance, datum.satellites, counted);
    rejectedWithoutPublication(negative);
    require(negative.reason == "COVARIANCE_NOT_POSITIVE_SEMIDEFINITE" && negative.trials.empty(),
        "negative ambient covariance accepted");
}

void testStrictRatioAndIncompleteContracts()
{
    Fixture input = fixture(5);
    auto equalRatio = [](const VectorXd& mean, const MatrixXd& covariance)
    {
        auto solution = acceptedStub(mean, covariance);
        solution.bestSquaredCost = 2;
        solution.secondSquaredCost = 4;
        return solution;
    };
    auto equality = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, equalRatio);
    rejectedWithoutPublication(equality);
    require(equality.trials[0].ratio == 2 && equality.trials[0].reason == "L1_RATIO_REJECTED", "ratio equality accepted");
    for (int failure = 0; failure < 6; ++failure)
    {
        auto invalid = [&](const VectorXd& mean, const MatrixXd& covariance)
        {
            auto solution = acceptedStub(mean, covariance);
            if (failure == 0) { solution.complete = false; solution.reason = "NODE_BUDGET_EXHAUSTED"; }
            if (failure == 1) solution.integers.conservativeResize(mean.size() - 1);
            if (failure == 2) solution.integers(0) = 0.5;
            if (failure == 3) solution.secondSquaredCost = 0.5;
            if (failure == 4) solution.bestSquaredCost = -1;
            if (failure == 5) solution.secondSquaredCost = std::numeric_limits<double>::infinity();
            return solution;
        };
        auto rejected = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, invalid);
        rejectedWithoutPublication(rejected);
        require(rejected.trials.size() == 1, "invalid full-rank contract spawned extra trial below minimum support");
        if (failure == 0) require(rejected.trials[0].reason == "L1_ILS_INCOMPLETE" &&
            rejected.trials[0].solverReason == "NODE_BUDGET_EXHAUSTED", "incomplete ILS not rejected");
    }
    auto zeroCost = [](const VectorXd& mean, const MatrixXd& covariance)
    {
        auto solution = acceptedStub(mean, covariance);
        solution.bestSquaredCost = 0;
        return solution;
    };
    const auto zero = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, zeroCost);
    rejectedWithoutPublication(zero);
    require(zero.trials[0].ratio == 0, "source zero-cost ratio semantics changed");
}

void testBudgetsAndInfeasibleWideLane()
{
    Fixture input = fixture(7);
    int calls = 0;
    auto failing = [&](const VectorXd& mean, const MatrixXd& covariance)
    {
        ++calls;
        auto solution = acceptedStub(mean, covariance);
        solution.complete = false;
        return solution;
    };
    ZhangUserGonsisOptions options;
    options.maximumTrials = 1;
    auto budget = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, failing, options);
    rejectedWithoutPublication(budget);
    require(budget.reason == "TRIAL_BUDGET_EXHAUSTED" && calls == 1 && budget.trials.size() == 1,
        "trial cap failed to reject incomplete search");
    calls = 0;
    auto exhausted = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, failing);
    rejectedWithoutPublication(exhausted);
    require(calls == 3 && exhausted.trials.size() == 3 && exhausted.removedSatelliteIndices.size() == 2,
        "zero cap did not terminate at minimum support");
    Fixture impossible = fixture(5);
    impossible.mean.setZero();
    for (int i = 1; i < 5; ++i) impossible.mean(2 * i) = 0.3;
    calls = 0;
    auto infeasible = zhangResolveUserGonsis(impossible.mean, impossible.covariance, impossible.satellites, failing);
    rejectedWithoutPublication(infeasible);
    require(calls == 0 && infeasible.trials[0].requestedL1Rank == 0, "infeasible WL reached ILS or publication");
    impossible.satellites[4].eligible = false;
    auto insufficient = zhangResolveUserGonsis(impossible.mean, impossible.covariance, impossible.satellites, acceptedStub);
    rejectedWithoutPublication(insufficient);
    require(insufficient.reason == "INSUFFICIENT_ELIGIBLE_SATELLITES" && insufficient.trials.empty(), "minimum four SD pairs weakened");
}

void testReferenceIdentityAndSampleThreshold()
{
    Fixture input = fixture(6);
    ZhangUserGonsisReferenceHint hint;
    hint.satellite = input.satellites[1].satellite;
    hint.identity = input.satellites[1].identity;
    hint.previousPublishedFix = true;
    input.satellites[0].effectiveSamples = 200;
    auto solve = [&]() { return zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, acceptedStub, {}, hint); };
    require(solve().referenceIndex == 1, "eligible prior reference not retained at >30 degrees");
    input.satellites[0].effectiveSamples = 201;
    require(solve().referenceIndex == 0, "higher elevation with >200 effective samples did not replace prior reference");
    input.satellites[0].effectiveSamples = 200;
    const std::string identity = input.satellites[1].identity;
    for (const std::string& changed : {"arc2", "state2", "product2"})
    {
        input.satellites[1].identity = identity + ":" + changed;
        require(solve().referenceIndex == 0, "changed physical/state/product identity retained old reference");
    }
    input.satellites[1].identity = identity;
    input.satellites[1].elevationRadians = std::acos(-1.0) / 6;
    require(solve().referenceIndex == 0, "prior reference at exactly 30 degrees retained");
    input.satellites[1].elevationRadians = 1.1;
    hint.previousPublishedFix = false;
    require(solve().referenceIndex == 0, "prior reference retained without previous published fix");
    hint.previousPublishedFix = true;
    input.satellites[1].eligible = false;
    require(solve().referenceIndex == 0, "ineligible prior reference retained");
}

void testInvalidInputs()
{
    Fixture input = fixture(5);
    auto invalid = input;
    invalid.covariance(0, 0) = -1;
    auto rejected = zhangResolveUserGonsis(invalid.mean, invalid.covariance, invalid.satellites, acceptedStub);
    rejectedWithoutPublication(rejected);
    require(rejected.reason == "COVARIANCE_NOT_POSITIVE_SEMIDEFINITE", "negative covariance accepted");
    invalid = input;
    invalid.covariance(0, 1) += 0.1;
    rejected = zhangResolveUserGonsis(invalid.mean, invalid.covariance, invalid.satellites, acceptedStub);
    require(!rejected.accepted && rejected.reason == "COVARIANCE_NOT_SYMMETRIC", "asymmetric covariance accepted");
    invalid = input;
    invalid.satellites[1].satellite = invalid.satellites[0].satellite;
    rejected = zhangResolveUserGonsis(invalid.mean, invalid.covariance, invalid.satellites, acceptedStub);
    require(!rejected.accepted && rejected.reason == "SATELLITE_METADATA_INVALID", "duplicate satellite made dependent SD pairs");
    ZhangUserGonsisOptions options;
    options.minimumWideLanePairs = 3;
    rejected = zhangResolveUserGonsis(input.mean, input.covariance, input.satellites, acceptedStub, options);
    require(!rejected.accepted && rejected.reason == "OPTIONS_INVALID", "support minimum relaxed below source minimum");
}
}

int main()
{
    try
    {
        testBoundariesAndNoSigmaGate();
        testCorrelatedEnumerationAndFullCovariance();
        testParRemovalsReferenceChangeAndFreshSnapshot();
        testDeterministicOmittedReferenceAndSearchedPositiveDefiniteness();
        testStrictRatioAndIncompleteContracts();
        testBudgetsAndInfeasibleWideLane();
        testReferenceIdentityAndSampleThreshold();
        testInvalidInputs();
        std::cout << "PASS GONSIS UDUC: frozen full covariance, independently bounded correlated ILS, "
            "inclusive WL/negative rounding, complete-rank strict ratio, distinct PAR removals, "
            "reference identity/deterministic PSD datum, strict searched PD, and rejection without WL-only publication\n";
    }
    catch (const std::exception& exception)
    {
        std::cerr << "FAIL GONSIS UDUC contract: " << exception.what() << '\n';
        return 1;
    }
    return 0;
}
