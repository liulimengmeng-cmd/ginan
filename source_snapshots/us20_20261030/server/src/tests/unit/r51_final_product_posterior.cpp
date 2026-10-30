#include "common/zhangFinalProductPosterior.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}

static KFKey keyFor(int i)
{
    KFKey key;
    key.type = i < 2 ? KF::PHASE_BIAS : KF::AMBIGUITY;
    key.Sat = SatSys(E_Sys::GPS, i % 30 + 2);
    key.num = i;
    key.str = "fixture";
    return key;
}

// Independent scalar-conditioning oracle. It uses the closed-form update for
// one scalar random variable at a time, without calling the production helper,
// LDLT, spectral pseudoinverse or any production conditioner.
static void scalarOracle(const Eigen::VectorXd& c, double z,
    Eigen::VectorXd& mean, Eigen::MatrixXd& covariance)
{
    Eigen::VectorXd cross = Eigen::VectorXd::Zero(mean.size());
    double current = 0;
    for (int i = 0; i < mean.size(); ++i)
    {
        current += c(i) * mean(i);
        for (int j = 0; j < mean.size(); ++j) cross(i) += covariance(i, j) * c(j);
    }
    double variance = 0;
    for (int i = 0; i < mean.size(); ++i) variance += c(i) * cross(i);
    require(variance > 0, "oracle scalar variance is not positive");
    const double delta = (z - current) / variance;
    Eigen::VectorXd nextMean(mean.size());
    Eigen::MatrixXd nextCovariance(mean.size(), mean.size());
    for (int i = 0; i < mean.size(); ++i)
    {
        nextMean(i) = mean(i) + cross(i) * delta;
        for (int j = 0; j < mean.size(); ++j)
            nextCovariance(i, j) = covariance(i, j) - cross(i) * cross(j) / variance;
    }
    mean = std::move(nextMean);
    covariance = std::move(nextCovariance);
}

int main()
{
    std::map<KFKey, int> indices;
    for (int i = 0; i < 4; ++i) indices.emplace(keyFor(i), i);
    Eigen::VectorXd sourceMean(4);
    sourceMean << 0.7, -0.3, 1.2, -0.8;
    Eigen::MatrixXd sourceCovariance(4, 4);
    sourceCovariance << 4, 1, 0.6, -0.2,
                        1, 3, 0.4, 0.3,
                        0.6, 0.4, 2, 0.5,
                        -0.2, 0.3, 0.5, 1.5;
    const auto sourceMeanBefore = sourceMean;
    const auto sourceCovarianceBefore = sourceCovariance;
    const ZhangCurrentEpochNumericConstraint mainConstraint{{{keyFor(2), 1}}, 2};
    const ZhangCurrentEpochNumericConstraint eventConstraint{
        {{keyFor(0), 1}, {keyFor(1), -1}}, 1};
    Eigen::VectorXd oracleMean = sourceMean;
    Eigen::MatrixXd oracleCovariance = sourceCovariance;
    Eigen::VectorXd first = Eigen::VectorXd::Zero(4), second = first;
    first(2) = 1;
    second(0) = 1;
    second(1) = -1;
    scalarOracle(first, 2, oracleMean, oracleCovariance);
    scalarOracle(second, 1, oracleMean, oracleCovariance);

    Eigen::VectorXd mean = sourceMean;
    Eigen::MatrixXd covariance = sourceCovariance;
    const auto joint = zhangConditionFinalProductPosterior(indices,
        {mainConstraint, eventConstraint}, mean, covariance);
    if (!joint.valid)
        std::cerr << "joint rejection=" << joint.reason << " rows=" << joint.independentRows
            << " residual=" << joint.maximumConstraintResidual << '\n';
    require(joint.valid && joint.independentRows == 2, "joint current posterior conditioning failed");
    require((mean - oracleMean).cwiseAbs().maxCoeff() < 1e-11,
        "joint mean differs from independent scalar oracle");
    require((covariance - oracleCovariance).cwiseAbs().maxCoeff() < 1e-11,
        "full joint cross covariance differs from scalar oracle");
    require((sourceMean - sourceMeanBefore).norm() == 0 &&
        (sourceCovariance - sourceCovarianceBefore).norm() == 0, "authoritative FLOAT mutated");
    require(std::abs(covariance(0, 3) - sourceCovariance(0, 3)) > 1e-3 &&
        std::abs(covariance(1, 3) - sourceCovariance(1, 3)) > 1e-3,
        "clock/phase/other cross covariance was not propagated");
    const auto onceMean = mean;
    const auto onceCovariance = covariance;
    auto repeated = zhangConditionFinalProductPosterior(indices,
        {mainConstraint, eventConstraint, eventConstraint}, mean, covariance);
    require(repeated.valid && (mean - onceMean).norm() < 1e-10 &&
        (covariance - onceCovariance).norm() < 1e-10,
        "repeated current constraints changed already-conditioned posterior");

    mean = sourceMean;
    covariance = sourceCovariance;
    require(zhangConditionFinalProductPosterior(indices, {mainConstraint}, mean, covariance).valid,
        "main conditioning failed");
    require(zhangConditionFinalProductPosterior(indices, {eventConstraint}, mean, covariance).valid,
        "sequential event conditioning failed");
    require((mean - onceMean).norm() < 1e-10 &&
        (covariance - onceCovariance).norm() < 1e-10,
        "unified current constraints differ from valid sequential conditioning");

    mean = sourceMean;
    covariance = sourceCovariance;
    const ZhangCurrentEpochNumericConstraint dependent{
        {{keyFor(0), 2}, {keyFor(1), -2}}, 2};
    const auto withDependent = zhangConditionFinalProductPosterior(indices,
        {mainConstraint, eventConstraint, dependent}, mean, covariance);
    require(withDependent.valid && withDependent.independentRows == 2 &&
        (mean - onceMean).norm() < 1e-10 && (covariance - onceCovariance).norm() < 1e-10,
        "dependent deterministic current row lost information or applied twice");

    auto rejectWithoutMutation = [&](const std::vector<ZhangCurrentEpochNumericConstraint>& constraints,
        const Eigen::VectorXd& suppliedMean, const Eigen::MatrixXd& suppliedCovariance,
        const char* why)
    {
        auto rejectedMean = suppliedMean;
        auto rejectedCovariance = suppliedCovariance;
        const auto rejected = zhangConditionFinalProductPosterior(indices, constraints,
            rejectedMean, rejectedCovariance);
        require(!rejected.valid, why);
        require((rejectedMean.array() == suppliedMean.array()).all() &&
            (rejectedCovariance.array() == suppliedCovariance.array()).all(),
            "failure partially committed private posterior");
    };
    auto conflict = eventConstraint;
    conflict.value = 2;
    rejectWithoutMutation({eventConstraint, conflict}, sourceMean, sourceCovariance,
        "exact duplicate integer conflict accepted");
    conflict = dependent;
    conflict.value = 3;
    rejectWithoutMutation({eventConstraint, conflict}, sourceMean, sourceCovariance,
        "dependent integer conflict accepted");
    rejectWithoutMutation({{{{keyFor(0), 0}}, 0}}, sourceMean, sourceCovariance,
        "structural zero admitted as integer conditioning evidence");
    rejectWithoutMutation({{{{keyFor(0), std::numeric_limits<double>::infinity()}}, 0}},
        sourceMean, sourceCovariance, "nonfinite numerical constraint accepted");
    auto unavailable = keyFor(8);
    rejectWithoutMutation({{{{unavailable, 1}}, 0}}, sourceMean, sourceCovariance,
        "unavailable KF column assumed zero");
    const ZhangExactInteger inexact = (ZhangExactInteger(1) << 70) + 1;
    rejectWithoutMutation({{{{keyFor(0), 1}}, inexact}}, sourceMean, sourceCovariance,
        "inexact integer RHS converted silently");

    Eigen::MatrixXd indefinite = sourceCovariance;
    indefinite(0, 1) = indefinite(1, 0) = 20;
    rejectWithoutMutation({mainConstraint}, sourceMean, indefinite,
        "non-PSD covariance in unconstrained directions accepted");
    rejectWithoutMutation({}, sourceMean, indefinite,
        "no-constraint path accepted non-PSD FLOAT");
    auto asymmetric = sourceCovariance;
    asymmetric(0, 1) += 0.1;
    rejectWithoutMutation({mainConstraint}, sourceMean, asymmetric,
        "asymmetric source silently reinterpreted");

    // A one-cycle dependency conflict must not disappear inside a tolerance
    // scaled by a billion-cycle RHS.
    Eigen::VectorXd largeMean = sourceMean;
    largeMean(0) = 1000000000;
    largeMean(1) = 0;
    auto largeFirst = eventConstraint;
    largeFirst.value = 1000000000;
    auto largeDependent = dependent;
    largeDependent.value = 2000000001;
    rejectWithoutMutation({largeFirst, largeDependent}, largeMean, sourceCovariance,
        "large RHS concealed one-cycle dependent conflict");

    mean = sourceMean;
    covariance = sourceCovariance;
    const auto trusted = zhangConditionFinalProductPosterior(indices,
        {mainConstraint, eventConstraint}, mean, covariance, true);
    require(trusted.valid && trusted.validationMode ==
        "UPSTREAM_SOURCE_PSD_VALIDATED_CONSTRAINT_AND_CLOSURE_CHECKS" &&
        (mean - onceMean).norm() < 1e-10 && (covariance - onceCovariance).norm() < 1e-10,
        "trusted PSD source mode changed conditioning mathematics");

    // Exercise sparse KF functional projection and more than one covariance
    // closure column block, while retaining all unconstrained cross terms.
    std::map<KFKey, int> largeIndices;
    const int dimensions = 100;
    Eigen::VectorXd largeSourceMean(dimensions), direction(dimensions);
    for (int i = 0; i < dimensions; ++i)
    {
        largeIndices.emplace(keyFor(i), i);
        largeSourceMean(i) = 0.1 * (i % 7 - 3);
        direction(i) = 0.03 * (i % 11 - 5);
    }
    const Eigen::MatrixXd largeSourceCovariance = Eigen::MatrixXd::Identity(dimensions, dimensions) +
        direction * direction.transpose();
    Eigen::VectorXd expectedMean = largeSourceMean;
    Eigen::MatrixXd expectedCovariance = largeSourceCovariance;
    Eigen::VectorXd functional = Eigen::VectorXd::Zero(dimensions);
    functional(2) = 1;
    scalarOracle(functional, 2, expectedMean, expectedCovariance);
    functional.setZero(); functional(70) = 1; functional(80) = -1;
    scalarOracle(functional, -1, expectedMean, expectedCovariance);
    functional.setZero(); functional(99) = 1;
    scalarOracle(functional, 4, expectedMean, expectedCovariance);
    auto largeActualMean = largeSourceMean;
    auto largeActualCovariance = largeSourceCovariance;
    const auto blocked = zhangConditionFinalProductPosterior(largeIndices,
        {{{{keyFor(2), 1}}, 2}, {{{keyFor(70), 1}, {keyFor(80), -1}}, -1},
         {{{keyFor(99), 1}}, 4}}, largeActualMean, largeActualCovariance);
    require(blocked.valid && (largeActualMean - expectedMean).cwiseAbs().maxCoeff() < 1e-10 &&
        (largeActualCovariance - expectedCovariance).cwiseAbs().maxCoeff() < 1e-10,
        "sparse/block conditioning differs from independent full-covariance oracle");
    auto aliasedIndices = indices;
    aliasedIndices[keyFor(1)] = aliasedIndices[keyFor(0)];
    mean = sourceMean;
    covariance = sourceCovariance;
    require(!zhangConditionFinalProductPosterior(aliasedIndices, {eventConstraint}, mean, covariance).valid &&
        (mean - sourceMean).norm() == 0 && (covariance - sourceCovariance).norm() == 0,
        "KF index alias cancellation manufactured nonzero numerical evidence");
    std::cout << "PASS final_product_posterior independent_scalar_joint_oracle_and_failures\n";
}
