#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

struct ZhangUserOutputValidation
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    double maximumMeanResidual = 0;
    double maximumMeanTolerance = 0;
    double maximumCovarianceResidual = 0;
    double maximumCovarianceTolerance = 0;
    double maximumSymmetryResidual = 0;
    double symmetryTolerance = 0;
    double minimumEigenvalue = std::numeric_limits<double>::quiet_NaN();
    double psdTolerance = std::numeric_limits<double>::quiet_NaN();
    double roundingFactor = std::numeric_limits<double>::quiet_NaN();
};

// A numerical commit gate for a disposable fixed posterior. This does not
// certify that an accepted integer is physically correct. No input is changed,
// no NIS is evaluated, and neither asymmetry nor a negative eigenvalue is
// repaired. Tolerances are floating-point error budgets, not statistical gates.
inline ZhangUserOutputValidation zhangValidateUserIntegerOutput(
    const Eigen::VectorXd& mean, const Eigen::MatrixXd& covariance,
    const Eigen::MatrixXd& integerRows, const Eigen::VectorXd& integers)
{
    ZhangUserOutputValidation out;
    const Eigen::Index dimension = mean.size();
    if (dimension <= 0 || covariance.rows() != dimension ||
        covariance.cols() != dimension || integerRows.cols() != dimension ||
        integerRows.rows() != integers.size())
    { out.reason = "OUTPUT_POSTERIOR_DIMENSION_MISMATCH"; return out; }
    if (!mean.allFinite() || !covariance.allFinite() ||
        !integerRows.allFinite() || !integers.allFinite())
    { out.reason = "OUTPUT_POSTERIOR_NONFINITE_INPUT"; return out; }

    // n accounts for dot-product lengths and eigensolver dimension; 64 leaves
    // room for the preceding finite conditioning operations. Absolute floors
    // are in the respective supplied state units, never arbitrary jitter.
    out.roundingFactor = 64 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max<Eigen::Index>(1, dimension));
    bool symmetryFailed = false;
    for (Eigen::Index row = 0; row < dimension; ++row)
    for (Eigen::Index column = row + 1; column < dimension; ++column)
    {
        const long double scale = std::max(1.0L,
            std::abs(static_cast<long double>(covariance(row,column))) +
            std::abs(static_cast<long double>(covariance(column,row))));
        const long double tolerance = out.roundingFactor * scale;
        const long double residual = std::abs(
            static_cast<long double>(covariance(row,column)) - covariance(column,row));
        if (residual > std::numeric_limits<double>::max() ||
            tolerance > std::numeric_limits<double>::max())
        { out.reason = "OUTPUT_POSTERIOR_SCALE_OVERFLOW"; return out; }
        out.maximumSymmetryResidual = std::max(out.maximumSymmetryResidual,
            static_cast<double>(residual));
        out.symmetryTolerance = std::max(out.symmetryTolerance,
            static_cast<double>(tolerance));
        symmetryFailed |= residual > tolerance;
    }
    if (symmetryFailed)
    { out.reason = "OUTPUT_POSTERIOR_ASYMMETRIC"; return out; }

    // Divide before symmetrising so finite very large entries cannot overflow
    // merely through P+P'. Symmetry above is a separate rejection gate.
    const double covarianceScale = std::max(1.0, covariance.cwiseAbs().maxCoeff());
    const Eigen::MatrixXd scaled = covariance / covarianceScale;
    const Eigen::MatrixXd symmetric = (0.5 * scaled + 0.5 * scaled.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetric);
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
    { out.reason = "OUTPUT_POSTERIOR_EIGEN_FAILED"; return out; }
    const double scaledTolerance = out.roundingFactor * std::max(1.0,
        eigen.eigenvalues().cwiseAbs().maxCoeff());
    const long double minimumEigenvalue =
        static_cast<long double>(eigen.eigenvalues().minCoeff()) * covarianceScale;
    const long double psdTolerance = static_cast<long double>(scaledTolerance) * covarianceScale;
    if (std::abs(minimumEigenvalue) > std::numeric_limits<double>::max() ||
        psdTolerance > std::numeric_limits<double>::max())
    { out.reason = "OUTPUT_POSTERIOR_SCALE_OVERFLOW"; return out; }
    out.minimumEigenvalue = static_cast<double>(minimumEigenvalue);
    out.psdTolerance = static_cast<double>(psdTolerance);
    if (eigen.eigenvalues().minCoeff() < -scaledTolerance)
    { out.reason = "OUTPUT_POSTERIOR_NOT_PSD"; return out; }

    // Empty row sets are defined without issuing an M=0 BLAS operation.
    if (integerRows.rows() == 0)
    { out.valid = true; out.reason = "NUMERICAL_OUTPUT_VALID_NO_INTEGER_ROWS"; return out; }
    const Eigen::VectorXd meanResidual = integerRows * mean - integers;
    const Eigen::MatrixXd covarianceResidual = integerRows * covariance;
    if (!meanResidual.allFinite() || !covarianceResidual.allFinite())
    { out.reason = "OUTPUT_POSTERIOR_NONFINITE_CLOSURE"; return out; }
    bool meanFailed = false;
    bool covarianceFailed = false;
    for (Eigen::Index row = 0; row < integerRows.rows(); ++row)
    {
        long double meanScale = std::abs(static_cast<long double>(integers(row)));
        for (Eigen::Index column = 0; column < dimension; ++column)
            meanScale += std::abs(static_cast<long double>(integerRows(row,column))) *
                std::abs(static_cast<long double>(mean(column)));
        const long double tolerance = out.roundingFactor * std::max(1.0L, meanScale);
        if (tolerance > std::numeric_limits<double>::max())
        { out.reason = "OUTPUT_POSTERIOR_SCALE_OVERFLOW"; return out; }
        const double residual = std::abs(meanResidual(row));
        out.maximumMeanResidual = std::max(out.maximumMeanResidual, residual);
        out.maximumMeanTolerance = std::max(out.maximumMeanTolerance,
            static_cast<double>(tolerance));
        meanFailed |= residual > tolerance;
        for (Eigen::Index column = 0; column < dimension; ++column)
        {
            long double covarianceScale = 0;
            for (Eigen::Index inner = 0; inner < dimension; ++inner)
                covarianceScale += std::abs(static_cast<long double>(integerRows(row,inner))) *
                    std::abs(static_cast<long double>(covariance(inner,column)));
            const long double covarianceTolerance =
                out.roundingFactor * std::max(1.0L, covarianceScale);
            if (covarianceTolerance > std::numeric_limits<double>::max())
            { out.reason = "OUTPUT_POSTERIOR_SCALE_OVERFLOW"; return out; }
            const double covarianceError = std::abs(covarianceResidual(row,column));
            out.maximumCovarianceResidual = std::max(out.maximumCovarianceResidual, covarianceError);
            out.maximumCovarianceTolerance = std::max(out.maximumCovarianceTolerance,
                static_cast<double>(covarianceTolerance));
            covarianceFailed |= covarianceError > covarianceTolerance;
        }
    }
    if (meanFailed)
    { out.reason = "OUTPUT_POSTERIOR_MEAN_CLOSURE_FAILED"; return out; }
    if (covarianceFailed)
    { out.reason = "OUTPUT_POSTERIOR_INTEGER_NULLSPACE_FAILED"; return out; }
    out.valid = true;
    out.reason = "NUMERICAL_FIXED_OUTPUT_VALID";
    return out;
}
