#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <map>
#include <limits>
#include <string>
#include <vector>
#include "common/zhangCurrentEpochEvidence.hpp"

struct ZhangFinalProductPosteriorResult
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    int inputRows = 0;
    int independentRows = 0;
    double maximumConstraintResidual = 0;
    std::string validationMode = "FULL_SOURCE_AND_FINAL_PSD";
};

/** Condition one PRIVATE copy of the authoritative FLOAT posterior on all
 * current-epoch numerical constraints. Evidence authorization is separate.
 * No historical constraint, diagonal approximation, NIS, or explicit matrix
 * inverse is introduced here. Dependent deterministic rows must agree. */
inline ZhangFinalProductPosteriorResult zhangConditionFinalProductPosterior(
    const std::map<KFKey, int>& indices,
    const std::vector<ZhangCurrentEpochNumericConstraint>& constraints,
    Eigen::VectorXd& mean, Eigen::MatrixXd& covariance,
    bool sourcePosteriorValidated = false)
{
    ZhangFinalProductPosteriorResult out;
    out.inputRows = static_cast<int>(constraints.size());
    out.validationMode = sourcePosteriorValidated
        ? "UPSTREAM_SOURCE_PSD_VALIDATED_CONSTRAINT_AND_CLOSURE_CHECKS"
        : "FULL_SOURCE_AND_FINAL_PSD";
    const int n = mean.size();
    if (n <= 0 || covariance.rows() != n || covariance.cols() != n ||
        !mean.allFinite() || !covariance.allFinite())
    {
        out.reason = "INVALID_FLOAT_POSTERIOR";
        return out;
    }
    const double covarianceScale = std::max(1.0, covariance.cwiseAbs().maxCoeff());
    if ((covariance - covariance.transpose()).cwiseAbs().maxCoeff() >
        1e-10 * covarianceScale)
    { out.reason = "ASYMMETRIC_FLOAT_POSTERIOR"; return out; }
    auto positiveSemidefinite = [](const Eigen::MatrixXd& matrix)
    {
        const double scale = std::max(1.0, matrix.cwiseAbs().maxCoeff());
        // A covariance after exact conditioning is normally singular. LDLT
        // is not rank revealing: near-zero pivots may produce a negative D
        // despite a PSD matrix. Strict validation therefore tests its actual
        // symmetric spectrum; the certified production mode skips this full
        // decomposition, rather than substituting a diagonal-only PSD claim.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(matrix,
            Eigen::EigenvaluesOnly);
        return spectrum.info() == Eigen::Success && spectrum.eigenvalues().allFinite() &&
            spectrum.eigenvalues().minCoeff() >= -1e-11 * scale;
    };
    // The opt-in fast mode requires a previously certified PSD FLOAT root;
    // it does not turn a local diagonal check into a whole-matrix PSD proof.
    // In this mode PSD preservation follows from conditioning that trusted
    // root with the positive-semidefinite S solve below. Strict/default mode
    // performs whole-matrix checks for standalone callers and tests.
    if (!sourcePosteriorValidated && !positiveSemidefinite(covariance))
    { out.reason = "NON_PSD_FLOAT_POSTERIOR"; return out; }
    if (constraints.empty())
    {
        out.valid = true;
        out.reason = "NO_NEW_NUMERICAL_CONSTRAINTS";
        return out;
    }
    // Collapse exact duplicate KF rows before numerical decomposition. A
    // contradictory duplicate is an exact conflict, not a small innovation.
    std::map<ZhangCoordinateRow, ZhangExactInteger> unique;
    for (const auto& constraint : constraints)
    {
        ZhangCoordinateRow normalized;
        for (const auto& [key, coefficient] : constraint.rowCycles)
        {
            if (!std::isfinite(coefficient))
            { out.reason = "NONFINITE_NUMERICAL_CONSTRAINT"; return out; }
            if (coefficient != 0) normalized.emplace(key, coefficient);
        }
        if (normalized.empty())
        {
            out.reason = "EMPTY_NUMERICAL_CONSTRAINT";
            return out;
        }
        auto [it, inserted] = unique.emplace(std::move(normalized), constraint.value);
        if (!inserted && it->second != constraint.value)
        {
            out.reason = "EXACT_NUMERICAL_CONSTRAINT_CONFLICT";
            return out;
        }
    }
    const int m = unique.size();
    // These are sparse KF functionals. Do not inflate them into m*n doubles
    // or evaluate the final C*P closure as a dense m*n*n multiplication.
    Eigen::SparseMatrix<double, Eigen::RowMajor> C(m, n);
    std::vector<Eigen::Triplet<double>> coefficients;
    std::size_t coefficientCount = 0;
    for (const auto& [terms, integer] : unique) coefficientCount += terms.size();
    coefficients.reserve(coefficientCount);
    Eigen::VectorXd z(m);
    int row = 0;
    for (const auto& [terms, integer] : unique)
    {
        try { z(row) = integer.convert_to<double>(); }
        catch (...) { out.reason = "INTEGER_RHS_OUT_OF_RANGE"; return out; }
        if (!std::isfinite(z(row)))
        { out.reason = "INTEGER_RHS_NONFINITE"; return out; }
        if (ZhangExactInteger(z(row)) != integer)
        { out.reason = "INTEGER_RHS_NOT_EXACTLY_REPRESENTABLE"; return out; }
        for (const auto& [key, coefficient] : terms)
        {
            const auto index = indices.find(key);
            if (index == indices.end() || index->second < 0 || index->second >= n ||
                !std::isfinite(coefficient))
            { out.reason = "CONSTRAINT_KF_COORDINATE_UNAVAILABLE"; return out; }
            coefficients.emplace_back(row, index->second, coefficient);
        }
        row++;
    }
    C.setFromTriplets(coefficients.begin(), coefficients.end());
    C.prune(0.0);
    for (int row = 0; row < m; ++row)
        if (!Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(C, row))
        { out.reason = "ZERO_NUMERICAL_CONSTRAINT"; return out; }
    const Eigen::MatrixXd cross = covariance * C.transpose();
    Eigen::MatrixXd S = C * cross;
    S = (0.5 * (S + S.transpose())).eval();
    const Eigen::VectorXd innovation = z - C * mean;
    const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(z.cwiseAbs().maxCoeff(),
            (C.cwiseAbs() * mean.cwiseAbs()).maxCoeff()));
    if (!std::isfinite(roundoff) || roundoff > 1e-5)
    { out.reason = "NUMERICAL_INTEGER_CLOSURE_PRECISION_INSUFFICIENT"; return out; }
    const double closureTolerance = std::max(1e-8, roundoff);
    const double scale = std::max(1.0, S.cwiseAbs().maxCoeff());
    const double tolerance = 1e-11 * scale;
    Eigen::LDLT<Eigen::MatrixXd> factor(S);
    Eigen::MatrixXd weights;
    Eigen::VectorXd delta;
    if (factor.info() == Eigen::Success && factor.vectorD().minCoeff() > tolerance)
    {
        out.independentRows = m;
        delta = factor.solve(innovation);
        weights = factor.solve(cross.transpose());
    }
    else
    {
        // Exact dependent rows and already-conditioned directions are legal.
        // Keep their consistency checks; never replace a singular S by jitter.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectral(S);
        if (spectral.info() != Eigen::Success ||
            spectral.eigenvalues().minCoeff() < -tolerance)
        { out.reason = "INVALID_CONSTRAINT_COVARIANCE"; return out; }
        Eigen::VectorXd reciprocal = Eigen::VectorXd::Zero(m);
        const Eigen::VectorXd coordinates = spectral.eigenvectors().transpose() * innovation;
        for (int i = 0; i < m; i++)
        {
            if (spectral.eigenvalues()(i) > tolerance)
            {
                reciprocal(i) = 1.0 / spectral.eigenvalues()(i);
                out.independentRows++;
            }
            else if (std::abs(coordinates(i)) > closureTolerance * std::sqrt(double(m)))
            { out.reason = "DETERMINISTIC_CONSTRAINT_CONFLICT"; return out; }
        }
        delta = spectral.eigenvectors() * reciprocal.asDiagonal() * coordinates;
        weights = spectral.eigenvectors() * reciprocal.asDiagonal() *
            (spectral.eigenvectors().transpose() * cross.transpose());
    }
    if (!delta.allFinite() || !weights.allFinite())
    { out.reason = "CONDITION_SOLVE_NONFINITE"; return out; }
    const Eigen::VectorXd candidateMean = mean + cross * delta;
    out.maximumConstraintResidual = (C * candidateMean - z).cwiseAbs().maxCoeff();
    if (!candidateMean.allFinite() ||
        out.maximumConstraintResidual > closureTolerance * std::sqrt(double(m)))
    { out.reason = "FINAL_CONSTRAINT_MEAN_CLOSURE_FAILED"; return out; }
    // Mutate only after validation. Preserve all clock/phase/other cross terms.
    Eigen::MatrixXd candidateCovariance = covariance;
    candidateCovariance.noalias() -= cross * weights;
    // Only one additional full P is needed to preserve atomic failure. Blend
    // each symmetric pair in place, not through another n*n transpose buffer.
    for (int column = 0; column < n; ++column)
        for (int row = column + 1; row < n; ++row)
        {
            const double symmetric = 0.5 * candidateCovariance(row, column) +
                0.5 * candidateCovariance(column, row);
            candidateCovariance(row, column) = candidateCovariance(column, row) = symmetric;
        }
    if (!candidateCovariance.allFinite() || candidateCovariance.diagonal().minCoeff() <
        -1e-8 * std::max(1.0, candidateCovariance.diagonal().cwiseAbs().maxCoeff()))
    { out.reason = "FINAL_POSTERIOR_NUMERIC_INVALID"; return out; }
    if (!sourcePosteriorValidated && !positiveSemidefinite(candidateCovariance))
    { out.reason = "NON_PSD_FINAL_POSTERIOR"; return out; }
    double constraintCovarianceResidual = 0;
    for (int column = 0; column < n; column += 64)
    {
        const Eigen::MatrixXd closure = C * candidateCovariance.middleCols(
            column, std::min(64, n - column));
        if (!closure.allFinite())
        { out.reason = "FINAL_CONSTRAINT_COVARIANCE_CLOSURE_FAILED"; return out; }
        constraintCovarianceResidual = std::max(constraintCovarianceResidual,
            closure.cwiseAbs().maxCoeff());
    }
    if (!std::isfinite(constraintCovarianceResidual) ||
        constraintCovarianceResidual > 1e-8 * std::max(1.0, cross.cwiseAbs().maxCoeff()))
    { out.reason = "FINAL_CONSTRAINT_COVARIANCE_CLOSURE_FAILED"; return out; }
    mean = candidateMean;
    covariance = std::move(candidateCovariance);
    out.valid = true;
    out.reason = "CURRENT_CONSTRAINTS_FULL_JOINT_POSTERIOR";
    return out;
}
