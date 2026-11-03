#pragma once

#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangProductRelationSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

// This workspace is call-scoped. The FLOAT posterior is never changed here;
// held integer equations condition a private copy, and only the exact free
// integer coordinates are offered to the next ILS/PAR search. An uncertified
// covariance null direction is an error, not an invitation to add jitter.
struct ZhangUserResidualQuotient
{
    bool valid = false;
    ZhangExactAffineIntegerQuotient exact;
    VectorXd conditionalMean;
    MatrixXd conditionalCovariance;
    VectorXd quotientMean;
    MatrixXd quotientCovariance;
    double minimumQuotientEigenvalue = std::numeric_limits<double>::quiet_NaN();
    double positiveVarianceTolerance = std::numeric_limits<double>::quiet_NaN();
    std::string reason = "NOT_EVALUATED";
};

struct ZhangUserIntegerUnion
{
    bool valid = false;
    int ambientDimension = 0;
    int rank = 0;
    int quotientDimension = 0;
    ZhangExactMatrix rows;
    ZhangExactVector values;
    std::string reason = "NOT_EVALUATED";
};

struct ZhangUserQuotientDecision
{
    bool valid = false;
    ZhangExactMatrix ambientRows;
    ZhangExactVector ambientValues;
    ZhangUserIntegerUnion united;
    int independentAddedRank = 0;
    std::string reason = "NOT_EVALUATED";
};

struct ZhangUserCandidateIntegerImage
{
    bool valid = false;
    int candidateDimension = 0;
    int fixedRank = 0;
    int freeImageRank = 0;
    // Complete candidate integer set: particularTarget + imageGenerators*Z^r.
    // The generators are NOT saturated: offset+2Z must remain offset+2Z.
    ZhangExactVector particularTarget;
    ZhangExactMatrix imageGenerators;
    // C*t = values are all integer equality consequences on candidate targets.
    // These are logical consequences of a feasible integer affine domain, NOT
    // necessarily integer-row-combination witnesses for a product certificate.
    ZhangExactMatrix fixedRows;
    ZhangExactVector fixedValues;
    ZhangExactMatrix ambientFixedRows;
    std::string reason = "NOT_EVALUATED";
};

inline bool zhangUserExactMatrixDimensions(
    const ZhangExactMatrix& rows, std::size_t columns)
{
    for (const auto& row : rows) if (row.size() != columns) return false;
    return true;
}

// Candidate masking is a TARGET selection, not permission to erase masked
// physical variables from H. Keep the full affine integer domain, including
// below-mask latent variables, and project T*(N0+K^T*q) exactly. In particular,
// 2*N_masked+N_live=1 leaves N_live odd even though no live equality is fixed.
inline ZhangUserCandidateIntegerImage zhangProjectUserCandidateIntegerImage(
    const ZhangExactAffineIntegerQuotient& domain,
    const ZhangExactMatrix& candidateSelectors)
{
    ZhangUserCandidateIntegerImage out;
    const int ambient = domain.ambientDimension;
    const int free = domain.quotientRank;
    out.candidateDimension = candidateSelectors.size();
    if (!domain.valid || ambient < 0 || free < 0 || domain.deterministicRank < 0 ||
        domain.deterministicRank + free != ambient ||
        domain.particularSolution.size() != static_cast<std::size_t>(ambient) ||
        domain.kernelBasis.size() != static_cast<std::size_t>(free) ||
        !zhangUserExactMatrixDimensions(domain.kernelBasis, ambient) ||
        !zhangUserExactMatrixDimensions(candidateSelectors, ambient))
    { out.reason = "CANDIDATE_INTEGER_IMAGE_INPUT_INVALID"; return out; }
    if (candidateSelectors.empty())
    {
        out.valid = true;
        out.reason = "EMPTY_CANDIDATE_INTEGER_IMAGE";
        return out;
    }
    const auto independent = zhangExactRowHermiteNormalForm(candidateSelectors);
    if (!independent.consistent || independent.basis.size() != candidateSelectors.size())
    { out.reason = "CANDIDATE_TARGET_ROWS_NOT_INDEPENDENT"; return out; }
    const int candidates = candidateSelectors.size();
    out.particularTarget.assign(candidates, 0);
    ZhangExactMatrix targetKernel(candidates, ZhangExactVector(free));
    for (int target = 0; target < candidates; ++target)
    for (int column = 0; column < ambient; ++column)
    {
        out.particularTarget[target] += candidateSelectors[target][column] *
            domain.particularSolution[column];
        for (int direction = 0; direction < free; ++direction)
            targetKernel[target][direction] += candidateSelectors[target][column] *
                domain.kernelBasis[direction][column];
    }
    if (free)
    {
        const auto image = zhangPrimitiveImageCoordinates(targetKernel);
        if (!image.valid)
        { out.reason = "CANDIDATE_INTEGER_IMAGE_GENERATORS_FAILED"; return out; }
        out.imageGenerators = image.imageGenerators;
        out.freeImageRank = image.primitiveRows.size();
    }
    else out.imageGenerators.assign(candidates, ZhangExactVector{});
    out.fixedRank = candidates - out.freeImageRank;
    if (out.fixedRank < 0)
    { out.reason = "CANDIDATE_INTEGER_IMAGE_OVER_RANK"; return out; }
    ZhangExactMatrix transpose(free, ZhangExactVector(candidates));
    for (int direction = 0; direction < free; ++direction)
    for (int target = 0; target < candidates; ++target)
        transpose[direction][target] = targetKernel[target][direction];
    const auto equalityKernel = zhangExactIntegerKernel(transpose, candidates);
    ZhangExactVector equalityValues(equalityKernel.size());
    for (std::size_t row = 0; row < equalityKernel.size(); ++row)
    for (int target = 0; target < candidates; ++target)
        equalityValues[row] += equalityKernel[row][target] * out.particularTarget[target];
    const auto equalities = zhangExactRowHermiteNormalForm(equalityKernel, equalityValues);
    if (!equalities.consistent || equalities.basis.size() != static_cast<std::size_t>(out.fixedRank))
    { out.reason = "CANDIDATE_FIXED_EQUALITY_RANK_MISMATCH"; return out; }
    out.fixedRows = equalities.basis;
    out.fixedValues = equalities.values;
    out.ambientFixedRows = zhangExactMultiply(out.fixedRows, candidateSelectors);
    out.valid = true;
    out.reason = "FULL_AFFINE_DOMAIN_CANDIDATE_INTEGER_IMAGE";
    return out;
}

// Integer coefficients/origins must survive the subsequent double operation
// exactly. Silently converting a huge exact witness into a rounded coefficient
// would invalidate the claimed correspondence between the two search charts.
inline bool zhangUserExactToDouble(
    const ZhangExactMatrix& source, int columns, MatrixXd& destination)
{
    if (columns < 0 || !zhangUserExactMatrixDimensions(source, columns)) return false;
    const ZhangExactInteger limit = ZhangExactInteger(1) << 53;
    destination.resize(source.size(), columns);
    for (int row = 0; row < destination.rows(); ++row)
    for (int column = 0; column < columns; ++column)
    {
        const auto& value = source[row][column];
        if (value < -limit || value > limit) return false;
        destination(row, column) = value.convert_to<double>();
    }
    return destination.allFinite();
}

inline bool zhangUserExactToDouble(
    const ZhangExactVector& source, VectorXd& destination)
{
    const ZhangExactInteger limit = ZhangExactInteger(1) << 53;
    destination.resize(source.size());
    for (int index = 0; index < destination.size(); ++index)
    {
        const auto& value = source[index];
        if (value < -limit || value > limit) return false;
        destination(index) = value.convert_to<double>();
    }
    return destination.allFinite();
}

inline ZhangUserIntegerUnion zhangAuditUserIntegerUnion(
    const ZhangExactMatrix& rows, const ZhangExactVector& values,
    int ambientDimension)
{
    ZhangUserIntegerUnion out;
    out.ambientDimension = ambientDimension;
    if (ambientDimension < 0 || rows.size() != values.size() ||
        !zhangUserExactMatrixDimensions(rows, ambientDimension))
    { out.reason = "INTEGER_UNION_DIMENSION_MISMATCH"; return out; }
    if (ambientDimension == 0)
    {
        for (const auto& value : values)
            if (value != 0)
            { out.reason = "INTEGER_UNION_ZERO_DIMENSION_CONFLICT"; return out; }
        out.valid = true;
        out.reason = "EMPTY_INTEGER_UNION";
        return out;
    }
    const auto exact = zhangExactAffineIntegerQuotient(rows, values,
        ambientDimension, ZhangExactQuotientWork::FEASIBILITY_ONLY);
    if (!exact.valid) { out.reason = exact.failureReason; return out; }
    out.rows = exact.deterministicBasis;
    out.values = exact.deterministicValues;
    out.rank = exact.deterministicRank;
    out.quotientDimension = ambientDimension - out.rank;
    out.valid = true;
    out.reason = "EXACT_INTEGER_UNION_FEASIBLE";
    return out;
}

inline ZhangUserResidualQuotient zhangBuildUserResidualQuotient(
    const VectorXd& mean, const MatrixXd& covariance,
    const ZhangExactMatrix& heldRows, const ZhangExactVector& heldValues)
{
    ZhangUserResidualQuotient out;
    const int ambient = mean.size();
    if (covariance.rows() != ambient || covariance.cols() != ambient ||
        !mean.allFinite() || !covariance.allFinite() ||
        heldRows.size() != heldValues.size() ||
        !zhangUserExactMatrixDimensions(heldRows, ambient))
    { out.reason = "RESIDUAL_QUOTIENT_INPUT_INVALID"; return out; }
    if (ambient == 0)
    {
        const auto empty = zhangAuditUserIntegerUnion(heldRows, heldValues, 0);
        if (!empty.valid) { out.reason = empty.reason; return out; }
        out.exact.valid = true;
        out.exact.failureReason = "NONE";
        out.conditionalMean = mean;
        out.conditionalCovariance = covariance;
        out.quotientMean.resize(0);
        out.quotientCovariance.resize(0, 0);
        out.valid = true;
        out.reason = "NO_RESIDUAL_INTEGER_DIRECTIONS";
        return out;
    }
    out.exact = zhangExactAffineIntegerQuotient(heldRows, heldValues, ambient);
    if (!out.exact.valid) { out.reason = out.exact.failureReason; return out; }
    MatrixXd constraints, projector;
    VectorXd integers, origin;
    if (!zhangUserExactToDouble(out.exact.deterministicBasis, ambient, constraints) ||
        !zhangUserExactToDouble(out.exact.deterministicValues, integers) ||
        !zhangUserExactToDouble(out.exact.quotientProjector, ambient, projector) ||
        !zhangUserExactToDouble(out.exact.particularSolution, origin))
    { out.reason = "RESIDUAL_EXACT_COEFFICIENT_NOT_DOUBLE_REPRESENTABLE"; return out; }
    // Validate PSD even when there are no held rows. A private conditioning
    // result cannot legalize an already invalid FLOAT covariance.
    const auto workspace = zhangBuildPosteriorEffectiveWorkspace(covariance);
    if (!workspace.valid) { out.reason = workspace.failureReason; return out; }
    out.conditionalMean = mean;
    out.conditionalCovariance = 0.5 * (covariance + covariance.transpose());
    if (out.exact.deterministicRank)
    {
        const auto conditioned = zhangConditionPosteriorEffectiveIntegers(
            mean, out.conditionalCovariance, constraints, integers);
        if (!conditioned.valid) { out.reason = conditioned.failureReason; return out; }
        out.conditionalMean = conditioned.mean;
        out.conditionalCovariance = conditioned.covariance;
    }
    if (!out.conditionalMean.allFinite() || !out.conditionalCovariance.allFinite())
    { out.reason = "RESIDUAL_CONDITIONED_MOMENTS_NONFINITE"; return out; }
    if (constraints.rows())
    {
        const VectorXd closure = constraints * out.conditionalMean - integers;
        const double valueScale = std::max(1.0, integers.cwiseAbs().maxCoeff());
        if (closure.cwiseAbs().maxCoeff() > 1e-7 * valueScale)
        { out.reason = "RESIDUAL_CONDITIONED_AFFINE_CLOSURE_FAILED"; return out; }
    }
    if (out.exact.quotientRank == 0)
    {
        // Production Eigen/BLAS must not receive an M=0 GEMV with LDA=0.
        // A completely held integer space has no searchable moments.
        out.quotientMean.resize(0);
        out.quotientCovariance.resize(0,0);
        out.valid = true;
        out.reason = "NO_RESIDUAL_INTEGER_DIRECTIONS";
        return out;
    }
    out.quotientMean = projector * (out.conditionalMean - origin);
    out.quotientCovariance = projector * out.conditionalCovariance * projector.transpose();
    out.quotientCovariance = (0.5 *
        (out.quotientCovariance + out.quotientCovariance.transpose())).eval();
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(out.quotientCovariance);
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite() ||
        !out.quotientMean.allFinite())
    { out.reason = "RESIDUAL_QUOTIENT_EIGEN_FAILED"; return out; }
    const double scale = std::max(1.0, eigen.eigenvalues().cwiseAbs().maxCoeff());
    out.positiveVarianceTolerance = std::max(1e-14,
        scale * std::numeric_limits<double>::epsilon() * out.exact.quotientRank * 64);
    out.minimumQuotientEigenvalue = eigen.eigenvalues().minCoeff();
    if (out.minimumQuotientEigenvalue <= out.positiveVarianceTolerance)
    { out.reason = "RESIDUAL_QUOTIENT_NOT_POSITIVE_DEFINITE"; return out; }
    out.valid = true;
    out.reason = "SEARCHABLE_EXACT_RESIDUAL_INTEGER_QUOTIENT";
    return out;
}

inline ZhangUserQuotientDecision zhangMapUserQuotientDecision(
    const ZhangUserResidualQuotient& quotient,
    const ZhangExactMatrix& rows, const ZhangExactVector& values)
{
    ZhangUserQuotientDecision out;
    if (!quotient.valid || !quotient.exact.valid || rows.size() != values.size() ||
        !zhangUserExactMatrixDimensions(rows, quotient.exact.quotientRank))
    { out.reason = "QUOTIENT_DECISION_INPUT_INVALID"; return out; }
    const int ambient = quotient.exact.ambientDimension;
    out.ambientRows.reserve(rows.size());
    out.ambientValues = values;
    // Explicit loops keep the zero-free-dimension case defined: an accepted
    // zero row with nonzero RHS must reach the exact union conflict check.
    for (std::size_t row = 0; row < rows.size(); ++row)
    {
        ZhangExactVector ambientRow(ambient);
        for (int free = 0; free < quotient.exact.quotientRank; ++free)
        for (int column = 0; column < ambient; ++column)
            ambientRow[column] += rows[row][free] *
                quotient.exact.quotientProjector[free][column];
        for (int column = 0; column < ambient; ++column)
            out.ambientValues[row] += ambientRow[column] *
                quotient.exact.particularSolution[column];
        out.ambientRows.push_back(std::move(ambientRow));
    }
    auto unitedRows = quotient.exact.deterministicBasis;
    auto unitedValues = quotient.exact.deterministicValues;
    unitedRows.insert(unitedRows.end(), out.ambientRows.begin(), out.ambientRows.end());
    unitedValues.insert(unitedValues.end(), out.ambientValues.begin(), out.ambientValues.end());
    out.united = zhangAuditUserIntegerUnion(unitedRows, unitedValues, ambient);
    if (!out.united.valid) { out.reason = out.united.reason; return out; }
    out.independentAddedRank = out.united.rank - quotient.exact.deterministicRank;
    out.valid = true;
    out.reason = "EXACT_QUOTIENT_DECISION_MAPPED_AND_UNITED";
    return out;
}

struct ZhangUserRestrictedIntegerLedger
{
    bool valid = false;
    // This function preserves all equality consequences obtainable by exact
    // integer combinations. It is not a full existential integer projection:
    // discarded latent variables can additionally imply modular constraints.
    bool equalityConsequencesOnly = true;
    // True only when full integer projection equivalence has been established
    // (all columns retained, or the lost coefficient image is primitive).
    bool knownFullProjection = false;
    int originalDimension = 0;
    int retainedDimension = 0;
    int lostColumns = 0;
    int discardedRank = 0;
    ZhangExactMatrix rows;
    ZhangExactVector values;
    ZhangExactMatrix sourceCombinations;
    std::string reason = "NOT_EVALUATED";
};

inline ZhangUserRestrictedIntegerLedger zhangRestrictUserIntegerLedger(
    const ZhangExactMatrix& rows, const ZhangExactVector& values,
    const std::vector<int>& retainedColumns, int emptyAmbientDimension = 0)
{
    ZhangUserRestrictedIntegerLedger out;
    const int ambient = rows.empty() ? emptyAmbientDimension : rows.front().size();
    out.originalDimension = ambient;
    out.retainedDimension = retainedColumns.size();
    if (ambient < 0 || rows.size() != values.size() ||
        !zhangUserExactMatrixDimensions(rows, ambient))
    { out.reason = "LEDGER_RESTRICTION_INPUT_INVALID"; return out; }
    std::set<int> selected;
    for (int column : retainedColumns)
    {
        if (column < 0 || column >= ambient || !selected.insert(column).second)
        { out.reason = "LEDGER_RESTRICTION_COLUMN_INVALID"; return out; }
    }
    const auto original = zhangAuditUserIntegerUnion(rows, values, ambient);
    if (!original.valid) { out.reason = original.reason; return out; }
    std::vector<int> lost;
    for (int column = 0; column < ambient; ++column)
        if (!selected.count(column)) lost.push_back(column);
    out.lostColumns = lost.size();
    // For combination c, every discarded coefficient must vanish:
    //     (H_lost)^T c = 0.
    // The saturated integer kernel enumerates ALL such integer combinations.
    ZhangExactMatrix lostTranspose(lost.size(), ZhangExactVector(rows.size()));
    for (std::size_t column = 0; column < lost.size(); ++column)
    for (std::size_t row = 0; row < rows.size(); ++row)
        lostTranspose[column][row] = rows[row][lost[column]];
    const auto combinations = zhangExactIntegerKernel(lostTranspose, rows.size());
    ZhangExactMatrix projected;
    ZhangExactVector projectedValues;
    for (const auto& combination : combinations)
    {
        ZhangExactVector projectedRow(retainedColumns.size());
        ZhangExactInteger rhs = 0;
        for (std::size_t row = 0; row < rows.size(); ++row)
        {
            rhs += combination[row] * values[row];
            for (std::size_t column = 0; column < retainedColumns.size(); ++column)
                projectedRow[column] += combination[row] * rows[row][retainedColumns[column]];
        }
        projected.push_back(std::move(projectedRow));
        projectedValues.push_back(std::move(rhs));
    }
    const auto hnf = zhangExactRowHermiteNormalForm(projected, projectedValues, true);
    if (!hnf.consistent)
    { out.reason = "LEDGER_RESTRICTION_AFFINE_CONFLICT"; return out; }
    const auto retained = zhangAuditUserIntegerUnion(hnf.basis, hnf.values,
        retainedColumns.size());
    if (!retained.valid) { out.reason = retained.reason; return out; }
    out.rows = retained.rows;
    out.values = retained.values;
    out.discardedRank = original.rank - retained.rank;
    if (lost.empty() || rows.empty()) out.knownFullProjection = true;
    else
    {
        // Primitive lost-column image has no hidden divisibility conditions:
        // integer membership is then exactly rational-image membership.
        const auto lostImage = zhangIntegerRowLatticeContains(
            lostTranspose, ZhangExactVector(rows.size()));
        out.knownFullProjection = lostImage.contained;
        for (const auto& invariant : lostImage.smithInvariants)
            if (invariant != 1 && invariant != -1) out.knownFullProjection = false;
    }
    // Store exact witnesses without another membership factorization.
    out.sourceCombinations = zhangExactMultiply(hnf.rowTransform, combinations);
    out.valid = true;
    out.reason = "ALL_RETAINED_INTEGER_EQUALITY_CONSEQUENCES";
    return out;
}
