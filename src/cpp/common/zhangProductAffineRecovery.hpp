#pragma once

#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangProductRelationSolver.hpp"

/** Exact recovery of a target on an integer-feasible history domain.
 * The certificate is
 *   target = posteriorRow * currentRows +
 *            (historyNumerators / historyDenominator) * historyRows,
 *   offset = (historyNumerators / historyDenominator) * historyValues.
 * posteriorRow is integral; rational history coefficients are allowed only
 * after integer feasibility has been established. No missing physical arc is
 * identified with zero by this operation. Before publication, the caller must
 * bind every nonzero historyNumerator to its versioned physical row and
 * accepted decision proof; this algebraic certificate alone is not a product
 * authorization. */
struct ZhangProductAffineRecovery
{
    bool valid = false;
    ZhangExactVector posteriorRow;
    ZhangExactInteger offset = 0;
    ZhangExactVector historyNumerators;
    ZhangExactInteger historyDenominator = 1;
    std::string reason = "NOT_EVALUATED";
};

inline bool zhangVerifyProductAffineRecovery(
    const ZhangExactMatrix& currentRows,
    const ZhangExactMatrix& historyRows,
    const ZhangExactVector& historyValues,
    const ZhangExactVector& target,
    const ZhangProductAffineRecovery& proof)
{
    const std::size_t dimension = target.size();
    if (!proof.valid || proof.historyDenominator <= 0 ||
        currentRows.size() != proof.posteriorRow.size() ||
        historyRows.size() != historyValues.size() ||
        historyRows.size() != proof.historyNumerators.size()) return false;
    ZhangExactVector residual = target;
    for (std::size_t row = 0; row < currentRows.size(); ++row)
    {
        if (currentRows[row].size() != dimension) return false;
        for (std::size_t column = 0; column < dimension; ++column)
            residual[column] -= proof.posteriorRow[row] * currentRows[row][column];
    }
    ZhangExactInteger rhs = 0;
    for (std::size_t row = 0; row < historyRows.size(); ++row)
    {
        if (historyRows[row].size() != dimension) return false;
        rhs += proof.historyNumerators[row] * historyValues[row];
    }
    if (historyRows.empty())
        for (auto value : residual) if (value != 0) return false;
    else
        for (std::size_t column = 0; column < dimension; ++column)
        {
            ZhangExactInteger represented = 0;
            for (std::size_t row = 0; row < historyRows.size(); ++row)
                represented += proof.historyNumerators[row] * historyRows[row][column];
            if (proof.historyDenominator * residual[column] != represented)
                return false;
        }
    return rhs == proof.historyDenominator * proof.offset;
}

inline ZhangProductAffineRecovery zhangRecoverProductAffineTarget(
    const ZhangExactMatrix& currentRows,
    const ZhangExactMatrix& historyRows,
    const ZhangExactVector& historyValues,
    const ZhangExactVector& target)
{
    ZhangProductAffineRecovery out;
    const int dimension = static_cast<int>(target.size());
    if (dimension <= 0 || historyRows.size() != historyValues.size())
    {
        out.reason = "AFFINE_RECOVERY_DIMENSION_MISMATCH";
        return out;
    }
    for (const auto& row : currentRows)
        if (row.size() != target.size())
        {
            out.reason = "AFFINE_RECOVERY_CURRENT_ROW_DIMENSION";
            return out;
        }
    for (const auto& row : historyRows)
        if (row.size() != target.size())
        {
            out.reason = "AFFINE_RECOVERY_HISTORY_ROW_DIMENSION";
            return out;
        }

    const auto domain = zhangExactAffineIntegerQuotient(
        historyRows, historyValues, dimension,
        ZhangExactQuotientWork::PARTICULAR_ONLY);
    if (!domain.valid)
    {
        out.reason = domain.failureReason;
        return out;
    }
    ZhangExactMatrix projectedTargets = currentRows;
    projectedTargets.push_back(target);
    const auto projected = zhangExactProjectedKernelColumns(
        historyRows, projectedTargets, dimension);
    if (projected.size() != projectedTargets.size())
    {
        out.reason = "AFFINE_RECOVERY_PROJECTED_KERNEL_INVALID";
        return out;
    }
    ZhangExactMatrix currentImage(projected.begin(), projected.end() - 1);
    if (projected.back().empty())
    {
        // A full-rank history domain has no free integer direction.  The
        // zero-dimensional membership API carries no combination vector.
        out.posteriorRow = ZhangExactVector(currentRows.size());
    }
    else
    {
        const auto membership = zhangIntegerRowLatticeContainsBatch(
            currentImage, {projected.back()});
        if (membership.size() != 1 || !membership.front().contained ||
            membership.front().combination.size() != currentRows.size())
        {
            out.reason = "TARGET_HAS_UNRESOLVED_INTEGER_DIRECTION";
            return out;
        }
        out.posteriorRow = membership.front().combination;
    }
    ZhangExactVector residual = target;
    for (std::size_t row = 0; row < currentRows.size(); ++row)
        for (int column = 0; column < dimension; ++column)
            residual[column] -= out.posteriorRow[row] * currentRows[row][column];
    for (int column = 0; column < dimension; ++column)
        out.offset += residual[column] * domain.particularSolution[column];

    const auto hnf = zhangExactRowHermiteNormalForm(
        historyRows, historyValues, true);
    if (!hnf.consistent)
    {
        out.reason = "AFFINE_RECOVERY_HISTORY_INCONSISTENT";
        return out;
    }
    out.historyNumerators.resize(historyRows.size());
    ZhangExactVector remaining = residual;
    for (std::size_t row = 0; row < hnf.basis.size(); ++row)
    {
        const auto& basis = hnf.basis[row];
        std::size_t pivot = 0;
        while (pivot < basis.size() && basis[pivot] == 0) ++pivot;
        if (pivot == basis.size() || remaining[pivot] == 0) continue;
        const auto oldPivot = remaining[pivot];
        const auto multiplier = basis[pivot] < 0 ? -oldPivot : oldPivot;
        const auto denominatorFactor = zhangExactAbs(basis[pivot]);
        out.historyDenominator *= denominatorFactor;
        for (auto& value : remaining) value *= denominatorFactor;
        for (auto& value : out.historyNumerators) value *= denominatorFactor;
        for (std::size_t column = 0; column < dimension; ++column)
            remaining[column] -= multiplier * basis[column];
        for (std::size_t source = 0; source < historyRows.size(); ++source)
            out.historyNumerators[source] +=
                multiplier * hnf.rowTransform[row][source];
    }
    for (const auto& value : remaining)
        if (value != 0)
        {
            out.reason = "AFFINE_RECOVERY_RATIONAL_WITNESS_FAILED";
            return out;
        }
    out.valid = true;
    out.reason = "EXACT_AFFINE_RECOVERY";
    if (!zhangVerifyProductAffineRecovery(
            currentRows, historyRows, historyValues, target, out))
    {
        out.valid = false;
        out.reason = "AFFINE_RECOVERY_WITNESS_MISMATCH";
    }
    return out;
}
