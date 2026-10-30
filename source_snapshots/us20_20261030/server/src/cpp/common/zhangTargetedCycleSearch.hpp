#pragma once

#include "common/zhangR48Bridge.hpp"

#include <functional>
#include <limits>
#include <memory>
#include <vector>

// This result intentionally retains the existing bridge receipt layout. rank
// is the rank of the requested INTEGER image; wlRank/l1Rank have no meaning
// for a general cycle block and are not populated by this helper.
using ZhangTargetedCycleSearchResult = ZhangR48BridgeResult;

/** Search a small current-network integer image, not necessarily a WL/L1 pair.
 *
 * All target and conditioner rows are in the SAME immutable root coordinates.
 * A borrowed workspace must belong to exactly the supplied root mean/full Q;
 * it is only a numerical cache, never a source of integer authorization. A
 * cached domain is reused only after its COMPLETE affine H/b is checked below.
 * Targets outside that domain, missing integer kernels, or stale domains are
 * compiled from the complete supplied H/b instead.
 *
 * The callback receives the full conditional image covariance, never its
 * diagonal. Its final bool is true, as in zhangR48SearchBridge: this is one
 * joint block request. allowPartial controls how many independent returned
 * rows may be admitted; the caller also configures the callback's PAR policy.
 * Ratio/ILS acceptance belongs to that callback. No row is authorized merely
 * because it is useful to the eventual product target.
 *
 * Exactly one callback can be made here. spent is charged immediately before
 * it, and is not refunded for an invalid proposal, infeasible integer union,
 * or failed numerical gate. The owner must also charge its shared attempt
 * counter before entering the actual ILS callback. This private helper does
 * not mutate a route, FLOAT state, ledger, registry, proof, or product graph.
 */
inline ZhangTargetedCycleSearchResult zhangSearchTargetedCycleBlock(
    const Eigen::VectorXd& rootMean, const Eigen::MatrixXd& rootCovariance,
    const ZhangExactMatrix& targets, const ZhangExactMatrix& held,
    const ZhangExactVector& heldValues, double allocation, double alpha,
    const std::function<ZhangSequentialShadowProposal(
        const Eigen::VectorXd&, const Eigen::MatrixXd&, double, bool)>& search,
    ZhangR48MarginalWorkspace* sharedWorkspace = nullptr,
    const ZhangR47ProductSearchFrame* domainFrame = nullptr,
    bool allowPartial = false)
{
    ZhangTargetedCycleSearchResult out;
    const int dimension = rootMean.size();
    if (zhangRatioOnly())
        out.schurIdentityError = std::numeric_limits<double>::quiet_NaN();
    if (!(allocation > 0) || !std::isfinite(allocation))
    {
        out.status = "NOT_EVALUATED_NO_SEARCH_BUDGET";
        return out;
    }
    if (dimension <= 0 || rootCovariance.rows() != dimension ||
        rootCovariance.cols() != dimension || !rootMean.allFinite() ||
        !rootCovariance.allFinite() || targets.empty() || targets.size() > 8 ||
        held.size() != heldValues.size() ||
        !zhangExactRectangularMatrix(targets, dimension) ||
        !zhangExactRectangularMatrix(held, dimension) ||
        !(alpha > 0 && alpha < 1) || !search)
    {
        out.status = "TARGETED_CYCLE_INVALID_INPUT";
        return out;
    }

    // HNF rank is exact integer algebra. In particular, a real row-space or
    // a matching number of rows cannot establish equality of affine lattices.
    const auto base = zhangExactRowHermiteNormalForm(held, heldValues);
    if (!base.consistent)
    {
        out.status = "AFFINE_INFEASIBLE";
        return out;
    }
    const auto cacheMatchesCompleteDomain = [&]()
    {
        if (!domainFrame || !domainFrame->valid || !domainFrame->affine ||
            !domainFrame->affine->valid)
            return false;
        const auto& affine = *domainFrame->affine;
        const auto width = domainFrame->columns.size();
        if (affine.ambientDimension != static_cast<int>(width) ||
            affine.deterministicBasis.size() != affine.deterministicValues.size() ||
            !zhangExactRectangularMatrix(affine.deterministicBasis, width) ||
            affine.quotientRank < 0 ||
            affine.kernelBasis.size() != static_cast<std::size_t>(affine.quotientRank) ||
            !zhangExactRectangularMatrix(affine.kernelBasis, width))
            return false;
        std::vector<bool> represented(dimension, false);
        for (int column : domainFrame->columns)
        {
            if (column < 0 || column >= dimension || represented[column])
                return false;
            represented[column] = true;
        }
        ZhangExactMatrix expanded(affine.deterministicBasis.size(),
            ZhangExactVector(dimension));
        for (std::size_t row = 0; row < expanded.size(); ++row)
            for (std::size_t column = 0; column < width; ++column)
                expanded[row][domainFrame->columns[column]] =
                    affine.deterministicBasis[row][column];
        const auto cached = zhangExactRowHermiteNormalForm(
            expanded, affine.deterministicValues);
        return cached.consistent && cached.basis == base.basis &&
            cached.values == base.values;
    };
    if (cacheMatchesCompleteDomain())
        out.frame = zhangR49CompileTargetOnDomain(targets, *domainFrame, dimension);
    else
        out.frame = zhangR47CompileProductSearchFrame(
            targets, held, heldValues, dimension);
    if (!out.frame.valid && out.frame.reason == "TARGET_OUTSIDE_CACHED_DOMAIN")
        out.frame = zhangR47CompileProductSearchFrame(
            targets, held, heldValues, dimension);
    if (!out.frame.valid)
    {
        out.status = out.frame.reason;
        return out;
    }
    out.rank = out.frame.searchRank;
    if (out.rank == 0)
    {
        out.status = "ALREADY_PROVEN";
        return out;
    }
    if (out.rank < 0 || out.rank > static_cast<int>(targets.size()) ||
        out.frame.projector.size() != static_cast<std::size_t>(out.rank) ||
        out.frame.offsets.size() != static_cast<std::size_t>(out.rank) ||
        !zhangExactRectangularMatrix(out.frame.projector, dimension))
    {
        out.status = "TARGETED_CYCLE_INVALID_INTEGER_IMAGE";
        return out;
    }

    std::unique_ptr<ZhangR48MarginalWorkspace> localWorkspace;
    if (!sharedWorkspace)
    {
        localWorkspace = std::make_unique<ZhangR48MarginalWorkspace>(
            rootMean, rootCovariance);
        sharedWorkspace = localWorkspace.get();
    }
    // The workspace API has no posterior identity token. Enforce the part of
    // its immutable-root contract that is directly observable here; its owner
    // remains responsible for supplying the matching full covariance.
    if (!sharedWorkspace->valid || sharedWorkspace->mean.size() != dimension ||
        !(sharedWorkspace->mean.array() == rootMean.array()).all())
    {
        out.status = "TARGETED_CYCLE_WORKSPACE_ROOT_MISMATCH";
        return out;
    }
    const auto projector = zhangR48Numeric(out.frame.projector, dimension);
    const auto heldNumeric = zhangR48Numeric(held, dimension);
    const auto heldNumericValues = zhangExactRowToDouble(heldValues);
    const auto marginal = sharedWorkspace->project(
        projector, heldNumeric, heldNumericValues);
    if (!marginal.valid)
    {
        out.status = marginal.failureReason;
        return out;
    }
    out.mean = marginal.mean + zhangExactRowToDouble(out.frame.offsets);
    out.covariance = marginal.covariance;
    if (!out.mean.allFinite() || !out.covariance.allFinite())
    {
        out.status = "TARGETED_CYCLE_NONFINITE_IMAGE";
        return out;
    }

    out.spent = allocation;
    const auto proposal = search(out.mean, out.covariance, allocation, true);
    const auto fixedRank = proposal.rows.size();
    if (!proposal.valid || fixedRank == 0 ||
        fixedRank > static_cast<std::size_t>(out.rank) ||
        (!allowPartial && fixedRank != static_cast<std::size_t>(out.rank)) ||
        proposal.values.size() != fixedRank ||
        !zhangExactRectangularMatrix(proposal.rows, out.rank) ||
        zhangExactRowHermiteNormalForm(proposal.rows).basis.size() != fixedRank ||
        !std::isfinite(proposal.failureProbability) ||
        proposal.failureProbability < 0 || proposal.failureProbability > 1 ||
        zhangRatioStatisticalReject(proposal.failureProbability > allocation))
    {
        out.status = "TARGETED_CYCLE_SEARCH_REJECTED";
        return out;
    }
    out.rows = zhangExactMultiply(proposal.rows, out.frame.projector);
    out.values = proposal.values;
    const auto offsets = zhangExactMatrixTimesColumn(
        proposal.rows, out.frame.offsets);
    for (std::size_t row = 0; row < out.values.size(); ++row)
        out.values[row] -= offsets[row];

    auto combined = held;
    auto combinedValues = heldValues;
    combined.insert(combined.end(), out.rows.begin(), out.rows.end());
    combinedValues.insert(combinedValues.end(), out.values.begin(), out.values.end());
    const auto united = zhangExactRowHermiteNormalForm(combined, combinedValues);
    if (!united.consistent || !zhangExactAffineIntegerQuotient(
            united.basis, united.values, dimension,
            ZhangExactQuotientWork::FEASIBILITY_ONLY).valid)
    {
        out.status = "AFFINE_INFEASIBLE";
        return out;
    }
    // Every accepted image row must add an independent constraint on the
    // complete current domain. Same-rank/entailed answers never advance it.
    if (united.basis.size() != base.basis.size() + fixedRank)
    {
        out.status = "TARGETED_CYCLE_NO_INDEPENDENT_RANK_GAIN";
        return out;
    }

    const auto jointNumeric = zhangR48Numeric(united.basis, dimension);
    out.nis = assessZhangIntegerCandidateNis(
        zhangExactRowToDouble(united.values) - jointNumeric * rootMean,
        jointNumeric * rootCovariance * jointNumeric.transpose(), alpha);
    // In ratio-only mode the assessor above performs only finite/PSD and
    // deterministic-null-residual validation: it computes neither a quadratic
    // NIS nor a chi-square threshold. Avoid the two diagnostic-only extra
    // covariance constructions/decompositions in that mode, and leave their
    // receipts explicitly unevaluated rather than manufacturing zero NIS.
    if (zhangRatioOnly())
    {
        out.baseNis.status = "NOT_EVALUATED_RATIO_ONLY";
        out.incrementNis.status = "NOT_EVALUATED_RATIO_ONLY";
    }
    else
    {
        if (held.empty())
        {
            out.baseNis.valid = true;
            out.baseNis.nis = 0;
            out.baseNis.rank = 0;
            out.baseNis.status = "EMPTY_CONDITIONER";
        }
        else
            out.baseNis = assessZhangIntegerCandidateNis(
                heldNumericValues - heldNumeric * rootMean,
                heldNumeric * rootCovariance * heldNumeric.transpose(), alpha);
        const auto increment = sharedWorkspace->project(
            zhangR48Numeric(out.rows, dimension), heldNumeric, heldNumericValues);
        if (increment.valid)
            out.incrementNis = assessZhangIntegerCandidateNis(
                zhangExactRowToDouble(out.values) - increment.mean,
                increment.covariance, alpha);
        if (out.nis.valid && out.baseNis.valid && out.incrementNis.valid)
            out.schurIdentityError = std::abs(
                out.nis.nis - out.baseNis.nis - out.incrementNis.nis);
    }
    out.accepted = zhangIntegerCandidateAdmissible(out.nis);
    out.status = out.accepted
        ? (allowPartial && fixedRank < static_cast<std::size_t>(out.rank)
            ? "PARTIAL_INTEGER_PROGRESS" : "TARGETED_CYCLE_ACCEPTED")
        : out.nis.status;
    return out;
}
