#pragma once
#include "common/zhangQuotientIntegerLattice.hpp"

/** The complete named catalogue and the independently searchable lattice are
 * separate objects. searchNamedRows*z = searchPosteriorRows*a + searchOffsets.
 * No value is assigned to a missing coordinate. */
struct ZhangWholeProductLattice
{
    bool valid = false;
    bool primitiveSearchImage = false;
    int targetIndependentRank = 0;
    int directlyAvailableRank = 0;
    int wholeAvailableRank = 0;
    ZhangExactMatrix namedExpressions;
    ZhangExactVector namedOffsets;
    std::vector<bool> availableColumns;
    ZhangExactMatrix searchNamedRows;
    ZhangExactMatrix searchPosteriorRows;
    ZhangExactVector searchOffsets;
    // Structural identities are equations between expressions, not AR proofs.
    ZhangExactMatrix structuralIdentityRows;
    ZhangExactVector structuralIdentityValues;
    std::vector<std::string> namedStatus;
    ZhangExactMatrix namedRecovery;
    ZhangExactVector namedRecoveryOffsets;
    std::vector<bool> namedRecoverable;
    ZhangExactVector searchSmithInvariants;
    std::string failureReason = "NOT_COMPILED";
};

inline ZhangWholeProductLattice zhangCompileWholeProductLattice(
    const ZhangExactMatrix& named, const ZhangExactVector& offsets,
    const std::vector<bool>& available)
{
    ZhangWholeProductLattice out;
    out.namedExpressions = named;
    out.namedOffsets = offsets;
    out.availableColumns = available;
    if (named.size() != offsets.size() ||
        !zhangExactRectangularMatrix(named, available.size()))
    {
        out.failureReason = "NAMED_CATALOGUE_DIMENSION_MISMATCH";
        return out;
    }
    ZhangExactMatrix transpose(available.size(), ZhangExactVector(named.size()));
    for (std::size_t n = 0; n < named.size(); ++n)
        for (std::size_t c = 0; c < available.size(); ++c) transpose[c][n] = named[n][c];
    out.structuralIdentityRows = zhangExactIntegerKernel(transpose, named.size());
    for (const auto& identity : out.structuralIdentityRows)
    {
        ZhangExactInteger value = 0;
        for (std::size_t n = 0; n < named.size(); ++n) value += identity[n]*offsets[n];
        out.structuralIdentityValues.push_back(value);
    }
    out.targetIndependentRank = zhangExactRowHermiteNormalForm(named).basis.size();
    ZhangExactMatrix direct;
    for (const auto& row : named)
    {
        bool supported = true;
        for (std::size_t c = 0; c < available.size(); ++c)
            supported &= available[c] || row[c] == 0;
        if (supported) direct.push_back(row);
    }
    out.directlyAvailableRank = zhangExactRowHermiteNormalForm(direct).basis.size();
    // Reuse the saturated integer kernel and tracked unimodular HNF. Values
    // here are zero because named targets are expressions, not equations.
    const auto surviving = zhangExactSurvivingLattice(
        named, ZhangExactVector(named.size()), available, true);
    if (!surviving.consistent)
    {
        out.failureReason = "WHOLE_LATTICE_ELIMINATION_FAILED";
        return out;
    }
    out.searchNamedRows = surviving.rowTransform;
    out.searchPosteriorRows = surviving.basis;
    out.wholeAvailableRank = surviving.basis.size();
    for (const auto& combination : out.searchNamedRows)
    {
        ZhangExactInteger offset = 0;
        for (std::size_t n = 0; n < named.size(); ++n) offset += combination[n]*offsets[n];
        out.searchOffsets.push_back(offset);
    }
    const std::size_t compactDimension = std::count(available.begin(), available.end(), true);
    const auto smith = zhangIntegerRowLatticeContains(
        out.searchPosteriorRows, ZhangExactVector(compactDimension));
    out.searchSmithInvariants = smith.smithInvariants;
    out.primitiveSearchImage = std::all_of(smith.smithInvariants.begin(),
        smith.smithInvariants.end(), [](const auto& d) { return zhangExactAbs(d) == 1; });
    // Verify every missing coefficient cancels in the actual selected basis.
    const auto rebuilt = zhangExactMultiply(out.searchNamedRows, named);
    for (std::size_t r = 0; r < rebuilt.size(); ++r)
    {
        std::size_t compact = 0;
        for (std::size_t c = 0; c < available.size(); ++c)
        {
            if (!available[c] && rebuilt[r][c] != 0)
            {
                out.failureReason = "MISSING_COLUMN_CANCELLATION_NOT_EXACT";
                return out;
            }
            if (available[c] && rebuilt[r][c] != out.searchPosteriorRows[r][compact++])
            {
                out.failureReason = "SEARCH_BASIS_PULLBACK_NOT_EXACT";
                return out;
            }
        }
    }
    ZhangExactMatrix selectedNamed;
    for (std::size_t n = 0; n < named.size(); ++n)
    {
        bool supported = true;
        ZhangExactVector compact;
        for (std::size_t c = 0; c < available.size(); ++c)
            if (available[c]) compact.push_back(named[n][c]);
            else supported &= named[n][c] == 0;
        auto membership = zhangIntegerRowLatticeContains(out.searchPosteriorRows, compact);
        const bool recoverable = supported && membership.contained;
        ZhangExactVector recovery(out.wholeAvailableRank);
        ZhangExactInteger shift = offsets[n];
        if (recoverable)
        {
            recovery = membership.combination;
            if (recovery.size() != out.searchOffsets.size())
            {
                // The empty search lattice still represents zero targets.
                if (!recovery.empty() || !out.searchOffsets.empty())
                { out.failureReason = "NAMED_RECOVERY_DIMENSION_MISMATCH"; return out; }
            }
            for (std::size_t r = 0; r < recovery.size(); ++r)
                shift -= recovery[r]*out.searchOffsets[r];
        }
        out.namedRecovery.push_back(std::move(recovery));
        out.namedRecoveryOffsets.push_back(shift);
        out.namedRecoverable.push_back(recoverable);
        auto augmented = selectedNamed;
        augmented.push_back(named[n]);
        const bool independent = zhangExactRowHermiteNormalForm(augmented).basis.size() > selectedNamed.size();
        if (recoverable && independent) selectedNamed = zhangExactRowHermiteNormalForm(augmented).basis;
        bool inCombination = false;
        for (const auto& combination : out.searchNamedRows) inCombination |= combination[n] != 0;
        out.namedStatus.push_back(recoverable
            ? (independent ? "REPRESENTABLE_INDEPENDENT" : "REPRESENTABLE_DEPENDENT")
            : (inCombination ? "AVAILABLE_ONLY_AS_COMBINATION" : "NOT_CURRENTLY_REPRESENTABLE"));
    }
    out.valid = true;
    out.failureReason = out.primitiveSearchImage
        ? "EXACT_WHOLE_PRODUCT_LATTICE" : "EXACT_PRODUCT_IMAGE_WITH_DIVISIBILITY_CONSTRAINTS";
    return out;
}
