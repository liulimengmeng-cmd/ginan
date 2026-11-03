#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "common/zhangProductRelationBasis.hpp"

/** An explicit request to project a genuinely conditioned native posterior.
 * It is NOT a pair-integer certificate and never turns 0=0 into an AR row.
 * The producer must separately validate the immutable current network domain,
 * its physical arc bindings and its decision parents. */
struct ZhangZeroDatumNativeRequest
{
    bool requested = false;
    std::vector<SatSys> satellites;
    int conditioningRank = 0;
    bool singleSBasisNative = false;
    std::string committedSBasisIdentity;
};

inline bool zhangZeroDatumNativeStructureValid(
    const ZhangProductRelationBasis& first,
    const ZhangProductRelationBasis& second,
    int dimension,
    std::vector<SatSys>& satellites,
    std::string& reason)
{
    satellites.clear();
    auto fail = [&](const char* why) { reason = why; return false; };
    if (dimension <= 0 || first.system == E_Sys::NONE ||
        first.system != second.system || first.observable == E_ObsCode::NONE ||
        second.observable == E_ObsCode::NONE || first.observable == second.observable ||
        first.referenceSatellite != second.referenceSatellite ||
        first.satellites != second.satellites || first.satellites.size() < 2 ||
        !std::is_sorted(first.satellites.begin(), first.satellites.end()) ||
        std::adjacent_find(first.satellites.begin(), first.satellites.end()) !=
            first.satellites.end() ||
        !zhangProductRelationSemanticOrderingMatches(first.namedRelations, second.namedRelations))
        return fail("ZERO_DATUM_NATIVE_CATALOGUE_OR_SIGNAL_MISMATCH");
    if (first.singleSBasisNative != second.singleSBasisNative ||
        first.committedSBasisIdentity != second.committedSBasisIdentity ||
        (first.singleSBasisNative && first.committedSBasisIdentity.empty()))
        return fail("ZERO_DATUM_NATIVE_COMMITTED_S_SNAPSHOT_MISMATCH");
    const std::set<SatSys> catalogue(first.satellites.begin(), first.satellites.end());
    if (!catalogue.contains(first.referenceSatellite))
        return fail("ZERO_DATUM_NATIVE_REFERENCE_OUTSIDE_CATALOGUE");
    for (const auto& satellite : catalogue)
        if (satellite.prn <= 0 || satellite.sys != first.system)
            return fail("ZERO_DATUM_NATIVE_INVALID_CATALOGUE_SATELLITE");
    const int namedCount = static_cast<int>(catalogue.size()) - 1;
    auto allZero = [](const auto& values)
    {
        return std::all_of(values.begin(), values.end(),
            [](const auto& value) { return value == 0; });
    };
    for (const auto* basis : {&first, &second})
    {
        if (!basis->valid || !basis->structuralMapValid || !basis->zeroDatumCorrection ||
            basis->exactRank != 0 || basis->datumCorrectionRank != 0 ||
            !basis->physicalExpansionValid || !basis->nuisanceOrthogonal ||
            !basis->networkClosureExactZero || !basis->networkLatticeContained ||
            !basis->admissibleCompletionProven || !basis->primitive ||
            basis->saturationIndex != 1 || !basis->canonicalTargetPrimitive ||
            basis->canonicalTargetRank != namedCount || basis->fullTargetRank != namedCount ||
            basis->namedRelationCount != namedCount ||
            basis->namedRelations.size() != static_cast<std::size_t>(namedCount) ||
            !basis->exactRowBasis.empty() || !basis->independentNamedIndices.empty() ||
            basis->mappableTargetRank != 0 || !basis->mappableNamedIndices.empty() ||
            !basis->unmappableNamedRelations.empty() || !basis->wholePosteriorRows.empty() ||
            basis->transform.rows() != 0 || basis->transform.cols() != dimension ||
            !basis->affineOffsets.empty() || !basis->wholeLattice.valid ||
            basis->wholeLattice.targetIndependentRank != 0 ||
            basis->wholeLattice.directlyAvailableRank != 0 ||
            basis->wholeLattice.wholeAvailableRank != 0 ||
            basis->wholeLattice.hasNonzeroSearchDirections ||
            !basis->wholeLattice.searchPosteriorRows.empty() ||
            !basis->wholeLattice.searchNamedRows.empty() ||
            !basis->wholeLattice.searchOffsets.empty())
            return fail("ZERO_DATUM_NATIVE_NOT_A_COMPLETE_ZERO_CORRECTION_MAP");
        const auto& whole = basis->wholeLattice;
        if (whole.namedExpressions.size() != static_cast<std::size_t>(namedCount) ||
            whole.namedOffsets.size() != static_cast<std::size_t>(namedCount) ||
            whole.namedRecoverable.size() != static_cast<std::size_t>(namedCount) ||
            whole.namedRecovery.size() != static_cast<std::size_t>(namedCount) ||
            whole.namedRecoveryOffsets.size() != static_cast<std::size_t>(namedCount) ||
            whole.availableColumns.size() != basis->currentChords.size() ||
            whole.structuralIdentityRows.size() != static_cast<std::size_t>(namedCount) ||
            whole.structuralIdentityValues.size() != static_cast<std::size_t>(namedCount) ||
            !allZero(whole.namedOffsets) || !allZero(whole.namedRecoveryOffsets) ||
            !allZero(whole.structuralIdentityValues) ||
            !std::all_of(whole.namedRecoverable.begin(), whole.namedRecoverable.end(),
                [](bool available) { return available; }))
            return fail("ZERO_DATUM_NATIVE_NAMED_CATALOGUE_INCOMPLETE");
        std::set<SatSys> namedSatellites;
        std::vector<int> expectedZeroIndices;
        for (int n = 0; n < namedCount; ++n)
        {
            const auto& named = basis->namedRelations[n];
            if (named.referenceSatellite != first.referenceSatellite ||
                named.satellite == first.referenceSatellite ||
                !catalogue.contains(named.satellite) ||
                !namedSatellites.insert(named.satellite).second ||
                named.currentCycleCoefficients.size() != basis->currentChords.size() ||
                !allZero(named.currentCycleCoefficients) || !allZero(named.nuisanceCoefficients) ||
                whole.namedExpressions[n].size() != basis->currentChords.size() ||
                !allZero(whole.namedExpressions[n]) || !whole.namedRecovery[n].empty() ||
                whole.structuralIdentityRows[n].size() != static_cast<std::size_t>(namedCount))
                return fail("ZERO_DATUM_NATIVE_NONZERO_OR_INCOMPLETE_NAMED_TARGET");
            for (const auto& [edge, coefficient] : named.physicalArcCoefficients)
                if (coefficient != 0)
                    return fail("ZERO_DATUM_NATIVE_NONZERO_PHYSICAL_TARGET");
            expectedZeroIndices.push_back(n);
        }
        auto actualZeroIndices = basis->structuralZeroNamedIndices;
        std::sort(actualZeroIndices.begin(), actualZeroIndices.end());
        if (actualZeroIndices != expectedZeroIndices)
            return fail("ZERO_DATUM_NATIVE_STRUCTURAL_ZERO_CATALOGUE_MISMATCH");
        // unmappableTargetRank is a legacy rank(G)-derived diagnostic, not a
        // count of absent native phase products. Do not reinterpret it here.
    }
    satellites = first.satellites;
    reason = "ZERO_DATUM_NATIVE_STRUCTURE_COMPLETE_NO_INTEGER_ROWS_CREATED";
    return true;
}

/** This verifies a receipt issued AFTER successful full-state conditioning;
 * it never issues one. Templates avoid a ProductRelationSolver/Evidence
 * include cycle. The caller supplies ACTUAL protected/result roots and state
 * order, not cached values copied from the receipt. This check does not
 * replace the producer's current physical/evidence/decision-parent contract. */
template<class Receipt, class Identity, class Recipe>
bool zhangZeroDatumNativePosteriorVerified(
    const ZhangZeroDatumNativeRequest& request,
    const Receipt& receipt,
    const Identity& actualIdentity,
    const Recipe& requiredRecipe,
    std::string& reason)
{
    auto fail = [&](const char* why) { reason = why; return false; };
    if (!request.requested || request.conditioningRank <= 0 ||
        request.satellites.size() < 2 || requiredRecipe.empty() ||
        !std::is_sorted(request.satellites.begin(), request.satellites.end()) ||
        std::adjacent_find(request.satellites.begin(), request.satellites.end()) !=
            request.satellites.end())
        return fail("ZERO_DATUM_NATIVE_NO_CURRENT_NONZERO_CONDITIONING_REQUEST");
    if (request.singleSBasisNative && request.committedSBasisIdentity.empty())
        return fail("SINGLE_S_BASIS_NATIVE_COMMITTED_SNAPSHOT_MISSING");
    const auto system = request.satellites.front().sys;
    for (const auto& satellite : request.satellites)
        if (satellite.prn <= 0 || satellite.sys != system || system == E_Sys::NONE)
            return fail("ZERO_DATUM_NATIVE_INVALID_REQUEST_SATELLITE");
    if (actualIdentity.runtimeId.empty() || actualIdentity.epochGpstSeconds <= 0 ||
        actualIdentity.posteriorId.empty() || actualIdentity.chartId.empty() ||
        actualIdentity.sourceRoot.empty() || actualIdentity.resultRoot.empty() ||
        actualIdentity.sourceStateOrder.empty() ||
        actualIdentity.sourceStateOrder != actualIdentity.resultStateOrder)
        return fail("ZERO_DATUM_NATIVE_INCOMPLETE_ACTUAL_POSTERIOR_IDENTITY");
    using Row = std::decay_t<decltype(requiredRecipe.front().rowCycles)>;
    std::map<Row, ZhangExactInteger> nonzeroRows;
    for (const auto& constraint : requiredRecipe)
    {
        Row row;
        for (const auto& [key, coefficient] : constraint.rowCycles)
        {
            if (!std::isfinite(coefficient))
                return fail("ZERO_DATUM_NATIVE_NONFINITE_CONDITIONING_ROW");
            if (coefficient != 0) row.emplace(key, coefficient);
        }
        if (row.empty())
        {
            if (constraint.value != 0)
                return fail("ZERO_DATUM_NATIVE_ZERO_ROW_CONFLICT");
            continue; // A structural identity contributes no AR authority.
        }
        const auto [where, inserted] = nonzeroRows.emplace(std::move(row), constraint.value);
        if (!inserted && where->second != constraint.value)
            return fail("ZERO_DATUM_NATIVE_DUPLICATE_RHS_CONFLICT");
    }
    if (nonzeroRows.empty() ||
        request.conditioningRank > static_cast<int>(nonzeroRows.size()))
        return fail("ZERO_DATUM_NATIVE_NONZERO_CONDITIONING_RECIPE_MISSING");
    if (!receipt.matches(actualIdentity, requiredRecipe, reason)) return false;
    reason = "ZERO_DATUM_NATIVE_ACTUAL_POSTERIOR_AND_COMPLETE_RECIPE_VERIFIED";
    return true;
}
