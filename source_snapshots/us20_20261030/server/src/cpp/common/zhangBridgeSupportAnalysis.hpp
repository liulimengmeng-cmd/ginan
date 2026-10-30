#pragma once

#include "common/zhangCurrentBridgeRepresentation.hpp"
#include <algorithm>

/** A representability diagnostic on an ALREADY VERIFIED feasible CURRENT
 * integer domain. It is neither a feasibility check nor an integer row-lattice
 * authorization. A rational row rewrite must not saturate the integer image. */
struct ZhangBridgeSupportAnalysis
{
    bool valid = false, representable = false;
    std::size_t missingRank = 0;
    std::vector<std::size_t> unavailableColumns;
    // Full physical-column integer vectors. Each is independently audited:
    // H*v=0, live(v)=0, and at least one target changes.
    std::vector<ZhangExactVector> counterexamples;
    std::string reason = "NOT_EVALUATED";
};

inline ZhangBridgeSupportAnalysis zhangAnalyzeBridgeSupport(
    const ZhangExactMatrix& currentRows, const ZhangExactMatrix& targets,
    const std::set<std::size_t>& unavailable, std::size_t dimension)
{
    ZhangPhaseTimer timer("BRIDGE_SUPPORT_PREFLIGHT");
    ZhangBridgeSupportAnalysis out;
    out.unavailableColumns.assign(unavailable.begin(), unavailable.end());
    if (!zhangExactRectangularMatrix(currentRows, dimension) ||
        !zhangExactRectangularMatrix(targets, dimension) ||
        (!unavailable.empty() && *unavailable.rbegin() >= dimension))
    { out.reason = "BRIDGE_SUPPORT_SHAPE_INVALID"; return out; }
    using Row = std::map<std::size_t, ZhangR51Rational>;
    std::map<std::size_t, Row> basis;
    auto add = [](Row& row, const Row& source, const ZhangR51Rational& factor) {
        for (const auto& [c, value] : source)
        {
            zhangExactWorkCheckpoint("BRIDGE_SUPPORT_SPARSE_REDUCTION");
            const ZhangR51Rational next = row[c] + factor * value;
            if (next == 0) row.erase(c); else row[c] = next;
        }
    };
    auto reduce = [&](Row row, const std::map<std::size_t, Row>& against) {
        for (const auto& [pivot, source] : against)
            if (const auto it = row.find(pivot); it != row.end())
            { const ZhangR51Rational factor = -it->second; add(row, source, factor); }
        return row;
    };
    auto restrict = [&](const ZhangExactVector& original) {
        Row row;
        for (const auto c : unavailable) if (original[c] != 0) row[c] = original[c];
        return row;
    };
    auto insert = [](std::map<std::size_t, Row>& into, Row row) {
        if (row.empty()) return;
        const auto pivot = row.begin()->first;
        const ZhangR51Rational divisor = row.begin()->second;
        for (auto& [c, value] : row) value /= divisor;
        into.emplace(pivot, std::move(row));
    };
    for (const auto& original : currentRows)
    {
        zhangExactWorkCheckpoint("BRIDGE_SUPPORT_CONDITION_ROW");
        insert(basis, reduce(restrict(original), basis));
    }
    auto augmented = basis;
    for (const auto& target : targets)
    {
        const auto reduced = reduce(restrict(target), basis);
        auto independent = reduce(restrict(target), augmented);
        if (!independent.empty()) { ++out.missingRank; insert(augmented, std::move(independent)); }
        if (reduced.empty()) continue;
        // A nonpivot column in the residual determines an exact null direction.
        Row vector;
        vector[reduced.begin()->first] = 1;
        for (auto it = basis.rbegin(); it != basis.rend(); ++it)
        {
            ZhangR51Rational value = 0;
            for (const auto& [c, coefficient] : it->second)
                if (c != it->first && vector.contains(c)) value -= coefficient * vector[c];
            if (value != 0) vector[it->first] = value;
        }
        ZhangExactInteger scale = 1;
        for (const auto& [c, value] : vector)
        {
            ZhangExactInteger a = scale, b = denominator(value);
            while (b != 0) { ZhangExactInteger r = a % b; a = b; b = r; }
            scale = scale / a * denominator(value);
        }
        ZhangExactVector witness(dimension);
        for (const auto& [c, value] : vector)
        { const ZhangR51Rational integerValue=value*scale; witness[c]=numerator(integerValue); }
        for (const auto& original : currentRows)
        {
            ZhangExactInteger value = 0;
            for (const auto c : unavailable) value += original[c] * witness[c];
            if (value != 0) { out.reason = "BRIDGE_SUPPORT_NULL_WITNESS_FAILED"; return out; }
        }
        ZhangExactInteger change = 0;
        for (const auto c : unavailable) change += target[c] * witness[c];
        if (change == 0) { out.reason = "BRIDGE_SUPPORT_TARGET_WITNESS_FAILED"; return out; }
        out.counterexamples.push_back(std::move(witness));
    }
    out.valid = true;
    out.representable = out.missingRank == 0;
    out.reason = out.representable ? "CURRENT_SUPPORT_REPRESENTABLE" :
        "CURRENT_CONSTRAINTS_CANNOT_ELIMINATE_UNAVAILABLE_ARC";
    return out;
}
