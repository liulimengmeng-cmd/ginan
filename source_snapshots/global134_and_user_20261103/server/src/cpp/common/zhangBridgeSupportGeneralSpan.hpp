#pragma once

#include "common/zhangIntegerAudit.hpp"
#include "common/zhangPhaseTimer.hpp"

using ZhangBridgeSpanRational = boost::multiprecision::cpp_rational;
using ZhangBridgeSpanRow = std::vector<ZhangBridgeSpanRational>;
using ZhangBridgeSpanMatrix = std::vector<ZhangBridgeSpanRow>;

// Mathematical representability only. The caller must separately prove that
// CURRENT H is integer feasible and each S row has its actual whole numeric
// binding on one owner/chart. Include nuisance columns explicitly, or prove
// their cancellation before supplying S. This does not certify an integer
// combination, change an integer image, or grant publication authority.
struct ZhangBridgeGeneralSpanAnalysis
{
    bool valid = false, representable = false;
    std::string reason = "NOT_EVALUATED";
    std::size_t missingRank = 0;
    std::vector<ZhangBridgeSpanRow> currentCombinations;
    std::vector<ZhangBridgeSpanRow> numericCombinations;
    std::vector<std::size_t> unresolvedTargets;
    // H*v=0 and S*v=0, but this target*v!=0. Denominators are cleared;
    // nuisance coordinates remain explicitly present, never set to zero.
    std::vector<ZhangExactVector> counterexamples;
};

inline ZhangBridgeGeneralSpanAnalysis analyzeZhangBridgeSupportOnActualSpan(
    const ZhangExactMatrix& currentH,
    const ZhangBridgeSpanMatrix& actualS,
    const ZhangBridgeSpanMatrix& targets,
    std::size_t dimension)
{
    ZhangPhaseTimer timer("BRIDGE_SUPPORT_GENERAL_SPAN");
    ZhangBridgeGeneralSpanAnalysis out;
    auto rectangular = [dimension](const auto& matrix)
    {
        return std::all_of(matrix.begin(), matrix.end(),
            [dimension](const auto& row) { return row.size() == dimension; });
    };
    if (!rectangular(currentH) || !rectangular(actualS) || !rectangular(targets))
    { out.reason = "BRIDGE_GENERAL_SPAN_SHAPE_INVALID"; return out; }

    using Rational = ZhangBridgeSpanRational;
    using Sparse = std::map<std::size_t, Rational>;
    struct BasisRow { Sparse coefficients, provenance; };
    using Basis = std::map<std::size_t, BasisRow>;
    Basis basis;
    auto sparse = [](const auto& row)
    {
        Sparse result;
        for (std::size_t column = 0; column < row.size(); ++column)
            if (row[column] != 0) result.emplace(column, Rational(row[column]));
        return result;
    };
    auto add = [](Sparse& destination, const Sparse& source, const Rational& scale)
    {
        if (scale == 0) return;
        for (const auto& [column, value] : source)
        {
            zhangExactWorkCheckpoint("BRIDGE_GENERAL_SPAN_SPARSE_ADD");
            const Rational next = destination[column] + scale * value;
            if (next == 0) destination.erase(column);
            else destination[column] = next;
        }
    };
    auto insert = [&](Basis& into, Sparse row, std::size_t sourceIndex)
    {
        Sparse provenance{{sourceIndex, 1}};
        for (const auto& [pivot, source] : into)
            if (const auto found = row.find(pivot); found != row.end())
            {
                const Rational scale = -found->second;
                add(row, source.coefficients, scale);
                add(provenance, source.provenance, scale);
            }
        if (row.empty()) return false;
        const auto pivot = row.begin()->first;
        const Rational divisor = row.begin()->second;
        for (auto& [column, value] : row) value /= divisor;
        for (auto& [index, value] : provenance) value /= divisor;
        into.emplace(pivot, BasisRow{std::move(row), std::move(provenance)});
        return true;
    };
    for (std::size_t row = 0; row < currentH.size(); ++row)
    {
        zhangExactWorkCheckpoint("BRIDGE_GENERAL_SPAN_CURRENT_ROW");
        insert(basis, sparse(currentH[row]), row);
    }
    for (std::size_t row = 0; row < actualS.size(); ++row)
    {
        zhangExactWorkCheckpoint("BRIDGE_GENERAL_SPAN_NUMERIC_ROW");
        insert(basis, sparse(actualS[row]), currentH.size() + row);
    }
    auto augmented = basis;
    out.currentCombinations.resize(targets.size());
    out.numericCombinations.resize(targets.size());
    for (std::size_t targetIndex = 0; targetIndex < targets.size(); ++targetIndex)
    {
        zhangExactWorkCheckpoint("BRIDGE_GENERAL_SPAN_TARGET");
        Sparse remaining = sparse(targets[targetIndex]), witness;
        for (const auto& [pivot, source] : basis)
            if (const auto found = remaining.find(pivot); found != remaining.end())
            {
                const Rational factor = found->second;
                add(remaining, source.coefficients, -factor);
                add(witness, source.provenance, factor);
            }
        if (remaining.empty())
        {
            auto& h = out.currentCombinations[targetIndex];
            auto& s = out.numericCombinations[targetIndex];
            h.assign(currentH.size(), 0); s.assign(actualS.size(), 0);
            for (const auto& [index, value] : witness)
                if (index < currentH.size()) h[index] = value;
                else s[index - currentH.size()] = value;
            // Independently audit the witness against original supplied rows.
            for (std::size_t column = 0; column < dimension; ++column)
            {
                Rational reconstructed = 0;
                for (std::size_t row = 0; row < currentH.size(); ++row)
                    reconstructed += h[row] * Rational(currentH[row][column]);
                for (std::size_t row = 0; row < actualS.size(); ++row)
                    reconstructed += s[row] * actualS[row][column];
                if (reconstructed != targets[targetIndex][column])
                { out.reason = "BRIDGE_GENERAL_SPAN_COMBINATION_AUDIT_FAILED"; return out; }
            }
            continue;
        }
        out.unresolvedTargets.push_back(targetIndex);
        if (insert(augmented, sparse(targets[targetIndex]),
                currentH.size() + actualS.size() + targetIndex)) ++out.missingRank;

        Sparse direction{{remaining.begin()->first, 1}};
        for (auto item = basis.rbegin(); item != basis.rend(); ++item)
        {
            Rational value = 0;
            for (const auto& [column, coefficient] : item->second.coefficients)
                if (column != item->first && direction.contains(column))
                    value -= coefficient * direction.at(column);
            if (value == 0) direction.erase(item->first);
            else direction[item->first] = value;
        }
        ZhangExactInteger scale = 1;
        for (const auto& [column, value] : direction)
        {
            ZhangExactInteger a = scale, b = denominator(value);
            while (b != 0) { ZhangExactInteger remainder = a % b; a = b; b = remainder; }
            scale = scale / a * denominator(value);
        }
        ZhangExactVector integerDirection(dimension);
        for (const auto& [column, value] : direction)
        {
            const Rational integerValue = value * scale;
            integerDirection[column] = numerator(integerValue);
        }
        auto dot = [&](const auto& row)
        {
            Rational value = 0;
            for (std::size_t column = 0; column < dimension; ++column)
                value += Rational(row[column]) * Rational(integerDirection[column]);
            return value;
        };
        if (std::any_of(currentH.begin(), currentH.end(),
                [&](const auto& row) { return dot(row) != 0; }) ||
            std::any_of(actualS.begin(), actualS.end(),
                [&](const auto& row) { return dot(row) != 0; }) ||
            dot(targets[targetIndex]) == 0)
        { out.reason = "BRIDGE_GENERAL_SPAN_NULL_WITNESS_FAILED"; return out; }
        out.counterexamples.push_back(std::move(integerDirection));
    }
    out.valid = true;
    out.representable = out.unresolvedTargets.empty();
    out.reason = out.representable ? "WHOLE_TARGET_IN_CURRENT_AND_ACTUAL_NUMERIC_SPAN" :
        "WHOLE_TARGET_OUTSIDE_CURRENT_AND_ACTUAL_NUMERIC_SPAN";
    return out;
}
