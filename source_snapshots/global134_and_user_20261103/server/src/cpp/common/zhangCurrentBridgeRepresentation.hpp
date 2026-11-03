#pragma once

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "common/zhangR51PhysicalImage.hpp"

/** Change only the numerical representation of an EXISTING integer-image
 * coordinate q=W*a+offset, using accepted CURRENT constraints H*a=z.
 * The original image generators, particular target, physical targets and
 * ordering are not changed. Rational W' is not an unconstrained integer
 * lattice: q remains the integer coordinate of the original exact image.
 * No missing physical variable is assigned zero or a historical integer. */
struct ZhangCurrentBridgeRepresentation
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    std::vector<ZhangR51RationalRow> projector;
    ZhangR51RationalRow offsets;
    std::set<std::size_t> unresolvedColumns;
};

inline ZhangCurrentBridgeRepresentation zhangRewriteCurrentBridgeRepresentation(
    const std::vector<ZhangR51RationalRow>& originalProjector,
    const ZhangR51RationalRow& originalOffsets,
    const ZhangExactMatrix& currentConditions,
    const ZhangExactVector& currentValues,
    const std::set<std::size_t>& unavailableColumns)
{
    ZhangCurrentBridgeRepresentation result;
    if (originalProjector.empty() || originalProjector.size() != originalOffsets.size() ||
        currentConditions.size() != currentValues.size())
    {
        result.reason = "CURRENT_BRIDGE_REWRITE_DIMENSION_INVALID";
        return result;
    }
    const std::size_t dimension = originalProjector.front().size();
    if (dimension == 0 || !zhangExactRectangularMatrix(currentConditions, dimension))
    {
        result.reason = "CURRENT_BRIDGE_REWRITE_DIMENSION_INVALID";
        return result;
    }
    for (const auto& row : originalProjector)
        if (row.size() != dimension)
        {
            result.reason = "CURRENT_BRIDGE_REWRITE_DIMENSION_INVALID";
            return result;
        }
    for (const auto column : unavailableColumns)
        if (column >= dimension)
        {
            result.reason = "CURRENT_BRIDGE_REWRITE_COLUMN_INVALID";
            return result;
        }

    using SparseRow = std::map<std::size_t, ZhangR51Rational>;
    struct Constraint
    {
        SparseRow coefficients;
        ZhangR51Rational value;
        SparseRow currentWitness;
    };
    auto addScaled = [](SparseRow& destination, const SparseRow& source,
                        const ZhangR51Rational& scale)
    {
        if (scale == 0) return;
        for (const auto& [column, coefficient] : source)
        {
            zhangExactWorkCheckpoint("CURRENT_BRIDGE_REWRITE_SPARSE_ADD");
            const ZhangR51Rational value = destination[column] + scale * coefficient;
            if (value == 0) destination.erase(column);
            else destination[column] = value;
        }
    };
    // Only rows touching unavailable columns can participate in their
    // elimination. Keep their complete live support, not a truncated H.
    std::vector<Constraint> conditions;
    for (std::size_t r = 0; r < currentConditions.size(); ++r)
    {
        zhangExactWorkCheckpoint("CURRENT_BRIDGE_REWRITE_SELECT_CONDITION");
        bool relevant = false;
        for (const auto column : unavailableColumns)
            relevant = relevant || currentConditions[r][column] != 0;
        if (!relevant) continue;
        Constraint row;
        for (std::size_t c = 0; c < dimension; ++c)
        {
            zhangExactWorkCheckpoint("CURRENT_BRIDGE_REWRITE_COPY_CONDITION");
            if (currentConditions[r][c] != 0)
                row.coefficients.emplace(c, currentConditions[r][c]);
        }
        row.value = currentValues[r];
        row.currentWitness.emplace(r, 1);
        conditions.push_back(std::move(row));
    }
    std::vector<std::pair<std::size_t, std::size_t>> pivots;
    std::size_t next = 0;
    for (const auto column : unavailableColumns)
    {
        zhangExactWorkCheckpoint("CURRENT_BRIDGE_REWRITE_PIVOT");
        std::size_t pivot = next;
        while (pivot < conditions.size() &&
               !conditions[pivot].coefficients.contains(column)) ++pivot;
        if (pivot == conditions.size()) continue;
        std::swap(conditions[next], conditions[pivot]);
        auto& row = conditions[next];
        const ZhangR51Rational divisor = row.coefficients.at(column);
        for (auto& [c, coefficient] : row.coefficients) coefficient /= divisor;
        for (auto& [r, coefficient] : row.currentWitness) coefficient /= divisor;
        row.value /= divisor;
        for (std::size_t r = next + 1; r < conditions.size(); ++r)
        {
            const auto found = conditions[r].coefficients.find(column);
            if (found == conditions[r].coefficients.end()) continue;
            const ZhangR51Rational factor = found->second;
            addScaled(conditions[r].coefficients, row.coefficients, -factor);
            addScaled(conditions[r].currentWitness, row.currentWitness, -factor);
            conditions[r].value -= factor * row.value;
        }
        pivots.emplace_back(column, next++);
    }

    result.projector = originalProjector;
    result.offsets = originalOffsets;
    for (std::size_t r = 0; r < originalProjector.size(); ++r)
    {
        SparseRow row, witness;
        for (std::size_t c = 0; c < dimension; ++c)
            if (originalProjector[r][c] != 0) row.emplace(c, originalProjector[r][c]);
        for (const auto& [column, index] : pivots)
        {
            const auto found = row.find(column);
            if (found == row.end()) continue;
            const ZhangR51Rational factor = found->second;
            addScaled(row, conditions[index].coefficients, -factor);
            addScaled(witness, conditions[index].currentWitness, factor);
            result.offsets[r] += factor * conditions[index].value;
        }
        for (const auto column : unavailableColumns)
            if (row.contains(column)) result.unresolvedColumns.insert(column);
        result.projector[r].assign(dimension, 0);
        for (const auto& [column, coefficient] : row)
            result.projector[r][column] = coefficient;

        // Independently audit the exact identity against the ORIGINAL current
        // rows: W-W'=alpha H and offset'-offset=alpha z. Fractional alpha is
        // permitted; it changes representation, never the integer image.
        ZhangR51RationalRow difference(dimension);
        ZhangR51Rational value = 0;
        for (const auto& [condition, coefficient] : witness)
        {
            value += coefficient * ZhangR51Rational(currentValues[condition]);
            for (std::size_t c = 0; c < dimension; ++c)
            {
                zhangExactWorkCheckpoint("CURRENT_BRIDGE_REWRITE_EXACT_AUDIT");
                if (currentConditions[condition][c] != 0)
                    difference[c] += coefficient * ZhangR51Rational(currentConditions[condition][c]);
            }
        }
        for (std::size_t c = 0; c < dimension; ++c)
            if (originalProjector[r][c] - result.projector[r][c] != difference[c])
            {
                result.reason = "CURRENT_BRIDGE_REWRITE_EXACT_IDENTITY_FAILED";
                return result;
            }
        if (result.offsets[r] - originalOffsets[r] != value)
        {
            result.reason = "CURRENT_BRIDGE_REWRITE_EXACT_OFFSET_FAILED";
            return result;
        }
    }
    if (!result.unresolvedColumns.empty())
    {
        result.reason = "CURRENT_CONSTRAINTS_CANNOT_ELIMINATE_UNAVAILABLE_ARC";
        return result;
    }
    result.valid = true;
    result.reason = unavailableColumns.empty() ? "CURRENT_BRIDGE_REPRESENTATION_UNCHANGED"
        : "CURRENT_CONSTRAINTS_ELIMINATED_UNAVAILABLE_ARCS";
    return result;
}
