#pragma once

#include "common/algebra.hpp"
#include "common/zhangBridgeSupportGeneralSpan.hpp"
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <tuple>

using ZhangFloatCarrierRational = ZhangBridgeSpanRational;
using ZhangFloatCarrierSparseRow = std::map<std::size_t, ZhangFloatCarrierRational>;

enum class ZhangProductFloatTargetSpace
{
    VERSIONED_PHASE_EDGE,
    INTEGER_CYCLE_WITH_EXPLICIT_NUISANCE,
    SOURCE_KF_COORDINATE
};

// A FLOAT carrier is not a receiver/satellite physical arc and is never an
// integer coordinate. Its metadata must use a distinct tagged record from
// the legacy RAW_ARC directory, even if KF::PRODUCT_HANDOFF is the KF category.
enum class ZhangProductFloatRecordKind { RAW_ARC, TARGET_CARRIER, TARGET_BINDING };

struct ZhangProductFloatTargetScope
{
    std::string runtimeId, chartId, stateOrderId;
    std::uint64_t coordinateGeneration = 0;
    auto fields() const
    { return std::tie(runtimeId, chartId, stateOrderId, coordinateGeneration); }
    bool operator==(const ZhangProductFloatTargetScope& other) const
    { return fields() == other.fields(); }
};

struct ZhangProductFloatVersionedTerm
{
    ZhangGraphEdge edge;
    E_ObsCode observable = E_ObsCode::NONE;
    int arcVersion = -1;
    bool operator<(const ZhangProductFloatVersionedTerm& other) const
    { return std::tie(edge, observable, arcVersion) <
        std::tie(other.edge, other.observable, other.arcVersion); }
    bool operator==(const ZhangProductFloatVersionedTerm& other) const
    { return edge==other.edge && observable==other.observable && arcVersion==other.arcVersion; }
};

// Registration contract. The physical integer signature and the actual
// complete numeric phase row are distinct fields. In particular, an integer
// functional alone does not supply its absolute phase anchor/nuisance row.
// A request without a complete bound numeric row prevents compression; it is
// not omitted from the active consumer directory.
struct ZhangProductFloatTargetRequest
{
    std::string targetId, consumer, lifecycleId, unit;
    ZhangProductFloatTargetSpace symbolicSpace =
        ZhangProductFloatTargetSpace::SOURCE_KF_COORDINATE;
    std::map<ZhangProductFloatVersionedTerm, ZhangExactInteger> physicalCoefficients;
    std::map<std::string, ZhangFloatCarrierRational> nuisanceCoefficients;
    ZhangFloatCarrierRational physicalScale = 1;
    ZhangFloatCarrierRational symbolicOffset = 0;
    ZhangProductFloatTargetScope scope;
    std::map<KFKey, double> numericRow;
    double numericOffset = 0;
    bool active = true, numericBindingComplete = false;
    std::string bindingReason = "WHOLE_FLOAT_TARGET_NOT_BOUND";
};

// Exact binary value of an IEEE-754 input coefficient, without decimal or
// tolerance-based rank decisions. Nonfinite coefficients are never admitted.
inline std::optional<ZhangFloatCarrierRational> zhangFloatCarrierExactDouble(double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t) &&
        std::numeric_limits<double>::is_iec559);
    if (!std::isfinite(value)) return std::nullopt;
    const auto bits = std::bit_cast<std::uint64_t>(value);
    const auto exponentBits = (bits >> 52) & 0x7ff;
    const auto fraction = bits & ((std::uint64_t(1) << 52) - 1);
    if (exponentBits == 0 && fraction == 0) return ZhangFloatCarrierRational(0);
    ZhangExactInteger numeratorValue = exponentBits
        ? (std::uint64_t(1) << 52) | fraction : fraction;
    const int exponent = exponentBits ? int(exponentBits) - 1023 - 52 : -1074;
    ZhangFloatCarrierRational result;
    if (exponent >= 0) result = ZhangFloatCarrierRational(numeratorValue << exponent);
    else result = ZhangFloatCarrierRational(numeratorValue,
        ZhangExactInteger(1) << -exponent);
    if (bits >> 63) result = -result;
    return result;
}

// Mathematical planner input on one immutable, fully indexed KF snapshot or
// virtual drafted graph destination. It grants no graph/integer authority.
struct ZhangFloatCarrierIndexedTarget
{
    std::string targetId, unit;
    ZhangFloatCarrierSparseRow sourceRow;
    ZhangFloatCarrierRational offset = 0;
    bool complete = false;
};

struct ZhangFloatCarrierIndexedResolution
{
    std::string targetId, unit;
    ZhangFloatCarrierSparseRow destinationRow;
    ZhangFloatCarrierRational offset = 0;
};

struct ZhangFloatCarrierRowPlan
{
    bool valid = false, reducesDimension = false;
    std::string reason = "NOT_EVALUATED";
    std::size_t sourceDimension = 0, protectedCount = 0, carrierRank = 0;
    // Protected live/cache rows first, followed by independent original raw
    // target rows. No covariance or n*n matrix is needed to select them.
    std::vector<ZhangFloatCarrierSparseRow> transformRows;
    std::vector<std::string> carrierSourceTargetIds, carrierUnits;
    std::vector<ZhangFloatCarrierIndexedResolution> targets;
};

inline ZhangFloatCarrierRowPlan planZhangProductFloatCarrierRows(
    std::size_t sourceDimension,
    const std::vector<std::size_t>& protectedIndices,
    const std::vector<ZhangFloatCarrierIndexedTarget>& completeTargets)
{
    ZhangPhaseTimer timer("FLOAT_TARGET_CARRIER_EXACT_ROW_PLAN");
    ZhangFloatCarrierRowPlan out;
    out.sourceDimension = sourceDimension;
    out.protectedCount = protectedIndices.size();
    std::map<std::size_t, std::size_t> protectedPosition;
    for (const auto index : protectedIndices)
    {
        if (index >= sourceDimension || protectedPosition.contains(index))
        { out.reason = "FLOAT_CARRIER_PROTECTED_INDEX_INVALID"; return out; }
        protectedPosition.emplace(index, out.transformRows.size());
        out.transformRows.push_back({{index, 1}});
    }
    std::map<std::string, const ZhangFloatCarrierIndexedTarget*> ordered;
    for (const auto& target : completeTargets)
    {
        if (!target.complete || target.targetId.empty() || target.unit.empty())
        { out.reason = "FLOAT_CARRIER_ACTIVE_WHOLE_TARGET_UNBOUND"; return out; }
        if (!ordered.emplace(target.targetId, &target).second)
        { out.reason = "FLOAT_CARRIER_DUPLICATE_TARGET_ID"; return out; }
        for (const auto& [index, value] : target.sourceRow)
            if (index >= sourceDimension)
            { out.reason = "FLOAT_CARRIER_TARGET_SOURCE_INDEX_INVALID"; return out; }
    }
    using Rational = ZhangFloatCarrierRational;
    using Row = ZhangFloatCarrierSparseRow;
    struct ReducedBasisRow { Row row, originalCarrierWeights; };
    std::map<std::size_t, ReducedBasisRow> basis;
    auto add = [](Row& destination, const Row& source, const Rational& scale)
    {
        for (const auto& [index, coefficient] : source)
        {
            zhangExactWorkCheckpoint("FLOAT_CARRIER_SPARSE_REDUCTION");
            const Rational next = destination[index] + scale * coefficient;
            if (next == 0) destination.erase(index); else destination[index] = next;
        }
    };
    auto rawPart = [&](const Row& source)
    {
        Row result;
        for (const auto& [index, coefficient] : source)
            if (coefficient != 0 && !protectedPosition.contains(index))
                result.emplace(index, coefficient);
        return result;
    };
    auto reduce = [&](Row row)
    {
        Row weights;
        for (const auto& [pivot, source] : basis)
            if (const auto found = row.find(pivot); found != row.end())
            {
                const Rational factor = found->second;
                add(row, source.row, -factor);
                add(weights, source.originalCarrierWeights, factor);
            }
        return std::pair<Row, Row>{std::move(row), std::move(weights)};
    };
    for (const auto& [id, target] : ordered)
    {
        zhangExactWorkCheckpoint("FLOAT_CARRIER_SELECT_ORIGINAL_TARGET_ROW");
        Row original = rawPart(target->sourceRow);
        auto [remaining, previousWeights] = reduce(original);
        if (remaining.empty()) continue;
        const auto carrier = out.carrierRank++;
        out.transformRows.push_back(std::move(original));
        out.carrierSourceTargetIds.push_back(id);
        out.carrierUnits.push_back(target->unit);
        const auto pivot = remaining.begin()->first;
        const Rational divisor = remaining.begin()->second;
        for (auto& [column, coefficient] : remaining) coefficient /= divisor;
        Row normalizedWeights{{carrier, Rational(1) / divisor}};
        for (const auto& [index, coefficient] : previousWeights)
            normalizedWeights[index] = -coefficient / divisor;
        basis.emplace(pivot, ReducedBasisRow{
            std::move(remaining), std::move(normalizedWeights)});
    }
    for (const auto& [id, target] : ordered)
    {
        auto [remaining, weights] = reduce(rawPart(target->sourceRow));
        if (!remaining.empty())
        { out.reason = "FLOAT_CARRIER_COMPLETE_TARGET_NOT_RECONSTRUCTABLE"; return out; }
        ZhangFloatCarrierIndexedResolution resolved;
        resolved.targetId = id; resolved.unit = target->unit; resolved.offset = target->offset;
        for (const auto& [index, coefficient] : target->sourceRow)
            if (coefficient != 0 && protectedPosition.contains(index))
                resolved.destinationRow[protectedPosition.at(index)] = coefficient;
        for (const auto& [carrier, coefficient] : weights)
            if (coefficient != 0)
                resolved.destinationRow[out.protectedCount + carrier] = coefficient;
        Row pulledBack;
        for (const auto& [destination, coefficient] : resolved.destinationRow)
            add(pulledBack, out.transformRows.at(destination), coefficient);
        Row original;
        for (const auto& [index, coefficient] : target->sourceRow)
            if (coefficient != 0) original.emplace(index, coefficient);
        if (pulledBack != original)
        { out.reason = "FLOAT_CARRIER_EXACT_PULLBACK_AUDIT_FAILED"; return out; }
        out.targets.push_back(std::move(resolved));
    }
    out.valid = true;
    out.reducesDimension = out.transformRows.size() < sourceDimension;
    out.reason = "COMPLETE_FLOAT_TARGET_ROWS_PRESERVED_EXACTLY";
    return out;
}
