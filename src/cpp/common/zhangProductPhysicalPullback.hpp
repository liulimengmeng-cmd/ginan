#pragma once
#include "common/zhangProductRelationBasis.hpp"
#include "common/zhangProductPhysicalCycleChart.hpp"

/** Compile intended physical targets, never private chord labels, into the
 * authoritative posterior. The complete chart and the available posterior
 * chart are distinct: a subgraph change can be rectangular. */
struct ZhangProductPhysicalPullback
{
    bool valid = false;
    ZhangExactMatrix posteriorRows;
    std::vector<bool> available;
    std::vector<std::string> reasons;
    std::string failureReason;
};

inline ZhangProductPhysicalPullback zhangPullbackProductPhysicalRelations(
    ZhangProductRelationBasis& target,
    const ZhangGraphBasis& authoritative,
    const std::map<ZhangGraphEdge, int>& versions,
    const std::vector<ZhangGraphEdge>& authoritativeChords,
    const std::map<ZhangGraphEdge, int>& posteriorColumns,
    int posteriorSize)
{
    ZhangProductPhysicalPullback out;
    ZhangProductPhysicalCycleChart complete, posterior;
    complete.columns = authoritativeChords.size();
    posterior.columns = posteriorSize;
    const std::string signal = std::to_string(static_cast<int>(target.observable));
    for (int c = 0; c < complete.columns; ++c)
    {
        if (!complete.add(c, signal, authoritativeChords[c], authoritative, versions))
        {
            out.failureReason = "AUTHORITATIVE_PHYSICAL_CHART_INVALID";
            return out;
        }
        auto column = posteriorColumns.find(authoritativeChords[c]);
        if (column != posteriorColumns.end() &&
            !posterior.add(column->second, signal, authoritativeChords[c], authoritative, versions))
        {
            out.failureReason = "POSTERIOR_PHYSICAL_CHART_INVALID";
            return out;
        }
    }
    auto rebound = target.namedRelations;
    for (auto& relation : rebound)
    {
        std::map<std::string, ZhangExactInteger> intended;
        for (const auto& [edge, value] : relation.physicalArcCoefficients)
        {
            if (value == 0) continue;
            auto version = versions.find(edge);
            if (version == versions.end())
            {
                out.failureReason = "INTENDED_PHYSICAL_ARC_UNVERSIONED";
                return out;
            }
            intended[signal + "|" + edge.receiver + "|" + edge.satellite.id() +
                "|V" + std::to_string(version->second)] += value;
        }
        // A zero target is a dependent, deterministic catalogue entry.
        ZhangExactVector full(complete.columns), actual(posterior.columns);
        std::string reason = "EXACT_ZERO_PHYSICAL_TARGET";
        if (!intended.empty() && !complete.project(intended, full, &reason))
        {
            out.failureReason = "INTENDED_TARGET_NOT_IN_AUTHORITATIVE_CHART:" + reason;
            return out;
        }
        const bool available = intended.empty() || posterior.project(intended, actual, &reason);
        relation.currentCycleCoefficients = std::move(full);
        out.posteriorRows.push_back(std::move(actual));
        out.available.push_back(available);
        out.reasons.push_back(reason);
    }
    target.namedRelations = std::move(rebound);
    target.currentChords = authoritativeChords;
    out.valid = true;
    out.failureReason = "EXACT_PHYSICAL_ROUNDTRIP";
    return out;
}

inline bool zhangExactPosteriorRowToDouble(
    const ZhangExactVector& exact, Eigen::VectorXd& numeric)
{
    numeric.resize(exact.size());
    const ZhangExactInteger limit = ZhangExactInteger(1) << 53;
    for (int c = 0; c < numeric.size(); ++c)
    {
        if (exact[c] > limit || exact[c] < -limit) return false;
        numeric(c) = exact[c].convert_to<double>();
    }
    return true;
}
