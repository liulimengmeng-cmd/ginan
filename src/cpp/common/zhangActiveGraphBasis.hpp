#pragma once
#include "common/zhangFullRank.hpp"

struct ZhangActiveGraphBasisResult
{
    bool valid = false;
    ZhangGraphBasis basis;
    std::set<ZhangGraphEdge> retainedDatumOnlyEdges;
    std::string failureReason;
};

// A transformed tree defines coordinates, even when one of its edges has
// no current observation/state eligibility. Preserve that chart without
// promoting datum-only edges to observationEdges or stateEdges. Inactive
// chords are not retained. Never silently choose a different integer basis.
inline ZhangActiveGraphBasisResult zhangActiveGraphBasisAfterTransform(
    const ZhangGraphBasis& represented,
    const std::set<ZhangGraphEdge>& eligibleEdges)
{
    ZhangActiveGraphBasisResult result;
    for (const auto& edge : represented.treeEdges)
    {
        if (!represented.edges.contains(edge))
        {
            result.failureReason = "REPRESENTED_TREE_EDGE_MISSING";
            return result;
        }
        if (!eligibleEdges.contains(edge))
            result.retainedDatumOnlyEdges.insert(edge);
    }
    auto edges = eligibleEdges;
    edges.insert(represented.treeEdges.begin(), represented.treeEdges.end());
    result.basis = zhangBuildSpanningTree(
        edges, represented.rootReceiver, represented.treeEdges);
    if (!result.basis.connected ||
        result.basis.treeEdges != represented.treeEdges)
    {
        result.failureReason = "ACTIVE_SUPPORT_REQUIRES_NEW_COORDINATE_TRANSFORM";
        return result;
    }
    result.valid = true;
    return result;
}
