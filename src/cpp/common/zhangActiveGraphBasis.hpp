#pragma once
#include "common/zhangFullRank.hpp"

struct ZhangActiveGraphBasisResult
{
    bool valid = false;
    ZhangGraphBasis basis;
    std::set<ZhangGraphEdge> retainedDatumOnlyEdges;
    std::string failureReason;
};

struct ZhangEligibilityPartition
{
    std::set<ZhangGraphEdge> withinRepresentedNodes;
    std::set<ZhangGraphEdge> requiresAugmentation;
};

inline ZhangEligibilityPartition partitionZhangEligibility(
    const ZhangGraphBasis& represented,
    const std::set<ZhangGraphEdge>& eligible)
{
    ZhangEligibilityPartition result;
    for (const auto& edge : eligible)
    {
        auto& destination = represented.receivers.contains(edge.receiver) &&
                            represented.satellites.contains(edge.satellite)
            ? result.withinRepresentedNodes : result.requiresAugmentation;
        destination.insert(edge);
    }
    return result;
}

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

// First validate the transported old chart. Only then extend its tree to
// new nodes. Preferred tree edges cannot be replaced by this augmentation;
// no new arc is used to reconstruct an old coordinate or certify an integer.
inline ZhangActiveGraphBasisResult zhangActiveGraphBasisWithAugmentation(
    const ZhangGraphBasis& represented,
    const std::set<ZhangGraphEdge>& eligible)
{
    const auto parts = partitionZhangEligibility(represented, eligible);
    auto result = zhangActiveGraphBasisAfterTransform(
        represented, parts.withinRepresentedNodes);
    if (!result.valid || parts.requiresAugmentation.empty()) return result;
    auto extendedEdges = result.basis.edges;
    extendedEdges.insert(parts.requiresAugmentation.begin(), parts.requiresAugmentation.end());
    auto extended = zhangBuildSpanningTree(
        extendedEdges, represented.rootReceiver, represented.treeEdges);
    if (!extended.connected || !std::includes(
            extended.treeEdges.begin(), extended.treeEdges.end(),
            represented.treeEdges.begin(), represented.treeEdges.end()))
    {
        result.valid = false;
        result.failureReason = "NEW_NODE_AUGMENTATION_NOT_A_TREE_EXTENSION";
        return result;
    }
    result.basis = std::move(extended);
    return result;
}

inline ZhangGraphBasis zhangRepresentedBasisAfterAugmentation(
    const ZhangGraphBasis& transported, const ZhangGraphBasis& active)
{
    auto represented = active;
    represented.edges.insert(transported.edges.begin(), transported.edges.end());
    return represented;
}
