#include "common/zhangEpochComponentContext.hpp"

#include <iostream>
#include <stdexcept>

static int checks = 0;
static void require(bool condition, const char* reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}
static ZhangGraphEdge edge(const char* receiver, int satellite)
{
    return {receiver, SatSys(E_Sys::GPS, satellite)};
}

int main()
{
    // Two independent bipartite squares. An independent graph-count oracle
    // gives one cycle in each, not a cross-component integer relationship.
    const std::set<ZhangGraphEdge> left{
        edge("A",1), edge("A",2), edge("B",1), edge("B",2)};
    const std::set<ZhangGraphEdge> right{
        edge("C",3), edge("C",4), edge("D",3), edge("D",4)};
    auto both = left;
    both.insert(right.begin(), right.end());
    const auto forest = zhangBuildSpanningTree(both, "A");
    const auto beforeEdges = forest.edges;
    const auto beforeTree = forest.treeEdges;
    const auto views = zhangEpochComponentGraphs(forest, {}, both);
    require(!forest.connected && forest.componentCount == 2, "fixture is not a forest");
    require(views.size() == 2, "independent components were discarded or merged");
    require(views[0].containsSourceRoot, "stable source root component not preferred");
    std::set<ZhangGraphEdge> reconstructed;
    for (const auto& view : views)
    {
        require(view.valid && view.basis.connected, "real connected component rejected");
        require(view.cycleRank == 1, "component cycle rank differs from independent E-V+1");
        require(view.basis.treeEdges.size() == 3, "S-tree was not restricted exactly");
        require(view.productBasis.connected, "current physical product target was not rebuilt");
        require(view.productBasis.rootReceiver == view.basis.rootReceiver &&
            view.productBasis.satellites == view.basis.satellites,
            "product target does not match its component snapshot");
        for (const auto& e : view.basis.treeEdges)
            require(beforeTree.contains(e), "view reselected authoritative S coordinates");
        reconstructed.insert(view.basis.edges.begin(), view.basis.edges.end());
    }
    require(reconstructed == beforeEdges && forest.treeEdges == beforeTree,
        "component extraction mutated/dropped authoritative FLOAT structure");
    require(views[0].componentId != views[1].componentId, "component identities alias");

    // An old target can have another root, extra satellite and dead receiver.
    // It must be rebuilt, not accepted by skipping the compatibility guard.
    auto staleEdges = left;
    staleEdges.insert(edge("OLD",5));
    const auto stale = zhangBuildSpanningTree(staleEdges, "OLD");
    const auto fresh = zhangEpochComponentGraphs(
        zhangBuildSpanningTree(left,"A"), stale, left);
    require(fresh.size() == 1 && fresh[0].productBasisRebuilt &&
        fresh[0].productBasis.connected, "node/root changes did not rebuild current target");
    require(fresh[0].productBasis.rootReceiver == "A" &&
        !fresh[0].productBasis.satellites.contains(SatSys(E_Sys::GPS,5)),
        "stale node or root leaked into the new product target");

    // Current observations which do not span the committed component may
    // not be silently completed with stale/grace edges.
    const auto unsupported = zhangEpochComponentGraphs(
        zhangBuildSpanningTree(left,"A"), fresh[0].productBasis, {edge("A",1)});
    require(unsupported[0].valid && !unsupported[0].productBasis.connected &&
        unsupported[0].failureReason == "CURRENT_PRODUCT_GRAPH_NOT_SPANNING_COMPONENT",
        "stale product support manufactured current connectivity");

    // A tree component is structurally valid but has no ILS cycle to search.
    auto treeEdges = left;
    treeEdges.erase(edge("B",2));
    const auto tree = zhangBuildSpanningTree(treeEdges,"A");
    const auto zero = zhangEpochComponentGraphs(tree, {}, treeEdges);
    require(zero.size() == 1 && zero[0].valid && zero[0].cycleRank == 0,
        "zero-cycle graph was confused with no observations");

    // Reconnection creates exactly one connected component; no previous
    // independent components or selected roots are inherited as certificates.
    auto joined = both;
    joined.insert(edge("C",2));
    const auto connected = zhangBuildSpanningTree(joined,"A");
    const auto merged = zhangEpochComponentGraphs(connected, views[1].productBasis, joined);
    require(merged.size() == 1 && merged[0].cycleRank == 2 &&
        merged[0].productBasisRebuilt, "component rejoin reused an incompatible product chart");

    // A malformed tree cannot be repaired by silently choosing a new S-tree.
    auto malformed = connected;
    malformed.treeEdges.clear();
    const auto rejected = zhangEpochComponentGraphs(malformed, {}, joined);
    require(rejected.size() == 1 && !rejected[0].valid &&
        rejected[0].failureReason == "COMMITTED_COMPONENT_S_TREE_INVALID",
        "malformed committed coordinates silently reparameterised");
    std::cout << "EPOCH_COMPONENT_CONTEXT_PASS checks=" << checks
              << " independent_components=2 source_S_unchanged=1 current_targets_only=1\n";
}
