#pragma once

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "common/zhangFullRank.hpp"

/** Read-only integer views of a committed FLOAT forest.  The existing S-tree
 * is restricted, never reselected: an integer query must not reparameterise
 * the KF merely to obtain a connected input. In explicit legacy mode product
 * trees may select targets; native mode never builds or reselects one. */
struct ZhangEpochComponentGraph
{
    ZhangGraphBasis basis;
    ZhangGraphBasis productBasis;
    std::string componentId;
    std::string failureReason = "NOT_EVALUATED";
    int sourceComponentCount = 0;
    int cycleRank = 0;
    bool containsSourceRoot = false;
    bool productBasisRebuilt = false;
    bool valid = false;
};

inline std::vector<ZhangEpochComponentGraph> zhangEpochComponentGraphs(
    const ZhangGraphBasis& source,
    const ZhangGraphBasis& previousProduct,
    const std::set<ZhangGraphEdge>& currentProductEdges,
    bool singleSBasisNative = false)
{
    using namespace zhang_graph_detail;
    DisjointSet components;
    for (const auto& edge : source.edges)
        components.unite(receiverNode(edge.receiver), satelliteNode(edge.satellite));
    std::map<std::string, std::set<ZhangGraphEdge>> edgeGroups;
    for (const auto& edge : source.edges)
        edgeGroups[components.find(receiverNode(edge.receiver))].insert(edge);
    std::vector<ZhangEpochComponentGraph> result;
    for (const auto& [identity, edges] : edgeGroups)
    {
        ZhangEpochComponentGraph out;
        out.sourceComponentCount = static_cast<int>(edgeGroups.size());
        std::set<ZhangGraphEdge> retainedTree;
        std::set_intersection(edges.begin(), edges.end(),
            source.treeEdges.begin(), source.treeEdges.end(),
            std::inserter(retainedTree, retainedTree.end()));
        std::set<std::string> receivers;
        for (const auto& edge : edges) receivers.insert(edge.receiver);
        out.containsSourceRoot = receivers.contains(source.rootReceiver);
        const auto root = out.containsSourceRoot
            ? source.rootReceiver : *receivers.begin();
        out.basis = zhangBuildSpanningTree(edges, root, retainedTree);
        out.cycleRank = static_cast<int>(edges.size()) -
            static_cast<int>(out.basis.receivers.size() + out.basis.satellites.size()) + 1;
        std::ostringstream componentId;
        componentId << root;
        // Identity includes every physical endpoint, not only a root name.
        for (const auto& edge : edges)
            componentId << '|' << edge.receiver << ':'
                << static_cast<int>(edge.satellite.sys) << ':' << edge.satellite.prn;
        out.componentId = componentId.str();
        if (!out.basis.connected || out.basis.treeEdges != retainedTree || out.cycleRank < 0)
        {
            out.failureReason = "COMMITTED_COMPONENT_S_TREE_INVALID";
            result.push_back(std::move(out));
            continue;
        }
        if (singleSBasisNative)
        {
            // Read-only restriction of ONE committed S forest. Current signal
            // support is checked separately by the owner and writer: missing
            // observations never authorize constructing a replacement datum.
            out.productBasis = out.basis;
            out.productBasisRebuilt = false;
            out.failureReason = "NONE";
            out.valid = true;
            result.push_back(std::move(out));
            continue;
        }
        std::set<ZhangGraphEdge> productEdges;
        std::set_intersection(edges.begin(), edges.end(),
            currentProductEdges.begin(), currentProductEdges.end(),
            std::inserter(productEdges, productEdges.end()));
        const bool oldTargetMatches = previousProduct.connected &&
            previousProduct.rootReceiver == root &&
            previousProduct.satellites == out.basis.satellites &&
            std::includes(out.basis.receivers.begin(), out.basis.receivers.end(),
                previousProduct.receivers.begin(), previousProduct.receivers.end()) &&
            std::includes(productEdges.begin(), productEdges.end(),
                previousProduct.treeEdges.begin(), previousProduct.treeEdges.end());
        if (oldTargetMatches)
        {
            out.productBasis = previousProduct;
            // Non-tree catalogue edges are not allowed to retain stale
            // physical support just because the old tree itself survived.
            out.productBasis.edges.clear();
            std::set_intersection(previousProduct.edges.begin(), previousProduct.edges.end(),
                productEdges.begin(), productEdges.end(),
                std::inserter(out.productBasis.edges, out.productBasis.edges.end()));
        }
        else
        {
            out.productBasis = zhangBuildRootedProductTree(
                productEdges, root, previousProduct.treeEdges, {}, productEdges);
            out.productBasisRebuilt = true;
        }
        if (!out.productBasis.connected ||
            out.productBasis.rootReceiver != root ||
            out.productBasis.satellites != out.basis.satellites ||
            !std::includes(out.basis.receivers.begin(), out.basis.receivers.end(),
                out.productBasis.receivers.begin(), out.productBasis.receivers.end()))
        {
            out.failureReason = "CURRENT_PRODUCT_GRAPH_NOT_SPANNING_COMPONENT";
            // The integer component itself is usable.  A missing product
            // selector remains explicit and must not borrow a stale target.
            out.productBasis = {};
        }
        else out.failureReason = "NONE";
        out.valid = true;
        result.push_back(std::move(out));
    }
    std::stable_sort(result.begin(), result.end(), [](const auto& a, const auto& b)
    {
        const auto score = [](const auto& x)
        {
            return std::make_tuple(x.valid && x.cycleRank > 0,
                x.containsSourceRoot, x.productBasis.connected,
                x.cycleRank, x.basis.satellites.size(), x.basis.edges.size());
        };
        if (score(a) != score(b)) return score(a) > score(b);
        return a.componentId < b.componentId;
    });
    return result;
}
