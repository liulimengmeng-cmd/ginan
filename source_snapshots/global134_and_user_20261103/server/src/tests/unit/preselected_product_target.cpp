#include "common/zhangPreselectedProductTarget.hpp"
#include <iostream>
#include <stdexcept>

static int checks = 0;
static void require(bool value, const char* why)
{
    ++checks;
    if (!value) throw std::runtime_error(why);
}
static ZhangGraphEdge edge(const char* receiver, int satellite)
{
    return {receiver, SatSys(E_Sys::GPS, satellite)};
}
int main()
{
    const std::set<ZhangGraphEdge> edges{edge("A",1),edge("A",2),edge("B",1),edge("B",2)};
    const auto committed = zhangBuildSpanningTree(edges,"A");
    const auto originalTree = committed.treeEdges;
    auto product = zhangBuildRootedProductTree(edges,"A");
    auto represented = product.treeEdges;
    std::map<ZhangGraphEdge,int> selected, current;
    std::set<std::pair<ZhangGraphEdge,int>> histories;
    std::map<ZhangGraphEdge,std::map<int,int>> selectedSignals, currentSignals;
    for (const auto& e : product.treeEdges)
    {
        selected[e] = 2; current[e] = 2; histories.insert({e,2});
        selectedSignals[e] = {{1,3},{2,4}};
        currentSignals[e] = {{1,3},{2,4}};
    }
    auto audit = [&]() { return zhangAuditPreselectedProductTarget(
        committed,product,represented,selected,current,histories,
        selectedSignals,currentSignals,{1,2}); };
    require(audit().valid,"well formed current selector rejected");
    const auto first = *product.treeEdges.begin();
    histories.erase({first,2});
    require(!audit().valid && audit().reason == "PRODUCT_TOPOLOGY_ARC_OBSERVATION_HISTORY_MISSING" &&
        audit().hasEdge && audit().edge == first, "missing exact version history not diagnosed");
    histories.insert({first,2});
    selectedSignals[first].erase(2);
    require(!audit().valid && audit().reason == "PRODUCT_BASELINE_SIGNAL_ARC_INCOMPLETE",
        "partial baseline identity accepted");
    selectedSignals[first][2] = 5;
    require(!audit().valid && audit().reason == "PRODUCT_SIGNAL_ARC_VERSION_INVALID",
        "future physical signal version accepted");
    selectedSignals[first][2] = 4;
    selected[first] = 3;
    require(!audit().valid && audit().reason == "PRODUCT_TOPOLOGY_ARC_VERSION_INVALID",
        "future legacy version accepted");
    selected[first] = 2;
    represented.erase(first);
    require(!audit().valid && audit().reason == "PRODUCT_TREE_EDGE_NOT_REPRESENTED_BY_COMMITTED_S",
        "posterior-missing edge accepted from topology alone");
    represented.insert(first);
    require(audit().valid,"restored valid selector failed");

    // Freeze the 05:07 causal shape: preselected P adds TGRI/new satellite,
    // while the final S graph did not represent that receiver/arc. No history
    // or physical state is manufactured to make the selector exportable.
    auto provisionalEdges = edges;
    provisionalEdges.insert(edge("TGRI",2)); provisionalEdges.insert(edge("TGRI",18));
    product = zhangBuildRootedProductTree(provisionalEdges,"A");
    const auto invalid = audit();
    require(!invalid.valid && invalid.hasEdge && invalid.edge.receiver == "TGRI" &&
        invalid.reason == "PRODUCT_TREE_EDGE_NOT_REPRESENTED_BY_COMMITTED_S",
        "pre-factor target-only receiver survived finalized S audit");
    require(committed.edges == edges && committed.treeEdges == originalTree &&
        !histories.contains({edge("TGRI",18),2}), "audit mutated source chart/history");
    std::cout << "PRESELECTED_PRODUCT_TARGET_PASS checks=" << checks
              << " unsupported_nodes_revoked=1 missing_history_rejected=1 source_unchanged=1\n";
}
