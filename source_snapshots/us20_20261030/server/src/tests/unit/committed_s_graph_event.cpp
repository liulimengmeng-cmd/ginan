// Regression SOURCE ONLY. Not compiled or executed as part of this change.
#include "common/zhangCommittedSGraphEvent.hpp"
#include "common/zhangEpochComponentContext.hpp"

#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}
static ZhangGraphEdge edge(const char* receiver, int prn)
{
    return {receiver, SatSys(E_Sys::GPS, prn)};
}

int main()
{
    const std::set<ZhangGraphEdge> treeEdges{
        edge("A",1), edge("A",2), edge("B",1)};
    const auto tree = zhangBuildSpanningTree(treeEdges, "A");
    auto squareEdges = treeEdges;
    squareEdges.insert(edge("B",2));
    const auto square = zhangBuildSpanningTree(squareEdges, "A", tree.treeEdges);
    ZhangGraphBasis augmented;
    require(zhangAugmentCommittedSChords(tree, square, tree.edges, augmented),
        "same-tree fresh chord refused");
    require(augmented.edges == squareEdges && augmented.treeEdges == tree.treeEdges,
        "fresh chord not committed or S-tree changed");
    const auto insertion = zhangClassifyCommittedSGraphEvent(tree, augmented, true, false);
    require(insertion.kind == "CHORD_ONLY_INSERTION" && insertion.chordCatalogueChanged &&
        !insertion.treeReplacement && !insertion.physicalArcBreak,
        "chord insertion confused with tree/physical event");
    require(zhangClassifyCommittedSGraphEvent(augmented, augmented, true, false).kind ==
        "TOPOLOGY_UNCHANGED", "unchanged graph reclassified");
    require(zhangClassifyCommittedSGraphEvent(augmented, augmented, true, true).kind ==
        "PHYSICAL_ARC_BREAK", "same topology hid a real physical signal break");

    auto leafEdges = squareEdges;
    leafEdges.insert(edge("C",2));
    const auto leaf = zhangBuildSpanningTree(leafEdges,"A",tree.treeEdges);
    require(zhangClassifyCommittedSGraphEvent(square,leaf,true,false).nodeAugmentation,
        "new node not classified separately");
    require(!zhangAugmentCommittedSChords(square,leaf,square.edges,augmented),
        "new node silently treated as chord-only insertion");

    const std::set<ZhangGraphEdge> alternateTree{
        edge("A",1),edge("B",1),edge("B",2)};
    const auto replacement = zhangBuildSpanningTree(squareEdges,"A",alternateTree);
    require(zhangClassifyCommittedSGraphEvent(square,replacement,true,false).treeReplacement,
        "tree replacement not classified");
    require(!zhangAugmentCommittedSChords(square,replacement,square.edges,augmented),
        "tree replacement bypassed coordinate transport");
    const auto rootChange = zhangBuildSpanningTree(squareEdges,"B",square.treeEdges);
    require(!zhangAugmentCommittedSChords(square,rootChange,square.edges,augmented),
        "root-only replacement bypassed coordinate transport");

    auto forestEdges = squareEdges;
    forestEdges.insert(edge("C",3));
    forestEdges.insert(edge("C",4));
    forestEdges.insert(edge("D",3));
    forestEdges.insert(edge("D",4));
    const auto forest = zhangBuildSpanningTree(forestEdges,"A");
    auto joinedEdges = forestEdges;
    joinedEdges.insert(edge("C",2));
    const auto joined = zhangBuildSpanningTree(joinedEdges,"A",forest.treeEdges);
    require(zhangClassifyCommittedSGraphEvent(joined,forest,true,false).componentSplit,
        "component split missing");
    require(zhangClassifyCommittedSGraphEvent(forest,joined,true,false).componentReconnect,
        "component rejoin missing");

    // A legacy stale P tree must be completely irrelevant in native mode.
    const auto staleP = zhangBuildSpanningTree({edge("OLD",9)},"OLD");
    const auto components = zhangEpochComponentGraphs(forest,staleP,{},true);
    require(components.size()==2,"native component restriction merged/lost a component");
    for (const auto& view : components)
    {
        require(view.valid && zhangSameCommittedGraphBasis(view.basis,view.productBasis),
            "native product view differs from the only S chart");
        require(!view.productBasisRebuilt,"native view invoked a P optimizer");
        for (const auto& e : view.basis.treeEdges)
            require(forest.treeEdges.contains(e),"component restriction reselected S-tree");
        require(view.cycleRank==1,"component E-V+1 oracle mismatch");
    }
    std::cout << "COMMITTED_S_GRAPH_EVENT_PASS\n";
}
