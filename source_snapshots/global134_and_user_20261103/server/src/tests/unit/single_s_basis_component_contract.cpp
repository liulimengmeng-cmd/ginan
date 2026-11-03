#include "common/zhangEpochComponentContext.hpp"

#include <stdexcept>

// Source-only regressions. These assertions have NOT been compiled or run.
// The reference catalogue below is stated explicitly, not built by calling
// the production spanning-tree selector a second time.
int main()
{
    const SatSys g1(E_Sys::GPS,1),g2(E_Sys::GPS,2),g3(E_Sys::GPS,3);
    const ZhangGraphEdge a1{"AAAA",g1},a2{"AAAA",g2},a3{"AAAA",g3};
    const ZhangGraphEdge b1{"BBBB",g1},b2{"BBBB",g2},b3{"BBBB",g3};
    auto require=[](bool value,const char* reason) {
        if(!value)throw std::runtime_error(reason);
    };
    ZhangGraphBasis old;
    old.rootReceiver="AAAA";
    old.edges={a1,a2,a3,b1,b2};
    old.treeEdges={a1,a2,a3,b1};
    old.receivers={"AAAA","BBBB"};old.satellites={g1,g2,g3};
    old.componentCount=1;old.connected=true;
    auto augmented=old;augmented.edges.insert(b3);
    require(old.treeEdges==augmented.treeEdges,"Fixture changed its S tree");
    require(!zhangSameCommittedGraphBasis(old,augmented),
        "A new chord was mistaken for the unchanged committed catalogue");
    require(zhangSameCommittedGraphBasis(augmented,augmented),
        "Identical complete S catalogue did not match");

    // No current product observations are supplied. This cannot manufacture
    // product permission, but may not cause an independent P-tree selection.
    const auto views=zhangEpochComponentGraphs(augmented,old,{},true);
    require(views.size()==1 && views.front().valid,"Committed component lost");
    const auto& view=views.front();
    require(view.basis.edges==std::set<ZhangGraphEdge>{a1,a2,a3,b1,b2,b3},
        "Current S chord catalogue truncated");
    require(view.basis.treeEdges==std::set<ZhangGraphEdge>{a1,a2,a3,b1},
        "Read-only component changed the S tree");
    require(view.cycleRank==2,"Wrong exact graph cycle count");
    require(zhangSameCommittedGraphBasis(view.basis,view.productBasis) &&
        !view.productBasisRebuilt,"Native view created an independent product tree");

    auto changedRoot=augmented;changedRoot.rootReceiver="BBBB";
    require(!zhangSameCommittedGraphBasis(augmented,changedRoot),
        "Root change omitted from complete S identity");
    auto changedComponents=augmented;changedComponents.componentCount=2;
    require(!zhangSameCommittedGraphBasis(augmented,changedComponents),
        "Component identity omitted from complete S identity");
}
