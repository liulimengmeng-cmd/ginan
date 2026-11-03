#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangR51PhysicalImage.hpp"

#include <chrono>
#include <iostream>
#include <stdexcept>

static void require(bool value,const char* reason)
{
    if(!value)throw std::runtime_error(reason);
}

int main()
{
    try {
        ZhangExactMatrix dense(24,ZhangExactVector(24));
        for(int r=0;r<24;++r)for(int c=0;c<24;++c)dense[r][c]=(r+3)*(c+5)+(r==c?1:0);
        {
            // Cancellation happens INSIDE the real HNF, not at an outer
            // post-hoc time check. The caller sees no partial exact result.
            ZhangExactWorkBudget budget;budget.maximumCheckpoints=130;
            ZhangExactRowHnf result;result.basis={{999}};
            require(!zhangRunWithinExactWorkBudget(budget,[&] {
                return zhangExactRowHermiteNormalForm(dense,{},true);
            },result),"mid-HNF operation budget unexpectedly succeeded");
            require(budget.exhausted && budget.checkpoints==130,"HNF did not checkpoint its work");
            require(result.basis.empty(),"cancelled HNF leaked partial basis");
            require(zhangActiveExactWorkBudget==nullptr,"HNF budget scope leaked");
        }
        {
            ZhangExactWorkBudget budget;budget.maximumCheckpoints=130;
            ZhangExactMatrix output{{999}};
            require(!zhangRunWithinExactWorkBudget(budget,[&] {
                return zhangExactProjectedKernelColumns(dense,{{1,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0}},24);
            },output),"projected kernel ignored operation budget");
            require(budget.exhausted && output.empty(),"cancelled projected kernel leaked image");
        }
        {
            // Unlike the wrapper-entry test, this expired budget first gets
            // observed by an arithmetic-loop checkpoint in the backend.
            ZhangExactWorkBudget budget;budget.deadline=std::chrono::steady_clock::now()-std::chrono::seconds(1);
            bool cancelled=false;
            { ZhangScopedExactWorkBudget scope(budget);
              try { (void)zhangExactRowHermiteNormalForm(dense,{},true); }
              catch(const ZhangExactWorkBudgetExceeded&) {cancelled=true;} }
            require(cancelled && budget.exhausted && std::string(budget.reason)=="WALL_TIMEOUT",
                "expired deadline was not caught within exact elimination");
            require(budget.checkpoints>=64,"exact deadline was checked only outside preparation");
        }
        {
            ZhangExactWorkBudget budget;budget.maximumCheckpoints=90;
            ZhangR51PhysicalImage output;output.valid=true;
            ZhangExactMatrix target(2,ZhangExactVector(24));target[0][0]=1;target[1][1]=1;
            require(!zhangRunWithinExactWorkBudget(budget,[&] {
                return zhangR51PhysicalImage(target,dense,ZhangExactVector(24),target,24);
            },output),"physical image preparation ignored exact budget");
            require(!output.valid && output.projector.empty() && output.generators.empty(),
                "cancelled physical image exposed an apparently valid certificate");
        }
        {
            ZhangExactWorkBudget outer,inner;
            outer.deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
            inner.deadline=outer.deadline+std::chrono::seconds(5);
            { ZhangScopedExactWorkBudget scope(outer);
              { ZhangScopedExactWorkBudget nested(inner);
                require(inner.deadline==outer.deadline,"nested scope extended parent deadline"); }
              require(zhangActiveExactWorkBudget==&outer,"nested cancellation scope did not restore parent"); }
            require(zhangActiveExactWorkBudget==nullptr,"budget TLS did not restore normal mode");
            require(zhangR47AffineIntegerFeasible({{2,1}},{1},2),"normal exact operation changed after timeout");
        }
        {
            // A structural-zero companion supplies no new integer row and
            // cannot authorise the unavailable old L2 segment. Only the
            // searched nonzero L1 condition enters the current registry.
            constexpr long epoch=1405216800;
            const std::string l1="L1C|R|G13|V2",oldL2="L2W|R|G13|V1";
            ZhangProductPhysicalCycleChart chart;chart.columns=2;
            chart.chordIdentities={{0,l1},{1,oldL2}};
            chart.expansions[0][l1]=1;chart.expansions[1][oldL2]=1;
            chart.activePhysicalIdentities={l1,oldL2};
            ZhangCurrentEpochEvidence evidence;
            require(evidence.bind("runtime",E_Sys::GPS,epoch,"posterior","chart",chart,"FLOAT"),"fixture bind failed");
            auto proof=std::make_shared<ZhangIntegerDecisionProof>();proof->id="current-L1";
            proof->observationProvenance="runtime|PPP_COMMITTED|E1405216800";
            proof->originalStatement="nonzero current L1 ratio plus exact image";proof->conditionalFailureBound=1e-6;
            require(evidence.addEventOverlap({{l1,1}},7,{proof},"complete-current-target"),"nonzero row admission failed");
            require(!evidence.addEventOverlap({},0,{proof},"structural-only"),"structural zero became new AR evidence");
            ZhangProductIntegerTransition target;target.valid=true;target.observable=E_ObsCode::L1C;
            target.physicalEdges={{"R",SatSys(E_Sys::GPS,13)}};
            target.physicalArcVersions={2};target.coefficients={1};
            const auto known=evidence.query(E_ObsCode::L1C,epoch,"posterior","chart",target,false);
            require(known.proven && known.value==7,"nonzero searched signal lost current proof");
            target.observable=E_ObsCode::L2W;target.physicalArcVersions={1};
            require(!evidence.query(E_ObsCode::L2W,epoch,"posterior","chart",target,false).proven,
                "identity companion authorised an unavailable old L2 integer");
            require(evidence.eventOverlapRows.size()==1,"zero companion added an integer fixed row");
        }
        {
            ZhangProductPhysicalCycleChart chart;chart.columns=2;
            chart.expansions[0]["L1C|R|G13|V2"]=1;
            std::map<std::string,ZhangExactInteger> output{{"old-output",99}};
            require(!chart.expand({1,1},output) && output.empty(),
                "unknown cycle column returned false with partial physical output");
            require(chart.expand({1,0},output) && output.size()==1 &&
                output.at("L1C|R|G13|V2")==1,"known physical expansion changed");
            require(!chart.expand({0,1},output) && output.empty(),
                "unknown physical column silently became zero");
            output["old-output"]=99;
            require(!chart.expand({1},output) && output.empty(),
                "wrong dimension left stale physical expansion output");
        }
        {
            // Independent, explicitly written target oracle. Column 1 is an
            // unknown physical variable, not zero; rows using it are merely
            // unavailable *proposals*. The untouched three-column numeric
            // domain must stay paired with physical rows 0 and 2 in order.
            ZhangProductPhysicalCycleChart chart;chart.columns=3;
            chart.expansions[0]={{"N0",1},{"N1",1}};
            chart.expansions[2]={{"N0",1},{"N2",1}};
            chart.chordIdentities={{0,"N1"},{2,"N2"}};
            chart.activePhysicalIdentities={"N0","N1","N2"};
            chart.unexpandedColumns[1]="UNKNOWN_SIGNAL_ARC_VERSION";
            const ZhangExactMatrix targets{{2,0,0},{0,1,0},{0,0,1},{1,1,0}};
            const auto subset=zhangPhysicalTargetProposalSubset(chart,targets,true);
            require(subset.valid && subset.retainedIndices==std::vector<int>({0,2}) &&
                subset.rejectedIndices==std::vector<int>({1,3}),"proposal subset changed original row order");
            require(subset.numericTargets==ZhangExactMatrix({{2,0,0},{0,0,1}}),
                "proposal filtering compressed or permuted the numeric columns");
            const std::vector<std::map<std::string,ZhangExactInteger>> expectedPhysical{
                {{"N0",2},{"N1",2}},{{"N0",1},{"N2",1}}};
            require(subset.physicalTargets==expectedPhysical,"physical and numeric target row order diverged");
            require(zhangProductPhysicalRowsRepresentable(chart,subset.numericTargets),
                "unused unknown column poisoned exact product candidate rows");
            require(!zhangProductPhysicalRowsRepresentable(chart,{{0,1,0}},true),
                "unknown nonzero product candidate or conditioner was authorized");
            require(!zhangProductPhysicalRowsRepresentable(chart,{{0,0,0}}) &&
                zhangProductPhysicalRowsRepresentable(chart,{{0,0,0}},true),
                "structural coordinate zero was confused with new integer evidence");
            // H*N=3 fixes N0 only. Hence targets are (6+2*N1,3+N2),
            // image diag(2,1), and q=(a0-3,a2-3) in the ORIGINAL x0,x1,x2
            // numerical chart. These expectations are algebraic, not a
            // second call to the same projection implementation.
            const auto image=zhangR51PhysicalImage({{2,2,0},{1,0,1}},{{1,0,0}},{3},
                subset.numericTargets,3);
            require(image.valid && image.generators==ZhangExactMatrix({{2,0},{0,1}}) &&
                image.particularTarget==ZhangExactVector({6,3}),"filtered physical integer image changed");
            require(image.projector==std::vector<ZhangR51RationalRow>({{1,0,0},{0,0,1}}) &&
                image.offsets==ZhangR51RationalRow({-3,-3}),"filtered numeric projector/offset differs from independent oracle");
            const auto strict=zhangPhysicalTargetProposalSubset(chart,targets,false);
            require(!strict.valid && strict.numericTargets.empty() && strict.physicalTargets.empty(),
                "legacy/mandatory target expansion silently omitted a row");
            std::map<std::string,ZhangExactInteger> history;
            require(!chart.expand({0,1,0},history) && history.empty(),
                "unknown conditioning variable was silently discarded");
            require(!zhangPhysicalTargetProposalSubset(chart,{{1,0}},true).valid,
                "shape error became a filtered current proposal");
            chart.failureReason="CURRENT_PHYSICAL_GRAPH_CONTEXT_MISSING";
            require(!zhangPhysicalTargetProposalSubset(chart,{{1,0,0}},true).valid,
                "global chart failure became a subset search");
        }
        {
            const SatSys a(E_Sys::GPS,1),b(E_Sys::GPS,2);
            const std::set<ZhangGraphEdge> edges{{"R1",a},{"R1",b},{"R2",a},{"R2",b}};
            const auto basis=zhangBuildSpanningTree(edges,"R1");
            ZhangGraphEdge chord;for(const auto& edge:edges)
                if(!basis.treeEdges.contains(edge))chord=edge;
            std::map<ZhangGraphEdge,int> versions;for(const auto& edge:edges)versions[edge]=4;
            const auto cycle=zhangFundamentalCycle(basis,chord);
            require(cycle.size()==4,"physical chart regression fixture has no full cycle");
            versions.erase(cycle.begin()->first);
            ZhangProductPhysicalCycleChart partial;partial.columns=1;
            require(partial.add(0,"L1C",chord,basis,versions,true) &&
                partial.expansions.empty() && partial.unexpandedColumns.contains(0),
                "unknown signal version invalidated the entire repair chart or became V0");
            std::map<std::string,ZhangExactInteger> physical;
            require(!partial.expand({1},physical) && physical.empty(),
                "unknown nonzero column was accepted as a physical zero");
            ZhangProductPhysicalCycleChart legacy;legacy.columns=1;
            require(!legacy.add(0,"L1C",chord,basis,versions),"legacy chart accepted an unknown signal version");
            versions.clear();for(const auto& edge:edges)versions[edge]=4;
            require(partial.add(0,"L1C",chord,basis,versions,true) &&
                !partial.unexpandedColumns.contains(0) && partial.expand({1},physical),
                "known signal version did not recover the original numeric coordinate");
            const auto knownPhysical=physical;ZhangExactVector projected;
            require(partial.project(knownPhysical,projected) && projected==ZhangExactVector({1}),
                "known physical cycle does not roundtrip before catalogue retirement");
            versions.erase(cycle.begin()->first);
            require(partial.add(0,"L1C",chord,basis,versions,true) &&
                partial.expansions.empty() && partial.chordIdentities.empty() &&
                !partial.expand({1},physical) && physical.empty(),
                "known-to-unknown column retained an old physical expansion");
            require(!partial.project(knownPhysical,projected),
                "stale active identity directory authorised an unknown column");
            versions.clear();for(const auto& edge:edges)versions[edge]=4;
            require(partial.add(0,"L1C",chord,basis,versions,true) && partial.expand({1},physical),
                "known-to-unknown-to-known column did not safely recover");
        }
        std::cout<<"Exact preparation deadlines, scoped cancellation, nonpartial results and structural companion tests passed\n";
        return 0;
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
