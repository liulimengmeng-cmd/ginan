#include <boost/test/unit_test.hpp>
#include "common/zhangActiveGraphBasis.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"

namespace {
KFState coordinateFixture(const ZhangGraphBasis& basis)
{
    KFState state;
    state.kfIndexMap.clear();
    auto add = [&](const KFKey& key) {
        if (!state.kfIndexMap.contains(key))
            state.kfIndexMap[key] = state.kfIndexMap.size();
        state.stateTransitionMap[key][key][0] = 1;
    };
    KFKey external; external.type=KF::IONO_STEC; external.str="external"; add(external);
    for (auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    {
        for (const auto& receiver : basis.receivers)
            if(receiver!=basis.rootReceiver) add(zhangTransportReceiverKey(E_Sys::GPS,code,receiver));
        for (const auto& satellite : basis.satellites) add(zhangTransportSatelliteKey(code,satellite));
        for (const auto& edge : basis.edges)
            if(!basis.treeEdges.contains(edge)) add(zhangTransportAmbiguityKey(code,edge));
    }
    int n=state.kfIndexMap.size();
    state.x=Eigen::VectorXd::LinSpaced(n,-3,4); state.dx=Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd lower=Eigen::MatrixXd::Identity(n,n);
    for(int i=0;i<n;i++) for(int j=0;j<i;j++) lower(i,j)=0.07*(1+(i+j)%5);
    state.P=lower*lower.transpose();
    return state;
}
ZhangGraphCoordinatePlan planFixture(const KFState& state,const ZhangGraphBasis& old,
    const ZhangGraphBasis& next,const std::set<ZhangGraphEdge>& retired={})
{
    return zhangPlanGraphCoordinateTransport(state,E_Sys::GPS,old,next,retired,
        {{E_ObsCode::L1C,0.1902936728},{E_ObsCode::L2W,0.2442102134}},
        [](const auto&,auto){return 100.;},[](const auto&,auto){return 25.;});
}
KFState evaluatePlan(const KFState& old,const ZhangGraphCoordinatePlan& plan,
    const ZhangGraphBasis& basis)
{
    BOOST_REQUIRE_MESSAGE(plan.valid,plan.failureReason);
    KFState next; next.kfIndexMap.clear(); next.stateTransitionMap.clear();
    const int n=plan.transform.size();
    Eigen::MatrixXd t=Eigen::MatrixXd::Zero(n,old.x.size());
    Eigen::MatrixXd g=Eigen::MatrixXd::Zero(n,plan.independentSourceVariances.size());
    std::map<KFKey,int> independent; for(const auto& [key,v]:plan.independentSourceVariances)
        independent[key]=independent.size();
    for(const auto& [key,row]:plan.transform) {
        int i=next.kfIndexMap.size(); next.kfIndexMap[key]=i; next.stateTransitionMap[key][key][0]=1;
        for(const auto& [source,c]:row) {
            if(old.kfIndexMap.contains(source)) t(i,old.kfIndexMap.at(source))=c;
            else g(i,independent.at(source))=c*std::sqrt(plan.independentSourceVariances.at(source));
        }
    }
    next.x=t*old.x; next.dx=t*old.dx; next.P=t*old.P*t.transpose()+g*g.transpose();
    int count=plan.survivingPhysicalRows.size();
    Eigen::MatrixXd oldH=Eigen::MatrixXd::Zero(count,old.x.size());
    Eigen::MatrixXd newH=Eigen::MatrixXd::Zero(count,n);
    int i=0;
    for(const auto& [id,row]:plan.survivingPhysicalRows) {
        const auto& [code,edge]=id;
        for(const auto& [key,c]:row) oldH(i,old.kfIndexMap.at(key))=c;
        if(edge.receiver!=basis.rootReceiver)
            newH(i,next.kfIndexMap.at(zhangTransportReceiverKey(E_Sys::GPS,code,edge.receiver)))=1;
        newH(i,next.kfIndexMap.at(zhangTransportSatelliteKey(code,edge.satellite)))=1;
        if(!basis.treeEdges.contains(edge)) newH(i,next.kfIndexMap.at(zhangTransportAmbiguityKey(code,edge)))=
            code==E_ObsCode::L1C?0.1902936728:0.2442102134;
        i++;
    }
    BOOST_CHECK_SMALL((newH*t-oldH).norm(),1e-10);
    BOOST_CHECK_SMALL((newH*g).norm(),1e-10);
    BOOST_CHECK_SMALL((newH*next.x-oldH*old.x).norm(),1e-10);
    BOOST_CHECK_SMALL((newH*next.P*newH.transpose()-oldH*old.P*oldH.transpose()).norm(),1e-9);
    // Full cross covariance with the independent non-phase state is conserved.
    KFKey external; external.type=KF::IONO_STEC; external.str="external";
    BOOST_CHECK_SMALL((newH*next.P.col(next.kfIndexMap.at(external))-
        oldH*old.P.col(old.kfIndexMap.at(external))).norm(),1e-10);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(next.P);
    BOOST_CHECK(eig.eigenvalues().minCoeff()>-1e-10);
    return next;
}
}

BOOST_AUTO_TEST_CASE(r51_component_transport_preserves_physical_rows_and_full_covariance)
{
    const ZhangGraphEdge a1{"A",SatSys("G01")},a2{"A",SatSys("G02")},
        b1{"B",SatSys("G01")},b2{"B",SatSys("G02")},
        c2{"C",SatSys("G02")},c3{"C",SatSys("G03")},
        d2{"D",SatSys("G02")},d3{"D",SatSys("G03")};
    auto old=zhangBuildSpanningTree({a1,a2,b1,b2,c2,c3,d2,d3},"A");
    auto state=coordinateFixture(old);
    auto pivot=zhangBuildSpanningTree(old.edges,"A",{a1,b1,b2,c2,c3,d3});
    auto plan=planFixture(state,old,pivot);
    BOOST_CHECK(plan.independentSourceVariances.empty());
    evaluatePlan(state,plan,pivot);
    // Physical removal of both edges to G02 splits A/B/G01 from C/D/G02/G03.
    auto split=zhangBuildSpanningTree({a1,b1,c2,c3,d2,d3},"A");
    BOOST_REQUIRE(!split.connected);
    plan=planFixture(state,old,split,{a2,b2});
    auto splitState=evaluatePlan(state,plan,split);
    BOOST_CHECK_EQUAL(plan.survivingPhysicalRows.size(),12);
    BOOST_CHECK(plan.independentSourceVariances.empty());
    auto rootKey=zhangTransportReceiverKey(E_Sys::GPS,E_ObsCode::L1C,"C");
    BOOST_CHECK(splitState.P(splitState.kfIndexMap.at(rootKey),splitState.kfIndexMap.at(rootKey))>0);
    BOOST_CHECK_EQUAL(split.edges.size()-split.treeEdges.size(),1); // internal cycle survives
    auto rejoin=zhangBuildSpanningTree(old.edges,"A");
    plan=planFixture(splitState,split,rejoin);
    BOOST_CHECK(plan.freshEdges==std::set<ZhangGraphEdge>({a2,b2}));
    BOOST_CHECK_EQUAL(plan.independentSourceVariances.size(),4);
    evaluatePlan(splitState,plan,rejoin);
}

BOOST_AUTO_TEST_CASE(r51_component_transport_new_node_retirement_and_dynamic_guard)
{
    const ZhangGraphEdge a1{"A",SatSys("G01")},a2{"A",SatSys("G02")},
        b1{"B",SatSys("G01")},b2{"B",SatSys("G02")},c2{"NEW",SatSys("G02")};
    auto old=zhangBuildSpanningTree({a1,a2,b1,b2},"A"); auto state=coordinateFixture(old);
    auto next=zhangBuildSpanningTree({a1,a2,b1,b2,c2},"A",{a1,b1,b2,c2});
    auto plan=planFixture(state,old,next,{a2});
    BOOST_CHECK_EQUAL(plan.survivingPhysicalRows.size(),6);
    BOOST_CHECK(plan.freshEdges==std::set<ZhangGraphEdge>({a2,c2}));
    evaluatePlan(state,plan,next);
    // A complete physical restart still has a valid destination chart. It
    // carries independent priors for every arc and certifies no old relation.
    auto restarted=planFixture(state,old,old,old.edges);
    BOOST_REQUIRE(restarted.valid);
    BOOST_CHECK(restarted.survivingPhysicalRows.empty());
    BOOST_CHECK_EQUAL(restarted.independentSourceVariances.size(),2*old.edges.size());
    evaluatePlan(state,restarted,old);
    auto chord=*std::find_if(old.edges.begin(),old.edges.end(),[&](auto e){return !old.treeEdges.contains(e);});
    auto key=zhangTransportAmbiguityKey(E_ObsCode::L1C,chord);
    state.stateTransitionMap.erase(key); // pending retirement must never be reused
    plan=planFixture(state,old,old);
    BOOST_CHECK(!plan.survivingPhysicalRows.contains({E_ObsCode::L1C,chord}));
    BOOST_CHECK(plan.freshEdges.contains(chord));
    state.procNoiseMap[key]=0.1;
    plan=planFixture(state,old,old);
    BOOST_CHECK(!plan.valid);
    BOOST_CHECK_EQUAL(plan.failureReason,"NONCONSTANT_PHASE_PROCESS_REQUIRES_FULL_MODEL_TRANSPORT");
}

BOOST_AUTO_TEST_CASE(r51_pivot_and_new_leaf_are_separate_chart_operations)
{
    const ZhangGraphEdge a1{"A", SatSys("G01")}, a2{"A", SatSys("G02")};
    const ZhangGraphEdge b1{"B", SatSys("G01")}, b2{"B", SatSys("G02")};
    const ZhangGraphEdge c1{"NEW", SatSys("G01")};
    const auto old = zhangBuildSpanningTree({a1,a2,b1,b2}, "A", {a1,b1,b2});
    const auto transported = zhangBuildSpanningTree({a1,a2,b2}, "A", {a1,a2,b2});
    BOOST_REQUIRE(old.connected);
    const std::set<ZhangGraphEdge> eligible{a1,a2,b2,c1};
    BOOST_CHECK(!zhangActiveGraphBasisAfterTransform(transported, eligible).valid);
    const auto parts = partitionZhangEligibility(transported, eligible);
    BOOST_CHECK(parts.requiresAugmentation == std::set<ZhangGraphEdge>{c1});
    const auto extended = zhangActiveGraphBasisWithAugmentation(transported, eligible);
    BOOST_REQUIRE(extended.valid);
    BOOST_CHECK(extended.basis.edges.contains(a2)); // replacement bridge survives
    BOOST_CHECK(extended.basis.treeEdges.contains(c1));
    BOOST_CHECK_EQUAL(extended.basis.treeEdges.size(), transported.treeEdges.size()+1);
    for (const auto& edge : transported.treeEdges)
        BOOST_CHECK(extended.basis.treeEdges.contains(edge));
}

BOOST_AUTO_TEST_CASE(r51_simultaneous_tree_and_chord_slips_cannot_supply_replacement_paths)
{
    const ZhangGraphEdge a1{"A",SatSys("G01")},a2{"A",SatSys("G02")},
        b1{"B",SatSys("G01")},b2{"B",SatSys("G02")};
    auto old=zhangBuildSpanningTree({a1,a2,b1,b2},"A",{a1,a2,b1});
    auto incompleteEvent=zhangPlanPivotBeforeRetire(old,old.edges,{b1});
    BOOST_REQUIRE(incompleteEvent.connected);
    BOOST_CHECK(incompleteEvent.replacementEdges.contains(b2));
    auto completeEvent=zhangPlanPivotBeforeRetire(old,old.edges,{b1,b2});
    BOOST_CHECK(!completeEvent.connected);
    BOOST_CHECK(!completeEvent.survivingRepresentedEdges.contains(b2));
}

BOOST_AUTO_TEST_CASE(r51_new_nodes_cannot_silently_exchange_an_old_tree)
{
    const ZhangGraphEdge a1{"A", SatSys("G01")}, a2{"A", SatSys("G02")};
    const ZhangGraphEdge b1{"B", SatSys("G01")}, b2{"B", SatSys("G02")};
    auto represented = zhangBuildSpanningTree({a1,a2,b1,b2}, "A", {a1,b1,b2});
    auto extended = zhangActiveGraphBasisWithAugmentation(
        represented, {a1,a2,b2,{"NEW",SatSys("G02")},{"NEW",SatSys("G03")}});
    BOOST_REQUIRE(extended.valid);
    BOOST_CHECK(extended.retainedDatumOnlyEdges.contains(b1));
    BOOST_CHECK(extended.basis.treeEdges.contains(b1));
    represented.edges.erase(b1);
    BOOST_CHECK(!zhangActiveGraphBasisWithAugmentation(represented,{a1,a2,b2}).valid);
}

#include "common/zhangIntegerCandidateNis.hpp"
BOOST_AUTO_TEST_CASE(r48_c0_covariance_safety) {
 Eigen::MatrixXd q=Eigen::MatrixXd::Identity(2,2);q(1,1)=-1;
 auto a=assessZhangIntegerCandidateNis(Eigen::VectorXd::Zero(2),q,0.001);
 BOOST_CHECK(!a.valid);BOOST_CHECK_EQUAL(a.status,"NON_PSD_COVARIANCE");
 q(1,1)=-1e-15;a=assessZhangIntegerCandidateNis(Eigen::VectorXd::Zero(2),q,0.001);
 BOOST_CHECK(a.valid);BOOST_CHECK_EQUAL(a.rank,1);
 Eigen::VectorXd v=Eigen::VectorXd::Zero(2);v(1)=0.1;
 a=assessZhangIntegerCandidateNis(v,q,0.001);
 BOOST_CHECK(!a.valid);BOOST_CHECK_EQUAL(a.status,"DETERMINISTIC_RESIDUAL_CONFLICT");
 BOOST_CHECK(!assessZhangIntegerCandidateNis(v,q,0).valid);
}

#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangProductRelationSolver.hpp"
BOOST_AUTO_TEST_CASE(r48_c1_batch_matches_scalar) {
 for(int n=2;n<=18;++n) {
  ZhangExactMatrix h(1,ZhangExactVector(n));h[0][0]=2;h[0][1]=3;
  auto f=zhangExactAffineIntegerQuotient(h,{7},n);BOOST_REQUIRE(f.valid);
  ZhangExactMatrix columns(n,ZhangExactVector(f.quotientRank));
  for(int i=0;i<n;++i)for(int j=0;j<f.quotientRank;++j)columns[i][j]=f.kernelBasis[j][i];
  for(int j=0;j<f.quotientRank;++j) {
   ZhangExactVector e(f.quotientRank);e[j]=1;
   auto scalar=zhangIntegerRowLatticeContains(columns,e);
   BOOST_REQUIRE(scalar.contained);BOOST_CHECK(scalar.combination==f.quotientProjector[j]);
  }
 }
 BOOST_CHECK(!zhangExactAffineIntegerQuotient({{2,0}},{5},2).valid);
}

#include "common/zhangR48ProductDatum.hpp"
BOOST_AUTO_TEST_CASE(r48_c2_healthy_tree_identity_and_hard_rebuild) {
 std::set<ZhangGraphEdge> edges={{"A",SatSys("G01")},{"A",SatSys("G02")},
   {"B",SatSys("G01")},{"B",SatSys("G02")}};
 auto old=zhangBuildSpanningTree(edges,"A");
 auto replacement=zhangBuildSpanningTree(edges,"A",{{"B",SatSys("G02")}});
 auto p=proposeProductDatum(old,replacement,{},true);auditProductDatumTransport(p,old);
 BOOST_REQUIRE(zhangR51Enabled()?p.exactBasisTransport:p.identityTransport);BOOST_CHECK(commitProductDatum(p).treeEdges==(zhangR51Enabled()?replacement.treeEdges:old.treeEdges));
 auto e=*old.treeEdges.begin();p=proposeProductDatum(old,replacement,{e},true);auditProductDatumTransport(p,old);
 BOOST_CHECK(!p.identityTransport);BOOST_CHECK_EQUAL(p.status,"TRANSPORT_UNPROVEN");
}
#include <boost/test/unit_test.hpp>
#include "common/zhangR48SafePrefix.hpp"
BOOST_AUTO_TEST_CASE(r48_c4_late_failure_keeps_safe_prefix) {
 Eigen::VectorXd mu=Eigen::VectorXd::Zero(3);
 Eigen::MatrixXd q=Eigen::MatrixXd::Identity(3,3)*0.001;
 int calls=0;std::vector<std::string> logs;
 auto result=zhangR48SafePrefix(mu,q,{}, {},0.001,1e-6,
  [&](const Eigen::VectorXd& m,const Eigen::MatrixXd& c,double risk,bool merge) {
   ZhangSequentialShadowProposal p;++calls;
   if(merge)return p;
   p.valid=true;p.failureProbability=0;
   p.rows=zhangExactIdentityMatrix(m.size());p.values=ZhangExactVector(m.size());
   if(calls>1)p.values.back()=20;
   return p;
  },[&](const std::string& s){logs.push_back(s);},1,2,3);
 BOOST_CHECK_EQUAL(result.rows.size(),2);
 BOOST_CHECK(result.values==ZhangExactVector(2));
 BOOST_CHECK(result.reservedRisk>0 && result.reservedRisk<=0.001);
 BOOST_CHECK_EQUAL(result.mergeAttempts,1);
 bool rejected=false;for(auto& s:logs)rejected|=s.find("safe_prefix_retained=1")!=std::string::npos;
 BOOST_CHECK(rejected);
}
BOOST_AUTO_TEST_CASE(r48_c4_rescue_uses_common_snapshot_and_bounded_risk) {
 Eigen::VectorXd mu=Eigen::VectorXd::Zero(3);
 Eigen::MatrixXd q=Eigen::MatrixXd::Identity(3,3)*0.001;
 int calls=0;bool common=false;
 auto result=zhangR48SafePrefix(mu,q,{}, {},0.001,1e-6,
  [&](const Eigen::VectorXd& m,const Eigen::MatrixXd& c,double risk,bool merge) {
   ZhangSequentialShadowProposal p;++calls;p.valid=true;p.failureProbability=0;
   p.rows=zhangExactIdentityMatrix(m.size());p.values=ZhangExactVector(m.size());
   if(!merge && calls==2)p.values[0]=1;
   if(merge)common=c.diagonal().minCoeff()>0.0009;
   return p;
  },[](const std::string&){},1,2,3);
 BOOST_CHECK(common);BOOST_CHECK_EQUAL(result.rows.size(),3);
 BOOST_CHECK_EQUAL(result.mergeAccepted,1);BOOST_CHECK(result.reservedRisk<=0.001);
}
BOOST_AUTO_TEST_CASE(r48_c4_joint_affine_not_pairwise) {
 BOOST_CHECK(zhangR47AffineIntegerFeasible({{1,0},{0,1}},{0,0},2));
 BOOST_CHECK(!zhangR47AffineIntegerFeasible({{1,0},{0,1},{1,1}},{0,0,1},2));
}


#include "common/zhangR48Bridge.hpp"
BOOST_AUTO_TEST_CASE(r48_c5_bridge_preserves_integer_image_and_joint_gate) {
 Eigen::VectorXd mu(2);mu<<3,2;
 Eigen::MatrixXd q=Eigen::MatrixXd::Identity(2,2)*0.001;
 auto search=[](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double,bool) {
  ZhangSequentialShadowProposal p;p.valid=true;p.failureProbability=0;
  p.rows=zhangExactIdentityMatrix(m.size());for(int i=0;i<m.size();++i)p.values.push_back(std::llround(m(i)));return p;
 };
 auto b=zhangR48SearchBridge(mu,q,{{2,0},{0,1}},{},{},1e-5,1e-6,search);
 BOOST_REQUIRE(b.accepted);BOOST_CHECK_EQUAL(b.rank,2);
 auto frame=zhangR47CompileProductSearchFrame({{2,0},{0,1}},b.rows,b.values,2);
 ZhangExactInteger value;BOOST_REQUIRE(zhangR47ProductConsequence(frame,{2,0},0,value));
 BOOST_CHECK(value==6);
 auto already=zhangR48SearchBridge(mu,q,{{2,0},{0,1}},b.rows,b.values,1e-5,1e-6,search);
 BOOST_CHECK_EQUAL(already.status,"ALREADY_PROVEN");BOOST_CHECK_EQUAL(already.spent,0);
 auto bad=[](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double,bool) {
  ZhangSequentialShadowProposal p;p.valid=true;p.failureProbability=0;
  p.rows=zhangExactIdentityMatrix(m.size());p.values=ZhangExactVector(m.size(),100);return p;
 };
 auto reject=zhangR48SearchBridge(mu,q,{{2,0},{0,1}},{},{},1e-5,1e-6,bad);
 BOOST_CHECK(!reject.accepted);BOOST_CHECK_EQUAL(reject.spent,1e-5);
}

BOOST_AUTO_TEST_CASE(r48_c6_product_image_lift_avoids_full_left_inverse) {
 ZhangExactMatrix t={{2,0,0},{0,3,0}},h={{1,1,1}};
 auto f=zhangR47CompileProductSearchFrame(t,h,{7},3);
 BOOST_REQUIRE(f.valid);BOOST_REQUIRE(f.affine);BOOST_CHECK(f.affine->quotientProjector.empty());
 BOOST_CHECK_EQUAL(f.searchRank,2);
 auto old=zhangExactAffineIntegerQuotient(h,{7},3);
 ZhangExactMatrix k(3,ZhangExactVector(old.quotientRank));
 for(int c=0;c<3;++c)for(int j=0;j<old.quotientRank;++j)k[c][j]=old.kernelBasis[j][c];
 auto image=zhangPrimitiveImageCoordinates(zhangExactMultiply(t,k));
 BOOST_CHECK(zhangExactMultiply(f.projector,k)==image.primitiveRows);
 auto feasible=zhangExactAffineIntegerQuotient({{2,0}},{6},2,ZhangExactQuotientWork::FEASIBILITY_ONLY);
 BOOST_CHECK(feasible.valid);BOOST_CHECK(feasible.kernelBasis.empty());
}
BOOST_AUTO_TEST_CASE(r48_c6_low_rank_marginal_matches_full_conditioning) {
 Eigen::VectorXd m(3);m<<0.1,0.2,0.3;
 Eigen::MatrixXd q(3,3);q<<2,.3,.1,.3,1,.2,.1,.2,3;
 Eigen::MatrixXd h(1,3);h<<1,1,0;Eigen::VectorXd v(1);v<<0;
 Eigen::MatrixXd j(2,3);j<<1,0,-1,0,2,0;
 ZhangR48MarginalWorkspace w(m,q);
 auto small=w.project(j,h,v);auto full=zhangConditionPosteriorEffectiveIntegers(m,q,h,v);
 BOOST_REQUIRE(small.valid && full.valid);
 BOOST_CHECK_SMALL((small.mean-j*full.mean).norm(),1e-10);
 BOOST_CHECK_SMALL((small.covariance-j*full.covariance*j.transpose()).norm(),1e-10);
 auto again=w.project(j,h,v);BOOST_CHECK_EQUAL(w.targetHits,1);BOOST_CHECK_EQUAL(w.decompositions,1);
 Eigen::VectorXd changed(1);changed<<1;auto different=w.project(j,h,changed);
 BOOST_CHECK_EQUAL(w.decompositions,1);BOOST_CHECK_EQUAL(w.hits,1);BOOST_CHECK((different.mean-small.mean).norm()>0.1);
}

#include "common/zhangActiveGraphBasis.hpp"
BOOST_AUTO_TEST_CASE(r48_active_basis_retains_transformed_datum_without_promoting_eligibility)
{
    const ZhangGraphEdge a9{"ROOT", SatSys("G09")}, a11{"ROOT", SatSys("G11")},
        b9{"DAV1", SatSys("G09")}, b11{"DAV1", SatSys("G11")};
    const std::set<ZhangGraphEdge> representedEdges{a9, a11, b9, b11};
    const auto transformed = zhangBuildSpanningTree(representedEdges, "ROOT", {a9,a11,b11});
    BOOST_REQUIRE(transformed.treeEdges.contains(b11));
    const std::set<ZhangGraphEdge> eligible{a9,a11,b9};
    auto broken = transformed; broken.edges = eligible;
    BOOST_CHECK(!std::includes(broken.edges.begin(), broken.edges.end(),
        broken.treeEdges.begin(), broken.treeEdges.end()));
    const auto repaired = zhangActiveGraphBasisAfterTransform(transformed, eligible);
    BOOST_REQUIRE_MESSAGE(repaired.valid, repaired.failureReason);
    BOOST_CHECK(repaired.basis.treeEdges == transformed.treeEdges);
    BOOST_CHECK(repaired.basis.edges == representedEdges);
    BOOST_CHECK(repaired.retainedDatumOnlyEdges == std::set<ZhangGraphEdge>{b11});
    BOOST_CHECK(!eligible.contains(b11));
    BOOST_CHECK(repaired.basis.receivers == transformed.receivers);
    BOOST_CHECK(repaired.basis.satellites == transformed.satellites);
    auto oldCycle = zhangFundamentalCycle(transformed,b9);
    BOOST_CHECK(zhangFundamentalCycle(repaired.basis,b9) == oldCycle);
    // Removing a non-tree edge must not retain its stale stochastic support.
    auto treeOnly = zhangActiveGraphBasisAfterTransform(transformed,{a9,a11});
    BOOST_REQUIRE(treeOnly.valid);
    BOOST_CHECK(!treeOnly.basis.edges.contains(b9));
    BOOST_CHECK(treeOnly.basis.treeEdges == transformed.treeEdges);
}

BOOST_AUTO_TEST_CASE(r48_active_basis_rejects_unrepresented_tree_and_new_chart)
{
    const ZhangGraphEdge a9{"ROOT", SatSys("G09")}, a11{"ROOT", SatSys("G11")},
        b9{"DAV1", SatSys("G09")}, b11{"DAV1", SatSys("G11")};
    auto basis = zhangBuildSpanningTree({a9,a11,b9,b11},"ROOT",{a9,a11,b11});
    auto invalid = basis;invalid.edges.erase(b11);
    auto result = zhangActiveGraphBasisAfterTransform(invalid,{a9,a11,b9});
    BOOST_CHECK(!result.valid);
    BOOST_CHECK_EQUAL(result.failureReason,"REPRESENTED_TREE_EDGE_MISSING");
    auto newNode=basis.edges;newNode.insert({"NEW",SatSys("G09")});
    result=zhangActiveGraphBasisAfterTransform(basis,newNode);
    BOOST_CHECK(!result.valid);
    BOOST_CHECK_EQUAL(result.failureReason,"ACTIVE_SUPPORT_REQUIRES_NEW_COORDINATE_TRANSFORM");
    result=zhangActiveGraphBasisAfterTransform(basis,basis.edges);
    BOOST_REQUIRE(result.valid);
    BOOST_CHECK(result.retainedDatumOnlyEdges.empty());
    BOOST_CHECK(result.basis.treeEdges==basis.treeEdges);
}

BOOST_AUTO_TEST_CASE(r49_partial_bridge_is_shadow_only_and_not_dual_certificate) {
 Eigen::VectorXd mu(2);mu<<3,2;
 Eigen::MatrixXd q=Eigen::MatrixXd::Identity(2,2)*0.001;
 auto partial=[](const Eigen::VectorXd&,const Eigen::MatrixXd&,double,bool) {
  ZhangSequentialShadowProposal p;p.valid=true;p.failureProbability=0;
  p.rows={{1,0}};p.values={3};return p;
 };
 auto formal=zhangR48SearchBridge(mu,q,{{1,0},{0,1}},{},{},1e-5,1e-6,partial);
 BOOST_CHECK(!formal.accepted);
 auto shadow=zhangR48SearchBridge(mu,q,{{1,0},{0,1}},{},{},1e-5,1e-6,partial,nullptr,nullptr,true);
 BOOST_REQUIRE(shadow.accepted);BOOST_CHECK_EQUAL(shadow.status,"PARTIAL_INTEGER_PROGRESS");
 auto domain=zhangR47CompileProductSearchFrame({{1,0},{0,1}},shadow.rows,shadow.values,2);
 ZhangExactInteger value;
 BOOST_CHECK(zhangR47ProductConsequence(domain,{1,0},0,value));
 BOOST_CHECK(!zhangR47ProductConsequence(domain,{0,1},0,value));
}

#include "common/zhangR49ConstraintNis.hpp"
BOOST_AUTO_TEST_CASE(r49_empty_product_domain_skips_blas_and_does_not_authorize) {
 const VectorXd mean=VectorXd::Zero(2588);
 const MatrixXd covariance=MatrixXd::Identity(2588,2588);
 const auto empty=zhangR49ConstraintNis(mean,covariance,MatrixXd(0,2588),VectorXd(0),1e-6);
 BOOST_CHECK(!empty.valid);
 BOOST_CHECK_EQUAL(empty.status,"NOT_EVALUATED_EMPTY_DOMAIN");
 MatrixXd row=MatrixXd::Zero(1,2588);row(0,0)=1;
 VectorXd value(1);value<<0.1;
 const auto actual=zhangR49ConstraintNis(mean,covariance,row,value,1e-6);
 const auto reference=assessZhangIntegerCandidateNis(value-row*mean,row*covariance*row.transpose(),1e-6);
 BOOST_REQUIRE(actual.valid);BOOST_CHECK_EQUAL(actual.nis,reference.nis);
 BOOST_CHECK_EQUAL(actual.threshold,reference.threshold);
}

#include "common/zhangPhaseArcDecision.hpp"
#include "common/zhangRatioGate.hpp"
BOOST_AUTO_TEST_CASE(r51_followup_lli_masks_half_cycle_and_recovery)
{
    for (unsigned lli : {0u,1u,2u,3u,4u}) {
        ZhangSignalTracking state;
        zhangDecidePhaseTracking(state,0,true,1,0,90);
        auto result=zhangDecidePhaseTracking(state,lli,true,1,30,90);
        BOOST_CHECK_EQUAL(result.breakArc,(lli&3)!=0);
        BOOST_CHECK_EQUAL(result.useCurrentPhase,(lli&2)==0);
        BOOST_CHECK(result.useCurrentCode);
        if (lli&2) {
            auto sequence=state.eventSequence;
            result=zhangDecidePhaseTracking(state,lli,true,1,60,90);
            BOOST_CHECK(!result.breakArc);BOOST_CHECK_EQUAL(state.eventSequence,sequence);
            result=zhangDecidePhaseTracking(state,0,true,1,90,90);
            BOOST_CHECK(result.breakArc && result.useCurrentPhase);
        }
    }
}
BOOST_AUTO_TEST_CASE(r51_followup_tracking_is_per_signal)
{
    ZhangSignalTracking l1,l2;
    zhangDecidePhaseTracking(l1,0,true,1,0,90);
    zhangDecidePhaseTracking(l2,0,true,2,0,90);
    BOOST_CHECK(zhangDecidePhaseTracking(l1,1,true,1,30,90).breakArc);
    BOOST_CHECK(!zhangDecidePhaseTracking(l2,0,true,2,30,90).breakArc);
    BOOST_CHECK(!zhangDecidePhaseTracking(l2,0,false,2,60,90).useCurrentPhase);
    BOOST_CHECK(!zhangDecidePhaseTracking(l2,0,true,2,90,90).breakArc);
    BOOST_CHECK(zhangDecidePhaseTracking(l2,0,true,2,210,90).gap);
}
BOOST_AUTO_TEST_CASE(r51_followup_predictor_isolated_spike_does_not_poison_history)
{
    ZhangCombinationDetector detector;
    for(int i=0;i<4;i++) BOOST_REQUIRE(zhangInspectCombination(detector,i*30,100,1e-4,1e-6,4,90).useSample);
    auto spike=zhangInspectCombination(detector,120,100.3,1e-4,1e-6,4,90);
    BOOST_CHECK(spike.suspect && !spike.breakArc && !spike.useSample);
    BOOST_CHECK_EQUAL(detector.accepted.size(),4);
    auto returned=zhangInspectCombination(detector,150,100,1e-4,1e-6,4,90);
    BOOST_CHECK(returned.useSample && !returned.breakArc && !returned.suspect);
    BOOST_CHECK_EQUAL(detector.accepted.size(),5);
}
BOOST_AUTO_TEST_CASE(r51_followup_predictor_persistent_step_and_irregular_sampling)
{
    ZhangCombinationDetector detector;
    for(double t : {0.,20.,55.,80.})
        BOOST_REQUIRE(zhangInspectCombination(detector,t,100+.0001*t,1e-4,1e-6,4,90).useSample);
    auto first=zhangInspectCombination(detector,110,100+.011+.4,1e-4,1e-6,4,90);
    auto second=zhangInspectCombination(detector,155,100+.0155+.4,1e-4,1e-6,4,90);
    BOOST_CHECK(first.suspect && !first.breakArc);
    BOOST_CHECK(second.breakArc && second.useSample);
    BOOST_CHECK_EQUAL(detector.accepted.size(),1);
    ZhangCombinationDetector noisy;
    zhangInspectCombination(noisy,0,100,.04,1e-6,4,90);
    auto low=zhangInspectCombination(noisy,30,100.08,.04,1e-6,4,90);
    BOOST_CHECK(low.useSample && !low.breakArc);
}
BOOST_AUTO_TEST_CASE(r51_followup_gf_mw_complementary_integer_slips)
{
    const double l1=.190293672798,l2=.244210213425;
    for (auto [n1,n2] : {std::pair{1,1},std::pair{9,7}}) {
        ZhangCombinationDetector gf,mw;
        for(int i=0;i<4;i++) {
            zhangInspectCombination(gf,i*30,0,1e-6,0,4,90);
            zhangInspectCombination(mw,i*30,0,.01,0,4,90);
        }
        const double dg=l1*n1-l2*n2,dm=n1-n2;
        zhangInspectCombination(gf,120,dg,1e-6,0,4,90);
        zhangInspectCombination(mw,120,dm,.01,0,4,90);
        auto a=zhangInspectCombination(gf,150,dg,1e-6,0,4,90);
        auto b=zhangInspectCombination(mw,150,dm,.01,0,4,90);
        BOOST_CHECK(a.breakArc || b.breakArc);
        if(n1==n2) BOOST_CHECK(!b.breakArc);
        if(n1==9) BOOST_CHECK(b.breakArc);
    }
}
BOOST_AUTO_TEST_CASE(r51_followup_ratio_boundary_contract)
{
    BOOST_CHECK(zhangEvaluateRatioGate(true,true,1,3,3).accepted);
    BOOST_CHECK(!zhangEvaluateRatioGate(true,true,1,2.99,3).accepted);
    BOOST_CHECK(!zhangEvaluateRatioGate(false,true,0,1,3).executed);
    BOOST_CHECK(!zhangEvaluateRatioGate(true,false,0,1,3).accepted);
    BOOST_CHECK(zhangEvaluateRatioGate(true,true,0,1,3).accepted);
    BOOST_CHECK(!zhangEvaluateRatioGate(true,true,0,0,3).accepted);
    BOOST_CHECK(!zhangEvaluateRatioGate(true,true,1,1,3).accepted);
    BOOST_CHECK(!zhangEvaluateRatioGate(true,true,1,std::numeric_limits<double>::infinity(),3).accepted);
    BOOST_CHECK(!zhangEvaluateRatioGate(true,true,1,3,1).accepted);
}
BOOST_AUTO_TEST_CASE(r51_followup_l1_retirement_preserves_l2_physical_covariance)
{
    const ZhangGraphEdge a1{"A",SatSys("G01")},a2{"A",SatSys("G02")},
        b1{"B",SatSys("G01")},b2{"B",SatSys("G02")};
    auto basis=zhangBuildSpanningTree({a1,a2,b1,b2},"A");
    auto state=coordinateFixture(basis);
    std::map<ZhangGraphEdge,std::set<E_ObsCode>> signals{{a1,{E_ObsCode::L1C}}};
    auto plan=zhangPlanGraphCoordinateTransport(state,E_Sys::GPS,basis,basis,{a1},
        {{E_ObsCode::L1C,.1902936728},{E_ObsCode::L2W,.2442102134}},
        [](const auto&,auto){return 100.;},[](const auto&,auto){return 25.;},&signals);
    BOOST_REQUIRE_MESSAGE(plan.valid,plan.failureReason);
    BOOST_CHECK_EQUAL(plan.survivingPhysicalRows.size(),7);
    BOOST_CHECK(plan.survivingPhysicalRows.contains({E_ObsCode::L2W,a1}));
    BOOST_CHECK(!plan.survivingPhysicalRows.contains({E_ObsCode::L1C,a1}));
    BOOST_CHECK_EQUAL(plan.independentSourceVariances.size(),1);
    evaluatePlan(state,plan,basis);
}

BOOST_AUTO_TEST_CASE(r51_diagnostic_switch_preserves_rejected_prefix_and_rescue) {
 auto run=[](bool diagnostics) {
  setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS",diagnostics?"1":"0",1);
  Eigen::VectorXd mu=Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd q=Eigen::MatrixXd::Identity(3,3)*0.001;
  int calls=0;std::string log;
  auto result=zhangR48SafePrefix(mu,q,{}, {},0.001,1e-6,
   [&](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double,bool merge) {
    ZhangSequentialShadowProposal p;++calls;if(merge)return p;
    p.valid=true;p.failureProbability=0;p.rows=zhangExactIdentityMatrix(m.size());
    p.values=ZhangExactVector(m.size());if(calls>1)p.values.back()=20;return p;
   },[&](const std::string& s){log+=s;},1,2,3);
  BOOST_CHECK_EQUAL(log.find("conditional_innovation=")!=std::string::npos,diagnostics);
  BOOST_CHECK_EQUAL(log.find("covariance=")!=std::string::npos,diagnostics);
  return result;
 };
 auto on=run(true),off=run(false);unsetenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS");
 BOOST_CHECK(on.rows==off.rows);BOOST_CHECK(on.values==off.values);
 BOOST_CHECK_EQUAL(on.attempts,off.attempts);BOOST_CHECK_EQUAL(on.mergeAttempts,off.mergeAttempts);
 BOOST_CHECK_EQUAL(on.reservedRisk,off.reservedRisk);BOOST_CHECK_EQUAL(on.status,off.status);
}
