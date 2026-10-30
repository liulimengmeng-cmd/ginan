#include <boost/test/unit_test.hpp>
#include "common/zhangActiveGraphBasis.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangProductPublicationGate.hpp"
#include "common/zhangProductConditionalMap.hpp"
#include "common/zhangProductChartIdentity.hpp"
#include "common/zhangPhaseContinuity.hpp"
#include "common/zhangR51CommittedBundle.hpp"
#include "pea/zhangPppAr.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>

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

BOOST_AUTO_TEST_CASE(r51_network_fixed_diagnostic_skips_held_membership)
{
    // No product-conditioned branch exists: the R51 FIXED row is only a
    // numerical control and cannot become a user-AR product.
    BOOST_CHECK(!zhangNeedsNamedProductIntegerMembership(
        true, false, "FIXED"));
    BOOST_CHECK(!zhangFormalPppArProductSolution("FIXED"));

    // Keep the pre-existing query in other modes.  PRODUCT_FIXED takes its
    // authority from the product certificate, not this legacy FIXED lookup.
    BOOST_CHECK(zhangNeedsNamedProductIntegerMembership(
        false, false, "FIXED"));
    BOOST_CHECK(zhangNeedsNamedProductIntegerMembership(
        true, true, "FIXED"));
    BOOST_CHECK(!zhangNeedsNamedProductIntegerMembership(
        true, true, "PRODUCT_FIXED"));
    BOOST_CHECK(!zhangNeedsNamedProductIntegerMembership(
        true, false, "FLOAT"));
    BOOST_CHECK(zhangFormalPppArProductSolution("PRODUCT_FIXED"));
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

BOOST_AUTO_TEST_CASE(r51_publication_gate_initial_and_unchanged_identity)
{
	ZhangProductPublicationEvidence evidence;
	evidence.nativeProductValid = true;
	evidence.currentSpatialCertificateValid = true;
	evidence.numericalStepValid = true;
	const auto initial = zhangProductPublicationGate(evidence);
	BOOST_CHECK(initial.route == ZhangProductPublicationRoute::NEW_SEGMENT);
	BOOST_CHECK(!initial.requireUserReset);
	BOOST_CHECK_EQUAL(initial.reason, "INITIAL_PRODUCT_SEGMENT");

	evidence.priorPublishedIdentity = true;
	evidence.samePhysicalSegment = true;
	evidence.sameProductFunctional = true;
	evidence.sameBackendChart = true;
	const auto unchanged = zhangProductPublicationGate(evidence);
	BOOST_CHECK(unchanged.route == ZhangProductPublicationRoute::SAME_SEGMENT);
	BOOST_CHECK(!unchanged.requireUserReset);
	// Neither a mode name nor an unnecessary fresh bridge is a prerequisite
	// when the current product identity is genuinely unchanged.
	BOOST_CHECK(!evidence.exactSameTimeMapProven);
	BOOST_CHECK(!evidence.currentEpochBridgeProven);
}

BOOST_AUTO_TEST_CASE(r51_new_observed_arc_is_current_support_without_old_integer_inheritance)
{
	const ZhangGraphEdge fresh{"ALBH", SatSys("G06")};
	const ZhangGraphEdge grace{"ALBH", SatSys("G07")};
	const ZhangGraphEdge qcBad{"ALBH", SatSys("G08")};
	const std::set<ZhangGraphEdge> observed{fresh, qcBad};
	BOOST_CHECK(zhangCurrentObservedProductEdge(
		fresh, observed, {}, {}));
	BOOST_CHECK(!zhangCurrentObservedProductEdge(
		grace, observed, {}, {}));
	BOOST_CHECK(!zhangCurrentObservedProductEdge(
		qcBad, observed, {qcBad}, {}));
	const auto old = zhangBuildSpanningTree({fresh, grace}, "ALBH");
	const auto current = zhangBuildSpanningTree({fresh, grace}, "ALBH");
	auto proposal = proposeProductDatum(old, current, {fresh}, true);
	auditProductDatumTransport(proposal, old);
	BOOST_CHECK(!proposal.transportProven());
}

BOOST_AUTO_TEST_CASE(r51_publication_gate_exact_chart_map_and_changed_functional_bridge)
{
	ZhangProductPublicationEvidence evidence;
	evidence.nativeProductValid = true;
	evidence.currentSpatialCertificateValid = true;
	evidence.numericalStepValid = true;
	evidence.priorPublishedIdentity = true;
	evidence.samePhysicalSegment = true;
	evidence.sameProductFunctional = true;
	evidence.sameBackendChart = false;

	const auto noProof = zhangProductPublicationGate(evidence);
	BOOST_CHECK(noProof.route == ZhangProductPublicationRoute::NEW_SEGMENT);
	BOOST_CHECK(noProof.requireUserReset);
	BOOST_CHECK_EQUAL(noProof.reason, "UNPROVEN_PRODUCT_HANDOFF");
	evidence.currentEpochBridgeProven = true;
	BOOST_CHECK(zhangProductPublicationGate(evidence).route ==
		ZhangProductPublicationRoute::NEW_SEGMENT);
	evidence.currentEpochBridgeProven = false;
	evidence.exactSameTimeMapProven = true;
	const auto sameFunctional = zhangProductPublicationGate(evidence);
	BOOST_CHECK(sameFunctional.route == ZhangProductPublicationRoute::SAME_SEGMENT);
	BOOST_CHECK_EQUAL(sameFunctional.reason,"EXACT_SAME_FUNCTIONAL_TRANSPORT");
	evidence.currentEpochBridgeProven = true;
	const auto proven = zhangProductPublicationGate(evidence);
	BOOST_CHECK(proven.route == ZhangProductPublicationRoute::SAME_SEGMENT);
	BOOST_CHECK(!proven.requireUserReset);

	// A changed functional needs a conditional product map and a bridge from
	// the current epoch. An unconditional map claim alone is insufficient.
	evidence.sameBackendChart = true;
	evidence.sameProductFunctional = false;
	BOOST_CHECK(zhangProductPublicationGate(evidence).route ==
		ZhangProductPublicationRoute::NEW_SEGMENT);
	evidence.conditionalProductMapProven = true;
	const auto conditional = zhangProductPublicationGate(evidence);
	BOOST_CHECK(conditional.route == ZhangProductPublicationRoute::SAME_SEGMENT);
	BOOST_CHECK_EQUAL(conditional.reason,"CURRENT_EPOCH_CONDITIONAL_HANDOFF");
	evidence.currentEpochBridgeProven = false;
	BOOST_CHECK(zhangProductPublicationGate(evidence).route ==
		ZhangProductPublicationRoute::NEW_SEGMENT);
}

BOOST_AUTO_TEST_CASE(r51_publication_gate_slip_reset_and_numeric_failure)
{
	ZhangProductPublicationEvidence evidence;
	evidence.nativeProductValid = true;
	evidence.currentSpatialCertificateValid = true;
	evidence.numericalStepValid = true;
	evidence.priorPublishedIdentity = true;
	evidence.samePhysicalSegment = false;
	evidence.sameProductFunctional = true;
	evidence.sameBackendChart = true;
	evidence.exactSameTimeMapProven = true;
	evidence.currentEpochBridgeProven = true;
	const auto slipped = zhangProductPublicationGate(evidence);
	BOOST_CHECK(slipped.route == ZhangProductPublicationRoute::NEW_SEGMENT);
	BOOST_CHECK(slipped.requireUserReset);
	BOOST_CHECK_EQUAL(slipped.reason, "PHYSICAL_SEGMENT_CHANGED");

	// Retiring an old graph arc is not necessarily a satellite product break:
	// a conditional FLOAT old-aux map plus this epoch's integer bridge may
	// certify the same published product across that event.
	evidence.currentEpochBridgeProven = false;
	evidence.conditionalProductMapProven = true;
	BOOST_CHECK(zhangProductPublicationGate(evidence).route ==
		ZhangProductPublicationRoute::NEW_SEGMENT);
	evidence.currentEpochBridgeProven = true;
	const auto retiredArcHandoff = zhangProductPublicationGate(evidence);
	BOOST_CHECK(retiredArcHandoff.route == ZhangProductPublicationRoute::SAME_SEGMENT);
	BOOST_CHECK(!retiredArcHandoff.requireUserReset);
	BOOST_CHECK_EQUAL(retiredArcHandoff.reason,
		"CURRENT_EPOCH_CONDITIONAL_HANDOFF");
	evidence.conditionalProductMapProven = false;

	// A large unexplained step is not legitimised by relabelling the segment.
	evidence.numericalStepValid = false;
	const auto badStep = zhangProductPublicationGate(evidence);
	BOOST_CHECK(badStep.route == ZhangProductPublicationRoute::SUSPENDED);
	BOOST_CHECK(!badStep.requireUserReset);
	BOOST_CHECK_EQUAL(badStep.reason, "NUMERICAL_STEP_INVALID");
	evidence.numericalStepValid = true;
	evidence.samePhysicalSegment = true;
	evidence.currentSpatialCertificateValid = false;
	BOOST_CHECK(zhangProductPublicationGate(evidence).route ==
		ZhangProductPublicationRoute::SUSPENDED);
	evidence.currentSpatialCertificateValid = true;
	evidence.nativeProductValid = false;
	BOOST_CHECK(zhangProductPublicationGate(evidence).route ==
		ZhangProductPublicationRoute::SUSPENDED);
}

BOOST_AUTO_TEST_CASE(r51_dual_frequency_publication_resets_both_user_arcs_on_either_new_segment)
{
	const ZhangProductPublicationDecision same{
		ZhangProductPublicationRoute::SAME_SEGMENT,false,"same"};
	const ZhangProductPublicationDecision fresh{
		ZhangProductPublicationRoute::NEW_SEGMENT,true,"new"};
	const ZhangProductPublicationDecision suspended{
		ZhangProductPublicationRoute::SUSPENDED,false,"suspended"};
	const auto bothSame = zhangAggregateDualFrequencyPublication(same,same);
	BOOST_CHECK(bothSame.dualFrequencyArUsable);
	BOOST_CHECK(!bothSame.resetBothUserAmbiguities);
	const auto firstNew = zhangAggregateDualFrequencyPublication(fresh,same);
	const auto secondNew = zhangAggregateDualFrequencyPublication(same,fresh);
	BOOST_CHECK(firstNew.dualFrequencyArUsable);
	BOOST_CHECK(secondNew.dualFrequencyArUsable);
	BOOST_CHECK(firstNew.resetBothUserAmbiguities);
	BOOST_CHECK(secondNew.resetBothUserAmbiguities);
	BOOST_CHECK_EQUAL(firstNew.reason,"DUAL_FREQUENCY_NEW_SEGMENT");
	// The other signal keeps its independent SAME_SEGMENT diagnostic route.
	BOOST_CHECK(same.route == ZhangProductPublicationRoute::SAME_SEGMENT);
	BOOST_CHECK(!zhangAggregateDualFrequencyPublication(same,suspended)
		.dualFrequencyArUsable);
	BOOST_CHECK(!zhangAggregateDualFrequencyPublication(suspended,same)
		.resetBothUserAmbiguities);
	const auto mixed = zhangAggregateDualFrequencyPublication(fresh,suspended);
	BOOST_CHECK(!mixed.dualFrequencyArUsable);
	BOOST_CHECK(mixed.resetBothUserAmbiguities);
	BOOST_CHECK_EQUAL(mixed.reason,
		"DUAL_FREQUENCY_SUSPENDED_WITH_NEW_SEGMENT");
	const auto missingSecond = zhangAggregateDualFrequencyPublication(
		std::optional<ZhangProductPublicationDecision>{same},
		std::optional<ZhangProductPublicationDecision>{});
	BOOST_CHECK(!missingSecond.dualFrequencyArUsable);
	BOOST_CHECK(!missingSecond.resetBothUserAmbiguities);
	BOOST_CHECK_EQUAL(missingSecond.reason,"DUAL_FREQUENCY_SIGNAL_MISSING");
	const auto missingFirstAfterNew = zhangAggregateDualFrequencyPublication(
		std::optional<ZhangProductPublicationDecision>{},
		std::optional<ZhangProductPublicationDecision>{fresh});
	BOOST_CHECK(!missingFirstAfterNew.dualFrequencyArUsable);
	BOOST_CHECK(missingFirstAfterNew.resetBothUserAmbiguities);
	BOOST_CHECK_EQUAL(missingFirstAfterNew.reason,
		"DUAL_FREQUENCY_SIGNAL_MISSING_WITH_NEW_SEGMENT");
}

BOOST_AUTO_TEST_CASE(r51_pure_graph_rebase_preserves_product_mean_and_full_covariance)
{
	const ZhangGraphEdge a1{"A",SatSys("G01")}, a2{"A",SatSys("G02")},
		b1{"B",SatSys("G01")}, b2{"B",SatSys("G02")};
	const std::set<ZhangGraphEdge> edges{a1,a2,b1,b2};
	const auto before = zhangBuildSpanningTree(edges,"A",{a1,b1,b2});
	const auto after = zhangBuildSpanningTree(edges,"A",{a1,a2,b2});
	BOOST_REQUIRE(before.connected && after.connected);
	const auto initial = coordinateFixture(before);
	const auto plan = planFixture(initial,before,after);
	BOOST_REQUIRE_MESSAGE(plan.valid,plan.failureReason);
	BOOST_CHECK(plan.independentSourceVariances.empty());
	const auto rebased = evaluatePlan(initial,plan,after);

	// Independent physical product model at root receiver A. The G02 row
	// includes its old chord ambiguity, while in the new tree it does not.
	// No row or covariance is taken from the graph transport's proof ledger.
	auto productRows = [&](const KFState& state,const ZhangGraphBasis& basis)
	{
		Eigen::MatrixXd rows = Eigen::MatrixXd::Zero(4,state.x.size());
		int row = 0;
		for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
			for (const auto sat : {SatSys("G01"),SatSys("G02")})
			{
				const int current = row++;
				KFKey key;
				key.type = KF::PHASE_BIAS;
				key.Sat = sat;
				key.num = static_cast<int>(code);
				const auto index = state.kfIndexMap.find(key);
				BOOST_REQUIRE(index != state.kfIndexMap.end());
				rows(current,index->second) = 1;
				const ZhangGraphEdge rootEdge{"A",sat};
				if (!basis.treeEdges.contains(rootEdge))
				{
					KFKey ambiguity;
					ambiguity.type = KF::AMBIGUITY;
					ambiguity.Sat = sat;
					ambiguity.str = "A";
					ambiguity.num = static_cast<int>(code);
					const auto ambiguityIndex = state.kfIndexMap.find(ambiguity);
					BOOST_REQUIRE(ambiguityIndex != state.kfIndexMap.end());
					const double wavelength = code == E_ObsCode::L1C
						? 0.1902936728 : 0.2442102134;
					rows(current,ambiguityIndex->second) = wavelength;
				}
			}
		return rows;
	};
	const auto oldRows = productRows(initial,before);
	const auto newRows = productRows(rebased,after);
	const Eigen::Vector4d offsets(0.125,-0.25,0.375,-0.5);
	const Eigen::VectorXd oldMean = oldRows * initial.x + offsets;
	const Eigen::VectorXd newMean = newRows * rebased.x;
	const Eigen::MatrixXd oldQ = oldRows * initial.P * oldRows.transpose();
	const Eigen::MatrixXd newQ = newRows * rebased.P * newRows.transpose();
	BOOST_CHECK(std::abs(oldQ(0,1)) > 1e-3); // Exercise off-diagonal propagation.
	BOOST_CHECK_SMALL((newMean+offsets-oldMean).norm(),1e-10);
	BOOST_CHECK_SMALL((newQ-oldQ).norm(),1e-9);
	auto namedRows = [](const KFState& state,const Eigen::MatrixXd& matrix)
	{
		std::vector<ZhangCoordinateRow> rows(matrix.rows());
		for (const auto& [key,index] : state.kfIndexMap)
			for (int i=0;i<matrix.rows();i++)
				if (matrix(i,index) != 0) rows[i][key] = matrix(i,index);
		return rows;
	};
	const auto oldNamed = namedRows(initial,oldRows);
	const auto newNamed = namedRows(rebased,newRows);
	std::vector<std::pair<ZhangCoordinateRow,double>> targets;
	for (int i=0;i<4;i++)
	{
		const auto proof = zhangVerifyGraphFunctionalTransport(
			plan,oldNamed[i],newNamed[i],offsets(i),offsets(i));
		BOOST_CHECK_MESSAGE(proof.sameTimeMapProven,proof.failureReason);
		targets.push_back({newNamed[i],offsets(i)});
	}
	const auto projected = zhangGraphFunctionalMoments(initial,plan,targets);
	BOOST_REQUIRE_MESSAGE(projected.valid,projected.failureReason);
	BOOST_CHECK_SMALL((projected.mean-oldMean).norm(),1e-10);
	BOOST_CHECK_SMALL((projected.covariance-oldQ).norm(),1e-9);
	const auto jointAgreement = zhangVerifyGraphJointMoments(
		oldMean,oldQ,projected.mean,projected.covariance);
	BOOST_REQUIRE_MESSAGE(jointAgreement.proven,
		jointAgreement.failureReason);
	// The production validator must check the whole joint Q, not just each
	// satellite's marginal sigma. Tamper with one symmetric cross term only.
	Eigen::MatrixXd badJointQ = projected.covariance;
	badJointQ(0,1) += 0.01;
	badJointQ(1,0) += 0.01;
	const auto badCross = zhangVerifyGraphJointMoments(
		oldMean,oldQ,projected.mean,badJointQ);
	BOOST_CHECK(!badCross.proven);
	BOOST_CHECK_EQUAL(badCross.failureReason,
		"JOINT_PRODUCT_COVARIANCE_TRANSPORT_MISMATCH");
	Eigen::VectorXd badJointMean = projected.mean;
	badJointMean(0) += 0.01;
	const auto badMean = zhangVerifyGraphJointMoments(
		oldMean,oldQ,badJointMean,projected.covariance);
	BOOST_CHECK(!badMean.proven);
	BOOST_CHECK_EQUAL(badMean.failureReason,
		"JOINT_PRODUCT_MEAN_TRANSPORT_MISMATCH");
	const auto badCoefficient = [&]()
	{
		auto row = newNamed[1];
		row.begin()->second += 0.1;
		return zhangVerifyGraphFunctionalTransport(
			plan,oldNamed[1],row,offsets(1),offsets(1));
	}();
	BOOST_CHECK(!badCoefficient.sameTimeMapProven);
	BOOST_CHECK_EQUAL(badCoefficient.failureReason,
		"PRODUCT_FUNCTIONAL_PULLBACK_MISMATCH");
	const auto badOffset = zhangVerifyGraphFunctionalTransport(
		plan,oldNamed[1],newNamed[1],offsets(1),offsets(1)+0.01);
	BOOST_CHECK(!badOffset.sameTimeMapProven);
	BOOST_CHECK_EQUAL(badOffset.failureReason,"AFFINE_OFFSET_MISMATCH");
	const auto exactMap = zhangVerifyGraphFunctionalTransport(
		plan,oldNamed[1],newNamed[1],offsets(1),offsets(1));
	BOOST_REQUIRE(exactMap.sameTimeMapProven);
	ZhangProductPublicationEvidence priorOnly;
	priorOnly.nativeProductValid = true;
	priorOnly.currentSpatialCertificateValid = true;
	priorOnly.priorPublishedIdentity = true;
	priorOnly.samePhysicalSegment = true;
	priorOnly.sameProductFunctional = false;
	priorOnly.sameBackendChart = false;
	priorOnly.exactSameTimeMapProven = exactMap.sameTimeMapProven;
	priorOnly.numericalStepValid = true;
	// A changed functional needs a current-epoch integer bridge. An exact
	// floating-state pullback and an old integer prior are insufficient.
	const auto oldIntegerOnly = zhangProductPublicationGate(priorOnly);
	BOOST_CHECK(oldIntegerOnly.route == ZhangProductPublicationRoute::NEW_SEGMENT);
	BOOST_CHECK(oldIntegerOnly.requireUserReset);
	priorOnly.currentEpochBridgeProven = true;
	BOOST_CHECK(zhangProductPublicationGate(priorOnly).route ==
		ZhangProductPublicationRoute::NEW_SEGMENT);
	priorOnly.conditionalProductMapProven = true;
	BOOST_CHECK(zhangProductPublicationGate(priorOnly).route ==
		ZhangProductPublicationRoute::SAME_SEGMENT);
	KFKey external; external.type=KF::IONO_STEC; external.str="external";
	const auto oldCross = oldRows * initial.P.col(initial.kfIndexMap.at(external));
	const auto newCross = newRows * rebased.P.col(rebased.kfIndexMap.at(external));
	BOOST_CHECK_SMALL((newCross-oldCross).norm(),1e-10);

	// Retiring both G02 arcs adds a fresh independent phase prior. Its FLOAT
	// mean and covariance still contain the old satellite phase-bias state;
	// that is not authority to reuse an old integer datum or segment identity.
	const auto slipped = planFixture(initial,before,after,{a2,b2});
	BOOST_REQUIRE_MESSAGE(slipped.valid,slipped.failureReason);
	const auto fresh = zhangPullbackGraphFunctional(slipped,newNamed[1]);
	BOOST_REQUIRE_MESSAGE(fresh.valid,fresh.failureReason);
	BOOST_CHECK(!fresh.independentPriorTerms.empty());
	const auto slipProof = zhangVerifyGraphFunctionalTransport(
		slipped,oldNamed[1],newNamed[1],offsets(1),offsets(1));
	BOOST_CHECK(!slipProof.sameTimeMapProven);
	BOOST_CHECK_EQUAL(slipProof.failureReason,
		"NEW_PRIOR_ENTERED_PUBLIC_FUNCTIONAL");
	ZhangCoordinateRow externalRow{{external,1}};
	const auto freshMoments = zhangGraphFunctionalMoments(initial,slipped,
		{{newNamed[1],0},{externalRow,0}});
	BOOST_REQUIRE_MESSAGE(freshMoments.valid,freshMoments.failureReason);
	const auto oldSatIndex = initial.kfIndexMap.at(
		zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G02")));
	const auto externalIndex = initial.kfIndexMap.at(external);
	BOOST_CHECK_SMALL(freshMoments.mean(0)-initial.x(oldSatIndex),1e-10);
	BOOST_CHECK_SMALL(freshMoments.covariance(0,0)-
		(initial.P(oldSatIndex,oldSatIndex)+100.0),1e-9);
	BOOST_CHECK_SMALL(freshMoments.covariance(0,1)-
		initial.P(oldSatIndex,externalIndex),1e-10);
}

BOOST_AUTO_TEST_CASE(r51_new_publication_segment_preserves_current_affine_map_but_slip_does_not)
{
	ZhangPhaseContinuityState continuity;
	continuity.counter = 4;
	continuity.datumVersion = 8;
	continuity.iod = 11;
	continuity.integerShiftCycles = 3;
	continuity.fractionalShiftCycles = 0.25;
	continuity.hasFixedDatum = true;
	GTime handoffTime; handoffTime.bigTime = 1000;
	continuity.beginPublicationSegment(handoffTime,"unproven_current_epoch_product_handoff");
	BOOST_CHECK_EQUAL(continuity.counter,5);
	BOOST_CHECK_EQUAL(continuity.datumVersion,9);
	BOOST_CHECK_EQUAL(continuity.iod,12);
	BOOST_CHECK(continuity.validFrom == handoffTime);
	BOOST_CHECK_EQUAL(continuity.integerShiftCycles,3);
	BOOST_CHECK_CLOSE(continuity.fractionalShiftCycles,0.25,1e-10);
	BOOST_CHECK(continuity.hasFixedDatum);

	GTime slipTime; slipTime.bigTime = 1030;
	continuity.reinitialise(slipTime,"physical_arc_slip",0);
	BOOST_CHECK_EQUAL(continuity.counter,6);
	BOOST_CHECK_EQUAL(continuity.datumVersion,10);
	BOOST_CHECK_EQUAL(continuity.integerShiftCycles,0);
	BOOST_CHECK_SMALL(continuity.fractionalShiftCycles,1e-15);
	BOOST_CHECK(!continuity.hasFixedDatum);
}

BOOST_AUTO_TEST_CASE(r51_current_epoch_integer_bridge_uses_only_frozen_posterior_decisions)
{
	const std::string arc = "L1C|R|G02|V1";
	ZhangProductPhysicalCycleChart chart;
	chart.columns = 1;
	chart.chordIdentities[0] = arc;
	chart.expansions[0][arc] = 1;
	chart.activePhysicalIdentities.insert(arc);
	auto proof = std::make_shared<ZhangIntegerDecisionProof>();
	proof->id = "current-epoch-decision";
	proof->originalStatement = "2*N=2";
	proof->observationProvenance = "runtime|PPP_COMMITTED|E1405216800|posterior-0";
	proof->conditionalFailureBound = 1e-4;
	ZhangCurrentEpochEvidence current;
	constexpr long int epoch = 1405216800;
	BOOST_REQUIRE(current.bind("runtime",E_Sys::GPS,epoch,
		"posterior-0","chart-0",chart));
	BOOST_CHECK(current.add({{arc,2}},ZhangExactInteger(2),{proof}));
	BOOST_CHECK(!current.add({{arc,2}},ZhangExactInteger(3),{proof}));
	ZhangProductIntegerTransition difference;
	difference.physicalEdges = {{"R",SatSys("G02")}};
	difference.physicalArcVersions = {1};
	difference.coefficients = {1};
	difference.affineOffsetCycles = 3;
	difference.valid = true;
	const auto exact = current.query(E_ObsCode::L1C,epoch,
		"posterior-0","chart-0",difference);
	BOOST_REQUIRE_MESSAGE(exact.proven,exact.reason);
	BOOST_CHECK(exact.value == ZhangExactInteger(4)); // 2N=2, then +3.
	BOOST_CHECK_EQUAL(exact.reason,
		"EXACT_CURRENT_AFFINE_INTEGER_CONSEQUENCE");
	BOOST_REQUIRE_EQUAL(exact.proofIds.size(),1);
	BOOST_CHECK_EQUAL(exact.proofIds.front(),proof->id);
	ZhangCurrentEpochEvidence nonIntegral;
	BOOST_REQUIRE(nonIntegral.bind("runtime",E_Sys::GPS,epoch,
		"posterior-0","chart-0",chart));
	BOOST_REQUIRE(nonIntegral.add({{arc,2}},ZhangExactInteger(3),{proof}));
	const auto parityConflict = nonIntegral.query(E_ObsCode::L1C,epoch,
		"posterior-0","chart-0",difference);
	BOOST_CHECK(!parityConflict.proven); // 2N=3 has no integer N.
	BOOST_CHECK_EQUAL(parityConflict.reason,
		"BRIDGE_NOT_EXACT_INTEGER_CONSEQUENCE");

	BOOST_CHECK(!current.query(E_ObsCode::L1C,epoch+30,
		"posterior-0","chart-0",difference).proven);
	BOOST_CHECK(!current.query(E_ObsCode::L1C,epoch,
		"posterior-1","chart-0",difference).proven);
	BOOST_CHECK(!current.query(E_ObsCode::L1C,epoch,
		"posterior-0","chart-1",difference).proven);
	BOOST_CHECK(!current.bind("runtime",E_Sys::GPS,epoch+30,
		"posterior-1","chart-1",chart));
	BOOST_CHECK(!current.query(E_ObsCode::L2W,epoch,
		"posterior-0","chart-0",difference).proven);
	const auto oldArc = difference;
	difference.physicalArcVersions = {2};
	BOOST_CHECK(!current.query(E_ObsCode::L1C,epoch,
		"posterior-0","chart-0",difference).proven);
	difference = oldArc;
	difference.physicalEdges.clear();
	difference.physicalArcVersions.clear();
	difference.coefficients.clear();
	const auto empty = current.query(E_ObsCode::L1C,epoch,
		"posterior-0","chart-0",difference);
	BOOST_CHECK(!empty.proven);
	BOOST_CHECK_EQUAL(empty.reason,"EMPTY_BRIDGE_HAS_NO_INTEGER_EVIDENCE");

	// A new epoch gets a fresh evidence object bound to its own posterior.
	// The accepted old integer is not silently copied into this new object.
	ZhangCurrentEpochEvidence next;
	BOOST_REQUIRE(next.bind("runtime",E_Sys::GPS,epoch+30,
		"posterior-1","chart-1",chart));
	const auto noPrior = next.query(E_ObsCode::L1C,epoch+30,
		"posterior-1","chart-1",oldArc);
	BOOST_CHECK(!noPrior.proven);
	BOOST_CHECK_EQUAL(noPrior.reason,"NO_CURRENT_ACCEPTED_INTEGER_ROWS");
}

BOOST_AUTO_TEST_CASE(r51_current_epoch_evidence_rejects_branch_or_prior_epoch_proof_ancestry)
{
	const std::string arc = "L1C|R|G02|V1";
	ZhangProductPhysicalCycleChart chart;
	chart.columns = 1;
	chart.chordIdentities[0] = arc;
	chart.expansions[0][arc] = 1;
	chart.activePhysicalIdentities.insert(arc);
	constexpr long int epoch = 1405216800;
	const std::string root = "runtime";
	auto makeProof = [](const std::string& id, const std::string& source,
		const ZhangDecisionProofs& parents = ZhangDecisionProofs{}) {
		auto proof = std::make_shared<ZhangIntegerDecisionProof>();
		proof->id = id;
		proof->originalStatement = "exact accepted integer row";
		proof->observationProvenance = source;
		proof->conditionalFailureBound = 1e-4;
		proof->parents = parents;
		return proof;
	};
	const auto currentParent = makeProof("network-parent",
		"runtime|PPP_COMMITTED|E1405216800");
	// The product decision is searched in a fixed branch, but its observation
	// authority is the protected current-epoch FLOAT root.
	const auto currentProduct = makeProof("r47-product-family",
		"runtime|PPP_COMMITTED|E1405216800|R51_RATIO_ONLY_NO_FFR_GUARANTEE",
		{currentParent});
	ZhangCurrentEpochEvidence accepted;
	BOOST_REQUIRE(accepted.bind(root, E_Sys::GPS, epoch,
		"float-posterior", "physical-chart", chart));
	BOOST_CHECK(accepted.add({{arc, 1}}, ZhangExactInteger(1),
		{currentProduct}));

	const auto branchProduct = makeProof("branch-product-family",
		"runtime/branch/fixed-transaction|PPP_COMMITTED|E1405216800",
		{currentParent});
	ZhangCurrentEpochEvidence wrongBranch;
	BOOST_REQUIRE(wrongBranch.bind(root, E_Sys::GPS, epoch,
		"float-posterior", "physical-chart", chart));
	BOOST_CHECK(!wrongBranch.add({{arc, 1}}, ZhangExactInteger(1),
		{branchProduct}));
	BOOST_CHECK(wrongBranch.rows.empty());

	const auto oldParent = makeProof("prior-epoch-parent",
		"runtime|PPP_COMMITTED|E1405216770");
	const auto currentWithOldParent = makeProof("current-with-old-parent",
		"runtime|PPP_COMMITTED|E1405216800", {oldParent});
	ZhangCurrentEpochEvidence wrongEpoch;
	BOOST_REQUIRE(wrongEpoch.bind(root, E_Sys::GPS, epoch,
		"float-posterior", "physical-chart", chart));
	BOOST_CHECK(!wrongEpoch.add({{arc, 1}}, ZhangExactInteger(1),
		{currentWithOldParent}));
	BOOST_CHECK(wrongEpoch.rows.empty());
}

BOOST_AUTO_TEST_CASE(r51_temporary_product_handoff_row_retains_cross_covariance_until_removal)
{
	KFKey first; first.type=KF::PHASE_BIAS; first.Sat=SatSys("G01");
	KFKey second; second.type=KF::PHASE_BIAS; second.Sat=SatSys("G02");
	KFKey external; external.type=KF::IONO_STEC; external.str="external";
	KFKey mixed; mixed.type=KF::PHASE_BIAS; mixed.Sat=SatSys("G03");
	KFKey auxiliary; auxiliary.type=KF::PRODUCT_HANDOFF;
	auxiliary.Sat=SatSys("G03"); auxiliary.str="old-tree-edge";
	KFState old;
	old.kfIndexMap = {{first,0},{second,1},{external,2}};
	old.x = Eigen::Vector3d(2,-1,0.5);
	old.P.resize(3,3);
	old.P << 4,1,0.2,
		1,9,-0.3,
		0.2,-0.3,1;
	ZhangGraphCoordinatePlan plan;
	plan.valid = true;
	plan.transform[mixed] = {{first,1},{second,1}};
	plan.transform[external] = {{external,1}};
	plan.transform[auxiliary] = {{first,1},{second,1}};
	const auto moments = zhangGraphFunctionalMoments(old,plan,
		{{{{mixed,1}},0},{{{external,1}},0},{{{auxiliary,1}},0}});
	BOOST_REQUIRE_MESSAGE(moments.valid,moments.failureReason);
	// Hand reference: s=x0+x1, Var(s)=4+9+2*1=15,
	// Cov(s,x2)=0.2-0.3=-0.1. The temporary row is an exact copy of s.
	Eigen::Vector3d referenceMean(1,0.5,1);
	Eigen::Matrix3d referenceQ;
	referenceQ << 15,-0.1,15,
		-0.1,1,-0.1,
		15,-0.1,15;
	BOOST_CHECK_SMALL((moments.mean-referenceMean).norm(),1e-12);
	BOOST_CHECK_SMALL((moments.covariance-referenceQ).norm(),1e-12);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(referenceQ);
	BOOST_CHECK(eigen.eigenvalues().minCoeff() > -1e-12);
	BOOST_CHECK_SMALL(eigen.eigenvalues().minCoeff(),1e-12);
	// Atomic retirement selects only live coordinates. No conditional
	// covariance or zero-variance pseudo-observation is substituted for it.
	Eigen::Matrix<double,2,3> keep;
	keep << 1,0,0,
		0,1,0;
	const Eigen::Vector2d liveMean = keep*moments.mean;
	const Eigen::Matrix2d liveQ =
		keep*moments.covariance*keep.transpose();
	Eigen::Matrix2d referenceLiveQ;
	referenceLiveQ << 15,-0.1,-0.1,1;
	BOOST_CHECK_SMALL((liveMean-Eigen::Vector2d(1,0.5)).norm(),1e-12);
	BOOST_CHECK_SMALL((liveQ-referenceLiveQ).norm(),1e-12);
}

BOOST_AUTO_TEST_CASE(r51_committed_epoch_bundle_requires_atomic_verified_pair)
{
	namespace fs = std::filesystem;
	const fs::path parent = fs::temp_directory_path();
	const fs::path root = parent / ("zhang-r51-bundle-unit-" +
		std::to_string(std::chrono::steady_clock::now()
			.time_since_epoch().count()));
	BOOST_REQUIRE(fs::create_directory(root));
	struct Cleanup
	{
		fs::path root;
		fs::path parent;
		~Cleanup()
		{
			if (root.parent_path() != parent ||
				root.filename().string().rfind("zhang-r51-bundle-unit-",0) != 0)
				return;
			std::error_code error;
			fs::remove_all(root,error);
		}
	} cleanup{root,parent};
	const fs::path pending = root / "1405216800_ENHANCED_1.pending";
	const fs::path committed = root / "1405216800_ENHANCED_1";
	BOOST_REQUIRE(fs::create_directory(pending));
	auto write = [](const fs::path& path,const std::string& contents)
	{
		std::ofstream output(path,std::ios::binary|std::ios::trunc);
		output << contents;
		output.flush();
		if (!output) throw std::runtime_error("R51_BUNDLE_UNIT_WRITE_FAILED");
	};
	write(pending/"products.csv","solution,satellite\nPRODUCT_FIXED,G02\n");
	write(pending/"covariance.csv","h\n0\n");
	write(pending/"certificate.txt","proof=current\n");
	boost::json::object manifest;
	manifest["epoch"] = "1405216800";
	manifest["route"] = "ENHANCED";
	manifest["root"] = "float-root";
	manifest["ar"] = true;
	manifest["products_bytes"] = static_cast<std::int64_t>(
		fs::file_size(pending/"products.csv"));
	manifest["covariance_bytes"] = static_cast<std::int64_t>(
		fs::file_size(pending/"covariance.csv"));
	write(pending/"COMMIT.json",boost::json::serialize(manifest)+"\n");
	std::ostringstream checksums;
	for (const auto* name : {"products.csv","covariance.csv",
			"certificate.txt","COMMIT.json"})
		checksums << zhangR51BundleFileHash(pending/name)
			<< "  " << name << '\n';
	write(pending/"SHA256SUMS",checksums.str());

	const auto hidden = zhangR51VerifyCommittedBundle(pending,true);
	BOOST_CHECK(!hidden.valid);
	BOOST_CHECK_EQUAL(hidden.failureReason,"R51_BUNDLE_NOT_COMMITTED");
	fs::rename(pending,committed);
	const auto valid = zhangR51VerifyCommittedBundle(committed,true);
	BOOST_REQUIRE_MESSAGE(valid.valid,valid.failureReason);
	BOOST_CHECK(valid.ar);
	BOOST_CHECK(valid.products == committed/"products.csv");
	BOOST_CHECK(valid.covariance == committed/"covariance.csv");
	BOOST_CHECK_EQUAL(valid.epoch,"1405216800");

	// Same-size corruption is detected by SHA, not only file size.
	write(committed/"covariance.csv","h\n1\n");
	const auto corrupt = zhangR51VerifyCommittedBundle(committed,true);
	BOOST_CHECK(!corrupt.valid);
	BOOST_CHECK_EQUAL(corrupt.failureReason,"R51_BUNDLE_HASH_MISMATCH");
	fs::remove(committed/"covariance.csv");
	const auto incomplete = zhangR51VerifyCommittedBundle(committed,true);
	BOOST_CHECK(!incomplete.valid);
	BOOST_CHECK_EQUAL(incomplete.failureReason,"R51_BUNDLE_COVARIANCE_MISSING");
}

BOOST_AUTO_TEST_CASE(r51_product_target_chart_is_not_the_kf_s_tree)
{
	const SatSys g01("G01");
	const SatSys g02("G02");
	const std::set<ZhangGraphEdge> edges{
		{"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}};
	const ZhangGraphBasis stateBasis = zhangBuildSpanningTree(
		edges, "R0", {{"R0", g01}, {"R0", g02}, {"R1", g01}});
	const ZhangGraphBasis targetChart = zhangBuildSpanningTree(
		edges, "R0", {{"R0", g01}, {"R1", g01}, {"R1", g02}});
	BOOST_REQUIRE(stateBasis.connected && targetChart.connected);
	BOOST_REQUIRE(stateBasis.treeEdges != targetChart.treeEdges);

	// The target chart selects product integer functionals; the estimator's
	// sole S-basis remains stateBasis and is never committed by this proposal.
	auto proposal = proposeProductDatum(
		targetChart, targetChart, {}, true);
	auditProductDatumTransport(proposal, targetChart);
	const auto committedTarget = commitProductDatum(proposal);
	BOOST_CHECK(committedTarget.treeEdges == targetChart.treeEdges);
	BOOST_CHECK(committedTarget.treeEdges != stateBasis.treeEdges);
	const auto target = zhangBuildSatelliteProductTarget(
		stateBasis, committedTarget, g01);
	BOOST_REQUIRE_MESSAGE(target.valid, target.failureReason);
	BOOST_REQUIRE_EQUAL(target.matrix.size(), 1);
	BOOST_CHECK(target.matrix[0] == ZhangExactVector({1}));
	const ZhangGraphEdge chord = target.currentChords.at(0);
	const std::map<ZhangGraphEdge, ZhangExactInteger> physical{
		{{"R0", g01}, 0}, {{"R0", g02}, 0},
		{{"R1", g01}, 0}, {{"R1", g02}, 5}};
	auto satelliteDifference = [&](const ZhangGraphBasis& basis)
		-> ZhangExactInteger
	{
		const auto audit = zhangCanonicalIntegerAudit(basis);
		BOOST_REQUIRE(audit.valid);
		ZhangExactVector treeValues;
		for (const auto& edge : audit.treeEdges)
			treeValues.push_back(physical.at(edge));
		const auto nodes = zhangExactMatrixTimesColumn(
			audit.treeInverse, treeValues);
		const std::size_t satelliteOffset = basis.receivers.size() - 1;
		return nodes[satelliteOffset + 1] - nodes[satelliteOffset];
	};
	ZhangExactInteger cycle = 0;
	for (const auto& [edge, coefficient] :
		 zhangFundamentalCycle(stateBasis, chord))
		cycle += coefficient * physical.at(edge);
	BOOST_CHECK(satelliteDifference(committedTarget) ==
		satelliteDifference(stateBasis) + target.matrix[0][0] * cycle);

	// The product phase is a linear functional of the single KF posterior;
	// its covariance must include the phase/ambiguity cross term.
	const KFState posterior = coordinateFixture(stateBasis);
	const int phaseIndex = posterior.kfIndexMap.at(
		zhangTransportSatelliteKey(E_ObsCode::L1C, g02));
	const int ambiguityIndex = posterior.kfIndexMap.at(
		zhangTransportAmbiguityKey(E_ObsCode::L1C, chord));
	const double lambda = 0.1902936728;
	Eigen::VectorXd projection = Eigen::VectorXd::Zero(posterior.x.size());
	projection(phaseIndex) = 1;
	projection(ambiguityIndex) = lambda;
	const double projectedMean = projection.dot(posterior.x);
	const double projectedVariance = projection.dot(posterior.P * projection);
	const double expandedVariance = posterior.P(phaseIndex, phaseIndex)
		+ lambda * lambda * posterior.P(ambiguityIndex, ambiguityIndex)
		+ 2 * lambda * posterior.P(phaseIndex, ambiguityIndex);
	BOOST_CHECK_SMALL(projectedMean -
		(posterior.x(phaseIndex) + lambda * posterior.x(ambiguityIndex)),
		1e-12);
	BOOST_CHECK_SMALL(projectedVariance - expandedVariance, 1e-12);
	BOOST_CHECK(std::abs(2 * lambda *
		posterior.P(phaseIndex, ambiguityIndex)) > 1e-6);
	const auto noChange = zhangBuildSatelliteProductTarget(
		stateBasis, stateBasis, g01);
	BOOST_REQUIRE_MESSAGE(noChange.valid, noChange.failureReason);
	BOOST_CHECK(noChange.matrix[0] == ZhangExactVector({0}));
}

BOOST_AUTO_TEST_CASE(r51_changed_product_tree_requires_exact_physical_conditional_map)
{
	const SatSys g01("G01"), g02("G02");
	const ZhangGraphEdge a1{"R0", g01}, a2{"R0", g02};
	const ZhangGraphEdge b1{"R1", g01}, b2{"R1", g02};
	const std::set<ZhangGraphEdge> edges{a1, a2, b1, b2};
	const auto stateBasis = zhangBuildSpanningTree(
		edges, "R0", {a1, a2, b1});
	const auto productBasis = zhangBuildSpanningTree(
		edges, "R0", {a1, b1, b2});
	BOOST_REQUIRE(stateBasis.connected && productBasis.connected);
	const std::map<ZhangGraphEdge, int> versions{
		{a1, 7}, {a2, 7}, {b1, 7}, {b2, 7}};
	const auto oldFunctions = zhangBuildProductIntegerFunctionals(
		stateBasis, versions, g01);
	const auto newFunctions = zhangBuildProductIntegerFunctionals(
		productBasis, versions, g01);
	BOOST_REQUIRE(oldFunctions.contains(g02) && newFunctions.contains(g02));
	const auto transition = zhangProductIntegerFunctionalDifference(
		oldFunctions.at(g02), newFunctions.at(g02));
	BOOST_REQUIRE(transition.valid);
	const double lambda = 0.1902936728;
	const auto phase = zhangTransportSatelliteKey(E_ObsCode::L1C, g02);
	const auto ambiguity = zhangTransportAmbiguityKey(E_ObsCode::L1C, b2);
	ZhangCoordinateRow oldRow{{phase, 1}};
	ZhangCoordinateRow writerRow{{phase, 1}, {ambiguity, lambda}};
	auto physical = [&](const ZhangGraphEdge& edge, int version,
		ZhangCoordinateRow& row, std::string& reason)
	{
		if (version != 7 || !edges.contains(edge))
		{
			reason = "VERSION_NOT_CURRENT";
			return false;
		}
		row.clear();
		if (edge.receiver != "R0")
			row[zhangTransportReceiverKey(E_Sys::GPS,
				E_ObsCode::L1C, edge.receiver)] = 1;
		row[zhangTransportSatelliteKey(E_ObsCode::L1C,
			edge.satellite)] = 1;
		if (!stateBasis.treeEdges.contains(edge))
			row[zhangTransportAmbiguityKey(E_ObsCode::L1C, edge)] = lambda;
		return true;
	};
	const auto proof = zhangProveProductConditionalMap(
		oldRow, writerRow, transition, lambda, true, physical);
	BOOST_REQUIRE_MESSAGE(proof.valid, proof.reason);
	BOOST_CHECK(proof.physicalEdges.size() == transition.coefficients.size());
	const auto posterior = coordinateFixture(stateBasis);
	const int p = posterior.kfIndexMap.at(phase);
	const int k = posterior.kfIndexMap.at(ambiguity);
	const double fullVariance = posterior.P(p, p)
		+ lambda * lambda * posterior.P(k, k)
		+ 2 * lambda * posterior.P(p, k);
	BOOST_CHECK(std::abs(posterior.P(p, k)) > 1e-6);
	Eigen::VectorXd l = Eigen::VectorXd::Zero(posterior.x.size());
	l(p) = 1;
	l(k) = lambda;
	BOOST_CHECK_SMALL(l.dot(posterior.P * l) - fullVariance, 1e-12);
}

BOOST_AUTO_TEST_CASE(r51_retired_product_chord_needs_correlated_auxiliary)
{
	const SatSys g01("G01"), g02("G02");
	const ZhangGraphEdge a1{"R0", g01}, a2{"R0", g02};
	const ZhangGraphEdge b1{"R1", g01}, b2{"R1", g02};
	const auto oldStateBasis = zhangBuildSpanningTree(
		{a1, a2, b1, b2}, "R0", {a1, a2, b1});
	const auto oldProductBasis = zhangBuildSpanningTree(
		{a1, a2, b1, b2}, "R0", {a1, b1, b2});
	const auto nextBasis = zhangBuildSpanningTree(
		{a1, a2, b1}, "R0", {a1, a2, b1});
	BOOST_REQUIRE(oldStateBasis.connected && oldProductBasis.connected
		&& nextBasis.connected);
	const std::map<ZhangGraphEdge, int> oldVersions{
		{a1, 7}, {a2, 7}, {b1, 7}, {b2, 7}};
	const std::map<ZhangGraphEdge, int> newVersions{
		{a1, 7}, {a2, 7}, {b1, 7}};
	const auto oldFunctions = zhangBuildProductIntegerFunctionals(
		oldProductBasis, oldVersions, g01);
	const auto newFunctions = zhangBuildProductIntegerFunctionals(
		nextBasis, newVersions, g01);
	BOOST_REQUIRE(oldFunctions.contains(g02) && newFunctions.contains(g02));
	const auto transition = zhangProductIntegerFunctionalDifference(
		oldFunctions.at(g02), newFunctions.at(g02));
	BOOST_REQUIRE(transition.valid);
	const double lambda = 0.1902936728;
	const auto frozen = coordinateFixture(oldStateBasis);
	const auto plan = zhangPlanGraphCoordinateTransport(
		frozen, E_Sys::GPS, oldStateBasis, nextBasis, {b2},
		{{E_ObsCode::L1C, lambda}, {E_ObsCode::L2W, 0.2442102134}},
		[](const auto&, auto) { return 100.; },
		[](const auto&, auto) { return 25.; },
		nullptr, true, &oldVersions, &oldProductBasis.treeEdges);
	BOOST_REQUIRE_MESSAGE(plan.valid, plan.failureReason);
	const auto auxiliary = plan.retainedOldArcAuxiliaryKeys.find(
		{E_ObsCode::L1C, b2});
	BOOST_REQUIRE(auxiliary != plan.retainedOldArcAuxiliaryKeys.end());
	const auto& auxRow = plan.transform.at(auxiliary->second);
	BOOST_CHECK(auxRow.contains(
		zhangTransportAmbiguityKey(E_ObsCode::L1C, b2)));
	const auto post = evaluatePlan(frozen, plan, nextBasis);
	const auto phase = zhangTransportSatelliteKey(E_ObsCode::L1C, g02);
	const auto receiver = zhangTransportReceiverKey(
		E_Sys::GPS, E_ObsCode::L1C, "R1");
	const auto ambiguity = zhangTransportAmbiguityKey(E_ObsCode::L1C, b2);
	const ZhangCoordinateRow oldProductRow{{phase, 1}, {ambiguity, lambda}};
	const ZhangCoordinateRow transportedOldRow{
		{auxiliary->second, 1}, {receiver, -1}};
	const ZhangCoordinateRow writerRow{{phase, 1}};
	const auto before = zhangGraphFunctionalMoments(
		frozen, plan, {{transportedOldRow, 0}});
	BOOST_REQUIRE_MESSAGE(before.valid, before.failureReason);
	const auto meanAfter = post.x(post.kfIndexMap.at(auxiliary->second))
		- post.x(post.kfIndexMap.at(receiver));
	const int aux = post.kfIndexMap.at(auxiliary->second);
	const int rec = post.kfIndexMap.at(receiver);
	const double covarianceAfter = post.P(aux, aux) + post.P(rec, rec)
		- 2 * post.P(aux, rec);
	const int oldPhase = frozen.kfIndexMap.at(phase);
	const int oldAmbiguity = frozen.kfIndexMap.at(ambiguity);
	const double oldProductVariance = frozen.P(oldPhase, oldPhase)
		+ lambda * lambda * frozen.P(oldAmbiguity, oldAmbiguity)
		+ 2 * lambda * frozen.P(oldPhase, oldAmbiguity);
	BOOST_CHECK_SMALL(before.mean(0) - meanAfter, 1e-10);
	BOOST_CHECK_SMALL(before.covariance(0, 0) - covarianceAfter, 1e-9);
	BOOST_CHECK_SMALL(covarianceAfter - oldProductVariance, 1e-9);
	BOOST_CHECK(std::abs(frozen.P(oldPhase, oldAmbiguity)) > 1e-6);
	BOOST_CHECK_SMALL(meanAfter -
		(frozen.x(oldPhase) + lambda * frozen.x(oldAmbiguity)), 1e-10);
	auto physical = [&](const ZhangGraphEdge& edge, int version,
		ZhangCoordinateRow& row, std::string& reason)
	{
		if (version != 7)
		{
			reason = "ARC_VERSION_MISMATCH";
			return false;
		}
		row.clear();
		if (edge == b2)
		{
			row[auxiliary->second] = 1;
			return true;
		}
		if (!nextBasis.edges.contains(edge))
		{
			reason = "RETIRED_PRODUCT_ARC_AUXILIARY_MISSING";
			return false;
		}
		if (edge.receiver != "R0")
			row[zhangTransportReceiverKey(E_Sys::GPS,
				E_ObsCode::L1C, edge.receiver)] = 1;
		row[zhangTransportSatelliteKey(E_ObsCode::L1C,
			edge.satellite)] = 1;
		return true;
	};
	const auto positive = zhangProveProductConditionalMap(
		transportedOldRow, writerRow, transition, lambda, true,
		physical);
	BOOST_REQUIRE_MESSAGE(positive.valid, positive.reason);
	auto missingAux = [&](const ZhangGraphEdge& edge, int version,
		ZhangCoordinateRow& row, std::string& reason)
	{
		if (edge == b2)
		{
			reason = "RETIRED_PRODUCT_ARC_AUXILIARY_MISSING";
			return false;
		}
		return physical(edge, version, row, reason);
	};
	const auto negative = zhangProveProductConditionalMap(
		transportedOldRow, writerRow, transition, lambda, true,
		missingAux);
	BOOST_CHECK(!negative.valid);
	BOOST_CHECK_EQUAL(negative.reason,
		"RETIRED_PRODUCT_ARC_AUXILIARY_MISSING");
	(void)oldProductRow;
}

BOOST_AUTO_TEST_CASE(r51_product_chart_only_identity_requires_versioned_unbroken_arcs)
{
	const SatSys g01("G01"), g02("G02");
	const ZhangGraphEdge a1{"R0", g01}, a2{"R0", g02};
	const ZhangGraphEdge b1{"R1", g01}, b2{"R1", g02};
	const std::set<ZhangGraphEdge> edges{a1, a2, b1, b2};
	const auto s = zhangBuildSpanningTree(
		edges, "R0", {a1, a2, b1});
	const auto p = zhangBuildSpanningTree(
		edges, "R0", {a1, b1, b2});
	BOOST_REQUIRE(s.connected && p.connected);
	std::map<ZhangGraphEdge, int> oldTreeVersions{
		{a1, 7}, {a2, 7}, {b1, 7}};
	std::map<ZhangGraphEdge, int> currentVersions{
		{a1, 7}, {a2, 7}, {b1, 7}, {b2, 7}};
	auto check = [&](const std::set<ZhangGraphEdge>& breaks = {})
	{
		return zhangCheckProductChartIdentityReceiptPreconditions(
			s, s, s, p, oldTreeVersions, currentVersions,
			edges, breaks);
	};
	BOOST_CHECK_MESSAGE(check().valid, check().reason);
	oldTreeVersions.erase(a2);
	BOOST_CHECK_EQUAL(check().reason,
		"OLD_PRODUCT_TREE_ARC_VERSION_MISSING");
	oldTreeVersions[a2] = 6;
	BOOST_CHECK_EQUAL(check().reason,
		"OLD_PRODUCT_TREE_ARC_VERSION_CHANGED");
	oldTreeVersions[a2] = 7;
	BOOST_CHECK_EQUAL(check({b2}).reason,
		"PHYSICAL_ARC_BREAK_REQUIRES_CORRELATED_TRANSFORM");
	currentVersions.erase(a1);
	BOOST_CHECK_EQUAL(check().reason, "OLD_S_TREE_ARC_NOT_CURRENT");
}

BOOST_AUTO_TEST_CASE(r51_old_aux_overlap_requires_two_admitted_signals_and_same_numeric_row)
{
	const SatSys g02("G02");
	const auto l1Phase = zhangTransportSatelliteKey(E_ObsCode::L1C, g02);
	const auto l2Phase = zhangTransportSatelliteKey(E_ObsCode::L2W, g02);
	KFKey l1Aux;
	l1Aux.type = KF::PRODUCT_HANDOFF;
	l1Aux.Sat = g02;
	l1Aux.str = "R1@4";
	l1Aux.num = static_cast<int>(E_ObsCode::L1C);
	KFKey l2Aux = l1Aux;
	l2Aux.num = static_cast<int>(E_ObsCode::L2W);
	std::array<ZhangCoordinateRow, 2> subject{
		ZhangCoordinateRow{{l1Phase, 1}, {l1Aux, -2}},
		ZhangCoordinateRow{{l2Phase, 1}, {l2Aux, -1}}};
	const std::array<ZhangCoordinateRow, 2> reference{};
	auto searched = subject;
	const auto accepted = zhangCheckProductOverlapNumericRows(
		subject, reference, searched, {true, true}, true);
	BOOST_REQUIRE_MESSAGE(accepted.valid, accepted.reason);
	const auto noL2 = zhangCheckProductOverlapNumericRows(
		subject, reference, searched, {true, false}, true);
	BOOST_CHECK(!noL2.valid);
	BOOST_CHECK_EQUAL(noL2.reason, "L2_CURRENT_ADMITTED_ARC_REQUIRED");
	searched[0][l1Phase] += 0.5;
	const auto wrongFunctional = zhangCheckProductOverlapNumericRows(
		subject, reference, searched, {true, true}, true);
	BOOST_CHECK(!wrongFunctional.valid);
	BOOST_CHECK_EQUAL(wrongFunctional.reason,
		"L1_PRODUCT_BRIDGE_NUMERIC_ROW_MISMATCH");
	const auto noAux = zhangCheckProductOverlapNumericRows(
		subject, reference, subject, {true, true}, false);
	BOOST_CHECK(!noAux.valid);
	BOOST_CHECK_EQUAL(noAux.reason, "RETIRED_ARC_AUXILIARY_REQUIRED");
}

BOOST_AUTO_TEST_CASE(r51_all_physical_arcs_retired_restarts_with_fresh_priors)
{
	const SatSys g01("G01"), g02("G02");
	const ZhangGraphEdge a1{"R0", g01}, a2{"R0", g02};
	const ZhangGraphEdge b1{"R1", g01}, b2{"R1", g02};
	const std::set<ZhangGraphEdge> all{a1, a2, b1, b2};
	const auto oldStateBasis = zhangBuildSpanningTree(
		all, "R0", {a1, a2, b1});
	const auto oldProductTarget = zhangBuildSpanningTree(
		all, "R0", {a1, b1, b2});
	BOOST_REQUIRE(oldStateBasis.connected && oldProductTarget.connected);
	const auto oldState = coordinateFixture(oldStateBasis);
	const std::map<ZhangGraphEdge, int> oldVersions{
		{a1, 3}, {a2, 3}, {b1, 3}, {b2, 3}};
	const auto plan = zhangPlanGraphCoordinateTransport(
		oldState, E_Sys::GPS, oldStateBasis, oldStateBasis, all,
		{{E_ObsCode::L1C, 0.1902936728},
		 {E_ObsCode::L2W, 0.2442102134}},
		[](const auto&, auto) { return 100.; },
		[](const auto&, auto) { return 25.; },
		nullptr, true, &oldVersions, &oldProductTarget.treeEdges);
	BOOST_REQUIRE_MESSAGE(plan.valid, plan.failureReason);
	BOOST_CHECK(plan.freshEdges == all);
	BOOST_CHECK(plan.survivingPhysicalRows.empty());
	BOOST_CHECK(plan.retainedOldArcAuxiliaryKeys.contains(
		{E_ObsCode::L1C, b2}));
	const auto auxiliary = plan.retainedOldArcAuxiliaryKeys.at(
		{E_ObsCode::L1C, b2});
	BOOST_CHECK(plan.transform.at(auxiliary).contains(
		zhangTransportAmbiguityKey(E_ObsCode::L1C, b2)));
	const auto next = evaluatePlan(oldState, plan, oldStateBasis);
	BOOST_CHECK(next.kfIndexMap.contains(auxiliary));
	BOOST_CHECK(plan.independentSourceVariances.size() > 0);
	// These observations are new versions of the same named edges; no old
	// versioned physical row is a surviving basis for an integer certificate.
	for (const auto& edge : all)
		BOOST_CHECK(!plan.survivingPhysicalRows.contains(
			{E_ObsCode::L1C, edge}));
}
