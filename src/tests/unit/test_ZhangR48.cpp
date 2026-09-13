#include <boost/test/unit_test.hpp>
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
 BOOST_REQUIRE(p.identityTransport);BOOST_CHECK(commitProductDatum(p).treeEdges==old.treeEdges);
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
 BOOST_REQUIRE(f.valid);BOOST_CHECK(f.affine.quotientProjector.empty());
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
