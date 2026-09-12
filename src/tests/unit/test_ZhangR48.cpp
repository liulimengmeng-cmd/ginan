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
