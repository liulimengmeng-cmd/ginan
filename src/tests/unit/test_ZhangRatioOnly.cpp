#include <boost/test/unit_test.hpp>
#include "common/zhangRatioOnly.hpp"
#include "common/zhangIntegerCandidateNis.hpp"
#include "common/zhangR47History.hpp"
#include <cstdlib>
#include "common/zhangProductIntegerLedger.hpp"
BOOST_AUTO_TEST_CASE(r51_ratio_policy_retains_structural_rejections)
{
    unsetenv("ZHANG_R51_RATIO_ONLY");
    BOOST_CHECK(!zhangRatioStatisticalAccept(false));
    setenv("ZHANG_R51_RATIO_ONLY","1",1);
    BOOST_CHECK(zhangRatioStatisticalAccept(false));
    BOOST_CHECK(!zhangRatioStatisticalReject(true));
    Eigen::VectorXd v=Eigen::VectorXd::Constant(1,10);
    Eigen::MatrixXd q=Eigen::MatrixXd::Identity(1,1);
    auto high=assessZhangIntegerCandidateNis(v,q,0.01);
    BOOST_CHECK(high.valid);BOOST_CHECK(high.nis>high.threshold);
    BOOST_CHECK(zhangRatioStatisticalAccept(high.nis<=high.threshold));
    q(0,0)=-1;BOOST_CHECK(!assessZhangIntegerCandidateNis(v,q,0.01).valid);
    q(0,0)=0;BOOST_CHECK(!assessZhangIntegerCandidateNis(v,q,0.01).valid);
    BOOST_CHECK(!zhangR47AffineIntegerFeasible({{1},{1}},{0,1},1));
    const auto recheck=zhangRecheckProductIntegerOnPosterior(0.1,1e-5,0,1e-3,0.01);
    BOOST_CHECK(recheck.reliable);BOOST_CHECK(recheck.nis>recheck.nisThreshold);
    auto proof=std::make_shared<ZhangIntegerDecisionProof>();
    proof->id="R51_RATIO_TEST";proof->originalStatement="x=0";proof->observationProvenance="TEST";proof->conditionalFailureBound=1;
    auto selected=zhangR47SelectHistorySubset({{1}},{0},{{proof}},{},{1},1,1e-3,2.5e-4);
    BOOST_CHECK(selected.valid);BOOST_CHECK_EQUAL(selected.selected.size(),1);BOOST_CHECK_EQUAL(selected.risk,1);
    auto missing=zhangR47SelectHistorySubset({{1}},{0},{{}},{},{1},1,1e-3,2.5e-4);
    BOOST_CHECK(missing.selected.empty());
    unsetenv("ZHANG_R51_RATIO_ONLY");
    auto baseline=zhangR47SelectHistorySubset({{1}},{0},{{proof}},{},{1},1,1e-3,2.5e-4);
    BOOST_CHECK(baseline.selected.empty());
    BOOST_CHECK(!zhangRecheckProductIntegerOnPosterior(0.1,1e-5,0,1e-3,0.01).reliable);
}
