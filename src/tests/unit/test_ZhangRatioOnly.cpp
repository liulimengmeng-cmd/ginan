#include <boost/test/unit_test.hpp>
#include "common/zhangRatioOnly.hpp"
#include "common/zhangIntegerCandidateNis.hpp"
#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangR47History.hpp"
#include <cstdlib>
#include "common/zhangProductIntegerLedger.hpp"
#include "common/zhangProductGaugeCertificateLedger.hpp"
#include "common/zhangProductRelationSolver.hpp"
BOOST_AUTO_TEST_CASE(r51_ratio_policy_retains_structural_rejections)
{
    unsetenv("ZHANG_R51_RATIO_ONLY");
    BOOST_CHECK(!zhangRatioStatisticalAccept(false));
    setenv("ZHANG_R51_RATIO_ONLY","1",1);
    BOOST_CHECK(zhangRatioStatisticalAccept(false));
    BOOST_CHECK(!zhangRatioStatisticalReject(true));
    setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS","1",1);
    BOOST_CHECK(!zhangExpensiveDiagnosticsEnabled());
    unsetenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS");
    Eigen::VectorXd v=Eigen::VectorXd::Constant(1,10);
    Eigen::MatrixXd q=Eigen::MatrixXd::Identity(1,1);
    int nisCalls=0;
    const auto skipped=zhangRatioOnlySkipNis([&]() {
        ++nisCalls;return assessZhangIntegerCandidateNis(v,q,0.01);
    });
    BOOST_CHECK_EQUAL(nisCalls,0);
    BOOST_CHECK(!skipped.valid);
    BOOST_CHECK(std::isnan(skipped.nis));
    BOOST_CHECK_EQUAL(skipped.status,"SKIPPED_RATIO_ONLY");
    // Skipping the statistical solve must not turn a deterministic conflict
    // into an admissible posterior update.
    const auto conflict=zhangConditionPosteriorEffectiveIntegers(
        Eigen::VectorXd::Zero(1),Eigen::MatrixXd::Zero(1,1),
        Eigen::MatrixXd::Identity(1,1),Eigen::VectorXd::Ones(1));
    BOOST_CHECK(!conflict.valid);
    auto high=assessZhangIntegerCandidateNis(v,q,0.01);
    BOOST_CHECK(high.valid);
    BOOST_CHECK(std::isnan(high.nis));
    BOOST_CHECK(std::isnan(high.threshold));
    BOOST_CHECK_EQUAL(high.status,"NIS_SKIPPED_RATIO_ONLY_GEOMETRY_VALID");
    BOOST_CHECK(zhangRatioStatisticalAccept(false));
    q(0,0)=-1;BOOST_CHECK(!assessZhangIntegerCandidateNis(v,q,0.01).valid);
    q(0,0)=0;BOOST_CHECK(!assessZhangIntegerCandidateNis(v,q,0.01).valid);
    v(0)=0;
    const auto deterministic=assessZhangIntegerCandidateNis(v,q,0.01);
    BOOST_CHECK(deterministic.valid);
    BOOST_CHECK(deterministic.deterministic);
    BOOST_CHECK(std::isnan(deterministic.nis));
    v(0)=10;
    Eigen::MatrixXd semidefinite=Eigen::MatrixXd::Zero(2,2);
    semidefinite(0,0)=1;
    Eigen::VectorXd allowed(2);allowed<<0.5,0;
    const auto validNull=assessZhangIntegerCandidateNis(
        allowed,semidefinite,0.01);
    BOOST_CHECK(validNull.valid);
    BOOST_CHECK(std::isnan(validNull.nis));
    allowed(1)=1;
    BOOST_CHECK(!assessZhangIntegerCandidateNis(
        allowed,semidefinite,0.01).valid);
    Eigen::VectorXd one=Eigen::VectorXd::Ones(1);
    const auto conditioned=zhangConditionExactProductRows(
        Eigen::VectorXd::Constant(1,0.2),
        Eigen::MatrixXd::Identity(1,1),
        Eigen::MatrixXd::Identity(1,1),one);
    BOOST_CHECK(conditioned.valid);
    BOOST_CHECK(std::isnan(conditioned.nis));
    Eigen::VectorXd measurements(2);measurements<<1,1.1;
    const auto bridge=zhangComponentBridgeGls(
        measurements,Eigen::MatrixXd::Identity(2,2));
    BOOST_CHECK(bridge.valid);
    BOOST_CHECK(std::isnan(bridge.residualNis));
    ProductGaugeCertificate certificate;
    certificate.wideLaneReliable=true;
    certificate.firstSignalReliable=true;
    certificate.exactProductLatticeMembership=true;
    certificate.cycleClosurePassed=true;
    certificate.temporalAlignmentCertified=true;
    certificate.jointNisPassed=false;
    BOOST_CHECK(zhangProductGaugeCertificateEvidenceComplete(certificate,1));
    BOOST_CHECK(!zhangR47AffineIntegerFeasible({{1},{1}},{0,1},1));
    const auto recheck=zhangRecheckProductIntegerOnPosterior(0.1,1e-5,0,1e-3,0.01);
    BOOST_CHECK(recheck.reliable);
    BOOST_CHECK(std::isnan(recheck.nis));
    BOOST_CHECK(std::isnan(recheck.nisThreshold));
    auto proof=std::make_shared<ZhangIntegerDecisionProof>();
    proof->id="R51_RATIO_TEST";proof->originalStatement="x=0";proof->observationProvenance="TEST";proof->conditionalFailureBound=1;
    auto selected=zhangR47SelectHistorySubset({{1}},{0},{{proof}},{},{1},1,1e-3,2.5e-4);
    BOOST_CHECK(selected.valid);BOOST_CHECK_EQUAL(selected.selected.size(),1);BOOST_CHECK_EQUAL(selected.risk,1);
    auto missing=zhangR47SelectHistorySubset({{1}},{0},{{}},{},{1},1,1e-3,2.5e-4);
    BOOST_CHECK(missing.selected.empty());
    unsetenv("ZHANG_R51_RATIO_ONLY");
    BOOST_CHECK(!zhangProductGaugeCertificateEvidenceComplete(certificate,1));
    q(0,0)=1;
    const auto evaluated=zhangRatioOnlySkipNis([&]() {
        ++nisCalls;return assessZhangIntegerCandidateNis(v,q,0.01);
    });
    BOOST_CHECK_EQUAL(nisCalls,1);
    BOOST_CHECK(evaluated.valid);
    BOOST_CHECK(evaluated.nis>evaluated.threshold);
    auto baseline=zhangR47SelectHistorySubset({{1}},{0},{{proof}},{},{1},1,1e-3,2.5e-4);
    BOOST_CHECK(baseline.selected.empty());
    BOOST_CHECK(!zhangRecheckProductIntegerOnPosterior(0.1,1e-5,0,1e-3,0.01).reliable);
}
