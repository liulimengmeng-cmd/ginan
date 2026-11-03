#include "common/zhangUserGonsisOutput.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

static void require(bool condition,const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
static bool close(double actual,double expected)
{
    return std::abs(actual-expected)<=1e-13;
}

int main()
{
    using Kind=ZhangUserGonsisMeasurementKind;
    Eigen::VectorXd x=Eigen::VectorXd::Zero(2),fixed(2);
    fixed<<0.1,-0.1;
    Eigen::MatrixXd h(3,2); h<<1,0,0,1,1,-1;
    Eigen::MatrixXd r(3,3);
    r<<0.01,0.0002,0.0001,0.0002,0.0004,0.00005,0.0001,0.00005,0.0004;
    const Eigen::VectorXd variance=r.diagonal();
    Eigen::VectorXd residual(3); residual<<0.11,-0.09,0.1;
    const std::vector<Kind> kinds{Kind::CODE,Kind::PHASE,Kind::REGIONAL_IONO_SD};

    // Independent literal oracle: H*(fixed-x)=[.1,-.1,.2]. The omitted
    // pseudo's fixed residual is -.1m=5 sigma while both physical rows pass.
    const auto oldPhysicalOnly=zhangValidateUserGonsisResiduals(x,fixed,h.topRows(2),
        residual.head(2),variance.head(2),{Kind::CODE,Kind::PHASE});
    require(oldPhysicalOnly.valid,"fixture must pass the old physical-only gate");
    const auto fail=zhangValidateUserGonsisResiduals(x,fixed,h,residual,variance,kinds,&r,1);
    require(!fail.valid && fail.reason=="GONSIS_FIXED_REGIONAL_IONO_4SIGMA_FAILED",
        "post-FIX regional inconsistency was not rejected");
    require(fail.codeRows==1 && fail.phaseRows==1 && fail.regionalIonoRows==1 &&
        fail.regionalIonoSummary.rowCount==1,"pseudo mixed into physical row counts");
    require(fail.fullMeasurementCovarianceValidated && fail.firstRejectedRow==2,
        "complete same-epoch R was not bound or rejection row was lost");
    require(close(fail.fixedResidual(2),-0.1) && close(fail.measurementSigmas(2),0.02) &&
        close(fail.normalizedResiduals(2),5) && close(fail.fourSigmaMarginsMetres(2),-0.02),
        "literal residual/sigma/margin oracle disagrees");
    require(close(fail.regionalIonoSummary.maximumNormalizedResidual,5) &&
        close(fail.regionalIonoSummary.minimumFourSigmaMarginMetres,-0.02),
        "regional summary does not expose the actual failing margin");
    std::ostringstream trace;
    zhangTraceUserGonsisResidualGate(trace,"time=EXAMPLE receiver=TEST attempt=COMBINED",fail);
    require(trace.str().find("kind=REGIONAL_IONO_SD rows=1")!=std::string::npos &&
        trace.str().find("quadratic_forms=0")!=std::string::npos &&
        trace.str().find("nis_evaluations=0")!=std::string::npos,
        "per-kind numeric/no-NIS audit receipt missing");

    // Do not substitute HPH' or whitened/marginal approximations for source R.
    require(!zhangValidateUserGonsisResiduals(x,fixed,h,residual,variance,kinds,nullptr,1).valid,
        "pseudo gate accepted a discarded full R");
    Eigen::MatrixXd wrongR=r; wrongR(0,2)=wrongR(2,0)=1;
    require(!zhangValidateUserGonsisResiduals(x,fixed,h,residual,variance,kinds,&wrongR,1).valid,
        "invalid CP-ionosphere cross block was ignored");
    wrongR=r; wrongR(2,2)*=2;
    require(!zhangValidateUserGonsisResiduals(x,fixed,h,residual,variance,kinds,&wrongR,1).valid,
        "mismatched covariance diagonal and sigma source were accepted");
    require(!zhangValidateUserGonsisResiduals(x,fixed,h,residual,variance,kinds,&r,2).valid,
        "regional builder/capture row count mismatch was ignored");
    const Eigen::MatrixXd cpR=r.topLeftCorner(2,2);
    require(!zhangValidateUserGonsisResiduals(x,fixed,h.topRows(2),residual.head(2),
        variance.head(2),{Kind::CODE,Kind::PHASE},&cpR,1).valid,
        "dropping a declared pseudo row silently became valid zero-pseudo geometry");
    require(zhangValidateUserGonsisResiduals(x,fixed,h.topRows(2),residual.head(2),
        variance.head(2),{Kind::CODE,Kind::PHASE},&cpR,0).valid,
        "actual final measurement with no pseudo rows was unconditionally rejected");

    // Four-sigma boundary uses the same inclusive comparison as CODE/PHASE.
    Eigen::VectorXd zero=Eigen::VectorXd::Zero(2);
    Eigen::VectorXd boundary(3); boundary<<0,0,4*std::sqrt(variance(2));
    require(zhangValidateUserGonsisResiduals(zero,zero,h,boundary,variance,kinds,&r,1).valid,
        "exactly four regional sigma must be accepted");
    boundary(2)=std::nextafter(boundary(2),1.0);
    require(!zhangValidateUserGonsisResiduals(zero,zero,h,boundary,variance,kinds,&r,1).valid,
        "above four regional sigma was accepted");

    bool row=false; std::string reason;
    require(zhangUserGonsisRegionalResidualRowContract(true,true,true,true,row,reason) && row,
        "actual independent builder marker and metre label did not bind");
    require(!zhangUserGonsisRegionalResidualRowContract(true,true,true,false,row,reason),
        "mistyped metre label silently dropped a constructed pseudo");
    require(!zhangUserGonsisRegionalResidualRowContract(true,true,false,true,row,reason),
        "missing builder marker silently dropped a labelled pseudo");
    require(!zhangUserGonsisRegionalResidualRowContract(true,true,false,false,row,reason),
        "regional IONOSPHERIC row was treated as zero because both tags disappeared");
    require(!zhangUserGonsisRegionalResidualRowContract(true,false,true,true,row,reason),
        "regional pseudo with a CODE/PHASE type was miscounted");
    require(zhangUserGonsisRegionalResidualRowContract(true,false,false,false,row,reason) && !row,
        "normal physical row was mistaken for regional pseudo");
    std::cout<<"USER_REGIONAL_GONSIS_RESIDUAL_CONTRACT_PASS\n";
}
