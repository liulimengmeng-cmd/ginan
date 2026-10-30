#include "common/zhangRegionalTropModel.hpp"
#include "common/zhangUserGonsisOutput.hpp"
#include <iostream>
#include <stdexcept>

void check(bool ok,const char* reason) { if(!ok) throw std::runtime_error(reason); }
int main()
{
    using namespace Eigen;
    MatrixXd a(20,4);
    for(int i=0;i<20;++i) a.row(i)<<1,std::cos(i*1.17),std::sin(i*1.17),0.05*(i%4)-0.08;
    Vector4d known;known<<0.04,0.006,-0.004,0.02;
    MatrixXd root=MatrixXd::Identity(20,20)*0.002;
    root.col(0).array()+=0.003;
    const MatrixXd q=root*root.transpose();
    auto f=zhangtrop::fit(a,a*known,q);
    check(f.valid && f.order==4,"full-rank fit");
    check((f.coefficients-known).norm()<1e-12 && f.looRms<1e-12,"coefficients and unbiased LOO closure");
    check(f.spatialVariance==0.001,"fixed extra user prediction variance in m^2");
    // Independent direct normal-equation oracle on a well-conditioned fixture.
    MatrixXd oracle=(a.transpose()*a).ldlt().solve(a.transpose());
    check((oracle-f.jacobian).norm()<1e-11,"SVD projection vs independent oracle");
    for(int omitted=0;omitted<20;++omitted) {
        MatrixXd reduced(19,4);for(int i=0,k=0;i<20;++i)if(i!=omitted)reduced.row(k++)=a.row(i);
        RowVectorXd w=a.row(omitted)*(reduced.transpose()*reduced).ldlt().solve(reduced.transpose());
        VectorXd contrast=VectorXd::Zero(20);contrast(omitted)=1;
        for(int i=0,k=0;i<20;++i)if(i!=omitted)contrast(i)=-w(k++);
        check(std::abs(contrast.dot(q*contrast)-f.leaveOneOutFormalVariances(omitted))<1e-14,"LOO includes negative source cross term");
    }
    const auto common=zhangtrop::fit(a,(a*known).array()+0.05,q);
    check(common.valid && common.looRms<1e-12 && std::abs(common.coefficients(0)-known(0)-0.05)<1e-12,"common bias is invisible to LOO; no absolute accuracy claim");
    MatrixXd flat=a;flat.col(3).setZero();
    check(zhangtrop::fit(flat,flat*known,q).order==3,"equal-height geometry downgrades to EN");
    flat.col(2)=flat.col(1);check(zhangtrop::fit(flat,VectorXd::Constant(20,.04),q).order==1,"collinear geometry downgrades to constant");
    VectorXd bad=a*known;bad(2)+=1;
    check(!zhangtrop::fit(a,bad,q).valid,"gross LOO mismatch rejected");
    MatrixXd corrupt=q;corrupt(0,0)=-1;
    check(!zhangtrop::fit(a,a*known,corrupt).valid,"indefinite source rejected");
    // Same full native P, correlated clock/phase/ZTD. No zero-filled cross block.
    MatrixXd sourceRoot=MatrixXd::Random(23,23)*0.01;
    MatrixXd sourceQ=sourceRoot*sourceRoot.transpose();
    MatrixXd projection=MatrixXd::Zero(7,23);projection.topLeftCorner(3,3).setIdentity();projection.bottomRightCorner(4,20)=f.jacobian;
    MatrixXd projected=projection*sourceQ*projection.transpose();
    check(projected.topRightCorner(3,4).norm()>1e-6,"nonzero CP-TROP covariance propagated");
    check((projected.topLeftCorner(3,3)-sourceQ.topLeftCorner(3,3)).norm()<1e-15,"OFF marginal is unchanged, not Schur conditioned");
    check(SelfAdjointEigenSolver<MatrixXd>(projected).eigenvalues().minCoeff()>-1e-13,"projected Q PSD");
    // Independent AR output check includes TROP with the existing exact 4-sigma rule.
    using K=ZhangUserGonsisMeasurementKind;
    VectorXd x=VectorXd::Zero(2), v=VectorXd::Zero(3), vars=VectorXd::Constant(3,.001);
    MatrixXd h=MatrixXd::Zero(3,2), r=vars.asDiagonal(); h(2,1)=1;r(0,2)=r(2,0)=.0001;
    std::vector<K> kinds{K::CODE,K::PHASE,K::REGIONAL_TROP_ZTD};
    auto gate=zhangValidateUserGonsisResiduals(x,x,h,v,vars,kinds,&r,0,1);
    check(gate.valid && gate.regionalTropRows==1,"TROP AR residual captured with joint R");
    v(2)=4*std::sqrt(.001)+1e-8;
    gate=zhangValidateUserGonsisResiduals(x,x,h,v,vars,kinds,&r,0,1);
    check(!gate.valid && gate.reason=="GONSIS_FIXED_REGIONAL_TROP_4SIGMA_FAILED","TROP residual over 4 sigma rejected");
    check(!zhangValidateUserGonsisResiduals(x,x,h,v,vars,kinds,nullptr,0,1).valid,"TROP requires full R");
    check(!zhangValidateUserGonsisResiduals(x,x,h,v,vars,kinds,&r,0,0).valid,"dropped TROP row detected");
    std::cout<<"TROP_MODEL_CONTRACT_PASS geometry_4_3_1=1 full_cross_Q=1 loo_common_bias_limit=1 exact_4sigma=1\n";
}
