#include "common/zhangServerNumericalSafety.hpp"
#include "common/zhangEffectiveClockDatum.hpp"
#include "common/zhangFloatThreading.hpp"
#include "common/observations.hpp"
#include <boost/json.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>

bool reacquireZhangReceiverClocks(Trace&,ReceiverMap&,KFState&);
static void require(bool b,const char* reason){if(!b)throw std::runtime_error(reason);}

static void resetContract()
{
    MatrixXd f,q;
    require(zhangBuildIndependentStateReset(4,0,{{1,123.,1e6}},f,q),"reset rejected");
    VectorXd x(4); x<<1,-80137,7,9;
    MatrixXd p=MatrixXd::Zero(4,4);
    p.bottomRightCorner(3,3)<<4,1,2,1,9,3,2,3,16;
    const VectorXd mean=f*x;
    const MatrixXd cov=f*p*f.transpose()+q;
    // Independent distribution oracle: retain the OTHER marginal, replace
    // exactly one variable by an independent N(123,1e6). Do not call reset again.
    require(mean[1]==123 && mean[2]==7 && mean[3]==9,"wrong marginal mean");
    require((cov.bottomRightCorner(2,2)-p.bottomRightCorner(2,2)).norm()==0,"other marginal changed");
    require(cov(1,1)==1e6 && cov(1,2)==0 && cov(1,3)==0,"old cross information inherited");
    require(zhangCheckFullCovariance(cov).valid,"reset not PSD");
    require(!zhangBuildIndependentStateReset(4,0,{{0,2,1}},f,q),"ONE reset accepted");
    require(!zhangBuildIndependentStateReset(4,0,{{1,2,1},{1,3,1}},f,q),"duplicate reset accepted");
    require(!zhangBuildIndependentStateReset(4,0,{{1,2,0}},f,q),"zero prior accepted");
    MatrixXd bad=MatrixXd::Identity(3,3);bad(0,1)=bad(1,0)=2;
    require(!zhangCheckFullCovariance(bad).valid,"positive diagonal hid indefinite matrix");
}

static void actualReacquisition()
{
    acsConfig.zhangFullRank.enable=true; acsConfig.zhangPppAr.user_adapter=false;
    acsConfig.process_sys.clear();acsConfig.process_sys[E_Sys::GPS]=true;
    acsConfig.zhangFullRank.sysOpts[E_Sys::GPS].baseline_observables={E_ObsCode::L1C,E_ObsCode::L2W};
    acsConfig.epoch_interval=30;
    ReceiverMap receivers;auto& rec=receivers["VITY"];rec.id="VITY";rec.ready=true;
    rec.obsList.push_back(std::make_shared<GObs>());
    rec.sol.status=E_Solution::SINGLE;
    auto& options=acsConfig.getRecOpts("VITY");options.clk.sigma={1000};options.clk.estimate={true};
    KFState state;tsync.bigTime=100000;state.time=tsync;
    const KFKey clock=zhangReceiverClockKey(acsConfig,"VITY",E_Sys::GPS,0,KF::REC_CLOCK,&rec);
    KFKey other;other.type=KF::TROP;other.str="OTHER";
    state.kfIndexMap.clear();state.kfIndexMap[KFState::oneKey]=0;
    state.kfIndexMap[clock]=1;state.kfIndexMap[other]=2;
    state.x=VectorXd::Zero(3);state.x<<1,-80137,3;
    state.P=MatrixXd::Zero(3,3);state.P(1,1)=483000;state.P(2,2)=4;state.P(1,2)=state.P(2,1)=1;
    state.dx=VectorXd::Zero(3);
    state.addKFState(clock,{.x=-80137,.P=483000});state.addKFState(other,{.x=3,.P=4});
    const auto offsetKey=zhangClockLifecycleMetadataKey("SPP_OFFSET",E_Sys::GPS,"VITY");
    auto receipt=[&](double gap){GTime previous=tsync;previous.bigTime-=gap;
        state.metaDataMap[offsetKey]=boost::json::serialize(boost::json::object{
            {"epoch",previous.to_string(9)},{"chart","CONTROLLED_PHYSICAL_CHART"},{"offset_m",-5.}});};
    // Controlled SPP/seed evidence is a unit fixture, not a real-data certificate.
    KFKey spp{KF::REC_SYS_BIAS,SatSys(E_Sys::GPS,0),"VITY"};
    rec.sol.sppState.kfIndexMap[spp]=1;
    rec.sol.sppState.x=VectorXd::Zero(2);rec.sol.sppState.x(1)=134361.;
    rec.sol.sppState.P=MatrixXd::Identity(2,2);
    rec.sol.sppState.metaDataMap[zhangClockLifecycleMetadataKey("SPP_CHART",E_Sys::GPS)]=
        boost::json::serialize(boost::json::object{{"epoch",tsync.to_string(9)},{"chart","CONTROLLED_PHYSICAL_CHART"}});
    state.metaDataMap[zhangClockLifecycleMetadataKey("SEED",E_Sys::GPS)]=
        boost::json::serialize(boost::json::object{{"epoch",tsync.to_string(9)},{"translation_m",-12.}});
    std::ostringstream trace;
    receipt(30);const auto before=state.P;
    require(reacquireZhangReceiverClocks(trace,receivers,state),"normal epoch failed");
    require((state.P-before).norm()==0 && state.x(1)==-80137,"continuous clock reset");
    receipt(4830);
    require(reacquireZhangReceiverClocks(trace,receivers,state),"actual reacquisition failed");
    require(state.x(1)==134349 && state.x(2)==3 && state.P(2,2)==4,"actual reset mean/marginal wrong");
    require(state.P(1,1)==1e6 && state.P(1,2)==0,"actual reset covariance wrong");
    require(trace.str().find("cross_information=RELEASED")!=std::string::npos,"reset evidence missing");
    require(state.metaDataMap[offsetKey].find("CONTROLLED_PHYSICAL_CHART")!=std::string::npos,
        "offset evidence unexpectedly removed");
}

static void actualKalman()
{
    KFState state;state.x=VectorXd::Zero(3);state.x(0)=1;state.P=MatrixXd::Zero(3,3);
    state.P(1,1)=4.6008e10;state.P(2,2)=.25;
    state.joseph_stabilisation=false;state.postfitOpts.omega_test=true;
    KFMeas meas;meas.H=MatrixXd::Zero(2,3);meas.H.block(0,1,2,2)<<1,1,1,-1;
    meas.R=MatrixXd::Identity(2,2)*4e-6;meas.V=VectorXd::Ones(2)*10;
    meas.H_star=MatrixXd::Zero(2,0);meas.uncorrelatedNoise=VectorXd::Zero(0);
    VectorXd xp=state.x,dx=VectorXd::Zero(3);MatrixXd pp=state.P,inv,invh;std::ostringstream trace;
    ZhangFloatThreading::inFloat=true;
    require(state.kFilter(trace,meas,xp,pp,dx,inv,invh,0,3,0,2,false),"actual Joseph solve failed");
    require(zhangCheckFullCovariance(pp).valid,"actual posterior not PSD");
    // Independent scalar information oracle for orthogonal sum/difference rows.
    const long double clock=1/(1/4.6008e10L+2/4e-6L);
    const long double other=1/(1/.25L+2/4e-6L);
    require(std::abs(pp(1,1)-clock)<2e-8 && std::abs(pp(2,2)-other)<2e-8,
        "Joseph covariance disagrees with independent information oracle");
    state.P(1,1)=state.P(2,2)=1;state.P(1,2)=state.P(2,1)=2;
    require(!state.kFilter(trace,meas,xp,pp,dx,inv,invh,0,3,0,2,false),"actual kernel accepted indefinite prior");
    ZhangFloatThreading::inFloat=false;
}
int main()
{
    try {setenv("ZHANG_SERVER_REACQ_FIX","1",1);resetContract();actualReacquisition();actualKalman();
        std::cout<<"SERVER_REACQ_NUMERICAL_CONTRACT_PASS actual_production_kernel=1 independent_marginal_oracle=1 independent_information_oracle=1\n";return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
