#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/receiver.hpp"
#include "common/zhangUserEpochStatus.hpp"
#include "pea/ppp.hpp"
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

extern "C" int __wrap_main(int, char**)
{
    unsigned checks = 0;
    const auto check = [&](bool pass, const char* message)
    {
        ++checks;
        if (!pass) throw std::runtime_error(message);
    };
    try
    {
        KFState state;
        state.time.bigTime = 100000;
        state.kfIndexMap.clear(); state.stateTransitionMap.clear();
        std::vector<KFKey> keys(5);
        keys[0].type=KF::ONE;
        for(int i=1;i<5;++i)
        {keys[i].type=i<4?KF::REC_POS:KF::TROP; keys[i].str="TEST";keys[i].num=i<4?i-1:0;}
        for(int i=0;i<5;++i)
        {state.kfIndexMap[keys[i]]=i;state.stateTransitionMap[keys[i]][keys[i]][0]=1;}
        state.x=VectorXd(5);state.x<<1,4200000,190000,4770000,.12;
        state.P=MatrixXd::Zero(5,5);
        MatrixXd a(4,4);a<<2,.2,.3,.1,0,3,.4,.2,0,0,4,.3,0,0,0,.5;
        state.P.bottomRightCorner(4,4)=a*a.transpose();
        state.lsqRequired=false;
        check(zhangUserMayContinueEmptyEpoch(true,0,state),"finite empty user epoch can continue");
        check(!zhangUserMayContinueEmptyEpoch(false,0,state),"server empty-epoch policy unchanged");
        check(!zhangUserMayContinueEmptyEpoch(true,1,state),"nonempty update cannot be skipped");
        auto invalid=state;invalid.P(1,2)=std::numeric_limits<double>::quiet_NaN();
        check(!zhangUserMayContinueEmptyEpoch(true,0,invalid),"nonfinite covariance cannot continue");
        const auto beforeX=state.x;const auto beforeP=state.P;
        check(!zhangUserHasMeasurementHistory(state),"no fake PPP initialisation");
        state.lastPppTransactionCommitted=true;state.lastPppQcConverged=true;
        zhangUserMarkEmptyEpoch(state);
        check(zhangUserHasMeasurementHistory(state),"legacy committed state retains initialisation");
        check(!zhangUserMayRunEpochAr(state),"previous successful epoch cannot authorize gap AR");
        check(state.x==beforeX && state.P==beforeP,"status never changes mean or full cross covariance");
        state.metaDataMap.erase(ZHANG_USER_LAST_UPDATE);
        acsConfig.process_ppp=true;acsConfig.zhangFullRank.enable=false;
        acsConfig.zhangPppAr.user_adapter=true;acsConfig.zhangPppAr.output_products=false;
        acsConfig.zhangPppAr.product_mode="HOU_OSB_LIKE";
        acsConfig.pppOpts.filter_reset_enable=false;
        acsConfig.pppOpts.merge_correlated_states=false;
        acsConfig.pppOpts.ionoOpts.use_if_combo=false;
        acsConfig.phaseClockOsb.enable=false;
        tsync.bigTime=100030;
        ReceiverMap receivers;KFState remote;std::ostringstream trace;
        check(ppp(trace,receivers,state,remote),"production empty PPP continues after valid prediction");
        check(state.x.isApprox(beforeX,1e-14) && state.P.isApprox(beforeP,1e-14),"production skip preserves full identity prediction and cross blocks");
        check(!state.lastPppTransactionCommitted && !state.lastPppQcConverged,"empty epoch is not a committed observation");
        check(state.metaDataMap.at(ZHANG_USER_EPOCH_STATUS)=="UNINITIALIZED","initial outage remains uninitialised");
        check(trace.str().find("boundary=B16")==std::string::npos,"empty epoch creates no measurement factor");
        state.metaDataMap[ZHANG_USER_LAST_UPDATE]="prior valid epoch";
        state.lastPppTransactionCommitted=true;state.lastPppQcConverged=true;
        tsync.bigTime=100060;
        check(ppp(trace,receivers,state,remote),"midstream gap continues");
        check(state.metaDataMap.at(ZHANG_USER_EPOCH_STATUS)=="PREDICTION_ONLY" && !zhangUserMayRunEpochAr(state),"midstream gap blocks stale AR without resetting state");
        acsConfig.zhangPppAr.user_adapter=false;
        check(!ppp(trace,receivers,state,remote),"service empty authoritative update still fails closed");
        std::cout<<"PASS user_empty_epoch_contract checks="<<checks<<'\n';return 0;
    }
    catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
