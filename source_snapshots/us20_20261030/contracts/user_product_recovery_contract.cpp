#include "user_chart_fixture.hpp"
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserProductChart.hpp"
#include "pea/zhangPppAr.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Real production reference/prediction/design/expiry hooks. Independent
// marginal/S-transform oracle, with a >900s lifecycle; no GNSS quality claim.
extern "C" int __wrap_main(int,char**)
{
    unsigned checks=0;
    auto check=[&](bool yes,const char* text)
    { ++checks; if(!yes) throw std::runtime_error(text); };
    try
    {
        unsetenv("ZHANG_USER_UNKNOWN_CHART");
        check(zhangUserUnknownChartPolicy()==ZhangUserUnknownChartPolicy::RECOVER,"safe recovery is default");
        setenv("ZHANG_USER_UNKNOWN_CHART","RECOVER",1);
        setenv("ZHANG_USER_QC","CAUSAL_TRACKING",1);
        setenv("ZHANG_USER_FLOAT_TRANSPORT","CONTRACT",1);
        setenv("ZHANG_USER_AR_SEARCH","GONSIS_UDUC",1);
        acsConfig.process_ppp=true;
        acsConfig.zhangPppAr.user_adapter=true;
        acsConfig.zhangPppAr.output_products=false;
        acsConfig.zhangPppAr.product_solution="PRODUCT_FIXED";
        acsConfig.zhangPppAr.integer_strategy="CANONICAL_USER_SD_WL_L1";
        acsConfig.zhangPppAr.user_use_full_product_covariance=false;
        acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow=false;
        acsConfig.zhangPppAr.targeted_besd_capture_shadow=false;
        acsConfig.zhangFullRank.enable=false;
        acsConfig.phaseClockOsb.enable=false;
        acsConfig.pppOpts.ionoOpts.use_if_combo=false;
        int serial=0;
        for (auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
        for (int mode=0;mode<6;++mode)
        {
            // 0 local unknown; 1 unknown reference + continuous replacement;
            // 2 all unknown; 3 unavailable product; 4 true physical slip;
            // 5 genuine observation outage and later return, not product loss.
            ++serial;
            const auto fixture=std::filesystem::current_path()/
                ("user_recovery_"+std::to_string(serial)+".csv");
            std::ofstream file(fixture);
            file << "synthetic CSV52 long lifecycle fixture\n";
            for (int epoch=0;epoch<=65;++epoch)
            for (int sat=1;sat<=4;++sat)
            {
                if (mode==3 && epoch==1 && sat==3) continue;
                const bool changed=epoch>0 && ((mode==0 && sat==3) || (mode==1 && sat==1) || mode==2);
                std::vector<std::string> r(52,"0");
                r[0]=std::to_string(100000+epoch*30); r[1]="PRODUCT_FIXED";
                r[2]=SatSys(E_Sys::GPS,sat).id(); r[3]=enum_to_string(code);
                r[5]=".01"; r[7]=".01"; r[10]=".02";
                r[6]=std::to_string(-.1*sat-.13*epoch); r[9]=std::to_string(.1*sat+.13*epoch);
                r[11]=changed?"1":"0"; r[14]="0"; r[15]=r[0]; r[17]="SYNTHETIC_RECOVERY";
                r[18]="1"; r[19]="CURRENT_ALIGNMENT_VALID";
                for (int f : {20,21,22,23,28,29,30,31,32,42,43,46}) r[f]="1";
                r[24]="RECOVERY-C0"; r[25]="RECOVERY-D0"; r[26]=r[0]; r[27]=r[0]; r[33]="NONE";
                r[34]=r[2]+"-PHYSICAL0"; r[35]="1"; r[38]="4"; r[39]="3"; r[40]="3";
                r[44]="AR_VALID"; r[47]="SYNTHETIC_FUNCTION";
                r[49]=changed && epoch==1?"0":"1"; r[50]=changed && epoch==1?"1":"0";
                r[51]=r[2]+"-PUBLIC"+(changed?"1":"0");
                for (int i=0;i<r.size();++i) file<<(i?",":"")<<r[i];
                file<<'\n';
            }
            file.close();
            acsConfig.zhangPppAr.product_filename=fixture.string();
            acsConfig.zhangPppAr.baseline_observables.clear();
            acsConfig.zhangPppAr.baseline_observables[E_Sys::GPS]={code};
            ReceiverMap receivers;
            auto& rec=receivers["RECOVERTEST"]; rec.id="RECOVERTEST"; rec.ready=true;
            const auto freq=code2Freq.at(E_Sys::GPS).at(code);
            for (int sat=1;sat<=4;++sat)
            {
                SatSys satellite(E_Sys::GPS,sat);
                auto& status=rec.satStatMap[satellite]; status.el=(70-10*sat)*D2R;
                auto& sig=status.sigStatMap[ft2string(freq)];
                sig.tracking.seen=true; sig.tracking.lastCode=static_cast<int>(code); sig.tracking.lastPhaseTime=100000;
                GObs obs; obs.time.bigTime=100000; obs.Sat=satellite; obs.mount=rec.id; obs.satStat_ptr=&status;
                obs.sigs[freq].code=code; obs.sigs[freq].P=2e7; obs.sigs[freq].L=1e8;
                rec.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
            }
            auto ambiguity=[&](int sat) { KFKey k; k.type=KF::AMBIGUITY; k.str=rec.id;
                k.Sat=SatSys(E_Sys::GPS,sat); k.num=static_cast<int>(code); return k; };
            KFKey clock; clock.type=KF::REC_CLOCK; clock.str=rec.id;
            KFKey position; position.type=KF::REC_POS; position.str=rec.id;
            KFKey phase; phase.type=KF::PHASE_BIAS; phase.str=rec.id; phase.Sat=SatSys(E_Sys::GPS,0); phase.num=static_cast<int>(code);
            KFState state; state.time.bigTime=100000;
            check(bindZhangCheckpointRuntimeId(state,"recovery-test-"+std::to_string(serial)),"bind runtime");
            std::vector<KFKey> keys{state.oneKey,clock,position,phase,ambiguity(2),ambiguity(3),ambiguity(4)};
            state.kfIndexMap.clear(); state.stateTransitionMap.clear();
            for (int i=0;i<keys.size();++i)
            {
                state.kfIndexMap[keys[i]]=i; state.stateTransitionMap[keys[i]][keys[i]][0]=1;
                if (keys[i].type==KF::AMBIGUITY)
                { auto key=keys[i]; key.estimatedTime=state.time; state.outageLimitMap[key]=900; }
            }
            state.x.resize(7); state.x<<1,7,1.2,20,3.25,-2.5,4.3;
            state.dx=VectorXd::Zero(7);
            MatrixXd factor=MatrixXd::Identity(7,7);
            factor(1,3)=.23; factor(2,4)=.31; factor(3,5)=.17; factor(4,1)=-.28; factor(5,2)=.19; factor(6,3)=-.11;
            state.P=factor*factor.transpose(); state.P.row(0).setZero(); state.P.col(0).setZero();
            std::ostringstream trace;
            initializeFixtureChart(trace,receivers,state);
            check(zhangPppArUserReference(state,rec.id,E_Sys::GPS,code)==SatSys("G01"),"initial known anchor");
            for (int epoch=1;epoch<=65;++epoch)
            {
                state.time.bigTime=100000+epoch*30;
                for (auto& obs : only<GObs>(rec.obsList))
                {
                    obs.time=state.time;
                    const bool absent=mode==5 && obs.Sat==SatSys("G04") && epoch<40;
                    obs.sigs[freq].L=absent?0:1e8;
                    if (!absent) obs.satStat_ptr->sigStatMap[ft2string(freq)].tracking.lastPhaseTime=state.time.bigTime;
                    if (mode==4 && obs.Sat==SatSys("G03"))
                    {
                        auto& sig=obs.satStat_ptr->sigStatMap[ft2string(freq)];
                        if (epoch==1) { sig.tracking.eventSequence=1; sig.userQcState=ZhangUserQcState::CONFIRMED_BREAK; }
                        else sig.userQcState=ZhangUserQcState::ACCEPTED;
                    }
                }
                if (mode==0 && epoch==1)
                    state.metaDataMap[zhangUserChartIsolationKey(rec.id,E_Sys::GPS,static_cast<int>(code))]="STALE_OLD_CHECKPOINT";
                const auto previousKeys=state.kfIndexMap; const auto previousX=state.x; const auto previousP=state.P;
                updateZhangPppArUserReferences(trace,receivers,state);
                check(!zhangUserChartIsolated(state,rec.id,E_Sys::GPS,static_cast<int>(code)),"no band-wide observation veto");
                if (epoch==1 && (mode==0 || mode==3 || mode==4))
                {
                    check(!state.kfIndexMap.contains(ambiguity(3)),"only unknown/broken coordinate retired");
                    for (const auto& [k,i] : state.kfIndexMap)
                    {
                        check(state.x(i)==previousX(previousKeys.at(k)),"independent surviving marginal mean");
                        for (const auto& [other,j] : state.kfIndexMap)
                            check(state.P(i,j)==previousP(previousKeys.at(k),previousKeys.at(other)),"independent full surviving covariance/cross blocks");
                    }
                }
                if (epoch==1 && mode==1)
                {
                    check(zhangPppArUserReference(state,rec.id,E_Sys::GPS,code)==SatSys("G02"),"continuous replacement anchor");
                    const double wavelength=299792458.0/(code==E_ObsCode::L1C?1575.42e6:1227.60e6);
                    MatrixXd oracle=MatrixXd::Zero(state.x.size(),previousX.size());
                    for (const auto& [k,i] : state.kfIndexMap)
                    {
                        if (k==phase) { oracle(i,previousKeys.at(phase))=1; oracle(i,previousKeys.at(ambiguity(2)))=wavelength; }
                        else { oracle(i,previousKeys.at(k))=1;
                            if (k.type==KF::AMBIGUITY) oracle(i,previousKeys.at(ambiguity(2)))=-1; }
                    }
                    check((state.x-oracle*previousX).norm()<1e-10,"independent S mean oracle");
                    check((state.P-oracle*previousP*oracle.transpose()).norm()<1e-10,"independent S full covariance oracle");
                }
                if (epoch==1 && mode==2)
                {
                    check(!state.kfIndexMap.contains(phase),"all-unknown chart releases unsupported phase coordinate");
                    for (const auto& [k,i] : state.kfIndexMap)
                        for (const auto& [other,j] : state.kfIndexMap)
                            check(state.P(i,j)==previousP(previousKeys.at(k),previousKeys.at(other)),"all-unknown preserves other FLOAT marginal");
                }
                const auto beforeDuplicateX=state.x; const auto beforeDuplicateP=state.P; const auto beforeDuplicateMetadata=state.metaDataMap;
                updateZhangPppArUserReferences(trace,receivers,state);
                check((state.x-beforeDuplicateX).norm()==0 && (state.P-beforeDuplicateP).norm()==0 &&
                    state.metaDataMap==beforeDuplicateMetadata,"duplicate event must not reset fresh chart twice");
                // Nonzero production design entries, not synthetic timer refresh.
                // The actual six-station replay additionally exercises full PPP.
                InitialState init; init.estimate=true; init.P=1000; init.outageLimit=900;
                for (const auto& obs : only<GObs>(rec.obsList))
                {
                    ZhangInternalProduct product;
                    if (!obs.sigs.at(freq).L || !queryZhangInternalProduct(state.time,obs.Sat,code,product)) continue;
                    KFMeasEntry measurement(&state); measurement.addDsgnEntry(phase,1,init);
                    if (!zhangPppArUserReferenceAmbiguity(state,rec.id,obs.Sat,code))
                        measurement.addDsgnEntry(ambiguity(obs.Sat.prn),1,init);
                }
                const auto now=state.time;
                check(state.stateTransition(trace,now),"real production prediction/consolidation");
                check(state.kfIndexMap.contains(phase),"fresh phase coordinate reinitializes");
                if (mode==3 && epoch==1)
                    check(!state.kfIndexMap.contains(ambiguity(3)),"missing product cannot initialize an unsupported coordinate");
                else
                    check(state.kfIndexMap.contains(ambiguity(3)),"observed coordinate recovers without indefinite quarantine");
                if (mode==5 && epoch>=33 && epoch<40)
                    check(!state.kfIndexMap.contains(ambiguity(4)),"genuine observation outage retains finite expiry");
                if (mode==5 && epoch>=40)
                    check(state.kfIndexMap.contains(ambiguity(4)),"returning observation reinitializes numerical state");
                for (int sat=1;sat<=4;++sat)
                    check(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,rec.id,
                        SatSys(E_Sys::GPS,sat).id(),static_cast<int>(code)))==unsigned(mode==4 && sat==3),
                        "product and numerical retirement never invent a physical slip");
                check(state.x.allFinite() && state.P.allFinite(),"finite long-lived state");
            }
            if (mode==5)
            {
                check(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(false,rec.id,"G04",static_cast<int>(code)))>0,
                    "numerical expiry advances state identity");
                check(trace.str().find("cause=NUMERICAL_STATE_OUTAGE")!=std::string::npos,"explicit nonphysical expiry receipt");
            }
            check(trace.str().find("action=ISOLATE_UNKNOWN_PHASE_BLOCK")==std::string::npos,"recovery never globally quarantines phase");
        }
        std::cout<<"USER_PRODUCT_RECOVERY_CONTRACT PASS checks="<<checks
                 <<" scenarios=12 production_reference_prediction_design=1 independent_full_P_oracle=1"
                 <<" long_lifecycle_seconds=1950 genuine_outage_expiry=1 physical_slip_locality=1"
                 <<" stale_isolation_recovery=1 no_fabricated_timer_refresh=1 NIS_computed=0\n";
        return 0;
    }
    catch (const std::exception& e) { std::cerr<<"USER_PRODUCT_RECOVERY_CONTRACT FAIL "<<e.what()<<'\n'; return 1; }
}
