#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangPublishedContinuity.hpp"
#include "common/zhangUserArcGeneration.hpp"
#include "common/zhangUserPhaseQc.hpp"
#include "common/zhangUserProductChart.hpp"
#include "pea/zhangPppAr.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// End-to-end production reference hook. The numerical oracle reconstructs
// physical phase predictions first, then independently defines the new chart;
// it never invokes the production map/reset/continuity decision for expected P.
extern "C" int __wrap_main(int, char**)
{
    unsigned checks=0;
    auto check=[&](bool ok,const char* what)
    { ++checks; if(!ok) throw std::runtime_error(std::string("FAIL ")+what); };
    try
    {
        setenv("ZHANG_USER_QC","CAUSAL_TRACKING",1);
        setenv("ZHANG_USER_FLOAT_TRANSPORT","CONTRACT",1);
        setenv("ZHANG_USER_AR_SEARCH","GONSIS_UDUC",1);
        setenv("ZHANG_USER_UNKNOWN_CHART","RELEASE",1); // explicit previous-policy regression
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

        unsigned serial=0;
        for(auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
        for(int trial=0;trial<26;++trial)
        {
            const bool isolationTrial=trial>=21;
            const int scenario=trial<21?trial:trial==22?6:3;
            setenv("ZHANG_USER_UNKNOWN_CHART",isolationTrial?"ISOLATE":"RELEASE",1);
            // 0 real-only affine; 1 pure re-reference + real affine;
            // 2 unknown old reference with known replacement; 3 all UNKNOWN;
            // 4 CSV49; 5 CSV52 same public segment; 6 CSV52 new public segment;
            // 7 predecessor epoch mismatch; 8 malformed nan declared mapping;
            // 9 predecessor identity mismatch; 10 LEGACY_RELEASE ablation;
            // 11 true explicit slip; 12 true omitted-reference slip;
            // 13 callback rejection; 14 certificate ordered at current epoch;
            // 15 actual receiver/product time reversal after a valid map;
            // 16 current integer components split with proved real transport;
            // 17 reference observation absent, product mapping still known;
            // 18 permuted old stochastic-coordinate insertion order;
            // 19 missing nonreference product while the known anchor shifts.
            // 20 UNKNOWN G03 falsely declares usable without complete joint Qp.
            setenv("ZHANG_USER_FLOAT_TRANSPORT",scenario==10?"LEGACY_RELEASE":"CONTRACT",1);
            const int columns=scenario==4?49:scenario==5 || scenario==6?52:65;
            const auto fixture=std::filesystem::current_path()/
                ("user_float_product_transport_"+std::to_string(++serial)+".csv");
            const std::array<double,4> shifts{.11,-.047,.29,.03};
            std::ofstream products(fixture);
            products << "synthetic frozen CSV" << columns << " production-entry fixture\n";
            for(int epoch=0;epoch<2;++epoch)
            for(int sat=1;sat<=4;++sat)
            {
                if(epoch && scenario==19 && sat==3) continue;
                std::vector<std::string> row(columns,"0");
                row[0]=epoch?"100030":"100000";
                row[1]="PRODUCT_FIXED";
                row[2]=SatSys(E_Sys::GPS,sat).id();
                row[3]=enum_to_string(code);
                row[5]=".01"; row[7]=".01"; row[10]=".02";
                // Temporal numerical values deliberately differ from the
                // same-physical-instant proof shift. No adjacent-delta oracle.
                row[6]=std::to_string(-.35*sat-.17*epoch);
                row[9]=std::to_string(.35*sat+.17*epoch);
                row[11]=epoch?"1":"0"; row[14]=epoch?"2":"1";
                row[15]=row[0]; row[16]="1"; row[17]="SYNTHETIC_PROVED_OR_UNKNOWN_HANDOFF";
                row[18]="1"; row[19]="CURRENT_ALIGNMENT_VALID";
                for(int f : {20,21,22,23,28,29,30,31,32,42,43,46}) row[f]="1";
                if (trial==25 && !epoch) { row[43]="0"; row[46]="0"; } // old FLOAT chart has no legal integer anchor
                row[24]="PRODUCT-GPS-C0"; row[25]="PRODUCT-GPS-C0-DATUM";
                if(epoch && scenario==16 && sat>=3)
                { row[24]="PRODUCT-GPS-C1"; row[25]="PRODUCT-GPS-C1-DATUM"; }
                row[26]=row[0]; row[27]=row[0]; row[33]="NONE";
                row[34]=row[2]+"-"+row[3]+"-PHYSICAL0";
                row[35]="1"; row[38]="4"; row[39]="3"; row[40]="3";
                row[44]="AR_VALID"; row[45]=epoch?"1":"0";
                row[47]="FROZEN_FUNCTIONAL";
                const std::string previousSegment=row[2]+"-"+row[3]+"-PUBLIC0";
                if(columns>=52)
                {
                    const bool sameSegment=!epoch || scenario==5;
                    row[49]=sameSegment?"1":"0";
                    row[50]=sameSegment?"0":"1";
                    row[51]=sameSegment?previousSegment:row[2]+"-"+row[3]+"-PUBLIC1";
                }
                if(columns==65)
                {
                    const bool known=epoch && scenario!=3 && !(scenario==2 && sat==1)
                        && !(scenario==20 && sat==3);
                    row[52]="1"; row[53]=known?"100000":"nan";
                    row[54]=known?"0":"-1"; row[55]=known?"1":"-1";
                    row[56]=known?previousSegment:"UNAVAILABLE";
                    row[57]=known?"1":"0";
                    row[58]=known?std::to_string(shifts[sat-1]):"nan";
                    row[59]="0"; row[60]="nan"; // wire-level unknown, never fabricated zero
                    row[61]=known?"FROZEN_REAL_MAP_PROOF":"UNAVAILABLE";
                    row[62]=known?"CERTIFIED_AFFINE":"UNKNOWN";
                    // UNKNOWN transport is not unknown product covariance:
                    // every legitimate fixture declares the complete joint Qp.
                    row[63]=epoch && scenario==20 && sat==3?"0":"1";
                    row[64]=epoch?"FROZEN_CURRENT_BUNDLE1":"FROZEN_CURRENT_BUNDLE0";
                    if(epoch && scenario==7) row[53]="100010"; // valid time, wrong runtime predecessor
                    if(epoch && scenario==8) row[58]="nan";
                    if(epoch && scenario==9) row[54]="99";
                    if(epoch && scenario==14) row[53]="100030"; // not strictly preceding
                }
                for(std::size_t f=0;f<row.size();++f) products << (f?",":"") << row[f];
                products << '\n';
            }
            products.close();
            acsConfig.zhangPppAr.product_filename=fixture.string();
            acsConfig.zhangPppAr.baseline_observables.clear();
            acsConfig.zhangPppAr.baseline_observables[E_Sys::GPS]={code};
            ReceiverMap receivers;
            auto& receiver=receivers["FLOATTEST"];
            receiver.id="FLOATTEST"; receiver.ready=true;
            const auto frequency=code2Freq.at(E_Sys::GPS).at(code);
            for(int sat=1;sat<=4;++sat)
            {
                const SatSys satellite(E_Sys::GPS,sat);
                auto& status=receiver.satStatMap[satellite];
                status.el=(70-sat*10)*D2R;
                auto& signal=status.sigStatMap[ft2string(frequency)];
                signal.tracking.seen=true;
                signal.tracking.lastCode=static_cast<int>(code);
                signal.tracking.lastPhaseTime=100000;
                GObs obs; obs.time.bigTime=100000; obs.Sat=satellite;
                obs.mount=receiver.id; obs.satStat_ptr=&status;
                obs.sigs[frequency].code=code; obs.sigs[frequency].P=2e7; obs.sigs[frequency].L=1e8;
                receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
            }
            auto ambiguity=[&](int sat)
            {
                KFKey key; key.type=KF::AMBIGUITY; key.str=receiver.id;
                key.Sat=SatSys(E_Sys::GPS,sat); key.num=static_cast<int>(code); return key;
            };
            KFKey clock; clock.type=KF::REC_CLOCK; clock.str=receiver.id;
            KFKey position; position.type=KF::REC_POS; position.str=receiver.id; position.num=0;
            KFKey phase; phase.type=KF::PHASE_BIAS; phase.str=receiver.id;
            phase.Sat=SatSys(E_Sys::GPS,0); phase.num=static_cast<int>(code);
            KFState state; state.time.bigTime=100000;
            check(bindZhangCheckpointRuntimeId(state,"float-map-"+std::to_string(serial)),"bind distinct production runtime");
            state.kfIndexMap.clear(); state.stateTransitionMap.clear();
            const std::vector<KFKey> canonicalKeys{state.oneKey,clock,position,phase,ambiguity(2),ambiguity(3),ambiguity(4)};
            auto keys=canonicalKeys;
            if(scenario==18) keys={state.oneKey,ambiguity(3),position,ambiguity(4),clock,ambiguity(2),phase};
            for(int i=0;i<7;++i)
            {
                state.kfIndexMap[keys[i]]=i;
                state.stateTransitionMap[keys[i]][keys[i]][0]=1;
            }
            VectorXd canonicalX(7); canonicalX << 1,7,1.2,20,3.25,-2.5,4.3;
            state.x.resize(7);
            state.dx=VectorXd::Zero(7);
            MatrixXd factor=MatrixXd::Identity(7,7);
            factor(1,3)=.23; factor(2,4)=.31; factor(3,5)=.17;
            factor(4,1)=-.28; factor(5,2)=.19; factor(6,3)=-.11;
            MatrixXd canonicalP=factor*factor.transpose();
            canonicalP.row(0).setZero(); canonicalP.col(0).setZero();
            state.P.resize(7,7);
            std::array<int,7> canonicalIndex{};
            for(int i=0;i<7;++i)
            {
                for(int j=0;j<7;++j) if(keys[i]==canonicalKeys[j]) canonicalIndex[i]=j;
                state.x(i)=canonicalX(canonicalIndex[i]);
            }
            for(int i=0;i<7;++i) for(int j=0;j<7;++j)
                state.P(i,j)=canonicalP(canonicalIndex[i],canonicalIndex[j]);
            std::ostringstream trace;
            updateZhangPppArUserReferences(trace,receivers,state);
            check(zhangPppArUserReference(state,receiver.id,E_Sys::GPS,code)==SatSys("G01"),"initial reference exists");
            check(!state.kfIndexMap.contains(ambiguity(1)),"initial S-basis omits old reference N");
            state.time.bigTime=100030;
            for(auto& obs : only<GObs>(receiver.obsList))
            {
                obs.time=state.time;
                obs.satStat_ptr->sigStatMap[ft2string(frequency)].tracking.lastPhaseTime=100030;
            }
            if(scenario==1) receiver.satStatMap[SatSys("G01")].sigStatMap[ft2string(frequency)].phaseQuarantined=true;
            if(scenario==17) for(auto& obs : only<GObs>(receiver.obsList))
                if(obs.Sat==SatSys("G01")) obs.sigs.clear();
            if(scenario==11 || scenario==12 || trial==23 || trial==24)
            {
                auto& signal=receiver.satStatMap[SatSys(E_Sys::GPS,scenario==11 || trial==23?3:1)].sigStatMap[ft2string(frequency)];
                signal.tracking.eventSequence=1;
                signal.userQcState=ZhangUserQcState::CONFIRMED_BREAK;
            }
            const auto beforeX=state.x; const auto beforeP=state.P; const auto beforeDx=state.dx;
            const auto beforeMetadata=state.metaDataMap;
            const auto beforeKeys=state.kfIndexMap; const auto beforeTransitions=state.stateTransitionMap;
            const auto beforeSequence=state.factorCommitSequence;
            if(scenario==20)
            {
                ZhangInternalProduct invalid;
                check(!queryZhangInternalProduct(state.time,SatSys("G03"),code,invalid),
                    "UNKNOWN transport plus incomplete joint Qp cannot declare a usable product");
                check(!invalid.ppp_usable && !invalid.pppar_usable,
                    "contradictory CSV65 withdraws both current PPP and AR eligibility");
            }
            unsigned rejectedCallbacks=0;
            if(scenario==13) state.exactStateTransformCallback=[&](auto&&...)
                { ++rejectedCallbacks; return false; };
            bool threw=false;
            try { updateZhangPppArUserReferences(trace,receivers,state); }
            catch(const std::runtime_error&) { threw=true; }
            if(scenario==13)
            {
                check(threw && rejectedCallbacks>0,"production transform callback rejection reaches atomic rollback");
                check((state.x-beforeX).norm()==0 && (state.P-beforeP).norm()==0 && (state.dx-beforeDx).norm()==0,"rollback complete numerical state");
                check(state.metaDataMap==beforeMetadata && state.kfIndexMap==beforeKeys && state.stateTransitionMap==beforeTransitions,"rollback identities indices and persistence");
                check(state.factorCommitSequence==beforeSequence,"rollback factor chronology");
                check(zhangPppArUserReference(state,receiver.id,E_Sys::GPS,code)==SatSys("G01"),"rollback reference registry");
                continue;
            }
            check(!threw,"production transport/fallback transaction accepted");
            if (isolationTrial)
            {
                const bool retain=trial!=24 && trial!=25;
                check(zhangUserChartIsolated(state,receiver.id,E_Sys::GPS,static_cast<int>(code))==retain,
                    "unknown chart isolation requires an unbroken legally identified old anchor");
                check(state.kfIndexMap.contains(phase)==retain,"unknown chart keeps legal prior or releases unsupported coordinates");
                if (retain)
                {
                    check(state.kfIndexMap.contains(ambiguity(3))==(trial!=23),"real explicit slip still retires its ambiguity during product isolation");
                    for (const auto& [key,index] : state.kfIndexMap)
                    {
                        check(state.x(index)==beforeX(beforeKeys.at(key)),"isolation does not relabel old FLOAT means");
                        for (const auto& [other,column] : state.kfIndexMap)
                            check(state.P(index,column)==beforeP(beforeKeys.at(key),beforeKeys.at(other)),"isolation retains full marginal and cross blocks");
                    }
                    const std::string receipt="USER_PRODUCT_PUBLIC_SEGMENT:"+receiver.id+":"+
                        std::to_string(static_cast<int>(E_Sys::GPS))+":"+std::to_string(static_cast<int>(code))+":G01";
                    check(state.metaDataMap.at(receipt+":EPOCH")=="100000" &&
                        state.metaDataMap.at(receipt+":RECEIVED_EPOCH")=="100030" &&
                        state.metaDataMap.at(receipt+":STATUS")=="UNRESOLVED","received event never becomes an applied chart receipt");
                    check(zhangPppArUserAmbiguityIntegerValid(state,receiver.id,SatSys("G02"),code),
                        "fresh AR can use the legally identified old integer prior without current chart reinterpretation");
                    if (trial==23) check(!zhangUserGonsisComponentCertificateValid(
                        zhangUserChartAnchor(state,receiver.id,SatSys("G03"),static_cast<int>(code))),"physical slip revokes old integer coordinate anchor");
                    const auto isolatedX=state.x; const auto isolatedP=state.P; const auto isolatedMetadata=state.metaDataMap;
                    updateZhangPppArUserReferences(trace,receivers,state);
                    check((state.x-isolatedX).norm()==0 && (state.P-isolatedP).norm()==0 &&
                        state.metaDataMap==isolatedMetadata && zhangUserChartIsolated(state,receiver.id,E_Sys::GPS,static_cast<int>(code)),
                        "duplicate unresolved receipt stays unresolved and is idempotent");
                }
                continue;
            }
            const bool mapped=scenario==0 || scenario==1 || scenario==2 || scenario==11 || scenario==12 || scenario>=15;
            const bool identity52=scenario==5;
            const int expectedReference=scenario==1 || scenario==2 || scenario==12 || scenario==17?2:scenario==8 || scenario==14?0:1;
            const auto actualReference=zhangPppArUserReference(state,receiver.id,E_Sys::GPS,code);
            check(expectedReference?actualReference==SatSys(E_Sys::GPS,expectedReference):actualReference.prn==0,"legal known anchor or fail-closed reference");
            const bool block=!mapped && !identity52;
            check(state.kfIndexMap.contains(phase)==!block,"unknown directions are truly released, not kept by relabelling");
            if(mapped)
            {
                check(state.kfIndexMap.contains(ambiguity(3))==(scenario!=11 && scenario!=19 && scenario!=20),"true-slip or missing-product unknown coordinate locally retired");
                check(state.kfIndexMap.contains(ambiguity(1))==(scenario==1 || scenario==17),"unknown/broken old reference never gains an inherited N");
                check(state.kfIndexMap.contains(ambiguity(4)),"known fourth physical direction retained");
            }
            // Independent physical constants, not the production wavelength
            // helper as an oracle for itself.
            const double lambda=299792458.0/(code==E_ObsCode::L1C?1575.42e6:1227.60e6);
            check(std::abs(zhangPppArUserPhaseCoordinateWavelength(E_Sys::GPS,code)-lambda)<1e-15,
                "production wavelength agrees with independent GPS carrier definition");
            // Independent frozen physical-space representation [clock,pos,phi1..4].
            MatrixXd physical=MatrixXd::Zero(6,7);
            physical(0,beforeKeys.at(clock))=1; physical(1,beforeKeys.at(position))=1;
            for(int sat=1;sat<=4;++sat)
            {
                physical(1+sat,beforeKeys.at(phase))=1;
                if(sat>1) physical(1+sat,beforeKeys.at(ambiguity(sat)))=lambda;
            }
            MatrixXd independent=MatrixXd::Zero(state.x.size(),7);
            VectorXd affine=VectorXd::Zero(state.x.size());
            for(const auto& [key,index] : state.kfIndexMap)
            {
                if(key==state.oneKey) independent(index,beforeKeys.at(state.oneKey))=1;
                else if(key==clock) independent.row(index)=physical.row(0);
                else if(key==position) independent.row(index)=physical.row(1);
                else if(key==phase)
                {
                    independent.row(index)=physical.row(1+expectedReference);
                    if(mapped) affine(index)=shifts[expectedReference-1];
                }
                else if(key.type==KF::AMBIGUITY)
                {
                    independent.row(index)=(physical.row(1+key.Sat.prn)-physical.row(1+expectedReference))/lambda;
                    if(mapped) affine(index)=(shifts[key.Sat.prn-1]-shifts[expectedReference-1])/lambda;
                }
                else throw std::runtime_error("unexpected destination key");
            }
            check((state.x-independent*beforeX-affine).cwiseAbs().maxCoeff()<1e-9,"independent complete physical-space affine mean oracle");
            check((state.P-independent*beforeP*independent.transpose()).cwiseAbs().maxCoeff()<1e-9,"independent full covariance and cross-covariance oracle");
            if (mapped) for (const auto& [key,index] : state.kfIndexMap)
                if (key.type==KF::AMBIGUITY)
                {
                    VectorXd integerRow=VectorXd::Zero(state.x.size()); integerRow(index)=1;
                    check(zhangUserCanonicalAmbiguityRow(state,receiver.id,key.Sat,static_cast<int>(code),integerRow),
                        "transported integer row is representable");
                    check(std::abs(integerRow.dot(state.x)-(independent*beforeX)(index))<1e-9,
                        "fractional chart shift is removed from actual integer search mean");
                }
            const int qi=state.kfIndexMap.at(position),ci=state.kfIndexMap.at(clock);
            check(std::abs(state.P(qi,qi)-beforeP(beforeKeys.at(position),beforeKeys.at(position)))<1e-12 && std::abs(state.P(qi,ci)-beforeP(beforeKeys.at(position),beforeKeys.at(clock)))<1e-12,"retained q is a marginal, not a conditional Schur complement");
            for(int sat=1;sat<=4;++sat)
            {
                const auto id=SatSys(E_Sys::GPS,sat).id();
                const auto physicalGeneration=zhangUserArcGeneration(state.metaDataMap,
                    zhangUserArcGenerationKey(true,receiver.id,id,static_cast<int>(code)));
                check(physicalGeneration==static_cast<unsigned>((scenario==11 && sat==3) || (scenario==12 && sat==1)),"only confirmed tracking break changes physical generation");
                if(scenario==0 || scenario==1 || (scenario>=15 && scenario!=19 && scenario!=20))
                {
                    check(zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(false,receiver.id,id,static_cast<int>(code)))==0,"known FLOAT map preserves stochastic state identity");
                    const std::string integerKey="USER_AR_INTEGER_AUTHORIZATION_GENERATION:"+receiver.id+":"+id+":"+std::to_string(static_cast<int>(code));
                    check(zhangUserArcGeneration(state.metaDataMap,integerKey)==1,"real-only map withdraws old integer authorization without deleting FLOAT");
                }
                if(scenario==19 || scenario==20)
                    check(zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(false,receiver.id,id,static_cast<int>(code)))==static_cast<unsigned>(sat==3),
                        "missing product retires only affected state identity, not all physical arcs");
            }
            const auto stableX=state.x; const auto stableP=state.P;
            const auto stableMetadata=state.metaDataMap; const auto stableSequence=state.factorCommitSequence;
            updateZhangPppArUserReferences(trace,receivers,state);
            check(state.x.size()==stableX.size() && (state.x-stableX).norm()==0 && (state.P-stableP).norm()==0,"same event cannot shift twice or release twice");
            check(state.metaDataMap==stableMetadata && state.factorCommitSequence==stableSequence,"same event authorization and factor chronology idempotence");
            if(scenario==15)
            {
                state.time.bigTime=100000;
                for(auto& obs : only<GObs>(receiver.obsList)) obs.time=state.time;
                bool reversed=false;
                try { updateZhangPppArUserReferences(trace,receivers,state); }
                catch(const std::runtime_error& e)
                { reversed=std::string(e.what()).find("USER_PRODUCT_TIME_REVERSED")!=std::string::npos; }
                check(reversed,"actual time reversal fails closed");
                check((state.x-stableX).norm()==0 && (state.P-stableP).norm()==0 && state.metaDataMap==stableMetadata,"time reversal restores mean full P and identities");
                check(state.factorCommitSequence==stableSequence && zhangPppArUserReference(state,receiver.id,E_Sys::GPS,code)==SatSys("G01"),"time reversal restores chronology and reference registry");
            }
        }
        unsetenv("ZHANG_USER_FLOAT_TRANSPORT");
        check(serial==52,"all 42 legacy and 10 new isolation scenarios executed");
        std::cout << "USER_FLOAT_PRODUCT_TRANSPORT_CONTRACT PASS checks=" << checks
                  << " scenarios=52 unknown_chart_isolation=1 received_applied_split=1 affine_integer_map=1"
                  << " production_reference_hook=1 physical_space_oracle=1 full_covariance=1"
                  << " real_only_integer_authorization_split=1 unknown_marginal_release=1"
                  << " repeated_event_idempotence=1 physical_slip_locality=1 callback_rollback=1"
                  << " schema49_52_65=1 malformed_ordered_identity_checked=1 time_reversal_atomic_checked=1"
                  << " component_split_float_checked=1 absent_reference_checked=1 insertion_order_checked=1"
                  << " absent_nonreference_product_checked=1 contradictory_joint_covariance_checked=1 NIS_computed=0\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e)
    {
        std::cerr << "USER_FLOAT_PRODUCT_TRANSPORT_CONTRACT FAIL " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
