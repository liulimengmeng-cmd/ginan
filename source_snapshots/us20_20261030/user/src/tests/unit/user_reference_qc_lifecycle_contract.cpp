#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/algebraTrace.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserArcGeneration.hpp"
#include "common/zhangUserPhaseQc.hpp"
#include "pea/zhangPppAr.hpp"

#include <cstdlib>
#include <cerrno>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/wait.h>
#include <unistd.h>

// Linked against production zhangPppAr/algebra/parser objects. No reference
// selector/transform/reset is reimplemented here. The oracle below is an
// independently specified dense linear map on the frozen input mean/full P.
extern "C" int __wrap_main(int argc,char** argv)
{
    unsigned checks = 0;
    int isolatedScenario=-1;
    int scenarioResultFd=-1;
    auto check = [&](bool value, const char* label)
    {
        ++checks;
        if (!value) throw std::runtime_error(std::string("FAIL ")+label);
    };
    try
    {
        for(int argument=1;argument<argc;++argument)
        {
            const std::string option=argv[argument];
            const std::string scenarioPrefix="--live-capture-scenario=";
            const std::string resultPrefix="--scenario-result-fd=";
            if(option.rfind(scenarioPrefix,0)==0)
                isolatedScenario=std::stoi(option.substr(scenarioPrefix.size()));
            else if(option.rfind(resultPrefix,0)==0)
                scenarioResultFd=std::stoi(option.substr(resultPrefix.size()));
            else throw std::runtime_error("unknown reference lifecycle contract option");
        }
        if(isolatedScenario>=0 && (isolatedScenario>=7 || scenarioResultFd<0))
            throw std::runtime_error("invalid isolated reference lifecycle scenario");
        setenv("ZHANG_USER_QC", "CAUSAL_TRACKING", 1);
        setenv("ZHANG_USER_UNKNOWN_CHART", "RELEASE", 1); // legacy reset/reference scenarios
        acsConfig.process_ppp = true;
        acsConfig.zhangPppAr.user_adapter = true;
        acsConfig.zhangPppAr.output_products = false;
        acsConfig.zhangPppAr.product_solution = "PRODUCT_FIXED";
        acsConfig.zhangPppAr.integer_strategy = "CANONICAL_USER_SD_WL_L1";
        acsConfig.zhangPppAr.user_use_full_product_covariance = false;
        acsConfig.zhangFullRank.enable = false;
        acsConfig.phaseClockOsb.enable = false;
        acsConfig.pppOpts.ionoOpts.use_if_combo = false;
        const auto fixture = std::filesystem::current_path()/"user_reference_qc_fixture.csv";
        std::ofstream products(fixture);
        products << "gpst_seconds,solution,satellite,observable,clock_m,clock_sigma_m,phase_m,phase_sigma_m,clock_phase_covariance_m2,correction_m,correction_sigma_m,discontinuity_counter,integer_shift_cycles,fractional_shift_cycles,datum_version,valid_from,product_iod,reset_reason,integer_valid\n";
        for (int satellite=1;satellite<=4;++satellite)
            for (auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                products << "100000,PRODUCT_FIXED,G0" << satellite << ',' << enum_to_string(code)
                         << ",0,.01,0,.01,0,0,.02,0,0,0,1,100000,1,NONE,1\n";
        products.close();
        acsConfig.zhangPppAr.product_filename = fixture.string();
        unsigned scenarioId = 0;
        if(isolatedScenario<0)
        {
        for (auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
            for (int scenario=0;scenario<4;++scenario)
            {
                acsConfig.zhangPppAr.baseline_observables.clear();
                acsConfig.zhangPppAr.baseline_observables[E_Sys::GPS] = {code};
                ReceiverMap receivers;
                auto& receiver = receivers["TEST"];
                receiver.id = "TEST";
                receiver.ready = true;
                const auto frequency = code2Freq.at(E_Sys::GPS).at(code);
                for (int satellite=1;satellite<=4;++satellite)
                {
                    const SatSys sat(E_Sys::GPS,satellite);
                    auto& status = receiver.satStatMap[sat];
                    status.el = (70-satellite*10)*D2R;
                    if (satellite==4) status.el=10*D2R;
                    auto& sig = status.sigStatMap[ft2string(frequency)];
                    sig.tracking.seen=true;
                    sig.tracking.lastCode=static_cast<int>(code);
                    sig.tracking.lastPhaseTime=100000;
                    GObs obs;
                    obs.time.bigTime=100000;
                    obs.Sat=sat;
                    obs.mount=receiver.id;
                    obs.satStat_ptr=&status;
                    obs.sigs[frequency].code=code;
                    obs.sigs[frequency].P=2e7;
                    obs.sigs[frequency].L=1e8;
                    receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
                }
                auto ambiguity = [&](int satellite)
                {
                    KFKey key;
                    key.type=KF::AMBIGUITY;
                    key.str=receiver.id;
                    key.Sat=SatSys(E_Sys::GPS,satellite);
                    key.num=static_cast<int>(code);
                    return key;
                };
                KFKey clock;
                clock.type=KF::REC_CLOCK;
                clock.str=receiver.id;
                KFKey phase;
                phase.type=KF::PHASE_BIAS;
                phase.str=receiver.id;
                phase.Sat=SatSys(E_Sys::GPS,0);
                phase.num=static_cast<int>(code);
                KFState state;
                state.time.bigTime=100000;
                check(bindZhangCheckpointRuntimeId(state,"reference-qc-"+std::to_string(++scenarioId)),"bind owner");
                state.kfIndexMap.clear();
                const std::vector<KFKey> beforeKeys{state.oneKey,clock,phase,ambiguity(2),ambiguity(3)};
                for (int i=0;i<5;++i)
                {
                    state.kfIndexMap[beforeKeys[i]]=i;
                    state.stateTransitionMap[beforeKeys[i]][beforeKeys[i]][0]=1;
                }
                state.x.resize(5);
                state.x << 1,7,20,3.25,-2.5;
                state.dx=VectorXd::Zero(5);
                MatrixXd factor=MatrixXd::Identity(5,5);
                factor(1,2)=.23; factor(2,3)=.31; factor(3,4)=.17;
                factor(4,1)=-.28; factor(4,2)=.19;
                state.P=factor*factor.transpose();
                state.P.row(0).setZero(); state.P.col(0).setZero();
                std::ostringstream trace;
                updateZhangPppArUserReferences(trace,receivers,state);
                check(zhangPppArUserReference(state,"TEST",E_Sys::GPS,code)==SatSys("G01"),"initial omitted ref");
                check(!state.kfIndexMap.contains(ambiguity(1)),"reference has no N key");
                const auto beforeX=state.x;
                const auto beforeP=state.P;
                auto& old=receiver.satStatMap[SatSys("G01")].sigStatMap[ft2string(frequency)];
                // Native cleanup clears transient bool but not tracking serial.
                old.userPhysicalBreakThisEpoch=false;
                old.tracking.pendingBreak=false;
                if(scenario==0 || scenario==2 || scenario==3)
                {
                    old.tracking.eventSequence=1;
                    old.userQcState=ZhangUserQcState::CONFIRMED_BREAK;
                }
                else
                {
                    old.phaseQuarantined=true;
                    old.userQcState=ZhangUserQcState::SUSPECT;
                }
                // Fresh highest-elevation star lacks a state; lower eligible
                // G02 must be selected if it can transport the existing chart.
                receiver.satStatMap[SatSys("G04")].el=80*D2R;
                if(scenario==1)
                    receiver.satStatMap[SatSys("G04")].sigStatMap[ft2string(frequency)].tracking.halfCycle=true;
                if(scenario==2)
                    for(int sat=2;sat<=4;++sat)
                        receiver.satStatMap[SatSys(E_Sys::GPS,sat)].sigStatMap[ft2string(frequency)].phaseQuarantined=true;
                if(scenario==3)
                {
                    // No state-backed alternative, despite a valid fresh star.
                    receiver.satStatMap[SatSys("G02")].sigStatMap[ft2string(frequency)].phaseQuarantined=true;
                    receiver.satStatMap[SatSys("G03")].sigStatMap[ft2string(frequency)].phaseQuarantined=true;
                }
                updateZhangPppArUserReferences(trace,receivers,state);
                if(scenario<2)
                {
                    check(zhangPppArUserReference(state,"TEST",E_Sys::GPS,code)==SatSys("G02"),"legal backed ref selected");
                    check(!state.kfIndexMap.contains(ambiguity(2)),"new reference omitted");
                    check(state.kfIndexMap.contains(ambiguity(3)),"surviving SD retained");
                    check(state.kfIndexMap.contains(ambiguity(1))==(scenario==1),"broken ref does not inherit new N");
                    const double lambda=zhangPppArUserPhaseCoordinateWavelength(E_Sys::GPS,code);
                    MatrixXd independent=MatrixXd::Zero(state.x.size(),5);
                    for(const auto& [key,index]:state.kfIndexMap)
                    {
                        if(key==phase) {independent(index,2)=1;independent(index,3)=lambda;}
                        else if(key==ambiguity(3)) {independent(index,4)=1;independent(index,3)=-1;}
                        else if(key==ambiguity(1)) independent(index,3)=-1;
                        else if(key==clock) independent(index,1)=1;
                        else if(key==state.oneKey) independent(index,0)=1;
                        else throw std::runtime_error("unexpected transformed state");
                    }
                    check((state.x-independent*beforeX).cwiseAbs().maxCoeff()<1e-10,"complete mean oracle");
                    check((state.P-independent*beforeP*independent.transpose()).cwiseAbs().maxCoeff()<1e-10,"complete covariance oracle");
                    check(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,"TEST","G01",static_cast<int>(code)))==(scenario==0?1:0),"only confirmed physical ref generation");
                    for(auto satellite : {"G02","G03"})
                        check(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,"TEST",satellite,static_cast<int>(code)))==0,"unbroken endpoints not retired");
                }
                else
                {
                    check(!state.kfIndexMap.contains(phase),"no alternative phase reset");
                    check(!state.kfIndexMap.contains(ambiguity(2)) && !state.kfIndexMap.contains(ambiguity(3)),"fallback stale ambiguities removed");
                    check(state.kfIndexMap.contains(clock),"fallback unrelated clock retained");
                    const auto i=state.kfIndexMap.at(clock);
                    check(state.x(i)==beforeX(1) && state.P(i,i)==beforeP(1,1),"fallback unrelated mean/P unchanged");
                }
                const auto stableX=state.x;
                const auto stableP=state.P;
                const auto stableMetadata=state.metaDataMap;
                updateZhangPppArUserReferences(trace,receivers,state);
                check(state.x.size()==stableX.size() && (state.x-stableX).norm()==0,"duplicate event mean idempotence");
                check((state.P-stableP).norm()==0,"duplicate event covariance idempotence");
                check(state.metaDataMap==stableMetadata,"duplicate event generation idempotence");
                if(scenario==0)
                {
                    InitialState fresh;
                    fresh.x=88.75; fresh.P=100; fresh.tau=0;
                    check(state.addKFState(ambiguity(1),fresh),"new arc adds fresh prior");
                    check(state.stateTransition(trace,state.time),"fresh prior consolidation");
                    const auto i=state.kfIndexMap.at(ambiguity(1));
                    check(std::abs(state.x(i)-88.75)<1e-10 && std::abs(state.P(i,i)-100)<1e-10,"fresh prior no old integer inheritance");
                    for(const auto& [key,j]:state.kfIndexMap) if(j!=i)
                        check(std::abs(state.P(i,j))<1e-10,"fresh prior independent of retained state");
                    // The original reference is now an explicit native N.
                    // Mirror native cleanup's exact receipt after the real
                    // production removeState/stateTransition operation; the
                    // reference updater must not retire/increment it twice.
                    old.tracking.eventSequence=2;
                    state.removeState(ambiguity(1));
                    state.metaDataMap[zhangUserQcAppliedSequenceKey("TEST","G01",static_cast<int>(code))]="2";
                    check(state.stateTransition(trace,state.time),"native explicit removal consolidated");
                    const auto nativeGeneration=zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(true,"TEST","G01",static_cast<int>(code)));
                    check(nativeGeneration==2,"native physical removal incremented once");
                    const auto nativeX=state.x;
                    const auto nativeP=state.P;
                    updateZhangPppArUserReferences(trace,receivers,state);
                    check((state.x-nativeX).norm()==0 && (state.P-nativeP).norm()==0,"native ack event no second transform");
                    check(zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(true,"TEST","G01",static_cast<int>(code)))==nativeGeneration,
                        "native ack event no second physical generation");
                }
            }
        // Product coordinates are not physical receiver arcs. These fixtures
        // exercise the production pre-measurement hook with independent dense
        // mean/full-covariance oracles, including clock/phase/N correlations.
        for (auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
            for (int scenario=0;scenario<11;++scenario)
            {
                const bool legacy49 = scenario==7;
                const auto fixture = std::filesystem::current_path()/
                    ("user_product_datum_fixture_"+std::to_string(++scenarioId)+".csv");
                std::ofstream products(fixture);
                products << "synthetic production CSV" << (legacy49 ? 49 : 52) << '\n';
                for(int epoch=0;epoch<2;++epoch)
                    for(int satellite=1;satellite<=4;++satellite)
                    {
                        if(epoch==1 && scenario==5 && satellite==1) continue;
                        const bool changed = epoch==1 &&
                            ((scenario==1) || (scenario==4 ? satellite==3 :
                              (scenario!=2 && satellite==1)));
                        std::vector<std::string> row(legacy49 ? 49 : 52,"0");
                        row[0]=epoch==0 ? "100000" : "100030";
                        row[1]="PRODUCT_FIXED";
                        row[2]="G0"+std::to_string(satellite);
                        row[3]=enum_to_string(code);
                        row[5]=".01"; row[7]=".01"; row[10]=".02";
                        row[11]=(changed && scenario!=8) || (epoch==1 && scenario==2) ? "1" : "0";
                        row[14]=(changed && scenario!=8) || (epoch==1 && scenario==2) ? "2" : "1";
                        row[15]=row[0]; row[16]="1"; row[17]="SYNTHETIC_PRODUCT_EVENT";
                        row[18]="1"; row[19]="CURRENT_ALIGNMENT_VALID";
                        for(int field : {20,21,22,23,28,29,30,31,32,42,43,46}) row[field]="1";
                        row[24]="PRODUCT-GPS-C0"; row[25]="PRODUCT-GPS-C0-DATUM";
                        row[26]=row[0]; row[27]=row[0]; row[33]="NONE";
                        row[34]=row[2]+"-"+row[3]+"-PHYSICAL_SEG0";
                        row[35]="1"; row[38]="4"; row[39]="3"; row[40]="3";
                        row[44]="AR_VALID"; row[45]=changed ? "1" : "0";
                        row[47]="SYNTHETIC_FROZEN_FUNCTIONAL";
                        if(!legacy49)
                        {
                            row[49]=changed ? "0" : "1";
                            row[50]=changed ? "1" : "0";
                            row[51]=row[2]+"-"+row[3]+(changed ? "-PUBLIC1" : "-PUBLIC0");
                            if(epoch==1 && satellite==1 && scenario==6) row[49]="1";
                            if(epoch==1 && satellite==1 && scenario==9) row[51]="UNRESOLVED";
                        }
                        for(std::size_t field=0;field<row.size();++field)
                            products << (field==0 ? "" : ",") << row[field];
                        products << '\n';
                    }
                products.close();
                acsConfig.zhangPppAr.product_filename=fixture.string();
                acsConfig.zhangPppAr.baseline_observables.clear();
                acsConfig.zhangPppAr.baseline_observables[E_Sys::GPS]={code};
                ReceiverMap receivers;
                auto& receiver=receivers["TEST"];
                receiver.id="TEST"; receiver.ready=true;
                const auto frequency=code2Freq.at(E_Sys::GPS).at(code);
                for(int satellite=1;satellite<=4;++satellite)
                {
                    const SatSys sat(E_Sys::GPS,satellite);
                    auto& status=receiver.satStatMap[sat];
                    status.el=(70-satellite*10)*D2R;
                    auto& signal=status.sigStatMap[ft2string(frequency)];
                    signal.tracking.seen=true;
                    signal.tracking.lastCode=static_cast<int>(code);
                    signal.tracking.lastPhaseTime=100000;
                    GObs obs; obs.time.bigTime=100000; obs.Sat=sat; obs.mount=receiver.id;
                    obs.satStat_ptr=&status;
                    obs.sigs[frequency].code=code;
                    obs.sigs[frequency].P=2e7; obs.sigs[frequency].L=1e8;
                    receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
                }
                auto ambiguity=[&](int satellite)
                {
                    KFKey key; key.type=KF::AMBIGUITY; key.str=receiver.id;
                    key.Sat=SatSys(E_Sys::GPS,satellite); key.num=static_cast<int>(code);
                    return key;
                };
                KFKey clock; clock.type=KF::REC_CLOCK; clock.str=receiver.id;
                KFKey phase; phase.type=KF::PHASE_BIAS; phase.str=receiver.id;
                phase.Sat=SatSys(E_Sys::GPS,0); phase.num=static_cast<int>(code);
                KFState state; state.time.bigTime=100000;
                check(bindZhangCheckpointRuntimeId(state,"product-datum-"+std::to_string(scenarioId)),"product bind owner");
                state.kfIndexMap.clear();
                const std::vector<KFKey> inputKeys{state.oneKey,clock,phase,ambiguity(2),ambiguity(3)};
                for(int index=0;index<5;++index)
                {
                    state.kfIndexMap[inputKeys[index]]=index;
                    state.stateTransitionMap[inputKeys[index]][inputKeys[index]][0]=1;
                }
                state.x.resize(5); state.x << 1,7,20,3.25,-2.5;
                state.dx=VectorXd::Zero(5);
                MatrixXd factor=MatrixXd::Identity(5,5);
                factor(1,2)=.23; factor(2,3)=.31; factor(3,4)=.17;
                factor(4,1)=-.28; factor(4,2)=.19;
                state.P=factor*factor.transpose();
                state.P.row(0).setZero(); state.P.col(0).setZero();
                std::ostringstream trace;
                updateZhangPppArUserReferences(trace,receivers,state);
                check(zhangPppArUserReference(state,"TEST",E_Sys::GPS,code)==SatSys("G01"),"product initial omitted reference");
                state.time.bigTime=100030;
                for(auto& obs : only<GObs>(receiver.obsList))
                {
                    obs.time.bigTime=100030;
                    obs.satStat_ptr->sigStatMap[ft2string(frequency)].tracking.lastPhaseTime=100030;
                }
                if(scenario==3)
                {
                    auto& old=receiver.satStatMap[SatSys("G01")].sigStatMap[ft2string(frequency)];
                    old.tracking.eventSequence=1;
                    old.userQcState=ZhangUserQcState::CONFIRMED_BREAK;
                }
                if(scenario==10) state.procNoiseMap[phase]=1;
                const auto beforeX=state.x;
                const auto beforeP=state.P;
                const auto beforeMetadata=state.metaDataMap;
                if(scenario==6 || scenario==9)
                {
                    ZhangInternalProduct malformed;
                    check(!queryZhangInternalProduct(state.time,SatSys("G01"),code,malformed),"malformed CSV52 fail closed");
                    check(malformed.invalid_reason=="PUBLIC_CONTINUITY_METADATA_INVALID","malformed public metadata reason");
                }
                bool threw=false;
                try { updateZhangPppArUserReferences(trace,receivers,state); }
                catch(const std::runtime_error&) { threw=true; }
                if(scenario==10)
                {
                    check(threw,"unsupported dynamic chart rejects atomically");
                    check((state.x-beforeX).norm()==0 && (state.P-beforeP).norm()==0,"failed product transaction mean full P unchanged");
                    check(state.metaDataMap==beforeMetadata,"failed product transaction metadata unchanged");
                    check(zhangPppArUserReference(state,"TEST",E_Sys::GPS,code)==SatSys("G01"),"failed transaction reference receipt restored");
                    continue;
                }
                check(!threw,"product transaction succeeds");
                const bool keepsReference=scenario==2 || scenario==4;
                const int expectedReference=keepsReference || scenario==1 ? 1 : 2;
                check(zhangPppArUserReference(state,"TEST",E_Sys::GPS,code)==SatSys(E_Sys::GPS,expectedReference),"continuous state-backed anchor selected");
                MatrixXd independent=MatrixXd::Zero(state.x.size(),5);
                const double lambda=zhangPppArUserPhaseCoordinateWavelength(E_Sys::GPS,code);
                for(const auto& [key,index] : state.kfIndexMap)
                {
                    if(key==state.oneKey) independent(index,0)=1;
                    else if(key==clock) independent(index,1)=1;
                    else if(key==phase)
                    {
                        independent(index,2)=1;
                        if(!keepsReference) independent(index,3)=lambda;
                    }
                    else if(key==ambiguity(2) && keepsReference) independent(index,3)=1;
                    else if(key==ambiguity(3))
                    {
                        independent(index,4)=1;
                        if(!keepsReference) independent(index,3)=-1;
                    }
                    else throw std::runtime_error("unexpected product destination coordinate");
                }
                check((state.x-independent*beforeX).cwiseAbs().maxCoeff()<1e-10,"product complete mean independent oracle");
                check((state.P-independent*beforeP*independent.transpose()).cwiseAbs().maxCoeff()<1e-10,"product complete covariance independent oracle");
                if(scenario==1)
                    check(!state.kfIndexMap.contains(phase) && !state.kfIndexMap.contains(ambiguity(2)) &&
                          !state.kfIndexMap.contains(ambiguity(3)),"unknown all products release block");
                else if(scenario==2)
                    check(state.x.size()==5 && (state.x-beforeX).norm()==0 && (state.P-beforeP).norm()==0,
                          "same public segment internal version change preserves complete posterior");
                else if(scenario==4)
                    check(state.kfIndexMap.contains(phase) && state.kfIndexMap.contains(ambiguity(2)) &&
                          !state.kfIndexMap.contains(ambiguity(3)),"unknown nonreference releases only its coordinate");
                else
                    check(state.kfIndexMap.contains(phase) && state.kfIndexMap.contains(ambiguity(3)) &&
                          !state.kfIndexMap.contains(ambiguity(1)) && !state.kfIndexMap.contains(ambiguity(2)),
                          "unknown omitted reference locally retired after transport");
                for(int satellite=1;satellite<=4;++satellite)
                {
                    const auto id=SatSys(E_Sys::GPS,satellite).id();
                    const auto physical=zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(true,"TEST",id,static_cast<int>(code)));
                    check(physical==static_cast<unsigned>(scenario==3 && satellite==1),"product event never fabricates physical break");
                }
                if(scenario==0 || scenario==5 || scenario==6 || scenario==7 || scenario==8 || scenario==9)
                    check(zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(false,"TEST","G01",static_cast<int>(code)))==1,
                        "unknown old product changes state identity exactly once");
                if(scenario==4)
                    check(zhangUserArcGeneration(state.metaDataMap,
                        zhangUserArcGenerationKey(false,"TEST","G03",static_cast<int>(code)))==1,
                        "unknown nonreference changes state identity");
                const auto stableX=state.x;
                const auto stableP=state.P;
                const auto stableMetadata=state.metaDataMap;
                updateZhangPppArUserReferences(trace,receivers,state);
                check(state.x.size()==stableX.size() && (state.x-stableX).norm()==0,"product same event mean idempotence");
                check((state.P-stableP).norm()==0,"product same event full covariance idempotence");
                check(state.metaDataMap==stableMetadata,"product same event receipts generations idempotence");
            }
        }
        // Exercise the installed production factor consumers, not callbacks
        // that merely return true. A real /PPP update anchors their raw replay
        // and authoritative commit chronology before either private KF copy.
        const auto liveFixtureRoot=std::filesystem::current_path()/
            ("user_live_capture_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(liveFixtureRoot);
        for (int scenario=0;scenario<7;++scenario)
        {
            if(isolatedScenario>=0 && scenario!=isolatedScenario) continue;
            if(isolatedScenario<0)
            {
                // Checkpoint export intentionally has one scientific runtime
                // owner per process. Keep that production safety boundary:
                // independent fixture lifetimes run as fresh exec processes,
                // not by clearing the owner guard or inheriting BLAS threads.
                int resultPipe[2];
                if(::pipe(resultPipe)!=0) throw std::runtime_error("live fixture result pipe failed");
                const auto scenarioArgument="--live-capture-scenario="+std::to_string(scenario);
                const auto resultArgument="--scenario-result-fd="+std::to_string(resultPipe[1]);
                const auto process=::fork();
                if(process<0) throw std::runtime_error("live fixture process fork failed");
                if(process==0)
                {
                    ::close(resultPipe[0]);
                    ::execl("/proc/self/exe","user_reference_qc_lifecycle_contract",
                        scenarioArgument.c_str(),resultArgument.c_str(),static_cast<char*>(nullptr));
                    ::_exit(EXIT_FAILURE);
                }
                ::close(resultPipe[1]);
                std::array<unsigned,4> receipt{};
                std::size_t received=0;
                while(received<sizeof(receipt))
                {
                    const auto count=::read(resultPipe[0],reinterpret_cast<char*>(&receipt)+received,
                        sizeof(receipt)-received);
                    if(count<0 && errno==EINTR) continue;
                    if(count<=0) break;
                    received+=count;
                }
                ::close(resultPipe[0]);
                int status=0;
                while(::waitpid(process,&status,0)<0)
                    if(errno!=EINTR) throw std::runtime_error("live fixture process wait failed");
                if(received!=sizeof(receipt) || receipt[0]!=0x52514643U || receipt[1]!=scenario ||
                    receipt[2]!=static_cast<unsigned>(process) || receipt[3]==0 ||
                    !WIFEXITED(status) || WEXITSTATUS(status)!=EXIT_SUCCESS)
                    throw std::runtime_error("isolated real capture scenario failed: "+std::to_string(scenario));
                checks+=receipt[3];
                check(true,"unique executed fresh-process scientific checkpoint owner receipt");
                continue;
            }
            auto completeIsolatedScenario=[&]
            {
                const std::array<unsigned,4> receipt{
                    0x52514643U,static_cast<unsigned>(scenario),static_cast<unsigned>(::getpid()),checks};
                std::cout << "USER_REFERENCE_LIVE_SCENARIO PASS scenario=" << scenario
                          << " pid=" << receipt[2] << " checks=" << checks << '\n' << std::flush;
                if(::write(scenarioResultFd,&receipt,sizeof(receipt))!=sizeof(receipt))
                    throw std::runtime_error("isolated live fixture result write failed");
                ::close(scenarioResultFd);
                ::_exit(EXIT_SUCCESS);
            };
            const bool sharedIf=scenario>=2;
            const bool rollback=scenario==1 || scenario==5;
            acsConfig.zhangPppAr.integer_strategy=sharedIf
                ? "CANONICAL_USER_IF_WL_L1" : "CANONICAL_USER_SD_WL_L1";
            acsConfig.zhangPppAr.baseline_observables[E_Sys::GPS]={E_ObsCode::L1C,E_ObsCode::L2W};
            acsConfig.zhangPppAr.user_use_full_product_covariance=false;
            acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow=false;
            acsConfig.zhangPppAr.targeted_besd_capture_shadow=false;
            acsConfig.pppOpts.queue_rts_outputs=false;
            const auto serial=++scenarioId;
            const auto fixture=liveFixtureRoot/
                ("user_live_capture_datum_"+std::to_string(serial)+".csv");
            const auto covarianceFixture=liveFixtureRoot/
                ("user_live_capture_covariance_"+std::to_string(serial)+".csv");
            std::ofstream products(fixture);
            products << "production CSV52 real hook fixture\n";
            for(int epoch=0;epoch<2;++epoch)
                for(int satellite=1;satellite<=4;++satellite)
                    for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                    {
                        const bool changed=epoch==1 && (satellite==1 ||
                            (scenario==3 && satellite==2 && code==E_ObsCode::L2W) ||
                            (scenario==5 && satellite==3) ||
                            (scenario==6 && code==E_ObsCode::L2W));
                        std::vector<std::string> row(52,"0");
                        row[0]=epoch==0?"100000":"100030";
                        row[1]="PRODUCT_FIXED"; row[2]="G0"+std::to_string(satellite);
                        row[3]=enum_to_string(code);
                        row[5]=".01"; row[7]=".01"; row[10]=".02";
                        row[11]=changed?"1":"0"; row[14]=changed?"2":"1";
                        row[15]=row[0]; row[16]="1"; row[17]="SYNTHETIC_PRODUCT_EVENT";
                        row[18]="1"; row[19]="CURRENT_ALIGNMENT_VALID";
                        for(int field : {20,21,22,23,28,29,30,31,32,42,43,46}) row[field]="1";
                        row[24]="PRODUCT-GPS-C0"; row[25]="PRODUCT-GPS-C0-DATUM";
                        row[26]=row[0]; row[27]=row[0]; row[33]="NONE";
                        row[34]=row[2]+"-"+row[3]+"-PHYSICAL_SEG0";
                        row[35]="1"; row[38]="4"; row[39]="3"; row[40]="3";
                        row[44]="AR_VALID"; row[45]=changed?"1":"0";
                        row[47]="SYNTHETIC_FROZEN_FUNCTIONAL";
                        row[49]=changed?"0":"1"; row[50]=changed?"1":"0";
                        row[51]=row[2]+"-"+row[3]+(changed?"-PUBLIC1":"-PUBLIC0");
                        for(std::size_t field=0;field<row.size();++field)
                            products << (field?",":"") << row[field];
                        products << '\n';
                    }
            products.close();
            std::ofstream covariance(covarianceFixture);
            covariance << "epoch,solution,row_sat,row_type,row_code,column_sat,column_type,column_code,covariance\n";
            for(int epoch : {100000,100030})
                for(int row=0;row<12;++row)
                    for(int column=row;column<12;++column)
                    {
                        const auto code=std::array{E_ObsCode::NONE,E_ObsCode::L1C,E_ObsCode::L2W}[row%3];
                        const auto other=std::array{E_ObsCode::NONE,E_ObsCode::L1C,E_ObsCode::L2W}[column%3];
                        const auto type=code==E_ObsCode::NONE?"CLOCK":"PHASE";
                        const auto otherType=other==E_ObsCode::NONE?"CLOCK":"PHASE";
                        covariance << epoch << ",PRODUCT_FIXED,G0" << row/3+1 << ',' << type << ','
                                   << enum_to_string(code) << ",G0" << column/3+1 << ',' << otherType << ','
                                   << enum_to_string(other) << ',' << (row==column?.0001:0) << '\n';
                    }
            covariance.close();
            acsConfig.zhangPppAr.product_filename=fixture.string();
            acsConfig.zhangPppAr.product_covariance_filename=covarianceFixture.string();
            ReceiverMap receivers;
            auto& receiver=receivers["LIVE"];
            receiver.id="LIVE"; receiver.ready=true;
            for(int satellite=1;satellite<=4;++satellite)
            {
                const SatSys sat(E_Sys::GPS,satellite);
                auto& status=receiver.satStatMap[sat]; status.el=(70-satellite*10)*D2R;
                GObs obs; obs.time.bigTime=100000; obs.Sat=sat; obs.mount=receiver.id;
                obs.satStat_ptr=&status;
                for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                {
                    const auto frequency=code2Freq.at(E_Sys::GPS).at(code);
                    auto& tracking=status.sigStatMap[ft2string(frequency)].tracking;
                    tracking.seen=true; tracking.lastCode=static_cast<int>(code);
                    tracking.lastPhaseTime=100000;
                    obs.sigs[frequency].code=code;
                    obs.sigs[frequency].P=2e7; obs.sigs[frequency].L=1e8;
                }
                receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
            }
            if(scenario==5)
            {
                auto& laterReceiver=receivers["ZZBLOCK"];
                laterReceiver.id="ZZBLOCK"; laterReceiver.ready=true;
                for(const auto& source : only<GObs>(receiver.obsList))
                {
                    auto& status=laterReceiver.satStatMap[source.Sat];
                    status=receiver.satStatMap.at(source.Sat);
                    GObs obs=source; obs.mount=laterReceiver.id; obs.satStat_ptr=&status;
                    laterReceiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
                }
            }
            auto ambiguity=[&](int satellite,E_ObsCode code)
            {
                KFKey key; key.type=KF::AMBIGUITY; key.str=receiver.id;
                key.Sat=SatSys(E_Sys::GPS,satellite);
                key.num=zhangPppArUserPhaseCoordinateNumber(E_Sys::GPS,code);
                return key;
            };
            auto phase=[&](E_ObsCode code)
            {
                KFKey key; key.type=KF::PHASE_BIAS; key.str=receiver.id;
                key.Sat=SatSys(E_Sys::GPS,0);
                key.num=zhangPppArUserPhaseCoordinateNumber(E_Sys::GPS,code);
                return key;
            };
            KFKey clock; clock.type=KF::REC_CLOCK; clock.str=receiver.id;
            KFState state; state.time.bigTime=100000;
            const auto runtimeId="live-reference-capture-"+std::to_string(serial);
            check(bindZhangCheckpointRuntimeId(state,runtimeId),"live hook bind authoritative owner");
            std::vector<KFKey> inputKeys{state.oneKey,clock,phase(E_ObsCode::L1C),
                ambiguity(2,E_ObsCode::L1C),ambiguity(3,E_ObsCode::L1C)};
            if(!sharedIf)
            {
                inputKeys.push_back(phase(E_ObsCode::L2W));
                inputKeys.push_back(ambiguity(2,E_ObsCode::L2W));
                inputKeys.push_back(ambiguity(3,E_ObsCode::L2W));
            }
            if(scenario==5)
            {
                KFKey laterClock=clock; laterClock.str="ZZBLOCK";
                KFKey laterPhase=phase(E_ObsCode::L1C); laterPhase.str="ZZBLOCK";
                KFKey laterSecond=ambiguity(2,E_ObsCode::L1C); laterSecond.str="ZZBLOCK";
                KFKey laterThird=ambiguity(3,E_ObsCode::L1C); laterThird.str="ZZBLOCK";
                inputKeys.insert(inputKeys.end(),{laterClock,laterPhase,laterSecond,laterThird});
            }
            const int dimension=inputKeys.size();
            state.kfIndexMap.clear(); state.stateTransitionMap.clear();
            state.x=VectorXd::Zero(dimension); state.x(0)=1; state.dx=VectorXd::Zero(dimension);
            MatrixXd factor=MatrixXd::Identity(dimension,dimension);
            for(int index=0;index<dimension;++index)
            {
                state.kfIndexMap[inputKeys[index]]=index;
                state.stateTransitionMap[inputKeys[index]][inputKeys[index]][0]=1;
                if(index) state.x(index)=.75+index*1.25;
                for(int other=1;other<dimension;++other)
                    if(index && index!=other) factor(index,other)=.013*(1+(index+other)%5);
            }
            state.P=factor*factor.transpose(); state.P.row(0).setZero(); state.P.col(0).setZero();
            state.prefitOpts.sigma_check=false; state.prefitOpts.omega_test=false;
            state.postfitOpts.sigma_check=false; state.postfitOpts.omega_test=false;
            state.chiSquareTest.enable=false;
            std::ostringstream trace;
            const auto rtsPath=liveFixtureRoot/
                ("user_live_reference_"+std::to_string(serial)+".rts");
            state.rts_basename=rtsPath.string();
            updateZhangPppArUserReferences(trace,receivers,state);
            for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                check(zhangPppArUserReference(state,"LIVE",E_Sys::GPS,code)==SatSys("G01"),"live initial common omitted anchor");
            check(state.rts_basename==rtsPath.string(),"no-event reference update preserves live RTS path");
            state.rts_basename.clear(); // The measurement itself is not the RTS test.
            acsConfig.zhangPppAr.user_use_full_product_covariance=sharedIf;
            acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow=true;
            acsConfig.zhangPppAr.fixed_lag_factor_capture_max_events=2000;
            acsConfig.zhangPppAr.fixed_lag_factor_capture_evaluation_stride=1000;
            configureZhangE18FactorCapture(state);
            check(bool(state.exactStateTransformCallback) && bool(state.authoritativeFactorTransactionCallback),
                "installed real factor shadow consumers");
            KFMeas measurement;
            measurement.time=state.time;
            measurement.H=MatrixXd::Zero(2,dimension);
            measurement.H(0,state.kfIndexMap.at(clock))=.25;
            measurement.H(1,state.kfIndexMap.at(phase(E_ObsCode::L1C)))=.5;
            measurement.R=MatrixXd::Identity(2,2)*.04;
            measurement.H_star=MatrixXd::Identity(2,2);
            measurement.uncorrelatedNoise=VectorXd::Constant(2,.04);
            measurement.Y=measurement.H*state.x;
            measurement.V=VectorXd::Zero(2); measurement.VV=VectorXd::Zero(2);
            measurement.prefitRatios=VectorXd::Zero(2); measurement.postfitRatios=VectorXd::Zero(2);
            measurement.W=VectorXd::Ones(2); measurement.metaDataMaps.resize(2);
            measurement.componentsMaps.resize(2);
            for(int row=0;row<2;++row)
            {
                KFKey key; key.type=row?KF::PHAS_MEAS:KF::CODE_MEAS;
                key.str=receiver.id; key.Sat=SatSys("G02"); key.num=static_cast<int>(E_ObsCode::L1C);
                measurement.obsKeys.push_back(key); measurement.noiseIndexMap[key]=row;
            }
            check(state.filterKalman(trace,measurement,"/PPP",true)==KFFilterResult::COMMITTED,
                "actual PPP measurement anchors real capture chronology");
            std::string anchorPayload;
            const auto anchored=exportZhangPppArCheckpointSection(state,runtimeId,anchorPayload);
            if(!anchored.valid || anchored.capturedFactorEvents!=1)
                std::cerr << "REAL_CAPTURE_CHECKPOINT scenario=" << scenario << " runtime=" << runtimeId
                          << " valid=" << anchored.valid << " events=" << anchored.capturedFactorEvents
                          << " failure_reason=" << anchored.failureReason << '\n';
            check(anchored.valid && anchored.capturedFactorEvents==1,"real anchored factor event checkpoint");
            if(sharedIf) check(anchored.e27SensitivityRows>0,"real IF measurement seeded joint-noise sensitivity sidecar");
            state.rts_basename=rtsPath.string();
            state.time.bigTime=100030;
            for(auto& [id,currentReceiver] : receivers)
                for(auto& obs : only<GObs>(currentReceiver.obsList))
                {
                    obs.time=state.time;
                    for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                        obs.satStat_ptr->sigStatMap[ft2string(code2Freq.at(E_Sys::GPS).at(code))].tracking.lastPhaseTime=100030;
                }
            if(scenario==4)
            {
                auto& signal=receiver.satStatMap[SatSys("G02")].sigStatMap[ft2string(E_FType::F2)];
                signal.tracking.eventSequence=1;
                signal.userQcState=ZhangUserQcState::CONFIRMED_BREAK;
            }
            if(scenario==1) state.procNoiseMap[phase(E_ObsCode::L2W)]=1;
            if(scenario==5)
            {
                // LIVE has already accepted its exact transform, local
                // releases and E27 sensitivity transport when the later
                // receiver hits a genuinely unsupported dynamic coordinate.
                KFKey laterPhase=phase(E_ObsCode::L1C); laterPhase.str="ZZBLOCK";
                state.procNoiseMap[laterPhase]=1;
            }
            if(sharedIf)
            {
                std::vector<double> factors;
                check(queryZhangInternalProductNoiseFactors(state.time,SatSys("G02"),E_ObsCode::L1C,true,factors),
                    "shared IF selector has actual covariance-backed first-band eligibility");
                check(queryZhangInternalProductNoiseFactors(state.time,SatSys("G03"),E_ObsCode::L2W,true,factors),
                    "shared IF selector has actual covariance-backed second-band eligibility");
            }
            const auto beforeX=state.x; const auto beforeP=state.P;
            const auto beforeMetadata=state.metaDataMap;
            const auto beforeSequence=state.factorCommitSequence;
            const auto beforeRts=state.rts_basename;
            std::string beforePayload;
            check(exportZhangPppArCheckpointSection(state,runtimeId,beforePayload).valid,"live before transaction checkpoint");
            std::ostringstream beforeAudit;
            traceZhangE18RawIntegerDatumWindow(beforeAudit,state,state.time);
            check(beforeAudit.str().find("measurement_accepted=1")!=std::string::npos &&
                  beforeAudit.str().find("owner_rejected=0")!=std::string::npos,"real installed hook accepted authoritative measurement");
            bool threw=false;
            const auto traceStart=trace.str().size();
            try {updateZhangPppArUserReferences(trace,receivers,state);}
            catch(const std::runtime_error&) {threw=true;}
            check(threw==rollback,"live production transport or guarded rollback outcome");
            check(state.factorCallbackOwnerOverride==nullptr,"candidate effective owner override never leaks to authoritative object");
            check(state.rts_basename==beforeRts,"live RTS path preserved after transaction or rollback");
            check(bool(state.exactStateTransformCallback) && bool(state.acceptedMeasurementFactorCallback),
                "real factor callbacks remain installed after candidate assignment");
            std::string afterPayload;
            const auto after=exportZhangPppArCheckpointSection(state,runtimeId,afterPayload);
            check(after.valid,"live after transaction checkpoint");
            std::ostringstream afterAudit;
            traceZhangE18RawIntegerDatumWindow(afterAudit,state,state.time);
            if(rollback)
            {
                check(trace.str().substr(traceStart).find("status=APPLIED")!=std::string::npos,
                    "later rejection follows a genuinely accepted production reference transform");
                check((state.x-beforeX).norm()==0 && (state.P-beforeP).norm()==0 && state.metaDataMap==beforeMetadata,
                    "later rejection restores complete authoritative mean covariance metadata");
                check(state.factorCommitSequence==beforeSequence,"later rejection restores authoritative factor sequence");
                check(afterPayload==beforePayload,"real capture and joint-noise sidecar checkpoint byte rollback");
                check(afterAudit.str()==beforeAudit.str(),"real callback counters chronology and rejection audit rollback");
                check(!std::filesystem::exists(beforeRts+FORWARD_SUFFIX),"rolled-back candidate writes no RTS transforms");
                for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                    check(zhangPppArUserReference(state,"LIVE",E_Sys::GPS,code)==SatSys("G01"),"later rejection restores both reference receipts");
                completeIsolatedScenario();
            }
            const int expectedReference=scenario==3 || scenario==4?3:scenario==6?1:2;
            for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
            {
                check(zhangPppArUserReference(state,"LIVE",E_Sys::GPS,code)==SatSys(E_Sys::GPS,expectedReference),
                    "both band getters follow actual coordinate owner reference");
                check(!state.kfIndexMap.contains(ambiguity(expectedReference,code)),"both bands identify actual omitted shared reference");
            }
            MatrixXd independent=MatrixXd::Zero(state.x.size(),dimension);
            for(const auto& [key,index] : state.kfIndexMap)
            {
                if(key==state.oneKey) independent(index,0)=1;
                else if(key==clock) independent(index,1)=1;
                else
                {
                    bool recognised=false;
                    for(const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                    {
                        if(sharedIf && code==E_ObsCode::L2W) continue;
                        const int phaseIndex=code==E_ObsCode::L1C?2:5;
                        const int anchorIndex=expectedReference==3?phaseIndex+2:phaseIndex+1;
                        if(key==phase(code))
                        {
                            independent(index,phaseIndex)=1;
                            independent(index,anchorIndex)=zhangPppArUserPhaseCoordinateWavelength(E_Sys::GPS,code);
                            recognised=true;
                        }
                        else if(key==ambiguity(3,code))
                        {
                            independent(index,phaseIndex+2)=1;
                            independent(index,anchorIndex)-=1;
                            recognised=true;
                        }
                    }
                    if(!recognised) throw std::runtime_error("unexpected live capture transported coordinate");
                }
            }
            check((state.x-independent*beforeX).norm()<1e-10,"real callback lifecycle independent dense complete mean oracle");
            check((state.P-independent*beforeP*independent.transpose()).norm()<1e-10,"real callback lifecycle independent dense full covariance oracle");
            if(scenario==3 || scenario==4 || scenario==6)
                check(!state.kfIndexMap.contains(ambiguity(2,E_ObsCode::L1C)),"other-frequency unknown or physical break cannot inherit shared N");
            if(scenario==6) check(!state.kfIndexMap.contains(phase(E_ObsCode::L1C)),"no joint continuous anchor releases IF block");
            check(afterAudit.str().find("owner_rejected=0")!=std::string::npos &&
                  afterAudit.str().find("transform_rejected=0")!=std::string::npos,"real installed shadow accepts private owner and nested copy transforms");
            check(after.capturedFactorEvents==1+state.factorCommitSequence-beforeSequence,
                "accepted transform and local release capture chronology matches live factor sequence");
            check(std::filesystem::exists(beforeRts+FORWARD_SUFFIX) && std::filesystem::file_size(beforeRts+FORWARD_SUFFIX)>0,
                "committed reference transforms written to original RTS path");
            const auto stablePayload=afterPayload; const auto stableSequence=state.factorCommitSequence;
            const auto stableSize=std::filesystem::file_size(beforeRts+FORWARD_SUFFIX);
            updateZhangPppArUserReferences(trace,receivers,state);
            check(state.factorCommitSequence==stableSequence,"duplicate datum event does not emit another factor");
            check(exportZhangPppArCheckpointSection(state,runtimeId,afterPayload).valid && afterPayload==stablePayload,
                "duplicate datum event preserves real capture sidecar and receipts");
            check(std::filesystem::file_size(beforeRts+FORWARD_SUFFIX)==stableSize,"duplicate event emits no RTS transform");
            if(scenario==0)
            {
                KFMeasurementTransaction activePacket;
                state.activeMeasurementTransaction=&activePacket;
                const auto activeX=state.x; const auto activeP=state.P;
                bool activeRejected=false;
                try {updateZhangPppArUserReferences(trace,receivers,state);}
                catch(const std::runtime_error& error)
                {activeRejected=std::string(error.what())=="USER_REFERENCE_ACTIVE_MEASUREMENT_TRANSACTION";}
                check(activeRejected,"active packet misuse stops current epoch before new-product consumption");
                check(state.activeMeasurementTransaction==&activePacket && state.rts_basename==beforeRts,
                    "reference hook refuses to detach active authoritative measurement packet");
                check((state.x-activeX).norm()==0 && (state.P-activeP).norm()==0 && state.factorCommitSequence==stableSequence,
                    "active packet rejection leaves complete state and sequence unchanged");
                state.activeMeasurementTransaction=nullptr;
            }
            if(scenario==4)
            {
                check(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,"LIVE","G02",static_cast<int>(E_ObsCode::L2W)))==1,
                    "shared IF physical second-frequency break increments its real physical identity");
                check(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,"LIVE","G02",static_cast<int>(E_ObsCode::L1C)))==0,
                    "shared IF second-frequency break does not fabricate first-frequency slip");
            }
            completeIsolatedScenario();
        }
        std::cout << "USER_REFERENCE_QC_LIFECYCLE_CONTRACT PASS checks=" << checks
                  << " mean_full_covariance_independent_oracle=1 NIS_computed=0"
                  << " product_event_split_checked=1 unknown_datum_relations_released=1"
                  << " transport_full_covariance_oracle=1 same_event_idempotence_checked=1"
                  << " malformed_public_metadata_fail_closed=1 reference_transaction_rollback_checked=1"
                  << " capture_callback_owner_checked=1 capture_sidecars_rollback_checked=1"
                  << " shared_if_reference_consistent_checked=1\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
