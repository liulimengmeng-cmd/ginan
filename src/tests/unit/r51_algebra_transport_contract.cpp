// Linked against the production PEA objects by the accompanying audit script.
// No KFState stubs: these checks exercise algebra.cpp's real transaction path.
#include "common/algebra.hpp"
#include "common/acsConfig.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/receiver.hpp"
#include "common/observations.hpp"
#include "pea/zhangReference.hpp"
#include "ambres/GNSSambres.hpp"
#include "pea/zhangReceiverCheckpoint.hpp"
void assignObservationValue(RawSig&,char,double,double);
#include <iostream>
#include <sstream>
#include <stdexcept>

static void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main()
{
    std::ostringstream trace;
    for(unsigned lli : {1u,2u,3u,4u,7u}) {
        RawSig raw;
        assignObservationValue(raw,'L',12345,lli);
        check(raw.LLI==lli,"RINEX assignment collapsed LLI bit mask");
    }
    for(auto mode : {E_ARmode::LAMBDA,E_ARmode::LAMBDA_ALT}) {
        for(double value : {0.,.01,.49,.5}) {
            GinAR_mtx matrix;matrix.aflt=VectorXd::Constant(1,value);
            matrix.Paflt=MatrixXd::Constant(1,1,.12*.12);
            GinAR_opt opt;opt.mode=mode;opt.min_lambda_fix_count=1;
            opt.lambda_candidate_nis_alpha=1e-6;
            const int count=GNSS_AR(trace,matrix,opt);
            check(matrix.searchDiagnostic.ilsComplete,"production ILS did not complete top two");
            check(matrix.searchDiagnostic.ratioExecuted,"production LAMBDA skipped ratio");
            check(matrix.searchDiagnostic.localNisExecuted,"production LAMBDA skipped local NIS");
            check((count==1)==(value<.1),"production ratio accepted ambiguous or rejected clear candidate");
        }
    }

    setenv("ZHANG_R51_RATIO_ONLY","1",1);
    for(auto mode : {E_ARmode::LAMBDA,E_ARmode::LAMBDA_ALT}) {
        for(auto test : {std::pair<double,double>{.1,1.0},{.1,.001},{.36,1.0},{.37,1.0},{.49,1.0},{.5,1.0}}) {
            GinAR_mtx matrix;matrix.aflt=VectorXd::Constant(1,test.first);
            matrix.Paflt=MatrixXd::Constant(1,1,test.second*test.second);
            GinAR_opt opt;opt.mode=mode;opt.min_lambda_fix_count=1;opt.sucthr=.9999;opt.ratthr=3;
            opt.lambda_candidate_nis_alpha=1e-6;
            int count=GNSS_AR(trace,matrix,opt);
            check(matrix.searchDiagnostic.ilsComplete,"ratio-only accepted incomplete ILS");
            check(matrix.searchDiagnostic.ratioExecuted,"ratio-only missed ratio gate");
            check((count==1)==(test.first<.365),"ratio-only failed bootstrap/NIS ablation or accepted ambiguous ratio");
            if(test.second==.001)check(matrix.lambda_candidate_nis>matrix.lambda_candidate_nis_threshold,"NIS-ablation fixture was not rejecting");
        }
    }
    // Exercise the real search twice; diagnostic eigensystems must not change
    // candidate enumeration, ratio, selected integer rows, or fixed values.
    for(auto mode : {E_ARmode::LAMBDA,E_ARmode::LAMBDA_ALT}) {
        for(double fraction : {.01,.36,.49}) {
            GinAR_mtx on,off;
            on.aflt=VectorXd::Constant(8,fraction);
            on.Paflt=MatrixXd::Identity(8,8)*.04+MatrixXd::Constant(8,8,.002);
            off=on;
            GinAR_opt opt;opt.mode=mode;opt.min_lambda_fix_count=8;opt.ratthr=3;
            opt.lambda_candidate_nis_alpha=1e-6;
            setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS","1",1);
            int onCount=GNSS_AR(trace,on,opt);
            setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS","0",1);
            int offCount=GNSS_AR(trace,off,opt);
            check(onCount==offCount,"diagnostic switch changed fixed count");
            check(on.Ztrs.rows()==off.Ztrs.rows() && (on.Ztrs-off.Ztrs).norm()==0,
                "diagnostic switch changed integer rows");
            check(on.zfix.size()==off.zfix.size() && (on.zfix-off.zfix).norm()==0,
                "diagnostic switch changed fixed integers");
            check(on.searchDiagnostic.ratio==off.searchDiagnostic.ratio,
                "diagnostic switch changed ratio");
            check(on.searchDiagnostic.ilsCalls==off.searchDiagnostic.ilsCalls,
                "diagnostic switch changed enumeration calls");
            check(std::isnan(off.lambda_candidate_min_sigma) && std::isnan(off.lambda_whitened_condition_number),
                "uncomputed diagnostics were reported as zeros");
            check(on.lambda_dominant_whitened_mode>=0,
                "diagnostic-on fixture did not exercise eigensystem");
            check(off.lambda_dominant_whitened_mode==-1 && off.lambda_dominant_original_loading.size()==0,
                "diagnostic-off still calculated whitened spectrum");
        }
    }
    unsetenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS");
    unsetenv("ZHANG_R51_RATIO_ONLY");

    KFKey a,b,z; a.type=KF::IONO_STEC; a.str="A";
    b.type=KF::IONO_STEC; b.str="B"; z.type=KF::AMBIGUITY; z.str="fresh";
    KFState source; source.kfIndexMap={{a,0},{b,1}};
    source.stateTransitionMap.clear();
    source.stateTransitionMap[a][a][0]=1; source.stateTransitionMap[b][b][0]=1;
    source.x=VectorXd(2); source.x<<2,3; source.dx=VectorXd::Zero(2);
    source.P=MatrixXd(2,2); source.P<<4,1,1,9;
    auto removed=source;
    removed.removeState(b);
    check(removed.kfIndexMap.contains(b),"fixture must still index pending removal");
    check(!removed.stateTransitionMap.contains(b),"removeState did not schedule retirement");
    check(removed.applyStateTransform(trace,{{a,{{a,1}}},{b,{{b,1}}}},"pending"),"identity transform failed");
    check(!removed.stateTransitionMap.contains(b),"transform resurrected retired STEC");
    check(removed.stateTransition(trace,removed.time),"pending marginalization failed");
    check(!removed.kfIndexMap.contains(b),"expired STEC persisted at next transition");
    check(std::abs(removed.P(removed.kfIndexMap.at(a),removed.kfIndexMap.at(a))-4)<1e-12,
        "retirement conditioned instead of marginalized surviving covariance");

    auto dynamic=source;
    check(dynamic.setProcessModel(a,.5,123.,600.),"explicit process model setter rejected valid model");
    check(dynamic.setProcessModel(a,0.,0.,-1.),"explicit zero model setter failed");
    check(dynamic.procNoiseMap.at(a)==0 && dynamic.gaussMarkovMuMap.at(a)==0 &&
        dynamic.gaussMarkovTauMap.at(a)==-1,"explicit zeros did not replace prior GM metadata");
    dynamic.stateTransitionMap[a][b][1]=1;
    check(dynamic.applyStateTransform(trace,{{a,{{a,1}}},{b,{{b,1}}}},"untouched_clock_rate"),
        "identity transport rejected unrelated clock/rate dynamics");
    check(dynamic.stateTransitionMap.at(a).at(b).at(1)==1,"transport erased clock/rate transition");
    auto beforeDynamic=dynamic;
    check(!dynamic.applyStateTransform(trace,{{a,{{a,1},{b,1}}},{b,{{b,1}}}},"coupled_reject"),
        "mixed dynamic coordinate was accepted without F/Q transport");
    check((dynamic.P-beforeDynamic.P).norm()==0,"dynamic rejection modified covariance");
    dynamic=source;dynamic.setProcessModel(a,.01,0.,-1.);
    check(!dynamic.applyStateTransform(trace,{{a,{{a,1},{b,1}}},{b,{{b,1}}}},"phase_q_reject"),
        "nonzero phase Q accepted by exact transport");

    auto state=source;
    bool exactCalled=false,transitionCalled=false;
    state.exactStateTransformCallback=[&](auto const&...){exactCalled=true;return true;};
    state.stateTransitionFactorCallback=[&](const KFState&,GTime,const auto&,const auto&,
        const auto&,const MatrixXd& q,const std::string&,const VectorXd& x,const MatrixXd& p,
        std::uint64_t,std::uint64_t) {
        transitionCalled=true;
        check(std::abs(q(0,1)-4)<1e-12,"fresh shared source lost correlated Q");
        check(std::abs(x(0)-5)<1e-12,"transformed mean is incorrect");
        check(std::abs(p(0,0)-17)<1e-12,"TPT plus Q is incorrect");
        return true;
    };
    std::map<KFKey,std::map<KFKey,double>> transform{{a,{{a,1},{b,1},{z,1}}},{b,{{b,1},{z,2}}}};
    check(state.applyStateTransform(trace,transform,"fresh_arc",{{z,2}}),"stochastic transform failed");
    check(transitionCalled&&!exactCalled,"fresh arc was published as an exact transform");
    check(std::abs(state.P(0,1)-14)<1e-12,"full transformed cross covariance is incorrect");
    check(state.factorCommitSequence==source.factorCommitSequence+1,"factor sequence not advanced once");

    auto rejected=source;
    rejected.stateTransitionFactorCallback=[](auto const&...){return false;};
    check(!rejected.applyStateTransform(trace,transform,"reject",{{z,2}}),"rejected factor committed");
    check((rejected.x-source.x).norm()==0 && (rejected.P-source.P).norm()==0,
        "factor rejection modified posterior");
    check(rejected.kfIndexMap==source.kfIndexMap && rejected.stateTransitionMap==source.stateTransitionMap,
        "factor rejection modified lifecycle metadata");
    check(rejected.factorCommitSequence==source.factorCommitSequence,"rejected sequence advanced");
    auto unsupported=source;
    unsupported.exactStateTransformCallback=[](auto const&...){return true;};
    check(!unsupported.applyStateTransform(trace,transform,"no_transition_consumer",{{z,2}}),
        "stochastic transform silently used exact-only consumer");
    const ZhangGraphEdge a1{"A",SatSys("G01")},a2{"A",SatSys("G02")},
        b1{"B",SatSys("G01")},b2{"B",SatSys("G02")};
    const auto graph=zhangBuildSpanningTree({a1,a2,b1,b2},"A",{a1,a2,b1});
    KFState phase; phase.kfIndexMap.clear(); phase.stateTransitionMap.clear();
    for(const auto& key : {zhangTransportReceiverKey(E_Sys::GPS,E_ObsCode::L1C,"B"),
        zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G01")),
        zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G02")),
        zhangTransportAmbiguityKey(E_ObsCode::L1C,b2)}) {
        phase.kfIndexMap[key]=phase.kfIndexMap.size();phase.stateTransitionMap[key][key][0]=1;
    }
    phase.x=VectorXd::LinSpaced(4,1,4);phase.dx=VectorXd::Zero(4);
    phase.P=MatrixXd::Identity(4,4);phase.P(0,1)=phase.P(1,0)=0.3;
    const auto original=phase;std::string reason;
    check(applyZhangGraphRetirementForAudit(trace,phase,E_Sys::GPS,{E_ObsCode::L1C},
        graph,graph,{b2},reason),"production audit demanded a deliberately retired chord");
    check(phase.x.size()==3,"retired chord was not projected out");
    auto invalid=original;
    check(!applyZhangGraphRetirementForAudit(trace,invalid,E_Sys::GPS,{E_ObsCode::L1C},
        graph,graph,{a1},reason),"retired tree arc was accepted as replacement support");
    check((invalid.x-original.x).norm()==0 && (invalid.P-original.P).norm()==0,
        "invalid graph changed state before rejection");
    for (int iterations : {1,2})
    {
        KFState scalar; scalar.kfIndexMap={{a,0}};
        scalar.stateTransitionMap.clear();scalar.stateTransitionMap[a][a][0]=1;
        scalar.initNoiseMap.clear();
        scalar.x=VectorXd::Zero(1);scalar.dx=VectorXd::Zero(1);scalar.P=MatrixXd::Identity(1,1);
        scalar.prefitOpts.sigma_check=false;scalar.prefitOpts.omega_test=false;
        scalar.postfitOpts.sigma_check=true;scalar.postfitOpts.omega_test=false;
        scalar.postfitOpts.state_sigma_threshold=1e9;scalar.postfitOpts.meas_sigma_threshold=1;
        scalar.postfitOpts.max_iterations=iterations;
        scalar.measRejectCallbacks.push_back(+[](RejectCallbackDetails detail) {
            detail.kfMeas.R(detail.measIndex,detail.measIndex)=1e6;
            return true;
        });
        KFMeas measurement;
        measurement.Y=VectorXd::Constant(1,100);measurement.H=MatrixXd::Identity(1,1);
        measurement.R=MatrixXd::Identity(1,1);measurement.V=VectorXd::Zero(1);
        measurement.VV=VectorXd::Zero(1);measurement.prefitRatios=VectorXd::Zero(1);
        measurement.postfitRatios=VectorXd::Zero(1);measurement.obsKeys={a};
        measurement.metaDataMaps.resize(1);measurement.componentsMaps.resize(1);
        std::ostringstream transactionTrace;
        check(scalar.filterKalman(transactionTrace,measurement,"GENERATION_TEST")==KFFilterResult::COMMITTED,
            "production generation regression did not commit");
        check(std::abs(scalar.x(0)-100./1000001)<1e-12,"final-noise posterior mean wrong");
        check(std::abs(scalar.P(0,0)-1000000./1000001)<1e-12,"final-noise posterior covariance wrong");
        const bool reconciled=transactionTrace.str().find("TERMINAL_RECONCILIATION_BEGIN")!=std::string::npos;
        check(reconciled==(iterations==1),"terminal solve was redundant or necessary reconciliation omitted");
    }
    // Exercise the production graph controller, not just the coordinate planner:
    // every represented arc breaks while valid new carrier observations exist.
    acsConfig.zhangFullRank.enable=true;acsConfig.process_sys[E_Sys::GPS]=true;
    acsConfig.exclude.LLI=false;acsConfig.ambErrors.resetOnSlip.LLI=true;
    acsConfig.zhangPppAr.tree_slip_pivot_before_retire=true;
    acsConfig.zhangPppAr.tree_slip_shadow_replay=false;
    auto& options=acsConfig.zhangFullRank.sysOpts[E_Sys::GPS];
    options.use_spanning_tree=true;options.reference_receiver="A";
    options.reference_satellite="G01";options.baseline_observables={E_ObsCode::L1C};
    options.product_core_min_satellite_support=0;
    ReceiverMap receivers;
    for (const auto& edge:graph.edges) {
        auto& rec=receivers[edge.receiver];rec.id=edge.receiver;rec.ready=true;
        auto observation=std::make_shared<GObs>();observation->Sat=edge.satellite;
        auto& stat=rec.satStatMap[edge.satellite];stat.el=1.5;observation->satStat_ptr=&stat;
        const auto frequency=code2Freq[E_Sys::GPS][E_ObsCode::L1C];
        auto& signal=observation->sigs[frequency];signal.code=E_ObsCode::L1C;
        signal.L=1e8;signal.P=2e7;rec.obsList.push_back(observation);
        auto& config=acsConfig.getRecOpts(edge.receiver,{"GPS","L1C"});
        config.ambiguity.sigma={1000};config.phase_bias.sigma={1000};
    }
    KFState network;
    network.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="R51_ALL_ARCS_RETIRED_TEST";
    updateZhangFullRankReferences(trace,receivers,network);
    ZhangGraphIntegerContext context;
    check(zhangGraphIntegerContext(network,E_Sys::GPS,context),"synthetic initial graph was not initialized");
    check(context.basis.treeEdges==graph.treeEdges,"synthetic KF chart does not match initialized graph");
    network=original;
    network.time=network.time+30;
    network.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="R51_ALL_ARCS_RETIRED_TEST";
    for(auto& [name,rec]:receivers) for(auto& [sat,stat]:rec.satStatMap) {
        auto& signal=stat.sigStatMap[ft2string(code2Freq[E_Sys::GPS][E_ObsCode::L1C])];
        signal.slip.LLI=true;signal.savedSlip.LLI=true;
    }
    updateZhangFullRankReferences(trace,receivers,network);
    for(const auto& edge:graph.edges)
        check(zhangGraphModelsObservation(network,edge.receiver,edge.satellite,E_ObsCode::L1C),
            "all-arcs-retired controller deadlocked with valid new observations");
    for(auto& [name,rec]:receivers) for(auto& [sat,stat]:rec.satStatMap)
        check(!stat.sigStatMap[ft2string(code2Freq[E_Sys::GPS][E_ObsCode::L1C])].savedSlip.any,
            "committed all-arc retirement left saved slip pending");
    check(trace.str().find("ZHANG_COMPONENT_TRANSPORT_COMMIT")!=std::string::npos,
        "all-arc restart did not commit component transport");
    eraseZhangGraphRuntime(network);
    auto initialiseNetwork=[&](const std::string& id) {
        for(auto& [name,rec]:receivers) for(auto& obs:only<GObs>(rec.obsList)) {
            obs.excludeElevation=false;
            auto& stat=*obs.satStat_ptr;
            stat.sigStatMap.clear();
        }
        KFState fresh;fresh.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=id;
        updateZhangFullRankReferences(trace,receivers,fresh);
        fresh=original;fresh.time=fresh.time+30;
        fresh.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=id;
        return fresh;
    };
    // Rejection must quarantine a replacement arc even when its endpoints
    // remain in the previous chart and the old tree could still model a row.
    auto failed=initialiseNetwork("R51_ADMISSION_REJECTION");
    auto failedBefore=failed;
    failed.exactStateTransformCallback=[](auto const&...){return false;};
    failed.stateTransitionFactorCallback=[](auto const&...){return false;};
    for(auto& [name,rec]:receivers) for(auto& [sat,stat]:rec.satStatMap) {
        auto& signal=stat.sigStatMap[ft2string(code2Freq[E_Sys::GPS][E_ObsCode::L1C])];
        signal.slip.LLI=true;signal.savedSlip.LLI=true;signal.tracking.pendingBreak=true;
    }
    updateZhangFullRankReferences(trace,receivers,failed);
    for(const auto& edge:graph.edges)
        check(!zhangGraphModelsObservation(failed,edge.receiver,edge.satellite,E_ObsCode::L1C),
            "rejected new arc was admitted through old endpoints");
    check((failed.x-failedBefore.x).norm()==0 && (failed.P-failedBefore.P).norm()==0,
        "rejected controller changed numerical posterior");
    for(auto& [name,rec]:receivers) for(auto& [sat,stat]:rec.satStatMap)
        check(stat.sigStatMap[ft2string(code2Freq[E_Sys::GPS][E_ObsCode::L1C])].tracking.pendingBreak,
            "rejected retirement consumed a pending event");
    eraseZhangGraphRuntime(failed);
    auto productFailure=initialiseNetwork("R51_PRODUCT_REJECTION_FLOAT_PRESERVED");
    options.product_core_min_satellite_support=1000;
    options.product_integer_support_core=true;
    for(auto& [name,rec]:receivers) for(auto& [sat,stat]:rec.satStatMap) {
        auto& status=stat.sigStatMap[ft2string(code2Freq[E_Sys::GPS][E_ObsCode::L1C])];
        status.slip.LLI=true;status.savedSlip.LLI=true;status.tracking.pendingBreak=true;
    }
    updateZhangFullRankReferences(trace,receivers,productFailure);
    for(const auto& edge:graph.edges) {
        check(zhangGraphModelsObservation(productFailure,edge.receiver,edge.satellite,E_ObsCode::L1C),
            "product rejection rolled back an accepted FLOAT arc");
        check(!zhangGraphProductSatelliteActive(productFailure,edge.satellite),
            "failed product transaction retained current authority");
    }
    std::string rejectedProductPayload,productWhy;
    check(exportZhangGraphCheckpointSection(productFailure,"R51_PRODUCT_REJECTION_FLOAT_PRESERVED",
        rejectedProductPayload,productWhy),productWhy.c_str());
    productFailure.time=productFailure.time+30;
    check(!zhangGraphModelsObservation(productFailure,"A",SatSys("G01"),E_ObsCode::L1C),
        "previous-epoch phase admission remained usable");
    eraseZhangGraphRuntime(productFailure);
    options.product_core_min_satellite_support=0;
    for(bool removeRoot : {false,true}) {
        const std::string id=removeRoot?"R51_ROOT_ABSENT_RESTORE":"R51_SPLIT_RESTORE";
        auto forest=initialiseNetwork(id);
        for(auto& [name,rec]:receivers) for(auto& obs:only<GObs>(rec.obsList)) {
            bool keep=(!removeRoot && name=="A" && obs.Sat==SatSys("G01")) ||
                (name=="B" && obs.Sat==SatSys("G02"));
            if(!keep) {
                obs.excludeElevation=true;
                auto& status=obs.satStat_ptr->sigStatMap[ft2string(code2Freq[E_Sys::GPS][E_ObsCode::L1C])];
                status.slip.LLI=true;status.savedSlip.LLI=true;status.tracking.pendingBreak=true;
            }
        }
        updateZhangFullRankReferences(trace,receivers,forest);
        ZhangGraphIntegerContext unauthorized;
        check(!zhangGraphIntegerContext(forest,E_Sys::GPS,unauthorized),
            "local FLOAT chart incorrectly authorized a global integer search");
        std::string payload,why;
        check(exportZhangGraphCheckpointSection(forest,id,payload,why),why.c_str());
        check(validateZhangGraphCheckpointSection(id,payload,why),why.c_str());
        KFState restored=forest;
        eraseZhangGraphRuntime(forest);
        check(importZhangGraphCheckpointSection(restored,id,payload,why),why.c_str());
        check(zhangGraphModelsObservation(restored,"B",SatSys("G02"),E_ObsCode::L1C),
            "surviving local FLOAT component was not restored");
        check(!zhangGraphProductSatelliteActive(restored,SatSys("G02")),
            "local FLOAT forest was promoted to current product authority");
        auto later=restored.time+30;
        check(restored.stateTransition(trace,later),"restored forest positive-time propagation failed");
        updateZhangFullRankReferences(trace,receivers,restored);
        check(restored.x.allFinite() && restored.P.allFinite(),"restored update produced nonfinite posterior");
        check(zhangGraphModelsObservation(restored,"B",SatSys("G02"),E_ObsCode::L1C),
            "restored component lost observation admission after positive-time update");
        KFMeas row;
        const auto rk=zhangTransportReceiverKey(E_Sys::GPS,E_ObsCode::L1C,"B");
        const auto sk=zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G02"));
        row.H=MatrixXd::Zero(1,restored.x.size());
        row.H(0,restored.kfIndexMap.at(rk))=1;row.H(0,restored.kfIndexMap.at(sk))=1;
        row.Y=row.H*restored.x+VectorXd::Constant(1,.01);row.R=MatrixXd::Constant(1,1,.001);
        row.V=VectorXd::Zero(1);row.VV=VectorXd::Zero(1);
        row.prefitRatios=VectorXd::Zero(1);row.postfitRatios=VectorXd::Zero(1);
        row.obsKeys={sk};row.metaDataMaps.resize(1);row.componentsMaps.resize(1);
        restored.prefitOpts.sigma_check=false;restored.prefitOpts.omega_test=false;
        restored.postfitOpts.sigma_check=false;restored.postfitOpts.omega_test=false;
        const double varianceBefore=(row.H*restored.P*row.H.transpose())(0,0);
        check(restored.filterKalman(trace,row,"RESTORED_COMPONENT_UPDATE")==KFFilterResult::COMMITTED,
            "restored component's real KF observation update failed");
        check((row.H*restored.P*row.H.transpose())(0,0)<varianceBefore,
            "restored local physical observation did not reduce uncertainty");
        eraseZhangGraphRuntime(restored);
    }
    std::cout<<trace.str()<<"\n";
    std::cout<<"PASS production algebra: pending retirement, marginal covariance, fresh correlated priors, "
        "callback routing, factor sequence, rejection rollback, exact-only rejection, "
        "retired-chord audit, retired-tree rejection before commit, converged-generation skip, "
        "last-iteration reconciliation, production controller all-arc restart, explicit zero process model, "
        "dynamic F/Q guards, LLI bit preservation, production LAMBDA ratio and NIS, failed arc admission, "
        "product invalidation without FLOAT rollback, split and root-absent checkpoint restore/update\n";
    return 0;
}
