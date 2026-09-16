// Linked against the production PEA objects by the accompanying audit script.
// No KFState stubs: these checks exercise algebra.cpp's real transaction path.
#include "common/algebra.hpp"
#include "common/acsConfig.hpp"
#include "common/constants.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/receiver.hpp"
#include "common/observations.hpp"
#include "pea/zhangReference.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

static void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main()
{
    std::ostringstream trace;
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
    std::cout<<"PASS production algebra: pending retirement, marginal covariance, fresh correlated priors, "
        "callback routing, factor sequence, rejection rollback, exact-only rejection, "
        "retired-chord audit, retired-tree rejection before commit, converged-generation skip, "
        "last-iteration reconciliation, production controller all-arc restart\n";
    return 0;
}
