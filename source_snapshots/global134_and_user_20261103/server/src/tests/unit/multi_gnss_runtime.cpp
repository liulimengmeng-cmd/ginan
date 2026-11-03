// Links the actual PEA graph/datum implementation, not a reimplementation.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/navigation.hpp"
#include "common/rinex.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangCurrentBridgeWorkIdentity.hpp"
#include "common/zhangProductConditioningReceipt.hpp"
#include "common/zhangSmallProductResponse.hpp"
#include "common/zhangProductIntegerLedger.hpp"
#include "common/zhangSystemClockDatum.hpp"
#include "ambres/GNSSambres.hpp"
#include "pea/zhangPppAr.hpp"
#include "pea/zhangReference.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace
{
void check(bool pass,const std::string& why)
{if(!pass)throw std::runtime_error(why);}

struct System
{
    E_Sys sys;
    const char* first;
    const char* second;
    std::vector<E_ObsCode> codes;
};
const std::vector<System> systems{
    {E_Sys::GPS,"G01","G02",{E_ObsCode::L1C,E_ObsCode::L2W}},
    {E_Sys::GAL,"E11","E12",{E_ObsCode::L1C,E_ObsCode::L5Q}},
    {E_Sys::BDS,"C19","C20",{E_ObsCode::L1P,E_ObsCode::L5P}}};

void configure(bool graph)
{
    setenv("ZHANG_R51_ENABLE","1",1);
    setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
    setenv("ZHANG_R51_PRODUCT_REPAIR","1",1);
    setenv("ZHANG_R51_RATIO_ONLY","1",1);
    setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);
    setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);
    acsConfig.zhangFullRank.enable=true;
    for(auto& [system,enabled]:acsConfig.process_sys)enabled=false;
    acsConfig.zhangFullRank.sysOpts.clear();
    acsConfig.zhangFullRank.bds3_only=true;
    acsConfig.zhangFullRank.bds3_satellite_whitelist={"C19","C20"};
    acsConfig.zhangPppAr.tree_slip_pivot_before_retire=true;
    acsConfig.zhangPppAr.tree_slip_shadow_replay=false;
    acsConfig.ambErrors.resetOnSlip.LLI=true;
    acsConfig.exclude.LLI=true;
    for(const auto& system:systems)
    {
        acsConfig.process_sys[system.sys]=true;
        auto& options=acsConfig.zhangFullRank.sysOpts[system.sys];
        options.use_spanning_tree=graph;
        options.auto_reference_switch=!graph;
        options.reference_receiver="A";
        options.reference_satellite=system.first;
        options.reference_receiver_candidates={"B","A"};
        options.reference_outage_epochs=1;
        options.baseline_observables=system.codes;
        options.product_core_min_satellite_support=0;
    }
}

void observations(ReceiverMap& receivers)
{
    for(const auto& system:systems)
    for(const std::string receiverId:{"A","B"})
    for(const char* id:{system.first,system.second})
    {
        auto& receiver=receivers[receiverId];
        receiver.id=receiverId;receiver.ready=true;
        auto obs=std::make_shared<GObs>();obs->Sat=SatSys(id);
        auto& status=receiver.satStatMap[obs->Sat];status.el=1.3;
        obs->satStat_ptr=&status;
        for(const auto code:system.codes)
        {
            const auto frequency=code2Freq.at(system.sys).at(code);
            auto& signal=obs->sigs[frequency];
            signal.code=code;signal.L=1e8;signal.P=2e7;
            auto& config=acsConfig.getRecOpts(receiverId,{enum_to_string(system.sys),enum_to_string(code)});
            config.ambiguity.sigma={1000};config.phase_bias.sigma={1000};
        }
        receiver.obsList.push_back(obs);
    }
}

void add(KFState& state,const KFKey& key)
{
    InitialState initial;initial.estimate=true;initial.x=2;initial.P=1;initial.Q=0;
    state.addKFState(key,initial);
}

ZhangGraphBasis star(const System& system,const std::string& root="A")
{
    ZhangGraphBasis basis;basis.rootReceiver=root;
    basis.receivers={"A","B"};basis.satellites={SatSys(system.first),SatSys(system.second)};
    for(const std::string receiver:{"A","B"})
    for(const auto& satellite:basis.satellites)
    {
        const ZhangGraphEdge edge{receiver,satellite};basis.edges.insert(edge);
        if(receiver==root || satellite==SatSys(system.first))basis.treeEdges.insert(edge);
    }
    basis.connected=true;basis.componentCount=1;return basis;
}

void states(KFState& state,std::ostringstream& trace,bool graph)
{
    for(const auto& system:systems)
    {
        ZhangGraphIntegerContext context;
        const auto basis=graph
            ? (check(zhangGraphIntegerContext(state,system.sys,context),"production graph missing"),context.basis)
            : star(system);
        check(basis.edges.size()==4 && basis.treeEdges.size()==3,"not a real 2-station/2-satellite graph");
        for(const auto code:system.codes)
        {
            for(const auto& receiver:basis.receivers)
                if(receiver!=basis.rootReceiver)add(state,zhangTransportReceiverKey(system.sys,code,receiver));
            for(const auto& satellite:basis.satellites)add(state,zhangTransportSatelliteKey(code,satellite));
            for(const auto& edge:basis.edges)
                if(!basis.treeEdges.contains(edge))add(state,zhangTransportAmbiguityKey(code,edge));
        }
        for(const auto type:{KF::REC_CLOCK,KF::REC_CLOCK_RATE})
            add(state,zhangReceiverClockKey(acsConfig,"B",system.sys,0,type));
        for(const auto& satellite:basis.satellites)
        for(const auto type:{KF::SAT_CLOCK,KF::SAT_CLOCK_RATE})
        {KFKey key;key.type=type;key.Sat=satellite;add(state,key);}
    }
    check(state.stateTransition(trace,state.time),"state creation failed");
    std::mt19937 random(20261021);std::normal_distribution<double> normal;
    MatrixXd seed(state.x.size(),state.x.size());
    for(int r=0;r<seed.rows();++r)
    for(int c=0;c<seed.cols();++c)seed(r,c)=normal(random);
    state.P=seed*seed.transpose()+MatrixXd::Identity(seed.rows(),seed.cols());
    for(const auto& [key,index]:state.kfIndexMap)
    {
        state.x(index)=0.2*index+1;
        if(key.type==KF::ONE)
        {state.x(index)=1;state.P.row(index).setZero();state.P.col(index).setZero();}
    }
    state.dx=VectorXd::Zero(state.x.size());
}

// Independent physical observation equation: r_phase+s_phase+lambda*N.
// This reference does not call graph coordinate expansion/transport routines.
MatrixXd phaseH(const KFState& state,const ZhangGraphBasis& basis,E_Sys system,E_ObsCode code)
{
    MatrixXd h=MatrixXd::Zero(basis.edges.size(),state.x.size());
    const double lambda=genericWavelength.at(code2Freq.at(system).at(code));
    int row=0;
    for(const auto& edge:basis.edges)
    {
        auto coefficient=[&](const KFKey& key,double value)
        {const auto found=state.kfIndexMap.find(key);check(found!=state.kfIndexMap.end(),"independent H state missing");h(row,found->second)=value;};
        if(edge.receiver!=basis.rootReceiver)
        {KFKey key;key.type=KF::PHASE_BIAS;key.Sat=SatSys(system,0);key.str=edge.receiver;key.num=int(code);coefficient(key,1);}
        KFKey sat;sat.type=KF::PHASE_BIAS;sat.Sat=edge.satellite;sat.num=int(code);coefficient(sat,1);
        if(!basis.treeEdges.contains(edge))
        {KFKey amb;amb.type=KF::AMBIGUITY;amb.Sat=edge.satellite;amb.str=edge.receiver;amb.num=int(code);coefficient(amb,lambda);}
        ++row;
    }
    return h;
}

MatrixXd independentTransform(const KFState& before,const KFState& after,
    const ZhangGraphBasis& oldBasis,const ZhangGraphBasis& newBasis,bool changeClock)
{
    MatrixXd transform=MatrixXd::Zero(after.x.size(),before.x.size());
    for(const auto& [key,index]:after.kfIndexMap)
    {
        const bool phase=key.Sat.sys==E_Sys::GAL &&
            (key.type==KF::PHASE_BIAS || key.type==KF::AMBIGUITY);
        const bool clock=changeClock && key.Sat.sys==E_Sys::GAL &&
            (key.type==KF::REC_CLOCK || key.type==KF::REC_CLOCK_RATE ||
             key.type==KF::SAT_CLOCK || key.type==KF::SAT_CLOCK_RATE);
        if(phase || clock)continue;
        const auto found=before.kfIndexMap.find(key);check(found!=before.kfIndexMap.end(),"identity complement removed");
        transform(index,found->second)=1;
    }
    for(const auto code:systems[1].codes)
    {
        const auto oldH=phaseH(before,oldBasis,E_Sys::GAL,code);
        const auto newH=phaseH(after,newBasis,E_Sys::GAL,code);
        std::vector<int> columns;
        for(const auto& [key,index]:after.kfIndexMap)
            if(key.Sat.sys==E_Sys::GAL && key.num==int(code) &&
                (key.type==KF::PHASE_BIAS || key.type==KF::AMBIGUITY))columns.push_back(index);
        check(columns.size()==4,"unexpected Galileo phase rank");
        MatrixXd design(4,4);for(int c=0;c<4;++c)design.col(c)=newH.col(columns[c]);
        Eigen::FullPivLU<MatrixXd> decomposition(design);check(decomposition.rank()==4,"independent new H rank defect");
        const MatrixXd solved=decomposition.solve(oldH);
        for(int r=0;r<4;++r)transform.row(columns[r])=solved.row(r);
    }
    if(changeClock)
    for(const auto& [key,index]:after.kfIndexMap)
    {
        if(key.Sat.sys!=E_Sys::GAL)continue;
        const bool receiver=key.type==KF::REC_CLOCK || key.type==KF::REC_CLOCK_RATE;
        const bool satellite=key.type==KF::SAT_CLOCK || key.type==KF::SAT_CLOCK_RATE;
        if(!receiver && !satellite)continue;
        KFKey pivot; pivot.type=(key.type==KF::REC_CLOCK_RATE || key.type==KF::SAT_CLOCK_RATE)
            ? KF::REC_CLOCK_RATE : KF::REC_CLOCK;pivot.Sat=SatSys(E_Sys::GAL,0);pivot.str="B";pivot.num=key.num;
        transform(index,before.kfIndexMap.at(pivot))=-1;
        if(satellite)transform(index,before.kfIndexMap.at(key))=1;
        else check(key.str=="A","unexpected post-reference receiver clock");
    }
    return transform;
}

void compare(const KFState& before,const KFState& after,const MatrixXd& transform,const char* label)
{
    const double meanError=(after.x-transform*before.x).cwiseAbs().maxCoeff();
    const double covarianceError=(after.P-transform*before.P*transform.transpose()).cwiseAbs().maxCoeff();
    check(meanError<1e-8,std::string(label)+" independent full mean differs");
    check(covarianceError<1e-7,std::string(label)+" independent full covariance differs");
    for(const auto& [key,index]:before.kfIndexMap)
        if(key.Sat.sys==E_Sys::GPS || key.Sat.sys==E_Sys::BDS)
        {
            const auto destination=after.kfIndexMap.find(key);check(destination!=after.kfIndexMap.end(),"G/C state identity lost");
            check(after.x(destination->second)==before.x(index),"E datum event changed G/C identity-state mean");
        }
    std::cout<<"PASS "<<label<<" full_mean_max_error="<<meanError
        <<" full_covariance_max_error="<<covarianceError<<" non_target_identity=1\n";
}

std::string datum(const System& system)
{
    std::ostringstream out;
    for(const auto code:system.codes)
    for(const char* id:{system.first,system.second})
    {
        const auto s=zhangSatelliteDatumStatus(system.sys,code,SatSys(id));
        out<<id<<":"<<int(code)<<":"<<s.phaseSegment<<":"<<s.discontinuityCounter
           <<":"<<s.datumVersion<<":"<<s.componentVersion<<":"<<s.alignmentGeneration
           <<":"<<s.alignmentCycles<<":"<<s.componentId<<";";
    }
    return out.str();
}

void jointBridgeScopeContract(KFState& state)
{
    state.lastPppTransactionCommitted=true;state.lastPppQcConverged=true;
    const auto runtime=zhangCheckpointRuntimeId(state);
    const auto epoch=static_cast<long int>(state.time.bigTime);
    const auto root=zhangR51NumericRoot(state.x,state.P);
    const auto order=zhangCurrentEvidenceStateOrderId(state.kfIndexMap);
    auto& registry=zhangCurrentEpochEvidenceRegistry();
    registry.clear();
    ZhangMultiGnssProductDelivery delivery;
    delivery.valid=true;delivery.runtimeId=runtime;delivery.epochGpstSeconds=epoch;
    delivery.floatRoot=root;delivery.floatStateOrderId=order;
    std::vector<ZhangCurrentEpochNumericConstraint> recipe;
    MatrixXd constraints=MatrixXd::Zero(3,state.x.size());VectorXd values(3);
    for(int i=0;i<systems.size();++i)
    {
        const auto& system=systems[i];
        ZhangGraphIntegerContext context;
        check(zhangGraphIntegerContext(state,system.sys,context),"scope fixture graph missing");
        const auto chord=std::find_if(context.basis.edges.begin(),context.basis.edges.end(),
            [&](const auto& edge){return !context.basis.treeEdges.contains(edge);});
        check(chord!=context.basis.edges.end(),"scope fixture cycle missing");
        const auto key=zhangTransportAmbiguityKey(system.codes[0],*chord);
        ZhangCurrentEpochEvidence evidence;evidence.valid=true;evidence.runtimeId=runtime;
        evidence.system=system.sys;evidence.epochGpstSeconds=epoch;
        evidence.posteriorId="DISTINCT_FLOAT_MARGINAL_"+enum_to_string(system.sys);
        evidence.chartId="DISTINCT_PHYSICAL_CHART_"+enum_to_string(system.sys);
        evidence.floatStateRoot=root;evidence.floatPosteriorPsdStateOrderId=order;
        evidence.constraintGeneration=1;
        evidence.numericConstraints.push_back({{{key,1}},i+2});
        registry[{runtime,system.sys}]=evidence;
        auto& scope=delivery.systems[system.sys];
        scope.constraints.system=system.sys;scope.constraints.reliable=true;
        scope.sourcePosteriorId=evidence.posteriorId;scope.authoritativeCycleChartId=evidence.chartId;
        recipe.push_back(evidence.numericConstraints[0]);
        constraints(i,state.kfIndexMap.at(key))=1;values(i)=i+2;
    }
    // Deliberately valid-looking but not actually delivered: never inherit it.
    auto failed=registry.at({runtime,E_Sys::BDS});failed.system=E_Sys::QZS;
    failed.posteriorId="FAILED_SYSTEM_DIAGNOSTIC";failed.numericConstraints[0].value=999;
    registry[{runtime,E_Sys::QZS}]=failed;
    check(zhangMultiGnssAuthorizeCurrentProductEvidence(delivery),"joint evidence scope registration failed");
    const auto& e=registry.at({runtime,E_Sys::GAL});
    auto key=zhangCurrentBridgeConditioningKey(registry,runtime,epoch,e.posteriorId,e.chartId,root,order);
    check(key.size()==3,"E bridge lost G/C current constraints or included failed-system evidence");

    // Independent full joint Schur reference creates the fixture posterior.
    const MatrixXd gram=constraints*state.P*constraints.transpose();
    const MatrixXd gain=state.P*constraints.transpose()*gram.ldlt().solve(MatrixXd::Identity(3,3));
    KFState conditioned=state;
    conditioned.x=state.x+gain*(values-constraints*state.x);
    conditioned.P=state.P-gain*constraints*state.P;
    ZhangProductConditioningIdentity identity;
    identity.runtimeId=runtime;identity.epochGpstSeconds=epoch;
    identity.posteriorId="MULTI_JOINT_FLOAT|"+root;identity.chartId="MULTI_JOINT_STATE_ORDER|"+order;
    identity.sourceRoot=root;identity.sourceStateOrder=order;
    identity.resultRoot=zhangR51NumericRoot(conditioned.x,conditioned.P);
    identity.resultStateOrder=order;
    ZhangProductConditioningReceipt receipt;std::string reason;
    check(receipt.bind(identity,recipe,reason),"joint numeric receipt fixture rejected: "+reason);
    zhangProductConditioningReceiptRegistry()[runtime]=receipt;
    const auto owner=zhangBindCurrentProductFloatOwner(state,E_Sys::GAL,epoch,
        e.posteriorId,e.chartId,reason,&conditioned);
    check(owner && owner->conditionedBaseOwner()==&conditioned && owner->appliedRecipe().size()==3,
        "production E owner factory failed to bind complete G/E/C numeric recipe: "+reason);
    const KFKey targetKey=zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("E11"));
    const std::vector<ZhangCoordinateRow> targets{{{targetKey,1}}};
    const auto response=zhangBuildSmallProductResponse(*owner->conditionedBase(order),targets,{});
    check(response.valid,"production complete-recipe small response rejected");
    const int target=state.kfIndexMap.at(targetKey);
    const VectorXd projected=response.targetProjector*response.mean;
    const MatrixXd projectedP=response.targetProjector*response.covariance*response.targetProjector.transpose();
    check(std::abs(projected[0]-conditioned.x[target])<1e-10 &&
        std::abs(projectedP(0,0)-conditioned.P(target,target))<1e-9,
        "E bridge numeric response lost cross-system conditional covariance");
    check(std::abs(conditioned.P(target,target)-state.P(target,target))>1e-5,
        "cross-system conditioning test fixture was degenerate");
    delivery.systems.erase(E_Sys::BDS);
    check(zhangMultiGnssAuthorizeCurrentProductEvidence(delivery),"partial successful scope rejected");
    key=zhangCurrentBridgeConditioningKey(registry,runtime,epoch,e.posteriorId,e.chartId,root,order);
    check(key.size()==2,"failed BDS registry residual entered joint bridge recipe");
    const auto stale=zhangBindCurrentProductFloatOwner(state,E_Sys::GAL,epoch,
        e.posteriorId,e.chartId,reason,&conditioned);
    check(stale && !stale->conditionedBaseOwner(),
        "conditioned base containing a failed-system row was reused");
    registry.clear();zhangMultiGnssCurrentEvidenceAuthorities().clear();
    zhangProductConditioningReceiptRegistry().erase(runtime);
    std::cout<<"PASS production_E_bridge_owner distinct_per_system_IDs full_joint_numeric_recipe=1 failed_registry_rows_excluded=1\n";
}

void seedLedgers(const KFState& state)
{
    const auto runtime=zhangCheckpointRuntimeId(state);
    for(const auto& system:systems)
    {
        ProductIntegerLedgerRow row;row.system=system.sys;
        row.firstObservable=system.codes[0];row.secondObservable=system.codes[1];
        row.productRow={1,0};row.integerValue=7;row.phaseSegmentFingerprint="SYNTHETIC_ISOLATION_FIXTURE";
        const auto code=enum_to_string(system.codes[0]);
        row.physicalExpansion={{code+"|A|"+system.first+"|V0",1},{code+"|B|"+system.second+"|V0",-1}};
        row.canonicalProductExpansion={{code+"|"+system.first,1},{code+"|"+system.second,-1}};
        check(zhangProductIntegerLedgerRegistry()[{runtime,system.sys}].observe(
            static_cast<long int>(state.time.bigTime),{row},1).valid,"synthetic real-ledger insertion failed");
        for(const auto signal:system.codes)
            check(promoteZhangSatelliteProductRelation(state.time,system.sys,signal,
                SatSys(system.first),SatSys(system.second),7,"SYNTHETIC_RUNTIME_ISOLATION"),
                "real datum-manager fixture relation failed");
    }
}

std::string ledger(const KFState& state,E_Sys system)
{
    std::ostringstream out;
    for(const auto& row:zhangProductIntegerLedgerRegistry().at({zhangCheckpointRuntimeId(state),system}).rows())
        out<<zhangProductLedgerIdentityFingerprint(row)<<":"<<row.integerValue
           <<":"<<row.backendBasisGeneration<<":"<<row.lastConfirmed<<":"<<row.confirmationEpochs
           <<":"<<row.certified<<":"<<zhangProductPhysicalRowFingerprint(row.physicalExpansion)<<";";
    return out.str();
}

void navigationParserContract(const std::string& path,int expectedVersion)
{
    std::ifstream input(path);check(bool(input),"NAV input not readable: "+path);
    Navigation navigation;ObsList observations;RinexStation station;
    char type=0;double version=0;E_Sys system=E_Sys::NONE;E_TimeSys timeSystem=E_TimeSys::GPST;
    std::map<E_Sys,std::map<int,CodeType>> codes;
    const int status=readRnx(input,type,observations,navigation,station,version,system,timeSystem,codes);
    std::cout<<"NAV_PARSER file="<<std::quoted(path)<<" version="<<version
        <<" header_type="<<type<<" status="<<status<<" parser=PRODUCTION_readRnx\n";
    check(status>0 && type=='N' && int(version)==expectedVersion,
        "NAV header/version/production parser result invalid: "+path);
    std::map<E_Sys,std::size_t> totals;
    std::map<E_NavMsgType,std::size_t> bdsCivil;
    auto count=[&](const auto& map,const char* name)
    {
        std::map<std::pair<E_Sys,E_NavMsgType>,std::size_t> messages;
        for(const auto& [satellite,types]:map)
        for(const auto& [message,epochs]:types)
        {
            messages[{satellite.sys,message}]+=epochs.size();
            totals[satellite.sys]+=epochs.size();
            if(std::string(name)=="cephMap" && satellite.sys==E_Sys::BDS)
                bdsCivil[message]+=epochs.size();
        }
        for(const auto& [identity,number]:messages)
            std::cout<<"NAV_COUNT map="<<name<<" system="<<enum_to_string(identity.first)
                <<" message="<<enum_to_string(identity.second)<<" records="<<number<<"\n";
    };
    count(navigation.ephMap,"ephMap");count(navigation.cephMap,"cephMap");
    for(const auto system:{E_Sys::GPS,E_Sys::GAL,E_Sys::BDS})
        check(totals[system]>0,"production NAV parser produced zero "+enum_to_string(system)+" ephemerides");
    if(expectedVersion>=4)
        for(const auto message:{E_NavMsgType::CNV1,E_NavMsgType::CNV2,E_NavMsgType::CNV3})
            check(bdsCivil[message]>0,"actual BDS "+enum_to_string(message)+" absent from production cephMap");
    std::cout<<"PASS production_NAV_parser version="<<expectedVersion
        <<" actual_ephemeris_records_nonzero_G_E_C=1 civil_CNV123_required="<<(expectedVersion>=4)
        <<" scientific_orbit_quality_accepted=0\n";
}
}

int main(int argc,char** argv)
{
    std::ostringstream trace;
    try
    {
        if(argc>1 && std::string(argv[1])=="--nav-only")
        {
            check(argc==4,"--nav-only requires RINEX3 and RINEX4 NAV paths");
            navigationParserContract(argv[2],3);navigationParserContract(argv[3],4);
            return 0;
        }
        configure(true);ReceiverMap receivers;observations(receivers);
        KFState state;state.time.bigTime=1405216800;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="MULTI_GNSS_RUNTIME_TREE";
        updateZhangFullRankReferences(trace,receivers,state);states(state,trace,true);
        jointBridgeScopeContract(state);
        ZhangGraphIntegerContext oldE;check(zhangGraphIntegerContext(state,E_Sys::GAL,oldE),"E graph missing");
        auto newBasis=oldE.basis;
        const auto oldChord=*std::find_if(newBasis.edges.begin(),newBasis.edges.end(),
            [&](const auto& edge){return !newBasis.treeEdges.contains(edge);});
        newBasis.treeEdges.insert(oldChord);
        // K2,2 has one four-edge cycle. Removing any OLD tree edge after
        // inserting the old chord is connected and guarantees a real exchange.
        const auto remove=*oldE.basis.treeEdges.begin();
        newBasis.treeEdges.erase(remove);
        check(newBasis.treeEdges!=oldE.basis.treeEdges,"tree exchange fixture was a no-op");
        const KFState before=state;SparseMatrix<double> productionTransform;std::string reason;
        check(applyZhangGraphBasisTransformForAudit(trace,state,E_Sys::GAL,systems[1].codes,
            oldE.basis,newBasis,productionTransform,reason),"production E-only tree transform failed: "+reason);
        const auto independent=independentTransform(before,state,oldE.basis,newBasis,false);
        compare(before,state,independent,"production_E_tree_exchange");
        const auto beforeH=phaseH(before,oldE.basis,E_Sys::GAL,E_ObsCode::L1C);
        const auto afterH=phaseH(state,newBasis,E_Sys::GAL,E_ObsCode::L1C);
        check((beforeH*before.x-afterH*state.x).norm()<1e-8,"frozen physical predictions changed on pure tree exchange");
        check((beforeH*before.P*beforeH.transpose()-afterH*state.P*afterH.transpose()).norm()<1e-7,
            "frozen physical prediction covariance changed on pure tree exchange");
        eraseZhangGraphRuntime(state);

        // Trigger a real E11 single-band retirement through the public graph
        // update. Repeat within this same global epoch: independently-AR
        // epoch clearing is intentional and is not an E-only side effect.
        configure(true);ReceiverMap slipReceivers;observations(slipReceivers);
        KFState slipState;slipState.time.bigTime=1405216830;
        slipState.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="MULTI_GNSS_RUNTIME_SLIP";
        updateZhangFullRankReferences(trace,slipReceivers,slipState);states(slipState,trace,true);
        seedLedgers(slipState);
        std::map<E_Sys,ZhangGraphIntegerContext> oldContexts;
        std::map<E_Sys,std::string> oldLedgers,oldDatums;
        for(const auto& system:systems)
        {
            check(zhangGraphIntegerContext(slipState,system.sys,oldContexts[system.sys]),"graph context missing before slip");
            oldLedgers[system.sys]=ledger(slipState,system.sys);oldDatums[system.sys]=datum(system);
        }
        const ZhangGraphEdge broken{"A",SatSys("E11")};
        auto& slipStatus=slipReceivers.at("A").satStatMap.at(SatSys("E11"));
        auto& signal=slipStatus.sigStatMap[ft2string(code2Freq.at(E_Sys::GAL).at(E_ObsCode::L1C))];
        signal.slip.LLI=1;signal.tracking.pendingBreak=true;
        updateZhangFullRankReferences(trace,slipReceivers,slipState);
        ZhangGraphIntegerContext newE;check(zhangGraphIntegerContext(slipState,E_Sys::GAL,newE),"E graph missing after slip");
        const auto oldVersion=zhangGraphSignalArcVersion(oldContexts.at(E_Sys::GAL),broken,E_ObsCode::L1C);
        const auto newVersion=zhangGraphSignalArcVersion(newE,broken,E_ObsCode::L1C);
        check(oldVersion && newVersion && *newVersion>*oldVersion,"real E11 single-band physical arc was inherited");
        for(const auto& system:systems)
        {
            if(system.sys==E_Sys::GAL)continue;
            ZhangGraphIntegerContext after;check(zhangGraphIntegerContext(slipState,system.sys,after),"G/C graph disappeared");
            const auto& old=oldContexts.at(system.sys);
            check(after.signalArcVersions==old.signalArcVersions && after.arcVersions==old.arcVersions &&
                after.basis.treeEdges==old.basis.treeEdges && after.eventId==old.eventId &&
                after.productDatumVersion==old.productDatumVersion &&
                after.integerComponentVersion==old.integerComponentVersion,
                "E11 slip changed unrelated G/C integer graph identities");
            check(ledger(slipState,system.sys)==oldLedgers.at(system.sys),"E11 slip changed G/C ledger fixture");
            check(datum(system)==oldDatums.at(system.sys),"E11 slip changed G/C public datum manager counters");
        }
        std::cout<<"PASS production_E11_L1C_slip E_version_before="<<*oldVersion
            <<" E_version_after="<<*newVersion<<" G_C_graph_ledger_public_datum_unchanged=1\n";
        eraseZhangGraphRuntime(slipState);

        // Public legacy-reference path, not merely the clock-key helper.
        configure(false);ReceiverMap legacyReceivers;observations(legacyReceivers);
        KFState legacy;legacy.time.bigTime=1405216860;
        legacy.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="MULTI_GNSS_RUNTIME_LEGACY_CLOCK";
        states(legacy,trace,false);
        const KFState legacyBefore=legacy;
        auto& aObservations=legacyReceivers.at("A").obsList;
        for(auto it=aObservations.begin();it!=aObservations.end();)
        {
            const auto observation=std::dynamic_pointer_cast<GObs>(*it);
            if(observation && observation->Sat.sys==E_Sys::GAL)it=aObservations.erase(it);else ++it;
        }
        updateZhangFullRankReferences(trace,legacyReceivers,legacy);
        check(acsConfig.zhangFullRank.sysOpts.at(E_Sys::GAL).reference_receiver=="B",
            "public update did not execute E-only receiver reference change");
        check(acsConfig.zhangFullRank.sysOpts.at(E_Sys::GPS).reference_receiver=="A" &&
            acsConfig.zhangFullRank.sysOpts.at(E_Sys::BDS).reference_receiver=="A",
            "E-only reference change modified G/C configured datums");
        const auto legacyT=independentTransform(legacyBefore,legacy,star(systems[1]),star(systems[1],"B"),true);
        compare(legacyBefore,legacy,legacyT,"public_E_receiver_clock_reference_change");
        if(argc>2)
        {
            check(argc==4,"optional NAV test requires both RINEX3 and RINEX4 paths");
            navigationParserContract(argv[2],3);navigationParserContract(argv[3],4);
        }
        std::cout<<"PASS multi-GNSS production runtime controlled tests; synthetic inputs only; no PPP filter run or full-hour acceptance\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<"FAIL "<<error.what()<<"\n";
        if(argc>1 && std::string(argv[1])!="--nav-only")
        {std::ofstream file(argv[1]);file<<trace.str();}
        return 1;
    }
    if(argc>1){std::ofstream file(argv[1]);file<<trace.str();}
    return 0;
}
