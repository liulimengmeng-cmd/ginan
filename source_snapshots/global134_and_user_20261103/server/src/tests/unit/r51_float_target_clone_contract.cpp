// Real graph clone overloads; no PEA/ILS/publication claim.
#include "common/acsConfig.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/receiver.hpp"
#include "common/observations.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangCertifiedTreeTransport.hpp"
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "pea/zhangReference.hpp"
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

static unsigned checks=0;
static void check(bool value,const std::string& why)
{++checks;if(!value)throw std::runtime_error(why);}
static KFKey phase(int prn)
{return zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys(E_Sys::GPS,prn));}
struct Fixture
{
    KFState state;ReceiverMap receivers;std::ostringstream trace;KFKey carrier;
    explicit Fixture(bool typed=true)
    {
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=typed?"FLOAT_CLONE_TYPED":"FLOAT_CLONE_LEGACY";
        for(const auto& edge:std::set<ZhangGraphEdge>{{"A",SatSys("G01")},{"A",SatSys("G02")},{"B",SatSys("G01")}})
        {
            auto& rec=receivers[edge.receiver];rec.id=edge.receiver;rec.ready=true;
            auto obs=std::make_shared<GObs>();obs->Sat=edge.satellite;
            auto& stat=rec.satStatMap[edge.satellite];stat.el=1.5;obs->satStat_ptr=&stat;
            for(const auto code:{E_ObsCode::L1C,E_ObsCode::L2W})
            {auto& sig=obs->sigs[code2Freq[E_Sys::GPS][code]];sig.code=code;sig.L=1e8;sig.P=2e7;}
            rec.obsList.push_back(obs);
        }
        updateZhangFullRankReferences(trace,receivers,state);
        auto add=[&](const KFKey& key,double mean){InitialState initial;initial.estimate=true;initial.x=mean;
            initial.P=1;initial.Q=0;state.addKFState(key,initial);};
        for(const auto code:{E_ObsCode::L1C,E_ObsCode::L2W})
        {add(zhangTransportReceiverKey(E_Sys::GPS,code,"B"),20);
            for(int sat:{1,2})add(zhangTransportSatelliteKey(code,SatSys(E_Sys::GPS,sat)),10+sat);}
        if(typed){carrier.type=KF::PRODUCT_HANDOFF;carrier.str="ZHANG_FLOAT_CARRIER:CLONE_SUM";add(carrier,23);}
        check(state.stateTransition(trace,state.time),"fixture state creation failed");
        const auto v=VectorXd::LinSpaced(state.x.size(),.01,.2).eval();
        state.P=MatrixXd::Identity(state.x.size(),state.x.size())+v*v.transpose();
        const int one=state.kfIndexMap.at(KFState::oneKey);state.x(one)=1;state.dx(one)=0;
        state.P.row(one).setZero();state.P.col(one).setZero();
        if(!typed)return;
        const int c=state.kfIndexMap.at(carrier),a=state.kfIndexMap.at(phase(1)),b=state.kfIndexMap.at(phase(2));
        state.x(c)=state.x(a)+state.x(b);state.dx(c)=state.dx(a)+state.dx(b);
        for(int i=0;i<state.P.rows();++i)if(i!=c)state.P(c,i)=state.P(a,i)+state.P(b,i);
        for(int i=0;i<state.P.rows();++i)if(i!=c)state.P(i,c)=state.P(c,i);
        state.P(c,c)=state.P(a,a)+state.P(b,b)+2*state.P(a,b);
        ZhangProductFloatTargetRequest target;target.targetId="live-whole-phase";target.consumer="CLONE_CONTRACT";
        target.lifecycleId="SEGMENT_1";target.unit="m";target.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
        target.physicalCoefficients[{{"A",SatSys("G01")},E_ObsCode::L1C,0}]=1;
        target.physicalCoefficients[{{"A",SatSys("G02")},E_ObsCode::L1C,0}]=1;
        target.scope={zhangProductFloatOwnerId(state),"CLONE_SOURCE_CHART",zhangFloatCarrierStateOrder(state.kfIndexMap),state.factorCommitSequence};
        target.numericRow={{carrier,1}};target.numericOffset=3.25;target.numericBindingComplete=true;
        target.bindingReason="LITERAL_COMPLETE_FLOAT_SOURCE";
        state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
        auto definition=target;definition.targetId="carrier-source-definition";definition.consumer="FLOAT_CARRIER_SOURCE_DEFINITION";
        definition.numericRow={{phase(1),1},{phase(2),1}};definition.numericOffset=0;
        state.metaDataMap[zhangProductFloatCarrierRecordKey(carrier)]=zhangEncodeFloatTarget(definition);
        state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]=target.scope.stateOrderId;
        state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(state.factorCommitSequence);
        state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(state.factorCommitSequence);
        // This is an owner-isolation fixture, not an accepted real AR fact.
        zhang_transport_certificate::Directory history;history.authoritativeOwner=zhangProductFloatOwnerId(state);
        history.system=int(E_Sys::GPS);history.lastCommittedEpoch=-1;
        history.sourceScopes["owner-contract"]={-1,"p","c","r","o"};
        history.proofAtoms["literal"]="LITERAL_OWNER_ISOLATION_NOT_REAL_AR";
        zhang_transport_certificate::Fact fact;fact.physical={{"L1C|A|G01|V0",1}};fact.value=3;
        fact.sourceScopeId="owner-contract";fact.proofIds={"literal"};history.facts={fact};
        const auto bytes=history.encode();state.metaDataMap[std::string(zhang_transport_certificate::metadataKey)+"/GPS"]=zhangCheckpointSha256(bytes)+"\n"+bytes;
    }
    ~Fixture(){eraseZhangGraphRuntime(state);}
};
static void bind(KFState& state,const KFState& source,const std::string& suffix)
{state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA]=zhangProductFloatOwnerId(source)+"/branch/"+suffix;}
static void positiveAndRejectedCopies()
{
    Fixture f;std::string reason;check(zhangProductFloatTargetMetadataValid(f.state,reason),"source metadata invalid:"+reason);
    const auto before=captureZhangCheckpointKfCore(f.state);
    KFState branch=f.state;bind(branch,f.state,"good");
    const auto oldHistory=branch.metaDataMap.at(std::string(zhang_transport_certificate::metadataKey)+"/GPS");
    check(cloneZhangGraphRuntime(f.state,branch,&reason),"real mutable V2 clone rejected:"+reason);
    check(zhangProductFloatTargetMetadataValid(branch,reason),"cloned V2 metadata invalid:"+reason);
    std::vector<ZhangProductFloatTargetRequest> targets;check(zhangReadFloatTargets(branch.metaDataMap,targets,reason) && targets.size()==1,"target missing");
    check(targets[0].scope.runtimeId==zhangProductFloatOwnerId(branch) && targets[0].numericRow==std::map<KFKey,double>{{f.carrier,1}} && targets[0].numericOffset==3.25,
        "clone changed complete phase row or failed owner rebind");
    ZhangProductFloatTargetRequest carrier;check(zhangDecodeFloatTarget(branch.metaDataMap.at(zhangProductFloatCarrierRecordKey(f.carrier)),carrier) &&
        carrier.scope.runtimeId==zhangProductFloatOwnerId(branch),"carrier definition owner not rebound");
    ZhangGraphIntegerContext graph;check(zhangGraphIntegerContext(branch,E_Sys::GPS,graph),"real branch graph was not cloned");
    check(branch.metaDataMap.at(std::string(zhang_transport_certificate::metadataKey)+"/GPS")==oldHistory,"clone promoted authoritative history owner");
    ZhangProductIntegerTransition old;old.valid=true;old.observable=E_ObsCode::L1C;
    old.physicalEdges={{"A",SatSys("G01")}};old.physicalArcVersions={0};old.coefficients={1};
    check(!zhangQueryCertifiedHistoricalTransport(branch,E_Sys::GPS,old).proven,"private V2 clone gained authoritative history permission");
    check(zhangCheckpointKfCoreBitwiseEqual(before,captureZhangCheckpointKfCore(f.state)),"clone changed source raw core");
    const auto cloned=captureZhangCheckpointKfCore(branch);const auto counter=branch.factorCommitSequence;
    const KFState& immutable=branch;
    check(!cloneZhangGraphRuntime(f.state,immutable,&reason),"const V2 destination escaped mutable clone gate");
    check(counter==branch.factorCommitSequence && zhangCheckpointKfCoreBitwiseEqual(cloned,captureZhangCheckpointKfCore(branch)),"const rejection mutated destination");
    eraseZhangGraphRuntime(branch);
    for(const auto& mutate:std::vector<std::function<void(KFState&)>>{
        [](auto&x){x.x(1)+=.125;},[](auto&x){x.dx(1)+=.125;},[](auto&x){x.P(1,2)+=.125;},
        [](auto&x){++x.factorCommitSequence;},
        [](auto&x){auto n=x.kfIndexMap.extract(x.kfIndexMap.begin());n.key().comment="DIFFERENT_KEY_ANNOTATION";x.kfIndexMap.insert(std::move(n));},
        [](auto&x){x.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]+="DIFFERENT_ORDER";},
        [](auto&x){x.metaDataMap.erase(ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2);}})
    {
        KFState bad=f.state;bind(bad,f.state,"rejected");mutate(bad);
        const auto prior=captureZhangCheckpointKfCore(bad);const auto seq=bad.factorCommitSequence;
        check(!cloneZhangGraphRuntime(f.state,bad,&reason),"changed posterior/order/generation was cloned as identical");
        check(seq==bad.factorCommitSequence && zhangCheckpointKfCoreBitwiseEqual(prior,captureZhangCheckpointKfCore(bad)),"rejected clone changed destination");
        check(!zhangGraphIntegerContext(bad,E_Sys::GPS,graph),"rejected clone published a graph runtime");
    }
}
static void legacyConstPath()
{
    Fixture f(false);KFState branch=f.state;bind(branch,f.state,"legacy");
    const auto metadata=branch.metaDataMap;const KFState& destination=branch;std::string reason;
    check(cloneZhangGraphRuntime(f.state,destination,&reason),"legacy const clone behavior changed:"+reason);
    check(branch.metaDataMap==metadata && !hasZhangFloatTargetV2State(branch.metaDataMap),"legacy clone injected empty V2 schema");
    eraseZhangGraphRuntime(branch);
}
static void staleCompleteRegistrationRejected()
{
    Fixture f;std::string reason;
    ZhangProductFloatTargetRequest request;
    request.targetId="stale-same-order-source";request.consumer="REGISTRATION_GENERATION_CONTRACT";
    request.lifecycleId="SEGMENT_STALE";request.unit="m";
    request.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
    request.physicalCoefficients[{{"A",SatSys("G01")},E_ObsCode::L1C,0}]=1;
    request.scope={zhangProductFloatOwnerId(f.state),"LITERAL_PRE_TRANSFORM_CHART",
        zhangFloatCarrierStateOrder(f.state.kfIndexMap),f.state.factorCommitSequence};
    request.numericRow={{phase(1),1}};request.numericBindingComplete=true;
    check(zhangRegisterGraphProductHandoffTarget(f.state,request,reason),"current complete source request rejected:"+reason);
    std::map<KFKey,std::map<KFKey,double>> transform;
    for(const auto&[key,index]:f.state.kfIndexMap)transform[key][key]=1;
    transform[phase(1)][phase(2)]=1;
    check(f.state.applyStateTransform(f.trace,transform,"TEST_SAME_ORDER_REGISTRATION_GENERATION"),"same-order real transform failed");
    const auto before=captureZhangCheckpointKfCore(f.state);const auto seq=f.state.factorCommitSequence;
    check(!zhangRegisterGraphProductHandoffTarget(f.state,request,reason),
        "old a row was silently relabelled current after a'=a+b");
    check(seq==f.state.factorCommitSequence && zhangCheckpointKfCoreBitwiseEqual(before,captureZhangCheckpointKfCore(f.state)),
        "stale complete registration rejection changed current state or directory");
}
int main()
{
    try {
        setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
        setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
        acsConfig.zhangFullRank.enable=true;acsConfig.zhangFullRank.product_basis_mode="SINGLE_S_BASIS_NATIVE";
        acsConfig.process_sys[E_Sys::GPS]=true;auto& opt=acsConfig.zhangFullRank.sysOpts[E_Sys::GPS];
        opt.use_spanning_tree=true;opt.reference_receiver="A";opt.baseline_observables={E_ObsCode::L1C,E_ObsCode::L2W};opt.product_core_min_satellite_support=0;
        positiveAndRejectedCopies();legacyConstPath();staleCompleteRegistrationRejected();
        std::cout<<"FLOAT_TARGET_CLONE_CONTRACT checks="<<checks<<" passed=1 real_graph_clone=1 raw_same_posterior_gate=1 target_carrier_owner_rebind=1 history_promotion=0 const_V2_rejected=1 legacy_no_tag=1 real_PEA_or_AR_executed=0\n";return 0;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
