// Exercise the real central KF transaction, not a second remap implementation.
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void check(bool ok,const char* why)
{if(!ok)throw std::runtime_error(why);}
struct Fixture
{
    KFState state;
    KFKey a,b,clock;
    Fixture()
    {
        a.type=KF::PHASE_BIAS;a.Sat=SatSys("G01");a.num=int(E_ObsCode::L1C);
        b=a;b.Sat=SatSys("G02");
        clock.type=KF::SAT_CLOCK;clock.Sat=SatSys("G01");
        state.time.bigTime=1000;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="FLOAT_META_TRANSFORM_TEST";
        state.kfIndexMap={{KFState::oneKey,0},{a,1},{b,2},{clock,3}};
        state.x.resize(4);state.x<<1,2,5,10;state.dx=VectorXd::Zero(4);
        state.P.resize(4,4);
        state.P<<0,0,0,0, 0,4,1,2, 0,1,9,3, 0,2,3,16;
        for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
        target({{a,1},{b,1}});
    }
    void target(const std::map<KFKey,double>& row,bool bound=true)
    {
        ZhangProductFloatTargetRequest value;
        value.targetId="active-whole-phase";value.consumer="CENTRAL_TRANSFORM_CONTRACT";
        value.lifecycleId="SEGMENT_1";value.unit="m";
        value.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
        value.physicalCoefficients[{{"A",SatSys("G01")},E_ObsCode::L1C,0}]=1;
        value.physicalCoefficients[{{"A",SatSys("G02")},E_ObsCode::L1C,0}]=1;
        value.scope={zhangProductFloatOwnerId(state),"LITERAL_SOURCE_CHART",
            zhangFloatCarrierStateOrder(state.kfIndexMap),state.factorCommitSequence};
        value.numericRow=row;value.numericOffset=3.25;
        value.active=true;value.numericBindingComplete=bound;
        value.bindingReason=bound?"LITERAL_COMPLETE_SOURCE_ROW":"EXPLICIT_UNBOUND_CONSUMER";
        state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+value.targetId]=zhangEncodeFloatTarget(value);
        state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]=value.scope.stateOrderId;
        state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(value.scope.coordinateGeneration);
        state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(state.factorCommitSequence);
    }
    ZhangProductFloatTargetRequest read() const
    {
        std::vector<ZhangProductFloatTargetRequest> targets;std::string reason;
        check(zhangReadFloatTargets(state.metaDataMap,targets,reason) && targets.size()==1,
            "target metadata did not survive its wire roundtrip");
        return targets.front();
    }
    std::map<KFKey,std::map<KFKey,double>> sumTransform() const
    {return {{KFState::oneKey,{{KFState::oneKey,1}}},{a,{{a,1},{b,1}}},{clock,{{clock,1}}}};}
};

static void unchanged(const KFState& actual,const KFState& before)
{
    check(actual.kfIndexMap==before.kfIndexMap && actual.metaDataMap==before.metaDataMap &&
        actual.stateTransitionMap==before.stateTransitionMap &&
        actual.factorCommitSequence==before.factorCommitSequence,
        "rejected transaction changed keys, metadata, models or sequence");
    check(actual.x.size()==before.x.size() && actual.P.size()==before.P.size() &&
        !std::memcmp(actual.x.data(),before.x.data(),actual.x.size()*sizeof(double)) &&
        !std::memcmp(actual.P.data(),before.P.data(),actual.P.size()*sizeof(double)) &&
        !std::memcmp(actual.dx.data(),before.dx.data(),actual.dx.size()*sizeof(double)),
        "rejected transaction changed raw posterior bytes");
}
static void realTransform()
{
    Fixture f;
    std::ostringstream trace;
    check(f.state.applyStateTransform(trace,f.sumTransform(),"TEST_WHOLE_ROW_SUM"),
        "real whole-target transform failed");
    const auto result=f.read();
    check(result.numericBindingComplete && result.numericRow==std::map<KFKey,double>{{f.a,1}} &&
        result.numericOffset==3.25 && result.scope.stateOrderId==zhangFloatCarrierStateOrder(f.state.kfIndexMap),
        "source row was left in old coordinates after the real transform");
    const int a=f.state.kfIndexMap.at(f.a),clock=f.state.kfIndexMap.at(f.clock);
    check(f.state.x(a)==7 && f.state.P(a,a)==15 && f.state.P(a,clock)==5,
        "whole sum mean/variance/clock covariance differs from literal oracle");
    std::string reason;
    check(zhangProductFloatTargetMetadataValid(f.state,reason),"new target binding is not valid on its actual state");
}
static void rejectedTransactions()
{
    std::ostringstream trace;
    Fixture f;f.target({{f.a,1},{f.b,-1}});
    const auto before=f.state;
    check(!f.state.applyStateTransform(trace,f.sumTransform(),ZHANG_FLOAT_TARGET_COMPRESSION_LABEL),
        "strict compression lost an active independent difference");
    unchanged(f.state,before);
    Fixture unbound;unbound.target({{unbound.a,1},{unbound.b,-1}},false);
    const auto unboundBefore=unbound.state;
    check(!unbound.state.applyStateTransform(trace,unbound.sumTransform(),ZHANG_FLOAT_TARGET_COMPRESSION_LABEL),
        "strict compression silently skipped an active unbound consumer");
    unchanged(unbound.state,unboundBefore);
    Fixture badOwner;auto value=badOwner.read();value.scope.runtimeId="UNRELATED_OWNER";
    badOwner.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+value.targetId]=zhangEncodeFloatTarget(value);
    const auto ownerBefore=badOwner.state;
    check(!badOwner.state.applyStateTransform(trace,badOwner.sumTransform(),"TEST_WRONG_OWNER"),
        "central transaction transported another owner's target metadata");
    unchanged(badOwner.state,ownerBefore);
    Fixture rejected;
    int callbacks=0;
    rejected.state.exactStateTransformCallback=[&](const KFState&,GTime,const std::map<KFKey,int>&,
        const std::map<KFKey,int>&,const SparseMatrix<double>&,const std::string&,
        const VectorXd&,const MatrixXd&,std::uint64_t,std::uint64_t){++callbacks;return false;};
    const auto callbackBefore=rejected.state;
    check(!rejected.state.applyStateTransform(trace,rejected.sumTransform(),"TEST_REJECTING_FACTOR_CONSUMER"),
        "callback rejection was ignored");
    unchanged(rejected.state,callbackBefore);
    check(callbacks==1,"callback was not reached by the valid staged remap");
}
static void metadataOverride()
{
    Fixture f;std::ostringstream trace;
    auto proposed=f.state.metaDataMap;
    proposed["ZHANG_HANDOFF_FLOAT_TEST_STAGED_RAW"]="EXPLICIT_NEW_RAW_METADATA";
    check(f.state.applyStateTransform(trace,f.sumTransform(),"TEST_METADATA_OVERRIDE",{},&proposed),
        "legitimate staged raw metadata override failed");
    check(f.state.metaDataMap.at("ZHANG_HANDOFF_FLOAT_TEST_STAGED_RAW")=="EXPLICIT_NEW_RAW_METADATA" &&
        f.read().numericRow==std::map<KFKey,double>{{f.a,1}},
        "old capture metadata overwrote the new whole-target coordinates");
    Fixture owner;auto forged=owner.state.metaDataMap;
    forged[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="FORGED_ROOT_OWNER";
    const auto before=owner.state;
    check(!owner.state.applyStateTransform(trace,owner.sumTransform(),"TEST_FORGED_OVERRIDE",{},&forged),
        "raw override replaced a core runtime identity");
    unchanged(owner.state,before);
}
static void timeTransition()
{
    Fixture f;std::ostringstream trace;
    KFKey extra;extra.type=KF::REC_CLOCK;extra.str="NEW_STATE_BEFORE_PHASE_KEYS";
    InitialState initial;initial.estimate=true;initial.x=42;initial.P=2;initial.Q=0;
    f.state.addKFState(extra,initial);
    check(f.state.stateTransition(trace,f.state.time),"new-state real initialization failed");
    const auto bound=f.read();
    check(bound.numericBindingComplete && bound.numericRow==std::map<KFKey,double>{{f.a,1},{f.b,1}} &&
        bound.scope.stateOrderId==zhangFloatCarrierStateOrder(f.state.kfIndexMap),
        "static target did not survive actual state-index insertion");
    // A random-walk rate is an actual process model only with tau < 0.
    // Hand-writing procNoiseMap alone leaves tau absent and production Q=0.
    check(f.state.setProcessModel(f.a,.25,0,-1),
        "fixture could not install the public random-walk process model");
    const auto beforeIndex=f.state.kfIndexMap;
    const double beforeVariance=f.state.P(beforeIndex.at(f.a),beforeIndex.at(f.a));
    MatrixXd actualF,actualQ;
    check(f.state.stateTransition(trace,f.state.time+30,&actualF,&actualQ),
        "legitimate noisy FLOAT prediction failed");
    const int currentA=f.state.kfIndexMap.at(f.a);
    check(actualF.rows()==f.state.x.size() && actualF.cols()==beforeIndex.size() &&
        actualF(currentA,beforeIndex.at(f.a))==1 && actualQ.rows()==f.state.x.size() &&
        actualQ.cols()==f.state.x.size() && actualQ(currentA,currentA)==7.5 &&
        f.state.P(currentA,currentA)==beforeVariance+7.5,
        "fixture did not inject the actual 0.25 m^2/s times 30 s process covariance");
    VectorXd destinationTarget=VectorXd::Zero(f.state.x.size());
    VectorXd sourceTarget=VectorXd::Zero(beforeIndex.size());
    for(const auto& [key,coefficient]:bound.numericRow)
    {
        destinationTarget(f.state.kfIndexMap.at(key))+=coefficient;
        sourceTarget(beforeIndex.at(key))+=coefficient;
    }
    const VectorXd pulledTarget=actualF.transpose()*destinationTarget;
    check((pulledTarget.array()==sourceTarget.array()).all(),
        "fixture changed the old whole target before injecting fresh process noise");
    check((destinationTarget.transpose()*actualQ*destinationTarget)(0,0)==7.5,
        "the old whole product row has no verified fresh-noise contribution");
    const auto unbound=f.read();
    check(!unbound.numericBindingComplete &&
        unbound.bindingReason=="UNBOUND_TIME_TRANSITION_NONIDENTITY_OR_RANDOM_SOURCE" &&
        unbound.physicalCoefficients==bound.physicalCoefficients,
        "new random source was certified as an exact old product transport");
}
static void tinyNonzeroIsNotExact()
{
    std::ostringstream trace;
    Fixture near;near.target({{near.a,1}});
    auto transform=near.sumTransform();transform[near.a][near.b]=1e-12;
    const auto before=near.state;
    check(!near.state.applyStateTransform(trace,transform,ZHANG_FLOAT_TARGET_COMPRESSION_LABEL),
        "tiny nonzero lost direction was rounded into exact target support");
    unchanged(near.state,before);
    check(near.state.applyStateTransform(trace,transform,"TEST_TINY_NONZERO_LEGITIMATE_KF_MAP") &&
        !near.read().numericBindingComplete,
        "legitimate non-exact map did not leave a typed unbound old target");
    Fixture noisy;noisy.target({{noisy.a,1}});
    KFKey source;source.type=KF::AMBIGUITY;source.Sat=SatSys("G03");source.str="INDEPENDENT_SOURCE";
    const std::map<KFKey,std::map<KFKey,double>> stochastic{
        {KFState::oneKey,{{KFState::oneKey,1}}},
        {noisy.a,{{noisy.a,1},{source,1e-12}}},{noisy.b,{{noisy.b,1}}},
        {noisy.clock,{{noisy.clock,1}}}};
    const auto noisyBefore=noisy.state;
    check(!noisy.state.applyStateTransform(trace,stochastic,ZHANG_FLOAT_TARGET_COMPRESSION_LABEL,{{source,1}}),
        "tiny nonzero fresh prior noise was certified as exact old target transport");
    unchanged(noisy.state,noisyBefore);
    check(noisy.state.applyStateTransform(trace,stochastic,"TEST_TINY_FRESH_PRIOR",{{source,1}}) &&
        !noisy.read().numericBindingComplete,
        "a legitimate fresh-noise map did not preserve typed non-exact status");
}
static void staleSameOrderRegistration()
{
    Fixture f;f.target({{f.a,1}});
    std::ostringstream trace;
    std::map<KFKey,std::map<KFKey,double>> canonical;
    for(const auto& [key,index]:f.state.kfIndexMap)canonical[key][key]=1;
    check(f.state.applyStateTransform(trace,canonical,"TEST_INITIAL_CANONICAL_KEY_ORDER"),
        "canonical identity coordinate ordering failed");
    const auto old=f.read();
    const std::map<KFKey,std::map<KFKey,double>> shear{
        {KFState::oneKey,{{KFState::oneKey,1}}},{f.a,{{f.a,1},{f.b,1}}},
        {f.b,{{f.b,1}}},{f.clock,{{f.clock,1}}}};
    check(f.state.applyStateTransform(trace,shear,"TEST_SAME_KEYS_NEW_COORDINATES"),
        "valid invertible same-key shear failed");
    const auto current=f.read();
    check(current.numericBindingComplete && current.scope.stateOrderId==old.scope.stateOrderId &&
        current.scope.coordinateGeneration!=old.scope.coordinateGeneration,
        "test failed to exercise a same-order different-generation target");
    const auto before=f.state.metaDataMap;
    std::string reason;
    check(!zhangRegisterGraphProductHandoffTarget(f.state,old,reason) && f.state.metaDataMap==before,
        "old complete numeric request was silently relabeled with the current coordinate generation");
}
static void destinationHintsAreCandidatesOnly()
{
    std::ostringstream trace;
    Fixture exact;
    const std::map<std::string,std::map<KFKey,double>> correct{
        {"active-whole-phase",{{exact.a,1}}}};
    check(exact.state.applyStateTransform(trace,exact.sumTransform(),ZHANG_FLOAT_TARGET_COMPRESSION_LABEL,
        {},nullptr,&correct) && exact.read().numericRow==std::map<KFKey,double>{{exact.a,1}},
        "correct destination candidate did not pass its actual exact source pullback");
    for(const int mode:{0,1,2})
    {
        Fixture f;
        std::map<std::string,std::map<KFKey,double>> bad;
        if(mode==0)bad["active-whole-phase"][f.a]=2;
        if(mode==1)bad["active-whole-phase"][f.b]=1; // b is absent from destination.
        if(mode==2)bad["UNREGISTERED_TARGET"][f.a]=1;
        const auto before=f.state;
        check(!f.state.applyStateTransform(trace,f.sumTransform(),ZHANG_FLOAT_TARGET_COMPRESSION_LABEL,
            {},nullptr,&bad),"false candidate coefficient/key/target ID became source evidence");
        unchanged(f.state,before);
    }
}
int main()
{
    realTransform();rejectedTransactions();metadataOverride();timeTransition();tinyNonzeroIsNotExact();
    staleSameOrderRegistration();
    destinationHintsAreCandidatesOnly();
    std::cout<<"PASS real_central_float_target_metadata full_row_sum=1 full_clock_crossQ=1 "
        "strict_independent_target_rejected=1 strict_unbound_rejected=1 wrong_owner_rejected=1 "
        "callback_rollback_raw_exact=1 staged_raw_and_V2_atomic=1 index_insertion=1 "
        "random_prediction_not_exact=1 actual_process_Q_verified=1 tiny_nonzero_not_rounded_to_exact=1 "
        "same_keys_stale_registration_rejected=1 candidate_hints_not_authority=1\n";
}
