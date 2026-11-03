// Exercise the real central KF transaction, not a second remap implementation.
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductConditionalMap.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "pea/zhangReference.hpp"
#include <cstring>
#include <cstdlib>
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
    check(!near.state.applyStateTransform(trace,transform,"TEST_TINY_NONZERO_LEGITIMATE_KF_MAP"),
        "ordinary real transform silently unbound an active old target");
    unchanged(near.state,before);
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
    check(!noisy.state.applyStateTransform(trace,stochastic,"TEST_TINY_FRESH_PRIOR",{{source,1}}),
        "ordinary fresh-noise map silently unbound an active old target");
    unchanged(noisy.state,noisyBefore);
}
static void graphDirectionProtection()
{
    Fixture f;
    f.a=zhangTransportReceiverKey(E_Sys::BDS,E_ObsCode::L5P,"FALK");
    f.b=zhangTransportReceiverKey(E_Sys::BDS,E_ObsCode::L5P,"AGGO");
    f.clock.Sat=SatSys("C27");
    f.state.kfIndexMap={{KFState::oneKey,0},{f.a,1},{f.b,2},{f.clock,3}};
    f.state.stateTransitionMap.clear();
    for(const auto& [key,index]:f.state.kfIndexMap)f.state.stateTransitionMap[key][key][0]=1;
    f.target({{f.a,1},{f.b,-1}});
    auto target=f.read();target.numericOffset=0;
    target.physicalCoefficients.clear();
    target.physicalCoefficients[{{"FALK",SatSys("C28")},E_ObsCode::L5P,0}]=1;
    target.physicalCoefficients[{{"AGGO",SatSys("C28")},E_ObsCode::L5P,0}]=-1;
    // Both old arcs share the satellite potential, so their complete
    // difference is r_FALK-r_AGGO. C28 L1's topological edge remains current;
    // no missing L5 individual coordinate is manufactured from it.
    const ZhangGraphEdge falk{"FALK",SatSys("C28")};
    ZhangGraphBasis currentTopology;currentTopology.edges.insert(falk);
    const std::map<int,std::set<ZhangGraphEdge>> currentSignals{{int(E_ObsCode::L1P),{falk}},{int(E_ObsCode::L5P),{}}};
    check(zhangGraphDestinationSignalPresent(currentTopology,E_ObsCode::L1P,falk,&currentSignals) &&
        !zhangGraphDestinationSignalPresent(currentTopology,E_ObsCode::L5P,falk,&currentSignals),
        "FALK fixture did not model L5 missing with the L1 edge still present");
    f.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
    auto duplicate=target;duplicate.targetId="dependent-target";
    duplicate.numericRow={{f.a,2},{f.b,-2}};
    for(auto& [term,value]:duplicate.physicalCoefficients)value*=2;
    f.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+duplicate.targetId]=zhangEncodeFloatTarget(duplicate);
    ZhangGraphCoordinatePlan plan;plan.valid=true;
    plan.transform={{KFState::oneKey,{{KFState::oneKey,1}}},{f.clock,{{f.clock,1}}}};
    const auto before=f.state;
    auto metadata=f.state.metaDataMap;
    std::map<std::string,ZhangCoordinateRow> hints;std::string reason;std::ostringstream trace;
    check(zhangProtectGraphFloatTargets(f.state,plan,metadata,hints,reason),reason.c_str());
    check(plan.transform.size()==3,"dependent old targets did not share one independent carrier");
    check(f.state.metaDataMap==before.metaDataMap && (f.state.P-before.P).norm()==0,
        "graph target preparation mutated source metadata or covariance");
    check(f.state.applyStateTransform(trace,plan.transform,"TEST_PRE_GRAPH_MINIMAL_CARRIER",{},&metadata,&hints),
        "protected graph transform failed the central exact [T G] remap");
    std::vector<ZhangProductFloatTargetRequest> targets;
    check(zhangReadFloatTargets(f.state.metaDataMap,targets,reason) && targets.size()==2,
        "protected target directory changed identities");
    const auto& mapped=targets.front().numericRow;
    check(mapped.size()==1 && zhangIsProductFloatCarrierKey(mapped.begin()->first),
        "lost old direction was not a typed real carrier");
    const int carrier=f.state.kfIndexMap.at(mapped.begin()->first);
    check(f.state.x(carrier)==-3 && f.state.P(carrier,carrier)==11 &&
        f.state.P(carrier,f.state.kfIndexMap.at(f.clock))==-1,
        "carrier lost literal full mean/variance/clock cross-Q");
    auto query=target;query.targetId="different-consumer-complete-query";
    ZhangCoordinateRow resolved;double offset;
    check(zhangResolveGraphProductHandoffTarget(f.state,query,resolved,offset,reason) && offset==0,
        "whole old target could not resolve without individual physical rows");
    ZhangProductIntegerTransition transition;transition.valid=true;
    transition.observable=E_ObsCode::L5P;
    transition.physicalEdges={{"FALK",SatSys("C28")},{"AGGO",SatSys("C28")}};
    transition.physicalArcVersions={0,0};transition.coefficients={1,-1};
    auto writer=resolved;zhangAddCoordinateRow(writer,resolved);
    const auto numericMap=zhangProveWholeProductConditionalMap(
        resolved,writer,transition,.25,true,true,resolved,0);
    check(numericMap.valid && numericMap.physicalArcVersions==std::vector<int>({0,0}),
        "complete typed-carrier bridge lost exact physical identity");
    check(!zhangProveWholeProductConditionalMap(resolved,writer,transition,.25,false,true,resolved,0).valid &&
        !zhangProveWholeProductConditionalMap(resolved,writer,transition,.25,true,false,resolved,0).valid &&
        !zhangProveWholeProductConditionalMap(resolved,resolved,transition,.25,true,true,resolved,0).valid,
        "whole bridge skipped joint-Q, exact span or independent writer-row identity gates");
    auto unavailable=query;unavailable.physicalCoefficients.erase({{"AGGO",SatSys("C28")},E_ObsCode::L5P,0});
    check(!zhangResolveGraphProductHandoffTarget(f.state,unavailable,resolved,offset,reason),
        "one compressed old arc was fabricated from its difference carrier");
    auto wrongVersion=query;
    wrongVersion.physicalCoefficients.clear();
    wrongVersion.physicalCoefficients[{{"FALK",SatSys("C28")},E_ObsCode::L5P,1}]=1;
    wrongVersion.physicalCoefficients[{{"AGGO",SatSys("C28")},E_ObsCode::L5P,0}]=-1;
    check(!zhangResolveGraphProductHandoffTarget(f.state,wrongVersion,resolved,offset,reason),
        "V1 physical slip inherited a V0 complete target");
    check(zhangProductFloatTargetMetadataValid(f.state,reason),reason.c_str());
    const auto integerOwners=zhangCurrentEpochEvidenceRegistry().size();
    const auto start=f.state.time;
    f.state.advanced_postfits=false;
    f.state.prefitOpts.sigma_check=f.state.prefitOpts.omega_test=false;
    f.state.postfitOpts.sigma_check=f.state.postfitOpts.omega_test=false;
    // 901 seconds is longer than the inactive RAW cache's hard 300 seconds.
    // The active sufficient statistic must keep receiving posterior updates;
    // a frozen disk mean/variance would fail the independent cross-Q oracle.
    for(const int seconds:{30,90,180,600,1,899,1800})
    {
        check(f.state.stateTransition(trace,f.state.time+seconds),"active carrier prediction after long gap failed");
        const int live=f.state.kfIndexMap.at(f.clock);
        const auto priorMean=f.state.x;const auto priorQ=f.state.P;
        const VectorXd cross=priorQ.col(live);
        const double innovation=.75,variance=priorQ(live,live)+2;
        const VectorXd oracleMean=priorMean+cross*(innovation/variance);
        const MatrixXd oracleQ=priorQ-cross*cross.transpose()/variance;
        KFMeas measurement;measurement.H=MatrixXd::Zero(1,f.state.x.size());measurement.H(0,live)=1;
        measurement.R=MatrixXd::Constant(1,1,2);measurement.V=VectorXd::Constant(1,innovation);
        measurement.Y=measurement.H*f.state.x+measurement.V;measurement.VV=measurement.V;
        measurement.H_star=MatrixXd::Zero(1,0);measurement.uncorrelatedNoise.resize(0);
        measurement.prefitRatios=VectorXd::Zero(1);measurement.postfitRatios=VectorXd::Zero(1);
        measurement.obsKeys={f.clock};measurement.metaDataMaps.resize(1);measurement.componentsMaps.resize(1);
        check(f.state.filterKalman(trace,measurement,"/PPP",true)==KFFilterResult::COMMITTED,
            "long-gap production FLOAT update failed");
        check((f.state.x-oracleMean).norm()<1e-11 && (f.state.P-oracleQ).norm()<1e-11,
            "long-gap carrier full posterior differs from rank-one independent oracle");
        check(retireZhangGraphProductEventAuxiliaries(trace,f.state,"ACTIVE_CARRIER_901_SECONDS"),
            "active target was charged to inactive expiry or cache128");
        check(zhangResolveGraphProductHandoffTarget(f.state,query,resolved,offset,reason),
            "long-gap active complete target lost its current joint posterior");
    }
    check((f.state.time-start).to_double()==3600 && zhangCurrentEpochEvidenceRegistry().size()==integerOwners,
        "long-gap FLOAT retention created current integer authority or shortened the gap");
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
static void measuredSameEpochReceiptRebase()
{
    Fixture f;std::ostringstream trace;std::string reason;
    auto published=f.read();published.consumer="PUBLISHED_PHASE";published.numericOffset=0;
    ZhangProductIntegerFunctional identity;
    identity.valid=true;identity.satellite=SatSys("G01");identity.referenceSatellite=SatSys("G02");
    identity.observable=E_ObsCode::L1C;identity.physicalEdges={{"A",SatSys("G01")},{"A",SatSys("G02")}};
    identity.physicalArcVersions={0,0};identity.networkCoefficients={1,1};
    const auto fingerprint=zhangProductPhysicalFunctionalFingerprint(identity);
    published.targetId="PUBLISHED_PHASE|G01|L1C|literal-initial-target";
    published.lifecycleId="published-segment-A|"+fingerprint;
    f.state.metaDataMap.erase(std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+f.read().targetId);
    f.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+published.targetId]=zhangEncodeFloatTarget(published);
    ZhangGraphProductEventReceipt first;
    first.time=f.state.time;first.system=E_Sys::GPS;first.observable=E_ObsCode::L1C;first.satellite=SatSys("G01");
    first.productReferenceSatellite=SatSys("G02");first.oldProductFunctionalEdges=identity.physicalEdges;
    first.oldProductFunctionalArcVersions=identity.physicalArcVersions;
    first.oldProductFunctionalCoefficients=identity.networkCoefficients;
    first.physicalPathEdges=identity.physicalEdges;first.oldPhysicalArcVersions={0,0};first.physicalPathCoefficients={1,1};
    first.sourceRuntimeId=zhangProductFloatOwnerId(f.state);
    first.sourceStateOrderId=zhangFloatCarrierStateOrder(f.state.kfIndexMap);
    first.sourceFactorSequence=f.state.factorCommitSequence;
    first.oldBackendGeneration=10;first.newBackendGeneration=11;
    first.oldChartGeneration=20;first.newChartGeneration=21;
    first.fullJointTransportId="actual-first-exact-map";
    first.initialPublishedTargetId=published.targetId;first.initialPublishedLifecycleId=published.lifecycleId;
    first.initialPublishedFunctionalFingerprint=fingerprint;
    first.oldPhaseRow=published.numericRow;
    check(f.state.applyStateTransform(trace,f.sumTransform(),"FIRST_PREGRAPH_EXACT_EVENT"),
        "first actual graph-coordinate transform failed");
    first.destinationFactorSequence=f.state.factorCommitSequence;
    first.destinationStateOrderId=zhangFloatCarrierStateOrder(f.state.kfIndexMap);
    first.oldProductProjectionCaptured=first.representableSameTimeMapProven=first.fullJointMomentTransportProven=true;
    first.newPhaseRow={{f.a,1}};
    check(f.state.x(f.state.kfIndexMap.at(f.a))==7 && f.state.P(f.state.kfIndexMap.at(f.a),f.state.kfIndexMap.at(f.clock))==5,
        "first event did not preserve independent full mean/cross-Q oracle");
    f.state.prefitOpts.sigma_check=f.state.prefitOpts.omega_test=false;
    f.state.postfitOpts.sigma_check=f.state.postfitOpts.omega_test=false;
    f.state.advanced_postfits=false;
    KFMeas measurement;measurement.H=MatrixXd::Zero(1,f.state.x.size());measurement.H(0,f.state.kfIndexMap.at(f.clock))=1;
    measurement.R=MatrixXd::Constant(1,1,2);measurement.V=VectorXd::Constant(1,1.5);
    measurement.Y=measurement.H*f.state.x+measurement.V;measurement.VV=measurement.V;
    measurement.H_star=MatrixXd::Zero(1,0);measurement.uncorrelatedNoise.resize(0);
    measurement.prefitRatios=measurement.postfitRatios=VectorXd::Zero(1);
    measurement.obsKeys={f.clock};measurement.metaDataMaps.resize(1);measurement.componentsMaps.resize(1);
    check(f.state.filterKalman(trace,measurement,"/PPP",true)==KFFilterResult::COMMITTED,
        "same-epoch real FLOAT QC transaction failed");
    const auto current=f.read();
    // A deliberately stale pre-QC moment must not be carried into the second
    // event. The registered complete target identifies its current posterior.
    first.oldPhaseMeanMetres=-1e8;first.oldPhaseVarianceMetres2=1e12;
    ZhangGraphProductEventReceipt rebased;
    check(zhangRebaseGraphReceiptOnCurrentFloatTarget(f.state,first,current,11,21,
        f.state.factorCommitSequence,21,rebased,reason),reason.c_str());
    const int phase=f.state.kfIndexMap.at(f.a);
    check(rebased.oldBackendGeneration==10 && rebased.oldChartGeneration==20 &&
        rebased.oldPhaseMeanMetres==f.state.x(phase) && rebased.oldPhaseVarianceMetres2==f.state.P(phase,phase) &&
        rebased.initialPublishedTargetId==first.initialPublishedTargetId && rebased.sameEpochSourceRebased,
        "same-epoch rebase reused stale moments or lost the initial published endpoint");
    auto bad=first;bad.initialPublishedLifecycleId="unrelated-segment";
    check(!zhangRebaseGraphReceiptOnCurrentFloatTarget(f.state,bad,current,11,21,
        f.state.factorCommitSequence,21,rebased,reason),"wrong initial target lifecycle passed source rebase");
    check(!zhangRebaseGraphReceiptOnCurrentFloatTarget(f.state,first,current,11,21,
        f.state.factorCommitSequence-1,21,rebased,reason),"stale measured factor sequence passed source rebase");
    bad=first;bad.newBackendGeneration=12;
    check(!zhangRebaseGraphReceiptOnCurrentFloatTarget(f.state,bad,current,11,21,
        f.state.factorCommitSequence,21,rebased,reason),"unrelated first destination backend passed composition");
    check(zhangRebaseGraphReceiptOnCurrentFloatTarget(f.state,first,current,11,21,
        f.state.factorCommitSequence,21,rebased,reason),reason.c_str());
    ZhangGraphCoordinatePlan second;second.valid=true;
    second.transform={{KFState::oneKey,{{KFState::oneKey,1}}},{f.a,{{f.a,-1}}},{f.clock,{{f.clock,1}}}};
    ZhangCoordinateRow mapped;
    check(zhangResolveGraphSourceFunctional(f.state,second,rebased.oldPhaseRow,mapped,reason) &&
        zhangVerifyGraphFunctionalTransport(second,rebased.oldPhaseRow,mapped,0,0,0).sameTimeMapProven,
        "second actual T/G could not carry initial target from the measured source");
    const auto measured=f.state;
    check(f.state.applyStateTransform(trace,second.transform,"SECOND_POSTQC_EXACT_EVENT"),
        "second exact graph-coordinate transaction failed");
    for(const auto& [key,index]:f.state.kfIndexMap)
    {
        double actual=0,expected=0;
        for(const auto& [left,a]:mapped)actual+=a*f.state.P(f.state.kfIndexMap.at(left),index);
        for(const auto& [left,a]:rebased.oldPhaseRow)for(const auto& [right,b]:second.transform.at(key))
            expected+=a*b*measured.P(measured.kfIndexMap.at(left),measured.kfIndexMap.at(right));
        check(std::abs(actual-expected)<1e-12,"composed target lost full destination cross-Q");
    }
}
int main()
{
    setenv("ZHANG_R51_ENABLE","1",1);
    setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
    realTransform();rejectedTransactions();metadataOverride();timeTransition();tinyNonzeroIsNotExact();
    staleSameOrderRegistration();
    destinationHintsAreCandidatesOnly();
    graphDirectionProtection();
    measuredSameEpochReceiptRebase();
    std::cout<<"PASS real_central_float_target_metadata full_row_sum=1 full_clock_crossQ=1 "
        "strict_independent_target_rejected=1 strict_unbound_rejected=1 wrong_owner_rejected=1 "
        "callback_rollback_raw_exact=1 staged_raw_and_V2_atomic=1 index_insertion=1 "
        "random_prediction_not_exact=1 actual_process_Q_verified=1 tiny_nonzero_not_rounded_to_exact=1 "
        "same_keys_stale_registration_rejected=1 candidate_hints_not_authority=1 "
        "ordinary_active_loss_rejected=1 minimal_pre_graph_carrier=1 typed_whole_target_resolved=1 "
        "compressed_individual_arc_unavailable=1 physical_V1_negative=1 "
        "active_carrier_30_120_300_900_901_1800_3600_seconds_fullQ=1 no_gap_AR_authority=1 "
        "same_epoch_pregraph_QC_second_transform=1 stale_moments_recomputed=1 wrong_factor_or_initial_target_rejected=1\n";
}
