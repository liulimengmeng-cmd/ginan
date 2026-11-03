// Real core archive encode/decode + real restore implementation, no PEA.
// This is a prerequisite to whole-process nine-section checkpoint replay.
#include "common/zhangCheckpoint.hpp"
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

static unsigned checks=0;
static void require(bool value,const std::string& reason)
{++checks;if(!value)throw std::runtime_error(reason);}
static ZhangCheckpointKfCore nativeWireRoundtrip(const ZhangCheckpointKfCore& core)
{
    const auto raw=serializeZhangCheckpointSectionPayload(core);
    std::istringstream in(raw,std::ios::binary|std::ios::in);
    ZhangCheckpointKfCore decoded;
    {boost::archive::binary_iarchive archive(in,boost::archive::no_header);archive>>decoded;}
    require(in.peek()==std::char_traits<char>::eof(),"core archive trailing bytes");
    require(serializeZhangCheckpointSectionPayload(decoded)==raw,"core native wire roundtrip changed bytes");
    return decoded;
}
struct Fixture
{
    KFState state;KFKey a,b,clock,carrier;
    Fixture()
    {
        a.type=KF::PHASE_BIAS;a.Sat=SatSys("G01");a.num=int(E_ObsCode::L1C);
        b=a;b.Sat=SatSys("G02");clock.type=KF::SAT_CLOCK;clock.Sat=SatSys("G01");
        carrier.type=KF::PRODUCT_HANDOFF;carrier.str="ZHANG_FLOAT_CARRIER:CHECKPOINT_SUM";
        state.time.bigTime=1000;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="FLOAT_TARGET_CHECKPOINT_COUNTER";
        state.kfIndexMap={{KFState::oneKey,0},{a,1},{b,2},{clock,3},{carrier,4}};
        state.x.resize(5);state.x<<1,2,5,10,7;state.dx=VectorXd::Zero(5);
        state.P.resize(5,5);
        state.P<<0,0,0,0,0, 0,4,1,2,5, 0,1,9,3,10, 0,2,3,16,5, 0,5,10,5,15;
        for(const auto&[key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
        state.procNoiseMap[carrier]=0;state.factorCommitSequence=42;
        ZhangProductFloatTargetRequest target;
        target.targetId="whole-phase-sum";target.consumer="CHECKPOINT_CONTRACT";
        target.lifecycleId="PUBLIC_SEGMENT_1";target.unit="m";
        target.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
        target.physicalCoefficients[{{"A",SatSys("G01")},E_ObsCode::L1C,0}]=1;
        target.physicalCoefficients[{{"A",SatSys("G02")},E_ObsCode::L1C,0}]=1;
        target.scope={zhangProductFloatOwnerId(state),"EXPLICIT_CHECKPOINT_CHART",
            zhangFloatCarrierStateOrder(state.kfIndexMap),37};
        target.numericRow={{carrier,1}};target.numericOffset=3.25;
        target.active=true;target.numericBindingComplete=true;target.bindingReason="LITERAL_CARRIER_BINDING";
        state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
        auto definition=target;definition.targetId="sum-carrier-definition";
        definition.consumer="FLOAT_CARRIER_SOURCE_DEFINITION";definition.numericRow={{a,1},{b,1}};
        definition.numericOffset=0;
        state.metaDataMap[zhangProductFloatCarrierRecordKey(carrier)]=zhangEncodeFloatTarget(definition);
        state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]=target.scope.stateOrderId;
        state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]="37";
        state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]="42";
    }
};
static void positiveAndAtomicNegatives()
{
    Fixture f;std::string reason;
    require(zhangProductFloatTargetMetadataValid(f.state,reason),"source V2 fixture invalid:"+reason);
    const auto source=nativeWireRoundtrip(captureZhangCheckpointKfCore(f.state));
    KFState restored;restored.factorCommitSequence=1;
    require(restoreZhangCheckpointKfCore(source,restored,&reason),"valid nonempty core restore rejected:"+reason);
    require(restored.factorCommitSequence==42,"saved true factor sequence not restored");
    require(restored.metaDataMap.at(ZHANG_FLOAT_TARGET_GENERATION_V2)=="37",
        "restore changed coordinate generation into the later measurement sequence");
    require(zhangProductFloatTargetMetadataValid(restored,reason),"restored target binding invalid:"+reason);
    require(zhangCheckpointKfCoreBitwiseEqual(source,captureZhangCheckpointKfCore(restored)),
        "nonempty V2 restore changed any of the original 33 core fields");
    auto reject=[&](const std::function<void(ZhangCheckpointKfCore&)>& mutate){
        auto bad=source;mutate(bad);bad=nativeWireRoundtrip(bad);
        KFState destination;destination.factorCommitSequence=918;
        destination.metaDataMap["DO_NOT_CHANGE"]="REJECTED_IMPORT_SENTINEL";
        const auto before=captureZhangCheckpointKfCore(destination);
        require(!restoreZhangCheckpointKfCore(bad,destination,&reason),"invalid V2 counter/binding accepted");
        require(destination.factorCommitSequence==918 && zhangCheckpointKfCoreBitwiseEqual(before,captureZhangCheckpointKfCore(destination)),
            "rejected V2 restore partially changed destination core or factor sequence");
    };
    reject([](auto&x){x.metaDataMap.erase(ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2);});
    for(const std::string value:{"","-1","+42"," 42","042","18446744073709551616","42junk"})
        reject([&](auto&x){x.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=value;});
    reject([](auto&x){x.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]="36";});
    reject([](auto&x){x.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]="FORGED_FULL_ORDER";});
    reject([](auto&x){x.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="OTHER_AUTHORITATIVE_OWNER";});
    reject([](auto&x){x.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]="43";});
    reject([](auto&x){
        for(auto it=x.metaDataMap.begin();it!=x.metaDataMap.end();)
            if(it->first.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || it->first.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX))
                it=x.metaDataMap.erase(it);else ++it;
    });
}
static void emptyLegacyAddsNothing()
{
    KFState legacy;legacy.factorCommitSequence=123;
    legacy.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="E29-GEC134-FIXEDCOORD-20261023-one_hour";
    const auto metadata=legacy.metaDataMap;
    const auto core=nativeWireRoundtrip(captureZhangCheckpointKfCore(legacy));
    const auto bytes=serializeZhangCheckpointSectionPayload(core);
    require(!hasZhangFloatTargetV2State(core.metaDataMap),"empty legacy fixture unexpectedly has V2 schema");
    KFState restored;restored.factorCommitSequence=777;std::string reason;
    require(restoreZhangCheckpointKfCore(core,restored,&reason),"legacy empty core rejected:"+reason);
    require(restored.factorCommitSequence==777,"legacy restore changed nonserialized runtime counter");
    require(restored.metaDataMap==metadata && !hasZhangFloatTargetV2State(restored.metaDataMap),"legacy restore added V2 tag");
    require(serializeZhangCheckpointSectionPayload(captureZhangCheckpointKfCore(restored))==bytes,
        "legacy empty 33-field core wire changed");
}
int main()
{
    try {positiveAndAtomicNegatives();emptyLegacyAddsNothing();
        std::cout<<"CHECKPOINT_FLOAT_TARGET_COUNTER_ROUNDTRIP checks="<<checks
            <<" passed=1 true_saved_sequence=1 coordinate_generation_unchanged=1 raw_core_fields=33 empty_legacy_no_tag=1 whole_process_restore_executed=0\n";return 0;}
    catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
