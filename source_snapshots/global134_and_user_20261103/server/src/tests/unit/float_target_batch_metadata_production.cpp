// Calls the unchanged production single Register/Resolve as the golden path.
// Synthetic metadata/kernel contract; never executes PEA main or a GNSS run.
#include "common/acsConfig.hpp"
#include "common/zhangProductFloatTargetBatch.hpp"
#include "common/zhangProductRelationAdmission.hpp"
#include "pea/zhangReference.hpp"
#include "ambres/GNSSambres.hpp"
#include <openssl/sha.h>
#include <chrono>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

static void check(bool ok,const std::string& why)
{if(!ok)throw std::runtime_error(why);}
static std::string metadataBatchRawByteSha(const void* data,std::size_t bytes)
{
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(static_cast<const unsigned char*>(data),bytes,digest);
    std::ostringstream out;for(auto c:digest)out<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c);
    return out.str();
}
static std::string hashKV(const std::map<std::string,std::string>& metadata)
{
    SHA256_CTX ctx;SHA256_Init(&ctx);
    for(const auto& [key,value]:metadata)
    {
        for(const auto* text:{&key,&value})
        {
            const std::uint64_t n=text->size();SHA256_Update(&ctx,&n,sizeof(n));
            SHA256_Update(&ctx,text->data(),text->size());
        }
    }
    unsigned char digest[SHA256_DIGEST_LENGTH];SHA256_Final(digest,&ctx);
    std::ostringstream out;for(auto c:digest)out<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c);
    return out.str();
}
struct Fixture
{
    KFState state;
    KFKey a,b,clock;
    ZhangGraphEdge ea{"META_A",SatSys("G01")},eb{"META_B",SatSys("G02")};
    explicit Fixture(int n=8)
    {
        check(n>=4,"fixture dimension");
        state.time.bigTime=1405216800;state.factorCommitSequence=17;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="METADATA_BATCH_SYNTHETIC_OWNER";
        a=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,ea,0);
        b=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,eb,0);
        clock.type=KF::SAT_CLOCK;clock.Sat=SatSys("G01");
        state.kfIndexMap={{KFState::oneKey,0},{a,1},{b,2},{clock,3}};
        for(int i=4;i<n;++i)
        {
            KFKey key;key.type=KF::AMBIGUITY;key.str="METADATA_STATION_"+std::to_string(i);
            key.Sat=SatSys(E_Sys::GPS,1+i%32);key.num=int(E_ObsCode::L1C);
            state.kfIndexMap.emplace(key,i);
        }
        state.x=VectorXd::LinSpaced(n,-.2,.4);state.x(0)=1;
        state.dx=VectorXd::LinSpaced(n,.001,.002);state.dx(0)=0;
        const VectorXd u=VectorXd::LinSpaced(n,.001,.002);
        state.P=MatrixXd::Identity(n,n)+u*u.transpose();
        state.P.row(0).setZero();state.P.col(0).setZero();
        for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
        check(zhangHandoffFloatRegisterRecord(state,a,ea,int(E_ObsCode::L1C),0),"RAW a registry");
        check(zhangHandoffFloatRegisterRecord(state,b,eb,int(E_ObsCode::L1C),0),"RAW b registry");
    }
    ZhangProductFloatTargetRequest target(const std::string& id,bool bound=true)const
    {
        ZhangProductFloatTargetRequest t;t.targetId=id;t.consumer="SYNTHETIC_COMPLETE_PHASE";
        t.lifecycleId="FIXED_LITERAL_LIFECYCLE";t.unit="m";
        t.scope={zhangProductFloatOwnerId(state),"LITERAL_COMPLETE_GAUSSIAN_CHART",
            zhangFloatCarrierStateOrder(state.kfIndexMap),generation()};
        t.numericRow={{a,2},{b,-1},{clock,.125}};t.numericOffset=.5;
        t.numericBindingComplete=bound;t.bindingReason="SYNTHETIC_SOURCE_ROW";
        return t;
    }
    std::uint64_t generation()const
    {
        const auto i=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);
        return i==state.metaDataMap.end()?state.factorCommitSequence:std::stoull(i->second);
    }
    void seed(std::size_t count)
    {
        const auto order=zhangFloatCarrierStateOrder(state.kfIndexMap);
        for(std::size_t i=0;i<count;++i)
        {
            std::ostringstream id;id<<"TARGET_"<<std::setw(5)<<std::setfill('0')<<i;
            auto t=target(id.str());t.active=(i%2==0);
            state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+t.targetId]=zhangEncodeFloatTarget(t);
        }
        if(count)
        {
            state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]=order;
            state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(state.factorCommitSequence);
            state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(state.factorCommitSequence);
        }
        std::string reason;check(zhangProductFloatTargetMetadataValid(state,reason),"seed:"+reason);
    }
    std::string moments()const
    {
        return metadataBatchRawByteSha(state.x.data(),state.x.size()*sizeof(double))+"|"+
            metadataBatchRawByteSha(state.P.data(),state.P.size()*sizeof(double))+"|"+
            metadataBatchRawByteSha(state.dx.data(),state.dx.size()*sizeof(double));
    }
};
static bool golden(Fixture& f,const std::vector<ZhangProductFloatTargetRequest>& requests,std::string& reason)
{
    try {zhangClearGraphProductHandoffTargets(f.state);}
    catch(const std::exception& error){reason=error.what();return false;}
    for(const auto& t:requests)if(!zhangRegisterGraphProductHandoffTarget(f.state,t,reason))return false;
    return true;
}
static std::string integerArchive(const Fixture& f)
{
    std::string wire;const auto result=exportZhangAmbresCheckpointSection(
        f.state,zhangProductFloatOwnerId(f.state),wire);
    check(result.valid,"integer archive export:"+result.failureReason);return wire;
}
static int contractCases=0;
static void equivalent(Fixture& f,const std::vector<ZhangProductFloatTargetRequest>& requests,const std::string& label)
{
    const auto before=f.state.metaDataMap;const auto moments=f.moments();const auto sequence=f.state.factorCommitSequence;
    const auto archive=integerArchive(f);std::string reason;
    check(golden(f,requests,reason),label+": old failed:"+reason);
    check(zhangProductFloatTargetMetadataValid(f.state,reason),label+": old final invalid:"+reason);
    const auto expected=f.state.metaDataMap;
    f.state.metaDataMap=before;ZhangFloatTargetBatchStats stats;
    check(zhangRefreshGraphProductHandoffTargetsBatch(f.state,requests,reason,&stats),label+": new failed:"+reason);
    check(f.state.metaDataMap==expected,label+": raw KV mismatch");
    check(zhangProductFloatTargetMetadataValid(f.state,reason),label+": original validator rejected new:"+reason);
    check(f.moments()==moments && f.state.factorCommitSequence==sequence,label+": x/P/dx/sequence changed");
    check(integerArchive(f)==archive,label+": nonempty integer archive changed");
    check(stats.atomicSwaps==1 && stats.sourceValidatorCalls==1 && stats.sourceDirectoryDecodeCalls==1 && stats.finalValidatorCalls==1,
        label+": nonconstant directory transaction count");
    ++contractCases;std::cout<<"CONTRACT_PASS "<<label<<" kv_sha="<<hashKV(expected)<<"\n";
}
static void rejectedAtomic(Fixture& f,const std::vector<ZhangProductFloatTargetRequest>& requests,const std::string& label)
{
    const auto before=f.state.metaDataMap;const auto moments=f.moments();const auto archive=integerArchive(f);
    std::string reason;ZhangFloatTargetBatchStats stats;
    check(!zhangRefreshGraphProductHandoffTargetsBatch(f.state,requests,reason,&stats),label+": wrongly accepted");
    check(f.state.metaDataMap==before && f.moments()==moments && integerArchive(f)==archive && stats.atomicSwaps==0,
        label+": failed transaction polluted original");
    ++contractCases;std::cout<<"CONTRACT_PASS "<<label<<" rejection="<<reason<<"\n";
}
static void contracts()
{
    Fixture f;
    equivalent(f,{},"legacy empty directory");
    equivalent(f,{f.target("B"),f.target("A")},"bound ordered raw metadata");
    equivalent(f,{f.target("B"),f.target("A")},"idempotent refresh");
    auto duplicate=f.target("B");duplicate.active=false;duplicate.numericOffset=1.25;
    equivalent(f,{f.target("B"),duplicate},"duplicate last write wins");
    auto unknown=f.target("UNKNOWN",false);unknown.symbolicSpace=ZhangProductFloatTargetSpace::SOURCE_KF_COORDINATE;
    KFKey missing;missing.type=KF::AMBIGUITY;missing.str="UNBOUND_UNKNOWN_COLUMN";missing.Sat=SatSys("G09");
    unknown.numericRow={{missing,2}};
    equivalent(f,{unknown},"unknown unbound column retained, not zero filled");
    auto bound=f.target("Z_PREVIOUS");bound.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
    ZhangProductFloatVersionedTerm fake{{"NOT_A_PHYSICAL_CURRENT_ARC",SatSys("G09")},E_ObsCode::L1C,3};
    bound.physicalCoefficients={{fake,1}};bound.numericRow={{f.a,2},{f.b,-1}};bound.numericOffset=.5;
    auto consequence=bound;consequence.targetId="A_AFTER_Z";consequence.numericBindingComplete=false;
    consequence.numericRow.clear();consequence.physicalCoefficients={{fake,3}};
    equivalent(f,{bound,consequence},"unbound sees prior ordered upsert");
    std::vector<ZhangProductFloatTargetRequest> parsed;std::string reason;
    check(zhangReadFloatTargets(f.state.metaDataMap,parsed,reason),reason);
    auto found=std::find_if(parsed.begin(),parsed.end(),[](const auto& t){return t.targetId=="A_AFTER_Z";});
    check(found!=parsed.end() && found->numericBindingComplete && found->numericRow.at(f.a)==6 &&
        found->numericRow.at(f.b)==-3 && found->numericOffset==1.5,"independent literal span/offset oracle");
    auto bad=f.target("FAIL");bad.scope.runtimeId="OTHER_OWNER";
    rejectedAtomic(f,{f.target("FIRST_VALID"),bad},"late stale owner atomic failure");
    bad=f.target("FAIL");bad.scope.stateOrderId="STALE_ORDER";rejectedAtomic(f,{bad},"stale full state order");
    bad=f.target("FAIL");bad.scope.coordinateGeneration=16;rejectedAtomic(f,{bad},"stale coordinate generation");
    bad=f.target("FAIL");bad.numericRow={{missing,1}};rejectedAtomic(f,{bad},"bound unavailable column");
    bad=f.target("FAIL");bad.numericRow={{f.clock,std::numeric_limits<double>::infinity()}};
    rejectedAtomic(f,{bad},"nonfinite bound full Gaussian row");
    bad=f.target("FAIL");bad.unit.clear();rejectedAtomic(f,{bad},"invalid unit");
    bad=f.target("FAIL",false);bad.numericRow={{missing,std::numeric_limits<double>::quiet_NaN()}};
    rejectedAtomic(f,{bad,f.target("FAIL")},"bad intermediate wire cannot be repaired by duplicate");
    {
        Fixture stale;stale.seed(2);stale.state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]="16";
        rejectedAtomic(stale,{stale.target("NEW")},"stale actual saved factor sequence");
    }
    {
        Fixture corrupt;corrupt.seed(2);corrupt.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+"TARGET_00000"]="BAD_WIRE";
        rejectedAtomic(corrupt,{},"invalid source wire");
    }
    {
        Fixture generation;generation.seed(2);std::vector<ZhangProductFloatTargetRequest> targets;
        check(zhangReadFloatTargets(generation.state.metaDataMap,targets,reason),reason);
        generation.state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]="9";
        for(auto& t:targets){t.scope.coordinateGeneration=9;generation.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+t.targetId]=zhangEncodeFloatTarget(t);}
        equivalent(generation,{generation.target("NEW")},"coordinate generation 9, actual factor counter 17");
    }
    {
        Fixture carrier;
        auto source=std::find_if(carrier.state.kfIndexMap.begin(),carrier.state.kfIndexMap.end(),[](const auto& item){return item.second==7;});
        const auto originalKey=source->first;carrier.state.kfIndexMap.erase(source);carrier.state.stateTransitionMap.erase(originalKey);
        KFKey key;key.type=KF::PRODUCT_HANDOFF;key.str="ZHANG_FLOAT_CARRIER:SYNTHETIC_COLUMN_ROLE";
        carrier.state.kfIndexMap.emplace(key,7);carrier.state.stateTransitionMap[key][key][0]=1;
        auto definition=carrier.target("CARRIER_SOURCE");definition.consumer="FLOAT_CARRIER_SOURCE_DEFINITION";
        definition.active=false;definition.numericBindingComplete=false;
        carrier.state.metaDataMap[zhangProductFloatCarrierRecordKey(key)]=zhangEncodeFloatTarget(definition);
        carrier.seed(1);
        auto target=carrier.target("CARRIER_CONSUMER");target.numericRow={{key,1},{carrier.clock,-1}};
        equivalent(carrier,{target},"typed FLOAT carrier role and complete clock covariance retained");
        carrier.state.stateTransitionMap[key][key][1]=1;
        rejectedAtomic(carrier,{target},"invalid nonconstant FLOAT carrier model");
    }
    {
        Fixture capacity;capacity.seed(8192);
        equivalent(capacity,{capacity.target("TARGET_00000")},"8192 legal boundary");
        const auto before=capacity.state.metaDataMap;auto extra=capacity.target("ZZZ_NEW_8193");
        check(golden(capacity,{extra},reason),"old lazy final insertion did not expose boundary fixture");
        check(!zhangProductFloatTargetMetadataValid(capacity.state,reason),"old 8193 unexpectedly passed final gate");
        capacity.state.metaDataMap=before;rejectedAtomic(capacity,{extra},"old final 8193 lazy hole rejected atomically");
    }
    std::cout<<"CONTRACT_SUMMARY passed="<<contractCases<<" old_register_resolve_validator_unchanged=1"
        <<" raw_KV_and_full_x_P_dx_equal=1 nonempty_integer_archive_unchanged=1 scientific_replay=0\n";
}
static void benchmark(const std::string& mode,int targets,int repeats)
{
    check(mode=="old" || mode=="batch","benchmark mode");check(targets>0 && targets<=308 && repeats==3,"bounded declared fixture sizes");
    Fixture f(7000);const int sourceTargets=targets*5/2;f.seed(sourceTargets);
    std::vector<ZhangProductFloatTargetRequest> requests;
    for(int i=0;i<targets;++i)
    {std::ostringstream id;id<<"TARGET_"<<std::setw(5)<<std::setfill('0')<<i;auto t=f.target(id.str());t.numericOffset=double(i)/16;requests.push_back(std::move(t));}
    const auto original=f.state.metaDataMap;const auto before=f.moments();const auto integerBefore=integerArchive(f);
    std::string expected;
    for(int repeat=0;repeat<repeats;++repeat)
    {
        f.state.metaDataMap=original;std::string reason;ZhangFloatTargetBatchStats stats;
        const auto begin=std::chrono::steady_clock::now();const auto cpu=std::clock();
        const bool ok=mode=="old"?golden(f,requests,reason):zhangRefreshGraphProductHandoffTargetsBatch(f.state,requests,reason,&stats);
        const double cpuSeconds=double(std::clock()-cpu)/CLOCKS_PER_SEC;
        const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        check(ok,"benchmark execution:"+reason);check(zhangProductFloatTargetMetadataValid(f.state,reason),"benchmark final gate:"+reason);
        const auto kv=hashKV(f.state.metaDataMap);if(expected.empty())expected=kv;else check(expected==kv,"repeat raw KV mismatch");
        check(f.moments()==before && integerArchive(f)==integerBefore,"benchmark changed moments or nonempty integer archive");
        std::cout<<std::setprecision(17)<<"BENCHMARK mode="<<mode<<" states=7000 source_targets="<<sourceTargets
            <<" requested_targets="<<targets<<" repeat="<<repeat<<" wall_seconds="<<wall<<" process_cpu_seconds="<<cpuSeconds
            <<" state_order_bytes="<<zhangFloatCarrierStateOrder(f.state.kfIndexMap).size()<<" kv_sha="<<kv
            <<" moment_raw_sha="<<before<<" integer_archive_sha="<<metadataBatchRawByteSha(integerBefore.data(),integerBefore.size())
            <<" atomic_swaps="<<stats.atomicSwaps<<" numerical_or_search_work=0"<<std::endl;
    }
}
int main(int argc,char** argv)
{
    try
    {
        setenv("OPENBLAS_NUM_THREADS","1",1);setenv("OMP_NUM_THREADS","1",1);
        setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
        setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);
        // Only this synthetic fixture's legacy archive export is enabled.
        // No real runner configuration or current AR authorization is changed.
        setenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT","0",1);
        auto& sentinel=zhangProductRelationAdmissionStateRegistry()[{"METADATA_BATCH_SYNTHETIC_OWNER",E_Sys::GPS}];
        sentinel.certifiedRows={{{"SYNTHETIC_INTEGER_NO_MUTATION_SENTINEL",1}}};sentinel.certifiedValues={7};
        sentinel.certifiedSatellites={"G13"};sentinel.certifiedRelationIds={"SYNTHETIC_ONLY"};
        if(argc==1)contracts();
        else {check(argc==4,"args: mode targets repeats");benchmark(argv[1],std::stoi(argv[2]),std::stoi(argv[3]));}
        return 0;
    }
    catch(const std::exception& e){std::cerr<<"METADATA_BATCH_DRIVER_FAILED: "<<e.what()<<"\n";return 1;}
}
