// The unchanged single-request production resolver is the golden implementation.
// This executable exercises coordinate/metadata contracts, not PEA scientific AR.
#include "common/acsConfig.hpp"
#include "common/zhangFloatTargetResolverBatch.hpp"
#include "common/zhangExactWorkBudget.hpp"
#include "common/zhangProductRelationAdmission.hpp"
#include "pea/zhangReference.hpp"
#include "pea/zhangPppAr.hpp"
#include "ambres/GNSSambres.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

static unsigned cases=0;
static void check(bool value,const std::string& why)
{if(!value)throw std::runtime_error(why);}
static bool sameDouble(double a,double b)
{return std::memcmp(&a,&b,sizeof(double))==0;}
static bool sameRow(const std::map<KFKey,double>& a,const std::map<KFKey,double>& b)
{
    if(a.size()!=b.size())return false;
    auto ai=a.begin(),bi=b.begin();
    for(;ai!=a.end();++ai,++bi)
        if(ai->first!=bi->first || !sameDouble(ai->second,bi->second))return false;
    return true;
}
struct Fixture
{
    KFState state;
    ZhangGraphEdge ea{"RESOLVE_A",SatSys("G01")},eb{"RESOLVE_B",SatSys("G02")};
    KFKey a,b,clock;
    Fixture()
    {
        state.time.bigTime=1405216800;state.factorCommitSequence=17;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="RESOLVER_BATCH_GOLDEN_OWNER";
        a=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,ea,0);
        b=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,eb,0);
        clock.type=KF::SAT_CLOCK;clock.Sat=SatSys("G01");
        state.kfIndexMap={{KFState::oneKey,0},{a,1},{b,2},{clock,3}};
        state.x.resize(4);state.x<<1,2,5,8;state.dx.resize(4);state.dx<<0,.25,-.125,.5;
        state.P.resize(4,4);state.P<<0,0,0,0,0,4,1,2,0,1,9,3,0,2,3,16;
        for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
        check(zhangHandoffFloatRegisterRecord(state,a,ea,int(E_ObsCode::L1C),0),"raw a registration");
        check(zhangHandoffFloatRegisterRecord(state,b,eb,int(E_ObsCode::L1C),0),"raw b registration");
    }
    ZhangProductFloatTargetRequest target(const std::string& id)const
    {
        ZhangProductFloatTargetRequest t;t.targetId=id;t.consumer="RESOLVER_CONTRACT";
        t.lifecycleId="LITERAL_SEGMENT";t.unit="m";t.numericBindingComplete=true;
        t.scope={zhangProductFloatOwnerId(state),"ACTUAL_LITERAL_CHART",
            zhangFloatCarrierStateOrder(state.kfIndexMap),9};
        t.numericRow={{a,2},{b,-1},{clock,.125}};t.numericOffset=.5;
        return t;
    }
    void seed(const std::vector<ZhangProductFloatTargetRequest>& targets)
    {
        for(const auto& t:targets)
            state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+t.targetId]=zhangEncodeFloatTarget(t);
        if(!targets.empty())
        {
            state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]=zhangFloatCarrierStateOrder(state.kfIndexMap);
            state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]="9";
            state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]="17";
        }
    }
};
static std::string integerArchive(const Fixture& f)
{
    std::string wire;const auto exported=exportZhangAmbresCheckpointSection(f.state,
        zhangProductFloatOwnerId(f.state),wire);
    check(exported.valid,"nonempty integer archive export:"+exported.failureReason);return wire;
}
static void compare(Fixture& f,const std::vector<ZhangProductFloatTargetRequest>& requests,
    const std::string& label,bool auditAuthoritativeIntegerArchive=true)
{
    const auto metadata=f.state.metaDataMap;
    const auto indices=f.state.kfIndexMap;
    const auto sequence=f.state.factorCommitSequence;
    const VectorXd x=f.state.x,dx=f.state.dx;const MatrixXd p=f.state.P;
    const auto archive=auditAuthoritativeIntegerArchive?integerArchive(f):std::string{};
    std::vector<ZhangFloatTargetResolution> golden;
    for(const auto& request:requests)
    {
        ZhangFloatTargetResolution answer;
        answer.resolved=zhangResolveGraphProductHandoffTarget(f.state,request,answer.row,answer.offset,answer.reason);
        golden.push_back(std::move(answer));
    }
    std::string validity;
    const bool valid=zhangProductFloatTargetMetadataValid(f.state,validity);
    std::vector<ZhangFloatTargetResolution> actual(1);
    std::vector<ZhangProductFloatTargetRequest> directory(1);
    ZhangFloatTargetResolverBatchStats stats;std::string reason;
    const bool accepted=zhangResolveGraphProductHandoffTargetsBatch(f.state,requests,actual,reason,&stats,&directory);
    check(accepted==valid,label+": batch/old directory validity differs:"+reason);
    if(valid)
    {
        check(actual.size()==golden.size(),label+": ordered output size changed");
        for(std::size_t i=0;i<golden.size();++i)
            check(actual[i].resolved==golden[i].resolved && actual[i].reason==golden[i].reason &&
                sameRow(actual[i].row,golden[i].row) && sameDouble(actual[i].offset,golden[i].offset),
                label+": per-target raw result differs at "+std::to_string(i));
        std::vector<ZhangProductFloatTargetRequest> oldDirectory;
        check(zhangReadFloatTargets(f.state.metaDataMap,oldDirectory,reason) && oldDirectory.size()==directory.size(),
            label+": validated directory differs");
        for(std::size_t i=0;i<directory.size();++i)
            check(zhangEncodeFloatTarget(directory[i])==zhangEncodeFloatTarget(oldDirectory[i]) &&
                directory[i].scope==oldDirectory[i].scope,label+": directory map order/raw meaning changed");
        check(stats.sourceValidatorCalls==1 && stats.explicitDirectoryDecodeCalls==1 && stats.stateOrderBuilds==1 &&
            stats.requestedTargets==requests.size(),label+": per-request directory rebuild returned");
    }
    else check(actual.empty() && directory.empty() && reason==validity,label+": fatal context was not atomic");
    check(f.state.metaDataMap==metadata && f.state.kfIndexMap==indices && f.state.factorCommitSequence==sequence &&
        !std::memcmp(f.state.x.data(),x.data(),x.size()*sizeof(double)) &&
        !std::memcmp(f.state.dx.data(),dx.data(),dx.size()*sizeof(double)) &&
        !std::memcmp(f.state.P.data(),p.data(),p.size()*sizeof(double)) &&
        (!auditAuthoritativeIntegerArchive || integerArchive(f)==archive),
        label+": source full x/P/dx/metadata/index/sequence/integer archive mutated");
    ++cases;std::cout<<"RESOLVER_BATCH_CONTRACT_PASS "<<label<<" targets="<<requests.size()
        <<" resolved="<<stats.resolvedTargets<<" unbound="<<stats.unboundTargets<<'\n';
}
static void tests()
{
    Fixture legacy;auto unknown=legacy.target("legacy-unsaved-anchor");
    compare(legacy,{},"empty legacy no V2 tag injection");
    compare(legacy,{unknown},"legacy absolute anchor stays UNKNOWN");
    Fixture f;auto bound=f.target("bound"),unbound=f.target("unbound");unbound.numericBindingComplete=false;
    f.seed({bound,unbound});
    auto alias=bound;alias.targetId="different-id-same-other-fields";
    auto changed=bound;changed.lifecycleId="OTHER_LIFECYCLE";
    compare(f,{bound,unbound,bound,alias,changed,bound},"registered unbound repeated requests alias full semantics");
    // Current physical version0 RAW directions, including a nonzero offset.
    auto raw=f.target("literal-physical");raw.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
    raw.physicalCoefficients[ {f.ea,E_ObsCode::L1C,0} ]=2;
    raw.physicalCoefficients[ {f.eb,E_ObsCode::L1C,0} ]=-3;raw.symbolicOffset=ZhangFloatCarrierRational(1)/4;
    auto absent=raw;absent.targetId="unknown-physical-version";
    absent.physicalCoefficients.clear();absent.physicalCoefficients[{f.ea,E_ObsCode::L1C,7}]=1;
    auto nuisance=raw;nuisance.targetId="nuisance-unbound";nuisance.nuisanceCoefficients["REAL_UNKNOWN_ANCHOR"]=1;
    compare(f,{raw,absent,nuisance,raw},"literal physical rows unknown versions and nuisance preserved");
    std::vector<ZhangFloatTargetResolution> answers;std::string reason;
    check(zhangResolveGraphProductHandoffTargetsBatch(f.state,{raw},answers,reason) && answers[0].resolved &&
        answers[0].row==std::map<KFKey,double>{{f.a,2},{f.b,-3}} && answers[0].offset==.25,
        "independent physical literal oracle failed");
    // Neither version7 individual arc exists. A retained complete sum can
    // represent the sum only; it must not invent either missing old endpoint.
    auto whole=f.target("retained-complete-sum");whole.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
    whole.physicalCoefficients[ {f.ea,E_ObsCode::L1C,7} ]=1;
    whole.physicalCoefficients[ {f.eb,E_ObsCode::L1C,7} ]=1;
    whole.numericRow={{f.a,1},{f.b,1}};whole.numericOffset=.5;whole.symbolicOffset=ZhangFloatCarrierRational(1)/2;
    f.seed({whole});auto sum=whole;sum.targetId="new-query-complete-sum";sum.physicalScale=2;sum.symbolicOffset=ZhangFloatCarrierRational(1)/4;
    compare(f,{sum,absent,sum},"whole exact span without individual RAW endpoints");
    check(zhangResolveGraphProductHandoffTargetsBatch(f.state,{sum,absent},answers,reason) && answers[0].resolved &&
        answers[0].row==std::map<KFKey,double>{{f.a,2},{f.b,2}} && answers[0].offset==.25 && !answers[1].resolved,
        "independent whole-span/unknown-direction oracle failed");
    {
        const auto metadata=f.state.metaDataMap;
        answers.resize(1);std::vector<ZhangProductFloatTargetRequest> directory(1);
        ZhangExactWorkBudget budget;budget.maximumCheckpoints=1;bool cancelled=false;
        try
        {
            ZhangScopedExactWorkBudget scope(budget);
            zhangResolveGraphProductHandoffTargetsBatch(f.state,{bound,sum},answers,reason,nullptr,&directory);
        }
        catch(const ZhangExactWorkBudgetExceeded&){cancelled=true;}
        check(cancelled && answers.empty() && directory.empty() && f.state.metaDataMap==metadata,
            "cancelled later query exposed a partial earlier batch result");
        ++cases;
    }
    const auto original=f.state.metaDataMap;
    auto corrupt=[&](const std::string& label,auto mutate) {mutate();compare(f,{bound,raw},label);f.state.metaDataMap=original;};
    corrupt("stale owner",[&]{auto t=bound;t.scope.runtimeId="FOREIGN";f.seed({t});});
    corrupt("stale exact order",[&]{f.state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]+="WRONG_ORDER";});
    corrupt("stale actual factor sequence",[&]{f.state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]="16";});
    corrupt("stale complete generation",[&]{f.state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]="8";});
    corrupt("bound absent KFKey",[&]{auto t=bound;KFKey missing;missing.type=KF::AMBIGUITY;missing.str="ABSENT";t.numericRow[missing]=1;f.seed({t});});
    corrupt("bound nonfinite numeric coefficient",[&]{auto t=bound;t.numericRow[f.a]=std::numeric_limits<double>::quiet_NaN();f.seed({t});});
    corrupt("wire target id disagreement",[&]{f.state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+"foreign-key"]=zhangEncodeFloatTarget(bound);});
    {
        Fixture typed;KFKey carrier;carrier.type=KF::PRODUCT_HANDOFF;carrier.str="ZHANG_FLOAT_CARRIER:RESOLVER_SUM";
        typed.state.kfIndexMap.emplace(carrier,4);typed.state.x.conservativeResize(5);typed.state.x(4)=7;
        typed.state.dx.conservativeResize(5);typed.state.dx(4)=.125;
        typed.state.P.conservativeResize(5,5);typed.state.P.row(4)<<0,5,10,5,15;
        typed.state.P.col(4)<<0,5,10,5,15;
        typed.state.stateTransitionMap[carrier][carrier][0]=1;
        auto target=typed.target("typed-whole-sum");target.numericRow={{carrier,1}};
        typed.seed({target});auto definition=target;definition.consumer="FLOAT_CARRIER_SOURCE_DEFINITION";
        definition.targetId="sum-definition";definition.numericRow={{typed.a,1},{typed.b,1}};
        typed.state.metaDataMap[zhangProductFloatCarrierRecordKey(carrier)]=zhangEncodeFloatTarget(definition);
        compare(typed,{target,target},"typed correlated carrier complete model");
        typed.state.procNoiseMap[carrier]=.01;
        compare(typed,{target},"typed carrier nonconstant process model rejected");
    }
    // A different coordinate chart requires a NEW context after real T.
    std::ostringstream trace;
    check(f.state.applyStateTransform(trace,{{KFState::oneKey,{{KFState::oneKey,1}}},
        {f.a,{{f.a,1}}},{f.b,{{f.b,1}}},{f.clock,{{f.clock,1},{f.a,1}}}},"RESOLVER_BATCH_REAL_T"),"real transform failed");
    compare(f,{bound,whole,sum,raw},"fresh batch AFTER actual coordinate transform");
    check(zhangResolveGraphProductHandoffTargetsBatch(f.state,{bound},answers,reason) && answers[0].resolved &&
        answers[0].row==std::map<KFKey,double>{{f.a,1.875},{f.b,-1},{f.clock,.125}} && answers[0].offset==.5,
        "independent actual T pullback failed");
    KFState restored=f.state;f.state=restored;
    compare(f,{bound,whole},"restored exact metadata snapshot new batch no serialized context");
    {
        const auto sourceMetadata=f.state.metaDataMap;
        const auto sourceIntegerArchive=integerArchive(f);
        Fixture copied;copied.state=f.state;
        copied.state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA]=
            zhangProductFloatOwnerId(f.state)+"/private-resolver-clone";
        // A disposable lease cannot export the authoritative integer archive.
        compare(copied,{bound},"unrebound private clone cannot reuse source-owner binding",false);
        check(cloneZhangGraphRuntime(f.state,copied.state,&reason),"public numeric graph-clone lease failed:"+reason);
        compare(copied,{bound,whole},"fresh private context after actual public clone owner rebind",false);
        check(f.state.metaDataMap==sourceMetadata && integerArchive(f)==sourceIntegerArchive,
            "private query/clone modified authoritative source metadata/integer archive");
        eraseZhangGraphRuntime(copied.state);
    }
    // The existing production fixture executes reserve -> real T
    // -> successful activation, rejection of lost reservation, later anchor
    // retention, collector and legacy RAW-preserving fallback. It now calls
    // BOTH newly wired batch consumers, rather than only the helper above.
    // Actual RAW/carrier retirement is covered by the existing independent
    // r51_product_float_carrier_production native fixture, not claimed here.
    std::string productionReport;
    check(zhangProductFloatConsumerProductionContractForAudit(productionReport),"actual consumer integration:"+productionReport);
    std::cout<<"PRODUCTION_CONSUMER_INTEGRATION "<<productionReport<<'\n';
}
int main()
{
    try
    {
        setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
        setenv("ZHANG_R51_PRODUCT_REPAIR","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
        setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);
        setenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT","0",1);
        acsConfig.zhangFullRank.enable=true;acsConfig.zhangPppAr.user_adapter=false;
        auto& sentinel=zhangProductRelationAdmissionStateRegistry()[{"RESOLVER_BATCH_GOLDEN_OWNER",E_Sys::GPS}];
        sentinel.certifiedRows={{{"CURRENT_INTEGER_ARCHIVE_SENTINEL",1}}};sentinel.certifiedValues={7};
        sentinel.certifiedSatellites={"G13"};sentinel.certifiedRelationIds={"ONLY_MUTATION_SENTINEL"};
        tests();std::cout<<"PASS immutable_resolver_batch cases="<<cases
            <<" old_single_resolver_is_golden=1 raw_results_equal=1 source_metadata_full_x_P_dx_unchanged=1"
            <<" ordinary_UNKNOWN_retained=1 fatal_context_atomic=1 current_integer_archive_unchanged=1 scientific_replay=0\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<"RESOLVER_BATCH_CONTRACT_FAILED "<<e.what()<<'\n';return 1;}
}
