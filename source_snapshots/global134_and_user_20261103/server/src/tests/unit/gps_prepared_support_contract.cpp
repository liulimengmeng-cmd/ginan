// Link the actual ppp_ambres/Reference/algebra production objects. The oracle
// is the literal four-edge cycle and full covariance below, not another call
// to the new preparation/decoding helpers.
#include "common/acsConfig.hpp"
#include "common/zhangPreparedWholeFloatSupport.hpp"
#include "common/zhangPendingCurrentBridge.hpp"
#include "common/zhangR51PhysicalImage.hpp"
#include "common/zhangR51Integrity.hpp"
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool ok,const char* why)
{if(!ok)throw std::runtime_error(why);}

struct Fixture
{
    KFState state;
    std::string runtime,posterior="LITERAL_CURRENT_POSTERIOR",chart="LITERAL_CURRENT_CHART";
    const std::array<E_ObsCode,2> codes={E_ObsCode::L1C,E_ObsCode::L2W};
    const std::vector<std::string> ids={"L1C|X|G03|V0","L1C|A|G01|V0","L1C|B|G02|V0",
        "L1C|A|G02|V0","L1C|B|G01|V0"};
    std::array<KFKey,4> keys;
    std::shared_ptr<const ZhangCurrentProductFloatOwnerToken> owner;
    std::shared_ptr<const ZhangCurrentEpochCompiledDomain> domain;
    int calls=0;
    bool wholeOnly=false;
    int anchorValue=0;
    Fixture(std::string name,bool whole=false,int anchor=0)
        :runtime(std::move(name)),wholeOnly(whole),anchorValue(anchor)
    {
        acsConfig.zhangFullRank.sysOpts.clear();
        state.time.bigTime=1405218360;state.factorCommitSequence=7;
        state.lastPppTransactionCommitted=true;state.lastPppQcConverged=true;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=runtime;
        state.kfIndexMap[KFState::oneKey]=0;
        for(int i=0;i<4;++i) {keys[i].type=KF::AMBIGUITY;keys[i].Sat=SatSys("G01");
            keys[i].num=i+1;keys[i].str="INDEPENDENT_LITERAL_COORDINATE";state.kfIndexMap[keys[i]]=i+1;}
        state.x.resize(5);state.x<<1,2,5,7,11;state.dx=VectorXd::Zero(5);
        state.P=MatrixXd::Identity(5,5);state.P(0,0)=0;state.P(1,2)=state.P(2,1)=0.125;
        for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
        const auto order=zhangFloatCarrierStateOrder(state.kfIndexMap);
        for(int i=0;i<(wholeOnly?0:200);++i) {
            ZhangProductFloatTargetRequest t;t.targetId="LITERAL_TARGET_"+std::to_string(i);
            t.consumer="PENDING_PUBLICATION_DEPENDENCY";t.lifecycleId="LITERAL_LIFECYCLE";
            t.unit="m";t.symbolicSpace=ZhangProductFloatTargetSpace::SOURCE_KF_COORDINATE;
            t.scope={runtime,chart,order,7};t.numericRow={{keys[0],1}};
            t.numericBindingComplete=true;
            state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+t.targetId]=zhangEncodeFloatTarget(t);
        }
        if(wholeOnly) {
            ZhangProductFloatTargetRequest t;t.targetId="LITERAL_COMPLETE_OLD_CYCLE";
            t.consumer="PENDING_PUBLICATION_DEPENDENCY";t.lifecycleId="LITERAL_COMPLETE_LIFECYCLE";
            t.unit="cycles";t.symbolicSpace=ZhangProductFloatTargetSpace::INTEGER_CYCLE_WITH_EXPLICIT_NUISANCE;
            t.scope={runtime,chart,order,7};t.numericBindingComplete=true;
            const std::array<ZhangGraphEdge,4> edges={ZhangGraphEdge{"A",SatSys("G01")},
                ZhangGraphEdge{"B",SatSys("G02")},ZhangGraphEdge{"A",SatSys("G02")},
                ZhangGraphEdge{"B",SatSys("G01")}};
            for(int i=0;i<4;++i) {
                t.physicalCoefficients[{edges[i],codes[0],0}]=i<2?1:-1;
                t.numericRow[keys[i]]=i<2?1:-1;
            }
            state.metaDataMap[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+t.targetId]=zhangEncodeFloatTarget(t);
        }
        state.metaDataMap[ZHANG_FLOAT_TARGET_SCOPE_V2]=order;
        state.metaDataMap[ZHANG_FLOAT_TARGET_GENERATION_V2]="7";
        state.metaDataMap[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]="7";
        auto proof=std::make_shared<ZhangIntegerDecisionProof>();
        proof->id=runtime+"|INDEPENDENT_LITERAL_ANCHOR";
        proof->observationProvenance=runtime+"|PPP_COMMITTED|E"+std::to_string(static_cast<long int>(state.time.bigTime));
        proof->conditionalFailureBound=1e-7;proof->originalStatement="LITERAL_ANCHOR_N_X_G03_V0";
        ZhangCurrentEpochEvidence e;e.valid=true;e.runtimeId=runtime;e.system=E_Sys::GPS;
        e.epochGpstSeconds=static_cast<long int>(state.time.bigTime);e.posteriorId=posterior;e.chartId=chart;
        e.floatStateRoot=zhangR51NumericRoot(state.x,state.P);e.floatPosteriorPsdStateOrderId=order;
        e.rows={{{ids.front(),1}}};e.values={anchorValue};e.rowProofs={{proof}};
        e.floatPosteriorPsdRootId=e.floatStateRoot;e.floatPosteriorPsdValidated=true;
        require(e.addNumericConstraint({{keys[0],1}},anchorValue,e.rows.front()),"literal current anchor numerical binding");
        auto& registry=zhangCurrentEpochEvidenceRegistry();registry[{runtime,E_Sys::GPS}]=std::move(e);
        auto& actual=registry.at({runtime,E_Sys::GPS});
        domain=actual.compiledCurrentDomain();require(domain && domain->feasible,"literal current integer domain infeasible");
        std::string why;
        owner=zhangBindCurrentProductFloatOwner(state,E_Sys::GPS,actual.epochGpstSeconds,posterior,chart,why);
        require(bool(owner),"actual immutable owner factory rejected fixture");
    }
    ~Fixture(){zhangCurrentEpochEvidenceRegistry().erase({runtime,E_Sys::GPS});}
    std::shared_ptr<ZhangPreparedWholeFloatSupport> context()
    {
        return zhangAcquirePreparedWholeFloatSupport(state,E_Sys::GPS,codes,owner,domain,ids,
            runtime,static_cast<long int>(state.time.bigTime),posterior,chart);
    }
    ZhangWholeFloatSupportResult query(ZhangPreparedWholeFloatSupport* p,
        const ZhangBridgeSpanMatrix& targets,const ZhangBridgeSpanRow& offsets)
    {
        auto arcs=[&](const std::string& id,ZhangCoordinateRow& row,std::string& why) {
            ++calls;
            if(wholeOnly){why="RETIRED_ARC_AUXILIARY_STATE_MISSING";return false;}
            for(int i=0;i<4;++i)if(id==ids[i+1]){row={{keys[i],1}};why="CURRENT_VERSIONED_ARC_ROW";return true;}
            why="RETIRED_ARC_AUXILIARY_STATE_MISSING";return false;
        };
        auto current=[](std::string&){return true;};
        return zhangResolveWholeFloatSupportProduction(state,E_Sys::GPS,codes,ids,
            {{1,0,0,0,0}},{anchorValue},targets,offsets,arcs,current,nullptr,"LITERAL_GPS_CYCLE",p);
    }
};
static void literal(const Fixture& f,const ZhangWholeFloatSupportResult& answer,double scale,int offset)
{
    require(answer.valid && answer.rows.size()==1 && answer.offsets[0]==offset,"literal target rejected");
    ZhangCoordinateRow expected={{f.keys[0],scale},{f.keys[1],scale},{f.keys[2],-scale},{f.keys[3],-scale}};
    require(answer.rows[0]==expected,"whole support differs from independent literal cycle");
    Eigen::RowVectorXd j=Eigen::RowVectorXd::Zero(5);
    for(const auto& [key,value]:answer.rows[0])j[f.state.kfIndexMap.at(key)]=value;
    require(std::abs((j*f.state.x)[0]+11*scale)<1e-13,"literal mean changed");
    require(std::abs((j*f.state.P*j.transpose())(0,0)-4.25*scale*scale)<1e-12,"complete cross-Q lost");
}
static ZhangPendingCurrentBridge pending(Fixture& f,bool affine=false)
{
    ZhangPendingCurrentBridge p;
    auto& base=zhangCurrentEpochEvidenceRegistry().at({f.runtime,E_Sys::GPS});
    p.transaction=std::make_shared<ZhangCurrentEpochEventTransaction>(base);
    p.exactAudited=true;p.system=E_Sys::GPS;p.codes=f.codes;p.satellite=SatSys("G03");
    p.physicalRows={{{f.ids[1],1},{f.ids[2],1},{f.ids[3],-1},{f.ids[4],-1}}};
    if(affine)p.physicalRows[0][f.ids[0]]=1;
    p.values={-11+(affine?f.anchorValue:0)};p.targetValues={p.values[0],p.values[0]};
    p.residualDimension=1;p.witness="LITERAL_PENDING_AUDITED_ORIGINAL_IMAGE";
    auto proof=std::make_shared<ZhangIntegerDecisionProof>();proof->id=f.runtime+"|LITERAL_CYCLE_DECISION";
    proof->observationProvenance=f.runtime+"|PPP_COMMITTED|E"+std::to_string(base.epochGpstSeconds);
    proof->conditionalFailureBound=1e-7;proof->originalStatement="INDEPENDENT_LITERAL_CYCLE_EQUALS_MINUS_ELEVEN";
    p.proofs={proof};return p;
}
int main()
{
    const ZhangBridgeSpanMatrix original={{0,1,1,-1,-1}},residual={{0,2,2,-2,-2}};
    Fixture f("GPS_PREPARED_SUPPORT_NATIVE");
    std::string why;require(zhangProductFloatTargetMetadataValid(f.state,why),"legacy validator fixture rejected");
    const auto before=f.state;
    auto context=f.context();
    literal(f,f.query(context.get(),original,{0}),1,0);
    const auto firstCalls=f.calls;
    literal(f,f.query(context.get(),residual,{3}),2,3);
    require(f.calls==firstCalls && context->queryCalls==2 && context->decodeEntries==200,
        "original and residual did not share one decoded/physical support context");
    require(context->registered.front().scope.stateOrderId.empty(),"shared order was copied into every target");
    require(f.context()==context,"same frozen owner/domain failed to reuse prepared context");
    auto cold=f.query(nullptr,original,{0});literal(f,cold,1,0);
    require(f.state.metaDataMap==before.metaDataMap && f.state.kfIndexMap==before.kfIndexMap &&
        !std::memcmp(f.state.x.data(),before.x.data(),5*sizeof(double)) &&
        !std::memcmp(f.state.P.data(),before.P.data(),25*sizeof(double)),"preparation mutated source posterior");
    f.state.metaDataMap["EXPLICIT_SAME_SEQUENCE_METADATA_MUTATION"]="1";
    require(!f.query(context.get(),original,{0}).valid,"same-sequence metadata mutation reused frozen support");
    f.state.metaDataMap.erase("EXPLICIT_SAME_SEQUENCE_METADATA_MUTATION");
    f.state.factorCommitSequence++;
    require(!f.query(context.get(),original,{0}).valid,"stale factor sequence reused frozen support");
    f.state.factorCommitSequence--;
    auto& evidence=zhangCurrentEpochEvidenceRegistry().at({f.runtime,E_Sys::GPS});
    evidence.invalidateCurrentDomain();
    require(!f.query(context.get(),original,{0}).valid,"stale current integer generation reused support");

    Fixture resumed("GPS_PREPARED_SUPPORT_CANCEL");
    auto paused=resumed.context();
    ZhangExactWorkBudget budget;budget.maximumCheckpoints=80;
    bool cancelled=false;
    try {ZhangScopedExactWorkBudget scope(budget);resumed.query(paused.get(),original,{0});}
    catch(const ZhangExactWorkBudgetExceeded&){cancelled=true;}
    require(cancelled && !paused->prepared,"controlled preparation cancellation was not observed");
    literal(resumed,resumed.query(paused.get(),original,{0}),1,0);
    require(paused->decodeEntries==200,"resumed metadata was decoded more than once");
    literal(resumed,resumed.query(paused.get(),residual,{3}),2,3);

    // End-to-end finalisation starts at an already exact-audited candidate.
    // No individual old arc is readable, yet its complete typed cycle is.
    // This exercises the real production finaliser, delta gates and commit;
    // it deliberately makes no claim of executing a GNSS ILS search here.
    for(bool affine:{false,true}) {
        Fixture old("GPS_PENDING_WHOLE_"+std::to_string(affine),true,affine?2:0);
        auto support=old.context();literal(old,old.query(support.get(),original,{0}),1,0);
        require(support->arcRows.empty() && support->arcResolved.size()==old.ids.size(),
            "missing individual old arcs unexpectedly acquired CURRENT rows");
        auto p=pending(old,affine);p.support=support;
        auto& base=zhangCurrentEpochEvidenceRegistry().at({old.runtime,E_Sys::GPS});
        std::ostringstream trace;std::string why;
        require(zhangFinishPendingCurrentBridge(trace,old.state,base,p,why),
            "whole-readable individual-missing lift did not finish/commit");
        auto expected=ZhangCoordinateRow{{old.keys[0],1},{old.keys[1],1},{old.keys[2],-1},{old.keys[3],-1}};
        if(affine)expected[KFState::oneKey]=2;
        require(p.numericAudited && p.numericRows.size()==1 && p.numericRows[0]==expected &&
            base.eventOverlapRows==p.physicalRows && base.eventOverlapValues==p.values &&
            base.validatePhysicalNumericBindings(why),"whole lifted numeric/affine/full-binding identity differs");
        const auto generation=base.constraintGeneration;
        require(zhangFinishPendingCurrentBridge(trace,old.state,base,p,why) &&
            base.constraintGeneration==generation,"committed finaliser was not idempotent");
        auto other=old.state;
        require(!zhangFinishPendingCurrentBridge(trace,other,base,p,why),"committed receipt accepted another KF owner");
        old.state.factorCommitSequence++;
        require(!zhangFinishPendingCurrentBridge(trace,old.state,base,p,why),"committed receipt accepted later factor sequence");
        old.state.factorCommitSequence--;
        old.state.x(1)+=1;
        require(!zhangFinishPendingCurrentBridge(trace,old.state,base,p,why),"committed receipt accepted changed full FLOAT root");
        old.state.x(1)-=1;
        old.state.metaDataMap["COMMITTED_SCOPE_FAULT"]="1";
        require(!zhangFinishPendingCurrentBridge(trace,old.state,base,p,why),"committed receipt accepted changed graph metadata");
        old.state.metaDataMap.erase("COMMITTED_SCOPE_FAULT");
    }
    for(int fault=0;fault<3;++fault) {
        Fixture old("GPS_PENDING_BAD_ONE_"+std::to_string(fault),true,2);
        auto support=old.context();require(old.query(support.get(),original,{0}).valid,"bad ONE setup whole query");
        auto p=pending(old,true);p.support=support;
        // Deliberate post-freeze fault injection: ordinary callers may not
        // mutate a borrowed owner, but the constant gate must fail closed.
        if(fault==0)old.state.x(0)=2;
        if(fault==1)old.state.P(0,1)=1e-300;
        if(fault==2)old.state.P(1,0)=1e-300;
        auto& base=zhangCurrentEpochEvidenceRegistry().at({old.runtime,E_Sys::GPS});
        std::ostringstream trace;
        require(!zhangFinishPendingCurrentBridge(trace,old.state,base,p,why) &&
            why=="CURRENT_EVENT_OVERLAP_PENDING_AFFINE_CONSTANT_NOT_EXACT" && base.eventOverlapRows.empty(),
            "nonunit/infinitesimal nonzero-Q ONE was accepted");
    }
    {
        Fixture old("GPS_PENDING_NONBINARY_AFFINE",true);
        auto& base=zhangCurrentEpochEvidenceRegistry().at({old.runtime,E_Sys::GPS});
        const ZhangCurrentEpochPhysicalRow relation={{old.ids[0],3},{old.ids[1],1},
            {old.ids[2],1},{old.ids[3],-1},{old.ids[4],-1}};
        base.rows={relation};base.values={1};base.numericConstraints.clear();base.physicalNumericBindings.clear();
        require(base.addNumericConstraint({{old.keys[0],4},{old.keys[1],1},{old.keys[2],-1},{old.keys[3],-1}},1,relation),
            "nonbinary fixture lacks current-H same-root numerical recipe");
        base.invalidateCurrentDomain();old.domain=base.compiledCurrentDomain();
        require(old.domain && old.domain->feasible && base.validatePhysicalNumericBindings(why),
            "3X+cycle=1 must be integer feasible with full numerical bindings");
        old.owner=zhangBindCurrentProductFloatOwner(old.state,E_Sys::GPS,base.epochGpstSeconds,
            old.posterior,old.chart,why);require(bool(old.owner),"nonbinary owner");
        std::vector<std::string> columns(old.domain->columns.size());
        for(const auto& [id,c]:old.domain->columns)columns[c]=id;
        auto support=zhangAcquirePreparedWholeFloatSupport(old.state,E_Sys::GPS,old.codes,old.owner,old.domain,
            columns,old.runtime,base.epochGpstSeconds,old.posterior,old.chart);
        ZhangBridgeSpanMatrix query(1,ZhangBridgeSpanRow(columns.size()));query[0][old.domain->columns.at(old.ids[0])]=1;
        const auto answer=zhangResolveWholeFloatSupportProduction(old.state,E_Sys::GPS,old.codes,columns,
            old.domain->basis,old.domain->basisValues,query,{0},
            [](const std::string&,ZhangCoordinateRow&,std::string& reason){reason="RETIRED_ARC_AUXILIARY_STATE_MISSING";return false;},
            [&](std::string& reason){return base.validatePhysicalNumericBindings(reason);},nullptr,"NONBINARY",support.get());
        require(answer.valid && answer.offsets[0]==ZhangBridgeSpanRational(1,3),
            "3X+cycle=1 whole oracle did not retain the exact one-third offset");
        auto p=pending(old);p.physicalRows={{{old.ids[0],1}}};p.values={0};p.targetValues={0,0};p.support=support;
        std::ostringstream trace;
        require(!zhangFinishPendingCurrentBridge(trace,old.state,base,p,why) &&
            why=="CURRENT_EVENT_OVERLAP_PENDING_AFFINE_CONSTANT_NOT_EXACT" && base.eventOverlapRows.empty(),
            "nonbinary current-H affine offset silently rounded into a binding");
        base.physicalNumericBindings.clear();
        auto missing=pending(old);missing.physicalRows=p.physicalRows;missing.values={0};missing.targetValues={0,0};missing.support=support;
        require(!zhangFinishPendingCurrentBridge(trace,old.state,base,missing,why) && base.eventOverlapRows.empty(),
            "missing current-H numerical binding was bypassed by integer feasibility");
    }
    // Sweep every proof-trace/reserve/commit cancellation boundary on the
    // actual finaliser. A cancelled private transaction must leave base empty;
    // a published registry must always correspond to a true return. Before
    // the fix, the post-commit proof closure threw at one of these limits.
    bool sawCancelled=false,sawPublished=false;std::set<std::string> cancelledPhases;int publishedCount=0;
    for(std::uint64_t limit=1;limit<=40;++limit) {
        Fixture f("GPS_PENDING_COMMIT_BOUNDARY_"+std::to_string(limit));
        auto p=pending(f);auto& base=zhangCurrentEpochEvidenceRegistry().at({f.runtime,E_Sys::GPS});
        p.numericRows={{{f.keys[0],1},{f.keys[1],1},{f.keys[2],-1},{f.keys[3],-1}}};
        p.numericAudited=true;
        require(p.transaction->stage(p.physicalRows[0],p.values[0],p.proofs,p.witness,p.numericRows[0],why) &&
            p.transaction->prepare(why),"private ready boundary fixture failed");
        p.prepared=true;p.stagedRows=1;
        const auto oldGeneration=base.constraintGeneration;
        ZhangExactWorkBudget b;b.maximumCheckpoints=limit;
        bool ok=false,threw=false;std::ostringstream trace;
        try {ZhangScopedExactWorkBudget scope(b);ok=zhangFinishPendingCurrentBridge(trace,f.state,base,p,why);}
        catch(const ZhangExactWorkBudgetExceeded&){threw=true;}
        const bool published=!base.eventOverlapRows.empty();
        require(published==ok && !(published && threw),"published decision was reported cancelled/rejected");
        if(!published)require(base.constraintGeneration==oldGeneration && base.numericConstraints.size()==1,
            "cancellation changed logical current evidence");
        sawCancelled|=threw;sawPublished|=published;
        if(threw)cancelledPhases.insert(b.phase);
        publishedCount+=published;
        if(threw)require(zhangFinishPendingCurrentBridge(trace,f.state,base,p,why),"private cancelled finaliser could not resume");
    }
    require(sawCancelled && sawPublished,"commit boundary sweep missed cancellation or publication");
    require((cancelledPhases.contains("DECISION_PROOF_CLOSURE_ROOT") ||
        cancelledPhases.contains("DECISION_PROOF_CLOSURE_VISIT")) &&
        cancelledPhases.contains("EVENT_DELTA_RESERVE_ROWS") &&
        cancelledPhases.contains("EVENT_DELTA_BEFORE_ATOMIC_PUBLICATION"),
        "boundary sweep missed private proof closure/reserve/last atomic gate");
    std::cout<<"commit_boundary_published="<<publishedCount<<" cancelled_phases=";
    for(const auto& phase:cancelledPhases)std::cout<<phase<<",";std::cout<<"\n";
    std::cout<<"PASS production frozen support original/residual shared decode=200 exact current H "
        "literal nuisance cancellation mean/fullQ unchanged stale metadata/sequence/generation rejected "
        "cancelled preparation resumed; pending whole lift/affine full binding committed; "
        "private cancellation never reported false after publication\n";
}

