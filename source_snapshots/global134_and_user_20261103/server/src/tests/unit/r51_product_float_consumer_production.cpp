#include "pea/zhangPppAr.hpp"
#include "pea/zhangReference.hpp"
#include "common/acsConfig.hpp"
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool value,const char* why)
{if(!value)throw std::runtime_error(why);}

static void currentReceiptRowsSurviveRepeatedCollection()
{
    KFState state;state.time.bigTime=1405216800;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="PRODUCT_FLOAT_RECEIPT_LIFETIME_AUDIT";
    const ZhangGraphEdge ea{"A",SatSys("G01")},eb{"B",SatSys("G02")};
    const auto a=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,ea,0);
    const auto b=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,eb,0);
    state.kfIndexMap={{KFState::oneKey,0},{a,1},{b,2}};
    state.x.resize(3);state.x<<1,2,5;state.dx=VectorXd::Zero(3);
    state.P=MatrixXd::Zero(3,3);state.P(1,1)=4;state.P(2,2)=9;
    for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
    require(zhangHandoffFloatRegisterRecord(state,a,ea,int(E_ObsCode::L1C),0) &&
        zhangHandoffFloatRegisterRecord(state,b,eb,int(E_ObsCode::L1C),0),"receipt RAW fixture registration failed");
    std::string reason;std::ostringstream trace;
    ZhangProductFloatTargetRequest literal;
    literal.scope={zhangProductFloatOwnerId(state),"POST_GRAPH_RECEIPT",
        zhangFloatCarrierStateOrder(state.kfIndexMap),state.factorCommitSequence};
    literal.lifecycleId=state.time.to_string(9);literal.numericBindingComplete=true;
    literal.numericRow={{a,1},{b,1}};literal.unit="m";
    for(const auto role:{"GRAPH_RECEIPT_OLD_MAPPED","GRAPH_RECEIPT_WRITER","GRAPH_RECEIPT_BRIDGE"})
    {
        auto target=literal;target.consumer=role;
        target.targetId=std::string(role)+"|"+target.lifecycleId+"|G01|L1C";
        require(zhangRegisterGraphProductHandoffTarget(state,target,reason),"same-epoch receipt registration failed");
    }
    auto collect=[&]() {
        std::vector<ZhangProductFloatTargetRequest> result;
        require(zhangCollectActiveProductFloatConsumerTargets(state,result,reason),"real collector rejected owned current receipt");
        require(result.size()==3,"owned current receipt was silently omitted");
        for(const auto& target:result)require(target.numericBindingComplete && target.numericRow==literal.numericRow &&
            target.numericOffset==0 && target.scope==literal.scope,"repeat collector changed receipt numerical meaning");
    };
    collect();collect();
    zhangRegisterPublishedProductHandoffDependencies(state);
    collect();
    require(zhangHandoffFloatDependencies(state).size()==2,"actual dependency refresh lost receipt RAW pins");
    const VectorXd originalX=state.x;const MatrixXd originalP=state.P;
    state.time=state.time+30;
    std::vector<ZhangProductFloatTargetRequest> expired;
    require(zhangCollectActiveProductFloatConsumerTargets(state,expired,reason) && expired.empty(),
        "canonical strictly old graph receipt did not expire");
    zhangRegisterPublishedProductHandoffDependencies(state);
    require(zhangHandoffFloatDependencies(state).empty(),"expired receipt retained obsolete RAW dependency");
    std::vector<ZhangProductFloatTargetRequest> stored;
    require(zhangReadFloatTargets(state.metaDataMap,stored,reason) && stored.size()==3 &&
        std::all_of(stored.begin(),stored.end(),[](const auto& t){return !t.active;}),
        "next epoch registration reactivated expired graph receipts");
    auto unknown=literal;unknown.consumer="GRAPH_RECEIPT_FOO";unknown.lifecycleId=state.time.to_string(9);
    unknown.targetId="GRAPH_RECEIPT_FOO|"+unknown.lifecycleId+"|G01|L1C";
    require(zhangRegisterGraphProductHandoffTarget(state,unknown,reason),"unknown-role fixture registration failed");
    std::vector<ZhangProductFloatTargetRequest> unbound;
    require(!zhangCollectActiveProductFloatConsumerTargets(state,unbound,reason) &&
        std::any_of(unbound.begin(),unbound.end(),[](const auto& t){return t.consumer=="GRAPH_RECEIPT_FOO" && !t.numericBindingComplete;}),
        "GRAPH_RECEIPT prefix granted an unknown role complete-consumer status");
    const auto before=state;
    require(zhangCompressGraphProductFloatTargets(trace,state,reason) && reason.starts_with("CONSUMER_DIRECTORY_INCOMPLETE:"),
        "unknown graph receipt role bypassed RAW-preserving compression fallback");
    require(state.metaDataMap==before.metaDataMap && state.kfIndexMap==before.kfIndexMap &&
        state.factorCommitSequence==before.factorCommitSequence &&
        !std::memcmp(state.x.data(),originalX.data(),originalX.size()*sizeof(double)) &&
        !std::memcmp(state.P.data(),originalP.data(),originalP.size()*sizeof(double)),
        "receipt lifecycle refresh or unknown role fallback changed FLOAT moments");
}

int main()
{
    setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
    setenv("ZHANG_R51_PRODUCT_REPAIR","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
    setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);
    acsConfig.zhangFullRank.enable=true;acsConfig.zhangPppAr.user_adapter=false;
    std::string report;
    if (!zhangProductFloatConsumerProductionContractForAudit(report)) throw std::runtime_error(report);
    currentReceiptRowsSurviveRepeatedCollection();
    std::cout << "PRODUCT_FLOAT_CONSUMER_PRODUCTION_CONTRACT " << report << '\n';
    std::cout << "PASS owned_current_graph_receipt_rows_preserved=1 repeated_collector=1 real_dependency_registration=1 "
        "canonical_old_receipt_expired=1 unknown_GRAPH_RECEIPT_FOO_unbound=1 RAW_preserved=1\n";
}
