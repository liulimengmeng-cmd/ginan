// Link with production PEA objects when testing is separately authorized.
// This source exercises the real update entry point, not a second invocation
// of the graph classifier as a numerical reference. Not executed here.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "pea/zhangReference.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}
static void addObservation(ReceiverMap& receivers, const ZhangGraphEdge& edge)
{
    auto& receiver=receivers[edge.receiver];
    receiver.id=edge.receiver; receiver.ready=true;
    auto observation=std::make_shared<GObs>();
    observation->Sat=edge.satellite;
    auto& status=receiver.satStatMap[edge.satellite];
    status.el=1.5; observation->satStat_ptr=&status;
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    {
        auto& signal=observation->sigs[code2Freq[E_Sys::GPS][code]];
        signal.code=code; signal.L=1e8; signal.P=2e7;
        auto& config=acsConfig.getRecOpts(edge.receiver,{"GPS",enum_to_string(code)});
        config.ambiguity.sigma={1000}; config.phase_bias.sigma={1000};
    }
    receiver.obsList.push_back(observation);
}
static void addPhase(KFState& state,const KFKey& key,double value)
{
    InitialState initial; initial.estimate=true;
    initial.x=value; initial.P=1; initial.Q=0;
    state.addKFState(key,initial);
}

int main()
{
    setenv("ZHANG_R51_ENABLE","1",1);
    setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
    setenv("ZHANG_R51_PRODUCT_REPAIR","1",1);
    acsConfig.zhangFullRank.product_basis_mode="SINGLE_S_BASIS_NATIVE";
    acsConfig.zhangFullRank.enable=true; acsConfig.process_sys[E_Sys::GPS]=true;
    auto& options=acsConfig.zhangFullRank.sysOpts[E_Sys::GPS];
    options.use_spanning_tree=true; options.reference_receiver="A";
    options.baseline_observables={E_ObsCode::L1C,E_ObsCode::L2W};
    options.product_core_min_satellite_support=0;
    ReceiverMap receivers;
    for (const auto& e : std::set<ZhangGraphEdge>{
        {"A",SatSys("G01")},{"A",SatSys("G02")},{"B",SatSys("G01")}})
        addObservation(receivers,e);
    KFState state;
    const std::string runtime="SINGLE_S_CHORD_CONTRACT";
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=runtime;
    std::ostringstream trace;
    updateZhangFullRankReferences(trace,receivers,state);
    ZhangGraphIntegerContext before;
    require(zhangGraphIntegerContext(state,E_Sys::GPS,before),"initial S chart missing");
    for (const auto code : options.baseline_observables)
    {
        for (const auto& receiver : before.basis.receivers)
            if (receiver!=before.basis.rootReceiver)
                addPhase(state,zhangTransportReceiverKey(E_Sys::GPS,code,receiver),20);
        for (const auto& satellite : before.basis.satellites)
            addPhase(state,zhangTransportSatelliteKey(code,satellite),10);
    }
    require(state.stateTransition(trace,state.time),"initial state creation failed");
    // Full dense correlated P, including all cross blocks, is the independent
    // oracle: a catalogue-only operation must preserve it BIT FOR BIT.
    VectorXd v=VectorXd::LinSpaced(state.x.size(),0.01,0.2);
    state.P=MatrixXd::Identity(state.x.size(),state.x.size())+v*v.transpose();
    const auto oldKeys=state.kfIndexMap;
    const VectorXd oldX=state.x; const MatrixXd oldP=state.P;
    const ZhangGraphEdge chord{"B",SatSys("G02")};
    addObservation(receivers,chord);
    state.time=state.time+30;
    updateZhangFullRankReferences(trace,receivers,state);
    ZhangGraphIntegerContext after;
    require(zhangGraphIntegerContext(state,E_Sys::GPS,after),"post-insertion S chart missing");
    require(after.basis.edges.contains(chord) && !after.basis.treeEdges.contains(chord),
        "production graph_update failed to commit the fresh chord");
    require(zhangSameCommittedGraphBasis(after.basis,after.productBasis),"independent P tree exists");
    require(after.representationVersion>before.representationVersion &&
        after.eventId==before.eventId &&
        after.productDatumVersion==before.productDatumVersion &&
        after.committedSGaugeVersion==before.committedSGaugeVersion,
        "pure chord insertion failed to version catalogue or spuriously changed S datum");
    require(after.committedSBasisIdentity!=before.committedSBasisIdentity,
        "complete snapshot identity ignored the fresh coordinate");
    require(state.kfIndexMap==oldKeys && state.x.size()==oldX.size() &&
        (state.x-oldX).norm()==0 && (state.P-oldP).norm()==0,
        "catalogue insertion reset/reordered old means, full P or cross blocks");
    for (const auto code : options.baseline_observables)
    {
        require(zhangGraphModelsObservation(state,"B",SatSys("G02"),code),
            "fresh signal admission missing before measurement construction");
        require(zhangGraphRetainsAmbiguity(state,"B",SatSys("G02"),code),
            "fresh cycle coordinate denied before measurement construction");
    }
    // The graph update above is deliberately pre-measurement. Complete the
    // same-epoch fresh-column creation which addDsgnEntry/addKFState performs
    // in the production measurement builder BEFORE taking a checkpoint.
    // Export below checks the graph section, not an entire KF checkpoint.
    for (const auto code : options.baseline_observables)
    {
        const auto key=zhangTransportAmbiguityKey(code,chord);
        require(!state.kfIndexMap.contains(key),"fresh chord already had a numerical column");
        addPhase(state,key,7);
    }
    require(state.stateTransition(trace,state.time),"same-epoch fresh-column creation failed");
    for (const auto code : options.baseline_observables)
        require(state.kfIndexMap.contains(zhangTransportAmbiguityKey(code,chord)),
            "same-epoch measurement preparation did not create fresh chord column");
    require(state.kfIndexMap.size()==oldKeys.size()+options.baseline_observables.size(),
        "fresh-column creation changed unrelated state dimensions");
    for (const auto& [key,oldIndex] : oldKeys)
    {
        const int newIndex=state.kfIndexMap.at(key);
        require(state.x[newIndex]==oldX[oldIndex],"fresh-column creation changed an old mean");
        for (const auto& [other,oldOtherIndex] : oldKeys)
            require(state.P(newIndex,state.kfIndexMap.at(other))==oldP(oldIndex,oldOtherIndex),
                "fresh-column creation lost an old covariance/cross block");
    }
    std::string bytes,reason;
    require(exportZhangGraphCheckpointSection(state,runtime,bytes,reason),reason.c_str());
    require(bytes.starts_with("GINAN_ZHANG_GRAPH_V5\n"),"mode-tagged checkpoint not exported");
    require(validateZhangGraphCheckpointSection(runtime,bytes,reason),reason.c_str());
    const auto fakeLegacy=std::string("GINAN_ZHANG_GRAPH_V4\n")+
        bytes.substr(std::string("GINAN_ZHANG_GRAPH_V5\n").size());
    require(!validateZhangGraphCheckpointSection(runtime,fakeLegacy,reason) &&
        reason.find("REQUIRES_EXPLICIT_MIGRATION")!=std::string::npos,
        "single-S loader accepted an archive without the mode contract");
    acsConfig.zhangFullRank.product_basis_mode="LEGACY_PRODUCT_TREE";
    ZhangGraphIntegerContext invalid;
    require(!zhangGraphIntegerContext(state,E_Sys::GPS,invalid),"live reload relabelled old S runtime");
    bool rejected=false;
    try { updateZhangFullRankReferences(trace,receivers,state); }
    catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"live mode change did not require restart/migration");
    acsConfig.zhangFullRank.product_basis_mode="SINGLE_S_BASIS_NATIVE";
    eraseZhangGraphRuntime(state);
    std::cout << "SINGLE_S_CHORD_PRODUCTION_CONTRACT_PASS\n";
}
