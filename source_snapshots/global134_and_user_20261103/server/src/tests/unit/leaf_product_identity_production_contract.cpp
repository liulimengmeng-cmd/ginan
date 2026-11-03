// Links the real PEA update entry point. The numerical oracle is the literal
// dense pre-event posterior, not another call to the transport implementation.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "pea/zhangReference.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

static void leaf_require(bool condition, const std::string& reason)
{
    if (!condition) throw std::runtime_error(reason);
}

static KFKey leaf_clock(const SatSys& satellite)
{
    KFKey key; key.type=KF::SAT_CLOCK; key.Sat=satellite; return key;
}

static void leaf_add_state(KFState& state,const KFKey& key,double value)
{
    InitialState initial; initial.estimate=true;
    initial.x=value; initial.P=1; initial.Q=0;
    state.addKFState(key,initial);
}

static void leaf_deterministic_one(KFState& state)
{
    // ONE is an affine constant, not a stochastic physical coordinate.
    // addKFState injects a new initial mean through ONE; assigning it a
    // variance would deliberately correlate every fresh prior with old states.
    const int one=state.kfIndexMap.at(KFState::oneKey);
    state.x(one)=1;
    state.dx(one)=0;
    state.P.row(one).setZero();
    state.P.col(one).setZero();
    leaf_require(state.P.row(one).isZero() && state.P.col(one).isZero(),
        "ONE must have zero variance and cross covariance");
}

static void leaf_add_observation(ReceiverMap& receivers,const ZhangGraphEdge& edge)
{
    auto& receiver=receivers[edge.receiver];
    receiver.id=edge.receiver; receiver.ready=true;
    auto observation=std::make_shared<GObs>(); observation->Sat=edge.satellite;
    auto& status=receiver.satStatMap[edge.satellite]; status.el=1.5;
    observation->satStat_ptr=&status;
    const auto system=edge.satellite.sys;
    const auto& codes=acsConfig.zhangFullRank.sysOpts[system].baseline_observables;
    leaf_require(!codes.empty(),"fixture baseline observables not configured");
    for (const auto code : codes)
    {
        auto& signal=observation->sigs[code2Freq[system][code]];
        signal.code=code; signal.L=1e8; signal.P=2e7;
        auto& config=acsConfig.getRecOpts(edge.receiver,{enum_to_string(system),enum_to_string(code)});
        config.ambiguity.sigma={1000}; config.phase_bias.sigma={1000};
    }
    receiver.obsList.push_back(observation);
}

struct LeafFixture
{
    KFState state;
    ReceiverMap receivers;
    std::ostringstream trace;
    ZhangGraphIntegerContext initial;
    explicit LeafFixture(const std::string& runtime)
    {
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=runtime;
        for (const auto& edge : std::set<ZhangGraphEdge>{
            {"A",SatSys("G01")},{"A",SatSys("G02")},{"B",SatSys("G01")}})
            leaf_add_observation(receivers,edge);
        updateZhangFullRankReferences(trace,receivers,state);
        leaf_require(zhangGraphIntegerContext(state,E_Sys::GPS,initial),"initial S missing");
        for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
        {
            for (const auto& receiver : initial.basis.receivers)
                if (receiver!=initial.basis.rootReceiver)
                    leaf_add_state(state,zhangTransportReceiverKey(E_Sys::GPS,code,receiver),20);
            for (const auto& satellite : initial.basis.satellites)
                leaf_add_state(state,zhangTransportSatelliteKey(code,satellite),10+int(code)*0.01);
        }
        for (const auto& satellite : initial.basis.satellites)
            leaf_add_state(state,leaf_clock(satellite),100+satellite.prn);
        leaf_require(state.stateTransition(trace,state.time),"initial state creation failed");
        VectorXd v=VectorXd::LinSpaced(state.x.size(),0.01,0.2);
        state.P=MatrixXd::Identity(state.x.size(),state.x.size())+v*v.transpose();
        leaf_deterministic_one(state);
    }
    ~LeafFixture() { eraseZhangGraphRuntime(state); }
    ZhangGraphIntegerContext extend(const std::set<ZhangGraphEdge>& edges,bool nextEpoch=true)
    {
        for (const auto& edge : edges) leaf_add_observation(receivers,edge);
        if (nextEpoch) state.time=state.time+30;
        updateZhangFullRankReferences(trace,receivers,state);
        ZhangGraphIntegerContext result;
        leaf_require(zhangGraphIntegerContext(state,E_Sys::GPS,result),"post-event S missing");
        return result;
    }
};

static ZhangGraphProductEventReceipt leaf_receipt(LeafFixture& fixture,
    E_ObsCode code,const SatSys& satellite)
{
    ZhangGraphProductEventReceipt receipt;
    leaf_require(zhangGraphProductEventReceipt(fixture.state,E_Sys::GPS,code,satellite,receipt),
        "old product receipt missing");
    return receipt;
}

static void leaf_check_posterior_unchanged(const LeafFixture& fixture,
    const std::map<KFKey,int>& keys,const VectorXd& means,const MatrixXd& covariance)
{
    leaf_require(fixture.state.kfIndexMap==keys && fixture.state.x.size()==means.size() &&
        (fixture.state.x-means).norm()==0 && (fixture.state.P-covariance).norm()==0,
        "leaf catalogue operation changed old means/full Q/cross blocks");
}

static void leaf_positive_and_same_epoch()
{
    LeafFixture f("LEAF_IDENTITY_DENSE_POSTERIOR");
    const auto keys=f.state.kfIndexMap;
    const VectorXd means=f.state.x; const MatrixXd covariance=f.state.P;
    const auto after=f.extend({{"C",SatSys("G02")},{"C",SatSys("G03")}});
    leaf_require(after.eventId==f.initial.eventId+1 &&
        after.representationVersion>f.initial.representationVersion,
        "leaf extension did not retain monotonic backend/catalogue generations");
    leaf_require(after.basis.rootReceiver==f.initial.basis.rootReceiver &&
        std::includes(after.basis.treeEdges.begin(),after.basis.treeEdges.end(),
            f.initial.basis.treeEdges.begin(),f.initial.basis.treeEdges.end()),
        "positive fixture failed to retain old root/tree");
    leaf_check_posterior_unchanged(f,keys,means,covariance);
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    for (const auto& satellite : f.initial.basis.satellites)
    {
        const auto receipt=leaf_receipt(f,code,satellite);
        leaf_require(receipt.eventKind=="LEAF_EXTENSION_EXISTING_PRODUCT_IDENTITY" &&
            receipt.oldBackendGeneration==f.initial.eventId &&
            receipt.newBackendGeneration==after.eventId &&
            receipt.oldProductProjectionCaptured && receipt.representableSameTimeMapProven &&
            receipt.fullJointMomentTransportProven && receipt.actualWriterProjectionBound &&
            receipt.actualWriterRowMatches && receipt.sameTimeMapProven,
            "leaf receipt lacks actual-writer/full-joint evidence: "+receipt.failureReason);
        const int phase=keys.at(zhangTransportSatelliteKey(code,satellite));
        const int clock=keys.at(leaf_clock(satellite));
        leaf_require(receipt.oldPhaseMeanMetres==means[phase] &&
            receipt.oldClockMeanMetres==means[clock] &&
            receipt.oldPhaseVarianceMetres2==covariance(phase,phase) &&
            receipt.oldClockVarianceMetres2==covariance(clock,clock),
            "old actual product/clock moments differ from independent posterior");
        for (const auto otherCode : {E_ObsCode::L1C,E_ObsCode::L2W})
        for (const auto& other : f.initial.basis.satellites)
        {
            const int q=keys.at(zhangTransportSatelliteKey(otherCode,other));
            const int c=keys.at(leaf_clock(other));
            leaf_require(receipt.oldPhaseCrossCovarianceMetres2.at({other,otherCode})==covariance(phase,q) &&
                receipt.oldPhaseClockCrossCovarianceMetres2.at(other)==covariance(phase,c) &&
                receipt.oldClockClockCrossCovarianceMetres2.at(other)==covariance(clock,c),
                "old full joint phase/clock covariance differs from dense oracle");
        }
        ZhangGraphProductEventReceipt fresh;
        leaf_require(!zhangGraphProductEventReceipt(f.state,E_Sys::GPS,code,SatSys("G03"),fresh),
            "new satellite received fabricated old endpoint");
    }
    // A second catalogue-only leaf in the same epoch must retain the first
    // published source. It is not an arbitrary second coordinate transform.
    const auto second=f.extend({{"D",SatSys("G01")}},false);
    leaf_require(second.eventId==after.eventId+1,"second leaf reused backend generation");
    leaf_check_posterior_unchanged(f,keys,means,covariance);
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    {
        const auto receipt=leaf_receipt(f,code,SatSys("G01"));
        leaf_require(receipt.oldBackendGeneration==f.initial.eventId &&
            receipt.newBackendGeneration==second.eventId && receipt.sameTimeMapProven,
            "same-epoch leaf replaced earliest source or invalidated unchanged destination");
    }
    // Complete the measurement-builder side of this real KF operation:
    // fresh coordinates have a zero cross block, while every old mean and
    // covariance entry is embedded without copying a made-up old endpoint.
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    {
        for (const auto& receiver : {std::string("C"),std::string("D")})
            leaf_add_state(f.state,zhangTransportReceiverKey(E_Sys::GPS,code,receiver),30);
        leaf_add_state(f.state,zhangTransportSatelliteKey(code,SatSys("G03")),13);
    }
    leaf_add_state(f.state,leaf_clock(SatSys("G03")),103);
    leaf_require(f.state.stateTransition(f.trace,f.state.time),"leaf coordinate augmentation failed");
    for (const auto& [key,index] : keys)
    {
        const int now=f.state.kfIndexMap.at(key);
        leaf_require(f.state.x[now]==means[index],"augmentation changed old mean");
        for (const auto& [other,oldOther] : keys)
            leaf_require(f.state.P(now,f.state.kfIndexMap.at(other))==covariance(index,oldOther),
                "augmentation changed old covariance/cross block");
        for (const auto& [fresh,freshIndex] : f.state.kfIndexMap)
            if (!keys.contains(fresh))
                leaf_require(f.state.P(now,freshIndex)==0 && f.state.P(freshIndex,now)==0,
                    "fresh leaf cross block was not a zero embedding");
    }
    updateZhangFullRankReferences(f.trace,f.receivers,f.state);
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    {
        const auto receipt=leaf_receipt(f,code,SatSys("G01"));
        leaf_require(receipt.oldBackendGeneration==f.initial.eventId &&
            receipt.newBackendGeneration==second.eventId && receipt.sameTimeMapProven,
            "augmented actual writer lost old-source binding");
        ZhangGraphProductEventReceipt fresh;
        leaf_require(!zhangGraphProductEventReceipt(f.state,E_Sys::GPS,code,SatSys("G03"),fresh),
            "fresh column creation fabricated a published old endpoint");
    }
}

static void leaf_missing_clock_or_phase(bool missingClock)
{
    LeafFixture f(missingClock?"LEAF_MISSING_CLOCK":"LEAF_MISSING_PHASE");
    const KFKey missing=missingClock?leaf_clock(SatSys("G02")):
        zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G02"));
    f.state.removeState(missing);
    leaf_require(f.state.stateTransition(f.trace,f.state.time),"negative state removal failed");
    f.extend({{"C",SatSys("G01")}});
    const auto receipt=leaf_receipt(f,E_ObsCode::L1C,SatSys("G02"));
    leaf_require(!receipt.sameTimeMapProven && !receipt.fullJointMomentTransportProven,
        "missing old clock/phase row was promoted to identity proof");
}

static void leaf_missing_covariance()
{
    LeafFixture f("LEAF_NONFINITE_CROSS_Q");
    const int p=f.state.kfIndexMap.at(zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G01")));
    const int c=f.state.kfIndexMap.at(leaf_clock(SatSys("G02")));
    f.state.P(p,c)=f.state.P(c,p)=std::numeric_limits<double>::quiet_NaN();
    const MatrixXd prior=f.state.P;
    const VectorXd mean=f.state.x;
    const auto keys=f.state.kfIndexMap;
    try { f.extend({{"C",SatSys("G01")}}); }
    catch (const std::runtime_error& error)
    {
        // The production global sweep can reject this deliberately corrupted
        // covariance before the leaf path is reached. That is a stronger
        // fail-closed outcome, not a reason to bypass the sweep in this test.
        leaf_require(std::string(error.what())==
            "ZHANG_HANDOFF_FLOAT_PRE_GRAPH_SWEEP_FAILED:HANDOFF_FLOAT_COVARIANCE_NONFINITE",
            "nonfinite negative failed for an unrelated reason");
        leaf_require(f.state.kfIndexMap==keys && (f.state.x-mean).isZero() &&
            f.state.P.rows()==prior.rows() && f.state.P.cols()==prior.cols() &&
            std::memcmp(f.state.P.data(),prior.data(),prior.size()*sizeof(double))==0,
            "nonfinite pre-graph rejection changed the numerical posterior");
        ZhangGraphProductEventReceipt rejected;
        leaf_require(!zhangGraphProductEventReceipt(f.state,E_Sys::GPS,
            E_ObsCode::L1C,SatSys("G01"),rejected),
            "nonfinite pre-graph rejection fabricated a current receipt");
        return;
    }
    const auto receipt=leaf_receipt(f,E_ObsCode::L1C,SatSys("G01"));
    leaf_require(!receipt.sameTimeMapProven && !receipt.fullJointMomentTransportProven,
        "missing/nonfinite joint covariance was filled or skipped");
}

static void leaf_chord_and_real_tree_change()
{
    LeafFixture f("LEAF_CHORD_NOT_LEAF");
    const auto keys=f.state.kfIndexMap;
    const VectorXd means=f.state.x; const MatrixXd covariance=f.state.P;
    const auto after=f.extend({{"B",SatSys("G02")}});
    leaf_require(after.eventId==f.initial.eventId,"pure chord created a backend event");
    leaf_check_posterior_unchanged(f,keys,means,covariance);
    leaf_require(f.trace.str().find("ZHANG_PRODUCT_LEAF_IDENTITY_RECEIPT")==std::string::npos,
        "pure chord was misclassified as leaf identity");
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
        leaf_add_state(f.state,zhangTransportAmbiguityKey(code,{"B",SatSys("G02")}),7);
    leaf_require(f.state.stateTransition(f.trace,f.state.time),"chord state creation failed");
    auto& observations=f.receivers.at("A").obsList;
    observations.erase(std::remove_if(observations.begin(),observations.end(),[](const auto& observation)
        { const auto obs=std::dynamic_pointer_cast<GObs>(observation);
          return obs && obs->Sat==SatSys("G02"); }),observations.end());
    f.state.time=f.state.time+30;
    updateZhangFullRankReferences(f.trace,f.receivers,f.state);
    leaf_require(f.trace.str().find("ZHANG_PRODUCT_LEAF_IDENTITY_RECEIPT")==std::string::npos,
        "real old-tree loss was authorized by leaf identity branch");
    ZhangGraphProductEventReceipt receipt;
    if (zhangGraphProductEventReceipt(f.state,E_Sys::GPS,E_ObsCode::L1C,SatSys("G02"),receipt))
        leaf_require(receipt.eventKind!="LEAF_EXTENSION_EXISTING_PRODUCT_IDENTITY",
            "real tree exchange reused identity leaf receipt");
}

static void leaf_after_real_transform_same_epoch()
{
    LeafFixture f("LEAF_AFTER_REAL_TREE_TRANSFORM");
    f.extend({{"B",SatSys("G02")}});
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
        leaf_add_state(f.state,zhangTransportAmbiguityKey(code,{"B",SatSys("G02")}),7);
    leaf_require(f.state.stateTransition(f.trace,f.state.time),"real-transform chord setup failed");
    VectorXd v=VectorXd::LinSpaced(f.state.x.size(),0.02,0.3);
    f.state.P=MatrixXd::Identity(f.state.x.size(),f.state.x.size())+v*v.transpose();
    leaf_deterministic_one(f.state);
    const auto sourceKeys=f.state.kfIndexMap;
    const VectorXd sourceMeans=f.state.x;
    const MatrixXd sourceCovariance=f.state.P;
    ZhangGraphIntegerContext source;
    leaf_require(zhangGraphIntegerContext(f.state,E_Sys::GPS,source),"source chart missing");
    auto& observations=f.receivers.at("A").obsList;
    observations.erase(std::remove_if(observations.begin(),observations.end(),[](const auto& observation)
        { const auto obs=std::dynamic_pointer_cast<GObs>(observation);
          return obs && obs->Sat==SatSys("G02"); }),observations.end());
    f.state.time=f.state.time+30;
    updateZhangFullRankReferences(f.trace,f.receivers,f.state);
    ZhangGraphIntegerContext transformed;
    leaf_require(zhangGraphIntegerContext(f.state,E_Sys::GPS,transformed) &&
        transformed.basis.treeEdges!=source.basis.treeEdges &&
        transformed.basis.rootReceiver==source.basis.rootReceiver &&
        transformed.eventId>source.eventId,
        "fixture did not execute a real production tree transform");
    std::map<std::pair<SatSys,E_ObsCode>,ZhangGraphProductEventReceipt> original;
    for (const auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
    for (const auto& satellite : source.basis.satellites)
    {
        const auto receipt=leaf_receipt(f,code,satellite);
        leaf_require(receipt.eventKind=="TREE_REBASE" || receipt.eventKind=="TREE_ARC_RETIREMENT" ||
            receipt.eventKind=="COMPONENT_REPARAMETERISATION",
            "fixture did not capture a real tree-event product source");
        leaf_require(receipt.oldProductProjectionCaptured &&
            receipt.representableSameTimeMapProven && receipt.fullJointMomentTransportProven &&
            receipt.actualWriterProjectionBound,
            "real tree receipt lacks complete source/destination evidence: "+receipt.failureReason);
        const int p=sourceKeys.at(zhangTransportSatelliteKey(code,satellite));
        const int c=sourceKeys.at(leaf_clock(satellite));
        leaf_require(receipt.oldPhaseMeanMetres==sourceMeans[p] &&
            receipt.oldClockMeanMetres==sourceMeans[c] &&
            receipt.oldPhaseVarianceMetres2==sourceCovariance(p,p) &&
            receipt.oldClockVarianceMetres2==sourceCovariance(c,c),
            "real transform captured a later source instead of independent pre-event posterior");
        for (const auto otherCode : {E_ObsCode::L1C,E_ObsCode::L2W})
        for (const auto& other : source.basis.satellites)
        {
            const int q=sourceKeys.at(zhangTransportSatelliteKey(otherCode,other));
            const int otherClock=sourceKeys.at(leaf_clock(other));
            leaf_require(receipt.oldPhaseCrossCovarianceMetres2.at({other,otherCode})==sourceCovariance(p,q) &&
                receipt.oldPhaseClockCrossCovarianceMetres2.at(other)==sourceCovariance(p,otherClock) &&
                receipt.oldClockClockCrossCovarianceMetres2.at(other)==sourceCovariance(c,otherClock),
                "real transform old joint cross-Q differs from literal pre-event oracle");
        }
        original[{satellite,code}]=receipt;
    }
    const auto destinationKeys=f.state.kfIndexMap;
    const VectorXd destinationMeans=f.state.x;
    const MatrixXd destinationCovariance=f.state.P;
    const auto extended=f.extend({{"C",SatSys("G01")}},false);
    leaf_require(extended.eventId==transformed.eventId+1,
        "same-epoch post-transform leaf failed to advance event generation");
    leaf_check_posterior_unchanged(f,destinationKeys,destinationMeans,destinationCovariance);
    for (const auto& [key,before] : original)
    {
        const auto after=leaf_receipt(f,key.second,key.first);
        leaf_require(after.eventKind==before.eventKind &&
            after.oldBackendGeneration==before.oldBackendGeneration &&
            after.newBackendGeneration==extended.eventId &&
            after.oldPhaseRow==before.oldPhaseRow && after.newPhaseRow==before.newPhaseRow &&
            after.actualWriterPhaseRow==before.actualWriterPhaseRow &&
            after.oldPhaseMeanMetres==before.oldPhaseMeanMetres &&
            after.oldClockMeanMetres==before.oldClockMeanMetres &&
            after.oldPhaseCrossCovarianceMetres2==before.oldPhaseCrossCovarianceMetres2 &&
            after.oldPhaseClockCrossCovarianceMetres2==before.oldPhaseClockCrossCovarianceMetres2 &&
            after.oldClockClockCrossCovarianceMetres2==before.oldClockClockCrossCovarianceMetres2 &&
            after.fullJointMomentTransportProven==before.fullJointMomentTransportProven &&
            after.sameTimeMapProven==before.sameTimeMapProven &&
            after.conditionalProductMapProven==before.conditionalProductMapProven,
            "same-epoch leaf overwrote source, changed rows/Q or promoted an unproven writer map");
    }
    leaf_require(f.trace.str().find("status=EXISTING_SOURCE_PRESERVED")!=std::string::npos,
        "real-transform-to-leaf path did not preserve the existing source receipt");
}

static void leaf_bds52_native_path(const std::vector<E_ObsCode>& codes,
    const std::string& caseName)
{
    // Same constellation and new-satellite category as the real C52 event;
    // correlated means/Q are intentionally synthetic and independently known.
    const auto system=E_Sys::BDS;
    acsConfig.process_sys[system]=true;
    auto& options=acsConfig.zhangFullRank.sysOpts[system];
    options.use_spanning_tree=true; options.reference_receiver="A";
    options.baseline_observables=codes;
    options.product_core_min_satellite_support=0;
    options.state_edge_grace_epochs=0;
    KFState state; ReceiverMap receivers; std::ostringstream trace;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="LEAF_BDS52_NATIVE_"+caseName;
    for (const auto& edge : std::set<ZhangGraphEdge>{
        {"A",SatSys("C19")},{"A",SatSys("C36")},{"B",SatSys("C19")}})
        leaf_add_observation(receivers,edge);
    updateZhangFullRankReferences(trace,receivers,state);
    ZhangGraphIntegerContext before;
    leaf_require(zhangGraphIntegerContext(state,system,before),"BDS initial chart missing");
    for (const auto code : codes)
    {
        leaf_add_state(state,zhangTransportReceiverKey(system,code,"B"),20);
        for (const auto& satellite : before.basis.satellites)
            leaf_add_state(state,zhangTransportSatelliteKey(code,satellite),satellite.prn);
    }
    for (const auto& satellite : before.basis.satellites)
        leaf_add_state(state,leaf_clock(satellite),100+satellite.prn);
    leaf_require(state.stateTransition(trace,state.time),"BDS state creation failed");
    VectorXd v=VectorXd::LinSpaced(state.x.size(),0.03,0.25);
    state.P=MatrixXd::Identity(state.x.size(),state.x.size())+v*v.transpose();
    leaf_deterministic_one(state);
    const auto keys=state.kfIndexMap;
    const VectorXd means=state.x; const MatrixXd covariance=state.P;
    leaf_add_observation(receivers,{"C",SatSys("C36")});
    leaf_add_observation(receivers,{"C",SatSys("C52")});
    state.time=state.time+30;
    updateZhangFullRankReferences(trace,receivers,state);
    ZhangGraphIntegerContext after;
    leaf_require(zhangGraphIntegerContext(state,system,after) &&
        after.eventId==before.eventId+1,"BDS leaf event missing or generation suppressed");
    leaf_require(state.kfIndexMap==keys && (state.x-means).norm()==0 &&
        (state.P-covariance).norm()==0,"BDS leaf changed old posterior");
    for (const auto code : codes)
    for (const auto& satellite : before.basis.satellites)
    {
        ZhangGraphProductEventReceipt receipt;
        leaf_require(zhangGraphProductEventReceipt(state,system,code,satellite,receipt) &&
            receipt.eventKind=="LEAF_EXTENSION_EXISTING_PRODUCT_IDENTITY" &&
            receipt.sameTimeMapProven && receipt.fullJointMomentTransportProven &&
            receipt.actualWriterRowMatches && receipt.oldBackendGeneration==before.eventId &&
            receipt.newBackendGeneration==after.eventId,
            "BDS C52 leaf stranded an existing product backend");
        const int p=keys.at(zhangTransportSatelliteKey(code,satellite));
        const int c=keys.at(leaf_clock(satellite));
        // This is the independent oracle: literal old dense x/P indexed by
        // the known native product and clock state keys. No second call to
        // the transport or product-projection implementation is made here.
        leaf_require(receipt.oldPhaseMeanMetres==means[p] &&
            receipt.oldClockMeanMetres==means[c] &&
            receipt.oldPhaseVarianceMetres2==covariance(p,p) &&
            receipt.oldClockVarianceMetres2==covariance(c,c) &&
            receipt.oldClockPhaseCovarianceMetres2==covariance(c,p),
            "BDS old product/clock joint mean/diagonal mismatch");
        leaf_require(receipt.oldPhaseRow==ZhangCoordinateRow{
            {zhangTransportSatelliteKey(code,satellite),1}} &&
            receipt.actualWriterPhaseRow==receipt.oldPhaseRow,
            "BDS native writer row is not the explicitly specified phase row");
        for (const auto& other : before.basis.satellites)
        {
            const int otherClock=keys.at(leaf_clock(other));
            leaf_require(receipt.oldPhaseClockCrossCovarianceMetres2.at(other)==
                covariance(p,otherClock) &&
                receipt.oldClockClockCrossCovarianceMetres2.at(other)==covariance(c,otherClock),
                "BDS old phase-clock/clock-clock cross Q mismatch");
            for (const auto otherCode : codes)
                leaf_require(receipt.oldPhaseCrossCovarianceMetres2.at({other,otherCode})==
                    covariance(p,keys.at(zhangTransportSatelliteKey(otherCode,other))),
                    "BDS old phase-phase cross Q mismatch");
        }
        ZhangGraphProductEventReceipt fresh;
        leaf_require(!zhangGraphProductEventReceipt(state,system,code,SatSys("C52"),fresh),
            "C52 incorrectly inherited a nonexistent old product endpoint");
    }
    eraseZhangGraphRuntime(state);
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
    options.state_edge_grace_epochs=0;
    leaf_positive_and_same_epoch();
    leaf_missing_clock_or_phase(true);
    leaf_missing_clock_or_phase(false);
    leaf_missing_covariance();
    leaf_chord_and_real_tree_change();
    leaf_after_real_transform_same_epoch();
    leaf_bds52_native_path({E_ObsCode::L2I,E_ObsCode::L6I},"L2I_L6I");
    leaf_bds52_native_path({E_ObsCode::L1P,E_ObsCode::L5P},"L1P_L5P_PRODUCTION");
    std::cout << "LEAF_PRODUCT_IDENTITY_PRODUCTION_CONTRACT_PASS\n";
}
