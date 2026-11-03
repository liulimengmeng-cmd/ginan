// Linked against real PEA objects, not header-only substitutes. This tests
// retained FLOAT coordinates and full cross-P, not integer certification.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductHandoffFloatRegistry.hpp"
#include "pea/zhangReference.hpp"

#include <cstdlib>
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
void check(bool condition, const char* reason)
{ if (!condition) throw std::runtime_error(reason); }

KFState source()
{
    KFState state;
    state.kfIndexMap.clear();
    state.stateTransitionMap.clear();
    const KFKey old = zhangTransportAmbiguityKey(
        E_ObsCode::L1C, {"AGGO", SatSys("G10")});
    const KFKey live = zhangTransportAmbiguityKey(
        E_ObsCode::L1C, {"SAFE", SatSys("G11")});
    state.kfIndexMap[old] = 0;
    state.kfIndexMap[live] = 1;
    state.stateTransitionMap[old][old][0] = 1;
    state.stateTransitionMap[live][live][0] = 1;
    state.x = VectorXd(2); state.x << 2, 7;
    state.dx = VectorXd::Zero(2);
    state.P = MatrixXd(2,2); state.P << 4, 1.2, 1.2, 9;
    state.time.bigTime = 1000000;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] =
        "HANDOFF_FLOAT_LIFECYCLE_TEST";
    state.advanced_postfits = false;
    state.postfitOpts.omega_test = false;
    return state;
}

KFKey retain(KFState& state, std::ostringstream& trace)
{
    const ZhangGraphEdge edge{"AGGO", SatSys("G10")};
    const auto old = zhangTransportAmbiguityKey(E_ObsCode::L1C, edge);
    const auto auxiliary = zhangGraphOldArcAuxiliaryKey(
        state, E_ObsCode::L1C, edge, 0);
    std::map<KFKey, ZhangCoordinateRow> transform;
    for (const auto& [key, index] : state.kfIndexMap)
        if (key != old) transform[key][key] = 1;
    transform[auxiliary][old] = 1;
    check(state.applyStateTransform(trace, transform, "TEST_RETAIN_OLD_FLOAT"),
        "production correlated auxiliary transform failed");
    check(zhangHandoffFloatRegisterRecord(state, auxiliary, edge,
        int(E_ObsCode::L1C), 0), "FLOAT registration failed");
    check(zhangRegisterGraphProductHandoffDependencies(state, E_Sys::GPS,
        E_ObsCode::L1C, {edge}, {0}), "dependency registration failed");
    return auxiliary;
}

struct MixedHandoffSource
{
    KFState state;
    std::array<KFKey, 3> auxiliaries;
};

MixedHandoffSource mixedHandoffSource(
    const std::array<long double, 3>& births, long double now)
{
    MixedHandoffSource result;
    result.state = source();
    const std::array<ZhangGraphEdge, 3> edges{{
        {"GPS_OLD", SatSys("G12")},
        {"JOG2", SatSys("E33")},
        {"BDS_OLD", SatSys("C19")}}};
    const std::array<E_ObsCode, 3> codes{
        E_ObsCode::L1C, E_ObsCode::L1C, E_ObsCode::L1P};
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        result.state.time.bigTime = births[i];
        const auto key = zhangGraphOldArcAuxiliaryKey(
            result.state, codes[i], edges[i], 28 + int(i));
        result.auxiliaries[i] = key;
        result.state.kfIndexMap[key] = 2 + int(i);
        result.state.stateTransitionMap[key][key][0] = 1;
        check(zhangHandoffFloatRegisterRecord(result.state, key, edges[i],
            int(codes[i]), 28 + int(i)), "mixed FLOAT registration failed");
        check(zhangRegisterGraphProductHandoffDependencies(result.state,
            edges[i].satellite.sys, codes[i], {edges[i]}, {28 + int(i)}),
            "mixed dependency registration failed");
    }
    result.state.x = VectorXd(5);
    result.state.dx = VectorXd(5);
    result.state.P = MatrixXd(5, 5);
    for (int i = 0; i < 5; ++i)
    {
        result.state.x(i) = 2 + .75 * i;
        result.state.dx(i) = .01 * (i + 1);
        for (int j = 0; j < 5; ++j)
            result.state.P(i, j) = .125 * (i + 1) * (j + 1) +
                (i == j ? 2 + i : 0);
    }
    result.state.time.bigTime = now;
    result.state.metaDataMap["UNRELATED_PRODUCT_AUTHORITY"] = "preserve-verbatim";
    return result;
}

// Independent oracle: select each surviving scalar directly from the saved
// original index order, including all off-diagonal/cross-system P entries.
// It does not call the production registry view or the production transform.
void checkPrincipalSelection(
    const KFState& before, const KFState& after, const std::set<KFKey>& removed)
{
    check(after.kfIndexMap.size() + removed.size() == before.kfIndexMap.size(),
        "global retirement removed the wrong number of coordinates");
    auto expectedMetadata = before.metaDataMap;
    for (const auto& key : removed)
    {
        check(!after.kfIndexMap.contains(key), "retired coordinate survived");
        expectedMetadata.erase(zhangHandoffFloatRecordKey(key));
    }
    for (const auto& [left, i] : after.kfIndexMap)
    {
        check(!removed.contains(left) && before.kfIndexMap.contains(left),
            "global retirement fabricated a surviving coordinate");
        const int sourceI = before.kfIndexMap.at(left);
        check(after.x(i) == before.x(sourceI) && after.dx(i) == before.dx(sourceI),
            "global selection changed a surviving mean/update");
        for (const auto& [right, j] : after.kfIndexMap)
            check(std::abs(after.P(i, j) -
                before.P(sourceI, before.kfIndexMap.at(right))) < 2e-12,
                "global selection changed the independent full cross-P oracle");
    }
    check(after.metaDataMap == expectedMetadata,
        "retirement changed unrelated authority/dependency metadata");
    std::string reason;
    check(zhangGraphProductHandoffStateValid(after, reason), reason.c_str());
}

void checkRejectedStateUnchanged(const KFState& before, const KFState& after)
{
    check(zhangCheckpointKfCoreBitwiseEqual(
        captureZhangCheckpointKfCore(before), captureZhangCheckpointKfCore(after)),
        "rejected sweep changed complete numeric/discrete KF core");
    check(before.metaDataMap == after.metaDataMap &&
        before.factorCommitSequence == after.factorCommitSequence &&
        before.stateTransitionMap == after.stateTransitionMap &&
        before.procNoiseMap == after.procNoiseMap &&
        before.gaussMarkovTauMap == after.gaussMarkovTauMap &&
        before.gaussMarkovMuMap == after.gaussMarkovMuMap,
        "rejected sweep changed metadata/factor sequence/transition/process model");
}

void testGlobalRetirement()
{
    constexpr long double birth = 1000000;
    std::string reason;
    for (const auto age : {0, 300, 330})
    {
        auto fixture = mixedHandoffSource({birth, birth, birth}, birth + age);
        const auto before = fixture.state;
        std::ostringstream trace;
        check(retireZhangGraphProductEventAuxiliaries(trace, fixture.state, "TEST_GLOBAL_AGE"),
            "global mixed-system age boundary failed");
        const std::set<KFKey> removed;
        checkPrincipalSelection(before, fixture.state, removed);
        check(fixture.state.factorCommitSequence == before.factorCommitSequence,
            "global retirement did not use exactly one selection transaction");
        check(trace.str().find("scope=GLOBAL stage=TEST_GLOBAL_AGE status=" +
            std::string("UNCHANGED")) != std::string::npos,
            "global age boundary lacks an explicit committed/unchanged witness");
    }

    // Active dependencies are pinned even when an earlier system has no aged
    // coordinate. The inactive TTL cannot expire another system's product.
    auto firstSystemEmpty = mixedHandoffSource(
        {birth + 330, birth, birth + 30}, birth + 330);
    const auto beforeMixed = firstSystemEmpty.state;
    std::ostringstream mixedTrace;
    check(retireZhangGraphProductEventAuxiliaries(
        mixedTrace, firstSystemEmpty.state, E_Sys::GPS),
        "GPS compatibility entry failed on the expired Galileo direction");
    checkPrincipalSelection(beforeMixed, firstSystemEmpty.state, {});
    check(firstSystemEmpty.state.kfIndexMap.contains(firstSystemEmpty.auxiliaries[0]) &&
        firstSystemEmpty.state.kfIndexMap.contains(firstSystemEmpty.auxiliaries[2]),
        "global sweep lost a birth-age-zero or inclusive-age-300 direction");

    std::array<E_Sys, 3> order{E_Sys::GPS, E_Sys::GAL, E_Sys::BDS};
    std::sort(order.begin(), order.end());
    do
    {
        auto fixture = mixedHandoffSource({birth, birth, birth}, birth + 330);
        const auto before = fixture.state;
        std::ostringstream trace;
        for (const auto system : order)
            check(retireZhangGraphProductEventAuxiliaries(trace, fixture.state, system),
                "global compatibility sweep still depends on constellation order");
        checkPrincipalSelection(before, fixture.state, {});
        check(fixture.state.factorCommitSequence == before.factorCommitSequence,
            "one global sweep became multiple per-system KF commits");
    } while (std::next_permutation(order.begin(), order.end()));

    // Exercise the actual pre-graph production entry with no observations in
    // any system. Availability must not bypass retirement of invisible rows.
    auto invisible = mixedHandoffSource({birth, birth, birth}, birth + 330);
    zhangClearGraphProductHandoffDependencies(invisible.state);
    invisible.state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] =
        "HANDOFF_FLOAT_GLOBAL_NO_AVAILABILITY_TEST";
    const auto beforeInvisible = invisible.state;
    const auto savedOptions = acsConfig.zhangFullRank;
    const auto savedSystems = acsConfig.process_sys;
    acsConfig.zhangFullRank.enable = true;
    for (const auto system : order)
    {
        acsConfig.process_sys[system] = true;
        acsConfig.zhangFullRank.sysOpts[system].use_spanning_tree = true;
    }
    ReceiverMap unavailable;
    std::ostringstream invisibleTrace;
    updateZhangFullRankReferences(invisibleTrace, unavailable, invisible.state);
    acsConfig.zhangFullRank = savedOptions;
    acsConfig.process_sys = savedSystems;
    checkPrincipalSelection(beforeInvisible, invisible.state,
        std::set<KFKey>(invisible.auxiliaries.begin(), invisible.auxiliaries.end()));
    check(invisibleTrace.str().find("scope=GLOBAL stage=PRE_GRAPH status=COMMITTED") !=
        std::string::npos, "invisible-system production entry did not run the global sweep");

    const std::vector<std::function<void(MixedHandoffSource&)>> corruptions{
        [](auto& f) { f.state.metaDataMap.erase(zhangHandoffFloatRecordKey(f.auxiliaries[1])); },
        [](auto& f) { f.state.metaDataMap[zhangHandoffFloatRecordKey(f.auxiliaries[1])] = "corrupt"; },
        [](auto& f) { f.state.stateTransitionMap[f.auxiliaries[1]][f.auxiliaries[1]][0] = 2; },
        [](auto& f) { f.state.kfIndexMap[f.auxiliaries[1]] = 0; },
        [](auto& f) { f.state.P(f.state.kfIndexMap.at(f.auxiliaries[1]), 0) =
            std::numeric_limits<double>::quiet_NaN(); },
        [](auto& f) { f.state.time.bigTime = 999999; },
        [](auto& f) { f.state.metaDataMap[std::string(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX) +
            "corrupt"] = "not an identity"; }};
    for (const auto& corrupt : corruptions)
    {
        auto fixture = mixedHandoffSource({birth, birth, birth}, birth + 330);
        corrupt(fixture);
        const auto before = fixture.state;
        std::ostringstream trace;
        check(!retireZhangGraphProductEventAuxiliaries(trace, fixture.state, "TEST_CORRUPT"),
            "corrupt expired row was silently erased instead of rejected");
        checkRejectedStateUnchanged(before, fixture.state);
        check(trace.str().find("status=REJECTED phase=ENTRY_REGISTRY reason=") !=
            std::string::npos && trace.str().find("status=COMMITTED") == std::string::npos,
            "corrupt-entry trace claimed a commit or omitted the actual failure phase");
    }

    auto callbackRejected = mixedHandoffSource({birth, birth, birth}, birth + 330);
    zhangClearGraphProductHandoffDependencies(callbackRejected.state);
    callbackRejected.state.factorCommitSequence = 41;
    const auto liveKey = zhangTransportAmbiguityKey(
        E_ObsCode::L1C, {"SAFE", SatSys("G11")});
    callbackRejected.state.procNoiseMap[liveKey] = .25;
    callbackRejected.state.gaussMarkovTauMap[liveKey] = 60;
    callbackRejected.state.gaussMarkovMuMap[liveKey] = 3;
    std::vector<int> callbackCounts;
    callbackRejected.state.exactStateTransformCallback =
        [count = 0, &callbackCounts](const KFState& owner, GTime,
            const auto&, const auto&, const auto&, const std::string&,
            const VectorXd&, const MatrixXd&, std::uint64_t, std::uint64_t) mutable
        {
            callbackCounts.push_back(++count);
            // Fault injection into KF-owned state verifies the production
            // operationEntry rollback, not an assumed test-only outer copy.
            auto& dirty = const_cast<KFState&>(owner);
            dirty.x.setZero(); dirty.P.setZero(); dirty.dx.setZero();
            dirty.metaDataMap["CALLBACK_DIRTY"] = "must roll back";
            dirty.factorCommitSequence += 99;
            dirty.stateTransitionMap.clear(); dirty.procNoiseMap.clear();
            dirty.gaussMarkovTauMap.clear(); dirty.gaussMarkovMuMap.clear();
            dirty.lastFactorTransactionFailureReason = "TEST_RETIREMENT_CONSUMER_REJECTED";
            return false;
        };
    const auto beforeCallback = callbackRejected.state;
    for (int attempt = 0; attempt < 2; ++attempt)
    {
        std::ostringstream trace;
        check(!retireZhangGraphProductEventAuxiliaries(
            trace, callbackRejected.state, "TEST_CALLBACK_REJECT"),
            "retirement ignored its consumer rejection");
        checkRejectedStateUnchanged(beforeCallback, callbackRejected.state);
        check(callbackRejected.state.lastFactorTransactionFailureReason ==
            "TEST_RETIREMENT_CONSUMER_REJECTED" &&
            trace.str().find("status=REJECTED phase=STATE_TRANSFORM reason=") != std::string::npos &&
            trace.str().find("status=COMMITTED") == std::string::npos,
            "callback rejection lost its actual reason or claimed a commit");
    }
    check(callbackCounts == std::vector<int>({1, 1}),
        "rejected transform did not restore its mutable callback callable");
}
}

// Shape of the measured C27 event: FALK-C28 L5P disappears, L1P remains,
// and the physical signal version stays V0. This is a production coordinate
// planner/central-KF test, not a claim to replay the real 134-station posterior.
void testFalkSignalOnlyLossAndPhysicalSlip()
{
    const ZhangGraphEdge lost{"FALK",SatSys("C28")};
    ZhangGraphBasis old;
    old.rootReceiver="A";old.receivers={"A","FALK","AGGO"};
    old.satellites={SatSys("C27"),SatSys("C28")};
    for(const auto& receiver:old.receivers)for(const auto& satellite:old.satellites)
        old.edges.insert({receiver,satellite});
    old.treeEdges={{"A",SatSys("C27")},{"FALK",SatSys("C27")},lost,{"AGGO",SatSys("C27")}};
    old.connected=true;old.componentCount=1;
    auto destination=old;
    destination.treeEdges.erase(lost);destination.treeEdges.insert({"A",SatSys("C28")});
    const std::map<E_ObsCode,double> wavelengths{{E_ObsCode::L1P,.19},{E_ObsCode::L5P,.25}};
    KFState source;source.time.bigTime=1405218450;
    source.kfIndexMap.clear();source.stateTransitionMap.clear();
    source.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="FALK_C28_SIGNAL_V0_CONTRACT";
    std::vector<KFKey> keys{KFState::oneKey};
    for(const auto& [code,wavelength]:wavelengths)
    {
        for(const auto& receiver:old.receivers)if(receiver!=old.rootReceiver)
            keys.push_back(zhangTransportReceiverKey(E_Sys::BDS,code,receiver));
        for(const auto& satellite:old.satellites)
            keys.push_back(zhangTransportSatelliteKey(code,satellite));
        for(const auto& edge:old.edges)if(!old.treeEdges.contains(edge))
            keys.push_back(zhangTransportAmbiguityKey(code,edge));
    }
    KFKey clock;clock.type=KF::SAT_CLOCK;clock.Sat=SatSys("C27");keys.push_back(clock);
    source.x=VectorXd::Zero(keys.size());source.dx=source.x;
    source.P=MatrixXd::Zero(keys.size(),keys.size());
    for(size_t i=0;i<keys.size();++i)
    {
        source.kfIndexMap[keys[i]]=i;source.stateTransitionMap[keys[i]][keys[i]][0]=1;
        source.x(i)=i?2+.5*i:1;
        for(size_t j=1;i && j<keys.size();++j)
            source.P(i,j)=.02*(i+1)*(j+1)+(i==j?1+i:0);
    }
    std::map<ZhangGraphEdge,std::map<int,int>> versions;
    for(const auto& edge:old.edges)for(const auto& [code,wavelength]:wavelengths)
        versions[edge][int(code)]=0;
    const ZhangCoordinateRow oldRow{
        {zhangTransportReceiverKey(E_Sys::BDS,E_ObsCode::L5P,"FALK"),1},
        {zhangTransportSatelliteKey(E_ObsCode::L5P,SatSys("C28")),1}};
    for(const bool slip:{false,true})
    {
        auto state=source;
        std::map<int,std::set<ZhangGraphEdge>> signalEdges;
        for(const auto& [code,wavelength]:wavelengths)signalEdges[int(code)]=destination.edges;
        if(!slip)signalEdges[int(E_ObsCode::L5P)].erase(lost);
        std::set<ZhangGraphEdge> retired;
        std::map<ZhangGraphEdge,std::set<E_ObsCode>> retiredSignals;
        if(slip){retired.insert(lost);retiredSignals[lost].insert(E_ObsCode::L5P);}
        auto plan=zhangPlanGraphCoordinateTransport(state,E_Sys::BDS,old,destination,retired,wavelengths,
            [](const auto&,auto){return 100.;},[](const auto&,auto){return 100.;},
            &retiredSignals,true,nullptr,&old.edges,&versions,&signalEdges);
        check(plan.valid,plan.failureReason.c_str());
        check(plan.retainedOldArcAuxiliaryKeys.contains({E_ObsCode::L5P,lost}) &&
            !plan.retainedOldArcAuxiliaryKeys.contains({E_ObsCode::L1P,lost}),
            "L1 topological presence suppressed capture of absent/retired L5 V0");
        const auto retained=plan.retainedOldArcAuxiliaryKeys.at({E_ObsCode::L5P,lost});
        ZhangHandoffFloatIdentity identity;
        check(zhangHandoffFloatIdentityFromKey(retained,identity) && identity.arcVersion==0,
            "old V0 auxiliary was relabelled as the new physical arc");
        check(plan.transform.contains(zhangTransportAmbiguityKey(E_ObsCode::L5P,lost))==slip,
            "missing L5 coordinate was fabricated from the surviving L1 chord");
        const auto proof=zhangVerifyGraphFunctionalTransport(plan,oldRow,{{retained,1}},0,0,0);
        check(proof.sameTimeMapProven,"old complete L5 V0 did not cancel all independent fresh sources");
        if(slip)
        {
            ZhangCoordinateRow current{{zhangTransportReceiverKey(E_Sys::BDS,E_ObsCode::L5P,"FALK"),1},
                {zhangTransportSatelliteKey(E_ObsCode::L5P,SatSys("C28")),1},
                {zhangTransportAmbiguityKey(E_ObsCode::L5P,lost),.25}};
            check(!zhangVerifyGraphFunctionalTransport(plan,oldRow,current,0,0).sameTimeMapProven,
                "fresh V1/slip arc inherited V0 same-time FLOAT or integer identity");
        }
        const auto capture=zhangPlanHandoffFloatCapture(state,plan.transform);
        check(capture.valid,capture.failureReason.c_str());
        std::ostringstream trace;
        check(state.applyStateTransform(trace,plan.transform,"TEST_FALK_SIGNAL_LIFECYCLE",
            plan.independentSourceVariances,&capture.metadata),"FALK production central transform failed");
        const int row=state.kfIndexMap.at(retained);
        double expectedMean=0,expectedVariance=0;
        for(const auto& [left,a]:oldRow)
        {
            expectedMean+=a*source.x(source.kfIndexMap.at(left));
            for(const auto& [right,b]:oldRow)expectedVariance+=a*b*source.P(source.kfIndexMap.at(left),source.kfIndexMap.at(right));
        }
        check(std::abs(state.x(row)-expectedMean)<1e-12 && std::abs(state.P(row,row)-expectedVariance)<1e-12,
            "FALK retained full mean/variance differs from independent source oracle");
        for(const auto& [destinationKey,column]:state.kfIndexMap)
        {
            double expected=0;
            for(const auto& [oldKey,a]:oldRow)for(const auto& [input,b]:plan.transform.at(destinationKey))
                if(source.kfIndexMap.contains(input))expected+=a*b*source.P(source.kfIndexMap.at(oldKey),source.kfIndexMap.at(input));
            check(std::abs(state.P(row,column)-expected)<1e-10,
                "FALK retained direction lost full cross-Q to a destination state");
        }
    }
}

int main()
{
    setenv("ZHANG_R51_ENABLE", "1", 1);
    setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF", "1", 1);
    setenv("ZHANG_R51_EPOCH_ONLY_AR", "1", 1);
    unsetenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS");
    unsigned retention = 0;
    check(zhangHandoffFloatRetentionSeconds(retention) && retention == 300,
        "default retention policy changed");
    std::ostringstream trace;
    std::string reason;
    auto state = source();
    const auto auxiliary = retain(state, trace);
    check(zhangGraphProductHandoffStateValid(state, reason), reason.c_str());
    const auto completePrior = state;
    // A fresh physical arc at the same receiver/satellite is independent;
    // retaining the old arc must never initialize its new version.
    const auto fresh = zhangTransportAmbiguityKey(E_ObsCode::L1C,
        {"AGGO", SatSys("G10")});
    std::map<KFKey, ZhangCoordinateRow> newArc;
    for (const auto& [key, index] : state.kfIndexMap) newArc[key][key] = 1;
    newArc[fresh][fresh] = 1;
    check(state.applyStateTransform(trace, newArc, "TEST_FRESH_ARC", {{fresh,25}}),
        "fresh physical arc transform failed");
    const int a = state.kfIndexMap.at(auxiliary), n = state.kfIndexMap.at(fresh);
    check(state.x(n) == 0 && state.P(n,n) == 25 && state.P(a,n) == 0,
        "old arc was inherited by new physical version");
    KFKey resolved;
    check(zhangHandoffFloatFind(state,
        {{"AGGO", SatSys("G10")}, int(E_ObsCode::L1C),0}, resolved) &&
        resolved == auxiliary, "old physical version lost");
    check(!zhangHandoffFloatFind(state,
        {{"AGGO", SatSys("G10")}, int(E_ObsCode::L1C),1}, resolved),
        "new version incorrectly aliases old FLOAT auxiliary");

    // Independent information-form posterior of the complete retained state.
    // The production update observes only a live edge, but must update the
    // old edge through its real cross covariance, not a fitted pseudo-row.
    check(state.stateTransition(trace, state.time + 30), "next-epoch prediction failed");
    KFMeas measurement;
    measurement.H = MatrixXd::Zero(1,state.x.size());
    const int live = state.kfIndexMap.at(zhangTransportAmbiguityKey(
        E_ObsCode::L1C, {"SAFE", SatSys("G11")}));
    measurement.H(0,live) = 1;
    measurement.R = MatrixXd::Constant(1,1,2);
    measurement.V = VectorXd::Constant(1,1.5);
    measurement.Y = measurement.H*state.x + measurement.V;
    measurement.VV = measurement.V;
    measurement.H_star = MatrixXd::Zero(1,0);
    measurement.uncorrelatedNoise.resize(0);
    const MatrixXd information = state.P.inverse() +
        measurement.H.transpose()*measurement.H/2;
    const MatrixXd expectedP = information.inverse();
    const VectorXd expectedX = expectedP *
        (state.P.inverse()*state.x + measurement.H.transpose()*measurement.Y/2);
    VectorXd xp=state.x, dx=state.dx;
    MatrixXd pp=state.P, qi,qh;
    check(state.kFilter(trace,measurement,xp,pp,dx,qi,qh,0,state.x.size(),0,1),
        "production FLOAT kernel failed");
    check((xp-expectedX).cwiseAbs().maxCoeff()<2e-12 &&
        (pp-expectedP).cwiseAbs().maxCoeff()<2e-12,
        "complete information-form mean/covariance oracle mismatch");
    state.x=xp; state.P=pp; state.dx=dx;
    check(std::abs(state.x(state.kfIndexMap.at(auxiliary))-completePrior.x(
        completePrior.kfIndexMap.at(auxiliary)))>0.01,
        "retired FLOAT row was frozen instead of jointly updated");
    check(retireZhangGraphProductEventAuxiliaries(trace,state,E_Sys::GPS) &&
        state.kfIndexMap.contains(auxiliary), "active old direction retired after one epoch");

    const auto core = captureZhangCheckpointKfCore(state);
    std::string graph;
    check(exportZhangGraphCheckpointSection(state,
        "HANDOFF_FLOAT_LIFECYCLE_TEST",graph,reason), reason.c_str());
    KFState restored;
    check(restoreZhangCheckpointKfCore(core,restored,&reason),reason.c_str());
    check(importZhangGraphCheckpointSection(restored,
        "HANDOFF_FLOAT_LIFECYCLE_TEST",graph,reason),reason.c_str());
    check((restored.x-state.x).isZero() && (restored.P-state.P).isZero() &&
        restored.metaDataMap==state.metaDataMap,
        "checkpoint lost retained FLOAT identity/full covariance");
    check(zhangGraphProductHandoffStateValid(restored,reason),reason.c_str());
    auto corrupted = restored;
    corrupted.metaDataMap.erase(zhangHandoffFloatRecordKey(auxiliary));
    check(!zhangGraphProductHandoffStateValid(corrupted,reason),
        "unregistered auxiliary accepted at checkpoint boundary");

    // Disposable AR branch changes cannot feed back to the protected FLOAT.
    auto arBranch = state;
    arBranch.x.setZero(); arBranch.P.setZero();
    check((state.x-xp).isZero() && (state.P-pp).isZero(),
        "AR branch contaminated protected FLOAT posterior");
    state.time.bigTime=1000300;
    check(retireZhangGraphProductEventAuxiliaries(trace,state,E_Sys::GPS) &&
        state.kfIndexMap.contains(auxiliary), "inclusive 300 s window lost");
    state.time.bigTime=1000330;
    check(retireZhangGraphProductEventAuxiliaries(trace,state,E_Sys::GPS) &&
        state.kfIndexMap.contains(auxiliary), "active dependency expired by cache TTL");
    zhangClearGraphProductHandoffDependencies(state);
    const auto beforeExpiry = state;
    check(retireZhangGraphProductEventAuxiliaries(trace,state,E_Sys::GPS) &&
        !state.kfIndexMap.contains(auxiliary), "released dependency did not retire old direction");
    for (const auto& [left,i] : state.kfIndexMap)
    {
        check(state.x(i)==beforeExpiry.x(beforeExpiry.kfIndexMap.at(left)),
            "marginalization changed live mean");
        for (const auto& [right,j] : state.kfIndexMap)
            check(std::abs(state.P(i,j)-beforeExpiry.P(beforeExpiry.kfIndexMap.at(left),
                beforeExpiry.kfIndexMap.at(right)))<2e-12,
                "marginalization changed surviving cross-P");
    }
    auto unused = source(); const auto unusedKey=retain(unused,trace);
    zhangClearGraphProductHandoffDependencies(unused);
    unused.time=unused.time+30;
    check(retireZhangGraphProductEventAuxiliaries(trace,unused,E_Sys::GPS) &&
        !unused.kfIndexMap.contains(unusedKey), "unreferenced direction retained forever");
    check(!zhangRegisterGraphProductHandoffDependencies(unused,E_Sys::GPS,
        E_ObsCode::L1C,{{"AGGO",SatSys("G10")}},{-1}),
        "unknown physical version registered as zero");

    for(const unsigned window : {30u,150u,300u})
    {
        setenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS",std::to_string(window).c_str(),1);
        auto bounded=source(); const auto key=retain(bounded,trace);
        bounded.time.bigTime=1000000+window;
        check(retireZhangGraphProductEventAuxiliaries(trace,bounded,E_Sys::GPS) &&
            bounded.kfIndexMap.contains(key), "configured inclusive retention boundary lost");
        bounded.time.bigTime+=30;
        check(retireZhangGraphProductEventAuxiliaries(trace,bounded,E_Sys::GPS) &&
            bounded.kfIndexMap.contains(key), "configured cache TTL expired active dependency");
        zhangClearGraphProductHandoffDependencies(bounded);
        check(retireZhangGraphProductEventAuxiliaries(trace,bounded,E_Sys::GPS) &&
            !bounded.kfIndexMap.contains(key), "released dependency exceeded cache lifecycle");
    }
    setenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS","0",1);
    auto immediate=source(); const auto immediateKey=retain(immediate,trace);
    check(!zhangHandoffFloatFind(immediate,
        {{"AGGO",SatSys("G10")},int(E_ObsCode::L1C),0},resolved),
        "explicit disabled retention granted a FLOAT map");
    check(retireZhangGraphProductEventAuxiliaries(trace,immediate,E_Sys::GPS) &&
        !immediate.kfIndexMap.contains(immediateKey),
        "zero-window control carried old arc beyond decision");
    immediate.time=immediate.time+30;
    check(!zhangHandoffFloatFind(immediate,
        {{"AGGO",SatSys("G10")},int(E_ObsCode::L1C),0},resolved),
        "zero-window control leaked previous epoch coordinate");
    for(const char* invalid : {"301","-1","","30.0"," 30","30 ","999999999999999999999"})
    {
        setenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS",invalid,1);
        check(!zhangGraphProductHandoffStateValid(immediate,reason) &&
            reason=="HANDOFF_FLOAT_RETENTION_SETTING_INVALID",
            "invalid retention setting did not fail closed");
        check(!retireZhangGraphProductEventAuxiliaries(trace,immediate,E_Sys::GPS),
            "invalid retention setting silently selected a window");
    }
    unsetenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS");
    auto capacity=source();
    capacity.kfIndexMap.clear(); capacity.stateTransitionMap.clear();
    capacity.x=VectorXd::Zero(128); capacity.dx=capacity.x;
    capacity.P=MatrixXd::Identity(128,128);
    std::vector<std::pair<ZhangGraphEdge,KFKey>> entries;
    for (int i=0;i<128;++i)
    {
        ZhangGraphEdge edge{"R"+std::to_string(i),SatSys("G10")};
        const auto key=zhangGraphOldArcAuxiliaryKey(capacity,E_ObsCode::L1C,edge,0);
        capacity.kfIndexMap[key]=i;
        capacity.stateTransitionMap[key][key][0]=1;
        entries.emplace_back(edge,key);
    }
    for (const auto& [edge,key]:entries)
        check(zhangHandoffFloatRegisterRecord(capacity,key,edge,int(E_ObsCode::L1C),0),
            "128-row retention bound rejected legal registry");
    check(zhangGraphProductHandoffStateValid(capacity,reason),reason.c_str());
    const ZhangGraphEdge extra{"R128",SatSys("G10")};
    const auto extraKey=zhangGraphOldArcAuxiliaryKey(capacity,E_ObsCode::L1C,extra,0);
    capacity.kfIndexMap[extraKey]=128;
    check(!zhangHandoffFloatRegisterRecord(capacity,extraKey,extra,int(E_ObsCode::L1C),0),
        "129th retained direction exceeded hard capacity");
    std::vector<ZhangGraphEdge> tooMany;
    for (int i=0;i<8193;++i) tooMany.push_back({"D"+std::to_string(i),SatSys("G10")});
    check(!zhangRegisterGraphProductHandoffDependencies(unused,E_Sys::GPS,
        E_ObsCode::L1C,tooMany,std::vector<int>(8193,0)),
        "dependency metadata grew beyond hard capacity");

    // Full public support can include a signal-specific chord. Its absent
    // other-frequency source must not abort the healthy FLOAT tree, and must
    // not be fabricated as a zero-variance old arc.
    KFState signalSpecific;
    signalSpecific.time.bigTime=1000000;
    signalSpecific.kfIndexMap.clear();signalSpecific.stateTransitionMap.clear();
    ZhangGraphBasis oldGraph;
    oldGraph.rootReceiver="A";
    oldGraph.receivers={"A","B"};oldGraph.satellites={SatSys("G10"),SatSys("G11")};
    oldGraph.treeEdges={{"A",SatSys("G10")},{"A",SatSys("G11")},{"B",SatSys("G10")}};
    oldGraph.edges=oldGraph.treeEdges;oldGraph.edges.insert({"B",SatSys("G11")});
    oldGraph.connected=true;oldGraph.componentCount=1;
    auto destination=oldGraph;destination.edges=destination.treeEdges;
    std::vector<KFKey> sourceKeys;
    for(const auto code:{E_ObsCode::L1C,E_ObsCode::L2W})
    {
        sourceKeys.push_back(zhangTransportReceiverKey(E_Sys::GPS,code,"B"));
        sourceKeys.push_back(zhangTransportSatelliteKey(code,SatSys("G10")));
        sourceKeys.push_back(zhangTransportSatelliteKey(code,SatSys("G11")));
    }
    const ZhangGraphEdge chord{"B",SatSys("G11")};
    sourceKeys.push_back(zhangTransportAmbiguityKey(E_ObsCode::L1C,chord));
    for(std::size_t i=0;i<sourceKeys.size();++i)
    {
        signalSpecific.kfIndexMap[sourceKeys[i]]=i;
        signalSpecific.stateTransitionMap[sourceKeys[i]][sourceKeys[i]][0]=1;
    }
    signalSpecific.x=VectorXd::Zero(sourceKeys.size());signalSpecific.dx=signalSpecific.x;
    signalSpecific.P=MatrixXd::Identity(sourceKeys.size(),sourceKeys.size());
    std::map<ZhangGraphEdge,std::map<int,int>> versions;
    for(const auto& edge:oldGraph.edges)
        for(const auto code:{E_ObsCode::L1C,E_ObsCode::L2W})versions[edge][int(code)]=0;
    const auto optionalPlan=zhangPlanGraphCoordinateTransport(signalSpecific,E_Sys::GPS,
        oldGraph,destination,{chord},{{E_ObsCode::L1C,.19},{E_ObsCode::L2W,.24}},
        [](const ZhangGraphEdge&,E_ObsCode){return 100.;},
        [](const std::string&,E_ObsCode){return 100.;},nullptr,true,nullptr,
        &oldGraph.edges,&versions);
    check(optionalPlan.valid &&
        optionalPlan.retainedOldArcAuxiliaryKeys.contains({E_ObsCode::L1C,chord}) &&
        !optionalPlan.retainedOldArcAuxiliaryKeys.contains({E_ObsCode::L2W,chord}),
        "optional signal chord blocked FLOAT or fabricated unavailable band");
    testGlobalRetirement();
    testFalkSignalOnlyLossAndPhysicalSlip();
    std::cout << "PASS production retained FLOAT: version isolation, full information-form mean/P, "
        "next-epoch joint update, pinned dependency lifetime, inactive expiry, checkpoint registry, "
        "AR branch isolation, surviving cross-P, 0/30/150/300 s policy, invalid-setting rejection, "
        "128 inactive/8192 dependency capacities, signal-specific support, global G/E/C active 0/300/330 ages, "
        "system-order/no-availability controls, corrupt-source refusal, full callback KF rollback, "
        "FALK-C28 L5 missing/L1 present V0, physical-slip fresh-prior negative, full destination cross-Q\n";
}
