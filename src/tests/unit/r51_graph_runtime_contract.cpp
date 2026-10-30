// Link against production PEA objects with run_r51_graph_runtime.py.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangPrivateGraphState.hpp"
#include "pea/zhangReference.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
void check(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}

void addPhaseState(KFState& state, const KFKey& key, double value)
{
    InitialState initial;
    initial.estimate = true;
    initial.x = value;
    initial.P = 1;
    initial.Q = 0;
    state.addKFState(key, initial);
}
}

int main()
{
    setenv("ZHANG_R51_ENABLE", "1", 1);
    setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF", "1", 1);
    setenv("ZHANG_R51_PRODUCT_REPAIR", "1", 1);
    acsConfig.zhangFullRank.enable = true;
    acsConfig.process_sys[E_Sys::GPS] = true;
    acsConfig.zhangPppAr.tree_slip_pivot_before_retire = true;
    acsConfig.zhangPppAr.tree_slip_shadow_replay = false;
    auto& options = acsConfig.zhangFullRank.sysOpts[E_Sys::GPS];
    options.use_spanning_tree = true;
    options.reference_receiver = "A";
    options.reference_satellite = "G01";
    options.baseline_observables = {E_ObsCode::L1C, E_ObsCode::L2W};
    options.product_core_min_satellite_support = 0;

    ReceiverMap receivers;
    const std::set<ZhangGraphEdge> physicalEdges{
        {"A", SatSys("G01")}, {"A", SatSys("G02")},
        {"B", SatSys("G01")}, {"B", SatSys("G02")}};
    for (const auto& edge : physicalEdges)
    {
        auto& receiver = receivers[edge.receiver];
        receiver.id = edge.receiver;
        receiver.ready = true;
        auto observation = std::make_shared<GObs>();
        observation->Sat = edge.satellite;
        auto& status = receiver.satStatMap[edge.satellite];
        status.el = 1.5;
        observation->satStat_ptr = &status;
        for (const auto code : options.baseline_observables)
        {
            const auto frequency = code2Freq[E_Sys::GPS][code];
            auto& signal = observation->sigs[frequency];
            signal.code = code;
            signal.L = 1e8;
            signal.P = 2e7;
            auto& config = acsConfig.getRecOpts(edge.receiver, {"GPS", enum_to_string(code)});
            config.ambiguity.sigma = {1000};
            config.phase_bias.sigma = {1000};
        }
        receiver.obsList.push_back(observation);
    }

    std::ostringstream trace;
    KFState authority;
    const std::string rootId = "R51_GRAPH_RUNTIME_CONTRACT";
    authority.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] = rootId;
    // Clean startup must remain valid before any phase coordinates exist.
    updateZhangFullRankReferences(trace, receivers, authority);
    ZhangGraphIntegerContext context;
    check(zhangGraphIntegerContext(authority, E_Sys::GPS, context),
        "clean graph initialization failed");
    for (const auto code : options.baseline_observables)
    {
        for (const auto& receiver : context.basis.receivers)
            if (receiver != context.basis.rootReceiver)
                addPhaseState(authority,
                    zhangTransportReceiverKey(E_Sys::GPS, code, receiver), 20);
        for (const auto& satellite : context.basis.satellites)
            addPhaseState(authority, zhangTransportSatelliteKey(code, satellite), 10);
        for (const auto& edge : context.basis.edges)
            if (!context.basis.treeEdges.contains(edge))
                addPhaseState(authority, zhangTransportAmbiguityKey(code, edge), 7);
    }
    check(authority.stateTransition(trace, authority.time), "phase state initialization failed");
    std::string reason;
    ZhangFloatGraphAuthorityToken before;
    check(captureZhangFloatGraphAuthorityToken(authority, before, reason), reason.c_str());

    const auto ownershipNumericalBefore = authority;
    auto cloneProduction = [](const KFState& source, const KFState& destination, std::string& why)
    {
        return cloneZhangGraphRuntime(source, destination, &why);
    };
    auto unchangedAuthority = [&]()
    {
        check(zhangFloatGraphAuthorityUnchanged(authority, before, reason),
            "prepared private state lifecycle changed authority graph");
        check((authority.x - ownershipNumericalBefore.x).norm() == 0 &&
            (authority.P - ownershipNumericalBefore.P).norm() == 0,
            "prepared private state lifecycle changed authority moments");
    };
    KFState releasedPrivate(authority);
    {
        ZhangPrivateGraphState prepared;
        check(zhangMakePrivateGraphState(authority, rootId, "final-product", cloneProduction,
            eraseZhangGraphRuntime, prepared, reason), reason.c_str());
        // This is the same owned pointer type used by prepared publication.
        // The source has a nonzero seven-cycle chord on both signals.
        check(prepared->x.isApprox(authority.x, 0), "private branch did not copy source moments");
        releasedPrivate = *prepared;
        prepared->x(1) += 3;
        prepared->P(1, 1) *= .5;
    }
    unchangedAuthority();
    ZhangFloatGraphAuthorityToken releasedToken;
    check(!captureZhangFloatGraphAuthorityToken(releasedPrivate, releasedToken, reason),
        "normal prepared destruction retained private graph");
    auto earlyReturn = [&]()
    {
        ZhangPrivateGraphState prepared;
        check(zhangMakePrivateGraphState(authority, rootId, "final-product", cloneProduction,
            eraseZhangGraphRuntime, prepared, reason), reason.c_str());
        releasedPrivate = *prepared;
        return false;
    };
    check(!earlyReturn(), "injected early return did not execute");
    unchangedAuthority();
    check(!captureZhangFloatGraphAuthorityToken(releasedPrivate, releasedToken, reason),
        "early return retained private graph");
    bool injectedException = false;
    try
    {
        ZhangPrivateGraphState prepared;
        check(zhangMakePrivateGraphState(authority, rootId, "final-product", cloneProduction,
            eraseZhangGraphRuntime, prepared, reason), reason.c_str());
        releasedPrivate = *prepared;
        throw std::runtime_error("INJECTED_PRIVATE_CONDITIONING_EXCEPTION");
    }
    catch (const std::runtime_error& error)
    {
        injectedException = std::string(error.what()) == "INJECTED_PRIVATE_CONDITIONING_EXCEPTION";
    }
    check(injectedException, "private conditioning exception did not execute");
    unchangedAuthority();
    check(!captureZhangFloatGraphAuthorityToken(releasedPrivate, releasedToken, reason),
        "conditioning exception retained private graph");
    {
        ZhangPrivateGraphState prepared;
        auto failedClone = [&](const KFState& source, const KFState& destination, std::string& why)
        {
            check(cloneProduction(source, destination, why), why.c_str());
            releasedPrivate = destination;
            why = "INJECTED_CLONE_FAILURE_AFTER_REGISTRY_INSERT";
            return false;
        };
        check(!zhangMakePrivateGraphState(authority, rootId, "final-product", failedClone,
            eraseZhangGraphRuntime, prepared, reason), "injected clone failure was accepted");
        check(!prepared, "failed clone published private state");
    }
    unchangedAuthority();
    check(!captureZhangFloatGraphAuthorityToken(releasedPrivate, releasedToken, reason),
        "failed clone retained allocated private graph");

    KFState alias(authority);
    check(!cloneZhangGraphRuntime(authority, alias, &reason), "same-ID clone was accepted");
    check(reason == "GRAPH_SOURCE_DESTINATION_RUNTIME_ALIAS", "same-ID rejection reason changed");
    check(zhangFloatGraphAuthorityUnchanged(authority, before, reason),
        "failed same-ID clone changed authority");

    KFState disposable(authority);
    disposable.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA] = rootId + "/branch/private";
    check(cloneZhangGraphRuntime(authority, disposable, &reason), reason.c_str());
    ZhangFloatGraphAuthorityToken privateBefore;
    check(captureZhangFloatGraphAuthorityToken(disposable, privateBefore, reason), reason.c_str());
    KFState missing(authority);
    missing.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] = rootId + "-missing";
    check(!cloneZhangGraphRuntime(missing, disposable, &reason), "missing source graph was accepted");
    check(zhangFloatGraphAuthorityUnchanged(disposable, privateBefore, reason),
        "missing source clone destroyed existing destination");
    KFState invalid(authority);
    invalid.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] = "invalid\nidentity";
    check(!cloneZhangGraphRuntime(invalid, disposable, &reason), "invalid source identity was accepted");
    check(zhangFloatGraphAuthorityUnchanged(disposable, privateBefore, reason),
        "invalid source clone changed destination");

    // An initialized GPS registry must not hide a missing Galileo chart.
    auto& galileo = acsConfig.zhangFullRank.sysOpts[E_Sys::GAL];
    galileo.use_spanning_tree = true;
    galileo.baseline_observables = {E_ObsCode::L1C};
    acsConfig.process_sys[E_Sys::GAL] = true;
    KFState incompleteMultiSystem(authority);
    addPhaseState(incompleteMultiSystem,
        zhangTransportSatelliteKey(E_ObsCode::L1C, SatSys("E01")), 12);
    check(incompleteMultiSystem.stateTransition(trace, incompleteMultiSystem.time),
        "multi-system state initialization failed");
    check(!cloneZhangGraphRuntime(incompleteMultiSystem, disposable, &reason),
        "GPS runtime hid missing Galileo runtime");
    check(zhangFloatGraphAuthorityUnchanged(disposable, privateBefore, reason),
        "multi-system rejection changed destination");
    check(!captureZhangFloatGraphAuthorityToken(incompleteMultiSystem, releasedToken, reason),
        "FLOAT token accepted estimated Galileo state without its graph");
    galileo.use_spanning_tree = false;
    acsConfig.process_sys[E_Sys::GAL] = false;
    eraseZhangGraphRuntime(disposable);
    check(zhangFloatGraphAuthorityUnchanged(authority, before, reason),
        "private graph cleanup erased FLOAT owner");

    // Non-graph full-rank configurations legitimately have no graph registry.
    options.use_spanning_tree = false;
    check(cloneZhangGraphRuntime(missing, disposable, &reason), "non-graph clone no-op was rejected");
    check(reason == "NO_GRAPH_CONFIGURED", "non-graph no-op was not explicit");
    check(captureZhangFloatGraphAuthorityToken(missing, releasedToken, reason),
        "non-graph token was rejected");
    check(zhangFloatGraphAuthorityUnchanged(missing, releasedToken, reason),
        "non-graph token did not remain valid");
    options.use_spanning_tree = true;

    std::string checkpoint;
    check(exportZhangGraphCheckpointSection(authority, rootId, checkpoint, reason), reason.c_str());
    const auto numericalBefore = authority;
    eraseZhangGraphRuntime(authority);
    check(!zhangFloatGraphAuthorityUnchanged(authority, before, reason),
        "FLOAT token did not detect erased graph");
    bool rejected = false;
    try
    {
        updateZhangFullRankReferences(trace, receivers, authority);
    }
    catch (const std::runtime_error& error)
    {
        rejected = std::string(error.what()) ==
            "ZHANG_FLOAT_GRAPH_RUNTIME_MISSING_WITH_ESTIMATED_PHASE_STATE";
    }
    check(rejected, "lost graph silently initialized a different S tree");
    check((authority.x - numericalBefore.x).norm() == 0 &&
        (authority.P - numericalBefore.P).norm() == 0,
        "lost-graph guard changed FLOAT moments");
    check(!captureZhangFloatGraphAuthorityToken(authority, privateBefore, reason),
        "lost-graph guard installed replacement runtime");

    check(importZhangGraphCheckpointSection(authority, rootId, checkpoint, reason), reason.c_str());
    check(zhangFloatGraphAuthorityUnchanged(authority, before, reason),
        "validated checkpoint did not restore exact FLOAT graph");
    authority.time = authority.time + 30;
    updateZhangFullRankReferences(trace, receivers, authority);
    check(captureZhangFloatGraphAuthorityToken(authority, privateBefore, reason),
        "validated checkpoint could not continue initialized graph");
    check(!zhangFloatGraphAuthorityUnchanged(authority, before, reason),
        "FLOAT token ignored physical epoch and admission changes");
    eraseZhangGraphRuntime(authority);
    std::cout << "PASS production graph runtime: clean startup, alias rejection, missing/invalid source, "
        "actual prepared cleanup/early return/conditioning exception/failed clone isolation, "
        "explicit non-graph no-op, lost-graph abort, exact checkpoint restore\n";
}
