// Native regression: real KF transform, global lifecycle and checkpoint code.
// The oracle selects saved scalar coordinates and every full cross-P entry;
// it does not invoke the production capture planner to compute expectations.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductHandoffFloatRegistry.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "pea/zhangReference.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
void check(bool condition, const std::string& reason)
{ if (!condition) throw std::runtime_error(reason); }

struct Candidate
{
    ZhangGraphEdge edge;
    E_ObsCode code = E_ObsCode::L1C;
    int version = 0;
    KFKey source;
    KFKey auxiliary;
    bool active = false;
};
struct Fixture
{
    KFState state;
    std::vector<Candidate> pending;
    std::map<KFKey, KFKey> origin;
};

Fixture fixture(int inactive, int active)
{
    Fixture result;
    auto& state = result.state;
    state.kfIndexMap.clear(); state.stateTransitionMap.clear();
    state.time.bigTime = 1000000;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] = "ACTIVE_FLOAT_RETENTION_TEST";
    const auto add = [&](const KFKey& key)
    {
        const int index = state.kfIndexMap.size();
        state.kfIndexMap[key] = index;
        state.stateTransitionMap[key][key][0] = 1;
        result.origin[key] = key;
    };
    for (int i = 0; i < inactive; ++i)
    {
        const ZhangGraphEdge edge{"INACTIVE_" + std::to_string(i), SatSys("G10")};
        const auto key = zhangGraphOldArcAuxiliaryKey(state, E_ObsCode::L1C, edge, 0);
        add(key);
    }
    const std::array<SatSys, 3> satellites{SatSys("G12"), SatSys("E11"), SatSys("C39")};
    for (int i = 0; i < active; ++i)
    {
        Candidate candidate;
        candidate.edge = {"ACTIVE_" + std::to_string(i), satellites[i % 3]};
        candidate.code = i % 2 ? E_ObsCode::L2W : E_ObsCode::L1C;
        candidate.version = 17 + i;
        candidate.active = true;
        candidate.source = zhangTransportAmbiguityKey(candidate.code, candidate.edge);
        candidate.auxiliary = zhangGraphOldArcAuxiliaryKey(
            state, candidate.code, candidate.edge, candidate.version);
        add(candidate.source);
        result.origin[candidate.auxiliary] = candidate.source;
        result.pending.push_back(candidate);
        check(zhangRegisterGraphProductHandoffDependencies(state,
            candidate.edge.satellite.sys, candidate.code, {candidate.edge}, {candidate.version}),
            "active dependency registration failed");
    }
    for (const auto& [key, index] : state.kfIndexMap)
        if (key.type == KF::PRODUCT_HANDOFF)
        {
            ZhangHandoffFloatIdentity identity;
            check(zhangHandoffFloatIdentityFromKey(key, identity), "fixture identity invalid");
            check(zhangHandoffFloatRegisterRecord(state, key, identity.edge,
                identity.observable, identity.arcVersion), "fixture inactive record rejected");
        }
    const int count = state.kfIndexMap.size();
    state.x.resize(count); state.dx.resize(count); state.P.resize(count, count);
    for (int i = 0; i < count; ++i)
    {
        state.x(i) = .25 * (i + 1); state.dx(i) = .001 * (i + 1);
        for (int j = 0; j < count; ++j)
            state.P(i,j) = .001 * (i + 1) * (j + 1) + (i == j ? 2 + i : 0);
    }
    return result;
}

void capture(Fixture& fixture, const std::vector<Candidate>& candidates)
{
    auto& state = fixture.state;
    std::map<KFKey, ZhangCoordinateRow> transform;
    for (const auto& [key, index] : state.kfIndexMap) transform[key][key] = 1;
    for (const auto& candidate : candidates)
    {
        transform.erase(candidate.source);
        transform[candidate.auxiliary][candidate.source] = 1;
    }
    const auto plan = zhangPlanHandoffFloatCapture(state, transform);
    check(plan.valid, plan.failureReason);
    for (const auto& key : plan.discarded) transform.erase(key);
    std::ostringstream trace;
    check(state.applyStateTransform(trace, transform, "TEST_GLOBAL_ACTIVE_CAPTURE"),
        "native capture transform failed");
    state.metaDataMap = plan.metadata;
    for (const auto& candidate : candidates)
        if (state.kfIndexMap.contains(candidate.auxiliary))
        {
            ZhangHandoffFloatRecord record;
            check(zhangHandoffFloatReadRecord(state, candidate.auxiliary, record) &&
                record.identity.edge == candidate.edge && record.identity.observable == int(candidate.code) &&
                record.identity.arcVersion == candidate.version,
                "capture plan did not prevalidate every selected new registry record");
        }
    std::string reason;
    check(zhangGraphProductHandoffStateValid(state, reason), reason);
}

void oracle(const Fixture& original, const Fixture& after)
{
    for (const auto& [left, i] : after.state.kfIndexMap)
    {
        const int oldI = original.state.kfIndexMap.at(after.origin.at(left));
        check(after.state.x(i) == original.state.x(oldI) &&
            after.state.dx(i) == original.state.dx(oldI), "capture changed retained mean/dx");
        for (const auto& [right, j] : after.state.kfIndexMap)
        {
            const int oldJ = original.state.kfIndexMap.at(after.origin.at(right));
            check(std::abs(after.state.P(i,j) - original.state.P(oldI, oldJ)) < 2e-12,
                "independent full cross-system P scalar oracle mismatch");
        }
    }
}

void noIntegerAuthority(const KFState& state)
{
    for (const auto& [name, value] : state.metaDataMap)
        check(name.starts_with(ZHANG_HANDOFF_FLOAT_RECORD_PREFIX) ||
            name.starts_with(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX) ||
            name == ZHANG_CHECKPOINT_RUNTIME_ID_METADATA,
            "FLOAT capture fabricated non-FLOAT metadata authority");
    for (const auto& [key, index] : state.kfIndexMap)
    {
        if (key.type != KF::PRODUCT_HANDOFF) continue;
        ZhangHandoffFloatRecord record;
        check(zhangHandoffFloatReadRecord(state, key, record), "integer negative fixture lacks FLOAT record");
        const auto sat = record.identity.edge.satellite;
        const auto code = static_cast<E_ObsCode>(record.identity.observable);
        ZhangGraphProductEventReceipt receipt;
        check(!zhangGraphProductEventReceipt(state, sat.sys, code, sat, receipt),
            "retaining a FLOAT key fabricated a product map receipt");
        ZhangProductIntegerTransition transition;
        transition.valid = true;
        transition.observable = code;
        transition.physicalEdges = {record.identity.edge};
        transition.physicalArcVersions = {record.identity.arcVersion};
        transition.coefficients = {ZhangExactInteger(1)};
        std::string reason;
        check(!zhangGraphProductEventTransitionWitness(state, sat.sys, code,
            sat, SatSys(sat.sys, 1), transition, reason) &&
            reason == "CURRENT_SIGNAL_FACTOR_ROWS_NOT_FINALIZED",
            "retained FLOAT key alone authorized an integer bridge");
    }
}

void roundTrip(const KFState& state)
{
    std::string graph, reason;
    check(exportZhangGraphCheckpointSection(state, "ACTIVE_FLOAT_RETENTION_TEST", graph, reason), reason);
    KFState restored;
    check(restoreZhangCheckpointKfCore(captureZhangCheckpointKfCore(state), restored, &reason), reason);
    check(importZhangGraphCheckpointSection(restored, "ACTIVE_FLOAT_RETENTION_TEST", graph, reason), reason);
    check(restored.metaDataMap == state.metaDataMap && restored.kfIndexMap == state.kfIndexMap &&
        (restored.x - state.x).isZero() && (restored.P - state.P).isZero(),
        "checkpoint changed active identities or complete correlated KF moments");
    check(zhangGraphProductHandoffStateValid(restored, reason), reason);
}

void duplicateBirthIdentity()
{
    auto value = fixture(0, 1);
    capture(value, value.pending);
    const auto candidate = value.pending.front();
    value.state.time.bigTime += 30;
    const auto laterBirth = zhangGraphOldArcAuxiliaryKey(value.state,
        candidate.code, candidate.edge, candidate.version);
    std::map<KFKey, ZhangCoordinateRow> transform;
    transform[candidate.auxiliary][candidate.auxiliary] = 1;
    transform[laterBirth][candidate.auxiliary] = 1;
    auto plan = zhangPlanHandoffFloatCapture(value.state, transform);
    check(plan.valid && plan.active == 1 && plan.aliases.at(laterBirth) == candidate.auxiliary &&
        plan.discarded.contains(laterBirth) &&
        plan.discardReasons.at(laterBirth) == "EXACT_DUPLICATE_COORDINATE",
        "duplicate physical identity acquired two active coordinates");
    transform[laterBirth][candidate.auxiliary] = 2;
    plan = zhangPlanHandoffFloatCapture(value.state, transform);
    check(!plan.valid && plan.failureReason == "HANDOFF_FLOAT_DUPLICATE_IDENTITY_UNPROVEN",
        "different random variables were merged on physical identity alone");
    // Metadata may not conceal a malformed birth as a successful new record.
    auto corrupt = value.state;
    corrupt.metaDataMap[zhangHandoffFloatRecordKey(candidate.auxiliary)] = "broken";
    const auto corruptCore = captureZhangCheckpointKfCore(corrupt);
    const auto corruptMetadata = corrupt.metaDataMap;
    check(!zhangHandoffFloatRegisterRecord(corrupt, candidate.auxiliary,
        candidate.edge, int(candidate.code), candidate.version),
        "registration silently repaired corrupt identity evidence");
    check(corrupt.metaDataMap == corruptMetadata && zhangCheckpointKfCoreBitwiseEqual(
        corruptCore, captureZhangCheckpointKfCore(corrupt)),
        "rejected registration mutated x/P or original metadata");
}

void malformedCaptureIsAtomic()
{
    const auto original = fixture(0, 1);
    for (const bool badBirth : {true, false})
    {
        auto value = original;
        const auto candidate = value.pending.front();
        KFKey invalid;
        if (badBirth)
        {
            invalid = candidate.auxiliary;
            invalid.str = "ZHANG_OLD_ARC:" + candidate.edge.receiver + "@" +
                std::to_string(candidate.version) + "#BAD";
        }
        else
        {
            auto otherEpoch = value.state;
            otherEpoch.time.bigTime += 30;
            invalid = zhangGraphOldArcAuxiliaryKey(otherEpoch,
                candidate.code, candidate.edge, candidate.version);
        }
        std::map<KFKey, ZhangCoordinateRow> transform{{invalid, {{candidate.source, 1}}}};
        const auto before = captureZhangCheckpointKfCore(value.state);
        const auto metadata = value.state.metaDataMap;
        const auto plan = zhangPlanHandoffFloatCapture(value.state, transform);
        check(!plan.valid && plan.failureReason == "HANDOFF_FLOAT_CAPTURE_NEW_IDENTITY_OR_EPOCH_INVALID",
            "malformed or wrong-epoch capture was accepted before KF commit");
        check(value.state.metaDataMap == metadata && zhangCheckpointKfCoreBitwiseEqual(
            before, captureZhangCheckpointKfCore(value.state)),
            "rejected capture changed original complete KF or metadata");
        // Independent direct registration negative: its coordinate exists,
        // so rejection must be identity/birth validation, not key absence.
        auto registration = value.state;
        registration.kfIndexMap.clear(); registration.stateTransitionMap.clear();
        registration.kfIndexMap[invalid] = 0;
        registration.stateTransitionMap[invalid][invalid][0] = 1;
        const auto registrationBefore = captureZhangCheckpointKfCore(registration);
        const auto registrationMetadata = registration.metaDataMap;
        check(!zhangHandoffFloatRegisterRecord(registration, invalid,
            candidate.edge, int(candidate.code), candidate.version),
            "invalid birth was accepted by direct registration");
        check(registration.metaDataMap == registrationMetadata && zhangCheckpointKfCoreBitwiseEqual(
            registrationBefore, captureZhangCheckpointKfCore(registration)),
            "rejected direct registration changed metadata or KF");
    }
}

void dependencyDirectoryBounds()
{
    auto value = fixture(0, 1);
    std::vector<ZhangGraphEdge> edges;
    for (std::size_t i = 1; i < ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES; ++i)
        edges.push_back({"DIRECTORY_" + std::to_string(i), SatSys("G12")});
    check(zhangRegisterGraphProductHandoffDependencies(value.state, E_Sys::GPS,
        E_ObsCode::L1C, edges, std::vector<int>(edges.size(), 0)),
        "production dependency registration rejected the inclusive 8192 identity bound");
    std::set<ZhangHandoffFloatIdentity> dependencies;
    std::string reason;
    check(zhangHandoffFloatReadDependencies(value.state.metaDataMap, dependencies, reason) &&
        dependencies.size() == ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES,
        "read/registration disagree at the 8192 unique dependency boundary");
    const auto before = captureZhangCheckpointKfCore(value.state);
    const auto metadata = value.state.metaDataMap;
    check(!zhangRegisterGraphProductHandoffDependencies(value.state, E_Sys::GPS,
        E_ObsCode::L1C, {{"DIRECTORY_OVERFLOW", SatSys("G12")}}, {0}),
        "production registration exceeded the 8192 identity directory");
    check(value.state.metaDataMap == metadata && zhangCheckpointKfCoreBitwiseEqual(
        before, captureZhangCheckpointKfCore(value.state)),
        "over-limit registration changed original KF or dependency metadata");
    // A corrupt imported 8193-entry directory must fail independently at
    // planner and strict checkpoint/AR-exit validation, never return a subset.
    const ZhangHandoffFloatIdentity overflow{{"IMPORTED_OVERFLOW", SatSys("G12")}, int(E_ObsCode::L1C), 0};
    const auto text = zhangHandoffFloatIdentityText(overflow);
    value.state.metaDataMap[std::string(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX) + text] = text;
    const auto plan = zhangPlanHandoffFloatCapture(value.state,
        {{value.pending.front().auxiliary, {{value.pending.front().source, 1}}}});
    check(!plan.valid && plan.failureReason == "HANDOFF_FLOAT_DEPENDENCY_CAPACITY_EXCEEDED",
        "capture accepted a corrupt over-limit dependency directory");
    check(!zhangGraphProductHandoffStateValid(value.state, reason) &&
        reason == "HANDOFF_FLOAT_DEPENDENCY_CAPACITY_EXCEEDED",
        "strict checkpoint/AR-exit validation accepted over-limit dependencies");
}

void captureDiscardReasons()
{
    auto expired = fixture(1, 0);
    expired.state.time.bigTime += 301;
    const auto old = expired.state.kfIndexMap.begin()->first;
    auto plan = zhangPlanHandoffFloatCapture(expired.state, {{old, {{old, 1}}}});
    check(plan.valid && plan.discarded.contains(old) &&
        plan.discardReasons.at(old) == "INACTIVE_CACHE_EXPIRED",
        "expired inactive capture was misreported as capacity exhaustion");
    auto full = fixture(128, 0);
    const ZhangGraphEdge extra{"NEW_INACTIVE_OVERFLOW", SatSys("G10")};
    const auto extraKey = zhangGraphOldArcAuxiliaryKey(full.state, E_ObsCode::L1C, extra, 0);
    std::map<KFKey, ZhangCoordinateRow> transform;
    for (const auto& [key, index] : full.state.kfIndexMap) transform[key][key] = 1;
    transform[extraKey][full.state.kfIndexMap.begin()->first] = 1;
    plan = zhangPlanHandoffFloatCapture(full.state, transform);
    check(plan.valid && plan.inactive == 128 && plan.discarded.size() == 1 &&
        plan.discardReasons.begin()->second == "INACTIVE_CACHE_CAPACITY",
        "inactive cache overflow lacks a distinct capacity diagnostic");
}

void inactiveBatchOrder()
{
    auto initial = fixture(128, 6);
    // Six additional genuinely distinct physical source coordinates. Inactive
    // events arrive first, in the same epoch, and cannot reserve active slots.
    std::vector<Candidate> extra;
    for (int i = 0; i < 6; ++i)
    {
        Candidate candidate = initial.pending[i];
        candidate.edge.receiver = "NEW_INACTIVE_" + std::to_string(i);
        candidate.version += 1000; candidate.active = false;
        candidate.source = zhangTransportAmbiguityKey(candidate.code, candidate.edge);
        candidate.auxiliary = zhangGraphOldArcAuxiliaryKey(initial.state,
            candidate.code, candidate.edge, candidate.version);
        const int index = initial.state.kfIndexMap.size();
        initial.state.kfIndexMap[candidate.source] = index;
        initial.state.stateTransitionMap[candidate.source][candidate.source][0] = 1;
        initial.origin[candidate.source] = candidate.source;
        initial.origin[candidate.auxiliary] = candidate.source;
        extra.push_back(candidate);
    }
    const int count = initial.state.kfIndexMap.size();
    initial.state.x.resize(count); initial.state.dx.resize(count); initial.state.P.resize(count, count);
    for (int i = 0; i < count; ++i)
    {
        initial.state.x(i) = .25 * (i + 1); initial.state.dx(i) = .001 * (i + 1);
        for (int j = 0; j < count; ++j)
            initial.state.P(i,j) = .001 * (i + 1) * (j + 1) + (i == j ? 2 + i : 0);
    }
    // Independent expected identities: every known requested auxiliary plus
    // the first 128 inactive keys, without invoking the production planner.
    std::set<KFKey> expectedKeys;
    std::vector<KFKey> inactiveKeys;
    for (const auto& [key, index] : initial.state.kfIndexMap)
        if (key.type == KF::PRODUCT_HANDOFF) inactiveKeys.push_back(key);
    for (const auto& candidate : initial.pending) expectedKeys.insert(candidate.auxiliary);
    for (const auto& candidate : extra) inactiveKeys.push_back(candidate.auxiliary);
    std::sort(inactiveKeys.begin(), inactiveKeys.end());
    expectedKeys.insert(inactiveKeys.begin(), inactiveKeys.begin() + 128);
    std::array<int, 3> order{0, 1, 2};
    do
    {
        auto batched = initial;
        for (const int system : order)
        {
            std::vector<Candidate> inactiveBatch, activeBatch;
            for (int i = system; i < 6; i += 3)
            { inactiveBatch.push_back(extra[i]); activeBatch.push_back(initial.pending[i]); }
            capture(batched, inactiveBatch);
            capture(batched, activeBatch);
        }
        std::set<KFKey> actualKeys;
        for (const auto& [key, index] : batched.state.kfIndexMap) actualKeys.insert(key);
        check(actualKeys == expectedKeys, "global cache depends on constellation processing order");
        oracle(initial, batched);
        noIntegerAuthority(batched.state);
    } while (std::next_permutation(order.begin(), order.end()));
}
}

int main()
{
    setenv("ZHANG_R51_ENABLE", "1", 1);
    setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF", "1", 1);
    setenv("ZHANG_R51_EPOCH_ONLY_AR", "1", 1);
    unsetenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS");
    for (const auto counts : {std::pair{128, 2}, std::pair{0, 130}})
    {
        auto value = fixture(counts.first, counts.second);
        const auto original = value;
        capture(value, value.pending);
        check(value.state.kfIndexMap.size() == std::size_t(counts.first + counts.second),
            "old inactive cache displaced newly requested active coordinates");
        oracle(original, value);
        noIntegerAuthority(value.state);
        value.state.time.bigTime += 900;
        std::ostringstream trace;
        check(retireZhangGraphProductEventAuxiliaries(trace, value.state, "TEST_ACTIVE_AGE_900"),
            "active age greater than cache TTL was rejected");
        check(value.state.kfIndexMap.size() == std::size_t(counts.second),
            "active pin or inactive lifecycle was violated");
        for (const auto& candidate : value.pending)
        {
            KFKey found;
            check(zhangHandoffFloatFind(value.state,
                {candidate.edge, int(candidate.code), candidate.version}, found) &&
                found == candidate.auxiliary, "Find disagrees with active lifecycle/checkpoint boundary");
        }
        oracle(original, value);
        roundTrip(value.state);
        // The projected validator may not accept an active-key marginalization.
        std::string reason;
        const auto key = value.pending.front().auxiliary;
        auto metadata = value.state.metaDataMap;
        metadata.erase(zhangHandoffFloatRecordKey(key));
        check(!zhangHandoffFloatRegistryViewValid(value.state, metadata, {key}, false, reason) &&
            reason == "HANDOFF_FLOAT_ACTIVE_DEPENDENCY_REMOVAL_FORBIDDEN",
            "projected registry allowed an active coordinate to disappear");
        // Keep an unrelated physical state so releasing all active auxiliaries
        // is a proper non-empty principal marginal rather than an empty KF.
        const auto live = zhangTransportAmbiguityKey(E_ObsCode::L1C, {"LIVE", SatSys("G03")});
        std::map<KFKey, ZhangCoordinateRow> keep;
        for (const auto& [saved, index] : value.state.kfIndexMap) keep[saved][saved] = 1;
        keep[live][live] = 1;
        check(value.state.applyStateTransform(trace, keep, "TEST_RELEASE_LIVE", {{live, 3}}),
            "release fixture augmentation failed");
        zhangClearGraphProductHandoffDependencies(value.state);
        check(retireZhangGraphProductEventAuxiliaries(trace, value.state, "TEST_RELEASE"),
            "released active dependency did not retire");
        check(value.state.kfIndexMap.size() == 1 && value.state.kfIndexMap.contains(live),
            "dependency deregistration left unbounded history");
    }
    inactiveBatchOrder();
    duplicateBirthIdentity();
    malformedCaptureIsAtomic();
    dependencyDirectoryBounds();
    captureDiscardReasons();
    auto disabled = fixture(0, 2);
    setenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS", "0", 1);
    std::map<KFKey, ZhangCoordinateRow> transform;
    for (const auto& [key, index] : disabled.state.kfIndexMap) transform[key][key] = 1;
    for (const auto& candidate : disabled.pending) transform[candidate.auxiliary][candidate.source] = 1;
    const auto plan = zhangPlanHandoffFloatCapture(disabled.state, transform);
    check(plan.valid && plan.explicitlyDisabled && plan.active == 0 && plan.inactive == 0 &&
        plan.discarded.size() == 2, "explicit zero retention granted active FLOAT capture");
    for (const auto& [key, reason] : plan.discardReasons)
        check(reason == "EXPLICIT_RETENTION_DISABLED", "explicit disable reason was obscured");
    unsetenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS");
    std::cout << "PASS native active FLOAT retention: existing128+active2, active130, cache order "
        "G/E/C all6 permutations, new same-epoch inactive, age900, release, full scalar cross-P "
        "oracle, checkpoint registry, active-drop rejection, structurally valid integer negative, "
        "BAD/wrong-epoch atomic refusal, production directory8192/8193, explicit0\n";
}
