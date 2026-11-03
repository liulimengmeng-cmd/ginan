#pragma once

#include "common/algebra.hpp"
#include "common/zhangFullRank.hpp"
#include "common/zhangProductFloatCarrierMetadata.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// These records contain FLOAT coordinate identities, never integer values,
// accepted AR rows, or pseudo-observations.  Keeping the actual KF coordinate
// is essential: an independent cached mean/variance would lose cross-P.
inline constexpr const char* ZHANG_HANDOFF_FLOAT_RECORD_PREFIX =
    "ZHANG_HANDOFF_FLOAT_V1/";
inline constexpr const char* ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX =
    "ZHANG_HANDOFF_FLOAT_DEP_V1/";
// This is an INACTIVE cache budget, not a budget for active public-product
// coordinates. Active identities are bounded by the dependency directory.
inline constexpr std::size_t ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES = 128;
inline constexpr std::size_t ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES = 8192;
inline constexpr long double ZHANG_HANDOFF_FLOAT_MAX_AGE_SECONDS = 300;

inline bool zhangHandoffFloatRetentionSeconds(
    unsigned& seconds, std::string* failureReason = nullptr)
{
    seconds = unsigned(ZHANG_HANDOFF_FLOAT_MAX_AGE_SECONDS);
    const char* setting = std::getenv("ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS");
    if (!setting) return true;
    const auto reject = [&]()
    {
        if (failureReason) *failureReason = "HANDOFF_FLOAT_RETENTION_SETTING_INVALID";
        return false;
    };
    if (*setting == '\0') return reject();
    seconds = 0;
    for (const char* p = setting; *p; ++p)
    {
        if (*p < '0' || *p > '9') return reject();
        seconds = seconds * 10 + unsigned(*p - '0');
        // Reject before an unsigned overflow can turn a huge configuration
        // into an apparently small allowance. The hard upper bound is fixed.
        if (seconds > unsigned(ZHANG_HANDOFF_FLOAT_MAX_AGE_SECONDS)) return reject();
    }
    return true;
}

struct ZhangHandoffFloatIdentity
{
    ZhangGraphEdge edge;
    int observable = 0;
    int arcVersion = -1;

    bool operator<(const ZhangHandoffFloatIdentity& other) const
    {
        return std::tie(edge, observable, arcVersion) <
            std::tie(other.edge, other.observable, other.arcVersion);
    }
};

struct ZhangHandoffFloatRecord
{
    ZhangHandoffFloatIdentity identity;
    long double birthTime = 0;
};

inline bool zhangHandoffFloatIdentityFromKey(
    const KFKey& key, ZhangHandoffFloatIdentity& identity)
{
    constexpr const char* prefix = "ZHANG_OLD_ARC:";
    if (key.type != KF::PRODUCT_HANDOFF || !key.str.starts_with(prefix))
        return false;
    const auto hash = key.str.rfind('#');
    const auto at = key.str.rfind('@', hash);
    if (hash == std::string::npos || at == std::string::npos ||
        at <= std::char_traits<char>::length(prefix)) return false;
    try
    {
        std::size_t end = 0;
        const auto value = key.str.substr(at + 1, hash - at - 1);
        identity.arcVersion = std::stoi(value, &end);
        if (end != value.size() || identity.arcVersion < 0) return false;
    }
    catch (...) { return false; }
    identity.edge = {key.str.substr(std::char_traits<char>::length(prefix),
        at - std::char_traits<char>::length(prefix)), key.Sat};
    identity.observable = key.num;
    return key.Sat.prn > 0;
}

inline std::string zhangHandoffFloatRecordKey(const KFKey& key)
{
    return std::string(ZHANG_HANDOFF_FLOAT_RECORD_PREFIX) +
        key.Sat.id() + "/" + std::to_string(key.num) + "/" + key.str;
}

inline std::string zhangHandoffFloatIdentityText(
    const ZhangHandoffFloatIdentity& identity)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::quoted(identity.edge.receiver) << ' '
        << std::quoted(identity.edge.satellite.id()) << ' '
        << identity.observable << ' ' << identity.arcVersion;
    return out.str();
}

inline bool zhangHandoffFloatParseIdentity(
    std::istream& in, ZhangHandoffFloatIdentity& identity)
{
    std::string satellite;
    if (!(in >> std::quoted(identity.edge.receiver) >> std::quoted(satellite)
              >> identity.observable >> identity.arcVersion) ||
        identity.edge.receiver.empty() || identity.arcVersion < 0)
        return false;
    identity.edge.satellite = SatSys(satellite.c_str());
    return identity.edge.satellite.prn > 0 &&
        identity.edge.satellite.sys != E_Sys::NONE &&
        identity.edge.satellite.id() == satellite;
}

inline bool zhangHandoffFloatReadRecord(
    const std::map<std::string, std::string>& metadata,
    const KFKey& key, ZhangHandoffFloatRecord& record)
{
    if (key.type != KF::PRODUCT_HANDOFF) return false;
    const auto item = metadata.find(zhangHandoffFloatRecordKey(key));
    if (item == metadata.end()) return false;
    std::istringstream in(item->second);
    in.imbue(std::locale::classic());
    std::string birth;
    if (!zhangHandoffFloatParseIdentity(in, record.identity) || !(in >> birth))
        return false;
    // stold accepts hexadecimal floating point, preserving long-double epoch
    // identity.  Iostream extraction of hexfloat is not portable on libstdc++.
    try
    {
        std::size_t end = 0;
        record.birthTime = std::stold(birth, &end);
        if (end != birth.size()) return false;
    }
    catch (...) { return false; }
    in >> std::ws;
    if (!in.eof() || !std::isfinite(record.birthTime) ||
        record.identity.edge.satellite != key.Sat ||
        record.identity.observable != key.num)
        return false;
    std::ostringstream birthIdentity;
    birthIdentity.imbue(std::locale::classic());
    birthIdentity << std::hexfloat
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << record.birthTime;
    return key.str == "ZHANG_OLD_ARC:" + record.identity.edge.receiver + "@" +
        std::to_string(record.identity.arcVersion) + "#" + birthIdentity.str();
}

inline bool zhangHandoffFloatReadRecord(
    const KFState& state, const KFKey& key, ZhangHandoffFloatRecord& record)
{
    return zhangHandoffFloatReadRecord(state.metaDataMap, key, record);
}

inline bool zhangHandoffFloatReadDependencies(
    const std::map<std::string, std::string>& metadata,
    std::set<ZhangHandoffFloatIdentity>& result, std::string& failureReason)
{
    result.clear();
    std::size_t count = 0;
    for (const auto& [name, value] : metadata)
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX))
        {
            if (++count > ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES)
            { failureReason = "HANDOFF_FLOAT_DEPENDENCY_CAPACITY_EXCEEDED"; return false; }
            std::istringstream input(value);
            input.imbue(std::locale::classic());
            ZhangHandoffFloatIdentity identity;
            if (!zhangHandoffFloatParseIdentity(input, identity))
            { failureReason = "HANDOFF_FLOAT_INVALID_DEPENDENCY"; return false; }
            input >> std::ws;
            if (!input.eof())
            { failureReason = "HANDOFF_FLOAT_INVALID_DEPENDENCY"; return false; }
            result.insert(std::move(identity));
        }
    // Whole-target consumers pin their actual correlated RAW coordinates as
    // well. A missing old anchor is not permission for cache eviction: until
    // that active request is resolved/retired, conservatively pin all RAW
    // records. These are FLOAT lifetime dependencies, never integer facts.
    std::vector<ZhangProductFloatTargetRequest> targets;
    if (!zhangReadFloatTargets(metadata, targets, failureReason)) return false;
    bool unresolvedActive = false;
    for (const auto& target : targets)
    {
        if (!target.active) continue;
        if (!target.numericBindingComplete) { unresolvedActive = true; continue; }
        for (const auto& [key, coefficient] : target.numericRow)
        {
            if (coefficient == 0 || key.type != KF::PRODUCT_HANDOFF ||
                zhangIsProductFloatCarrierKey(key)) continue;
            ZhangHandoffFloatIdentity identity;
            if (!zhangHandoffFloatIdentityFromKey(key, identity))
            { failureReason = "HANDOFF_FLOAT_TARGET_RAW_IDENTITY_INVALID"; return false; }
            result.insert(std::move(identity));
        }
    }
    if (unresolvedActive)
        for (const auto& [name, value] : metadata)
            if (name.starts_with(ZHANG_HANDOFF_FLOAT_RECORD_PREFIX))
            {
                std::istringstream input(value);input.imbue(std::locale::classic());
                ZhangHandoffFloatIdentity identity;
                if (!zhangHandoffFloatParseIdentity(input, identity))
                { failureReason = "HANDOFF_FLOAT_UNBOUND_TARGET_PIN_RECORD_INVALID"; return false; }
                result.insert(std::move(identity));
            }
    if (result.size() > ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES)
    { failureReason = "HANDOFF_FLOAT_DEPENDENCY_CAPACITY_EXCEEDED"; return false; }
    return true;
}

inline std::set<ZhangHandoffFloatIdentity> zhangHandoffFloatDependencies(
    const KFState& state)
{
    std::set<ZhangHandoffFloatIdentity> result;
    std::string reason;
    // Invalid directories fail closed, rather than granting a subset of the
    // active identities while silently ignoring a corrupt dependency.
    if (!zhangHandoffFloatReadDependencies(state.metaDataMap, result, reason))
        result.clear();
    return result;
}

inline bool zhangHandoffFloatRecordCurrent(
    const KFState& state, const ZhangHandoffFloatRecord& record,
    unsigned retention, const std::set<ZhangHandoffFloatIdentity>& dependencies)
{
    const long double age = state.time.bigTime - record.birthTime;
    return retention > 0 && std::isfinite(age) && age >= 0 &&
        (dependencies.contains(record.identity) || age <= retention);
}

inline bool zhangHandoffFloatRecordCurrent(
    const KFState& state, const ZhangHandoffFloatRecord& record)
{
    unsigned retention = 0;
    if (!zhangHandoffFloatRetentionSeconds(retention) || retention == 0) return false;
    const long double age = state.time.bigTime - record.birthTime;
    if (!std::isfinite(age) || age < 0) return false;
    std::set<ZhangHandoffFloatIdentity> dependencies;
    std::string reason;
    if (!zhangHandoffFloatReadDependencies(state.metaDataMap, dependencies, reason))
        return false;
    return zhangHandoffFloatRecordCurrent(state, record, retention, dependencies);
}

inline bool zhangHandoffFloatMakeCurrentRecord(
    const KFState& state, const KFKey& key,
    const ZhangHandoffFloatIdentity& identity, std::string& value)
{
    if (!std::isfinite(state.time.bigTime) || key.type != KF::PRODUCT_HANDOFF ||
        key.Sat != identity.edge.satellite || key.num != identity.observable ||
        identity.arcVersion < 0) return false;
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << zhangHandoffFloatIdentityText(identity) << ' ' << std::hexfloat
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << state.time.bigTime;
    value = output.str();
    const std::map<std::string, std::string> proposed{
        {zhangHandoffFloatRecordKey(key), value}};
    ZhangHandoffFloatRecord record;
    return zhangHandoffFloatReadRecord(proposed, key, record) &&
        record.birthTime == state.time.bigTime;
}

inline bool zhangHandoffFloatRegisterRecord(
    KFState& state, const KFKey& key,
    const ZhangGraphEdge& edge, int observable, int version)
{
    if (key.type != KF::PRODUCT_HANDOFF || key.Sat != edge.satellite ||
        key.num != observable || version < 0 ||
        !state.kfIndexMap.contains(key)) return false;
    ZhangHandoffFloatRecord existing;
    if (zhangHandoffFloatReadRecord(state, key, existing))
        return existing.identity.edge == edge &&
            existing.identity.observable == observable &&
            existing.identity.arcVersion == version;
    if (state.metaDataMap.contains(zhangHandoffFloatRecordKey(key))) return false;
    std::set<ZhangHandoffFloatIdentity> dependencies;
    std::string reason;
    if (!zhangHandoffFloatReadDependencies(state.metaDataMap, dependencies, reason))
        return false;
    std::size_t inactive = 0;
    std::set<ZhangHandoffFloatIdentity> identities;
    for (const auto& [candidate, index] : state.kfIndexMap)
    {
        if (candidate.type != KF::PRODUCT_HANDOFF) continue;
        if (zhangIsProductFloatCarrierKey(candidate)) continue;
        ZhangHandoffFloatIdentity physical;
        if (!zhangHandoffFloatIdentityFromKey(candidate, physical) ||
            !identities.insert(physical).second) return false;
        inactive += !dependencies.contains(physical);
    }
    if (inactive > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES) return false;
    std::string value;
    if (!zhangHandoffFloatMakeCurrentRecord(
            state, key, {edge, observable, version}, value)) return false;
    state.metaDataMap[zhangHandoffFloatRecordKey(key)] = std::move(value);
    return true;
}

inline bool zhangHandoffFloatFind(
    const KFState& state, const ZhangHandoffFloatIdentity& identity, KFKey& key)
{
    unsigned retention = 0;
    std::set<ZhangHandoffFloatIdentity> dependencies;
    std::string reason;
    if (!zhangHandoffFloatRetentionSeconds(retention) || retention == 0 ||
        !zhangHandoffFloatReadDependencies(state.metaDataMap, dependencies, reason)) return false;
    for (const auto& [candidate, index] : state.kfIndexMap)
    {
        if (candidate.type != KF::PRODUCT_HANDOFF) continue;
        ZhangHandoffFloatRecord record;
        if (zhangHandoffFloatReadRecord(state, candidate, record) &&
            zhangHandoffFloatRecordCurrent(state, record, retention, dependencies) &&
            !(identity < record.identity) && !(record.identity < identity))
        {
            key = candidate;
            return true;
        }
    }
    return false;
}

struct ZhangHandoffFloatCapturePlan
{
    bool valid = false;
    std::string failureReason;
    std::set<KFKey> discarded;
    std::map<KFKey, std::string> discardReasons;
    std::map<KFKey, KFKey> aliases;
    std::map<std::string, std::string> metadata;
    std::size_t active = 0;
    std::size_t inactive = 0;
    bool explicitlyDisabled = false;
};

// Select from ALL existing and new destinations together. Keeping the lowest
// inactive KFKeys makes repeated constellation batches equivalent to one
// global selection. Active physical identities are never charged to cache128.
// This is retention of actual correlated FLOAT coordinates, not an integer
// certificate or a target-space rank compression.
inline ZhangHandoffFloatCapturePlan zhangPlanHandoffFloatCapture(
    const KFState& state,
    const std::map<KFKey, std::map<KFKey, double>>& transform)
{
    ZhangHandoffFloatCapturePlan plan;
    unsigned retention = 0;
    if (!zhangHandoffFloatRetentionSeconds(retention, &plan.failureReason)) return plan;
    std::set<ZhangHandoffFloatIdentity> dependencies;
    if (!zhangHandoffFloatReadDependencies(state.metaDataMap, dependencies,
            plan.failureReason)) return plan;
    plan.explicitlyDisabled = retention == 0;
    plan.metadata = state.metaDataMap;
    std::set<std::string> sourceRecords;
    std::set<ZhangHandoffFloatIdentity> sourceIdentities;
    for (const auto& [key, index] : state.kfIndexMap)
    {
        if (key.type != KF::PRODUCT_HANDOFF) continue;
        if (zhangIsProductFloatCarrierKey(key))
        {
            if (!zhangProductFloatCarrierCoordinateValid(state, state.metaDataMap, key, plan.failureReason) ||
                !transform.contains(key) || transform.at(key) != std::map<KFKey,double>{{key,1}})
            { plan.failureReason = "FLOAT_CARRIER_CAPTURE_MUST_PRESERVE_TYPED_COORDINATE"; return plan; }
            continue;
        }
        ZhangHandoffFloatRecord record;
        if (!zhangHandoffFloatReadRecord(state, key, record))
        { plan.failureReason = "HANDOFF_FLOAT_CAPTURE_SOURCE_REGISTRY_INVALID"; return plan; }
        if (!sourceIdentities.insert(record.identity).second)
        { plan.failureReason = "HANDOFF_FLOAT_CAPTURE_SOURCE_DUPLICATE_IDENTITY"; return plan; }
        const auto age = state.time.bigTime - record.birthTime;
        if (!std::isfinite(age) || age < 0)
        { plan.failureReason = "HANDOFF_FLOAT_FUTURE_OR_NONFINITE_BIRTH"; return plan; }
        sourceRecords.insert(zhangHandoffFloatRecordKey(key));
        if (transform.contains(key)) continue;
        if (!plan.explicitlyDisabled && dependencies.contains(record.identity))
        { plan.failureReason = "HANDOFF_FLOAT_ACTIVE_CAPTURE_DESTINATION_MISSING"; return plan; }
        plan.discarded.insert(key);
        plan.discardReasons[key] = plan.explicitlyDisabled ? "EXPLICIT_RETENTION_DISABLED" :
            age > retention ? "INACTIVE_CACHE_EXPIRED" : "INACTIVE_DESTINATION_NOT_REQUESTED";
        plan.metadata.erase(zhangHandoffFloatRecordKey(key));
    }
    for (const auto& [name, value] : state.metaDataMap)
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_RECORD_PREFIX) &&
            !sourceRecords.contains(name))
        { plan.failureReason = "HANDOFF_FLOAT_ORPHAN_REGISTRY_ENTRY"; return plan; }
    std::map<ZhangHandoffFloatIdentity, KFKey> canonical;
    for (const auto& [key, row] : transform)
    {
        if (key.type != KF::PRODUCT_HANDOFF) continue;
        if (zhangIsProductFloatCarrierKey(key)) continue;
        ZhangHandoffFloatIdentity identity;
        if (!zhangHandoffFloatIdentityFromKey(key, identity))
        { plan.failureReason = "HANDOFF_FLOAT_CAPTURE_IDENTITY_INVALID"; return plan; }
        if (state.kfIndexMap.contains(key))
        {
            ZhangHandoffFloatRecord record;
            if (!zhangHandoffFloatReadRecord(state, key, record))
            { plan.failureReason = "HANDOFF_FLOAT_CAPTURE_SOURCE_REGISTRY_INVALID"; return plan; }
            const auto age = state.time.bigTime - record.birthTime;
            if (!std::isfinite(age) || age < 0)
            { plan.failureReason = "HANDOFF_FLOAT_FUTURE_OR_NONFINITE_BIRTH"; return plan; }
        }
        else
        {
            std::string value;
            if (!zhangHandoffFloatMakeCurrentRecord(state, key, identity, value))
            { plan.failureReason = "HANDOFF_FLOAT_CAPTURE_NEW_IDENTITY_OR_EPOCH_INVALID"; return plan; }
            plan.metadata[zhangHandoffFloatRecordKey(key)] = std::move(value);
        }
        auto [entry, inserted] = canonical.emplace(identity, key);
        if (!inserted)
        {
            const auto old = entry->second;
            if (transform.at(old) != row)
            {
                // Same receiver/satellite/code/version is not sufficient to
                // merge two different random variables or two physical arcs.
                plan.failureReason = "HANDOFF_FLOAT_DUPLICATE_IDENTITY_UNPROVEN";
                return plan;
            }
            // Prefer an already registered coordinate, then lexical key.
            if (!state.kfIndexMap.contains(old) && state.kfIndexMap.contains(key))
                entry->second = key;
        }
    }
    std::set<KFKey> retained;
    std::set<KFKey> inactive;
    for (const auto& [identity, key] : canonical)
        if (!plan.explicitlyDisabled && dependencies.contains(identity))
        { retained.insert(key); ++plan.active; }
        else
        {
            if (state.kfIndexMap.contains(key))
            {
                ZhangHandoffFloatRecord record;
                zhangHandoffFloatReadRecord(state, key, record);
                const auto age = state.time.bigTime - record.birthTime;
                if (!std::isfinite(age) || age < 0)
                { plan.failureReason = "HANDOFF_FLOAT_FUTURE_OR_NONFINITE_BIRTH"; return plan; }
                if (age > retention) continue;
            }
            inactive.insert(key);
        }
    if (!plan.explicitlyDisabled)
        for (const auto& key : inactive)
            if (plan.inactive < ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES)
            { retained.insert(key); ++plan.inactive; }
    for (const auto& [key, row] : transform)
    {
        if (key.type != KF::PRODUCT_HANDOFF) continue;
        if (zhangIsProductFloatCarrierKey(key)) continue;
        ZhangHandoffFloatIdentity identity;
        zhangHandoffFloatIdentityFromKey(key, identity); // validated above
        const auto& chosen = canonical.at(identity);
        if (retained.contains(chosen) && key != chosen) plan.aliases[key] = chosen;
        if (!retained.contains(key))
        {
            plan.discarded.insert(key);
            std::string reason = plan.explicitlyDisabled ? "EXPLICIT_RETENTION_DISABLED" :
                plan.aliases.contains(key) ? "EXACT_DUPLICATE_COORDINATE" : "INACTIVE_CACHE_CAPACITY";
            if (!plan.explicitlyDisabled && !plan.aliases.contains(key) &&
                state.kfIndexMap.contains(key))
            {
                ZhangHandoffFloatRecord record;
                zhangHandoffFloatReadRecord(state, key, record);
                if (state.time.bigTime - record.birthTime > retention)
                    reason = "INACTIVE_CACHE_EXPIRED";
            }
            plan.discardReasons[key] = std::move(reason);
            plan.metadata.erase(zhangHandoffFloatRecordKey(key));
        }
    }
    plan.valid = true;
    plan.failureReason = "DEPENDENCY_PINNED_FLOAT_CAPTURE_PLANNED";
    return plan;
}

// A read-only principal-submatrix view: no second KF or dense covariance is
// constructed. The maintenance entry may admit an otherwise valid expired
// record only so the global sweep can remove it. Future births and corrupt
// records are never retirement candidates. Public validation stays strict.
inline bool zhangHandoffFloatRegistryViewValid(
    const KFState& state,
    const std::map<std::string, std::string>& metadata,
    const std::set<KFKey>& removedKeys,
    bool allowRetirableExpired,
    std::string& failureReason)
{
    unsigned retention = 0;
    if (!zhangHandoffFloatRetentionSeconds(retention, &failureReason)) return false;
    std::set<ZhangHandoffFloatIdentity> dependencies;
    if (!zhangHandoffFloatReadDependencies(metadata, dependencies, failureReason))
        return false;
    std::size_t inactive = 0;
    std::set<std::string> expected;
    std::set<ZhangHandoffFloatIdentity> identities;
    if (state.x.rows() != state.P.rows() || state.P.rows() != state.P.cols() ||
        state.dx.rows() != state.x.rows() ||
        state.kfIndexMap.size() != static_cast<std::size_t>(state.x.rows()))
    { failureReason = "HANDOFF_FLOAT_MOMENT_DIMENSION_MISMATCH"; return false; }
    if (!std::isfinite(state.time.bigTime))
    { failureReason = "HANDOFF_FLOAT_EPOCH_INVALID"; return false; }
    for (const auto& key : removedKeys)
    {
        if (key.type != KF::PRODUCT_HANDOFF || !state.kfIndexMap.contains(key))
        { failureReason = "HANDOFF_FLOAT_RETIREMENT_SELECTION_INVALID"; return false; }
        ZhangHandoffFloatRecord record;
        if (retention > 0 && zhangHandoffFloatReadRecord(state, key, record) &&
            dependencies.contains(record.identity))
        { failureReason = "HANDOFF_FLOAT_ACTIVE_DEPENDENCY_REMOVAL_FORBIDDEN"; return false; }
    }
    std::set<int> sourceIndices;
    std::vector<int> survivingIndices;
    for (const auto& [key, index] : state.kfIndexMap)
    {
        if (index < 0 || index >= state.x.rows() || !sourceIndices.insert(index).second)
        { failureReason = "HANDOFF_FLOAT_STATE_INDEX_INVALID"; return false; }
        if (removedKeys.contains(key)) continue;
        survivingIndices.push_back(index);
        if (!std::isfinite(state.x(index)) || !std::isfinite(state.dx(index)))
        { failureReason = "HANDOFF_FLOAT_STATE_MOMENT_NONFINITE"; return false; }
        if (key.type != KF::PRODUCT_HANDOFF) continue;
        if (zhangIsProductFloatCarrierKey(key))
        {
            if (!zhangProductFloatCarrierCoordinateValid(state, metadata, key, failureReason)) return false;
            continue;
        }
        ZhangHandoffFloatRecord record;
        if (!zhangHandoffFloatReadRecord(metadata, key, record) ||
            !identities.insert(record.identity).second ||
            !state.stateTransitionMap.contains(key))
        { failureReason = "HANDOFF_FLOAT_IDENTITY_OR_MOMENT_INVALID"; return false; }
        const long double age = state.time.bigTime - record.birthTime;
        if (!std::isfinite(age) || age < 0)
        { failureReason = "HANDOFF_FLOAT_FUTURE_OR_NONFINITE_BIRTH"; return false; }
        const bool active = dependencies.contains(record.identity);
        inactive += !active;
        if (!allowRetirableExpired &&
            (retention == 0 || (!active && age > retention)))
        { failureReason = "HANDOFF_FLOAT_IDENTITY_OR_MOMENT_INVALID"; return false; }
        const auto& transition = state.stateTransitionMap.at(key);
        if (transition.size() != 1 || !transition.contains(key) ||
            transition.at(key).size() != 1 ||
            !transition.at(key).contains(0) || transition.at(key).at(0) != 1 ||
            (state.procNoiseMap.contains(key) && state.procNoiseMap.at(key) != 0) ||
            (state.gaussMarkovTauMap.contains(key) && state.gaussMarkovTauMap.at(key) > 0))
        { failureReason = "HANDOFF_FLOAT_NONCONSTANT_AUXILIARY_MODEL"; return false; }
        expected.insert(zhangHandoffFloatRecordKey(key));
    }
    // Inspect the actual retained principal submatrix, including every
    // cross-system/off-diagonal covariance. A diagonal-only view could hide
    // a corrupt moment behind an otherwise legal retirement selection.
    for (const int right : survivingIndices)
        for (const int left : survivingIndices)
            if (!std::isfinite(state.P(left, right)))
            { failureReason = "HANDOFF_FLOAT_COVARIANCE_NONFINITE"; return false; }
    for (const auto& [name, value] : metadata)
    {
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_RECORD_PREFIX) &&
            !expected.contains(name))
        { failureReason = "HANDOFF_FLOAT_ORPHAN_REGISTRY_ENTRY"; return false; }
    }
    if (!allowRetirableExpired && inactive > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES)
    { failureReason = "HANDOFF_FLOAT_RETENTION_CAPACITY_EXCEEDED"; return false; }
    failureReason = "DEPENDENCY_PINNED_CORRELATED_FLOAT_REGISTRY_VALID";
    return true;
}

inline bool zhangHandoffFloatRegistryValid(
    const KFState& state, std::string& failureReason)
{
    return zhangHandoffFloatRegistryViewValid(
        state, state.metaDataMap, {}, false, failureReason);
}
