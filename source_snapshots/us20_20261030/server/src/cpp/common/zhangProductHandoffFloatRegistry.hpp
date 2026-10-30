#pragma once

#include "common/algebra.hpp"
#include "common/zhangFullRank.hpp"

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

inline bool zhangHandoffFloatRecordCurrent(
    const KFState& state, const ZhangHandoffFloatRecord& record)
{
    unsigned retention = 0;
    if (!zhangHandoffFloatRetentionSeconds(retention)) return false;
    const long double age = state.time.bigTime - record.birthTime;
    return std::isfinite(age) && age >= 0 &&
        age <= retention;
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
    std::size_t count = 0;
    for (const auto& [candidate, index] : state.kfIndexMap)
        count += candidate.type == KF::PRODUCT_HANDOFF;
    if (count > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES) return false;
    ZhangHandoffFloatIdentity identity{edge, observable, version};
    std::ostringstream value;
    value.imbue(std::locale::classic());
    value << zhangHandoffFloatIdentityText(identity) << ' ' << std::hexfloat
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << state.time.bigTime;
    state.metaDataMap[zhangHandoffFloatRecordKey(key)] = value.str();
    return true;
}

inline std::set<ZhangHandoffFloatIdentity> zhangHandoffFloatDependencies(
    const KFState& state)
{
    std::set<ZhangHandoffFloatIdentity> result;
    for (const auto& [name, value] : state.metaDataMap)
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX))
        {
            std::istringstream input(value);
            input.imbue(std::locale::classic());
            ZhangHandoffFloatIdentity identity;
            if (zhangHandoffFloatParseIdentity(input, identity))
                result.insert(std::move(identity));
        }
    return result;
}

inline bool zhangHandoffFloatFind(
    const KFState& state, const ZhangHandoffFloatIdentity& identity, KFKey& key)
{
    for (const auto& [candidate, index] : state.kfIndexMap)
    {
        if (candidate.type != KF::PRODUCT_HANDOFF) continue;
        ZhangHandoffFloatRecord record;
        if (zhangHandoffFloatReadRecord(state, candidate, record) &&
            zhangHandoffFloatRecordCurrent(state, record) &&
            !(identity < record.identity) && !(record.identity < identity))
        {
            key = candidate;
            return true;
        }
    }
    return false;
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
    std::set<std::string> expected;
    std::set<ZhangHandoffFloatIdentity> identities;
    if (state.x.rows() != state.P.rows() || state.P.rows() != state.P.cols() ||
        state.dx.rows() != state.x.rows() ||
        state.kfIndexMap.size() != static_cast<std::size_t>(state.x.rows()))
    { failureReason = "HANDOFF_FLOAT_MOMENT_DIMENSION_MISMATCH"; return false; }
    if (!std::isfinite(state.time.bigTime))
    { failureReason = "HANDOFF_FLOAT_EPOCH_INVALID"; return false; }
    for (const auto& key : removedKeys)
        if (key.type != KF::PRODUCT_HANDOFF || !state.kfIndexMap.contains(key))
        { failureReason = "HANDOFF_FLOAT_RETIREMENT_SELECTION_INVALID"; return false; }
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
        ZhangHandoffFloatRecord record;
        if (!zhangHandoffFloatReadRecord(metadata, key, record) ||
            !identities.insert(record.identity).second ||
            !state.stateTransitionMap.contains(key))
        { failureReason = "HANDOFF_FLOAT_IDENTITY_OR_MOMENT_INVALID"; return false; }
        const long double age = state.time.bigTime - record.birthTime;
        if (!std::isfinite(age) || age < 0)
        { failureReason = "HANDOFF_FLOAT_FUTURE_OR_NONFINITE_BIRTH"; return false; }
        if (!allowRetirableExpired && age > retention)
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
    std::size_t dependencies = 0;
    for (const auto& [name, value] : metadata)
    {
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_RECORD_PREFIX) &&
            !expected.contains(name))
        { failureReason = "HANDOFF_FLOAT_ORPHAN_REGISTRY_ENTRY"; return false; }
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX))
        {
            ++dependencies;
            std::istringstream input(value);
            input.imbue(std::locale::classic());
            ZhangHandoffFloatIdentity identity;
            if (!zhangHandoffFloatParseIdentity(input, identity))
            { failureReason = "HANDOFF_FLOAT_INVALID_DEPENDENCY"; return false; }
            input >> std::ws;
            if (!input.eof())
            { failureReason = "HANDOFF_FLOAT_INVALID_DEPENDENCY"; return false; }
        }
    }
    if (expected.size() > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES ||
        dependencies > ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES)
    { failureReason = "HANDOFF_FLOAT_RETENTION_CAPACITY_EXCEEDED"; return false; }
    failureReason = "BOUNDED_CORRELATED_FLOAT_REGISTRY_VALID";
    return true;
}

inline bool zhangHandoffFloatRegistryValid(
    const KFState& state, std::string& failureReason)
{
    return zhangHandoffFloatRegistryViewValid(
        state, state.metaDataMap, {}, false, failureReason);
}
