#pragma once

#include <boost/json.hpp>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "common/receiver.hpp"
#include "common/zhangSystemClockDatum.hpp"

/** The Zhang network states are effective clocks in a deterministic zero
 * reference-receiver gauge. Their configured P/Q are effective-coordinate
 * priors, not independent physical-clock priors with an unmodelled random
 * reference. A seed translation changes no existing state or covariance.
 */
inline bool zhangUsesEffectiveNetworkClocks(const ACSConfig& config)
{
    return config.zhangFullRank.enable && !config.zhangPppAr.user_adapter;
}

inline std::vector<E_Source> zhangPhysicalClockSources(
    const std::vector<E_Source>& sources)
{
    std::vector<E_Source> physical;
    for (const auto source : sources)
        if (source != E_Source::KALMAN && source != E_Source::REMOTE)
            physical.push_back(source);
    return physical;
}

inline std::string zhangClockLifecycleMetadataKey(
    const char* kind, E_Sys system, const std::string& receiver = {})
{
    return std::string("ZHANG_CLOCK_LIFECYCLE_V2|") + kind + "|" +
        enum_to_string(system) + "|" + receiver;
}

inline bool zhangReadClockMetadata(
    const KFState& state, const std::string& key, boost::json::object& result)
{
    const auto found = state.metaDataMap.find(key);
    if (found == state.metaDataMap.end()) return false;
    boost::system::error_code error;
    const auto value = boost::json::parse(found->second, error);
    if (error || !value.is_object()) return false;
    result = value.as_object();
    return true;
}

inline std::string zhangSppPhysicalClockChart(const Solution& solution, E_Sys system)
{
    boost::json::object metadata;
    if (!zhangReadClockMetadata(solution.sppState,
            zhangClockLifecycleMetadataKey("SPP_CHART", system), metadata))
        return {};
    const auto* epoch = metadata.if_contains("epoch");
    const auto* chart = metadata.if_contains("chart");
    if (!epoch || !epoch->is_string() || epoch->as_string() != tsync.to_string(9) ||
        !chart || !chart->is_string()) return {};
    return std::string(chart->as_string());
}

inline bool zhangEffectiveClockSeed(
    const KFState& state, E_Sys system, double& translation, std::string* anchor = nullptr)
{
    boost::json::object metadata;
    if (!zhangReadClockMetadata(state,
            zhangClockLifecycleMetadataKey("SEED", system), metadata)) return false;
    const auto* epoch = metadata.if_contains("epoch");
    const auto* value = metadata.if_contains("translation_m");
    if (!epoch || !epoch->is_string() || epoch->as_string() != tsync.to_string(9) ||
        !value || !value->is_number()) return false;
    translation = value->to_number<double>();
    if (!std::isfinite(translation)) return false;
    if (anchor)
    {
        const auto* identity = metadata.if_contains("anchor");
        *anchor = identity && identity->is_string() ? std::string(identity->as_string()) : "";
    }
    return true;
}

inline bool zhangSppClockForSystem(const Receiver& receiver, E_Sys system, double& clock)
{
    if (receiver.sol.status != E_Solution::SINGLE) return false;
    const KFKey key{KF::REC_SYS_BIAS, SatSys(system, 0), receiver.id};
    return receiver.sol.sppState.getKFValue(key, clock) != E_Source::NONE &&
        std::isfinite(clock) && !zhangSppPhysicalClockChart(receiver.sol, system).empty();
}

inline bool zhangSavedSppClockOffset(
    const KFState& state, const Receiver& receiver, E_Sys system, double& offset)
{
    boost::json::object metadata;
    if (!zhangReadClockMetadata(state,
            zhangClockLifecycleMetadataKey("SPP_OFFSET", system, receiver.id), metadata))
        return false;
    const auto* epoch = metadata.if_contains("epoch");
    const auto* chart = metadata.if_contains("chart");
    const auto* value = metadata.if_contains("offset_m");
    if (!epoch || !epoch->is_string() || epoch->as_string() != state.time.to_string(9) ||
        !chart || !chart->is_string() ||
        chart->as_string() != zhangSppPhysicalClockChart(receiver.sol, system) ||
        chart->as_string().empty() || !value || !value->is_number()) return false;
    offset = value->to_number<double>();
    return std::isfinite(offset);
}

/** Serial before measurement workers. Entries are epoch-bound, numerical
 * initial guesses only. Existing means/P/Q are deliberately untouched.
 */
void prepareZhangEffectiveClockSeeds(
    Trace& trace, ReceiverMap& receivers, KFState& state, KFState& remoteState);
