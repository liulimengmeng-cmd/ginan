#pragma once

#include <boost/json.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "common/receiver.hpp"
#include "common/constants.hpp"
#include "common/zhangSystemClockDatum.hpp"

/** Diagnostic only. No state/offset/source choice is changed by this switch. */
inline bool zhangClockLifecycleDiagnosticsEnabled()
{
    const char* value = std::getenv("ZHANG_CLOCK_LIFECYCLE_DIAGNOSTICS");
    return value && std::string(value) == "1";
}

inline boost::json::value zhangClockLifecycleNumber(double value)
{
    return std::isfinite(value) ? boost::json::value(value) : boost::json::value(nullptr);
}

/** Machine-readable physical clock boundaries, also available without a
 * receiver TRACE. KF clock, SPP coordinate and observation source are separate
 * fields: matching constellation names do not prove matching clock charts.
 * No diagnostic metadata is written into the authoritative KF or checkpoint.
 */
inline void zhangClockLifecycleRecord(
    const char* phase, Receiver& receiver, const KFState* state,
    boost::json::object details = {})
{
    if (!zhangClockLifecycleDiagnosticsEnabled()) return;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    auto& solution = receiver.sol;
    boost::json::object record;
    record["schema"] = "ZHANG_CLOCK_LIFECYCLE_V1";
    record["phase"] = phase;
    record["epoch"] = tsync.to_string(0);
    record["kf_epoch"] = state ? boost::json::value(state->time.to_string(0))
                               : boost::json::value(nullptr);
    record["receiver"] = receiver.id;
    record["spp_status"] = enum_to_string(solution.status);
    record["spp_system"] = enum_to_string(solution.clkRefSys);
    record["spp_clock_m"] = zhangClockLifecycleNumber(solution.sppClk);
    record["saved_offset_m"] = zhangClockLifecycleNumber(solution.sppPppClkOffset);
    record["saved_offset_system"] = enum_to_string(solution.sppPppClkOffsetSys);
    record["adjust_ready"] = solution.clkAdjustReady;
    record["scoped_keys"] = zhangUsesSystemReceiverClocks(acsConfig);
    boost::json::array clocks;
    for (const auto& [system, enabled] : acsConfig.process_sys)
    {
        if (!enabled) continue;
        auto key = zhangReceiverClockKey(acsConfig, receiver.id, system);
        double clock = missing;
        double rate = missing;
        E_Source found = E_Source::NONE;
        if (state) found = state->getKFValue(key, clock);
        key.type = KF::REC_CLOCK_RATE;
        if (state) state->getKFValue(key, rate);
        boost::json::object row;
        row["system"] = enum_to_string(system);
        row["key_sat"] = key.Sat.id();
        row["kf_clock_m"] = zhangClockLifecycleNumber(clock);
        row["kf_rate_mps"] = zhangClockLifecycleNumber(rate);
        row["kf_found_source"] = enum_to_string(found);
        auto options = acsConfig.zhangFullRank.sysOpts.find(system);
        row["configured_reference"] = options == acsConfig.zhangFullRank.sysOpts.end()
            ? boost::json::value(nullptr) : boost::json::value(options->second.reference_receiver);
        clocks.push_back(std::move(row));
    }
    record["model_clocks"] = std::move(clocks);
    boost::json::array observations;
    for (const auto& pointer : receiver.obsList)
    {
        auto* observation = dynamic_cast<GObs*>(pointer.get());
        if (!observation) continue;
        const auto enabled = acsConfig.process_sys.find(observation->Sat.sys);
        if (enabled == acsConfig.process_sys.end() || !enabled->second) continue;
        boost::json::object row;
        row["satellite"] = observation->Sat.id();
        row["system"] = enum_to_string(observation->Sat.sys);
        row["source"] = enum_to_string(observation->clkSource);
        row["clock_m"] = zhangClockLifecycleNumber(observation->satClk * CLIGHT);
        row["eph_clock_valid"] = observation->ephClkValid;
        row["excluded"] = observation->exclude != 0;
        row["spp_valid"] = observation->sppValid;
        double clock = missing;
        if (state)
        {
            KFKey key;
            key.type = KF::SAT_CLOCK;
            key.Sat = observation->Sat;
            state->getKFValue(key, clock);
        }
        row["kf_clock_m"] = zhangClockLifecycleNumber(clock);
        observations.push_back(std::move(row));
    }
    record["observations"] = std::move(observations);
    for (auto& item : details) record[item.key()] = std::move(item.value());
    BOOST_LOG_TRIVIAL(info) << "ZHANG_CLOCK_LIFECYCLE " << boost::json::serialize(record);
}
