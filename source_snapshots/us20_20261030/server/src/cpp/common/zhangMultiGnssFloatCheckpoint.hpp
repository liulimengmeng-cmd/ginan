#pragma once

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "common/acsConfig.hpp"
#include "common/zhangR51Warmup.hpp"
#include "common/zhangGec134Migration.hpp"
#include "common/zhangGec134ArContinuation.hpp"

// Deliberately separate from the legacy GPS/AR checkpoint contract. This
// opt-in authorises only cold FLOAT schedules, never multi-GNSS integer state.
inline bool zhangMultiGnssFloatCheckpointActive(const ACSConfig& config)
{
    const auto* setting = std::getenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT");
    if (!setting || std::strcmp(setting, "1") != 0) return false;
    for (const auto& [system, enabled] : config.process_sys)
        if (enabled && system != E_Sys::GPS) return true;
    return false;
}

inline bool zhangCheckpointFlagEqualsOne(const char* name)
{
    const auto* setting = std::getenv(name);
    return setting && std::strcmp(setting, "1") == 0;
}

// A separate, versioned capability; never reinterpret a FLOAT archive as an
// AR archive merely by moving the AR start time inside its schedule.
inline bool zhangMultiGnssArCheckpointActive(const ACSConfig& config)
{
    if (!zhangCheckpointFlagEqualsOne("ZHANG_MULTI_GNSS_AR_CHECKPOINT") ||
        zhangCheckpointFlagEqualsOne("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT")) return false;
    std::size_t active = 0;
    for (const auto& [system, enabled] : config.process_sys)
        if (enabled && system != E_Sys::GPS && system != E_Sys::GAL && system != E_Sys::BDS)
            return false;
        else if (enabled) active++;
    if (active != 3) return false;
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        const auto found = config.process_sys.find(system);
        if (found == config.process_sys.end() || !found->second) return false;
    }
    return true;
}

inline std::string zhangCheckpointModeForConfig(
    const ACSConfig& config, std::string& failureReason)
{
    failureReason.clear();
    std::vector<E_Sys> systems;
    for (const auto& [system, enabled] : config.process_sys)
        if (enabled) systems.push_back(system);
    auto fail = [&](const char* reason) -> std::string
    { failureReason = reason; return {}; };
    if (systems.empty()) return fail("CHECKPOINT_NO_ACTIVE_SYSTEM");
    if (std::abs(config.epoch_interval - 30) > 1e-9 &&
        std::abs(config.epoch_interval - 60) > 1e-9)
        return fail("CHECKPOINT_REQUIRES_30_OR_60_SECOND_INTERVAL");
    const bool legacyGps = systems.size() == 1 && systems.front() == E_Sys::GPS;
    const bool arMode = zhangMultiGnssArCheckpointActive(config);
    if (zhangCheckpointFlagEqualsOne("ZHANG_MULTI_GNSS_AR_CHECKPOINT") &&
        zhangCheckpointFlagEqualsOne("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT"))
        return fail("MULTI_GNSS_CHECKPOINT_FLOAT_AR_MODE_CONFLICT");
    if (zhangCheckpointFlagEqualsOne("ZHANG_MULTI_GNSS_AR_CHECKPOINT") && !arMode)
        return fail("GEC134_AR_CHECKPOINT_REQUIRES_ALL_THREE_SYSTEMS");
    if (arMode && zhangMultiGnssFloatCheckpointActive(config))
        return fail("MULTI_GNSS_CHECKPOINT_FLOAT_AR_MODE_CONFLICT");
    if (!legacyGps && !zhangMultiGnssFloatCheckpointActive(config) && !arMode)
        return fail("MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_EXPLICIT_OPT_IN");
    for (const auto system : systems)
    {
        std::vector<E_ObsCode> expected;
        switch (system)
        {
            case E_Sys::GPS: expected = {E_ObsCode::L1C, E_ObsCode::L2W}; break;
            case E_Sys::GAL: expected = {E_ObsCode::L1C, E_ObsCode::L5Q}; break;
            case E_Sys::BDS: expected = {E_ObsCode::L1P, E_ObsCode::L5P}; break;
            default: return fail("FLOAT_CHECKPOINT_UNSUPPORTED_SYSTEM");
        }
        const auto options = config.zhangFullRank.sysOpts.find(system);
        if (options == config.zhangFullRank.sysOpts.end() ||
            options->second.baseline_observables != expected)
            return fail("CHECKPOINT_ORDERED_BASELINES_MISMATCH");
        if (!legacyGps)
        {
            const auto priorities = config.code_priorities.find(system);
            const auto products = config.zhangPppAr.baseline_observables.find(system);
            if (priorities == config.code_priorities.end() || priorities->second != expected ||
                products == config.zhangPppAr.baseline_observables.end() || products->second != expected)
                return fail("FLOAT_CHECKPOINT_EXACT_SIGNAL_PRIORITIES_MISMATCH");
        }
    }
    if (legacyGps) return "E29_GPS_L1C_L2W_ZHANG_FULL_RANK";
    if (!std::isfinite(config.epoch_interval))
        return fail("CHECKPOINT_REQUIRES_30_OR_60_SECOND_INTERVAL");
    if (!config.zhangFullRank.enable || !config.zhangPppAr.output_products ||
        config.zhangPppAr.user_adapter ||
        !config.zhangPppAr.transactional_integer_fixing ||
        config.zhangPppAr.product_relation_feedback ||
        !zhangCheckpointFlagEqualsOne("ZHANG_R51_ENABLE") ||
        !zhangCheckpointFlagEqualsOne("ZHANG_R51_RATIO_ONLY") ||
        !zhangCheckpointFlagEqualsOne("ZHANG_R51_EPOCH_ONLY_AR"))
        return fail("MULTI_GNSS_CHECKPOINT_REQUIRES_ISOLATED_FLOAT_AUTHORITY");
    if (config.start_epoch.is_special() ||
        config.end_epoch.is_special() ||
        config.end_epoch < config.start_epoch)
        return fail("MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_BOUNDED_EPOCHS");
    long double arStart = 0;
    try { arStart = zhangR51ParseArStart(std::getenv("ZHANG_R51_AR_START_GPST_SECONDS")); }
    catch (const std::exception&)
    { return fail("MULTI_GNSS_FLOAT_CHECKPOINT_INVALID_AR_SCHEDULE"); }
    // GTime::bigTime / experiment epoch_seconds use elapsed GPST from
    // 1980-01-06, not Unix seconds and not UTC with a leap-second correction.
    const boost::posix_time::ptime origin(boost::gregorian::date(1980, 1, 6));
    const long double end = (config.end_epoch - origin).total_microseconds() / 1000000.0L;
    if (arMode)
    {
        const bool continuation = zhangCheckpointFlagEqualsOne(ZHANG_GEC134_AR_CONTINUATION_FLAG);
        if (systems.size() != 3 || !zhangCheckpointFlagEqualsOne("ZHANG_GEC134_FLOAT_TO_AR_MIGRATION") ||
            std::getenv("ZHANG_R51_RATIO_THRESHOLD") == nullptr ||
            std::strcmp(std::getenv("ZHANG_R51_RATIO_THRESHOLD"), "2") != 0 ||
            arStart != 1405213230.0L || end < arStart ||
            end > (continuation ? 1405216830.0L : 1405216800.0L) ||
            config.epoch_interval != 30 ||
            config.zhangPppAr.checkpoint_runtime_id != ZHANG_GEC134_SOURCE_RUNTIME_ID)
            return fail("GEC134_AR_CHECKPOINT_BOUNDED_POLICY_REJECTED");
        // Do not widen the original server/FLOAT-migration classifier. Only
        // this explicit local continuation profile admits the extra 30 s;
        // its source/target authority is checked before any restore commit.
        const auto* continuationProfile = continuation ?
            zhangGec134ArContinuationProfileForOutputs(config.outputs_root) : nullptr;
        if (continuation &&
            (!continuationProfile || !zhangCheckpointFlagEqualsOne(ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG) ||
             end != continuationProfile->endBigTime || config.max_epochs != continuationProfile->maxEpochs ||
             config.zhangPppAr.checkpoint_restore_path != ZHANG_GEC134_AR_E122_PATH ||
             config.zhangPppAr.checkpoint_capture_epochs != continuationProfile->captures))
            return fail("GEC134_AR_CONTINUATION_EXACT_SCHEDULE_REQUIRED");
        return ZHANG_GEC134_AR_MODE;
    }
    if (!(arStart > end))
        return fail("MULTI_GNSS_FLOAT_CHECKPOINT_END_MUST_PRECEDE_AR");
    std::string mode = "E29_MULTI_GNSS_FLOAT_V2";
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        const auto active = config.process_sys.find(system);
        if (active != config.process_sys.end() && active->second)
            mode += "_" + enum_to_string(system);
    }
    return mode;
}
