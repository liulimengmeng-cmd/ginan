#pragma once

#include <cmath>
#include <set>
#include <string>

#include "common/acsConfig.hpp"
#include "common/algebra.hpp"

/** Systems participating in the opt-in Zhang clock model.
 *
 * The decision is configuration based, not epoch-visibility based: losing all
 * observations of one constellation must not rename the surviving states.
 * Traditional Ginan and a single-system Zhang run retain their original NONE
 * receiver-clock key, including existing GPS-only checkpoints.
 */
inline std::set<E_Sys> zhangSystemReceiverClockSystems(const ACSConfig& config)
{
    std::set<E_Sys> systems;
    if (!config.zhangFullRank.enable && !config.zhangPppAr.user_adapter)
        return systems;
    for (const auto& [system, enabled] : config.process_sys)
    {
        if (!enabled) continue;
        if (config.zhangFullRank.enable)
        {
            const auto it = config.zhangFullRank.sysOpts.find(system);
            if (it != config.zhangFullRank.sysOpts.end() &&
                !it->second.baseline_observables.empty())
                systems.insert(system);
        }
        else
        {
            const auto it = config.zhangPppAr.baseline_observables.find(system);
            if (it != config.zhangPppAr.baseline_observables.end() &&
                !it->second.empty())
                systems.insert(system);
        }
    }
    return systems;
}

inline bool zhangUsesSystemReceiverClocks(const ACSConfig& config)
{
    if (!config.zhangFullRank.enable && !config.zhangPppAr.user_adapter)
        return false;
    int count = 0;
    for (const auto& [system, enabled] : config.process_sys)
    {
        if (!enabled) continue;
        bool configured = false;
        if (config.zhangFullRank.enable)
        {
            const auto it = config.zhangFullRank.sysOpts.find(system);
            configured = it != config.zhangFullRank.sysOpts.end() &&
                !it->second.baseline_observables.empty();
        }
        else
        {
            const auto it = config.zhangPppAr.baseline_observables.find(system);
            configured = it != config.zhangPppAr.baseline_observables.end() &&
                !it->second.empty();
        }
        if (configured && ++count > 1) return true;
    }
    return false;
}

/** num remains the original clock component (polynomial/drift component).
 * rec_ptr is metadata, but must survive construction for existing callbacks.
 */
inline KFKey zhangReceiverClockKey(
    const ACSConfig& config,
    const std::string& receiver,
    E_Sys system,
    int component = 0,
    KF type = KF::REC_CLOCK,
    Receiver* receiverPointer = nullptr)
{
    KFKey key;
    key.type = type;
    key.str = receiver;
    key.num = component;
    key.rec_ptr = receiverPointer;
    if (zhangUsesSystemReceiverClocks(config))
        key.Sat = SatSys(system, 0);
    return key;
}

/** A legacy NONE receiver key belongs only to the sole configured system.
 * In particular it is never implicitly imported as G/E/C in a joint run.
 */
inline bool zhangReceiverClockMatchesSystem(
    const ACSConfig& config, const KFKey& key, E_Sys system)
{
    if (key.type != KF::REC_CLOCK && key.type != KF::REC_CLOCK_RATE)
        return false;
    if (key.Sat.prn != 0) return false;
    const auto systems = zhangSystemReceiverClockSystems(config);
    if (systems.find(system) == systems.end()) return false;
    return systems.size() > 1
        ? key.Sat.sys == system
        : key.Sat.sys == E_Sys::NONE;
}

/** Reference-clock subtraction must never mix constellation datums. */
inline bool zhangClockDatumPivotApplies(
    const ACSConfig& config, const KFKey& target, const KFKey& pivot)
{
    if (!zhangUsesSystemReceiverClocks(config)) return true;
    return pivot.Sat.prn == 0 && pivot.Sat.sys == target.Sat.sys;
}

/** One square-root factor belongs to the complete epoch/solution Q, not to
 * a constellation. Different receivers and G/E/C measurements intentionally
 * share these unit-variance random variables. The caller must first validate
 * that the loaded covariance is complete; absent cross terms are not zeros.
 */
inline KFKey zhangProductJointNoiseKey(
    long long epoch, const std::string& solution, int factor)
{
    KFKey key;
    key.type = KF::SAT_CLOCK;
    key.str = "ZHANG_PRODUCT_COV_" + solution + "_" + std::to_string(epoch);
    key.num = factor;
    // Sat deliberately remains NONE: the columns are global Q factors.
    return key;
}
