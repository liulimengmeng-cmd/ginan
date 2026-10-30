#pragma once

#include <algorithm>
#include "common/acsConfig.hpp"

/** The configured whitelist is produced from dated metadata by the input
 * audit.  This predicate does not claim that a PRN number identifies a BDS
 * generation, and an empty requested whitelist is deliberately fail-closed. */
inline bool zhangMultiGnssSatelliteAllowed(
    const ACSConfig& config, const SatSys& satellite)
{
    if (satellite.sys != E_Sys::BDS || !config.zhangFullRank.bds3_only)
        return true;
    const auto& allowed = config.zhangFullRank.bds3_satellite_whitelist;
    return satellite.prn > 0 && std::find(allowed.begin(), allowed.end(),
        satellite.id()) != allowed.end();
}
