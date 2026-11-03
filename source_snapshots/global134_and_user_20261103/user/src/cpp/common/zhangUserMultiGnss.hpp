#pragma once

#include <utility>
#include "common/enums.h"

// This experiment has exactly one dated dual-frequency product convention
// per constellation. Do not infer Galileo/BDS observables from GPS labels.
inline bool zhangUserGecSystem(E_Sys system)
{
    return system == E_Sys::GPS || system == E_Sys::GAL || system == E_Sys::BDS;
}

inline std::pair<E_ObsCode, E_ObsCode> zhangUserGecCodes(E_Sys system)
{
    if (system == E_Sys::GPS) return {E_ObsCode::L1C, E_ObsCode::L2W};
    if (system == E_Sys::GAL) return {E_ObsCode::L1C, E_ObsCode::L5Q};
    if (system == E_Sys::BDS) return {E_ObsCode::L1P, E_ObsCode::L5P};
    return {E_ObsCode::NONE, E_ObsCode::NONE};
}

inline bool zhangUserGecPair(E_Sys system, E_ObsCode first, E_ObsCode second)
{
    const auto expected = zhangUserGecCodes(system);
    return zhangUserGecSystem(system) && first == expected.first && second == expected.second;
}

inline E_Sys zhangUserClockSystem(bool internalUser, int activeSystems, E_Sys system)
{
    // A common clock would constrain the independent network code datums.
    // Preserve the legacy single-system state identity for old experiments.
    return internalUser && activeSystems > 1 ? system : E_Sys::NONE;
}
