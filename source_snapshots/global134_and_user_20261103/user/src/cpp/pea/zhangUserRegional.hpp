#pragma once

#include <string>
#include <vector>
#include "common/algebra.hpp"
#include "common/gTime.hpp"
#include "common/trace.hpp"

// The default is OFF. Unknown modes and incompatible enabled configurations
// fail closed; this switch does not change integer or HOLD authorisation.
bool zhangUserRegionalEnabled();

std::string exportZhangUs20RegionalCheckpoint(double completed);
void importZhangUs20RegionalCheckpoint(const std::string& payload,double completed,bool dryRun);

// Value-based metadata name, not an address/token that can differ across TUs.
// The marker follows each actual pseudo row through final QC row selection.
inline constexpr const char* ZHANG_USER_REGIONAL_SD_RESIDUAL_MARKER = "regionalSdIonoL1Metres";
inline constexpr const char* ZHANG_USER_REGIONAL_SD_METRE_LABEL = "REGIONAL_SD_L1_METRES:";

// Factors use ONE complete current-epoch CLOCK/PHASE/IONO covariance root.
// Code selects CLOCK, phase selects CLOCK-PHASE. The same noiseId and column
// numbering are also used by the regional ionosphere pseudo-observations.
bool zhangUserRegionalNoiseFactors(
    const GTime& time,
    const std::string& receiver,
    const SatSys& satellite,
    E_ObsCode observable,
    bool phase,
    std::vector<double>& factors,
    std::string& noiseId,
    int* rank = nullptr,
    std::string* reason = nullptr);

// Append at most one independent satellite-difference ionosphere row per
// non-reference satellite, using ONLY states already supported by constructed
// current physical observations. No absolute receiver GF datum is invented.
void zhangUserRegionalIonoPseudoObs(
    Trace& trace, KFState& state, KFMeasEntryList& entries);

inline constexpr const char* ZHANG_USER_REGIONAL_TROP_MARKER="regionalTropZtdMetres";
inline constexpr const char* ZHANG_USER_REGIONAL_TROP_LABEL="REGIONAL_TROP_ZTD_METRES:";
void zhangUserRegionalTropPseudoObs(Trace& trace, KFState& state, KFMeasEntryList& entries);
