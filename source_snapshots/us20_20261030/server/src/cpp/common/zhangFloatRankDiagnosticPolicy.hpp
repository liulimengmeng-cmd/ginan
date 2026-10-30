#pragma once

#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <string_view>

// This controls only the added pure-observation rank report. It must never
// control state initialisation, LSQ, innovation factorisation, KF or QC.
// Normal GINAN estimation does not need this extra measurement copy/SVD.
namespace ZhangFloatRankDiagnostic
{
inline constexpr const char* environmentVariable = "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC";
inline constexpr std::size_t boundedMaximumRows = 1024;
inline constexpr std::size_t boundedMaximumColumns = 1024;

enum class Mode { Off, Bounded, Full };

inline Mode parseMode(std::string_view value)
{
    if (value.empty() || value == "off") return Mode::Off;
    if (value == "bounded") return Mode::Bounded;
    if (value == "full") return Mode::Full;
    throw std::invalid_argument(
        "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC must be off, bounded, or full; "
        "full explicitly opts into the legacy unbounded diagnostic SVD");
}

inline Mode modeFromEnvironment()
{
    const char* value = std::getenv(environmentVariable);
    return parseMode(value ? std::string_view(value) : std::string_view{});
}

inline const char* modeName(Mode mode)
{
    switch (mode)
    {
        case Mode::Off: return "off";
        case Mode::Bounded: return "bounded";
        case Mode::Full: return "full";
    }
    return "invalid";
}

inline bool requested(Mode mode, bool parentDiagnosticsEnabled)
{
    return parentDiagnosticsEnabled && mode != Mode::Off;
}

inline bool allowsRows(Mode mode, std::size_t rows)
{
    return mode == Mode::Full || (mode == Mode::Bounded && rows <= boundedMaximumRows);
}

inline bool allowsColumns(Mode mode, std::size_t columns)
{
    return mode == Mode::Full || (mode == Mode::Bounded && columns <= boundedMaximumColumns);
}
}  // namespace ZhangFloatRankDiagnostic
