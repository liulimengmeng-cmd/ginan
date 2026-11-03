#pragma once

#include "common/zhangProductDatumTransition.hpp"
#include <stdexcept>

// Candidate numerical representation only. This is not an AR authorization,
// a stochastic state, or evidence for transport across an unknown boundary.
// The production adapter must consume already validated endpoint contracts.
struct ZhangUserPredictionChartEndpoint
{
    ZhangProductDatumEndpoint anchor;
    ZhangProductDatumEndpoint applied;
    double wavelength = 0;
    long long cycles = 0;
    std::string lastProof;
};

inline bool zhangPredictionChartCyclesValid(long long value)
{
    return std::abs(static_cast<long double>(value)) <= 4503599627370496.L;
}

inline bool zhangPredictionChartEndpointValid(const ZhangUserPredictionChartEndpoint& x)
{
    return zhangDatumEndpointValid(x.anchor) && zhangDatumEndpointValid(x.applied)
        && x.anchor.epoch <= x.applied.epoch
        && x.anchor.satellite == x.applied.satellite
        && x.anchor.observable == x.applied.observable
        && x.wavelength > 0 && std::isfinite(x.wavelength)
        && zhangPredictionChartCyclesValid(x.cycles);
}

// Called explicitly for first observation of a product endpoint or for a new
// representation root after UNKNOWN. It asserts no connection to an old root.
inline ZhangUserPredictionChartEndpoint zhangPredictionChartNewRoot(
    const ZhangProductDatumEndpoint& endpoint, double wavelength)
{
    ZhangUserPredictionChartEndpoint out{endpoint, endpoint, wavelength, 0, {}};
    if (!zhangPredictionChartEndpointValid(out))
        throw std::runtime_error("PREDICTION_CHART_NEW_ROOT_INVALID");
    return out;
}

// The shift remains an exact integer. No division/rounding of metres and no
// integer permission is inferred from a current phase difference or FLOAT.
// Return false without modifying the destination on any unsupported map.
inline bool zhangPredictionChartTransport(
    const ZhangUserPredictionChartEndpoint& before,
    const ZhangProductDatumTransition& transition,
    bool validatedFloatContract,
    ZhangUserPredictionChartEndpoint& after)
{
    if (!validatedFloatContract || !zhangPredictionChartEndpointValid(before)
        || !zhangDatumTransitionPureIntegerValid(transition)
        || transition.wavelength != before.wavelength)
        return false;
    const auto proof = zhangEncodeDatumTransition(transition);
    if (zhangDatumEndpointIdentity(before.applied) == zhangDatumEndpointIdentity(transition.target))
    {
        if (before.lastProof != proof) return false;
        after = before;
        return true;
    }
    if (zhangDatumEndpointIdentity(before.applied) != zhangDatumEndpointIdentity(transition.source))
        return false;
    // Both summands are at most 2^52: addition is within int64, then enforce
    // the declared exactly representable range before committing.
    const auto sum = before.cycles + transition.correctionShiftCycles;
    if (!zhangPredictionChartCyclesValid(sum)) return false;
    auto candidate = before;
    candidate.applied = transition.target;
    candidate.cycles = sum;
    candidate.lastProof = proof;
    if (!zhangPredictionChartEndpointValid(candidate)) return false;
    after = std::move(candidate);
    return true;
}

inline std::string zhangPredictionChartEncode(const ZhangUserPredictionChartEndpoint& value)
{
    if (!zhangPredictionChartEndpointValid(value))
        throw std::runtime_error("PREDICTION_CHART_SERIALIZE_INVALID");
    std::ostringstream out;
    out << "US20_PREDICTION_ENDPOINT_V1 ";
    zhangWriteDatumEndpoint(out, value.anchor);
    out << ' ';
    zhangWriteDatumEndpoint(out, value.applied);
    out << ' ' << std::setprecision(17) << value.wavelength << ' '
        << value.cycles << ' ' << std::quoted(value.lastProof);
    return out.str();
}

inline bool zhangPredictionChartDecode(const std::string& text,
    ZhangUserPredictionChartEndpoint& destination)
{
    std::istringstream input(text);
    std::string schema;
    ZhangUserPredictionChartEndpoint candidate;
    if (!(input >> schema) || schema != "US20_PREDICTION_ENDPOINT_V1"
        || !zhangReadDatumEndpoint(input, candidate.anchor)
        || !zhangReadDatumEndpoint(input, candidate.applied)
        || !(input >> candidate.wavelength >> candidate.cycles >> std::quoted(candidate.lastProof)))
        return false;
    input >> std::ws;
    if (!input.eof() || !zhangPredictionChartEndpointValid(candidate)) return false;
    if (!candidate.lastProof.empty())
    {
        ZhangProductDatumTransition last;
        if (!zhangDecodeDatumTransition(candidate.lastProof, last)
            || !zhangDatumTransitionPureIntegerValid(last)
            || zhangDatumEndpointIdentity(last.target) != zhangDatumEndpointIdentity(candidate.applied)
            || last.wavelength != candidate.wavelength)
            return false;
    }
    else if (candidate.cycles != 0
        || zhangDatumEndpointIdentity(candidate.anchor) != zhangDatumEndpointIdentity(candidate.applied))
        return false;
    if (zhangPredictionChartEncode(candidate) != text) return false;
    destination = std::move(candidate);
    return true;
}

struct ZhangUserPredictionValues
{
    double correctionMetres = 0;
    double receiverPhaseMetres = 0;
    double ambiguityCycles = 0;
};

// One authoritative per-satellite cumulative product shift g. The receiver
// common shift is g_reference; relative numerical offset is g_s-g_reference.
// This remains independent of the AR offset, which may be retired locally.
// Callers must verify current endpoint ownership before asking for this view.
// False means leave all three original components untouched, not partially
// canonicalize a row. The KF state, derivatives, Q and noise IDs are untouched.
inline bool zhangPredictionChartView(
    double clock, double phase, double rawCorrection,
    double receiverPhase, double ambiguity, double wavelength,
    long long satelliteCycles, long long referenceCycles,
    ZhangUserPredictionValues& destination)
{
    if (satelliteCycles == 0 && referenceCycles == 0) return false;
    if (!zhangPredictionChartCyclesValid(satelliteCycles)
        || !zhangPredictionChartCyclesValid(referenceCycles)
        || !(wavelength > 0) || !std::isfinite(wavelength)
        || !std::isfinite(clock) || !std::isfinite(phase)
        || !std::isfinite(rawCorrection) || !std::isfinite(receiverPhase)
        || !std::isfinite(ambiguity)) return false;
    // An independently represented nominal correction must not be silently
    // replaced by clock-phase just because an older input gate allows a small
    // difference. Unsupported decompositions keep the complete legacy row.
    if (rawCorrection != clock - phase) return false;
    const auto relative = satelliteCycles - referenceCycles;
    if (!zhangPredictionChartCyclesValid(relative)) return false;
    ZhangUserPredictionValues candidate;
    candidate.correctionMetres = clock - (phase + wavelength * static_cast<double>(satelliteCycles));
    candidate.receiverPhaseMetres = receiverPhase - wavelength * static_cast<double>(referenceCycles);
    candidate.ambiguityCycles = ambiguity - static_cast<double>(relative);
    if (!std::isfinite(candidate.correctionMetres)
        || !std::isfinite(candidate.receiverPhaseMetres)
        || !std::isfinite(candidate.ambiguityCycles)) return false;
    destination = candidate;
    return true;
}
