#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

// This is an independent receiver QC route. It never enables the network graph
// and never calls the legacy SCDIA quadratic-form/statistical test.
inline bool zhangUserCausalQcRequested()
{
    const char* value = std::getenv("ZHANG_USER_QC");
    if (!value || std::string(value) == "LEGACY") return false;
    if (std::string(value) == "CAUSAL_TRACKING") return true;
    throw std::runtime_error("INVALID_ZHANG_USER_QC: use LEGACY or CAUSAL_TRACKING");
}

inline std::string zhangUserQcAppliedSequenceKey(
    const std::string& receiver, const std::string& satellite, int code)
{
    return "USER_AR_QC_APPLIED_SEQUENCE:" + receiver + ":" + satellite + ":" + std::to_string(code);
}

enum class ZhangUserQcState { ACCEPTED, SUSPECT, CONFIRMED_BREAK };

inline const char* zhangUserQcStateName(ZhangUserQcState state)
{
    switch (state)
    {
        case ZhangUserQcState::ACCEPTED: return "ACCEPTED";
        case ZhangUserQcState::SUSPECT: return "SUSPECT";
        case ZhangUserQcState::CONFIRMED_BREAK: return "CONFIRMED_BREAK";
    }
    return "SUSPECT";
}

struct ZhangUserTrendSample
{
    double time = 0, value = 0, sigma = 0;
};

struct ZhangUserQcTrend
{
    std::vector<ZhangUserTrendSample> accepted;
    unsigned consecutiveSuspects = 0;
    bool rawSeen = false;
    double lastRawTime = 0;
    double firstSuspectThreshold = 0;
};

struct ZhangUserTrendDecision
{
    bool accepted = false;
    bool warmup = false;
    bool gap = false;
    double prediction = 0;
    double innovation = 0;
    double threshold = 0;
    const char* reason = "INVALID_SAMPLE";
};

// A causal local quadratic trend, with actual per-sample standard deviations.
// No squared-normalised innovation or tail probability is calculated. Gates
// are engineering bounds, not a false-alarm certification. Bad samples never
// enter accepted history. A persistent jump remains unresolved, rather than
// being relabelled as a physical cycle slip merely because it persists.
inline ZhangUserTrendDecision zhangInspectUserTrend(
    ZhangUserQcTrend& state, double time, double value, double sigma,
    double minimumGate, double sigmaGate, double varianceRate,
    double maximumGap, bool reset = false)
{
    ZhangUserTrendDecision out;
    if (!std::isfinite(time) || !std::isfinite(value) || !std::isfinite(sigma) ||
        sigma <= 0 || !std::isfinite(minimumGate) || minimumGate <= 0 ||
        !std::isfinite(sigmaGate) || sigmaGate <= 0 ||
        !std::isfinite(varianceRate) || varianceRate < 0 ||
        !std::isfinite(maximumGap) || maximumGap <= 0) return out;
    out.gap = state.rawSeen &&
        (time <= state.lastRawTime || time - state.lastRawTime > maximumGap);
    if (reset || out.gap)
    {
        state.accepted.clear();
        state.consecutiveSuspects = 0;
        state.firstSuspectThreshold = 0;
    }
    state.rawSeen = true;
    state.lastRawTime = time;
    if (state.accepted.empty())
    {
        state.accepted.push_back({time, value, sigma});
        out.accepted = out.warmup = true;
        out.prediction = value;
        out.reason = "TREND_WARMUP";
        return out;
    }
    const auto& last = state.accepted.back();
    const double dt = time - last.time;
    // Missing *accepted* samples are not a physical gap. Never erase a known
    // unresolved innovation and silently learn its shifted value after timeout.
    if (state.consecutiveSuspects && dt > maximumGap)
    {
        out.reason = "STALE_ACCEPTED_HISTORY_UNRESOLVED";
        return out;
    }
    out.prediction = last.value;
    double predictedNoise = last.sigma;
    if (state.accepted.size() >= 2)
    {
        // Work in units of the actual last observation separation. Scaling is
        // important when GPST seconds are large and when sampling is irregular.
        const double scale = last.time - state.accepted[state.accepted.size()-2].time;
        if (!(scale > 0)) return out;
        const int degree = state.accepted.size() >= 5 ? 2 : 1;
        double normal[3][7] = {};
        for (const auto& sample : state.accepted)
        {
            const double x = (sample.time - time) / scale;
            const double basis[3] = {1, x, x*x};
            const double weight = 1 / sample.sigma / sample.sigma;
            for (int i=0; i<=degree; ++i)
            {
                for (int j=0; j<=degree; ++j) normal[i][j] += weight*basis[i]*basis[j];
                normal[i][3] += weight*basis[i]*sample.value;
            }
        }
        for (int i=0; i<=degree; ++i) normal[i][4+i] = 1;
        bool valid = true;
        for (int i=0; i<=degree; ++i)
        {
            int pivot = i;
            for (int j=i+1; j<=degree; ++j)
                if (std::abs(normal[j][i]) > std::abs(normal[pivot][i])) pivot = j;
            for (int k=0; k<7; ++k) std::swap(normal[pivot][k], normal[i][k]);
            if (!std::isfinite(normal[i][i]) || std::abs(normal[i][i]) < 1e-18)
            { valid = false; break; }
            const double divisor = normal[i][i];
            for (int k=i; k<7; ++k) normal[i][k] /= divisor;
            for (int j=0; j<=degree; ++j) if (j != i)
            {
                const double multiple = normal[j][i];
                for (int k=i; k<7; ++k) normal[j][k] -= multiple*normal[i][k];
            }
        }
        if (!valid) return out;
        out.prediction = normal[0][3];
        // Formal trend extrapolation noise under independent raw sample errors.
        // Model discrepancy is budgeted separately; this is not a false-alarm
        // probability and is not a squared-normalised innovation calculation.
        if (!(normal[0][4] > 0)) return out;
        predictedNoise = std::sqrt(normal[0][4]);
    }
    out.warmup = state.accepted.size() < 3;
    out.innovation = value - out.prediction;
    const double noiseAllowance = sigmaGate * std::hypot(sigma, predictedNoise);
    const double modelAllowance = std::sqrt(varianceRate * dt);
    // The second sample cannot distinguish initial ionospheric slope from a
    // jump. Bound acquisition conservatively; subsequent samples use trend.
    out.threshold = std::max(minimumGate * (state.accepted.size()==1 ? 4 : 1),
                             noiseAllowance + modelAllowance);
    // Uncertainty growth during quarantine must not turn the same unresolved
    // level shift into an accepted sample without new continuity evidence.
    if (state.consecutiveSuspects && state.firstSuspectThreshold > 0)
        out.threshold = std::min(out.threshold, state.firstSuspectThreshold);
    if (!std::isfinite(out.prediction) || !std::isfinite(out.threshold)) return out;
    if (std::abs(out.innovation) > out.threshold)
    {
        if (state.consecutiveSuspects == 0) state.firstSuspectThreshold = out.threshold;
        if (state.consecutiveSuspects < 1000000) ++state.consecutiveSuspects;
        out.reason = state.consecutiveSuspects > 1 ? "PERSISTENT_UNRESOLVED_INNOVATION"
                                                  : "UNRESOLVED_INNOVATION";
        return out;
    }
    state.consecutiveSuspects = 0;
    state.firstSuspectThreshold = 0;
    state.accepted.push_back({time, value, sigma});
    if (state.accepted.size() > 8) state.accepted.erase(state.accepted.begin());
    out.accepted = true;
    out.reason = out.warmup ? "TREND_WARMUP" : "TREND_ACCEPTED";
    return out;
}

struct ZhangUserQcEligibility
{
    bool usePhase = true, useCode = true;
    ZhangUserQcState state = ZhangUserQcState::ACCEPTED;
    const char* reason = "OBSERVATION_ACCEPTED";
};

inline ZhangUserQcEligibility zhangUserCombinationEligibility(
    bool physicalBreak, bool halfCycle, bool gfAvailable, bool gfAccepted,
    bool mwAvailable, bool mwAccepted, bool independentCodeEvidenceBad = false)
{
    ZhangUserQcEligibility out;
    if (physicalBreak)
    {
        out.state = ZhangUserQcState::CONFIRMED_BREAK;
        out.reason = "TRACKING_DISCONTINUITY";
    }
    if (halfCycle)
    {
        out.usePhase = false;
        out.state = physicalBreak ? ZhangUserQcState::CONFIRMED_BREAK : ZhangUserQcState::SUSPECT;
        out.reason = "HALF_CYCLE_UNRESOLVED";
        return out;
    }
    if (mwAvailable && !mwAccepted)
    {
        // MW contains code and phase. Without independent code evidence it
        // cannot justify deleting good pseudorange. Quarantine the phase
        // continuity claim; no physical integer is retired.
        out.usePhase = false;
        if (!physicalBreak) out.state = ZhangUserQcState::SUSPECT;
        out.reason = "MW_SOURCE_UNIDENTIFIED";
    }
    else if (gfAvailable && !gfAccepted)
    {
        out.usePhase = false; // GF has no code: preserve independently good code.
        if (!physicalBreak) out.state = ZhangUserQcState::SUSPECT;
        out.reason = "GF_PHASE_INNOVATION_UNRESOLVED";
    }
    if (independentCodeEvidenceBad)
    {
        out.useCode = false;
        if (!physicalBreak) out.state = ZhangUserQcState::SUSPECT;
        out.reason = "INDEPENDENT_CODE_QC_EVIDENCE";
    }
    return out;
}
