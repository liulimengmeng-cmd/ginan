#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

// Detector evidence and observation eligibility are separate from a committed
// graph arc serial. Only the graph transaction may consume pendingBreak.
struct ZhangSignalTracking
{
    bool seen = false;
    bool halfCycle = false;
    bool pendingBreak = false;
    double lastPhaseTime = 0;
    int lastCode = 0;
    unsigned long long eventSequence = 0;
    template<class A> void serialize(A& a, const unsigned int&)
    { a & seen & halfCycle & pendingBreak & lastPhaseTime & lastCode & eventSequence; }
};

struct ZhangPhaseArcDecision
{
    bool breakArc = false;
    bool useCurrentPhase = true;
    bool useCurrentCode = true;
    bool halfCycleResolved = true;
    bool gap = false;
};

inline ZhangPhaseArcDecision zhangDecidePhaseTracking(
    ZhangSignalTracking& state, unsigned int lli, bool phasePresent,
    int code, double time, double maximumContinuousGap)
{
    ZhangPhaseArcDecision out;
    if (!phasePresent) { out.useCurrentPhase=false; return out; }
    const bool half = lli & 2;
    out.halfCycleResolved = !half;
    out.useCurrentPhase = !half;
    out.gap = state.seen && (time <= state.lastPhaseTime ||
        time-state.lastPhaseTime > maximumContinuousGap || code != state.lastCode);
    // A persistent half-cycle interval is one quarantine, not one retirement
    // per epoch. Recovery starts a fresh whole-cycle arc without guessing parity.
    out.breakArc = (half != state.halfCycle) || (!half && ((lli & 1) || out.gap));
    if (out.breakArc) { state.pendingBreak=true; ++state.eventSequence; }
    state.halfCycle=half; state.seen=true; state.lastPhaseTime=time; state.lastCode=code;
    return out;
}

struct ZhangAcceptedSample
{
    double time=0, value=0, variance=0;
    template<class A> void serialize(A& a, const unsigned int&)
    { a & time & value & variance; }
};

struct ZhangCombinationDetector
{
    std::vector<ZhangAcceptedSample> accepted;
    bool suspect=false;
    double suspectTime=0, suspectInnovation=0, suspectVariance=0;
    double lastRawValue=0;
    template<class A> void serialize(A& a, const unsigned int&)
    { a & accepted & suspect & suspectTime & suspectInnovation & suspectVariance & lastRawValue; }
};

struct ZhangCombinationDecision
{
    bool useSample=false, breakArc=false, suspect=false;
    double prediction=0, innovation=0, variance=0, statistic=0;
};

// Weighted local-linear prediction on accepted samples only. R is in the
// combination's own units. q is an explicit detector-model variance rate,
// never the main filter's STEC process noise. Independent temporal errors are
// an assumption of this baseline; its nominal tail probability is not a
// measured false-alarm rate.
inline ZhangCombinationDecision zhangInspectCombination(
    ZhangCombinationDetector& state, double time, double value, double variance,
    double varianceRate, double sigmaGate, double maximumGap, bool reset=false)
{
    ZhangCombinationDecision out;
    if (std::isfinite(value)) state.lastRawValue=value;
    if (!std::isfinite(time) || !std::isfinite(value) || !std::isfinite(variance) ||
        variance <= 0 || !std::isfinite(varianceRate) || varianceRate < 0 || !std::isfinite(sigmaGate) || sigmaGate <= 0 ||
        !std::isfinite(maximumGap) || maximumGap <= 0)
    { out.suspect=true; return out; }
    if (reset || (!state.accepted.empty() &&
        (time <= state.accepted.back().time || time-state.accepted.back().time > maximumGap)))
    { state.accepted.clear(); state.suspect=false; }
    if (state.accepted.empty())
    {
        state.accepted.push_back({time,value,variance});
        out.useSample=true;out.prediction=value;out.variance=variance;return out;
    }
    double sw=0, st=0, stt=0, sy=0, sty=0;
    for (const auto& v:state.accepted)
    {
        const double t=v.time-time, w=1/v.variance;
        sw+=w;st+=w*t;stt+=w*t*t;sy+=w*v.value;sty+=w*t*v.value;
    }
    const double det=sw*stt-st*st;
    double predictionVariance=1/sw;
    out.prediction=sy/sw;
    if (state.accepted.size()>=3 && det>1e-12*sw*stt)
    {
        out.prediction=(stt*sy-st*sty)/det;
        predictionVariance=stt/det;
    }
    out.innovation=value-out.prediction;
    out.variance=variance+predictionVariance+varianceRate*(time-state.accepted.back().time);
    out.statistic=out.innovation*out.innovation/out.variance;
    if (out.statistic > sigmaGate*sigmaGate)
    {
        const bool persistent=state.suspect && time>state.suspectTime &&
            time-state.suspectTime<=maximumGap &&
            out.innovation*state.suspectInnovation>0 &&
            std::abs(out.innovation-state.suspectInnovation)<=
                sigmaGate*std::sqrt(out.variance+state.suspectVariance);
        if (!persistent)
        {
            state.suspect=true;state.suspectTime=time;
            state.suspectInnovation=out.innovation;state.suspectVariance=out.variance;
            out.suspect=true;return out; // do not poison accepted history
        }
        out.breakArc=true;state.accepted.clear();
    }
    state.suspect=false;out.useSample=true;
    state.accepted.push_back({time,value,variance});
    if (state.accepted.size()>8) state.accepted.erase(state.accepted.begin());
    return out;
}
