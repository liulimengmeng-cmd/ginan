#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Candidate-only policy. All inputs are values: applying this gate cannot
// remove observations, mutate their validity, reset an ambiguity or condition
// a filter. Without explicit opt-in the historical candidate set is retained;
// when enabled the caller supplies the final resolved YAML threshold.
struct ZhangUserArElevationDecision
{
    bool enabled = false;
    bool eligible = true;
    bool geometryValid = false;
    bool observationTimeValid = false;
    const char* reason = "AR_ELEVATION_NOT_ENABLED";
};

inline ZhangUserArElevationDecision zhangUserArElevationDecision(
    bool explicitlyConfigured,
    double minimumDegrees,
    double elevationRadians,
    bool observationGeometryValid,
    double observationAgeSeconds = 0,
    double geometryToleranceSeconds = 0.5)
{
    constexpr double degreesToRadians = 3.14159265358979323846 / 180.0;
    ZhangUserArElevationDecision out;
    out.enabled = explicitlyConfigured;
    out.geometryValid = observationGeometryValid && std::isfinite(elevationRadians);
    out.observationTimeValid = std::isfinite(observationAgeSeconds) &&
        std::isfinite(geometryToleranceSeconds) && geometryToleranceSeconds >= 0 &&
        std::abs(observationAgeSeconds) <= geometryToleranceSeconds;
    if (!explicitlyConfigured) return out;
    out.eligible = false;
    if (!std::isfinite(minimumDegrees) || minimumDegrees < 0 || minimumDegrees > 90)
        out.reason = "INVALID_AR_ELEVATION_CONFIGURATION";
    else if (!out.geometryValid)
        out.reason = "INVALID_OBSERVATION_GEOMETRY";
    else if (!std::isfinite(geometryToleranceSeconds) || geometryToleranceSeconds < 0)
        out.reason = "INVALID_GEOMETRY_TIME_TOLERANCE";
    else if (!out.observationTimeValid)
        out.reason = "STALE_OBSERVATION_GEOMETRY";
    else if (elevationRadians < -90 * degreesToRadians ||
             elevationRadians > 90 * degreesToRadians)
        out.reason = "INVALID_ELEVATION_RANGE";
    else if (elevationRadians < minimumDegrees * degreesToRadians)
        out.reason = "BELOW_AR_ELEVATION_MASK";
    else
    {
        out.eligible = true;
        out.reason = "AR_ELEVATION_ACCEPTED";
    }
    return out;
}

// Experiment policy only. These starting thresholds are not calibrated
// reliability bounds, and repeated candidates are correlated observations.
struct ZhangUserArExperimentPolicy
{
    double minimumArcSeconds = 120;
    std::size_t minimumArcSamples = 5;
    double maximumGapSeconds = 45;
    std::size_t confirmationEpochs = 3;
    double confirmationIntervalSeconds = 30;
    double timeToleranceSeconds = 1e-6;
};

struct ZhangUserArExperimentKey
{
    std::string runtime;
    std::string receiver;
    // Satellite for the arc tracker; constellation/group for the first gate.
    std::string object;

    bool operator<(const ZhangUserArExperimentKey& other) const
    {
        return std::tie(runtime, receiver, object) <
            std::tie(other.runtime, other.receiver, other.object);
    }
    bool valid() const
    {
        return !runtime.empty() && !receiver.empty() && !object.empty();
    }
};

inline std::string zhangUserArEncodedField(const std::string& field)
{
    return std::to_string(field.size()) + ":" + field;
}

// The caller supplies actual physical-arc and product identities, not merely
// satellite names. Both signals must be available for a common dual-band arc.
inline std::string zhangUserArArcIdentity(
    const std::string& firstPhysical,
    const std::string& firstProduct,
    const std::string& secondPhysical,
    const std::string& secondProduct)
{
    if (firstPhysical.empty() || firstProduct.empty() ||
        secondPhysical.empty() || secondProduct.empty()) return {};
    return zhangUserArEncodedField(firstPhysical) +
        zhangUserArEncodedField(firstProduct) +
        zhangUserArEncodedField(secondPhysical) +
        zhangUserArEncodedField(secondProduct);
}

struct ZhangUserArExperimentArc
{
    bool valid = false;
    bool eligible = false;
    double first = 0;
    double last = 0;
    double age = 0;
    std::size_t count = 0;
    std::string identity;
    std::string reason = "NOT_OBSERVED";
};

struct ZhangUserArcTrackerCheckpoint
{
    ZhangUserArExperimentPolicy policy;
    std::map<ZhangUserArExperimentKey,ZhangUserArExperimentArc> arcs;
};

inline bool zhangUserArcRuntimeFamily(const std::string& value,const std::string& root)
{return value==root || value.rfind(root+"/branch/",0)==0;}

inline bool zhangUserArcPolicyEqual(const ZhangUserArExperimentPolicy& a,const ZhangUserArExperimentPolicy& b)
{
    return std::tie(a.minimumArcSeconds,a.minimumArcSamples,a.maximumGapSeconds,a.confirmationEpochs,
        a.confirmationIntervalSeconds,a.timeToleranceSeconds)==
        std::tie(b.minimumArcSeconds,b.minimumArcSamples,b.maximumGapSeconds,b.confirmationEpochs,
        b.confirmationIntervalSeconds,b.timeToleranceSeconds);
}

class ZhangUserArExperimentArcTracker
{
public:
    explicit ZhangUserArExperimentArcTracker(
        ZhangUserArExperimentPolicy policy = {}) : policy_(policy) {}

    ZhangUserArExperimentArc update(
        const ZhangUserArExperimentKey& key,
        double time,
        const std::string& identity,
        bool valid = true)
    {
        auto existing = arcs_.find(key);
        if (!key.valid() || !valid || identity.empty() || !std::isfinite(time))
        {
            if (existing != arcs_.end()) arcs_.erase(existing);
            ZhangUserArExperimentArc out;
            out.reason = "INVALID_OR_DISCONTINUOUS";
            return out;
        }
        if (existing != arcs_.end() &&
            time < existing->second.last - policy_.timeToleranceSeconds)
        {
            arcs_.erase(existing);
            ZhangUserArExperimentArc out;
            out.reason = "TIME_REVERSED";
            return out;
        }

        std::string resetReason;
        if (existing == arcs_.end()) resetReason = "NEW_ARC";
        else if (identity != existing->second.identity) resetReason = "IDENTITY_CHANGED";
        else if (time - existing->second.last > policy_.maximumGapSeconds)
            resetReason = "OBSERVATION_GAP";

        if (!resetReason.empty())
        {
            ZhangUserArExperimentArc fresh;
            fresh.valid = true;
            fresh.first = fresh.last = time;
            fresh.count = 1;
            fresh.identity = identity;
            fresh.reason = std::move(resetReason);
            arcs_[key] = std::move(fresh);
        }
        else
        {
            auto& arc = existing->second;
            if (time - arc.last > policy_.timeToleranceSeconds)
            {
                arc.last = time;
                ++arc.count;
                arc.reason = "CONTINUOUS_SAMPLE";
            }
            else arc.reason = "DUPLICATE_EPOCH";
        }
        auto& arc = arcs_.at(key);
        arc.age = arc.last - arc.first;
        arc.eligible = arc.valid &&
            arc.age + policy_.timeToleranceSeconds >= policy_.minimumArcSeconds &&
            arc.count >= policy_.minimumArcSamples;
        return arc;
    }

    void erase(const ZhangUserArExperimentKey& key) { arcs_.erase(key); }
    void clear() { arcs_.clear(); }
    void clearRuntime(const std::string& runtime)
    {
        for (auto it = arcs_.begin(); it != arcs_.end();)
            if (it->first.runtime == runtime) it = arcs_.erase(it);
            else ++it;
    }
    std::size_t size() const { return arcs_.size(); }
    ZhangUserArcTrackerCheckpoint checkpoint(const std::string& runtime) const
    {
        ZhangUserArcTrackerCheckpoint value;value.policy=policy_;
        for(const auto& [key,arc]:arcs_)if(zhangUserArcRuntimeFamily(key.runtime,runtime))value.arcs.emplace(key,arc);
        return value;
    }
    void restoreCheckpoint(const ZhangUserArcTrackerCheckpoint& value,const std::string& runtime,double completed)
    {
        if(runtime.empty() || !std::isfinite(completed) || value.arcs.size()>100000 ||
            !zhangUserArcPolicyEqual(policy_,value.policy))throw std::runtime_error("US20_AR_TRACKER_POLICY_OR_BOUNDARY");
        for(const auto& [key,arc]:value.arcs)
        {
            const bool knownReason=arc.reason=="NEW_ARC" || arc.reason=="IDENTITY_CHANGED" ||
                arc.reason=="OBSERVATION_GAP" || arc.reason=="CONTINUOUS_SAMPLE" || arc.reason=="DUPLICATE_EPOCH";
            if(!key.valid() || !zhangUserArcRuntimeFamily(key.runtime,runtime) || !arc.valid ||
                !std::isfinite(arc.first) || !std::isfinite(arc.last) || !std::isfinite(arc.age) ||
                arc.first>arc.last || arc.last>completed || arc.age!=arc.last-arc.first || !arc.count ||
                arc.identity.empty() || !knownReason || arc.eligible!=(arc.age+policy_.timeToleranceSeconds>=policy_.minimumArcSeconds &&
                    arc.count>=policy_.minimumArcSamples))throw std::runtime_error("US20_AR_TRACKER_INVALID_RECORD");
        }
        auto candidate=arcs_;
        for(auto it=candidate.begin();it!=candidate.end();)
            if(zhangUserArcRuntimeFamily(it->first.runtime,runtime))it=candidate.erase(it);else ++it;
        candidate.insert(value.arcs.begin(),value.arcs.end());arcs_.swap(candidate);
    }

private:
    ZhangUserArExperimentPolicy policy_;
    std::map<ZhangUserArExperimentKey, ZhangUserArExperimentArc> arcs_;
};

struct ZhangUserArSeparatedArcAges
{
    ZhangUserArExperimentArc physical;
    ZhangUserArExperimentArc state;
    ZhangUserArExperimentArc product;
    ZhangUserArExperimentArc effective;
};

// A phase quarantine is an interruption of usable evidence, not proof of a
// new physical ambiguity. Product datum/segment changes similarly reset the
// product/effective age but cannot rewrite the physical age. The caller must
// supply a true tracking/state generation for the physical identity.
class ZhangUserArSeparatedArcTracker
{
public:
    explicit ZhangUserArSeparatedArcTracker(ZhangUserArExperimentPolicy policy)
        : physical_(ageOnlyPolicy(policy)), state_(ageOnlyPolicy(policy)),
          product_(ageOnlyPolicy(policy)), effective_(policy) {}

    ZhangUserArSeparatedArcAges update(const ZhangUserArExperimentKey& key,
        double time, const std::string& physicalIdentity, bool physicallyObserved,
        const std::string& productIdentity, bool productAvailable, bool effectiveValid,
        const std::string& stateIdentity = {}, bool stateAvailable = true)
    {
        ZhangUserArSeparatedArcAges out;
        out.physical = physical_.update(key, time, physicalIdentity, physicallyObserved);
        const std::string resolvedStateIdentity = stateIdentity.empty() ? physicalIdentity : stateIdentity;
        out.state = state_.update(key, time, resolvedStateIdentity,
            stateAvailable && physicallyObserved);
        out.product = product_.update(key, time, productIdentity, productAvailable);
        const std::string effectiveIdentity = physicallyObserved && productAvailable
            ? zhangUserArEncodedField(physicalIdentity) + zhangUserArEncodedField(resolvedStateIdentity) +
                zhangUserArEncodedField(productIdentity)
            : std::string{};
        out.effective = effective_.update(key, time, effectiveIdentity,
            effectiveValid && out.physical.valid && out.state.valid && out.product.valid);
        return out;
    }
    std::array<ZhangUserArcTrackerCheckpoint,4> checkpoint(const std::string& runtime) const
    {return {physical_.checkpoint(runtime),state_.checkpoint(runtime),product_.checkpoint(runtime),effective_.checkpoint(runtime)};}
    void restoreCheckpoint(const std::array<ZhangUserArcTrackerCheckpoint,4>& value,const std::string& runtime,double completed)
    {
        auto candidate=*this;
        candidate.physical_.restoreCheckpoint(value[0],runtime,completed);
        candidate.state_.restoreCheckpoint(value[1],runtime,completed);
        candidate.product_.restoreCheckpoint(value[2],runtime,completed);
        candidate.effective_.restoreCheckpoint(value[3],runtime,completed);
        *this=std::move(candidate);
    }

private:
    static ZhangUserArExperimentPolicy ageOnlyPolicy(ZhangUserArExperimentPolicy policy)
    { policy.minimumArcSeconds = 0; policy.minimumArcSamples = 1; return policy; }
    ZhangUserArExperimentArcTracker physical_, state_, product_, effective_;
};

// Integers are strings to preserve arbitrary precision. Canonical decimal
// conversion rejects whitespace, decimal points, exponents and sign-only text.
inline bool zhangUserArCanonicalInteger(std::string input, std::string& output)
{
    if (input.empty()) return false;
    std::size_t begin = 0;
    const bool negative = input.front() == '-';
    if (negative || input.front() == '+') ++begin;
    if (begin == input.size()) return false;
    for (std::size_t i = begin; i < input.size(); ++i)
        if (input[i] < '0' || input[i] > '9') return false;
    while (begin + 1 < input.size() && input[begin] == '0') ++begin;
    output = input.substr(begin);
    if (negative && output != "0") output.insert(output.begin(), '-');
    return true;
}

struct ZhangUserArNamedInteger
{
    std::string target;
    std::string identity;
    std::string wideLane;
    std::string first;
};

// The full named WL/L1 set must be passed; an omitted member deliberately
// changes the signature. Reference changes also reset confirmation (no implicit
// basis transport). Sorting removes irrelevant input insertion order.
inline std::string zhangUserArCandidateSignature(
    const std::string& reference,
    const std::string& referenceIdentity,
    std::vector<ZhangUserArNamedInteger> targets)
{
    if (reference.empty() || referenceIdentity.empty() || targets.empty()) return {};
    std::sort(targets.begin(), targets.end(), [](const auto& a, const auto& b)
    {
        return a.target < b.target;
    });
    std::string out = zhangUserArEncodedField(reference) +
        zhangUserArEncodedField(referenceIdentity) +
        zhangUserArEncodedField(std::to_string(targets.size()));
    std::string previous;
    for (const auto& target : targets)
    {
        std::string wideLane, first;
        if (target.target.empty() || target.target == reference ||
            target.target == previous || target.identity.empty() ||
            !zhangUserArCanonicalInteger(target.wideLane, wideLane) ||
            !zhangUserArCanonicalInteger(target.first, first)) return {};
        previous = target.target;
        out += zhangUserArEncodedField(target.target) +
            zhangUserArEncodedField(target.identity) +
            zhangUserArEncodedField(wideLane) +
            zhangUserArEncodedField(first);
    }
    return out;
}

struct ZhangUserArFirstFeedbackDecision
{
    bool eligible = false;
    bool released = false;
    double first = 0;
    double last = 0;
    std::size_t count = 0;
    std::string signature;
    std::string reason = "NO_CANDIDATE";
};

// This gate stores candidate metadata only, never estimator state. Before its
// first successful commit, callers MUST search disposable copies of the FLOAT
// posterior and suppress all factor-history/held-registry side effects.
class ZhangUserArFirstFeedbackGate
{
public:
    explicit ZhangUserArFirstFeedbackGate(
        ZhangUserArExperimentPolicy policy = {}) : policy_(policy) {}

    ZhangUserArFirstFeedbackDecision observe(
        const ZhangUserArExperimentKey& key,
        double time,
        const std::string& signature,
        bool valid = true)
    {
        if (!key.valid())
        {
            ZhangUserArFirstFeedbackDecision out;
            out.reason = "INVALID_KEY";
            return out;
        }
        auto& state = states_[key];
        const bool usable = valid && !signature.empty() && std::isfinite(time);
        if (state.released)
        {
            state.eligible = usable;
            state.reason = usable ? "FIRST_FEEDBACK_ALREADY_COMMITTED" : "INVALID_CANDIDATE";
            return state;
        }
        if (!usable)
        {
            state = {};
            state.reason = "INVALID_OR_EMPTY_CANDIDATE";
            return state;
        }
        if (state.count && time < state.last - policy_.timeToleranceSeconds)
        {
            state = {};
            state.reason = "TIME_REVERSED";
            return state;
        }

        std::string resetReason;
        if (!state.count) resetReason = "FIRST_CANDIDATE";
        else if (state.signature != signature) resetReason = "CANDIDATE_CHANGED";
        else
        {
            const double delta = time - state.last;
            if (delta > policy_.maximumGapSeconds) resetReason = "CANDIDATE_GAP";
            else if (delta > policy_.timeToleranceSeconds &&
                std::abs(delta - policy_.confirmationIntervalSeconds) > policy_.timeToleranceSeconds)
                resetReason = "CANDIDATE_CADENCE_CHANGED";
        }
        if (!resetReason.empty())
        {
            state = {};
            state.signature = signature;
            state.first = state.last = time;
            state.count = 1;
            state.reason = std::move(resetReason);
        }
        else if (time - state.last > policy_.timeToleranceSeconds)
        {
            state.last = time;
            ++state.count;
            state.reason = "CONSISTENT_CANDIDATE";
        }
        else state.reason = "DUPLICATE_EPOCH";
        state.eligible = state.count >= std::max<std::size_t>(1, policy_.confirmationEpochs);
        return state;
    }

    // Only call after the actual estimator transaction has succeeded. A false
    // receipt or mismatched candidate resets pending evidence, never releases.
    bool commit(
        const ZhangUserArExperimentKey& key,
        double time,
        const std::string& signature,
        bool success)
    {
        auto it = states_.find(key);
        if (it == states_.end()) return false;
        auto& state = it->second;
        if (state.released) return success;
        if (!success || !std::isfinite(time) || !state.eligible ||
            std::abs(time - state.last) > policy_.timeToleranceSeconds ||
            signature != state.signature)
        {
            state = {};
            state.reason = success ? "COMMIT_RECEIPT_MISMATCH" : "COMMIT_FAILED";
            return false;
        }
        state.released = true;
        state.reason = "FIRST_FEEDBACK_COMMITTED";
        return true;
    }

    bool released(const ZhangUserArExperimentKey& key) const
    {
        auto it = states_.find(key);
        return it != states_.end() && it->second.released;
    }
    bool isReleased(const ZhangUserArExperimentKey& key) const { return released(key); }
    void erase(const ZhangUserArExperimentKey& key) { states_.erase(key); }
    void clear() { states_.clear(); }
    void clearRuntime(const std::string& runtime)
    {
        for (auto it = states_.begin(); it != states_.end();)
            if (it->first.runtime == runtime) it = states_.erase(it);
            else ++it;
    }
    std::size_t size() const { return states_.size(); }

private:
    ZhangUserArExperimentPolicy policy_;
    std::map<ZhangUserArExperimentKey, ZhangUserArFirstFeedbackDecision> states_;
};
