#pragma once
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <stdexcept>
#include <string>

// Explicit research profile. These values are screening rules under the
// assumed FLOAT stochastic model, NOT a calibrated fixed-failure-rate test.
// Legacy experiments remain reproducible; safety changes require an explicit
// profile receipt and must never globally disable ratio-only (which would
// resurrect unrelated historical NIS paths).
struct ZhangUserArReliabilityPolicy
{
    bool safe = false;
    bool gonsis = false;
    double minimumArcSeconds = 0;
    std::size_t minimumArcSamples = 1;
    std::size_t confirmationEpochs = 1;
    bool enforceBootstrapSuccess = false;
    double maximumConditionalSigmaCycles = 0;
    int minimumParRank = 1;
    int maximumParTrials = 0;
    int minimumSupportSatellites = 2;
    bool contradictionIsolation = false;
    const char* name = "LEGACY_RATIO";
};

inline double zhangUserArPolicyNumber(const char* name, double fallback,
    double minimum, double maximum)
{
    const char* input = std::getenv(name);
    if (!input) return fallback;
    char* end = nullptr;
    const double value = std::strtod(input, &end);
    if (end == input || *end || !std::isfinite(value) || value < minimum || value > maximum)
        throw std::runtime_error(std::string("INVALID_USER_AR_POLICY:") + name);
    return value;
}

inline int zhangUserArPolicyCount(const char* name, int fallback, int minimum, int maximum)
{
    const double value = zhangUserArPolicyNumber(name, fallback, minimum, maximum);
    if (value != std::floor(value))
        throw std::runtime_error(std::string("NONINTEGER_USER_AR_POLICY:") + name);
    return static_cast<int>(value);
}

inline ZhangUserArReliabilityPolicy zhangUserArReliabilityPolicy()
{
    ZhangUserArReliabilityPolicy out;
    const char* profile = std::getenv("ZHANG_USER_AR_RELIABILITY_PROFILE");
    if (profile && std::string(profile) != "LEGACY_RATIO" &&
        std::string(profile) != "SAFE_RESEARCH" && std::string(profile) != "GONSIS_UDUC")
        throw std::runtime_error("INVALID_USER_AR_RELIABILITY_PROFILE");
    out.safe = profile && std::string(profile) == "SAFE_RESEARCH";
    out.gonsis = profile && std::string(profile) == "GONSIS_UDUC";
    if (out.safe)
    {
        out.name = "SAFE_RESEARCH";
        out.minimumArcSeconds = 120;
        out.minimumArcSamples = 5;
        out.confirmationEpochs = 3;
        out.enforceBootstrapSuccess = true;
        out.maximumConditionalSigmaCycles = .15;
        out.minimumParRank = 2;
        out.maximumParTrials = 4;
        out.contradictionIsolation = true;
    }
    if (out.gonsis)
    {
        // Reimplementation of the uploaded GONSIS UCDF path, not a weaker
        // SAFE_RESEARCH suffix. Source lockhold is zero; position precision,
        // WL rounding, complete L1 ILS and output residual gates are separate.
        out.name = "GONSIS_UDUC";
        out.minimumParRank = 4;
    }
    out.minimumArcSeconds = zhangUserArPolicyNumber("ZHANG_USER_AR_MIN_ARC_SECONDS",
        out.minimumArcSeconds, 0, 86400);
    out.minimumArcSamples = zhangUserArPolicyCount("ZHANG_USER_AR_MIN_ARC_SAMPLES",
        out.safe ? 5 : (out.minimumArcSeconds > 0 ? 5 : 1), 1, 100000);
    out.confirmationEpochs = zhangUserArPolicyCount("ZHANG_USER_AR_CONFIRM_EPOCHS",
        out.confirmationEpochs, 1, 100);
    out.maximumConditionalSigmaCycles = zhangUserArPolicyNumber(
        "ZHANG_USER_AR_MAX_CONDITIONAL_SIGMA_CYCLES", out.maximumConditionalSigmaCycles, 0, 100);
    out.minimumParRank = zhangUserArPolicyCount("ZHANG_USER_AR_PAR_MIN_RANK", out.minimumParRank, 1, 10000);
    out.maximumParTrials = zhangUserArPolicyCount("ZHANG_USER_AR_PAR_MAX_TRIALS", out.maximumParTrials, 0, 10000);
    out.minimumSupportSatellites = zhangUserArPolicyCount("ZHANG_USER_AR_MIN_SUPPORT_SATELLITES",
        out.minimumSupportSatellites, 2, 1000);
    // A safe profile cannot silently fall back to its old permissive policy.
    if (out.safe && (out.minimumArcSeconds < 120 || out.minimumArcSamples < 5 ||
        out.confirmationEpochs < 2 || out.minimumParRank < 2 || out.maximumParTrials < 1 ||
        !(out.maximumConditionalSigmaCycles > 0)))
        throw std::runtime_error("SAFE_RESEARCH_POLICY_WEAKENING_REJECTED");
    if (out.gonsis && (out.confirmationEpochs != 1 ||
        out.maximumConditionalSigmaCycles != 0 || out.minimumParRank < 4))
        throw std::runtime_error("GONSIS_UDUC_INCOMPATIBLE_SHARED_RELIABILITY_GATES");
    return out;
}

// HOLD retains exact integer evidence separately from the authoritative FLOAT
// filter. Both modes output a disposable conditional branch each epoch.
inline bool zhangUserConditionalHoldRequested()
{
    const char* policy = std::getenv("ZHANG_USER_AR_FEEDBACK");
    if (!policy || std::string(policy) != "CONDITIONAL_HOLD") return false;
    if (!zhangUserArReliabilityPolicy().gonsis)
        throw std::runtime_error("CONDITIONAL_HOLD_REQUIRES_GONSIS_UDUC");
    return true;
}

inline bool zhangUserQuotientHoldRequested()
{
    const char* policy = std::getenv("ZHANG_USER_AR_FEEDBACK");
    if (!policy || std::string(policy) == "EPOCH_ONLY") return false;
    if (zhangUserConditionalHoldRequested()) return false; // never enable the legacy HOLD implementation
    if (std::string(policy) == "QUOTIENT_HOLD")
    {
        const char* research = std::getenv("ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH");
        const auto reliability = zhangUserArReliabilityPolicy();
        if (reliability.safe || reliability.gonsis || !research || std::string(research) != "1")
            throw std::runtime_error("USER_QUOTIENT_HOLD_NOT_DEPLOYABLE; legacy research requires explicit opt-in");
        if (reliability.confirmationEpochs != 1)
            throw std::runtime_error("LEGACY_HOLD_RESEARCH_DOES_NOT_SUPPORT_CANDIDATE_CONFIRMATION");
        return true;
    }
    throw std::runtime_error("INVALID_USER_AR_FEEDBACK: use EPOCH_ONLY or CONDITIONAL_HOLD (legacy QUOTIENT_HOLD is research-only)");
}

inline bool zhangUserSharedSearchRequested()
{
    const char* search = std::getenv("ZHANG_USER_AR_SEARCH");
    const auto reliability = zhangUserArReliabilityPolicy();
    if (reliability.gonsis)
    {
        const char* qc = std::getenv("ZHANG_USER_QC");
        if (!qc || std::string(qc) != "CAUSAL_TRACKING")
            throw std::runtime_error("GONSIS_UDUC_REQUIRES_CAUSAL_TRACKING_QC");
        if (search && std::string(search) != "GONSIS_UDUC")
            throw std::runtime_error("GONSIS_UDUC_REQUIRES_MATCHED_SEARCH_PROFILE");
        if (zhangUserQuotientHoldRequested())
            throw std::runtime_error("GONSIS_UDUC_REQUIRES_EPOCH_ONLY");
        return false;
    }
    if (reliability.safe)
    {
        const char* qc = std::getenv("ZHANG_USER_QC");
        if (!qc || std::string(qc) != "CAUSAL_TRACKING")
            throw std::runtime_error("SAFE_RESEARCH_REQUIRES_CAUSAL_TRACKING_QC");
        if (search && std::string(search) != "SHARED_WL_L1_PAR")
            throw std::runtime_error("SAFE_RESEARCH_REQUIRES_SHARED_WL_L1_PAR");
    }
    if (!search || std::string(search) == "SHARED_WL_L1_PAR") return true;
    if (std::string(search) == "LEGACY_EPOCH" && !zhangUserQuotientHoldRequested()) return false;
    throw std::runtime_error("INVALID_USER_AR_SEARCH_OR_LEGACY_HOLD");
}

inline bool zhangUserGonsisSearchRequested()
{
    if (!zhangUserArReliabilityPolicy().gonsis) return false;
    // Validate the same explicit QC/search/feedback boundary before dispatch.
    zhangUserSharedSearchRequested();
    return true;
}
