#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "common/algebra.hpp"
#include "common/biases.hpp"

// External products often supply only marginal standard deviations. A zero
// off-diagonal in our working model is an independence ASSUMPTION, never proof
// that missing satellite/frequency/clock correlations are zero. The current
// internal files likewise have no cross-epoch joint covariance.
enum class ZhangUserProductErrorPolicy { WHITE, SEGMENT_CONSTANT, EXTERNAL_FIXED_CONDITIONAL };
enum class ZhangUserProductErrorUnsupported { WHITE_APPROX, FAIL };

inline ZhangUserProductErrorPolicy zhangParseUserProductErrorPolicy(const std::string& text)
{
    if (text.empty() || text == "WHITE") return ZhangUserProductErrorPolicy::WHITE;
    if (text == "SEGMENT_CONSTANT") return ZhangUserProductErrorPolicy::SEGMENT_CONSTANT;
    if (text == "EXTERNAL_FIXED_CONDITIONAL") return ZhangUserProductErrorPolicy::EXTERNAL_FIXED_CONDITIONAL;
    throw std::runtime_error("INVALID_ZHANG_USER_PRODUCT_ERROR_POLICY");
}

inline ZhangUserProductErrorUnsupported zhangParseUserProductErrorUnsupported(const std::string& text)
{
    if (text.empty() || text == "WHITE_APPROX") return ZhangUserProductErrorUnsupported::WHITE_APPROX;
    if (text == "FAIL") return ZhangUserProductErrorUnsupported::FAIL;
    throw std::runtime_error("INVALID_ZHANG_USER_PRODUCT_ERROR_UNSUPPORTED");
}

inline ZhangUserProductErrorPolicy zhangUserProductErrorPolicy()
{
    static const auto policy = []()
    {
        const char* value = std::getenv("ZHANG_USER_PRODUCT_ERROR_POLICY");
        return zhangParseUserProductErrorPolicy(value ? value : "WHITE");
    }();
    return policy;
}

inline ZhangUserProductErrorUnsupported zhangUserProductErrorUnsupported()
{
    static const auto policy = []()
    {
        const char* value = std::getenv("ZHANG_USER_PRODUCT_ERROR_UNSUPPORTED");
        return zhangParseUserProductErrorUnsupported(value ? value : "WHITE_APPROX");
    }();
    return policy;
}

inline const char* zhangUserProductErrorPolicyName(ZhangUserProductErrorPolicy policy)
{
    switch (policy)
    {
        case ZhangUserProductErrorPolicy::WHITE: return "WHITE";
        case ZhangUserProductErrorPolicy::SEGMENT_CONSTANT: return "SEGMENT_CONSTANT";
        case ZhangUserProductErrorPolicy::EXTERNAL_FIXED_CONDITIONAL: return "EXTERNAL_FIXED_CONDITIONAL";
    }
    throw std::runtime_error("INVALID_ZHANG_USER_PRODUCT_ERROR_POLICY_ENUM");
}

inline void zhangUserValidateProductErrorProvider(
    ZhangUserProductErrorPolicy policy, bool internalProvider)
{
    if (internalProvider && policy == ZhangUserProductErrorPolicy::EXTERNAL_FIXED_CONDITIONAL)
        throw std::runtime_error("EXTERNAL_FIXED_CONDITIONAL_CANNOT_IGNORE_INTERNAL_JOINT_PRODUCT_COVARIANCE");
}

struct ZhangUserExternalBiasVarianceDecision
{
    double rawInputVariance = 0;
    double modeledVariance = 0;
    bool fixedConditional = false;
};

// This is a conditional experiment with nominal external OSBs held fixed,
// NOT an assertion that their physical error or missing covariance is zero.
// Only the propagated CODE/PHAS OSB variance is omitted. Satellite clock and
// orbit uncertainty retain their separately configured observation models.
// No nominal bias is changed, and an explicitly estimated bias's -1 noise
// sentinel/continuous prior must remain owned by that original KF state.
inline ZhangUserExternalBiasVarianceDecision zhangUserExternalBiasVariance(
    ZhangUserProductErrorPolicy policy, bool internalProvider,
    bool explicitlyEstimated, bool nominalFound, bool directNominalValid,
    double appliedBias, double appliedVariance)
{
    zhangUserValidateProductErrorProvider(policy, internalProvider);
    ZhangUserExternalBiasVarianceDecision result;
    result.rawInputVariance = result.modeledVariance = appliedVariance;
    if (policy != ZhangUserProductErrorPolicy::EXTERNAL_FIXED_CONDITIONAL || explicitlyEstimated)
        return result;
    if (!nominalFound || !directNominalValid || !std::isfinite(appliedBias) ||
        !std::isfinite(appliedVariance) || appliedVariance < 0)
        throw std::runtime_error("EXTERNAL_FIXED_CONDITIONAL_REQUIRES_VALID_DIRECT_NOMINAL_BSX");
    result.modeledVariance = 0;
    result.fixedConditional = true;
    return result;
}

inline void zhangUserValidateProductErrorEpochAlignment(
    ZhangUserProductErrorPolicy policy, GTime observationTime, GTime epochTime)
{
    // The first segment-state implementation retires immutable records in
    // nominal-epoch measurement assembly. A mixed-time validity-boundary
    // batch needs referenced-key-aware retirement, not per-row expiry.
    // Reject that unsupported model instead of dropping a referenced prior.
    if (policy == ZhangUserProductErrorPolicy::SEGMENT_CONSTANT &&
        (observationTime == GTime::noTime() || epochTime == GTime::noTime() ||
         observationTime != epochTime))
        throw std::runtime_error("SEGMENT_CONSTANT_PRODUCT_ERRORS_REQUIRE_EXACT_EPOCH_ALIGNMENT");
}

inline void zhangUserMarkProductUncertaintyIncomplete(KFState& state, const std::string& reason)
{
    const std::lock_guard<std::recursive_mutex> guard(state.kfStateMutex);
    state.metaDataMap["ZHANG_USER_PRODUCT_UNCERTAINTY_ASSESSED"] = "1";
    // Set-like metadata is bounded by the model's reason vocabulary, not by
    // elapsed epochs or historic records. KF copies retain these receipts.
    state.metaDataMap["ZHANG_USER_PRODUCT_UNCERTAINTY_MISSING:" + reason] = "1";
}

inline bool zhangUserProductUncertaintyComplete(const KFState& state)
{
    const std::lock_guard<std::recursive_mutex> guard(state.kfStateMutex);
    if (!state.metaDataMap.count("ZHANG_USER_PRODUCT_UNCERTAINTY_ASSESSED")) return false;
    for (const auto& [key, value] : state.metaDataMap)
        if (key.find("ZHANG_USER_PRODUCT_UNCERTAINTY_MISSING:") == 0 && value == "1") return false;
    return true;
}

struct ZhangUserConstantProductErrorRecord
{
    bool supported = false;
    std::string reason = "NO_DIRECT_VALID_BIAS_RECORD";
    std::string owner;
    std::string identity;
    SatSys satellite;
    int code = 0;
    long double validFrom = 0;
    long double validTo = 0;
    double variance = 0;
};

// Only a direct, single constant BSX record is supported. A chained/interpolated
// OSB requires the constituent random variables and their transformation;
// freezing its epochwise aggregate variance would invent a temporal model.
inline ZhangUserConstantProductErrorRecord zhangUserDirectConstantBiasRecord(
    const SatSys& satellite, E_ObsCode code, E_MeasType measurement, GTime time,
    double appliedBias, double appliedVariance, bool oneFrequencyPhaseBias,
    bool requireKnownNonzeroVariance = true)
{
    ZhangUserConstantProductErrorRecord result;
    result.satellite = satellite;
    E_ObsCode lookupCode = code;
    if (oneFrequencyPhaseBias && measurement == PHAS)
    {
        const auto frequencies = code2Freq.find(satellite.sys);
        if (frequencies == code2Freq.end() || !frequencies->second.count(code)) return result;
        lookupCode = freq2CodeHax(satellite.sys, frequencies->second.at(code));
    }
    // When one-frequency lookup aliases several physical codes to the same
    // source OSB, they share ONE error variable, not independent duplicates.
    result.code = static_cast<int>(lookupCode);
    result.owner = "ZHANG_USER_PRODUCT_ERROR:" + satellite.id() + ":" +
        std::string(1, satellite.sysChar()) + ':' +
        std::to_string(static_cast<int>(measurement)) + ':' + std::to_string(static_cast<int>(lookupCode));
    // GLO's device lookup differs from GPS. Do not accidentally authorize a
    // different record; the currently validated profile is GPS direct OSBs.
    if (satellite.sys != E_Sys::GPS)
    {
        result.reason = "UNSUPPORTED_NON_GPS_DIRECT_BIAS_LOOKUP";
        return result;
    }
    const std::string deviceId = satellite.id() + ":" + std::string(1, satellite.sysChar());
    const auto device = biasMaps[measurement].find(deviceId);
    if (device == biasMaps[measurement].end()) return result;
    const auto first = device->second.find(lookupCode);
    if (first == device->second.end()) return result;
    const auto zero = first->second.find(E_ObsCode::NONE);
    if (zero == first->second.end()) return result;
    const auto selected = zero->second.lower_bound(time);
    if (selected == zero->second.end()) return result;
    const auto& entry = selected->second;
    if (entry.source != "bsx" || entry.measType != measurement ||
        entry.Sat != satellite || !entry.name.empty() || entry.cod1 != lookupCode ||
        entry.cod2 != E_ObsCode::NONE || entry.tini == GTime::noTime() ||
        entry.tfin == GTime::noTime() || time < entry.tini || !(time < entry.tfin)) return result;
    if (!std::isfinite(entry.bias) || !std::isfinite(entry.var) || entry.var < 0 ||
        !std::isfinite(entry.slop) || !std::isfinite(entry.slpv) || entry.slpv < 0)
    {
        result.reason = "BIAS_VALUE_OR_STD_UNAVAILABLE";
        return result;
    }
    if (entry.var == 0 && requireKnownNonzeroVariance)
    {
        // The legacy BSX parser maps a missing/unparseable std field to zero.
        // Without an availability bit, zero cannot authorize an exact prior.
        // WHITE compatibility is untouched, but this is explicitly UNSUPPORTED.
        result.reason = "BIAS_STD_ZERO_OR_MISSING_NOT_DISTINGUISHABLE";
        return result;
    }
    if (entry.slop != 0 || entry.slpv != 0)
    {
        result.reason = "BIAS_SLOPE_RANDOM_VARIABLE_NOT_IMPLEMENTED";
        return result;
    }
    if (selected != zero->second.begin() && std::prev(selected)->second.refTime != entry.refTime)
    {
        result.reason = "BIAS_INTERPOLATION_RANDOM_VARIABLES_NOT_IMPLEMENTED";
        return result;
    }
    if (!std::isfinite(appliedBias) || !std::isfinite(appliedVariance) || appliedVariance < 0 ||
        std::abs(appliedBias - entry.bias) > 1e-12 * std::max(1.0, std::abs(entry.bias)) ||
        std::abs(appliedVariance - entry.var) > 1e-12 * std::max(1.0, std::abs(entry.var)))
    {
        result.reason = "APPLIED_BIAS_DIFFERS_FROM_IMMUTABLE_RECORD";
        return result;
    }
    std::ostringstream token;
    token << std::setprecision(24) << deviceId << ':' << static_cast<int>(measurement)
          << ':' << static_cast<int>(lookupCode) << ':' << entry.tini.bigTime
          << ':' << entry.tfin.bigTime << ':' << entry.refTime.bigTime
          << ':' << std::setprecision(17) << entry.bias << ':' << entry.var;
    result.identity = token.str();
    result.validFrom = entry.tini.bigTime;
    result.validTo = entry.tfin.bigTime;
    result.variance = entry.var;
    result.supported = true;
    result.reason = requireKnownNonzeroVariance ? "DIRECT_CONSTANT_BSX_RECORD" :
        "DIRECT_CONSTANT_BSX_NOMINAL_CONDITIONAL_STD_NOT_AUTHORIZED_AS_TRUTH";
    return result;
}

inline bool zhangUserProductErrorKey(const KFKey& key)
{
    return key.type == KF::CALC && key.str.find("ZHANG_USER_PRODUCT_ERROR:") == 0;
}

inline bool zhangUserRequiresJointProductErrorFilter(const KFState& state)
{
    const std::lock_guard<std::recursive_mutex> guard(state.kfStateMutex);
    // Sharing a latent product error can also correlate different receivers.
    // Its later retirement does NOT erase those retained cross-covariances.
    // Keep the joint-filter requirement for the lifetime of this KF instance.
    if (state.metaDataMap.count("ZHANG_USER_PRODUCT_ERROR_JOINT_FILTER_REQUIRED")) return true;
    for (const auto& [key, index] : state.kfIndexMap)
        if (zhangUserProductErrorKey(key)) return true;
    for (const auto& [key, transition] : state.stateTransitionMap)
        if (zhangUserProductErrorKey(key)) return true;
    return false;
}

inline void zhangUserValidateConditionalProductErrorState(
    ZhangUserProductErrorPolicy policy, const KFState& state)
{
    // A posterior that learned latent errors (including after retirement)
    // cannot be relabelled as conditional-fixed by removing noise at AR time.
    // Start a fresh FLOAT reprocessing run under the new stochastic model.
    if (policy == ZhangUserProductErrorPolicy::EXTERNAL_FIXED_CONDITIONAL &&
        zhangUserRequiresJointProductErrorFilter(state))
        throw std::runtime_error("EXTERNAL_FIXED_CONDITIONAL_REQUIRES_FRESH_FLOAT_NOT_SEGMENT_ERROR_POSTERIOR");
}

inline KFKey zhangUserProductErrorStateKey(const ZhangUserConstantProductErrorRecord& record)
{
    KFKey key;
    // CALC is a continuous nuisance parameter, not an integer ambiguity or an
    // estimated physical satellite bias. Namespaced identity prevents aliasing.
    key.type = KF::CALC;
    key.Sat = record.satellite;
    key.num = record.code;
    key.str = record.owner + '|' + record.identity;
    std::ostringstream expiry;
    expiry << "ZHANG_PRODUCT_ERROR_EXPIRES:" << std::setprecision(24) << record.validTo;
    key.comment = expiry.str();
    return key;
}

inline long double zhangUserProductErrorExpiry(const KFKey& key)
{
    constexpr const char* prefix = "ZHANG_PRODUCT_ERROR_EXPIRES:";
    if (key.comment.find(prefix) != 0)
        throw std::runtime_error("PRODUCT_ERROR_STATE_EXPIRY_MISSING");
    const std::string text = key.comment.substr(std::char_traits<char>::length(prefix));
    std::size_t consumed = 0;
    const long double result = std::stold(text, &consumed);
    if (consumed != text.size() || !std::isfinite(result))
        throw std::runtime_error("PRODUCT_ERROR_STATE_EXPIRY_INVALID");
    return result;
}

inline unsigned zhangUserRetireExpiredProductErrors(KFState& state, GTime time, Trace& trace)
{
    const std::lock_guard<std::recursive_mutex> guard(state.kfStateMutex);
    std::vector<KFKey> retire;
    for (const auto& [key, transition] : state.stateTransitionMap)
        if (zhangUserProductErrorKey(key) && time.bigTime >= zhangUserProductErrorExpiry(key))
            retire.push_back(key);
    for (const auto& key : retire)
    {
        state.removeState(key, false);
        trace << "\nZHANG_USER_PRODUCT_ERROR_RETIRE time=" << time.to_string(0)
              << " key=" << std::quoted(std::string(key))
              << " reason=VALIDITY_INTERVAL_ENDED operation=STATE_TRANSITION_MARGINALIZATION"
              << " hard_fix=0 integer_state_removed=0";
    }
    return static_cast<unsigned>(retire.size());
}

struct ZhangUserProductErrorAttachment
{
    KFKey key;
    bool priorInjected = false;
    unsigned retired = 0;
    double mean = 0;
    double posteriorVariance = 0;
};

inline ZhangUserProductErrorAttachment zhangUserAttachConstantProductError(
    KFState& state, KFMeasEntry& measurement,
    const ZhangUserConstantProductErrorRecord& record, GTime time,
    double coefficient, Trace& trace)
{
    if (!record.supported || !std::isfinite(record.variance) || record.variance < 0 ||
        !(record.validFrom <= time.bigTime && time.bigTime < record.validTo) ||
        !std::isfinite(coefficient) || coefficient == 0)
        throw std::runtime_error("INVALID_CONSTANT_PRODUCT_ERROR_ATTACHMENT");
    const std::lock_guard<std::recursive_mutex> guard(state.kfStateMutex);
    ZhangUserProductErrorAttachment out;
    out.key = zhangUserProductErrorStateKey(record);
    std::vector<KFKey> retire;
    const std::string ownerPrefix = record.owner + '|';
    for (const auto& [key, transition] : state.stateTransitionMap)
        if (zhangUserProductErrorKey(key) && key.str.find(ownerPrefix) == 0 && key != out.key)
            retire.push_back(key);
    for (const auto& key : retire)
    {
        state.removeState(key, false);
        trace << "\nZHANG_USER_PRODUCT_ERROR_RETIRE time=" << time.to_string(0)
              << " key=" << std::quoted(std::string(key))
              << " reason=IMMUTABLE_RECORD_REPLACED operation=STATE_TRANSITION_MARGINALIZATION"
              << " hard_fix=0 integer_state_removed=0";
    }
    out.retired = static_cast<unsigned>(retire.size());
    if (record.variance == 0) return out;
    state.metaDataMap["ZHANG_USER_PRODUCT_ERROR_JOINT_FILTER_REQUIRED"] = "1";
    const bool exists = state.stateTransitionMap.count(out.key) != 0;
    InitialState initial;
    initial.estimate = true;
    initial.x = 0;
    initial.P = record.variance;
    initial.Q = 0;
    initial.tau = -1;
    initial.comment = out.key.comment;
    out.posteriorVariance = record.variance;
    if (exists)
        state.getKFValue(out.key, initial.x, &out.posteriorVariance, nullptr, false);
    out.mean = initial.x;
    // addKFState is idempotent: existing P is NOT reinitialized. No repeated
    // zero-valued pseudo-observation/prior is injected each epoch.
    measurement.addDsgnEntry(out.key, coefficient, initial);
    out.priorInjected = !exists;
    return out;
}

inline void zhangUserProductErrorUnsupportedCheck(const std::string& reason)
{
    if (zhangUserProductErrorPolicy() == ZhangUserProductErrorPolicy::SEGMENT_CONSTANT &&
        zhangUserProductErrorUnsupported() == ZhangUserProductErrorUnsupported::FAIL)
        throw std::runtime_error("SEGMENT_CONSTANT_PRODUCT_ERROR_UNAVAILABLE:" + reason);
}
