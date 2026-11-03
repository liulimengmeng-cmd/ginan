#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "common/zhangProductPublicationGate.hpp"

/** Numerical evidence is separate from exact physical/chart transport.
 * In particular an expired reference is NOT a measured zero correction step.
 * This helper does not establish any integer or structural certificate. */
enum class ZhangProductNumericalReferenceKind
{
    INITIAL_BASELINE,
    NEW_SEGMENT_BASELINE,
    INCOMPARABLE_NEW_DATUM,
    COMPARABLE,
    GAP_EXCEEDED,
    UNASSESSED,
    INVALID
};

enum class ZhangProductComparisonDomain
{
    LEGACY_PUBLISHED_IDENTITY,
    SAME_VERIFIED_DOMAIN,
    DIFFERENT_CURRENT_DOMAIN,
    UNASSESSED
};

struct ZhangProductNumericalReferenceInput
{
    bool currentNumericValid = false;
    bool hasPrevious = false;
    bool samePublishedIdentity = false;
    ZhangProductComparisonDomain comparisonDomain =
        ZhangProductComparisonDomain::LEGACY_PUBLISHED_IDENTITY;
    ZhangProductPublicationRoute nominalRoute =
        ZhangProductPublicationRoute::SUSPENDED;
    double currentCorrection = std::numeric_limits<double>::quiet_NaN();
    double previousCorrection = std::numeric_limits<double>::quiet_NaN();
    double currentSigma = std::numeric_limits<double>::quiet_NaN();
    double previousSigma = std::numeric_limits<double>::quiet_NaN();
    double gapSeconds = std::numeric_limits<double>::quiet_NaN();
    double maximumGapSeconds = std::numeric_limits<double>::quiet_NaN();
    // Only undo the deterministic offset staged by THIS writer invocation
    // while testing an otherwise unchanged native/physical comparison chart.
    // This is a numerical outlier comparison, never an old integer RHS or a
    // certificate of temporal continuity. The emitted product keeps its NEW
    // datum and is still required to open a new public segment.
    bool removeStagedIntegerDatumShift = false;
    bool stagedIntegerDatumComparisonProven = false;
    double correctionBeforeStagedDatumChange = std::numeric_limits<double>::quiet_NaN();
    double stagedCorrectionShiftMetres = std::numeric_limits<double>::quiet_NaN();
};

struct ZhangProductNumericalReference
{
    ZhangProductNumericalReferenceKind kind =
        ZhangProductNumericalReferenceKind::UNASSESSED;
    bool hasDelta = false;
    double delta = std::numeric_limits<double>::quiet_NaN();
    const char* reason = "NUMERICAL_REFERENCE_NOT_EVALUATED";
};

inline ZhangProductNumericalReference zhangProductNumericalReference(
    const ZhangProductNumericalReferenceInput& input)
{
    using Kind = ZhangProductNumericalReferenceKind;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    if (!input.currentNumericValid || !std::isfinite(input.currentCorrection) ||
        !std::isfinite(input.currentSigma) || input.currentSigma < 0)
        return {Kind::INVALID, false, missing, "NUMERICAL_CURRENT_ENDPOINT_INVALID"};
    if (!std::isfinite(input.maximumGapSeconds) || input.maximumGapSeconds <= 0)
        return {Kind::INVALID, false, missing, "NUMERICAL_REFERENCE_GAP_LIMIT_INVALID"};
    double comparisonCorrection = input.currentCorrection;
    if (input.removeStagedIntegerDatumShift)
    {
        if (!input.stagedIntegerDatumComparisonProven || !input.hasPrevious ||
            !input.samePublishedIdentity ||
            input.nominalRoute != ZhangProductPublicationRoute::NEW_SEGMENT)
            return {Kind::INVALID, false, missing, "NUMERICAL_STAGED_DATUM_COMPARISON_UNPROVEN"};
        const double before = input.correctionBeforeStagedDatumChange;
        const double shift = input.stagedCorrectionShiftMetres;
        const double reconstructed = before + shift;
        if (!std::isfinite(before) || !std::isfinite(shift) || !std::isfinite(reconstructed))
            return {Kind::INVALID, false, missing, "NUMERICAL_STAGED_DATUM_SHIFT_INVALID"};
        const double tolerance = 64 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(before), std::abs(shift), std::abs(input.currentCorrection)});
        if (std::abs(input.currentCorrection - reconstructed) > tolerance)
            return {Kind::INVALID, false, missing, "NUMERICAL_STAGED_DATUM_SHIFT_MISMATCH"};
        comparisonCorrection = before;
    }
    if (!input.hasPrevious)
        return {Kind::INITIAL_BASELINE, false, missing, "INITIAL_NUMERICAL_BASELINE"};
    // A changed target/chart is not comparable merely because its public
    // counter has not yet advanced. Counter staging happens AFTER this test.
    // Conversely a nominal NEW route in an unchanged domain must still pass
    // the real outlier test below: NEW is not an escape from a measured jump.
    if (input.comparisonDomain == ZhangProductComparisonDomain::DIFFERENT_CURRENT_DOMAIN)
    {
        if (input.nominalRoute == ZhangProductPublicationRoute::NEW_SEGMENT)
            return {Kind::INCOMPARABLE_NEW_DATUM, false, missing,
                "CURRENT_TARGET_CHART_CHANGED_NEW_BASELINE"};
        return {Kind::UNASSESSED, false, missing,
            "CURRENT_TARGET_CHART_CHANGE_WITHOUT_NEW_SEGMENT"};
    }
    if (input.comparisonDomain == ZhangProductComparisonDomain::UNASSESSED)
        return {Kind::UNASSESSED, false, missing, "NUMERICAL_COMPARISON_DOMAIN_UNASSESSED"};
    // An already-required physical/publication segment change starts a new
    // baseline. It cannot claim numerical continuity with a different datum.
    if (!input.samePublishedIdentity && input.comparisonDomain !=
        ZhangProductComparisonDomain::SAME_VERIFIED_DOMAIN)
    {
        if (input.nominalRoute == ZhangProductPublicationRoute::NEW_SEGMENT)
            return {Kind::NEW_SEGMENT_BASELINE, false, missing,
                "EXPLICIT_NEW_SEGMENT_NUMERICAL_BASELINE"};
        return {Kind::UNASSESSED, false, missing,
            "NUMERICAL_REFERENCE_IDENTITY_MISMATCH"};
    }
    if (!std::isfinite(input.previousCorrection) ||
        !std::isfinite(input.previousSigma) || input.previousSigma < 0)
        return {Kind::INVALID, false, missing, "NUMERICAL_PREVIOUS_ENDPOINT_INVALID"};
    const double sigmaSum = input.currentSigma + input.previousSigma;
    if (!std::isfinite(sigmaSum) ||
        sigmaSum > std::numeric_limits<double>::max() / 10.0)
        return {Kind::INVALID, false, missing, "NUMERICAL_UNCERTAINTY_GUARD_OVERFLOW"};
    // These finite-arithmetic checks apply even to an expired reference: in
    // that case the missing comparison must be the ONLY missing evidence.
    const double rawDelta = input.currentCorrection - input.previousCorrection;
    const double delta = comparisonCorrection - input.previousCorrection;
    if (!std::isfinite(rawDelta) || !std::isfinite(delta))
        return {Kind::INVALID, false, missing, "NUMERICAL_REFERENCE_DELTA_NONFINITE"};
    if (!std::isfinite(input.gapSeconds) || input.gapSeconds <= 0)
        return {Kind::UNASSESSED, false, missing, "NUMERICAL_REFERENCE_TIME_INVALID"};
    if (input.gapSeconds > input.maximumGapSeconds)
        return {Kind::GAP_EXCEEDED, false, missing, "NUMERICAL_REFERENCE_GAP_EXCEEDED"};
    return {Kind::COMPARABLE, true, delta, input.removeStagedIntegerDatumShift
        ? "COMPARABLE_BEFORE_STAGED_INTEGER_DATUM_CHANGE" : "COMPARABLE_NUMERICAL_REFERENCE"};
}

/** Numerical comparability is NOT old-segment integer authorization. A NEW
 * product may anchor the clock common mode only when its independently
 * frozen current evidence and its actual comparison coordinates are valid.
 * For a NEW writer datum, any nonzero late gamma must have been removed by
 * the already-validated narrow normalization, never by assuming zero. */
inline bool zhangProductCommonModeAnchorEligible(
    const ZhangCurrentProductSpatialEvidence& current,
    const ZhangProductNumericalReferenceInput& input,
    const ZhangProductNumericalReference& reference)
{
    if (!current.nativeValid() || !current.currentIntegerCertificateValid ||
        !input.currentNumericValid || !input.hasPrevious ||
        input.comparisonDomain != ZhangProductComparisonDomain::SAME_VERIFIED_DOMAIN ||
        input.nominalRoute == ZhangProductPublicationRoute::SUSPENDED ||
        reference.kind != ZhangProductNumericalReferenceKind::COMPARABLE ||
        !reference.hasDelta || !std::isfinite(reference.delta))
        return false;
    if (input.nominalRoute == ZhangProductPublicationRoute::NEW_SEGMENT)
    {
        if (!std::isfinite(input.stagedCorrectionShiftMetres)) return false;
        if (input.stagedCorrectionShiftMetres != 0 &&
            (!input.removeStagedIntegerDatumShift ||
             !input.stagedIntegerDatumComparisonProven)) return false;
    }
    return true;
}

struct ZhangProductCommonModeAssessment
{
    bool assessed = false;
    double value = std::numeric_limits<double>::quiet_NaN();
    std::size_t distinctAnchors = 0;
    const char* reason = "INSUFFICIENT_COMPARABLE_COMMON_MODE_ANCHORS";
};

/** Three DISTINCT satellites are the minimum for detecting one 100 m
 * outlier with a median. Two samples (0,100) must not select 100 and let
 * the corrupted product cancel itself. This is a robust diagnostic, NOT a
 * statistical certificate that a majority of anchors is correct. */
inline ZhangProductCommonModeAssessment zhangEstimateProductCommonMode(
    const std::vector<std::pair<std::string, double>>& anchors)
{
    ZhangProductCommonModeAssessment out;
    std::map<std::string, double> distinct;
    for (const auto& [satellite, delta] : anchors)
    {
        if (satellite.empty() || !std::isfinite(delta))
        {
            out.reason = "INVALID_COMMON_MODE_ANCHOR";
            return out;
        }
        const auto [where, inserted] = distinct.emplace(satellite, delta);
        if (!inserted && where->second != delta)
        {
            out.reason = "CONFLICTING_DUPLICATE_COMMON_MODE_ANCHOR";
            return out;
        }
    }
    out.distinctAnchors = distinct.size();
    if (distinct.size() < 3) return out;
    std::vector<double> deltas;
    for (const auto& [satellite, delta] : distinct) deltas.push_back(delta);
    const auto middle = deltas.begin() + deltas.size() / 2;
    std::nth_element(deltas.begin(), middle, deltas.end());
    out.value = *middle;
    out.assessed = true;
    out.reason = "ROBUST_MEDIAN_VERIFIED_COMPARISON_DOMAINS";
    return out;
}

struct ZhangProductNumericalAssessment
{
    bool evaluated = false;
    bool valid = false;
    bool gapOnlyUnassessed = false;
    double residualStep = std::numeric_limits<double>::quiet_NaN();
    const char* state = "UNASSESSED";
    const char* reason = "NUMERICAL_REFERENCE_NOT_EVALUATED";
};

inline ZhangProductNumericalAssessment zhangAssessProductNumericalContinuity(
    const ZhangProductNumericalReference& reference,
    double commonMode,
    double residualLimit,
    bool commonModeAssessed = true,
    const char* unassessedCommonModeReason = "NO_COMPARABLE_COMMON_MODE_ANCHORS")
{
    using Kind = ZhangProductNumericalReferenceKind;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    if (reference.kind == Kind::INITIAL_BASELINE ||
        reference.kind == Kind::NEW_SEGMENT_BASELINE ||
        reference.kind == Kind::INCOMPARABLE_NEW_DATUM)
        return {false, true, false, missing, "BASELINE_INITIALIZATION", reference.reason};
    if (reference.kind != Kind::COMPARABLE)
        return {false, false, reference.kind == Kind::GAP_EXCEEDED, missing,
            reference.kind == Kind::INVALID ? "INVALID_INPUT" : "UNASSESSED",
            reference.reason};
    if (!commonModeAssessed)
        return {false, false, false, missing, "UNASSESSED",
            unassessedCommonModeReason};
    if (!reference.hasDelta || !std::isfinite(reference.delta) ||
        !std::isfinite(commonMode) || !std::isfinite(residualLimit) ||
        residualLimit <= 0)
        return {false, false, false, missing, "INVALID_INPUT",
            "NUMERICAL_COMPARISON_INPUT_INVALID"};
    const double residual = std::abs(reference.delta - commonMode);
    // The repaired PRODUCT_FIXED writer has an unchanged >=1 m uncertainty
    // guard. Unlike the legacy generic guard, zero/negative comparison limits
    // cannot silently disable this fixed-product evidence requirement.
    if (!std::isfinite(residual) || residual > residualLimit)
        return {true, false, false, residual, "EVALUATED_INVALID",
            "COMMON_MODE_REMOVED_STEP_EXCEEDED"};
    return {true, true, false, residual, "EVALUATED_VALID",
        "COMMON_MODE_REMOVED_STEP_WITHIN_LIMIT"};
}

/** A finite expired reference cannot certify continuity, but must not prevent
 * an independently valid current spatial product from starting a NEW segment.
 * This includes an already-required NEW segment whose old transport chain is
 * unavailable. Native/spatial failures and measured jumps remain suspended. */
inline ZhangProductPublicationDecision zhangProductPublicationWithNumericalAssessment(
    ZhangProductPublicationEvidence evidence,
    const ZhangProductNumericalAssessment& assessment)
{
    evidence.numericalStepValid = true; // Evaluate independent structural facts.
    const auto nominal = zhangProductPublicationGate(evidence);
    if (nominal.route == ZhangProductPublicationRoute::SUSPENDED)
        return nominal;
    if (assessment.valid)
    {
        // A missing reference must never be mistaken for a previous identity.
        if (!assessment.evaluated &&
            nominal.route == ZhangProductPublicationRoute::SAME_SEGMENT)
            return {ZhangProductPublicationRoute::SUSPENDED, false,
                "NUMERICAL_BASELINE_CANNOT_AUTHORIZE_SAME_SEGMENT"};
        return nominal;
    }
    if (!assessment.evaluated && assessment.gapOnlyUnassessed &&
        evidence.priorPublishedIdentity &&
        (nominal.route == ZhangProductPublicationRoute::SAME_SEGMENT ||
         nominal.route == ZhangProductPublicationRoute::NEW_SEGMENT))
        return {ZhangProductPublicationRoute::NEW_SEGMENT, true,
            "NUMERICAL_REFERENCE_GAP_UNASSESSED_NEW_SEGMENT"};
    return {ZhangProductPublicationRoute::SUSPENDED, false,
        assessment.evaluated ? "NUMERICAL_STEP_INVALID" : assessment.reason};
}

/** Bind current evidence explicitly instead of rereading row usability flags
 * which the temporal or dual-signal gate may already have changed. */
inline ZhangProductPublicationDecision zhangProductPublicationWithCurrentEvidence(
    const ZhangCurrentProductSpatialEvidence& current,
    ZhangProductPublicationEvidence temporal,
    const ZhangProductNumericalAssessment& assessment)
{
    temporal.nativeProductValid = current.nativeValid();
    temporal.currentSpatialCertificateValid = current.currentIntegerCertificateValid;
    return zhangProductPublicationWithNumericalAssessment(temporal, assessment);
}
