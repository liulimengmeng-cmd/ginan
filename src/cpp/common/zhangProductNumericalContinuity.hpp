#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "common/zhangProductPublicationGate.hpp"

/** Numerical evidence is separate from exact physical/chart transport.
 * In particular an expired reference is NOT a measured zero correction step.
 * This helper does not establish any integer or structural certificate. */
enum class ZhangProductNumericalReferenceKind
{
    INITIAL_BASELINE,
    NEW_SEGMENT_BASELINE,
    COMPARABLE,
    GAP_EXCEEDED,
    UNASSESSED,
    INVALID
};

struct ZhangProductNumericalReferenceInput
{
    bool currentNumericValid = false;
    bool hasPrevious = false;
    bool samePublishedIdentity = false;
    ZhangProductPublicationRoute nominalRoute =
        ZhangProductPublicationRoute::SUSPENDED;
    double currentCorrection = std::numeric_limits<double>::quiet_NaN();
    double previousCorrection = std::numeric_limits<double>::quiet_NaN();
    double currentSigma = std::numeric_limits<double>::quiet_NaN();
    double previousSigma = std::numeric_limits<double>::quiet_NaN();
    double gapSeconds = std::numeric_limits<double>::quiet_NaN();
    double maximumGapSeconds = std::numeric_limits<double>::quiet_NaN();
    // A proposed NEW public segment can explicitly reset its integer chart
    // constant. That deterministic writer operation is not a measured native
    // product step. Compare the same native chart before/after the operation,
    // but only when its complete FLOAT mean/Q map has independently passed.
    // These fields never certify the old integer bridge or SAME_SEGMENT.
    bool publicationGammaResetApplied = false;
    bool nativeFullMomentMapProven = false;
    long long gammaBeforePublicationReset = 0;
    long long gammaAfterPublicationReset = 0;
    double publicationWavelengthMetres = std::numeric_limits<double>::quiet_NaN();
    double appliedPhaseAdjustmentMetres = std::numeric_limits<double>::quiet_NaN();
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
    if (input.publicationGammaResetApplied)
    {
        if (input.nominalRoute != ZhangProductPublicationRoute::NEW_SEGMENT || !input.hasPrevious)
            return {Kind::INVALID, false, missing, "NUMERICAL_GAMMA_RESET_NOT_NEW_SEGMENT"};
        // An explicitly changed counter/datum already starts a new baseline;
        // it is not a comparison in the previous public chart. Do not impose
        // the old-chart FLOAT map on that pre-existing physical reset path.
        if (!input.samePublishedIdentity)
            return {Kind::NEW_SEGMENT_BASELINE, false, missing,
                "EXPLICIT_NEW_SEGMENT_NUMERICAL_BASELINE"};
        if (!input.nativeFullMomentMapProven)
            return {Kind::INVALID, false, missing, "NUMERICAL_GAMMA_RESET_NATIVE_MAP_UNPROVEN"};
        // Restrict compensation to the actual old-gamma -> zero operation.
        // Never subtract an arbitrary fitted jump or a newly searched SAME
        // segment gamma, which could conceal a physical/product anomaly.
        constexpr long long exactDoubleInteger = 9007199254740991LL;
        if (input.gammaAfterPublicationReset != 0 || input.gammaBeforePublicationReset == 0 ||
            input.gammaBeforePublicationReset < -exactDoubleInteger ||
            input.gammaBeforePublicationReset > exactDoubleInteger ||
            !std::isfinite(input.publicationWavelengthMetres) || input.publicationWavelengthMetres <= 0 ||
            !std::isfinite(input.appliedPhaseAdjustmentMetres))
            return {Kind::INVALID, false, missing, "NUMERICAL_GAMMA_RESET_IDENTITY_INVALID"};
        const double expectedPhaseAdjustment = input.publicationWavelengthMetres *
            (static_cast<double>(input.gammaAfterPublicationReset) -
                static_cast<double>(input.gammaBeforePublicationReset));
        // The writer stores this exact expression, not an independently
        // estimated value. A relative tolerance would admit large arbitrary
        // offsets when gamma is large, so require the actual operation exactly.
        if (!std::isfinite(expectedPhaseAdjustment) ||
            expectedPhaseAdjustment != input.appliedPhaseAdjustmentMetres)
            return {Kind::INVALID, false, missing, "NUMERICAL_GAMMA_RESET_APPLIED_SHIFT_MISMATCH"};
        // q = C-B; undo ONLY the explicit candidate B chart shift for the
        // comparison. The emitted q/B/gamma remain in their proposed NEW chart.
        comparisonCorrection += input.appliedPhaseAdjustmentMetres;
        if (!std::isfinite(comparisonCorrection))
            return {Kind::INVALID, false, missing, "NUMERICAL_GAMMA_RESET_COMPARISON_NONFINITE"};
    }
    if (!input.hasPrevious)
        return {Kind::INITIAL_BASELINE, false, missing, "INITIAL_NUMERICAL_BASELINE"};
    // An already-required physical/publication segment change starts a new
    // baseline. It cannot claim numerical continuity with a different datum.
    if (!input.samePublishedIdentity)
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
    const double delta = comparisonCorrection - input.previousCorrection;
    if (!std::isfinite(delta))
        return {Kind::INVALID, false, missing, "NUMERICAL_REFERENCE_DELTA_NONFINITE"};
    if (!std::isfinite(input.gapSeconds) || input.gapSeconds <= 0)
        return {Kind::UNASSESSED, false, missing, "NUMERICAL_REFERENCE_TIME_INVALID"};
    if (input.gapSeconds > input.maximumGapSeconds)
        return {Kind::GAP_EXCEEDED, false, missing, "NUMERICAL_REFERENCE_GAP_EXCEEDED"};
    return {Kind::COMPARABLE, true, delta, "COMPARABLE_NUMERICAL_REFERENCE"};
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
    double residualLimit)
{
    using Kind = ZhangProductNumericalReferenceKind;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    if (reference.kind == Kind::INITIAL_BASELINE ||
        reference.kind == Kind::NEW_SEGMENT_BASELINE)
        return {false, true, false, missing, "BASELINE_INITIALIZATION", reference.reason};
    if (reference.kind != Kind::COMPARABLE)
        return {false, false, reference.kind == Kind::GAP_EXCEEDED, missing,
            reference.kind == Kind::INVALID ? "INVALID_INPUT" : "UNASSESSED",
            reference.reason};
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

/** An expired finite reference may start a clearly labelled NEW segment only
 * when every original current-epoch structural/spatial gate would otherwise
 * prove SAME_SEGMENT. It is never presented as a successful comparison.
 * All measured failures and all other unknown evidence remain suspended. */
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
        nominal.route == ZhangProductPublicationRoute::SAME_SEGMENT)
        return {ZhangProductPublicationRoute::NEW_SEGMENT, true,
            "NUMERICAL_REFERENCE_GAP_UNASSESSED_NEW_SEGMENT"};
    return {ZhangProductPublicationRoute::SUSPENDED, false,
        assessment.evaluated ? "NUMERICAL_STEP_INVALID" : assessment.reason};
}
