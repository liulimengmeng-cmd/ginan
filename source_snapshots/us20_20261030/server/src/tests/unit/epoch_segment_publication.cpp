#include "common/zhangProductNumericalContinuity.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
int checks = 0;
void require(bool result, const char* message)
{
    ++checks;
    if (!result) throw std::runtime_error(message);
}

ZhangCurrentProductSpatialEvidence current()
{ return {true, true, true}; }

ZhangProductPublicationEvidence temporal()
{
    ZhangProductPublicationEvidence out;
    out.nativeProductValid = true;
    out.currentSpatialCertificateValid = true;
    out.priorPublishedIdentity = true;
    out.samePhysicalSegment = true;
    out.sameProductFunctional = true;
    out.sameBackendChart = true;
    out.requireFullJointHandoff = true;
    out.fullJointHandoffProven = true;
    out.structuralProductMapProven = true;
    out.numericalStepValid = true;
    return out;
}

ZhangProductNumericalReferenceInput endpoint()
{
    ZhangProductNumericalReferenceInput in;
    in.currentNumericValid = true;
    in.hasPrevious = true;
    in.samePublishedIdentity = true;
    in.comparisonDomain = ZhangProductComparisonDomain::SAME_VERIFIED_DOMAIN;
    in.nominalRoute = ZhangProductPublicationRoute::SAME_SEGMENT;
    in.currentCorrection = 10.02;
    in.previousCorrection = 10;
    in.currentSigma = in.previousSigma = 0.01;
    in.gapSeconds = 30;
    in.maximumGapSeconds = 120;
    return in;
}
}

int main()
{
    using Route = ZhangProductPublicationRoute;
    using Domain = ZhangProductComparisonDomain;
    using Kind = ZhangProductNumericalReferenceKind;
    auto in = endpoint();
    auto evidence = temporal();
    auto assess = [](const auto& input, bool anchors = true)
    {
        return zhangAssessProductNumericalContinuity(
            zhangProductNumericalReference(input), 0, 1, anchors);
    };
    require(zhangProductPublicationWithCurrentEvidence(current(), evidence, assess(in)).route ==
        Route::SAME_SEGMENT, "ordinary comparable product lost continuity");

    // Public counters have not advanced yet. The target/chart identity, not
    // those staged counters, makes this candidate a new numerical baseline.
    in.comparisonDomain = Domain::DIFFERENT_CURRENT_DOMAIN;
    in.nominalRoute = Route::NEW_SEGMENT;
    in.currentCorrection += 100;
    evidence.fullJointHandoffProven = false;
    evidence.sameBackendChart = false;
    const auto changedReference = zhangProductNumericalReference(in);
    const auto changed = assess(in, false);
    require(changedReference.kind == Kind::INCOMPARABLE_NEW_DATUM,
        "staged old counter overrode the changed target/chart");
    require(!changedReference.hasDelta && std::isnan(changedReference.delta),
        "incomparable target fabricated a delta");
    require(changed.valid && !changed.evaluated && std::isnan(changed.residualStep),
        "new baseline claimed a measured zero residual");
    const auto opened = zhangProductPublicationWithCurrentEvidence(current(), evidence, changed);
    require(opened.route == Route::NEW_SEGMENT && opened.requireUserReset,
        "valid current product was blocked only by missing old transport");

    // Current evidence is independent and cannot be manufactured by temporal
    // eligibility, old mutable flags, or permission to start a new segment.
    for (int missing = 0; missing < 3; ++missing)
    {
        auto bad = current();
        if (missing == 0) bad.nativeProjectionValid = false;
        if (missing == 1) bad.currentIntegerCertificateValid = false;
        if (missing == 2) bad.currentJointMomentsValid = false;
        require(zhangProductPublicationWithCurrentEvidence(bad, evidence, changed).route ==
            Route::SUSPENDED, "new segment bypassed current evidence");
    }
    auto staleFlags = evidence;
    staleFlags.nativeProductValid = staleFlags.currentSpatialCertificateValid = false;
    require(zhangProductPublicationWithCurrentEvidence(current(), staleFlags, changed).route ==
        Route::NEW_SEGMENT, "mutable temporal flags contaminated frozen current evidence");

    // A missing bridge alone does not change the comparison domain. A genuine
    // 100 m same-domain corruption MUST still fail even with nominal NEW.
    in = endpoint();
    in.nominalRoute = Route::NEW_SEGMENT;
    in.currentCorrection += 100;
    const auto corruption = assess(in);
    require(corruption.evaluated && !corruption.valid && corruption.residualStep > 100,
        "100 m same-domain corruption escaped the numerical test");
    require(zhangProductPublicationWithCurrentEvidence(current(), evidence, corruption).route ==
        Route::SUSPENDED, "NEW_SEGMENT relabelled a measured numerical failure");
    in.comparisonDomain = Domain::UNASSESSED;
    require(zhangProductPublicationWithCurrentEvidence(current(), evidence, assess(in)).route ==
        Route::SUSPENDED, "unknown domain was misrepresented as a known new datum");

    in = endpoint();
    const auto noAnchors = assess(in, false);
    require(!noAnchors.valid && !noAnchors.evaluated && std::isnan(noAnchors.residualStep),
        "missing common-mode anchors were silently replaced by a measured zero");
    require(std::string(noAnchors.reason) == "NO_COMPARABLE_COMMON_MODE_ANCHORS",
        "missing anchors lost their diagnostic reason");
    require(zhangProductPublicationWithCurrentEvidence(current(), temporal(), noAnchors).route ==
        Route::SUSPENDED, "unassessed comparable product was called continuous");

    // One changed signal is legal dual-frequency AR; an actually invalid
    // partner is not. The reason string itself is not a rejection status.
    auto pair = zhangAggregateDualFrequencyPublication(opened,
        {Route::SAME_SEGMENT, false, "FULL_JOINT_STRUCTURAL_IDENTITY_NO_AR_DECISION"});
    require(pair.dualFrequencyArUsable && pair.resetBothUserAmbiguities,
        "NEW plus SAME incorrectly removed current dual AR");
    pair = zhangAggregateDualFrequencyPublication(opened,
        {Route::SUSPENDED, false, "NUMERICAL_STEP_INVALID"});
    require(!pair.dualFrequencyArUsable,
        "an invalid signal authorized its dual-frequency partner");

    // Nonfinite CURRENT products cannot be salvaged as incomparable baselines.
    in = endpoint();
    in.nominalRoute = Route::NEW_SEGMENT;
    in.comparisonDomain = Domain::DIFFERENT_CURRENT_DOMAIN;
    in.currentCorrection = std::numeric_limits<double>::quiet_NaN();
    require(!assess(in).valid, "nonfinite current endpoint entered new baseline");

    // A missing old integer bridge does not make an unchanged numerical
    // domain incomparable. Three current NEW candidates can anchor a robust
    // common mode, while an individual 100 m corruption remains rejected.
    auto newTemporal = temporal();
    newTemporal.fullJointHandoffProven = false;
    std::vector<ZhangProductNumericalReferenceInput> comparableNew;
    std::vector<std::pair<std::string, double>> anchors;
    for (int i = 0; i < 3; ++i)
    {
        auto candidate = endpoint();
        candidate.nominalRoute = Route::NEW_SEGMENT;
        candidate.stagedCorrectionShiftMetres = 0;
        candidate.currentCorrection = candidate.previousCorrection + (i == 2 ? 100 : 0);
        const auto reference = zhangProductNumericalReference(candidate);
        require(zhangProductCommonModeAnchorEligible(current(), candidate, reference),
            "same-domain NEW product with current evidence could not anchor");
        comparableNew.push_back(candidate);
        anchors.emplace_back("G0" + std::to_string(i + 1), reference.delta);
    }
    const auto robust = zhangEstimateProductCommonMode(anchors);
    require(robust.assessed && robust.distinctAnchors == 3 && robust.value == 0,
        "median of distinct (0,0,100) anchors was not zero");
    for (int i = 0; i < 3; ++i)
    {
        const auto comparison = zhangAssessProductNumericalContinuity(
            zhangProductNumericalReference(comparableNew[i]), robust.value, 1,
            robust.assessed, robust.reason);
        const auto result = zhangProductPublicationWithCurrentEvidence(
            current(), newTemporal, comparison);
        require(result.route == (i == 2 ? Route::SUSPENDED : Route::NEW_SEGMENT),
            "NEW anchor median accepted the bad star or refused the good stars");
    }
    std::reverse(anchors.begin(), anchors.end());
    require(zhangEstimateProductCommonMode(anchors).value == robust.value,
        "common mode changed with insertion order");
    const auto two = zhangEstimateProductCommonMode({{"G01", 0}, {"G03", 100}});
    require(!two.assessed && two.distinctAnchors == 2 && std::isnan(two.value),
        "two anchors selected the corrupted star's delta");
    const auto twoAssessment = zhangAssessProductNumericalContinuity(
        zhangProductNumericalReference(comparableNew[2]), two.value, 1,
        two.assessed, two.reason);
    require(!twoAssessment.evaluated && !twoAssessment.valid &&
        std::isnan(twoAssessment.residualStep) &&
        std::string(twoAssessment.reason) == "INSUFFICIENT_COMPARABLE_COMMON_MODE_ANCHORS",
        "insufficient anchors were silently treated as a measured zero common mode");
    require(zhangProductPublicationWithCurrentEvidence(current(), newTemporal, twoAssessment).route ==
        Route::SUSPENDED, "two anchors bypassed the numerical gate via NEW");
    const auto duplicates = zhangEstimateProductCommonMode({{"G01", 0}, {"G01", 0}, {"G02", 100}});
    require(!duplicates.assessed && duplicates.distinctAnchors == 2,
        "duplicate satellite rows manufactured the third independent anchor");
    require(!zhangEstimateProductCommonMode({{"G01", 0}, {"G01", 100}, {"G02", 0}}).assessed,
        "conflicting duplicate anchor values were accepted");

    auto normalized = endpoint();
    normalized.nominalRoute = Route::NEW_SEGMENT;
    normalized.stagedCorrectionShiftMetres = 5 * 0.190293672798365;
    normalized.correctionBeforeStagedDatumChange = normalized.currentCorrection;
    normalized.currentCorrection += normalized.stagedCorrectionShiftMetres;
    require(!zhangProductCommonModeAnchorEligible(current(), normalized,
        zhangProductNumericalReference(normalized)), "unproved nonzero staged gamma became an anchor");
    normalized.removeStagedIntegerDatumShift = true;
    normalized.stagedIntegerDatumComparisonProven = true;
    auto normalizedReference = zhangProductNumericalReference(normalized);
    require(zhangProductCommonModeAnchorEligible(current(), normalized, normalizedReference) &&
        std::abs(normalizedReference.delta - 0.02) < 1e-12,
        "proved staged gamma normalization did not use the pre-gamma delta");
    auto noCurrentInteger = current(); noCurrentInteger.currentIntegerCertificateValid = false;
    require(!zhangProductCommonModeAnchorEligible(noCurrentInteger, normalized, normalizedReference),
        "missing CURRENT integer certificate became a numerical anchor");
    normalized.comparisonDomain = Domain::DIFFERENT_CURRENT_DOMAIN;
    require(!zhangProductCommonModeAnchorEligible(current(), normalized,
        zhangProductNumericalReference(normalized)), "different-domain NEW product polluted the common mode");
    std::cout << "epoch_segment_publication checks=" << checks << " PASS\n";
}
