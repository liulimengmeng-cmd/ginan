#include "common/zhangProductNumericalContinuity.hpp"
#include "common/zhangPublishedContinuity.hpp"

#include <iostream>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
int checks = 0;
void require(bool value, const char* message)
{
    ++checks;
    if (!value) throw std::runtime_error(message);
}

ZhangProductPublicationEvidence unchanged()
{
    ZhangProductPublicationEvidence evidence;
    evidence.nativeProductValid = true;
    evidence.currentSpatialCertificateValid = true;
    evidence.priorPublishedIdentity = true;
    evidence.samePhysicalSegment = true;
    evidence.sameProductFunctional = true;
    evidence.sameBackendChart = true;
    evidence.exactSameTimeMapProven = true;
    evidence.conditionalProductMapProven = true;
    evidence.currentEpochBridgeProven = true;
    evidence.requireFullJointHandoff = true;
    evidence.fullJointHandoffProven = true;
    evidence.structuralProductMapProven = true;
    evidence.numericalStepValid = true;
    return evidence;
}

ZhangProductNumericalReferenceInput comparable(double gap = 30)
{
    ZhangProductNumericalReferenceInput input;
    input.currentNumericValid = true;
    input.hasPrevious = true;
    input.samePublishedIdentity = true;
    input.nominalRoute = ZhangProductPublicationRoute::SAME_SEGMENT;
    input.currentCorrection = 100.02;
    input.previousCorrection = 100;
    input.currentSigma = input.previousSigma = 0.01;
    input.gapSeconds = gap;
    input.maximumGapSeconds = 120;
    return input;
}

ZhangProductNumericalAssessment assess(
    const ZhangProductNumericalReferenceInput& input, double common = 0, double limit = 1)
{
    return zhangAssessProductNumericalContinuity(
        zhangProductNumericalReference(input), common, limit);
}

void mustSuspend(const ZhangProductNumericalReferenceInput& input, const char* message)
{
    const auto result = assess(input);
    const auto decision = zhangProductPublicationWithNumericalAssessment(unchanged(), result);
    require(!result.valid && !result.gapOnlyUnassessed, message);
    require(decision.route == ZhangProductPublicationRoute::SUSPENDED, message);
    require(!decision.requireUserReset, "invalid comparison must not relabel a new segment");
}
}

int main()
{
    using Route = ZhangProductPublicationRoute;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    auto initial = comparable();
    initial.hasPrevious = false;
    initial.samePublishedIdentity = false;
    initial.previousCorrection = initial.previousSigma = initial.gapSeconds = nan;
    initial.nominalRoute = Route::NEW_SEGMENT;
    auto initialEvidence = unchanged();
    initialEvidence.priorPublishedIdentity = false;
    const auto initialAssessment = assess(initial);
    const auto first = zhangProductPublicationWithNumericalAssessment(
        initialEvidence, initialAssessment);
    require(initialAssessment.valid && !initialAssessment.evaluated,
        "finite first product must initialize without fictitious historical comparison");
    require(std::isnan(initialAssessment.residualStep), "first residual must be unavailable, not zero");
    require(first.route == Route::NEW_SEGMENT && !first.requireUserReset,
        "first product must retain INITIAL_PRODUCT_SEGMENT with no user reset");
    require(std::string(first.reason) == "INITIAL_PRODUCT_SEGMENT", "first product reason changed");
    require(zhangProductPublicationWithNumericalAssessment(unchanged(), initialAssessment).route ==
        Route::SUSPENDED, "missing reference must not assert an existing SAME_SEGMENT");

    auto atLimit = comparable(120);
    const auto atLimitReference = zhangProductNumericalReference(atLimit);
    const auto atLimitAssessment = assess(atLimit);
    require(atLimitReference.hasDelta && atLimitAssessment.evaluated && atLimitAssessment.valid,
        "120 second boundary must actually compare the finite correction");
    require(std::abs(atLimitAssessment.residualStep - 0.02) < 1e-12,
        "120 second comparison changed the direct arithmetic oracle");
    require(zhangProductPublicationWithNumericalAssessment(unchanged(), atLimitAssessment).route ==
        Route::SAME_SEGMENT, "valid boundary comparison should remain same segment");

    auto expired = comparable(150);
    // Real E11 01:01:00 -> 01:03:30 endpoints.  B moved by about 8 mm while
    // C-B moved by about -5.118 m; this fixture is not a physical slip oracle.
    expired.previousCorrection = 703808.87845273 - 22.23767320095;
    expired.currentCorrection = 703803.76813968 - 22.24578612355;
    const auto expiredReference = zhangProductNumericalReference(expired);
    const auto expiredAssessment = assess(expired, -0.02298464099294506);
    const auto recovered = zhangProductPublicationWithNumericalAssessment(
        unchanged(), expiredAssessment);
    require(!expiredReference.hasDelta && std::isnan(expiredReference.delta),
        "expired delta must remain explicitly unavailable");
    require(!expiredAssessment.evaluated && !expiredAssessment.valid &&
        expiredAssessment.gapOnlyUnassessed && std::isnan(expiredAssessment.residualStep),
        "150 second gap must not produce zero or a valid numerical comparison");
    require(recovered.route == Route::NEW_SEGMENT && recovered.requireUserReset,
        "150 second reference may recover only as explicit new segment plus reset");
    require(std::string(recovered.reason) == "NUMERICAL_REFERENCE_GAP_UNASSESSED_NEW_SEGMENT",
        "gap recovery must carry an explicit unassessed reason");
    const auto recoveredPair = zhangAggregateDualFrequencyPublication(recovered, recovered);
    require(recoveredPair.dualFrequencyArUsable && recoveredPair.resetBothUserAmbiguities,
        "two certified NEW_SEGMENT rows must remain usable with a dual ambiguity reset");
    // Current integer authorization and Qp do not require a claim of transport
    // from the old public segment. UNKNOWN/nan is honest for this new datum.
    ZhangPublishedContinuityContract newContract;
    newContract.bundleId = "1405213410:PRODUCT_FIXED";
    newContract.previousEpochGpstSeconds = 1405213260;
    newContract.previousCounter = newContract.previousDatumVersion = 0;
    newContract.previousPublicSegment = "E11-old";
    newContract.jointCovarianceComplete = true;
    require(zhangPublishedContinuityValid(newContract, 1405213410, "E11-new"),
        "new segment with current Qp must not require fictitious temporal identity");
    require(!zhangAuthorizePublishedIdentity(newContract, 1405213410, "E11-new",
        false, true, true, "E11-current-proof"),
        "new segment incorrectly manufactured a same-public-chart certificate");
    require(newContract.kind == "UNKNOWN" && std::isnan(newContract.floatCorrectionShiftMetres),
        "new-segment certificate fabricated a zero historical transport");

    // The same measured endpoint change at the inclusive 120 s boundary must
    // be rejected. NEW_SEGMENT is not an escape hatch from a measured failure.
    expired.gapSeconds = 120;
    const auto measuredFailure = assess(expired, -0.02298464099294506);
    require(measuredFailure.evaluated && !measuredFailure.valid &&
        measuredFailure.residualStep > 5, "large measured E11 correction change was not rejected");
    require(zhangProductPublicationWithNumericalAssessment(unchanged(), measuredFailure).route ==
        Route::SUSPENDED, "measured failure must not be relabelled NEW_SEGMENT");
    auto nominalNew = unchanged();
    nominalNew.fullJointHandoffProven = false;
    require(zhangProductPublicationGate(nominalNew).route == Route::NEW_SEGMENT,
        "missing handoff fixture should nominally open a segment");
    require(zhangProductPublicationWithNumericalAssessment(nominalNew, measuredFailure).route ==
        Route::SUSPENDED, "nominal new segment bypassed a measured failure");
    require(zhangProductPublicationWithNumericalAssessment(nominalNew, expiredAssessment).route ==
        Route::SUSPENDED, "gap alone must not bypass missing full joint handoff evidence");
    auto noNative = unchanged(); noNative.nativeProductValid = false;
    require(zhangProductPublicationWithNumericalAssessment(noNative, expiredAssessment).route ==
        Route::SUSPENDED, "gap recovery bypassed native numeric/physical gate");
    auto noInteger = unchanged(); noInteger.currentSpatialCertificateValid = false;
    require(zhangProductPublicationWithNumericalAssessment(noInteger, expiredAssessment).route ==
        Route::SUSPENDED, "gap recovery bypassed current integer certificate");

    // Simulate only the existing post-COMMIT history update: a new segment at
    // t=150 becomes the finite reference. Subsequent 30 s epochs must not
    // remain suspended or recursively open a new segment on every call.
    auto next = comparable();
    next.previousCorrection = expired.currentCorrection;
    int resets = recovered.requireUserReset ? 1 : 0;
    for (int epoch = 0; epoch < 10; ++epoch)
    {
        next.currentCorrection = next.previousCorrection + 0.01;
        const auto decision = zhangProductPublicationWithNumericalAssessment(unchanged(), assess(next));
        require(decision.route == Route::SAME_SEGMENT && !decision.requireUserReset,
            "committed gap recovery did not return to tested same-segment publication");
        resets += decision.requireUserReset ? 1 : 0;
        next.previousCorrection = next.currentCorrection;
    }
    require(resets == 1, "recovery repeatedly reset instead of advancing committed reference");

    for (double invalid : {nan, inf, -inf})
    {
        auto input = comparable(); input.currentCorrection = invalid;
        mustSuspend(input, "nonfinite current correction published");
        input = comparable(150); input.currentSigma = invalid;
        mustSuspend(input, "nonfinite current sigma bypassed gap recovery");
        input = comparable(150); input.previousCorrection = invalid;
        mustSuspend(input, "nonfinite historical correction bypassed gap recovery");
        input = comparable(150); input.previousSigma = invalid;
        mustSuspend(input, "nonfinite historical sigma bypassed gap recovery");
        input = comparable(); input.gapSeconds = invalid;
        mustSuspend(input, "nonfinite elapsed time published");
        input = comparable(); input.maximumGapSeconds = invalid;
        mustSuspend(input, "nonfinite maximum gap published");
        auto failed = assess(comparable(), invalid);
        require(!failed.valid && zhangProductPublicationWithNumericalAssessment(
            unchanged(), failed).route == Route::SUSPENDED, "nonfinite common mode published");
        failed = assess(comparable(), 0, invalid);
        require(!failed.valid && zhangProductPublicationWithNumericalAssessment(
            unchanged(), failed).route == Route::SUSPENDED, "nonfinite residual limit published");
    }
    for (double invalidTime : {0.0, -30.0})
    {
        auto input = comparable(invalidTime);
        mustSuspend(input, "nonpositive reference time published");
    }
    auto invalidCurrent = comparable(150); invalidCurrent.currentNumericValid = false;
    mustSuspend(invalidCurrent, "current numeric_valid=false bypassed gap recovery");
    auto negativeSigma = comparable(150); negativeSigma.currentSigma = -0.01;
    mustSuspend(negativeSigma, "negative current sigma bypassed gap recovery");
    negativeSigma = comparable(150); negativeSigma.previousSigma = -0.01;
    mustSuspend(negativeSigma, "negative historical sigma bypassed gap recovery");
    for (double invalidLimit : {0.0, -1.0})
    {
        const auto result = assess(comparable(), 0, invalidLimit);
        require(!result.valid && !result.evaluated &&
            zhangProductPublicationWithNumericalAssessment(unchanged(), result).route ==
                Route::SUSPENDED, "nonpositive fixed-product limit silently disabled the gate");
    }
    auto mismatched = comparable(150); mismatched.samePublishedIdentity = false;
    mustSuspend(mismatched, "identity mismatch mislabeled as gap-only recovery");
    auto overflow = comparable();
    overflow.currentCorrection = std::numeric_limits<double>::max();
    overflow.previousCorrection = -std::numeric_limits<double>::max();
    mustSuspend(overflow, "overflowing finite-endpoint delta published");
    overflow.gapSeconds = 150;
    mustSuspend(overflow, "expired reference relabelled arithmetic overflow as gap-only evidence");
    overflow = comparable(150);
    overflow.currentSigma = overflow.previousSigma = std::numeric_limits<double>::max();
    mustSuspend(overflow, "overflowing endpoint sigma sum published");
    overflow = comparable(150);
    overflow.currentSigma = std::numeric_limits<double>::max() / 5;
    mustSuspend(overflow, "overflowing ten-sigma uncertainty guard published");
    overflow = comparable();
    overflow.currentCorrection = std::numeric_limits<double>::max();
    overflow.previousCorrection = 0;
    const auto overflowResidual = assess(overflow, -std::numeric_limits<double>::max());
    require(!overflowResidual.valid && overflowResidual.evaluated &&
        zhangProductPublicationWithNumericalAssessment(unchanged(), overflowResidual).route ==
            Route::SUSPENDED, "finite delta/common-mode residual overflow published");
    auto equality = comparable(); equality.currentCorrection = 101;
    require(assess(equality).valid, "residual limit equality changed from inclusive to exclusive");
    equality.currentCorrection = std::nextafter(101.0, inf);
    require(!assess(equality).valid, "residual above unchanged limit was accepted");
    auto explicitSegment = comparable();
    explicitSegment.nominalRoute = Route::NEW_SEGMENT;
    explicitSegment.samePublishedIdentity = false;
    require(zhangProductPublicationWithNumericalAssessment(nominalNew,
        assess(explicitSegment)).route == Route::NEW_SEGMENT,
        "pre-existing explicit datum change must retain new-baseline semantics");

    // Independent scalar oracle for the observed candidate gamma-reset bug.
    // q_candidate = q_native - lambda*gamma_new, whereas the preceding public
    // q used gamma_old. Only the explicit old->zero operation is removed from
    // the NUMERICAL comparison, never from the emitted candidate or AR proof.
    const std::array<std::array<long long, 2>, 4> observedGammas{{
        {{31, 61}}, {{88, 125}}, {{44, 275}}, {{44, 275}}
    }}; // C39, C36, C27, C33; BDS B1C/B2a, not GPS L2W.
    const std::array<double, 2> lambdas{{
        299792458.0 / 1575420000.0, 299792458.0 / 1176450000.0
    }};
    for (const auto& pair : observedGammas)
    {
        std::array<ZhangProductPublicationDecision, 2> decisions;
        for (int signal = 0; signal < 2; ++signal)
        {
            auto reset = comparable();
            reset.nominalRoute = Route::NEW_SEGMENT;
            reset.gammaBeforePublicationReset = pair[signal];
            reset.gammaAfterPublicationReset = 0;
            reset.publicationWavelengthMetres = lambdas[signal];
            reset.appliedPhaseAdjustmentMetres = -lambdas[signal] * pair[signal];
            const double nativeCorrectionStep = 0.02;
            reset.currentCorrection = reset.previousCorrection + nativeCorrectionStep -
                reset.appliedPhaseAdjustmentMetres;
            const double candidateCorrection = reset.currentCorrection;
            require(!assess(reset).valid, "raw new-chart jump fixture no longer exposes the old bug");
            reset.publicationGammaResetApplied = true;
            reset.nativeFullMomentMapProven = true;
            const auto fixedAssessment = assess(reset);
            require(fixedAssessment.evaluated && fixedAssessment.valid &&
                std::abs(fixedAssessment.residualStep - nativeCorrectionStep) < 1e-12,
                "known gamma reset did not preserve the independently calculated native step");
            decisions[signal] = zhangProductPublicationWithNumericalAssessment(nominalNew, fixedAssessment);
            require(decisions[signal].route == Route::NEW_SEGMENT && decisions[signal].requireUserReset,
                "gamma reset correction fabricated old-segment integer continuity");
            require(reset.currentCorrection == candidateCorrection && reset.gammaAfterPublicationReset == 0,
                "numerical comparison changed the candidate publication chart");
            auto bad = reset;
            bad.currentCorrection += 5;
            const auto realJump = assess(bad);
            require(realJump.evaluated && !realJump.valid && realJump.residualStep > 5,
                "gamma correction absorbed a real extra native product anomaly");
            require(zhangProductPublicationWithNumericalAssessment(nominalNew, realJump).route == Route::SUSPENDED,
                "real native anomaly bypassed the unchanged numerical limit");
            bad = reset; bad.nativeFullMomentMapProven = false;
            require(!assess(bad).valid && !assess(bad).evaluated,
                "unknown native mean/Q map authorized gamma compensation");
            bad = reset; bad.nominalRoute = Route::SAME_SEGMENT;
            require(!assess(bad).valid, "SAME_SEGMENT silently removed a searched bridge gamma");
            bad = reset; bad.gammaAfterPublicationReset = 1;
            require(!assess(bad).valid, "arbitrary nonzero new gamma was compensated");
            bad = reset; bad.appliedPhaseAdjustmentMetres += 5;
            require(!assess(bad).valid, "a fitted arbitrary phase offset was accepted as the writer operation");
            bad = reset; bad.gapSeconds = 150;
            const auto gapReset = assess(bad);
            require(!gapReset.valid && gapReset.gapOnlyUnassessed &&
                zhangProductPublicationWithNumericalAssessment(nominalNew, gapReset).route == Route::SUSPENDED,
                "gamma reset manufactured a numerical comparison across an expired reference");
            auto absentNative = nominalNew; absentNative.nativeProductValid = false;
            auto absentSpatial = nominalNew; absentSpatial.currentSpatialCertificateValid = false;
            require(zhangProductPublicationWithNumericalAssessment(absentNative, fixedAssessment).route == Route::SUSPENDED &&
                zhangProductPublicationWithNumericalAssessment(absentSpatial, fixedAssessment).route == Route::SUSPENDED,
                "comparison compensation bypassed native or current spatial certificate gates");
        }
        const auto dualNew = zhangAggregateDualFrequencyPublication(decisions[0], decisions[1]);
        require(dualNew.dualFrequencyArUsable && dualNew.resetBothUserAmbiguities,
            "two valid rebaselined signals did not request a dual user reset");
        decisions[1] = {Route::SUSPENDED, false, "UNKNOWN_NATIVE_MAP"};
        require(!zhangAggregateDualFrequencyPublication(decisions[0], decisions[1]).dualFrequencyArUsable,
            "one valid rebaselined signal authorized an invalid dual-frequency partner");
    }
    auto resetIdentity = comparable();
    resetIdentity.nominalRoute = Route::NEW_SEGMENT;
    resetIdentity.publicationGammaResetApplied = resetIdentity.nativeFullMomentMapProven = true;
    resetIdentity.gammaBeforePublicationReset = -31;
    resetIdentity.gammaAfterPublicationReset = 0;
    resetIdentity.publicationWavelengthMetres = lambdas[0];
    resetIdentity.appliedPhaseAdjustmentMetres = 31*lambdas[0];
    resetIdentity.currentCorrection = resetIdentity.previousCorrection + 0.02 - resetIdentity.appliedPhaseAdjustmentMetres;
    require(assess(resetIdentity).valid, "negative old gamma sign was reversed");
    for (const auto largeGamma : {1000000000000LL, 9007199254740991LL})
    {
        auto largeReset = resetIdentity;
        largeReset.gammaBeforePublicationReset = largeGamma;
        largeReset.publicationWavelengthMetres = 1;
        largeReset.appliedPhaseAdjustmentMetres = -static_cast<double>(largeGamma);
        largeReset.currentCorrection = largeReset.previousCorrection - largeReset.appliedPhaseAdjustmentMetres;
        require(assess(largeReset).valid, "exact large-gamma writer operation was rejected");
        for (double forgedExtra : {16.0, 1024.0})
        {
            auto forged = largeReset;
            forged.appliedPhaseAdjustmentMetres += forgedExtra;
            forged.currentCorrection = forged.previousCorrection - forged.appliedPhaseAdjustmentMetres;
            const auto reference = zhangProductNumericalReference(forged);
            require(reference.kind == ZhangProductNumericalReferenceKind::INVALID &&
                std::string(reference.reason) == "NUMERICAL_GAMMA_RESET_APPLIED_SHIFT_MISMATCH",
                "large gamma relative tolerance concealed a forged non-writer shift");
            require(!assess(forged).valid, "large gamma forged shift acquired numerical authorization");
        }
    }
    auto physicalReset = resetIdentity;
    physicalReset.samePublishedIdentity = false;
    physicalReset.nativeFullMomentMapProven = false;
    physicalReset.publicationGammaResetApplied = false;
    physicalReset.previousCorrection = physicalReset.previousSigma = physicalReset.gapSeconds = nan;
    const auto physicalBaseline = assess(physicalReset);
    require(physicalBaseline.valid && !physicalBaseline.evaluated &&
        std::isnan(physicalBaseline.residualStep),
        "explicitly changed physical counter/datum incorrectly required an old-chart map");
    auto physicalEvidence = unchanged();
    physicalEvidence.samePhysicalSegment = false;
    physicalEvidence.fullJointHandoffProven = false;
    const auto physicalDecision = zhangProductPublicationWithNumericalAssessment(
        physicalEvidence, physicalBaseline);
    require(physicalDecision.route == Route::NEW_SEGMENT && physicalDecision.requireUserReset,
        "explicit physical reset did not preserve its honest new-segment/user-reset route");
    physicalReset.publicationGammaResetApplied = true;
    const auto flaggedPhysicalBaseline = assess(physicalReset);
    require(flaggedPhysicalBaseline.valid && !flaggedPhysicalBaseline.evaluated &&
        std::isnan(flaggedPhysicalBaseline.residualStep),
        "stale reset-operation metadata blocked an already changed physical identity");
    require(zhangProductPublicationWithNumericalAssessment(unchanged(), flaggedPhysicalBaseline).route ==
        Route::SUSPENDED, "new physical baseline falsely authorized SAME_SEGMENT");
    auto absentPhysicalNative = physicalEvidence;
    absentPhysicalNative.nativeProductValid = false;
    require(zhangProductPublicationWithNumericalAssessment(absentPhysicalNative, physicalBaseline).route ==
        Route::SUSPENDED, "physical new baseline bypassed the native product gate");
    auto invalidPhysicalReset = physicalReset;
    invalidPhysicalReset.currentNumericValid = false;
    require(!assess(invalidPhysicalReset).valid,
        "physical new baseline bypassed invalid current numerical endpoint");
    auto badReset = resetIdentity; badReset.gammaBeforePublicationReset = 0;
    require(!assess(badReset).valid, "no-op gamma was labelled an applied reset");
    badReset = resetIdentity; badReset.gammaBeforePublicationReset = std::numeric_limits<long long>::max();
    require(!assess(badReset).valid, "non-exact-double gamma identity was accepted");
    badReset = resetIdentity; badReset.hasPrevious = false;
    require(!assess(badReset).valid, "initial product fabricated an old gamma reference");
    for (double badValue : {nan, inf, -inf, 0.0, -0.1})
    {
        badReset = resetIdentity; badReset.publicationWavelengthMetres = badValue;
        require(!assess(badReset).valid, "invalid rebase wavelength was accepted");
    }
    for (double badValue : {nan, inf, -inf})
    {
        badReset = resetIdentity; badReset.appliedPhaseAdjustmentMetres = badValue;
        require(!assess(badReset).valid, "nonfinite applied reset was accepted");
    }
    std::cout << "PRODUCT_NUMERICAL_CONTINUITY_PASS checks=" << checks
        << " boundary120=assessed gap150=unassessed_new_segment"
        << " initial_reset=0 recovery_resets=1 measured_failure=suspended"
        << " gamma_reset_comparison=verified_native_map_only explicit_new_segment=1"
        << " same_segment_no_compensation=1 real_anomaly_still_rejected=1"
        << " ratio_unchanged=1 nis_unchanged=1 no_solver=1\n";
}
