#include "common/zhangProductNumericalContinuity.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

struct RealCase {
    const char* name;
    long long oldGamma;
    double currentQ, previousQ, currentSigma, previousSigma;
    double lambda, adjustment, expectedDelta, commonMode, limit, expectedResidual;
};
#include "four_star_cases.inc"

static int assertions = 0;
static void check(bool valid, const std::string& message) {
    ++assertions;
    if (!valid) throw std::runtime_error(message);
}
static bool close(double actual, double expected) {
    return std::isfinite(actual) && std::abs(actual - expected) < 1e-8;
}
static ZhangProductNumericalReferenceInput input(const RealCase& r) {
    ZhangProductNumericalReferenceInput in;
    in.currentNumericValid = true;
    in.hasPrevious = true;
    in.samePublishedIdentity = true;
    in.nominalRoute = ZhangProductPublicationRoute::NEW_SEGMENT;
    in.currentCorrection = r.currentQ;
    in.previousCorrection = r.previousQ;
    in.currentSigma = r.currentSigma;
    in.previousSigma = r.previousSigma;
    in.gapSeconds = 30;
    in.maximumGapSeconds = 120;
    in.publicationGammaResetApplied = true;
    in.nativeFullMomentMapProven = true;
    in.gammaBeforePublicationReset = r.oldGamma;
    in.gammaAfterPublicationReset = 0;
    in.publicationWavelengthMetres = r.lambda;
    // This input records the ACTUAL double operation performed by the writer,
    // not a separately rounded exact-decimal product. The oracle delta below
    // is independently reconstructed from frozen clock/phase operands.
    in.appliedPhaseAdjustmentMetres = r.lambda * (0.0 - static_cast<double>(r.oldGamma));
    return in;
}
static ZhangProductPublicationEvidence evidence() {
    ZhangProductPublicationEvidence e;
    e.nativeProductValid = e.currentSpatialCertificateValid = true;
    e.priorPublishedIdentity = e.samePhysicalSegment = true;
    e.sameProductFunctional = e.exactSameTimeMapProven = true;
    e.conditionalProductMapProven = true;
    e.requireFullJointHandoff = true;
    // The actual cases have a complete FLOAT map but no current integer
    // authorization for the OLD segment. These must remain false.
    e.fullJointHandoffProven = e.currentEpochBridgeProven = false;
    e.numericalStepValid = true;
    return e;
}
static bool assessed(const ZhangProductNumericalReferenceInput& in, const RealCase& r) {
    return zhangAssessProductNumericalContinuity(zhangProductNumericalReference(in),
                                               r.commonMode, r.limit).valid;
}

int main(int argc, char** argv) {
    try {
        for (const auto& r : realCases) {
            auto in = input(r);
            const auto original = in;
            auto old = in;
            old.publicationGammaResetApplied = false;
            check(!assessed(old, r), std::string(r.name) + " old raw comparison no longer reproduces failure");
            const auto ref = zhangProductNumericalReference(in);
            check(ref.hasDelta && close(ref.delta, r.expectedDelta), std::string(r.name) + " independent real CSV delta mismatch");
            const auto assessment = zhangAssessProductNumericalContinuity(ref, r.commonMode, r.limit);
            check(assessment.evaluated && assessment.valid && close(assessment.residualStep, r.expectedResidual),
                  std::string(r.name) + " independent real CSV residual mismatch");
            const auto decision = zhangProductPublicationWithNumericalAssessment(evidence(), assessment);
            check(decision.route == ZhangProductPublicationRoute::NEW_SEGMENT && decision.requireUserReset,
                  std::string(r.name) + " FLOAT map incorrectly authorized old integer continuity");
            check(in.currentCorrection == original.currentCorrection && in.gammaBeforePublicationReset == original.gammaBeforePublicationReset,
                  "Numerical assessment mutated publication input");

            auto bad = in; bad.nativeFullMomentMapProven = false;
            check(!assessed(bad, r), "Missing native full-moment map admitted");
            bad = in; bad.nominalRoute = ZhangProductPublicationRoute::SAME_SEGMENT;
            check(!assessed(bad, r), "Gamma reset authorized SAME_SEGMENT");
            bad = in; bad.gammaAfterPublicationReset = 1;
            check(!assessed(bad, r), "Nonzero new gamma admitted");
            bad = in; bad.appliedPhaseAdjustmentMetres += 0.25;
            check(!assessed(bad, r), "Arbitrary adjustment admitted");
            bad = in; bad.currentCorrection += 2 * r.limit + std::abs(r.expectedDelta - r.commonMode) + 1;
            check(!assessed(bad, r), "Actual additional native jump hidden by gamma cancellation");
            bad = in; bad.gapSeconds = 150;
            check(!assessed(bad, r), "Expired numerical reference admitted by gamma reset");
            bad = in; bad.currentSigma = std::numeric_limits<double>::quiet_NaN();
            check(!assessed(bad, r), "Invalid sigma admitted");
            bad = in; bad.appliedPhaseAdjustmentMetres = std::numeric_limits<double>::infinity();
            check(!assessed(bad, r), "Nonfinite adjustment admitted");
            auto missingSpatial = evidence(); missingSpatial.currentSpatialCertificateValid = false;
            check(zhangProductPublicationWithNumericalAssessment(missingSpatial, assessment).route ==
                  ZhangProductPublicationRoute::SUSPENDED, "Numerical map bypassed current integer spatial certificate");
        }
        const auto newSegment = ZhangProductPublicationDecision{ZhangProductPublicationRoute::NEW_SEGMENT, true, "TEST"};
        const auto suspended = ZhangProductPublicationDecision{ZhangProductPublicationRoute::SUSPENDED, false, "TEST"};
        auto pair = zhangAggregateDualFrequencyPublication(newSegment, newSegment);
        check(pair.dualFrequencyArUsable && pair.resetBothUserAmbiguities, "New pair did not require both user resets");
        pair = zhangAggregateDualFrequencyPublication(newSegment, suspended);
        check(!pair.dualFrequencyArUsable, "Suspended partner admitted dual-frequency AR");
        pair = zhangAggregateDualFrequencyPublication(std::optional<ZhangProductPublicationDecision>(newSegment), std::nullopt);
        check(!pair.dualFrequencyArUsable, "Missing partner admitted dual-frequency AR");
        // An explicitly claimed arithmetic reset must match the actual known
        // operation, not merely a relative tolerance that grows to kilometres.
        // Construct the q endpoint independently so a forged shift could
        // otherwise hide in an apparently harmless numerical residual.
        auto huge = input(realCases.front());
        huge.gammaBeforePublicationReset = 9007199254740991LL;
        const double exactHugeAdjustment = -huge.publicationWavelengthMetres *
            static_cast<double>(huge.gammaBeforePublicationReset);
        huge.appliedPhaseAdjustmentMetres = exactHugeAdjustment + 16.0;
        huge.currentCorrection = huge.previousCorrection + realCases.front().commonMode -
            huge.appliedPhaseAdjustmentMetres;
        check(zhangProductNumericalReference(huge).kind == ZhangProductNumericalReferenceKind::INVALID,
              "Large-gamma relative tolerance admitted a forged 16m reset shift");
        if (argc == 2) {
            std::ofstream out(argv[1]);
            out << "{\"schema\":\"INDEPENDENT_REAL_CSV_NUMERICAL_CORE_BLACKBOX_V1\",\"passed\":true,\"real_signal_cases\":8,\"assertions\":"
                << assertions << ",\"pea_integration_replay\":false,\"atomic_runtime_rollback_tested\":false,\"full_q_audit_tested\":false}\n";
            if (!out) throw std::runtime_error("Cannot save independent result");
        }
        std::cout << "PASS independent real CSV helper blackbox: " << assertions << " assertions, 8 real signal cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
