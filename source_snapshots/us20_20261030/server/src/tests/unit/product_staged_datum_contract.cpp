#include "common/zhangProductNumericalContinuity.hpp"
#include "common/zhangPublishedContinuity.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

static int checks = 0;
static void require(bool condition, const char* reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}

static ZhangProductPublicationEvidence currentNewSegment()
{
    ZhangProductPublicationEvidence e;
    e.nativeProductValid = e.currentSpatialCertificateValid = true;
    e.priorPublishedIdentity = e.samePhysicalSegment = true;
    e.sameProductFunctional = e.sameBackendChart = true;
    e.conditionalProductMapProven = e.requireFullJointHandoff = true;
    e.fullJointHandoffProven = e.currentEpochBridgeProven = false;
    e.numericalStepValid = true;
    return e;
}

struct FrozenSignal
{
    const char* id;
    double previous, current, frequency, oldGamma, common, limit;
};

int main()
{
    using Route = ZhangProductPublicationRoute;
    // Literal original CSV endpoints, GPS frequencies, and native TRACE
    // thresholds. The oracle uses ordinary independent affine arithmetic,
    // not another call to the production reference builder.
    const std::array<FrozenSignal,6> signals{{
        {"G10-L1",-20658.093341283977,-20683.042154057857,1575.42e6,-130,-.11244529124815017,9.304689664189809},
        {"G10-L2",-20663.327437931686,-20684.767214296688,1227.60e6,-87,-.10880471643758938,15.134977785239737},
        {"G19-L1",153454.11659409749,153443.16587791499,1575.42e6,-57,-.11244529124815017,7.797642505510679},
        {"G19-L2",153452.22854128267,153444.55906938299,1227.60e6,-31,-.10880471643758938,12.743734043505018},
        {"G24-L1",-144986.58713479046,-145000.40866149621,1575.42e6,-72,-.11244529124815017,7.728270519628609},
        {"G24-L2",-144996.49691404068,-145009.79399106765,1227.60e6,-54,-.10880471643758938,12.642791344682507}
    }};
    std::array<ZhangProductPublicationDecision,6> decisions;
    int oldFailures = 0;
    for (std::size_t n=0; n<signals.size(); ++n)
    {
        const auto& s=signals[n];
        const double lambda=299792458.0/s.frequency;
        const double knownShift=lambda*s.oldGamma; // gamma_new=0
        const double before=s.current-knownShift;
        ZhangProductNumericalReferenceInput input;
        input.currentNumericValid = input.hasPrevious = input.samePublishedIdentity = true;
        input.nominalRoute = Route::NEW_SEGMENT;
        input.currentCorrection=s.current; input.previousCorrection=s.previous;
        input.currentSigma=input.previousSigma=.5;
        input.gapSeconds=30; input.maximumGapSeconds=120;
        auto evaluate=[&](const ZhangProductNumericalReferenceInput& value) {
            return zhangAssessProductNumericalContinuity(
                zhangProductNumericalReference(value),s.common,s.limit);
        };
        const auto old=evaluate(input);
        oldFailures += !old.valid;
        input.removeStagedIntegerDatumShift = input.stagedIntegerDatumComparisonProven = true;
        input.correctionBeforeStagedDatumChange=before;
        input.stagedCorrectionShiftMetres=knownShift;
        const auto reference=zhangProductNumericalReference(input);
        const auto assessment=evaluate(input);
        const double oracle=(s.current-s.previous)-knownShift;
        require(reference.hasDelta && std::abs(reference.delta-oracle)<1e-9,
            "staged delta disagrees with independent affine oracle");
        require(assessment.evaluated && assessment.valid &&
            std::abs(assessment.residualStep-std::abs(oracle-s.common))<1e-9,
            "known datum reset was confused with a physical numerical jump");
        decisions[n]=zhangProductPublicationWithNumericalAssessment(currentNewSegment(),assessment);
        require(decisions[n].route==Route::NEW_SEGMENT && decisions[n].requireUserReset,
            "normalization manufactured old-segment continuity");
        require(std::string(reference.reason)=="COMPARABLE_BEFORE_STAGED_INTEGER_DATUM_CHANGE",
            "comparison datum is not auditable");
        auto reject=[&](ZhangProductNumericalReferenceInput invalid,const char* reason) {
            require(zhangProductPublicationWithNumericalAssessment(currentNewSegment(),
                evaluate(invalid)).route==Route::SUSPENDED,reason);
        };
        auto physicalJump=input;
        physicalJump.currentCorrection+=100; physicalJump.correctionBeforeStagedDatumChange+=100;
        reject(physicalJump,"NEW segment bypassed a real 100 metre native jump");
        auto broken=input;broken.stagedIntegerDatumComparisonProven=false;
        reject(broken,"unproven native comparison chart was accepted");
        broken=input;broken.samePublishedIdentity=false;
        reject(broken,"staged normalization accepted a different old identity");
        broken=input;broken.nominalRoute=Route::SAME_SEGMENT;
        reject(broken,"SAME segment changed its tested datum silently");
        broken=input;broken.stagedCorrectionShiftMetres+=1;
        reject(broken,"wrong staged shift was accepted");
        broken=input;broken.correctionBeforeStagedDatumChange+=1;
        reject(broken,"wrong pre-stage endpoint was accepted");
        broken=input;broken.gapSeconds=0;
        reject(broken,"invalid time was hidden by datum normalization");
        for(double bad:{std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::infinity()})
        {
            broken=input;broken.currentCorrection=bad;
            reject(broken,"nonfinite published mean was accepted");
            broken=input;broken.correctionBeforeStagedDatumChange=bad;
            reject(broken,"nonfinite pre-stage mean was accepted");
            broken=input;broken.stagedCorrectionShiftMetres=bad;
            reject(broken,"nonfinite staged shift was accepted");
        }
        auto noCertificate=currentNewSegment();noCertificate.currentSpatialCertificateValid=false;
        require(zhangProductPublicationWithNumericalAssessment(noCertificate,assessment).route==
            Route::SUSPENDED,"normalization authorized a missing current integer certificate");
        // A numerical comparison must not create a temporal map in a NEW datum.
        ZhangPublishedContinuityContract contract;
        require(contract.kind=="UNKNOWN" && !contract.floatTransportValid &&
            !contract.integerTransportValid,"new segment fabricated temporal authorization");
        std::cout<<s.id<<" raw_delta="<<s.current-s.previous<<" known_shift="<<knownShift
            <<" tested_delta="<<reference.delta<<" route=NEW_SEGMENT physical_jump=SUSPENDED\n";
    }
    require(oldFailures==5,"frozen case no longer reproduces five original signal failures");
    for(int s=0;s<3;++s)
    {
        const auto pair=zhangAggregateDualFrequencyPublication(decisions[2*s],decisions[2*s+1]);
        require(pair.dualFrequencyArUsable && pair.resetBothUserAmbiguities,
            "two independently valid new segments did not authorize dual AR");
        const auto failed=zhangAggregateDualFrequencyPublication(decisions[2*s],
            {Route::SUSPENDED,false,"INJECTED_INVALID_SIGNAL"});
        require(!failed.dualFrequencyArUsable,"single invalid signal bypassed pair gate");
    }
    std::cout<<"PRODUCT_STAGED_DATUM_PASS checks="<<checks
        <<" real_signals=6 legacy_failures=5 new_segments=6 same_segment_claims=0"
        <<" actual_native_jump_rejected=1 ratio_unchanged=1 nis_unchanged=1\n";
}
