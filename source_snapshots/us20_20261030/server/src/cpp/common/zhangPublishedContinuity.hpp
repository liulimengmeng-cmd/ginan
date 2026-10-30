#pragma once

#include "common/zhangProductDatumTransition.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

/** Version 1 is an endpoint-bound deterministic chart map, NOT a difference
 * between product values at adjacent epochs. The 49/52 legacy columns remain
 * unchanged. Unknown transport is written as nan, never as a fabricated zero.
 * Current Qp is in the emitted chart; this contract supplies no temporal Qp.
 */
struct ZhangPublishedContinuityContract
{
    int version = 1;
    double previousEpochGpstSeconds = std::numeric_limits<double>::quiet_NaN();
    int previousCounter = -1;
    int previousDatumVersion = -1;
    std::string previousPublicSegment = "UNAVAILABLE";
    bool floatTransportValid = false;
    double floatCorrectionShiftMetres = std::numeric_limits<double>::quiet_NaN();
    bool integerTransportValid = false;
    long long integerTransportShiftCycles = 0;
    std::string proofId = "UNAVAILABLE";
    std::string kind = "UNKNOWN";
    bool jointCovarianceComplete = false;
    std::string bundleId = "UNAVAILABLE";
};

inline constexpr const char* zhangPublishedContinuityCsvHeader =
    "continuity_contract_version,previous_epoch_gpst_seconds,"
    "previous_discontinuity_counter,previous_datum_version,"
    "previous_public_phase_segment_id,float_chart_transport_valid,"
    "float_correction_shift_m,integer_transport_valid,"
    "integer_transport_shift_cycles,continuity_proof_id,"
    "continuity_transport_kind,joint_covariance_complete,continuity_bundle_id";

inline bool zhangPublishedContinuityToken(const std::string& value)
{
    return !value.empty() && value != "UNAVAILABLE" && value != "UNRESOLVED"
        && value.find_first_of(",\r\n") == std::string::npos;
}

inline bool zhangPublishedContinuityValid(
    const ZhangPublishedContinuityContract& value, double currentEpoch,
    const std::string& currentPublicSegment)
{
    if ((value.version != 1 && value.version != 2) || !zhangPublishedContinuityToken(value.bundleId))
        return false;
    if(value.version==2)
    {
        ZhangProductDatumTransition t;
        if(!zhangReadDatumWireProof(value.proofId,t) || !zhangDatumEndpointValid(t.target) ||
            t.target.epoch!=currentEpoch || t.target.segment!=currentPublicSegment ||
            t.target.bundle!=value.bundleId)return false;
        if(value.kind!="UNKNOWN" && (!zhangDatumTransitionPureIntegerValid(t) ||
            t.source.epoch!=value.previousEpochGpstSeconds || t.source.counter!=value.previousCounter ||
            t.source.datum!=value.previousDatumVersion || t.source.segment!=value.previousPublicSegment ||
            t.correctionShiftMetres!=value.floatCorrectionShiftMetres ||
            (value.integerTransportValid && t.correctionShiftCycles!=value.integerTransportShiftCycles)))return false;
    }
    if (value.kind == "UNKNOWN")
        return !value.floatTransportValid && !value.integerTransportValid
            && std::isnan(value.floatCorrectionShiftMetres);
    if (value.kind != "IDENTITY" && value.kind != "CERTIFIED_AFFINE")
        return false;
    if (!value.floatTransportValid || !value.jointCovarianceComplete
        || !std::isfinite(value.floatCorrectionShiftMetres)
        || !std::isfinite(value.previousEpochGpstSeconds)
        || !std::isfinite(currentEpoch)
        || value.previousEpochGpstSeconds >= currentEpoch
        || value.previousCounter < 0 || value.previousDatumVersion < 0
        || !zhangPublishedContinuityToken(value.previousPublicSegment)
        || !zhangPublishedContinuityToken(currentPublicSegment)
        || !zhangPublishedContinuityToken(value.proofId))
        return false;
    if (value.kind == "IDENTITY")
        return value.previousPublicSegment == currentPublicSegment
            && value.floatCorrectionShiftMetres == 0
            && (!value.integerTransportValid || value.integerTransportShiftCycles == 0);
    return true;
}

// All booleans are independently established producer evidence. This factory
// never searches for an integer and never changes a product mean/covariance.
inline bool zhangAuthorizePublishedAffine(
    ZhangPublishedContinuityContract& out,double currentEpoch,
    const std::string& currentPublicSegment,const ZhangProductDatumTransition& transition,
    bool fullJointMapProven,bool currentIntegerValueProven,bool actualWriterProjectionBound,
    bool sameTimeNumericSignProven,bool integerAuthorizationProven)
{
    if(!fullJointMapProven || !currentIntegerValueProven || !actualWriterProjectionBound ||
        !sameTimeNumericSignProven || !out.jointCovarianceComplete ||
        !zhangDatumTransitionPureIntegerValid(transition))return false;
    auto candidate=out;
    candidate.version=2;
    candidate.kind="CERTIFIED_AFFINE";
    candidate.floatTransportValid=true;
    candidate.floatCorrectionShiftMetres=transition.correctionShiftMetres;
    candidate.integerTransportValid=integerAuthorizationProven;
    candidate.integerTransportShiftCycles=transition.correctionShiftCycles;
    candidate.proofId=zhangDatumWireProof(transition);
    if(!zhangPublishedContinuityValid(candidate,currentEpoch,currentPublicSegment))return false;
    out=std::move(candidate);return true;
}

/** Only call after an actual writer map, full joint moment map, current
 * integer authorization and final Qp coverage have independently passed.
 * The server has already absorbed its native affine/gamma offset into the
 * persistent public chart. Thus the PUBLIC-to-PUBLIC map is identity; gamma
 * is deliberately not exported as an additional client correction.
 */
inline bool zhangAuthorizePublishedIdentity(
    ZhangPublishedContinuityContract& out, double currentEpoch,
    const std::string& currentPublicSegment, bool samePublicChartProven,
    bool fullMomentMapProven, bool integerAuthorizationProven,
    const std::string& proofId)
{
    if (!samePublicChartProven || !fullMomentMapProven
        || !out.jointCovarianceComplete || !zhangPublishedContinuityToken(proofId))
        return false;
    auto candidate = out;
    candidate.floatTransportValid = true;
    candidate.floatCorrectionShiftMetres = 0;
    candidate.integerTransportValid = integerAuthorizationProven;
    candidate.integerTransportShiftCycles = 0;
    candidate.proofId = proofId;
    candidate.kind = "IDENTITY";
    if (!zhangPublishedContinuityValid(candidate, currentEpoch, currentPublicSegment))
        return false;
    out = std::move(candidate);
    return true;
}

inline void zhangWritePublishedContinuityCsv(
    std::ostream& out, const ZhangPublishedContinuityContract& value)
{
    out << std::setprecision(17) << value.version << ','
        << value.previousEpochGpstSeconds << ',' << value.previousCounter << ','
        << value.previousDatumVersion << ',' << value.previousPublicSegment << ','
        << value.floatTransportValid << ',';
    if (value.floatTransportValid) out << value.floatCorrectionShiftMetres;
    else out << "nan";
    out << ',' << value.integerTransportValid << ',';
    if (value.integerTransportValid) out << value.integerTransportShiftCycles;
    else out << "nan";
    out << ',' << value.proofId << ',' << value.kind << ','
        << value.jointCovarianceComplete << ',' << value.bundleId;
}

/** Publishing a PPP-only row does not permit filtering its Qp columns away.
 * Missing clock OR phase coverage is fail-closed for both PPP and AR.
 */
inline bool zhangPublishedProductHasNoiseParameters(
    bool clockPresent, bool phasePresent)
{
    return clockPresent && phasePresent;
}

inline bool zhangPublishedProductRequiresJointCovariance(
    const std::string& solution, bool pppUsable)
{
    // AR authorization is deliberately absent from this API. A PPP-only
    // row needs its noise parameters and cross covariance just as an AR row.
    return solution == "PRODUCT_FIXED" && pppUsable;
}

