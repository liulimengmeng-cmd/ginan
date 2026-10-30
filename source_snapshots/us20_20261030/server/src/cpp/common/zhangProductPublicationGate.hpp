#pragma once

#include <cstdlib>
#include <optional>
#include <string_view>

inline bool zhangR51SingleBasisHandoffEnabled()
{
	const char* value = std::getenv("ZHANG_R51_SINGLE_BASIS_HANDOFF");
	return value != nullptr && std::string_view(value) == "1";
}

inline bool zhangR51ProductRepairEnabled()
{
    const char* value = std::getenv("ZHANG_R51_PRODUCT_REPAIR");
    const char* epoch = std::getenv("ZHANG_R51_EPOCH_ONLY_AR");
    const char* enabled = std::getenv("ZHANG_R51_ENABLE");
    return value && epoch && enabled && std::string_view(value) == "1" &&
        std::string_view(epoch) == "1" && std::string_view(enabled) == "1";
}

/** Publication evidence for one satellite/signal product at one epoch.
 *
 * These are independently established facts, not configuration-mode flags.
 * In particular, a current spatial certificate does not prove transport of a
 * previously published product through a backend coordinate change.
 */
struct ZhangProductPublicationEvidence
{
	bool nativeProductValid = false;
	bool currentSpatialCertificateValid = false;
	bool priorPublishedIdentity = false;
	bool samePhysicalSegment = false;
	bool sameProductFunctional = false;
	bool sameBackendChart = false;
	bool exactSameTimeMapProven = false;
	bool conditionalProductMapProven = false;
	bool currentEpochBridgeProven = false;
	// Under the repaired writer even an unchanged relative functional needs
	// complete clock/L1/L2 chart evidence, and a nonzero segment difference
	// needs a freshly authorized TOTAL bridge, not a previous epoch's gamma.
	bool requireFullJointHandoff = false;
	bool fullJointHandoffProven = false;
	bool structuralProductMapProven = false;
	bool numericalStepValid = false;
    // A valid native nonzero bridge changes the PUBLIC endpoint. Certification
    // and reset authority are finalized only after complete Qp/dual publication.
    bool nativePublicChartChanges = false;
};

/** Current-epoch facts, frozen BEFORE any temporal/publication flags change.
 * A failed old-segment comparison cannot invalidate these facts, and neither
 * opening a segment nor clearing a mutable row flag can manufacture them. */
struct ZhangCurrentProductSpatialEvidence
{
    bool nativeProjectionValid = false;
    bool currentIntegerCertificateValid = false;
    bool currentJointMomentsValid = false;

    constexpr bool nativeValid() const
    { return nativeProjectionValid && currentJointMomentsValid; }
};

enum class ZhangProductPublicationRoute
{
	SAME_SEGMENT,
	NEW_SEGMENT,
	SUSPENDED
};

struct ZhangProductPublicationDecision
{
	ZhangProductPublicationRoute route = ZhangProductPublicationRoute::SUSPENDED;
	bool requireUserReset = false;
	const char* reason = "NOT_EVALUATED";
};

struct ZhangDualFrequencyPublicationDecision
{
	bool dualFrequencyArUsable = false;
	bool resetBothUserAmbiguities = false;
	const char* reason = "NOT_EVALUATED";
};

/** Aggregate two required phase signals for one satellite.  The per-signal
 * decisions are retained for diagnostics: a continuous L1 row does not lose
 * that fact merely because L2 opens a new segment.  User WL/L1 history,
 * however, spans both signals and must be reset together on either new
 * segment. A suspended signal cannot authorize a dual-frequency AR product. */
inline constexpr ZhangDualFrequencyPublicationDecision
zhangAggregateDualFrequencyPublication(
	const ZhangProductPublicationDecision& first,
	const ZhangProductPublicationDecision& second)
{
	const bool newSegment =
		first.route == ZhangProductPublicationRoute::NEW_SEGMENT ||
		second.route == ZhangProductPublicationRoute::NEW_SEGMENT;
	const bool suspended =
		first.route == ZhangProductPublicationRoute::SUSPENDED ||
		second.route == ZhangProductPublicationRoute::SUSPENDED;
	if (suspended)
		return {false,newSegment,
			newSegment ? "DUAL_FREQUENCY_SUSPENDED_WITH_NEW_SEGMENT"
				: "DUAL_FREQUENCY_SIGNAL_SUSPENDED"};
	if (newSegment)
		return {true,true,"DUAL_FREQUENCY_NEW_SEGMENT"};
	return {true,false,"DUAL_FREQUENCY_SAME_SEGMENT"};
}

/** A missing required phase row is suspended, even if the present row alone
 * carries a spatial certificate. No single-signal row can commit dual AR. */
inline constexpr ZhangDualFrequencyPublicationDecision
zhangAggregateDualFrequencyPublication(
	const std::optional<ZhangProductPublicationDecision>& first,
	const std::optional<ZhangProductPublicationDecision>& second)
{
	if (!first || !second)
	{
		const bool newSegment =
			(first && first->route == ZhangProductPublicationRoute::NEW_SEGMENT) ||
			(second && second->route == ZhangProductPublicationRoute::NEW_SEGMENT);
		return {false,newSegment,
			newSegment ? "DUAL_FREQUENCY_SIGNAL_MISSING_WITH_NEW_SEGMENT"
				: "DUAL_FREQUENCY_SIGNAL_MISSING"};
	}
	return zhangAggregateDualFrequencyPublication(*first,*second);
}

/** Fail closed on invalid current evidence. A retired physical arc can remain
 * continuous only when a conditional product map and a current-epoch integer
 * bridge prove the handoff; an unproved physical break opens a new segment.
 * A changed chart with the identical versioned physical product functional
 * needs an exact same-time map. A changed functional needs its conditional
 * product map and an integer bridge certified from the current epoch; past
 * integers alone are not a substitute. This gate establishes neither proof.
 */
inline constexpr ZhangProductPublicationDecision zhangProductPublicationGate(
	const ZhangProductPublicationEvidence& evidence)
{
	if (!evidence.nativeProductValid)
		return {ZhangProductPublicationRoute::SUSPENDED, false,
			"NATIVE_PRODUCT_INVALID"};
	if (!evidence.currentSpatialCertificateValid)
		return {ZhangProductPublicationRoute::SUSPENDED, false,
			"CURRENT_SPATIAL_CERTIFICATE_MISSING"};
	if (!evidence.numericalStepValid)
		return {ZhangProductPublicationRoute::SUSPENDED, false,
			"NUMERICAL_STEP_INVALID"};
	if (!evidence.priorPublishedIdentity)
		return {ZhangProductPublicationRoute::NEW_SEGMENT, false,
			"INITIAL_PRODUCT_SEGMENT"};
    if(evidence.nativePublicChartChanges)
        return {ZhangProductPublicationRoute::NEW_SEGMENT,true,
            "NATIVE_NEW_SEGMENT_PENDING_ENDPOINT_CERTIFICATE"};
	if (evidence.requireFullJointHandoff && !evidence.fullJointHandoffProven)
		return {ZhangProductPublicationRoute::NEW_SEGMENT, true,
			"FULL_JOINT_CURRENT_SEGMENT_HANDOFF_UNPROVEN"};
	if (evidence.conditionalProductMapProven &&
		evidence.currentEpochBridgeProven)
		return {ZhangProductPublicationRoute::SAME_SEGMENT, false,
			"CURRENT_EPOCH_CONDITIONAL_HANDOFF"};
	if (evidence.fullJointHandoffProven && evidence.structuralProductMapProven)
		return {ZhangProductPublicationRoute::SAME_SEGMENT, false,
			"FULL_JOINT_STRUCTURAL_IDENTITY_NO_AR_DECISION"};
	if (!evidence.samePhysicalSegment)
		return {ZhangProductPublicationRoute::NEW_SEGMENT, true,
			"PHYSICAL_SEGMENT_CHANGED"};
	if (evidence.sameProductFunctional && evidence.sameBackendChart)
		return {ZhangProductPublicationRoute::SAME_SEGMENT, false,
			"UNCHANGED_PRODUCT_IDENTITY"};
	if (evidence.sameProductFunctional && evidence.exactSameTimeMapProven)
		return {ZhangProductPublicationRoute::SAME_SEGMENT, false,
			"EXACT_SAME_FUNCTIONAL_TRANSPORT"};
	return {ZhangProductPublicationRoute::NEW_SEGMENT, true,
		"UNPROVEN_PRODUCT_HANDOFF"};
}
