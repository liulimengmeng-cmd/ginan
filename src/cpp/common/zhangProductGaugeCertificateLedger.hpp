#pragma once
#include "common/zhangRatioOnly.hpp"

#include "common/zhangIntegerDecisionProof.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "common/enums.h"
#include "common/satSys.hpp"
#include "common/zhangIntegerCandidateNis.hpp"
#include "common/zhangIntegerAudit.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"

enum class ProductGaugeCertificateState
{
	PROVISIONAL,
	CONFIRMING,
	ACTIVE,
	SUSPENDED,
	REVOKED
};

enum class ProductGaugeCertificateInactiveReason
{
	NONE,
	SUPERSEDED,
	CONFLICT,
	PHASE_SEGMENT_RETIRED,
	ALIGNMENT_SUSPENDED,
	PRODUCT_COMPONENT_REVOKED
};

/** Immutable statistical lineage for one admitted integer family.
 *
 * currentCertificateRisk is the only risk used to authorise a certificate in
 * the current epoch.  lifetimeBroadcastRisk is a monotone audit quantity and
 * must never be fed back into a current writer gate. */
struct ZhangProductGaugeRiskReceipt
{
	std::string certificateFamilyId;
	std::string exactLatticeHash;
	std::string affineRhsHash;
	std::string factorProvenanceHash;
	std::vector<std::string> parentCertificateIds;
	bool isNewIntegerDecision = false;
	double conditionalFailureProbability = 1;
	double currentCertificateRisk = 1;
	double lifetimeBroadcastRisk = 0;
};

inline const char* zhangProductGaugeCertificateInactiveReasonName(
	ProductGaugeCertificateInactiveReason reason)
{
	switch (reason)
	{
		case ProductGaugeCertificateInactiveReason::NONE:
			return "NONE";
		case ProductGaugeCertificateInactiveReason::SUPERSEDED:
			return "SUPERSEDED";
		case ProductGaugeCertificateInactiveReason::CONFLICT:
			return "CONFLICT";
		case ProductGaugeCertificateInactiveReason::PHASE_SEGMENT_RETIRED:
			return "PHASE_SEGMENT_RETIRED";
		case ProductGaugeCertificateInactiveReason::ALIGNMENT_SUSPENDED:
			return "ALIGNMENT_SUSPENDED";
		case ProductGaugeCertificateInactiveReason::PRODUCT_COMPONENT_REVOKED:
			return "PRODUCT_COMPONENT_REVOKED";
	}
	return "UNKNOWN";
}

inline const char* zhangProductGaugeCertificateStateName(
	ProductGaugeCertificateState state)
{
	switch (state)
	{
		case ProductGaugeCertificateState::PROVISIONAL: return "PROVISIONAL";
		case ProductGaugeCertificateState::CONFIRMING: return "CONFIRMING";
		case ProductGaugeCertificateState::ACTIVE:      return "ACTIVE";
		case ProductGaugeCertificateState::SUSPENDED:   return "SUSPENDED";
		case ProductGaugeCertificateState::REVOKED:     return "REVOKED";
	}
	return "UNKNOWN";
}

/** The product-integer ledger and the product-gauge certificate ledger share
 * one stable service runtime namespace.  Either feature therefore requires
 * the caller to propagate the runtime ID from producer to scheduler/solver. */
inline bool zhangPersistentProductEvidenceNeedsRuntime(
	bool productIntegerLedgerEnabled,
	bool productGaugeCertificateLedgerEnabled)
{
	return productIntegerLedgerEnabled || productGaugeCertificateLedgerEnabled;
}

/** Persistent frontend knowledge of one dual-frequency satellite product
 * relation.  This deliberately does not retain receiver ambiguity arcs or a
 * backend S-basis generation: those belong to ProductIntegerLedger.  The
 * identity is the satellite pair plus the four physical phase-product
 * segments.  A tree/basis change therefore cannot retire this certificate,
 * while a satellite phase-segment change does. */
struct ProductGaugeCertificate
{
	ZhangDecisionProofs decisionProofs;
	E_Sys system = E_Sys::NONE;
	E_ObsCode firstObservable = E_ObsCode::NONE;
	E_ObsCode secondObservable = E_ObsCode::NONE;
	SatSys satellite;
	SatSys reference;
	ZhangExactInteger wideLaneInteger = 0;
	ZhangExactInteger firstSignalInteger = 0;
	// Canonical satellite-space affine relation q.  This survives receiver
	// tree/chord reparameterisation because it contains no receiver arc.
	std::map<SatSys, ZhangExactInteger> productRelationQ;
	std::array<std::string, 2> satellitePhaseSegments;
	std::array<std::string, 2> referencePhaseSegments;
	// Complete physical support for an arbitrary canonical q row.  The legacy
	// endpoint arrays remain populated for pair certificates and compatibility.
	std::map<SatSys, std::array<std::string, 2>> phaseSegmentsBySatellite;
	std::string componentUuid;
	long validFrom = 0;
	long validTo = 0; // zero means open-ended
	std::vector<std::string> sourceCertificateIds;
	double failureProbabilityBudget = 1;
	// The first joint lattice admission is the statistical family that made
	// this certificate eligible for later confirmation.  Several pair rows
	// emitted by one square/component block share the same family and therefore
	// share one union-bound cost.  Keeping this immutable lineage prevents the
	// final Ledger union from either charging the same joint event once per row
	// or, more dangerously, treating historical Ledger evidence as Perr=0.
	std::string admissionFamilyId;
	double admissionFamilyFailureProbabilityBudget = 1;
	ZhangProductGaugeRiskReceipt riskReceipt;
	int componentId = -1;
	int componentVersion = 0;
	long firstCertified = 0;
	long lastConfirmed = 0;
	int confirmationEpochs = 0;
	// Admission evidence is stored with the certificate rather than inferred
	// from a current backend generation.  A gauge product therefore survives a
	// tree/S-basis reparameterisation, but cannot silently survive a failed
	// temporal alignment or a changed physical phase segment.
	bool wideLaneReliable = false;
	bool firstSignalReliable = false;
	bool exactProductLatticeMembership = false;
	bool jointNisPassed = false;
	bool cycleClosurePassed = false;
	bool temporalAlignmentCertified = false;
	bool noResidualDof = false;
	int independentSupportPaths = 0;
	// A fixed-lag candidate is born from a cumulative raw-factor window.  Two
	// overlapping window posteriors are not independent confirmations.  The
	// exact accepted-measurement prefix and the candidate-specific joint
	// [WL,L1] covariance are therefore retained until a strictly positive new
	// information increment has been demonstrated.
	bool requiresIndependentTemporalEvidence = false;
	std::string factorWindowIdentity;
	std::vector<std::size_t> factorSequences;
	std::size_t factorMeasurementRows = 0;
	Matrix2d temporalEvidenceCovariance = Matrix2d::Constant(
		std::numeric_limits<double>::quiet_NaN());
	ProductGaugeCertificateState state = ProductGaugeCertificateState::PROVISIONAL;
	ProductGaugeCertificateInactiveReason inactiveReason =
		ProductGaugeCertificateInactiveReason::NONE;
	bool currentAlignmentValid = false;
	bool active = false;
	std::string source = "PRODUCT_COMPONENT_GAUGE";
};

struct ProductGaugeCertificateLedgerUpdate
{
	bool valid = false;
	int inputCertificates = 0;
	int freshCertificates = 0;
	int confirmedCertificates = 0;
	int conflicts = 0;
	int retiredSegmentCertificates = 0;
	int suspendedCertificates = 0;
	int revokedCertificates = 0;
	int rejectedCandidates = 0;
	int duplicateEvidenceCandidates = 0;
	int nonNestedEvidenceCandidates = 0;
	int insufficientInformationGainCandidates = 0;
	int independentEvidenceConfirmations = 0;
	int maximumNewFactorCount = 0;
	int maximumConditionalInformationRank = 0;
	double maximumOverlapFraction = 0;
	double maximumConditionalInformationGain = 0;
	std::string failureReason = "NOT_EVALUATED";
};

struct ZhangProductGaugeTemporalEvidenceAudit
{
	bool valid = false;
	bool nestedFactorPrefix = false;
	bool independentInformation = false;
	int previousFactorCount = 0;
	int currentFactorCount = 0;
	int overlappingFactorCount = 0;
	int newFactorCount = 0;
	int previousInformationRank = 0;
	int currentInformationRank = 0;
	int conditionalInformationRank = 0;
	double overlapFraction = 0;
	double conditionalInformationGain = 0;
	std::string failureReason = "NOT_EVALUATED";
};

/** Audit whether a later cumulative fixed-lag window contains a genuinely new
 * information increment for the same dual-frequency product relation.
 *
 * The previous accepted-measurement IDs must be an exact subset of the later
 * prefix.  Precision is compared only in the two-dimensional [WL,L1] relation
 * space.  A new epoch or a repeated posterior is not evidence: at least one
 * new raw measurement factor and a positive-semidefinite, positive-rank
 * conditional precision increment are both mandatory. */
inline ZhangProductGaugeTemporalEvidenceAudit
zhangProductGaugeTemporalEvidenceIncrement(
	const ProductGaugeCertificate& previous,
	const ProductGaugeCertificate& current,
	double relativeTolerance = 1e-10)
{
	ZhangProductGaugeTemporalEvidenceAudit result;
	result.previousFactorCount = previous.factorSequences.size();
	result.currentFactorCount = current.factorSequences.size();
	if (!previous.requiresIndependentTemporalEvidence ||
		!current.requiresIndependentTemporalEvidence ||
		previous.factorWindowIdentity.empty() ||
		previous.factorWindowIdentity != current.factorWindowIdentity ||
		previous.factorSequences.empty() || current.factorSequences.empty() ||
		!std::is_sorted(previous.factorSequences.begin(),
			previous.factorSequences.end()) ||
		!std::is_sorted(current.factorSequences.begin(),
			current.factorSequences.end()) ||
		std::adjacent_find(previous.factorSequences.begin(),
			previous.factorSequences.end()) != previous.factorSequences.end() ||
		std::adjacent_find(current.factorSequences.begin(),
			current.factorSequences.end()) != current.factorSequences.end() ||
		!previous.temporalEvidenceCovariance.allFinite() ||
		!current.temporalEvidenceCovariance.allFinite())
	{
		result.failureReason = "TEMPORAL_EVIDENCE_METADATA_INVALID";
		return result;
	}
	std::vector<std::size_t> overlap;
	std::set_intersection(previous.factorSequences.begin(),
		previous.factorSequences.end(), current.factorSequences.begin(),
		current.factorSequences.end(), std::back_inserter(overlap));
	result.overlappingFactorCount = overlap.size();
	result.newFactorCount = std::max(0,
		result.currentFactorCount - result.overlappingFactorCount);
	result.overlapFraction = result.currentFactorCount > 0
		? static_cast<double>(result.overlappingFactorCount) /
			result.currentFactorCount : 0;
	result.nestedFactorPrefix = std::includes(
		current.factorSequences.begin(), current.factorSequences.end(),
		previous.factorSequences.begin(), previous.factorSequences.end());
	if (!result.nestedFactorPrefix)
	{
		result.valid = true;
		result.failureReason = "TEMPORAL_FACTOR_WINDOW_NOT_NESTED";
		return result;
	}
	if (result.newFactorCount <= 0)
	{
		result.valid = true;
		result.failureReason = "NO_NEW_RAW_MEASUREMENT_FACTOR";
		return result;
	}

	auto precision = [&](const Matrix2d& covariance, Matrix2d& information,
		int& rank)
	{
		const Matrix2d symmetric = 0.5 * (
			covariance + covariance.transpose());
		Eigen::SelfAdjointEigenSolver<Matrix2d> eigen(symmetric);
		if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
			return false;
		const double scale = std::max(1.0,
			eigen.eigenvalues().cwiseAbs().maxCoeff());
		const double tolerance = std::max(1e-14,
			relativeTolerance * scale);
		if (eigen.eigenvalues().minCoeff() < -tolerance) return false;
		Vector2d inverse = Vector2d::Zero();
		rank = 0;
		for (int index = 0; index < 2; index++)
		{
			if (eigen.eigenvalues()(index) <= tolerance) continue;
			inverse(index) = 1.0 / eigen.eigenvalues()(index);
			rank++;
		}
		information = eigen.eigenvectors() * inverse.asDiagonal() *
			eigen.eigenvectors().transpose();
		information = 0.5 * (information + information.transpose());
		return information.allFinite();
	};
	Matrix2d previousInformation;
	Matrix2d currentInformation;
	if (!precision(previous.temporalEvidenceCovariance, previousInformation,
			result.previousInformationRank) ||
		!precision(current.temporalEvidenceCovariance, currentInformation,
			result.currentInformationRank))
	{
		result.failureReason = "TEMPORAL_EVIDENCE_COVARIANCE_INVALID";
		return result;
	}
	const Matrix2d increment = 0.5 * (
		(currentInformation - previousInformation) +
		(currentInformation - previousInformation).transpose());
	Eigen::SelfAdjointEigenSolver<Matrix2d> incrementEigen(increment);
	if (incrementEigen.info() != Eigen::Success ||
		!incrementEigen.eigenvalues().allFinite())
	{
		result.failureReason = "TEMPORAL_INFORMATION_INCREMENT_INVALID";
		return result;
	}
	const double informationScale = std::max({1.0,
		previousInformation.cwiseAbs().maxCoeff(),
		currentInformation.cwiseAbs().maxCoeff()});
	const double informationTolerance = std::max(1e-12,
		relativeTolerance * informationScale);
	if (incrementEigen.eigenvalues().minCoeff() < -informationTolerance)
	{
		result.valid = true;
		result.failureReason = "NON_POSITIVE_CONDITIONAL_INFORMATION_INCREMENT";
		return result;
	}
	for (int index = 0; index < 2; index++)
	{
		if (incrementEigen.eigenvalues()(index) <= informationTolerance) continue;
		result.conditionalInformationRank++;
		result.conditionalInformationGain += incrementEigen.eigenvalues()(index);
	}
	result.valid = true;
	result.independentInformation = result.conditionalInformationRank > 0 &&
		result.conditionalInformationGain > informationTolerance;
	result.failureReason = result.independentInformation ? "NONE" :
		"NO_POSITIVE_CONDITIONAL_INFORMATION_GAIN";
	return result;
}

/** Preserve the original family-wise strength of a gauge certificate during
 * a later-posterior consistency check.  Pending certificates are correlated
 * historical hypotheses, not a new family whose budget should shrink with
 * the number of redundant rows retained by the ledger. */
inline double zhangProductGaugePosteriorFailureBudget(
	const ProductGaugeCertificate& certificate,
	double configuredFamilyFailureBudget)
{
	const double family = std::clamp(
		std::isfinite(configuredFamilyFailureBudget)
			? configuredFamilyFailureBudget : 1e-3,
		1e-12, 1.0);
	if (!std::isfinite(certificate.failureProbabilityBudget) ||
		certificate.failureProbabilityBudget < 0 ||
		zhangRatioStatisticalReject(certificate.failureProbabilityBudget > 1))
	{
		return family;
	}
	return std::min(family,
		std::max(1e-12, certificate.failureProbabilityBudget));
}

struct ZhangProductGaugePosteriorRecheck
{
	bool valid = false;
	bool sameInteger = false;
	bool deterministic = false;
	bool reliable = false;
	bool alternativeReliable = false;
	ZhangExactInteger nearestWideLaneInteger = 0;
	ZhangExactInteger nearestFirstSignalInteger = 0;
	double wideLaneFailureProbability = 1;
	double firstSignalFailureProbability = 1;
	double jointFailureProbability = 1;
	double nis = std::numeric_limits<double>::quiet_NaN();
	double nisThreshold = std::numeric_limits<double>::quiet_NaN();
	double alternativeNis = std::numeric_limits<double>::quiet_NaN();
	double alternativeNisThreshold = std::numeric_limits<double>::quiet_NaN();
	std::string failureReason = "NOT_EVALUATED";
};

struct ZhangProductGaugePosteriorProjection
{
	bool valid = false;
	ZhangExactVector productRow;
	Vector2d floatingValue = Vector2d::Zero();
	Matrix2d covariance = Matrix2d::Zero();
	std::string failureReason = "NOT_EVALUATED";
};

/** Map a reference-invariant canonical satellite relation q to the current
 * product chart and form the correlated [WL, L1] posterior.  The chart stores
 * z_s = phase_s - phase_reference, hence a valid q must sum to zero and the
 * reference coefficient is implicit. */
inline ZhangProductGaugePosteriorProjection
zhangProjectProductGaugeToPosterior(
	const ProductGaugeCertificate& certificate,
	const SatSys& referenceSatellite,
	const std::map<SatSys, int>& coordinate,
	int productDimension,
	const VectorXd& fullJointMean,
	const MatrixXd& fullJointCovariance)
{
	ZhangProductGaugePosteriorProjection result;
	if (productDimension <= 0 ||
		fullJointMean.size() != 2 * productDimension ||
		fullJointCovariance.rows() != 2 * productDimension ||
		fullJointCovariance.cols() != 2 * productDimension ||
		!fullJointMean.allFinite() || !fullJointCovariance.allFinite() ||
		certificate.productRelationQ.empty())
	{
		result.failureReason = "PRODUCT_GAUGE_POSTERIOR_MAPPING_INPUT_INVALID";
		return result;
	}
	result.productRow = ZhangExactVector(productDimension);
	ZhangExactInteger coefficientSum = 0;
	for (const auto& [satellite, coefficient] : certificate.productRelationQ)
	{
		if (coefficient == 0) continue;
		coefficientSum += coefficient;
		if (satellite == referenceSatellite) continue;
		auto found = coordinate.find(satellite);
		if (found == coordinate.end() || found->second < 0 ||
			found->second >= productDimension)
		{
			result.failureReason = "PRODUCT_GAUGE_SATELLITE_NOT_IN_CURRENT_BASIS";
			return result;
		}
		result.productRow[found->second] += coefficient;
	}
	if (coefficientSum != 0)
	{
		result.failureReason = "PRODUCT_GAUGE_RELATION_NOT_REFERENCE_INVARIANT";
		return result;
	}
	if (std::all_of(result.productRow.begin(), result.productRow.end(),
		[](const auto& value) { return value == 0; }))
	{
		result.failureReason = "PRODUCT_GAUGE_ZERO_PRODUCT_ROW";
		return result;
	}
	const VectorXd numeric = zhangExactRowToDouble(result.productRow);
	if (!numeric.allFinite())
	{
		result.failureReason = "PRODUCT_GAUGE_COEFFICIENT_OUT_OF_NUMERIC_RANGE";
		return result;
	}
	MatrixXd rows = MatrixXd::Zero(2, 2 * productDimension);
	rows.row(0).head(productDimension) = numeric.transpose();
	rows.row(0).tail(productDimension) = -numeric.transpose();
	rows.row(1).head(productDimension) = numeric.transpose();
	result.floatingValue = rows * fullJointMean;
	result.covariance = rows *
		(0.5 * (fullJointCovariance + fullJointCovariance.transpose())) *
		rows.transpose();
	if (!result.floatingValue.allFinite() || !result.covariance.allFinite())
	{
		result.failureReason = "PRODUCT_GAUGE_POSTERIOR_PROJECTION_NOT_FINITE";
		return result;
	}
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

/** Re-evaluate one pending dual-frequency satellite-product relation on a
 * later posterior without conditioning that posterior.  The two coordinates
 * are [WL, L1].  Marginal rounding errors are combined conservatively, while
 * the correlated joint innovation is tested with the full 2x2 covariance.
 * A repeated posterior is only a consistency confirmation; probabilities are
 * never multiplied across epochs. */
inline ZhangProductGaugePosteriorRecheck
zhangRecheckProductGaugeOnPosterior(
	const Vector2d& floatingValue,
	Matrix2d covariance,
	const ZhangExactInteger& expectedWideLaneInteger,
	const ZhangExactInteger& expectedFirstSignalInteger,
	double maximumFailureProbability,
	double nisAlpha)
{
	ZhangProductGaugePosteriorRecheck result;
	if (!floatingValue.allFinite() || !covariance.allFinite() ||
		maximumFailureProbability <= 0 || maximumFailureProbability >= 1 ||
		nisAlpha <= 0 || nisAlpha >= 1)
	{
		result.failureReason = "PRODUCT_GAUGE_RECHECK_INPUT_INVALID";
		return result;
	}
	covariance = 0.5 * (covariance + covariance.transpose());
	Eigen::SelfAdjointEigenSolver<Matrix2d> eigenSolver(covariance);
	if (eigenSolver.info() != Eigen::Success ||
		!eigenSolver.eigenvalues().allFinite())
	{
		result.failureReason = "PRODUCT_GAUGE_RECHECK_COVARIANCE_INVALID";
		return result;
	}
	const double largestEigenvalue = eigenSolver.eigenvalues().maxCoeff();
	const double rankTolerance = std::max(
		1e-14, 1e-12 * std::max(0.0, largestEigenvalue));
	if (eigenSolver.eigenvalues().minCoeff() < -rankTolerance ||
		eigenSolver.eigenvalues().minCoeff() < -rankTolerance)
	{
		result.failureReason = "PRODUCT_GAUGE_RECHECK_COVARIANCE_NOT_PSD";
		return result;
	}
	if (floatingValue(0) <
			static_cast<double>(std::numeric_limits<long long>::min()) ||
		floatingValue(0) >
			static_cast<double>(std::numeric_limits<long long>::max()) ||
		floatingValue(1) <
			static_cast<double>(std::numeric_limits<long long>::min()) ||
		floatingValue(1) >
			static_cast<double>(std::numeric_limits<long long>::max()))
	{
		result.failureReason = "PRODUCT_GAUGE_RECHECK_INTEGER_RANGE_INVALID";
		return result;
	}
	result.nearestWideLaneInteger = std::llround(floatingValue(0));
	result.nearestFirstSignalInteger = std::llround(floatingValue(1));
	result.sameInteger =
		result.nearestWideLaneInteger == expectedWideLaneInteger &&
		result.nearestFirstSignalInteger == expectedFirstSignalInteger;
	const Vector2d nearest(
		result.nearestWideLaneInteger.convert_to<double>(),
		result.nearestFirstSignalInteger.convert_to<double>());
	const Vector2d expected(
		expectedWideLaneInteger.convert_to<double>(),
		expectedFirstSignalInteger.convert_to<double>());
	result.wideLaneFailureProbability = zhangIntegerRoundFailureProbability(
		floatingValue(0) - nearest(0), covariance(0, 0));
	result.firstSignalFailureProbability = zhangIntegerRoundFailureProbability(
		floatingValue(1) - nearest(1), covariance(1, 1));
	result.jointFailureProbability = std::min(1.0,
		result.wideLaneFailureProbability +
		result.firstSignalFailureProbability);

	const Vector2d expectedInnovation = expected - floatingValue;
	const Vector2d alternativeInnovation = nearest - floatingValue;
	if (largestEigenvalue <= rankTolerance)
	{
		result.valid = true;
		result.deterministic = true;
		result.nis = 0;
		result.nisThreshold = std::numeric_limits<double>::infinity();
		result.alternativeNis = 0;
		result.alternativeNisThreshold =
			std::numeric_limits<double>::infinity();
		result.reliable = result.sameInteger &&
			expectedInnovation.lpNorm<Eigen::Infinity>() <= 1e-7;
		result.alternativeReliable = !result.sameInteger &&
			alternativeInnovation.lpNorm<Eigen::Infinity>() <= 1e-7 &&
			zhangRatioStatisticalAccept(result.jointFailureProbability <=
				maximumFailureProbability + 1e-12);
		result.failureReason = result.reliable ? "NONE" :
			(result.alternativeReliable ? "RELIABLE_INTEGER_CHANGED" :
			 "DETERMINISTIC_AFFINE_CONFLICT");
		return result;
	}

	const auto expectedNis = assessZhangIntegerCandidateNis(
		expectedInnovation, covariance, nisAlpha);
	const auto alternativeNis = assessZhangIntegerCandidateNis(
		alternativeInnovation, covariance, nisAlpha);
	result.valid = expectedNis.valid || alternativeNis.valid;
	result.nis = expectedNis.nis;
	result.nisThreshold = expectedNis.threshold;
	result.alternativeNis = alternativeNis.nis;
	result.alternativeNisThreshold = alternativeNis.threshold;
	result.reliable = result.sameInteger && expectedNis.valid &&
		zhangRatioStatisticalAccept(expectedNis.nis <= expectedNis.threshold) &&
		zhangRatioStatisticalAccept(result.jointFailureProbability <= maximumFailureProbability + 1e-12);
	result.alternativeReliable = !result.sameInteger && alternativeNis.valid &&
		zhangRatioStatisticalAccept(alternativeNis.nis <= alternativeNis.threshold) &&
		zhangRatioStatisticalAccept(result.jointFailureProbability <= maximumFailureProbability + 1e-12);
	result.failureReason = result.reliable ? "NONE" :
		(result.alternativeReliable ? "RELIABLE_INTEGER_CHANGED" :
		 (!result.sameInteger ? "INTEGER_CHANGED_UNRELIABLE" :
		  (!expectedNis.valid ? "NIS_INVALID" :
		   (zhangRatioStatisticalReject(expectedNis.nis > expectedNis.threshold) ? "NIS_REJECTED" :
		    "FAILURE_PROBABILITY_EXCEEDED"))));
	return result;
}

inline void zhangCanonicaliseProductGaugeCertificate(
	ProductGaugeCertificate& certificate)
{
	if (certificate.productRelationQ.empty() &&
		certificate.satellite != certificate.reference)
	{
		certificate.productRelationQ[certificate.satellite] = 1;
		certificate.productRelationQ[certificate.reference] = -1;
	}
	if (certificate.phaseSegmentsBySatellite.empty())
	{
		certificate.phaseSegmentsBySatellite[certificate.satellite] =
			certificate.satellitePhaseSegments;
		certificate.phaseSegmentsBySatellite[certificate.reference] =
			certificate.referencePhaseSegments;
	}
	if (certificate.reference < certificate.satellite)
	{
		std::swap(certificate.satellite, certificate.reference);
		std::swap(certificate.satellitePhaseSegments,
			certificate.referencePhaseSegments);
		certificate.wideLaneInteger = -certificate.wideLaneInteger;
		certificate.firstSignalInteger = -certificate.firstSignalInteger;
		for (auto& [_, coefficient] : certificate.productRelationQ)
			coefficient = -coefficient;
	}
}

/** Only a literal primitive e_s-e_r relation may carry writer edge metadata.
 * A dense gauge row can imply pair relations through exact lattice membership,
 * but its arbitrarily selected positive/negative representatives are not an
 * edge until that pair consequence is materialised explicitly. */
inline bool zhangProductGaugeCertificateIsPrimitivePair(
	const ProductGaugeCertificate& certificate)
{
	if (certificate.satellite == certificate.reference ||
		certificate.productRelationQ.size() != 2) return false;
	const auto satellite = certificate.productRelationQ.find(
		certificate.satellite);
	const auto reference = certificate.productRelationQ.find(
		certificate.reference);
	return satellite != certificate.productRelationQ.end() &&
		reference != certificate.productRelationQ.end() &&
		satellite->second == 1 && reference->second == -1;
}

inline bool zhangResolvedPhaseProductSegment(const std::string& id)
{
	return !id.empty() &&
		id != "UNRESOLVED" &&
		id != "NONE" &&
		id.find("-SEG") != std::string::npos;
}

inline bool zhangProductGaugeCertificateSegmentsResolved(
	const ProductGaugeCertificate& certificate)
{
	if (certificate.productRelationQ.empty() ||
		certificate.phaseSegmentsBySatellite.empty()) return false;
	for (const auto& [satellite, coefficient] : certificate.productRelationQ)
	{
		if (coefficient == 0) continue;
		auto segments = certificate.phaseSegmentsBySatellite.find(satellite);
		if (segments == certificate.phaseSegmentsBySatellite.end() ||
			!zhangResolvedPhaseProductSegment(segments->second[0]) ||
			!zhangResolvedPhaseProductSegment(segments->second[1])) return false;
	}
	return true;
}

inline std::string zhangProductGaugeCertificatePairKey(
	const ProductGaugeCertificate& certificate)
{
	std::ostringstream stream;
	stream << enum_to_string(certificate.system) << "|"
		<< static_cast<int>(certificate.firstObservable) << "|"
		<< static_cast<int>(certificate.secondObservable) << "|q=";
	for (const auto& [satellite, coefficient] : certificate.productRelationQ)
		if (coefficient != 0)
			stream << static_cast<int>(satellite) << ":" << coefficient << ";";
	return stream.str();
}

inline std::string zhangProductGaugeCertificateIdentity(
	const ProductGaugeCertificate& certificate)
{
	std::ostringstream stream;
	stream << zhangProductGaugeCertificatePairKey(certificate) << "|segments=";
	for (const auto& [satellite, segments] :
		certificate.phaseSegmentsBySatellite)
		stream << static_cast<int>(satellite) << ":"
			<< segments[0] << "," << segments[1] << ";";
	return stream.str();
}

struct ZhangProductGaugeAdmissionFamilyBudgetAudit
{
	bool valid = false;
	int certificates = 0;
	int families = 0;
	int legacyFallbackFamilies = 0;
	double failureProbability = 1;
	std::string failureReason = "NOT_EVALUATED";
};

/** Structural score used when several independently admitted gauge families
 * compete for one global error-probability budget.  The largest common integer
 * component is primary because disconnected pair islands do not share a usable
 * product datum.  Within an equal largest component, exact dual graph rank and
 * total certified coverage remain the secondary completion objectives. */
struct ZhangProductGaugeFamilySubsetScore
{
	bool valid = false;
	// PPP-AR users can resolve only within one common integer datum component.
	// Prefer consolidating a large component before accumulating rank in
	// disconnected pair islands.
	int largestCertifiedComponentSize = 0;
	int certifiedComponentCount = 0;
	int dualGraphRank = 0;
	int certifiedSatelliteCount = 0;
	int conditioningRank = 0;
	int wideLaneRank = 0;
	int firstSignalRank = 0;
	int continuityPriority = 0;
	double referenceInvariantProductGain = 0;
};

struct ZhangProductGaugeFamilySubsetSelection
{
	bool valid = false;
	bool exactEnumeration = false;
	bool fallbackAvailable = false;
	bool fallbackPreserved = false;
	int familyCount = 0;
	int individuallyBudgetEligibleFamilies = 0;
	int evaluatedSubsets = 0;
	int budgetFeasibleSubsets = 0;
	std::uint64_t searchableNonEmptySubsets = 0;
	std::uint64_t evaluatedUniqueSubsets = 0;
	std::uint64_t prunedNotTestedSubsets = 0;
	double failureProbability = 0;
	std::vector<int> selectedFamilyIndices;
	ZhangProductGaugeFamilySubsetScore score;
};

inline bool zhangProductGaugeFamilySubsetScoreBetter(
	const ZhangProductGaugeFamilySubsetScore& left,
	double leftFailureProbability,
	const std::vector<int>& leftFamilies,
	const ZhangProductGaugeFamilySubsetScore& right,
	double rightFailureProbability,
	const std::vector<int>& rightFamilies,
	double tolerance = 1e-15)
{
	if (left.valid != right.valid) return left.valid;
	if (!left.valid)
	{
		const int leftBalanced = std::min(
			left.wideLaneRank, left.firstSignalRank);
		const int rightBalanced = std::min(
			right.wideLaneRank, right.firstSignalRank);
		if (leftBalanced != rightBalanced)
			return leftBalanced > rightBalanced;
		const int leftPartial = left.wideLaneRank + left.firstSignalRank;
		const int rightPartial = right.wideLaneRank + right.firstSignalRank;
		if (leftPartial != rightPartial) return leftPartial > rightPartial;
	}
	if (left.largestCertifiedComponentSize !=
		right.largestCertifiedComponentSize)
		return left.largestCertifiedComponentSize >
			right.largestCertifiedComponentSize;
	if (left.certifiedSatelliteCount != right.certifiedSatelliteCount)
		return left.certifiedSatelliteCount > right.certifiedSatelliteCount;
	if (left.dualGraphRank != right.dualGraphRank)
		return left.dualGraphRank > right.dualGraphRank;
	if (left.continuityPriority != right.continuityPriority)
		return left.continuityPriority > right.continuityPriority;
	if (left.certifiedComponentCount != right.certifiedComponentCount)
		return left.certifiedComponentCount < right.certifiedComponentCount;
	if (left.conditioningRank != right.conditioningRank)
		return left.conditioningRank > right.conditioningRank;
	const double leftGain = std::isfinite(left.referenceInvariantProductGain)
		? left.referenceInvariantProductGain : -std::numeric_limits<double>::infinity();
	const double rightGain = std::isfinite(right.referenceInvariantProductGain)
		? right.referenceInvariantProductGain : -std::numeric_limits<double>::infinity();
	if (std::abs(leftGain - rightGain) > tolerance) return leftGain > rightGain;
	if (std::abs(leftFailureProbability - rightFailureProbability) > tolerance)
		return leftFailureProbability < rightFailureProbability;
	return leftFamilies < rightFamilies;
}

/** Select a subset of jointly admitted evidence families without exceeding a
 * single global Perr budget.
 *
 * The evaluator performs the exact affine-union, posterior-NIS and product
 * graph audit for one subset.  Up to exactFamilyLimit families all feasible
 * subsets are enumerated.  Larger problems use a deterministic bounded beam;
 * this keeps runtime bounded while preserving the same safety gates.  Invalid
 * subsets are retained in the beam after valid ones so a later complementary
 * family can still turn an individually non-broadcastable lattice into a valid
 * dual-frequency product lattice. */
template<typename Evaluator>
inline ZhangProductGaugeFamilySubsetSelection
zhangSelectProductGaugeAdmissionFamilies(
	const std::vector<double>& familyFailureProbabilities,
	double maximumFamilyFailureProbability,
	Evaluator&& evaluate,
	int exactFamilyLimit = 8,
	int maximumBeamStates = 256,
	const std::vector<int>& fallbackFamilyIndices = {},
	const std::vector<std::string>& parentFamilyIds = {},
	const std::vector<ZhangDecisionProofs>& decisionProofsByFamily = {})
{
	ZhangProductGaugeFamilySubsetSelection result;
	result.familyCount = static_cast<int>(familyFailureProbabilities.size());
	if (!std::isfinite(maximumFamilyFailureProbability) ||
		maximumFamilyFailureProbability <= 0 || exactFamilyLimit < 0 ||
		maximumBeamStates <= 0) return result;
	const double budget = std::min(1.0, maximumFamilyFailureProbability);
	const double tolerance = std::max(1e-15, 1e-12 * budget);
	if (!parentFamilyIds.empty() && parentFamilyIds.size() != familyFailureProbabilities.size())
		return result;
	if (!decisionProofsByFamily.empty() && decisionProofsByFamily.size() != familyFailureProbabilities.size())
		return result;
	// Subblocks of one joint decision retain its entire original upper bound.
	// Selecting two such blocks charges the parent once, never zero or twice.
	auto subsetRisk = [&](const std::vector<int>& selected)
	{
		if (!decisionProofsByFamily.empty())
		{
			ZhangDecisionProofs roots;
			for (const int index : selected)
			{
				if (index < 0 || index >= result.familyCount || decisionProofsByFamily[index].empty())
					return std::numeric_limits<double>::infinity();
				roots = zhangMergeDecisionProofs(roots, decisionProofsByFamily[index]);
			}
			const auto closure = zhangDecisionRiskClosure(roots);
			return closure.valid ? closure.bound : std::numeric_limits<double>::infinity();
		}
		std::map<std::string, double> risks;
		for (int index : selected)
		{
			if (index < 0 || index >= result.familyCount ||
				!std::isfinite(familyFailureProbabilities[index]) || familyFailureProbabilities[index] < 0)
				return std::numeric_limits<double>::infinity();
			const auto id = parentFamilyIds.empty() ? std::to_string(index) : parentFamilyIds[index];
			if (id.empty()) return std::numeric_limits<double>::infinity();
			risks[id] = std::max(risks[id], familyFailureProbabilities[index]);
		}
		double total = 0;
		for (const auto& [id, risk] : risks) total += risk;
		return total;
	};
	for (const double probability : familyFailureProbabilities)
	{
		if (std::isfinite(probability) && probability >= 0 &&
			zhangRatioStatisticalAccept(probability <= budget + tolerance))
			result.individuallyBudgetEligibleFamilies++;
	}

	struct State
	{
		std::vector<int> selected;
		double failureProbability = 0;
		ZhangProductGaugeFamilySubsetScore score;
	};
	std::map<std::vector<int>, ZhangProductGaugeFamilySubsetScore> scoreCache;
	auto evaluateState = [&](State& state)
	{
		auto selected = state.selected;
		std::sort(selected.begin(), selected.end());
		selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
		state.selected = std::move(selected);
		if (const auto cached = scoreCache.find(state.selected);
			cached != scoreCache.end())
		{
			state.score = cached->second;
			return;
		}
		state.score = evaluate(state.selected);
		scoreCache.emplace(state.selected, state.score);
		result.evaluatedSubsets++;
		result.evaluatedUniqueSubsets++;
	};
	std::optional<ZhangProductGaugeFamilySubsetScore> fallbackFloor;
	auto consider = [&](const State& state)
	{
		if (!state.score.valid) return;
		if (fallbackFloor && (state.score.dualGraphRank < fallbackFloor->dualGraphRank ||
			state.score.certifiedSatelliteCount < fallbackFloor->certifiedSatelliteCount)) return;
		if (!result.valid || zhangProductGaugeFamilySubsetScoreBetter(
			state.score, state.failureProbability, state.selected,
			result.score, result.failureProbability,
			result.selectedFamilyIndices, tolerance))
		{
			result.valid = true;
			result.score = state.score;
			result.failureProbability = state.failureProbability;
			result.selectedFamilyIndices = state.selected;
		}
	};
	// A previously feasible baseline is an immutable fallback, not a mandatory
	// member of every trial.  Better safe historical subsets may replace it, but
	// a later local conflict or beam truncation can never erase it.
	if (!fallbackFamilyIndices.empty())
	{
		State fallback;
		fallback.selected = fallbackFamilyIndices;
		std::sort(fallback.selected.begin(), fallback.selected.end());
		fallback.selected.erase(std::unique(
			fallback.selected.begin(), fallback.selected.end()),
			fallback.selected.end());
		bool fallbackBudgetValid = true;
		for (const int family : fallback.selected)
		{
			if (family < 0 || family >= result.familyCount)
			{
				fallbackBudgetValid = false;
				break;
			}
			const double probability = familyFailureProbabilities[family];
			if (!std::isfinite(probability) || probability < 0)
			{
				fallbackBudgetValid = false;
				break;
			}
		}
		fallback.failureProbability = subsetRisk(fallback.selected);
		if (fallbackBudgetValid &&
			zhangRatioStatisticalAccept(fallback.failureProbability <= budget + tolerance))
		{
			evaluateState(fallback);
			result.fallbackAvailable = fallback.score.valid;
			if (fallback.score.valid) fallbackFloor = fallback.score;
			consider(fallback);
		}
	}

	const int familyCount = result.familyCount;
	result.searchableNonEmptySubsets = familyCount < 64
		? ((std::uint64_t{1} << familyCount) - 1)
		: std::numeric_limits<std::uint64_t>::max();
	result.exactEnumeration = familyCount <= exactFamilyLimit && familyCount < 63;
	if (result.exactEnumeration)
	{
		const std::uint64_t subsetCount = std::uint64_t{1} << familyCount;
		for (std::uint64_t mask = 1; mask < subsetCount; mask++)
		{
			State state;
			bool validBudget = true;
			for (int family = 0; family < familyCount; family++)
			{
				if (!(mask & (std::uint64_t{1} << family))) continue;
				const double probability = familyFailureProbabilities[family];
				if (!std::isfinite(probability) || probability < 0)
				{
					validBudget = false;
					break;
				}
				state.selected.push_back(family);
			}
			state.failureProbability = subsetRisk(state.selected);
			if (!validBudget || zhangRatioStatisticalReject(state.failureProbability > budget + tolerance))
				continue;
			result.budgetFeasibleSubsets++;
			evaluateState(state);
			consider(state);
		}
		result.fallbackPreserved = result.fallbackAvailable && result.valid;
		result.prunedNotTestedSubsets = 0;
		return result;
	}

	std::vector<State> frontier(1);
	for (int family = 0; family < familyCount; family++)
	{
		std::vector<State> next = frontier;
		const double probability = familyFailureProbabilities[family];
		if (std::isfinite(probability) && probability >= 0)
		{
			for (const auto& state : frontier)
			{
				State candidate = state;
				candidate.selected.push_back(family);
				candidate.failureProbability = subsetRisk(candidate.selected);
				if (zhangRatioStatisticalReject(candidate.failureProbability > budget + tolerance)) continue;
				evaluateState(candidate);
				result.budgetFeasibleSubsets++;
				next.push_back(std::move(candidate));
			}
		}
		std::vector<State> validStates;
		std::vector<State> invalidStates;
		std::optional<State> emptySeed;
		for (auto& state : next)
		{
			if (state.selected.empty())
			{
				emptySeed = std::move(state);
				continue;
			}
			(state.score.valid ? validStates : invalidStates)
				.push_back(std::move(state));
		}
		auto beamBetter = [&](const State& left, const State& right)
		{
			return zhangProductGaugeFamilySubsetScoreBetter(
				left.score, left.failureProbability, left.selected,
				right.score, right.failureProbability, right.selected, tolerance);
		};
		std::stable_sort(validStates.begin(), validStates.end(), beamBetter);
		std::stable_sort(invalidStates.begin(), invalidStates.end(), beamBetter);

		frontier.clear();
		if (emptySeed) frontier.push_back(std::move(*emptySeed));
		// The empty state is an expansion seed, not a scored beam candidate.
		// Do not let it consume the configured non-empty state capacity.
		const int remainingCapacity = maximumBeamStates;
		const int invalidReserve = invalidStates.empty() ? 0 :
			(validStates.empty() ? remainingCapacity :
			 (remainingCapacity >= 2 ? std::min(
				 remainingCapacity - 1,
				 std::max(1, remainingCapacity / 4)) : 0));
		const int validReserve = std::min(
			static_cast<int>(validStates.size()),
			remainingCapacity - invalidReserve);
		frontier.insert(frontier.end(),
			std::make_move_iterator(validStates.begin()),
			std::make_move_iterator(validStates.begin() + validReserve));
		const int retainedInvalid = std::min(
			static_cast<int>(invalidStates.size()), invalidReserve);
		frontier.insert(frontier.end(),
			std::make_move_iterator(invalidStates.begin()),
			std::make_move_iterator(invalidStates.begin() + retainedInvalid));

		const int retainedSeed = emptySeed ? 1 : 0;
		int freeCapacity = maximumBeamStates -
			(static_cast<int>(frontier.size()) - retainedSeed);
		for (int index = validReserve;
			index < static_cast<int>(validStates.size()) && freeCapacity > 0;
			index++, freeCapacity--)
			frontier.push_back(std::move(validStates[index]));
		for (int index = retainedInvalid;
			index < static_cast<int>(invalidStates.size()) && freeCapacity > 0;
			index++, freeCapacity--)
			frontier.push_back(std::move(invalidStates[index]));
	}
	for (const auto& state : frontier) consider(state);
	result.fallbackPreserved = result.fallbackAvailable && result.valid;
	result.prunedNotTestedSubsets =
		result.searchableNonEmptySubsets ==
			std::numeric_limits<std::uint64_t>::max()
		? std::numeric_limits<std::uint64_t>::max()
		: result.searchableNonEmptySubsets - std::min(
			result.searchableNonEmptySubsets, result.evaluatedUniqueSubsets);
	return result;
}

/** Conservatively account for a set of persistent gauge certificates.
 *
 * Rows born from one jointly admitted lattice share a family id and are
 * charged once at the maximum recorded bound for that family.  Older callers
 * without family metadata fail safely to one independent family per exact
 * certificate identity. */
inline ZhangProductGaugeAdmissionFamilyBudgetAudit
zhangProductGaugeAdmissionFamilyBudgetAudit(
	const std::vector<const ProductGaugeCertificate*>& certificates)
{
	ZhangProductGaugeAdmissionFamilyBudgetAudit result;
	result.certificates = certificates.size();
	std::map<std::string, double> familyBudgets;
	for (const auto* certificate : certificates)
	{
		if (!certificate)
		{
			result.failureReason = "PRODUCT_GAUGE_ADMISSION_FAMILY_NULL";
			return result;
		}
		const bool hasReceipt =
			!certificate->riskReceipt.certificateFamilyId.empty();
		const bool legacy = !hasReceipt && certificate->admissionFamilyId.empty();
		const std::string family = hasReceipt
			? certificate->riskReceipt.certificateFamilyId
			: (legacy
				? "LEGACY|" + zhangProductGaugeCertificateIdentity(*certificate)
				: certificate->admissionFamilyId);
		const double budget = hasReceipt
			? certificate->riskReceipt.currentCertificateRisk
			: (legacy
				? certificate->failureProbabilityBudget
				: certificate->admissionFamilyFailureProbabilityBudget);
		if (!std::isfinite(budget) || budget < 0 || budget > 1)
		{
			result.failureReason =
				"PRODUCT_GAUGE_ADMISSION_FAMILY_BUDGET_INVALID";
			return result;
		}
		if (legacy) result.legacyFallbackFamilies++;
		auto [existing, inserted] = familyBudgets.try_emplace(family, budget);
		if (!inserted) existing->second = std::max(existing->second, budget);
	}
	double total = 0;
	for (const auto& [_, budget] : familyBudgets)
		total = std::min(1.0, total + budget);
	result.valid = true;
	result.families = familyBudgets.size();
	result.failureProbability = total;
	result.failureReason = "NONE";
	return result;
}

inline std::string zhangProductGaugeCertificateFunctionalFingerprint(
	const ProductGaugeCertificate& certificate)
{
	std::ostringstream stream;
	stream << enum_to_string(certificate.system) << "|"
		<< static_cast<int>(certificate.firstObservable) << "|"
		<< static_cast<int>(certificate.secondObservable) << "|q=";
	for (const auto& [satellite, coefficient] : certificate.productRelationQ)
		if (coefficient != 0)
			stream << static_cast<int>(satellite) << ":" << coefficient << ";";
	stream << "|segments=";
	for (const auto& [satellite, segments] :
		certificate.phaseSegmentsBySatellite)
		stream << static_cast<int>(satellite) << ":"
			<< segments[0] << "," << segments[1] << ";";
	stream << "|affine=" << certificate.wideLaneInteger
		<< "," << certificate.firstSignalInteger;
	return stream.str();
}

inline std::string zhangProductGaugeRiskHash(const std::string& text)
{
	std::ostringstream stream;
	stream << std::hex << std::setw(16) << std::setfill('0')
		<< zhangAuditFnv1a(1469598103934665603ULL, text);
	return stream.str();
}

/** Populate the immutable receipt at first admission.  Reconfirmations keep
 * the original receipt: a different view or a later epoch is not a new
 * integer decision and must not be charged again. */
inline void zhangInitialiseProductGaugeRiskReceipt(
	ProductGaugeCertificate& certificate)
{
	auto& receipt = certificate.riskReceipt;
	if (!receipt.certificateFamilyId.empty()) return;
	const bool legacy = certificate.admissionFamilyId.empty();
	receipt.certificateFamilyId = legacy
		? "LEGACY|" + zhangProductGaugeCertificateIdentity(certificate)
		: certificate.admissionFamilyId;
	receipt.exactLatticeHash = zhangProductGaugeRiskHash(
		zhangProductGaugeCertificatePairKey(certificate));
	receipt.affineRhsHash = zhangProductGaugeRiskHash(
		certificate.wideLaneInteger.convert_to<std::string>() + "|" +
		certificate.firstSignalInteger.convert_to<std::string>());
	std::ostringstream provenance;
	provenance << certificate.factorWindowIdentity << "|rows="
		<< certificate.factorMeasurementRows << "|seq=";
	for (const auto sequence : certificate.factorSequences)
		provenance << sequence << ",";
	std::vector<std::string> parents = certificate.sourceCertificateIds;
	std::sort(parents.begin(), parents.end());
	parents.erase(std::unique(parents.begin(), parents.end()), parents.end());
	for (const auto& parent : parents) provenance << "|parent=" << parent;
	receipt.factorProvenanceHash = zhangProductGaugeRiskHash(provenance.str());
	receipt.parentCertificateIds = std::move(parents);
	receipt.isNewIntegerDecision = true;
	receipt.conditionalFailureProbability = legacy
		? certificate.failureProbabilityBudget
		: certificate.admissionFamilyFailureProbabilityBudget;
	receipt.currentCertificateRisk = receipt.conditionalFailureProbability;
}

inline bool zhangProductGaugeCertificateEvidenceComplete(
	const ProductGaugeCertificate& certificate,
	int requiredConfirmations)
{
	if (!certificate.wideLaneReliable || !certificate.firstSignalReliable ||
		!certificate.exactProductLatticeMembership || !certificate.jointNisPassed ||
		!certificate.cycleClosurePassed || !certificate.temporalAlignmentCertified)
		return false;
	if (certificate.requiresIndependentTemporalEvidence &&
		(certificate.factorWindowIdentity.empty() ||
		 certificate.factorSequences.empty() ||
		 certificate.factorMeasurementRows == 0 ||
		 !certificate.temporalEvidenceCovariance.allFinite()))
		return false;
	// A square bridge has no residual degrees of freedom.  It needs either a
	// separate support path or a multi-epoch confirmation; a one-epoch rounded
	// relation is never a strict dual-frequency certificate.
	if (certificate.noResidualDof &&
		certificate.independentSupportPaths <= 0 && requiredConfirmations < 2)
		return false;
	return true;
}

/** A frontend certificate is usable only for the exact pair of current
 * physical phase-product segments that created it.  This is intentionally a
 * segment test, not a backend-generation or product-tree test: changing an
 * S-basis leaves the physical relation intact, whereas replacing either
 * satellite phase segment invalidates it. */
inline bool zhangProductGaugeCertificateMatchesCurrentSegments(
	const ProductGaugeCertificate& certificate,
	E_ObsCode firstObservable,
	E_ObsCode secondObservable,
	const std::array<std::string, 2>& satelliteSegments,
	const std::array<std::string, 2>& referenceSegments)
{
	return certificate.state == ProductGaugeCertificateState::ACTIVE
		&& certificate.active
		&& certificate.currentAlignmentValid
		&& certificate.firstObservable == firstObservable
		&& certificate.secondObservable == secondObservable
		&& certificate.satellitePhaseSegments == satelliteSegments
		&& certificate.referencePhaseSegments == referenceSegments;
}

inline bool zhangProductGaugeCertificateMatchesCurrentSegments(
	const ProductGaugeCertificate& certificate,
	E_ObsCode firstObservable,
	E_ObsCode secondObservable,
	const std::map<SatSys, std::array<std::string, 2>>& currentSegments)
{
	if (certificate.state != ProductGaugeCertificateState::ACTIVE ||
		!certificate.active || !certificate.currentAlignmentValid ||
		certificate.firstObservable != firstObservable ||
		certificate.secondObservable != secondObservable)
		return false;
	for (const auto& [satellite, coefficient] : certificate.productRelationQ)
	{
		if (coefficient == 0) continue;
		auto stored = certificate.phaseSegmentsBySatellite.find(satellite);
		auto current = currentSegments.find(satellite);
		if (stored == certificate.phaseSegmentsBySatellite.end() ||
			current == currentSegments.end() || stored->second != current->second)
			return false;
	}
	return true;
}

/** Segment identity check without requiring ACTIVE state.  Pending posterior
 * confirmation needs this exact physical guard before a CONFIRMING certificate
 * can be evaluated on a later posterior. */
inline bool zhangProductGaugeCertificateHasCurrentSegments(
	const ProductGaugeCertificate& certificate,
	E_ObsCode firstObservable,
	E_ObsCode secondObservable,
	const std::map<SatSys, std::array<std::string, 2>>& currentSegments)
{
	if (certificate.firstObservable != firstObservable ||
		certificate.secondObservable != secondObservable)
		return false;
	for (const auto& [satellite, coefficient] : certificate.productRelationQ)
	{
		if (coefficient == 0) continue;
		auto stored = certificate.phaseSegmentsBySatellite.find(satellite);
		auto current = currentSegments.find(satellite);
		if (stored == certificate.phaseSegmentsBySatellite.end() ||
			current == currentSegments.end() || stored->second != current->second)
			return false;
	}
	return true;
}

inline void zhangSetProductGaugeCertificateState(
	ProductGaugeCertificate& certificate,
	ProductGaugeCertificateState state,
	ProductGaugeCertificateInactiveReason inactiveReason =
		ProductGaugeCertificateInactiveReason::NONE)
{
	certificate.state = state;
	certificate.active = state == ProductGaugeCertificateState::ACTIVE;
	certificate.currentAlignmentValid = certificate.active &&
		certificate.temporalAlignmentCertified;
	certificate.inactiveReason = certificate.active
		? ProductGaugeCertificateInactiveReason::NONE : inactiveReason;
}

class ProductGaugeCertificateLedger
{
public:
	ProductGaugeCertificateLedgerUpdate observe(
		long epoch,
		const std::vector<ProductGaugeCertificate>& candidates,
		int requiredConfirmations)
	{
		ProductGaugeCertificateLedgerUpdate result;
		result.inputCertificates = candidates.size();
		if (epoch <= 0 || requiredConfirmations < 1)
		{
			result.failureReason = "PRODUCT_GAUGE_LEDGER_INPUT_INVALID";
			return result;
		}
		for (auto candidate : candidates)
		{
			zhangCanonicaliseProductGaugeCertificate(candidate);
			zhangInitialiseProductGaugeRiskReceipt(candidate);
			if (candidate.system == E_Sys::NONE ||
				candidate.satellite == candidate.reference)
			{
				result.failureReason = "PRODUCT_GAUGE_CERTIFICATE_SEGMENT_MISSING";
				return result;
			}
			if (!zhangProductGaugeCertificateSegmentsResolved(candidate))
			{
				result.rejectedCandidates++;
				result.failureReason =
					"PRODUCT_GAUGE_CERTIFICATE_SEGMENT_UNRESOLVED";
				return result;
			}
			const auto pairKey = zhangProductGaugeCertificatePairKey(candidate);
			const auto identity = zhangProductGaugeCertificateIdentity(candidate);
			for (auto& existing : certificates_)
			{
				if (zhangProductGaugeCertificatePairKey(existing) == pairKey &&
					zhangProductGaugeCertificateIdentity(existing) != identity)
				{
					zhangSetProductGaugeCertificateState(existing,
						ProductGaugeCertificateState::REVOKED,
						ProductGaugeCertificateInactiveReason::PHASE_SEGMENT_RETIRED);
					result.retiredSegmentCertificates++;
					result.revokedCertificates++;
					existing.validTo = epoch;
				}
			}
			if (!zhangProductGaugeCertificateEvidenceComplete(
				candidate, requiredConfirmations))
			{
				for (auto& existing : certificates_)
				if (zhangProductGaugeCertificateIdentity(existing) == identity &&
					existing.state != ProductGaugeCertificateState::REVOKED)
				{
					zhangSetProductGaugeCertificateState(existing,
						ProductGaugeCertificateState::SUSPENDED,
						ProductGaugeCertificateInactiveReason::ALIGNMENT_SUSPENDED);
					result.suspendedCertificates++;
				}
				result.rejectedCandidates++;
				continue;
			}
			auto existing = std::find_if(certificates_.rbegin(), certificates_.rend(),
				[&](const auto& certificate)
				{
					return zhangProductGaugeCertificateIdentity(certificate) == identity &&
						certificate.wideLaneInteger == candidate.wideLaneInteger &&
						certificate.firstSignalInteger == candidate.firstSignalInteger;
				});
			if (existing == certificates_.rend())
			{
				for (auto& historical : certificates_)
				if (zhangProductGaugeCertificateIdentity(historical) == identity &&
					historical.state != ProductGaugeCertificateState::REVOKED)
				{
					zhangSetProductGaugeCertificateState(historical,
						ProductGaugeCertificateState::SUSPENDED,
						ProductGaugeCertificateInactiveReason::CONFLICT);
					result.suspendedCertificates++;
					result.conflicts++;
				}
			candidate.firstCertified = epoch;
			candidate.validFrom = epoch;
			candidate.validTo = 0;
				candidate.lastConfirmed = epoch;
				candidate.confirmationEpochs = 1;
				candidate.componentVersion = ++version_;
				zhangSetProductGaugeCertificateState(candidate,
					requiredConfirmations == 1
						? ProductGaugeCertificateState::ACTIVE
						: ProductGaugeCertificateState::CONFIRMING);
				certificates_.push_back(std::move(candidate));
				if (certificates_.back().active)
					recordRiskActivation(certificates_.back());
				result.freshCertificates++;
				result.confirmedCertificates += certificates_.back().active;
				continue;
			}
			// reverse_iterator lets us retain the conflicting historical
			// certificate as SUSPENDED instead of overwriting its provenance.
			auto& confirmed = *existing;
			std::optional<ZhangProductGaugeTemporalEvidenceAudit> temporalAudit;
			if (confirmed.state != ProductGaugeCertificateState::ACTIVE &&
				confirmed.requiresIndependentTemporalEvidence &&
				candidate.requiresIndependentTemporalEvidence)
			{
				temporalAudit = zhangProductGaugeTemporalEvidenceIncrement(
					confirmed, candidate);
				result.maximumNewFactorCount = std::max(
					result.maximumNewFactorCount, temporalAudit->newFactorCount);
				result.maximumConditionalInformationRank = std::max(
					result.maximumConditionalInformationRank,
					temporalAudit->conditionalInformationRank);
				result.maximumOverlapFraction = std::max(
					result.maximumOverlapFraction, temporalAudit->overlapFraction);
				result.maximumConditionalInformationGain = std::max(
					result.maximumConditionalInformationGain,
					temporalAudit->conditionalInformationGain);
				if (!temporalAudit->valid)
				{
					result.insufficientInformationGainCandidates++;
					continue;
				}
				if (!temporalAudit->nestedFactorPrefix)
				{
					result.nonNestedEvidenceCandidates++;
					continue;
				}
				if (temporalAudit->newFactorCount <= 0)
				{
					result.duplicateEvidenceCandidates++;
					continue;
				}
				if (!temporalAudit->independentInformation)
				{
					result.insufficientInformationGainCandidates++;
					continue;
				}
				result.independentEvidenceConfirmations++;
			}
			for (const auto& sourceId : candidate.sourceCertificateIds)
				if (std::find(confirmed.sourceCertificateIds.begin(),
					confirmed.sourceCertificateIds.end(), sourceId) ==
					confirmed.sourceCertificateIds.end())
					confirmed.sourceCertificateIds.push_back(sourceId);
			confirmed.failureProbabilityBudget = std::min(
				confirmed.failureProbabilityBudget,
				candidate.failureProbabilityBudget);
			if (temporalAudit && temporalAudit->independentInformation)
			{
				confirmed.factorWindowIdentity = candidate.factorWindowIdentity;
				confirmed.factorSequences = candidate.factorSequences;
				confirmed.factorMeasurementRows = candidate.factorMeasurementRows;
				confirmed.temporalEvidenceCovariance =
					candidate.temporalEvidenceCovariance;
			}
			if (confirmed.lastConfirmed != epoch)
			{
				confirmed.lastConfirmed = epoch;
				confirmed.confirmationEpochs++;
			}
			const bool wasActive = confirmed.active;
			zhangSetProductGaugeCertificateState(confirmed,
				confirmed.confirmationEpochs >= requiredConfirmations
					? ProductGaugeCertificateState::ACTIVE
					: ProductGaugeCertificateState::CONFIRMING);
			confirmed.riskReceipt.isNewIntegerDecision = false;
			if (!wasActive && confirmed.active) recordRiskActivation(confirmed);
			result.confirmedCertificates += confirmed.active;
		}
		refreshComponentLineage();
		result.valid = true;
		result.failureReason = "NONE";
		return result;
	}

	const std::vector<ProductGaugeCertificate>& certificates() const
	{
		return certificates_;
	}

	double lifetimeBroadcastRisk() const
	{
		return lifetimeBroadcastRisk_;
	}

	bool invalidate(
		const ProductGaugeCertificate& rejected,
		long epoch,
		ProductGaugeCertificateInactiveReason reason =
			ProductGaugeCertificateInactiveReason::CONFLICT)
	{
		ProductGaugeCertificate canonical = rejected;
		zhangCanonicaliseProductGaugeCertificate(canonical);
		const std::string fingerprint =
			zhangProductGaugeCertificateFunctionalFingerprint(canonical);
		auto certificate = std::find_if(
			certificates_.rbegin(), certificates_.rend(), [&](const auto& stored)
			{
				return zhangProductGaugeCertificateFunctionalFingerprint(stored) ==
					fingerprint &&
					stored.state != ProductGaugeCertificateState::REVOKED;
			});
		if (certificate == certificates_.rend()) return false;
		zhangSetProductGaugeCertificateState(*certificate,
			ProductGaugeCertificateState::SUSPENDED, reason);
		certificate->validTo = epoch;
		certificate->confirmationEpochs = 0;
		refreshComponentLineage();
		return true;
	}

	/** Exact persistent satellite-product lattice rank.  Counting active rows
	 * overstates knowledge whenever redundant pair certificates are present. */
	int activeRank() const
	{
		std::set<SatSys> satellites;
		for (const auto& certificate : certificates_)
			if (certificate.state == ProductGaugeCertificateState::ACTIVE &&
				certificate.active && certificate.currentAlignmentValid)
				for (const auto& [satellite, coefficient] :
					certificate.productRelationQ)
					if (coefficient != 0) satellites.insert(satellite);
		std::map<SatSys, int> column;
		int index = 0;
		for (const auto& satellite : satellites) column[satellite] = index++;
		ZhangExactMatrix rows;
		for (const auto& certificate : certificates_)
		{
			if (certificate.state != ProductGaugeCertificateState::ACTIVE ||
				!certificate.active || !certificate.currentAlignmentValid) continue;
			ZhangExactVector row(satellites.size());
			for (const auto& [satellite, coefficient] :
				certificate.productRelationQ)
				if (coefficient != 0) row[column.at(satellite)] += coefficient;
			rows.push_back(std::move(row));
		}
		return static_cast<int>(
			zhangExactRowHermiteNormalForm(rows).basis.size());
	}

	int rawActiveCount() const
	{
		return std::count_if(certificates_.begin(), certificates_.end(),
			[](const auto& certificate)
			{
				return certificate.state == ProductGaugeCertificateState::ACTIVE &&
					certificate.active;
			});
	}

	std::string componentUuidFor(const SatSys& satellite) const
	{
		auto found = componentLineage_.find(satellite);
		return found == componentLineage_.end() ? std::string{} : found->second;
	}

private:
	void recordRiskActivation(ProductGaugeCertificate& certificate)
	{
		const auto& receipt = certificate.riskReceipt;
		if (receipt.certificateFamilyId.empty() ||
			!std::isfinite(receipt.currentCertificateRisk) ||
			receipt.currentCertificateRisk < 0 ||
			receipt.currentCertificateRisk > 1) return;
		auto [family, inserted] = activatedFamilyRisk_.try_emplace(
			receipt.certificateFamilyId, receipt.currentCertificateRisk);
		if (!inserted) family->second = std::max(
			family->second, receipt.currentCertificateRisk);
		lifetimeBroadcastRisk_ = 0;
		for (const auto& [_, probability] : activatedFamilyRisk_)
			lifetimeBroadcastRisk_ = std::min(
				1.0, lifetimeBroadcastRisk_ + probability);
		certificate.riskReceipt.lifetimeBroadcastRisk = lifetimeBroadcastRisk_;
	}

	void refreshComponentLineage()
	{
		std::map<SatSys, SatSys> parent;
		auto find = [&](auto&& self, SatSys satellite) -> SatSys
		{
			auto found = parent.find(satellite);
			if (found == parent.end()) return parent[satellite] = satellite;
			if (found->second == satellite) return satellite;
			return found->second = self(self, found->second);
		};
		std::set<SatSys> satellites;
		for (const auto& certificate : certificates_)
			if (certificate.active && certificate.currentAlignmentValid)
				for (const auto& [satellite, coefficient] :
					certificate.productRelationQ)
					if (coefficient != 0) satellites.insert(satellite);
		std::map<SatSys, int> column;
		int columnIndex = 0;
		for (const auto& satellite : satellites)
		{
			column[satellite] = columnIndex++;
			parent[satellite] = satellite;
		}
		ZhangExactMatrix rows;
		for (const auto& certificate : certificates_)
		{
			if (!certificate.active || !certificate.currentAlignmentValid) continue;
			ZhangExactVector row(satellites.size());
			for (const auto& [satellite, coefficient] :
				certificate.productRelationQ)
				if (coefficient != 0) row[column.at(satellite)] += coefficient;
			rows.push_back(std::move(row));
		}
		// Component membership is exact lattice membership of every primitive
		// e_s-e_p direction.  Merely co-occurring in a multi-satellite row does
		// not prove that the pair difference is integer-known.
		std::vector<std::pair<SatSys, SatSys>> satellitePairs;
		ZhangExactMatrix pairDifferences;
		for (auto left = satellites.begin(); left != satellites.end(); left++)
		for (auto right = std::next(left); right != satellites.end(); right++)
		{
			ZhangExactVector difference(satellites.size());
			difference[column.at(*left)] = 1;
			difference[column.at(*right)] = -1;
			satellitePairs.emplace_back(*left, *right);
			pairDifferences.push_back(std::move(difference));
		}
		const auto pairMemberships = zhangIntegerRowLatticeContainsBatch(
			rows, pairDifferences);
		for (std::size_t index = 0; index < satellitePairs.size(); index++)
		{
			if (!pairMemberships[index].contained) continue;
			const SatSys leftRoot = find(find, satellitePairs[index].first);
			const SatSys rightRoot = find(find, satellitePairs[index].second);
			if (leftRoot != rightRoot) parent[rightRoot] = leftRoot;
		}
		std::map<SatSys, std::vector<SatSys>> members;
		for (const auto& [satellite, ignored] : parent)
			members[find(find, satellite)].push_back(satellite);
		std::map<SatSys, std::string> lineage;
		for (auto& [root, component] : members)
		{
			std::sort(component.begin(), component.end());
			std::ostringstream stream;
			stream << "SATELLITE_LATTICE_COMPONENT:";
			for (std::size_t index = 0; index < component.size(); index++)
				stream << (index ? "," : "")
					<< static_cast<int>(component[index]);
			lineage[root] = stream.str();
		}
		componentLineage_.clear();
		for (const auto& satellite : satellites)
			componentLineage_[satellite] = lineage[find(find, satellite)];
		for (auto& certificate : certificates_)
		{
			if (!certificate.active || certificate.productRelationQ.empty()) continue;
			std::set<std::string> relationComponents;
			for (const auto& [satellite, coefficient] :
				certificate.productRelationQ)
				if (coefficient != 0)
					relationComponents.insert(componentLineage_[satellite]);
			std::ostringstream relationLineage;
			for (const auto& component : relationComponents)
				relationLineage << (relationLineage.tellp() > 0 ? "+" : "")
					<< component;
			certificate.componentUuid = relationLineage.str();
		}
	}

	std::vector<ProductGaugeCertificate> certificates_;
	std::map<SatSys, std::string> componentLineage_;
	std::map<std::string, double> activatedFamilyRisk_;
	double lifetimeBroadcastRisk_ = 0;
	int version_ = 0;
};

inline std::map<std::pair<std::string, E_Sys>, ProductGaugeCertificateLedger>&
zhangProductGaugeCertificateLedgerRegistry()
{
	static std::map<std::pair<std::string, E_Sys>,
		ProductGaugeCertificateLedger> registry;
	return registry;
}
