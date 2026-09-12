#pragma once


#include <algorithm>
#include <functional>
#include <initializer_list>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "common/eigenIncluder.hpp"
#include "common/zhangIntegerAudit.hpp"
#include "common/zhangIarGainAudit.hpp"
#include "common/zhangProductRelationBasis.hpp"
#include "common/zhangQuotientIntegerLattice.hpp"

/** Product closure is a private inference operation.  A fresh network-WL
 * fix is one admissible source of evidence, but it is not the only one: a
 * held physical lattice, a product-integer ledger, a frontend gauge
 * certificate, pending certified temporal work, or an available persistent
 * fixed-lag product marginal must also keep the closure branch alive.  The
 * latter is required for bootstrap: otherwise the raw marginal can never
 * create its first provisional Ledger certificate. */
// Exact preimage of the APPLIED physical affine lattice in product coordinates.
// Intersect whole lattices first: testing individual HNF rows loses combinations
// that cancel receiver coordinates. Never saturate or round a covariance mode.
inline ZhangExactRowHnf zhangAppliedLatticeProductImage(
	const ZhangExactMatrix& appliedRows, const ZhangExactVector& appliedValues,
	const ZhangExactMatrix& productTransform, const ZhangExactVector& offsets)
{
	ZhangExactRowHnf result;
	result.consistent = false;
	if (appliedRows.empty()) return zhangExactRowHermiteNormalForm({}, {});
	if (productTransform.empty() || productTransform.size() != offsets.size()) return result;
	const auto intersection = zhangExactAffineRowLatticeIntersection(
		appliedRows, appliedValues, productTransform,
		ZhangExactVector(productTransform.size()));
	if (!intersection.valid) return result;
	ZhangExactMatrix rows;
	ZhangExactVector values;
	const auto memberships = zhangIntegerRowLatticeContainsBatch(
		productTransform, intersection.rows);
	if (memberships.size() != intersection.rows.size()) return result;
	for (std::size_t i = 0; i < memberships.size(); ++i)
	{
		if (!memberships[i].contained || memberships[i].combination.size() != offsets.size())
			return result;
		auto value = intersection.firstValues[i];
		for (std::size_t j = 0; j < offsets.size(); ++j)
			value += memberships[i].combination[j] * offsets[j];
		rows.push_back(memberships[i].combination);
		values.push_back(value);
	}
	return zhangExactRowHermiteNormalForm(rows, values);
}

inline bool zhangShouldRunProductClosure(
	bool hasNewNetworkWl,
	bool hasHeldNetworkIntegers,
	bool hasActiveProductLedger,
	bool hasActiveGaugeCertificates,
	bool hasPendingTemporalWork,
	bool hasFixedLagProductWork = false)
{
	return hasNewNetworkWl
		|| hasHeldNetworkIntegers
		|| hasActiveProductLedger
		|| hasActiveGaugeCertificates
		|| hasPendingTemporalWork
		|| hasFixedLagProductWork;
}

/** Ordering used only to quarantine an inconsistent duplicate/cycle.  It is
 * intentionally not a reliability waiver: a lower-trust group is removed and
 * the remaining integer incidence block is solved again under the same gates. */
inline int zhangComponentEvidenceTrust(const std::string& source)
{
	if (source == "PRODUCT_GAUGE_CERTIFICATE") return 5;
	if (source == "CURRENT_GENERATION_EXACT_PAIR") return 4;
	if (source == "CROSS_GENERATION_PROJECTED_LEDGER") return 3;
	if (source == "TEMPORAL_RECERTIFIED") return 2;
	if (source == "BESD_CANDIDATE") return 1;
	return 0;
}

/** One exact satellite-pair coordinate in a named star ambient lattice.
 * Nodes [0,namedCount) are named satellite-minus-reference coordinates and
 * node namedCount is the canonical reference satellite. */
struct ZhangCertifiedPairRelation
{
	int firstNode = -1;
	int secondNode = -1;
	ZhangExactInteger value = 0;
	ZhangExactVector parentCombination;
	std::string coordinate = "UNKNOWN";
	// True only when at least one of the WL/L1 coordinates forming this
	// dual-frequency edge came from a previously certified ledger row that
	// survived the current private-branch joint-NIS admission.
	bool fromTemporalLedger = false;
	// A frontend product-gauge certificate has a different lifetime from a
	// physical-arc ledger row.  It is keyed by satellite phase segments and
	// remains valid across a backend S-basis change.
	bool fromProductGaugeLedger = false;
	// Set only by the final certificate-authorisation pass, after every exact
	// union and pair-catalogue rebuild has completed.  A historical ACTIVE flag
	// or a matching transient node index is never sufficient writer authority.
	bool currentPosteriorReauthorized = false;
	bool exactHistoricalTransportWitness = false;
	int productGaugeComponentVersion = 0;
	std::uint64_t authorizedBackendBasisGeneration = 0;
	std::string productGaugeSnapshotId;
	std::string validatedDeliveryMomentId;
	std::string authorizedPhaseSegmentFingerprint;
	std::vector<std::string> authorizedFamilyIds;
	std::vector<std::string> authorizedDecisionIds;
	double authorizedFailureProbability = 1;
};

/** Exact product-lattice constraints returned by the product IAR solver.
 *
 * Product rows are expressed in the complete mappable product coordinate,
 * not in a transient named-subset or LAMBDA-reduced suffix.  networkRows and
 * networkIntegers are the exact affine pull-back to the ambiguity coordinates
 * consumed by GinAR_mtx.  Mixed rows remain valid conditioning information
 * even when they contain no individually recoverable satellite-pair edge. */
#include "common/zhangIntegerDecisionProof.hpp"

// Bind each exact target row to the source decisions actually used in its
// integer representation. Never assign the union of unrelated source proofs
// to every HNF/projection row. Conditional parents of a used decision remain.
struct ZhangExactRowProofBinding
{
	bool valid = false;
	std::vector<ZhangDecisionProofs> rowProofs;
	std::string reason = "ROW_PROOF_INPUT_INVALID";
};

inline ZhangExactRowProofBinding zhangBindExactRowDecisionProofs(
	const ZhangExactMatrix& sourceRows, const ZhangExactVector& sourceValues,
	const std::vector<ZhangDecisionProofs>& sourceProofs,
	const ZhangExactMatrix& targetRows, const ZhangExactVector& targetValues,
	const ZhangExactMatrix* exactTransform = nullptr)
{
	ZhangExactRowProofBinding result;
	if (sourceRows.size() != sourceValues.size() || sourceRows.size() != sourceProofs.size() ||
		targetRows.size() != targetValues.size()) return result;
	if (targetRows.empty()) { result.valid=true; result.reason="EMPTY_TARGET"; return result; }
	if (sourceRows.empty()) return result;
	const auto dimension=sourceRows.front().size();
	for (const auto& row:sourceRows) if(row.size()!=dimension) return result;
	for (const auto& row:targetRows) if(row.size()!=dimension) return result;
	// HNF already provides its exact transform; reuse it rather than factor twice.
	std::vector<ZhangIntegerLatticeMembership> membership;
	if (exactTransform)
	{
		if (exactTransform->size()!=targetRows.size()) return result;
		for (const auto& combination:*exactTransform)
		{
			ZhangIntegerLatticeMembership item;
			item.contained=true; item.combination=combination;
			membership.push_back(std::move(item));
		}
	}
	else membership=zhangIntegerRowLatticeContainsBatch(sourceRows,targetRows);
	if (membership.size()!=targetRows.size()) return result;
	for(std::size_t r=0;r<targetRows.size();++r)
	{
		const auto& item=membership[r];
		if(!item.contained || item.combination.size()!=sourceRows.size())
		{ result.reason="ROW_PROOF_NOT_INTEGER_CONSEQUENCE"; return result; }
		ZhangExactInteger rhs=0;
		ZhangExactVector reconstructed(dimension);
		ZhangDecisionProofs proofs;
		for(std::size_t k=0;k<sourceRows.size();++k)
		{
			if(item.combination[k]==0) continue;
			if(sourceProofs[k].empty() || !zhangDecisionRiskClosure(sourceProofs[k]).valid)
			{ result.reason="ROW_PROOF_SOURCE_MISSING"; return result; }
			rhs+=item.combination[k]*sourceValues[k];
			for(std::size_t c=0;c<dimension;++c) reconstructed[c]+=item.combination[k]*sourceRows[k][c];
			proofs=zhangMergeDecisionProofs(proofs,sourceProofs[k]);
		}
		if(reconstructed!=targetRows[r] || rhs!=targetValues[r])
		{ result.reason="ROW_PROOF_AFFINE_MISMATCH"; return result; }
		if(!zhangDecisionRiskClosure(proofs).valid)
		{ result.reason="ROW_PROOF_DEPENDENCY_INVALID"; return result; }
		result.rowProofs.push_back(std::move(proofs));
	}
	result.valid=true;
	result.reason="EXACT_ROW_DEPENDENCIES_BOUND";
	return result;
}

struct ZhangProductIntegerConstraintSet
{
	ZhangDecisionProofs decisionProofs;
	bool reliable = false;
	bool exactNetworkMapping = false;
	E_Sys system = E_Sys::NONE;
	E_ObsCode firstObservable = E_ObsCode::NONE;
	E_ObsCode secondObservable = E_ObsCode::NONE;
	SatSys referenceSatellite;
	// Node i of every product row denotes coordinateSatellites[i] minus
	// referenceSatellite.  The reference itself is graph node N.
	std::vector<SatSys> coordinateSatellites;
	int productCoordinateDimension = 0;
	int networkAmbiguityDimension = 0;

	ZhangExactMatrix wideLaneProductRows;
	ZhangExactVector wideLaneIntegers;
	ZhangExactMatrix firstSignalProductRows;
	ZhangExactVector firstSignalIntegers;
	// A no-residual-DOF component-gauge solution is deliberately not a
	// constraint.  These rows are only observations for the independent
	// ProductGaugeCertificateLedger multi-epoch confirmation path.
	ZhangExactMatrix provisionalWideLaneProductRows;
	ZhangExactVector provisionalWideLaneIntegers;
	ZhangExactMatrix provisionalFirstSignalProductRows;
	ZhangExactVector provisionalFirstSignalIntegers;
	// Statistical admission belongs to the provisional parent lattice itself,
	// not to whichever direct product lattice happened to be accepted in the
	// same epoch.  The writer uses these fields when constructing square-block
	// observations for the multi-epoch gauge ledger.
	bool provisionalGaugeEvidenceValidated = false;
	double provisionalGaugeFailureProbability = 1;

	ZhangExactMatrix networkRows;
	ZhangExactVector networkIntegers;
	// Complete current [z1,z2] product coordinate pulled back to the network
	// ambiguity state.  This is intentionally separate from jointProductRows:
	// the latter contains only accepted direct constraints, whereas this matrix
	// is the exact coordinate chart used to re-express a physical Ledger row in
	// the current product lattice after a backend basis-generation change.
	ZhangExactMatrix fullJointProductNetworkRows;
	ZhangExactVector fullJointProductAffineOffsets;
	bool fullJointProductMappingExact = false;
	// Joint [z1,z2] rows corresponding one-to-one with networkRows.
	ZhangExactMatrix jointProductRows;
	// Immutable current physical-arc identities, populated by AMBRES after the
	// solver returns because only it owns the KF ambiguity map and arc versions.
	std::vector<std::map<std::string, ZhangExactInteger>> physicalNetworkRows;
	std::string phaseSegmentFingerprint;
	std::uint64_t backendBasisGeneration = 0;
	// Exact historical affine lattice selected by the gauge-ledger presearch.
	// It remains in the canonical [WL,L1] product chart until the final pair
	// catalogue exists.  The finaliser below recovers *all* exact pair
	// consequences from the whole lattice, so a basis row that is not itself a
	// pair does not incorrectly imply that a pair combination is unmappable.
	ZhangExactMatrix currentReauthorizedHistoricalDualRows;
	ZhangExactVector currentReauthorizedHistoricalDualValues;
	ZhangDecisionProofs currentReauthorizedHistoricalDecisionProofs;
	std::vector<std::string> currentReauthorizedHistoricalFamilyIds;
	std::string currentReauthorizedHistoricalSnapshotId;
	std::string currentReauthorizedSourceMomentId;
	double currentReauthorizedHistoricalFailureProbability = 1;
	int currentReauthorizedHistoricalComponentVersion = 0;
	bool currentReauthorizedHistoricalSegmentsValidated = false;
	bool currentReauthorizedHistoricalExactTransport = false;

	std::vector<ZhangCertifiedPairRelation> certifiedPairs;
	// Pair relations proved by both WL and conditional-L1.  These, and only
	// these, are allowed to form the dual-frequency broadcast certificate graph.
	std::vector<ZhangCertifiedPairRelation> dualFrequencyCertifiedPairs;
	ZhangExactMatrix conditioningOnlyRows;
	ZhangExactVector conditioningOnlyIntegers;

	int conditioningRank = 0;
	int certifiedPairRank = 0;
	int pairRecoveryBatchTargets = 0;
	int pairRecoveryScalarDecompositionsAvoided = 0;
	double jointNis = std::numeric_limits<double>::quiet_NaN();
	double jointNisThreshold = std::numeric_limits<double>::quiet_NaN();
	double failureProbability = 1;
	double referenceInvariantProductGain = 0;
	std::string failureReason = "NOT_EVALUATED";
};

enum class ZhangProductCandidateFirstFailure
{
	NONE,
	EXACT_AFFINE_CONFLICT,
	CURRENT_MAPPING_FAIL,
	NONPRIMITIVE_OR_CONGRUENCE_FAIL,
	RISK_PARENT_INCOMPLETE,
	RISK_BUDGET_EXCEEDED,
	DETERMINISTIC_RESIDUAL_FAIL,
	CURRENT_NIS_FAIL,
	ZERO_PAIR_IMAGE,
	SEARCH_PRUNED_NOT_TESTED
};

inline const char* zhangProductCandidateFirstFailureName(
	ZhangProductCandidateFirstFailure reason)
{
	switch (reason)
	{
		case ZhangProductCandidateFirstFailure::NONE:
			return "NONE";
		case ZhangProductCandidateFirstFailure::EXACT_AFFINE_CONFLICT:
			return "EXACT_AFFINE_CONFLICT";
		case ZhangProductCandidateFirstFailure::CURRENT_MAPPING_FAIL:
			return "CURRENT_MAPPING_FAIL";
		case ZhangProductCandidateFirstFailure::NONPRIMITIVE_OR_CONGRUENCE_FAIL:
			return "NONPRIMITIVE_OR_CONGRUENCE_FAIL";
		case ZhangProductCandidateFirstFailure::RISK_PARENT_INCOMPLETE:
			return "RISK_PARENT_INCOMPLETE";
		case ZhangProductCandidateFirstFailure::RISK_BUDGET_EXCEEDED:
			return "RISK_BUDGET_EXCEEDED";
		case ZhangProductCandidateFirstFailure::DETERMINISTIC_RESIDUAL_FAIL:
			return "DETERMINISTIC_RESIDUAL_FAIL";
		case ZhangProductCandidateFirstFailure::CURRENT_NIS_FAIL:
			return "CURRENT_NIS_FAIL";
		case ZhangProductCandidateFirstFailure::ZERO_PAIR_IMAGE:
			return "ZERO_PAIR_IMAGE";
		case ZhangProductCandidateFirstFailure::SEARCH_PRUNED_NOT_TESTED:
			return "SEARCH_PRUNED_NOT_TESTED";
	}
	return "UNKNOWN";
}

/** Preserve the admission provenance of exact dual-frequency edges when an
 * HNF union rebuilds the certified pair catalogue.  Rebuilding the catalogue
 * is algebraically correct but otherwise erases whether an identical edge was
 * supplied by a current-posterior-admitted physical or gauge ledger.  The
 * writer uses these flags only for lineage and may never infer them again from
 * ledger ACTIVE state. */
inline void zhangPropagateProductPairProvenance(
	ZhangProductIntegerConstraintSet& output,
	std::initializer_list<const ZhangProductIntegerConstraintSet*> sources)
{
	using Edge = std::pair<int, int>;
	std::map<Edge, std::pair<bool, bool>> provenance;
	for (const auto* source : sources)
	{
		if (!source) continue;
		for (const auto& pair : source->dualFrequencyCertifiedPairs)
		{
			if (pair.firstNode < 0 || pair.secondNode < 0 ||
				pair.firstNode == pair.secondNode) continue;
			const Edge edge{std::min(pair.firstNode, pair.secondNode),
				std::max(pair.firstNode, pair.secondNode)};
			auto& flags = provenance[edge];
			flags.first |= pair.fromTemporalLedger;
			flags.second |= pair.fromProductGaugeLedger;
		}
	}
	for (auto& pair : output.dualFrequencyCertifiedPairs)
	{
		if (pair.firstNode < 0 || pair.secondNode < 0 ||
			pair.firstNode == pair.secondNode) continue;
		const Edge edge{std::min(pair.firstNode, pair.secondNode),
			std::max(pair.firstNode, pair.secondNode)};
		auto source = provenance.find(edge);
		if (source == provenance.end()) continue;
		pair.fromTemporalLedger |= source->second.first;
		pair.fromProductGaugeLedger |= source->second.second;
	}
}

/** Result boundary for the direct satellite-product integer solver.
 *
 * This object deliberately separates structural estimability, statistical
 * reliability and frontend admission.  A positive network fixed rank alone
 * can never set certifiedForProduct. */
struct ZhangR47Candidate;
struct ZhangProductRelationFixResult
{
	std::shared_ptr<const ZhangR47Candidate> r47Candidate;
    bool r47ControlledSearch=false;
    ZhangExactMatrix r47AdmittedRows,r47NewRows;
    ZhangExactVector r47AdmittedValues,r47NewValues;
	bool basisValid = false;
	bool mappingValid = false;
	bool wideLaneReliable = false;
	bool firstSignalReliable = false;
	bool certifiedForProduct = false;
	// The ProductIntegerLedger posterior and its exact/risk receipts form one
	// transaction.  When the receipts cannot be authorised, the caller must
	// discard that posterior and retry from the pre-Ledger search state.
	bool ledgerPosteriorRollbackRequired = false;
	int fullTargetRank = 0;
	int mappableTargetRank = 0;
	int wideLaneFixedRank = 0;
	int firstSignalFixedRank = 0;
	int namedFirstSignalFixed = 0;
	int namedSecondSignalFixed = 0;
	int evaluatedBranches = 0;
	int namedRoundWideLaneCandidates = 0;
	int namedRoundWideLaneRetained = 0;
	int selectedRawPartialFixedRank = 0;
	int selectedRecoveredNamedRank = 0;
	int selectedParentBranchRank = 0;
	double selectedPartialFixFraction = 0;
	int componentCoverageGain = 0;
	int certifiedJointIntegerRank = 0;
	// Private component-gauge diagnostics.  These are copied out of the
	// shadow/private solver so a frozen-posterior closure audit can report the
	// dual graph rank and component count at every iteration.
	int componentGaugeComponentsBefore = 0;
	int componentGaugeComponentsAfter = 0;
	int componentGaugeNewDualGraphRank = 0;
	double maximumWideLanePerr = 1;
	double maximumWideLaneMarginalRoundPerr = 1;
	double wideLaneParentFailureProbabilityBound = 1;
	double maximumFirstSignalPerr = 1;
	double productInformationGain = 0;
	double realSubspaceUpperBoundAtSelectedRank =
		std::numeric_limits<double>::quiet_NaN();
	// Ratio to the unconstrained real-subspace relaxation.  This is not an
	// integer-search efficiency until an integer-constrained frontier exists.
	double relaxedRealUpperBoundCapture =
		std::numeric_limits<double>::quiet_NaN();
	double realIntegerGainGap =
		std::numeric_limits<double>::quiet_NaN();
	double relationRho5 = std::numeric_limits<double>::quiet_NaN();
	double relationRho10 = std::numeric_limits<double>::quiet_NaN();
	double relationRho20 = std::numeric_limits<double>::quiet_NaN();
	double relationRho40 = std::numeric_limits<double>::quiet_NaN();
	double relationRho80 = std::numeric_limits<double>::quiet_NaN();
	int realSubspaceRank80 = 0;
	int realSubspaceRank90 = 0;
	int realSubspaceRank95 = 0;
	std::string gainSpectrumDiagnosis = "NOT_EVALUATED";
	double jointNis = std::numeric_limits<double>::quiet_NaN();
	double jointNisThreshold = std::numeric_limits<double>::quiet_NaN();
	std::map<std::size_t, ZhangExactInteger> namedWideLane;
	std::map<std::size_t, ZhangExactInteger> namedFirstSignal;
	std::map<std::size_t, ZhangExactInteger> namedSecondSignal;
	std::vector<int> selectedNamedRelationIndices;
	std::string selectedCanonicalHnf;
	std::string wideLaneCertificateSource = "NONE";
	// Provenance of accepted fixed-lag raw-factor rows.  These flags describe
	// rows that survived exact union and the current-posterior joint-NIS gate;
	// merely proposing a raw-factor integer never sets either flag.
	bool fixedLagWideLaneContribution = false;
	bool fixedLagFirstSignalContribution = false;
	bool namedOrderingValid = false;
	bool namedSubsetCertificate = false;
	ZhangProductIntegerConstraintSet constraints;
	std::string status = "NOT_EVALUATED";
	std::string failureReason = "NONE";
};

/** Private dual-frequency closure of the remaining component gauges.
 *
 * Direct product IAR first certifies within-component satellite differences.
 * This result represents only the K-1 datum-free gauges between those
 * components.  WL is resolved first; L1 is resolved conditional on the
 * accepted WL lattice.  The rows remain product-lattice rows until the caller
 * completes the final exact union and one-shot PRODUCT_FIXED conditioning. */
struct ZhangDualComponentGaugeFixResult
{
	bool valid = false;
	bool reliable = false;
	int componentsBefore = 0;
	int componentsAfter = 0;
	int gaugeTargetRank = 0;
	int wlGaugeFixedRank = 0;
	int firstSignalGaugeFixedRank = 0;
	int combinedCertifiedRank = 0;
	int newDualGraphRank = 0;
	int measurementRank = 0;
	int estimableGaugeRank = 0;
	std::vector<int> estimableGaugeColumns;
	int residualDof = 0;
	bool confirmationRequired = false;
	double residualNis = std::numeric_limits<double>::quiet_NaN();
	double residualNisThreshold = std::numeric_limits<double>::quiet_NaN();
	double maximumNullResidual = std::numeric_limits<double>::quiet_NaN();
	double wlFailureProbability = 1;
	double firstSignalFailureProbability = 1;
	ZhangExactMatrix wideLaneProductRows;
	ZhangExactVector wideLaneIntegers;
	ZhangExactMatrix firstSignalProductRows;
	ZhangExactVector firstSignalIntegers;
	std::vector<ZhangCertifiedPairRelation> newDualPairCertificates;
	std::string failureReason = "NOT_EVALUATED";
};

/** One independently admitted connected component-gauge block.
 *
 * A block never borrows an integer decision from another block.  Its rows are
 * exact product-lattice rows and may therefore be unioned by the caller
 * without converting a real SVD/QR direction into an integer constraint.
 */
struct ZhangComponentGaugeBlockResult
{
	std::vector<int> componentIds;
	int targetRank = 0;
	int estimableRank = 0;
	int measurementRank = 0;
	double residualNis = std::numeric_limits<double>::quiet_NaN();
	double residualNisThreshold = std::numeric_limits<double>::quiet_NaN();
	double maximumNullResidual = std::numeric_limits<double>::quiet_NaN();
	ZhangExactMatrix wlRows;
	ZhangExactVector wlIntegers;
	ZhangExactMatrix l1Rows;
	ZhangExactVector l1Integers;
	int newDualRank = 0;
	// valid is distinct from reliable.  A square gauge block has no residual
	// degrees of freedom and therefore cannot be authoritative in its first
	// epoch, but its exact dual-frequency rows are still valid provisional
	// evidence for the multi-epoch ProductGaugeCertificateLedger.
	bool valid = false;
	bool reliable = false;
	bool confirmationRequired = false;
	double wlFailureProbability = 1;
	double l1FailureProbability = 1;
	std::string failureReason = "NOT_EVALUATED";
	// Only rows admitted by this block.  This is intentionally distinct from
	// conditioning-only ledger rows and from the caller's final global union.
	ZhangProductIntegerConstraintSet constraints;
};

inline bool zhangComponentGaugeBlockRowsComplete(
	const ZhangComponentGaugeBlockResult& block)
{
	return block.targetRank > 0 && block.newDualRank > 0 &&
		!block.wlRows.empty() &&
		block.wlRows.size() == block.wlIntegers.size() &&
		!block.l1Rows.empty() &&
		block.l1Rows.size() == block.l1Integers.size();
}

inline bool zhangComponentGaugeBlockProvidesAuthoritativeEvidence(
	const ZhangComponentGaugeBlockResult& block)
{
	return block.valid && block.reliable && !block.confirmationRequired &&
		zhangComponentGaugeBlockRowsComplete(block);
}

inline bool zhangComponentGaugeBlockProvidesProvisionalEvidence(
	const ZhangComponentGaugeBlockResult& block)
{
	return block.valid && !block.reliable && block.confirmationRequired &&
		zhangComponentGaugeBlockRowsComplete(block);
}

/** Select independently admitted component blocks under one global family
 * error budget.
 *
 * Each block solver already guards its own integer lattice.  Treating every
 * block as though it owned the complete configured Perr budget would however
 * exceed that budget when several disconnected blocks succeed in one epoch.
 * This small exact Pareto dynamic program maximises the eventual dual graph
 * rank, then the rank that is authoritative immediately, and finally chooses
 * the smallest union-bound probability.  It also gives the caller one list
 * with which to keep authoritative and provisional rows strictly separated. */
struct ZhangComponentGaugeEvidenceSelection
{
	std::vector<int> selectedBlockIndices;
	int evidenceBlocks = 0;
	int budgetEligibleBlocks = 0;
	int selectedAuthoritativeBlocks = 0;
	int selectedProvisionalBlocks = 0;
	int selectedTotalRank = 0;
	int selectedAuthoritativeRank = 0;
	double familyFailureProbability = 0;
};

inline ZhangComponentGaugeEvidenceSelection
zhangSelectComponentGaugeEvidenceBlocks(
	const std::vector<ZhangComponentGaugeBlockResult>& blocks,
	double maximumFamilyFailureProbability)
{
	ZhangComponentGaugeEvidenceSelection result;
	if (!std::isfinite(maximumFamilyFailureProbability) ||
		maximumFamilyFailureProbability <= 0) return result;
	const double budget = std::min(1.0, maximumFamilyFailureProbability);
	const double tolerance = std::max(1e-15, 1e-12 * budget);

	struct State
	{
		double probability = 0;
		std::vector<int> selected;
	};
	// score = {eventual dual graph rank, immediately authoritative rank}.
	std::map<std::pair<int, int>, State> frontier;
	frontier[{0, 0}] = {};
	for (int index = 0; index < static_cast<int>(blocks.size()); index++)
	{
		const auto& block = blocks[index];
		const bool authoritative =
			zhangComponentGaugeBlockProvidesAuthoritativeEvidence(block);
		const bool provisional =
			zhangComponentGaugeBlockProvidesProvisionalEvidence(block);
		if (!authoritative && !provisional) continue;
		result.evidenceBlocks++;
		if (!std::isfinite(block.wlFailureProbability) ||
			!std::isfinite(block.l1FailureProbability) ||
			block.wlFailureProbability < 0 ||
			block.l1FailureProbability < 0) continue;
		const double probability =
			block.wlFailureProbability + block.l1FailureProbability;
		if (!std::isfinite(probability) || probability > budget + tolerance)
			continue;
		result.budgetEligibleBlocks++;

		auto next = frontier;
		for (const auto& [score, state] : frontier)
		{
			const double combinedProbability = state.probability + probability;
			if (combinedProbability > budget + tolerance) continue;
			const std::pair<int, int> combinedScore{
				score.first + block.newDualRank,
				score.second + (authoritative ? block.newDualRank : 0)};
			State candidate = state;
			candidate.probability = combinedProbability;
			candidate.selected.push_back(index);
			auto existing = next.find(combinedScore);
			if (existing == next.end() ||
				candidate.probability < existing->second.probability - tolerance ||
				(std::abs(candidate.probability -
					existing->second.probability) <= tolerance &&
				 candidate.selected < existing->second.selected))
				next[combinedScore] = std::move(candidate);
		}
		frontier = std::move(next);
	}

	std::pair<int, int> bestScore{0, 0};
	State best;
	for (const auto& [score, state] : frontier)
	{
		if (score > bestScore ||
			(score == bestScore && state.probability < best.probability - tolerance))
		{
			bestScore = score;
			best = state;
		}
	}
	result.selectedBlockIndices = std::move(best.selected);
	result.selectedTotalRank = bestScore.first;
	result.selectedAuthoritativeRank = bestScore.second;
	result.familyFailureProbability = best.probability;
	for (const int index : result.selectedBlockIndices)
	{
		if (zhangComponentGaugeBlockProvidesAuthoritativeEvidence(blocks[index]))
			result.selectedAuthoritativeBlocks++;
		else
			result.selectedProvisionalBlocks++;
	}
	return result;
}

/** Provenance for one component support edge.  Integer rows are always built
 * from this incidence provenance; real null vectors are diagnostic only. */
struct ZhangComponentEdgeId
{
	int firstComponent = -1;
	int secondComponent = -1;
	int firstNode = -1;
	int secondNode = -1;
	// Physical support is optional, but when the exact product-to-network map
	// is available it lets the forest prefer genuinely independent bridges
	// instead of repeatedly using the same receiver/ambiguity arc backbone.
	std::set<std::string> receiverSupport;
	std::set<std::string> arcSupport;
};

/** A failed covariance-nullspace consistency equation.  The coefficients are
 * a locator for quarantine, never an integer row to be fixed. */
struct ZhangNullSpaceConflict
{
	double nullResidual = std::numeric_limits<double>::quiet_NaN();
	std::vector<ZhangComponentEdgeId> dominantEdges;
	// c_k,i = u_0,k,i * r_i.  These are the quantities used for local
	// quarantine, rather than the bare null-vector coefficients.
	std::vector<double> contributions;
	// Kept for trace compatibility only; never use this field for ranking.
	std::vector<double> coefficients;
	std::vector<E_ObsCode> signals;
	std::vector<uint64_t> backendGenerations;
	std::vector<std::string> phaseSegments;
	std::vector<std::string> sources;
	std::vector<ZhangExactInteger> affineOffsets;
};

using ZhangNullConflict = ZhangNullSpaceConflict;

/** One product-aware PAR beam state.
 *
 * The candidate coordinates are always original named satellite relations.
 * LAMBDA may decorrelate internally, but a branch is certifiable only when
 * its fixed integer row lattice recovers every named coordinate in this list.
 */
struct ProductParBranch
{
	std::vector<int> namedRelationIndices;
	int integerRank = 0;
	bool reliabilityPassed = false;
	int rawPartialFixedRank = 0;
	int recoveredNamedRank = 0;
	int parentBranchRank = 0;
	bool inheritedFromParentFixedLattice = false;
	double partialFixFraction = 0;
	double maximumNamedPerr = 1;
	double normalizedCandidateNis = std::numeric_limits<double>::infinity();
	double maxPerr = 1;
	double jointNis = std::numeric_limits<double>::quiet_NaN();
	double jointNisThreshold = std::numeric_limits<double>::quiet_NaN();
	int componentCoverageGain = 0;
	double productInformationGain = 0;
	std::string canonicalHnf;
};

/** Lexicographic product-PAR score.  Reliability is an admission gate, not a
 * weighted reward: an unreliable candidate always ranks below every reliable
 * one, regardless of coverage or gain. */
struct ZhangProductRelationLexicographicScore
{
	bool reliabilityPassed = false;
	double partialFixFraction = 0;
	int rawPartialFixedRank = 0;
	int recoveredNamedRank = 0;
	double normalizedCandidateNis = std::numeric_limits<double>::infinity();
	double maximumNamedPerr = 1;
	int componentCoverageGain = 0;
	double productInformationGain = 0;

	bool operator<(const ZhangProductRelationLexicographicScore& other) const
	{
		if (reliabilityPassed != other.reliabilityPassed)
		{
			return reliabilityPassed < other.reliabilityPassed;
		}
		if (!reliabilityPassed)
		{
			if (recoveredNamedRank != other.recoveredNamedRank)
			{
				return recoveredNamedRank < other.recoveredNamedRank;
			}
			if (partialFixFraction != other.partialFixFraction)
			{
				return partialFixFraction < other.partialFixFraction;
			}
			if (rawPartialFixedRank != other.rawPartialFixedRank)
			{
				return rawPartialFixedRank < other.rawPartialFixedRank;
			}
			const double leftNis = std::isfinite(normalizedCandidateNis)
				? normalizedCandidateNis : std::numeric_limits<double>::infinity();
			const double rightNis = std::isfinite(other.normalizedCandidateNis)
				? other.normalizedCandidateNis : std::numeric_limits<double>::infinity();
			if (leftNis != rightNis)
			{
				return leftNis > rightNis;
			}
			if (maximumNamedPerr != other.maximumNamedPerr)
			{
				return maximumNamedPerr > other.maximumNamedPerr;
			}
		}
		if (componentCoverageGain != other.componentCoverageGain)
		{
			return componentCoverageGain < other.componentCoverageGain;
		}
		return productInformationGain < other.productInformationGain;
	}
};

inline ZhangProductRelationLexicographicScore zhangProductParScore(
	const ProductParBranch& branch)
{
	return {
		branch.reliabilityPassed,
		branch.partialFixFraction,
		branch.rawPartialFixedRank,
		branch.recoveredNamedRank,
		branch.normalizedCandidateNis,
		branch.maximumNamedPerr,
		branch.componentCoverageGain,
		branch.productInformationGain
	};
}

/** One named satellite edge proposed to the product-aware forest beam. */
struct ZhangNamedPairBeamCandidate
{
	ZhangExactVector row;
	double perr = 1;
	double gain = 0;
	double variance = std::numeric_limits<double>::infinity();
	std::vector<int> nodes;
	int dualGraphRankGain = 0;
	int productLatticeRankGain = 0;
};

/** One primitive component-incidence integer candidate in a datum-free gauge
 * coordinate.  Component zero is the implicit datum; coordinates 0..K-2 are
 * components 1..K-1.  Dense decorrelated directions are deliberately absent
 * from this catalogue, so callers can use them for conditioning diagnostics
 * without ever promoting them to graph certificates. */
struct ZhangPrimitiveComponentGaugeCandidate
{
	ZhangExactVector row;
	int firstComponent = -1;
	int secondComponent = -1;
	double floatingValue = std::numeric_limits<double>::quiet_NaN();
	double variance = std::numeric_limits<double>::quiet_NaN();
	double fractional = std::numeric_limits<double>::quiet_NaN();
};

/** Enumerate every finite, stochastic primitive e_a-e_b gauge edge.
 *
 * When allowedPairLattice is non-empty, only primitive edges contained in that
 * exact lattice are returned.  This is used for L1-given-WL: an independently
 * fixed L1 edge cannot become a dual-frequency graph certificate unless the
 * same named component edge is already an exact consequence of the accepted
 * WL lattice. */
inline std::vector<ZhangPrimitiveComponentGaugeCandidate>
zhangPrimitiveComponentGaugeCandidates(
	const Eigen::VectorXd& mean,
	const Eigen::MatrixXd& covariance,
	int componentCount,
	const ZhangExactMatrix& allowedPairLattice = {})
{
	std::vector<ZhangPrimitiveComponentGaugeCandidate> result;
	const int gauges = componentCount - 1;
	if (componentCount <= 1 || mean.size() != gauges ||
		covariance.rows() != gauges || covariance.cols() != gauges ||
		!mean.allFinite() || !covariance.allFinite()) return result;

	std::vector<std::pair<int, int>> nodes;
	ZhangExactMatrix rows;
	for (int first = 0; first < componentCount; first++)
	for (int second = first + 1; second < componentCount; second++)
	{
		ZhangExactVector row(gauges);
		if (first > 0) row[first - 1] += 1;
		if (second > 0) row[second - 1] -= 1;
		nodes.push_back({first, second});
		rows.push_back(std::move(row));
	}
	std::vector<ZhangIntegerLatticeMembership> memberships;
	if (!allowedPairLattice.empty())
		memberships = zhangIntegerRowLatticeContainsBatch(
			allowedPairLattice, rows);
	const double varianceTolerance = std::max(1e-20,
		1e-14 * std::max(1.0, covariance.norm()));
	for (int index = 0; index < static_cast<int>(rows.size()); index++)
	{
		if (!allowedPairLattice.empty() &&
			(memberships.size() != rows.size() ||
			 !memberships[index].contained)) continue;
		const Eigen::VectorXd numeric = zhangExactRowToDouble(rows[index]);
		const double value = numeric.dot(mean);
		const double variance =
			(numeric.transpose() * covariance * numeric)(0, 0);
		if (!std::isfinite(value) || !std::isfinite(variance) ||
			variance <= varianceTolerance) continue;
		result.push_back({rows[index], nodes[index].first, nodes[index].second,
			value, variance, value - std::round(value)});
	}
	return result;
}

struct ZhangNamedPairBeamBranch
{
	std::vector<int> selected;
	double maximumPerr = 0;
	double conditionalProductGain = 0;
	double summedGain = 0;
	double summedVariance = 0;
	int coveredNodes = 0;
	int dualGraphRankGain = 0;
	int productLatticeRankGain = 0;
};

/** Count only graph-rank gained beyond an already certified dual-frequency
 * component partition.  Summing standalone edge flags is incorrect: two
 * individually useful edges can close the same component cut and contribute
 * only one independent graph relation when selected together. */
inline int zhangNamedEdgeGraphRankGain(
	int nodeCount,
	const std::vector<int>& existingComponentLabels,
	const std::vector<std::vector<int>>& edges)
{
	if (nodeCount <= 0) return 0;
	std::vector<int> parent(nodeCount);
	for (int node = 0; node < nodeCount; node++) parent[node] = node;
	auto root = [&](int node)
	{
		while (parent[node] != node)
		{
			parent[node] = parent[parent[node]];
			node = parent[node];
		}
		return node;
	};
	auto unite = [&](int first, int second)
	{
		if (first < 0 || second < 0 || first >= nodeCount || second >= nodeCount)
			return false;
		first = root(first);
		second = root(second);
		if (first == second) return false;
		parent[second] = first;
		return true;
	};
	if (existingComponentLabels.size() == static_cast<std::size_t>(nodeCount))
	{
		std::map<int, int> representative;
		for (int node = 0; node < nodeCount; node++)
		{
			const int label = existingComponentLabels[node];
			auto [iterator, inserted] = representative.try_emplace(label, node);
			if (!inserted) unite(iterator->second, node);
		}
	}
	int gain = 0;
	for (const auto& edge : edges)
	{
		if (edge.size() != 2) continue;
		gain += unite(edge[0], edge[1]);
	}
	return gain;
}

/** Retain alternate primitive named-edge forests at every rank.
 *
 * This is deliberately an expansion beam, not deletion from one greedy tree:
 * an edge excluded from the first high-gain forest remains reachable through
 * another branch.  Ordering is the requested reliability -> satellite
 * dual-frequency graph gain -> persistent-quotient rank gain -> satellite
 * coverage -> product gain
 * lexicographic order. */
inline std::vector<std::vector<ZhangNamedPairBeamBranch>>
zhangNamedPairForestBeamLevels(
	const std::vector<ZhangNamedPairBeamCandidate>& candidates,
	int dimension,
	int beamWidth,
	const std::function<double(const std::vector<int>&)>& conditionalGain = {},
	int maximumExpansionsPerLevel = std::numeric_limits<int>::max(),
	bool* expansionCapped = nullptr,
	const ZhangExactMatrix& existingSignalRows = {},
	const std::vector<int>& existingDualComponentLabels = {})
{
	std::vector<std::vector<ZhangNamedPairBeamBranch>> levels;
	if (expansionCapped) *expansionCapped = false;
	if (dimension <= 0 || beamWidth <= 0 || candidates.empty()) return levels;
	const int existingRank = static_cast<int>(
		zhangExactRowHermiteNormalForm(existingSignalRows).basis.size());
	auto quality = [&](const std::vector<int>& selected)
	{
		ZhangNamedPairBeamBranch branch;
		branch.selected = selected;
		std::set<int> nodes;
		ZhangExactMatrix unitedRows = existingSignalRows;
		std::vector<std::vector<int>> dualEdges;
		for (int index : selected)
		{
			if (index < 0 || index >= static_cast<int>(candidates.size()))
				return ZhangNamedPairBeamBranch{};
			branch.maximumPerr = std::max(
				branch.maximumPerr, candidates[index].perr);
			branch.summedGain += candidates[index].gain;
			branch.summedVariance += candidates[index].variance;
			if (candidates[index].dualGraphRankGain > 0)
				dualEdges.push_back(candidates[index].nodes);
			unitedRows.push_back(candidates[index].row);
			nodes.insert(candidates[index].nodes.begin(),
				candidates[index].nodes.end());
		}
		branch.coveredNodes = nodes.size();
		branch.dualGraphRankGain = zhangNamedEdgeGraphRankGain(
			dimension + 1, existingDualComponentLabels, dualEdges);
		branch.productLatticeRankGain = std::max(0,
			static_cast<int>(zhangExactRowHermiteNormalForm(
				unitedRows).basis.size()) - existingRank);
		branch.conditionalProductGain = conditionalGain
			? conditionalGain(selected)
			: branch.summedGain;
		return branch;
	};
	auto better = [](const auto& left, const auto& right)
	{
		// Every edge entering this beam has already passed the scalar Perr/NIS
		// admission gate.  Ranking again by the continuous Perr value before the
		// exact objectives starves a slightly less overconfident edge that closes
		// a dual-frequency or persistent-quotient direction.  Exact completion
		// therefore leads; Perr remains the safety-margin tie breaker.
		if (left.dualGraphRankGain != right.dualGraphRankGain)
			return left.dualGraphRankGain > right.dualGraphRankGain;
		if (left.productLatticeRankGain != right.productLatticeRankGain)
			return left.productLatticeRankGain >
				right.productLatticeRankGain;
		if (left.coveredNodes != right.coveredNodes)
			return left.coveredNodes > right.coveredNodes;
		if (left.maximumPerr != right.maximumPerr)
			return left.maximumPerr < right.maximumPerr;
		if (left.conditionalProductGain != right.conditionalProductGain)
			return left.conditionalProductGain > right.conditionalProductGain;
		if (left.summedVariance != right.summedVariance)
			return left.summedVariance < right.summedVariance;
		return left.selected < right.selected;
	};

	std::vector<ZhangNamedPairBeamBranch> frontier = {
		ZhangNamedPairBeamBranch{}};
	for (int targetRank = 1;
		targetRank <= dimension && !frontier.empty(); targetRank++)
	{
		std::map<std::vector<int>, ZhangNamedPairBeamBranch> unique;
		int expanded = 0;
		bool cappedThisLevel = false;
		for (const auto& parent : frontier)
		for (int index = 0; index < static_cast<int>(candidates.size()); index++)
		{
			// The LAMBDA evaluation cap must also bound the preceding exact
			// forest construction.  Otherwise a large all-pair dictionary can
			// allocate and score thousands of temporary HNF/gain candidates before
			// the later evaluation loop ever sees its cap.
			if (expanded >= maximumExpansionsPerLevel)
			{
				cappedThisLevel = true;
				break;
			}
			expanded++;
			if (std::binary_search(
				parent.selected.begin(), parent.selected.end(), index)) continue;
			auto selected = parent.selected;
			selected.push_back(index);
			std::sort(selected.begin(), selected.end());
			ZhangExactMatrix rows;
			for (int selectedIndex : selected)
				rows.push_back(candidates[selectedIndex].row);
			int exactRank = 0;
			if (!zhangExactPrimitiveRowLattice(rows, dimension, &exactRank) ||
				exactRank != static_cast<int>(rows.size())) continue;
			unique.try_emplace(selected, quality(selected));
		}
		if (cappedThisLevel && expansionCapped) *expansionCapped = true;
		frontier.clear();
		for (auto& [selected, branch] : unique)
			frontier.push_back(std::move(branch));
		std::sort(frontier.begin(), frontier.end(), better);
		if (frontier.size() > static_cast<std::size_t>(beamWidth))
			frontier.resize(beamWidth);
		if (!frontier.empty()) levels.push_back(frontier);
	}
	return levels;
}

/** Prove that L1 and L2 use the same named satellite-minus-reference rows.
 * Numeric row indices alone are insufficient because two independently
 * compiled bases can assign the same index to different satellites. */
inline bool zhangProductNamedOrderingMatches(
	const ZhangProductRelationBasis& first,
	const ZhangProductRelationBasis& second)
{
	if (first.mappableNamedIndices.size() !=
		second.mappableNamedIndices.size())
	{
		return false;
	}
	for (std::size_t local = 0;
		 local < first.mappableNamedIndices.size(); local++)
	{
		const int firstIndex = first.mappableNamedIndices[local];
		const int secondIndex = second.mappableNamedIndices[local];
		if (firstIndex < 0 || secondIndex < 0 ||
			firstIndex >= static_cast<int>(first.namedRelations.size()) ||
			secondIndex >= static_cast<int>(second.namedRelations.size()))
		{
			return false;
		}
		const auto& firstRow = first.namedRelations[firstIndex];
		const auto& secondRow = second.namedRelations[secondIndex];
		if (firstRow.satellite != secondRow.satellite ||
			firstRow.referenceSatellite != secondRow.referenceSatellite)
		{
			return false;
		}
	}
	return true;
}

/** Exact affine pull-back of product WL and first-signal integer rows.
 *
 * For Fw(z1-z2)=nw and F1 z1=n1 this produces
 * [Fw(G1-G2); F1 G1] a =
 * [nw-Fw(c1-c2); n1-F1 c1].  Numeric basis matrices are accepted only when
 * every coefficient is exactly integral within the representation tolerance. */
inline bool zhangPullBackProductIntegerConstraints(
	const ZhangProductRelationBasis& firstBasis,
	const ZhangProductRelationBasis& secondBasis,
	const ZhangExactMatrix& wideLaneRows,
	const ZhangExactVector& wideLaneIntegers,
	const ZhangExactMatrix& firstSignalRows,
	const ZhangExactVector& firstSignalIntegers,
	ZhangExactMatrix& networkRows,
	ZhangExactVector& networkIntegers,
	std::string& failureReason)
{
	networkRows.clear();
	networkIntegers.clear();
	const int productRank = firstBasis.mappableTargetRank;
	const int networkDimension = firstBasis.transform.cols();
	if (productRank <= 0 || networkDimension <= 0 ||
		secondBasis.mappableTargetRank != productRank ||
		secondBasis.transform.cols() != networkDimension ||
		firstBasis.transform.rows() != productRank ||
		secondBasis.transform.rows() != productRank ||
		firstBasis.affineOffsets.size() != static_cast<std::size_t>(productRank) ||
		secondBasis.affineOffsets.size() != static_cast<std::size_t>(productRank) ||
		wideLaneRows.size() != wideLaneIntegers.size() ||
		firstSignalRows.size() != firstSignalIntegers.size())
	{
		failureReason = "PRODUCT_CONSTRAINT_DIMENSION_MISMATCH";
		return false;
	}
	auto exactBasisCoefficient = [](const MatrixXd& transform,
		int row, int column, ZhangExactInteger& coefficient)
	{
		const long long rounded = std::llround(transform(row, column));
		if (std::abs(transform(row, column) - rounded) > 1e-8) return false;
		coefficient = rounded;
		return true;
	};
	auto append = [&](const ZhangExactVector& productRow,
		const ZhangExactInteger& integer, bool wideLane)
	{
		if (productRow.size() != static_cast<std::size_t>(productRank))
			return false;
		ZhangExactVector networkRow(networkDimension);
		ZhangExactInteger rhs = integer;
		for (int product = 0; product < productRank; product++)
		{
			const ZhangExactInteger multiplier = productRow[product];
			if (multiplier == 0) continue;
			for (int column = 0; column < networkDimension; column++)
			{
				ZhangExactInteger first = 0;
				ZhangExactInteger second = 0;
				if (!exactBasisCoefficient(firstBasis.transform,
					product, column, first) ||
					!exactBasisCoefficient(secondBasis.transform,
					product, column, second)) return false;
				networkRow[column] += multiplier *
					(wideLane ? first - second : first);
			}
			const ZhangExactInteger offset = wideLane
				? firstBasis.affineOffsets.at(product) -
					secondBasis.affineOffsets.at(product)
				: firstBasis.affineOffsets.at(product);
			rhs -= multiplier * offset;
		}
		networkRows.push_back(std::move(networkRow));
		networkIntegers.push_back(std::move(rhs));
		return true;
	};
	for (std::size_t row = 0; row < wideLaneRows.size(); row++)
	{
		if (!append(wideLaneRows[row], wideLaneIntegers[row], true))
		{
			failureReason = "WL_PRODUCT_TO_NETWORK_MAPPING_FAILED";
			return false;
		}
	}
	for (std::size_t row = 0; row < firstSignalRows.size(); row++)
	{
		if (!append(firstSignalRows[row], firstSignalIntegers[row], false))
		{
			failureReason = "L1_PRODUCT_TO_NETWORK_MAPPING_FAILED";
			return false;
		}
	}
	failureReason = "NONE";
	return true;
}

/** Split exact rows expressed in [L1,L2] product coordinates into the two
 * single-signal lattices used by the dual-frequency completion search.
 *
 * A row [q,-q] is a WL integer, while [q,0] is an L1 integer.  Every other
 * mixed row remains useful as a conditioner but is deliberately excluded
 * from graph support.  Affine HNF is applied together with the RHS so a
 * contradictory historical value fails closed. */
struct ZhangDualSignalIntegerSupport
{
	ZhangExactMatrix wideLaneRows;
	ZhangExactVector wideLaneValues;
	ZhangExactMatrix firstSignalRows;
	ZhangExactVector firstSignalValues;
	int inputRows = 0;
	int mixedConditioningOnlyRows = 0;
	int invalidRows = 0;
	bool valid = false;
	std::string failureReason = "NOT_EVALUATED";
};

/** Pull exact rows from the unimodular [WL,L1] chart back to [L1,L2].
 *
 * With w=a1-a2 and y=[w,a1]^T, y=T[a1,a2]^T for
 * T=[[I,-I],[I,0]].  A fixed row r_y therefore becomes r_x=r_y T, i.e.
 * [r_w+r_1,-r_w].  Keeping this operation exact is essential: rounding a
 * dense LAMBDA row here can manufacture a satellite-product integer that was
 * never fixed by the parent joint search. */
inline ZhangExactMatrix zhangWideLaneFirstRowsToSignalRows(
	const ZhangExactMatrix& wideLaneFirstRows,
	int signalDimension,
	bool* dimensionsValid = nullptr)
{
	if (dimensionsValid) *dimensionsValid = false;
	if (signalDimension <= 0) return {};
	ZhangExactMatrix signalRows;
	signalRows.reserve(wideLaneFirstRows.size());
	for (const auto& row : wideLaneFirstRows)
	{
		if (row.size() != static_cast<std::size_t>(2 * signalDimension))
			return {};
		ZhangExactVector signalRow(2 * signalDimension);
		for (int column = 0; column < signalDimension; column++)
		{
			signalRow[column] = row[column] + row[signalDimension + column];
			signalRow[signalDimension + column] = -row[column];
		}
		signalRows.push_back(std::move(signalRow));
	}
	if (dimensionsValid) *dimensionsValid = true;
	return signalRows;
}

/** Exact affine coordinates for the stochastic quotient of an already-fixed
 * integer lattice D z = d.
 *
 * The rows of kernelBasis form the saturated integer kernel of D.  Hence every
 * admissible integer is represented uniquely as
 *
 *     z = particularSolution + kernelBasis^T q,  q in Z^quotientRank.
 *
 * quotientProjector is an exact integer left inverse of kernelBasis^T, so it
 * maps the (possibly singular) ambient posterior to the only coordinates that
 * still require ILS.  This is the affine-kernel form of the same primitive
 * HNF/SNF quotient; no floating null-space vector or rounded inverse enters
 * the construction. */
struct ZhangExactAffineIntegerQuotient
{
	bool valid = false;
	int ambientDimension = 0;
	int deterministicRank = 0;
	int quotientRank = 0;
	ZhangExactMatrix deterministicBasis;
	ZhangExactVector deterministicValues;
	ZhangExactVector particularSolution;
	ZhangExactMatrix kernelBasis;
	ZhangExactMatrix quotientProjector;
	std::vector<ZhangExactInteger> smithInvariants;
	std::string failureReason = "NOT_EVALUATED";
};

inline ZhangExactAffineIntegerQuotient zhangExactAffineIntegerQuotient(
	const ZhangExactMatrix& deterministicRows,
	const ZhangExactVector& deterministicValues,
	int ambientDimension)
{
	ZhangExactAffineIntegerQuotient result;
	result.ambientDimension = ambientDimension;
	if (ambientDimension <= 0 ||
		deterministicRows.size() != deterministicValues.size())
	{
		result.failureReason = "AFFINE_QUOTIENT_DIMENSION_MISMATCH";
		return result;
	}
	for (const auto& row : deterministicRows)
	{
		if (row.size() != static_cast<std::size_t>(ambientDimension))
		{
			result.failureReason = "AFFINE_QUOTIENT_INVALID_ROW";
			return result;
		}
	}
	const auto hnf = zhangExactRowHermiteNormalForm(
		deterministicRows, deterministicValues);
	if (!hnf.consistent)
	{
		result.failureReason = "AFFINE_QUOTIENT_INTEGER_CONFLICT";
		return result;
	}
	result.deterministicBasis = hnf.basis;
	result.deterministicValues = hnf.values;
	result.deterministicRank = static_cast<int>(hnf.basis.size());
	if (result.deterministicRank > ambientDimension)
	{
		result.failureReason = "AFFINE_QUOTIENT_OVER_RANK";
		return result;
	}

	const auto primitive = zhangIntegerRowLatticeContains(
		result.deterministicBasis, ZhangExactVector(ambientDimension));
	result.smithInvariants = primitive.smithInvariants;
	// A nonprimitive equation is not infeasible: 2*a=6 is a=3. The exact
	// column-lattice membership below enforces divisibility (2*a=5 fails).
	if (primitive.rank != result.deterministicRank)
	{
		result.failureReason = "AFFINE_QUOTIENT_RANK_MISMATCH";
		return result;
	}
	if (result.deterministicRank == 0)
	{
		result.particularSolution = ZhangExactVector(ambientDimension);
		result.kernelBasis = zhangExactIdentityMatrix(ambientDimension);
		result.quotientProjector = zhangExactIdentityMatrix(ambientDimension);
		result.quotientRank = ambientDimension;
		result.valid = true;
		result.failureReason = "NONE";
		return result;
	}

	// Solve D z0=d exactly by treating the columns of D as generators of Z^r.
	ZhangExactMatrix columnGenerators(ambientDimension,
		ZhangExactVector(result.deterministicRank));
	for (int column = 0; column < ambientDimension; column++)
	for (int row = 0; row < result.deterministicRank; row++)
		columnGenerators[column][row] = result.deterministicBasis[row][column];
	const auto particular = zhangIntegerRowLatticeContains(
		columnGenerators, result.deterministicValues);
	if (!particular.contained ||
		particular.combination.size() !=
			static_cast<std::size_t>(ambientDimension))
	{
		result.failureReason = "AFFINE_QUOTIENT_NO_INTEGER_ORIGIN";
		return result;
	}
	result.particularSolution = particular.combination;

	result.kernelBasis = zhangExactIntegerKernel(
		result.deterministicBasis, ambientDimension);
	result.quotientRank = static_cast<int>(result.kernelBasis.size());
	if (result.deterministicRank + result.quotientRank != ambientDimension)
	{
		result.failureReason = "AFFINE_QUOTIENT_KERNEL_RANK_MISMATCH";
		return result;
	}

	// Find an exact integer left inverse L K^T=I.  The saturated kernel is
	// primitive, so every quotient unit vector must lie in the column lattice
	// of K.
	ZhangExactMatrix kernelColumns(ambientDimension,
		ZhangExactVector(result.quotientRank));
	for (int column = 0; column < ambientDimension; column++)
	for (int row = 0; row < result.quotientRank; row++)
		kernelColumns[column][row] = result.kernelBasis[row][column];
    const auto inverseRows=zhangIntegerRowLatticeContainsBatch(
        kernelColumns,zhangExactIdentityMatrix(result.quotientRank));
	result.quotientProjector.reserve(result.quotientRank);
	for (int quotient = 0; quotient < result.quotientRank; quotient++)
	{
		const auto& inverseRow=inverseRows[quotient];
		if (!inverseRow.contained ||
			inverseRow.combination.size() !=
				static_cast<std::size_t>(ambientDimension))
		{
			result.failureReason = "AFFINE_QUOTIENT_LEFT_INVERSE_FAILED";
			return result;
		}
		result.quotientProjector.push_back(inverseRow.combination);
	}
	const auto identityAudit = zhangExactMultiply(
		result.quotientProjector,
		[&]()
		{
			ZhangExactMatrix transpose(ambientDimension,
				ZhangExactVector(result.quotientRank));
			for (int row = 0; row < result.quotientRank; row++)
			for (int column = 0; column < ambientDimension; column++)
				transpose[column][row] = result.kernelBasis[row][column];
			return transpose;
		}());
	if (identityAudit != zhangExactIdentityMatrix(result.quotientRank))
	{
		result.failureReason = "AFFINE_QUOTIENT_LEFT_INVERSE_AUDIT_FAILED";
		return result;
	}
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

/** One immutable branch input captured at the same posterior moment as P.
 * The snapshot contains only exact, currently reauthorised persistent rows.
 * Posterior nullity may be larger because of external hard priors; that fact
 * is recorded separately and is not mislabelled as missing persistent state. */
struct ZhangPersistentCertificateSnapshot
{
	bool valid = false;
	ZhangExactMatrix lattice;
	ZhangExactVector affineRhs;
	std::vector<std::string> certificateFamilyIds;
	std::string latticeHash;
	std::string affineHash;
	std::string provenanceHash;
	std::string posteriorMomentId;
	std::string snapshotId;
	int deterministicRank = 0;
	int expectedPosteriorNullity = 0;
	int externalPosteriorNullity = 0;
	double currentCertificateRisk = 0;
	double lifetimeBroadcastRisk = 0;
	std::string failureReason = "NOT_EVALUATED";
};

/** Coordinates of the exact image A Z^n.  Never replace a nonprimitive
 * image by Z^r: B contains its congruences and A=B U with U surjective.
 * Solving integers in U coordinates therefore preserves, e.g., 2Z. */
struct ZhangPrimitiveImageCoordinates
{
	bool valid = false;
	ZhangExactMatrix primitiveRows;
	ZhangExactMatrix imageGenerators;
};

inline ZhangPrimitiveImageCoordinates zhangPrimitiveImageCoordinates(const ZhangExactMatrix& rows)
{
	ZhangPrimitiveImageCoordinates result;
	if (rows.empty() || rows.front().empty()) return result;
	const int m = rows.size(), n = rows.front().size();
	ZhangExactMatrix columns(n, ZhangExactVector(m));
	for (int i=0; i<m; ++i)
	{
		if (rows[i].size() != n) return result;
		for (int j=0; j<n; ++j) columns[j][i]=rows[i][j];
	}
	const auto image = zhangExactRowHermiteNormalForm(columns);
	if (!image.consistent) return result;
	const int r = image.basis.size();
	result.primitiveRows.assign(r, ZhangExactVector(n));
	result.imageGenerators.assign(m, ZhangExactVector(r));
	for (int i=0;i<m;++i) for(int k=0;k<r;++k)
		result.imageGenerators[i][k]=image.basis[k][i];
	for (int j=0;j<n;++j)
	{
		const auto membership = zhangIntegerRowLatticeContains(image.basis, columns[j]);
		if (!membership.contained || membership.combination.size() != r) return result;
		for (int k=0;k<r;++k) result.primitiveRows[k][j]=membership.combination[k];
	}
	if (r && (!zhangExactPrimitiveRowLattice(result.primitiveRows, n) ||
		zhangExactMultiply(result.imageGenerators, result.primitiveRows) != rows)) return result;
	result.valid=true;
	return result;
}

struct ZhangCertificateFamilyRisk
{
	std::string familyId;
	double currentRisk = 1;
	double lifetimeRisk = 0;
	bool currentReauthorised = false;
	ZhangDecisionProofs decisionProofs;
};

struct ZhangCertificateRiskAudit
{
	bool valid = false;
	int currentFamilies = 0;
	int lifetimeFamilies = 0;
	double currentCertificateRisk = 0;
	double lifetimeBroadcastRisk = 0;
	std::string failureReason = "NOT_EVALUATED";
};

/** Deduplicate statistical views by immutable family id.  Current risk and
 * lifetime service risk are deliberately accumulated into separate ledgers. */
inline ZhangCertificateRiskAudit zhangAuditCertificateFamilyRisk(
	const std::vector<ZhangCertificateFamilyRisk>& receipts)
{
	ZhangCertificateRiskAudit result;
	std::map<std::string, double> current;
	std::map<std::string, double> lifetime;
	ZhangDecisionProofs allProofs;
	for (const auto& receipt : receipts)
		allProofs = zhangMergeDecisionProofs(allProofs, receipt.decisionProofs);
	const auto identityAudit = zhangDecisionRiskClosure(allProofs);
	if (!identityAudit.valid)
	{
		result.failureReason = identityAudit.reason;
		return result;
	}
	for (const auto& receipt : receipts)
	{
		if (!receipt.decisionProofs.empty())
		{
			const auto closure = zhangDecisionRiskClosure(receipt.decisionProofs);
			if (!closure.valid)
			{
				result.failureReason = closure.reason;
				return result;
			}
			for (const auto& [id, atom] : closure.atoms)
			{
				if (receipt.currentReauthorised)
					current[id] = std::max(current[id], atom->conditionalFailureBound);
				lifetime[id] = std::max(lifetime[id], atom->conditionalFailureBound);
			}
			continue;
		}
		if (receipt.familyId.empty() || !std::isfinite(receipt.currentRisk) ||
			receipt.currentRisk < 0 || receipt.currentRisk > 1 ||
			!std::isfinite(receipt.lifetimeRisk) || receipt.lifetimeRisk < 0 ||
			receipt.lifetimeRisk > 1)
		{
			result.failureReason = "CERTIFICATE_FAMILY_RISK_INVALID";
			return result;
		}
		if (receipt.currentReauthorised)
			current[receipt.familyId] = std::max(
				current[receipt.familyId], receipt.currentRisk);
		lifetime[receipt.familyId] = std::max(
			lifetime[receipt.familyId], receipt.lifetimeRisk);
	}
	for (const auto& [_, risk] : current)
		result.currentCertificateRisk = std::min(
			1.0, result.currentCertificateRisk + risk);
	for (const auto& [_, risk] : lifetime)
		result.lifetimeBroadcastRisk = std::min(
			1.0, result.lifetimeBroadcastRisk + risk);
	result.currentFamilies = static_cast<int>(current.size());
	result.lifetimeFamilies = static_cast<int>(lifetime.size());
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangPersistentCertificateSnapshot
zhangBuildPersistentCertificateSnapshot(
	const ZhangExactMatrix& rows,
	const ZhangExactVector& values,
	int ambientDimension,
	int posteriorNullity,
	std::vector<std::string> certificateFamilyIds,
	double currentCertificateRisk,
	double lifetimeBroadcastRisk,
	std::string posteriorMomentId)
{
	ZhangPersistentCertificateSnapshot result;
	result.expectedPosteriorNullity = posteriorNullity;
	result.currentCertificateRisk = currentCertificateRisk;
	result.lifetimeBroadcastRisk = lifetimeBroadcastRisk;
	result.posteriorMomentId = std::move(posteriorMomentId);
	if (ambientDimension <= 0 || posteriorNullity < 0 ||
		posteriorNullity > ambientDimension || rows.size() != values.size() ||
		result.posteriorMomentId.empty() ||
		!std::isfinite(currentCertificateRisk) || currentCertificateRisk < 0 ||
		currentCertificateRisk > 1 || !std::isfinite(lifetimeBroadcastRisk) ||
		lifetimeBroadcastRisk < 0 || lifetimeBroadcastRisk > 1)
	{
		result.failureReason = "PERSISTENT_SNAPSHOT_INPUT_INVALID";
		return result;
	}
	for (const auto& row : rows)
		if (row.size() != static_cast<std::size_t>(ambientDimension))
		{
			result.failureReason = "PERSISTENT_SNAPSHOT_DIMENSION_MISMATCH";
			return result;
		}
	const auto hnf = zhangExactRowHermiteNormalForm(rows, values);
	if (!hnf.consistent)
	{
		result.failureReason = "PERSISTENT_SNAPSHOT_AFFINE_CONFLICT";
		return result;
	}
	result.lattice = hnf.basis;
	result.affineRhs = hnf.values;
	result.deterministicRank = static_cast<int>(result.lattice.size());
	if (result.deterministicRank > posteriorNullity)
	{
		result.failureReason = "INTERNAL_PERSISTENT_SNAPSHOT_DIVERGENCE";
		return result;
	}
	result.externalPosteriorNullity = posteriorNullity - result.deterministicRank;
	std::sort(certificateFamilyIds.begin(), certificateFamilyIds.end());
	certificateFamilyIds.erase(std::unique(certificateFamilyIds.begin(),
		certificateFamilyIds.end()), certificateFamilyIds.end());
	result.certificateFamilyIds = std::move(certificateFamilyIds);
	result.latticeHash = zhangExactMatrixFingerprint(result.lattice);
	result.affineHash = zhangExactMatrixFingerprint({result.affineRhs});
	std::uint64_t provenance = 1469598103934665603ULL;
	for (const auto& family : result.certificateFamilyIds)
		provenance = zhangAuditFnv1a(provenance, family + "|");
	std::ostringstream provenanceText;
	provenanceText << std::hex << provenance;
	result.provenanceHash = provenanceText.str();
	result.snapshotId = result.posteriorMomentId + ":" + result.latticeHash +
		":" + result.affineHash + ":" + result.provenanceHash;
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

/** Union-bound accounting for a current reauthorization epoch.  Reusing an
 * already admitted persistent family retains its historical risk.  An epoch
 * that accepts no new integer claim contributes exactly zero incremental
 * risk; it is neither a failure with probability one nor fresh evidence that
 * resets the old risk to zero. */
inline double zhangPersistentIncrementalFailureProbability(
	double persistentFailureProbability,
	bool newIntegerClaimAccepted,
	double newFailureProbability)
{
	if (!std::isfinite(persistentFailureProbability) ||
		persistentFailureProbability < 0 || persistentFailureProbability > 1)
		return 1;
	if (!newIntegerClaimAccepted) return persistentFailureProbability;
	if (!std::isfinite(newFailureProbability) ||
		newFailureProbability < 0 || newFailureProbability > 1)
		return 1;
	return std::min(1.0,
		persistentFailureProbability + newFailureProbability);
}

enum class ZhangIncrementalRiskBudgetStatus
{
	INVALID,
	PARENT_INFEASIBLE,
	BUDGET_EXHAUSTED,
	SEARCH_WITH_POSITIVE_BUDGET
};

struct ZhangIncrementalRiskBudget
{
	ZhangIncrementalRiskBudgetStatus status =
		ZhangIncrementalRiskBudgetStatus::INVALID;
	double parentRisk = 1;
	double totalBudget = 0;
	double remainingBudget = 0;
	bool parentFeasible = false;
	bool searchAuthorized = false;
};

struct ZhangPosteriorReceiptAtomicDecision
{
	bool useConditionedPosterior = false;
	bool rollbackRequired = false;
	std::string status = "PRECONDITION_NOT_APPLIED";
};

/** A private posterior may be used only together with every exact constraint
 * and risk receipt that conditioned it. */
inline ZhangPosteriorReceiptAtomicDecision
zhangSelectPosteriorReceiptTransaction(
	bool posteriorConditioned,
	bool receiptAuthorized)
{
	ZhangPosteriorReceiptAtomicDecision result;
	if (!posteriorConditioned) return result;
	if (!receiptAuthorized)
	{
		result.rollbackRequired = true;
		result.status = "ROLLBACK_CONDITIONED_POSTERIOR_WITHOUT_RECEIPT";
		return result;
	}
	result.useConditionedPosterior = true;
	result.status = "CONDITIONED_POSTERIOR_AND_RECEIPT_COMMITTED";
	return result;
}

/** Allocate risk to a genuinely new integer decision.  A negative remainder
 * is never converted into a tiny positive search budget: an infeasible parent
 * or an exhausted budget is a terminal state for the incremental search. */
inline ZhangIncrementalRiskBudget zhangAllocateIncrementalRiskBudget(
	double parentRisk,
	double totalBudget)
{
	ZhangIncrementalRiskBudget result;
	result.parentRisk = parentRisk;
	result.totalBudget = totalBudget;
	if (!std::isfinite(parentRisk) || parentRisk < 0 || parentRisk > 1 ||
		!std::isfinite(totalBudget) || totalBudget < 0 || totalBudget > 1)
		return result;
	if (parentRisk > totalBudget)
	{
		result.status = ZhangIncrementalRiskBudgetStatus::PARENT_INFEASIBLE;
		return result;
	}
	result.parentFeasible = true;
	result.remainingBudget = totalBudget - parentRisk;
	if (result.remainingBudget <= 0)
	{
		result.status = ZhangIncrementalRiskBudgetStatus::BUDGET_EXHAUSTED;
		return result;
	}
	result.searchAuthorized = true;
	result.status = ZhangIncrementalRiskBudgetStatus::SEARCH_WITH_POSITIVE_BUDGET;
	return result;
}

struct ZhangOptionalLedgerRiskDecision
{
	bool baselineFeasible = false;
	bool optionalRequested = false;
	bool optionalAuthorized = false;
	bool baselinePreserved = false;
	double baselineFailureProbability = 1;
	double optionalFailureProbability = 1;
	double proposedUnionBound = 1;
	double selectedFailureProbability = 1;
	double remainingRiskBudget = 0;
	std::string status = "INVALID_RISK_INPUT";
};

/** Keep an already feasible authoritative baseline immutable while evaluating
 * an optional historical Ledger family.  The optional family is charged by a
 * conservative union bound and is admitted only inside the remaining risk
 * budget.  A rejected or malformed optional family therefore cannot turn a
 * feasible baseline into an infeasible product certificate. */
inline ZhangOptionalLedgerRiskDecision
zhangSelectOptionalLedgerWithoutPoisoningBaseline(
	bool baselineAvailable,
	double baselineFailureProbability,
	bool optionalRequested,
	double optionalFailureProbability,
	double totalRiskBudget)
{
	ZhangOptionalLedgerRiskDecision result;
	result.optionalRequested = optionalRequested;
	const bool budgetValid = std::isfinite(totalRiskBudget) &&
		totalRiskBudget >= 0 && totalRiskBudget <= 1;
	const bool baselineRiskValid = std::isfinite(baselineFailureProbability) &&
		baselineFailureProbability >= 0 && baselineFailureProbability <= 1;
	result.baselineFeasible = baselineAvailable && budgetValid &&
		baselineRiskValid && baselineFailureProbability <= totalRiskBudget;
	// baselineAvailable describes a publishable product, NOT whether the
	// posterior already depends on integer decisions. Charge parents either way.
	if (budgetValid && baselineRiskValid && baselineFailureProbability <= totalRiskBudget)
	{
		result.baselinePreserved = result.baselineFeasible;
		result.baselineFailureProbability = baselineFailureProbability;
		result.selectedFailureProbability = baselineFailureProbability;
		result.remainingRiskBudget = totalRiskBudget - baselineFailureProbability;
		if (!optionalRequested)
		{
			result.optionalFailureProbability = 0;
			result.proposedUnionBound = baselineFailureProbability;
			result.status = baselineAvailable ? "BASELINE_ONLY_NO_OPTIONAL_LEDGER"
				: "POSTERIOR_DEPENDENCIES_ONLY_NO_OPTIONAL_LEDGER";
			return result;
		}
		const bool optionalRiskValid = std::isfinite(optionalFailureProbability) &&
			optionalFailureProbability >= 0 && optionalFailureProbability <= 1;
		if (!optionalRiskValid)
		{
			result.status = "BASELINE_PRESERVED_INVALID_OPTIONAL_RISK";
			return result;
		}
		result.optionalFailureProbability = optionalFailureProbability;
		result.proposedUnionBound = std::min(1.0,
			baselineFailureProbability + optionalFailureProbability);
		if (result.proposedUnionBound <= totalRiskBudget)
		{
			result.optionalAuthorized = true;
			result.selectedFailureProbability = result.proposedUnionBound;
			result.status = baselineAvailable ? "BASELINE_PLUS_OPTIONAL_WITHIN_RISK_BUDGET"
				: "POSTERIOR_PLUS_OPTIONAL_WITHIN_RISK_BUDGET";
		}
		else
		{
			result.status = baselineAvailable ? "BASELINE_PRESERVED_OPTIONAL_EXCEEDS_RISK_BUDGET"
				: "POSTERIOR_PRESERVED_OPTIONAL_EXCEEDS_RISK_BUDGET";
		}
		return result;
	}

	result.baselineFailureProbability = baselineRiskValid
		? baselineFailureProbability : 1;
	result.status = budgetValid && baselineRiskValid
		? "BASELINE_EXCEEDS_RISK_BUDGET"
		: "INVALID_BASELINE_RISK_INPUT";
	return result;
}

inline ZhangDualSignalIntegerSupport zhangSplitJointProductIntegerSupport(
	const ZhangExactMatrix& jointRows,
	const ZhangExactVector& jointValues,
	int productDimension)
{
	ZhangDualSignalIntegerSupport result;
	result.inputRows = jointRows.size();
	if (productDimension <= 0 || jointRows.size() != jointValues.size())
	{
		result.failureReason = "JOINT_PRODUCT_SUPPORT_DIMENSION_MISMATCH";
		return result;
	}
	for (const auto& row : jointRows)
	{
		if (row.size() != static_cast<std::size_t>(2 * productDimension))
		{
			result.invalidRows++;
			result.failureReason = "JOINT_PRODUCT_SUPPORT_INVALID_ROW";
			return result;
		}
	}
	const auto jointHnf = zhangExactRowHermiteNormalForm(jointRows, jointValues);
	if (!jointHnf.consistent)
	{
		result.failureReason = "JOINT_PRODUCT_SUPPORT_AFFINE_CONFLICT";
		return result;
	}

	// Intersect the *whole lattice*, rather than classifying only its current
	// basis rows.  HNF can mix two pure signal constraints, and conversely two
	// mixed historical rows can have a pure WL or L1 integer combination.
	// For WL use the unimodular chart [a,b] -> [a,a+b]; removing its second
	// block gives exactly b=-a.  For L1, removing the original b block gives
	// exactly b=0.  RHS values follow the same exact integer combinations.
	ZhangExactMatrix wideLaneChart = jointHnf.basis;
	for (auto& row : wideLaneChart)
	for (int column = 0; column < productDimension; column++)
		row[productDimension + column] += row[column];
	std::vector<bool> keepFirst(productDimension * 2, false);
	std::fill(keepFirst.begin(), keepFirst.begin() + productDimension, true);
	const auto wideLaneIntersection = zhangExactSurvivingLattice(
		wideLaneChart, jointHnf.values, keepFirst);
	const auto firstSignalIntersection = zhangExactSurvivingLattice(
		jointHnf.basis, jointHnf.values, keepFirst);
	if (!wideLaneIntersection.consistent)
	{
		result.failureReason = "JOINT_PRODUCT_SUPPORT_WL_INTERSECTION_FAILED";
		return result;
	}
	if (!firstSignalIntersection.consistent)
	{
		result.failureReason = "JOINT_PRODUCT_SUPPORT_L1_INTERSECTION_FAILED";
		return result;
	}
	result.wideLaneRows = wideLaneIntersection.basis;
	result.wideLaneValues = wideLaneIntersection.values;
	result.firstSignalRows = firstSignalIntersection.basis;
	result.firstSignalValues = firstSignalIntersection.values;
	result.mixedConditioningOnlyRows = std::max(0,
		static_cast<int>(jointHnf.basis.size()) -
		static_cast<int>(result.wideLaneRows.size()) -
		static_cast<int>(result.firstSignalRows.size()));
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

/** Recover every individually proven named target from an exact fixed row
 * lattice.  Each returned unit row has passed exact integer HNF/Smith lattice
 * containment; unrecoverable named rows are omitted rather than fabricated. */
inline std::map<std::size_t, ZhangExactInteger>
zhangRecoverCertifiedNamedProductSubset(
	const ZhangExactMatrix& fixedRows,
	const ZhangExactVector& fixedValues,
	int namedCount)
{
	return ProductConstraintPromotion::recoverNamedTargets(
		fixedRows, fixedValues, namedCount);
}

/** Full-subset wrapper used by outer named PAR.  Inner LAMBDA/PAR may expose a
 * smaller exact named seed through zhangRecoverCertifiedNamedProductSubset(),
 * but a branch is a complete certificate only when every row is recovered. */
inline std::map<std::size_t, ZhangExactInteger>
zhangRecoverCompleteNamedProductSubset(
	const ZhangExactMatrix& fixedRows,
	const ZhangExactVector& fixedValues,
	int namedCount)
{
	auto recovered = zhangRecoverCertifiedNamedProductSubset(
		fixedRows, fixedValues, namedCount);
	if (recovered.size() != static_cast<std::size_t>(namedCount))
	{
		return {};
	}
	for (int index = 0; index < namedCount; index++)
	{
		if (!recovered.contains(index)) return {};
	}
	return recovered;
}

/** Exact promotion of named coordinates from an already accepted parent
 * integer lattice.
 *
 * If the parent mixed LAMBDA/PAR lattice has passed its bootstrap-success and
 * absolute NIS gates, every unit row contained in that lattice is a
 * deterministic integer function of the accepted parent solution.  Such a
 * row must inherit the parent certificate; re-running LAMBDA on its marginal
 * covariance is mathematically different because it discards the joint
 * constraints that made the row identifiable.  This helper performs only the
 * exact HNF membership/value proof.  The caller remains responsible for
 * proving that the parent search was statistically accepted. */
struct ZhangInheritedNamedCertificate
{
	bool exact = false;
	int parentFixedRank = 0;
	std::map<std::size_t, ZhangExactInteger> values;
};

/** Recover every directly named satellite-pair edge contained in an accepted
 * mixed fixed lattice.  Unlike named-star promotion this tests both unit rows
 * e_s and pair rows e_s-e_t.  Higher-order fixed integers that contain no pair
 * edge are deliberately retained only as conditioning evidence by callers. */
inline std::vector<ZhangCertifiedPairRelation>
zhangRecoverCertifiedPairRelations(
	const ZhangExactMatrix& fixedRows,
	const ZhangExactVector& fixedValues,
	int namedCount,
	bool parentStatisticallyAccepted,
	int* batchTargets = nullptr,
	int* scalarDecompositionsAvoided = nullptr)
{
	std::vector<ZhangCertifiedPairRelation> result;
	if (batchTargets) *batchTargets = 0;
	if (scalarDecompositionsAvoided) *scalarDecompositionsAvoided = 0;
	if (!parentStatisticallyAccepted || namedCount <= 0 ||
		fixedRows.empty() || fixedRows.size() != fixedValues.size())
	{
		return result;
	}
	std::vector<std::pair<int, int>> pairNodes;
	ZhangExactMatrix pairRows;
	for (int first = 0; first <= namedCount; first++)
	for (int second = first + 1; second <= namedCount; second++)
	{
		ZhangExactVector pairRow(namedCount);
		if (first < namedCount) pairRow[first] += 1;
		if (second < namedCount) pairRow[second] -= 1;
		pairNodes.push_back({first, second});
		pairRows.push_back(std::move(pairRow));
	}
	const auto memberships = zhangIntegerRowLatticeContainsBatch(
		fixedRows, pairRows);
	if (batchTargets) *batchTargets = pairRows.size();
	if (scalarDecompositionsAvoided)
		*scalarDecompositionsAvoided = std::max(0,
			static_cast<int>(pairRows.size()) - 1);
	if (memberships.size() != pairRows.size()) return {};
	for (int pair = 0; pair < static_cast<int>(pairRows.size()); pair++)
	{
		const auto& membership = memberships[pair];
		if (!membership.contained ||
			membership.combination.size() != fixedValues.size())
		{
			continue;
		}
		ZhangExactInteger value = 0;
		for (std::size_t row = 0; row < fixedValues.size(); row++)
		{
			value += membership.combination[row] * fixedValues[row];
		}
		result.push_back({pairNodes[pair].first, pairNodes[pair].second, value,
			membership.combination});
	}
	return result;
}

struct ZhangPairReliabilityEdge
{
	int firstNode = -1;
	int secondNode = -1;
	double perr = 1;
	double variance = std::numeric_limits<double>::infinity();
};

struct ZhangExactDualGraphScore
{
	bool valid = false;
	int graphRank = 0, satellites = 0, components = 0, largest = 0;
};

struct ZhangExactDualGraphProfile
{
	ZhangExactDualGraphScore score;
	// Every canonical satellite node is represented.  A node absent from the
	// current strict dual graph is a singleton component and therefore remains
	// an explicit attachment target rather than disappearing from the search.
	std::vector<int> componentByNode;
	std::vector<int> componentSizeByNode;
	std::set<std::pair<int, int>> dualFrequencyEdges;
};

// Recover exact pair consequences common to BOTH frequency lattices, including
// singleton nodes needed by bridge/attachment-directed PAR. Dense integer rank
// by itself is not a usable satellite component.
inline ZhangExactDualGraphProfile zhangExactDualGraphProfile(
	const ZhangExactMatrix& signalRows, const ZhangExactVector& values, int namedCount)
{
	ZhangExactDualGraphProfile result;
	if (namedCount < 0) return result;
	const auto split = zhangSplitJointProductIntegerSupport(signalRows, values, namedCount);
	if (!split.valid) return result;
	const auto wl = zhangRecoverCertifiedPairRelations(split.wideLaneRows,
		split.wideLaneValues, namedCount, true);
	const auto l1 = zhangRecoverCertifiedPairRelations(split.firstSignalRows,
		split.firstSignalValues, namedCount, true);
	std::set<std::pair<int, int>> firstPairs;
	for (const auto& pair : l1) firstPairs.insert({pair.firstNode, pair.secondNode});
	std::vector<int> parent(namedCount + 1);
	std::iota(parent.begin(), parent.end(), 0);
	auto root = [&](int i) { while (i != parent[i]) i = parent[i]; return i; };
	std::set<int> nodes;
	for (const auto& pair : wl)
	{
		if (!firstPairs.count({pair.firstNode, pair.secondNode})) continue;
		if (pair.firstNode < 0 || pair.firstNode > namedCount ||
			pair.secondNode < 0 || pair.secondNode > namedCount) continue;
		result.dualFrequencyEdges.insert({pair.firstNode, pair.secondNode});
		nodes.insert(pair.firstNode); nodes.insert(pair.secondNode);
		int a = root(pair.firstNode), b = root(pair.secondNode);
		if (a != b) { parent[b] = a; ++result.score.graphRank; }
	}
	std::map<int, int> sizes;
	for (int node = 0; node <= namedCount; node++) ++sizes[root(node)];
	result.componentByNode.resize(namedCount + 1);
	result.componentSizeByNode.resize(namedCount + 1);
	for (int node = 0; node <= namedCount; node++)
	{
		result.componentByNode[node] = root(node);
		result.componentSizeByNode[node] = sizes[root(node)];
	}
	std::set<int> activeComponents;
	for (int node : nodes) activeComponents.insert(root(node));
	result.score.satellites = nodes.size();
	result.score.components = activeComponents.size();
	for (int component : activeComponents)
		result.score.largest = std::max(result.score.largest, sizes[component]);
	result.score.valid = true;
	return result;
}

/** Finalise writer authority for historical product-gauge evidence.
 *
 * This is deliberately run only after the final exact union and pair-catalogue
 * rebuild.  It recovers pair consequences from the complete selected affine
 * lattice and matches both WL and L1 integer values.  Node coincidence alone,
 * an ACTIVE ledger record, or a presearch NIS pass cannot authorise a writer
 * edge.  The caller supplies the exact posterior moment that will be cloned
 * and conditioned by the writer path.
 */
inline int zhangFinalizeHistoricalProductPairAuthorization(
	ZhangProductIntegerConstraintSet& constraints,
	const std::string& validatedDeliveryMomentId,
	double totalFailureProbabilityBudget = 1e-3)
{
	for (auto& pair : constraints.dualFrequencyCertifiedPairs)
	{
		pair.currentPosteriorReauthorized = false;
		pair.exactHistoricalTransportWitness = false;
		pair.productGaugeComponentVersion = 0;
		pair.authorizedBackendBasisGeneration = 0;
		pair.productGaugeSnapshotId.clear();
		pair.validatedDeliveryMomentId.clear();
		pair.authorizedPhaseSegmentFingerprint.clear();
		pair.authorizedFamilyIds.clear();
		pair.authorizedDecisionIds.clear();
		pair.authorizedFailureProbability = 1;
	}
	const int dimension = constraints.productCoordinateDimension;
	const auto& historicalRows =
		constraints.currentReauthorizedHistoricalDualRows;
	const auto& historicalValues =
		constraints.currentReauthorizedHistoricalDualValues;
	const auto proofClosure = zhangDecisionRiskClosure(
		constraints.currentReauthorizedHistoricalDecisionProofs);
	const double chargedRisk = std::max(
		constraints.currentReauthorizedHistoricalFailureProbability,
		proofClosure.valid ? proofClosure.bound : 1.0);
	if (dimension <= 0 || historicalRows.empty() ||
		historicalRows.size() != historicalValues.size() ||
		constraints.currentReauthorizedHistoricalFamilyIds.empty() ||
		constraints.currentReauthorizedHistoricalSnapshotId.empty() ||
		constraints.currentReauthorizedSourceMomentId.empty() ||
		validatedDeliveryMomentId.empty() ||
		constraints.phaseSegmentFingerprint.empty() ||
		!constraints.currentReauthorizedHistoricalSegmentsValidated ||
		!constraints.currentReauthorizedHistoricalExactTransport ||
		constraints.currentReauthorizedHistoricalDecisionProofs.empty() ||
		!proofClosure.valid || !std::isfinite(chargedRisk) || chargedRisk < 0 ||
		chargedRisk > totalFailureProbabilityBudget + 1e-12)
		return 0;
	std::vector<std::string> familyIds =
		constraints.currentReauthorizedHistoricalFamilyIds;
	std::sort(familyIds.begin(), familyIds.end());
	familyIds.erase(std::unique(familyIds.begin(), familyIds.end()),
		familyIds.end());
	if (familyIds.empty() || std::any_of(familyIds.begin(), familyIds.end(),
		[](const std::string& id) { return id.empty(); })) return 0;
	std::vector<std::string> decisionIds;
	decisionIds.reserve(proofClosure.atoms.size());
	for (const auto& [id, proof] : proofClosure.atoms)
	{
		(void)proof;
		if (id.empty()) return 0;
		decisionIds.push_back(id);
	}
	if (decisionIds.empty()) return 0;

	bool dimensionsValid = false;
	const auto signalRows = zhangWideLaneFirstRowsToSignalRows(
		historicalRows, dimension, &dimensionsValid);
	if (!dimensionsValid) return 0;
	const auto support = zhangSplitJointProductIntegerSupport(
		signalRows, historicalValues, dimension);
	if (!support.valid) return 0;
	const auto historicalWide = zhangRecoverCertifiedPairRelations(
		support.wideLaneRows, support.wideLaneValues, dimension, true);
	const auto historicalFirst = zhangRecoverCertifiedPairRelations(
		support.firstSignalRows, support.firstSignalValues, dimension, true);
	using Edge = std::pair<int, int>;
	std::map<Edge, ZhangExactInteger> wideValues;
	std::map<Edge, ZhangExactInteger> firstValues;
	for (const auto& pair : historicalWide)
		wideValues[{pair.firstNode, pair.secondNode}] = pair.value;
	for (const auto& pair : historicalFirst)
		firstValues[{pair.firstNode, pair.secondNode}] = pair.value;

	int authorised = 0;
	for (auto& pair : constraints.dualFrequencyCertifiedPairs)
	{
		int first = pair.firstNode;
		int second = pair.secondNode;
		ZhangExactInteger value = pair.value;
		if (second < first)
		{
			std::swap(first, second);
			value = -value;
		}
		const Edge edge{first, second};
		const auto wide = wideValues.find(edge);
		const auto l1 = firstValues.find(edge);
		if (wide == wideValues.end() || l1 == firstValues.end() ||
			l1->second != value)
			continue;
		// The WL integer is not stored on the dual edge itself.  Require the
		// final certified-pair catalogue to contain the identical exact value.
		const bool finalWideValueMatches = std::any_of(
			constraints.certifiedPairs.begin(), constraints.certifiedPairs.end(),
			[&](const ZhangCertifiedPairRelation& candidate)
			{
				if (candidate.coordinate != "WL") return false;
				int candidateFirst = candidate.firstNode;
				int candidateSecond = candidate.secondNode;
				ZhangExactInteger candidateValue = candidate.value;
				if (candidateSecond < candidateFirst)
				{
					std::swap(candidateFirst, candidateSecond);
					candidateValue = -candidateValue;
				}
				return candidateFirst == first && candidateSecond == second &&
					candidateValue == wide->second;
			});
		if (!finalWideValueMatches) continue;
		pair.fromProductGaugeLedger = true;
		pair.currentPosteriorReauthorized = true;
		pair.exactHistoricalTransportWitness = true;
		pair.productGaugeComponentVersion =
			constraints.currentReauthorizedHistoricalComponentVersion;
		pair.authorizedBackendBasisGeneration =
			constraints.backendBasisGeneration;
		pair.productGaugeSnapshotId =
			constraints.currentReauthorizedHistoricalSnapshotId;
		pair.validatedDeliveryMomentId = validatedDeliveryMomentId;
		pair.authorizedPhaseSegmentFingerprint =
			constraints.phaseSegmentFingerprint;
		pair.authorizedFamilyIds = familyIds;
		pair.authorizedDecisionIds = decisionIds;
		pair.authorizedFailureProbability = chargedRisk;
		authorised++;
	}
	return authorised;
}

inline ZhangExactDualGraphScore zhangExactDualGraphScore(
	const ZhangExactMatrix& signalRows, const ZhangExactVector& values, int namedCount)
{
	return zhangExactDualGraphProfile(signalRows, values, namedCount).score;
}

/** Reference-invariant all-pair incidence for a star-coordinate ambient
 * lattice.  Columns are K_s-K_ref for s<namedCount and the implicit final
 * node is the canonical reference.  Every unordered satellite pair appears
 * once, so D Q D' and its trace are invariant to the chosen star reference. */
inline MatrixXd zhangAllPairIncidence(int namedCount)
{
	if (namedCount <= 0) return MatrixXd(0, 0);
	const int nodeCount = namedCount + 1;
	MatrixXd result = MatrixXd::Zero(
		nodeCount * (nodeCount - 1) / 2, namedCount);
	int row = 0;
	for (int first = 0; first < nodeCount; first++)
	for (int second = first + 1; second < nodeCount; second++)
	{
		if (first < namedCount) result(row, first) += 1;
		if (second < namedCount) result(row, second) -= 1;
		row++;
	}
	return result;
}

struct ZhangExactConditioningAudit
{
	bool valid = false;
	int effectiveRank = 0;
	double maximumNullInnovation = 0;
	double nis = std::numeric_limits<double>::quiet_NaN();
	VectorXd mean;
	MatrixXd covariance;
	MatrixXd reduction;
};

/** Condition a Gaussian product coordinate on exact integer rows A z=n.
 * Rank-deficient constraint covariance is handled with an eigen pseudo-
 * inverse; a nonzero innovation in its null space fails closed. */
inline ZhangExactConditioningAudit zhangConditionExactProductRows(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const MatrixXd& rows,
	const VectorXd& integers)
{
	ZhangExactConditioningAudit result;
	result.mean = mean;
	result.covariance = covariance;
	result.reduction = MatrixXd::Zero(covariance.rows(), covariance.cols());
	if (covariance.rows() != covariance.cols() || mean.size() != covariance.rows() ||
		rows.cols() != mean.size() || rows.rows() != integers.size()) return result;
	if (rows.rows() == 0)
	{
		result.valid = true;
		return result;
	}
	const MatrixXd symmetric = 0.5 * (covariance + covariance.transpose());
	MatrixXd constraint = rows * symmetric * rows.transpose();
	constraint = 0.5 * (constraint + constraint.transpose());
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(constraint);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
		return result;
	const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
	const double tolerance = std::max(1e-14, 1e-12 * largest);
	VectorXd inverse = VectorXd::Zero(rows.rows());
	const VectorXd innovation = integers - rows * mean;
	const VectorXd coordinates = eigen.eigenvectors().transpose() * innovation;
	for (int index = 0; index < rows.rows(); index++)
	{
		if (eigen.eigenvalues()(index) > tolerance)
		{
			inverse(index) = 1 / eigen.eigenvalues()(index);
			result.effectiveRank++;
		}
		else result.maximumNullInnovation = std::max(
			result.maximumNullInnovation, std::abs(coordinates(index)));
	}
	if (result.maximumNullInnovation > 1e-7) return result;
	const MatrixXd pseudoInverse = eigen.eigenvectors() * inverse.asDiagonal() *
		eigen.eigenvectors().transpose();
	const MatrixXd cross = symmetric * rows.transpose();
	result.reduction = cross * pseudoInverse * cross.transpose();
	result.mean = mean + cross * pseudoInverse * innovation;
	result.covariance = symmetric - result.reduction;
	result.covariance = 0.5 * (result.covariance + result.covariance.transpose());
	result.nis = innovation.dot(pseudoInverse * innovation);
	result.valid = result.mean.allFinite() && result.covariance.allFinite() &&
		result.reduction.allFinite() && std::isfinite(result.nis);
	return result;
}

inline double zhangReferenceInvariantPairTrace(const MatrixXd& covariance)
{
	if (covariance.rows() <= 0 || covariance.rows() != covariance.cols()) return 0;
	const MatrixXd pairs = zhangAllPairIncidence(covariance.rows());
	return (pairs * covariance * pairs.transpose()).trace();
}

struct ZhangComponentBridgeGls
{
	bool valid = false;
	int effectiveRank = 0;
	double mean = std::numeric_limits<double>::quiet_NaN();
	double variance = std::numeric_limits<double>::quiet_NaN();
	double residualNis = std::numeric_limits<double>::quiet_NaN();
	double maximumNullResidual = 0;
};

struct ZhangComponentGaugeGls
{
	bool valid = false;
	int measurementRank = 0;
	int gaugeRank = 0;
	Eigen::VectorXd mean;
	Eigen::MatrixXd covariance;
	Eigen::MatrixXd information;
	Eigen::VectorXd residual;
	std::vector<Eigen::VectorXd> nullModes;
	double residualNis = std::numeric_limits<double>::quiet_NaN();
	double maximumNullResidual = 0;
	std::string failureReason = "NOT_EVALUATED";
};

struct ZhangComponentGaugeProductRow
{
	bool valid = false;
	ZhangExactVector row;
	ZhangExactInteger value = 0;
};

/** Map one fixed datum-free component-gauge row back to the named satellite
 * product lattice, including the already-certified within-component offsets.
 * An anchor equal to namedDimension denotes the implicit canonical reference
 * and therefore contributes no explicit coordinate column. */
inline ZhangComponentGaugeProductRow zhangComponentGaugeToProductRow(
	const ZhangExactVector& gaugeCombination,
	const std::vector<int>& componentAnchors,
	const ZhangExactVector& componentAnchorPotentials,
	int namedDimension,
	const ZhangExactInteger& gaugeIntegerValue)
{
	ZhangComponentGaugeProductRow result;
	if (namedDimension <= 0 || componentAnchors.size() < 2 ||
		componentAnchors.size() != componentAnchorPotentials.size() ||
		gaugeCombination.size() + 1 != componentAnchors.size()) return result;
	if (std::any_of(componentAnchors.begin(), componentAnchors.end(),
		[namedDimension](int anchor)
		{ return anchor < 0 || anchor > namedDimension; })) return result;
	result.row = ZhangExactVector(namedDimension);
	result.value = gaugeIntegerValue;
	const int datumAnchor = componentAnchors.front();
	for (int gauge = 0; gauge < static_cast<int>(gaugeCombination.size()); gauge++)
	{
		const auto& coefficient = gaugeCombination[gauge];
		const int componentAnchor = componentAnchors[gauge + 1];
		if (componentAnchor < namedDimension)
			result.row[componentAnchor] += coefficient;
		if (datumAnchor < namedDimension)
			result.row[datumAnchor] -= coefficient;
		result.value += coefficient *
			(componentAnchorPotentials[gauge + 1] -
			 componentAnchorPotentials.front());
	}
	result.valid = std::any_of(result.row.begin(), result.row.end(),
		[](const auto& value) { return value != 0; });
	return result;
}

/** Joint GLS reduction of all correlated cross-component observations
 *
 *     y = D_C c + e,  c in Z^(K-1).
 *
 * The measurement covariance may be singular because the complete set of
 * satellite-pair edges is deliberately retained.  The returned gauge
 * covariance is usable for ILS/PAR only when every datum-free component gauge
 * is estimable.  Null-space disagreement fails closed instead of being
 * discarded by the pseudo inverse.
 */
inline ZhangComponentGaugeGls zhangComponentGaugeGls(
	const Eigen::VectorXd& measurements,
	const Eigen::MatrixXd& covariance,
	const Eigen::MatrixXd& design)
{
	ZhangComponentGaugeGls result;
	const int count = measurements.size();
	const int gauges = design.cols();
	if (count == 0 || gauges == 0 || design.rows() != count ||
		covariance.rows() != count || covariance.cols() != count ||
		!measurements.allFinite() || !design.allFinite() ||
		!covariance.allFinite())
	{
		result.failureReason = "COMPONENT_GAUGE_INPUT_INVALID";
		return result;
	}

	const Eigen::MatrixXd symmetric =
		0.5 * (covariance + covariance.transpose());
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetric);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
	{
		result.failureReason = "COMPONENT_GAUGE_COVARIANCE_EIGEN_FAILED";
		return result;
	}
	const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
	const double tolerance = std::max(1e-14, 1e-12 * largest);
	if (eigen.eigenvalues().minCoeff() < -tolerance)
	{
		result.failureReason = "COMPONENT_GAUGE_COVARIANCE_NOT_PSD";
		return result;
	}
	Eigen::VectorXd inverse = Eigen::VectorXd::Zero(count);
	for (int index = 0; index < count; index++)
	{
		if (eigen.eigenvalues()(index) <= tolerance) continue;
		inverse(index) = 1 / eigen.eigenvalues()(index);
		result.measurementRank++;
	}
	const Eigen::MatrixXd pseudoInverse =
		eigen.eigenvectors() * inverse.asDiagonal() *
		eigen.eigenvectors().transpose();
	const Eigen::MatrixXd information =
		design.transpose() * pseudoInverse * design;
	result.information = information;
	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> informationSolver(
		information);
	result.gaugeRank = informationSolver.rank();
	if (result.gaugeRank != gauges)
	{
		result.failureReason = "COMPONENT_GAUGE_DESIGN_RANK_DEFICIENT";
		return result;
	}
	result.covariance = informationSolver.solve(
		Eigen::MatrixXd::Identity(gauges, gauges));
	result.covariance = 0.5 *
		(result.covariance + result.covariance.transpose());
	result.mean = result.covariance * design.transpose() *
		pseudoInverse * measurements;
	const Eigen::VectorXd residual = measurements - design * result.mean;
	result.residual = residual;
	result.residualNis = residual.dot(pseudoInverse * residual);
	const Eigen::VectorXd nullResidual =
		residual - symmetric * pseudoInverse * residual;
	result.maximumNullResidual = nullResidual.lpNorm<Eigen::Infinity>();
	for (int mode = 0; mode < eigen.eigenvalues().size(); mode++)
	{
		if (eigen.eigenvalues()(mode) <= tolerance)
		{
			const auto& vector = eigen.eigenvectors().col(mode);
			if (std::abs(vector.dot(residual)) > tolerance)
				result.nullModes.push_back(vector);
		}
	}
	result.valid = result.mean.allFinite() && result.covariance.allFinite() &&
		std::isfinite(result.residualNis) &&
		result.maximumNullResidual <= 1e-7;
	result.failureReason = result.valid
		? "NONE"
		: (result.maximumNullResidual > 1e-7
			? "COMPONENT_GAUGE_COVARIANCE_NULLSPACE_CONFLICT"
			: "COMPONENT_GAUGE_NUMERICAL_FAILURE");
	return result;
}

/** Select a largest paired L1/L2 coordinate sublattice supported by the
 * information matrix.  QR/COD is used only to decide whether the *named
 * incidence coordinates* add rank; the returned coordinates are unit rows of
 * the component forest and hence primitive integer functions. */
inline std::vector<int> zhangMaxEstimableDualGaugeForest(
	const Eigen::MatrixXd& information,
	int gaugeCount)
{
	std::vector<int> selected;
	if (gaugeCount <= 0 || information.rows() != 2 * gaugeCount ||
		information.cols() != 2 * gaugeCount || !information.allFinite())
		return selected;
	auto matrixRank = [](const Eigen::MatrixXd& matrix)
	{
		if (matrix.size() == 0) return 0;
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(matrix);
		return static_cast<int>(cod.rank());
	};
	int rank = 0;
	for (int gauge = 0; gauge < gaugeCount; gauge++)
	{
		std::vector<int> columns;
		for (const int existing : selected)
		{
			columns.push_back(existing);
			columns.push_back(gaugeCount + existing);
		}
		columns.push_back(gauge);
		columns.push_back(gaugeCount + gauge);
		Eigen::MatrixXd candidate(columns.size(), columns.size());
		for (int row = 0; row < static_cast<int>(columns.size()); row++)
		for (int column = 0; column < static_cast<int>(columns.size()); column++)
			candidate(row, column) = information(columns[row], columns[column]);
		const int candidateRank = matrixRank(candidate);
		if (candidateRank >= rank + 2)
		{
			selected.push_back(gauge);
			rank = candidateRank;
		}
	}
	return selected;
}

/** Test whether one named datum-free component-incidence functional belongs
 * to the real estimable space of I.  This is an estimability diagnostic only:
 * callers still use the original integer incidence row as their constraint.
 */
inline bool zhangComponentGaugeFunctionalEstimable(
	const Eigen::MatrixXd& information,
	const Eigen::VectorXd& functional,
	double tolerance = 1e-8)
{
	if (information.rows() == 0 || information.rows() != information.cols() ||
		functional.size() != information.cols() || !information.allFinite() ||
		!functional.allFinite()) return false;
	const Eigen::MatrixXd symmetric =
		0.5 * (information + information.transpose());
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetric);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
		return false;
	const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
	const double cutoff = std::max(1e-14, 1e-12 * largest);
	Eigen::VectorXd inverse = Eigen::VectorXd::Zero(information.rows());
	for (int index = 0; index < inverse.size(); index++)
		if (eigen.eigenvalues()(index) > cutoff)
			inverse(index) = 1 / eigen.eigenvalues()(index);
	const Eigen::MatrixXd pseudoInverse = eigen.eigenvectors() *
		inverse.asDiagonal() * eigen.eigenvectors().transpose();
	const Eigen::VectorXd unestimable = functional -
		symmetric * pseudoInverse * functional;
	return unestimable.norm() <= tolerance * std::max(1.0, functional.norm());
}

struct ZhangComponentGaugeLeaveOneGroupOut
{
	bool valid = false;
	bool estimable = false;
	int removedCount = 0;
	int retainedCount = 0;
	double mean = std::numeric_limits<double>::quiet_NaN();
	double variance = std::numeric_limits<double>::quiet_NaN();
	ZhangComponentGaugeGls gauge;
};

/** Refit the component-gauge GLS after deleting one physical observation group.
 *
 * Indices refer to rows of the supplied measurement/design system.  Both the
 * observation covariance and the design matrix are subset before the solve;
 * this is therefore a genuine leave-one-group-out refit at the component-gauge
 * GLS level, not a zeroed posterior functional.
 */
inline ZhangComponentGaugeLeaveOneGroupOut
zhangComponentGaugeLeaveOneGroupOut(
	const Eigen::VectorXd& measurements,
	const Eigen::MatrixXd& covariance,
	const Eigen::MatrixXd& design,
	const std::vector<int>& removedIndices,
	const Eigen::VectorXd& functional)
{
	ZhangComponentGaugeLeaveOneGroupOut result;
	const int count = measurements.size();
	if (count <= 0 || design.rows() != count || covariance.rows() != count ||
		covariance.cols() != count || functional.size() != design.cols())
		return result;
	std::vector<bool> removed(count, false);
	for (const int index : removedIndices)
	{
		if (index < 0 || index >= count || removed[index]) continue;
		removed[index] = true;
		result.removedCount++;
	}
	result.retainedCount = count - result.removedCount;
	if (result.removedCount == 0 || result.retainedCount == 0) return result;
	Eigen::VectorXd retainedMeasurements(result.retainedCount);
	Eigen::MatrixXd retainedDesign(result.retainedCount, design.cols());
	Eigen::MatrixXd retainedCovariance(
		result.retainedCount, result.retainedCount);
	std::vector<int> retained;
	retained.reserve(result.retainedCount);
	for (int index = 0; index < count; index++)
		if (!removed[index]) retained.push_back(index);
	for (int row = 0; row < result.retainedCount; row++)
	{
		retainedMeasurements(row) = measurements(retained[row]);
		retainedDesign.row(row) = design.row(retained[row]);
		for (int column = 0; column < result.retainedCount; column++)
			retainedCovariance(row, column) =
				covariance(retained[row], retained[column]);
	}
	result.gauge = zhangComponentGaugeGls(
		retainedMeasurements, retainedCovariance, retainedDesign);
	if (!result.gauge.valid) return result;
	result.estimable = zhangComponentGaugeFunctionalEstimable(
		result.gauge.information, functional);
	if (!result.estimable) return result;
	result.mean = functional.dot(result.gauge.mean);
	result.variance = (functional.transpose() * result.gauge.covariance *
		functional)(0, 0);
	result.valid = std::isfinite(result.mean) &&
		std::isfinite(result.variance) && result.variance >= 0;
	return result;
}

/** Return a maximum-information forest of dual-frequency estimable component
 * differences.  QR/eigen analysis here is never converted to an integer row:
 * every returned edge is the primitive incidence e_a-e_b, for both signals.
 */
inline std::vector<ZhangComponentEdgeId>
zhangMaximumEstimableDualComponentForest(
	const Eigen::MatrixXd& information,
	int componentCount,
	const std::vector<ZhangComponentEdgeId>& candidates,
	double estimabilityTolerance = 1e-8)
{
	std::vector<ZhangComponentEdgeId> forest;
	const int gauges = componentCount - 1;
	if (componentCount <= 1 || information.rows() != 2 * gauges ||
		information.cols() != 2 * gauges) return forest;

	const Eigen::MatrixXd symmetric =
		0.5 * (information + information.transpose());
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetric);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
		return forest;
	const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
	const double cutoff = std::max(1e-14, 1e-12 * largest);
	Eigen::VectorXd inverse = Eigen::VectorXd::Zero(2 * gauges);
	for (int index = 0; index < inverse.size(); index++)
		if (eigen.eigenvalues()(index) > cutoff)
			inverse(index) = 1 / eigen.eigenvalues()(index);
	const Eigen::MatrixXd pseudoInverse = eigen.eigenvectors() *
		inverse.asDiagonal() * eigen.eigenvectors().transpose();
	const double varianceTolerance = std::max(1e-14,
		1e-12 * std::max(1.0, pseudoInverse.norm()));

	struct WeightedEdge
	{
		ZhangComponentEdgeId edge;
		Eigen::VectorXd wideLaneFunctional;
	};
	std::vector<WeightedEdge> eligible;
	for (const auto& candidate : candidates)
	{
		const int first = candidate.firstComponent;
		const int second = candidate.secondComponent;
		if (first < 0 || second < 0 || first >= componentCount ||
			second >= componentCount || first == second) continue;
		Eigen::VectorXd firstSignal = Eigen::VectorXd::Zero(2 * gauges);
		Eigen::VectorXd secondSignal = Eigen::VectorXd::Zero(2 * gauges);
		if (first > 0)
		{
			firstSignal(first - 1) += 1;
			secondSignal(gauges + first - 1) += 1;
		}
		if (second > 0)
		{
			firstSignal(second - 1) -= 1;
			secondSignal(gauges + second - 1) -= 1;
		}
		if (!zhangComponentGaugeFunctionalEstimable(
			information, firstSignal, estimabilityTolerance) ||
			!zhangComponentGaugeFunctionalEstimable(
				information, secondSignal, estimabilityTolerance)) continue;
		// WL is g1-g2.  Its variance therefore contains the full
		// Q11+Q22-Q12-Q21 expression; summing two marginal variances would
		// silently discard the cross-frequency covariance.
		const Eigen::VectorXd wideLane = firstSignal - secondSignal;
		const double wideLaneVariance =
			(wideLane.transpose() * pseudoInverse * wideLane)(0, 0);
		if (!std::isfinite(wideLaneVariance) ||
			wideLaneVariance <= varianceTolerance)
			continue;
		eligible.push_back({candidate, wideLane});
	}
	std::vector<int> parent(componentCount);
	std::iota(parent.begin(), parent.end(), 0);
	auto root = [&parent](int node)
	{
		int value = node;
		while (parent[value] != value)
		{
			parent[value] = parent[parent[value]];
			value = parent[value];
		}
		return value;
	};
	std::vector<Eigen::VectorXd> selectedFunctionals;
	std::set<std::string> usedReceivers;
	std::set<std::string> usedArcs;
	while (forest.size() + 1 < static_cast<std::size_t>(componentCount))
	{
		int bestIndex = -1;
		double bestWeakestInformation = -1;
		double bestConditionalInformation = -1;
		int bestReceiverNovelty = -1;
		int bestArcNovelty = -1;
		for (int index = 0; index < static_cast<int>(eligible.size()); index++)
		{
			const auto& weighted = eligible[index];
			if (root(weighted.edge.firstComponent) ==
				root(weighted.edge.secondComponent)) continue;
			std::vector<Eigen::VectorXd> trial = selectedFunctionals;
			trial.push_back(weighted.wideLaneFunctional);
			Eigen::MatrixXd trialCovariance(trial.size(), trial.size());
			for (int row = 0; row < static_cast<int>(trial.size()); row++)
			for (int column = 0; column < static_cast<int>(trial.size()); column++)
				trialCovariance(row, column) = (trial[row].transpose() *
					pseudoInverse * trial[column])(0, 0);
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> trialEigen(
				0.5 * (trialCovariance + trialCovariance.transpose()));
			if (trialEigen.info() != Eigen::Success ||
				!trialEigen.eigenvalues().allFinite() ||
				trialEigen.eigenvalues().minCoeff() <= varianceTolerance) continue;
			const double weakestInformation =
				1 / trialEigen.eigenvalues().maxCoeff();
			const int previous = static_cast<int>(selectedFunctionals.size());
			double conditionalVariance = trialCovariance(previous, previous);
			if (previous > 0)
			{
				const Eigen::MatrixXd held =
					trialCovariance.topLeftCorner(previous, previous);
				Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> solver(held);
				const Eigen::VectorXd cross =
					trialCovariance.topRightCorner(previous, 1);
				conditionalVariance -= cross.dot(solver.solve(cross));
			}
			if (!std::isfinite(conditionalVariance) ||
				conditionalVariance <= varianceTolerance)
				continue;
			const double conditionalInformation = 1 / conditionalVariance;
			int receiverNovelty = 0;
			for (const auto& receiver : weighted.edge.receiverSupport)
				receiverNovelty += !usedReceivers.contains(receiver);
			int arcNovelty = 0;
			for (const auto& arc : weighted.edge.arcSupport)
				arcNovelty += !usedArcs.contains(arc);
			const bool better = weakestInformation > bestWeakestInformation + 1e-14 ||
				(std::abs(weakestInformation - bestWeakestInformation) <= 1e-14 &&
				 (receiverNovelty > bestReceiverNovelty ||
				  (receiverNovelty == bestReceiverNovelty &&
				   (arcNovelty > bestArcNovelty ||
				    (arcNovelty == bestArcNovelty &&
				     conditionalInformation > bestConditionalInformation)))));
			if (better)
			{
				bestIndex = index;
				bestWeakestInformation = weakestInformation;
				bestConditionalInformation = conditionalInformation;
				bestReceiverNovelty = receiverNovelty;
				bestArcNovelty = arcNovelty;
			}
		}
		if (bestIndex < 0) break;
		const auto selected = eligible[bestIndex];
		const int left = root(selected.edge.firstComponent);
		const int right = root(selected.edge.secondComponent);
		parent[right] = left;
		forest.push_back(selected.edge);
		selectedFunctionals.push_back(selected.wideLaneFunctional);
		usedReceivers.insert(selected.edge.receiverSupport.begin(),
			selected.edge.receiverSupport.end());
		usedArcs.insert(selected.edge.arcSupport.begin(),
			selected.edge.arcSupport.end());
	}
	return forest;
}

/** Partition a dual-frequency component support graph.  A support edge is
 * admitted only when it has simultaneous L1 and L2 evidence; a missing edge
 * cannot make a different block fail. */
inline std::vector<std::vector<int>> zhangDualComponentSupportBlocks(
	int componentCount,
	const std::vector<ZhangComponentEdgeId>& dualSupportEdges)
{
	std::vector<std::vector<int>> blocks;
	if (componentCount <= 0) return blocks;
	std::vector<int> parent(componentCount);
	std::iota(parent.begin(), parent.end(), 0);
	auto root = [&](int node)
	{
		int value = node;
		while (parent[value] != value)
		{
			parent[value] = parent[parent[value]];
			value = parent[value];
		}
		return value;
	};
	for (const auto& edge : dualSupportEdges)
	{
		if (edge.firstComponent < 0 || edge.secondComponent < 0 ||
			edge.firstComponent >= componentCount ||
			edge.secondComponent >= componentCount ||
			edge.firstComponent == edge.secondComponent) continue;
		const int left = root(edge.firstComponent);
		const int right = root(edge.secondComponent);
		if (left != right) parent[right] = left;
	}
	std::map<int, std::vector<int>> grouped;
	for (int component = 0; component < componentCount; component++)
		grouped[root(component)].push_back(component);
	for (auto& [_, block] : grouped)
	{
		std::sort(block.begin(), block.end());
		blocks.push_back(std::move(block));
	}
	return blocks;
}

/** Map an exact covariance-nullspace contradiction back to the observations
 * that dominate it.  This performs localization only; callers must quarantine
 * the reported edge group rather than relax a numerical tolerance. */
inline std::vector<ZhangNullConflict> zhangLocalizeNullConflicts(
	const ZhangComponentGaugeGls& gls,
	const std::vector<ZhangComponentEdgeId>& edges,
	const std::vector<E_ObsCode>& signals = {},
	const std::vector<uint64_t>& generations = {},
	const std::vector<std::string>& segments = {},
	const std::vector<std::string>& sources = {},
	const std::vector<ZhangExactInteger>& affineOffsets = {})
{
	std::vector<ZhangNullConflict> conflicts;
	if (gls.residual.size() != static_cast<int>(edges.size())) return conflicts;
	for (const auto& mode : gls.nullModes)
	{
		if (mode.size() != gls.residual.size()) continue;
		ZhangNullConflict conflict;
		conflict.nullResidual = mode.dot(gls.residual);
		std::vector<int> order(mode.size());
		std::iota(order.begin(), order.end(), 0);
		std::sort(order.begin(), order.end(), [&mode, &gls](int left, int right)
			{ return std::abs(mode(left) * gls.residual(left)) >
				std::abs(mode(right) * gls.residual(right)); });
		for (const int index : order)
		{
			const double contribution = mode(index) * gls.residual(index);
			if (std::abs(contribution) <= 1e-10) break;
			conflict.dominantEdges.push_back(edges[index]);
			conflict.contributions.push_back(contribution);
			conflict.coefficients.push_back(mode(index));
			if (index < static_cast<int>(signals.size()))
				conflict.signals.push_back(signals[index]);
			if (index < static_cast<int>(generations.size()))
				conflict.backendGenerations.push_back(generations[index]);
			if (index < static_cast<int>(segments.size()))
				conflict.phaseSegments.push_back(segments[index]);
			if (index < static_cast<int>(sources.size()))
				conflict.sources.push_back(sources[index]);
			if (index < static_cast<int>(affineOffsets.size()))
				conflict.affineOffsets.push_back(affineOffsets[index]);
			if (conflict.dominantEdges.size() == 8) break;
		}
		if (!conflict.dominantEdges.empty()) conflicts.push_back(std::move(conflict));
	}
	return conflicts;
}

/** GLS estimate of one shared relative integer gauge from correlated cross-
 * component edge measurements y=1*d+e. */
inline ZhangComponentBridgeGls zhangComponentBridgeGls(
	const VectorXd& measurements, const MatrixXd& covariance)
{
	ZhangComponentBridgeGls result;
	const int count = measurements.size();
	if (count == 0 || covariance.rows() != count || covariance.cols() != count)
		return result;
	const MatrixXd symmetric = 0.5 * (covariance + covariance.transpose());
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(symmetric);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
		return result;
	const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
	const double tolerance = std::max(1e-14, 1e-12 * largest);
	VectorXd inverse = VectorXd::Zero(count);
	for (int index = 0; index < count; index++)
	{
		if (eigen.eigenvalues()(index) > tolerance)
		{
			inverse(index) = 1 / eigen.eigenvalues()(index);
			result.effectiveRank++;
		}
	}
	const MatrixXd pseudoInverse = eigen.eigenvectors() * inverse.asDiagonal() *
		eigen.eigenvectors().transpose();
	const VectorXd ones = VectorXd::Ones(count);
	const double information = ones.dot(pseudoInverse * ones);
	if (!(information > 0) || !std::isfinite(information)) return result;
	result.variance = 1 / information;
	result.mean = result.variance * ones.dot(pseudoInverse * measurements);
	const VectorXd residual = measurements - ones * result.mean;
	result.residualNis = residual.dot(pseudoInverse * residual);
	const VectorXd nullResidual = residual - symmetric * pseudoInverse * residual;
	result.maximumNullResidual = nullResidual.lpNorm<Eigen::Infinity>();
	result.valid = std::isfinite(result.mean) && std::isfinite(result.variance) &&
		result.variance > 0 && std::isfinite(result.residualNis) &&
		result.maximumNullResidual <= 1e-7;
	return result;
}

/** Deterministic Kruskal forest over statistically admissible pair edges.
 * Perr is the primary order and variance only breaks ties.  The routine never
 * relaxes the reliability ceiling and therefore cannot turn a weak edge into
 * an integer certificate. */
inline std::vector<ZhangPairReliabilityEdge> zhangPairReliabilityForest(
	int nodeCount,
	std::vector<ZhangPairReliabilityEdge> edges,
	double maximumPerr)
{
	std::vector<ZhangPairReliabilityEdge> forest;
	if (nodeCount < 2 || maximumPerr < 0) return forest;
	std::sort(edges.begin(), edges.end(), [](const auto& left, const auto& right)
	{
		if (left.perr != right.perr) return left.perr < right.perr;
		if (left.variance != right.variance)
			return left.variance < right.variance;
		if (left.firstNode != right.firstNode)
			return left.firstNode < right.firstNode;
		return left.secondNode < right.secondNode;
	});
	std::vector<int> parent(nodeCount);
	std::iota(parent.begin(), parent.end(), 0);
	auto root = [&](int node)
	{
		int value = node;
		while (parent[value] != value) value = parent[value];
		while (parent[node] != node)
		{
			const int next = parent[node];
			parent[node] = value;
			node = next;
		}
		return value;
	};
	for (const auto& edge : edges)
	{
		if (!std::isfinite(edge.perr) || edge.perr > maximumPerr ||
			edge.firstNode < 0 || edge.secondNode < 0 ||
			edge.firstNode >= nodeCount || edge.secondNode >= nodeCount)
		{
			continue;
		}
		int first = root(edge.firstNode);
		int second = root(edge.secondNode);
		if (first == second) continue;
		parent[second] = first;
		forest.push_back(edge);
	}
	return forest;
}

/** Completion value of a statistically accepted product constraint set.
 *
 * Candidate generation and the named-pair beam already prefer rows that add
 * missing dual-frequency or persistent-quotient rank.  The outer solver must
 * use the same ordering when it compares whole L1 alternatives; otherwise a
 * larger, but already represented, pair forest can overwrite the smaller set
 * that actually connects a new satellite component.
 *
 * The score is exact and affine-aware.  It first forms the certified HNF union
 * with the pre-existing WL/L1 lattices, then measures (a) the additional rank
 * of the common named-pair graph and (b) the additional rank of the joint
 * [WL,L1] lattice.  No reliability decision is made here: callers may score
 * only candidates that have already passed their Perr and NIS gates. */
struct ZhangProductConstraintCompletionScore
{
	bool valid = false;
	int existingDualGraphRank = 0;
	int completedDualGraphRank = 0;
	int dualGraphRankGain = 0;
	int existingJointLatticeRank = 0;
	int completedJointLatticeRank = 0;
	int quotientRankGain = 0;
	std::string failureReason = "NOT_EVALUATED";
};

inline ZhangProductConstraintCompletionScore
zhangScoreProductConstraintCompletion(
	const ZhangExactMatrix& existingWideLaneRows,
	const ZhangExactVector& existingWideLaneValues,
	const ZhangExactMatrix& existingFirstSignalRows,
	const ZhangExactVector& existingFirstSignalValues,
	const ZhangExactMatrix& candidateWideLaneRows,
	const ZhangExactVector& candidateWideLaneValues,
	const ZhangExactMatrix& candidateFirstSignalRows,
	const ZhangExactVector& candidateFirstSignalValues,
	int productDimension,
	int precomputedExistingDualGraphRank = -1)
{
	ZhangProductConstraintCompletionScore result;
	auto dimensionsValid = [&](const ZhangExactMatrix& rows,
		const ZhangExactVector& values)
	{
		if (rows.size() != values.size()) return false;
		return std::all_of(rows.begin(), rows.end(), [&](const auto& row)
		{
			return row.size() == static_cast<std::size_t>(productDimension);
		});
	};
	if (productDimension <= 0 ||
		!dimensionsValid(existingWideLaneRows, existingWideLaneValues) ||
		!dimensionsValid(existingFirstSignalRows, existingFirstSignalValues) ||
		!dimensionsValid(candidateWideLaneRows, candidateWideLaneValues) ||
		!dimensionsValid(candidateFirstSignalRows, candidateFirstSignalValues))
	{
		result.failureReason = "PRODUCT_COMPLETION_SCORE_DIMENSION_MISMATCH";
		return result;
	}

	const auto target = zhangExactIdentityMatrix(productDimension);
	const auto wideLaneUnion = zhangExactCertifiedUnionAudit(
		target,
		existingWideLaneRows, existingWideLaneValues,
		candidateWideLaneRows, candidateWideLaneValues);
	const auto firstSignalUnion = zhangExactCertifiedUnionAudit(
		target,
		existingFirstSignalRows, existingFirstSignalValues,
		candidateFirstSignalRows, candidateFirstSignalValues);
	if (!wideLaneUnion.consistent ||
		!wideLaneUnion.certifiedContainedInTarget ||
		!firstSignalUnion.consistent ||
		!firstSignalUnion.certifiedContainedInTarget)
	{
		result.failureReason = "PRODUCT_COMPLETION_SCORE_AFFINE_UNION_FAILED";
		return result;
	}

	auto dualGraphRank = [&](const ZhangExactMatrix& wideRows,
		const ZhangExactVector& wideValues,
		const ZhangExactMatrix& firstRows,
		const ZhangExactVector& firstValues)
	{
		const auto widePairs = zhangRecoverCertifiedPairRelations(
			wideRows, wideValues, productDimension, true);
		const auto firstPairs = zhangRecoverCertifiedPairRelations(
			firstRows, firstValues, productDimension, true);
		std::set<std::pair<int, int>> wideIds;
		for (const auto& pair : widePairs)
			wideIds.insert({pair.firstNode, pair.secondNode});
		std::vector<ZhangPairReliabilityEdge> common;
		for (const auto& pair : firstPairs)
		{
			if (wideIds.contains({pair.firstNode, pair.secondNode}))
				common.push_back({pair.firstNode, pair.secondNode, 0, 0});
		}
		return static_cast<int>(zhangPairReliabilityForest(
			productDimension + 1, std::move(common), 0).size());
	};

	result.existingDualGraphRank = precomputedExistingDualGraphRank >= 0
		? precomputedExistingDualGraphRank
		: dualGraphRank(
			existingWideLaneRows, existingWideLaneValues,
			existingFirstSignalRows, existingFirstSignalValues);
	const bool candidateEmpty = candidateWideLaneRows.empty() &&
		candidateFirstSignalRows.empty();
	result.completedDualGraphRank = candidateEmpty
		? result.existingDualGraphRank
		: dualGraphRank(
			wideLaneUnion.certifiedBasis, wideLaneUnion.certifiedValues,
			firstSignalUnion.certifiedBasis, firstSignalUnion.certifiedValues);
	result.dualGraphRankGain = std::max(
		0, result.completedDualGraphRank - result.existingDualGraphRank);
	// WL [q,-q] and L1 [r,0] subspaces intersect only at zero, so their
	// exact ranks add in the joint two-frequency coordinate.
	result.existingJointLatticeRank =
		wideLaneUnion.heldRank + firstSignalUnion.heldRank;
	result.completedJointLatticeRank = static_cast<int>(
		wideLaneUnion.certifiedBasis.size() +
		firstSignalUnion.certifiedBasis.size());
	result.quotientRankGain = std::max(
		0, result.completedJointLatticeRank -
		result.existingJointLatticeRank);
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline bool zhangProductCompletionScoreBetter(
	const ZhangProductConstraintCompletionScore& left,
	const ZhangProductConstraintCompletionScore& right)
{
	if (left.valid != right.valid) return left.valid;
	if (left.dualGraphRankGain != right.dualGraphRankGain)
		return left.dualGraphRankGain > right.dualGraphRankGain;
	if (left.quotientRankGain != right.quotientRankGain)
		return left.quotientRankGain > right.quotientRankGain;
	return false;
}

inline ZhangInheritedNamedCertificate
zhangPromoteNamedCertificateFromAcceptedParent(
	const ZhangExactMatrix& parentFixedRows,
	const ZhangExactVector& parentFixedValues,
	int namedCount,
	bool parentStatisticallyAccepted)
{
	ZhangInheritedNamedCertificate result;
	result.parentFixedRank = static_cast<int>(parentFixedRows.size());
	if (!parentStatisticallyAccepted || namedCount <= 0 ||
		parentFixedRows.empty() ||
		parentFixedRows.size() != parentFixedValues.size())
	{
		return result;
	}
	result.values = zhangRecoverCertifiedNamedProductSubset(
		parentFixedRows, parentFixedValues, namedCount);
	result.exact = !result.values.empty();
	return result;
}

/** Primitive satellite-difference basis that removes one component-common
 * real gauge.  Different anchors are coordinate choices only: their exact row
 * HNFs must agree. */
inline ZhangExactMatrix zhangComponentRelativeGaugeBasis(
	int satelliteCount, int anchorIndex = 0)
{
	ZhangExactMatrix result;
	if (satelliteCount < 2 || anchorIndex < 0 ||
		anchorIndex >= satelliteCount)
	{
		return result;
	}
	for (int satellite = 0; satellite < satelliteCount; satellite++)
	{
		if (satellite == anchorIndex) continue;
		ZhangExactVector row(satelliteCount);
		row[satellite] = 1;
		row[anchorIndex] = -1;
		result.push_back(std::move(row));
	}
	return result;
}

/** Generate one bounded outer named-PAR child without enumerating all n
 * deletions.  An exact recoverable mixed-LAMBDA seed is preferred; otherwise
 * one deterministic backward-elimination child is returned.  Starting from a
 * single full branch therefore reaches minimumRank in at most fullRank steps. */
inline std::vector<std::vector<int>> zhangProductNamedBackwardChildren(
	const std::vector<int>& selected,
	const std::vector<int>& recoverableLocalIndices,
	int worstNamedLocalPosition,
	int minimumRank)
{
	std::vector<std::vector<int>> children;
	if (static_cast<int>(selected.size()) <= minimumRank)
	{
		return children;
	}
	if (!recoverableLocalIndices.empty() &&
		recoverableLocalIndices.size() < selected.size())
	{
		std::vector<int> exactSeed;
		for (int local : recoverableLocalIndices)
		{
			if (local >= 0 && local < static_cast<int>(selected.size()))
			{
				exactSeed.push_back(selected.at(local));
			}
		}
		std::sort(exactSeed.begin(), exactSeed.end());
		exactSeed.erase(std::unique(exactSeed.begin(), exactSeed.end()),
			exactSeed.end());
		if (static_cast<int>(exactSeed.size()) >= minimumRank)
		{
			children.push_back(std::move(exactSeed));
			return children;
		}
	}
	if (worstNamedLocalPosition >= 0 &&
		worstNamedLocalPosition < static_cast<int>(selected.size()))
	{
		auto child = selected;
		child.erase(child.begin() + worstNamedLocalPosition);
		if (static_cast<int>(child.size()) >= minimumRank)
		{
			children.push_back(std::move(child));
		}
	}
	return children;
}

/** Identity-weighted product covariance gain from exact constraints Bq=k.
 *
 * With q=Lx and candidate rows A_S x, exact fixing gives
 *   P_fix = P - P A_S' (A_S P A_S')^-1 A_S P,
 *   Delta Q_q = L P A_S' (A_S P A_S')^-1 A_S P L'.
 * A general product weighting would use
 *   g(S)=tr(W_q Delta Q_q)/tr(W_q Q_q).
 * Here q is already the joint named product coordinate and B maps q to the
 * candidate WL/L1 relations.  The returned scalar is
 *   tr(Delta Q_q) / tr(Q_q).
 * No user geometry is introduced in this first implementation. */
inline double zhangNamedProductInformationGain(
	const MatrixXd& productCovariance,
	const ZhangIarFunctional& candidateRows)
{
	const double denominator = productCovariance.trace();
	const ZhangIarCovarianceCondition condition =
		zhangIarCovarianceCondition(productCovariance, candidateRows);
	if (!condition.valid || !(denominator > 0) ||
		!std::isfinite(denominator))
	{
		return 0;
	}
	const double reduction = condition.reductionFactor.squaredNorm();
	return std::clamp(reduction / denominator, 0.0, 1.0);
}

/** Separate a rank ceiling problem from an integer candidate alignment
 * problem.  The thresholds only classify diagnostics; they never authorize
 * fixing or alter the beam. */
inline std::string zhangProductGainSpectrumDiagnosis(
	double realSubspaceUpperBound,
	double namedIntegerSubsetGain,
	double desiredCoverage = 0.80,
	double minimumEfficiency = 0.25)
{
	if (!std::isfinite(realSubspaceUpperBound) ||
		!std::isfinite(namedIntegerSubsetGain) ||
		realSubspaceUpperBound < 0 || namedIntegerSubsetGain < 0)
	{
		return "GAIN_COMPARISON_INVALID";
	}
	if (realSubspaceUpperBound < desiredCoverage)
	{
		return "REAL_RANK_CEILING_LOW_INCREASE_RANK";
	}
	const double efficiency = realSubspaceUpperBound > 0
		? namedIntegerSubsetGain / realSubspaceUpperBound : 0;
	if (efficiency < minimumEfficiency)
	{
		return "INTEGER_CANDIDATE_SUBSPACE_MISALIGNED";
	}
	return "INTEGER_SUBSET_USES_REAL_CEILING_EFFICIENTLY";
}

/** Exact commuting-square check for a derived pair.  Product integers include
 * the affine product offset; the physical Ledger RHS must exclude it. */
inline bool zhangCertifiedPairPhysicalRhs(
    const ZhangExactVector& pairRow,
    const ZhangExactInteger& pairValue,
    const ZhangExactVector& parentCombination,
    int parentOffset,
    const ZhangExactMatrix& parentProductRows,
    const ZhangExactVector& parentPhysicalValues,
    const ZhangExactVector& affineOffsets,
    ZhangExactInteger& physicalValue)
{
    if (pairRow.empty() || pairRow.size() != affineOffsets.size() ||
        parentCombination.empty() || parentOffset < 0 ||
        parentProductRows.size() != parentPhysicalValues.size() ||
        static_cast<std::size_t>(parentOffset) + parentCombination.size() > parentProductRows.size())
        return false;
    ZhangExactVector reconstructed(pairRow.size());
    ZhangExactInteger rhs = 0;
    for (std::size_t i = 0; i < parentCombination.size(); ++i)
    {
        const auto& row = parentProductRows[parentOffset + i];
        if (row.size() != pairRow.size()) return false;
        for (std::size_t j = 0; j < pairRow.size(); ++j)
            reconstructed[j] += parentCombination[i] * row[j];
        rhs += parentCombination[i] * parentPhysicalValues[parentOffset + i];
    }
    if (reconstructed != pairRow) return false;
    ZhangExactInteger productValue = rhs;
    for (std::size_t j = 0; j < pairRow.size(); ++j)
        productValue += pairRow[j] * affineOffsets[j];
    if (productValue != pairValue) return false;
    physicalValue = std::move(rhs);
    return true;
}

/** Reserve the same total risk as the previous four fallback attempts.  The
 * attempt list is fixed before search; failed allocations are not recycled. */
inline double zhangConnectivityTargetFailureAllocation(double budget, std::size_t targetCount)
{
    if (!std::isfinite(budget) || budget <= 0 || targetCount == 0) return 0;
    return budget / (8.0 * static_cast<double>(std::max<std::size_t>(4, targetCount)));
}
