#pragma once

#include "common/zhangIntegerDecisionProof.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "common/enums.h"
#include "common/zhangIntegerCandidateNis.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangIntegerAudit.hpp"
#include "common/zhangR51Integrity.hpp"

enum class ZhangProductIntegerLedgerSource
{
	PRODUCT_ILS,
	DERIVED_PAIR,
	TEMPORAL_RECERTIFIED,
	BESD
};

inline const char* zhangProductIntegerLedgerSourceName(
	ZhangProductIntegerLedgerSource source)
{
	switch (source)
	{
		case ZhangProductIntegerLedgerSource::PRODUCT_ILS:
			return "PRODUCT_ILS";
		case ZhangProductIntegerLedgerSource::DERIVED_PAIR:
			return "DERIVED_PAIR";
		case ZhangProductIntegerLedgerSource::TEMPORAL_RECERTIFIED:
			return "TEMPORAL_RECERTIFIED";
		case ZhangProductIntegerLedgerSource::BESD:
			return "BESD";
	}
	return "UNKNOWN";
}

/** One product-lattice integer retained independently of the generic network
 * held lattice.  physicalExpansion is authoritative and includes observable,
 * receiver-satellite arc and arc version in every string key.  productRow is
 * only the coordinate used when the row was discovered. canonicalProductExpansion
 * is the immutable dual-frequency satellite-incidence functional used for
 * exact transport across product-tree and product-image changes. */
struct ProductIntegerLedgerRow
{
	ZhangDecisionProofs decisionProofs;
	E_Sys system = E_Sys::NONE;
	E_ObsCode firstObservable = E_ObsCode::NONE;
	E_ObsCode secondObservable = E_ObsCode::NONE;
	ZhangExactVector productRow;
	ZhangExactInteger integerValue = 0;
	std::map<std::string, ZhangExactInteger> physicalExpansion;
    bool physicalExpansionExact = false;
	std::map<std::string, ZhangExactInteger> canonicalProductExpansion;
	// Named edge metadata is populated only for exact pair certificates.  It is
	// retained separately from productRow so an offline full-oracle builder can
	// reconstruct a reference-invariant satellite graph after basis changes.
	std::string coordinate;
	std::string firstSatellite;
	std::string secondSatellite;
	std::string phaseSegmentFingerprint;
	std::uint64_t backendBasisGeneration = 0;
	// Family-wise failure bound of the jointly admitted product lattice from
	// which this row was derived.  A later-posterior consistency check must not
	// replace this fixed admission bound with a progressively smaller divisor
	// based on the number of redundant rows accumulated in the Ledger.
	double admissionFailureProbabilityBound = 1;
	long int firstCertified = 0;
	long int lastConfirmed = 0;
	int confirmationEpochs = 0;
	ZhangProductIntegerLedgerSource source =
		ZhangProductIntegerLedgerSource::PRODUCT_ILS;
	bool conditioningOnly = true;
	bool pairCertificate = false;
	bool certified = false;
};

/** Preserve the statistical strength of the original joint-lattice admission
 * when checking temporal consistency.  Repeated posteriors are correlated and
 * are not new hypothesis tests; the recheck may tighten the original bound but
 * must neither relax it nor divide it again by a growing redundant registry. */
inline double zhangProductLedgerPosteriorFailureBudget(
	const ProductIntegerLedgerRow& row,
	double configuredFamilyFailureBudget)
{
	const double family = std::clamp(
		std::isfinite(configuredFamilyFailureBudget)
			? configuredFamilyFailureBudget : 1e-3,
		1e-12, 1.0);
	if (!std::isfinite(row.admissionFailureProbabilityBound) ||
		row.admissionFailureProbabilityBound < 0 ||
		row.admissionFailureProbabilityBound > 1)
	{
		// An invalid admission receipt cannot be repaired with the current
		// configured budget: that would silently replace historical evidence
		// with a different statistical authority.
		return std::numeric_limits<double>::quiet_NaN();
	}
	return std::min(family,
		std::max(1e-12, row.admissionFailureProbabilityBound));
}

/** A mature Ledger row from another phase segment or outside the current
 * exact product image is not a competing integer for the current identity.
 * Those absence states must not poison a safe fresh baseline.  Every
 * compatibility, transport, posterior, decision-proof or affine-conflict
 * failure remains fail-closed by default. */
inline bool zhangProductLedgerPresearchAllowsFreshWithoutHistory(
	const std::string& status)
{
	return status == "NO_CERTIFIED_LEDGER_ROWS" ||
		status == "NO_CURRENT_SEGMENT_PROJECTABLE_ROWS" ||
		status == "NO_LEDGER_ROWS_IN_CURRENT_PRODUCT_IMAGE";
}

struct ZhangProductIntegerLedgerUpdate
{
    ZhangR51ConflictWitness exactConflict;
    std::vector<std::string> conflictRows;
	bool valid = false;
	int inputRows = 0;
	int freshRows = 0;
	int confirmedRows = 0;
	int conflictingRows = 0;
	int activeRankBefore = 0;
	int activeRankAfter = 0;
	std::string failureReason = "NOT_EVALUATED";
    int conflictingCandidateIndex = -1;
    ZhangExactInteger incumbentInteger = 0;
    ZhangExactInteger proposedInteger = 0;
    bool conflictPhysicalRowsEqual = false;
    std::uint64_t incumbentBackendGeneration = 0;
    std::uint64_t proposedBackendGeneration = 0;
    std::string conflictIdentity;
    std::string incumbentPhysicalRow;
    std::string proposedPhysicalRow;
};

/** The product writer is the final publication boundary.  A transaction that
 * was rejected, or one that had to restart an unconfirmed integer hypothesis,
 * cannot authorise PRODUCT_FIXED in the same epoch.  This keeps publication
 * atomic with respect to the complete current proof family. */
inline bool zhangProductLedgerWriterCommitAuthorized(
	const ZhangProductIntegerLedgerUpdate& update)
{
	return update.valid && update.conflictingRows == 0;
}

struct ZhangProductIntegerPosteriorRecheck
{
	bool valid = false;
	bool sameInteger = false;
	bool deterministic = false;
	bool reliable = false;
	ZhangExactInteger nearestInteger = 0;
	double innovation = std::numeric_limits<double>::quiet_NaN();
	double variance = std::numeric_limits<double>::quiet_NaN();
	double failureProbability = 1;
	double nis = std::numeric_limits<double>::quiet_NaN();
	double nisThreshold = std::numeric_limits<double>::quiet_NaN();
	std::string failureReason = "NOT_EVALUATED";
};

/** Recheck one previously accepted exact physical integer on a later
 * posterior without conditioning that posterior.  This is a confirmation
 * gate, not repeated-data probability multiplication: the original ILS
 * certificate remains mandatory and the later epoch must independently keep
 * the same nearest integer under scalar Perr and NIS limits. */
inline ZhangProductIntegerPosteriorRecheck
zhangRecheckProductIntegerOnPosterior(
	double floatingValue,
	double variance,
	const ZhangExactInteger& expectedInteger,
	double maximumFailureProbability,
	double nisAlpha)
{
	ZhangProductIntegerPosteriorRecheck result;
	result.variance = variance;
	if (!std::isfinite(floatingValue) || !std::isfinite(variance) ||
		variance < -1e-12 || maximumFailureProbability <= 0 ||
		maximumFailureProbability >= 1 || nisAlpha <= 0 || nisAlpha >= 1 ||
		floatingValue < static_cast<double>(std::numeric_limits<long long>::min()) ||
		floatingValue > static_cast<double>(std::numeric_limits<long long>::max()))
	{
		result.failureReason = "PRODUCT_LEDGER_RECHECK_INPUT_INVALID";
		return result;
	}
	result.nearestInteger = std::llround(floatingValue);
	result.sameInteger = result.nearestInteger == expectedInteger;
	const double expected = expectedInteger.convert_to<double>();
	result.innovation = expected - floatingValue;
	const double fractional = floatingValue -
		result.nearestInteger.convert_to<double>();
	result.failureProbability = zhangIntegerRoundFailureProbability(
		fractional, variance);
	if (variance <= 1e-12)
	{
		result.valid = true;
		result.deterministic = true;
		result.nis = 0;
		result.nisThreshold = std::numeric_limits<double>::infinity();
		result.reliable = result.sameInteger &&
			std::abs(result.innovation) <= 1e-7;
		result.failureReason = result.reliable
			? "NONE" : "DETERMINISTIC_AFFINE_CONFLICT";
		return result;
	}
	VectorXd innovation(1);
	innovation(0) = result.innovation;
	MatrixXd covariance(1, 1);
	covariance(0, 0) = variance;
	const auto nis = assessZhangIntegerCandidateNis(
		innovation, covariance, nisAlpha);
	result.valid = nis.valid;
	result.nis = nis.nis;
	result.nisThreshold = nis.threshold;
	result.reliable = result.sameInteger && nis.valid &&
		nis.nis <= nis.threshold &&
		result.failureProbability <= maximumFailureProbability + 1e-12;
	result.failureReason = result.reliable ? "NONE" :
		(!result.sameInteger ? "INTEGER_CHANGED" :
		 (!nis.valid ? "NIS_INVALID" :
		  (nis.nis > nis.threshold ? "NIS_REJECTED" :
		   "FAILURE_PROBABILITY_EXCEEDED")));
	return result;
}

enum class ZhangExactAffineMembershipStatus
{
	INVALID,
	EXACT_NEW,
	EXACT_REDUNDANT,
	AFFINE_CONFLICT
};

struct ZhangExactAffineMembership
{
	ZhangExactAffineMembershipStatus status =
		ZhangExactAffineMembershipStatus::INVALID;
	ZhangExactInteger impliedValue = 0;
	ZhangExactVector combination;
};

/** Classify q'x=n against an existing exact affine integer lattice.
 *
 * Row membership alone is insufficient: a row already generated by Q may
 * carry an inconsistent right-hand side.  Such an affine conflict must never
 * be routed to a numerical conditioner or hidden by an HNF rank check. */
inline ZhangExactAffineMembership zhangExactAffineMembership(
	const ZhangExactMatrix& rows,
	const ZhangExactVector& values,
	const ZhangExactVector& candidateRow,
	const ZhangExactInteger& candidateValue)
{
	ZhangExactAffineMembership result;
	if (rows.size() != values.size() || candidateRow.empty()) return result;
	for (const auto& row : rows)
		if (row.size() != candidateRow.size()) return result;
	if (rows.empty())
	{
		result.status = ZhangExactAffineMembershipStatus::EXACT_NEW;
		return result;
	}
	const auto membership = zhangIntegerRowLatticeContains(rows, candidateRow);
	if (!membership.contained)
	{
		result.status = ZhangExactAffineMembershipStatus::EXACT_NEW;
		return result;
	}
	if (membership.combination.size() != values.size()) return result;
	result.combination = membership.combination;
	for (int row = 0; row < static_cast<int>(values.size()); row++)
		result.impliedValue += result.combination[row] * values[row];
	result.status = result.impliedValue == candidateValue
		? ZhangExactAffineMembershipStatus::EXACT_REDUNDANT
		: ZhangExactAffineMembershipStatus::AFFINE_CONFLICT;
	return result;
}

inline std::string zhangProductPhysicalRowFingerprint(
	const std::map<std::string, ZhangExactInteger>& row)
{
	std::ostringstream stream;
	for (const auto& [identity, coefficient] : row)
	{
		if (coefficient == 0) continue;
		stream << identity << "=" << coefficient << ";";
	}
	return stream.str();
}

inline std::string zhangProductCanonicalCoordinateKey(
	E_ObsCode observable, const std::string& satellite)
{
	return std::to_string(static_cast<int>(observable)) + "|" + satellite;
}

inline bool zhangBuildCanonicalProductExpansion(
	const ZhangExactVector& discoveryRow,
	const std::vector<std::string>& coordinateSatellites,
	const std::string& referenceSatellite,
	E_ObsCode firstObservable,
	E_ObsCode secondObservable,
	std::map<std::string, ZhangExactInteger>& expansion)
{
	expansion.clear();
	const int rank = static_cast<int>(coordinateSatellites.size());
	if (rank <= 0 || discoveryRow.size() != static_cast<std::size_t>(2 * rank) ||
		referenceSatellite.empty() || firstObservable == E_ObsCode::NONE ||
		secondObservable == E_ObsCode::NONE || firstObservable == secondObservable)
		return false;
	for (int signal = 0; signal < 2; signal++)
	{
		const E_ObsCode observable = signal == 0
			? firstObservable : secondObservable;
		ZhangExactInteger sum = 0;
		for (int coordinate = 0; coordinate < rank; coordinate++)
		{
			const auto coefficient = discoveryRow[signal * rank + coordinate];
			if (coefficient == 0) continue;
			if (coordinateSatellites[coordinate].empty() ||
				coordinateSatellites[coordinate] == referenceSatellite)
			{
				expansion.clear();
				return false;
			}
			expansion[zhangProductCanonicalCoordinateKey(
				observable, coordinateSatellites[coordinate])] += coefficient;
			sum += coefficient;
		}
		if (sum != 0)
			expansion[zhangProductCanonicalCoordinateKey(
				observable, referenceSatellite)] -= sum;
	}
	for (auto iterator = expansion.begin(); iterator != expansion.end();)
		if (iterator->second == 0) iterator = expansion.erase(iterator);
		else ++iterator;
	return !expansion.empty();
}

inline ZhangExactVector zhangDenseCanonicalProductExpansion(
	const std::map<std::string, ZhangExactInteger>& expansion,
	const std::map<std::string, int>& columns)
{
	ZhangExactVector dense(columns.size());
	for (const auto& [identity, coefficient] : expansion)
	{
		auto column = columns.find(identity);
		if (column == columns.end()) return {};
		dense[column->second] += coefficient;
	}
	return dense;
}

/** Canonicalise the sign of one exact product relation.  A row and its
 * simultaneous integer negation are the same lattice certificate and must not
 * occupy two ledger identities merely because a compiler chose the opposite
 * satellite/reference orientation.  The immutable canonical product image,
 * when present, is the sole sign authority; the backend physical chart is only
 * a transport view and may change independently. */
inline void zhangCanonicaliseProductLedgerRow(ProductIntegerLedgerRow& row)
{
	for (auto iterator = row.physicalExpansion.begin();
		 iterator != row.physicalExpansion.end();)
	{
		if (iterator->second == 0)
			iterator = row.physicalExpansion.erase(iterator);
		else ++iterator;
	}
	for (auto iterator = row.canonicalProductExpansion.begin();
		 iterator != row.canonicalProductExpansion.end();)
	{
		if (iterator->second == 0)
			iterator = row.canonicalProductExpansion.erase(iterator);
		else ++iterator;
	}
	const auto& signAuthority = row.physicalExpansionExact || row.canonicalProductExpansion.empty()
		? row.physicalExpansion : row.canonicalProductExpansion;
	if (signAuthority.empty() || signAuthority.begin()->second > 0) return;
	for (auto& [identity, coefficient] : row.physicalExpansion)
		coefficient = -coefficient;
	for (auto& [identity, coefficient] : row.canonicalProductExpansion)
		coefficient = -coefficient;
	for (auto& coefficient : row.productRow) coefficient = -coefficient;
	row.integerValue = -row.integerValue;
	if (row.pairCertificate && !row.firstSatellite.empty() &&
		!row.secondSatellite.empty())
		std::swap(row.firstSatellite, row.secondSatellite);
}

/** Immutable ledger identity.  Once an exact canonical satellite-product row
 * exists, backend generation and temporary physical-arc coordinates are views,
 * not identities.  The exact canonical row plus row-local phase segments form
 * the stable key; legacy rows without that certificate retain the old strict
 * generation/physical identity and therefore cannot be silently transported. */
inline std::string zhangProductLedgerIdentityFingerprint(
	const ProductIntegerLedgerRow& row)
{
    std::ostringstream stream;
    if(row.physicalExpansionExact) {
        stream << "PHYSICAL_CYCLE_V1|SYS" << static_cast<int>(row.system)
               << "|ROW{" << zhangProductPhysicalRowFingerprint(row.physicalExpansion) << "}";
        return stream.str();
    }
    if (!row.canonicalProductExpansion.empty())
	{
		stream << "CANONICAL|SEG{" << row.phaseSegmentFingerprint << "}|ROW{"
			<< zhangProductPhysicalRowFingerprint(
				row.canonicalProductExpansion) << "}";
	}
	else
	{
		stream << "LEGACY|GEN" << row.backendBasisGeneration
			<< "|SEG{" << row.phaseSegmentFingerprint << "}|ROW{"
			<< zhangProductPhysicalRowFingerprint(row.physicalExpansion) << "}";
	}
	return stream.str();
}

/** Re-express a retained physical integer row in a current ambiguity-column
 * map.  Backend generation is deliberately not part of this algebraic
 * operation: safety comes from exact arc/version identities, with phase
 * segment and statistical gates applied by the caller. */
inline bool zhangProjectProductLedgerPhysicalRow(
	const ProductIntegerLedgerRow& row,
	const std::map<std::string, int>& currentIdentityColumns,
	int columnCount,
	ZhangExactVector& projected)
{
	projected = ZhangExactVector(std::max(0, columnCount));
	if (columnCount <= 0 || row.physicalExpansion.empty()) return false;
	for (const auto& [identity, coefficient] : row.physicalExpansion)
	{
		if (coefficient == 0) continue;
		auto column = currentIdentityColumns.find(identity);
		if (column == currentIdentityColumns.end() || column->second < 0 ||
			column->second >= columnCount) return false;
		projected[column->second] += coefficient;
	}
	return std::any_of(projected.begin(), projected.end(),
		[](const auto& coefficient) { return coefficient != 0; });
}

inline int zhangProductLedgerExactRank(
	const std::vector<ProductIntegerLedgerRow>& rows,
	bool certifiedOnly = true)
{
	std::set<std::string> columnSet;
	for (const auto& row : rows)
	{
		if (certifiedOnly && !row.certified) continue;
		for (const auto& [identity, coefficient] : row.physicalExpansion)
			if (coefficient != 0) columnSet.insert(identity);
	}
	std::vector<std::string> columns(columnSet.begin(), columnSet.end());
	std::map<std::string, int> columnIndex;
	for (int column = 0; column < static_cast<int>(columns.size()); column++)
		columnIndex[columns[column]] = column;
	ZhangExactMatrix matrix;
	for (const auto& row : rows)
	{
		if (certifiedOnly && !row.certified) continue;
		ZhangExactVector dense(columns.size());
		for (const auto& [identity, coefficient] : row.physicalExpansion)
			if (coefficient != 0) dense[columnIndex.at(identity)] = coefficient;
		matrix.push_back(std::move(dense));
	}
	return static_cast<int>(
		zhangExactRowHermiteNormalForm(matrix).basis.size());
}

class ProductIntegerLedger
{
public:
	ZhangProductIntegerLedgerUpdate observe(
		long int epoch,
		const std::vector<ProductIntegerLedgerRow>& candidates,
		int requiredConfirmations = 1)
	{
		ZhangProductIntegerLedgerUpdate result;
		result.inputRows = candidates.size();
		result.activeRankBefore = zhangProductLedgerExactRank(rows_);
		if (epoch <= 0 || requiredConfirmations < 1)
		{
			result.failureReason = "PRODUCT_LEDGER_INPUT_INVALID";
			return result;
		}
		if (std::any_of(rows_.begin(), rows_.end(), [epoch](const auto& row)
			{
				return epoch < row.lastConfirmed;
			}))
		{
			// Ledger time is global for one runtime/system registry.  Accepting an
			// older observation for any identity would make confirmation history
			// depend on callback order.  Reject the complete transaction before a
			// proposed copy can be committed.
			result.failureReason = "PRODUCT_LEDGER_EPOCH_REGRESSION";
			return result;
		}
		auto proposedRows = rows_;
        for (std::size_t candidateIndex = 0; candidateIndex < candidates.size(); ++candidateIndex)
        {
            auto candidate = candidates[candidateIndex];
			zhangCanonicaliseProductLedgerRow(candidate);
			if (candidate.system == E_Sys::NONE ||
				candidate.physicalExpansion.empty() ||
				candidate.canonicalProductExpansion.empty() ||
				candidate.phaseSegmentFingerprint.empty())
			{
				result.failureReason = "PRODUCT_LEDGER_IDENTITY_INCOMPLETE";
				return result;
			}
			if (!std::isfinite(candidate.admissionFailureProbabilityBound) ||
				candidate.admissionFailureProbabilityBound < 0 ||
				candidate.admissionFailureProbabilityBound > 1)
			{
				result.failureReason =
					"PRODUCT_LEDGER_ADMISSION_FAILURE_BOUND_INVALID";
				return result;
			}
			auto exactGcd = [](ZhangExactInteger left, ZhangExactInteger right)
			{
				left = left < 0 ? -left : left;
				right = right < 0 ? -right : right;
				while (right != 0)
				{
					const ZhangExactInteger remainder = left % right;
					left = right;
					right = remainder;
				}
				return left;
			};
			ZhangExactInteger coefficientGcd = 0;
			for (const auto& [identity, coefficient] : candidate.physicalExpansion)
				coefficientGcd = exactGcd(
					coefficientGcd, coefficient < 0 ? -coefficient : coefficient);
			if (coefficientGcd != 1)
			{
				result.failureReason = "PRODUCT_LEDGER_ROW_NOT_PRIMITIVE";
				return result;
			}
			const std::string identity =
				zhangProductLedgerIdentityFingerprint(candidate);
			auto existing = std::find_if(proposedRows.begin(), proposedRows.end(),
				[&](const auto& row)
				{
					return zhangProductLedgerIdentityFingerprint(row) == identity;
				});
			if (existing == proposedRows.end())
			{
				candidate.firstCertified = epoch;
				candidate.lastConfirmed = epoch;
				candidate.confirmationEpochs = 1;
				candidate.certified = requiredConfirmations == 1;
				proposedRows.push_back(std::move(candidate));
				result.freshRows++;
				result.confirmedRows += proposedRows.back().certified;
				continue;
			}
			if (existing->integerValue != candidate.integerValue)
			{
				// The immutable identity already contains the phase-segment
				// fingerprint.  A different integer observed for an established
				// certificate in that same identity is therefore a competing
				// hypothesis, not a temporal continuation.  Preserve the confirmed
				// certificate and fail closed for this observation; replacing it
				// here would let one fresh epoch erase both its value and its
				// accumulated confirmation history before posterior arbitration.
				// Unconfirmed rows retain the legacy restart behaviour below.
                result.conflictingRows++;
                result.conflictingCandidateIndex = static_cast<int>(candidateIndex);
                result.incumbentInteger = existing->integerValue;
                result.proposedInteger = candidate.integerValue;
                result.conflictPhysicalRowsEqual = existing->physicalExpansion == candidate.physicalExpansion;
                result.incumbentBackendGeneration = existing->backendBasisGeneration;
                result.proposedBackendGeneration = candidate.backendBasisGeneration;
                result.conflictIdentity = identity;
                result.incumbentPhysicalRow = zhangProductPhysicalRowFingerprint(existing->physicalExpansion);
                result.proposedPhysicalRow = zhangProductPhysicalRowFingerprint(candidate.physicalExpansion);
				if (existing->certified)
				{
					// A mature conflict rejects the entire transaction.  proposedRows
					// is a private copy, so neither earlier non-conflicting candidates
					// in this batch nor the competing value can reach rows_.
					result.activeRankAfter = result.activeRankBefore;
					result.failureReason =
						"PRODUCT_LEDGER_MATURE_INTEGER_CONFLICT";
					return result;
				}
				candidate.firstCertified = epoch;
				candidate.lastConfirmed = epoch;
				candidate.confirmationEpochs = 1;
				candidate.certified = requiredConfirmations == 1;
				*existing = std::move(candidate);
				result.confirmedRows += existing->certified;
				continue;
			}
            if (!candidate.physicalExpansionExact && existing->canonicalProductExpansion !=
                candidate.canonicalProductExpansion)
			{
				result.failureReason =
					"PRODUCT_LEDGER_CANONICAL_IDENTITY_CONFLICT";
				return result;
			}
			if (existing->lastConfirmed != epoch)
			{
				existing->lastConfirmed = epoch;
				existing->confirmationEpochs++;
			}
            if(candidate.physicalExpansionExact) {
                existing->canonicalProductExpansion=candidate.canonicalProductExpansion;
                existing->phaseSegmentFingerprint=candidate.phaseSegmentFingerprint;
            }
            existing->productRow = std::move(candidate.productRow);
			existing->physicalExpansion = std::move(candidate.physicalExpansion);
			existing->backendBasisGeneration = candidate.backendBasisGeneration;
			// Exact pair membership is a stronger semantic certificate than a
			// generic conditioning row.  Upgrade the retained physical integer if
			// it was first seen as a mixed PRODUCT_ILS row.
			if (candidate.pairCertificate)
			{
				existing->pairCertificate = true;
				existing->conditioningOnly = false;
				existing->source = candidate.source;
				existing->coordinate = candidate.coordinate;
				existing->firstSatellite = candidate.firstSatellite;
				existing->secondSatellite = candidate.secondSatellite;
			}
			existing->certified =
				existing->confirmationEpochs >= requiredConfirmations;
			result.confirmedRows += existing->certified;
		}
        if(std::any_of(candidates.begin(),candidates.end(),[](const auto& row){return row.physicalExpansionExact;})) {
            std::map<std::string,int> columns;
            for(const auto& row:proposedRows) if(row.physicalExpansionExact)
                for(const auto& [id,value]:row.physicalExpansion)
                    if(value!=0 && !columns.contains(id)) columns[id]=static_cast<int>(columns.size());
            ZhangExactMatrix physical; ZhangExactVector values;
            for(const auto& row:proposedRows) if(row.physicalExpansionExact) {
                ZhangExactVector dense(columns.size());
                for(const auto& [id,value]:row.physicalExpansion) if(value!=0) dense[columns.at(id)]=value;
                physical.push_back(std::move(dense));values.push_back(row.integerValue);
            }
            const auto joint=zhangExactRowHermiteNormalForm(physical,values);
            // F*N=n has an integer solution iff n belongs to the column
            // lattice of F. This also rejects parity/divisibility conflicts
            // that rational affine consistency alone cannot detect.
            ZhangExactMatrix columnLattice(columns.size(),ZhangExactVector(physical.size()));
            for(int r=0;r<physical.size();++r) for(int c=0;c<columns.size();++c)
                columnLattice[c][r]=physical[r][c];
            const bool integerFeasible=values.empty() ||
                zhangIntegerRowLatticeContains(columnLattice,values).contained;
            if(!joint.consistent || !integerFeasible) {
                result.activeRankAfter=result.activeRankBefore;
                result.conflictingRows++;
                result.failureReason="PRODUCT_LEDGER_TRUE_PHYSICAL_AFFINE_CONFLICT";
                if(zhangR51Enabled()) {
                    result.exactConflict=zhangR51ExactConflictWitness(physical,values);
                    for(const auto& row:proposedRows) if(row.physicalExpansionExact)
                        result.conflictRows.push_back(zhangProductPhysicalRowFingerprint(row.physicalExpansion));
                }
                return result;
            }
        }
		result.activeRankAfter = zhangProductLedgerExactRank(proposedRows);
		rows_ = std::move(proposedRows);
		result.valid = true;
		result.failureReason = "NONE";
		return result;
	}

    struct Preflight {
        const ZhangProductIntegerLedgerUpdate update;
        const std::string root,physicalEpoch,baseIdentity,proposedIdentity;
        const long int epoch;
        const int requiredConfirmations;
        const std::vector<ProductIntegerLedgerRow> proposed;
    private:
        friend class ProductIntegerLedger;
        Preflight(ZhangProductIntegerLedgerUpdate result,std::string rootId,std::string physicalId,
            std::string baseId,std::string proposalId,long int time,int confirmations,
            std::vector<ProductIntegerLedgerRow> rows)
            :update(std::move(result)),root(std::move(rootId)),physicalEpoch(std::move(physicalId)),
             baseIdentity(std::move(baseId)),proposedIdentity(std::move(proposalId)),epoch(time),
             requiredConfirmations(confirmations),proposed(std::move(rows)){}
    };
    static std::string snapshotIdentity(const std::vector<ProductIntegerLedgerRow>& rows) {
        std::ostringstream out;
        for(const auto& row:rows) {
            out<<zhangProductLedgerIdentityFingerprint(row)<<"="<<row.integerValue
               <<":"<<row.lastConfirmed<<":"<<row.confirmationEpochs<<":"<<row.certified
               <<":"<<row.backendBasisGeneration<<":"<<row.phaseSegmentFingerprint
               <<":"<<zhangProductPhysicalRowFingerprint(row.physicalExpansion)<<"|";
            out<<std::setprecision(17)<<":"<<int(row.system)<<":"<<int(row.firstObservable)<<":"<<int(row.secondObservable)
               <<":"<<row.physicalExpansionExact<<":"<<row.coordinate<<":"<<row.firstSatellite<<":"<<row.secondSatellite
               <<":"<<row.admissionFailureProbabilityBound<<":"<<row.firstCertified<<":"<<int(row.source)
               <<":"<<row.conditioningOnly<<":"<<row.pairCertificate<<":";
            for(const auto& value:row.productRow)out<<value<<",";
            out<<":"<<zhangProductPhysicalRowFingerprint(row.canonicalProductExpansion)<<":";
            const auto closure=zhangDecisionRiskClosure(row.decisionProofs);
            out<<closure.valid<<":"<<closure.reason<<":";
            for(const auto& [id,proof]:closure.atoms) {
                out<<std::quoted(id)<<std::quoted(proof->originalStatement)<<std::quoted(proof->observationProvenance)
                   <<proof->conditionalFailureBound;
                for(const auto& parent:proof->parents)out<<std::quoted(parent?parent->id:"NULL");
            }
        }
        return out.str();
    }
    Preflight preflight(long int epoch,const std::vector<ProductIntegerLedgerRow>& candidates,
        int confirmations,const std::string& root,const std::string& physicalEpoch) const {
        ProductIntegerLedger trial=*this;
        auto result=trial.observe(epoch,candidates,confirmations);
        if(root.empty() || physicalEpoch.empty()) {
            result.valid=false;result.failureReason="R51_PREFLIGHT_IDENTITY_EMPTY";
        }
        const auto base=snapshotIdentity(rows_),proposal=snapshotIdentity(trial.rows_);
        return Preflight(std::move(result),root,physicalEpoch,base,proposal,epoch,confirmations,std::move(trial.rows_));
    }
    ZhangProductIntegerLedgerUpdate commit(const Preflight& receipt,
        const std::string& root,const std::string& physicalEpoch) {
        auto result=receipt.update;
        if(!zhangProductLedgerWriterCommitAuthorized(result)) return result;
        if(root!=receipt.root || physicalEpoch!=receipt.physicalEpoch ||
           snapshotIdentity(rows_)!=receipt.baseIdentity ||
           snapshotIdentity(receipt.proposed)!=receipt.proposedIdentity) {
            result.valid=false;result.failureReason="R51_PREFLIGHT_RECEIPT_STALE";return result;
        }
        rows_=receipt.proposed;return result;
    }
	const std::vector<ProductIntegerLedgerRow>& rows() const { return rows_; }

	std::vector<ProductIntegerLedgerRow> rowsForGeneration(
		std::uint64_t backendBasisGeneration,
		bool certifiedOnly = true) const
	{
		std::vector<ProductIntegerLedgerRow> selected;
		for (const auto& row : rows_)
		{
			if (row.backendBasisGeneration != backendBasisGeneration) continue;
			if (certifiedOnly && !row.certified) continue;
			selected.push_back(row);
		}
		return selected;
	}

private:
	std::vector<ProductIntegerLedgerRow> rows_;
};

inline std::map<std::pair<std::string, E_Sys>, ProductIntegerLedger>&
zhangProductIntegerLedgerRegistry()
{
	static std::map<std::pair<std::string, E_Sys>, ProductIntegerLedger> registry;
	return registry;
}
