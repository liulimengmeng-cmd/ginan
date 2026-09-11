#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "common/zhangIntegerAudit.hpp"

/** Exact target/held quotient audit in a shared physical-ambiguity ambient
 * coordinate.
 *
 * targetRows must be an integer basis of the product target lattice. heldRows
 * are already certified integer facts in the same physical-arc coordinate.
 * No covariance or small-variance decision is accepted as integer evidence.
 */
struct ZhangHeldQuotientAudit
{
	ZhangExactMatrix targetBasis;
	ZhangExactMatrix heldIntersectionPhysicalBasis;
	ZhangExactMatrix heldIntersectionTargetCoordinates;
	ZhangExactVector heldIntersectionValues;
	ZhangExactMatrix quotientTargetCoordinates;

	int targetRank = 0;
	int heldIntersectionRank = 0;
	int quotientRank = 0;
	bool heldIntersectionPrimitiveInTarget = false;
	bool exactClosure = false;
	bool valid = false;
	std::string failureReason;
};

struct ZhangCertifiedUnionAudit
{
	ZhangExactMatrix certifiedBasis;
	ZhangExactVector certifiedValues;
	int targetRank = 0;
	int heldRank = 0;
	int newlyFixedRank = 0;
	int combinedCertifiedRank = 0;
	bool targetContainedInCertified = false;
	bool certifiedContainedInTarget = false;
	bool exactTargetEquality = false;
	bool consistent = false;
	std::string failureReason;
};

struct ZhangDeterministicQuotientAudit
{
	int covarianceRank = 0;
	int nullity = 0;
	double maximumNullFractionalInteger = 0;
	bool covarianceValid = false;
	bool integerConsistent = true;
	std::string status;
};

/** Maximum exact dual-frequency satellite-product quotient lattice.
 *
 * Network rows A and product rows P1/P2 share one physical-ambiguity ambient
 * coordinate.  The returned q rows are the complete integer projection of
 * the saturated kernel of
 *
 *   [ A^T  0   -P1^T ] [u1]   [0]
 *   [ 0    A^T -(P1-P2)^T ] [uw] = [0].
 *
 * No named-pair order, reference-tree choice, covariance, or rounded float
 * candidate participates in this calculation.  Integer values are recovered
 * independently by exact membership in A for every canonical q row. */
struct ZhangDualProductIntersectionAudit
{
	ZhangExactMatrix productRows;
	ZhangExactVector firstSignalIntegers;
	ZhangExactVector wideLaneIntegers;
	std::vector<ZhangExactInteger> productSmithInvariants;
	int networkRank = 0;
	int fullProductRank = 0;
	int firstSignalIntersectionRank = 0;
	int wideLaneIntersectionRank = 0;
	int dualIntersectionRank = 0;
	bool exactClosure = false;
	bool productLatticePrimitive = false;
	bool valid = false;
	std::string failureReason;
};

/** Exact intersection of two affine row lattices in one product coordinate.
 *
 * The coefficient lattices may use unrelated primitive bases (for example a
 * raw fixed-lag WL solve and a conditional L1 solve).  The saturated integer
 * kernel of [A^T,-B^T] is projected back to the common physical row q.  Both
 * integer right-hand sides are carried through independent exact HNF
 * reductions, so a hidden zero-row affine contradiction is rejected instead
 * of being turned into a dual-frequency certificate. */
struct ZhangExactAffineLatticeIntersection
{
	ZhangExactMatrix rows;
	ZhangExactVector firstValues;
	ZhangExactVector secondValues;
	std::vector<ZhangExactInteger> smithInvariants;
	int firstRank = 0;
	int secondRank = 0;
	int intersectionRank = 0;
	bool primitive = false;
	bool exactClosure = false;
	bool valid = false;
	std::string failureReason = "NOT_EVALUATED";
};

/** Common integer functionals on which two affine row lattices also agree on
 * the integer right-hand side.
 *
 * zhangExactAffineRowLatticeIntersection() deliberately carries the two RHS
 * vectors independently.  This second step computes the exact integer kernel
 * of their difference.  It therefore keeps combinations whose individual
 * basis rows disagree but whose combined functional agrees; comparing named
 * pair scalars before this operation is mathematically invalid.
 */
struct ZhangExactAffineAgreementLattice
{
	ZhangExactMatrix rows;
	ZhangExactVector values;
	int commonRank = 0;
	int agreementRank = 0;
	int conflictRank = 0;
	bool exactClosure = false;
	bool valid = false;
	std::string failureReason = "NOT_EVALUATED";
};

/** Inclusion-minimal candidate-side conflict certificate relative to an
 * immutable baseline affine lattice.  The result identifies only candidate
 * rows that are jointly needed for one contradiction; callers can quarantine
 * that certificate family/component and retain every other consistent row.
 */
struct ZhangMinimalAffineConflictCertificate
{
	std::vector<int> candidateIndices;
	bool baselineConsistent = false;
	bool combinedConsistent = false;
	bool minimal = false;
	bool valid = false;
	std::string failureReason = "NOT_EVALUATED";
};

inline bool zhangExactRectangularMatrix(
	const ZhangExactMatrix& rows,
	std::size_t dimension);

inline ZhangExactVector zhangExactRowCombination(
	const ZhangExactVector& coefficients,
	const ZhangExactMatrix& rows);

/** Embed exact rows from a currently representable canonical subset into the
 * complete canonical catalogue.  This is deliberately exact and rejects
 * duplicate/out-of-range coordinates; an unavailable direction is represented
 * by a zero column in the embedded q, never by a zero physical target row. */
inline bool zhangEmbedCanonicalSubsetRows(
	ZhangExactMatrix& rows,
	const std::vector<int>& canonicalIndices,
	std::size_t canonicalDimension)
{
	std::vector<bool> seen(canonicalDimension, false);
	for (const int index : canonicalIndices)
	{
		if (index < 0 || static_cast<std::size_t>(index) >= canonicalDimension ||
			seen[index]) return false;
		seen[index] = true;
	}
	for (auto& row : rows)
	{
		if (row.size() != canonicalIndices.size()) return false;
		ZhangExactVector embedded(canonicalDimension);
		for (std::size_t local = 0; local < canonicalIndices.size(); local++)
			embedded[canonicalIndices[local]] = row[local];
		row = std::move(embedded);
	}
	return true;
}

inline ZhangDualProductIntersectionAudit zhangExactDualProductIntersection(
	const ZhangExactMatrix& networkRows,
	const ZhangExactVector& networkValues,
	const ZhangExactMatrix& firstProductRows,
	const ZhangExactMatrix& secondProductRows)
{
	ZhangDualProductIntersectionAudit result;
	if (networkRows.empty() || firstProductRows.empty() ||
		firstProductRows.size() != secondProductRows.size() ||
		networkValues.size() != networkRows.size())
	{
		result.failureReason = "DUAL_PRODUCT_INTERSECTION_INPUT_EMPTY_OR_MISMATCHED";
		return result;
	}
	const std::size_t ambient = networkRows.front().size();
	const std::size_t productDimension = firstProductRows.size();
	if (!zhangExactRectangularMatrix(networkRows, ambient) ||
		!zhangExactRectangularMatrix(firstProductRows, ambient) ||
		!zhangExactRectangularMatrix(secondProductRows, ambient))
	{
		result.failureReason = "DUAL_PRODUCT_INTERSECTION_AMBIENT_MISMATCH";
		return result;
	}
	const auto networkHnf = zhangExactRowHermiteNormalForm(
		networkRows, networkValues);
	if (!networkHnf.consistent)
	{
		result.failureReason = "DUAL_PRODUCT_NETWORK_AFFINE_INCONSISTENCY";
		return result;
	}
	result.networkRank = static_cast<int>(networkHnf.basis.size());
	result.fullProductRank = static_cast<int>(
		zhangExactRowHermiteNormalForm(firstProductRows).basis.size());
	ZhangExactMatrix wideLaneRows = firstProductRows;
	for (std::size_t row = 0; row < productDimension; row++)
	for (std::size_t column = 0; column < ambient; column++)
		wideLaneRows[row][column] -= secondProductRows[row][column];
	const std::size_t networkRank = networkHnf.basis.size();
	auto singleIntersectionRank = [&](const ZhangExactMatrix& productRows)
	{
		ZhangExactMatrix singleRelation(ambient,
			ZhangExactVector(networkRank + productDimension));
		for (std::size_t column = 0; column < ambient; column++)
		{
			for (std::size_t row = 0; row < networkRank; row++)
				singleRelation[column][row] = networkHnf.basis[row][column];
			for (std::size_t row = 0; row < productDimension; row++)
				singleRelation[column][networkRank + row] =
					-productRows[row][column];
		}
		ZhangExactMatrix coordinates;
		for (const auto& kernelRow : zhangExactIntegerKernel(singleRelation))
		{
			ZhangExactVector q(kernelRow.begin() + networkRank, kernelRow.end());
			if (std::any_of(q.begin(), q.end(),
				[](const auto& value) { return value != 0; }))
				coordinates.push_back(std::move(q));
		}
		return static_cast<int>(
			zhangExactRowHermiteNormalForm(coordinates).basis.size());
	};
	result.firstSignalIntersectionRank =
		singleIntersectionRank(firstProductRows);
	result.wideLaneIntersectionRank = singleIntersectionRank(wideLaneRows);

	ZhangExactMatrix relation(2 * ambient,
		ZhangExactVector(2 * networkRank + productDimension));
	for (std::size_t column = 0; column < ambient; column++)
	{
		for (std::size_t row = 0; row < networkRank; row++)
		{
			relation[column][row] = networkHnf.basis[row][column];
			relation[ambient + column][networkRank + row] =
				networkHnf.basis[row][column];
		}
		for (std::size_t row = 0; row < productDimension; row++)
		{
			relation[column][2 * networkRank + row] =
				-firstProductRows[row][column];
			relation[ambient + column][2 * networkRank + row] =
				-wideLaneRows[row][column];
		}
	}
	ZhangExactMatrix projected;
	for (const auto& kernelRow : zhangExactIntegerKernel(relation))
	{
		ZhangExactVector q(kernelRow.begin() + 2 * networkRank,
			kernelRow.end());
		if (std::any_of(q.begin(), q.end(),
			[](const auto& value) { return value != 0; }))
			projected.push_back(std::move(q));
	}
	const auto projectedHnf = zhangExactRowHermiteNormalForm(projected);
	if (!projectedHnf.consistent)
	{
		result.failureReason = "DUAL_PRODUCT_PROJECTED_HNF_FAILED";
		return result;
	}
	result.productRows = projectedHnf.basis;
	result.dualIntersectionRank = static_cast<int>(result.productRows.size());
	const auto productSmith = zhangIntegerRowLatticeContains(
		result.productRows, ZhangExactVector(productDimension));
	result.productSmithInvariants = productSmith.smithInvariants;
	result.productLatticePrimitive =
		static_cast<int>(result.productSmithInvariants.size()) ==
			result.dualIntersectionRank &&
		std::all_of(result.productSmithInvariants.begin(),
			result.productSmithInvariants.end(), [](const auto& invariant)
			{
				return zhangExactAbs(invariant) == 1;
			});
	for (const auto& q : result.productRows)
	{
		const auto firstPhysical = zhangExactRowCombination(q, firstProductRows);
		const auto widePhysical = zhangExactRowCombination(q, wideLaneRows);
		const auto firstMembership = zhangIntegerRowLatticeContains(
			networkHnf.basis, firstPhysical);
		const auto wideMembership = zhangIntegerRowLatticeContains(
			networkHnf.basis, widePhysical);
		if (!firstMembership.contained || !wideMembership.contained)
		{
			result.failureReason = "DUAL_PRODUCT_PROJECTED_ROW_NOT_IN_NETWORK";
			return result;
		}
		ZhangExactInteger firstValue = 0;
		ZhangExactInteger wideValue = 0;
		for (std::size_t row = 0; row < networkHnf.values.size(); row++)
		{
			firstValue += firstMembership.combination[row] * networkHnf.values[row];
			wideValue += wideMembership.combination[row] * networkHnf.values[row];
		}
		result.firstSignalIntegers.push_back(firstValue);
		result.wideLaneIntegers.push_back(wideValue);
	}
	result.exactClosure = true;
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangExactAffineLatticeIntersection
zhangExactAffineRowLatticeIntersection(
	const ZhangExactMatrix& firstRows,
	const ZhangExactVector& firstValues,
	const ZhangExactMatrix& secondRows,
	const ZhangExactVector& secondValues)
{
	ZhangExactAffineLatticeIntersection result;
	if (firstRows.empty() || secondRows.empty() ||
		firstRows.size() != firstValues.size() ||
		secondRows.size() != secondValues.size())
	{
		result.failureReason =
			"AFFINE_LATTICE_INTERSECTION_INPUT_EMPTY_OR_MISMATCHED";
		return result;
	}
	const std::size_t dimension = firstRows.front().size();
	if (dimension == 0 ||
		!zhangExactRectangularMatrix(firstRows, dimension) ||
		!zhangExactRectangularMatrix(secondRows, dimension))
	{
		result.failureReason =
			"AFFINE_LATTICE_INTERSECTION_DIMENSION_MISMATCH";
		return result;
	}
	const ZhangExactRowHnf firstHnf =
		zhangExactRowHermiteNormalForm(firstRows, firstValues);
	const ZhangExactRowHnf secondHnf =
		zhangExactRowHermiteNormalForm(secondRows, secondValues);
	if (!firstHnf.consistent || !secondHnf.consistent)
	{
		result.failureReason =
			"AFFINE_LATTICE_INTERSECTION_INPUT_AFFINE_CONFLICT";
		return result;
	}
	result.firstRank = static_cast<int>(firstHnf.basis.size());
	result.secondRank = static_cast<int>(secondHnf.basis.size());
	if (result.firstRank == 0 || result.secondRank == 0)
	{
		result.valid = true;
		result.primitive = true;
		result.exactClosure = true;
		result.failureReason = "NONE";
		return result;
	}

	ZhangExactMatrix relation(
		dimension,
		ZhangExactVector(firstHnf.basis.size() + secondHnf.basis.size()));
	for (std::size_t column = 0; column < dimension; column++)
	{
		for (std::size_t row = 0; row < firstHnf.basis.size(); row++)
			relation[column][row] = firstHnf.basis[row][column];
		for (std::size_t row = 0; row < secondHnf.basis.size(); row++)
			relation[column][firstHnf.basis.size() + row] =
				-secondHnf.basis[row][column];
	}
	ZhangExactMatrix projectedRows;
	ZhangExactVector projectedFirstValues;
	ZhangExactVector projectedSecondValues;
	for (const auto& kernelRow : zhangExactIntegerKernel(relation))
	{
		ZhangExactVector firstCoefficients(
			kernelRow.begin(),
			kernelRow.begin() + firstHnf.basis.size());
		ZhangExactVector secondCoefficients(
			kernelRow.begin() + firstHnf.basis.size(), kernelRow.end());
		ZhangExactVector common = zhangExactRowCombination(
			firstCoefficients, firstHnf.basis);
		if (common.empty() || std::all_of(common.begin(), common.end(),
			[](const auto& value) { return value == 0; }))
			continue;
		projectedRows.push_back(std::move(common));
		ZhangExactInteger firstValue = 0;
		ZhangExactInteger secondValue = 0;
		for (std::size_t row = 0; row < firstCoefficients.size(); row++)
			firstValue += firstCoefficients[row] * firstHnf.values[row];
		for (std::size_t row = 0; row < secondCoefficients.size(); row++)
			secondValue += secondCoefficients[row] * secondHnf.values[row];
		projectedFirstValues.push_back(std::move(firstValue));
		projectedSecondValues.push_back(std::move(secondValue));
	}
	if (projectedRows.empty())
	{
		result.valid = true;
		result.primitive = true;
		result.exactClosure = true;
		result.failureReason = "NONE";
		return result;
	}

	const ZhangExactRowHnf firstProjection =
		zhangExactRowHermiteNormalForm(projectedRows, projectedFirstValues);
	const ZhangExactRowHnf secondProjection =
		zhangExactRowHermiteNormalForm(projectedRows, projectedSecondValues);
	if (!firstProjection.consistent || !secondProjection.consistent)
	{
		result.failureReason =
			"AFFINE_LATTICE_INTERSECTION_PROJECTED_AFFINE_CONFLICT";
		return result;
	}
	if (firstProjection.basis != secondProjection.basis)
	{
		result.failureReason =
			"AFFINE_LATTICE_INTERSECTION_HNF_DISAGREEMENT";
		return result;
	}
	result.rows = firstProjection.basis;
	result.firstValues = firstProjection.values;
	result.secondValues = secondProjection.values;
	result.intersectionRank = static_cast<int>(result.rows.size());
	if (result.rows.size() != result.firstValues.size() ||
		result.rows.size() != result.secondValues.size())
	{
		result.failureReason =
			"AFFINE_LATTICE_INTERSECTION_OUTPUT_DIMENSION_MISMATCH";
		return result;
	}
	for (std::size_t row = 0; row < result.rows.size(); row++)
	{
		const auto firstMembership = zhangIntegerRowLatticeContains(
			firstHnf.basis, result.rows[row]);
		const auto secondMembership = zhangIntegerRowLatticeContains(
			secondHnf.basis, result.rows[row]);
		if (!firstMembership.contained || !secondMembership.contained)
		{
			result.failureReason =
				"AFFINE_LATTICE_INTERSECTION_MEMBERSHIP_FAILURE";
			return result;
		}
		ZhangExactInteger recoveredFirst = 0;
		ZhangExactInteger recoveredSecond = 0;
		for (std::size_t parent = 0;
			 parent < firstMembership.combination.size(); parent++)
			recoveredFirst += firstMembership.combination[parent] *
				firstHnf.values[parent];
		for (std::size_t parent = 0;
			 parent < secondMembership.combination.size(); parent++)
			recoveredSecond += secondMembership.combination[parent] *
				secondHnf.values[parent];
		if (recoveredFirst != result.firstValues[row] ||
			recoveredSecond != result.secondValues[row])
		{
			result.failureReason =
				"AFFINE_LATTICE_INTERSECTION_VALUE_RECOVERY_FAILURE";
			return result;
		}
	}
	const auto smith = zhangIntegerRowLatticeContains(
		result.rows, ZhangExactVector(dimension));
	result.smithInvariants = smith.smithInvariants;
	result.primitive = static_cast<int>(result.smithInvariants.size()) ==
		result.intersectionRank &&
		std::all_of(result.smithInvariants.begin(),
			result.smithInvariants.end(), [](const auto& invariant)
			{
				return zhangExactAbs(invariant) == 1;
			});
	result.exactClosure = true;
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangExactAffineAgreementLattice
zhangExactAffineAgreementLattice(
	const ZhangExactMatrix& firstRows,
	const ZhangExactVector& firstValues,
	const ZhangExactMatrix& secondRows,
	const ZhangExactVector& secondValues)
{
	ZhangExactAffineAgreementLattice result;
	const auto common = zhangExactAffineRowLatticeIntersection(
		firstRows, firstValues, secondRows, secondValues);
	if (!common.valid)
	{
		result.failureReason = common.failureReason;
		return result;
	}
	result.commonRank = common.intersectionRank;
	if (common.rows.empty())
	{
		result.valid = true;
		result.exactClosure = true;
		result.failureReason = "NONE";
		return result;
	}
	if (common.rows.size() != common.firstValues.size() ||
		common.rows.size() != common.secondValues.size())
	{
		result.failureReason = "AFFINE_AGREEMENT_COMMON_DIMENSION_MISMATCH";
		return result;
	}

	ZhangExactVector delta(common.rows.size());
	for (std::size_t row = 0; row < common.rows.size(); row++)
		delta[row] = common.firstValues[row] - common.secondValues[row];
	if (std::all_of(delta.begin(), delta.end(),
		[](const auto& value) { return value == 0; }))
	{
		result.rows = common.rows;
		result.values = common.firstValues;
		result.agreementRank = common.intersectionRank;
		result.conflictRank = 0;
		result.valid = true;
		result.exactClosure = true;
		result.failureReason = "NONE";
		return result;
	}

	// delta * k^T = 0.  The one equation is a matrix row and its columns are
	// the coefficients of the common-lattice basis.
	const ZhangExactMatrix deltaEquation{delta};
	const auto agreementCoefficients = zhangExactIntegerKernel(deltaEquation);
	ZhangExactMatrix agreementRows;
	ZhangExactVector agreementValues;
	for (const auto& coefficients : agreementCoefficients)
	{
		if (coefficients.size() != common.rows.size())
		{
			result.failureReason =
				"AFFINE_AGREEMENT_KERNEL_DIMENSION_MISMATCH";
			return result;
		}
		auto row = zhangExactRowCombination(coefficients, common.rows);
		if (row.empty() || std::all_of(row.begin(), row.end(),
			[](const auto& value) { return value == 0; }))
			continue;
		ZhangExactInteger value = 0;
		for (std::size_t parent = 0; parent < coefficients.size(); parent++)
			value += coefficients[parent] * common.firstValues[parent];
		agreementRows.push_back(std::move(row));
		agreementValues.push_back(std::move(value));
	}
	if (!agreementRows.empty())
	{
		const auto hnf = zhangExactRowHermiteNormalForm(
			agreementRows, agreementValues);
		if (!hnf.consistent || hnf.basis.size() != hnf.values.size())
		{
			result.failureReason = "AFFINE_AGREEMENT_HNF_CONFLICT";
			return result;
		}
		result.rows = hnf.basis;
		result.values = hnf.values;
	}
	result.agreementRank = static_cast<int>(result.rows.size());
	result.conflictRank = result.commonRank - result.agreementRank;
	for (std::size_t row = 0; row < result.rows.size(); row++)
	{
		const auto firstMembership = zhangIntegerRowLatticeContains(
			firstRows, result.rows[row]);
		const auto secondMembership = zhangIntegerRowLatticeContains(
			secondRows, result.rows[row]);
		if (!firstMembership.contained || !secondMembership.contained)
		{
			result.failureReason = "AFFINE_AGREEMENT_MEMBERSHIP_FAILURE";
			return result;
		}
		ZhangExactInteger firstValue = 0;
		ZhangExactInteger secondValue = 0;
		for (std::size_t parent = 0;
			 parent < firstMembership.combination.size(); parent++)
			firstValue += firstMembership.combination[parent] *
				firstValues[parent];
		for (std::size_t parent = 0;
			 parent < secondMembership.combination.size(); parent++)
			secondValue += secondMembership.combination[parent] *
				secondValues[parent];
		if (firstValue != secondValue || firstValue != result.values[row])
		{
			result.failureReason = "AFFINE_AGREEMENT_VALUE_RECOVERY_FAILURE";
			return result;
		}
	}
	result.valid = true;
	result.exactClosure = true;
	result.failureReason = "NONE";
	return result;
}

inline bool zhangExactAffineSystemConsistent(
	const ZhangExactMatrix& rows,
	const ZhangExactVector& values)
{
	if (rows.size() != values.size()) return false;
	if (rows.empty()) return true;
	return zhangExactRowHermiteNormalForm(rows, values).consistent;
}

inline ZhangMinimalAffineConflictCertificate
zhangMinimalAffineConflictCertificate(
	const ZhangExactMatrix& baselineRows,
	const ZhangExactVector& baselineValues,
	const ZhangExactMatrix& candidateRows,
	const ZhangExactVector& candidateValues)
{
	ZhangMinimalAffineConflictCertificate result;
	if (baselineRows.size() != baselineValues.size() ||
		candidateRows.size() != candidateValues.size())
	{
		result.failureReason = "AFFINE_CONFLICT_INPUT_DIMENSION_MISMATCH";
		return result;
	}
	result.baselineConsistent = zhangExactAffineSystemConsistent(
		baselineRows, baselineValues);
	if (!result.baselineConsistent)
	{
		result.failureReason = "AFFINE_CONFLICT_BASELINE_INCONSISTENT";
		return result;
	}
	ZhangExactMatrix combinedRows = baselineRows;
	ZhangExactVector combinedValues = baselineValues;
	combinedRows.insert(combinedRows.end(), candidateRows.begin(), candidateRows.end());
	combinedValues.insert(
		combinedValues.end(), candidateValues.begin(), candidateValues.end());
	result.combinedConsistent = zhangExactAffineSystemConsistent(
		combinedRows, combinedValues);
	result.valid = true;
	if (result.combinedConsistent)
	{
		result.minimal = true;
		result.failureReason = "NONE";
		return result;
	}

	result.candidateIndices.resize(candidateRows.size());
	for (int index = 0; index < static_cast<int>(candidateRows.size()); index++)
		result.candidateIndices[index] = index;
	for (std::size_t position = 0;
		 position < result.candidateIndices.size();)
	{
		std::vector<int> trial = result.candidateIndices;
		trial.erase(trial.begin() + position);
		ZhangExactMatrix trialRows = baselineRows;
		ZhangExactVector trialValues = baselineValues;
		for (const int index : trial)
		{
			trialRows.push_back(candidateRows[index]);
			trialValues.push_back(candidateValues[index]);
		}
		if (!zhangExactAffineSystemConsistent(trialRows, trialValues))
			result.candidateIndices = std::move(trial);
		else
			position++;
	}
	result.minimal = !result.candidateIndices.empty();
	result.failureReason = result.minimal
		? "MINIMAL_AFFINE_CONFLICT_CERTIFICATE" :
		  "AFFINE_CONFLICT_HAS_NO_CANDIDATE_WITNESS";
	return result;
}

/** Distinguish an untracked deterministic integer direction from an affine
 * contradiction.  Eigenvectors are used only to locate the real nullspace;
 * authorization still fails closed and never turns the mode into an integer
 * certificate. */
inline ZhangDeterministicQuotientAudit zhangAuditDeterministicQuotientModes(
	const Eigen::VectorXd& mean,
	const Eigen::MatrixXd& covariance,
	double relativeTolerance = 1e-12,
	double integerTolerance = 1e-8)
{
	ZhangDeterministicQuotientAudit result;
	if (mean.size() == 0 || covariance.rows() != mean.size() ||
		covariance.cols() != mean.size() || !mean.allFinite() ||
		!covariance.allFinite())
	{
		result.status = "INVALID_QUOTIENT_COVARIANCE";
		return result;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(
		0.5 * (covariance + covariance.transpose()));
	if (eigen.info() != Eigen::Success)
	{
		result.status = "QUOTIENT_EIGENSOLVER_FAILED";
		return result;
	}
	const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
	const double tolerance = std::max(1e-14, relativeTolerance * largest);
	for (int mode = 0; mode < eigen.eigenvalues().size(); mode++)
	{
		if (eigen.eigenvalues()(mode) > tolerance)
		{
			result.covarianceRank++;
			continue;
		}
		result.nullity++;
	}
	// Eigenvectors have arbitrary real scale and cannot be treated as integer
	// functions.  A contradiction is asserted only for a canonical quotient
	// coordinate whose own variance is zero.  Non-axis-aligned null modes remain
	// UNTRACKED until an exact integer null row is recovered.
	bool canonicalContradiction = false;
	for (int coordinate = 0; coordinate < mean.size(); coordinate++)
	{
		if (std::abs(covariance(coordinate, coordinate)) > tolerance) continue;
		const double fractional = std::abs(
			mean(coordinate) - std::round(mean(coordinate)));
		result.maximumNullFractionalInteger = std::max(
			result.maximumNullFractionalInteger, fractional);
		canonicalContradiction |= fractional > integerTolerance;
	}
	result.covarianceValid = true;
	result.integerConsistent = !canonicalContradiction;
	result.status = result.nullity == 0 ? "FULL_RANK" :
		(result.integerConsistent ? "UNTRACKED_DETERMINISTIC_RELATION" :
			"DETERMINISTIC_INTEGER_INCONSISTENCY");
	return result;
}

inline bool zhangExactRectangularMatrix(
	const ZhangExactMatrix& matrix,
	std::size_t columns)
{
	return std::all_of(matrix.begin(), matrix.end(),
		[columns](const auto& row) { return row.size() == columns; });
}

inline ZhangExactVector zhangExactRowCombination(
	const ZhangExactVector& coefficients,
	const ZhangExactMatrix& rows)
{
	if (coefficients.size() != rows.size() || rows.empty()) return {};
	ZhangExactVector result(rows.front().size());
	for (std::size_t row = 0; row < rows.size(); row++)
	for (std::size_t column = 0; column < result.size(); column++)
	{
		result[column] += coefficients[row] * rows[row][column];
	}
	return result;
}

inline bool zhangExactLatticeContainsAll(
	const ZhangExactMatrix& lattice,
	const ZhangExactMatrix& targets)
{
	return std::all_of(targets.begin(), targets.end(),
		[&](const auto& row)
		{
			return zhangIntegerRowLatticeContains(lattice, row).contained;
		});
}

inline bool zhangExactPrimitiveRowLattice(
	const ZhangExactMatrix& rows,
	std::size_t dimension,
	int* rank = nullptr)
{
	if (!zhangExactRectangularMatrix(rows, dimension)) return false;
	if (rows.empty())
	{
		if (rank) *rank = 0;
		return true;
	}
	const auto smith = zhangIntegerRowLatticeContains(
		rows, ZhangExactVector(dimension));
	if (rank) *rank = smith.rank;
	return std::all_of(smith.smithInvariants.begin(),
		smith.smithInvariants.end(), [](const auto& invariant)
		{
			return zhangExactAbs(invariant) == 1;
		});
}

inline ZhangHeldQuotientAudit zhangExactHeldQuotientAudit(
	const ZhangExactMatrix& targetRows,
	const ZhangExactMatrix& heldRows,
	const ZhangExactVector& heldValues = {})
{
	ZhangHeldQuotientAudit result;
	if (targetRows.empty())
	{
		result.failureReason = "EMPTY_TARGET_LATTICE";
		return result;
	}
	const std::size_t dimension = targetRows.front().size();
	if (!zhangExactRectangularMatrix(targetRows, dimension) ||
		!zhangExactRectangularMatrix(heldRows, dimension) ||
		(!heldValues.empty() && heldValues.size() != heldRows.size()))
	{
		result.failureReason = "TARGET_HELD_DIMENSION_MISMATCH";
		return result;
	}
	const ZhangExactRowHnf targetHnf =
		zhangExactRowHermiteNormalForm(targetRows);
	result.targetBasis = targetHnf.basis;
	result.targetRank = static_cast<int>(result.targetBasis.size());
	if (!targetHnf.consistent || result.targetRank !=
		static_cast<int>(targetRows.size()))
	{
		result.failureReason = "TARGET_ROWS_NOT_AN_INTEGER_BASIS";
		return result;
	}
	if (heldRows.empty())
	{
		result.heldIntersectionPrimitiveInTarget = true;
		for (int coordinate = 0; coordinate < result.targetRank; coordinate++)
		{
			ZhangExactVector unit(result.targetRank);
			unit[coordinate] = 1;
			result.quotientTargetCoordinates.push_back(std::move(unit));
		}
		result.quotientRank = result.targetRank;
		result.exactClosure = true;
		result.valid = true;
		return result;
	}

	// Solve a*T = b*H exactly.  Kernel vectors of [T^T,-H^T] contain
	// target coefficients a followed by held coefficients b.
	ZhangExactMatrix relation(dimension,
		ZhangExactVector(targetRows.size() + heldRows.size()));
	for (std::size_t column = 0; column < dimension; column++)
	{
		for (std::size_t row = 0; row < targetRows.size(); row++)
			relation[column][row] = targetRows[row][column];
		for (std::size_t row = 0; row < heldRows.size(); row++)
			relation[column][targetRows.size() + row] = -heldRows[row][column];
	}
	const ZhangExactMatrix relationKernel = zhangExactIntegerKernel(relation);
	ZhangExactMatrix intersectionCoordinates;
	ZhangExactVector intersectionValues;
	for (const auto& kernelRow : relationKernel)
	{
		ZhangExactVector targetCoefficients(
			kernelRow.begin(), kernelRow.begin() + targetRows.size());
		if (std::all_of(targetCoefficients.begin(), targetCoefficients.end(),
			[](const auto& value) { return value == 0; })) continue;
		intersectionCoordinates.push_back(std::move(targetCoefficients));
		ZhangExactInteger value = 0;
		if (!heldValues.empty())
		{
			for (std::size_t row = 0; row < heldRows.size(); row++)
				value += kernelRow[targetRows.size() + row] * heldValues[row];
		}
		intersectionValues.push_back(value);
	}
	const ZhangExactRowHnf intersectionHnf = zhangExactRowHermiteNormalForm(
		intersectionCoordinates, intersectionValues);
	if (!intersectionHnf.consistent)
	{
		result.failureReason = "HELD_INTERSECTION_AFFINE_INCONSISTENCY";
		return result;
	}
	result.heldIntersectionTargetCoordinates = intersectionHnf.basis;
	result.heldIntersectionValues = intersectionHnf.values;
	result.heldIntersectionRank = static_cast<int>(intersectionHnf.basis.size());
	for (const auto& coordinateRow : intersectionHnf.basis)
	{
		result.heldIntersectionPhysicalBasis.push_back(
			zhangExactRowCombination(coordinateRow, targetRows));
	}
	result.heldIntersectionPrimitiveInTarget = zhangExactPrimitiveRowLattice(
		result.heldIntersectionTargetCoordinates, result.targetRank);
	if (!result.heldIntersectionPrimitiveInTarget)
	{
		result.failureReason = "HELD_INTERSECTION_NOT_PRIMITIVE_IN_TARGET";
		return result;
	}

	// Choose an explicit primitive complement from canonical target unit rows.
	// Every accepted addition must increase rank while retaining index one.
	ZhangExactMatrix completed = result.heldIntersectionTargetCoordinates;
	int completedRank = result.heldIntersectionRank;
	for (int coordinate = 0;
		 coordinate < result.targetRank && completedRank < result.targetRank;
		 coordinate++)
	{
		ZhangExactVector unit(result.targetRank);
		unit[coordinate] = 1;
		ZhangExactMatrix candidate = completed;
		candidate.push_back(unit);
		int candidateRank = 0;
		if (zhangExactPrimitiveRowLattice(
				candidate, result.targetRank, &candidateRank) &&
			candidateRank == completedRank + 1)
		{
			result.quotientTargetCoordinates.push_back(std::move(unit));
			completed = std::move(candidate);
			completedRank = candidateRank;
		}
	}
	result.quotientRank = static_cast<int>(
		result.quotientTargetCoordinates.size());
	result.exactClosure = completedRank == result.targetRank &&
		result.heldIntersectionRank + result.quotientRank == result.targetRank &&
		zhangExactLatticeContainsAll(completed,
			zhangExactIdentityMatrix(result.targetRank));
	if (!result.exactClosure)
	{
		result.failureReason = "PRIMITIVE_QUOTIENT_COMPLETION_FAILED";
		return result;
	}
	result.valid = true;
	return result;
}

inline ZhangCertifiedUnionAudit zhangExactCertifiedUnionAudit(
	const ZhangExactMatrix& targetRows,
	const ZhangExactMatrix& heldIntersectionRows,
	const ZhangExactVector& heldIntersectionValues,
	const ZhangExactMatrix& newlyFixedRows,
	const ZhangExactVector& newlyFixedValues)
{
	ZhangCertifiedUnionAudit result;
	if (targetRows.empty() ||
		heldIntersectionRows.size() != heldIntersectionValues.size() ||
		newlyFixedRows.size() != newlyFixedValues.size())
	{
		result.failureReason = "CERTIFIED_UNION_DIMENSION_MISMATCH";
		return result;
	}
	ZhangExactMatrix rows = heldIntersectionRows;
	rows.insert(rows.end(), newlyFixedRows.begin(), newlyFixedRows.end());
	ZhangExactVector values = heldIntersectionValues;
	values.insert(values.end(), newlyFixedValues.begin(), newlyFixedValues.end());
	const auto targetHnf = zhangExactRowHermiteNormalForm(targetRows);
	const auto heldHnf = zhangExactRowHermiteNormalForm(
		heldIntersectionRows, heldIntersectionValues);
	const auto fixedHnf = zhangExactRowHermiteNormalForm(
		newlyFixedRows, newlyFixedValues);
	const auto unionHnf = zhangExactRowHermiteNormalForm(rows, values);
	result.targetRank = targetHnf.basis.size();
	result.heldRank = heldHnf.basis.size();
	result.combinedCertifiedRank = unionHnf.basis.size();
	// This is the *incremental* certified rank, not the standalone rank of
	// the proposed batch.  Reporting fixedHnf here made a duplicate component
	// edge look like a new certificate and could inflate a closure audit.
	result.newlyFixedRank = std::max(
		0, result.combinedCertifiedRank - result.heldRank);
	result.consistent = targetHnf.consistent && heldHnf.consistent &&
		fixedHnf.consistent && unionHnf.consistent;
	if (!result.consistent)
	{
		result.failureReason = "CERTIFIED_UNION_AFFINE_INCONSISTENCY";
		return result;
	}
	result.certifiedBasis = unionHnf.basis;
	result.certifiedValues = unionHnf.values;
	result.targetContainedInCertified = zhangExactLatticeContainsAll(
		result.certifiedBasis, targetHnf.basis);
	result.certifiedContainedInTarget = zhangExactLatticeContainsAll(
		targetHnf.basis, result.certifiedBasis);
	result.exactTargetEquality = result.targetContainedInCertified &&
		result.certifiedContainedInTarget;
	if (!result.certifiedContainedInTarget)
		result.failureReason = "CERTIFIED_ROW_OUTSIDE_TARGET_LATTICE";
	else if (!result.targetContainedInCertified)
		result.failureReason = "TARGET_LATTICE_NOT_FULLY_CERTIFIED";
	return result;
}
