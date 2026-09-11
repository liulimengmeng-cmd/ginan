#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>

#include "common/zhangIntegerProductGainFrontier.hpp"

struct ProductIntegerCandidate
{
	ZhangExactVector row;
	double variance = std::numeric_limits<double>::quiet_NaN();
	double fractional = std::numeric_limits<double>::quiet_NaN();
	double perr = 1;
	double incrementalProductGain = 0;
	int dualGraphRankGain = 0;
	int signalGraphRankGain = 0;
	int productLatticeRankGain = 1;
	std::string source = "UNKNOWN";
	bool reliabilityPassed = false;
};

struct ProductIntegerCandidateGenerationResult
{
	std::vector<ProductIntegerCandidate> candidates;
	int allPairRows = 0;
	int reducedRows = 0;
	int realModeApproximations = 0;
	int reliableRows = 0;
	int dualSupportedReliableRowsBeforeTruncation = 0;
	int dualSupportedReliableRowsRetained = 0;
	int dualGraphRankGainBeforeTruncation = 0;
	int dualGraphRankGainRetained = 0;
	int quotientGainingReliableRowsBeforeTruncation = 0;
	int quotientGainingReliableRowsRetained = 0;
	int quotientRankGainBeforeTruncation = 0;
	int quotientRankGainRetained = 0;
	int dualMembershipBatchTargets = 0;
	int dualScalarDecompositionsAvoided = 0;
	int existingDualGraphRank = 0;
	std::vector<int> existingDualComponentLabels;
	int quotientMembershipBatchTargets = 0;
	int quotientScalarDecompositionsAvoided = 0;
	bool valid = false;
	std::string failureReason = "NOT_EVALUATED";
};

inline bool zhangProductCandidateIsNamedPairRow(
	const ZhangExactVector& row);

/** Recover the two graph nodes represented by an exact named-pair row.
 *
 * Product coordinates omit the local reference satellite, so a single unit
 * coefficient is an edge to the implicit node row.size(). */
inline bool zhangProductNamedPairNodes(
	const ZhangExactVector& row,
	int& firstNode,
	int& secondNode)
{
	if (!zhangProductCandidateIsNamedPairRow(row)) return false;
	firstNode = static_cast<int>(row.size());
	secondNode = static_cast<int>(row.size());
	for (int index = 0; index < static_cast<int>(row.size()); index++)
	{
		if (row[index] == 1) firstNode = index;
		if (row[index] == -1) secondNode = index;
	}
	if (firstNode == secondNode) return false;
	if (secondNode < firstNode) std::swap(firstNode, secondNode);
	return true;
}

/** Exact graph-rank represented by the currently retained completion
 * candidates relative to the pre-existing dual-frequency components. */
inline int zhangProductCandidateDualGraphRankGain(
	const std::vector<ProductIntegerCandidate>& candidates,
	const std::vector<int>& existingComponentLabels)
{
	if (candidates.empty()) return 0;
	const int nodeCount = static_cast<int>(candidates.front().row.size()) + 1;
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
			auto [iterator, inserted] = representative.try_emplace(
				existingComponentLabels[node], node);
			if (!inserted) unite(iterator->second, node);
		}
	}
	int gain = 0;
	for (const auto& candidate : candidates)
	{
		if (!candidate.reliabilityPassed || candidate.dualGraphRankGain <= 0)
			continue;
		int first = 0;
		int second = 0;
		if (zhangProductNamedPairNodes(candidate.row, first, second))
			gain += unite(first, second);
	}
	return gain;
}

/** Exact quotient rank represented by a candidate dictionary.  Counting
 * individually novel rows overstates capacity whenever several rows generate
 * the same quotient direction. */
inline int zhangProductCandidateLatticeRankGain(
	const std::vector<ProductIntegerCandidate>& candidates,
	const ZhangExactMatrix& existingSignalRows)
{
	ZhangExactMatrix united = existingSignalRows;
	for (const auto& candidate : candidates)
	{
		if (candidate.reliabilityPassed && candidate.productLatticeRankGain > 0)
			united.push_back(candidate.row);
	}
	return std::max(0,
		static_cast<int>(zhangExactRowHermiteNormalForm(united).basis.size()) -
		static_cast<int>(zhangExactRowHermiteNormalForm(
			existingSignalRows).basis.size()));
}

/** Reorder an already generated single-signal candidate dictionary so rows
 * that close a certified relation in the complementary signal are visited
 * first.
 *
 * This is an exact lattice-membership preference only.  It never changes a
 * candidate's scalar reliability decision, integer value, Perr, NIS, or
 * product gain.  A dense mixed row cannot become a graph certificate merely
 * because it belongs to the complementary lattice: only named satellite-pair
 * rows are marked as dual-graph rank gaining. */
inline int zhangPrioritiseProductCandidatesForDualLattice(
	std::vector<ProductIntegerCandidate>& candidates,
	const ZhangExactMatrix& complementarySignalRows,
	int* batchTargets = nullptr,
	int* scalarDecompositionsAvoided = nullptr,
	const ZhangExactMatrix& existingSignalRows = {},
	int* existingDualGraphRank = nullptr,
	std::vector<int>* existingDualComponentLabels = nullptr)
{
	if (batchTargets) *batchTargets = 0;
	if (scalarDecompositionsAvoided) *scalarDecompositionsAvoided = 0;
	if (existingDualGraphRank) *existingDualGraphRank = 0;
	if (existingDualComponentLabels) existingDualComponentLabels->clear();
	int supported = 0;
	std::vector<int> named;
	ZhangExactMatrix targets;
	for (int index = 0; index < static_cast<int>(candidates.size()); index++)
	{
		auto& candidate = candidates[index];
		candidate.dualGraphRankGain = 0;
		if (!zhangProductCandidateIsNamedPairRow(candidate.row) ||
			complementarySignalRows.empty())
		{
			continue;
		}
		named.push_back(index);
		targets.push_back(candidate.row);
	}
	const auto complementaryMemberships = zhangIntegerRowLatticeContainsBatch(
		complementarySignalRows, targets);

	const int dimension = candidates.empty()
		? 0 : static_cast<int>(candidates.front().row.size());
	std::vector<int> parent(std::max(0, dimension + 1));
	for (int node = 0; node < static_cast<int>(parent.size()); node++)
		parent[node] = node;
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
		first = root(first);
		second = root(second);
		if (first == second) return false;
		parent[second] = first;
		return true;
	};
	int baseRank = 0;
	for (int target = 0; target < static_cast<int>(named.size()); target++)
	{
		const bool inComplementary = target < static_cast<int>(
			complementaryMemberships.size()) &&
			complementaryMemberships[target].contained;
		// productLatticeRankGain was assigned by the preceding batched exact
		// membership query against existingSignalRows.  Reuse it here instead of
		// repeating a second multiprecision decomposition for every named pair.
		const bool inExisting = !existingSignalRows.empty() &&
			candidates[named[target]].productLatticeRankGain == 0;
		if (!inComplementary || !inExisting) continue;
		int first = 0;
		int second = 0;
		if (zhangProductNamedPairNodes(
				candidates[named[target]].row, first, second) &&
			unite(first, second)) baseRank++;
	}
	if (existingDualGraphRank) *existingDualGraphRank = baseRank;
	if (existingDualComponentLabels)
	{
		existingDualComponentLabels->resize(parent.size());
		for (int node = 0; node < static_cast<int>(parent.size()); node++)
			(*existingDualComponentLabels)[node] = root(node);
	}
	for (int target = 0; target < static_cast<int>(named.size()); target++)
	{
		auto& candidate = candidates[named[target]];
		const bool inComplementary = target < static_cast<int>(
			complementaryMemberships.size()) &&
			complementaryMemberships[target].contained;
		int first = 0;
		int second = 0;
		candidate.dualGraphRankGain = candidate.reliabilityPassed &&
			inComplementary &&
			zhangProductNamedPairNodes(candidate.row, first, second) &&
			root(first) != root(second) ? 1 : 0;
		supported += candidate.dualGraphRankGain;
	}
	if (batchTargets) *batchTargets = targets.size();
	if (scalarDecompositionsAvoided)
		*scalarDecompositionsAvoided = std::max(0,
			static_cast<int>(targets.size()) - 1);
	std::sort(candidates.begin(), candidates.end(),
		[](const auto& left, const auto& right)
		{
			if (left.reliabilityPassed != right.reliabilityPassed)
				return left.reliabilityPassed > right.reliabilityPassed;
			if (left.dualGraphRankGain != right.dualGraphRankGain)
				return left.dualGraphRankGain > right.dualGraphRankGain;
			if (left.productLatticeRankGain != right.productLatticeRankGain)
				return left.productLatticeRankGain > right.productLatticeRankGain;
			if (left.signalGraphRankGain != right.signalGraphRankGain)
				return left.signalGraphRankGain > right.signalGraphRankGain;
			if (left.incrementalProductGain != right.incrementalProductGain)
				return left.incrementalProductGain > right.incrementalProductGain;
			if (left.variance != right.variance)
				return left.variance < right.variance;
			return left.row < right.row;
		});
	return supported;
}

/** Convert a bootstrapped success probability into a conservative failure
 * bound without letting binary representation of a configured decimal
 * threshold (for example 1 - 0.999) spuriously exceed that same budget. */
inline double zhangProductFailureProbabilityBound(
	double bootstrapSuccess,
	double failureProbabilityBudget)
{
	if (!std::isfinite(bootstrapSuccess) || bootstrapSuccess <= 0)
		return 1;
	const double budget = std::clamp(failureProbabilityBudget, 0.0, 1.0);
	const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
		std::max({1.0, std::abs(bootstrapSuccess), std::abs(budget)});
	if (bootstrapSuccess > 1 + roundoff) return 1;
	const double failure = std::max(0.0, 1 - bootstrapSuccess);
	// Only collapse the sub-ulp decimal tail at the configured boundary.  The
	// previous unconditional min(budget, failure) converted every genuinely
	// unsafe candidate (for example success=0.9) into an apparent pass.
	return failure <= budget + roundoff
		? std::min(failure, budget)
		: failure;
}

inline bool zhangProductFailureProbabilityPassed(
	double failureProbability,
	double failureProbabilityBudget)
{
	return std::isfinite(failureProbability) &&
		failureProbability <= failureProbabilityBudget + 1e-12;
}

inline bool zhangNormalisePrimitiveIntegerCandidate(ZhangExactVector& row)
{
	auto gcd = [](ZhangExactInteger left, ZhangExactInteger right)
	{
		left = zhangExactAbs(left);
		right = zhangExactAbs(right);
		while (right != 0)
		{
			const auto remainder = left % right;
			left = right;
			right = remainder;
		}
		return left;
	};
	ZhangExactInteger divisor = 0;
	int first = -1;
	for (int index = 0; index < static_cast<int>(row.size()); index++)
	{
		if (row[index] != 0 && first < 0) first = index;
		divisor = gcd(divisor, row[index]);
	}
	if (first < 0 || divisor == 0) return false;
	for (auto& coefficient : row) coefficient /= divisor;
	if (row[first] < 0)
		for (auto& coefficient : row) coefficient = -coefficient;
	return true;
}

/** A product-coordinate row is a named satellite pair only when it is one
 * reference edge (+/- one unit coordinate) or one exact +/-1 difference.
 * Merely having support two is insufficient: 2*z1-z2 is a mixed lattice row,
 * not a graph certificate. */
inline bool zhangProductCandidateIsNamedPairRow(
	const ZhangExactVector& row)
{
	int positive = 0;
	int negative = 0;
	int nonzero = 0;
	for (const auto& coefficient : row)
	{
		if (coefficient == 0) continue;
		nonzero++;
		if (coefficient == 1) positive++;
		else if (coefficient == -1) negative++;
		else return false;
	}
	return (nonzero == 1 && positive + negative == 1) ||
		(nonzero == 2 && positive == 1 && negative == 1);
}

/** Generate legal primitive rows in the product coordinate itself.
 *
 * Real product-gain modes are candidate-generation guides only.  Every
 * returned row is an exact primitive integer vector and is independently
 * evaluated with scalar Perr/NIS and Mahalanobis length.  No support-count
 * gate is applied to real-mode approximations. */
inline ProductIntegerCandidateGenerationResult
generateProductIntegerCandidates(
	const Eigen::VectorXd& mean,
	const Eigen::MatrixXd& covariance,
	const Eigen::MatrixXd& productCrossCovariance,
	double maximumPerr,
	double nisAlpha,
	const ZhangExactMatrix& reducedRows = {},
	int realModeCount = 8,
	int maximumApproximationScale = 12,
	std::size_t maximumCandidates = 512,
	const ZhangExactMatrix& complementarySignalRows = {},
	const ZhangExactMatrix& existingSignalRows = {})
{
	ProductIntegerCandidateGenerationResult result;
	const int dimension = mean.size();
	if (dimension < 1 || covariance.rows() != dimension ||
		covariance.cols() != dimension ||
		productCrossCovariance.cols() != dimension || maximumPerr <= 0 ||
		maximumPerr >= 1 || nisAlpha <= 0 || nisAlpha >= 1 ||
		realModeCount < 0 || maximumApproximationScale < 1 ||
		maximumCandidates < 1)
	{
		result.failureReason = "PRODUCT_CANDIDATE_INPUT_INVALID";
		return result;
	}
	for (const auto& row : complementarySignalRows)
	if (row.size() != static_cast<std::size_t>(dimension))
	{
		result.failureReason =
			"PRODUCT_CANDIDATE_COMPLEMENTARY_LATTICE_DIMENSION_MISMATCH";
		return result;
	}
	for (const auto& row : existingSignalRows)
	if (row.size() != static_cast<std::size_t>(dimension))
	{
		result.failureReason =
			"PRODUCT_CANDIDATE_EXISTING_LATTICE_DIMENSION_MISMATCH";
		return result;
	}
	const Eigen::MatrixXd symmetric = 0.5 *
		(covariance + covariance.transpose());
	if (!symmetric.allFinite())
	{
		result.failureReason = "PRODUCT_CANDIDATE_COVARIANCE_NONFINITE";
		return result;
	}
	std::map<ZhangExactVector, std::string> rows;
	auto add = [&](ZhangExactVector row, const std::string& source)
	{
		if (row.size() != static_cast<std::size_t>(dimension) ||
			!zhangNormalisePrimitiveIntegerCandidate(row)) return;
		rows.try_emplace(std::move(row), source);
	};
	// Named star rows and every pair difference form the complete graphic seed
	// set, including the implicit reference node represented by the zero vector.
	for (int first = 0; first < dimension; first++)
	{
		ZhangExactVector unit(dimension);
		unit[first] = 1;
		add(unit, "ALL_PAIR_ROWS");
		result.allPairRows++;
		for (int second = first + 1; second < dimension; second++)
		{
			ZhangExactVector pair(dimension);
			pair[first] = 1;
			pair[second] = -1;
			add(pair, "ALL_PAIR_ROWS");
			result.allPairRows++;
		}
	}
	for (auto row : reducedRows)
	{
		add(std::move(row), "LAMBDA_REDUCED_ROWS");
		result.reducedRows++;
	}

	// M describes product information captured by one product integer row.
	const Eigen::MatrixXd information = productCrossCovariance.transpose() *
		productCrossCovariance;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qSolver(symmetric);
	if (qSolver.info() != Eigen::Success)
	{
		result.failureReason = "PRODUCT_CANDIDATE_COVARIANCE_EIGEN_FAILED";
		return result;
	}
	const double largestQ = std::max(1.0,
		qSolver.eigenvalues().cwiseAbs().maxCoeff());
	const double floor = 1e-12 * largestQ;
	Eigen::MatrixXd regularised = symmetric;
	regularised.diagonal().array() += floor;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gainSolver(
		0.5 * (information + information.transpose()), regularised);
	if (gainSolver.info() == Eigen::Success)
	{
		const int modes = std::min(realModeCount, dimension);
		for (int order = 0; order < modes; order++)
		{
			Eigen::VectorXd mode = gainSolver.eigenvectors().col(
				dimension - 1 - order);
			const double maximum = mode.cwiseAbs().maxCoeff();
			if (!(maximum > 0) || !std::isfinite(maximum)) continue;
			mode /= maximum;
			for (int scale = 1; scale <= maximumApproximationScale; scale++)
			{
				ZhangExactVector approximation(dimension);
				for (int index = 0; index < dimension; index++)
					approximation[index] = std::llround(scale * mode(index));
				const auto before = rows.size();
				add(std::move(approximation),
					"PRODUCT_GAIN_REAL_MODE_APPROXIMATION");
				result.realModeApproximations += rows.size() > before;
			}
		}
	}

	boost::math::chi_squared scalarDistribution(1);
	const double scalarThreshold = boost::math::quantile(
		boost::math::complement(scalarDistribution, nisAlpha));
	for (const auto& [row, source] : rows)
	{
		const Eigen::VectorXd numeric = zhangExactRowToDouble(row);
		const double variance =
			(numeric.transpose() * symmetric * numeric)(0, 0);
		const double floating = numeric.dot(mean);
		const double fractional = floating - std::round(floating);
		const double perr = zhangIntegerRoundFailureProbability(
			fractional, variance);
		const double nis = variance > 0
			? fractional * fractional / variance
			: std::numeric_limits<double>::infinity();
		Eigen::MatrixXd oneRow(1, dimension);
		oneRow.row(0) = numeric.transpose();
		const double gain = zhangIntegerConstraintProductGain(
			oneRow, symmetric, productCrossCovariance);
		ProductIntegerCandidate candidate;
		candidate.row = row;
		candidate.variance = variance;
		candidate.fractional = fractional;
		candidate.perr = perr;
		candidate.incrementalProductGain = std::max(0.0, gain);
		candidate.signalGraphRankGain =
			zhangProductCandidateIsNamedPairRow(row) ? 1 : 0;
		candidate.productLatticeRankGain = 1;
		// A candidate cannot create a dual-frequency graph edge on its own;
		// that rank is evaluated after its WL/L1 partner has been selected.
		candidate.dualGraphRankGain = 0;
		candidate.source = source;
		candidate.reliabilityPassed = std::isfinite(variance) && variance > 0 &&
			std::isfinite(perr) && perr <= maximumPerr && nis <= scalarThreshold;
		result.reliableRows += candidate.reliabilityPassed;
		result.candidates.push_back(std::move(candidate));
	}
	if (!existingSignalRows.empty())
	{
		ZhangExactMatrix targets;
		targets.reserve(result.candidates.size());
		for (const auto& candidate : result.candidates)
			targets.push_back(candidate.row);
		const auto memberships = zhangIntegerRowLatticeContainsBatch(
			existingSignalRows, targets);
		if (memberships.size() != result.candidates.size())
		{
			result.failureReason =
				"PRODUCT_CANDIDATE_QUOTIENT_MEMBERSHIP_BATCH_FAILED";
			return result;
		}
		for (int index = 0; index < static_cast<int>(result.candidates.size());
			 index++)
		{
			result.candidates[index].productLatticeRankGain =
				memberships[index].contained ? 0 : 1;
		}
		result.quotientMembershipBatchTargets = targets.size();
		result.quotientScalarDecompositionsAvoided = std::max(0,
			static_cast<int>(targets.size()) - 1);
	}
	result.quotientGainingReliableRowsBeforeTruncation = std::count_if(
		result.candidates.begin(), result.candidates.end(), [](const auto& candidate)
		{
			return candidate.reliabilityPassed &&
				candidate.productLatticeRankGain > 0;
		});
	if (!complementarySignalRows.empty())
	{
		// Apply the dual-frequency objective before the bounded dictionary is
		// truncated.  Reordering only after resize silently loses a low-gain pair
		// that would close an already certified relation in the other signal.
		result.dualSupportedReliableRowsBeforeTruncation =
			zhangPrioritiseProductCandidatesForDualLattice(
				result.candidates, complementarySignalRows,
				&result.dualMembershipBatchTargets,
				&result.dualScalarDecompositionsAvoided,
				existingSignalRows,
				&result.existingDualGraphRank,
				&result.existingDualComponentLabels);
	}
	else
	{
		std::sort(result.candidates.begin(), result.candidates.end(),
			[](const auto& left, const auto& right)
			{
				if (left.reliabilityPassed != right.reliabilityPassed)
					return left.reliabilityPassed > right.reliabilityPassed;
				if (left.dualGraphRankGain != right.dualGraphRankGain)
					return left.dualGraphRankGain > right.dualGraphRankGain;
				if (left.productLatticeRankGain != right.productLatticeRankGain)
					return left.productLatticeRankGain > right.productLatticeRankGain;
				if (left.signalGraphRankGain != right.signalGraphRankGain)
					return left.signalGraphRankGain > right.signalGraphRankGain;
				if (left.incrementalProductGain != right.incrementalProductGain)
					return left.incrementalProductGain > right.incrementalProductGain;
				if (left.variance != right.variance)
					return left.variance < right.variance;
				return left.row < right.row;
			});
	}
	result.dualGraphRankGainBeforeTruncation =
		zhangProductCandidateDualGraphRankGain(
			result.candidates, result.existingDualComponentLabels);
	result.quotientRankGainBeforeTruncation =
		zhangProductCandidateLatticeRankGain(
			result.candidates, existingSignalRows);
	if (result.candidates.size() > maximumCandidates)
		result.candidates.resize(maximumCandidates);
	result.dualSupportedReliableRowsRetained = std::count_if(
		result.candidates.begin(), result.candidates.end(), [](const auto& candidate)
		{
			return candidate.reliabilityPassed && candidate.dualGraphRankGain > 0;
		});
	result.quotientGainingReliableRowsRetained = std::count_if(
		result.candidates.begin(), result.candidates.end(), [](const auto& candidate)
		{
			return candidate.reliabilityPassed &&
				candidate.productLatticeRankGain > 0;
		});
	result.dualGraphRankGainRetained = zhangProductCandidateDualGraphRankGain(
		result.candidates, result.existingDualComponentLabels);
	result.quotientRankGainRetained = zhangProductCandidateLatticeRankGain(
		result.candidates, existingSignalRows);
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}
