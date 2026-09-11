#pragma once

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "common/eigenIncluder.hpp"
#include "common/zhangIncrementalRawSquareRoot.hpp"

struct ZhangPersistentRawTargetMarginal
{
	bool valid = false;
	int stateDimension = 0;
	int targetCount = 0;
	int exactConstraintsApplied = 0;
	int coordinateRepresentableTargets = 0;
	int coordinateUnrepresentableTargets = 0;
	std::size_t skippedUnrepresentableRebinds = 0;
	VectorXd mean;
	MatrixXd covariance;
	std::vector<std::string> identities;
	std::vector<std::string> physicalVersions;
	std::string failureReason;
};

struct ZhangPersistentTargetDifferenceMarginal
{
	bool valid = false;
	int availableCount = 0;
	int stochasticRank = 0;
	VectorXd mean;
	MatrixXd covariance;
	std::vector<bool> availableRows;
	std::string failureReason;
};

/** Project one persistent-target marginal into exact endpoint differences.
 * Cross-target covariance is retained, so a common stochastic gauge cancels
 * instead of being counted once per endpoint.  Identity/version mismatches are
 * missing rows; callers that authorise AR must require complete coverage. */
inline ZhangPersistentTargetDifferenceMarginal
zhangPersistentTargetDifferenceMarginal(
	const ZhangPersistentRawTargetMarginal& marginal,
	const std::vector<std::pair<std::string, std::string>>& endpointPairs,
	bool requireComplete = false,
	double relativeRankTolerance = 1e-11)
{
	ZhangPersistentTargetDifferenceMarginal result;
	result.availableRows.assign(endpointPairs.size(), false);
	if (!marginal.valid
	 || marginal.identities.size() != static_cast<std::size_t>(marginal.mean.size())
	 || marginal.physicalVersions.size() != marginal.identities.size()
	 || marginal.covariance.rows() != marginal.mean.size()
	 || marginal.covariance.cols() != marginal.mean.size()
	 || endpointPairs.empty() || !marginal.mean.allFinite()
	 || !marginal.covariance.allFinite())
	{
		result.failureReason = marginal.valid
			? "INVALID_PERSISTENT_TARGET_DIFFERENCE_INPUT"
			: (marginal.failureReason.empty()
				? "PERSISTENT_TARGET_MARGINAL_INVALID"
				: marginal.failureReason);
		return result;
	}
	std::map<std::string, int> indices;
	for (int index = 0; index < static_cast<int>(marginal.identities.size()); index++)
	{
		if (marginal.identities[index].empty()
		 || marginal.physicalVersions[index] != marginal.identities[index]
		 || !indices.emplace(marginal.identities[index], index).second)
		{
			result.failureReason =
				"DUPLICATE_OR_MISMATCHED_PERSISTENT_TARGET_IDENTITY";
			return result;
		}
	}
	MatrixXd transform = MatrixXd::Zero(
		endpointPairs.size(), marginal.mean.size());
	for (int row = 0; row < static_cast<int>(endpointPairs.size()); row++)
	{
		const auto first = indices.find(endpointPairs[row].first);
		const auto second = indices.find(endpointPairs[row].second);
		if (first == indices.end() || second == indices.end()
		 || first->second == second->second)
		{
			continue;
		}
		transform(row, first->second) = -1;
		transform(row, second->second) = +1;
		result.availableRows[row] = true;
		result.availableCount++;
	}
	if (result.availableCount == 0
	 || (requireComplete && result.availableCount != transform.rows()))
	{
		result.failureReason = result.availableCount == 0
			? "PERSISTENT_TARGET_DIFFERENCE_IDENTITIES_MISSING"
			: "INCOMPLETE_PERSISTENT_TARGET_DIFFERENCE_IDENTITIES";
		return result;
	}
	result.mean = transform * marginal.mean;
	result.covariance = transform * marginal.covariance * transform.transpose();
	result.covariance = 0.5
		* (result.covariance + result.covariance.transpose());
	if (!result.mean.allFinite() || !result.covariance.allFinite())
	{
		result.failureReason = "NONFINITE_PERSISTENT_TARGET_DIFFERENCE_MARGINAL";
		return result;
	}
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(result.covariance);
	if (eigen.info() != Eigen::Success)
	{
		result.failureReason = "PERSISTENT_TARGET_DIFFERENCE_EIGEN_FAILED";
		return result;
	}
	const double scale = std::max(
		1.0, eigen.eigenvalues().cwiseAbs().maxCoeff());
	if (eigen.eigenvalues().minCoeff() < -100 * relativeRankTolerance * scale)
	{
		result.failureReason = "PERSISTENT_TARGET_DIFFERENCE_COVARIANCE_NOT_PSD";
		return result;
	}
	result.stochasticRank =
		(eigen.eigenvalues().array() > relativeRankTolerance * scale).count();
	result.failureReason = result.availableCount == transform.rows()
		? "NONE" : "PARTIAL_PERSISTENT_TARGET_DIFFERENCE_MARGINAL";
	result.valid = true;
	return result;
}

/** Raw-factor square-root boundary augmented with persistent physical target
 * variables.  The augmented order is [current network state, persistent
 * targets].  Targets are introduced by an exact affine transform, carried
 * with identity transition and zero process noise, and reconnected to each
 * epoch's current state coordinates by zero-variance exact constraints. */
class ZhangPersistentRawTargetWindow
{
public:
	void clear()
	{
		initialised = false;
		stateDimension = 0;
		targets.clear();
		targetIndex.clear();
		window.clear();
		failureReason.clear();
		skippedUnrepresentableRebinds = 0;
	}

	bool initialise(const VectorXd& mean, const MatrixXd& covariance)
	{
		clear();
		if (!window.initialise(mean, covariance))
		{
			return fail(window.summary().failureReason);
		}
		stateDimension = mean.size();
		initialised = true;
		return true;
	}

	bool addAcceptedMeasurement(
		const MatrixXd& design,
		const MatrixXd& covariance,
		const VectorXd& observation)
	{
		if (!initialised || design.cols() != stateDimension)
		{
			return fail("INVALID_PERSISTENT_TARGET_MEASUREMENT");
		}
		MatrixXd augmented = MatrixXd::Zero(
			design.rows(), stateDimension + targets.size());
		augmented.leftCols(stateDimension) = design;
		if (!window.addAcceptedMeasurement(
				augmented, covariance, observation))
		{
			return fail(window.summary().failureReason);
		}
		return true;
	}

	bool advance(
		const MatrixXd& transition,
		const MatrixXd& processCovariance,
		const VectorXd& affineOffset = {})
	{
		if (!initialised || transition.cols() != stateDimension
		 || processCovariance.rows() != transition.rows()
		 || processCovariance.cols() != transition.rows())
		{
			return fail("INVALID_PERSISTENT_TARGET_TRANSITION");
		}
		if (affineOffset.size() != 0
		 && (affineOffset.size() != transition.rows()
			 || !affineOffset.allFinite()))
		{
			return fail("INVALID_PERSISTENT_TARGET_TRANSITION_OFFSET");
		}
		const int oldStateDimension = stateDimension;
		const int newStateDimension = transition.rows();
		const int targetCount = targets.size();
		MatrixXd augmentedTransition = MatrixXd::Zero(
			newStateDimension + targetCount,
			oldStateDimension + targetCount);
		augmentedTransition.topLeftCorner(
			newStateDimension, oldStateDimension) = transition;
		augmentedTransition.bottomRightCorner(targetCount, targetCount)
			= MatrixXd::Identity(targetCount, targetCount);
		MatrixXd augmentedProcess = MatrixXd::Zero(
			newStateDimension + targetCount,
			newStateDimension + targetCount);
		augmentedProcess.topLeftCorner(
			newStateDimension, newStateDimension) = processCovariance;
		VectorXd augmentedOffset = VectorXd::Zero(
			newStateDimension + targetCount);
		if (affineOffset.size() != 0)
		{
			augmentedOffset.head(newStateDimension) = affineOffset;
		}
		if (!window.advance(
				augmentedTransition, augmentedProcess, augmentedOffset))
		{
			return fail(window.summary().failureReason);
		}
		stateDimension = newStateDimension;
		return true;
	}

	bool applyExactCoordinateTransform(
		const MatrixXd& transform,
		const VectorXd& translation = {})
	{
		if (!initialised || transform.cols() != stateDimension)
		{
			return fail("INVALID_PERSISTENT_TARGET_COORDINATE_TRANSFORM");
		}
		const int oldStateDimension = stateDimension;
		const int newStateDimension = transform.rows();
		const int targetCount = targets.size();
		MatrixXd augmentedTransform = MatrixXd::Zero(
			newStateDimension + targetCount,
			oldStateDimension + targetCount);
		augmentedTransform.topLeftCorner(
			newStateDimension, oldStateDimension) = transform;
		augmentedTransform.bottomRightCorner(targetCount, targetCount)
			= MatrixXd::Identity(targetCount, targetCount);
		VectorXd augmentedTranslation = VectorXd::Zero(
			newStateDimension + targetCount);
		VectorXd coordinateTranslation = VectorXd::Zero(newStateDimension);
		if (translation.size() != 0)
		{
			if (translation.size() != newStateDimension)
			{
				return fail("INVALID_PERSISTENT_TARGET_COORDINATE_TRANSLATION");
			}
			augmentedTranslation.head(newStateDimension) = translation;
			coordinateTranslation = translation;
		}
		std::vector<VectorXd> transportedRows;
		std::vector<double> transportedOffsets;
		std::vector<bool> representationAvailable;
		transportedRows.reserve(targets.size());
		transportedOffsets.reserve(targets.size());
		representationAvailable.reserve(targets.size());

		// Every target solves the same equation
		//
		//     transform^T r_new = r_old.
		//
		// Factorising transform^T once per target is both unnecessary and, for
		// a global network, can retain several gigabytes of decomposition
		// workspace in one state-transition call.  A single rank-revealing COD
		// with all target covectors as right-hand sides is algebraically
		// identical and gives each target the same representability test.
		MatrixXd transportedMatrix = MatrixXd::Zero(
			newStateDimension, targetCount);
		if (targetCount > 0)
		{
			MatrixXd oldRows = MatrixXd::Zero(oldStateDimension, targetCount);
			for (int index = 0; index < targetCount; index++)
			{
				if (targets[index].coordinateRepresentationAvailable)
				{
					oldRows.col(index) = targets[index].currentStateRow;
				}
			}
			MatrixXd transformTranspose = transform.transpose();
			Eigen::CompleteOrthogonalDecomposition<MatrixXd> decomposition(
				transformTranspose);
			transportedMatrix = decomposition.solve(oldRows);
		}
		for (int index = 0; index < targetCount; index++)
		{
			const auto& target = targets[index];
			VectorXd transported = transportedMatrix.col(index);
			const double scale = std::max(1.0, target.currentStateRow.norm());
			const bool representable =
				target.coordinateRepresentationAvailable &&
				transported.allFinite() &&
				(transform.transpose() * transported
					- target.currentStateRow).norm() <= 1e-10 * scale;
			transportedRows.push_back(representable
				? transported : VectorXd::Zero(newStateDimension));
			transportedOffsets.push_back(representable
				? target.currentOffset - transported.dot(coordinateTranslation)
				: target.currentOffset);
			representationAvailable.push_back(representable);
		}
		if (!window.applyExactCoordinateTransform(
				augmentedTransform, augmentedTranslation))
		{
			return fail(window.summary().failureReason);
		}
		stateDimension = newStateDimension;
		for (int index = 0; index < static_cast<int>(targets.size()); index++)
		{
			targets[index].currentStateRow = std::move(transportedRows[index]);
			targets[index].currentOffset = transportedOffsets[index];
			targets[index].coordinateRepresentationAvailable =
				representationAvailable[index];
			targets[index].coordinateTransformPendingValidation = true;
		}
		return true;
	}

	/** Bind or rebind one canonical target at an event sequence.  The
	 * physicalVersion must remain unchanged; a changed version is a hard reset
	 * boundary and is never mixed into the current window. */
	bool bindTarget(
		const std::string& identity,
		const std::string& physicalVersion,
		const VectorXd& row,
		double offset,
		std::size_t eventSequence)
	{
		if (!initialised || identity.empty() || physicalVersion.empty()
		 || row.size() != stateDimension || !row.allFinite()
		 || !std::isfinite(offset))
		{
			return fail("INVALID_PERSISTENT_RAW_TARGET");
		}
		auto found = targetIndex.find(identity);
		if (found == targetIndex.end())
		{
			MatrixXd rows(1, stateDimension);
			rows.row(0) = row.transpose();
			return bindNewTargets(
				{identity}, {physicalVersion}, rows,
				VectorXd::Constant(1, offset), eventSequence);
		}
		auto& target = targets[found->second];
		if (target.physicalVersion != physicalVersion)
		{
			return fail("PERSISTENT_RAW_TARGET_PHYSICAL_VERSION_CHANGED");
		}
		if (target.lastBoundSequence == eventSequence)
		{
			return true;
		}
		if (!target.coordinateRepresentationAvailable)
		{
			// Once a rectangular projection has removed part of the exact
			// covector, later current-state rows are not proof that the lost
			// historical functional has reappeared.  The explicit target
			// variable remains queryable, but it must not receive a new
			// zero-noise constraint under the old identity.
			target.coordinateTransformPendingValidation = false;
			target.lastBoundSequence = eventSequence;
			skippedUnrepresentableRebinds++;
			return true;
		}
		if (target.coordinateTransformPendingValidation)
		{
			target.coordinateTransformPendingValidation = false;
			const double rowScale = std::max({
				1.0, row.norm(), target.currentStateRow.norm()});
			if ((row - target.currentStateRow).norm() > 1e-10 * rowScale ||
				std::abs(offset - target.currentOffset) > 1e-10 *
					std::max({1.0, std::abs(offset),
						std::abs(target.currentOffset)}))
			{
				return fail("PERSISTENT_RAW_TARGET_REBIND_MISMATCH");
			}
		}
		MatrixXd constraint = MatrixXd::Zero(
			1, stateDimension + targets.size());
		constraint.leftCols(stateDimension) = -row.transpose();
		constraint(0, stateDimension + found->second) = 1;
		if (!window.applyExactConstraint(
				constraint, VectorXd::Constant(1, offset)))
		{
			return fail(window.summary().failureReason);
		}
		target.lastBoundSequence = eventSequence;
		target.currentStateRow = row;
		target.currentOffset = offset;
		target.coordinateRepresentationAvailable = true;
		return true;
	}

	/** Introduce a same-epoch block of new immutable targets with one exact
	 * affine augmentation.  Applying K one-row transforms is algebraically
	 * identical, but repeats the full square-root boundary decomposition K
	 * times.  This block form performs that decomposition once. */
	bool bindNewTargets(
		const std::vector<std::string>& identities,
		const std::vector<std::string>& physicalVersions,
		const MatrixXd& rows,
		const VectorXd& offsets,
		std::size_t eventSequence)
	{
		const int count = identities.size();
		if (!initialised || count == 0
		 || physicalVersions.size() != identities.size()
		 || rows.rows() != count || rows.cols() != stateDimension
		 || offsets.size() != count || !rows.allFinite()
		 || !offsets.allFinite())
		{
			return fail("INVALID_PERSISTENT_RAW_TARGET_BATCH");
		}
		std::map<std::string, int> batchIndex;
		for (int index = 0; index < count; index++)
		{
			if (identities[index].empty() || physicalVersions[index].empty())
			{
				return fail("INVALID_PERSISTENT_RAW_TARGET_BATCH_IDENTITY");
			}
			if (targetIndex.find(identities[index]) != targetIndex.end())
			{
				return fail("PERSISTENT_RAW_TARGET_BATCH_ALREADY_BOUND");
			}
			if (!batchIndex.emplace(identities[index], index).second)
			{
				return fail("DUPLICATE_PERSISTENT_RAW_TARGET_BATCH_IDENTITY");
			}
		}

		const int oldDimension = stateDimension + targets.size();
		MatrixXd augment = MatrixXd::Zero(
			oldDimension + count, oldDimension);
		augment.topRows(oldDimension) =
			MatrixXd::Identity(oldDimension, oldDimension);
		augment.bottomLeftCorner(count, stateDimension) = rows;
		VectorXd translation = VectorXd::Zero(oldDimension + count);
		translation.tail(count) = offsets;
		if (!window.applyExactCoordinateTransform(augment, translation))
		{
			return fail(window.summary().failureReason);
		}
		for (int index = 0; index < count; index++)
		{
			targetIndex[identities[index]] = targets.size();
			targets.push_back({
				identities[index], physicalVersions[index], eventSequence,
				rows.row(index).transpose(), offsets(index), true, false});
		}
		return true;
	}

	ZhangPersistentRawTargetMarginal targetMarginal() const
	{
		ZhangPersistentRawTargetMarginal result;
		result.stateDimension = stateDimension;
		result.targetCount = targets.size();
		result.exactConstraintsApplied =
			window.summary().exactConstraintsApplied;
		result.skippedUnrepresentableRebinds =
			skippedUnrepresentableRebinds;
		for (const auto& target : targets)
		{
			result.identities.push_back(target.identity);
			result.physicalVersions.push_back(target.physicalVersion);
			if (target.coordinateRepresentationAvailable)
			{
				result.coordinateRepresentableTargets++;
			}
			else
			{
				result.coordinateUnrepresentableTargets++;
			}
		}
		if (!initialised || targets.empty())
		{
			result.failureReason = "NO_PERSISTENT_RAW_TARGETS";
			return result;
		}
		const auto selected = window.coordinateBlockMarginal(
			stateDimension,
			static_cast<int>(targets.size()),
			VectorXd::Zero(targets.size()));
		if (!selected.valid)
		{
			result.failureReason = selected.failureReason.empty()
				? "PERSISTENT_RAW_TARGET_MARGINAL_FAILED"
				: selected.failureReason;
			return result;
		}
		result.mean = selected.mean;
		result.covariance = selected.covariance;
		result.covariance = 0.5
			* (result.covariance + result.covariance.transpose());
		result.valid = result.mean.allFinite() && result.covariance.allFinite();
		if (!result.valid)
		{
			result.failureReason = "NONFINITE_PERSISTENT_RAW_TARGET_MARGINAL";
		}
		return result;
	}

	/** Direct joint marginal of requested physical target differences.  This is
	 * the production path for current satellite-minus-reference L1/L2 product
	 * relations: it preserves all cross-target covariance but never constructs
	 * the covariance of unrelated retained history. */
	ZhangPersistentTargetDifferenceMarginal targetDifferenceMarginal(
		const std::vector<std::pair<std::string, std::string>>& endpointPairs,
		bool requireComplete = false) const
	{
		ZhangPersistentTargetDifferenceMarginal result;
		result.availableRows.assign(endpointPairs.size(), false);
		if (!initialised || endpointPairs.empty())
		{
			result.failureReason = "INVALID_PERSISTENT_TARGET_DIFFERENCE_QUERY";
			return result;
		}
		std::map<std::string, int> indices;
		for (int index = 0; index < static_cast<int>(targets.size()); index++)
		{
			if (targets[index].identity.empty()
			 || targets[index].physicalVersion != targets[index].identity
			 || !indices.emplace(targets[index].identity, index).second)
			{
				result.failureReason =
					"DUPLICATE_OR_MISMATCHED_PERSISTENT_TARGET_IDENTITY";
				return result;
			}
		}
		std::vector<std::pair<int, int>> availablePairs;
		std::vector<int> availableOutputRows;
		for (int row = 0; row < static_cast<int>(endpointPairs.size()); row++)
		{
			const auto negative = indices.find(endpointPairs[row].first);
			const auto positive = indices.find(endpointPairs[row].second);
			if (negative == indices.end() || positive == indices.end()
			 || negative->second == positive->second)
			{
				continue;
			}
			availablePairs.emplace_back(negative->second, positive->second);
			availableOutputRows.push_back(row);
			result.availableRows[row] = true;
		}
		result.availableCount = static_cast<int>(availablePairs.size());
		if (availablePairs.empty()
		 || (requireComplete && availablePairs.size() != endpointPairs.size()))
		{
			result.failureReason = availablePairs.empty()
				? "PERSISTENT_TARGET_DIFFERENCE_IDENTITIES_MISSING"
				: "INCOMPLETE_PERSISTENT_TARGET_DIFFERENCE_IDENTITIES";
			return result;
		}
		const auto selected = window.coordinateDifferenceMarginal(
			stateDimension, static_cast<int>(targets.size()), availablePairs,
			VectorXd::Zero(availablePairs.size()));
		if (!selected.valid)
		{
			result.failureReason = selected.failureReason.empty()
				? "PERSISTENT_TARGET_DIFFERENCE_MARGINAL_FAILED"
				: selected.failureReason;
			return result;
		}
		result.mean = VectorXd::Zero(endpointPairs.size());
		result.covariance = MatrixXd::Zero(
			endpointPairs.size(), endpointPairs.size());
		for (int row = 0; row < static_cast<int>(availablePairs.size()); row++)
		{
			result.mean(availableOutputRows[row]) = selected.mean(row);
			for (int column = 0;
				 column < static_cast<int>(availablePairs.size()); column++)
			{
				result.covariance(
					availableOutputRows[row], availableOutputRows[column]) =
					selected.covariance(row, column);
			}
		}
		result.stochasticRank = selected.targetRank;
		result.failureReason = availablePairs.size() == endpointPairs.size()
			? "NONE" : "PARTIAL_PERSISTENT_TARGET_DIFFERENCE_MARGINAL";
		result.valid = true;
		return result;
	}

	/** Marginalise released snapshot variables while preserving the exact
	 * current Gaussian marginal of the network state and retained targets.
	 * This is the lifecycle operation corresponding to reference-count expiry;
	 * it never drops information carried by a still-referenced target. */
	bool retainTargets(const std::set<std::string>& retainedIdentities)
	{
		if (!initialised)
		{
			return fail("PERSISTENT_RAW_TARGET_WINDOW_NOT_INITIALISED");
		}
		std::vector<int> retainedTargetIndices;
		for (int index = 0; index < static_cast<int>(targets.size()); index++)
		{
			if (retainedIdentities.count(targets[index].identity) > 0)
			{
				retainedTargetIndices.push_back(index);
			}
		}
		if (retainedTargetIndices.size() == targets.size())
		{
			return true;
		}

		VectorXd fullMean;
		MatrixXd fullCovariance;
		if (!window.currentMarginal(fullMean, fullCovariance))
		{
			return fail(window.summary().failureReason.empty()
				? "PERSISTENT_RAW_TARGET_LIFECYCLE_MARGINAL_FAILED"
				: window.summary().failureReason);
		}
		const int retainedDimension =
			stateDimension + retainedTargetIndices.size();
		VectorXd retainedMean(retainedDimension);
		MatrixXd retainedCovariance(retainedDimension, retainedDimension);
		std::vector<int> sourceIndices;
		sourceIndices.reserve(retainedDimension);
		for (int index = 0; index < stateDimension; index++)
		{
			sourceIndices.push_back(index);
		}
		for (int targetIndex : retainedTargetIndices)
		{
			sourceIndices.push_back(stateDimension + targetIndex);
		}
		for (int row = 0; row < retainedDimension; row++)
		{
			retainedMean(row) = fullMean(sourceIndices[row]);
			for (int column = 0; column < retainedDimension; column++)
			{
				retainedCovariance(row, column) = fullCovariance(
					sourceIndices[row], sourceIndices[column]);
			}
		}
		retainedCovariance = 0.5 *
			(retainedCovariance + retainedCovariance.transpose());
		if (!window.initialise(retainedMean, retainedCovariance))
		{
			return fail(window.summary().failureReason.empty()
				? "PERSISTENT_RAW_TARGET_LIFECYCLE_REANCHOR_FAILED"
				: window.summary().failureReason);
		}

		std::vector<Target> retainedTargets;
		retainedTargets.reserve(retainedTargetIndices.size());
		for (int index : retainedTargetIndices)
		{
			retainedTargets.push_back(targets[index]);
		}
		targets = std::move(retainedTargets);
		targetIndex.clear();
		for (int index = 0; index < static_cast<int>(targets.size()); index++)
		{
			targetIndex[targets[index].identity] = index;
		}
		failureReason.clear();
		return true;
	}

	int currentStateDimension() const { return stateDimension; }
	int targetCount() const { return targets.size(); }
	const std::string& lastFailureReason() const { return failureReason; }
	const ZhangIncrementalRawSquareRootSummary summary() const
	{
		return window.summary();
	}

private:
	struct Target
	{
		std::string identity;
		std::string physicalVersion;
		std::size_t lastBoundSequence = 0;
		VectorXd currentStateRow;
		double currentOffset = 0;
		bool coordinateRepresentationAvailable = true;
		bool coordinateTransformPendingValidation = false;
	};

	bool fail(const std::string& reason)
	{
		failureReason = reason.empty()
			? "UNKNOWN_PERSISTENT_RAW_TARGET_FAILURE" : reason;
		return false;
	}

	bool initialised = false;
	int stateDimension = 0;
	std::vector<Target> targets;
	std::map<std::string, int> targetIndex;
	ZhangIncrementalRawSquareRoot window;
	std::string failureReason;
	std::size_t skippedUnrepresentableRebinds = 0;
};
