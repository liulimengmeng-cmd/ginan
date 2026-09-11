#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/**
 * Small, fail-closed hypothesis set used only after an integer conflict.
 *
 * A hypothesis is identified by the exact canonical affine-lattice
 * fingerprint supplied by the caller.  Confirmation always requires the same
 * hypothesis on a strictly later epoch with a different evidence fingerprint;
 * multiple observations from one epoch can never confirm a hypothesis.
 */
enum class ZhangConflictHypothesisStatus
{
	INVALID_INPUT,
	REJECTED_NOT_RELIABLE,
	INSERTED_PENDING,
	SAME_EPOCH_PENDING,
	DEPENDENT_EVIDENCE_PENDING,
	CONFIRMED,
	ALREADY_CONFIRMED,
	NON_MONOTONIC_EPOCH,
	CAPACITY_REJECTED
};

inline const char* zhangConflictHypothesisStatusName(
	ZhangConflictHypothesisStatus status)
{
	switch (status)
	{
		case ZhangConflictHypothesisStatus::INVALID_INPUT:
			return "INVALID_INPUT";
		case ZhangConflictHypothesisStatus::REJECTED_NOT_RELIABLE:
			return "REJECTED_NOT_RELIABLE";
		case ZhangConflictHypothesisStatus::INSERTED_PENDING:
			return "INSERTED_PENDING";
		case ZhangConflictHypothesisStatus::SAME_EPOCH_PENDING:
			return "SAME_EPOCH_PENDING";
		case ZhangConflictHypothesisStatus::DEPENDENT_EVIDENCE_PENDING:
			return "DEPENDENT_EVIDENCE_PENDING";
		case ZhangConflictHypothesisStatus::CONFIRMED:
			return "CONFIRMED";
		case ZhangConflictHypothesisStatus::ALREADY_CONFIRMED:
			return "ALREADY_CONFIRMED";
		case ZhangConflictHypothesisStatus::NON_MONOTONIC_EPOCH:
			return "NON_MONOTONIC_EPOCH";
		case ZhangConflictHypothesisStatus::CAPACITY_REJECTED:
			return "CAPACITY_REJECTED";
	}
	return "UNKNOWN";
}

struct ZhangConflictPendingHypothesis
{
	std::string hypothesisFingerprint;
	long firstEpoch = 0;
	/** Last observed epoch, used only to reject non-monotonic input. */
	long lastEpoch = 0;
	/**
	 * Last epoch that supplied a previously unseen evidence fingerprint.
	 * Same-epoch and repeated evidence must never extend this lifetime.
	 */
	long lastIndependentEvidenceEpoch = 0;
	long confirmationEpoch = 0;
	int independentEpochs = 0;
	bool confirmed = false;
	std::set<std::string> evidenceFingerprints;
};

struct ZhangConflictHypothesisUpdate
{
	ZhangConflictHypothesisStatus status =
		ZhangConflictHypothesisStatus::INVALID_INPUT;
	bool acceptedObservation = false;
	bool inserted = false;
	bool sameEpoch = false;
	bool dependentEvidence = false;
	bool independentEvidence = false;
	bool confirmed = false;
	int independentEpochs = 0;
	std::size_t candidateCount = 0;
	long firstEpoch = 0;
	long lastIndependentEvidenceEpoch = 0;
	long confirmationEpoch = 0;
};

class ZhangConflictAwareHypothesisSet
{
public:
	static constexpr std::size_t maximumCandidates = 4;

	ZhangConflictHypothesisUpdate observe(
		const std::string& hypothesisFingerprint,
		long epoch,
		const std::string& evidenceFingerprint,
		bool reliable)
	{
		ZhangConflictHypothesisUpdate result;
		auto finish = [&]()
		{
			result.candidateCount = candidates_.size();
			return result;
		};
		if (hypothesisFingerprint.empty() || evidenceFingerprint.empty() ||
			epoch <= 0)
		{
			result.status = ZhangConflictHypothesisStatus::INVALID_INPUT;
			return finish();
		}
		if (!reliable)
		{
			result.status =
				ZhangConflictHypothesisStatus::REJECTED_NOT_RELIABLE;
			return finish();
		}

		auto found = std::find_if(candidates_.begin(), candidates_.end(),
			[&](const auto& candidate)
			{
				return candidate.hypothesisFingerprint ==
					hypothesisFingerprint;
			});
		if (found == candidates_.end())
		{
			if (candidates_.size() >= maximumCandidates)
			{
				result.status =
					ZhangConflictHypothesisStatus::CAPACITY_REJECTED;
				return finish();
			}
			ZhangConflictPendingHypothesis candidate;
			candidate.hypothesisFingerprint = hypothesisFingerprint;
			candidate.firstEpoch = epoch;
			candidate.lastEpoch = epoch;
			candidate.lastIndependentEvidenceEpoch = epoch;
			candidate.independentEpochs = 1;
			candidate.evidenceFingerprints.insert(evidenceFingerprint);
			candidates_.push_back(std::move(candidate));
			std::sort(candidates_.begin(), candidates_.end(),
				[](const auto& left, const auto& right)
				{
					return left.hypothesisFingerprint <
						right.hypothesisFingerprint;
				});
			result.status =
				ZhangConflictHypothesisStatus::INSERTED_PENDING;
			result.acceptedObservation = true;
			result.inserted = true;
			result.independentEpochs = 1;
			result.firstEpoch = epoch;
			result.lastIndependentEvidenceEpoch = epoch;
			return finish();
		}

		result.firstEpoch = found->firstEpoch;
		result.lastIndependentEvidenceEpoch =
			found->lastIndependentEvidenceEpoch;
		result.confirmationEpoch = found->confirmationEpoch;
		result.independentEpochs = found->independentEpochs;
		result.confirmed = found->confirmed;
		if (epoch < found->lastEpoch)
		{
			result.status =
				ZhangConflictHypothesisStatus::NON_MONOTONIC_EPOCH;
			return finish();
		}
		if (epoch == found->lastEpoch)
		{
			// Keep the provenance for later dependency checks, but never count
			// two observations from the same epoch as independent confirmation
			// or refresh lastIndependentEvidenceEpoch.
			found->evidenceFingerprints.insert(evidenceFingerprint);
			result.status = found->confirmed
				? ZhangConflictHypothesisStatus::ALREADY_CONFIRMED
				: ZhangConflictHypothesisStatus::SAME_EPOCH_PENDING;
			result.acceptedObservation = true;
			result.sameEpoch = true;
			return finish();
		}

		found->lastEpoch = epoch;
		const bool freshEvidence =
			found->evidenceFingerprints.insert(evidenceFingerprint).second;
		if (!freshEvidence)
		{
			// lastEpoch records ordering only.  In particular, dependent
			// evidence cannot keep a pending or confirmed slot alive.
			result.status = found->confirmed
				? ZhangConflictHypothesisStatus::ALREADY_CONFIRMED
				: ZhangConflictHypothesisStatus::DEPENDENT_EVIDENCE_PENDING;
			result.acceptedObservation = true;
			result.dependentEvidence = true;
			return finish();
		}

		found->lastIndependentEvidenceEpoch = epoch;
		found->independentEpochs++;
		result.independentEvidence = true;
		result.independentEpochs = found->independentEpochs;
		result.lastIndependentEvidenceEpoch = epoch;
		if (found->confirmed)
		{
			result.status =
				ZhangConflictHypothesisStatus::ALREADY_CONFIRMED;
			result.acceptedObservation = true;
			result.confirmed = true;
			return finish();
		}

		found->confirmed = true;
		found->confirmationEpoch = epoch;
		result.status = ZhangConflictHypothesisStatus::CONFIRMED;
		result.acceptedObservation = true;
		result.confirmed = true;
		result.confirmationEpoch = epoch;
		return finish();
	}

	/**
	 * Expire pending and confirmed slots using only the last independent
	 * evidence epoch.  The caller supplies an inclusive lower bound to retain.
	 */
	std::size_t eraseExpiredBefore(
		long minimumLastIndependentEvidenceEpoch)
	{
		const auto previous = candidates_.size();
		candidates_.erase(std::remove_if(
			candidates_.begin(), candidates_.end(),
			[&](const auto& candidate)
			{
				return candidate.lastIndependentEvidenceEpoch <
					minimumLastIndependentEvidenceEpoch;
			}), candidates_.end());
		return previous - candidates_.size();
	}

	/** Compatibility API for callers that intentionally retain confirmations. */
	std::size_t eraseUnconfirmedBefore(long minimumLastEpoch)
	{
		const auto previous = candidates_.size();
		candidates_.erase(std::remove_if(
			candidates_.begin(), candidates_.end(),
			[&](const auto& candidate)
			{
				return !candidate.confirmed &&
					candidate.lastIndependentEvidenceEpoch < minimumLastEpoch;
			}), candidates_.end());
		return previous - candidates_.size();
	}

	/** Explicitly reclaim a pending or confirmed slot after consume/reject. */
	bool eraseHypothesis(const std::string& hypothesisFingerprint)
	{
		const auto found = std::find_if(candidates_.begin(), candidates_.end(),
			[&](const auto& candidate)
			{
				return candidate.hypothesisFingerprint ==
					hypothesisFingerprint;
			});
		if (found == candidates_.end()) return false;
		candidates_.erase(found);
		return true;
	}

	/** Consume only a confirmed hypothesis; pending entries fail closed. */
	bool consumeConfirmed(const std::string& hypothesisFingerprint)
	{
		const auto found = std::find_if(candidates_.begin(), candidates_.end(),
			[&](const auto& candidate)
			{
				return candidate.hypothesisFingerprint ==
					hypothesisFingerprint;
			});
		if (found == candidates_.end() || !found->confirmed) return false;
		candidates_.erase(found);
		return true;
	}

	void clear()
	{
		candidates_.clear();
	}

	const std::vector<ZhangConflictPendingHypothesis>& candidates() const
	{
		return candidates_;
	}

private:
	std::vector<ZhangConflictPendingHypothesis> candidates_;
};

/** One potential dual-frequency bridge between two current components. */
struct ZhangConflictAwareForestEdge
{
	std::string id;
	int firstNode = -1;
	int secondNode = -1;
	int firstComponent = -1;
	int secondComponent = -1;
	int mergedSize = 0;
	double covarianceTrace = std::numeric_limits<double>::infinity();
};

/** A deterministic, acyclic family of two or three component bridges. */
struct ZhangConflictAwareForestFamily
{
	bool valid = false;
	std::string fingerprint;
	std::vector<std::size_t> inputEdgeIndices;
	std::vector<std::string> edgeFingerprints;
	int trueMergeCount = 0;
	int largestMergedComponentCount = 0;
	int touchedComponentCount = 0;
	int coveredNodeCount = 0;
	int maximumMergedSize = 0;
	long long summedMergedSize = 0;
	double covarianceTrace = 0;
};

namespace zhang_conflict_aware_detail
{
struct NormalizedForestEdge
{
	ZhangConflictAwareForestEdge edge;
	std::size_t inputIndex = 0;
	std::string fingerprint;
};

inline std::string forestEdgeFingerprint(
	const ZhangConflictAwareForestEdge& edge)
{
	std::ostringstream stream;
	stream << "N" << edge.firstNode << '-' << edge.secondNode
		<< "|C" << edge.firstComponent << '-' << edge.secondComponent;
	if (!edge.id.empty()) stream << "|I" << edge.id;
	return stream.str();
}

inline bool forestEdgeBetter(
	const NormalizedForestEdge& left,
	const NormalizedForestEdge& right)
{
	if (left.edge.mergedSize != right.edge.mergedSize)
		return left.edge.mergedSize > right.edge.mergedSize;
	if (left.edge.covarianceTrace != right.edge.covarianceTrace)
		return left.edge.covarianceTrace < right.edge.covarianceTrace;
	if (left.edge.firstComponent != right.edge.firstComponent)
		return left.edge.firstComponent < right.edge.firstComponent;
	if (left.edge.secondComponent != right.edge.secondComponent)
		return left.edge.secondComponent < right.edge.secondComponent;
	if (left.edge.firstNode != right.edge.firstNode)
		return left.edge.firstNode < right.edge.firstNode;
	if (left.edge.secondNode != right.edge.secondNode)
		return left.edge.secondNode < right.edge.secondNode;
	if (left.edge.id != right.edge.id) return left.edge.id < right.edge.id;
	return left.inputIndex < right.inputIndex;
}

inline bool forestFamilyBetter(
	const ZhangConflictAwareForestFamily& left,
	const ZhangConflictAwareForestFamily& right)
{
	if (left.valid != right.valid) return left.valid;
	if (left.largestMergedComponentCount !=
		right.largestMergedComponentCount)
		return left.largestMergedComponentCount >
			right.largestMergedComponentCount;
	if (left.maximumMergedSize != right.maximumMergedSize)
		return left.maximumMergedSize > right.maximumMergedSize;
	if (left.trueMergeCount != right.trueMergeCount)
		return left.trueMergeCount > right.trueMergeCount;
	if (left.touchedComponentCount != right.touchedComponentCount)
		return left.touchedComponentCount > right.touchedComponentCount;
	if (left.coveredNodeCount != right.coveredNodeCount)
		return left.coveredNodeCount > right.coveredNodeCount;
	if (left.summedMergedSize != right.summedMergedSize)
		return left.summedMergedSize > right.summedMergedSize;
	if (left.covarianceTrace != right.covarianceTrace)
		return left.covarianceTrace < right.covarianceTrace;
	return left.fingerprint < right.fingerprint;
}

inline ZhangConflictAwareForestFamily evaluateForestFamily(
	const std::vector<NormalizedForestEdge>& edges,
	const std::vector<std::size_t>& selected)
{
	ZhangConflictAwareForestFamily result;
	if (selected.size() < 2 || selected.size() > 3) return result;
	std::map<int, int> parent;
	std::map<int, int> componentSize;
	auto addComponent = [&](int component)
	{
		parent.emplace(component, component);
		componentSize.emplace(component, 1);
	};
	for (const auto selectedIndex : selected)
	{
		if (selectedIndex >= edges.size()) return result;
		addComponent(edges[selectedIndex].edge.firstComponent);
		addComponent(edges[selectedIndex].edge.secondComponent);
	}
	auto root = [&](int component)
	{
		int value = component;
		while (parent.at(value) != value) value = parent.at(value);
		while (parent.at(component) != component)
		{
			const int next = parent.at(component);
			parent[component] = value;
			component = next;
		}
		return value;
	};
	std::set<int> touchedComponents;
	std::set<int> coveredNodes;
	for (const auto selectedIndex : selected)
	{
		const auto& normalized = edges[selectedIndex];
		const auto& edge = normalized.edge;
		int first = root(edge.firstComponent);
		int second = root(edge.secondComponent);
		// A component cycle is not a forest and cannot claim another merge.
		if (first == second) return ZhangConflictAwareForestFamily{};
		if (componentSize.at(first) < componentSize.at(second))
			std::swap(first, second);
		parent[second] = first;
		componentSize[first] += componentSize.at(second);
		result.trueMergeCount++;
		result.maximumMergedSize = std::max(
			result.maximumMergedSize, edge.mergedSize);
		result.summedMergedSize += edge.mergedSize;
		result.covarianceTrace += edge.covarianceTrace;
		result.inputEdgeIndices.push_back(normalized.inputIndex);
		result.edgeFingerprints.push_back(normalized.fingerprint);
		touchedComponents.insert(edge.firstComponent);
		touchedComponents.insert(edge.secondComponent);
		coveredNodes.insert(edge.firstNode);
		coveredNodes.insert(edge.secondNode);
	}
	for (const auto& item : parent)
		result.largestMergedComponentCount = std::max(
			result.largestMergedComponentCount,
			componentSize.at(root(item.first)));
	result.touchedComponentCount = static_cast<int>(touchedComponents.size());
	result.coveredNodeCount = static_cast<int>(coveredNodes.size());
	std::sort(result.inputEdgeIndices.begin(), result.inputEdgeIndices.end());
	std::sort(result.edgeFingerprints.begin(), result.edgeFingerprints.end());
	std::ostringstream fingerprint;
	bool firstFingerprint = true;
	for (const auto& edge : result.edgeFingerprints)
	{
		if (!firstFingerprint) fingerprint << ';';
		fingerprint << edge;
		firstFingerprint = false;
	}
	result.fingerprint = fingerprint.str();
	result.valid = result.trueMergeCount == static_cast<int>(selected.size()) &&
		std::isfinite(result.covarianceTrace);
	return result;
}
} // namespace zhang_conflict_aware_detail

/**
 * Generate at most eight deterministic two/three-edge component forests.
 *
 * The construction uses a bounded expansion beam, rejects component cycles,
 * and deduplicates both input edges and complete families.  One candidate from
 * each available family size is retained before filling the remaining slots,
 * so a difficult three-edge ILS still has a two-edge fallback.
 */
inline std::vector<ZhangConflictAwareForestFamily>
zhangGenerateConflictAwareForestFamilies(
	const std::vector<ZhangConflictAwareForestEdge>& inputEdges,
	std::size_t requestedMaximumFamilies = 8,
	std::size_t beamWidth = 64)
{
	using namespace zhang_conflict_aware_detail;
	const std::size_t maximumFamilies =
		std::min<std::size_t>(8, requestedMaximumFamilies);
	if (maximumFamilies == 0 || beamWidth == 0) return {};

	std::map<std::tuple<int, int, int, int>, NormalizedForestEdge> uniqueEdges;
	for (std::size_t inputIndex = 0; inputIndex < inputEdges.size(); inputIndex++)
	{
		auto edge = inputEdges[inputIndex];
		if (edge.firstNode < 0 || edge.secondNode < 0 ||
			edge.firstNode == edge.secondNode || edge.firstComponent < 0 ||
			edge.secondComponent < 0 ||
			edge.firstComponent == edge.secondComponent ||
			edge.mergedSize <= 0 || !std::isfinite(edge.covarianceTrace) ||
			edge.covarianceTrace < 0)
			continue;
		if (edge.secondNode < edge.firstNode)
		{
			std::swap(edge.firstNode, edge.secondNode);
			std::swap(edge.firstComponent, edge.secondComponent);
		}
		const auto firstComponent = std::min(
			edge.firstComponent, edge.secondComponent);
		const auto secondComponent = std::max(
			edge.firstComponent, edge.secondComponent);
		const auto key = std::make_tuple(
			edge.firstNode, edge.secondNode, firstComponent, secondComponent);
		NormalizedForestEdge normalized{
			edge, inputIndex, forestEdgeFingerprint(edge)};
		auto found = uniqueEdges.find(key);
		if (found == uniqueEdges.end() || forestEdgeBetter(
			normalized, found->second))
			uniqueEdges[key] = std::move(normalized);
	}
	std::vector<NormalizedForestEdge> edges;
	for (auto& item : uniqueEdges)
		edges.push_back(std::move(item.second));
	std::sort(edges.begin(), edges.end(), forestEdgeBetter);
	if (edges.size() < 2) return {};

	struct PartialFamily
	{
		std::vector<std::size_t> selected;
		ZhangConflictAwareForestFamily score;
	};
	std::vector<PartialFamily> frontier(1);
	std::vector<ZhangConflictAwareForestFamily> generated;
	for (std::size_t depth = 1; depth <= 3 && !frontier.empty(); depth++)
	{
		std::map<std::string, PartialFamily> unique;
		for (const auto& parent : frontier)
		{
			const std::size_t begin = parent.selected.empty()
				? 0 : parent.selected.back() + 1;
			for (std::size_t edge = begin; edge < edges.size(); edge++)
			{
				auto selected = parent.selected;
				selected.push_back(edge);
				ZhangConflictAwareForestFamily score;
				if (selected.size() >= 2)
				{
					score = evaluateForestFamily(edges, selected);
					if (!score.valid) continue;
				}
				else
				{
					// A one-edge partial branch cannot yet be emitted, but its
					// deterministic key keeps the expansion beam bounded.
					score.valid = true;
					score.fingerprint = edges[edge].fingerprint;
					score.trueMergeCount = 1;
					score.largestMergedComponentCount = 2;
					score.touchedComponentCount = 2;
					score.coveredNodeCount = 2;
					score.maximumMergedSize = edges[edge].edge.mergedSize;
					score.summedMergedSize = edges[edge].edge.mergedSize;
					score.covarianceTrace =
						edges[edge].edge.covarianceTrace;
				}
				unique.try_emplace(score.fingerprint,
					PartialFamily{std::move(selected), std::move(score)});
			}
		}
		frontier.clear();
		for (auto& item : unique)
			frontier.push_back(std::move(item.second));
		std::sort(frontier.begin(), frontier.end(),
			[](const auto& left, const auto& right)
			{
				return forestFamilyBetter(left.score, right.score);
			});
		if (frontier.size() > beamWidth) frontier.resize(beamWidth);
		if (depth >= 2)
			for (const auto& family : frontier)
				generated.push_back(family.score);
	}
	std::map<std::string, ZhangConflictAwareForestFamily> uniqueFamilies;
	for (auto& family : generated)
		uniqueFamilies.try_emplace(family.fingerprint, std::move(family));
	generated.clear();
	for (auto& item : uniqueFamilies)
		generated.push_back(std::move(item.second));
	std::sort(generated.begin(), generated.end(), forestFamilyBetter);

	std::vector<ZhangConflictAwareForestFamily> selected;
	std::set<std::string> selectedFingerprints;
	auto retain = [&](const ZhangConflictAwareForestFamily& family)
	{
		if (selected.size() >= maximumFamilies ||
			!selectedFingerprints.insert(family.fingerprint).second)
			return;
		selected.push_back(family);
	};
	if (maximumFamilies >= 2)
	{
		for (const auto& family : generated)
			if (family.inputEdgeIndices.size() == 3)
			{
				retain(family);
				break;
			}
		for (const auto& family : generated)
			if (family.inputEdgeIndices.size() == 2)
			{
				retain(family);
				break;
			}
	}
	for (const auto& family : generated) retain(family);
	std::sort(selected.begin(), selected.end(), forestFamilyBetter);
	return selected;
}

/** Exact post-evaluation score; all fields are supplied by the caller's exact
 * affine-union, graph and rank audits. */
struct ZhangConflictAwareProductBranchScore
{
	bool valid = false;
	bool safe = false;
	int exactLargestComponentSize = 0;
	int exactSatelliteCount = 0;
	int exactGraphRank = 0;
	int exactLatticeRank = 0;
	int exactComponentCount = 0;
	int temporalOverlap = 0;
	double failureProbability = 1;
	std::string fingerprint;
};

struct ZhangConflictAwareProductBranchSelection
{
	bool valid = false;
	bool replacedBaseline = false;
	int selectedCandidateIndex = -1;
	ZhangConflictAwareProductBranchScore selected;
	std::string reason = "NO_SAFE_BRANCH";
};

inline bool zhangConflictAwareProductBranchScoreWellFormed(
	const ZhangConflictAwareProductBranchScore& score)
{
	return score.valid && score.safe && score.exactLargestComponentSize >= 0 &&
		score.exactSatelliteCount >= 0 && score.exactGraphRank >= 0 &&
		score.exactLatticeRank >= 0 && score.exactComponentCount >= 0 &&
		score.temporalOverlap >= 0 &&
		std::isfinite(score.failureProbability) &&
		score.failureProbability >= 0 && score.failureProbability <= 1;
}

/** An existing empty baseline remains a legal fail-closed retained result. */
inline bool zhangConflictAwareProductBaselineRetainable(
	const ZhangConflictAwareProductBranchScore& score)
{
	return zhangConflictAwareProductBranchScoreWellFormed(score);
}

/** A new branch is adoptable only if its exact product graph is non-degenerate. */
inline bool zhangConflictAwareProductCandidateAdmissible(
	const ZhangConflictAwareProductBranchScore& score)
{
	return zhangConflictAwareProductBranchScoreWellFormed(score) &&
		score.exactSatelliteCount >= 2 && score.exactGraphRank > 0 &&
		score.exactLatticeRank > 0;
}

/** Backward-compatible name denotes validity for candidate adoption. */
inline bool zhangConflictAwareProductBranchScoreValid(
	const ZhangConflictAwareProductBranchScore& score)
{
	return zhangConflictAwareProductCandidateAdmissible(score);
}

inline bool zhangConflictAwareProductBranchNonRegressive(
	const ZhangConflictAwareProductBranchScore& baseline,
	const ZhangConflictAwareProductBranchScore& candidate)
{
	if (!zhangConflictAwareProductCandidateAdmissible(candidate)) return false;
	if (!zhangConflictAwareProductBaselineRetainable(baseline)) return true;
	return candidate.exactLargestComponentSize >=
			baseline.exactLargestComponentSize &&
		candidate.exactSatelliteCount >= baseline.exactSatelliteCount &&
		candidate.exactGraphRank >= baseline.exactGraphRank &&
		candidate.exactLatticeRank >= baseline.exactLatticeRank;
}

inline bool zhangConflictAwareProductBranchStrictlyBetter(
	const ZhangConflictAwareProductBranchScore& candidate,
	const ZhangConflictAwareProductBranchScore& reference)
{
	if (candidate.exactLargestComponentSize !=
		reference.exactLargestComponentSize)
		return candidate.exactLargestComponentSize >
			reference.exactLargestComponentSize;
	if (candidate.exactSatelliteCount != reference.exactSatelliteCount)
		return candidate.exactSatelliteCount > reference.exactSatelliteCount;
	if (candidate.exactGraphRank != reference.exactGraphRank)
		return candidate.exactGraphRank > reference.exactGraphRank;
	if (candidate.exactLatticeRank != reference.exactLatticeRank)
		return candidate.exactLatticeRank > reference.exactLatticeRank;
	if (candidate.temporalOverlap != reference.temporalOverlap)
		return candidate.temporalOverlap > reference.temporalOverlap;
	if (candidate.exactComponentCount != reference.exactComponentCount)
		return candidate.exactComponentCount < reference.exactComponentCount;
	return candidate.failureProbability < reference.failureProbability;
}

inline ZhangConflictAwareProductBranchSelection
zhangSelectConflictAwareProductBranch(
	const ZhangConflictAwareProductBranchScore& baseline,
	const std::vector<ZhangConflictAwareProductBranchScore>& candidates)
{
	ZhangConflictAwareProductBranchSelection result;
	const bool baselineValid =
		zhangConflictAwareProductBaselineRetainable(baseline);
	if (baselineValid)
	{
		result.valid = true;
		result.selected = baseline;
		result.reason = "BASELINE_RETAINED";
	}
	for (int index = 0; index < static_cast<int>(candidates.size()); index++)
	{
		const auto& candidate = candidates[index];
		if (!zhangConflictAwareProductBranchNonRegressive(
			baseline, candidate)) continue;
		const bool betterThanBaseline = !baselineValid ||
			zhangConflictAwareProductBranchStrictlyBetter(candidate, baseline);
		if (!betterThanBaseline) continue;
		const bool betterThanSelected = !result.valid ||
			zhangConflictAwareProductBranchStrictlyBetter(
				candidate, result.selected) ||
			(!zhangConflictAwareProductBranchStrictlyBetter(
				result.selected, candidate) &&
			 candidate.fingerprint < result.selected.fingerprint);
		if (!betterThanSelected) continue;
		result.valid = true;
		result.replacedBaseline = baselineValid;
		result.selectedCandidateIndex = index;
		result.selected = candidate;
		result.reason = baselineValid
			? "SAFE_NON_REGRESSIVE_REPLACEMENT"
			: "SAFE_CANDIDATE_SELECTED";
	}
	return result;
}
