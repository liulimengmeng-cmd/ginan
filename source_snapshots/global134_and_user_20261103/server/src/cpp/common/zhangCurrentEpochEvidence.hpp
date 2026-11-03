#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "common/zhangIntegerDecisionProof.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangRatioOnly.hpp"
#include "common/zhangProductPhysicalCycleChart.hpp"
#include "common/zhangR47ProductDomain.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "common/zhangCurrentProductBridgePlan.hpp"
#include "common/zhangPhaseTimer.hpp"

struct KFState;

using ZhangCurrentEpochPhysicalRow =
	std::map<std::string, ZhangExactInteger>;

struct ZhangCurrentEpochNumericConstraint
{
	ZhangCoordinateRow rowCycles;
	ZhangExactInteger value;
};

/** Admission-time link between one CURRENT physical equation and the
 * numerical row actually added to the immutable FLOAT conditioning recipe.
 * This is not old-arc storage or historical integer authority. */
struct ZhangCurrentEpochPhysicalNumericBinding
{
	ZhangCurrentEpochPhysicalRow physical;
	ZhangCoordinateRow rowCycles;
	ZhangExactInteger value;
	long int epochGpstSeconds = 0;
	E_Sys system = E_Sys::NONE;
	std::string runtimeId, posteriorId, rootId, chartId, stateOrderId;
};

struct ZhangCurrentEpochDomainIdentity
{
	long int epochGpstSeconds = 0;
	std::string runtimeId, posteriorId, rootId, chartId, physicalColumnsHash;
	E_Sys system = E_Sys::NONE;
	std::uint64_t constraintGeneration = 0;
	bool operator==(const ZhangCurrentEpochDomainIdentity&) const = default;
};

/** Immutable current H*a=z. The feasibility decomposition also recovers one
 * exact origin for later target images, without a square transform, physical
 * kernel or integer lift. In particular,
 * 2*a=6 entails a=3 here, but that is NOT an integer row-combination witness.
 * Product pair authorization must keep using its row-lattice certificate. */
struct ZhangCurrentEpochCompiledDomain
{
	ZhangCurrentEpochDomainIdentity identity;
	std::vector<ZhangCurrentEpochPhysicalRow> acceptedRows;
	ZhangExactVector acceptedValues;
	std::vector<ZhangDecisionProofs> acceptedProofGroups;
	std::vector<std::string> eventOverlapWitnesses;
	std::map<std::string, std::size_t> columns;
	ZhangExactAffineIntegerQuotient affine;
	// Read-only compatibility views: the large immutable HNF is owned once.
	const ZhangExactMatrix& basis = affine.deterministicBasis;
	const ZhangExactVector& basisValues = affine.deterministicValues;
	std::vector<std::size_t> pivots;
	std::shared_ptr<const ZhangDecisionProofs> sharedProofs;
	ZhangDecisionRiskClosure proofClosure;
	bool feasible = false;
	double buildSeconds = 0;
	std::string reason = "CURRENT_DOMAIN_NOT_COMPILED";
	ZhangCurrentEpochCompiledDomain() = default;
	ZhangCurrentEpochCompiledDomain(const ZhangCurrentEpochCompiledDomain&) = delete;
	ZhangCurrentEpochCompiledDomain& operator=(const ZhangCurrentEpochCompiledDomain&) = delete;

	bool consequence(const ZhangCurrentEpochPhysicalRow& target,
		ZhangExactInteger& value) const
	{
		if (!feasible) return false;
		using Rational = boost::multiprecision::cpp_rational;
		std::vector<Rational> remaining(columns.size());
		for (const auto& [id, coefficient] : target)
		{
			if (coefficient == 0) continue;
			const auto column = columns.find(id);
			if (column == columns.end()) return false;
			remaining[column->second] = Rational(coefficient);
		}
		Rational result = 0;
		for (std::size_t row = 0; row < basis.size(); ++row)
		{
			const auto pivot = pivots[row];
			if (remaining[pivot] == 0) continue;
			const Rational factor = remaining[pivot] / Rational(basis[row][pivot]);
			for (std::size_t column = pivot; column < remaining.size(); ++column)
				if (basis[row][column] != 0)
					remaining[column] -= factor * Rational(basis[row][column]);
			result += factor * Rational(basisValues[row]);
		}
		if (std::any_of(remaining.begin(), remaining.end(),
			[](const Rational& coefficient) { return coefficient != 0; }) ||
			denominator(result) != 1) return false;
		value = numerator(result);
		return true;
	}
};

struct ZhangCurrentEpochDomainCache
{
	std::mutex mutex;
	std::shared_ptr<const ZhangCurrentEpochCompiledDomain> domain;
	std::size_t builds = 0, hits = 0;
};

/** Lightweight immutable affine view. Matching the FULL current identity is
 * required at the ownership boundary; pointer identity or epoch alone cannot
 * authorize reuse after an admitted event changes the constraint generation. */
inline std::shared_ptr<const ZhangExactAffineIntegerQuotient>
zhangCurrentEpochAffineDomainView(
	const std::shared_ptr<const ZhangCurrentEpochCompiledDomain>& domain,
	const ZhangCurrentEpochDomainIdentity& expectedIdentity)
{
	if (!domain || !domain->feasible || !domain->affine.valid ||
		domain->identity != expectedIdentity) return {};
	return {domain, &domain->affine};
}

struct ZhangCurrentEpochBridgeProof
{
	bool proven = false;
	// A ratio-only decision is an integer search result, not a calibrated
	// false-fix probability statement.
	bool statisticallyCalibrated = false;
	ZhangExactInteger value = 0;
	ZhangDecisionProofs proofs;
	std::shared_ptr<const ZhangDecisionProofs> sharedProofs;
	std::vector<std::string> proofIds;
	std::string reason = "NO_CURRENT_EPOCH_EVIDENCE";
	int evidenceRows = 0;
	const ZhangDecisionProofs& proofsView() const
	{
		return sharedProofs ? *sharedProofs : proofs;
	}
};

/** One immutable FLOAT posterior/chart may admit several integer decisions,
 * but a later epoch must start with a new object. No held integer store or
 * historical product ledger is consulted by this class. */
struct ZhangCurrentEpochEvidence
{
	std::string runtimeId;
	E_Sys system = E_Sys::NONE;
	long int epochGpstSeconds = 0;
	std::string posteriorId;
	std::string chartId;
	std::string floatStateRoot;
	// Only a successful whole-FLOAT covariance gate may set these fields.
	// Ambiguity-submatrix PSD checks cannot authorize the trusted-source mode
	// of final joint product conditioning.
	bool floatPosteriorPsdValidated = false;
	std::string floatPosteriorPsdRootId;
	std::string floatPosteriorPsdStateOrderId;
	ZhangProductPhysicalCycleChart chart;
	std::vector<ZhangCurrentEpochPhysicalRow> rows;
	ZhangExactVector values;
	std::vector<ZhangDecisionProofs> rowProofs;
	std::vector<ZhangCurrentEpochPhysicalRow> eventOverlapRows;
	ZhangExactVector eventOverlapValues;
	std::vector<ZhangDecisionProofs> eventOverlapProofs;
	std::vector<std::string> eventOverlapWitnesses;
	std::vector<ZhangCurrentEpochNumericConstraint> numericConstraints;
	std::vector<ZhangCurrentEpochPhysicalNumericBinding> physicalNumericBindings;
	std::uint64_t constraintGeneration = 0;
	// Copies made for transactional admission share the immutable domain.
	// A subsequent mutation replaces the slot before modifying the snapshot.
	mutable std::shared_ptr<ZhangCurrentEpochDomainCache> domainCache =
		std::make_shared<ZhangCurrentEpochDomainCache>();
	bool valid = false;

	void invalidateCurrentDomain()
	{
		++constraintGeneration;
		domainCache = std::make_shared<ZhangCurrentEpochDomainCache>();
	}

	bool addNumericConstraint(const ZhangCoordinateRow& rowCycles,
		const ZhangExactInteger& value,
		const ZhangCurrentEpochPhysicalRow& physical = {})
	{
		if (!valid || floatStateRoot.empty() || rowCycles.empty()) return false;
		ZhangCoordinateRow normalized;
		for (const auto& [key, coefficient] : rowCycles)
		{
			if (!std::isfinite(coefficient)) return false;
			if (coefficient != 0) normalized.emplace(key, coefficient);
		}
		if (normalized.empty()) return false;
		bool numericExists = false;
		for (const auto& constraint : numericConstraints)
			if (constraint.rowCycles == normalized)
			{
				if (constraint.value != value) return false;
				numericExists = true;
				break;
			}
		bool bindingExists = false;
		if (!physical.empty())
		{
			if (floatPosteriorPsdStateOrderId.empty()) return false;
			bool physicalExists = false;
			for (const auto& [source, sourceValues] :
				{std::pair{&rows, &values}, std::pair{&eventOverlapRows, &eventOverlapValues}})
			{
				if (source->size() != sourceValues->size()) return false;
				for (std::size_t r = 0; r < source->size(); ++r)
					if ((*source)[r] == physical)
					{
						if ((*sourceValues)[r] != value) return false;
						physicalExists = true;
					}
			}
			if (!physicalExists) return false;
			for (const auto& binding : physicalNumericBindings)
				if (binding.physical == physical)
				{
					if (binding.rowCycles != normalized || binding.value != value ||
						binding.epochGpstSeconds != epochGpstSeconds || binding.system != system ||
						binding.runtimeId != runtimeId || binding.posteriorId != posteriorId ||
						binding.rootId != floatStateRoot || binding.chartId != chartId ||
						binding.stateOrderId != floatPosteriorPsdStateOrderId) return false;
					bindingExists = true;
					break;
				}
		}
		// Check all conflicts before mutating either recipe. Production callers
		// still build a provisional evidence copy before the registry commit.
		if (!numericExists) numericConstraints.push_back({normalized, value});
		if (!physical.empty() && !bindingExists)
			physicalNumericBindings.push_back({physical, normalized, value,
				epochGpstSeconds, system, runtimeId, posteriorId, floatStateRoot,
				chartId, floatPosteriorPsdStateOrderId});
		return true;
	}

	/** Required before a missing-arc rewrite. Every raw H row (therefore its
	 * HNF basis) must be linked to an actual same-root numerical condition.
	 * Old callers without bindings remain valid for pure consequence queries,
	 * but cannot authorize this stronger numerical representation operation. */
	bool validatePhysicalNumericBindings(std::string& reason) const
	{
		reason = "CURRENT_PHYSICAL_NUMERIC_BINDINGS_NOT_EVALUATED";
		if (!valid || floatStateRoot.empty() || floatPosteriorPsdStateOrderId.empty() ||
			rows.size() != values.size() || eventOverlapRows.size() != eventOverlapValues.size())
		{
			reason = "CURRENT_PHYSICAL_NUMERIC_BINDING_SCOPE_INVALID";
			return false;
		}
		std::map<ZhangCoordinateRow, ZhangExactInteger> numerical;
		for (const auto& constraint : numericConstraints)
		{
			zhangExactWorkCheckpoint("CURRENT_PHYSICAL_NUMERIC_RECIPE_AUDIT");
			if (constraint.rowCycles.empty())
			{
				reason = "CURRENT_PHYSICAL_NUMERIC_RECIPE_EMPTY";
				return false;
			}
			for (const auto& [key, coefficient] : constraint.rowCycles)
				if (!std::isfinite(coefficient) || coefficient == 0)
				{
					reason = "CURRENT_PHYSICAL_NUMERIC_RECIPE_INVALID";
					return false;
				}
			const auto [entry, inserted] = numerical.emplace(constraint.rowCycles, constraint.value);
			if (!inserted && entry->second != constraint.value)
			{
				reason = "CURRENT_PHYSICAL_NUMERIC_RECIPE_CONFLICT";
				return false;
			}
		}
		std::map<ZhangCurrentEpochPhysicalRow, const ZhangCurrentEpochPhysicalNumericBinding*> linked;
		for (const auto& binding : physicalNumericBindings)
		{
			zhangExactWorkCheckpoint("CURRENT_PHYSICAL_NUMERIC_BINDING_AUDIT");
			if (binding.physical.empty() || binding.epochGpstSeconds != epochGpstSeconds ||
				binding.system != system || binding.runtimeId != runtimeId ||
				binding.posteriorId != posteriorId || binding.rootId != floatStateRoot ||
				binding.chartId != chartId || binding.stateOrderId != floatPosteriorPsdStateOrderId)
			{
				reason = "CURRENT_PHYSICAL_NUMERIC_BINDING_SCOPE_MISMATCH";
				return false;
			}
			const auto numeric = numerical.find(binding.rowCycles);
			if (numeric == numerical.end() || numeric->second != binding.value)
			{
				reason = "CURRENT_PHYSICAL_NUMERIC_BINDING_RECIPE_MISSING_OR_RHS_MISMATCH";
				return false;
			}
			const auto [entry, inserted] = linked.emplace(binding.physical, &binding);
			if (!inserted && (entry->second->value != binding.value ||
				entry->second->rowCycles != binding.rowCycles))
			{
				reason = "CURRENT_PHYSICAL_NUMERIC_BINDING_CONFLICT";
				return false;
			}
		}
		for (const auto& [source, sourceValues] :
			{std::pair{&rows, &values}, std::pair{&eventOverlapRows, &eventOverlapValues}})
			for (std::size_t r = 0; r < source->size(); ++r)
			{
				zhangExactWorkCheckpoint("CURRENT_PHYSICAL_NUMERIC_RAW_DOMAIN_AUDIT");
				const auto binding = linked.find((*source)[r]);
				if (binding == linked.end() || binding->second->value != (*sourceValues)[r])
				{
					reason = "CURRENT_PHYSICAL_NUMERIC_RAW_ROW_UNBOUND_OR_RHS_MISMATCH";
					return false;
				}
			}
		reason = "CURRENT_PHYSICAL_NUMERIC_BINDINGS_VERIFIED";
		return true;
	}

	std::string physicalColumnsHash() const
	{
		// Stable content identity, independent of object addresses. Full chart
		// content is checked by bind; hashes are diagnostic/cache keys only.
		std::uint64_t hash = 14695981039346656037ULL;
		auto append = [&](const std::string& text)
		{
			const auto length = std::to_string(text.size()) + ":";
			for (const unsigned char byte : length) { hash ^= byte; hash *= 1099511628211ULL; }
			for (const unsigned char byte : text) { hash ^= byte; hash *= 1099511628211ULL; }
		};
		append(std::to_string(chart.columns));
		for (const auto& signal : chart.cataloguedSignals) append(signal);
		for (const auto& id : chart.activePhysicalIdentities) append(id);
		for (const auto& [column, id] : chart.chordIdentities)
		{ append(std::to_string(column)); append(id); }
		for (const auto& [column, row] : chart.expansions)
		{
			append(std::to_string(column));
			for (const auto& [id, coefficient] : row)
			{ append(id); append(coefficient.convert_to<std::string>()); }
		}
		std::ostringstream text;
		text << std::hex << std::setw(16) << std::setfill('0') << hash;
		return text.str();
	}

	ZhangCurrentEpochDomainIdentity currentDomainIdentity() const
	{
		return {epochGpstSeconds, runtimeId, posteriorId, floatStateRoot, chartId,
			physicalColumnsHash(), system, constraintGeneration};
	}

	std::shared_ptr<const ZhangCurrentEpochCompiledDomain> compiledCurrentDomain() const
	{
		const auto identity = currentDomainIdentity();
		std::lock_guard<std::mutex> lock(domainCache->mutex);
		const auto& prior = domainCache->domain;
		auto sameRows = [&](const ZhangCurrentEpochCompiledDomain& domain)
		{
			if (domain.acceptedRows.size() != rows.size() + eventOverlapRows.size() ||
				domain.acceptedValues.size() != values.size() + eventOverlapValues.size() ||
				domain.acceptedProofGroups.size() != rowProofs.size() + eventOverlapProofs.size())
				return false;
			return std::equal(rows.begin(), rows.end(), domain.acceptedRows.begin()) &&
				std::equal(eventOverlapRows.begin(), eventOverlapRows.end(),
					domain.acceptedRows.begin() + rows.size()) &&
				std::equal(values.begin(), values.end(), domain.acceptedValues.begin()) &&
				std::equal(eventOverlapValues.begin(), eventOverlapValues.end(),
					domain.acceptedValues.begin() + values.size()) &&
				std::equal(rowProofs.begin(), rowProofs.end(),
					domain.acceptedProofGroups.begin()) &&
				std::equal(eventOverlapProofs.begin(), eventOverlapProofs.end(),
					domain.acceptedProofGroups.begin() + rowProofs.size()) &&
				domain.eventOverlapWitnesses == eventOverlapWitnesses;
		};
		if (prior && prior->identity == identity && sameRows(*prior))
		{ ++domainCache->hits; return prior; }
		auto domain = std::make_shared<ZhangCurrentEpochCompiledDomain>();
		domain->identity = identity;
		const auto start = std::chrono::steady_clock::now();
		domain->acceptedRows = rows;
		domain->acceptedRows.insert(domain->acceptedRows.end(),
			eventOverlapRows.begin(), eventOverlapRows.end());
		domain->acceptedValues = values;
		domain->acceptedValues.insert(domain->acceptedValues.end(),
			eventOverlapValues.begin(), eventOverlapValues.end());
		domain->acceptedProofGroups = rowProofs;
		domain->acceptedProofGroups.insert(domain->acceptedProofGroups.end(),
			eventOverlapProofs.begin(), eventOverlapProofs.end());
		domain->eventOverlapWitnesses = eventOverlapWitnesses;
		ZhangDecisionProofs proofs;
		std::set<const ZhangIntegerDecisionProof*> seen;
		for (const auto* source : {&rowProofs, &eventOverlapProofs})
			for (const auto& group : *source)
				for (const auto& proof : group)
					if (seen.insert(proof.get()).second) proofs.push_back(proof);
		domain->sharedProofs = std::make_shared<const ZhangDecisionProofs>(std::move(proofs));
		domain->proofClosure = zhangDecisionRiskClosure(*domain->sharedProofs);
		if (valid && domain->acceptedRows.size() == domain->acceptedValues.size() &&
			rows.size() == rowProofs.size() &&
			eventOverlapRows.size() == eventOverlapProofs.size() &&
			eventOverlapRows.size() == eventOverlapWitnesses.size())
		{
			for (const auto& row : domain->acceptedRows)
				for (const auto& [id, coefficient] : row)
					if (coefficient != 0 && !domain->columns.contains(id))
						domain->columns.emplace(id, domain->columns.size());
			ZhangExactMatrix h(domain->acceptedRows.size(),
				ZhangExactVector(domain->columns.size()));
			for (std::size_t row = 0; row < domain->acceptedRows.size(); ++row)
				for (const auto& [id, coefficient] : domain->acceptedRows[row])
					if (coefficient != 0) h[row][domain->columns.at(id)] = coefficient;
			auto hnf = zhangExactRowHermiteNormalForm(std::move(h), domain->acceptedValues);
			ZhangIntegerLatticeMembership origin;
			if (hnf.consistent)
			{
				ZhangR48ExactTimer timer("R51_CURRENT_DOMAIN_INTEGER_ORIGIN",
					hnf.basis, domain->columns.size());
				if (hnf.basis.empty())
				{
					origin.contained = true;
					origin.combination.resize(domain->columns.size());
				}
				else
				{
					ZhangExactMatrix generators(domain->columns.size(),
						ZhangExactVector(hnf.basis.size()));
					for (std::size_t c = 0; c < domain->columns.size(); ++c)
					for (std::size_t r = 0; r < hnf.basis.size(); ++r)
					{
						zhangExactWorkCheckpoint("CURRENT_DOMAIN_ORIGIN_TRANSPOSE");
						generators[c][r] = hnf.basis[r][c];
					}
					// Identical Smith operations and origin as the old product
					// quotient, but reuse this one feasibility decomposition.
					origin = zhangIntegerRowLatticeContains(generators, hnf.values, true, true);
				}
			}
			if (hnf.consistent && origin.contained &&
				origin.combination.size() == domain->columns.size() &&
				zhangExactMatrixTimesColumn(hnf.basis, origin.combination) == hnf.values)
			{
				domain->affine.ambientDimension = domain->columns.size();
				domain->affine.deterministicRank = hnf.basis.size();
				domain->affine.quotientRank = domain->columns.size() - hnf.basis.size();
				domain->affine.deterministicBasis = std::move(hnf.basis);
				domain->affine.deterministicValues = std::move(hnf.values);
				domain->affine.particularSolution = std::move(origin.combination);
				domain->affine.smithInvariants = std::move(origin.smithInvariants);
				domain->affine.valid = true;
				domain->affine.failureReason = "NONE";
				for (const auto& row : domain->basis)
					domain->pivots.push_back(std::distance(row.begin(),
						std::find_if(row.begin(), row.end(),
							[](const auto& coefficient) { return coefficient != 0; })));
				domain->feasible = true;
				domain->reason = "CURRENT_INTEGER_DOMAIN_FEASIBLE";
			}
			else domain->reason = "CURRENT_INTEGER_DOMAIN_INFEASIBLE";
		}
		else domain->reason = "CURRENT_INTEGER_DOMAIN_INVALID";
		domain->buildSeconds = std::chrono::duration<double>(
			std::chrono::steady_clock::now() - start).count();
		domainCache->domain = domain;
		++domainCache->builds;
		return domain;
	}

	bool bind(
		const std::string& runtime,
		E_Sys sys,
		long int epoch,
		const std::string& posterior,
		const std::string& chartIdentity,
		const ZhangProductPhysicalCycleChart& frozenChart,
		const std::string& frozenFloatStateRoot = {})
	{
		if (runtime.empty() || sys == E_Sys::NONE || epoch <= 0 ||
			posterior.empty() || chartIdentity.empty() ||
			frozenChart.columns <= 0 || frozenChart.expansions.empty())
			return false;
		if (valid)
			return runtimeId == runtime && system == sys &&
				epochGpstSeconds == epoch && posteriorId == posterior &&
				chartId == chartIdentity &&
				floatStateRoot == frozenFloatStateRoot &&
				chart.columns == frozenChart.columns &&
				chart.activePhysicalIdentities == frozenChart.activePhysicalIdentities &&
				chart.cataloguedSignals == frozenChart.cataloguedSignals &&
				chart.chordIdentities == frozenChart.chordIdentities &&
				chart.expansions == frozenChart.expansions;
		runtimeId = runtime;
		system = sys;
		epochGpstSeconds = epoch;
		posteriorId = posterior;
		chartId = chartIdentity;
		floatStateRoot = frozenFloatStateRoot;
		chart = frozenChart;
		valid = true;
		return true;
	}

	bool matches(long int epoch, const std::string& posterior,
		const std::string& chartIdentity) const
	{
		return valid && epochGpstSeconds == epoch &&
			posteriorId == posterior && chartId == chartIdentity;
	}

	bool add(const ZhangCurrentEpochPhysicalRow& physical,
		const ZhangExactInteger& value, const ZhangDecisionProofs& proofs)
	{
		if (!valid || physical.empty() || proofs.empty() ||
			rows.size() != values.size() || rows.size() != rowProofs.size())
			return false;
		const auto risk = zhangDecisionRiskClosure(proofs);
		if (!risk.valid ||
			zhangRatioStatisticalReject(risk.bound > 1e-3)) return false;
		const std::string currentSource = runtimeId + "|PPP_COMMITTED|E" +
			std::to_string(epochGpstSeconds);
		for (const auto& [id, proof] : risk.atoms)
			if (proof->observationProvenance != currentSource &&
				proof->observationProvenance.rfind(currentSource + "|", 0) != 0)
				return false;
		ZhangExactVector projected;
		if (!chart.project(physical, projected)) return false;
		for (std::size_t i = 0; i < rows.size(); ++i)
		{
			if (rows[i] != physical) continue;
			if (values[i] != value) return false;
			auto merged = zhangMergeDecisionProofs(rowProofs[i], proofs);
			const auto mergedRisk = zhangDecisionRiskClosure(merged);
			if (!mergedRisk.valid ||
				zhangRatioStatisticalReject(mergedRisk.bound > 1e-3))
				return false;
			if (merged == rowProofs[i]) return true;
			rowProofs[i] = std::move(merged);
			invalidateCurrentDomain();
			return true;
		}
		rows.push_back(physical);
		values.push_back(value);
		rowProofs.push_back(proofs);
		invalidateCurrentDomain();
		return true;
	}

	bool addEventOverlap(const ZhangCurrentEpochPhysicalRow& physical,
		const ZhangExactInteger& value, const ZhangDecisionProofs& proofs,
		const std::string& graphEventWitness)
	{
		if (!valid || floatStateRoot.empty() || graphEventWitness.empty() ||
			physical.empty() || proofs.empty() ||
			eventOverlapRows.size() != eventOverlapValues.size() ||
			eventOverlapRows.size() != eventOverlapProofs.size() ||
			eventOverlapRows.size() != eventOverlapWitnesses.size())
			return false;
		const auto risk = zhangDecisionRiskClosure(proofs);
		if (!risk.valid ||
			zhangRatioStatisticalReject(risk.bound > 1e-3)) return false;
		const std::string currentSource = runtimeId + "|PPP_COMMITTED|E" +
			std::to_string(epochGpstSeconds);
		for (const auto& [id, proof] : risk.atoms)
			if (proof->observationProvenance != currentSource &&
				proof->observationProvenance.rfind(currentSource + "|", 0) != 0)
				return false;
		for (std::size_t i = 0; i < eventOverlapRows.size(); ++i)
		{
			if (eventOverlapRows[i] != physical) continue;
			return eventOverlapValues[i] == value &&
				eventOverlapWitnesses[i] == graphEventWitness;
		}
		eventOverlapRows.push_back(physical);
		eventOverlapValues.push_back(value);
		eventOverlapProofs.push_back(proofs);
		eventOverlapWitnesses.push_back(graphEventWitness);
		invalidateCurrentDomain();
		return true;
	}

	ZhangCurrentEpochBridgeProof query(
		E_ObsCode observable,
		long int epoch,
		const std::string& posterior,
		const std::string& chartIdentity,
		const ZhangProductIntegerTransition& difference,
		bool materializeProofs = true) const
	{
		ZhangCurrentEpochBridgeProof out;
		out.evidenceRows = static_cast<int>(rows.size() + eventOverlapRows.size());
		if (!matches(epoch, posterior, chartIdentity))
		{
			out.reason = "CURRENT_EPOCH_CONTEXT_MISMATCH";
			return out;
		}
		if (observable == E_ObsCode::NONE || !difference.valid ||
			difference.physicalEdges.size() != difference.coefficients.size() ||
			difference.physicalEdges.size() !=
				difference.physicalArcVersions.size())
		{
			out.reason = "CURRENT_EPOCH_BRIDGE_INVALID";
			return out;
		}
		// An empty difference is a structural identity check performed by the
		// writer. It cannot create an AR certificate or rescue an old relation.
		if (difference.coefficients.empty())
		{
			out.reason = "EMPTY_BRIDGE_HAS_NO_INTEGER_EVIDENCE";
			return out;
		}
		ZhangCurrentEpochPhysicalRow target;
		for (std::size_t i = 0; i < difference.coefficients.size(); ++i)
		{
			if (difference.coefficients[i] == 0) continue;
			const auto& edge = difference.physicalEdges[i];
			if (edge.receiver.empty() || edge.satellite.sys != system ||
				difference.physicalArcVersions[i] < 0)
			{
				out.reason = "BRIDGE_ARC_IDENTITY_INVALID";
				return out;
			}
			const std::string id = enum_to_string(observable) +
				std::string("|") + edge.receiver + "|" + edge.satellite.id() +
				"|V" + std::to_string(difference.physicalArcVersions[i]);
			target[id] += difference.coefficients[i];
		}
		for (auto it = target.begin(); it != target.end();)
			if (it->second == 0) it = target.erase(it);
			else ++it;
		if (target.empty())
		{
			out.reason = "ZERO_BRIDGE_HAS_NO_INTEGER_EVIDENCE";
			return out;
		}
		ZhangExactVector projected;
		const bool inCurrentChart = chart.project(target, projected);
		if (!inCurrentChart && eventOverlapRows.empty())
		{
			out.reason = "BRIDGE_OUTSIDE_FROZEN_FLOAT_CHART";
			return out;
		}
		if (rows.size() != values.size() || rows.size() != rowProofs.size() ||
			eventOverlapRows.size() != eventOverlapValues.size() ||
			eventOverlapRows.size() != eventOverlapProofs.size() ||
			eventOverlapRows.size() != eventOverlapWitnesses.size() ||
			(rows.empty() && eventOverlapRows.empty()))
		{
			out.reason = "NO_CURRENT_ACCEPTED_INTEGER_ROWS";
			return out;
		}
		ZhangDecisionProofs usedProofs;
		ZhangExactInteger physicalValue = 0;
		bool direct = false;
		const auto domain = compiledCurrentDomain();
		if (!domain->feasible)
		{
			out.reason = "BRIDGE_NOT_EXACT_INTEGER_CONSEQUENCE";
			return out;
		}
		const auto& acceptedRows = domain->acceptedRows;
		const auto& acceptedValues = domain->acceptedValues;
		for (std::size_t i = 0; i < acceptedRows.size(); ++i)
		{
			const auto& acceptedProofs = i < rowProofs.size()
				? rowProofs[i] : eventOverlapProofs[i - rowProofs.size()];
			if (acceptedRows[i] == target)
			{
				physicalValue = acceptedValues[i];
				usedProofs = acceptedProofs;
				direct = true;
				break;
			}
			bool negative = acceptedRows[i].size() == target.size();
			if (negative)
				for (const auto& [id, coefficient] : acceptedRows[i])
				{
					const auto term = target.find(id);
					if (term == target.end() || term->second != -coefficient)
					{ negative = false; break; }
				}
			if (negative)
			{
				physicalValue = -acceptedValues[i];
				usedProofs = acceptedProofs;
				direct = true;
				break;
			}
		}
		if (!direct)
		{
			if (!domain->consequence(target, physicalValue))
			{
				out.reason = "BRIDGE_NOT_EXACT_INTEGER_CONSEQUENCE";
				return out;
			}
			// Conservatively include every accepted current decision if exact
			// elimination used multiple rows. Merge/validate once per immutable
			// domain rather than repeatedly copying all roots for every target.
			out.sharedProofs = domain->sharedProofs;
		}
		ZhangDecisionRiskClosure directRisk;
		if (direct) directRisk = zhangDecisionRiskClosure(usedProofs);
		const auto& risk = direct ? directRisk : domain->proofClosure;
		if ((direct ? usedProofs.empty() : domain->sharedProofs->empty()) || !risk.valid ||
			zhangRatioStatisticalReject(risk.bound > 1e-3))
		{
			out.reason = "CURRENT_EPOCH_PROOF_RISK_INVALID";
			return out;
		}
		out.proven = true;
		out.statisticallyCalibrated = !zhangRatioOnly();
		out.value = physicalValue + difference.affineOffsetCycles;
		if (direct)
		{
			out.sharedProofs = std::make_shared<const ZhangDecisionProofs>(std::move(usedProofs));
			if (materializeProofs) out.proofs = *out.sharedProofs;
		}
		else if (materializeProofs) out.proofs = *out.sharedProofs;
		for (const auto& [id, proof] : risk.atoms) out.proofIds.push_back(id);
		out.reason = direct
			? (inCurrentChart ? "EXACT_CURRENT_ROW" :
				"EXACT_CURRENT_EVENT_OVERLAP_ROW")
			: (inCurrentChart ? "EXACT_CURRENT_AFFINE_INTEGER_CONSEQUENCE"
				: "EXACT_CURRENT_EVENT_OVERLAP_CONSEQUENCE");
		return out;
	}
};

/** Same-parent event admission without cloning the complete numerical recipe.
 * The caller still performs exact residual-image admission and whole-FLOAT
 * verification. This class only stages the already verified rows and applies
 * the existing provenance/conflict gates to base UNION delta. No candidate is
 * authoritative before commit. A transaction is local to one registry owner
 * and one immutable FLOAT/chart/constraint generation; it cannot cross epochs.
 *
 * Registry mutation is single-writer. All supported mutations append rows or
 * replace the domain generation/cache; in-place editing of admitted records is
 * forbidden. Reserves may change capacity on cancellation, never evidence.
 */
class ZhangCurrentEpochEventTransaction
{
	const ZhangCurrentEpochEvidence* parent_;
	ZhangCurrentEpochDomainIdentity identity_;
	std::string stateOrder_, psdRoot_;
	bool psdValidated_, valid_;
	std::shared_ptr<ZhangCurrentEpochDomainCache> cache_;
	std::vector<std::size_t> sizes_;
	ZhangCurrentEpochEvidence delta_;
	std::shared_ptr<ZhangCurrentEpochDomainCache> nextCache_;
	bool prepared_ = false, committed_ = false;
	std::uint64_t committedGeneration_ = 0;

	static std::vector<std::size_t> sizes(const ZhangCurrentEpochEvidence& e)
	{
		return {e.rows.size(), e.values.size(), e.rowProofs.size(),
			e.eventOverlapRows.size(), e.eventOverlapValues.size(),
			e.eventOverlapProofs.size(), e.eventOverlapWitnesses.size(),
			e.numericConstraints.size(), e.physicalNumericBindings.size()};
	}
	bool parentMatches(const ZhangCurrentEpochEvidence& e) const
	{
		return &e == parent_ && e.valid == valid_ &&
			e.floatPosteriorPsdValidated == psdValidated_ &&
			e.floatPosteriorPsdRootId == psdRoot_ &&
			e.floatPosteriorPsdStateOrderId == stateOrder_ &&
			e.currentDomainIdentity() == identity_ && e.domainCache == cache_ &&
			sizes(e) == sizes_;
	}
	void setScope(ZhangCurrentEpochEvidence& e) const
	{
		e.valid = valid_; e.runtimeId = identity_.runtimeId;
		e.system = identity_.system; e.epochGpstSeconds = identity_.epochGpstSeconds;
		e.posteriorId = identity_.posteriorId; e.floatStateRoot = identity_.rootId;
		e.chartId = identity_.chartId;
		e.floatPosteriorPsdValidated = psdValidated_;
		e.floatPosteriorPsdRootId = psdRoot_;
		e.floatPosteriorPsdStateOrderId = stateOrder_;
	}
	bool baseConflicts(const ZhangCurrentEpochPhysicalRow& physical,
		const ZhangExactInteger& value, const std::string& witness,
		const ZhangCoordinateRow& numeric, bool& eventExists,
		bool& numericExists, bool& bindingExists, std::string& reason) const
	{
		const auto& base = *parent_;
		for (const auto& [rows, values] : {std::pair{&base.rows, &base.values},
			std::pair{&base.eventOverlapRows, &base.eventOverlapValues}})
		{
			if (rows->size() != values->size())
			{ reason = "EVENT_DELTA_BASE_SHAPE_INVALID"; return true; }
			for (std::size_t i = 0; i < rows->size(); ++i)
			{
				zhangExactWorkCheckpoint("EVENT_DELTA_PHYSICAL_CONFLICT_SCAN");
				if ((*rows)[i] != physical) continue;
				if ((*values)[i] != value)
				{ reason = "EVENT_DELTA_PHYSICAL_RHS_CONFLICT"; return true; }
				if (rows == &base.eventOverlapRows)
				{
					if (i >= base.eventOverlapWitnesses.size() ||
						base.eventOverlapWitnesses[i] != witness)
					{ reason = "EVENT_DELTA_EVENT_WITNESS_CONFLICT"; return true; }
					eventExists = true;
				}
			}
		}
		for (const auto& item : base.numericConstraints)
		{
			zhangExactWorkCheckpoint("EVENT_DELTA_NUMERIC_CONFLICT_SCAN");
			if (item.rowCycles != numeric) continue;
			if (item.value != value)
			{ reason = "EVENT_DELTA_NUMERIC_RHS_CONFLICT"; return true; }
			numericExists = true;
		}
		for (const auto& item : base.physicalNumericBindings)
		{
			zhangExactWorkCheckpoint("EVENT_DELTA_BINDING_CONFLICT_SCAN");
			if (item.physical != physical) continue;
			if (item.rowCycles != numeric || item.value != value ||
				item.epochGpstSeconds != identity_.epochGpstSeconds ||
				item.system != identity_.system || item.runtimeId != identity_.runtimeId ||
				item.posteriorId != identity_.posteriorId || item.rootId != identity_.rootId ||
				item.chartId != identity_.chartId || item.stateOrderId != stateOrder_)
			{ reason = "EVENT_DELTA_PHYSICAL_NUMERIC_BINDING_CONFLICT"; return true; }
			bindingExists = true;
		}
		return false;
	}
	// Capacities are reserved before publication. With noexcept element moves,
	// no allocation, budget callback or validation can split the logical commit.
	template<class T> static void append(std::vector<T>& into, std::vector<T>& from) noexcept
	{
		static_assert(std::is_nothrow_move_constructible_v<T>);
		for (auto& item : from) into.emplace_back(std::move(item));
	}
public:
	explicit ZhangCurrentEpochEventTransaction(const ZhangCurrentEpochEvidence& base)
		: parent_(&base), identity_(base.currentDomainIdentity()),
		  stateOrder_(base.floatPosteriorPsdStateOrderId), psdRoot_(base.floatPosteriorPsdRootId),
		  psdValidated_(base.floatPosteriorPsdValidated), valid_(base.valid),
		  cache_(base.domainCache), sizes_(sizes(base))
	{ setScope(delta_); }
	bool committed() const { return committed_; }
	bool matches(const ZhangCurrentEpochEvidence& base) const
	{ return !committed_ && parentMatches(base); }
	bool stage(const ZhangCurrentEpochPhysicalRow& physical,
		const ZhangExactInteger& value, const ZhangDecisionProofs& proofs,
		const std::string& witness, const ZhangCoordinateRow& numeric,
		std::string& reason)
	{
		if (committed_ || prepared_ || !parentMatches(*parent_))
		{ reason = "EVENT_DELTA_PARENT_CHANGED_OR_CLOSED"; return false; }
		zhangExactLibraryBoundary("EVENT_DELTA_STAGE_ENTRY");
		ZhangCoordinateRow normalized;
		for (const auto& [key, coefficient] : numeric)
		{
			if (!std::isfinite(coefficient))
			{ reason = "EVENT_DELTA_NUMERIC_NONFINITE"; return false; }
			if (coefficient != 0) normalized.emplace(key, coefficient);
		}
		bool eventExists = false, numericExists = false, bindingExists = false;
		if (baseConflicts(physical, value, witness, normalized,
			eventExists, numericExists, bindingExists, reason)) return false;
		// Only copy the small event delta, never the base chart, numerical
		// recipe or thousands of repeated state-order strings in old bindings.
		zhangExactLibraryBoundary("EVENT_DELTA_PRIVATE_COPY_BEGIN");
		auto candidate = delta_;
		zhangExactLibraryBoundary("EVENT_DELTA_PRIVATE_COPY_END");
		if (!candidate.addEventOverlap(physical, value, proofs, witness) ||
			!candidate.addNumericConstraint(normalized, value, physical))
		{ reason = "EVENT_DELTA_EXISTING_ADMISSION_GATE_REJECTED"; return false; }
		zhangExactLibraryBoundary("EVENT_DELTA_STAGE_VALIDATED");
		delta_ = std::move(candidate);
		reason = "EVENT_DELTA_STAGED_PRIVATE"; return true;
	}
	bool prepare(std::string& reason)
	{
		if (committed_ || !parentMatches(*parent_))
		{ reason = "EVENT_DELTA_PARENT_CHANGED_OR_CLOSED"; return false; }
		if (prepared_) return true;
		if (parent_->rows.size() != parent_->rowProofs.size() ||
			parent_->eventOverlapRows.size() != parent_->eventOverlapProofs.size() ||
			parent_->eventOverlapRows.size() != parent_->eventOverlapWitnesses.size())
		{ reason = "EVENT_DELTA_PARENT_PROOF_SHAPE_INVALID"; return false; }
		// Base and delta can refer to one proof ID through different parents.
		// Check their joint closure before publication, retaining pointer sharing.
		ZhangDecisionProofs unionProofs;
		for (const auto* groups : {&parent_->rowProofs, &parent_->eventOverlapProofs,
			static_cast<const std::vector<ZhangDecisionProofs>*>(&delta_.eventOverlapProofs)})
			for (const auto& group : *groups)
			{
				zhangExactWorkCheckpoint("EVENT_DELTA_UNION_PROOF_GROUP");
				unionProofs.insert(unionProofs.end(), group.begin(), group.end());
			}
		zhangExactLibraryBoundary("EVENT_DELTA_UNION_PROOF_CLOSURE_BEGIN");
		const auto closure = [&]()
		{
			ZhangPhaseTimer timer("BRIDGE_DELTA_UNION_PROOF_CLOSURE");
			return zhangDecisionRiskClosure(unionProofs);
		}();
		zhangExactLibraryBoundary("EVENT_DELTA_UNION_PROOF_CLOSURE_END");
		if (!closure.valid || zhangRatioStatisticalReject(closure.bound > 1e-3))
		{ reason = "EVENT_DELTA_BASE_UNION_DELTA_PROOF_CONFLICT"; return false; }
		ZhangCurrentEpochEvidence ready; setScope(ready);
		for (std::size_t i = 0; i < delta_.eventOverlapRows.size(); ++i)
		{
			zhangExactWorkCheckpoint("EVENT_DELTA_PREPARE");
			const auto& row = delta_.eventOverlapRows[i];
			const auto& value = delta_.eventOverlapValues[i];
			const auto binding = std::find_if(delta_.physicalNumericBindings.begin(),
				delta_.physicalNumericBindings.end(), [&](const auto& b) { return b.physical == row; });
			if (binding == delta_.physicalNumericBindings.end())
			{ reason = "EVENT_DELTA_NUMERIC_BINDING_MISSING"; return false; }
			bool eventExists = false, numericExists = false, bindingExists = false;
			if (baseConflicts(row, value, delta_.eventOverlapWitnesses[i], binding->rowCycles,
				eventExists, numericExists, bindingExists, reason)) return false;
			if (!eventExists)
			{
				ready.eventOverlapRows.push_back(row); ready.eventOverlapValues.push_back(value);
				ready.eventOverlapProofs.push_back(delta_.eventOverlapProofs[i]);
				ready.eventOverlapWitnesses.push_back(delta_.eventOverlapWitnesses[i]);
			}
			if (!numericExists && std::none_of(ready.numericConstraints.begin(),
				ready.numericConstraints.end(), [&](const auto& c) { return c.rowCycles == binding->rowCycles; }))
				ready.numericConstraints.push_back({binding->rowCycles, value});
			if (!bindingExists) ready.physicalNumericBindings.push_back(*binding);
		}
		nextCache_ = std::make_shared<ZhangCurrentEpochDomainCache>();
		delta_ = std::move(ready); prepared_ = true;
		reason = "EVENT_DELTA_READY_PRIVATE"; return true;
	}
	bool commit(ZhangCurrentEpochEvidence& target, std::string& reason)
	{
		if (committed_)
		{
			auto currentIdentity = target.currentDomainIdentity();
			currentIdentity.constraintGeneration = identity_.constraintGeneration;
			if (&target != parent_ || target.constraintGeneration != committedGeneration_ ||
				currentIdentity != identity_ || target.valid != valid_ ||
				target.floatPosteriorPsdValidated != psdValidated_ ||
				target.floatPosteriorPsdRootId != psdRoot_ ||
				target.floatPosteriorPsdStateOrderId != stateOrder_)
			{ reason = "EVENT_DELTA_COMMITTED_OWNER_OR_GENERATION_CHANGED"; return false; }
			reason = "EVENT_DELTA_ALREADY_COMMITTED"; return true;
		}
		if (!parentMatches(target))
		{ reason = "EVENT_DELTA_COMMIT_PARENT_CHANGED"; return false; }
		if (!prepare(reason)) return false;
		const auto count = delta_.eventOverlapRows.size();
		const bool changed = count || !delta_.numericConstraints.empty() || !delta_.physicalNumericBindings.empty();
		const auto increment = changed ? (count ? count : 1) : 0;
		if (increment > std::numeric_limits<std::uint64_t>::max() - target.constraintGeneration)
		{ reason = "EVENT_DELTA_CONSTRAINT_GENERATION_OVERFLOW"; return false; }
		zhangExactLibraryBoundary("EVENT_DELTA_BEFORE_RESERVE");
		target.eventOverlapRows.reserve(target.eventOverlapRows.size() + delta_.eventOverlapRows.size());
		zhangExactLibraryBoundary("EVENT_DELTA_RESERVE_ROWS");
		target.eventOverlapValues.reserve(target.eventOverlapValues.size() + delta_.eventOverlapValues.size());
		target.eventOverlapProofs.reserve(target.eventOverlapProofs.size() + delta_.eventOverlapProofs.size());
		target.eventOverlapWitnesses.reserve(target.eventOverlapWitnesses.size() + delta_.eventOverlapWitnesses.size());
		zhangExactLibraryBoundary("EVENT_DELTA_RESERVE_EVENTS");
		target.numericConstraints.reserve(target.numericConstraints.size() + delta_.numericConstraints.size());
		target.physicalNumericBindings.reserve(target.physicalNumericBindings.size() + delta_.physicalNumericBindings.size());
		zhangExactLibraryBoundary("EVENT_DELTA_BEFORE_ATOMIC_PUBLICATION");
		if (!parentMatches(target))
		{ reason = "EVENT_DELTA_COMMIT_PARENT_CHANGED"; return false; }
		// No throwing operation below until all six containers and generation
		// have been published. Keep the old per-added-row generation semantics.
		append(target.eventOverlapRows, delta_.eventOverlapRows);
		append(target.eventOverlapValues, delta_.eventOverlapValues);
		append(target.eventOverlapProofs, delta_.eventOverlapProofs);
		append(target.eventOverlapWitnesses, delta_.eventOverlapWitnesses);
		append(target.numericConstraints, delta_.numericConstraints);
		append(target.physicalNumericBindings, delta_.physicalNumericBindings);
		if (changed)
		{
			target.constraintGeneration += increment;
			target.domainCache = std::move(nextCache_);
		}
		committedGeneration_ = target.constraintGeneration; committed_ = true;
		reason = "EVENT_DELTA_ATOMICALLY_COMMITTED"; return true;
	}
};

using ZhangCurrentEpochEvidenceKey = std::pair<std::string, E_Sys>;

inline std::map<ZhangCurrentEpochEvidenceKey, ZhangCurrentEpochEvidence>&
zhangCurrentEpochEvidenceRegistry()
{
	static std::map<ZhangCurrentEpochEvidenceKey, ZhangCurrentEpochEvidence> registry;
	return registry;
}

inline bool zhangBindCurrentEpochEvidence(
	const std::string& runtimeId, E_Sys system, long int epoch,
	const std::string& posteriorId, const std::string& chartId,
	const ZhangProductPhysicalCycleChart& chart,
	const std::string& floatStateRoot = {})
{
	auto& registry = zhangCurrentEpochEvidenceRegistry();
	auto& evidence = registry[{runtimeId, system}];
	return evidence.bind(runtimeId, system, epoch, posteriorId, chartId,
		chart, floatStateRoot);
}

inline ZhangCurrentEpochBridgeProof zhangQueryCurrentEpochBridge(
	const std::string& runtimeId, E_Sys system, E_ObsCode observable,
	long int epoch, const std::string& posteriorId,
	const std::string& chartId,
	const ZhangProductIntegerTransition& difference,
	bool materializeProofs = true)
{
	const auto& registry = zhangCurrentEpochEvidenceRegistry();
	const auto found = registry.find({runtimeId, system});
	if (found == registry.end()) return {};
	return found->second.query(observable, epoch, posteriorId, chartId, difference,
		materializeProofs);
}

// Implemented by ppp_ambres.cpp; graph update and AR entrance both call it.
// Repeated calls for the same runtime/epoch are idempotent.
void zhangBeginEpochOnlyIntegerTransaction(const KFState& state);

// The writer may ask this only for an actual changed product functional that
// lacks a current-chart consequence. The implementation verifies a current
// graph-event witness and searches the joint WL/L1 FLOAT posterior with the
// existing LAMBDA ratio gate; it never rounds a scalar FLOAT into a proof.
bool zhangTryAdmitCurrentEventOverlap(
	std::ostream& trace, const KFState& state, E_Sys system,
	E_ObsCode observable, const SatSys& satellite,
	const SatSys& referenceSatellite, long int epoch,
	const std::string& posteriorId, const std::string& chartId,
	const ZhangProductIntegerTransition& difference,
	double ratioThreshold, std::string& reason,
	const ZhangCurrentProductBridgePlan* completeTargetPlan=nullptr);
