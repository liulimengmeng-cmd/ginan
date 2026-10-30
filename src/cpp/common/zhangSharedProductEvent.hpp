#pragma once

#include "pea/zhangReference.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "common/zhangHandoffSearchBudget.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// This key groups work, not authority. A connected FLOAT/product graph can
// contain several independently certified integer components; publication
// must still check each product's exact current certificate and risk closure.
struct ZhangSharedProductEventKey
{
    std::string runtimeId;
    long int epochGpstSeconds = 0;
    E_Sys system = E_Sys::NONE;
    std::string graphComponentId;
    std::string posteriorId;
    std::string chartId;
    std::string fullJointTransportId;
    int oldBackendGeneration = -1;
    int newBackendGeneration = -1;
    int oldChartGeneration = -1;
    int newChartGeneration = -1;

    bool operator<(const ZhangSharedProductEventKey& other) const
    {
        return std::tie(runtimeId, epochGpstSeconds, system,
                   graphComponentId, posteriorId, chartId,
                   fullJointTransportId, oldBackendGeneration,
                   newBackendGeneration, oldChartGeneration,
                   newChartGeneration) <
            std::tie(other.runtimeId, other.epochGpstSeconds,
                   other.system, other.graphComponentId,
                   other.posteriorId, other.chartId,
                   other.fullJointTransportId,
                   other.oldBackendGeneration,
                   other.newBackendGeneration,
                   other.oldChartGeneration,
                   other.newChartGeneration);
    }
};

inline bool zhangMakeSharedProductEventKey(
    const ZhangGraphProductEventReceipt& receipt,
    const std::string& runtimeId,
    long int epochGpstSeconds,
    const std::string& graphComponentId,
    const std::string& posteriorId,
    const std::string& chartId,
    int currentBackendGeneration,
    int currentChartGeneration,
    ZhangSharedProductEventKey& key,
    std::string& reason)
{
    key = {};
    reason = "NOT_EVALUATED";
    if (runtimeId.empty() || graphComponentId.empty() ||
        posteriorId.empty() || chartId.empty() ||
        receipt.system == E_Sys::NONE ||
        receipt.fullJointTransportId.empty() ||
        receipt.oldBackendGeneration < 0 ||
        receipt.newBackendGeneration != currentBackendGeneration ||
        receipt.oldChartGeneration < 0 ||
        receipt.newChartGeneration != currentChartGeneration ||
        !receipt.oldProductProjectionCaptured ||
        !receipt.actualWriterProjectionBound ||
        !receipt.representableSameTimeMapProven ||
        !receipt.fullJointMomentTransportProven ||
        !(receipt.sameTimeMapProven ||
          receipt.conditionalProductMapProven))
    {
        reason = "SHARED_EVENT_RECEIPT_OR_GENERATION_INVALID";
        return false;
    }
    key = {runtimeId, epochGpstSeconds, receipt.system,
        graphComponentId, posteriorId, chartId,
        receipt.fullJointTransportId,
        receipt.oldBackendGeneration,
        receipt.newBackendGeneration,
        receipt.oldChartGeneration,
        receipt.newChartGeneration};
    reason = "SHARED_EVENT_IDENTITY_BOUND";
    return true;
}

struct ZhangSharedProductTargetKey
{
    E_ObsCode observable = E_ObsCode::NONE;
    std::vector<std::tuple<ZhangGraphEdge, int, ZhangExactInteger>> terms;
    ZhangExactInteger affineOffsetCycles = 0;

    bool operator<(const ZhangSharedProductTargetKey& other) const
    {
        if (observable != other.observable)
            return observable < other.observable;
        if (terms != other.terms)
            return terms < other.terms;
        return affineOffsetCycles < other.affineOffsetCycles;
    }
};

inline bool zhangCanonicalSharedProductTarget(
    const ZhangProductIntegerTransition& transition,
    ZhangSharedProductTargetKey& key,
    std::string& reason)
{
    key = {};
    reason = "NOT_EVALUATED";
    if (!transition.valid ||
        transition.observable == E_ObsCode::NONE ||
        transition.physicalEdges.size() !=
            transition.physicalArcVersions.size() ||
        transition.physicalEdges.size() !=
            transition.coefficients.size())
    {
        reason = "SHARED_EVENT_TARGET_INVALID";
        return false;
    }
    std::map<std::pair<ZhangGraphEdge, int>, ZhangExactInteger> combined;
    for (std::size_t index = 0;
         index < transition.coefficients.size(); index++)
    {
        const int version = transition.physicalArcVersions[index];
        if (version < 0)
        {
            reason = "SHARED_EVENT_SIGNAL_VERSION_INVALID";
            return false;
        }
        combined[{transition.physicalEdges[index], version}] +=
            transition.coefficients[index];
    }
    key.observable = transition.observable;
    key.affineOffsetCycles = transition.affineOffsetCycles;
    for (const auto& [identity, coefficient] : combined)
        if (coefficient != 0)
            key.terms.emplace_back(
                identity.first, identity.second, coefficient);
    reason = "SHARED_EVENT_EXACT_TARGET_CANONICAL";
    return true;
}

struct ZhangSharedProductEventSlot
{
    std::string id;
    SatSys satellite;
    SatSys referenceSatellite;
    // Search admission can be a one-step event increment while the query is
    // the complete segment-B-to-current target. They are not interchangeable.
    ZhangProductIntegerTransition eventIncrement;
};

// Scheduling diagnostics are not integer evidence. In particular, completing
// this phase (or reusing a failed plan) does not certify a target. The caller
// owns the canonical paired-plan cache and the actual pair-search budget.
struct ZhangSharedProductEventSearchVisit
{
    std::string phase = "NOT_VISITED";
    std::string reason = "NOT_SEARCHED_EVENT_NOT_VISITED";
    bool attempted = false;
    bool reused = false;
    bool stopEvent = false;
    std::string unvisitedReason = "NOT_SEARCHED_EVENT_STOPPED";
    bool currentQueryProven = false;
    bool searchVisited = false;
    // Only an actual incomplete attempt's representative may request a
    // retry. An alias reusing that outcome must not add another queue entry.
    bool retryPending = false;
    bool incomplete = false;
    std::size_t attemptCount = 0;
    std::size_t retryCount = 0;
    std::string lastAttemptReason;
};

struct ZhangSharedProductEventSearchReport
{
    std::map<std::string, ZhangSharedProductEventSearchVisit> visits;
    std::size_t initialQueries = 0;
    std::size_t searchCallbacks = 0;
    std::size_t attempted = 0;
    std::size_t reused = 0;
    std::size_t unvisited = 0;
    std::size_t actualAttempts = 0;
    std::size_t retries = 0;
    std::size_t retryCallbacks = 0;
    std::size_t incomplete = 0;
    std::size_t pendingRetries = 0;

    void refreshCounts()
    {
        attempted = reused = unvisited = actualAttempts = retries =
            incomplete = pendingRetries = 0;
        for (const auto& [id, visit] : visits)
        {
            attempted += visit.attempted;
            reused += visit.reused;
            unvisited += !visit.searchVisited && !visit.currentQueryProven;
            actualAttempts += visit.attemptCount;
            retries += visit.retryCount;
            incomplete += visit.incomplete && !visit.currentQueryProven;
            pendingRetries += visit.retryPending && !visit.currentQueryProven;
        }
    }
};

// Every unique complete versioned target is initially queried once. Slots
// with the same exact target receive the same final proof only after all
// candidate searches, from a stable current-evidence generation. This is a
// real work-sharing boundary, not a replacement for the current integer gate.
class ZhangSharedProductEventBatch
{
public:
    struct Target
    {
        ZhangProductIntegerTransition transition;
        std::vector<ZhangSharedProductEventSlot> slots;
    };

    explicit ZhangSharedProductEventBatch(ZhangSharedProductEventKey identity)
        : identity_(std::move(identity)) {}

    const ZhangSharedProductEventKey& identity() const { return identity_; }
    std::size_t uniqueTargetCount() const { return targets_.size(); }
    std::size_t slotCount() const { return slotKeys_.size(); }

    // Phase one runs before any final product posterior is constructed.
    // Complete the initial batch queries before entering any pair search:
    // expensive shared-domain preparation must not consume a pair deadline.
    // The search target is the COMPLETE native-to-segment difference, not a
    // latest-event increment. Slots sharing this single-signal target can
    // have different companion-signal targets. Try each slot until the
    // complete target is proven; the caller deduplicates actual searches by
    // the canonical pair of complete signal targets and evidence generation.
    template <class QUERY, class SEARCH, class IS_PROVEN>
    bool executeSearchPhase(QUERY&& query, SEARCH&& search,
                            IS_PROVEN&& isProven,
                            std::string& reason) const
    {
        ZhangSharedProductEventSearchReport report;
        // Compatibility for the original void callback. It cannot report
        // whether caller-side pair deduplication actually performed a search.
        // Do not infer attempted=true from merely returning from the callback.
        auto visitSearch = [&](const auto& slot, const auto& transition)
        {
            search(slot, transition);
            ZhangSharedProductEventSearchVisit visit;
            visit.phase = "CALLER_SEARCH_CALLBACK";
            visit.reason = "CALLER_SEARCH_CALLBACK_COMPLETED_OUTCOME_UNREPORTED";
            return visit;
        };
        return executeSearchPhase(std::forward<QUERY>(query), visitSearch,
            std::forward<IS_PROVEN>(isProven), reason, report);
    }

    // The structured callback must catch a recoverable per-pair timeout and
    // return its actual failure reason. Only exhausted aggregate search time
    // should set stopEvent. Unexpected exceptions remain exceptions: annotate
    // all unvisited slots and rethrow rather than declaring a completed search.
    // A three-argument callback additionally receives retryRound (zero for
    // every initial visit). It may request bounded retries only for actual
    // incomplete wall-time attempts. The caller still owns pair/generation
    // deduplication and must charge every retry to the same aggregate budget.
    template <class QUERY, class SEARCH, class IS_PROVEN>
    bool executeSearchPhase(QUERY&& query, SEARCH&& search,
                            IS_PROVEN&& isProven,
                            std::string& reason,
                            ZhangSharedProductEventSearchReport& report) const
    {
        report = {};
        reason = "NOT_EVALUATED";
        if (targets_.empty())
        {
            reason = "SHARED_EVENT_NO_TARGETS";
            return false;
        }
        for (const auto& [targetKey, target] : targets_)
            for (const auto& slot : target.slots)
                report.visits.emplace(slot.id, ZhangSharedProductEventSearchVisit{});

        using Proof = decltype(query(targets_.begin()->second.transition));
        std::map<ZhangSharedProductTargetKey, Proof> initialProofs;
        using Pending = std::pair<ZhangSharedProductTargetKey,
            const ZhangSharedProductEventSlot*>;
        std::vector<Pending> pending;
        auto callSearch = [&](const ZhangSharedProductEventSlot& slot,
                              const ZhangProductIntegerTransition& transition,
                              unsigned retryRound)
        {
            if constexpr (requires { search(slot, transition, retryRound); })
                return search(slot, transition, retryRound);
            else
                return search(slot, transition);
        };
        auto recordSearch = [&](ZhangSharedProductEventSearchVisit& recorded,
                                ZhangSharedProductEventSearchVisit next,
                                unsigned retryRound)
        {
            const bool actualAttempt = next.attempted;
            next.attemptCount = recorded.attemptCount + actualAttempt;
            next.retryCount = recorded.retryCount + (retryRound != 0 && actualAttempt);
            if (actualAttempt) next.lastAttemptReason = next.reason;
            else if (next.lastAttemptReason.empty())
                next.lastAttemptReason = recorded.lastAttemptReason;
            next.attempted = recorded.attempted || actualAttempt;
            next.reused = recorded.reused || next.reused;
            next.searchVisited = true;
            // Do not turn a reused alias into another actual-plan retry.
            next.retryPending = next.retryPending && actualAttempt && next.incomplete;
            recorded = std::move(next);
            return recorded.retryPending;
        };
        auto markProven = [&](const Target& target, const Proof& proof)
        {
            std::string proofReason = "CURRENT_QUERY_PROVEN_NO_SEARCH_REQUIRED";
            if constexpr (requires { std::string(proof.reason); })
                if (!std::string(proof.reason).empty()) proofReason = proof.reason;
            for (const auto& slot : target.slots)
            {
                auto& visit = report.visits.at(slot.id);
                visit.currentQueryProven = true;
                visit.retryPending = false;
                // Preserve actual failed/search/reuse history for a visited
                // slot even if a later companion supplies a valid consequence.
                if (visit.searchVisited) continue;
                visit.phase = "CURRENT_QUERY_PROVEN";
                visit.reason = proofReason;
            }
        };
        auto markUnvisited = [&](const std::string& phase,
                                 const std::string& notSearchedReason)
        {
            for (auto& [id, visit] : report.visits)
                if (!visit.searchVisited && !visit.currentQueryProven)
                {
                    visit.phase = phase;
                    visit.reason = notSearchedReason;
                }
            report.refreshCounts();
        };
        try
        {
            for (const auto& [targetKey, target] : targets_)
            {
                ++report.initialQueries;
                Proof current = query(target.transition);
                if (isProven(current)) markProven(target, current);
                initialProofs.emplace(targetKey, std::move(current));
            }
        }
        catch (...)
        {
            reason = "SHARED_EVENT_INITIAL_QUERY_EXCEPTION";
            markUnvisited("INITIAL_QUERY_ABORTED", "NOT_SEARCHED_INITIAL_QUERY_EXCEPTION");
            throw;
        }
        for (const auto& [targetKey, target] : targets_)
        {
            Proof current = initialProofs.at(targetKey);
            // Earlier targets may have admitted constraints after this
            // target's initial query. Requery before an otherwise needless
            // companion search; final authority is still checked in phase two.
            if (report.searchCallbacks != 0)
            {
                try { current = query(target.transition); }
                catch (...)
                {
                    reason = "SHARED_EVENT_CURRENT_QUERY_EXCEPTION";
                    markUnvisited("CURRENT_QUERY_ABORTED", "NOT_SEARCHED_CURRENT_QUERY_EXCEPTION");
                    throw;
                }
            }
            if (isProven(current))
            {
                markProven(target, current);
                continue;
            }
            for (const auto& slot : target.slots)
            {
                auto& visit = report.visits.at(slot.id);
                if (visit.searchVisited) continue;
                visit.currentQueryProven = false;
                visit.phase = "PAIR_SEARCH_NOT_STARTED";
                visit.reason = "NOT_SEARCHED_CURRENT_TARGET_UNPROVEN";
            }
            for (const auto& slot : target.slots)
            {
                if (isProven(current)) break;
                auto& recorded = report.visits.at(slot.id);
                ++report.searchCallbacks;
                try
                {
                    if (recordSearch(recorded, callSearch(slot, target.transition, 0), 0))
                        pending.emplace_back(targetKey, &slot);
                }
                catch (...)
                {
                    recorded.phase = "SEARCH_CALLBACK_EXCEPTION";
                    recorded.reason = "SEARCH_OUTCOME_UNAVAILABLE_CALLBACK_EXCEPTION";
                    recorded.searchVisited = true;
                    reason = "SHARED_EVENT_SEARCH_CALLBACK_EXCEPTION";
                    markUnvisited("EVENT_ABORTED", "NOT_SEARCHED_SEARCH_CALLBACK_EXCEPTION");
                    throw;
                }
                if (recorded.stopEvent)
                {
                    reason = recorded.reason;
                    markUnvisited("PAIR_SEARCH_NOT_STARTED", recorded.unvisitedReason);
                    return true;
                }
                try { current = query(target.transition); }
                catch (...)
                {
                    reason = "SHARED_EVENT_POST_SEARCH_QUERY_EXCEPTION";
                    markUnvisited("POST_SEARCH_QUERY_ABORTED", "NOT_SEARCHED_POST_SEARCH_QUERY_EXCEPTION");
                    throw;
                }
            }
            if (isProven(current)) markProven(target, current);
        }
        // Fair deterministic rounds begin only after ALL first visits.
        // The callback revalidates its canonical pair and current generation;
        // no stale outcome or unused slice acquires mathematical authority.
        for (unsigned round = 1;
             !pending.empty() && round <= ZhangHandoffSearchBudget::maximumRetryRounds;
             ++round)
        {
            std::vector<Pending> nextPending;
            for (const auto& [targetKey, slot] : pending)
            {
                const auto& target = targets_.at(targetKey);
                auto& recorded = report.visits.at(slot->id);
                if (!recorded.retryPending) continue;
                Proof current = initialProofs.at(targetKey);
                try { current = query(target.transition); }
                catch (...)
                {
                    reason = "SHARED_EVENT_RETRY_QUERY_EXCEPTION";
                    markUnvisited("RETRY_QUERY_ABORTED", "NOT_SEARCHED_RETRY_QUERY_EXCEPTION");
                    throw;
                }
                if (isProven(current))
                {
                    markProven(target, current);
                    continue;
                }
                ++report.searchCallbacks;
                ++report.retryCallbacks;
                try
                {
                    if (recordSearch(recorded,
                        callSearch(*slot, target.transition, round), round))
                        nextPending.emplace_back(targetKey, slot);
                }
                catch (...)
                {
                    recorded.phase = "SEARCH_CALLBACK_EXCEPTION";
                    recorded.reason = "SEARCH_OUTCOME_UNAVAILABLE_CALLBACK_EXCEPTION";
                    recorded.searchVisited = true;
                    recorded.retryPending = false;
                    recorded.incomplete = true;
                    reason = "SHARED_EVENT_RETRY_CALLBACK_EXCEPTION";
                    markUnvisited("EVENT_ABORTED", "NOT_SEARCHED_RETRY_CALLBACK_EXCEPTION");
                    throw;
                }
                if (recorded.stopEvent)
                {
                    reason = recorded.reason;
                    markUnvisited("PAIR_SEARCH_NOT_STARTED", recorded.unvisitedReason);
                    return true;
                }
                try { current = query(target.transition); }
                catch (...)
                {
                    reason = "SHARED_EVENT_RETRY_POST_SEARCH_QUERY_EXCEPTION";
                    markUnvisited("POST_SEARCH_QUERY_ABORTED", "NOT_SEARCHED_RETRY_POST_SEARCH_QUERY_EXCEPTION");
                    throw;
                }
                if (isProven(current)) markProven(target, current);
            }
            pending = std::move(nextPending);
        }
        report.refreshCounts();
        reason = report.pendingRetries != 0
            ? "SHARED_EVENT_RETRY_ROUND_LIMIT_INCOMPLETE"
            : "SHARED_EVENT_COMPANION_VARIANTS_TRIED_CALLER_PAIR_DEDUP";
        return true;
    }

    // Phase two is called only after ALL graph-event groups have finished
    // admitting candidates. The caller takes one final global generation
    // snapshot and requires it to remain unchanged across every group.
    template <class PROOF, class QUERY, class GENERATION>
    bool finalizeQueryPhase(QUERY&& query, GENERATION&& generation,
                            std::uint64_t expectedGeneration,
                            std::map<std::string, PROOF>& slotProofs,
                            std::string& reason) const
    {
        slotProofs.clear();
        reason = "NOT_EVALUATED";
        if (targets_.empty())
        {
            reason = "SHARED_EVENT_NO_TARGETS";
            return false;
        }
        if (generation() != expectedGeneration)
        {
            reason = "SHARED_EVENT_EVIDENCE_GENERATION_CHANGED";
            return false;
        }
        std::map<ZhangSharedProductTargetKey, PROOF> finalProofs;
        for (const auto& [targetKey, target] : targets_)
        {
            if (generation() != expectedGeneration)
            {
                reason = "SHARED_EVENT_EVIDENCE_GENERATION_CHANGED";
                return false;
            }
            finalProofs.emplace(targetKey, query(target.transition));
            if (generation() != expectedGeneration)
            {
                reason = "SHARED_EVENT_EVIDENCE_GENERATION_CHANGED";
                return false;
            }
        }
        if (generation() != expectedGeneration)
        {
            reason = "SHARED_EVENT_EVIDENCE_GENERATION_CHANGED";
            return false;
        }
        for (const auto& [targetKey, target] : targets_)
            for (const auto& slot : target.slots)
                slotProofs.emplace(slot.id, finalProofs.at(targetKey));
        reason = "SHARED_EVENT_FINAL_GENERATION_REQUERIED";
        return true;
    }

    bool add(const ZhangSharedProductEventSlot& slot,
             const ZhangProductIntegerTransition& transition,
             std::string& reason)
    {
        reason = "NOT_EVALUATED";
        if (slot.id.empty() ||
            slot.satellite.sys != identity_.system ||
            slot.referenceSatellite.sys != identity_.system ||
            slotKeys_.contains(slot.id))
        {
            reason = "SHARED_EVENT_SLOT_IDENTITY_INVALID";
            return false;
        }
        ZhangSharedProductTargetKey targetKey;
        if (!zhangCanonicalSharedProductTarget(
                transition, targetKey, reason))
            return false;
        auto [target, inserted] = targets_.try_emplace(targetKey);
        if (inserted)
            target->second.transition = transition;
        target->second.slots.push_back(slot);
        slotKeys_.emplace(slot.id, targetKey);
        reason = inserted ? "SHARED_EVENT_UNIQUE_TARGET"
                          : "SHARED_EVENT_TARGET_REUSED";
        return true;
    }

    template <class PROOF, class QUERY, class SEARCH,
              class GENERATION, class IS_PROVEN>
    bool execute(QUERY&& query, SEARCH&& search,
                 GENERATION&& generation, IS_PROVEN&& isProven,
                 std::map<std::string, PROOF>& slotProofs,
                 std::string& reason) const
    {
        slotProofs.clear();
        reason = "NOT_EVALUATED";
        if (targets_.empty())
        {
            reason = "SHARED_EVENT_NO_TARGETS";
            return false;
        }
        auto incrementalSearch = [&](const auto& slot, const auto& transition)
        {
            search(slot, slot.eventIncrement.valid ? slot.eventIncrement : transition);
        };
        if (!executeSearchPhase(query, incrementalSearch, isProven, reason))
            return false;
        const std::uint64_t finalGeneration = generation();
        if (!finalizeQueryPhase<PROOF>(query, generation,
                finalGeneration, slotProofs, reason)) return false;
        reason = "SHARED_EVENT_ALL_TARGETS_REQUERIED";
        return true;
    }

private:
    ZhangSharedProductEventKey identity_;
    std::map<ZhangSharedProductTargetKey, Target> targets_;
    std::map<std::string, ZhangSharedProductTargetKey> slotKeys_;
};
