#include <boost/test/unit_test.hpp>

#include "common/zhangSharedProductEvent.hpp"
#include "common/zhangExactWorkBudget.hpp"
#include "common/zhangHandoffSearchBudget.hpp"

#include <set>
#include <stdexcept>

namespace
{
ZhangProductIntegerTransition sharedEventTestTarget(const std::string& receiver,
    int prn, E_ObsCode code = E_ObsCode::L1C)
{
    ZhangProductIntegerTransition target;
    target.valid = true;
    target.observable = code;
    target.physicalEdges = {{receiver, SatSys(E_Sys::GPS, prn)}};
    target.physicalArcVersions = {2};
    target.coefficients = {1};
    return target;
}

ZhangSharedProductEventBatch sharedEventTestBatch()
{
    ZhangSharedProductEventKey key;
    key.system = E_Sys::GPS;
    return ZhangSharedProductEventBatch(key);
}
}

BOOST_AUTO_TEST_CASE(shared_product_event_groups_one_exact_target_once)
{
    ZhangGraphProductEventReceipt receipt;
    receipt.system = E_Sys::GPS;
    receipt.fullJointTransportId = "same-graph-transaction";
    receipt.oldBackendGeneration = 5;
    receipt.newBackendGeneration = 6;
    receipt.oldChartGeneration = 9;
    receipt.newChartGeneration = 10;
    receipt.oldProductProjectionCaptured = true;
    receipt.actualWriterProjectionBound = true;
    receipt.representableSameTimeMapProven = true;
    receipt.fullJointMomentTransportProven = true;
    receipt.conditionalProductMapProven = true;
    ZhangSharedProductEventKey key;
    std::string reason;
    BOOST_REQUIRE(zhangMakeSharedProductEventKey(
        receipt, "runtime", 1000, "common-component", "posterior",
        "chart", 6, 10, key, reason));
    ZhangSharedProductEventBatch batch(key);
    const ZhangGraphEdge a01{"A", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge b13{"B", SatSys(E_Sys::GPS, 13)};
    ZhangProductIntegerTransition repeated;
    repeated.valid = true;
    repeated.observable = E_ObsCode::L1C;
    repeated.physicalEdges = {b13, a01, b13};
    repeated.physicalArcVersions = {1, 0, 1};
    repeated.coefficients = {1, 2, -1};
    ZhangProductIntegerTransition canonical;
    canonical.valid = true;
    canonical.observable = E_ObsCode::L1C;
    canonical.physicalEdges = {a01};
    canonical.physicalArcVersions = {0};
    canonical.coefficients = {2};
    BOOST_REQUIRE(batch.add({"G13-L1", SatSys(E_Sys::GPS, 13),
        SatSys(E_Sys::GPS, 1)}, repeated, reason));
    BOOST_REQUIRE(batch.add({"G14-L1", SatSys(E_Sys::GPS, 14),
        SatSys(E_Sys::GPS, 1)}, canonical, reason));
    BOOST_CHECK_EQUAL(batch.slotCount(), 2);
    BOOST_CHECK_EQUAL(batch.uniqueTargetCount(), 1);

    struct Proof { bool proven = false; int value = 0; };
    std::uint64_t generation = 0;
    int queries = 0;
    int searches = 0;
    std::map<std::string, Proof> proofs;
    const bool executed = batch.execute<Proof>(
        [&](const auto&) { ++queries; return Proof{generation != 0, 7}; },
        [&](const auto&, const auto&) { ++searches; ++generation; },
        [&]() { return generation; },
        [](const Proof& proof) { return proof.proven; },
        proofs, reason);
    BOOST_REQUIRE(executed);
    BOOST_CHECK_EQUAL(searches, 1);
    BOOST_CHECK_EQUAL(queries, 3); // initial + post-search + final generation
    BOOST_REQUIRE_EQUAL(proofs.size(), 2);
    BOOST_CHECK(proofs.at("G13-L1").proven);
    BOOST_CHECK(proofs.at("G14-L1").proven);
    BOOST_CHECK_EQUAL(proofs.at("G13-L1").value,
        proofs.at("G14-L1").value);

    // Production uses a two-phase batch: ALL groups search first, then every
    // unique complete target is queried on one immutable final generation.
    generation = 0;
    queries = 0;
    searches = 0;
    proofs.clear();
    BOOST_REQUIRE(batch.executeSearchPhase(
        [&](const auto&) { ++queries; return Proof{generation != 0, 7}; },
        [&](const auto&, const auto&) { ++searches; ++generation; },
        [](const Proof& proof) { return proof.proven; }, reason));
    BOOST_CHECK_EQUAL(searches, 1);
    BOOST_CHECK_EQUAL(queries, 2); // initial + post-search consequence
    const std::uint64_t finalGeneration = generation;
    BOOST_REQUIRE(batch.finalizeQueryPhase<Proof>(
        [&](const auto&) { ++queries; return Proof{true, 7}; },
        [&]() { return generation; }, finalGeneration, proofs, reason));
    BOOST_CHECK_EQUAL(queries, 3);
    BOOST_CHECK_EQUAL(proofs.size(), 2);
}

BOOST_AUTO_TEST_CASE(shared_product_event_tries_distinct_companion_pairs_until_proven)
{
    ZhangSharedProductEventKey key;
    key.system = E_Sys::GPS;
    ZhangSharedProductEventBatch batch(key);
    ZhangProductIntegerTransition completeL1;
    completeL1.valid = true;
    completeL1.observable = E_ObsCode::L1C;
    completeL1.physicalEdges = {{"A", SatSys(E_Sys::GPS, 13)}};
    completeL1.physicalArcVersions = {2};
    completeL1.coefficients = {1};
    std::string reason;
    // All three slots have the same complete L1 target. The first two have
    // the same L2 companion; only the third offers a distinct paired search.
    for (const auto& [id, prn] :
         std::vector<std::pair<std::string, int>>{{"pair-a", 13},
                                                  {"pair-a-duplicate", 14},
                                                  {"pair-b", 15}})
        BOOST_REQUIRE(batch.add({id, SatSys(E_Sys::GPS, prn),
            SatSys(E_Sys::GPS, 1)}, completeL1, reason));
    BOOST_REQUIRE_EQUAL(batch.uniqueTargetCount(), 1);
    struct Proof { bool proven = false; };
    std::set<std::string> attemptedPairs;
    bool completeTargetProven = false;
    int actualSearches = 0;
    int queries = 0;
    BOOST_REQUIRE(batch.executeSearchPhase(
        [&](const auto&) { ++queries; return Proof{completeTargetProven}; },
        [&](const ZhangSharedProductEventSlot& slot, const auto&)
        {
            const std::string pair = slot.id == "pair-b" ? "L1+L2-B"
                                                          : "L1+L2-A";
            if (!attemptedPairs.insert(pair).second) return;
            ++actualSearches;
            if (pair == "L1+L2-B") completeTargetProven = true;
        },
        [](const Proof& proof) { return proof.proven; }, reason));
    BOOST_CHECK_EQUAL(actualSearches, 2);
    BOOST_CHECK_EQUAL(attemptedPairs.size(), 2);
    BOOST_CHECK_EQUAL(queries, 4); // initial and after each slot variant
    BOOST_CHECK(completeTargetProven);
    BOOST_CHECK_EQUAL(reason,
        "SHARED_EVENT_COMPANION_VARIANTS_TRIED_CALLER_PAIR_DEDUP");
}

BOOST_AUTO_TEST_CASE(shared_product_event_rejects_wrong_generation_and_mutating_final_query)
{
    ZhangGraphProductEventReceipt receipt;
    receipt.system = E_Sys::GPS;
    receipt.fullJointTransportId = "one-event";
    receipt.oldBackendGeneration = 4;
    receipt.newBackendGeneration = 5;
    receipt.oldChartGeneration = 7;
    receipt.newChartGeneration = 8;
    receipt.oldProductProjectionCaptured = true;
    receipt.actualWriterProjectionBound = true;
    receipt.representableSameTimeMapProven = true;
    receipt.fullJointMomentTransportProven = true;
    receipt.sameTimeMapProven = true;
    ZhangSharedProductEventKey key;
    std::string reason;
    BOOST_CHECK(!zhangMakeSharedProductEventKey(
        receipt, "runtime", 1000, "component", "posterior",
        "chart", 6, 8, key, reason));
    BOOST_REQUIRE(zhangMakeSharedProductEventKey(
        receipt, "runtime", 1000, "component", "posterior",
        "chart", 5, 8, key, reason));
    ZhangSharedProductEventBatch batch(key);
    ZhangProductIntegerTransition transition;
    transition.valid = true;
    transition.observable = E_ObsCode::L2W;
    transition.physicalEdges = {{"A", SatSys(E_Sys::GPS, 13)}};
    transition.physicalArcVersions = {2};
    transition.coefficients = {1};
    BOOST_REQUIRE(batch.add({"G13-L2", SatSys(E_Sys::GPS, 13),
        SatSys(E_Sys::GPS, 1)}, transition, reason));
    struct Proof { bool proven = true; };
    int calls = 0;
    std::uint64_t generation = 0;
    std::map<std::string, Proof> proofs;
    BOOST_CHECK(!batch.execute<Proof>(
        [&](const auto&) {
            if (++calls == 2) ++generation;
            return Proof{};
        },
        [&](const auto&, const auto&) {},
        [&]() { return generation; },
        [](const Proof& proof) { return proof.proven; },
        proofs, reason));
    BOOST_CHECK_EQUAL(reason, "SHARED_EVENT_EVIDENCE_GENERATION_CHANGED");
    BOOST_CHECK(proofs.empty());
}

BOOST_AUTO_TEST_CASE(shared_product_event_prequeries_do_not_charge_pair_budget)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    BOOST_REQUIRE(batch.add({"A", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("A", 13), reason));
    BOOST_REQUIRE(batch.add({"B", SatSys(E_Sys::GPS, 14), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("B", 14), reason));
    struct Proof { bool proven = false; };
    ZhangHandoffSearchBudget budget(2000);
    auto now = ZhangHandoffSearchBudget::Clock::time_point{};
    int queries = 0, actualSearches = 0;
    bool targetBProven = false;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase(
        [&](const auto& target)
        {
            ++queries;
            // Deterministic slow shared-domain work: no wall-clock sleeps.
            // Advancing the clock alone must not charge the search allowance.
            now += std::chrono::milliseconds(5000);
            return Proof{target.physicalEdges.front().receiver == "B" && targetBProven};
        },
        [&](const auto& slot, const auto&)
        {
            BOOST_CHECK_EQUAL(report.initialQueries, 2);
            BOOST_CHECK(queries >= 2); // all unique targets queried first
            BOOST_REQUIRE(!budget.exhausted());
            const auto started = now;
            BOOST_CHECK(budget.pairDeadline(started) == started + std::chrono::milliseconds(1000));
            ++actualSearches;
            now += std::chrono::milliseconds(slot.id == "A" ? 1000 : 20);
            budget.charge(started, now);
            targetBProven = slot.id == "B";
            ZhangSharedProductEventSearchVisit visit;
            visit.phase = "PAIR_SEARCH";
            visit.reason = slot.id == "A" ? "CURRENT_EVENT_COOPERATIVE_BUDGET_EXHAUSTED"
                                         : "CURRENT_EVENT_OVERLAP_ADMITTED";
            visit.attempted = true;
            return visit;
        },
        [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(actualSearches, 2);
    BOOST_CHECK_EQUAL(report.attempted, 2);
    BOOST_CHECK_EQUAL(report.unvisited, 0);
    BOOST_CHECK_CLOSE(budget.consumedSeconds(), 1.02, 1e-8);
    BOOST_CHECK_EQUAL(report.visits.at("A").reason,
        "CURRENT_EVENT_COOPERATIVE_BUDGET_EXHAUSTED");
    BOOST_CHECK(report.visits.at("B").currentQueryProven);
}

BOOST_AUTO_TEST_CASE(shared_product_event_pair_timeout_does_not_stop_next_companion)
{
    auto batch = sharedEventTestBatch();
    const auto target = sharedEventTestTarget("A", 13);
    std::string reason;
    BOOST_REQUIRE(batch.add({"timeout", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)}, target, reason));
    BOOST_REQUIRE(batch.add({"healthy", SatSys(E_Sys::GPS, 14), SatSys(E_Sys::GPS, 1)}, target, reason));
    struct Proof { bool proven = false; };
    bool proven = false;
    int attempts = 0;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase(
        [&](const auto&) { return Proof{proven}; },
        [&](const auto& slot, const auto&)
        {
            ++attempts;
            ZhangSharedProductEventSearchVisit visit;
            visit.phase = "PAIR_SEARCH";
            visit.attempted = true;
            if (slot.id == "timeout")
            {
                ZhangExactWorkBudget pairBudget;
                pairBudget.maximumCheckpoints = 1;
                int candidate = 7;
                BOOST_CHECK(!zhangRunWithinExactWorkBudget<int>(pairBudget,
                    [] { return 99; }, candidate));
                BOOST_CHECK_EQUAL(candidate, 0); // interrupted work is not evidence
                visit.reason = "CURRENT_EVENT_COOPERATIVE_BUDGET_EXHAUSTED";
            }
            else
            {
                proven = true;
                visit.reason = "CURRENT_EVENT_OVERLAP_ADMITTED";
            }
            return visit;
        },
        [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(attempts, 2);
    BOOST_CHECK_EQUAL(report.attempted, 2);
    BOOST_CHECK_EQUAL(report.visits.at("timeout").reason,
        "CURRENT_EVENT_COOPERATIVE_BUDGET_EXHAUSTED");
    BOOST_CHECK(!report.visits.at("timeout").stopEvent);
    BOOST_CHECK(report.visits.at("healthy").currentQueryProven);
    std::map<std::string, Proof> finalProofs;
    BOOST_REQUIRE(batch.finalizeQueryPhase<Proof>(
        [&](const auto&) { return Proof{proven}; }, [] { return std::uint64_t{1}; },
        1, finalProofs, reason));
    BOOST_REQUIRE_EQUAL(finalProofs.size(), 2);
    BOOST_CHECK(finalProofs.at("timeout").proven);
    BOOST_CHECK(finalProofs.at("healthy").proven);
}

BOOST_AUTO_TEST_CASE(shared_product_event_alias_reuses_actual_failure_without_authority)
{
    auto batch = sharedEventTestBatch();
    const auto target = sharedEventTestTarget("A", 13);
    std::string reason;
    for (const auto& [id, prn] : std::vector<std::pair<std::string, int>>{
        {"original", 13}, {"alias", 14}, {"other-pair", 15}})
        BOOST_REQUIRE(batch.add({id, SatSys(E_Sys::GPS, prn), SatSys(E_Sys::GPS, 1)}, target, reason));
    struct Proof { bool proven = false; };
    std::map<std::pair<std::string, std::uint64_t>, ZhangHandoffPairOutcome> outcomes;
    int actualSearches = 0;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase(
        [](const auto&) { return Proof{}; },
        [&](const auto& slot, const auto&)
        {
            const std::string pair = slot.id == "other-pair" ? "pair-B" : "pair-A";
            const auto plan = std::make_pair(pair, std::uint64_t{4});
            ZhangSharedProductEventSearchVisit visit;
            const auto found = outcomes.find(plan);
            if (found != outcomes.end())
            {
                visit.phase = "PAIR_RESULT_REUSED";
                visit.reason = found->second.reason;
                visit.reused = true;
                return visit;
            }
            ++actualSearches;
            ZhangHandoffPairOutcome outcome;
            outcome.reason = pair == "pair-A"
                ? "CURRENT_EVENT_OVERLAP_VERSIONED_PHYSICAL_ARC_NOT_IN_CURRENT_POSTERIOR"
                : "CURRENT_EVENT_OVERLAP_RATIO_REJECTED";
            outcome.sourceSlot = slot.id;
            outcomes.emplace(plan, outcome);
            visit.phase = "PAIR_SEARCH";
            visit.reason = outcome.reason;
            visit.attempted = true;
            return visit;
        },
        [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(actualSearches, 2);
    BOOST_CHECK_EQUAL(report.attempted, 2);
    BOOST_CHECK_EQUAL(report.reused, 1);
    BOOST_CHECK_EQUAL(report.unvisited, 0);
    BOOST_CHECK_EQUAL(report.visits.at("original").reason, report.visits.at("alias").reason);
    BOOST_CHECK_EQUAL(report.visits.at("alias").phase, "PAIR_RESULT_REUSED");
    BOOST_CHECK(!report.visits.at("alias").attempted);
    std::map<std::string, Proof> proofs;
    BOOST_REQUIRE(batch.finalizeQueryPhase<Proof>([](const auto&) { return Proof{}; },
        [] { return std::uint64_t{4}; }, 4, proofs, reason));
    for (const auto& [id, proof] : proofs) BOOST_CHECK(!proof.proven);

    // Production continues cheap plan visits after the aggregate pair budget
    // is exhausted: aliases must still receive the original real outcome.
    // A new evidence generation must not reuse an old generation's failure.
    ZhangHandoffSearchBudget budget(1);
    const auto started = ZhangHandoffSearchBudget::Clock::time_point{};
    actualSearches = 0;
    BOOST_REQUIRE(batch.executeSearchPhase([](const auto&) { return Proof{}; },
        [&](const auto& slot, const auto&)
        {
            const std::string pair = slot.id == "other-pair" ? "pair-B" : "pair-A";
            const auto plan = std::make_pair(pair, std::uint64_t{5});
            ZhangSharedProductEventSearchVisit visit;
            if (const auto prior = outcomes.find(plan); prior != outcomes.end())
            {
                visit.phase = "PAIR_RESULT_REUSED";
                visit.reason = prior->second.reason;
                visit.reused = true;
                return visit;
            }
            if (budget.exhausted())
            {
                visit.phase = "NOT_SEARCHED";
                visit.reason = "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED";
                return visit; // do not stop later cheap alias visits
            }
            ++actualSearches;
            budget.charge(started, started + std::chrono::milliseconds(1));
            outcomes.emplace(plan, ZhangHandoffPairOutcome{
                "CURRENT_EVENT_OVERLAP_VERSIONED_PHYSICAL_ARC_NOT_IN_CURRENT_POSTERIOR", slot.id, false});
            visit.phase = "PAIR_SEARCH";
            visit.reason = outcomes.at(plan).reason;
            visit.attempted = true;
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(actualSearches, 1);
    BOOST_CHECK_EQUAL(report.attempted, 1);
    BOOST_CHECK_EQUAL(report.reused, 1);
    BOOST_CHECK_EQUAL(report.searchCallbacks, 3);
    BOOST_CHECK_EQUAL(report.unvisited, 0);
    BOOST_CHECK_EQUAL(report.visits.at("alias").reason, report.visits.at("original").reason);
    BOOST_CHECK_EQUAL(report.visits.at("other-pair").reason, "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED");
    BOOST_REQUIRE(batch.finalizeQueryPhase<Proof>([](const auto&) { return Proof{}; },
        [] { return std::uint64_t{5}; }, 5, proofs, reason));
    for (const auto& [id, proof] : proofs) BOOST_CHECK(!proof.proven);
}

BOOST_AUTO_TEST_CASE(shared_product_event_exhausted_budget_records_zero_attempts)
{
    auto batch = sharedEventTestBatch();
    const auto target = sharedEventTestTarget("A", 13);
    std::string reason;
    for (int prn = 13; prn <= 15; ++prn)
        BOOST_REQUIRE(batch.add({"G" + std::to_string(prn), SatSys(E_Sys::GPS, prn),
            SatSys(E_Sys::GPS, 1)}, target, reason));
    struct Proof { bool proven = false; };
    ZhangHandoffSearchBudget budget(1);
    const auto started = ZhangHandoffSearchBudget::Clock::time_point{};
    budget.charge(started, started + std::chrono::milliseconds(1));
    BOOST_REQUIRE(budget.exhausted());
    int attempts = 0;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase([](const auto&) { return Proof{}; },
        [&](const auto&, const auto&)
        {
            ZhangSharedProductEventSearchVisit visit;
            if (budget.exhausted())
            {
                visit.phase = "PAIR_SEARCH_NOT_STARTED";
                visit.reason = "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED";
                visit.stopEvent = true;
                visit.unvisitedReason = visit.reason;
                return visit;
            }
            ++attempts;
            visit.attempted = true;
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(attempts, 0);
    BOOST_CHECK_EQUAL(report.attempted, 0);
    BOOST_CHECK_EQUAL(report.searchCallbacks, 1);
    BOOST_CHECK_EQUAL(report.unvisited, 2); // the stop callback itself was visited
    BOOST_CHECK_EQUAL(reason, "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED");
    for (const auto& [id, visit] : report.visits)
    {
        BOOST_CHECK_EQUAL(visit.phase, "PAIR_SEARCH_NOT_STARTED");
        BOOST_CHECK_EQUAL(visit.reason, "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED");
        BOOST_CHECK(!visit.attempted);
        BOOST_CHECK(!visit.currentQueryProven);
    }
    // The production policy visits every cheap plan even with zero remaining
    // pair time, so aliases can be resolved; none may turn into a real attempt.
    BOOST_REQUIRE(batch.executeSearchPhase([](const auto&) { return Proof{}; },
        [&](const auto&, const auto&)
        {
            BOOST_REQUIRE(budget.exhausted());
            ZhangSharedProductEventSearchVisit visit;
            visit.phase = "NOT_SEARCHED";
            visit.reason = "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED";
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(report.attempted, 0);
    BOOST_CHECK_EQUAL(report.reused, 0);
    BOOST_CHECK_EQUAL(report.searchCallbacks, 3);
    BOOST_CHECK_EQUAL(report.unvisited, 0);
    for (const auto& [id, visit] : report.visits)
    {
        BOOST_CHECK_EQUAL(visit.reason, "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED");
        BOOST_CHECK(!visit.attempted);
        BOOST_CHECK(!visit.currentQueryProven);
    }
}

BOOST_AUTO_TEST_CASE(shared_product_event_requeries_later_targets_after_admission)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    BOOST_REQUIRE(batch.add({"A", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("A", 13), reason));
    BOOST_REQUIRE(batch.add({"B", SatSys(E_Sys::GPS, 14), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("B", 14), reason));
    struct Proof { bool proven = false; std::string reason = "EXACT_CURRENT_CONSEQUENCE"; };
    bool proven = false;
    int attempts = 0;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase([&](const auto&) { return Proof{proven}; },
        [&](const auto&, const auto&)
        {
            ++attempts;
            proven = true; // one search proves both complete targets
            ZhangSharedProductEventSearchVisit visit;
            visit.phase = "PAIR_SEARCH";
            visit.reason = "CURRENT_EVENT_OVERLAP_ADMITTED";
            visit.attempted = true;
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(report.initialQueries, 2);
    BOOST_CHECK_EQUAL(attempts, 1);
    BOOST_CHECK_EQUAL(report.visits.at("B").phase, "CURRENT_QUERY_PROVEN");
    BOOST_CHECK_EQUAL(report.visits.at("B").reason, "EXACT_CURRENT_CONSEQUENCE");
    BOOST_CHECK(!report.visits.at("B").attempted);
    BOOST_CHECK_EQUAL(report.unvisited, 0);
}

BOOST_AUTO_TEST_CASE(shared_product_event_query_and_search_exceptions_are_not_success)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    const auto target = sharedEventTestTarget("A", 13);
    BOOST_REQUIRE(batch.add({"A", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)}, target, reason));
    BOOST_REQUIRE(batch.add({"B", SatSys(E_Sys::GPS, 14), SatSys(E_Sys::GPS, 1)}, target, reason));
    struct Proof { bool proven = false; };
    int callbacks = 0;
    ZhangSharedProductEventSearchReport report;
    BOOST_CHECK_THROW(batch.executeSearchPhase(
        [](const auto&) -> Proof { throw std::runtime_error("controlled initial query failure"); },
        [&](const auto&, const auto&) { ++callbacks; return ZhangSharedProductEventSearchVisit{}; },
        [](const Proof& proof) { return proof.proven; }, reason, report), std::runtime_error);
    BOOST_CHECK_EQUAL(callbacks, 0);
    BOOST_CHECK_EQUAL(report.unvisited, 2);
    BOOST_CHECK_EQUAL(reason, "SHARED_EVENT_INITIAL_QUERY_EXCEPTION");
    for (const auto& [id, visit] : report.visits)
        BOOST_CHECK_EQUAL(visit.reason, "NOT_SEARCHED_INITIAL_QUERY_EXCEPTION");

    BOOST_CHECK_THROW(batch.executeSearchPhase(
        [](const auto&) { return Proof{}; },
        [](const auto&, const auto&) -> ZhangSharedProductEventSearchVisit { throw ZhangExactWorkBudgetExceeded{}; },
        [](const Proof& proof) { return proof.proven; }, reason, report), ZhangExactWorkBudgetExceeded);
    BOOST_CHECK_EQUAL(reason, "SHARED_EVENT_SEARCH_CALLBACK_EXCEPTION");
    BOOST_CHECK_EQUAL(report.visits.at("A").reason, "SEARCH_OUTCOME_UNAVAILABLE_CALLBACK_EXCEPTION");
    BOOST_CHECK_EQUAL(report.visits.at("B").reason, "NOT_SEARCHED_SEARCH_CALLBACK_EXCEPTION");
    BOOST_CHECK_EQUAL(report.unvisited, 1);
    BOOST_CHECK_EQUAL(report.attempted, 0); // unknown outcome is not invented
}

BOOST_AUTO_TEST_CASE(shared_product_event_final_generation_change_stops_before_next_query)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    BOOST_REQUIRE(batch.add({"A", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("A", 13), reason));
    BOOST_REQUIRE(batch.add({"B", SatSys(E_Sys::GPS, 14), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("B", 14), reason));
    struct Proof { bool proven = true; };
    std::uint64_t generation = 4;
    int queries = 0;
    std::map<std::string, Proof> proofs{{"stale", Proof{}}};
    BOOST_CHECK(!batch.finalizeQueryPhase<Proof>(
        [&](const auto&) { ++queries; ++generation; return Proof{}; },
        [&] { return generation; }, 4, proofs, reason));
    BOOST_CHECK_EQUAL(queries, 1);
    BOOST_CHECK(proofs.empty());
    BOOST_CHECK_EQUAL(reason, "SHARED_EVENT_EVIDENCE_GENERATION_CHANGED");
}

BOOST_AUTO_TEST_CASE(shared_product_event_retries_incomplete_plans_after_fair_first_round)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    const auto targetA = sharedEventTestTarget("A", 13);
    BOOST_REQUIRE(batch.add({"A", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)}, targetA, reason));
    BOOST_REQUIRE(batch.add({"A-alias", SatSys(E_Sys::GPS, 14), SatSys(E_Sys::GPS, 1)}, targetA, reason));
    BOOST_REQUIRE(batch.add({"B", SatSys(E_Sys::GPS, 15), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("B", 15), reason));
    struct Proof { bool proven = false; };
    ZhangHandoffSearchBudget budget(10000);
    auto now = ZhangHandoffSearchBudget::Clock::time_point{};
    std::map<std::string, ZhangHandoffPairOutcome> outcomes;
    std::set<std::string> proven;
    std::vector<std::string> actualOrder;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase(
        [&](const auto& target)
        {
            now += std::chrono::seconds(5); // shared queries stay separately timed even on retries
            return Proof{proven.contains(target.physicalEdges.front().receiver)};
        },
        [&](const auto& slot, const auto&, unsigned round)
        {
            const std::string pair = slot.id == "B" ? "B" : "A";
            auto prior = outcomes.find(pair);
            ZhangSharedProductEventSearchVisit visit;
            if (prior != outcomes.end() && (round == 0 || !prior->second.retryable()))
            {
                visit.reason = prior->second.reason;
                visit.reused = true;
                visit.incomplete = prior->second.incomplete();
                visit.lastAttemptReason = prior->second.reason;
                return visit;
            }
            BOOST_REQUIRE(!budget.exhausted());
            const unsigned completedAttempts = prior == outcomes.end() ? 0 : prior->second.attempts;
            const auto started = now;
            BOOST_CHECK(budget.pairDeadline(started, completedAttempts) == started +
                std::chrono::milliseconds(completedAttempts == 0 ? 1000 : 2000));
            actualOrder.push_back(pair + std::to_string(round));
            now += std::chrono::milliseconds(completedAttempts == 0 ? 1000 : 1500);
            budget.charge(started, now);
            const bool admitted = completedAttempts != 0;
            if (admitted) proven.insert(pair);
            auto& outcome = outcomes[pair];
            outcome = {admitted ? "CURRENT_EVENT_OVERLAP_ADMITTED"
                                : "CURRENT_EVENT_OVERLAP_END_TO_END_WALL_TIMEOUT",
                slot.id, admitted, completedAttempts + 1};
            visit.attempted = true;
            visit.reason = outcome.reason;
            visit.incomplete = outcome.incomplete();
            visit.retryPending = outcome.retryable() && !budget.exhausted();
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    const std::vector<std::string> expected{"A0", "B0", "A1", "B1"};
    BOOST_CHECK(actualOrder == expected);
    BOOST_CHECK_EQUAL(report.actualAttempts, 4);
    BOOST_CHECK_EQUAL(report.retries, 2);
    BOOST_CHECK_EQUAL(report.reused, 1);
    BOOST_CHECK_EQUAL(report.incomplete, 0);
    BOOST_CHECK_EQUAL(report.pendingRetries, 0);
    BOOST_CHECK_EQUAL(report.visits.at("A-alias").attemptCount, 0);
    BOOST_CHECK(report.visits.at("A-alias").currentQueryProven);
    // Even a later proven alias retains the real reused attempt reason. It
    // cannot rely on production's final unproven-alias refresh to restore it.
    BOOST_CHECK_EQUAL(report.visits.at("A-alias").lastAttemptReason,
        "CURRENT_EVENT_OVERLAP_END_TO_END_WALL_TIMEOUT");
    BOOST_CHECK_EQUAL(budget.consumedSeconds(), 5.0);
    BOOST_CHECK_EQUAL(budget.remainingSeconds(), 5.0);
}

BOOST_AUTO_TEST_CASE(shared_product_event_retries_charge_one_aggregate_budget_without_renewal)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    for (const auto& [id, prn] : std::vector<std::pair<std::string, int>>{{"A", 13}, {"B", 14}})
        BOOST_REQUIRE(batch.add({id, SatSys(E_Sys::GPS, prn), SatSys(E_Sys::GPS, 1)},
            sharedEventTestTarget(id, prn), reason));
    struct Proof { bool proven = false; };
    ZhangHandoffSearchBudget budget(5000);
    auto now = ZhangHandoffSearchBudget::Clock::time_point{};
    std::map<std::string, ZhangHandoffPairOutcome> outcomes;
    std::vector<std::string> order;
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase([](const auto&) { return Proof{}; },
        [&](const auto& slot, const auto&, unsigned round)
        {
            ZhangSharedProductEventSearchVisit visit;
            if (budget.exhausted())
            {
                visit.reason = "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED";
                visit.incomplete = true;
                return visit;
            }
            const auto prior = outcomes.find(slot.id);
            const unsigned completedAttempts = prior == outcomes.end() ? 0 : prior->second.attempts;
            const auto started = now;
            now = budget.pairDeadline(started, completedAttempts);
            budget.charge(started, now);
            order.push_back(slot.id + std::to_string(round));
            outcomes[slot.id] = {"CURRENT_EVENT_OVERLAP_WALL_TIMEOUT", slot.id,
                false, completedAttempts + 1};
            visit.reason = outcomes.at(slot.id).reason;
            visit.attempted = true;
            visit.incomplete = true;
            // Queue an incomplete representative even on exact exhaustion:
            // next callback records why no retry could actually start.
            visit.retryPending = true;
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    const std::vector<std::string> expected{"A0", "B0", "A1", "B1"};
    BOOST_CHECK(order == expected);
    BOOST_CHECK_EQUAL(report.actualAttempts, 4);
    BOOST_CHECK_EQUAL(report.retries, 2);
    BOOST_CHECK_EQUAL(report.incomplete, 2);
    BOOST_CHECK_EQUAL(report.pendingRetries, 0);
    BOOST_CHECK_EQUAL(budget.consumedSeconds(), 5.0);
    BOOST_CHECK_EQUAL(budget.remainingSeconds(), 0.0);
    BOOST_CHECK_EQUAL(report.visits.at("A").reason, "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED");
    BOOST_CHECK_EQUAL(report.visits.at("B").reason, "NOT_SEARCHED_SEARCH_BUDGET_EXHAUSTED");
    std::map<std::string, Proof> proofs;
    BOOST_REQUIRE(batch.finalizeQueryPhase<Proof>([](const auto&) { return Proof{}; },
        [] { return std::uint64_t{5}; }, 5, proofs, reason));
    for (const auto& [id, proof] : proofs) BOOST_CHECK(!proof.proven);
}

BOOST_AUTO_TEST_CASE(shared_product_event_retry_cache_is_generation_scoped)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    for (const auto& [id, prn] : std::vector<std::pair<std::string, int>>{{"A", 13}, {"B", 14}})
        BOOST_REQUIRE(batch.add({id, SatSys(E_Sys::GPS, prn), SatSys(E_Sys::GPS, 1)},
            sharedEventTestTarget(id, prn), reason));
    struct Proof { bool proven = false; };
    std::uint64_t generation = 0;
    std::set<std::string> proven;
    std::map<std::pair<std::string, std::uint64_t>, ZhangHandoffPairOutcome> outcomes;
    std::vector<std::pair<std::string, std::uint64_t>> attempts;
    ZhangHandoffSearchBudget budget(10000);
    auto now = ZhangHandoffSearchBudget::Clock::time_point{};
    ZhangSharedProductEventSearchReport report;
    BOOST_REQUIRE(batch.executeSearchPhase(
        [&](const auto& target) { return Proof{proven.contains(target.physicalEdges.front().receiver)}; },
        [&](const auto& slot, const auto&, unsigned round)
        {
            const auto key = std::make_pair(slot.id, generation);
            const auto prior = outcomes.find(key);
            ZhangSharedProductEventSearchVisit visit;
            if (prior != outcomes.end() && (round == 0 || !prior->second.retryable()))
            {
                visit.reused = true;
                visit.reason = prior->second.reason;
                return visit;
            }
            const unsigned completed = prior == outcomes.end() ? 0 : prior->second.attempts;
            // A stale generation's timeout must not grow the fresh domain's
            // first attempt or supply a cached failure on the new domain.
            BOOST_CHECK_EQUAL(completed, 0);
            BOOST_CHECK(budget.pairAllowance(completed) == std::chrono::milliseconds(1000));
            const auto started = now;
            now += std::chrono::milliseconds(100);
            budget.charge(started, now);
            attempts.push_back(key);
            const bool admitted = slot.id == "B" || generation != 0;
            outcomes[key] = {admitted ? "CURRENT_EVENT_OVERLAP_ADMITTED"
                                      : "CURRENT_EVENT_OVERLAP_WALL_TIMEOUT", slot.id, admitted, 1};
            if (admitted) { proven.insert(slot.id); ++generation; }
            visit.reason = outcomes.at(key).reason;
            visit.attempted = true;
            visit.incomplete = outcomes.at(key).incomplete();
            visit.retryPending = outcomes.at(key).retryable();
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    const std::vector<std::pair<std::string, std::uint64_t>> expected{{"A", 0}, {"B", 0}, {"A", 1}};
    BOOST_CHECK(attempts == expected);
    BOOST_CHECK_EQUAL(report.actualAttempts, 3);
    BOOST_CHECK_EQUAL(report.retries, 1);
    BOOST_CHECK_EQUAL(report.incomplete, 0);
    BOOST_CHECK_EQUAL(outcomes.size(), 3);
}

BOOST_AUTO_TEST_CASE(shared_product_event_retry_round_bound_and_terminal_results)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    for (const auto& [id, prn] : std::vector<std::pair<std::string, int>>{{"A", 13}, {"B", 14}})
        BOOST_REQUIRE(batch.add({id, SatSys(E_Sys::GPS, prn), SatSys(E_Sys::GPS, 1)},
            sharedEventTestTarget(id, prn), reason));
    struct Proof { bool proven = false; };
    ZhangSharedProductEventSearchReport report;
    std::vector<std::string> calls;
    BOOST_REQUIRE(batch.executeSearchPhase([](const auto&) { return Proof{}; },
        [&](const auto& slot, const auto&, unsigned round)
        {
            calls.push_back(slot.id + std::to_string(round));
            ZhangSharedProductEventSearchVisit visit;
            visit.attempted = true;
            visit.reason = slot.id == "A" ? "CURRENT_EVENT_OVERLAP_WALL_TIMEOUT"
                                           : "CURRENT_EVENT_OVERLAP_RATIO_REJECTED";
            visit.incomplete = zhangHandoffOutcomeIncomplete(visit.reason);
            visit.retryPending = zhangHandoffTimeoutRetryable(visit.reason);
            return visit; // intentionally zero fake clock cost: bound still holds
        }, [](const Proof& proof) { return proof.proven; }, reason, report));
    BOOST_CHECK_EQUAL(report.actualAttempts, ZhangHandoffSearchBudget::maximumRetryRounds + 2);
    BOOST_CHECK_EQUAL(report.retries, ZhangHandoffSearchBudget::maximumRetryRounds);
    BOOST_CHECK_EQUAL(report.visits.at("B").attemptCount, 1);
    BOOST_CHECK_EQUAL(report.pendingRetries, 1);
    BOOST_CHECK_EQUAL(report.incomplete, 1);
    BOOST_CHECK_EQUAL(reason, "SHARED_EVENT_RETRY_ROUND_LIMIT_INCOMPLETE");
    BOOST_CHECK_EQUAL(calls.front(), "A0");
    BOOST_CHECK_EQUAL(calls[1], "B0");
}

BOOST_AUTO_TEST_CASE(shared_product_event_retry_exception_is_not_a_completed_search)
{
    auto batch = sharedEventTestBatch();
    std::string reason;
    BOOST_REQUIRE(batch.add({"A", SatSys(E_Sys::GPS, 13), SatSys(E_Sys::GPS, 1)},
        sharedEventTestTarget("A", 13), reason));
    struct Proof { bool proven = false; };
    ZhangSharedProductEventSearchReport report;
    BOOST_CHECK_THROW(batch.executeSearchPhase([](const auto&) { return Proof{}; },
        [](const auto&, const auto&, unsigned round)
        {
            if (round != 0) throw std::runtime_error("controlled retry failure");
            ZhangSharedProductEventSearchVisit visit;
            visit.reason = "CURRENT_EVENT_OVERLAP_WALL_TIMEOUT";
            visit.attempted = visit.incomplete = visit.retryPending = true;
            return visit;
        }, [](const Proof& proof) { return proof.proven; }, reason, report), std::runtime_error);
    BOOST_CHECK_EQUAL(reason, "SHARED_EVENT_RETRY_CALLBACK_EXCEPTION");
    BOOST_CHECK_EQUAL(report.visits.at("A").reason, "SEARCH_OUTCOME_UNAVAILABLE_CALLBACK_EXCEPTION");
    BOOST_CHECK_EQUAL(report.visits.at("A").attemptCount, 1);
    BOOST_CHECK(report.visits.at("A").incomplete);
    BOOST_CHECK(!report.visits.at("A").retryPending);
}
