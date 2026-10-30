#include "common/zhangHandoffSearchBudget.hpp"
#include "common/zhangExactWorkBudget.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message)
{ if (!condition) throw std::runtime_error(message); }

int main()
{
    using Event = ZhangHandoffEventBudget;
    using Search = ZhangHandoffSearchBudget;
    using namespace std::chrono;
    const auto start = Event::Clock::time_point{};
    Event event(60000, start);
    Search pairs(10000);
    require(event.deadline() == start + seconds(60), "event deadline not absolute");
    require(!event.exhausted(start + seconds(59)) && event.exhausted(start + seconds(60)),
        "event end boundary incorrect");
    require(event.pairDeadline(pairs, start + seconds(59)) == start + seconds(60),
        "pair extended the event end-to-end deadline");
    require(event.pairDeadline(pairs, start + milliseconds(59900), 4) == start + seconds(60),
        "larger retry renewed the event deadline");
    require(event.remainingSeconds(start + seconds(100)) == 0,
        "expired event reported negative remaining time");
    // Expensive shared preparation does not steal the separate pair allowance,
    // but the event cannot run indefinitely before its first actual ILS call.
    const auto afterPreparation = start + seconds(435);
    require(pairs.remainingSeconds() == 10 && event.exhausted(afterPreparation),
        "shared preparation incorrectly free in event accounting");
    require(event.pairDeadline(pairs, afterPreparation) == event.deadline(),
        "overrun preparation authorized a new pair slice");
    long setting = 0;
    require(Event::parseTotalMilliseconds(nullptr, setting) && setting == 60000,
        "default event budget changed");
    require(Event::parseTotalMilliseconds("600000", setting) && setting == 600000,
        "maximum configured event budget rejected");
    for (const auto invalid : {"0", "-5", "600001", "60000junk", " 60000", "999999999999999999999"})
        require(!Event::parseTotalMilliseconds(invalid, setting) && setting == 60000,
            "invalid event budget accepted or default corrupted");
    require(Event(-1, start).deadline() == start + milliseconds(1), "negative event allowance escaped clamp");
    require(Event(900000, start).deadline() == start + seconds(600), "event allowance escaped upper cap");
    // Nested budgets inherit the earlier EVENT deadline, even if a retry asks
    // for a longer allowance. No sleeps or wall-time-sensitive expectations.
    ZhangExactWorkBudget outer;
    outer.deadline = Event::Clock::now() - seconds(1);
    ZhangExactWorkBudget inner;
    inner.deadline = Event::Clock::time_point::max();
    bool caught = false;
    {
        ZhangScopedExactWorkBudget outerScope(outer);
        try
        {
            ZhangScopedExactWorkBudget innerScope(inner);
            require(inner.deadline == outer.deadline, "nested exact budget overrode outer event deadline");
            zhangExactWorkCheckpoint("TEST_SHARED_PREPARATION_ENTRY", true);
        }
        catch (const ZhangExactWorkBudgetExceeded&) { caught = true; }
        require(zhangActiveExactWorkBudget == &outer, "unwinding lost outer budget ownership");
    }
    require(caught && inner.exhausted && !zhangActiveExactWorkBudget,
        "expired shared preparation was not cooperatively cancelled");
    for (const auto reason : {"SHARED_EVENT_PREPARATION_COOPERATIVE_BUDGET_EXHAUSTED",
         "SHARED_EVENT_TOTAL_COOPERATIVE_BUDGET_EXHAUSTED",
         "SHARED_EVENT_FINAL_QUERY_COOPERATIVE_BUDGET_EXHAUSTED"})
    {
        require(zhangHandoffOutcomeIncomplete(reason), "unfinished event labeled completed mathematical failure");
        require(!zhangHandoffTimeoutRetryable(reason), "exhausted event renewed by per-pair retry");
    }
    std::cout << "HANDOFF_EVENT_BUDGET_PASS event_default_ms=60000 pair_total_ms=10000"
        << " preparation_charged_to_event=1 nested_deadline=1 unfinished_is_not_proof=1 hard_wall_bound=0\n";
}
