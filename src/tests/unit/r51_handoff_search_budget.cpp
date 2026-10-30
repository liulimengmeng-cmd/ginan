#include "common/zhangHandoffSearchBudget.hpp"

#include <iostream>
#include <map>
#include <stdexcept>

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

int main()
{
    using Budget = ZhangHandoffSearchBudget;
    using namespace std::chrono;
    const auto epoch = Budget::Clock::time_point{};
    // An arbitrarily slow shared query does not charge the pair allowance.
    Budget budget(2500);
    const auto afterQuery = epoch + seconds(435);
    require(budget.pairDeadline(afterQuery) == afterQuery + seconds(1), "query stole pair allowance");
    budget.charge(afterQuery, afterQuery + milliseconds(1100));
    // One cooperative timeout does not suppress the next independent pair.
    const auto second = afterQuery + seconds(20);
    require(budget.pairDeadline(second) == second + seconds(1), "timeout suppressed next pair");
    budget.charge(second, second + milliseconds(1000));
    const auto third = second + seconds(20);
    require(budget.pairDeadline(third) == third + milliseconds(400), "remaining allowance renewed");
    budget.charge(third, third + milliseconds(450));
    require(budget.exhausted() && budget.pairDeadline(third) == third, "exhaustion accepted another pair");
    budget.charge(third, third - seconds(1));
    require(budget.exhausted(), "negative duration refunded budget");

    long value = 0;
    require(Budget::parseTotalMilliseconds(nullptr, value) && value == 10000, "default budget");
    require(Budget::parseTotalMilliseconds("20000", value) && value == 20000, "configured budget");
    for (const auto invalid : {"0", "-1", "1junk", "600001", "999999999999999999999999", " 1000"})
        require(!Budget::parseTotalMilliseconds(invalid, value) && value == 10000, "invalid budget accepted");
    std::map<int, ZhangHandoffPairOutcome> attempts;
    attempts.emplace(7, ZhangHandoffPairOutcome{"CURRENT_EVENT_OVERLAP_RATIO_REJECTED", "G18|L1C", false});
    require(attempts.at(7).reason == "CURRENT_EVENT_OVERLAP_RATIO_REJECTED" &&
        !attempts.at(7).admitted, "alias lost actual failure result");
    require(Budget(100).pairDeadline(epoch) == epoch + milliseconds(100), "short total budget exceeded");
    Budget retries(10000);
    require(retries.pairDeadline(epoch, 0) == epoch + seconds(1), "first slice changed");
    retries.charge(epoch, epoch + seconds(1));
    require(retries.pairDeadline(epoch, 1) == epoch + seconds(2), "retry did not grow allowance");
    retries.charge(epoch, epoch + seconds(2));
    require(retries.pairDeadline(epoch, 2) == epoch + seconds(4), "second retry did not grow allowance");
    retries.charge(epoch, epoch + seconds(4));
    require(retries.pairDeadline(epoch, 3) == epoch + seconds(3), "retry renewed aggregate budget");
    require(retries.pairDeadline(epoch, 1000000) == epoch + seconds(3), "retry shift overflow");
    require(retries.remainingSeconds() == 3.0, "remaining time audit incorrect");
    for (const auto retryable : {"CURRENT_EVENT_OVERLAP_WALL_TIMEOUT",
         "CURRENT_EVENT_OVERLAP_END_TO_END_WALL_TIMEOUT",
         "CURRENT_EVENT_PAIR_COOPERATIVE_BUDGET_EXHAUSTED",
         "CURRENT_EVENT_COOPERATIVE_BUDGET_EXHAUSTED"})
    {
        ZhangHandoffPairOutcome outcome{retryable, "A", false};
        require(outcome.incomplete() && outcome.retryable(), "wall timeout treated as terminal proof");
        outcome.admitted = true;
        require(!outcome.incomplete() && !outcome.retryable(), "admitted outcome retried");
    }
    for (const auto terminal : {"CURRENT_EVENT_OVERLAP_RATIO_REJECTED",
         "CURRENT_EVENT_OVERLAP_INCOMPLETE_TOP_TWO",
         "CURRENT_EVENT_OVERLAP_NODE_LIMIT",
         "CURRENT_EVENT_OVERLAP_EXACT_OPERATION_LIMIT",
         "CURRENT_EVENT_OVERLAP_VERSIONED_PHYSICAL_ARC_NOT_IN_CURRENT_POSTERIOR"})
        require(!zhangHandoffTimeoutRetryable(terminal), "non-wall outcome retried");
    require(zhangHandoffOutcomeIncomplete("CURRENT_EVENT_OVERLAP_EXACT_OPERATION_LIMIT"),
        "operation-limit result mistaken for completed evidence");
    require(zhangHandoffOutcomeIncomplete("CURRENT_EVENT_OVERLAP_INCOMPLETE_TOP_TWO") &&
        zhangHandoffOutcomeIncomplete("CURRENT_EVENT_OVERLAP_NODE_LIMIT"),
        "incomplete ILS result mistaken for completed evidence");
    std::cout << "handoff search budget and outcome contracts passed\n";
}
