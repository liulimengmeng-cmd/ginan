#pragma once

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstring>
#include <string>
#include <string_view>

// Charge only actual pair attempts. Shared exact-domain compilation and
// consequence queries must not consume the allowance before a pair starts.
// This is a cooperative limit, not a hard wall-clock bound on Eigen/Boost.
class ZhangHandoffSearchBudget
{
public:
    using Clock = std::chrono::steady_clock;
    static constexpr long defaultTotalMilliseconds = 10000;
    static constexpr long pairMilliseconds = 1000;
    // An incomplete attempt is not a negative mathematical result. Give
    // distinct plans their first slice before revisiting any of them, and
    // bound retries even if a callback reports a timeout without using time.
    static constexpr unsigned maximumRetryRounds = 8;

    static bool parseTotalMilliseconds(const char* setting, long& milliseconds)
    {
        milliseconds = defaultTotalMilliseconds;
        if (!setting || !*setting) return true;
        long parsed = 0;
        const auto end = setting + std::strlen(setting);
        const auto result = std::from_chars(setting, end, parsed);
        if (result.ec != std::errc{} || result.ptr != end || parsed < 1 || parsed > 600000)
            return false;
        milliseconds = parsed;
        return true;
    }

    explicit ZhangHandoffSearchBudget(long milliseconds = defaultTotalMilliseconds)
        : allowance_(std::chrono::milliseconds(std::clamp(milliseconds, 1L, 600000L))) {}

    Clock::duration remaining() const
    {
        return consumed_ >= allowance_ ? Clock::duration::zero() : allowance_ - consumed_;
    }
    bool exhausted() const { return remaining() == Clock::duration::zero(); }
    Clock::duration pairAllowance(unsigned completedAttempts = 0) const
    {
        const long multiplier = 1L << std::min(completedAttempts, maximumRetryRounds);
        const Clock::duration slice = std::chrono::milliseconds(pairMilliseconds * multiplier);
        return std::min(remaining(), slice);
    }
    Clock::time_point pairDeadline(Clock::time_point start,
                                   unsigned completedAttempts = 0) const
    {
        return start + pairAllowance(completedAttempts);
    }
    void charge(Clock::time_point start, Clock::time_point finish)
    {
        if (finish > start) consumed_ += finish - start;
    }
    double consumedSeconds() const
    {
        return std::chrono::duration<double>(consumed_).count();
    }
    double remainingSeconds() const
    {
        return std::chrono::duration<double>(remaining()).count();
    }

private:
    Clock::duration allowance_;
    Clock::duration consumed_ = Clock::duration::zero();
};

// End-to-end EVENT allowance is separate from the unchanged 10-second sum
// of actual pair attempts. Shared domain preparation/consequence queries and
// final stable-generation verification now have a cooperative deadline too.
// Nested exact budgets must use min(pairDeadline, eventDeadline). This is NOT
// a hard wall bound on a single Boost/Eigen operation or a worker thread.
class ZhangHandoffEventBudget
{
public:
    using Clock = ZhangHandoffSearchBudget::Clock;
    static constexpr long defaultTotalMilliseconds = 60000;
    static constexpr long maximumTotalMilliseconds = 600000;

    static bool parseTotalMilliseconds(const char* setting, long& milliseconds)
    {
        milliseconds = defaultTotalMilliseconds;
        if (!setting || !*setting) return true;
        long parsed = 0;
        const auto end = setting + std::strlen(setting);
        const auto result = std::from_chars(setting, end, parsed);
        if (result.ec != std::errc{} || result.ptr != end || parsed < 1 ||
            parsed > maximumTotalMilliseconds) return false;
        milliseconds = parsed;
        return true;
    }

    explicit ZhangHandoffEventBudget(long milliseconds = defaultTotalMilliseconds,
        Clock::time_point started = Clock::now())
        : started_(started), deadline_(started + std::chrono::milliseconds(
            std::clamp(milliseconds, 1L, maximumTotalMilliseconds))) {}

    Clock::time_point deadline() const { return deadline_; }
    bool exhausted(Clock::time_point now = Clock::now()) const { return now >= deadline_; }
    Clock::time_point pairDeadline(const ZhangHandoffSearchBudget& search,
        Clock::time_point started, unsigned completedAttempts = 0) const
    { return std::min(deadline_, search.pairDeadline(started, completedAttempts)); }
    double elapsedSeconds(Clock::time_point now = Clock::now()) const
    { return std::chrono::duration<double>(now - started_).count(); }
    double remainingSeconds(Clock::time_point now = Clock::now()) const
    { return std::chrono::duration<double>(std::max(Clock::duration::zero(), deadline_ - now)).count(); }

private:
    Clock::time_point started_, deadline_;
};

// Deliberately use an explicit whitelist. A low ratio, incomplete top-two
// enumeration, node limit, invalid binding, or exact operation-count limit
// is not made successful merely by extending the wall-clock allowance.
inline bool zhangHandoffTimeoutRetryable(std::string_view reason)
{
    return reason == "CURRENT_EVENT_OVERLAP_WALL_TIMEOUT" ||
        reason == "CURRENT_EVENT_OVERLAP_END_TO_END_WALL_TIMEOUT" ||
        reason == "CURRENT_EVENT_PAIR_COOPERATIVE_BUDGET_EXHAUSTED" ||
        reason == "CURRENT_EVENT_COOPERATIVE_BUDGET_EXHAUSTED";
}

inline bool zhangHandoffOutcomeIncomplete(std::string_view reason)
{
    return zhangHandoffTimeoutRetryable(reason) ||
        reason == "SHARED_EVENT_TOTAL_COOPERATIVE_BUDGET_EXHAUSTED" ||
        reason == "SHARED_EVENT_PREPARATION_COOPERATIVE_BUDGET_EXHAUSTED" ||
        reason == "SHARED_EVENT_FINAL_QUERY_COOPERATIVE_BUDGET_EXHAUSTED" ||
        reason == "CURRENT_EVENT_OVERLAP_EXACT_OPERATION_LIMIT" ||
        reason == "CURRENT_EVENT_OVERLAP_INCOMPLETE_TOP_TWO" ||
        reason == "CURRENT_EVENT_OVERLAP_NODE_LIMIT" ||
        reason == "CURRENT_EVENT_SEARCH_NOT_COMPLETED";
}

// Deduplication is a work-sharing decision, not a mathematical outcome.
// Preserve the original attempt's result when another slot uses that plan.
struct ZhangHandoffPairOutcome
{
    std::string reason = "CURRENT_EVENT_SEARCH_NOT_COMPLETED";
    std::string sourceSlot;
    bool admitted = false;
    unsigned attempts = 1;

    bool incomplete() const { return !admitted && zhangHandoffOutcomeIncomplete(reason); }
    bool retryable() const { return !admitted && zhangHandoffTimeoutRetryable(reason); }
};
