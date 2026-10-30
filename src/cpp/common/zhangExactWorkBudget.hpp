#pragma once

#include <chrono>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>

/** Internal cooperative cancellation signal. It is deliberately not a
 * std::exception, so unrelated parser/arithmetic catch blocks cannot swallow
 * it and continue with a partial exact result. Every budget-owning production
 * entry must catch it before publication/admission. */
struct ZhangExactWorkBudgetExceeded {};

struct ZhangExactWorkBudget
{
    using Clock=std::chrono::steady_clock;
    Clock::time_point deadline=Clock::time_point::max();
    // Zero means unlimited. A finite operation budget also supports a
    // deterministic mid-elimination cancellation regression test.
    std::uint64_t maximumCheckpoints=0;
    std::uint64_t checkpoints=0,clockChecks=0,libraryBoundaries=0;
    bool exhausted=false;
    const char* phase="NOT_STARTED";
    const char* reason="NONE";
};

inline thread_local ZhangExactWorkBudget* zhangActiveExactWorkBudget=nullptr;

class ZhangScopedExactWorkBudget
{
    ZhangExactWorkBudget* previous;
public:
    explicit ZhangScopedExactWorkBudget(ZhangExactWorkBudget& budget)
        :previous(zhangActiveExactWorkBudget)
    {
        if(previous && previous->deadline<budget.deadline)budget.deadline=previous->deadline;
        zhangActiveExactWorkBudget=&budget;
    }
    ~ZhangScopedExactWorkBudget() {zhangActiveExactWorkBudget=previous;}
    ZhangScopedExactWorkBudget(const ZhangScopedExactWorkBudget&)=delete;
    ZhangScopedExactWorkBudget& operator=(const ZhangScopedExactWorkBudget&)=delete;
};

inline void zhangExactWorkCheckpoint(const char* phase,bool forceClockCheck=false)
{
    auto* budget=zhangActiveExactWorkBudget;
    if(!budget)return; // No deadline and no clock/counter cost in normal paths.
    budget->phase=phase;++budget->checkpoints;
    if(budget->maximumCheckpoints && budget->checkpoints>=budget->maximumCheckpoints) {
        budget->exhausted=true;budget->reason="EXACT_OPERATION_LIMIT";
        throw ZhangExactWorkBudgetExceeded{};
    }
    if(forceClockCheck || (budget->checkpoints&63)==0 || budget->exhausted) {
        ++budget->clockChecks;
        if(budget->exhausted || ZhangExactWorkBudget::Clock::now()>=budget->deadline) {
            budget->exhausted=true;budget->reason="WALL_TIMEOUT";
            throw ZhangExactWorkBudgetExceeded{};
        }
    }
}

/** Eigen factorizations, allocation/copies, and one multiprecision operation
 * are not hard-interruptible here. Check both sides of such a call; do not
 * advertise the whole overlap as a hard wall-clock upper bound. */
inline void zhangExactLibraryBoundary(const char* phase)
{
    if(zhangActiveExactWorkBudget)++zhangActiveExactWorkBudget->libraryBoundaries;
    zhangExactWorkCheckpoint(phase,true);
}

template<class Result,class Work>
inline bool zhangRunWithinExactWorkBudget(ZhangExactWorkBudget& budget,
    Work&& work,Result& output)
{
    ZhangScopedExactWorkBudget scope(budget);
    try {
        zhangExactWorkCheckpoint("EXACT_BUDGET_ENTRY",true);
        Result candidate=std::forward<Work>(work)();
        zhangExactWorkCheckpoint("EXACT_BUDGET_RESULT",true);
        output=std::move(candidate);return true;
    } catch(const ZhangExactWorkBudgetExceeded&) {
        output=Result{};return false;
    }
}
