#pragma once

#include "common/zhangExactWorkBudget.hpp"
#include "common/zhangHandoffSearchBudget.hpp"

#include <utility>

// The production search inserts its residual ticket before enumeration. A
// wall timeout must release it, including cooperative stack unwinding before
// the entry point's catch assigns a user-facing reason. Terminal decisions
// and operation-count failures remain deduplicated.
template<class ERASE>
class ZhangHandoffIncompleteTicket
{
    ERASE erase_;
    const std::string& reason_;
    const ZhangExactWorkBudget& budget_;
public:
    ZhangHandoffIncompleteTicket(ERASE erase, const std::string& reason,
        const ZhangExactWorkBudget& budget)
        : erase_(std::move(erase)), reason_(reason), budget_(budget) {}
    ~ZhangHandoffIncompleteTicket()
    {
        if (zhangHandoffTimeoutRetryable(reason_) ||
            (budget_.exhausted && std::string(budget_.reason) == "WALL_TIMEOUT")) erase_();
    }
    ZhangHandoffIncompleteTicket(const ZhangHandoffIncompleteTicket&) = delete;
    ZhangHandoffIncompleteTicket& operator=(const ZhangHandoffIncompleteTicket&) = delete;
};
