// Links the completed production zhangPppAr.cpp object. The isolated audit seam
// installs and calls its actual exactStateTransformCallback and restores its
// own reserved runtime; no network input, checkpoint, or solver is consumed.
#include <iostream>
#include <string>
#include <type_traits>
#include "common/zhangFactorCapture.hpp"

static_assert(std::is_nothrow_move_assignable_v<ZhangFactorCaptureBuffer>,
    "Final capture commit must remain a nonthrowing move");

bool zhangExactTransformConsumerAtomicityContractForAudit(std::string& report);

int main()
{
    std::string report;
    const bool passed = zhangExactTransformConsumerAtomicityContractForAudit(report);
    std::cout << report << '\n';
    return passed ? 0 : 1;
}
