#include "common/zhangNativeProductTargetAdmission.hpp"

#include <cassert>
#include <iostream>

int main()
{
    // Zero supported GPS targets must not invalidate independently certified
    // GPS rows or the same complete joint posterior serving Galileo and BDS.
    const auto gpsPending = zhangNativeProductTargetAdmission(true, true, true, 0);
    const auto galileo = zhangNativeProductTargetAdmission(true, true, true, 22);
    const auto bds = zhangNativeProductTargetAdmission(true, true, true, 26);
    assert(!gpsPending.rejectSharedBranch && gpsPending.targetsPending);
    assert(!galileo.rejectSharedBranch && !galileo.targetsPending);
    assert(!bds.rejectSharedBranch && !bds.targetsPending);

    // The single-target case is also pending. No missing second signal or
    // missing target is converted into a zero observation/covariance row.
    const auto one = zhangNativeProductTargetAdmission(true, true, true, 1);
    assert(!one.rejectSharedBranch && one.targetsPending);
    const auto two = zhangNativeProductTargetAdmission(true, true, true, 2);
    assert(!two.rejectSharedBranch && !two.targetsPending);

    // Invalid shared evidence or full covariance remains a whole-branch
    // rejection regardless of how many apparently supported targets exist.
    for (const unsigned targets : {0u, 1u, 2u, 29u})
    {
        assert(zhangNativeProductTargetAdmission(false, true, true, targets).rejectSharedBranch);
        assert(zhangNativeProductTargetAdmission(true, false, true, targets).rejectSharedBranch);
        assert(zhangNativeProductTargetAdmission(true, true, false, targets).rejectSharedBranch);
    }
    std::cout << "NATIVE_PRODUCT_TARGET_ADMISSION_PASS" << std::endl;
}
