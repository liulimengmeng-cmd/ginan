#pragma once

#include <cstddef>

struct ZhangNativeProductTargetAdmission
{
    bool rejectSharedBranch = true;
    bool targetsPending = false;
    const char* reason = "ZERO_DATUM_NATIVE_CURRENT_EVIDENCE_INVALID";
};

// A native system contributes independently certified current integer rows to
// one full joint posterior. An unavailable product target does not invalidate
// those rows. Missing/invalid shared evidence or covariance does invalidate the
// branch, while fewer than two supported dual-frequency targets suspends this
// system's products only. No target receives authority from this decision.
inline ZhangNativeProductTargetAdmission zhangNativeProductTargetAdmission(
    bool currentEvidenceVerified,
    bool jointCovarianceStaged,
    bool regionalJointReady,
    std::size_t supportedDualTargets)
{
    if (!currentEvidenceVerified)
        return {};
    if (!jointCovarianceStaged)
        return {true, false, "ZERO_DATUM_NATIVE_SHARED_JOINT_COVARIANCE_MISSING"};
    if (!regionalJointReady)
        return {true, false, "ZERO_DATUM_NATIVE_REGIONAL_JOINT_COVARIANCE_INVALID"};
    if (supportedDualTargets < 2)
        return {false, true, "ZERO_DATUM_NATIVE_DUAL_SIGNAL_TARGET_SUBSET_PENDING"};
    return {false, false, "CURRENT_NATIVE_POSTERIOR_AND_TARGET_SUBSET_VERIFIED"};
}
