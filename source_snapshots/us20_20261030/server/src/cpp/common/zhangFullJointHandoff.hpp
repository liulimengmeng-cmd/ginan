#pragma once
#include <map>
#include <string>
#include <memory>
#include <tuple>
#include <cstdint>
#include <limits>
#include "common/zhangUserIntegerFunctional.hpp"

/** A segment target is structural information, never a historical integer
 * decision. The complete native-to-segment difference must be re-certified
 * against the CURRENT epoch after every change, including a reference change.
 * Summation cancels exact versioned terms; it never saturates an integer image. */
inline bool zhangFullJointTransportTailMatches(
    const std::string& publishedNativeFingerprint,
    const std::string& receiptOldNativeFingerprint,
    int publishedBackend, int receiptOldBackend,
    int receiptNewBackend, int currentBackend)
{
    return !publishedNativeFingerprint.empty() && publishedNativeFingerprint != "INVALID" &&
        publishedNativeFingerprint == receiptOldNativeFingerprint && publishedBackend >= 0 &&
        receiptOldBackend == publishedBackend && receiptNewBackend >= 0 &&
        receiptNewBackend == currentBackend;
}

inline ZhangProductIntegerTransition zhangComposeFullJointDifference(
    const ZhangProductIntegerTransition& prefix,
    const ZhangProductIntegerTransition& increment)
{
    ZhangProductIntegerTransition out;
    out.observable = increment.observable;
    if (!prefix.valid || !increment.valid ||
        prefix.observable != increment.observable ||
        prefix.physicalEdges.size() != prefix.coefficients.size() ||
        prefix.physicalEdges.size() != prefix.physicalArcVersions.size() ||
        increment.physicalEdges.size() != increment.coefficients.size() ||
        increment.physicalEdges.size() != increment.physicalArcVersions.size())
    { out.failureReason = "FULL_JOINT_DIFFERENCE_IDENTITY_MISMATCH"; return out; }
    std::map<std::pair<ZhangGraphEdge, int>, ZhangExactInteger> terms;
    for (const auto* difference : {&prefix, &increment})
        for (std::size_t i = 0; i < difference->coefficients.size(); i++)
            terms[{difference->physicalEdges[i], difference->physicalArcVersions[i]}]
                += difference->coefficients[i];
    for (const auto& [identity, coefficient] : terms)
    {
        if (coefficient == 0) continue;
        if (identity.second < 0)
        { out.failureReason = "INVALID_FULL_JOINT_SIGNAL_ARC"; return out; }
        out.physicalEdges.push_back(identity.first);
        out.physicalArcVersions.push_back(identity.second);
        out.coefficients.push_back(coefficient);
    }
    out.affineOffsetCycles = prefix.affineOffsetCycles + increment.affineOffsetCycles;
    out.valid = true;
    out.failureReason = "NONE";
    return out;
}

/** A structural chain is anchored to the last COMMITTED public endpoint,
 * but its tail follows every proven same-time chart change, including epochs
 * at which no AR product is published.  It stores no integer values, means,
 * covariances or acceptance flags.  The cumulative target is always queried
 * anew in the current integer domain before it can authorize publication. */
struct ZhangFullJointStructuralChain
{
    static constexpr std::uint64_t maximumLinks = 4096;
    static constexpr std::size_t maximumTerms = 65536;
    std::string publishedNativeFingerprint;
    std::string publishedSegment;
    int publishedBackend = -1;
    std::string tailNativeFingerprint;
    int tailBackend = -1;
    std::uint64_t linkCount = 0;
    ZhangProductIntegerTransition nativeToSegment;

    bool anchoredTo(const std::string& fingerprint, int backend,
        const std::string& segment) const
    {
        return nativeToSegment.valid && !fingerprint.empty() && fingerprint != "INVALID" &&
            !segment.empty() && publishedNativeFingerprint == fingerprint &&
            publishedBackend == backend && backend >= 0 && publishedSegment == segment;
    }

    bool tailMatches(const std::string& fingerprint, int backend) const
    {
        return nativeToSegment.valid && !fingerprint.empty() && fingerprint != "INVALID" &&
            tailNativeFingerprint == fingerprint && tailBackend == backend && backend >= 0;
    }
};

inline ZhangFullJointStructuralChain zhangAnchorFullJointStructuralChain(
    const std::string& fingerprint, int backend, const std::string& segment,
    const ZhangProductIntegerTransition& nativeToSegment)
{
    ZhangFullJointStructuralChain out;
    if (!nativeToSegment.valid || fingerprint.empty() || fingerprint == "INVALID" ||
        backend < 0 || segment.empty()) return out;
    out.publishedNativeFingerprint = out.tailNativeFingerprint = fingerprint;
    out.publishedBackend = out.tailBackend = backend;
    out.publishedSegment = segment;
    out.nativeToSegment = nativeToSegment;
    return out;
}

inline bool zhangAdvanceFullJointStructuralChain(
    ZhangFullJointStructuralChain& chain,
    const std::string& receiptOldFingerprint, int receiptOldBackend,
    const std::string& currentFingerprint, int receiptNewBackend, int currentBackend,
    const ZhangProductIntegerTransition& increment, bool fullSameTimeMapProven,
    std::string& reason)
{
    if (!fullSameTimeMapProven)
    { reason = "STRUCTURAL_CHAIN_FULL_MOMENT_MAP_UNPROVEN"; return false; }
    if (!zhangFullJointTransportTailMatches(chain.tailNativeFingerprint,
        receiptOldFingerprint, chain.tailBackend, receiptOldBackend,
        receiptNewBackend, currentBackend) || currentFingerprint.empty() ||
        currentFingerprint == "INVALID")
    { reason = "STRUCTURAL_CHAIN_ENDPOINT_MISMATCH"; return false; }
    // A receipt represents a directed change. Reject stale/reverse receipts;
    // replay of the same current tail is handled by the caller, not composed.
    if (receiptNewBackend < receiptOldBackend ||
        chain.linkCount >= ZhangFullJointStructuralChain::maximumLinks)
    { reason = "STRUCTURAL_CHAIN_GENERATION_OR_COUNT_INVALID"; return false; }
    auto complete = zhangComposeFullJointDifference(chain.nativeToSegment, increment);
    if (!complete.valid)
    { reason = complete.failureReason; return false; }
    if (complete.coefficients.size() > ZhangFullJointStructuralChain::maximumTerms)
    { reason = "STRUCTURAL_CHAIN_RESOURCE_LIMIT"; return false; }
    chain.nativeToSegment = std::move(complete);
    chain.tailNativeFingerprint = currentFingerprint;
    chain.tailBackend = currentBackend;
    ++chain.linkCount;
    reason = "STRUCTURAL_CHAIN_CURRENT_NATIVE_TAIL_PROVEN";
    return true;
}

struct ZhangFullJointHandoffTask
{
    std::string eventIdentity;
    bool mapProven = false;
    bool structuralZero = false;
    bool currentValueProven = false;
    std::uint64_t structuralLinkCount = 0;
    std::string mapReason = "NOT_EVALUATED";
    std::string bridgeReason = "NOT_EVALUATED";
    ZhangProductIntegerTransition nativeToSegment;
    ZhangExactInteger currentValue = 0;
    std::vector<std::string> proofIds;
    std::vector<std::string> proofStatements;
    std::shared_ptr<const std::vector<std::string>> sharedProofStatements;
};
