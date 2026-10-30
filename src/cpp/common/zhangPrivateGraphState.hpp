#pragma once

#include <atomic>
#include <memory>
#include <string>
#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"

// Runtime maps are keyed by identity, not by KFState address. A plain KF copy
// therefore MUST NOT acquire cleanup rights over its inherited owner identity.
inline std::string zhangPrivateGraphIdentity(const KFState& state)
{
    const auto branch = state.metaDataMap.find(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA);
    if (branch != state.metaDataMap.end()) return branch->second;
    const auto owner = state.metaDataMap.find(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);
    return owner == state.metaDataMap.end() ? std::string{} : owner->second;
}

struct ZhangPrivateGraphStateDeleter
{
    std::string ownerId;
    std::string ownedId;
    void (*erase)(const KFState&) = nullptr;

    void operator()(KFState* state) const noexcept
    {
        if (!state) return;
        // Cleanup is authorized by the independently minted identity, never
        // by the current metadata alone. A tampered/misbound copy cannot erase
        // its authoritative parent or another live private branch.
        if (erase && !ownerId.empty() && !ownedId.empty() && ownedId != ownerId &&
            ownedId.rfind(ownerId + "/branch/", 0) == 0 &&
            zhangPrivateGraphIdentity(*state) == ownedId)
            erase(*state);
        delete state;
    }
};

using ZhangPrivateGraphState = std::unique_ptr<KFState, ZhangPrivateGraphStateDeleter>;
inline std::atomic<std::uint64_t> zhangPrivateGraphNextLease{0};

/** Bind before cloning or conditioning. Failed cloning, failed conditioning,
 * early return and normal destruction all release only this private lease.
 * The clone callback must leave its source unchanged and report failure. */
template<class Clone>
bool zhangMakePrivateGraphState(
    const KFState& source, const std::string& expectedOwnerId,
    const std::string& role, Clone&& clone, void (*erase)(const KFState&),
    ZhangPrivateGraphState& destination, std::string& reason)
{
    const auto validId = [](const std::string& id)
    {
        if (id.empty() || id.size() > 512) return false;
        for (unsigned char ch : id) if (ch < 0x21 || ch > 0x7e) return false;
        return true;
    };
    if (!validId(expectedOwnerId) || zhangPrivateGraphIdentity(source) != expectedOwnerId ||
        !validId(role) || role.find('/') != std::string::npos || !erase)
    {
        reason = "PRIVATE_GRAPH_OWNER_OR_ROLE_INVALID";
        return false;
    }
    // Never recycle a role-only branch ID while another prepared publication
    // may still own it. This identity is process-local, not checkpoint state.
    const std::string privateId = expectedOwnerId + "/branch/" + role + "-" +
        std::to_string(++zhangPrivateGraphNextLease);
    if (!validId(privateId))
    {
        reason = "PRIVATE_GRAPH_ID_INVALID";
        return false;
    }
    ZhangPrivateGraphState candidate(new KFState(source),
        ZhangPrivateGraphStateDeleter{expectedOwnerId, privateId, erase});
    candidate->metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA] = privateId;
    if (!clone(source, *candidate, reason)) return false;
    destination = std::move(candidate);
    reason = "PRIVATE_GRAPH_BOUND_AND_CLONED";
    return true;
}
