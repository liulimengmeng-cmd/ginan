#pragma once

#include "common/zhangProductFloatCarrierMetadata.hpp"
#include <cstddef>

// Observations of this one directory transaction, not a process resource guard.
struct ZhangFloatTargetBatchStats
{
    std::size_t requestedTargets=0, sourceTargets=0, finalTargets=0;
    std::size_t sourceMetadataEntries=0, finalMetadataEntries=0;
    std::size_t stateOrderBytes=0, duplicateRequests=0;
    std::size_t sourceValidatorCalls=0, sourceDirectoryDecodeCalls=0, finalValidatorCalls=0;
    std::size_t registeredBound=0, registeredUnbound=0, resolverCalls=0;
    std::size_t requestWireValidationCalls=0;
    std::size_t serializedTargets=0, serializedTargetBytes=0, atomicSwaps=0;
};

// Transactional equivalent of ClearGraphProductHandoffTargets + the ordered
// single-target Register loop, on the domain whose final directory is valid.
// The old single-target Register/Resolve/MetadataValid are unchanged oracles.
// Invalid final wire/schema states (including the old last-insert 8193 hole)
// are rejected before the only metadata swap. No x/P/dx/integer owner changes.
bool zhangRefreshGraphProductHandoffTargetsBatch(KFState& state,
    const std::vector<ZhangProductFloatTargetRequest>& requested,
    std::string& reason, ZhangFloatTargetBatchStats* statistics=nullptr);
