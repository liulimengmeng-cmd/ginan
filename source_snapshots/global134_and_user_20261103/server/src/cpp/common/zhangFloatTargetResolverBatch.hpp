#pragma once

#include "common/zhangProductFloatCarrierMetadata.hpp"
#include <cstddef>

// Observational work counters, not acceptance/resource thresholds.
struct ZhangFloatTargetResolverBatchStats
{
    std::size_t requestedTargets=0, sourceTargets=0;
    std::size_t sourceValidatorCalls=0, explicitDirectoryDecodeCalls=0;
    std::size_t stateOrderBuilds=0, indexedTargets=0, sourceGuardChecks=0;
    std::size_t registeredMatches=0, spanQueries=0, physicalLookups=0;
    std::size_t resolvedTargets=0, unboundTargets=0;
};

struct ZhangFloatTargetResolution
{
    bool resolved=false;
    std::map<KFKey,double> row;
    double offset=0;
    std::string reason;
};

// The validated immutable context is deliberately PRIVATE to this synchronous
// call. It cannot escape into a caller, survive a transform/retirement/metadata
// swap, or be reused by epoch/root name. No callbacks are accepted. Source is
// const throughout; exact source/graph checks precede the only output commit.
//
// Return false is a fatal directory/scope failure: both outputs are empty.
// Return true supplies one result for EVERY ordered request, including ordinary
// UNKNOWN results with the unchanged single-resolver reason/row/offset. Repeated
// requests are not deduplicated. registeredDirectory is the validated ORIGINAL
// map-order directory, not a directory changed by queries. No x/P/dx, metadata,
// current integer evidence or historical integer RHS is written.
// sourceValidatorCalls counts the unchanged old validator once; its INTERNAL
// codec passes are not claimed to be one parse. explicitDirectoryDecodeCalls
// counts the one additional directory materialisation used by this context.
bool zhangResolveGraphProductHandoffTargetsBatch(const KFState& state,
    const std::vector<ZhangProductFloatTargetRequest>& requests,
    std::vector<ZhangFloatTargetResolution>& output,std::string& reason,
    ZhangFloatTargetResolverBatchStats* statistics=nullptr,
    std::vector<ZhangProductFloatTargetRequest>* registeredDirectory=nullptr);
