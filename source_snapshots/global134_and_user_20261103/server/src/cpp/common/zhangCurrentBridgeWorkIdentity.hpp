#pragma once

#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangMultiGnssInteger.hpp"

#include <tuple>

// Work/result deduplication must bind the COMPLETE numerical conditioning
// recipe, not just this constellation's integer generation or a row count.
// Equality is exact content equality (no floating tolerance or hash-only key).
using ZhangCurrentBridgeConditioningKey = std::vector<std::tuple<E_Sys,
    std::uint64_t, std::vector<std::pair<ZhangCoordinateRow, ZhangExactInteger>>,
    std::string,std::string>>;

inline bool zhangCurrentBridgeEvidenceMatchesScope(
    const ZhangCurrentEpochEvidenceKey& scope,const ZhangCurrentEpochEvidence& evidence,
    const std::string& runtime,long int epoch,const std::string& posterior,
    const std::string& chart,const std::string& root,const std::string& order)
{
    if(scope.first!=runtime || evidence.runtimeId!=runtime || evidence.system!=scope.second ||
        evidence.floatStateRoot!=root || evidence.floatPosteriorPsdStateOrderId!=order)
        return false;
    const auto& authorized=zhangMultiGnssCurrentEvidenceAuthorities();
    const auto joint=authorized.find(runtime);
    if(joint==authorized.end())return evidence.matches(epoch,posterior,chart);
    const auto& authority=joint->second;
    if(authority.epochGpstSeconds!=epoch || authority.floatRoot!=root ||
        authority.floatStateOrderId!=order)return false;
    const auto system=authority.systemIdentities.find(scope.second);
    return system!=authority.systemIdentities.end() && evidence.matches(epoch,
        system->second.first,system->second.second);
}

template<class PAIR_KEY>
inline PAIR_KEY zhangCurrentBridgeRekey(PAIR_KEY key,std::uint64_t generation,
    const ZhangCurrentBridgeConditioningKey& conditioner)
{
    std::get<3>(key)=generation;
    std::get<4>(key)=conditioner;
    return key;
}

inline ZhangCurrentBridgeConditioningKey zhangCurrentBridgeConditioningKey(
    const std::map<ZhangCurrentEpochEvidenceKey, ZhangCurrentEpochEvidence>& registry,
    const std::string& runtime, long int epoch, const std::string& posterior,
    const std::string& chart, const std::string& root, const std::string& order)
{
    ZhangCurrentBridgeConditioningKey result;
    for (const auto& [scope, evidence] : registry)
    {
        if (!zhangCurrentBridgeEvidenceMatchesScope(scope,evidence,runtime,epoch,
            posterior,chart,root,order))
            continue;
        std::vector<std::pair<ZhangCoordinateRow, ZhangExactInteger>> recipe;
        recipe.reserve(evidence.numericConstraints.size());
        for (const auto& condition : evidence.numericConstraints)
            recipe.emplace_back(condition.rowCycles, condition.value);
        result.emplace_back(scope.second,evidence.constraintGeneration,std::move(recipe),
            evidence.posteriorId,evidence.chartId);
    }
    return result;
}
