#pragma once
#include <map>
#include <set>
#include <string>

// Only these Gaussian-target/transport bookkeeping fields may advance at a
// product publication boundary. A directory is not a CURRENT integer recipe.
inline bool zhangFloatPublicationMetadataKey(const std::string& key)
{
    return key.starts_with("ZHANG_FLOAT_TARGET_V2/") ||
        key=="ZHANG_FLOAT_TARGET_SCOPE_V2" || key=="ZHANG_FLOAT_TARGET_GENERATION_V2" ||
        key=="ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2" ||
        key=="ZHANG_TRANSPORT_ONLY_COMMITTED_V1/GPS" ||
        key=="ZHANG_TRANSPORT_ONLY_COMMITTED_V1/GAL" ||
        key=="ZHANG_TRANSPORT_ONLY_COMMITTED_V1/BDS";
}

inline bool zhangValidateFloatPublicationMetadataDelta(
    const std::map<std::string,std::string>& before,
    const std::map<std::string,std::string>& after,
    const std::set<std::string>& authorizedKeys,
    const std::map<std::string,std::string>& expectedChangedEntries,
    std::string& reason)
{
    std::map<std::string,std::string> actual;
    for(const auto& [key,value]:before)
        if(!after.contains(key))
        {reason="FLOAT_METADATA_COMMIT_DELETED_EXISTING_KEY:"+key;return false;}
    for(const auto& [key,value]:after)
    {
        const auto old=before.find(key);
        if(old!=before.end() && old->second==value)continue;
        if(!zhangFloatPublicationMetadataKey(key) || !authorizedKeys.contains(key))
        {reason="FLOAT_METADATA_COMMIT_UNAUTHORISED_KEY:"+key;return false;}
        actual.emplace(key,value);
    }
    if(actual!=expectedChangedEntries)
    {reason="FLOAT_METADATA_COMMIT_CHANGED_VALUES_OR_MEMBERSHIP_DIFFER";return false;}
    reason="EXACT_PUBLICATION_METADATA_DELTA_VALID";return true;
}
