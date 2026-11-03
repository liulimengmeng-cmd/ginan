#pragma once

#include "common/zhangProductFloatCarriers.hpp"
#include "common/zhangExactFloatTargetSpan.hpp"
#include "common/zhangCheckpoint.hpp"
#include <Eigen/SparseQR>
#include <iomanip>
#include <locale>
#include <sstream>

inline constexpr const char* ZHANG_FLOAT_TARGET_V2_PREFIX = "ZHANG_FLOAT_TARGET_V2/";
inline constexpr const char* ZHANG_FLOAT_CARRIER_V2_PREFIX = "ZHANG_FLOAT_CARRIER_V2/";
inline constexpr const char* ZHANG_FLOAT_TARGET_SCOPE_V2 = "ZHANG_FLOAT_TARGET_SCOPE_V2";
inline constexpr const char* ZHANG_FLOAT_TARGET_GENERATION_V2 = "ZHANG_FLOAT_TARGET_GENERATION_V2";
inline constexpr const char* ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2 = "ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2";
inline constexpr const char* ZHANG_FLOAT_TARGET_COMPRESSION_LABEL =
    "ZHANG_PRODUCT_HANDOFF_TARGET_COMPRESSION";

inline std::string zhangProductFloatOwnerId(const KFState& state)
{
    const auto branch=state.metaDataMap.find(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA);
    if(branch!=state.metaDataMap.end() && !branch->second.empty())return branch->second;
    const auto root=state.metaDataMap.find(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);
    return root==state.metaDataMap.end()?std::string{}:root->second;
}
inline bool zhangIsProductFloatCarrierKey(const KFKey& key)
{return key.type==KF::PRODUCT_HANDOFF && key.str.starts_with("ZHANG_FLOAT_CARRIER:");}
inline std::string zhangProductFloatCarrierRecordKey(const KFKey& key)
{return std::string(ZHANG_FLOAT_CARRIER_V2_PREFIX)+static_cast<std::string>(key);}

inline std::string zhangFloatCarrierStateOrder(const std::map<KFKey,int>& indices)
{
    std::vector<std::string> order(indices.size());
    std::vector<bool> seen(indices.size());
    for (const auto& [key,index] : indices)
    {
        if (index < 0 || index >= order.size() || seen[index]) return {};
        order[index] = static_cast<std::string>(key); seen[index] = true;
    }
    std::ostringstream out; out << "R51_FULL_STATE_ORDER_V1|";
    for (std::size_t index=0;index<order.size();++index)
        out << index << '=' << order[index].size() << ':' << order[index] << ';';
    return out.str();
}

inline void zhangWriteFloatTargetKey(std::ostream& out,const KFKey& key)
{
    out << int(key.type) << ' ' << std::quoted(key.Sat.id()) << ' '
        << std::quoted(key.str) << ' ' << key.num << ' ' << std::quoted(key.comment);
}
inline bool zhangReadFloatTargetKey(std::istream& in,KFKey& key)
{
    int type;std::string satellite;
    if (!(in >> type >> std::quoted(satellite) >> std::quoted(key.str)
          >> key.num >> std::quoted(key.comment))) return false;
    key.type=static_cast<KF>(type);key.Sat=SatSys(satellite.c_str());
    return key.Sat.id()==satellite;
}
inline std::string zhangEncodeFloatTarget(const ZhangProductFloatTargetRequest& target)
{
    std::ostringstream out;out.imbue(std::locale::classic());
    out << "FLOAT_TARGET_V2 " << std::quoted(target.targetId) << ' '
        << std::quoted(target.consumer) << ' ' << std::quoted(target.lifecycleId) << ' '
        << std::quoted(target.unit) << ' ' << int(target.symbolicSpace) << ' '
        << std::quoted(target.scope.runtimeId) << ' ' << std::quoted(target.scope.chartId)
        << ' ' << target.scope.coordinateGeneration << ' ' << target.active << ' '
        << target.numericBindingComplete << ' ' << std::quoted(target.bindingReason) << ' '
        << std::quoted(target.symbolicOffset.str()) << ' ' << std::quoted(target.physicalScale.str()) << ' ' << std::hexfloat
        << std::setprecision(std::numeric_limits<double>::max_digits10)
        << target.numericOffset << std::defaultfloat << ' ' << target.physicalCoefficients.size();
    for (const auto& [term,coefficient] : target.physicalCoefficients)
        out << ' ' << std::quoted(term.edge.receiver) << ' '
            << std::quoted(term.edge.satellite.id()) << ' ' << int(term.observable) << ' '
            << term.arcVersion << ' ' << std::quoted(coefficient.str());
    out << ' ' << target.nuisanceCoefficients.size();
    for (const auto& [id,value] : target.nuisanceCoefficients)
        out << ' ' << std::quoted(id) << ' ' << std::quoted(value.str());
    out << ' ' << target.numericRow.size();
    for (const auto& [key,value] : target.numericRow)
    {
        out << ' ';zhangWriteFloatTargetKey(out,key);
        out << ' ' << std::hexfloat << std::setprecision(std::numeric_limits<double>::max_digits10)
            << value << std::defaultfloat;
    }
    return out.str();
}
inline bool zhangDecodeFloatTarget(const std::string& text,ZhangProductFloatTargetRequest& target)
{
    std::istringstream in(text);in.imbue(std::locale::classic());
    std::string magic,offset,scale,numeric;int space;std::size_t count;
    if (!(in >> magic >> std::quoted(target.targetId) >> std::quoted(target.consumer)
        >> std::quoted(target.lifecycleId) >> std::quoted(target.unit) >> space
        >> std::quoted(target.scope.runtimeId) >> std::quoted(target.scope.chartId)
        >> target.scope.coordinateGeneration >> target.active >> target.numericBindingComplete
        >> std::quoted(target.bindingReason) >> std::quoted(offset) >> std::quoted(scale) >> numeric) ||
        magic!="FLOAT_TARGET_V2" || target.targetId.empty() || target.unit.empty() ||
        space<0 || space>int(ZhangProductFloatTargetSpace::SOURCE_KF_COORDINATE)) return false;
    target.symbolicSpace=static_cast<ZhangProductFloatTargetSpace>(space);
    auto number=[](const std::string& token,double& value) {
        try {std::size_t used=0;value=std::stod(token,&used);
            return used==token.size() && std::isfinite(value);} catch (...) {return false;}};
    try {target.symbolicOffset=ZhangFloatCarrierRational(offset);target.physicalScale=ZhangFloatCarrierRational(scale);}
    catch (...) {return false;}
    if (!number(numeric,target.numericOffset) || !(in>>count) || count>65536) return false;
    for (std::size_t i=0;i<count;++i)
    {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_DECODE_PHYSICAL_TERM");
        ZhangProductFloatVersionedTerm term;std::string satellite,value;int observable;
        if (!(in >> std::quoted(term.edge.receiver) >> std::quoted(satellite) >> observable
            >> term.arcVersion >> std::quoted(value)) || term.arcVersion<0) return false;
        term.edge.satellite=SatSys(satellite.c_str());term.observable=static_cast<E_ObsCode>(observable);
        if (term.edge.satellite.prn<=0 || term.edge.satellite.id()!=satellite) return false;
        try {if(!target.physicalCoefficients.emplace(term,ZhangExactInteger(value)).second)return false;}
        catch (...) {return false;}
    }
    if (!(in>>count) || count>65536) return false;
    for (std::size_t i=0;i<count;++i)
    {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_DECODE_NUISANCE_TERM");
        std::string id,value;if(!(in>>std::quoted(id)>>std::quoted(value)))return false;
        try {if(!target.nuisanceCoefficients.emplace(id,ZhangFloatCarrierRational(value)).second)return false;}
        catch (...) {return false;}
    }
    if (!(in>>count) || count>65536) return false;
    for (std::size_t i=0;i<count;++i)
    {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_DECODE_NUMERIC_TERM");
        KFKey key;double value;
        if(!zhangReadFloatTargetKey(in,key) || !(in>>numeric) || !number(numeric,value) ||
            !target.numericRow.emplace(key,value).second)return false;
    }
    in>>std::ws;return in.eof();
}
inline bool zhangReadFloatTargets(const std::map<std::string,std::string>& metadata,
    std::vector<ZhangProductFloatTargetRequest>& targets,std::string& reason)
{
    targets.clear();
    for(const auto& [name,value]:metadata)if(name.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX))
    {
        if(targets.size()>=8192){reason="FLOAT_TARGET_DIRECTORY_CAPACITY";return false;}
        ZhangProductFloatTargetRequest target;
        if(!zhangDecodeFloatTarget(value,target) ||
            name!=std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId)
        {reason="FLOAT_TARGET_METADATA_INVALID";return false;}
        targets.push_back(std::move(target));
    }
    const auto scope=metadata.find(ZHANG_FLOAT_TARGET_SCOPE_V2);
    for(auto& target:targets)
        if(scope!=metadata.end())target.scope.stateOrderId=scope->second;
    reason="FLOAT_TARGET_DIRECTORY_PARSED";return true;
}

inline bool hasZhangFloatTargetV2State(const std::map<std::string,std::string>& metadata)
{
    for(const auto& [name,value]:metadata)
        if(name.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX) ||
            name==ZHANG_FLOAT_TARGET_SCOPE_V2 || name==ZHANG_FLOAT_TARGET_GENERATION_V2 ||
            name==ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2)return true;
    return false;
}
inline bool hasZhangFloatTargetV2Records(const std::map<std::string,std::string>& metadata)
{
    for(const auto& [name,value]:metadata)
        if(name.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX))return true;
    return false;
}
inline bool readZhangFloatTargetSavedFactorSequence(const std::map<std::string,std::string>& metadata,
    const std::string& expectedOwner,const std::string& expectedOrder,
    std::uint64_t& savedSequence,std::string& reason)
{
    savedSequence=0;
    bool records=false;
    for(const auto& [name,value]:metadata)
        records |= name.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX);
    if(!records){reason=hasZhangFloatTargetV2State(metadata)?"ORPHAN_V2_FLOAT_TARGET_SCHEMA":"NO_V2_FLOAT_TARGET_STATE";return false;}
    auto decimal=[&](const char* name,std::uint64_t& value) {
        const auto entry=metadata.find(name);if(entry==metadata.end() || entry->second.empty())return false;
        try {std::size_t used=0;value=std::stoull(entry->second,&used);
            return used==entry->second.size() && std::to_string(value)==entry->second;}catch(...){return false;}};
    std::uint64_t generation;
    const auto order=metadata.find(ZHANG_FLOAT_TARGET_SCOPE_V2);
    if(expectedOwner.empty() || expectedOrder.empty() || order==metadata.end() || order->second!=expectedOrder ||
        !decimal(ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2,savedSequence) ||
        !decimal(ZHANG_FLOAT_TARGET_GENERATION_V2,generation) || generation>savedSequence)
    {reason="FLOAT_TARGET_SAVED_SEQUENCE_SCOPE_OR_DECIMAL_INVALID";return false;}
    std::vector<ZhangProductFloatTargetRequest> targets;
    if(!zhangReadFloatTargets(metadata,targets,reason))return false;
    for(const auto& target:targets)
        if(target.scope.runtimeId!=expectedOwner || target.scope.coordinateGeneration>savedSequence ||
            (target.numericBindingComplete && (target.scope.coordinateGeneration!=generation || target.scope.chartId.empty())))
        {reason="FLOAT_TARGET_SAVED_SEQUENCE_OWNER_OR_GENERATION_MISMATCH";return false;}
    for(const auto& [name,value]:metadata)if(name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX))
    {
        ZhangProductFloatTargetRequest definition;
        if(!zhangDecodeFloatTarget(value,definition) || definition.consumer!="FLOAT_CARRIER_SOURCE_DEFINITION" ||
            definition.scope.runtimeId!=expectedOwner || definition.scope.coordinateGeneration>savedSequence)
        {reason="FLOAT_CARRIER_SAVED_SEQUENCE_DEFINITION_INVALID";return false;}
    }
    reason="FLOAT_TARGET_ACTUAL_FACTOR_SEQUENCE_VALIDATED";return true;
}

inline bool zhangProductFloatCarrierCoordinateValid(const KFState& state,
    const std::map<std::string,std::string>& metadata,const KFKey& key,std::string& reason)
{
    if(!zhangIsProductFloatCarrierKey(key) || !state.kfIndexMap.contains(key))
    {reason="FLOAT_CARRIER_KEY_OR_STATE_MISSING";return false;}
    const auto record=metadata.find(zhangProductFloatCarrierRecordKey(key));
    ZhangProductFloatTargetRequest definition;
    if(record==metadata.end() || !zhangDecodeFloatTarget(record->second,definition) ||
        definition.consumer!="FLOAT_CARRIER_SOURCE_DEFINITION" ||
        definition.scope.runtimeId!=zhangProductFloatOwnerId(state))
    {reason="FLOAT_CARRIER_TAGGED_DEFINITION_OR_OWNER_INVALID";return false;}
    const auto model=state.stateTransitionMap.find(key);
    if(model==state.stateTransitionMap.end() || model->second.size()!=1 ||
        !model->second.contains(key) || model->second.at(key).size()!=1 ||
        !model->second.at(key).contains(0) || model->second.at(key).at(0)!=1 ||
        (state.procNoiseMap.contains(key) && state.procNoiseMap.at(key)!=0) ||
        (state.gaussMarkovTauMap.contains(key) && state.gaussMarkovTauMap.at(key)>0))
    {reason="FLOAT_CARRIER_NONCONSTANT_MODEL";return false;}
    reason="TYPED_CORRELATED_FLOAT_CARRIER_VALID";return true;
}

inline bool zhangProductFloatTargetMetadataValid(const KFState& state,std::string& reason)
{
    const auto owner=zhangProductFloatOwnerId(state);
    const auto stateOrder=zhangFloatCarrierStateOrder(state.kfIndexMap);
    if(hasZhangFloatTargetV2State(state.metaDataMap))
    {
        std::uint64_t saved;
        if(!readZhangFloatTargetSavedFactorSequence(state.metaDataMap,owner,
            stateOrder,saved,reason) || saved!=state.factorCommitSequence)
        {if(reason=="FLOAT_TARGET_ACTUAL_FACTOR_SEQUENCE_VALIDATED")reason="FLOAT_TARGET_SAVED_SEQUENCE_NOT_ACTUAL";return false;}
    }
    std::vector<ZhangProductFloatTargetRequest> targets;
    if(!zhangReadFloatTargets(state.metaDataMap,targets,reason))return false;
    std::optional<std::uint64_t> generation;
    if(const auto entry=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);entry!=state.metaDataMap.end())
        try {std::size_t used=0;const auto value=std::stoull(entry->second,&used);
            if(used==entry->second.size())generation=value;} catch (...) {}
    for(const auto& target:targets)
    {
        if(target.scope.runtimeId!=owner)
        {reason="FLOAT_TARGET_OWNER_MISMATCH";return false;}
        if(!target.numericBindingComplete)continue;
        if(!generation || target.scope.coordinateGeneration!=*generation || *generation>state.factorCommitSequence)
        {reason="FLOAT_TARGET_COMPLETE_GENERATION_MISMATCH";return false;}
        if(target.scope.stateOrderId!=stateOrder)
        {reason="FLOAT_TARGET_COMPLETE_ORDER_MISMATCH";return false;}
        for(const auto& [key,value]:target.numericRow)
            if(!std::isfinite(value) || !state.kfIndexMap.contains(key))
            {reason="FLOAT_TARGET_BOUND_KEY_OR_COEFFICIENT_INVALID";return false;}
    }
    std::set<std::string> records;
    for(const auto& [key,index]:state.kfIndexMap)if(zhangIsProductFloatCarrierKey(key))
    {
        if(!zhangProductFloatCarrierCoordinateValid(state,state.metaDataMap,key,reason))return false;
        records.insert(zhangProductFloatCarrierRecordKey(key));
    }
    for(const auto& [name,value]:state.metaDataMap)
        if(name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX) && !records.contains(name))
        {reason="FLOAT_CARRIER_ORPHAN_DEFINITION";return false;}
    reason="WHOLE_FLOAT_TARGET_AND_CARRIER_METADATA_VALID";return true;
}

/** Validate one already-decoded directory without materialising the shared
 * (potentially very large) state-order string in every target. This has the
 * same owner/generation/model/key gates as the legacy reader above. It is
 * preparation only, never an integer/publication authority. The caller must
 * decode the literal directory entries once and keep their immutable owner. */
inline bool zhangValidateDecodedProductFloatTargetMetadata(
    const KFState& state,const std::vector<ZhangProductFloatTargetRequest>& targets,
    const std::map<std::string,ZhangProductFloatTargetRequest>& carriers,
    const std::string& stateOrder,std::string& reason)
{
    ZhangPhaseTimer timer("BRIDGE_METADATA_SCOPE_AND_MODEL_AUDIT");
    const auto owner=zhangProductFloatOwnerId(state);
    std::size_t targetCount=0,carrierCount=0;
    for(const auto& [name,value]:state.metaDataMap) {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_DIRECTORY_COVERAGE");
        targetCount+=name.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX);
        carrierCount+=name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX);
    }
    if(targetCount!=targets.size() || carrierCount!=carriers.size())
    {reason="FLOAT_TARGET_DECODED_DIRECTORY_COVERAGE_MISMATCH";return false;}
    std::optional<std::uint64_t> generation;
    if(hasZhangFloatTargetV2State(state.metaDataMap)) {
        if(!targetCount && !carrierCount){reason="ORPHAN_V2_FLOAT_TARGET_SCHEMA";return false;}
        auto decimal=[&](const char* name,std::uint64_t& result) {
            const auto entry=state.metaDataMap.find(name);
            if(entry==state.metaDataMap.end() || entry->second.empty())return false;
            try {std::size_t used=0;result=std::stoull(entry->second,&used);
                return used==entry->second.size() && std::to_string(result)==entry->second;
            }catch(...){return false;}
        };
        std::uint64_t saved=0,current=0;
        const auto order=state.metaDataMap.find(ZHANG_FLOAT_TARGET_SCOPE_V2);
        if(owner.empty() || stateOrder.empty() || order==state.metaDataMap.end() ||
            order->second!=stateOrder || !decimal(ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2,saved) ||
            !decimal(ZHANG_FLOAT_TARGET_GENERATION_V2,current) || current>saved)
        {reason="FLOAT_TARGET_SAVED_SEQUENCE_SCOPE_OR_DECIMAL_INVALID";return false;}
        if(saved!=state.factorCommitSequence)
        {reason="FLOAT_TARGET_SAVED_SEQUENCE_NOT_ACTUAL";return false;}
        generation=current;
    }
    std::set<std::string> targetIds;
    for(const auto& target:targets) {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_TARGET_SCOPE_AUDIT");
        if(!targetIds.insert(target.targetId).second ||
            !state.metaDataMap.contains(std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId))
        {reason="FLOAT_TARGET_DECODED_IDENTITY_MISMATCH";return false;}
        if(target.scope.runtimeId!=owner ||
            target.scope.coordinateGeneration>state.factorCommitSequence ||
            (target.numericBindingComplete && (!generation ||
                target.scope.coordinateGeneration!=*generation || target.scope.chartId.empty())))
        {reason="FLOAT_TARGET_SAVED_SEQUENCE_OWNER_OR_GENERATION_MISMATCH";return false;}
        if(!target.numericBindingComplete)continue;
        if(!target.scope.stateOrderId.empty() && target.scope.stateOrderId!=stateOrder)
        {reason="FLOAT_TARGET_COMPLETE_ORDER_MISMATCH";return false;}
        for(const auto& [key,value]:target.numericRow) {
            zhangExactWorkCheckpoint("BRIDGE_METADATA_NUMERIC_KEY_AUDIT");
            if(!std::isfinite(value) || !state.kfIndexMap.contains(key))
            {reason="FLOAT_TARGET_BOUND_KEY_OR_COEFFICIENT_INVALID";return false;}
        }
    }
    std::set<std::string> records;
    for(const auto& [key,index]:state.kfIndexMap)if(zhangIsProductFloatCarrierKey(key)) {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_CARRIER_MODEL_AUDIT");
        const auto name=zhangProductFloatCarrierRecordKey(key);
        const auto found=carriers.find(name);
        if(found==carriers.end()) {reason="FLOAT_CARRIER_TAGGED_DEFINITION_OR_OWNER_INVALID";return false;}
        const auto& definition=found->second;
        if(definition.consumer!="FLOAT_CARRIER_SOURCE_DEFINITION" ||
            definition.scope.runtimeId!=owner ||
            definition.scope.coordinateGeneration>state.factorCommitSequence)
        {reason="FLOAT_CARRIER_SAVED_SEQUENCE_DEFINITION_INVALID";return false;}
        const auto model=state.stateTransitionMap.find(key);
        if(model==state.stateTransitionMap.end() || model->second.size()!=1 ||
            !model->second.contains(key) || model->second.at(key).size()!=1 ||
            !model->second.at(key).contains(0) || model->second.at(key).at(0)!=1 ||
            (state.procNoiseMap.contains(key) && state.procNoiseMap.at(key)!=0) ||
            (state.gaussMarkovTauMap.contains(key) && state.gaussMarkovTauMap.at(key)>0))
        {reason="FLOAT_CARRIER_NONCONSTANT_MODEL";return false;}
        records.insert(name);
    }
    for(const auto& [name,definition]:carriers) {
        zhangExactWorkCheckpoint("BRIDGE_METADATA_CARRIER_DIRECTORY_AUDIT");
        if(!records.contains(name)) {reason="FLOAT_CARRIER_ORPHAN_DEFINITION";return false;}
    }
    reason="WHOLE_FLOAT_TARGET_AND_CARRIER_METADATA_VALID";return true;
}


// Staged remap used by EVERY real KF coordinate transform. The source row is
// resolved as a whole against T; no individual missing source is assigned 0.
// An active already-bound functional may never silently become UNBOUND.
// Graph callers preserve its direction before projection; an unprotected
// destructive map is rejected here before any x/P or callback mutation.
inline bool transportZhangProductFloatTargetMetadata(const KFState& state,
    const std::map<KFKey,int>& destination,const SparseMatrix<double>& transform,
    const SparseMatrix<double>& noiseInjection,const std::string& label,
    std::uint64_t afterSequence,std::map<std::string,std::string>& metadata,
    std::string& reason,const std::map<std::string,std::string>* metadataOverride=nullptr,
    const std::map<std::string,std::map<KFKey,double>>* targetDestinationHints=nullptr)
{
    if(hasZhangFloatTargetV2State(state.metaDataMap) && !zhangProductFloatTargetMetadataValid(state,reason))return false;
    if(metadataOverride)
    {
        auto mutableHandoff=[](const std::string& name) {
            return name.starts_with("ZHANG_HANDOFF_FLOAT_") ||
                name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX);};
        for(const auto& [name,value]:state.metaDataMap)
            if(!mutableHandoff(name))
            {
                const auto proposed=metadataOverride->find(name);
                if(proposed==metadataOverride->end() || proposed->second!=value)
                {reason="FLOAT_TARGET_METADATA_OVERRIDE_CHANGED_SOURCE_OWNER_OR_CORE";return false;}
            }
        for(const auto& [name,value]:*metadataOverride)
            if(!mutableHandoff(name) && !state.metaDataMap.contains(name))
            {reason="FLOAT_TARGET_METADATA_OVERRIDE_INJECTED_NONHANDOFF_DATA";return false;}
    }
    metadata=metadataOverride?*metadataOverride:state.metaDataMap;
    std::vector<ZhangProductFloatTargetRequest> targets;
    // A caller cannot inject or replace the source numeric rows via override.
    if(!zhangReadFloatTargets(state.metaDataMap,targets,reason))return false;
    if(targetDestinationHints)
        for(const auto& [id,hint]:*targetDestinationHints)
        {
            if(std::none_of(targets.begin(),targets.end(),[&](const auto& target){return target.targetId==id;}))
            {reason="FLOAT_TARGET_DESTINATION_HINT_UNKNOWN_TARGET";return false;}
            for(const auto& [key,value]:hint)
                if(!destination.contains(key) || !std::isfinite(value))
                {reason="FLOAT_TARGET_DESTINATION_HINT_INVALID_KEY_OR_COEFFICIENT";return false;}
        }
    if(!hasZhangFloatTargetV2State(state.metaDataMap))return true;
    metadata[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(afterSequence);
    if(targets.empty())
    {
        if(!hasZhangFloatTargetV2Records(metadata))
        {
            metadata.erase(ZHANG_FLOAT_TARGET_SCOPE_V2);metadata.erase(ZHANG_FLOAT_TARGET_GENERATION_V2);
            metadata.erase(ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2);return true;
        }
        metadata[ZHANG_FLOAT_TARGET_SCOPE_V2]=zhangFloatCarrierStateOrder(destination);
        metadata[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(afterSequence);
        return true;
    }
    const bool strict=label==ZHANG_FLOAT_TARGET_COMPRESSION_LABEL;
    SparseMatrix<double> augmented=transform;
    if(noiseInjection.cols())
    {
        std::vector<Eigen::Triplet<double>> entries;
        entries.reserve(transform.nonZeros()+noiseInjection.nonZeros());
        for(int outer=0;outer<transform.outerSize();++outer)
            for(SparseMatrix<double>::InnerIterator e(transform,outer);e;++e)
                entries.emplace_back(e.row(),e.col(),e.value());
        for(int outer=0;outer<noiseInjection.outerSize();++outer)
            for(SparseMatrix<double>::InnerIterator e(noiseInjection,outer);e;++e)
                entries.emplace_back(e.row(),transform.cols()+e.col(),e.value());
        augmented.resize(transform.rows(),transform.cols()+noiseInjection.cols());
        augmented.setFromTriplets(entries.begin(),entries.end());
    }
    SparseMatrix<double> transposed=augmented.transpose();transposed.makeCompressed();
    std::vector<std::vector<std::pair<std::size_t,ZhangFloatCarrierRational>>> exactRows(augmented.rows());
    for(int outer=0;outer<augmented.outerSize();++outer)
        for(SparseMatrix<double>::InnerIterator entry(augmented,outer);entry;++entry)
        {
            const auto exact=zhangFloatCarrierExactDouble(entry.value());
            if(!exact){reason="FLOAT_TARGET_TRANSFORM_COEFFICIENT_NONFINITE";return false;}
            if(*exact!=0)exactRows[entry.row()].emplace_back(entry.col(),*exact);
        }
    std::optional<Eigen::SparseQR<SparseMatrix<double>,Eigen::COLAMDOrdering<int>>> qr;
    // Lifetime is exactly this immutable actual [T G] call. Components are
    // indexed/factored only after QR failed an existing exact pullback gate.
    std::optional<zhang_float_target_q_span::ImmutableComponentSpan> exactSpan;
    const auto order=zhangFloatCarrierStateOrder(destination);
    const auto generationEntry=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);
    std::optional<std::uint64_t> generation;
    if(generationEntry!=state.metaDataMap.end())try {
        std::size_t used=0;const auto value=std::stoull(generationEntry->second,&used);
        if(used==generationEntry->second.size())generation=value;
    } catch (...) {}
    for(auto& target:targets)
    {
        if(target.scope.runtimeId!=zhangProductFloatOwnerId(state))
        {reason="FLOAT_TARGET_SOURCE_OWNER_MISMATCH";return false;}
        if(strict && target.active && !target.numericBindingComplete)
        {reason="FLOAT_COMPRESSION_ACTIVE_TARGET_UNBOUND:"+target.targetId;return false;}
        if(!target.numericBindingComplete)continue;
        if(!generation || target.scope.coordinateGeneration!=*generation ||
            *generation>state.factorCommitSequence)
        {reason="FLOAT_TARGET_SOURCE_COORDINATE_GENERATION_INVALID";return false;}
        VectorXd old=VectorXd::Zero(transposed.rows()),next=VectorXd::Zero(destination.size());
        ZhangFloatCarrierSparseRow exactOld;
        if(target.scope.stateOrderId!=zhangFloatCarrierStateOrder(state.kfIndexMap))
        {reason="FLOAT_TARGET_SOURCE_ORDER_INVALID";return false;}
        bool keysPresent=true;
        for(const auto& [key,value]:target.numericRow)
        {
            const auto source=state.kfIndexMap.find(key);
            if(source==state.kfIndexMap.end())
            {reason="FLOAT_TARGET_SOURCE_KEY_INVALID";return false;}
            old(source->second)+=value;
            const auto exact=zhangFloatCarrierExactDouble(value);
            if(!exact){reason="FLOAT_TARGET_SOURCE_COEFFICIENT_NONFINITE";return false;}
            if(*exact!=0)exactOld[source->second]+=*exact;
            if(const auto kept=destination.find(key);kept!=destination.end())next(kept->second)+=value;
        }
        const bool hinted=targetDestinationHints && targetDestinationHints->contains(target.targetId);
        if(hinted)
        {
            next.setZero();
            for(const auto& [key,value]:targetDestinationHints->at(target.targetId))next(destination.at(key))+=value;
        }
        auto residual=[&]() {return (transposed*next-old).norm()/std::max(1.0,old.norm());};
        auto exactPullback=[&]()
        {
            ZhangFloatCarrierSparseRow pulled;
            for(int output=0;output<next.size();++output)
            {
                if(next(output)==0)continue;
                const auto weight=zhangFloatCarrierExactDouble(next(output));
                if(!weight)return false;
                for(const auto& [source,coefficient]:exactRows[output])
                {
                    const ZhangFloatCarrierRational value=pulled[source]+*weight*coefficient;
                    if(value==0)pulled.erase(source);else pulled[source]=value;
                }
            }
            // Fresh independent sources occupy appended columns in [T G];
            // exactOld has zero there. Any nonzero amount, however small,
            // therefore fails this equality rather than passing a tolerance.
            return pulled==exactOld;
        };
        double error=keysPresent?residual():std::numeric_limits<double>::infinity();
        std::string exactFallbackReason;
        if(hinted && (!std::isfinite(error) || error>1e-10 || !exactPullback()))
        {reason="FLOAT_TARGET_DESTINATION_HINT_FAILED_EXACT_PULLBACK";return false;}
        if(!hinted && keysPresent && (!std::isfinite(error) || error>1e-10 || !exactPullback()))
        {
            if(!qr){qr.emplace();qr->compute(transposed);}
            if(qr->info()==Eigen::Success){next=qr->solve(old);error=residual();}
        }
        if(!hinted && keysPresent && (!next.allFinite() || !std::isfinite(error) || error>1e-10 || !exactPullback()))
        {
            ZhangPhaseTimer exactTimer("FLOAT_TARGET_EXACT_SPAN_FALLBACK");
            if(!exactSpan)exactSpan.emplace(exactRows,augmented.cols());
            const auto recovered=exactSpan->resolve(exactOld);
            exactFallbackReason=recovered.reason;
            if(recovered.losslesslyEncodable)
            {
                next.setZero();
                for(const auto&[row,value]:recovered.coefficients)next(row)=value;
                error=residual();
            }
        }
        if(!keysPresent || !next.allFinite() || !std::isfinite(error) || error>1e-10 ||
            !exactPullback())
        {
            if(target.active){reason=(strict?"FLOAT_COMPRESSION_ACTIVE_TARGET_NOT_TRANSPORTABLE:":
                "FLOAT_ACTIVE_BOUND_TARGET_NOT_TRANSPORTABLE:")+target.targetId;return false;}
            target.numericBindingComplete=false;
            target.bindingReason="UNBOUND_AFTER_REAL_TRANSFORM:"+label+
                (exactFallbackReason.empty()?std::string():":"+exactFallbackReason);
        }
        else
        {
            target.numericRow.clear();
            for(const auto& [key,index]:destination)if(next(index)!=0)target.numericRow.emplace(key,next(index));
            target.scope.coordinateGeneration=afterSequence;
            target.scope.chartId="KF_REAL_TRANSFORM:"+label+":"+std::to_string(afterSequence);
            target.scope.stateOrderId=order;
            target.bindingReason="WHOLE_FLOAT_ROW_TRANSPORTED_WITHOUT_FRESH_NOISE";
        }
        metadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
    }
    metadata[ZHANG_FLOAT_TARGET_SCOPE_V2]=order;
    metadata[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(afterSequence);
    reason="FLOAT_TARGET_METADATA_REMAP_STAGED";return true;
}

// Time prediction is not an exact historical transport when a dependency
// has dynamics, reinitialisation or process noise. Only proven static identity
// rows may acquire the new index/order binding; otherwise retain a typed
// UNBOUND semantic request so optional compression keeps the raw state.
inline bool stageZhangFloatTargetTimeTransitionMetadata(const KFState& state,
    const std::map<KFKey,int>& destination,const SparseMatrix<double>& transition,
    const VectorXd& processDiagonal,std::uint64_t afterSequence,
    std::map<std::string,std::string>& metadata,std::string& reason)
{
    if(hasZhangFloatTargetV2State(state.metaDataMap) && !zhangProductFloatTargetMetadataValid(state,reason))return false;
    metadata=state.metaDataMap;
    std::vector<ZhangProductFloatTargetRequest> targets;
    if(!zhangReadFloatTargets(metadata,targets,reason))return false;
    if(!hasZhangFloatTargetV2State(state.metaDataMap))return true;
    metadata[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(afterSequence);
    std::vector<std::map<int,double>> rows(transition.rows());
    for(int outer=0;outer<transition.outerSize();++outer)
        for(SparseMatrix<double>::InnerIterator entry(transition,outer);entry;++entry)
            if(entry.value()!=0)rows[entry.row()].emplace(entry.col(),entry.value());
    const auto sourceOrder=zhangFloatCarrierStateOrder(state.kfIndexMap);
    const auto destinationOrder=zhangFloatCarrierStateOrder(destination);
    const auto generationEntry=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);
    std::optional<std::uint64_t> generation;
    if(generationEntry!=state.metaDataMap.end())try {
        std::size_t used=0;const auto value=std::stoull(generationEntry->second,&used);
        if(used==generationEntry->second.size())generation=value;
    } catch (...) {}
    for(auto& target:targets)
    {
        if(target.scope.runtimeId!=zhangProductFloatOwnerId(state))
        {reason="FLOAT_TARGET_TIME_SOURCE_OWNER_MISMATCH";return false;}
        if(!target.numericBindingComplete)continue;
        if(!generation || target.scope.coordinateGeneration!=*generation ||
            *generation>state.factorCommitSequence || target.scope.stateOrderId!=sourceOrder)
        {reason="FLOAT_TARGET_TIME_SOURCE_BINDING_INVALID";return false;}
        bool identity=true;
        for(const auto& [key,coefficient]:target.numericRow)
        {
            if(coefficient==0)continue;
            const auto old=state.kfIndexMap.find(key),next=destination.find(key);
            if(old==state.kfIndexMap.end() || next==destination.end()) {identity=false;break;}
            const auto& row=rows.at(next->second);
            if(row.size()!=1 || !row.contains(old->second) || row.at(old->second)!=1 ||
                next->second>=processDiagonal.size() || processDiagonal(next->second)!=0)
            {identity=false;break;}
        }
        if(!identity)
        {
            target.numericBindingComplete=false;
            target.bindingReason="UNBOUND_TIME_TRANSITION_NONIDENTITY_OR_RANDOM_SOURCE";
        }
        else
        {
            target.scope.coordinateGeneration=afterSequence;
            target.scope.stateOrderId=destinationOrder;
            target.bindingReason="STATIC_IDENTITY_TIME_TRANSITION_REBOUND";
        }
        metadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
    }
    metadata[ZHANG_FLOAT_TARGET_SCOPE_V2]=destinationOrder;
    metadata[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(afterSequence);
    reason="FLOAT_TARGET_TIME_TRANSITION_METADATA_STAGED";return true;
}

inline bool stageZhangFloatTargetStochasticMetadata(const KFState& state,
    const std::map<KFKey,int>& destination,const SparseMatrix<double>& transition,
    const MatrixXd& processCovariance,std::uint64_t afterSequence,
    std::map<std::string,std::string>& metadata,std::string& reason)
{
    if(!hasZhangFloatTargetV2State(state.metaDataMap)){metadata=state.metaDataMap;return true;}
    VectorXd touched=VectorXd::Zero(processCovariance.rows());
    for(int column=0;column<processCovariance.cols();++column)
    for(int row=0;row<processCovariance.rows();++row)
        if(processCovariance(row,column)!=0)touched(row)=1;
    return stageZhangFloatTargetTimeTransitionMetadata(state,destination,transition,
        touched,afterSequence,metadata,reason);
}

inline bool stageZhangFloatTargetMeasurementMetadata(const KFState& state,
    const std::vector<KFStagedFactor>* factors,std::uint64_t afterSequence,
    std::map<std::string,std::string>& metadata,std::string& reason)
{
    metadata=state.metaDataMap;
    if(!hasZhangFloatTargetV2State(metadata))return true;
    if(!zhangProductFloatTargetMetadataValid(state,reason))return false;
    std::vector<ZhangProductFloatTargetRequest> targets;
    if(!zhangReadFloatTargets(metadata,targets,reason))return false;
    if(factors)for(const auto& factor:*factors)
    {
        if(factor.kind!=KFStagedFactorKind::STATE_TRANSITION)continue;
        std::vector<std::map<int,double>> rows(factor.transition.rows());
        for(int outer=0;outer<factor.transition.outerSize();++outer)
            for(SparseMatrix<double>::InnerIterator e(factor.transition,outer);e;++e)
                if(e.value()!=0)rows[e.row()].emplace(e.col(),e.value());
        std::vector<bool> noisy(factor.processCovariance.rows());
        for(int outer=0;outer<factor.processCovariance.outerSize();++outer)
            for(SparseMatrix<double>::InnerIterator e(factor.processCovariance,outer);e;++e)
                if(e.value()!=0)noisy[e.row()]=true;
        for(auto& target:targets)
        {
            if(!target.numericBindingComplete)continue;
            bool identity=true;
            for(const auto& [key,coefficient]:target.numericRow)
            {
                if(coefficient==0)continue;
                const auto source=factor.sourceIndexMap.find(key),destination=factor.destinationIndexMap.find(key);
                if(source==factor.sourceIndexMap.end() || destination==factor.destinationIndexMap.end() ||
                    destination->second>=rows.size() || destination->second>=noisy.size())
                {identity=false;break;}
                const auto& row=rows[destination->second];
                if(row.size()!=1 || !row.contains(source->second) || row.at(source->second)!=1 || noisy[destination->second])
                {identity=false;break;}
            }
            if(!identity)
            {
                target.numericBindingComplete=false;
                target.bindingReason="UNBOUND_BY_STAGED_QC_TRANSITION_NOT_AN_EXACT_OLD_TARGET";
            }
        }
    }
    for(const auto& target:targets)
        metadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
    // An accepted measurement changes posterior moments and the REAL factor
    // sequence, but leaves the coordinate binding generation untouched.
    metadata[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(afterSequence);
    reason="FLOAT_TARGET_MEASUREMENT_SEQUENCE_STAGED";return true;
}

inline void isolateZhangFloatTargetsForNonAuthoritativeReplay(KFState& state)
{
    for(auto item=state.metaDataMap.begin();item!=state.metaDataMap.end();)
        if(item->first.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || item->first.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX) ||
            item->first==ZHANG_FLOAT_TARGET_SCOPE_V2 || item->first==ZHANG_FLOAT_TARGET_GENERATION_V2 ||
            item->first==ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2)item=state.metaDataMap.erase(item);
        else ++item;
    state.lastPppTransactionCommitted=false;state.lastPppQcConverged=false;
}

bool zhangRegisterGraphProductHandoffTarget(KFState&,const ZhangProductFloatTargetRequest&,std::string&);
bool zhangResolveGraphProductHandoffTarget(const KFState&,const ZhangProductFloatTargetRequest&,
    std::map<KFKey,double>&,double&,std::string&);
void zhangClearGraphProductHandoffTargets(KFState&);
bool zhangCollectActiveProductFloatConsumerTargets(const KFState&,
    std::vector<ZhangProductFloatTargetRequest>&,std::string&);
bool zhangCompressGraphProductFloatTargets(std::ostream&,KFState&,std::string&);
