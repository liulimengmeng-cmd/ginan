// Actual-data COMPONENT test, not a checkpoint provenance override and not a
// nine-section process restore. All persisted FLOAT values come from the real
// independently exported core.payload; no synthetic x/P or zero filling.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangMultiGnssInteger.hpp"
#include "pea/zhangPppAr.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/json.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace j=boost::json;
namespace
{
constexpr const char* scope="REAL_GEC_SECTION_WITH_ACTUAL_KF_CORE_NOT_FULL_NINE_SECTION_RESTORE";
constexpr std::uint64_t realGecBytes=504018523;
constexpr const char* realCoreSha="6163ae25308fec882ad2f63c967e53cce7a2e52b6d90d0ea5d3c1901bcb43cdf";
constexpr const char* realGecSha="fb07adc9f827bd5d311dd50b35629d5cbdfaa16f76ef5b00df7c98530369cbde";
void require(bool value,const std::string& why)
{ if(!value)throw std::runtime_error(why); }
bool hashText(const std::string& value)
{ return value.size()==64 && value.find_first_not_of("0123456789abcdef")==std::string::npos; }
void fileHash(const std::string& path,const std::string& expected)
{
    std::string reason;
    const auto actual=zhangCheckpointFileSha256(path,&reason);
    require(hashText(expected) && actual==expected,"REAL_INPUT_HASH_MISMATCH:"+path+":"+reason);
}
void mode()
{
    setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
    setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);
    setenv("ZHANG_MULTI_GNSS_AR_CHECKPOINT","1",1);unsetenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT");
    acsConfig.zhangFullRank.enable=true;
    for(auto& [system,enabled]:acsConfig.process_sys)enabled=false;
    acsConfig.zhangFullRank.sysOpts.clear();
    acsConfig.process_sys[E_Sys::GPS]=acsConfig.process_sys[E_Sys::GAL]=acsConfig.process_sys[E_Sys::BDS]=true;
    acsConfig.zhangFullRank.sysOpts[E_Sys::GPS].baseline_observables={E_ObsCode::L1C,E_ObsCode::L2W};
    acsConfig.zhangFullRank.sysOpts[E_Sys::GAL].baseline_observables={E_ObsCode::L1C,E_ObsCode::L5Q};
    acsConfig.zhangFullRank.sysOpts[E_Sys::BDS].baseline_observables={E_ObsCode::L1P,E_ObsCode::L5P};
}
std::set<std::string> receiverIds(const ZhangCheckpointKfCore& core)
{
    std::set<std::string> ids;
    auto key=[&](const ZhangCheckpointKfKey& k)
    { if(!k.receiverPointerId.empty())ids.insert(k.receiverPointerId); };
    auto scalar=[&](const auto& rows){for(const auto& [k,value]:rows)key(k);};
    scalar(core.kfIndexMap);
    for(const auto& [destination,sources]:core.stateTransitionMap)
    { key(destination);for(const auto& [source,orders]:sources)key(source); }
    scalar(core.gaussMarkovTauMap);scalar(core.gaussMarkovMuMap);scalar(core.procNoiseMap);
    scalar(core.initNoiseMap);scalar(core.sigmaMaxMap);scalar(core.outageLimitMap);
    scalar(core.exponentialNoiseMap);scalar(core.errorCountMap);
    for(const auto& [destination,sources]:core.pseudoStateMap)
    { key(destination);for(const auto& [source,value]:sources)key(source); }
    for(const auto& [child,parent]:core.pseudoParentMap){key(child);key(parent);}
    return ids;
}
}

int main(int argc,char** argv)
{
    j::object report{{"scope",scope},{"passed",false},{"full_nine_section_restore",false},
        {"checkpoint_binary_sha_gate_bypassed",false},{"scientific_precision_accepted",false},
        {"stage","ARGUMENTS"},{"numerical_epochs_consumed",0}};
    try
    {
        require(argc==8,"usage: real_gec core.payload core_sha256 GEC.payload GEC_sha256 after.payload dimension expected_runtime");
        const std::string corePath=argv[1],coreSha=argv[2],gecPath=argv[3],gecSha=argv[4],outputPath=argv[5];
        const int dimension=std::stoi(argv[6]);const std::string expectedRuntime=argv[7];
        require(dimension==7262 && coreSha==realCoreSha && gecSha==realGecSha,
            "UNAPPROVED_REAL_E135_COMPONENT_IDENTITY");
        require(!std::filesystem::exists(outputPath),"FRESH_OUTPUT_FILE_REQUIRED");
        require(std::filesystem::file_size(gecPath)==realGecBytes,"WRONG_REAL_GEC_INPUT_LENGTH");
        fileHash(corePath,coreSha);fileHash(gecPath,gecSha);
        report["core_sha256"]=coreSha;report["source_GEC_sha256"]=gecSha;
        report["source_GEC_bytes"]=realGecBytes;report["stage"]="ACTUAL_CORE_DECODE";
        ZhangCheckpointKfCore core;
        {
            std::ifstream input(corePath,std::ios::binary);
            require(bool(input),"ACTUAL_CORE_OPEN_FAILED");
            boost::archive::binary_iarchive archive(input,boost::archive::no_header);
            archive>>core;
            require(input.peek()==std::char_traits<char>::eof(),"ACTUAL_CORE_TRAILING_BYTES");
        }
        require(core.x.size()==dimension && core.P.rows()==dimension && core.P.cols()==dimension &&
            core.dx.size()==dimension && core.kfIndexMap.size()==std::size_t(dimension),"ACTUAL_CORE_DIMENSION_DIFFERS");
        const auto runtimeEntry=core.metaDataMap.find(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);
        require(runtimeEntry!=core.metaDataMap.end() && runtimeEntry->second==expectedRuntime &&
            !expectedRuntime.empty(),"ACTUAL_CORE_RUNTIME_ID_DIFFERS");
        require(!core.metaDataMap.contains(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA),"ACTUAL_CORE_IS_DISPOSABLE_BRANCH");
        std::map<std::string,Receiver> receivers;
        for(const auto& id:receiverIds(core))receivers.try_emplace(id).first->second.id=id;
        KFState owner;std::string reason;
        const auto resolver=[&](const std::string& id)->Receiver*
        {const auto found=receivers.find(id);return found==receivers.end()?nullptr:&found->second;};
        report["stage"]="ACTUAL_CORE_RESTORE";
        require(restoreZhangCheckpointKfCoreWithReceiverResolver(core,owner,resolver,&reason),
            "ACTUAL_CORE_NATIVE_RESTORE_REJECTED:"+reason);
        require(zhangCheckpointKfCoreMatchesState(core,owner),"ACTUAL_CORE_NOT_BITWISE_EQUAL_AFTER_RESTORE");
        require(zhangCheckpointRuntimeId(owner)==expectedRuntime,"ACTUAL_OWNER_RUNTIME_REBOUND_INCORRECTLY");
        report["actual_core_dimension"]=dimension;report["stable_receiver_ids_from_core_only"]=receivers.size();
        report["actual_x_dx_full_P_and_all_persisted_core_equal"]=true;
        mode();
        require(zhangGecArRuntimeRegistriesEmpty(),"REAL_TEST_PROCESS_ALREADY_HAS_GEC_AUTHORITY");
        std::string source(realGecBytes,'\0');
        {
            std::ifstream input(gecPath,std::ios::binary);
            require(bool(input),"ACTUAL_GEC_OPEN_FAILED");
            input.read(source.data(),static_cast<std::streamsize>(source.size()));
            require(input.gcount()==static_cast<std::streamsize>(source.size()) &&
                input.peek()==std::char_traits<char>::eof(),"ACTUAL_GEC_TRUNCATED_OR_TRAILING_BYTES");
        }
        report["stage"]="ACTUAL_GEC_NATIVE_PREFLIGHT";
        require(importZhangGecArRuntimeSection(owner,expectedRuntime,source,true,dimension,reason),
            "REAL_GEC_PREFLIGHT_REJECTED:"+reason);
        require(zhangGecArRuntimeRegistriesEmpty(),"REAL_GEC_PREFLIGHT_MUTATED_REGISTRIES");
        require(zhangCheckpointKfCoreMatchesState(core,owner),"REAL_GEC_PREFLIGHT_CHANGED_ACTUAL_CORE");
        for(int round=1;round<=2;++round)
        {
            report["stage"]="ACTUAL_GEC_NATIVE_IMPORT";report["round"]=round;
            require(importZhangGecArRuntimeSection(owner,expectedRuntime,source,false,dimension,reason),
                "REAL_GEC_IMPORT_REJECTED:"+reason);
            require(zhangCheckpointKfCoreMatchesState(core,owner),"REAL_GEC_IMPORT_CHANGED_ACTUAL_CORE");
            report["stage"]="ACTUAL_GEC_NATIVE_EXPORT";
            std::string after;
            require(exportZhangGecArRuntimeSection(owner,expectedRuntime,after,reason),"REAL_GEC_EXPORT_REJECTED:"+reason);
            require(zhangCheckpointKfCoreMatchesState(core,owner),"REAL_GEC_EXPORT_CHANGED_ACTUAL_CORE");
            report["actual_core_unchanged_after_import_and_export"]=true;
            if(round==1)
            {
                std::ofstream output(outputPath,std::ios::binary|std::ios::out);
                require(bool(output),"FRESH_ACTUAL_GEC_OUTPUT_OPEN_FAILED");
                output.write(after.data(),static_cast<std::streamsize>(after.size()));output.close();
                require(bool(output),"ACTUAL_GEC_OUTPUT_WRITE_FAILED");
                report["raw_reexport_written"]=outputPath;report["reexport_bytes"]=after.size();
                report["reexport_sha256"]=zhangCheckpointSha256(after);
            }
            if(after!=source)
            {
                const auto commonBytes=std::min(after.size(),source.size());std::size_t offset=0;
                while(offset<commonBytes && source[offset]==after[offset])++offset;
                report["first_different_byte_offset"]=offset;report["stage"]="COMPLETE_WIRE_COMPARE";
                throw std::runtime_error("GEC_AR_REAL_COMPONENT_WIRE_CHANGED_NO_TOLERANCE");
            }
        }
        require(zhangCheckpointKfCoreMatchesState(core,owner),"REAL_GEC_FINAL_ACTUAL_CORE_DIFFERS");
        fileHash(corePath,coreSha);fileHash(gecPath,gecSha);
        report["stage"]="COMPONENT_PASS";report["complete_raw_GEC_bytes_equal"]=true;
        report["source_core_and_GEC_hashes_unchanged"]=true;report["actual_production_roundtrips"]=2;
        report["no_synthetic_core_or_zero_covariance"]=true;report["passed"]=true;
        std::cout<<"REAL_GEC_ROUNDTRIP_RESULT "<<j::serialize(report)<<'\n';return 0;
    }
    catch(const std::exception& error)
    {
        report["error"]=error.what();
        std::cout<<"REAL_GEC_ROUNDTRIP_RESULT "<<j::serialize(report)<<'\n';return 1;
    }
}
