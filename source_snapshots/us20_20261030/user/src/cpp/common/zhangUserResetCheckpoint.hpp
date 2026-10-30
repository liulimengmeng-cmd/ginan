#pragma once

#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"
#include <boost/json.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#if !defined(_WIN32)
extern char** environ;
#endif

// This is a separate, explicitly frozen experiment fork mode. It never
// normalizes a legacy E29 configuration or accepts an arbitrary new policy.
inline bool zhangUs20CheckpointEnabled()
{
    const char* value=std::getenv("ZHANG_US20_RESET_CHECKPOINT");
    return value && std::string(value)=="1";
}
inline std::string zhangUs20Environment(const char* name)
{
    const char* value=std::getenv(name);return value?value:"";
}
inline std::string zhangUs20ReadFile(const std::filesystem::path& path)
{
    if(!std::filesystem::is_regular_file(path))throw std::runtime_error("US20_FORK_FILE_MISSING:"+path.string());
    std::ifstream input(path,std::ios::binary);
    std::string text{std::istreambuf_iterator<char>(input),{}};
    if(input.bad())throw std::runtime_error("US20_FORK_FILE_READ:"+path.string());
    return text;
}
inline std::string zhangUs20JsonString(const boost::json::object& value,const char* key)
{
    return std::string(value.at(key).as_string());
}
inline std::pair<std::string,boost::json::object> zhangUs20ForkPlan()
{
    const auto path=zhangUs20Environment("ZHANG_US20_FORK_PLAN");
    const auto text=zhangUs20ReadFile(path);
    const auto digest=zhangCheckpointSha256(text);
    if(digest!=zhangUs20Environment("ZHANG_US20_FORK_PLAN_SHA256"))throw std::runtime_error("US20_FORK_PLAN_HASH_MISMATCH");
    auto plan=boost::json::parse(text).as_object();
    if(zhangUs20JsonString(plan,"schema")!="US20_GEC_RESET_FORK_PLAN_V1")throw std::runtime_error("US20_FORK_PLAN_SCHEMA");
    return {text,std::move(plan)};
}
inline boost::json::object zhangUs20ForkBranch(const boost::json::object& plan)
{
    const auto id=zhangUs20Environment("ZHANG_US20_FORK_BRANCH");
    if(id.empty())throw std::runtime_error("US20_FORK_BRANCH_MISSING");
    return plan.at("branches").as_object().at(id).as_object();
}

inline void zhangUs20AppendExactEnvironment(std::string& config)
{
    std::map<std::string,std::string> values;
#if !defined(_WIN32)
    for(char** current=environ;current && *current;++current)
    {
        const std::string field=*current;const auto split=field.find('=');
        if(split==std::string::npos)continue;
        const auto key=field.substr(0,split);
        if(!(key.starts_with("ZHANG_") || key=="OPENBLAS_NUM_THREADS" || key=="OMP_NUM_THREADS" || key=="MKL_NUM_THREADS"))continue;
        if(key=="ZHANG_US20_FORK_PLAN" || key=="ZHANG_US20_FORK_PLAN_SHA256" ||
           key=="ZHANG_US20_FORK_BRANCH" || key=="ZHANG_US20_DUMP_PROVENANCE" ||
           key=="ZHANG_US20_RESTORE_RECEIPT" || key=="ZHANG_US20_RESTORE_RECEIPT_SHA256" ||
           key=="ZHANG_E29_RESTORE_BUNDLE")continue;
        values[key]=field.substr(split+1);
    }
#else
    throw std::runtime_error("US20_RESET_CHECKPOINT_REQUIRES_TESTED_LINUX_RUNTIME");
#endif
    std::ostringstream out;out<<"US20_EXACT_ENVIRONMENT_V1\n";
    for(const auto& [key,value]:values)out<<std::quoted(key)<<' '<<std::quoted(value)<<'\n';
    config+=out.str();
}

inline void zhangUs20AppendInput(std::string& inputs,const std::string& role,const std::string& filename)
{
    if(filename.empty())throw std::runtime_error("US20_CHECKPOINT_INPUT_PATH_MISSING:"+role);
    const auto path=std::filesystem::canonical(filename);
    const auto digest=zhangCheckpointFileSha256(path.string());
    if(digest.empty())throw std::runtime_error("US20_CHECKPOINT_INPUT_HASH_FAILED:"+role);
    std::ostringstream out;out<<std::quoted(role)<<' '<<std::quoted(path.string())<<' '
        <<std::filesystem::file_size(path)<<' '<<digest<<'\n';inputs+=out.str();
}

inline bool zhangUs20AuthoriseForkProvenance(std::string& config,std::string& inputs,const std::string& binary)
{
    // A dump is an input-inspection operation. main exits before any epoch.
    const auto dump=zhangUs20Environment("ZHANG_US20_DUMP_PROVENANCE");
    if(!dump.empty())
    {
        const auto root=std::filesystem::path(dump);
        if(std::filesystem::exists(root))throw std::runtime_error("US20_PROVENANCE_DUMP_ALREADY_EXISTS");
        std::filesystem::create_directories(root);
        std::ofstream(root/"canonical_config.txt",std::ios::binary)<<config;
        std::ofstream(root/"input_manifest.txt",std::ios::binary)<<inputs;
        boost::json::object receipt{{"binary_sha256",zhangCheckpointFileSha256(binary)},
            {"canonical_config_sha256",zhangCheckpointSha256(config)},
            {"input_manifest_sha256",zhangCheckpointSha256(inputs)},
            {"policy",zhangUs20Environment("ZHANG_USER_DATUM_POLICY")}};
        std::ofstream(root/"receipt.json")<<boost::json::serialize(receipt)<<'\n';
        return false;
    }
    const auto [text,plan]=zhangUs20ForkPlan();const auto branch=zhangUs20ForkBranch(plan);
    if(zhangUs20JsonString(plan,"binary_sha256")!=zhangCheckpointFileSha256(binary) ||
       zhangUs20JsonString(branch,"canonical_config_sha256")!=zhangCheckpointSha256(config) ||
       zhangUs20JsonString(branch,"input_manifest_sha256")!=zhangCheckpointSha256(inputs) ||
       zhangUs20JsonString(branch,"policy")!=zhangUs20Environment("ZHANG_USER_DATUM_POLICY"))
        throw std::runtime_error("US20_UNAPPROVED_FORK_CONFIGURATION_OR_INPUT");
    // Every permitted branch has its exact bytes enumerated in the immutable
    // common plan. Its construction separately verifies the narrow scientific
    // delta (policy or certified future gauge fixture), never a fuzzy path match.
    config="US20_GEC_USER_CHECKPOINT_V1_APPROVED_FORKS\n"+text;
    inputs="US20_GEC_USER_CHECKPOINT_V1_APPROVED_INPUTS\n"+text;
    return true;
}

inline void zhangUs20ValidateLedgerItem(const std::string& name,const std::string& text,long double completed)
{
    std::istringstream input(text);
    std::string product,receiver,epoch,extra;std::getline(input,product);std::getline(input,receiver);std::getline(input,epoch);
    std::size_t used=0;const auto seconds=std::stoll(epoch,&used);
    if(product.empty() || receiver.empty() || used!=epoch.size() || seconds>completed+1e-6L ||
       name!=zhangCheckpointSha256(receiver+":"+epoch) || std::getline(input,extra))
        throw std::runtime_error("US20_TROP_LEDGER_IDENTITY_OR_TIME_MISMATCH");
}
inline boost::json::object zhangUs20LedgerEntries(const std::filesystem::path& root,long double completed)
{
    boost::json::object entries;
    if(!std::filesystem::exists(root))return entries;
    for(const auto& entry:std::filesystem::directory_iterator(root))
    {
        if(entry.is_symlink() || !entry.is_directory())throw std::runtime_error("US20_TROP_LEDGER_UNEXPECTED_ENTRY");
        const auto name=entry.path().filename().string();const auto source=entry.path()/"source.txt";
        if(std::filesystem::is_symlink(source))throw std::runtime_error("US20_TROP_LEDGER_SOURCE_SYMLINK");
        std::size_t files=0;for(const auto& ignored:std::filesystem::directory_iterator(entry.path()))++files;
        if(files!=1)throw std::runtime_error("US20_TROP_LEDGER_DIRECTORY_NOT_EXACT");
        const auto text=zhangUs20ReadFile(source);zhangUs20ValidateLedgerItem(name,text,completed);
        entries[name]=text;
    }
    return entries;
}

inline void zhangUs20CaptureLedger(KFState& state)
{
    if(!zhangUs20CheckpointEnabled())return;
    const auto path=zhangUs20Environment("ZHANG_USER_TROP_CONSUMPTION_LEDGER");
    if(path.empty())throw std::runtime_error("US20_TROP_LEDGER_PATH_REQUIRED");
    const auto [text,plan]=zhangUs20ForkPlan();
    state.metaDataMap["US20_RESET_CHECKPOINT_LEDGER_V1"]=boost::json::serialize(zhangUs20LedgerEntries(path,state.time.bigTime));
    state.metaDataMap["US20_RESET_CHECKPOINT_PLAN_SHA256"]=zhangCheckpointSha256(text);
    state.metaDataMap["US20_RESET_CHECKPOINT_SOURCE_BRANCH"]=zhangUs20Environment("ZHANG_US20_FORK_BRANCH");
}

inline void zhangUs20RestoreLedger(const KFState& state,bool dryRun=false)
{
    if(!zhangUs20CheckpointEnabled())return;
    const auto [text,plan]=zhangUs20ForkPlan();const auto branch=zhangUs20ForkBranch(plan);
    if(state.metaDataMap.at("US20_RESET_CHECKPOINT_PLAN_SHA256")!=zhangCheckpointSha256(text) ||
       state.metaDataMap.at("US20_RESET_CHECKPOINT_SOURCE_BRANCH")!=zhangUs20JsonString(branch,"restore_from") ||
       state.time.to_string(0)!=zhangUs20JsonString(branch,"restore_epoch"))
        throw std::runtime_error("US20_FORK_CHECKPOINT_SOURCE_OR_EPOCH_REJECTED");
    const auto ledger=boost::json::parse(state.metaDataMap.at("US20_RESET_CHECKPOINT_LEDGER_V1")).as_object();
    for(const auto& entry:ledger)
        zhangUs20ValidateLedgerItem(std::string(entry.key()),std::string(entry.value().as_string()),state.time.bigTime);
    const auto root=std::filesystem::path(zhangUs20Environment("ZHANG_USER_TROP_CONSUMPTION_LEDGER"));
    if(root.empty())throw std::runtime_error("US20_TROP_LEDGER_PATH_REQUIRED");
    const auto existing=zhangUs20LedgerEntries(root,state.time.bigTime);
    for(const auto& entry:existing)
        if(!ledger.contains(entry.key()) || ledger.at(entry.key())!=entry.value())throw std::runtime_error("US20_FORK_TARGET_LEDGER_CONFLICT");
    if(dryRun)return;
    std::filesystem::create_directories(root);
    for(const auto& entry:ledger)
    {
        if(existing.contains(entry.key()))continue;
        const auto directory=root/std::string(entry.key());
        if(!std::filesystem::create_directory(directory))throw std::runtime_error("US20_LEDGER_RESTORE_CREATE_FAILED");
        std::ofstream output(directory/"source.txt",std::ios::binary);output<<std::string(entry.value().as_string());output.flush();
        if(!output)throw std::runtime_error("US20_LEDGER_RESTORE_WRITE_FAILED");
    }
    if(zhangUs20LedgerEntries(root,state.time.bigTime)!=ledger)throw std::runtime_error("US20_LEDGER_RESTORE_NOT_EXACT");
}

inline std::string zhangUs20CaptureLedgerPayload(const KFState& state)
{
    KFState temporary;temporary.time=state.time;zhangUs20CaptureLedger(temporary);
    boost::json::object payload;
    for(const auto& [key,value]:temporary.metaDataMap)payload[key]=value;
    return boost::json::serialize(payload);
}
inline void zhangUs20RestoreLedgerPayload(const std::string& payload,const GTime& completed,bool dryRun)
{
    KFState temporary;temporary.time=completed;
    const auto parsed=boost::json::parse(payload).as_object();
    if(parsed.size()!=3)throw std::runtime_error("US20_LEDGER_SECTION_FIELDS");
    for(const auto& entry:parsed)temporary.metaDataMap[std::string(entry.key())]=std::string(entry.value().as_string());
    zhangUs20RestoreLedger(temporary,dryRun);
}

inline void zhangUs20ValidateParentCheckpoint(const std::filesystem::path& directory)
{
    if(!zhangUs20CheckpointEnabled())return;
    const auto receiptText=zhangUs20ReadFile(zhangUs20Environment("ZHANG_US20_RESTORE_RECEIPT"));
    if(zhangCheckpointSha256(receiptText)!=zhangUs20Environment("ZHANG_US20_RESTORE_RECEIPT_SHA256"))
        throw std::runtime_error("US20_RESTORE_RECEIPT_HASH_MISMATCH");
    const auto receipt=boost::json::parse(receiptText).as_object();const auto [text,plan]=zhangUs20ForkPlan();
    const auto branch=zhangUs20ForkBranch(plan);
    if(zhangUs20JsonString(receipt,"schema")!="US20_EXACT_PARENT_RECEIPT_V1" ||
       zhangUs20JsonString(receipt,"branch")!=zhangUs20Environment("ZHANG_US20_FORK_BRANCH") ||
       zhangUs20JsonString(receipt,"source_branch")!=zhangUs20JsonString(branch,"restore_from") ||
       zhangUs20JsonString(receipt,"plan_sha256")!=zhangCheckpointSha256(text) ||
       zhangUs20JsonString(receipt,"parent_checkpoint_sha256")!=zhangCheckpointFileSha256((directory/"checkpoint.bundle").string()))
        throw std::runtime_error("US20_RESTORE_PARENT_IDENTITY_MISMATCH");
}
