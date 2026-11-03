#pragma once

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>

// Same SHA-256 primitive as the strict checkpoint reader. Forward-declared to
// keep this policy header independent of KF/Boost archive implementation.
std::string zhangCheckpointSha256(const std::string& bytes);

// This is one versioned, audited experiment transition, not an option to
// ignore checkpoint provenance.  In particular, the old archive is never
// relabelled or rewritten to claim that the new executable produced it.
inline constexpr const char* ZHANG_GEC134_MIGRATION_POLICY = "GEC134_FLOAT_E121_TO_AR_V1";
inline constexpr const char* ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG =
    "ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6";
inline constexpr const char* ZHANG_GEC134_LOCAL_ONE_EPOCH_POLICY =
    "GEC134_FLOAT_E121_TO_AR_LOCAL_12_6_ONE_EPOCH_V1";
inline constexpr const char* ZHANG_GEC134_LOCAL_ONE_EPOCH_ROOT =
    "/home/rx/GINAN/multi_gnss_float_profile_20261022/"
    "server134_fixedcoord_local_memory_20261024_01";
inline constexpr const char* ZHANG_GEC134_SOURCE_MODE = "E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS";
inline constexpr const char* ZHANG_GEC134_AR_MODE = "E29_MULTI_GNSS_AR_V1_GPS_GAL_BDS";
inline constexpr const char* ZHANG_GEC134_SOURCE_BUNDLE_SHA256 =
    "1310d0f6543e552ac318208b79491139d93886de87d8331ee14cb14dc3bf283b";
inline constexpr const char* ZHANG_GEC134_SOURCE_BINARY_SHA256 =
    "86f69f9aca5c6e5ba238ba193e1127bb57b8ba5b092f99bb8a13ff102316b288";
inline constexpr const char* ZHANG_GEC134_SOURCE_CONFIG_SHA256 =
    "6ec5bbb76c9c0cc977b7c2bcc5197a44b3ab70ff239cb6a418ac5fd9a085b21a";
inline constexpr const char* ZHANG_GEC134_SOURCE_INPUT_SHA256 =
    "e325b47f85b6061b0859bf9af7e9256c7f6bd88a1cf2830099cd4faca06f5da4";
inline constexpr const char* ZHANG_GEC134_SOURCE_RUNTIME_ID =
    "E29-GEC134-FIXEDCOORD-20261023-one_hour";
inline constexpr const char* ZHANG_GEC134_SOURCE_CHECKPOINT_ID =
    "E29-20240717010000-e121-r32f626de9463";
inline constexpr const char* ZHANG_GEC134_SOURCE_BINARY_PATH =
    "/home/rx/GINAN/multi_gnss_float_profile_20261022/"
    "server180_matrix14_parallel_validation01_staging/producer_matrix_parallel_v1/pea";
inline constexpr const char* ZHANG_GEC134_SOURCE_CHECKPOINT_PATH =
    "/home/rx/GINAN/multi_gnss_float_profile_20261022/"
    "server134_fixedcoord_float1h_offline_v3_20261023_01/checkpoints/"
    "E29-20240717010000-e121-r32f626de9463";
inline constexpr const char* ZHANG_GEC134_AR_START = "1405213230";

class ZhangGec134MigrationAuthority;
inline bool zhangGec134LocalOneEpochRequested(
    const std::map<std::string, std::string>& environment)
{
    const auto local = environment.find(ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG);
    return local != environment.end() && local->second == "1";
}

inline bool validateZhangGec134Environment(
    const std::map<std::string, std::string>& environment, std::string& reason)
{
    std::map<std::string, std::string> required = {
        {"ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", "1"},
        {"ZHANG_MULTI_GNSS_AR_CHECKPOINT", "1"},
        {"ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "UNSET"},
        {"ZHANG_R51_REUSE_FLOAT_20240717", "UNSET"},
        {"ZHANG_R51_REUSE_RATIO_20240717", "UNSET"},
        {"ZHANG_R51_ENABLE", "1"}, {"ZHANG_R51_RATIO_ONLY", "1"},
        {"ZHANG_R51_RATIO_THRESHOLD", "2"}, {"ZHANG_R51_EPOCH_ONLY_AR", "1"},
        {"ZHANG_R51_AR_START_GPST_SECONDS", ZHANG_GEC134_AR_START},
        {"ZHANG_R51_FLOAT_THREADS", "28"}, {"OMP_NUM_THREADS", "28"},
        {"OPENBLAS_NUM_THREADS", "1"}, {"MKL_NUM_THREADS", "1"},
        {"ZHANG_FLOAT_MATRIX_THREADS", "14"}, {"ZHANG_FLOAT_MEMORY_MODE", "parallel"},
        {"ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC", "off"},
        {"ZHANG_FLOAT_PROFILE", "1"}, {"ZHANG_P0_RESOURCE_PROBE", "1"},
        {"ZHANG_R49_FUSION_WEIGHT", "0.10"}, {"ZHANG_R51_SINGLE_BASIS_HANDOFF", "1"},
        {"ZHANG_R51_PRODUCT_REPAIR", "1"}, {"ZHANG_R51_INTERN_HANDOFF_PROOFS", "1"}};
    // This is not an available-CPU fallback for the server experiment. Only
    // the explicitly named local profile changes these three environment
    // values, and the canonical validator below also requires actual 12/6
    // CPU evidence, the one exact output root and a single terminal epoch.
    if (zhangGec134LocalOneEpochRequested(environment))
    {
        required.emplace(ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG, "1");
        required.at("ZHANG_R51_FLOAT_THREADS") = "12";
        required.at("OMP_NUM_THREADS") = "12";
        required.at("ZHANG_FLOAT_MATRIX_THREADS") = "6";
    }
    for (const auto& [key, expected] : required)
    {
        const auto found = environment.find(key);
        if (found == environment.end() || found->second != expected)
        { reason = "GEC134_MIGRATION_ENVIRONMENT_REJECTED:" + key; return false; }
    }
    for (const auto& [key, value] : environment)
    {
        if (required.count(key)) continue;
        if (key == "ZHANG_E29_RESTORE_BUNDLE" && !value.empty() && value != "UNSET") continue;
        if (key == "ZHANG_GEC134_RESTORE_ONLY" && (value == "UNSET" || value == "1")) continue;
        reason = "GEC134_MIGRATION_UNLISTED_ENVIRONMENT:" + key;
        return false;
    }
    reason.clear();
    return true;
}
inline bool validateZhangGec134Migration(
    const std::string&, const std::string&, const std::string&,
    const std::string&, const std::string&, const std::string&,
    const std::string&, const std::map<std::string, std::string>&,
    ZhangGec134MigrationAuthority&, std::string&);

/** Only the complete whitelist validator can create a valid authority. The
 * token binds the immutable source to the *exact* target canonical text; it
 * is not a general "different binary is OK" flag. */
class ZhangGec134MigrationAuthority
{
public:
    bool valid() const { return valid_; }
    const std::string& targetConfigText() const { return targetConfigText_; }
    bool matchesSource(const std::string& mode, const std::string& binarySha,
        const std::string& configSha, const std::string& inputSha) const
    {
        return valid_ && mode == ZHANG_GEC134_SOURCE_MODE &&
            binarySha == ZHANG_GEC134_SOURCE_BINARY_SHA256 &&
            configSha == ZHANG_GEC134_SOURCE_CONFIG_SHA256 &&
            inputSha == ZHANG_GEC134_SOURCE_INPUT_SHA256;
    }
private:
    bool valid_ = false;
    std::string targetConfigText_;
    friend bool validateZhangGec134Migration(
        const std::string&, const std::string&, const std::string&,
        const std::string&, const std::string&, const std::string&,
        const std::string&, const std::map<std::string, std::string>&,
        ZhangGec134MigrationAuthority&, std::string&);
};

namespace ZhangGec134MigrationDetail
{
inline bool unsignedValue(const std::string& text, std::size_t& value)
{
    if (text.empty() || text.size() > 18 ||
        !std::all_of(text.begin(), text.end(), [](unsigned char c){ return std::isdigit(c); }))
        return false;
    value = 0;
    for (const char c : text)
    {
        if (value > (std::numeric_limits<std::size_t>::max() - (c - '0')) / 10)
            return false;
        value = value * 10 + (c - '0');
    }
    return true;
}

inline bool canonicalFields(const std::string& text,
    std::map<std::string, std::string>& fields, std::string& reason)
{
    fields.clear();
    const std::string header = "E29_CANONICAL_CONFIG_V1\n";
    if (text.size() > 16 * 1024 * 1024 || text.compare(0, header.size(), header) != 0)
    { reason = "GEC134_MIGRATION_CANONICAL_HEADER"; return false; }
    std::size_t cursor = header.size();
    while (cursor < text.size())
    {
        const auto newline = text.find('\n', cursor);
        const auto equal = text.find('=', cursor);
        if (newline == std::string::npos || equal == std::string::npos || equal >= newline)
        { reason = "GEC134_MIGRATION_CANONICAL_RECORD"; return false; }
        const auto key = text.substr(cursor, equal - cursor);
        const auto sizeText = text.substr(equal + 1, newline - equal - 1);
        std::size_t size = 0;
        if (!unsignedValue(sizeText, size) || fields.size() >= 4096)
        { reason = "GEC134_MIGRATION_CANONICAL_LENGTH"; return false; }
        cursor = newline + 1;
        const bool scalar = key == "command_argument_count" || key == "config_file_count" ||
            key == "effective_capture_epoch_count" ||
            (key.rfind("config_file_", 0) == 0 && key.size() >= 5 &&
             key.compare(key.size() - 5, 5, "_size") == 0);
        std::string value = sizeText;
        if (!scalar)
        {
            if (size > text.size() - cursor || cursor + size >= text.size() ||
                text[cursor + size] != '\n')
            { reason = "GEC134_MIGRATION_CANONICAL_TRUNCATED"; return false; }
            value = text.substr(cursor, size);
            cursor += size + 1;
        }
        if (key.empty() || !fields.emplace(key, std::move(value)).second)
        { reason = "GEC134_MIGRATION_CANONICAL_DUPLICATE"; return false; }
    }
    return true;
}

inline bool epochOffset(const std::string& time, int& seconds)
{
    if (time.size() != 19 || time.compare(0, 11, "2024-07-17 ") != 0 ||
        time[13] != ':' || time[16] != ':') return false;
    for (const auto index : {11, 12, 14, 15, 17, 18})
        if (!std::isdigit(static_cast<unsigned char>(time[index]))) return false;
    const int hour = (time[11] - '0') * 10 + time[12] - '0';
    const int minute = (time[14] - '0') * 10 + time[15] - '0';
    const int second = (time[17] - '0') * 10 + time[18] - '0';
    if (hour > 2 || minute >= 60 || (second != 0 && second != 30)) return false;
    seconds = hour * 3600 + minute * 60 + second;
    return true;
}

inline bool flattenOverlay(const YAML::Node& node, const std::string& path,
    std::map<std::string, YAML::Node>& leaves, std::string& reason, std::size_t depth = 0)
{
    if (depth > 8 || leaves.size() > 16)
    { reason = "GEC134_MIGRATION_OVERLAY_DEPTH_OR_SIZE"; return false; }
    if (node.IsMap())
    {
        std::set<std::string> keys;
        if (node.size() == 0) { reason = "GEC134_MIGRATION_EMPTY_OVERLAY_MAP"; return false; }
        for (const auto& entry : node)
        {
            if (!entry.first.IsScalar()) { reason = "GEC134_MIGRATION_OVERLAY_KEY"; return false; }
            const auto key = entry.first.as<std::string>();
            if (!keys.insert(key).second || key.empty() || key.find('.') != std::string::npos)
            { reason = "GEC134_MIGRATION_OVERLAY_DUPLICATE_KEY"; return false; }
            if (!flattenOverlay(entry.second, path.empty() ? key : path + "." + key,
                                leaves, reason, depth + 1)) return false;
        }
        return true;
    }
    if ((!node.IsScalar() && !node.IsSequence()) || !leaves.emplace(path, node).second)
    { reason = "GEC134_MIGRATION_OVERLAY_NODE"; return false; }
    return true;
}
} // namespace ZhangGec134MigrationDetail

/** The source hashes must be computed over actual bytes by the caller. All
 * first ten included files (boxwing + nine original CLI layers) remain exact;
 * the sole extra layer cannot contain models, station coordinates or threads.
 * Unknown environment policy and unknown canonical records fail closed. */
inline bool validateZhangGec134Migration(
    const std::string& sourceBundleSha, const std::string& sourceBinarySha,
    const std::string& sourceConfigSha, const std::string& sourceInputSha,
    const std::string& sourceConfigText, const std::string& targetConfigText,
    const std::string& targetInputSha, const std::map<std::string, std::string>& environment,
    ZhangGec134MigrationAuthority& authority, std::string& reason)
{
    using namespace ZhangGec134MigrationDetail;
    authority = {};
    auto fail = [&](const std::string& code){ reason = code; return false; };
    if (sourceBundleSha != ZHANG_GEC134_SOURCE_BUNDLE_SHA256 ||
        sourceBinarySha != ZHANG_GEC134_SOURCE_BINARY_SHA256 ||
        sourceConfigSha != ZHANG_GEC134_SOURCE_CONFIG_SHA256 ||
        sourceInputSha != ZHANG_GEC134_SOURCE_INPUT_SHA256 || targetInputSha != sourceInputSha)
        return fail("GEC134_MIGRATION_PINNED_SOURCE_IDENTITY_REJECTED");
    if (zhangCheckpointSha256(sourceConfigText) != sourceConfigSha)
        return fail("GEC134_MIGRATION_SOURCE_CANONICAL_HASH_REJECTED");
    if (!validateZhangGec134Environment(environment, reason)) return false;
    const bool localOneEpoch = zhangGec134LocalOneEpochRequested(environment);
    std::map<std::string, std::string> source, target;
    if (!canonicalFields(sourceConfigText, source, reason) ||
        !canonicalFields(targetConfigText, target, reason)) return false;
    auto take = [](std::map<std::string, std::string>& values, const std::string& key,
                   std::string& value)
    {
        const auto found = values.find(key);
        if (found == values.end()) return false;
        value = found->second; values.erase(found); return true;
    };
    std::string value;
    if (!take(source, "environment_ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", value) || value != "1" ||
        !take(target, "environment_ZHANG_MULTI_GNSS_AR_CHECKPOINT", value) || value != "1" ||
        !take(target, "environment_ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", value) || value != "1")
        return fail("GEC134_MIGRATION_EXPLICIT_MODE_REJECTED");
    if (localOneEpoch)
    {
        if (!take(target, std::string("environment_") + ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG, value) ||
            value != "1" ||
            !take(target, "gec134_local_migration_policy", value) ||
            value != ZHANG_GEC134_LOCAL_ONE_EPOCH_POLICY)
            return fail("GEC134_LOCAL_MIGRATION_EXPLICIT_CANONICAL_POLICY_REJECTED");
        // Values on the left are taken from the SHA-locked REAL source
        // canonical archive: comput exposed 28 logical AND 28 physical CPUs.
        // Never overwrite/normalise target canonical bytes to pretend the
        // local Ryzen machine has those CPUs. Consume only these four exact
        // source/target pairs; every other original canonical field must still
        // compare byte-for-byte in the final map comparison below.
        const std::map<std::string, std::pair<std::string, std::string>> localCpu = {
            {"environment_OMP_NUM_THREADS", {"28", "12"}},
            {"environment_ZHANG_R51_FLOAT_THREADS", {"28", "12"}},
            {"float_available_cpus", {"28", "12"}},
            {"float_physical_cores", {"28", "6"}}};
        for (const auto& [key, expected] : localCpu)
        {
            std::string original, actual;
            if (!take(source, key, original) || original != expected.first ||
                !take(target, key, actual) || actual != expected.second)
                return fail("GEC134_LOCAL_MIGRATION_CPU_POLICY_REJECTED:" + key);
        }
    }
    for (const auto& key : {"ZHANG_FLOAT_MATRIX_THREADS", "ZHANG_FLOAT_MEMORY_MODE",
                            "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC", "ZHANG_FLOAT_PROFILE"})
        if (!take(target, std::string("environment_") + key, value) ||
            value != environment.at(key))
            return fail(std::string("GEC134_MIGRATION_NUMERICAL_POLICY_REJECTED:") + key);
    if (!take(source, "environment_ZHANG_R51_AR_START_GPST_SECONDS", value) || value != "1405295970" ||
        !take(target, "environment_ZHANG_R51_AR_START_GPST_SECONDS", value) || value != ZHANG_GEC134_AR_START)
        return fail("GEC134_MIGRATION_AR_SCHEDULE_REJECTED");
    std::string sourceCountText, targetCountText;
    std::size_t sourceCount = 0, targetCount = 0;
    if (!take(source, "config_file_count", sourceCountText) ||
        !take(target, "config_file_count", targetCountText) ||
        !unsignedValue(sourceCountText, sourceCount) ||
        !unsignedValue(targetCountText, targetCount) || sourceCount != 10 || targetCount != 11)
        return fail("GEC134_MIGRATION_CONFIG_LAYER_COUNT_REJECTED");
    const auto prefix = "config_file_" + std::to_string(sourceCount) + "_";
    std::string overlayPath, overlayBytes, overlaySize, overlayHash;
    if (!take(target, prefix + "path", overlayPath) ||
        !take(target, prefix + "bytes", overlayBytes) ||
        !take(target, prefix + "size", overlaySize) ||
        !take(target, prefix + "sha256", overlayHash) || overlayHash.size() != 64 ||
        !unsignedValue(overlaySize, targetCount) || targetCount != overlayBytes.size())
        return fail("GEC134_MIGRATION_EXTRA_OVERLAY_RECORD_REJECTED");
    if (zhangCheckpointSha256(overlayBytes) != overlayHash)
        return fail("GEC134_MIGRATION_EXTRA_OVERLAY_HASH_REJECTED");
    if (!take(source, "command_argument_count", sourceCountText) ||
        !take(target, "command_argument_count", targetCountText) ||
        !unsignedValue(sourceCountText, sourceCount) ||
        !unsignedValue(targetCountText, targetCount) || sourceCount != 149 || targetCount != 150)
        return fail("GEC134_MIGRATION_ARGUMENT_COUNT_REJECTED");
    for (std::size_t i = 1; i <= sourceCount; i++)
    {
        std::string oldArgument, newArgument;
        if (!take(source, "command_argument_" + std::to_string(i), oldArgument) ||
            !take(target, "command_argument_" + std::to_string(i <= 10 ? i : i + 1), newArgument) ||
            oldArgument != newArgument)
            return fail("GEC134_MIGRATION_ARGUMENT_CHANGED:" + std::to_string(i));
    }
    if (!take(target, "command_argument_11", value) || value != overlayPath)
        return fail("GEC134_MIGRATION_OVERLAY_ARGUMENT_REJECTED");

    std::map<std::string, YAML::Node> leaves;
    try
    {
        const auto documents = YAML::LoadAll(overlayBytes);
        if (documents.size() != 1 || !flattenOverlay(documents.front(), "", leaves, reason))
            return fail(reason.empty() ? "GEC134_MIGRATION_OVERLAY_DOCUMENTS" : reason);
    }
    catch (const std::exception&) { return fail("GEC134_MIGRATION_OVERLAY_PARSE_FAILED"); }
    const std::string epoch = "processing_options.epoch_control.";
    const std::string ar = "processing_options.gnss_general.zhang_pppar.";
    const std::set<std::string> allowed = {
        epoch + "end_epoch", epoch + "max_epochs", ar + "checkpoint_capture_epochs",
        ar + "checkpoint_output_directory", ar + "checkpoint_restore_path",
        ar + "product_filename", ar + "product_covariance_filename",
        "outputs.outputs_root", "outputs.metadata.config_description"};
    if (leaves.size() != allowed.size()) return fail("GEC134_MIGRATION_OVERLAY_LEAF_COUNT");
    for (const auto& [key, node] : leaves)
        if (!allowed.count(key) || (key != ar + "checkpoint_capture_epochs" && !node.IsScalar()))
            return fail("GEC134_MIGRATION_OVERLAY_FORBIDDEN_KEY:" + key);
    try
    {
        const auto outputRoot = leaves.at("outputs.outputs_root").as<std::string>();
        const auto runRoot = std::filesystem::path(outputRoot).parent_path();
        const std::filesystem::path base("/home/rx/GINAN/multi_gnss_float_profile_20261022");
        const auto runName = runRoot.filename().string();
        if (runRoot.parent_path() != base || runRoot.filename().empty() ||
            runRoot != runRoot.lexically_normal() ||
            runName.rfind("server134_fixedcoord_", 0) != 0 ||
            !std::all_of(runName.begin(), runName.end(), [](unsigned char c)
                { return std::isalnum(c) || c == '_' || c == '-'; }) ||
            runRoot.filename() == "server134_fixedcoord_float1h_offline_v3_20261023_01" ||
            std::filesystem::path(outputRoot) != runRoot / "outputs" ||
            std::filesystem::path(overlayPath).parent_path() != runRoot / "config" ||
            leaves.at(ar + "checkpoint_output_directory").as<std::string>() != (runRoot / "checkpoints").string() ||
            leaves.at(ar + "checkpoint_restore_path").as<std::string>() != ZHANG_GEC134_SOURCE_CHECKPOINT_PATH ||
            leaves.at(ar + "product_filename").as<std::string>() != (runRoot / "outputs/zhang_internal_products.csv").string() ||
            leaves.at(ar + "product_covariance_filename").as<std::string>() != (runRoot / "outputs/zhang_internal_product_covariance.csv").string())
            return fail("GEC134_MIGRATION_ISOLATED_OUTPUT_PATHS_REJECTED");
        if (localOneEpoch && runRoot != std::filesystem::path(ZHANG_GEC134_LOCAL_ONE_EPOCH_ROOT))
            return fail("GEC134_LOCAL_MIGRATION_EXACT_OUTPUT_ROOT_REJECTED");
        // Existing parents may not redirect a nominally fresh output tree
        // through a symlink into the immutable source run or another tree.
        std::error_code pathError;
        for (const auto& path : {runRoot, std::filesystem::path(outputRoot),
                                runRoot / "checkpoints", runRoot / "restore_audit"})
            if (std::filesystem::weakly_canonical(path, pathError) != path || pathError)
                return fail("GEC134_MIGRATION_OUTPUT_SYMLINK_REJECTED");
        const auto end = leaves.at(epoch + "end_epoch").as<std::string>();
        int endSeconds = 0;
        if (!epochOffset(end, endSeconds) || endSeconds < 3630 || endSeconds > 7200 ||
            leaves.at(epoch + "max_epochs").as<int>() != 121 + (endSeconds - 3600) / 30)
            return fail("GEC134_MIGRATION_BOUNDED_END_REJECTED");
        if (localOneEpoch && (end != "2024-07-17 01:00:30" ||
                leaves.at(epoch + "max_epochs").as<int>() != 122))
            return fail("GEC134_LOCAL_MIGRATION_ONE_EPOCH_REQUIRED");
        const auto captures = leaves.at(ar + "checkpoint_capture_epochs");
        if (!captures.IsSequence() || captures.size() == 0 || captures.size() > 4)
            return fail("GEC134_MIGRATION_CAPTURE_SCHEDULE_REJECTED");
        if (localOneEpoch && (captures.size() != 1 || !captures[0].IsScalar() ||
                captures[0].as<std::string>() != "2024-07-17 01:00:30"))
            return fail("GEC134_LOCAL_MIGRATION_SINGLE_TERMINAL_CAPTURE_REQUIRED");
        std::string previous = "2024-07-17 01:00:00";
        for (const auto& capture : captures)
        {
            if (!capture.IsScalar()) return fail("GEC134_MIGRATION_CAPTURE_VALUE_REJECTED");
            const auto time = capture.as<std::string>();
            int captureSeconds = 0;
            if (!epochOffset(time, captureSeconds) || time <= previous || time > end)
                return fail("GEC134_MIGRATION_CAPTURE_ORDER_REJECTED");
            previous = time;
        }
        if (previous != end) return fail("GEC134_MIGRATION_FINAL_CHECKPOINT_REQUIRED");
        const std::map<std::string, std::string> effective = {
            {"effective_end_epoch", std::string(end).replace(10, 1, "T")},
            {"effective_max_epochs", leaves.at(epoch + "max_epochs").as<std::string>()},
            {"effective_outputs_root", outputRoot},
            {"effective_product_filename", leaves.at(ar + "product_filename").as<std::string>()},
            {"effective_product_covariance_filename", leaves.at(ar + "product_covariance_filename").as<std::string>()},
            {"effective_checkpoint_output_directory", leaves.at(ar + "checkpoint_output_directory").as<std::string>()},
            {"effective_checkpoint_restore_path", ZHANG_GEC134_SOURCE_CHECKPOINT_PATH}};
        for (const auto& [key, expected] : effective)
        {
            if (!take(source, key, value) || !take(target, key, value) || value != expected)
                return fail("GEC134_MIGRATION_EFFECTIVE_VALUE_REJECTED:" + key);
        }
        if (!take(source, "effective_capture_epoch_count", sourceCountText) ||
            !unsignedValue(sourceCountText, sourceCount) || sourceCount != 2 ||
            !take(target, "effective_capture_epoch_count", value) || value != std::to_string(captures.size()))
            return fail("GEC134_MIGRATION_EFFECTIVE_CAPTURE_COUNT_REJECTED");
        for (std::size_t i = 0; i < sourceCount; i++)
            if (!take(source, "effective_capture_epoch_" + std::to_string(i), value))
                return fail("GEC134_MIGRATION_SOURCE_CAPTURE_MISSING");
        for (std::size_t i = 0; i < captures.size(); i++)
            if (!take(target, "effective_capture_epoch_" + std::to_string(i), value) ||
                value != captures[i].as<std::string>())
                return fail("GEC134_MIGRATION_EFFECTIVE_CAPTURE_REJECTED");
    }
    catch (const std::exception&) { return fail("GEC134_MIGRATION_OVERLAY_VALUE_REJECTED"); }
    // This catches changes to *any* original input, physical option, included
    // file bytes/hash/path, CPU/thread setting, runtime owner or station XYZ.
    if (source != target) return fail("GEC134_MIGRATION_NONWHITELISTED_CONFIG_CHANGE");
    authority.valid_ = true;
    authority.targetConfigText_ = targetConfigText;
    reason.clear();
    return true;
}
