#pragma once

#include "common/zhangGec134Migration.hpp"

// This capability is deliberately disjoint from FLOAT->AR migration. It may
// change only the output/schedule envelope of one already-AR e122 archive;
// every AR runtime section must still be preflighted and actually imported.
inline constexpr const char* ZHANG_GEC134_AR_CONTINUATION_FLAG =
    "ZHANG_GEC134_AR_CONTINUATION_12_6";
inline constexpr const char* ZHANG_GEC134_AR_CONTINUATION_POLICY =
    "GEC134_AR_E122_LOCAL_12_6_ONE_HOUR_V1";
inline constexpr const char* ZHANG_GEC134_AR_CONTINUATION_ROOT =
    "/home/rx/GINAN/multi_gnss_float_profile_20261022/"
    "server134_fixedcoord_local_ar_e122_hour_memory_20261024_01";
inline constexpr const char* ZHANG_GEC134_AR_CONTINUATION_OVERLAY =
    "/home/rx/GINAN/multi_gnss_float_profile_20261022/"
    "server134_fixedcoord_local_ar_e122_hour_memory_20261024_01/config/"
    "100_gec134_ar_e122_hour_continuation_overlay.yaml";
inline constexpr const char* ZHANG_GEC134_AR_E122_PATH =
    "/home/rx/GINAN/multi_gnss_float_profile_20261022/"
    "server134_fixedcoord_local_memory_20261024_01/checkpoints/"
    "E29-20240717010030-e122-r32f626de9463";
inline constexpr const char* ZHANG_GEC134_AR_E122_ID =
    "E29-20240717010030-e122-r32f626de9463";
inline constexpr const char* ZHANG_GEC134_AR_E122_BINARY_PATH =
    "/home/rx/GINAN/local_real_epoch_memory_20261024_v4/candidate/pea";
inline constexpr const char* ZHANG_GEC134_AR_E122_BINARY_SHA256 =
    "835a818e5b4f880ad7607271546d6fef1a3837a97dfd9f70d5c006fefd48e3cc";
inline constexpr const char* ZHANG_GEC134_AR_E122_BUNDLE_SHA256 =
    "e251a3ca23c36d082c690011b421bcb6f79318cfadc76ef5d41fdb3565bdcdca";
inline constexpr const char* ZHANG_GEC134_AR_E122_PAYLOAD_SHA256 =
    "5fd2125f63381392267b80a8761353f122cc7c78a60381c0cd5737426b0f1c6c";
inline constexpr const char* ZHANG_GEC134_AR_E122_CONFIG_SHA256 =
    "2b803cffa50e5bcc75e196f7172a095ed3b6356992c62433d63255cbd3701258";

// Literal capabilities, not a user-supplied output/schedule override.
// The existing v5 capability is unchanged. The two v6 capabilities have one
// pinned candidate path and disjoint outputs; neither may consume the other
// run's products or substitute its schedule for the approved e122 source.
struct ZhangGec134ArContinuationProfile
{
    std::string root, overlay, endEpoch, endCanonical, candidatePath;
    int maxEpochs;
    long double endBigTime;
    std::vector<std::string> captures;
};

inline const std::vector<ZhangGec134ArContinuationProfile>& zhangGec134ArContinuationProfiles()
{
    static const std::vector<ZhangGec134ArContinuationProfile> profiles = {
        {ZHANG_GEC134_AR_CONTINUATION_ROOT, ZHANG_GEC134_AR_CONTINUATION_OVERLAY,
         "2024-07-17 02:00:30", "2024-07-17T02:00:30", "", 242, 1405216830.0L,
         {"2024-07-17 01:30:30", "2024-07-17 02:00:30"}},
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01/config/100_gec134_ar_e122_v6_pilot_overlay.yaml",
         "2024-07-17 01:07:00", "2024-07-17T01:07:00",
         "/home/rx/GINAN/local_real_epoch_memory_20261025_v6/candidate/pea", 135, 1405213620.0L,
         {"2024-07-17 01:07:00"}},
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_hour_20261025_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_hour_20261025_01/config/100_gec134_ar_e122_v6_hour_overlay.yaml",
         "2024-07-17 02:00:30", "2024-07-17T02:00:30",
         "/home/rx/GINAN/local_real_epoch_memory_20261025_v6/candidate/pea", 242, 1405216830.0L,
         {"2024-07-17 01:30:30", "2024-07-17 02:00:30"}},
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v8_pilot_20261025_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v8_pilot_20261025_01/config/100_gec134_ar_e122_v8_pilot_overlay.yaml",
         "2024-07-17 01:07:00", "2024-07-17T01:07:00",
         "/home/rx/GINAN/local_real_epoch_memory_20261025_v8/candidate/pea", 135, 1405213620.0L,
         {"2024-07-17 01:07:00"}},
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v8_hour_20261025_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v8_hour_20261025_01/config/100_gec134_ar_e122_v8_hour_overlay.yaml",
         "2024-07-17 02:00:30", "2024-07-17T02:00:30",
         "/home/rx/GINAN/local_real_epoch_memory_20261025_v8/candidate/pea", 242, 1405216830.0L,
         {"2024-07-17 01:30:30", "2024-07-17 02:00:30"}},
        // Matched one-epoch A/B capabilities. These change only the literal
        // output envelope and candidate identity; both retain the original
        // e122 payload, numerical environment and full strict restore path.
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_baseline_20261025_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_baseline_20261025_01/config/100_gec134_ar_e122_dual_opt_baseline_20261025_01.yaml",
         "2024-07-17 01:01:00", "2024-07-17T01:01:00",
         "/home/rx/GINAN/local_float_ar_dual_opt_replay_baseline_20261025_v1/candidate/pea", 123, 1405213260.0L,
         {"2024-07-17 01:01:00"}},
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_candidate_20261025_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_candidate_20261025_01/config/100_gec134_ar_e122_dual_opt_candidate_20261025_01.yaml",
         "2024-07-17 01:01:00", "2024-07-17T01:01:00",
         "/home/rx/GINAN/local_float_ar_dual_opt_replay_candidate_20261025_v1/candidate/pea", 123, 1405213260.0L,
         {"2024-07-17 01:01:00"}},
        // Fresh literal one-hour optimization validation. This capability
        // does not extend either paired one-epoch profile or alter old rules.
        {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_hour_20261026_01",
         "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_hour_20261026_01/config/100_gec134_ar_e122_dual_opt_hour_20261026_01.yaml",
         "2024-07-17 02:00:30", "2024-07-17T02:00:30",
         "/home/rx/GINAN/local_float_ar_dual_opt_hour_20261026_v1/candidate/pea", 242, 1405216830.0L,
         {"2024-07-17 01:30:30", "2024-07-17 02:00:30"}}};
    return profiles;
}

inline const ZhangGec134ArContinuationProfile* zhangGec134ArContinuationProfileForOverlay(
    const std::string& overlay)
{
    for (const auto& profile : zhangGec134ArContinuationProfiles())
        if (overlay == profile.overlay) return &profile;
    return nullptr;
}

inline const ZhangGec134ArContinuationProfile* zhangGec134ArContinuationProfileForOutputs(
    const std::string& outputs)
{
    for (const auto& profile : zhangGec134ArContinuationProfiles())
        if (outputs == profile.root + "/outputs") return &profile;
    return nullptr;
}

inline bool zhangGec134ArContinuationRequested(
    const std::map<std::string, std::string>& environment)
{
    const auto found = environment.find(ZHANG_GEC134_AR_CONTINUATION_FLAG);
    return found != environment.end() && found->second == "1";
}

inline bool validateZhangGec134ArContinuationEnvironment(
    const std::map<std::string, std::string>& environment, std::string& reason)
{
    if (!zhangGec134ArContinuationRequested(environment) ||
        !zhangGec134LocalOneEpochRequested(environment))
    { reason = "GEC134_AR_CONTINUATION_EXPLICIT_LOCAL_PROFILE_REQUIRED"; return false; }
    // Retain all source numerical environment values and reject every other
    // new ZHANG/BLAS/OpenMP setting. The old FLOAT validator is unchanged.
    auto sourceCompatible = environment;
    sourceCompatible.erase(ZHANG_GEC134_AR_CONTINUATION_FLAG);
    if (!validateZhangGec134Environment(sourceCompatible, reason)) return false;
    const auto restore = environment.find("ZHANG_E29_RESTORE_BUNDLE");
    bool approvedRestore = restore != environment.end() && restore->second == ZHANG_GEC134_AR_E122_PATH;
    for (const auto& profile : zhangGec134ArContinuationProfiles())
    {
        const std::string checkpoints = profile.root + "/checkpoints/";
        if (restore != environment.end() &&
            ((profile.maxEpochs == 123 && restore->second == checkpoints + "E29-20240717010100-e123-r32f626de9463") ||
             (profile.maxEpochs == 135 && restore->second == checkpoints + "E29-20240717010700-e135-r32f626de9463") ||
             (profile.maxEpochs == 242 &&
              (restore->second == checkpoints + "E29-20240717013030-e182-r32f626de9463" ||
               restore->second == checkpoints + "E29-20240717020030-e242-r32f626de9463"))))
            approvedRestore = true;
    }
    if (!approvedRestore)
    { reason = "GEC134_AR_CONTINUATION_RESTORE_PATH_REJECTED"; return false; }
    reason.clear();
    return true;
}

// Routing only, never restore authority: a checkpoint produced by this exact
// new binary/config/input uses the ordinary strict AR reader/importer. This
// keeps the new profile's own e182/e242 archives restartable without granting
// the e122 cross-provenance token to another source. The reader independently
// checks the complete manifest, ABI and all payload bytes afterwards.
inline bool zhangGec134ArContinuationUsesOrdinaryRestore(
    const std::string& sourceMode, const std::string& sourceBinarySha,
    const std::string& sourceConfigSha, const std::string& sourceInputSha,
    const std::string& sourceConfigText, const std::string& sourceInputText,
    const std::string& targetMode, const std::string& targetBinarySha,
    const std::string& targetConfigText, const std::string& targetInputText)
{
    return sourceMode == ZHANG_GEC134_AR_MODE && targetMode == ZHANG_GEC134_AR_MODE &&
        sourceBinarySha == targetBinarySha && sourceConfigText == targetConfigText &&
        sourceInputText == targetInputText &&
        sourceConfigSha == zhangCheckpointSha256(targetConfigText) &&
        sourceInputSha == zhangCheckpointSha256(targetInputText);
}

class ZhangGec134ArContinuationAuthority;
inline bool validateZhangGec134ArContinuation(
    const std::string&, const std::string&, const std::string&, const std::string&,
    const std::string&, const std::string&, const std::string&, const std::string&,
    const std::string&, const std::map<std::string, std::string>&,
    ZhangGec134ArContinuationAuthority&, std::string&);

class ZhangGec134ArContinuationAuthority
{
public:
    bool valid() const { return valid_; }
    const std::string& targetConfigText() const { return targetConfigText_; }
    bool matchesSource(const std::string& mode, const std::string& binarySha,
        const std::string& configSha, const std::string& inputSha) const
    {
        return valid_ && mode == ZHANG_GEC134_AR_MODE &&
            binarySha == ZHANG_GEC134_AR_E122_BINARY_SHA256 &&
            configSha == ZHANG_GEC134_AR_E122_CONFIG_SHA256 &&
            inputSha == ZHANG_GEC134_SOURCE_INPUT_SHA256;
    }
    bool matchesTarget(const std::string& mode, const std::string& binaryPath,
        const std::string& binarySha, const std::string& configText,
        const std::string& inputSha) const
    {
        return valid_ && mode == ZHANG_GEC134_AR_MODE &&
            binaryPath == targetBinaryPath_ && binarySha == targetBinarySha_ &&
            configText == targetConfigText_ && inputSha == ZHANG_GEC134_SOURCE_INPUT_SHA256;
    }
private:
    bool valid_ = false;
    std::string targetConfigText_, targetBinaryPath_, targetBinarySha_;
    friend bool validateZhangGec134ArContinuation(
        const std::string&, const std::string&, const std::string&, const std::string&,
        const std::string&, const std::string&, const std::string&, const std::string&,
        const std::string&, const std::map<std::string, std::string>&,
        ZhangGec134ArContinuationAuthority&, std::string&);
};

// Only one new final YAML layer and its corresponding argv element may be
// added. The original eleven included files and all 150 original argv values
// remain byte-exact; no option is re-inferred or normalised to match.
inline bool validateZhangGec134ArContinuation(
    const std::string& sourceBundleSha, const std::string& sourceBinarySha,
    const std::string& sourceConfigSha, const std::string& sourceInputSha,
    const std::string& sourceConfigText, const std::string& targetConfigText,
    const std::string& targetInputSha, const std::string& targetBinaryPath,
    const std::string& targetBinarySha,
    const std::map<std::string, std::string>& environment,
    ZhangGec134ArContinuationAuthority& authority, std::string& reason)
{
    using namespace ZhangGec134MigrationDetail;
    authority = {};
    auto fail = [&](const std::string& code){ reason = code; return false; };
    const auto binaryPath = std::filesystem::path(targetBinaryPath);
    if (sourceBundleSha != ZHANG_GEC134_AR_E122_BUNDLE_SHA256 ||
        sourceBinarySha != ZHANG_GEC134_AR_E122_BINARY_SHA256 ||
        sourceConfigSha != ZHANG_GEC134_AR_E122_CONFIG_SHA256 ||
        sourceInputSha != ZHANG_GEC134_SOURCE_INPUT_SHA256 || targetInputSha != sourceInputSha ||
        zhangCheckpointSha256(sourceConfigText) != sourceConfigSha ||
        !binaryPath.is_absolute() || binaryPath != binaryPath.lexically_normal() ||
        targetBinaryPath == ZHANG_GEC134_AR_E122_BINARY_PATH || targetBinarySha.size() != 64 ||
        !std::all_of(targetBinarySha.begin(), targetBinarySha.end(), [](unsigned char c)
            { return std::isdigit(c) || (c >= 'a' && c <= 'f'); }))
        return fail("GEC134_AR_CONTINUATION_PINNED_IDENTITY_REJECTED");
    if (!validateZhangGec134ArContinuationEnvironment(environment, reason)) return false;
    if (environment.at("ZHANG_E29_RESTORE_BUNDLE") != ZHANG_GEC134_AR_E122_PATH)
        return fail("GEC134_AR_CONTINUATION_CROSS_PROVENANCE_REQUIRES_E122");
    std::map<std::string, std::string> source, target;
    if (!canonicalFields(sourceConfigText, source, reason) ||
        !canonicalFields(targetConfigText, target, reason)) return false;
    const auto overlayRecord = target.find("config_file_11_path");
    const auto* profile = overlayRecord == target.end() ? nullptr :
        zhangGec134ArContinuationProfileForOverlay(overlayRecord->second);
    if (!profile || (!profile->candidatePath.empty() && targetBinaryPath != profile->candidatePath))
        return fail("GEC134_AR_CONTINUATION_EXACT_PROFILE_CANDIDATE_REQUIRED");
    auto take = [](std::map<std::string, std::string>& fields, const std::string& key,
                   std::string& value)
    {
        const auto found = fields.find(key);
        if (found == fields.end()) return false;
        value = found->second; fields.erase(found); return true;
    };
    std::string value;
    if (!take(target, std::string("environment_") + ZHANG_GEC134_AR_CONTINUATION_FLAG, value) ||
        value != "1" || !take(target, "gec134_ar_continuation_policy", value) ||
        value != ZHANG_GEC134_AR_CONTINUATION_POLICY)
        return fail("GEC134_AR_CONTINUATION_CANONICAL_POLICY_REJECTED");
    // Explicitly check these science/CPU invariants as well as the final exact
    // remainder comparison; neither the AR start nor the source owner changes.
    const std::map<std::string, std::string> fixed = {
        {"environment_ZHANG_MULTI_GNSS_AR_CHECKPOINT", "1"},
        {"environment_ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", "1"},
        {std::string("environment_") + ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG, "1"},
        {"gec134_local_migration_policy", ZHANG_GEC134_LOCAL_ONE_EPOCH_POLICY},
        {"environment_ZHANG_R51_RATIO_THRESHOLD", "2"},
        {"effective_R51_ratio_threshold", "2"}, {"environment_ZHANG_R51_RATIO_ONLY", "1"},
        {"environment_ZHANG_R51_EPOCH_ONLY_AR", "1"},
        {"environment_ZHANG_R51_AR_START_GPST_SECONDS", ZHANG_GEC134_AR_START},
        {"environment_OMP_NUM_THREADS", "12"}, {"environment_ZHANG_R51_FLOAT_THREADS", "12"},
        {"environment_ZHANG_FLOAT_MATRIX_THREADS", "6"},
        {"float_available_cpus", "12"}, {"float_physical_cores", "6"},
        {"effective_checkpoint_runtime_id", ZHANG_GEC134_SOURCE_RUNTIME_ID}};
    for (const auto& [key, expected] : fixed)
        if (!source.count(key) || !target.count(key) ||
            source.at(key) != expected || target.at(key) != expected)
            return fail("GEC134_AR_CONTINUATION_FIXED_POLICY_REJECTED:" + key);
    if (!take(source, "config_file_count", value) || value != "11" ||
        !take(target, "config_file_count", value) || value != "12")
        return fail("GEC134_AR_CONTINUATION_CONFIG_COUNT_REJECTED");
    std::string overlayPath, overlayBytes, overlaySize, overlaySha;
    std::size_t byteCount = 0;
    if (!take(target, "config_file_11_path", overlayPath) ||
        overlayPath != profile->overlay ||
        !take(target, "config_file_11_bytes", overlayBytes) ||
        !take(target, "config_file_11_size", overlaySize) ||
        !unsignedValue(overlaySize, byteCount) || byteCount != overlayBytes.size() ||
        !take(target, "config_file_11_sha256", overlaySha) ||
        overlaySha != zhangCheckpointSha256(overlayBytes))
        return fail("GEC134_AR_CONTINUATION_OVERLAY_RECORD_REJECTED");
    if (!take(source, "command_argument_count", value) || value != "150" ||
        !take(target, "command_argument_count", value) || value != "151")
        return fail("GEC134_AR_CONTINUATION_ARGUMENT_COUNT_REJECTED");
    for (std::size_t i = 1; i <= 150; ++i)
    {
        std::string previous, current;
        if (!take(source, "command_argument_" + std::to_string(i), previous) ||
            !take(target, "command_argument_" + std::to_string(i < 12 ? i : i + 1), current) ||
            previous != current)
            return fail("GEC134_AR_CONTINUATION_ARGUMENT_CHANGED:" + std::to_string(i));
    }
    if (!take(target, "command_argument_12", value) || value != overlayPath)
        return fail("GEC134_AR_CONTINUATION_OVERLAY_ARGUMENT_REJECTED");
    const std::string epoch = "processing_options.epoch_control.";
    const std::string ar = "processing_options.gnss_general.zhang_pppar.";
    const std::string root = profile->root;
    const std::map<std::string, std::string> leavesExpected = {
        {epoch + "end_epoch", profile->endEpoch}, {epoch + "max_epochs", std::to_string(profile->maxEpochs)},
        {ar + "checkpoint_output_directory", root + "/checkpoints"},
        {ar + "checkpoint_restore_path", ZHANG_GEC134_AR_E122_PATH},
        {ar + "product_filename", root + "/outputs/zhang_internal_products.csv"},
        {ar + "product_covariance_filename", root + "/outputs/zhang_internal_product_covariance.csv"},
        {"outputs.outputs_root", root + "/outputs"}};
    std::map<std::string, YAML::Node> leaves;
    try
    {
        const auto documents = YAML::LoadAll(overlayBytes);
        if (documents.size() != 1 || !flattenOverlay(documents.front(), "", leaves, reason) ||
            leaves.size() != leavesExpected.size() + 2)
            return fail("GEC134_AR_CONTINUATION_OVERLAY_LEAVES_REJECTED");
        for (const auto& [key, expected] : leavesExpected)
            if (!leaves.count(key) || !leaves.at(key).IsScalar() ||
                leaves.at(key).as<std::string>() != expected)
                return fail("GEC134_AR_CONTINUATION_OVERLAY_VALUE_REJECTED:" + key);
        const auto description = leaves.find("outputs.metadata.config_description");
        const auto captures = leaves.find(ar + "checkpoint_capture_epochs");
        if (description == leaves.end() || !description->second.IsScalar() ||
            description->second.as<std::string>().empty() || captures == leaves.end() ||
            !captures->second.IsSequence() || captures->second.size() != profile->captures.size())
            return fail("GEC134_AR_CONTINUATION_CAPTURE_OR_DESCRIPTION_REJECTED");
        for (std::size_t i = 0; i < profile->captures.size(); ++i)
            if (!captures->second[i].IsScalar() || captures->second[i].as<std::string>() != profile->captures[i])
                return fail("GEC134_AR_CONTINUATION_CAPTURE_OR_DESCRIPTION_REJECTED");
    }
    catch (const std::exception&) { return fail("GEC134_AR_CONTINUATION_OVERLAY_PARSE_REJECTED"); }
    const std::map<std::string, std::pair<std::string, std::string>> changed = {
        {"effective_end_epoch", {"2024-07-17T01:00:30", profile->endCanonical}},
        {"effective_max_epochs", {"122", std::to_string(profile->maxEpochs)}},
        {"effective_outputs_root", {std::string(ZHANG_GEC134_LOCAL_ONE_EPOCH_ROOT) + "/outputs", root + "/outputs"}},
        {"effective_product_filename", {std::string(ZHANG_GEC134_LOCAL_ONE_EPOCH_ROOT) + "/outputs/zhang_internal_products.csv", root + "/outputs/zhang_internal_products.csv"}},
        {"effective_product_covariance_filename", {std::string(ZHANG_GEC134_LOCAL_ONE_EPOCH_ROOT) + "/outputs/zhang_internal_product_covariance.csv", root + "/outputs/zhang_internal_product_covariance.csv"}},
        {"effective_checkpoint_output_directory", {std::string(ZHANG_GEC134_LOCAL_ONE_EPOCH_ROOT) + "/checkpoints", root + "/checkpoints"}},
        {"effective_checkpoint_restore_path", {ZHANG_GEC134_SOURCE_CHECKPOINT_PATH, ZHANG_GEC134_AR_E122_PATH}},
        {"effective_capture_epoch_count", {"1", std::to_string(profile->captures.size())}},
        {"effective_capture_epoch_0", {"2024-07-17 01:00:30", profile->captures.front()}}};
    for (const auto& [key, expected] : changed)
    {
        std::string previous, current;
        if (!take(source, key, previous) || previous != expected.first ||
            !take(target, key, current) || current != expected.second)
            return fail("GEC134_AR_CONTINUATION_EFFECTIVE_VALUE_REJECTED:" + key);
    }
    if (profile->captures.size() == 2 &&
        (!take(target, "effective_capture_epoch_1", value) || value != profile->captures[1]))
        return fail("GEC134_AR_CONTINUATION_FINAL_CAPTURE_REJECTED");
    std::error_code error;
    for (const auto& path : {std::filesystem::path(root), std::filesystem::path(root + "/outputs"),
                            std::filesystem::path(root + "/checkpoints"),
                            std::filesystem::path(root + "/restore_audit")})
        if (std::filesystem::weakly_canonical(path, error) != path || error)
            return fail("GEC134_AR_CONTINUATION_OUTPUT_SYMLINK_REJECTED");
    if (source != target) return fail("GEC134_AR_CONTINUATION_NONWHITELISTED_CHANGE");
    authority.valid_ = true;
    authority.targetConfigText_ = targetConfigText;
    authority.targetBinaryPath_ = targetBinaryPath;
    authority.targetBinarySha_ = targetBinarySha;
    reason.clear();
    return true;
}
