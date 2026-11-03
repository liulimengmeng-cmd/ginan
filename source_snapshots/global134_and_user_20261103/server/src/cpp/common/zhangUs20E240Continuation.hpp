#pragma once

// One finite source-to-target capability. This never re-labels an archive or
// relaxes the original US20/e122 restore validators.
#include "common/zhangGec134Migration.hpp"
#include "common/zhangUs20Global134Checkpoint.hpp"

inline constexpr const char* ZHANG_US20_E240_CONTINUATION_FLAG =
    "ZHANG_US20_E240_TRANSPORT_CONTINUATION";
inline constexpr const char* ZHANG_US20_E240_CONTINUATION_POLICY =
    "US20_GLOBAL134_E240_CODEC_FIX_V2";
inline constexpr const char* ZHANG_US20_E240_TARGET_ROOT =
    "/home/rx/GINAN/ar_checkpoint_codec_fix_20261103_01/run";
inline constexpr const char* ZHANG_US20_E240_TARGET_BINARY =
    "/home/rx/GINAN/ar_checkpoint_codec_fix_20261103_01/candidate/pea";
inline constexpr const char* ZHANG_US20_E240_OVERLAY =
    "/home/rx/GINAN/ar_checkpoint_codec_fix_20261103_01/run/config/"
    "100_us20_e240_codec_fix_five_hour_overlay.yaml";
inline constexpr const char* ZHANG_US20_E240_VALIDATION_ROOT =
    "/home/rx/GINAN/ar_checkpoint_codec_fix_20261103_01/validation_e241/run";
inline constexpr const char* ZHANG_US20_E240_VALIDATION_OVERLAY =
    "/home/rx/GINAN/ar_checkpoint_codec_fix_20261103_01/validation_e241/run/config/"
    "100_us20_e240_codec_e241_validation_overlay.yaml";
inline constexpr const char* ZHANG_US20_E240_SOURCE_PATH =
    "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/checkpoints/"
    "E29-20240717015930-e240-r32f626de9463";
inline constexpr const char* ZHANG_US20_E240_SOURCE_BINARY_SHA =
    "263db653de36e433d2458de7a34064ee3090fb6aa834324b5fb5e5e78ba5da16";
inline constexpr const char* ZHANG_US20_E240_SOURCE_BUNDLE_SHA =
    "abe8b436a814cd2f9dab87b979675c1328ba85aa5a1cc797cbdcc91a28189db5";
inline constexpr const char* ZHANG_US20_E240_SOURCE_MANIFEST_SHA =
    "690ca1d4559dfb312d3b30138081d8830fb9ed1c4506489d3cd540f86d34f517";
inline constexpr const char* ZHANG_US20_E240_SOURCE_CONFIG_SHA =
    "752fa1a4ef93d4f6529767d111b80bcfbe826207d698b063322f124feb23946c";
inline constexpr const char* ZHANG_US20_E240_SOURCE_INPUT_SHA =
    "e325b47f85b6061b0859bf9af7e9256c7f6bd88a1cf2830099cd4faca06f5da4";
inline constexpr const char* ZHANG_US20_E240_TARGET_DESCRIPTION =
    "US20 e240 codec repair FLOAT30m then AR270m";
// Two additional literal capabilities for the separately authorised one-hour
// experiment. No arbitrary root, binary, duration or schedule is accepted.
inline constexpr const char* ZHANG_US20_E240_DATUM_AR1H_ROOT =
    "/home/rx/GINAN/datum_unknown_ar1h_20261103_01/run";
inline constexpr const char* ZHANG_US20_E240_DATUM_AR1H_BINARY =
    "/home/rx/GINAN/datum_unknown_ar1h_20261103_01/candidate/pea";
inline constexpr const char* ZHANG_US20_E240_DATUM_AR1H_OVERLAY =
    "/home/rx/GINAN/datum_unknown_ar1h_20261103_01/run/config/100_datum_unknown_ar1h_overlay.yaml";
inline constexpr const char* ZHANG_US20_E240_DATUM_VALIDATION_ROOT =
    "/home/rx/GINAN/datum_unknown_ar1h_20261103_01/validation_e241/run";
inline constexpr const char* ZHANG_US20_E240_DATUM_VALIDATION_OVERLAY =
    "/home/rx/GINAN/datum_unknown_ar1h_20261103_01/validation_e241/run/config/"
    "100_datum_unknown_e241_validation_overlay.yaml";

struct ZhangUs20E240FiniteProfile
{
    std::string root, endEpoch, overlay, description;
    int maxEpochs;
    std::int64_t arStart;
    std::string binary = ZHANG_US20_E240_TARGET_BINARY;
};

inline const std::vector<ZhangUs20E240FiniteProfile>& zhangUs20E240FiniteProfiles()
{
    static const std::vector<ZhangUs20E240FiniteProfile> profiles = {
        {ZHANG_US20_E240_TARGET_ROOT,"2024-07-17 06:59:30",ZHANG_US20_E240_OVERLAY,
            ZHANG_US20_E240_TARGET_DESCRIPTION,840,1405218600},
        {ZHANG_US20_E240_VALIDATION_ROOT,"2024-07-17 02:00:30",ZHANG_US20_E240_VALIDATION_OVERLAY,
            "US20 e240 codec repair e241-e242 validation",242,1405216800},
        {ZHANG_US20_E240_DATUM_AR1H_ROOT,"2024-07-17 02:59:30",ZHANG_US20_E240_DATUM_AR1H_OVERLAY,
            "Datum UNKNOWN repair direct AR1h",360,1405216800,ZHANG_US20_E240_DATUM_AR1H_BINARY},
        {ZHANG_US20_E240_DATUM_VALIDATION_ROOT,"2024-07-17 02:00:30",ZHANG_US20_E240_DATUM_VALIDATION_OVERLAY,
            "Datum UNKNOWN repair e241-e242 validation",242,1405216800,ZHANG_US20_E240_DATUM_AR1H_BINARY}};
    return profiles;
}

inline const ZhangUs20E240FiniteProfile* zhangUs20E240ProfileForRoot(const std::string& root)
{
    for (const auto& profile:zhangUs20E240FiniteProfiles())
        if (root==profile.root) return &profile;
    return nullptr;
}

inline const ZhangUs20E240FiniteProfile* zhangUs20E240ProfileForOutputs(const std::string& outputs)
{
    for (const auto& profile:zhangUs20E240FiniteProfiles())
        if (outputs==profile.root+"/outputs") return &profile;
    return nullptr;
}

inline bool zhangUs20E240ContinuationRequested(
    const std::map<std::string, std::string>& environment)
{
    return zhangUs20Global134EnvironmentValue(
        environment, ZHANG_US20_E240_CONTINUATION_FLAG) == "1";
}

inline std::vector<std::string> zhangUs20E240TargetCaptures(
    const std::string& root=ZHANG_US20_E240_TARGET_ROOT)
{
    if (root==ZHANG_US20_E240_VALIDATION_ROOT || root==ZHANG_US20_E240_DATUM_VALIDATION_ROOT)
        return {"2024-07-17 02:00:00","2024-07-17 02:00:30"};
    if (root==ZHANG_US20_E240_DATUM_AR1H_ROOT)
        return {"2024-07-17 02:29:30","2024-07-17 02:59:30"};
    if (root!=ZHANG_US20_E240_TARGET_ROOT) return {};
    return {"2024-07-17 02:29:30", "2024-07-17 02:59:30",
        "2024-07-17 03:29:30", "2024-07-17 03:59:30",
        "2024-07-17 04:29:30", "2024-07-17 04:59:30",
        "2024-07-17 05:29:30", "2024-07-17 05:59:30",
        "2024-07-17 06:29:30", "2024-07-17 06:59:30"};
}

inline std::string zhangUs20E240ValidationSnapshotPath(
    const std::string& root=ZHANG_US20_E240_TARGET_ROOT)
{
    if (!zhangUs20E240ProfileForRoot(root)) return {};
    return root + "/restore_audit/ar_roundtrip_" +
        ZHANG_US20_GLOBAL134_WARMUP_CAPTURE + "-validation-only";
}

// Only the explicitly listed target checkpoints can decode as own archives.
// The source archive and the validation snapshot retain separate authorities.
inline const std::map<std::string,std::pair<int,std::string>>& zhangUs20E240TargetCheckpoints(
    const std::string& root=ZHANG_US20_E240_TARGET_ROOT)
{
    static const std::map<std::string,std::pair<int,std::string>> probe = {
        {"E29-20240717020000-e241-r32f626de9463",{241,"2024-07-17 02:00:00"}},
        {"E29-20240717020030-e242-r32f626de9463",{242,"2024-07-17 02:00:30"}}};
    static const std::map<std::string,std::pair<int,std::string>> empty;
    if (root==ZHANG_US20_E240_VALIDATION_ROOT || root==ZHANG_US20_E240_DATUM_VALIDATION_ROOT) return probe;
    static const std::map<std::string,std::pair<int,std::string>> datumHour = {
        {"E29-20240717022930-e300-r32f626de9463",{300,"2024-07-17 02:29:30"}},
        {"E29-20240717025930-e360-r32f626de9463",{360,"2024-07-17 02:59:30"}}};
    if (root==ZHANG_US20_E240_DATUM_AR1H_ROOT) return datumHour;
    if (root!=ZHANG_US20_E240_TARGET_ROOT) return empty;
    static const std::map<std::string,std::pair<int,std::string>> values = {
        {"E29-20240717022930-e300-r32f626de9463",{300,"2024-07-17 02:29:30"}},
        {"E29-20240717025930-e360-r32f626de9463",{360,"2024-07-17 02:59:30"}},
        {"E29-20240717032930-e420-r32f626de9463",{420,"2024-07-17 03:29:30"}},
        {"E29-20240717035930-e480-r32f626de9463",{480,"2024-07-17 03:59:30"}},
        {"E29-20240717042930-e540-r32f626de9463",{540,"2024-07-17 04:29:30"}},
        {"E29-20240717045930-e600-r32f626de9463",{600,"2024-07-17 04:59:30"}},
        {"E29-20240717052930-e660-r32f626de9463",{660,"2024-07-17 05:29:30"}},
        {"E29-20240717055930-e720-r32f626de9463",{720,"2024-07-17 05:59:30"}},
        {"E29-20240717062930-e780-r32f626de9463",{780,"2024-07-17 06:29:30"}},
        {"E29-20240717065930-e840-r32f626de9463",{840,"2024-07-17 06:59:30"}}};
    return values;
}

inline bool zhangUs20E240OwnCheckpointIdentity(const std::string& path,
    int& completed, std::string& epoch, std::string& id,const std::string& boundRoot={})
{
    completed=0;epoch.clear();id.clear();
    for (const auto& profile:zhangUs20E240FiniteProfiles())
    {
        if (!boundRoot.empty() && profile.root!=boundRoot) continue;
        for (const auto& [candidate,identity]:zhangUs20E240TargetCheckpoints(profile.root))
            if (path==profile.root+"/checkpoints/"+candidate)
            {completed=identity.first;epoch=identity.second;id=candidate;return true;}
    }
    return false;
}

inline bool zhangUs20E240RestorePathAllowed(const std::string& path,const std::string& boundRoot={})
{
    if (!boundRoot.empty() && !zhangUs20E240ProfileForRoot(boundRoot)) return false;
    int completed=0;std::string epoch,id;
    if (path == ZHANG_US20_E240_SOURCE_PATH) return true;
    for (const auto& profile:zhangUs20E240FiniteProfiles())
        if ((boundRoot.empty() || boundRoot==profile.root) &&
            path==zhangUs20E240ValidationSnapshotPath(profile.root)) return true;
    return zhangUs20E240OwnCheckpointIdentity(path,completed,epoch,id,boundRoot);
}

inline bool validateZhangUs20E240TargetPolicy(
    const ZhangUs20Global134CheckpointPolicy& value, std::string& reason)
{
    const auto env = [&](const char* name) {
        return zhangUs20Global134EnvironmentValue(value.environment, name);
    };
    auto fail = [&](const std::string& code) { reason = code; return false; };
    if (!zhangUs20E240ContinuationRequested(value.environment) ||
        env(ZHANG_US20_GLOBAL134_CHECKPOINT_FLAG) != "1" ||
        env(ZHANG_US20_GLOBAL134_STAGE_FLAG) != "AR")
        return fail("US20_E240_EXPLICIT_AR_CAPABILITY_REQUIRED");
    const auto* profile=zhangUs20E240ProfileForOutputs(value.outputsRoot);
    if (!profile) return fail("US20_E240_FINITE_PROFILE_REQUIRED");
    const std::string& root = profile->root;
    if (value.outputsRoot != root + "/outputs" ||
        value.binaryPath != profile->binary ||
        value.productFilename != root + "/outputs/zhang_internal_products.csv" ||
        value.covarianceFilename != root + "/outputs/zhang_internal_product_covariance.csv" ||
        value.checkpointOutputDirectory != root + "/checkpoints" ||
        value.checkpointRuntimeId != ZHANG_US20_GLOBAL134_RUNTIME_ID ||
        !zhangUs20E240RestorePathAllowed(value.checkpointRestorePath,root) ||
        env("ZHANG_E29_RESTORE_BUNDLE") != value.checkpointRestorePath)
        return fail("US20_E240_EXACT_PATH_REQUIRED");
    if (value.startEpoch != "2024-07-17 00:00:00" ||
        value.endEpoch != profile->endEpoch || value.maxEpochs != profile->maxEpochs ||
        value.epochIntervalSeconds != 30 ||
        value.arStartGpstSeconds != profile->arStart ||
        value.captureEpochs != zhangUs20E240TargetCaptures(root))
        return fail("US20_E240_EXACT_FINITE_SCHEDULE_REQUIRED");
    if (value.receiverCount != 134 ||
        value.systems != std::set<std::string>{"GPS","GAL","BDS"} ||
        !value.fixedCoordinates || value.productBasisMode != "SINGLE_S_BASIS_NATIVE" ||
        !value.physicalIonosphereEnabled || !value.physicalTroposphereEnabled ||
        !value.stecEstimated || !value.ztdEstimated || value.regionalJointEnabled ||
        value.regionalTropEnabled || value.ratio != 2 || value.jointNisEnabled ||
        !value.independentEpochAr)
        return fail("US20_E240_SCIENCE_POLICY_REQUIRED");
    const std::map<std::string,std::string> required = {
        {"ZHANG_R51_ENABLE","1"},{"ZHANG_R51_RATIO_ONLY","1"},
        {"ZHANG_R51_RATIO_THRESHOLD","2"},{"ZHANG_R51_EPOCH_ONLY_AR","1"},
        {"ZHANG_R51_AR_START_GPST_SECONDS",std::to_string(profile->arStart)},
        {"OMP_NUM_THREADS","12"},{"OPENBLAS_NUM_THREADS","1"},{"MKL_NUM_THREADS","1"},
        {"ZHANG_R51_FLOAT_THREADS","12"},{"ZHANG_FLOAT_MATRIX_THREADS","6"},
        {"ZHANG_FLOAT_MEMORY_MODE","parallel"},{"ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC","off"},
        {"ZHANG_REGIONAL_JOINT_EXPORT","0"},{"ZHANG_REGIONAL_TROPOSPHERE","OFF"}};
    for (const auto& [key, expected] : required)
        if (zhangUs20Global134EnvironmentValue(value.environment,key) != expected)
            return fail("US20_E240_ENVIRONMENT_REJECTED:" + key);
    for (const auto* name : {"ZHANG_GEC134_FLOAT_TO_AR_MIGRATION",
            "ZHANG_GEC134_AR_CONTINUATION_12_6","ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6",
            "ZHANG_GEC134_LOCAL_ONE_EPOCH_MIGRATION","ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT",
            "ZHANG_MULTI_GNSS_AR_CHECKPOINT","ZHANG_R51_REUSE_FLOAT_20240717",
            "ZHANG_R51_REUSE_RATIO_20240717"})
        if (!env(name).empty() && env(name) != "UNSET")
            return fail(std::string("US20_E240_LEGACY_CAPABILITY_FORBIDDEN:") + name);
    reason.clear(); return true;
}

inline std::string zhangUs20E240OverlayText(const std::string& root=ZHANG_US20_E240_TARGET_ROOT)
{
    const auto* profile=zhangUs20E240ProfileForRoot(root);
    if (!profile) return {};
    std::string captureText;
    for (const auto& epoch:zhangUs20E240TargetCaptures(root))
        captureText+="        - '"+epoch+"'\n";
    return "processing_options:\n  epoch_control:\n"
        "    end_epoch: '"+profile->endEpoch+"'\n    max_epochs: "+std::to_string(profile->maxEpochs)+"\n"
        "  gnss_general:\n    zhang_pppar:\n"
        "      checkpoint_output_directory: '" + root + "/checkpoints'\n"
        "      checkpoint_capture_epochs:\n"+captureText+
        "      product_filename: '" + root + "/outputs/zhang_internal_products.csv'\n"
        "      product_covariance_filename: '" + root + "/outputs/zhang_internal_product_covariance.csv'\n"
        "outputs:\n  outputs_root: '" + root + "/outputs'\n"
        "  metadata:\n    config_description: '" + profile->description + "'\n";
}

namespace ZhangUs20E240Detail
{
inline bool take(std::map<std::string,std::string>& values,
    const std::string& key, std::string& value)
{
    const auto item = values.find(key);
    if (item == values.end()) return false;
    value = item->second; values.erase(item); return true;
}

inline bool environmentDirectory(std::map<std::string,std::string>& fields,
    std::map<std::string,std::string>& result, std::string& reason)
{
    std::string countText; std::size_t count = 0;
    if (!take(fields,"us20_global134_science_environment_count",countText) ||
        !ZhangGec134MigrationDetail::unsignedValue(countText,count) || count > 512)
    { reason="US20_E240_SCIENCE_DIRECTORY_COUNT"; return false; }
    for (std::size_t i=0;i<count;++i)
    {
        const auto prefix="us20_global134_science_environment_"+std::to_string(i)+"_";
        std::string key,value;
        if (!take(fields,prefix+"key",key) || !take(fields,prefix+"value",value) ||
            key.empty() || !result.emplace(key,value).second)
        { reason="US20_E240_SCIENCE_DIRECTORY_RECORD"; return false; }
    }
    return true;
}
} // namespace ZhangUs20E240Detail

// Pure canonical-delta check; this by itself mints no restore authority.
inline bool validateZhangUs20E240CanonicalDelta(
    const std::string& sourceText, const std::string& targetText,
    const std::map<std::string,std::string>& actualEnvironment, std::string& reason)
{
    using namespace ZhangUs20E240Detail;
    std::map<std::string,std::string> source,target;
    auto fail=[&](const std::string& code){reason=code;return false;};
    if (!ZhangGec134MigrationDetail::canonicalFields(sourceText,source,reason) ||
        !ZhangGec134MigrationDetail::canonicalFields(targetText,target,reason)) return false;
    const auto working=target.find("working_directory");
    const auto* profile=working==target.end()?nullptr:zhangUs20E240ProfileForRoot(working->second);
    if (!profile) return fail("US20_E240_FINITE_CANONICAL_PROFILE_REQUIRED");
    std::string endIso=profile->endEpoch;endIso[10]='T';
    std::map<std::string,std::string> oldEnvironment,newEnvironment;
    if (!environmentDirectory(source,oldEnvironment,reason) ||
        !environmentDirectory(target,newEnvironment,reason)) return false;
    if (oldEnvironment.contains(ZHANG_US20_E240_CONTINUATION_FLAG) ||
        !newEnvironment.contains(ZHANG_US20_E240_CONTINUATION_FLAG) ||
        newEnvironment.at(ZHANG_US20_E240_CONTINUATION_FLAG)!="1" ||
        newEnvironment != zhangUs20Global134ScienceEnvironment(actualEnvironment))
        return fail("US20_E240_TARGET_SCIENCE_DIRECTORY_REJECTED");
    newEnvironment.erase(ZHANG_US20_E240_CONTINUATION_FLAG);
    // The formal extension alone postpones actual AR from 02:00 to 02:30.
    // Both copies of the explicit key are bound to the selected finite
    // profile; every other science setting remains byte-for-byte identical.
    const auto oldStart=oldEnvironment.find("ZHANG_R51_AR_START_GPST_SECONDS");
    const auto newStart=newEnvironment.find("ZHANG_R51_AR_START_GPST_SECONDS");
    if (oldStart==oldEnvironment.end() || newStart==newEnvironment.end() ||
        oldStart->second!="1405216800" || newStart->second!=std::to_string(profile->arStart))
        return fail("US20_E240_EXACT_AR_START_DELTA_REQUIRED");
    newStart->second=oldStart->second;
    if (oldEnvironment!=newEnvironment)
        return fail("US20_E240_NONWHITELISTED_SCIENCE_ENVIRONMENT");
    const std::map<std::string,std::pair<std::string,std::string>> delta = {
        {"working_directory",{"/mnt/d/GINAN_R20/inputData",profile->root}},
        {"us20_global134_producer_path",{ZHANG_US20_GLOBAL134_BINARY,profile->binary}},
        {"us20_global134_checkpoint_policy",{"US20_GLOBAL134_COLD_FLOAT2H_AR5H_V2",ZHANG_US20_E240_CONTINUATION_POLICY}},
        {"effective_end_epoch",{"2024-07-17T06:59:30",endIso}},
        {"effective_max_epochs",{"840",std::to_string(profile->maxEpochs)}},
        {"environment_ZHANG_R51_AR_START_GPST_SECONDS",{"1405216800",std::to_string(profile->arStart)}},
        {"effective_outputs_root",{std::string(ZHANG_US20_GLOBAL134_ROOT)+"/outputs",profile->root+"/outputs"}},
        {"effective_product_filename",{std::string(ZHANG_US20_GLOBAL134_ROOT)+"/outputs/zhang_internal_products.csv",profile->root+"/outputs/zhang_internal_products.csv"}},
        {"effective_product_covariance_filename",{std::string(ZHANG_US20_GLOBAL134_ROOT)+"/outputs/zhang_internal_product_covariance.csv",profile->root+"/outputs/zhang_internal_product_covariance.csv"}},
        {"effective_checkpoint_output_directory",{std::string(ZHANG_US20_GLOBAL134_ROOT)+"/checkpoints",profile->root+"/checkpoints"}}};
    for (const auto& [key,expected] : delta)
    {
        std::string before,after;
        if (!take(source,key,before) || before!=expected.first ||
            !take(target,key,after) || after!=expected.second)
            return fail("US20_E240_CANONICAL_DELTA_REJECTED:"+key);
    }
    std::string value;
    if (!take(source,"config_file_count",value) || value!="10" ||
        !take(target,"config_file_count",value) || value!="11")
        return fail("US20_E240_CONFIG_COUNT_REJECTED");
    const auto overlay=zhangUs20E240OverlayText(profile->root);
    const std::map<std::string,std::string> overlayFields = {
        {"path",profile->overlay},{"bytes",overlay},
        {"size",std::to_string(overlay.size())},{"sha256",zhangCheckpointSha256(overlay)}};
    for (const auto& [suffix,expected] : overlayFields)
        if (!take(target,"config_file_10_"+suffix,value) || value!=expected)
            return fail("US20_E240_OVERLAY_RECORD_REJECTED:"+suffix);
    if (!take(source,"command_argument_count",value) || value!="149" ||
        !take(target,"command_argument_count",value) || value!="150")
        return fail("US20_E240_ARGUMENT_COUNT_REJECTED");
    for (std::size_t i=1;i<=149;++i)
    {
        std::string before,after;
        if (!take(source,"command_argument_"+std::to_string(i),before) ||
            !take(target,"command_argument_"+std::to_string(i<11?i:i+1),after) || before!=after)
            return fail("US20_E240_ORIGINAL_ARGUMENT_CHANGED:"+std::to_string(i));
    }
    if (!take(target,"command_argument_11",value) || value!=profile->overlay)
        return fail("US20_E240_OVERLAY_ARGUMENT_REJECTED");
    const auto newCaptures=zhangUs20E240TargetCaptures(profile->root);
    if (!take(source,"effective_capture_epoch_count",value) || value!="6" ||
        !take(target,"effective_capture_epoch_count",value) || value!=std::to_string(newCaptures.size()))
        return fail("US20_E240_CAPTURE_COUNT_REJECTED");
    const auto oldCaptures=zhangUs20Global134CheckpointProfiles().front().captures;
    for (std::size_t i=0;i<oldCaptures.size();++i)
        if (!take(source,"effective_capture_epoch_"+std::to_string(i),value) || value!=oldCaptures[i])
            return fail("US20_E240_SOURCE_CAPTURE_REJECTED");
    for (std::size_t i=0;i<newCaptures.size();++i)
        if (!take(target,"effective_capture_epoch_"+std::to_string(i),value) || value!=newCaptures[i])
            return fail("US20_E240_TARGET_CAPTURE_REJECTED");
    if (source!=target) return fail("US20_E240_NONWHITELISTED_CANONICAL_CHANGE");
    reason.clear();return true;
}

struct ZhangUs20E240SourceIdentity : ZhangUs20Global134CheckpointIdentity
{
    std::string checkpointId;
    int stateDimension=0;
};

class ZhangUs20E240ContinuationAuthority;
inline bool validateZhangUs20E240Continuation(
    const std::string&,const std::string&,const ZhangUs20E240SourceIdentity&,
    const std::string&,const std::string&,const std::string&,const std::string&,
    const std::string&,const std::string&,const std::map<std::string,std::string>&,
    ZhangUs20E240ContinuationAuthority&,std::string&);

class ZhangUs20E240ContinuationAuthority
{
public:
    bool valid() const {return valid_;}
    bool matchesSource(const std::string& mode,const std::string& binarySha,
        const std::string& configSha,const std::string& inputSha) const
    {
        return valid_ && mode==ZHANG_US20_GLOBAL134_MODE &&
            binarySha==ZHANG_US20_E240_SOURCE_BINARY_SHA &&
            configSha==sourceConfigSha_ && inputSha==inputSha_;
    }
    bool matchesTarget(const std::string& mode,const std::string& binaryPath,
        const std::string& binarySha,const std::string& configText,
        const std::string& inputSha) const
    {
        return valid_ && mode==ZHANG_US20_GLOBAL134_MODE &&
            binaryPath==targetBinaryPath_ && binarySha==targetBinarySha_ &&
            configText==targetConfigText_ && inputSha==inputSha_;
    }
private:
    bool valid_=false;
    std::string sourceConfigSha_,inputSha_,targetConfigText_,targetBinarySha_,targetBinaryPath_;
    friend bool validateZhangUs20E240Continuation(
        const std::string&,const std::string&,const ZhangUs20E240SourceIdentity&,
        const std::string&,const std::string&,const std::string&,const std::string&,
        const std::string&,const std::string&,const std::map<std::string,std::string>&,
        ZhangUs20E240ContinuationAuthority&,std::string&);
};

inline bool validateZhangUs20E240Continuation(
    const std::string& bundleSha,const std::string& sidecarSha,
    const ZhangUs20E240SourceIdentity& source,const std::string& sourceConfig,
    const std::string& sourceInputs,const std::string& targetConfig,
    const std::string& targetInputs,const std::string& targetBinaryPath,
    const std::string& targetBinarySha,const std::map<std::string,std::string>& environment,
    ZhangUs20E240ContinuationAuthority& authority,std::string& reason)
{
    authority={};
    auto fail=[&](const std::string& code){reason=code;return false;};
    if (bundleSha!=ZHANG_US20_E240_SOURCE_BUNDLE_SHA ||
        sidecarSha!=ZHANG_US20_E240_SOURCE_MANIFEST_SHA ||
        source.binarySha256!=ZHANG_US20_E240_SOURCE_BINARY_SHA ||
        source.configSha256!=ZHANG_US20_E240_SOURCE_CONFIG_SHA ||
        source.inputSha256!=ZHANG_US20_E240_SOURCE_INPUT_SHA ||
        source.binaryPath!=ZHANG_US20_GLOBAL134_BINARY ||
        source.experimentMode!=ZHANG_US20_GLOBAL134_MODE ||
        source.runtimeId!=ZHANG_US20_GLOBAL134_RUNTIME_ID ||
        source.checkpointId!=ZHANG_US20_GLOBAL134_WARMUP_CAPTURE ||
        source.epoch!="2024-07-17 01:59:30" || source.completedEpoch!=240 ||
        source.nextEpoch!=241 || source.stateDimension!=7274 ||
        source.graphSchema!=5 || source.sectionCount!=9 ||
        zhangCheckpointSha256(sourceConfig)!=source.configSha256 ||
        zhangCheckpointSha256(sourceInputs)!=source.inputSha256 ||
        targetInputs!=sourceInputs ||
        !zhangUs20Global134Sha256WellFormed(targetBinarySha) ||
        targetBinarySha==source.binarySha256)
        return fail("US20_E240_PINNED_SOURCE_OR_TARGET_REJECTED");
    if (!zhangUs20E240ContinuationRequested(environment) ||
        zhangUs20Global134EnvironmentValue(environment,"ZHANG_E29_RESTORE_BUNDLE")!=ZHANG_US20_E240_SOURCE_PATH ||
        !validateZhangUs20E240CanonicalDelta(sourceConfig,targetConfig,environment,reason)) return false;
    std::map<std::string,std::string> targetFields;
    if (!ZhangGec134MigrationDetail::canonicalFields(targetConfig,targetFields,reason)) return false;
    const auto working=targetFields.find("working_directory");
    const auto* profile=working==targetFields.end()?nullptr:zhangUs20E240ProfileForRoot(working->second);
    if (!profile) return fail("US20_E240_FINITE_TARGET_ROOT_REQUIRED");
    if (targetBinaryPath!=profile->binary)
        return fail("US20_E240_SELECTED_PROFILE_BINARY_PATH_REQUIRED");
    std::error_code error;
    const std::filesystem::path root(profile->root);
    for (const auto& path : {root,root/"outputs",root/"checkpoints",root/"restore_audit"})
    {
        error.clear();
        if (std::filesystem::weakly_canonical(path,error)!=path || error)
            return fail("US20_E240_TARGET_PATH_SYMLINK_REJECTED");
    }
    authority.valid_=true;authority.sourceConfigSha_=source.configSha256;
    authority.inputSha_=source.inputSha256;authority.targetConfigText_=targetConfig;
    authority.targetBinarySha_=targetBinarySha;authority.targetBinaryPath_=targetBinaryPath;
    reason.clear();return true;
}
