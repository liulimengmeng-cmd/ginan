#pragma once

// A finite, cold-start experiment policy, independent of ACSConfig and KFState.
// Returning true authorises only the stated CONFIGURATION. It is not evidence
// that a checkpoint was decoded, its nine sections imported, or AR was correct.
// This policy never grants a legacy/e122 cross-binary migration capability.
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

inline constexpr const char* ZHANG_US20_GLOBAL134_CHECKPOINT_FLAG =
    "ZHANG_US20_GLOBAL134_CHECKPOINT";
inline constexpr const char* ZHANG_US20_GLOBAL134_STAGE_FLAG =
    "ZHANG_US20_GLOBAL134_STAGE";
inline constexpr const char* ZHANG_US20_GLOBAL134_ROOT =
    "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run";
inline constexpr const char* ZHANG_US20_GLOBAL134_PROBE_ROOT =
    "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/probe/run";
inline constexpr const char* ZHANG_US20_GLOBAL134_BINARY =
    "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/candidate/pea";
// Existing nine-section runtime schema has this literal internal identity.
// It does not substitute for a producer's binary/config/input hashes.
inline constexpr const char* ZHANG_US20_GLOBAL134_RUNTIME_ID =
    "E29-GEC134-FIXEDCOORD-20261023-one_hour";
inline constexpr const char* ZHANG_US20_GLOBAL134_MODE =
    "E29_MULTI_GNSS_AR_V1_GPS_GAL_BDS";
inline constexpr const char* ZHANG_US20_GLOBAL134_WARMUP_CAPTURE =
    "E29-20240717015930-e240-r32f626de9463";
inline constexpr const char* ZHANG_US20_GLOBAL134_PROBE_CAPTURE =
    "E29-20240717000000-e1-r32f626de9463";
inline constexpr std::int64_t ZHANG_US20_GLOBAL134_AR_START = 1405216800;

struct ZhangUs20Global134CheckpointProfile
{
    std::string root;
    std::string endEpoch;
    int maxEpochs;
    std::vector<std::string> captures;
};

inline const std::vector<ZhangUs20Global134CheckpointProfile>&
zhangUs20Global134CheckpointProfiles()
{
    static const std::vector<ZhangUs20Global134CheckpointProfile> profiles = {
        {ZHANG_US20_GLOBAL134_ROOT, "2024-07-17 06:59:30", 840,
         {"2024-07-17 01:59:30", "2024-07-17 02:59:30",
          "2024-07-17 03:59:30", "2024-07-17 04:59:30",
          "2024-07-17 05:59:30", "2024-07-17 06:59:30"}},
        {ZHANG_US20_GLOBAL134_PROBE_ROOT, "2024-07-17 00:00:00", 1,
         {"2024-07-17 00:00:00"}}
    };
    return profiles;
}

inline const ZhangUs20Global134CheckpointProfile*
zhangUs20Global134CheckpointProfileForOutputs(const std::string& outputsRoot)
{
    for (const auto& profile : zhangUs20Global134CheckpointProfiles())
        if (outputsRoot == profile.root + "/outputs") return &profile;
    return nullptr;
}

inline std::string zhangUs20Global134EnvironmentValue(
    const std::map<std::string, std::string>& environment, const std::string& key)
{
    const auto found = environment.find(key);
    return found == environment.end() ? std::string{} : found->second;
}

inline bool zhangUs20Global134CheckpointRequested(
    const std::map<std::string, std::string>& environment)
{
    return zhangUs20Global134EnvironmentValue(
        environment, ZHANG_US20_GLOBAL134_CHECKPOINT_FLAG) == "1";
}

inline bool zhangUs20Global134ScienceEnvironmentKey(const std::string& key)
{
    return key.rfind("ZHANG_", 0) == 0 || key.rfind("OMP_", 0) == 0 ||
        key.rfind("OPENBLAS_", 0) == 0 || key.rfind("MKL_", 0) == 0 ||
        key == "GOTO_NUM_THREADS" || key == "BLIS_NUM_THREADS" ||
        key == "VECLIB_MAXIMUM_THREADS";
}

// The caller supplies the ACTUAL process environment. Preserve every scoped
// value, including unrecognised future switches and present-but-empty values.
// A missing entry is different from an entry containing "UNSET". The three
// exact lifecycle controls alone may differ between cold warmup and restoration;
// their values are still separately constrained by validate...Policy below.
inline std::map<std::string, std::string> zhangUs20Global134ScienceEnvironment(
    const std::map<std::string, std::string>& actualEnvironment)
{
    std::map<std::string, std::string> science;
    for (const auto& [key, value] : actualEnvironment)
        if (zhangUs20Global134ScienceEnvironmentKey(key) &&
            key != ZHANG_US20_GLOBAL134_STAGE_FLAG &&
            key != "ZHANG_E29_RESTORE_BUNDLE" && key != "ZHANG_GEC134_RESTORE_ONLY")
            science.emplace(key, value);
    return science;
}

// Root integration must fill these fields from actual effective configuration
// and the actual process environment; do not fill them from intended settings.
struct ZhangUs20Global134CheckpointPolicy
{
    std::string outputsRoot, binaryPath, startEpoch, endEpoch;
    std::string productFilename, covarianceFilename, checkpointOutputDirectory;
    // Actual startup restore directory, resolved from the stage environment.
    // The frozen YAML restore field remains empty in both stages; never change
    // the canonical model merely to move WARMUP -> AR.
    std::string checkpointRestorePath, checkpointRuntimeId, productBasisMode;
    std::set<std::string> systems;
    std::vector<std::string> captureEpochs;
    std::map<std::string, std::string> environment;
    int receiverCount = 0, epochIntervalSeconds = 0, maxEpochs = 0;
    std::int64_t arStartGpstSeconds = 0;
    double ratio = 0;
    bool jointNisEnabled = true;
    bool independentEpochAr = false, fixedCoordinates = false;
    bool physicalIonosphereEnabled = false, physicalTroposphereEnabled = false;
    bool stecEstimated = false, ztdEstimated = false;
    bool regionalJointEnabled = true, regionalTropEnabled = true;
};

inline bool validateZhangUs20Global134CheckpointPolicy(
    const ZhangUs20Global134CheckpointPolicy& value, std::string& reason)
{
    auto fail = [&](const std::string& code) { reason = code; return false; };
    const auto env = [&](const char* key) {
        return zhangUs20Global134EnvironmentValue(value.environment, key);
    };
    if (!zhangUs20Global134CheckpointRequested(value.environment))
        return fail("US20_GLOBAL134_EXPLICIT_CHECKPOINT_FLAG_REQUIRED");
    const auto* profile = zhangUs20Global134CheckpointProfileForOutputs(value.outputsRoot);
    if (!profile) return fail("US20_GLOBAL134_OUTPUT_PROFILE_REJECTED");
    if (value.binaryPath != ZHANG_US20_GLOBAL134_BINARY)
        return fail("US20_GLOBAL134_BINARY_PATH_REJECTED");
    if (value.startEpoch != "2024-07-17 00:00:00" ||
        value.endEpoch != profile->endEpoch || value.epochIntervalSeconds != 30 ||
        value.maxEpochs != profile->maxEpochs ||
        value.arStartGpstSeconds != ZHANG_US20_GLOBAL134_AR_START ||
        value.captureEpochs != profile->captures)
        return fail("US20_GLOBAL134_EXACT_SCHEDULE_REQUIRED");
    if (value.productFilename != profile->root + "/outputs/zhang_internal_products.csv" ||
        value.covarianceFilename != profile->root + "/outputs/zhang_internal_product_covariance.csv" ||
        value.checkpointOutputDirectory != profile->root + "/checkpoints" ||
        value.checkpointRuntimeId != ZHANG_US20_GLOBAL134_RUNTIME_ID)
        return fail("US20_GLOBAL134_PRODUCT_OR_CHECKPOINT_PATH_REJECTED");
    if (value.receiverCount != 134 || value.systems != std::set<std::string>{"GPS", "GAL", "BDS"} ||
        !value.fixedCoordinates || value.productBasisMode != "SINGLE_S_BASIS_NATIVE")
        return fail("US20_GLOBAL134_NETWORK_OR_BASIS_REJECTED");
    if (!std::isfinite(value.ratio) || value.ratio != 2 || value.jointNisEnabled ||
        !value.independentEpochAr || env("ZHANG_R51_ENABLE") != "1" ||
        env("ZHANG_R51_RATIO_ONLY") != "1" ||
        env("ZHANG_R51_RATIO_THRESHOLD") != "2" || env("ZHANG_R51_EPOCH_ONLY_AR") != "1" ||
        env("ZHANG_R51_AR_START_GPST_SECONDS") != "1405216800")
        return fail("US20_GLOBAL134_AR_POLICY_REJECTED");
    if (!value.physicalIonosphereEnabled || !value.physicalTroposphereEnabled ||
        !value.stecEstimated || !value.ztdEstimated ||
        value.regionalJointEnabled || value.regionalTropEnabled)
        return fail("US20_GLOBAL134_ATMOSPHERE_POLICY_REJECTED");
    if ((!env("ZHANG_REGIONAL_JOINT_EXPORT").empty() &&
         env("ZHANG_REGIONAL_JOINT_EXPORT") != "0" && env("ZHANG_REGIONAL_JOINT_EXPORT") != "UNSET") ||
        (!env("ZHANG_REGIONAL_TROPOSPHERE").empty() &&
         env("ZHANG_REGIONAL_TROPOSPHERE") != "OFF" && env("ZHANG_REGIONAL_TROPOSPHERE") != "UNSET"))
        return fail("US20_GLOBAL134_REGIONAL_EXPORT_REJECTED");
    for (const auto* forbidden : {"ZHANG_GEC134_FLOAT_TO_AR_MIGRATION",
            "ZHANG_GEC134_AR_CONTINUATION_12_6", "ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6",
            "ZHANG_GEC134_LOCAL_ONE_EPOCH_MIGRATION",
            "ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "ZHANG_MULTI_GNSS_AR_CHECKPOINT",
            "ZHANG_R51_REUSE_FLOAT_20240717", "ZHANG_R51_REUSE_RATIO_20240717"})
    {
        const auto setting = env(forbidden);
        if (!setting.empty() && setting != "UNSET")
            return fail(std::string("US20_GLOBAL134_LEGACY_CAPABILITY_FORBIDDEN:") + forbidden);
    }
    const auto stage = env(ZHANG_US20_GLOBAL134_STAGE_FLAG);
    const auto restoreEnvironment = env("ZHANG_E29_RESTORE_BUNDLE");
    const bool environmentHasRestore = !restoreEnvironment.empty() && restoreEnvironment != "UNSET";
    if (stage == "WARMUP")
    {
        if (profile->root == ZHANG_US20_GLOBAL134_ROOT)
        {
            if (!value.checkpointRestorePath.empty() || environmentHasRestore)
                return fail("US20_GLOBAL134_WARMUP_MUST_BE_COLD");
        }
        else
        {
            const std::string probeCheckpoint = profile->root + "/checkpoints/" +
                ZHANG_US20_GLOBAL134_PROBE_CAPTURE;
            if ((!value.checkpointRestorePath.empty() && value.checkpointRestorePath != probeCheckpoint) ||
                (environmentHasRestore && restoreEnvironment != probeCheckpoint) ||
                (environmentHasRestore != !value.checkpointRestorePath.empty()))
                return fail("US20_GLOBAL134_PROBE_RESTORE_PATH_REJECTED");
        }
    }
    else if (stage == "AR")
    {
        const std::string expected = std::string(ZHANG_US20_GLOBAL134_ROOT) + "/checkpoints/" +
            ZHANG_US20_GLOBAL134_WARMUP_CAPTURE;
        if (profile->root != ZHANG_US20_GLOBAL134_ROOT ||
            value.checkpointRestorePath != expected || !environmentHasRestore ||
            restoreEnvironment != expected)
            return fail("US20_GLOBAL134_AR_REQUIRES_OWN_E240");
    }
    else return fail("US20_GLOBAL134_STAGE_REJECTED");
    reason.clear();
    return true;
}

inline bool zhangUs20Global134WarmupTerminal(
    const ZhangUs20Global134CheckpointPolicy& policy, const std::string& completedEpoch,
    int completedEpochNumber)
{
    std::string reason;
    if (!validateZhangUs20Global134CheckpointPolicy(policy, reason) ||
        zhangUs20Global134EnvironmentValue(policy.environment, ZHANG_US20_GLOBAL134_STAGE_FLAG) != "WARMUP")
        return false;
    const bool probe = policy.outputsRoot == std::string(ZHANG_US20_GLOBAL134_PROBE_ROOT) + "/outputs";
    return probe ? completedEpoch == "2024-07-17 00:00:00" && completedEpochNumber == 1 :
        completedEpoch == "2024-07-17 01:59:30" && completedEpochNumber == 240;
}

struct ZhangUs20Global134CheckpointIdentity
{
    std::string binarySha256, configSha256, inputSha256, binaryPath;
    std::string experimentMode, runtimeId, epoch;
    int completedEpoch = 0, nextEpoch = 0, graphSchema = 0, sectionCount = 0;
};

inline bool zhangUs20Global134Sha256WellFormed(const std::string& value)
{
    if (value.size() != 64) return false;
    for (const char c : value)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}

// Metadata admission ONLY. The caller must still authenticate the bundle,
// compare the canonical config/input bytes and import every section. Passing
// this function cannot waive the existing full-Q/graph/runtime exact audits.
inline bool validateZhangUs20Global134StrictRestoreIdentity(
    const ZhangUs20Global134CheckpointPolicy& target,
    const ZhangUs20Global134CheckpointIdentity& source,
    const ZhangUs20Global134CheckpointIdentity& expectedProducer, std::string& reason)
{
    if (!validateZhangUs20Global134CheckpointPolicy(target, reason)) return false;
    auto fail = [&](const std::string& code) { reason = code; return false; };
    const bool probe = target.outputsRoot == std::string(ZHANG_US20_GLOBAL134_PROBE_ROOT) + "/outputs";
    const auto expectedEpoch = probe ? "2024-07-17 00:00:00" : "2024-07-17 01:59:30";
    const int number = probe ? 1 : 240;
    if (target.checkpointRestorePath.empty() ||
        source.experimentMode != ZHANG_US20_GLOBAL134_MODE ||
        source.runtimeId != ZHANG_US20_GLOBAL134_RUNTIME_ID ||
        source.epoch != expectedEpoch || source.completedEpoch != number ||
        source.nextEpoch != number + 1 || source.graphSchema != 5 || source.sectionCount != 9)
        return fail("US20_GLOBAL134_SOURCE_CHECKPOINT_IDENTITY_REJECTED");
    if (!zhangUs20Global134Sha256WellFormed(source.binarySha256) ||
        !zhangUs20Global134Sha256WellFormed(source.configSha256) ||
        !zhangUs20Global134Sha256WellFormed(source.inputSha256) ||
        source.binaryPath != ZHANG_US20_GLOBAL134_BINARY ||
        expectedProducer.binaryPath != ZHANG_US20_GLOBAL134_BINARY ||
        source.binarySha256 != expectedProducer.binarySha256 ||
        source.configSha256 != expectedProducer.configSha256 ||
        source.inputSha256 != expectedProducer.inputSha256)
        return fail("US20_GLOBAL134_SAME_PRODUCER_REQUIRED");
    reason.clear();
    return true;
}
