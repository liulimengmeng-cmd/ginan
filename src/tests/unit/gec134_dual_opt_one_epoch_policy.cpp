// Independent policy fixtures for one real e122 -> e123 A/B replay. Read only
// checkpoint_manifest.json; do not load the bundle, PEA, KF or runtime importer.
// This is a capability test, not evidence that a real epoch/restoration passed.
#include "common/zhangGec134ArContinuation.hpp"
#include <openssl/evp.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

std::string zhangCheckpointSha256(const std::string& input)
{
    unsigned char digest[EVP_MAX_MD_SIZE]; unsigned int size = 0;
    if (!EVP_Digest(input.data(), input.size(), digest, &size, EVP_sha256(), nullptr) || size != 32)
        throw std::runtime_error("Independent standard SHA-256 failed");
    std::ostringstream out;
    for (unsigned int i = 0; i < size; ++i)
        out << std::hex << std::setw(2) << std::setfill('0') << unsigned(digest[i]);
    return out.str();
}

namespace
{
using Fields = std::map<std::string, std::string>;
std::size_t checks = 0;
void require(bool value, const std::string& message)
{
    ++checks;
    if (!value) throw std::runtime_error(message);
}
bool scalar(const std::string& key)
{
    return key == "command_argument_count" || key == "config_file_count" ||
        key == "effective_capture_epoch_count" ||
        (key.rfind("config_file_", 0) == 0 && key.ends_with("_size"));
}
// Independent length-prefix codec: deliberately not the production parser.
Fields decode(const std::string& text)
{
    std::istringstream input(text); std::string line;
    require(bool(std::getline(input, line)) && line == "E29_CANONICAL_CONFIG_V1", "Fixture header");
    Fields result;
    while (std::getline(input, line))
    {
        const auto at = line.find('=');
        require(at != std::string::npos, "Fixture record");
        const auto key = line.substr(0, at), length = line.substr(at + 1);
        std::size_t consumed = 0;
        const auto size = std::stoull(length, &consumed);
        require(consumed == length.size() && size < 16 * 1024 * 1024, "Fixture length");
        std::string value = length;
        if (!scalar(key))
        {
            value.resize(size); input.read(value.data(), size);
            require(std::size_t(input.gcount()) == size && input.get() == '\n', "Fixture bytes");
        }
        require(result.emplace(key, value).second, "Fixture duplicate");
    }
    return result;
}
std::string encode(const Fields& fields)
{
    std::ostringstream result; result << "E29_CANONICAL_CONFIG_V1\n";
    for (const auto& [key, value] : fields)
    {
        result << key << '=' << (scalar(key) ? value : std::to_string(value.size())) << '\n';
        if (!scalar(key)) result << value << '\n';
    }
    return result.str();
}

struct Fixture
{
    std::string root, overlay, end, canonicalEnd, binary;
    int counter;
    long double seconds;
    std::vector<std::string> captures;
};
const std::vector<Fixture> oldFixtures = {
    {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_hour_memory_20261024_01",
     "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_hour_memory_20261024_01/config/100_gec134_ar_e122_hour_continuation_overlay.yaml",
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
     {"2024-07-17 01:30:30", "2024-07-17 02:00:30"}}};
const std::vector<Fixture> replayFixtures = {
    {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_baseline_20261025_01",
     "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_baseline_20261025_01/config/100_gec134_ar_e122_dual_opt_baseline_20261025_01.yaml",
     "2024-07-17 01:01:00", "2024-07-17T01:01:00",
     "/home/rx/GINAN/local_float_ar_dual_opt_replay_baseline_20261025_v1/candidate/pea", 123, 1405213260.0L,
     {"2024-07-17 01:01:00"}},
    {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_candidate_20261025_01",
     "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_dual_opt_candidate_20261025_01/config/100_gec134_ar_e122_dual_opt_candidate_20261025_01.yaml",
     "2024-07-17 01:01:00", "2024-07-17T01:01:00",
     "/home/rx/GINAN/local_float_ar_dual_opt_replay_candidate_20261025_v1/candidate/pea", 123, 1405213260.0L,
     {"2024-07-17 01:01:00"}}};

std::string overlay(const Fixture& f)
{
    std::string captures;
    for (const auto& stamp : f.captures) captures += "        - \"" + stamp + "\"\n";
    return "processing_options:\n  epoch_control:\n    end_epoch: \"" + f.end +
        "\"\n    max_epochs: " + std::to_string(f.counter) + "\n  gnss_general:\n    zhang_pppar:\n"
        "      checkpoint_capture_epochs:\n" + captures +
        "      checkpoint_restore_path: " + ZHANG_GEC134_AR_E122_PATH +
        "\n      checkpoint_output_directory: " + f.root + "/checkpoints\n"
        "      product_filename: " + f.root + "/outputs/zhang_internal_products.csv\n"
        "      product_covariance_filename: " + f.root + "/outputs/zhang_internal_product_covariance.csv\n"
        "outputs:\n  outputs_root: " + f.root + "/outputs\n"
        "  metadata:\n    config_description: independent literal one epoch capability\n";
}
void setOverlay(Fields& fields, const std::string& bytes)
{
    fields["config_file_11_bytes"] = bytes;
    fields["config_file_11_size"] = std::to_string(bytes.size());
    fields["config_file_11_sha256"] = zhangCheckpointSha256(bytes);
}
Fields target(const Fields& source, const Fixture& f)
{
    Fields result = source;
    result[std::string("environment_") + ZHANG_GEC134_AR_CONTINUATION_FLAG] = "1";
    result["gec134_ar_continuation_policy"] = ZHANG_GEC134_AR_CONTINUATION_POLICY;
    result["config_file_count"] = "12";
    result["config_file_11_path"] = f.overlay;
    result["command_argument_count"] = "151";
    for (int i = 150; i >= 12; --i)
        result["command_argument_" + std::to_string(i + 1)] = source.at("command_argument_" + std::to_string(i));
    result["command_argument_12"] = f.overlay;
    result["effective_end_epoch"] = f.canonicalEnd;
    result["effective_max_epochs"] = std::to_string(f.counter);
    result["effective_outputs_root"] = f.root + "/outputs";
    result["effective_product_filename"] = f.root + "/outputs/zhang_internal_products.csv";
    result["effective_product_covariance_filename"] = f.root + "/outputs/zhang_internal_product_covariance.csv";
    result["effective_checkpoint_output_directory"] = f.root + "/checkpoints";
    result["effective_checkpoint_restore_path"] = ZHANG_GEC134_AR_E122_PATH;
    result["effective_capture_epoch_count"] = std::to_string(f.captures.size());
    for (std::size_t i = 0; i < f.captures.size(); ++i)
        result["effective_capture_epoch_" + std::to_string(i)] = f.captures[i];
    setOverlay(result, overlay(f));
    return result;
}
Fields environment()
{
    return {
        {"ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", "1"}, {"ZHANG_MULTI_GNSS_AR_CHECKPOINT", "1"},
        {"ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "UNSET"}, {"ZHANG_R51_REUSE_FLOAT_20240717", "UNSET"},
        {"ZHANG_R51_REUSE_RATIO_20240717", "UNSET"}, {"ZHANG_R51_ENABLE", "1"},
        {"ZHANG_R51_RATIO_ONLY", "1"}, {"ZHANG_R51_RATIO_THRESHOLD", "2"},
        {"ZHANG_R51_EPOCH_ONLY_AR", "1"}, {"ZHANG_R51_AR_START_GPST_SECONDS", "1405213230"},
        {"ZHANG_R51_FLOAT_THREADS", "12"}, {"OMP_NUM_THREADS", "12"}, {"OPENBLAS_NUM_THREADS", "1"},
        {"MKL_NUM_THREADS", "1"}, {"ZHANG_FLOAT_MATRIX_THREADS", "6"}, {"ZHANG_FLOAT_MEMORY_MODE", "parallel"},
        {"ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC", "off"}, {"ZHANG_FLOAT_PROFILE", "1"},
        {"ZHANG_P0_RESOURCE_PROBE", "1"}, {"ZHANG_R49_FUSION_WEIGHT", "0.10"},
        {"ZHANG_R51_SINGLE_BASIS_HANDOFF", "1"}, {"ZHANG_R51_PRODUCT_REPAIR", "1"},
        {"ZHANG_R51_INTERN_HANDOFF_PROOFS", "1"}, {"ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6", "1"},
        {"ZHANG_GEC134_AR_CONTINUATION_12_6", "1"}, {"ZHANG_E29_RESTORE_BUNDLE", ZHANG_GEC134_AR_E122_PATH}};
}
void exactProfile(const Fixture& f, std::size_t index)
{
    const auto& p = zhangGec134ArContinuationProfiles().at(index);
    require(p.root == f.root && p.overlay == f.overlay && p.endEpoch == f.end &&
        p.endCanonical == f.canonicalEnd && p.candidatePath == f.binary && p.maxEpochs == f.counter &&
        p.endBigTime == f.seconds && p.captures == f.captures, "Literal profile field changed");
    require(zhangGec134ArContinuationProfileForOverlay(f.overlay) == &p &&
        zhangGec134ArContinuationProfileForOutputs(f.root + "/outputs") == &p,
        "Overlay/output route differs");
}
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected real e122 checkpoint_manifest.json path");
        const auto manifest = YAML::LoadFile(argv[1]);
        const auto sourceText = manifest["config_text"].as<std::string>();
        const auto inputText = manifest["input_manifest_text"].as<std::string>();
        require(zhangCheckpointSha256(sourceText) == ZHANG_GEC134_AR_E122_CONFIG_SHA256 &&
            zhangCheckpointSha256(inputText) == ZHANG_GEC134_SOURCE_INPUT_SHA256 &&
            manifest["experiment_mode"].as<std::string>() == ZHANG_GEC134_AR_MODE &&
            manifest["checkpoint_id"].as<std::string>() == ZHANG_GEC134_AR_E122_ID,
            "Fixture is not exact original AR e122");
        const auto source = decode(sourceText);
        const auto env = environment();
        const std::string candidateSha(64, 'a'), otherSha(64, 'b');
        ZhangGec134ArContinuationAuthority authority; std::string reason;
        auto validate = [&](const Fields& fields, const Fields& settings, const std::string& binary,
            const std::string& bundle = ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
            const std::string& sourceBinary = ZHANG_GEC134_AR_E122_BINARY_SHA256,
            const std::string& sourceConfig = ZHANG_GEC134_AR_E122_CONFIG_SHA256,
            const std::string& sourceInput = ZHANG_GEC134_SOURCE_INPUT_SHA256,
            const std::string& targetInput = ZHANG_GEC134_SOURCE_INPUT_SHA256)
        {
            return validateZhangGec134ArContinuation(bundle, sourceBinary, sourceConfig, sourceInput,
                sourceText, encode(fields), targetInput, binary, candidateSha, settings, authority, reason);
        };
        require(zhangGec134ArContinuationProfiles().size() == 8, "Unexpected authority surface");
        for (std::size_t i = 0; i < oldFixtures.size(); ++i)
        {
            const auto& f = oldFixtures[i]; exactProfile(f, i);
            const auto fields = target(source, f);
            const auto binary = f.binary.empty() ? "/home/rx/GINAN/local_real_epoch_memory_20261024_v5/candidate/pea" : f.binary;
            require(validate(fields, env, binary), "Old positive changed: " + reason);
            auto wrong = env; wrong["ZHANG_R51_RATIO_THRESHOLD"] = "3";
            require(!validate(fields, wrong, binary) && !authority.valid(), "Old ratio gate weakened");
            wrong = env; wrong["ZHANG_R51_JOINT_NIS"] = "1";
            require(!validate(fields, wrong, binary) && !authority.valid(), "Old NIS gate weakened");
            auto changed = fields; changed["effective_max_epochs"] += "0";
            require(!validate(changed, env, binary) && !authority.valid(), "Old schedule gate weakened");
        }
        for (std::size_t i = 0; i < replayFixtures.size(); ++i)
        {
            const auto& f = replayFixtures[i]; const auto& other = replayFixtures[1 - i];
            exactProfile(f, oldFixtures.size() + i);
            const auto fields = target(source, f); const auto text = encode(fields);
            require(validate(fields, env, f.binary), "New exact positive rejected: " + reason);
            require(authority.matchesSource(ZHANG_GEC134_AR_MODE, ZHANG_GEC134_AR_E122_BINARY_SHA256,
                ZHANG_GEC134_AR_E122_CONFIG_SHA256, ZHANG_GEC134_SOURCE_INPUT_SHA256), "Exact e122 token lost");
            require(authority.matchesTarget(ZHANG_GEC134_AR_MODE, f.binary, candidateSha, text,
                ZHANG_GEC134_SOURCE_INPUT_SHA256), "Exact target token lost");
            require(!authority.matchesTarget(ZHANG_GEC134_AR_MODE, f.binary, otherSha, text,
                ZHANG_GEC134_SOURCE_INPUT_SHA256), "Candidate SHA token unbound");
            require(!authority.matchesTarget(ZHANG_GEC134_AR_MODE, other.binary, candidateSha, text,
                ZHANG_GEC134_SOURCE_INPUT_SHA256), "Cross-role token accepted");
            auto reject = [&](const Fields& changed, const Fields& settings, const std::string& label)
            {
                require(validate(fields, env, f.binary), "Warm token failed");
                require(!validate(changed, settings, f.binary) && !authority.valid() &&
                    authority.targetConfigText().empty() && !reason.empty(), "Invalid request accepted: " + label);
            };
            for (const auto& [key, value] : fields)
            {
                auto changed = fields; changed[key] = value + "!";
                reject(changed, env, "changed canonical " + key);
                changed = fields; changed.erase(key); reject(changed, env, "missing canonical " + key);
            }
            for (const auto& [key, value] : env)
            {
                auto changed = env; changed[key] = value + "!";
                reject(fields, changed, "changed environment " + key);
                changed = env; changed.erase(key); reject(fields, changed, "missing environment " + key);
            }
            for (const auto& wrongBinary : {other.binary, oldFixtures[3].binary, f.binary + ".other",
                    std::string("/tmp/pea"), std::string(ZHANG_GEC134_AR_E122_BINARY_PATH)})
                require(!validate(fields, env, wrongBinary) && !authority.valid(), "Other candidate path accepted");
            for (const auto& stamp : {std::string("2024-07-17 01:00:30"), std::string("2024-07-17 01:01:30")})
            {
                auto changedFixture = f; changedFixture.end = stamp; changedFixture.captures = {stamp};
                auto changed = fields; setOverlay(changed, overlay(changedFixture));
                reject(changed, env, "other epoch/capture");
            }
            auto changedFixture = f; changedFixture.counter = 124;
            auto changed = fields; setOverlay(changed, overlay(changedFixture)); reject(changed, env, "two epochs");
            changedFixture = f; changedFixture.root += "_other";
            changed = fields; setOverlay(changed, overlay(changedFixture)); reject(changed, env, "other output root");
            changed = fields; setOverlay(changed, overlay(f) + "forbidden: true\n"); reject(changed, env, "extra overlay leaf");
            changed = fields; changed["effective_capture_epoch_1"] = "2024-07-17 01:01:30";
            reject(changed, env, "extra canonical capture");
            require(!validate(fields, env, f.binary, otherSha) && !authority.valid(), "Other bundle accepted");
            require(!validate(fields, env, f.binary, ZHANG_GEC134_AR_E122_BUNDLE_SHA256, otherSha) &&
                !authority.valid(), "Other source binary accepted");
            require(!validate(fields, env, f.binary, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
                ZHANG_GEC134_AR_E122_BINARY_SHA256, otherSha) && !authority.valid(), "Other source config accepted");
            require(!validate(fields, env, f.binary, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
                ZHANG_GEC134_AR_E122_BINARY_SHA256, ZHANG_GEC134_AR_E122_CONFIG_SHA256, otherSha) &&
                !authority.valid(), "Other source input accepted");
            require(!validate(fields, env, f.binary, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
                ZHANG_GEC134_AR_E122_BINARY_SHA256, ZHANG_GEC134_AR_E122_CONFIG_SHA256,
                ZHANG_GEC134_SOURCE_INPUT_SHA256, otherSha) && !authority.valid(), "Other target input accepted");
            const std::string captureId = "E29-20240717010100-e123-r32f626de9463";
            auto own = env; own["ZHANG_GEC134_RESTORE_ONLY"] = "1";
            own["ZHANG_E29_RESTORE_BUNDLE"] = f.root + "/checkpoints/" + captureId;
            require(validateZhangGec134ArContinuationEnvironment(own, reason), "Own exact e123 environment rejected");
            reject(fields, own, "own capture must not acquire e122 cross-provenance authority");
            require(zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, candidateSha,
                zhangCheckpointSha256(text), zhangCheckpointSha256(inputText), text, inputText,
                ZHANG_GEC134_AR_MODE, candidateSha, text, inputText), "Own ordinary strict routing rejected");
            require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, candidateSha,
                zhangCheckpointSha256(text), zhangCheckpointSha256(inputText), text, inputText,
                ZHANG_GEC134_AR_MODE, otherSha, text, inputText), "Cross-binary ordinary restore accepted");
            const auto otherText = encode(target(source, other));
            require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, candidateSha,
                zhangCheckpointSha256(text), zhangCheckpointSha256(inputText), text, inputText,
                ZHANG_GEC134_AR_MODE, candidateSha, otherText, inputText), "Cross-role config restore accepted");
            require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, candidateSha,
                otherSha, zhangCheckpointSha256(inputText), text, inputText,
                ZHANG_GEC134_AR_MODE, candidateSha, text, inputText), "Wrong source config SHA accepted as strict");
            require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, candidateSha,
                zhangCheckpointSha256(text), otherSha, text, inputText,
                ZHANG_GEC134_AR_MODE, candidateSha, text, inputText), "Wrong input SHA accepted as strict");
            require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_SOURCE_MODE, candidateSha,
                zhangCheckpointSha256(text), zhangCheckpointSha256(inputText), text, inputText,
                ZHANG_GEC134_AR_MODE, candidateSha, text, inputText), "FLOAT source accepted as strict AR");
            auto crossRole = own; crossRole["ZHANG_E29_RESTORE_BUNDLE"] = other.root + "/checkpoints/" + captureId;
            reject(fields, crossRole, "cross-role own capture");
            for (const auto& id : {std::string("E29-20240717010130-e124-r32f626de9463"),
                    std::string("E29-20240717010700-e135-r32f626de9463"), std::string("unknown")})
            {
                auto wrong = own; wrong["ZHANG_E29_RESTORE_BUNDLE"] = f.root + "/checkpoints/" + id;
                require(!validateZhangGec134ArContinuationEnvironment(wrong, reason), "Other own-capture endpoint accepted");
            }
            auto wrongSource = env; wrongSource["ZHANG_E29_RESTORE_BUNDLE"] = ZHANG_GEC134_SOURCE_CHECKPOINT_PATH;
            reject(fields, wrongSource, "FLOAT source path");
            require(validate(fields, env, f.binary), "Positive failed after rejection matrix: " + reason);
        }
        std::cout << "GEC134_DUAL_OPT_ONE_EPOCH_POLICY_PASSED checks=" << checks
            << " old_profiles_unchanged=5 exact_roles=2 e122_to_e123=1 no_native_epoch=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GEC134_DUAL_OPT_ONE_EPOCH_POLICY_FAIL " << error.what() << '\n';
        return 1;
    }
}
