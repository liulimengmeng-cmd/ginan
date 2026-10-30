// Independent small policy test. Reads only the real e122 JSON manifest, never
// checkpoint.bundle; no PEA main/KF/importer is linked or executed. Fixtures are
// constructed by an independent length-prefix codec, not the production parser.
#include "common/zhangGec134ArContinuation.hpp"
#include <openssl/evp.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

std::string zhangCheckpointSha256(const std::string& input)
{
    unsigned char output[EVP_MAX_MD_SIZE]; unsigned int size = 0;
    if (!EVP_Digest(input.data(), input.size(), output, &size, EVP_sha256(), nullptr) || size != 32)
        throw std::runtime_error("Independent standard SHA-256 failed");
    std::ostringstream text;
    for (unsigned int i = 0; i < size; ++i)
        text << std::hex << std::setw(2) << std::setfill('0') << unsigned(output[i]);
    return text.str();
}

namespace
{
std::size_t checks = 0;
void require(bool value, const std::string& reason)
{
    ++checks;
    if (!value) throw std::runtime_error(reason);
}
bool scalar(const std::string& key)
{
    return key == "command_argument_count" || key == "config_file_count" ||
        key == "effective_capture_epoch_count" ||
        (key.rfind("config_file_", 0) == 0 && key.ends_with("_size"));
}
std::map<std::string, std::string> decode(const std::string& text)
{
    std::istringstream input(text);
    std::string line;
    require(bool(std::getline(input, line)) && line == "E29_CANONICAL_CONFIG_V1", "fixture header");
    std::map<std::string, std::string> fields;
    while (std::getline(input, line))
    {
        const auto equal = line.find('=');
        require(equal != std::string::npos, "fixture record");
        const auto key = line.substr(0, equal), length = line.substr(equal + 1);
        std::size_t consumed = 0;
        const auto size = std::stoull(length, &consumed);
        require(consumed == length.size() && size < 16 * 1024 * 1024, "fixture length");
        std::string value = length;
        if (!scalar(key))
        {
            value.resize(size); input.read(value.data(), size);
            require(std::size_t(input.gcount()) == size && input.get() == '\n', "fixture bytes");
        }
        require(fields.emplace(key, value).second, "fixture duplicate");
    }
    return fields;
}
std::string encode(const std::map<std::string, std::string>& fields)
{
    std::ostringstream output; output << "E29_CANONICAL_CONFIG_V1\n";
    for (const auto& [key, value] : fields)
    {
        output << key << '=' << (scalar(key) ? value : std::to_string(value.size())) << '\n';
        if (!scalar(key)) output << value << '\n';
    }
    return output.str();
}
std::string overlay()
{
    const std::string root = ZHANG_GEC134_AR_CONTINUATION_ROOT;
    return "processing_options:\n  epoch_control:\n    end_epoch: 2024-07-17 02:00:30\n"
        "    max_epochs: 242\n  gnss_general:\n    zhang_pppar:\n"
        "      checkpoint_capture_epochs: [2024-07-17 01:30:30, 2024-07-17 02:00:30]\n"
        "      checkpoint_output_directory: " + root + "/checkpoints\n"
        "      checkpoint_restore_path: " + ZHANG_GEC134_AR_E122_PATH + "\n"
        "      product_filename: " + root + "/outputs/zhang_internal_products.csv\n"
        "      product_covariance_filename: " + root + "/outputs/zhang_internal_product_covariance.csv\n"
        "outputs:\n  outputs_root: " + root + "/outputs\n"
        "  metadata:\n    config_description: independent AR continuation contract\n";
}
void setOverlay(std::map<std::string, std::string>& target, const std::string& bytes)
{
    target["config_file_11_bytes"] = bytes;
    target["config_file_11_size"] = std::to_string(bytes.size());
    target["config_file_11_sha256"] = zhangCheckpointSha256(bytes);
}
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected real e122 checkpoint_manifest.json path");
        const auto manifest = YAML::LoadFile(argv[1]);
        const std::string sourceText = manifest["config_text"].as<std::string>();
        require(zhangCheckpointSha256(sourceText) == ZHANG_GEC134_AR_E122_CONFIG_SHA256,
            "Fixture is not the pinned e122 source canonical");
        require(manifest["experiment_mode"].as<std::string>() == ZHANG_GEC134_AR_MODE &&
            manifest["checkpoint_id"].as<std::string>() == ZHANG_GEC134_AR_E122_ID,
            "Fixture is not an AR e122 archive");
        const auto source = decode(sourceText);
        auto target = source;
        target[std::string("environment_") + ZHANG_GEC134_AR_CONTINUATION_FLAG] = "1";
        target["gec134_ar_continuation_policy"] = ZHANG_GEC134_AR_CONTINUATION_POLICY;
        target["config_file_count"] = "12";
        target["config_file_11_path"] = ZHANG_GEC134_AR_CONTINUATION_OVERLAY;
        setOverlay(target, overlay());
        target["command_argument_count"] = "151";
        for (int i = 150; i >= 12; --i)
            target["command_argument_" + std::to_string(i + 1)] = source.at("command_argument_" + std::to_string(i));
        target["command_argument_12"] = ZHANG_GEC134_AR_CONTINUATION_OVERLAY;
        const std::string root = ZHANG_GEC134_AR_CONTINUATION_ROOT;
        target["effective_end_epoch"] = "2024-07-17T02:00:30";
        target["effective_max_epochs"] = "242";
        target["effective_outputs_root"] = root + "/outputs";
        target["effective_product_filename"] = root + "/outputs/zhang_internal_products.csv";
        target["effective_product_covariance_filename"] = root + "/outputs/zhang_internal_product_covariance.csv";
        target["effective_checkpoint_output_directory"] = root + "/checkpoints";
        target["effective_checkpoint_restore_path"] = ZHANG_GEC134_AR_E122_PATH;
        target["effective_capture_epoch_count"] = "2";
        target["effective_capture_epoch_0"] = "2024-07-17 01:30:30";
        target["effective_capture_epoch_1"] = "2024-07-17 02:00:30";
        std::map<std::string, std::string> environment = {
            {"ZHANG_GEC134_FLOAT_TO_AR_MIGRATION","1"}, {"ZHANG_MULTI_GNSS_AR_CHECKPOINT","1"},
            {"ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT","UNSET"}, {"ZHANG_R51_REUSE_FLOAT_20240717","UNSET"},
            {"ZHANG_R51_REUSE_RATIO_20240717","UNSET"}, {"ZHANG_R51_ENABLE","1"},
            {"ZHANG_R51_RATIO_ONLY","1"}, {"ZHANG_R51_RATIO_THRESHOLD","2"},
            {"ZHANG_R51_EPOCH_ONLY_AR","1"}, {"ZHANG_R51_AR_START_GPST_SECONDS",ZHANG_GEC134_AR_START},
            {"ZHANG_R51_FLOAT_THREADS","12"}, {"OMP_NUM_THREADS","12"}, {"OPENBLAS_NUM_THREADS","1"},
            {"MKL_NUM_THREADS","1"}, {"ZHANG_FLOAT_MATRIX_THREADS","6"},
            {"ZHANG_FLOAT_MEMORY_MODE","parallel"}, {"ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC","off"},
            {"ZHANG_FLOAT_PROFILE","1"}, {"ZHANG_P0_RESOURCE_PROBE","1"}, {"ZHANG_R49_FUSION_WEIGHT","0.10"},
            {"ZHANG_R51_SINGLE_BASIS_HANDOFF","1"}, {"ZHANG_R51_PRODUCT_REPAIR","1"},
            {"ZHANG_R51_INTERN_HANDOFF_PROOFS","1"}, {ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG,"1"},
            {ZHANG_GEC134_AR_CONTINUATION_FLAG,"1"}, {"ZHANG_E29_RESTORE_BUNDLE",ZHANG_GEC134_AR_E122_PATH}};
        const std::string binary = "/home/rx/GINAN/local_real_epoch_memory_20261024_v5/candidate/pea";
        const std::string binarySha(64, 'a');
        const auto acceptedText = encode(target);
        ZhangGec134ArContinuationAuthority authority; std::string reason;
        auto validate = [&](const auto& fields, const auto& env,
                            const std::string& bundle = ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
                            const std::string& input = ZHANG_GEC134_SOURCE_INPUT_SHA256,
                            const std::string& binaryPath = std::string(),
                            const std::string& candidateSha = std::string())
        {
            return validateZhangGec134ArContinuation(bundle, ZHANG_GEC134_AR_E122_BINARY_SHA256,
                ZHANG_GEC134_AR_E122_CONFIG_SHA256, ZHANG_GEC134_SOURCE_INPUT_SHA256,
                sourceText, encode(fields), input, binaryPath.empty() ? binary : binaryPath,
                candidateSha.empty() ? binarySha : candidateSha, env, authority, reason);
        };
        require(validate(target, environment), "Positive rejected: " + reason);
        require(authority.valid() && authority.targetConfigText() == acceptedText, "Token target lost");
        require(authority.matchesSource(ZHANG_GEC134_AR_MODE, ZHANG_GEC134_AR_E122_BINARY_SHA256,
            ZHANG_GEC134_AR_E122_CONFIG_SHA256, ZHANG_GEC134_SOURCE_INPUT_SHA256), "Source lost");
        require(!authority.matchesSource(ZHANG_GEC134_SOURCE_MODE, ZHANG_GEC134_SOURCE_BINARY_SHA256,
            ZHANG_GEC134_SOURCE_CONFIG_SHA256, ZHANG_GEC134_SOURCE_INPUT_SHA256), "FLOAT source accepted");
        require(authority.matchesTarget(ZHANG_GEC134_AR_MODE, binary, binarySha,
            acceptedText, ZHANG_GEC134_SOURCE_INPUT_SHA256), "Candidate binding lost");
        require(!authority.matchesTarget(ZHANG_GEC134_AR_MODE, binary, std::string(64,'b'),
            acceptedText, ZHANG_GEC134_SOURCE_INPUT_SHA256), "Other candidate accepted");
        auto reject = [&](const auto& fields, const auto& env, const std::string& label)
        {
            require(validate(target, environment), "Warm token failed");
            require(!validate(fields, env), "Invalid fixture accepted: " + label);
            require(!authority.valid() && authority.targetConfigText().empty() && !reason.empty(),
                "Rejected request left authority: " + label);
        };
        for (const auto& [key, value] : target)
        {
            auto changed = target; changed[key] = value + "!";
            reject(changed, environment, "changed canonical " + key);
            changed = target; changed.erase(key);
            reject(changed, environment, "missing canonical " + key);
        }
        for (const auto& [key, value] : environment)
        {
            auto changed = environment; changed[key] = value + "!";
            reject(target, changed, "changed environment " + key);
            changed = environment; changed.erase(key);
            reject(target, changed, "missing environment " + key);
        }
        for (const auto& injection : {"extra: forbidden\n", "processing_options:\n  epoch_control:\n    start_epoch: 2024-07-17 01:00:30\n"})
        {
            auto changed = target; setOverlay(changed, overlay() + injection);
            reject(changed, environment, "forbidden overlay");
        }
        for (const auto& [from, to] : std::map<std::string,std::string>{
                {"max_epochs: 242", "max_epochs: 243"}, {"02:00:30", "02:00:00"},
                {"01:30:30", "01:30:00"}, {"config_description: independent", "ratio: 3\n    config_description: independent"}})
        {
            auto bytes = overlay(); const auto at = bytes.find(from);
            require(at != std::string::npos, "Missing mutation fixture");
            bytes.replace(at, from.size(), to); auto changed = target; setOverlay(changed, bytes);
            reject(changed, environment, "overlay leaf " + from);
        }
        require(validate(target, environment), "Warm before SHA rejection");
        require(!validate(target, environment, std::string(64,'0')), "Wrong raw source bundle accepted");
        require(!authority.valid(), "Invalid raw bundle left authority");
        require(!validate(target, environment, ZHANG_GEC134_AR_E122_BUNDLE_SHA256, std::string(64,'0')),
            "Changed scientific inputs accepted");
        require(!validate(target, environment, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
            ZHANG_GEC134_SOURCE_INPUT_SHA256, ZHANG_GEC134_AR_E122_BINARY_PATH), "Old producer accepted as new continuation");
        require(!validate(target, environment, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
            ZHANG_GEC134_SOURCE_INPUT_SHA256, "relative/pea"), "Relative candidate accepted");
        require(!validate(target, environment, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
            ZHANG_GEC134_SOURCE_INPUT_SHA256, binary, std::string(64,'X')), "Invalid candidate SHA accepted");
        // The new profile's own checkpoints must use ordinary strict restore,
        // not acquire an e122 cross-provenance capability. No bundle is read
        // by this routing contract; native full-section validation is separate.
        const std::string ownInputText = manifest["input_manifest_text"].as<std::string>();
        for (const auto& id : {"E29-20240717013030-e182-r32f626de9463",
                               "E29-20240717020030-e242-r32f626de9463"})
        {
            auto ownEnvironment = environment;
            ownEnvironment["ZHANG_E29_RESTORE_BUNDLE"] = root + "/checkpoints/" + id;
            require(validateZhangGec134ArContinuationEnvironment(ownEnvironment, reason),
                "Own strict restore environment rejected: " + reason);
            require(zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, binarySha,
                zhangCheckpointSha256(acceptedText), zhangCheckpointSha256(ownInputText),
                acceptedText, ownInputText, ZHANG_GEC134_AR_MODE, binarySha, acceptedText, ownInputText),
                "Own exact checkpoint cannot take ordinary strict path");
            reject(target, ownEnvironment, "own checkpoint must not get cross-provenance token");
        }
        require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE,
            ZHANG_GEC134_AR_E122_BINARY_SHA256, ZHANG_GEC134_AR_E122_CONFIG_SHA256,
            ZHANG_GEC134_SOURCE_INPUT_SHA256, sourceText, ownInputText, ZHANG_GEC134_AR_MODE,
            binarySha, acceptedText, ownInputText), "Old source routed as new strict source");
        require(!zhangGec134ArContinuationUsesOrdinaryRestore(ZHANG_GEC134_AR_MODE, binarySha,
            zhangCheckpointSha256(acceptedText), zhangCheckpointSha256(ownInputText),
            acceptedText, ownInputText, ZHANG_GEC134_AR_MODE, binarySha, acceptedText + "!", ownInputText),
            "Changed config routed as ordinary strict source");
        auto unknownEnvironment = environment;
        unknownEnvironment["ZHANG_E29_RESTORE_BUNDLE"] = root + "/checkpoints/unknown";
        require(!validateZhangGec134ArContinuationEnvironment(unknownEnvironment, reason),
            "Unbounded own restore source accepted");
        // Independent literal v6 fixtures: do not derive expected roots or
        // endpoint/capture counts from the production profile table.
        struct V6Fixture
        {
            std::string root, overlayPath, end, canonicalEnd, counter;
            std::vector<std::string> captures;
        };
        const std::string v6Binary = "/home/rx/GINAN/local_real_epoch_memory_20261025_v6/candidate/pea";
        const std::vector<V6Fixture> v6Fixtures = {
            {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01",
             "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01/config/100_gec134_ar_e122_v6_pilot_overlay.yaml",
             "2024-07-17 01:07:00", "2024-07-17T01:07:00", "135", {"2024-07-17 01:07:00"}},
            {"/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_hour_20261025_01",
             "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_hour_20261025_01/config/100_gec134_ar_e122_v6_hour_overlay.yaml",
             "2024-07-17 02:00:30", "2024-07-17T02:00:30", "242", {"2024-07-17 01:30:30", "2024-07-17 02:00:30"}}};
        for (const auto& fixture : v6Fixtures)
        {
            auto v6 = target;
            v6["config_file_11_path"] = fixture.overlayPath;
            v6["command_argument_12"] = fixture.overlayPath;
            v6["effective_end_epoch"] = fixture.canonicalEnd;
            v6["effective_max_epochs"] = fixture.counter;
            v6["effective_outputs_root"] = fixture.root + "/outputs";
            v6["effective_product_filename"] = fixture.root + "/outputs/zhang_internal_products.csv";
            v6["effective_product_covariance_filename"] = fixture.root + "/outputs/zhang_internal_product_covariance.csv";
            v6["effective_checkpoint_output_directory"] = fixture.root + "/checkpoints";
            v6["effective_capture_epoch_count"] = std::to_string(fixture.captures.size());
            v6.erase("effective_capture_epoch_1");
            std::string captureText;
            for (std::size_t i = 0; i < fixture.captures.size(); ++i)
            {
                v6["effective_capture_epoch_" + std::to_string(i)] = fixture.captures[i];
                captureText += "        - \"" + fixture.captures[i] + "\"\n";
            }
            const std::string bytes = "processing_options:\n  epoch_control:\n    end_epoch: \"" + fixture.end +
                "\"\n    max_epochs: " + fixture.counter + "\n  gnss_general:\n    zhang_pppar:\n      checkpoint_capture_epochs:\n" + captureText +
                "      checkpoint_output_directory: " + fixture.root + "/checkpoints\n      checkpoint_restore_path: " + ZHANG_GEC134_AR_E122_PATH +
                "\n      product_filename: " + fixture.root + "/outputs/zhang_internal_products.csv\n      product_covariance_filename: " + fixture.root +
                "/outputs/zhang_internal_product_covariance.csv\noutputs:\n  outputs_root: " + fixture.root +
                "/outputs\n  metadata:\n    config_description: independent fixed v6 fixture\n";
            setOverlay(v6, bytes);
            auto v6Validate = [&](const auto& fields, const auto& env, const std::string& candidate = std::string())
            {
                return validate(fields, env, ZHANG_GEC134_AR_E122_BUNDLE_SHA256,
                    ZHANG_GEC134_SOURCE_INPUT_SHA256, candidate.empty() ? v6Binary : candidate);
            };
            require(v6Validate(v6, environment), "V6 bounded positive rejected: " + reason);
            require(authority.matchesTarget(ZHANG_GEC134_AR_MODE, v6Binary, binarySha, encode(v6),
                ZHANG_GEC134_SOURCE_INPUT_SHA256), "V6 candidate/target token binding lost");
            for (const auto key : {"effective_end_epoch", "effective_max_epochs", "effective_outputs_root",
                    "effective_capture_epoch_count", "effective_capture_epoch_0", "config_file_11_path", "command_argument_12"})
            {
                auto changed = v6; changed[key] += "!";
                require(!v6Validate(changed, environment) && !authority.valid(), "V6 changed bound accepted: " + std::string(key));
            }
            auto extra = v6; extra["effective_capture_epoch_2"] = "2024-07-17 01:07:30";
            require(!v6Validate(extra, environment), "V6 extra capture accepted");
            require(!v6Validate(v6, environment, binary), "V6 old candidate path accepted");
            require(!v6Validate(v6, environment, "/tmp/pea"), "V6 unpinned candidate path accepted");
            auto wrongEnvironment = environment; wrongEnvironment["ZHANG_R51_RATIO_THRESHOLD"] = "3";
            require(!v6Validate(v6, wrongEnvironment), "V6 changed ratio accepted");
            wrongEnvironment = environment; wrongEnvironment["ZHANG_R51_EPOCH_ONLY_AR"] = "0";
            require(!v6Validate(v6, wrongEnvironment), "V6 historical integer mode accepted");
            wrongEnvironment = environment; wrongEnvironment["ZHANG_E29_RESTORE_BUNDLE"] = fixture.root + "/checkpoints/unknown";
            require(!v6Validate(v6, wrongEnvironment), "V6 unknown checkpoint accepted");
            wrongEnvironment = environment;
            wrongEnvironment["ZHANG_E29_RESTORE_BUNDLE"] = fixture.root + "/checkpoints/" +
                (fixture.counter == "135" ? "E29-20240717010700-e135-r32f626de9463" : "E29-20240717020030-e242-r32f626de9463");
            require(validateZhangGec134ArContinuationEnvironment(wrongEnvironment, reason), "V6 own strict restart path rejected");
            require(!v6Validate(v6, wrongEnvironment), "V6 own checkpoint acquired source-e122 authority");
        }
        std::cout << "GEC134_AR_CONTINUATION_POLICY_PASSED checks=" << checks
            << " actual_source_e122=1 independent_sha256=1 no_PEA_or_KF=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GEC134_AR_CONTINUATION_POLICY_FAIL " << error.what() << '\n';
        return 1;
    }
}
