// Read-only, authenticated CORE metadata inventory. No solver/restore mutation.
#include "common/zhangCheckpoint.hpp"
#include <boost/json.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

static void require(bool pass, const std::string& why)
{ if (!pass) throw std::runtime_error(why); }

int main(int argc, char** argv)
try {
    require(argc == 5, "Usage: bundle producer-sha256 expected-mode output.json");
    const std::filesystem::path path(argv[1]);
    std::ifstream manifest(path.parent_path() / "checkpoint_manifest.json");
    std::ostringstream text; text << manifest.rdbuf();
    const auto parsed = boost::json::parse(text.str());
    const auto& doc = parsed.as_object();
    auto field = [&](const char* key) { return std::string(doc.at(key).as_string()); };
    ZhangCheckpointExpectations expected;
    expected.experimentMode = argv[3]; expected.binarySha256 = argv[2];
    expected.configSha256 = zhangCheckpointSha256(field("config_text"));
    expected.inputManifestSha256 = zhangCheckpointSha256(field("input_manifest_text"));
    expected.platformFingerprint = field("platform");
    expected.compilerFingerprint = field("compiler");
    expected.linearAlgebraFingerprint = field("linear_algebra");
    expected.endianness = field("endianness");
    require(field("experiment_mode") == argv[3] && field("binary_sha256") == argv[2],
        "WRONG_CHECKPOINT_PRODUCER_OR_MODE");
    require(expected.configSha256 == field("config_sha256") &&
        expected.inputManifestSha256 == field("input_manifest_sha256"), "MANIFEST_HASH_MISMATCH");
    const auto envelope = verifyZhangCheckpointBundleEnvelope(path.string());
    require(envelope.valid, envelope.failureReason);
    ZhangCheckpointBundle bundle;
    const auto read = readZhangCheckpointBundle(path.string(), expected, bundle);
    require(read.valid, read.failureReason);
    const auto epoch = restoreZhangCheckpointTime(bundle.kfCore.time).to_string(9);
    boost::json::object metadata, offsetsBySystem, seedsBySystem;
    std::set<std::string> recClocks, savedOffsets;
    for (const auto& [key, index] : bundle.kfCore.kfIndexMap)
        if (key.type == KF::REC_CLOCK)
            recClocks.insert(enum_to_string(key.satellite.sys) + "|" + key.receiver);
    bool valid = true;
    unsigned stale = 0, invalid = 0;
    for (const auto& [key, bytes] : bundle.kfCore.metaDataMap)
    {
        const std::string prefix = "ZHANG_CLOCK_LIFECYCLE_V2|";
        if (!key.starts_with(prefix)) continue;
        boost::system::error_code error;
        const auto value = boost::json::parse(bytes, error);
        if (error || !value.is_object()) { ++invalid; valid = false; continue; }
        const auto& obj = value.as_object();
        metadata[key] = value;
        const auto* savedEpoch = obj.if_contains("epoch");
        if (!savedEpoch || !savedEpoch->is_string() || savedEpoch->as_string() != epoch)
        { ++stale; valid = false; }
        const auto offsetPrefix = prefix + "SPP_OFFSET|";
        const auto seedPrefix = prefix + "SEED|";
        const auto numberName = key.starts_with(offsetPrefix) ? "offset_m" : "translation_m";
        const auto* number = obj.if_contains(numberName);
        if (!number || !number->is_number() || !std::isfinite(number->to_number<double>()))
        { ++invalid; valid = false; }
        const auto parts = key.substr(key.starts_with(offsetPrefix) ? offsetPrefix.size() : seedPrefix.size());
        const auto split = parts.find('|');
        const auto system = parts.substr(0, split);
        auto& counts = key.starts_with(offsetPrefix) ? offsetsBySystem : seedsBySystem;
        const auto* count = counts.if_contains(system);
        counts[system] = count ? count->as_int64() + 1 : 1;
        if (key.starts_with(offsetPrefix)) {
            savedOffsets.insert(parts);
            const auto* chart = obj.if_contains("chart");
            if (!chart || !chart->is_string() || chart->as_string().empty()) { ++invalid; valid = false; }
        }
    }
    boost::json::array missing;
    for (const auto& key : recClocks) if (!savedOffsets.contains(key)) missing.push_back(boost::json::value(key));
    boost::json::object report{{"schema", "CHECKPOINT_CLOCK_METADATA_INVENTORY_V1"},
        {"bundle_sha256", zhangCheckpointFileSha256(path.string())}, {"producer_sha256", argv[2]},
        {"experiment_mode", argv[3]}, {"checkpoint_epoch", bundle.manifest.epoch}, {"core_epoch", epoch},
        {"v2_metadata", metadata}, {"offset_count_by_system", offsetsBySystem},
        {"seed_count_by_system", seedsBySystem}, {"receiver_clock_count", recClocks.size()},
        {"missing_v2_offset_receiver_clocks", missing}, {"stale_epoch_count", stale},
        {"invalid_metadata_count", invalid}, {"all_recorded_v2_values_finite_current_epoch", valid},
        {"scope", "READ_ONLY_ACTUAL_SHORT_REPLAY_CORE_NOT_180_STATION_OR_SCIENTIFIC_ACCURACY"}};
    std::ofstream output(argv[4]); output << boost::json::serialize(report) << '\n';
    require(bool(output), "REPORT_WRITE_FAILED");
    std::cout << "CHECKPOINT_CLOCK_METADATA_INVENTORY records=" << metadata.size()
              << " stale=" << stale << " invalid=" << invalid << " missing=" << missing.size()
              << " no_estimation=1\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << "CHECKPOINT_CLOCK_METADATA_AUDIT_FAILED " << e.what() << '\n'; return 1;
}
