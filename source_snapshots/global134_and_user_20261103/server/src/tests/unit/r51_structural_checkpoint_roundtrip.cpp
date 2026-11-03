// Link this DRIVER against the already built production PEA objects. It does
// not include/recompile zhangPppAr.cpp and never seeds a private fixture map.
#include "common/zhangCheckpoint.hpp"
#include "common/algebra.hpp"
#include "common/acsConfig.hpp"
#include "common/receiver.hpp"
#include "pea/zhangPppAr.hpp"
#include "pea/zhangReference.hpp"
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

static void require(bool condition, const std::string& reason)
{ if (!condition) throw std::runtime_error(reason); }

static ZhangCheckpointBundle load(const std::filesystem::path& path, const std::string& binary)
{
    std::ifstream file(path.parent_path() / "checkpoint_manifest.json");
    require(bool(file), "CHECKPOINT_MANIFEST_MISSING");
    std::ostringstream bytes;
    bytes << file.rdbuf();
    const auto document = boost::json::parse(bytes.str());
    const auto& object = document.as_object();
    auto field = [&](const char* name)
    {
        const auto& value = object.at(name).as_string();
        return std::string(value.data(), value.size());
    };
    ZhangCheckpointExpectations expected;
    expected.experimentMode = field("experiment_mode");
    expected.binarySha256 = binary;
    expected.configSha256 = zhangCheckpointSha256(field("config_text"));
    expected.inputManifestSha256 = zhangCheckpointSha256(field("input_manifest_text"));
    expected.platformFingerprint = field("platform");
    expected.compilerFingerprint = field("compiler");
    expected.linearAlgebraFingerprint = field("linear_algebra");
    expected.endianness = field("endianness");
    require(binary == field("binary_sha256") && expected.configSha256 == field("config_sha256") &&
        expected.inputManifestSha256 == field("input_manifest_sha256"), "CHECKPOINT_PROVENANCE_MISMATCH");
    const auto envelope = verifyZhangCheckpointBundleEnvelope(path.string());
    require(envelope.valid, envelope.failureReason);
    ZhangCheckpointBundle bundle;
    const auto decoded = readZhangCheckpointBundle(path.string(), expected, bundle);
    require(decoded.valid, decoded.failureReason);
    return bundle;
}

int main(int argc, char** argv)
{
    require(argc == 3 || (argc == 4 && std::string(argv[3]) == "--upgrade-real-schema5"),
        "Usage: real.bundle producer-binary-sha256 [--upgrade-real-schema5]");
    const bool upgradeRealSchema5 = argc == 4;
    const auto bundle = load(argv[1], argv[2]);
    const auto& section = bundle.sections.at(ZHANG_PPP_AR_CHECKPOINT_SECTION_NAME);
    require(section.schemaVersion == (upgradeRealSchema5 ? 5 : 6),
        "INPUT_SCHEMA_DOES_NOT_MATCH_EXPLICIT_TEST_MODE");
    std::string input = section.payload;
    require(input.size() >= 4, "PPP_AR_PAYLOAD_TOO_SHORT");
    // The schema6 envelope's LAST field is structuralTransportChains. Every
    // live committed/composed transition is serialized with failureReason
    // NONE as its last member. An empty final map terminates in its binary
    // zero count/item-version, not NONE. This narrow wire assertion is checked
    // only after the production importer fully decodes/validates the payload.
    // It proves nonempty storage, NOT a positive number of unpublished links.

    KFState owner;
    ReceiverMap receiverBindings;
    std::string reason;
    require(restoreZhangCheckpointKfCoreWithReceiverResolver(bundle.kfCore, owner,
        [&](const std::string& id) -> Receiver*
        {
            auto& receiver = receiverBindings[id];
            receiver.id = id;
            return &receiver;
        }, &reason), reason);
    require(zhangCheckpointRuntimeId(owner) == bundle.manifest.runtimeId, "CORE_RUNTIME_ID_MISMATCH");
    require(importZhangGraphCheckpointSection(owner, bundle.manifest.runtimeId,
        bundle.sections.at(ZHANG_GRAPH_CHECKPOINT_SECTION_NAME).payload, reason), reason);
    const auto preflight = importZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId,
        input, true, owner.x.size());
    require(preflight.valid, "PRODUCTION_PREFLIGHT_FAILED:" + preflight.failureReason);
    if (upgradeRealSchema5)
    {
        const auto imported = importZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId,
            input, false, owner.x.size());
        require(imported.valid, "REAL_SCHEMA5_IMPORT_FAILED:" + imported.failureReason);
        // These two transient commit flags are not part of the core wire
        // format. Restore them ONLY in this explicitly named checkpoint-barrier
        // test, using an authenticated accepted real POST_EPOCH_COMMITTED core.
        // This invokes no search/conditioning and creates no historical RHS.
        owner.lastPppTransactionCommitted = true;
        owner.lastPppQcConverged = true;
        acsConfig.zhangPppAr.baseline_observables =
            {{E_Sys::GPS, {E_ObsCode::L1C, E_ObsCode::L2W}}};
        ::setenv("ZHANG_R51_ENABLE", "1", 1);
        ::setenv("ZHANG_R51_EPOCH_ONLY_AR", "1", 1);
        ::setenv("ZHANG_R51_PRODUCT_REPAIR", "1", 1);
        const auto beforeMean = owner.x;
        const auto beforeCovariance = owner.P;
        zhangAdvanceCurrentProductStructuralChains(std::cout, owner, owner);
        require((owner.x.array() == beforeMean.array()).all() &&
            (owner.P.array() == beforeCovariance.array()).all(), "STRUCTURAL_ANCHOR_UPGRADE_MUTATED_FLOAT");
        const auto upgraded = exportZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId, input);
        require(upgraded.valid, "SCHEMA6_ANCHOR_EXPORT_FAILED:" + upgraded.failureReason);
        std::cout << "REAL_SCHEMA5_ANCHOR_UPGRADE source_producer=" << argv[2]
            << " source_schema=5 exported_schema=6 source_epoch=" << bundle.manifest.epoch
            << " production_advance_hook=1 manual_map_seed=0 transient_commit_flags_restored=1"
            << " tested_positive_unpublished_link_count=0\n";
    }
    const bool nonemptyChainTailMarker = input.size() >= 4 &&
        input.compare(input.size() - 4, 4, "NONE") == 0;
    require(nonemptyChainTailMarker, "REAL_SCHEMA6_STRUCTURAL_MAP_EMPTY_OR_NONCANONICAL_TAIL");

    const auto originalMean = owner.x;
    const auto originalCovariance = owner.P;
    for (int round = 0; round < 2; ++round)
    {
        const auto imported = importZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId,
            input, false, owner.x.size());
        require(imported.valid, "PRODUCTION_IMPORT_FAILED:" + imported.failureReason);
        std::string output;
        const auto exported = exportZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId, output);
        require(exported.valid, "PRODUCTION_EXPORT_FAILED:" + exported.failureReason);
        if (input != output)
        {
            std::size_t first = 0, suffix = 0;
            while (first < input.size() && first < output.size() && input[first] == output[first]) ++first;
            while (suffix < input.size() - first && suffix < output.size() - first &&
                input[input.size() - 1 - suffix] == output[output.size() - 1 - suffix]) ++suffix;
            std::cerr << "ROUNDTRIP_DIFFERENCE input_bytes=" << input.size()
                << " output_bytes=" << output.size() << " first_offset=" << first
                << " equal_suffix_bytes=" << suffix << '\n';
        }
        require(input == output, "PRODUCTION_PPP_AR_PAYLOAD_ROUNDTRIP_NOT_BYTE_IDENTICAL");
        require((owner.x.array() == originalMean.array()).all() &&
            (owner.P.array() == originalCovariance.array()).all(), "RUNTIME_ROUNDTRIP_MUTATED_FLOAT_X_OR_P");
        std::cout << "STRUCTURAL_CHECKPOINT_ROUNDTRIP round=" << round + 1
            << " bytes=" << output.size() << " sha256=" << zhangCheckpointSha256(output)
            << " entire_ppp_ar_payload_byte_identical=1 float_x_p_identical=1\n";
    }
    const auto truncated = importZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId,
        input.substr(0, input.size() - 1), true, owner.x.size());
    require(!truncated.valid, "TRUNCATED_SCHEMA6_CHAIN_TAIL_ACCEPTED");
    std::string afterInvalidPreflight;
    const auto after = exportZhangPppArCheckpointSection(owner, bundle.manifest.runtimeId, afterInvalidPreflight);
    require(after.valid && afterInvalidPreflight == input, "FAILED_PREFLIGHT_MUTATED_LIVE_STRUCTURAL_CHAIN");
    std::cout << "REAL_STRUCTURAL_CHECKPOINT_PASS schema=6 nonempty_chain_map_tail_marker=1"
        << " actual_production_import_export=1 full_payload_roundtrips=2"
        << " truncated_tail_rejected=1 failed_preflight_side_effect_free=1"
        << " positive_unpublished_link_count_not_claimed=1"
        << " input_mode=" << (upgradeRealSchema5 ? "REAL_SCHEMA5_UPGRADE" : "REAL_SCHEMA6")
        << " epoch=" << bundle.manifest.epoch << '\n';
}
