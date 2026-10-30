// Compare two independently executed real PEA checkpoints, not two
// reconstructions from the same saved epoch. This reader does no estimation.
#include "common/zhangCheckpoint.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <stdexcept>

namespace fs = std::filesystem;

static void require(bool condition, const std::string& why)
{ if (!condition) throw std::runtime_error(why); }

static ZhangCheckpointBundle load(const fs::path& path, const std::string& producer,
                                 const std::string& expectedMode)
{
    std::ifstream input(path.parent_path() / "checkpoint_manifest.json");
    std::ostringstream text; text << input.rdbuf();
    const auto value = boost::json::parse(text.str());
    const auto& doc = value.as_object();
    auto field = [&](const char* name) { return std::string(doc.at(name).as_string().c_str()); };
    ZhangCheckpointExpectations expected;
    expected.experimentMode = expectedMode;
    expected.binarySha256 = producer;
    expected.configSha256 = zhangCheckpointSha256(field("config_text"));
    expected.inputManifestSha256 = zhangCheckpointSha256(field("input_manifest_text"));
    expected.platformFingerprint = field("platform");
    expected.compilerFingerprint = field("compiler");
    expected.linearAlgebraFingerprint = field("linear_algebra");
    expected.endianness = field("endianness");
    require(field("experiment_mode") == expectedMode, "UNEXPECTED_EXPERIMENT_MODE");
    require(producer == field("binary_sha256"), "WRONG_PRODUCER_BINARY");
    require(expected.configSha256 == field("config_sha256") &&
            expected.inputManifestSha256 == field("input_manifest_sha256"), "MANIFEST_PROVENANCE_HASH_FAILED");
    const auto envelope = verifyZhangCheckpointBundleEnvelope(path.string());
    require(envelope.valid, "CHECKPOINT_ENVELOPE:" + envelope.failureReason);
    ZhangCheckpointBundle result;
    const auto decoded = readZhangCheckpointBundle(path.string(), expected, result);
    require(decoded.valid, "CHECKPOINT_DECODE:" + decoded.failureReason);
    return result;
}

template<class T> static std::string wire(const T& value)
{
    std::ostringstream stream(std::ios::binary);
    boost::archive::binary_oarchive archive(stream);
    archive << value;
    return stream.str();
}

static boost::json::object compare(const ZhangCheckpointBundle& left, const ZhangCheckpointBundle& right)
{
    boost::json::object out, fields;
    boost::json::array failures;
    auto check = [&](const std::string& name, bool equal)
    { fields[name] = equal; if (!equal) failures.push_back(boost::json::value(name)); };
    // Both independent processes deliberately use the exact same config,
    // input/output paths and logical runtime ID. No numerical normalization
    // or binary string substitution is allowed/needed.
    const auto& lm = left.manifest; const auto& rm = right.manifest;
    check("manifest.experimentMode", lm.experimentMode == rm.experimentMode);
    check("manifest.runtimeId", lm.runtimeId == rm.runtimeId);
    check("manifest.checkpointId", lm.checkpointId == rm.checkpointId);
    check("manifest.parentCheckpointId", lm.parentCheckpointId == rm.parentCheckpointId);
    check("manifest.epoch", lm.epoch == rm.epoch);
    check("manifest.binarySha256", lm.binarySha256 == rm.binarySha256);
    check("manifest.configText", lm.configText == rm.configText);
    check("manifest.inputManifestText", lm.inputManifestText == rm.inputManifestText);
    check("manifest.platform", lm.platformFingerprint == rm.platformFingerprint);
    check("manifest.compiler", lm.compilerFingerprint == rm.compilerFingerprint);
    check("manifest.linearAlgebra", lm.linearAlgebraFingerprint == rm.linearAlgebraFingerprint);
    check("manifest.endianness", lm.endianness == rm.endianness);
    const auto& a = left.kfCore; const auto& b = right.kfCore;
    auto numeric = [&](const char* name, const auto& first, const auto& second)
    {
        const bool shape = first.rows() == second.rows() && first.cols() == second.cols();
        const bool finite = shape && first.allFinite() && second.allFinite();
        const bool exact = finite && (first.array() == second.array()).all();
        check(std::string("core.") + name, exact);
        boost::json::object metric{{"shape_equal", shape}, {"all_finite", finite}, {"exact_equal", exact},
            {"left_rows", first.rows()}, {"left_columns", first.cols()},
            {"right_rows", second.rows()}, {"right_columns", second.cols()}};
        if (finite)
        {
            const double maximum = first.size() ? (first - second).cwiseAbs().maxCoeff() : 0.;
            metric["maximum_absolute_difference_native"] = maximum;
            metric["relative_frobenius_difference"] = (first - second).norm() / std::max(1., second.norm());
        }
        out[std::string("numeric_") + name] = std::move(metric);
    };
    numeric("x", a.x, b.x); numeric("P", a.P, b.P); numeric("dx", a.dx, b.dx);
    numeric("prefitRatios", a.prefitRatios, b.prefitRatios);
    numeric("postfitRatios", a.postfitRatios, b.postfitRatios);
#define FIELD(name) check("core." #name, wire(a.name) == wire(b.name))
    FIELD(time); FIELD(kfIndexMap); FIELD(stateTransitionMap); FIELD(gaussMarkovTauMap);
    FIELD(gaussMarkovMuMap); FIELD(procNoiseMap); FIELD(initNoiseMap); FIELD(sigmaMaxMap);
    FIELD(outageLimitMap); FIELD(exponentialNoiseMap); FIELD(pseudoStateMap); FIELD(pseudoParentMap);
    FIELD(errorCountMap); FIELD(filterChunkMap); FIELD(metaDataMap); FIELD(lsqRequired);
    FIELD(sigmaPass); FIELD(chiQCPass); FIELD(chi2); FIELD(dof); FIELD(chi2PerDof); FIELD(qc);
    FIELD(id); FIELD(rtsBasename); FIELD(outputResiduals); FIELD(outputMongoMeasurements);
    FIELD(statisticsMap); FIELD(statisticsMapSum);
#undef FIELD
    boost::json::array clockMetadata;
    std::set<std::string> clockSystems;
    for (const auto& [key, value] : a.metaDataMap)
        if (key.starts_with("ZHANG_CLOCK_LIFECYCLE_V2|SPP_OFFSET|"))
        {
            clockMetadata.push_back(boost::json::object{{"key", key}, {"value", value},
                {"equal", b.metaDataMap.contains(key) && b.metaDataMap.at(key) == value}});
            for (const std::string system : {"GPS", "GAL", "BDS"})
                if (key.starts_with("ZHANG_CLOCK_LIFECYCLE_V2|SPP_OFFSET|" + system + "|"))
                    clockSystems.insert(system);
        }
    check("core.authoritative_spp_offset_all_GEC_present", clockSystems == std::set<std::string>{"GPS", "GAL", "BDS"});
    out["authoritative_spp_metadata"] = std::move(clockMetadata);
    std::map<std::string, std::set<int>> positions;
    bool noRate = true, positivePositionVariance = true;
    for (const auto& [key, index] : a.kfIndexMap)
    {
        if (key.type == KF::REC_POS)
        {
            positions[key.receiver].insert(key.number);
            positivePositionVariance = positivePositionVariance && index >= 0 && index < a.P.rows() &&
                std::isfinite(a.P(index, index)) && a.P(index, index) > 0;
        }
        if (key.type == KF::REC_POS_RATE) noRate = false;
    }
    boost::json::object coordinates;
    bool xyz = true;
    for (const auto& [receiver, components] : positions)
    {
        xyz = xyz && components == std::set<int>{0, 1, 2};
        boost::json::array numbers;
        for (const int number : components) numbers.push_back(number);
        coordinates[receiver] = std::move(numbers);
    }
    out["coordinate_components_by_receiver"] = std::move(coordinates);
    check("core.static_positions_all_XYZ_positive_variance", xyz && positivePositionVariance && positions.size() == 12);
    check("core.no_receiver_position_rate", noRate);
    boost::json::array sections;
    check("sections.inventory", left.sections.size() == right.sections.size());
    for (const auto& [name, first] : left.sections)
    {
        const auto item = right.sections.find(name);
        const bool present = item != right.sections.end();
        const bool equal = present && first.schemaVersion == item->second.schemaVersion &&
                           first.payload == item->second.payload && first.sha256 == item->second.sha256;
        check("section." + name, equal);
        sections.push_back(boost::json::object{{"name", name}, {"schema", first.schemaVersion},
            {"left_bytes", first.payload.size()}, {"right_present", present},
            {"entire_payload_byte_identical", equal}, {"left_sha256", first.sha256},
            {"right_sha256", present ? item->second.sha256 : "MISSING"}});
    }
    out["sections"] = std::move(sections);
    out["checks"] = std::move(fields);
    out["failed_fields"] = std::move(failures);
    out["passed"] = out.at("failed_fields").as_array().empty();
    out["checkpoint_epoch"] = lm.epoch;
    out["state_dimension"] = a.x.size();
    out["covariance_entries"] = a.P.size();
    out["identity_policy"] = "Same original runtime/config/input/output paths; only manifest.createdUtc is not compared";
    out["numeric_policy"] = "Strict entrywise exact equality; nonzero differences are reported, never rounded or ignored";
    out["scope"] = "Independent continuous process versus actual fresh-process checkpoint continuation; no numerical reconstruction oracle";
    return out;
}

int main(int argc, char** argv)
{
    require(argc == 6, "Usage: continuous.bundle resumed.bundle producer-sha256 expected-mode report.json");
    const auto left = load(argv[1], argv[3], argv[4]);
    const auto right = load(argv[2], argv[3], argv[4]);
    const auto report = compare(left, right);
    std::ofstream output(argv[5]); output << boost::json::serialize(report) << '\n';
    require(bool(output), "REPORT_WRITE_FAILED");
    std::cout << "REAL_MULTI_GNSS_CHECKPOINT_COMPARE passed=" << report.at("passed")
              << " dimension=" << left.kfCore.x.size() << " sections=" << left.sections.size()
              << " entire_P_entries=" << left.kfCore.P.size() << '\n';
    return report.at("passed").as_bool() ? 0 : 1;
}
