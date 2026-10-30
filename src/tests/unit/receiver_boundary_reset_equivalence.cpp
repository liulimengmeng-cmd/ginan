// Test-only independent wire mutations, calling the production read-only
// boundary comparator. Never imports a receiver or executes a PEA epoch.
#include "multi_gnss_receiver_runtime_audit_schema.hpp"
#include <boost/json.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <functional>
#include <sstream>
#include <stdexcept>

using namespace ReceiverAuditSchema;
namespace json = boost::json;

static void require(bool value, const std::string& message)
{ if (!value) throw std::runtime_error(message); }

template<class T> static std::string wire(const T& value)
{
    std::ostringstream stream(std::ios::binary | std::ios::out);
    { boost::archive::binary_oarchive archive(stream, boost::archive::no_header); archive << value; }
    return stream.str();
}

static ZhangReceiverRuntimeEnvelope decode(const std::string& payload)
{
    ZhangReceiverRuntimeEnvelope result;
    std::istringstream stream(payload, std::ios::binary | std::ios::in);
    { boost::archive::binary_iarchive archive(stream, boost::archive::no_header); archive >> result; }
    require(stream.peek() == std::char_traits<char>::eof(), "MIRROR_TRAILING_BYTES");
    require(wire(result) == payload, "MIRROR_DECODE_REENCODE_NOT_EXACT");
    return result;
}

static std::string read(const char* path)
{
    std::ifstream stream(path, std::ios::binary);
    require(bool(stream), "RAW_PAYLOAD_OPEN_FAILED");
    std::string result((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
    require(!stream.bad() && !result.empty(), "RAW_PAYLOAD_READ_FAILED");
    return result;
}

int main(int argc, char** argv)
try
{
    require(argc == 2 || argc == 3, "Usage: SOURCE_RECEIVER.payload [ACTUAL_AFTER_RECEIVER.payload]");
    const auto sourceBytes = read(argv[1]);
    const auto source = decode(sourceBytes);
    require(source.schemaVersion == 3 && !source.receivers.empty(), "NEED_REAL_SCHEMA3_RECEIVER_SOURCE");
    auto after = source;
    std::size_t ready = 0, worksets = 0, observations = 0;
    std::vector<std::string> readyIds, worksetIds;
    for (auto& receiver : after.receivers)
    {
        if (receiver.ready) { ++ready; readyIds.push_back(receiver.mapKey); }
        if (receiver.resetObservationCount)
        { ++worksets; observations += receiver.resetObservationCount; worksetIds.push_back(receiver.mapKey); }
        receiver.ready = false;
        receiver.resetObservationCount = 0;
    }
    const auto afterBytes = wire(after);
    auto accepted = auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, sourceBytes, afterBytes);
    require(accepted.valid && accepted.failureReason == "NONE", "INDEPENDENT_NORMALIZATION_REJECTED:" + accepted.failureReason);
    require(accepted.receiverCount == source.receivers.size() && accepted.readyResetCount == ready &&
        accepted.observationWorksetResetReceiverCount == worksets && accepted.discardedObservationWorksetCount == observations &&
        accepted.readyResetReceiverIds == readyIds && accepted.observationWorksetResetReceiverIds == worksetIds,
        "AUDIT_RESET_COUNTS_OR_IDENTITIES_WRONG");
    require(accepted.sourcePayloadSha256 == zhangCheckpointSha256(sourceBytes) &&
        accepted.afterPayloadSha256 == zhangCheckpointSha256(afterBytes) &&
        accepted.normalizedSourcePayloadSha256 == accepted.afterPayloadSha256 &&
        accepted.rawPayloadsExact == (sourceBytes == afterBytes), "AUDIT_RAW_OR_NORMALIZED_HASH_WRONG");
    const auto exact = auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, afterBytes, afterBytes);
    require(exact.valid && exact.rawPayloadsExact && exact.readyResetCount == 0 &&
        exact.observationWorksetResetReceiverCount == 0, "ALREADY_RESET_EXACT_PAYLOAD_REJECTED");
    json::array negatives;
    auto reject = [&](const char* name, const std::function<void(ZhangReceiverRuntimeEnvelope&)>& mutate)
    {
        auto changed = after;
        mutate(changed);
        const auto result = auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, sourceBytes, wire(changed));
        require(!result.valid && !result.failureReason.empty() && result.failureReason != "NONE",
            std::string("MUTATION_WAS_ACCEPTED:") + name);
        negatives.emplace_back(name);
    };
    reject("apriori_coordinate", [](auto& x){ x.receivers.front().aprioriPosition(0) += 0.125; });
    reject("SNX_coordinate", [](auto& x){ x.receivers.front().sinex.position(1) += 0.125; });
    reject("SPP_full_covariance", [](auto& x){
        auto& p = x.receivers.front().solution.sppState.P;
        require(p.rows() > 1 && p.cols() > 1, "NEED_REAL_SPP_P"); p(1, 1) += 0.125;
    });
    reject("SPP_state_order", [](auto& x){
        auto& keys = x.receivers.front().solution.sppState.kfIndexMap;
        require(keys.size() > 1, "NEED_REAL_SPP_KEYS"); auto a = keys.begin(), b = std::next(a); std::swap(a->second, b->second);
    });
    reject("cumulative_observation_count", [](auto& x){ ++x.receivers.front().observationCount; });
    reject("completed_receiver_epoch_count", [](auto& x){ ++x.receivers.front().epochCount; });
    reject("last_epoch", [](auto& x){
        auto& value = x.receivers.front().lastEpoch;
        value = captureZhangCheckpointTime(restoreZhangCheckpointTime(value) + 30.0);
    });
    reject("cumulative_slip_count", [](auto& x){ ++x.receivers.front().slipCount; });
    reject("satellite_GF_history", [](auto& x){
        auto& sats = x.receivers.front().satelliteStatuses;
        require(!sats.empty(), "NEED_REAL_SATELLITE_HISTORY"); sats.begin()->second.geometryFree += 0.125;
    });
    reject("signal_QC_history", [](auto& x){
        for (auto& [sat, status] : x.receivers.front().satelliteStatuses)
            if (!status.signalStatuses.empty()) { ++status.signalStatuses.begin()->second.phaseRejectCount; return; }
        throw std::runtime_error("NEED_REAL_SIGNAL_HISTORY");
    });
    reject("retained_ready_true", [](auto& x){ x.receivers.front().ready = true; });
    reject("retained_observation_workset", [](auto& x){ x.receivers.front().resetObservationCount = 1; });
    reject("unknown_ready_policy", [](auto& x){ x.readyPolicy = "KEEP_READY"; });
    reject("unknown_observation_policy", [](auto& x){ x.observationPolicy = "KEEP_WORKSET"; });
    reject("wrong_boundary", [](auto& x){ x.boundary = "MID_EPOCH"; });
    reject("old_schema", [](auto& x){ x.schemaVersion = 2; });
    reject("runtime_owner", [](auto& x){ x.runtimeId += "-foreign"; });
    reject("dynamic_alias_inventory", [](auto& x){ x.customAliasesMap["MAL2"].insert("INJECTED_ALIAS"); });
    reject("trace_path", [](auto& x){ x.receivers.front().traceFilename += ".changed"; });
    auto unknownSource = source; unknownSource.readyPolicy = "UNKNOWN";
    require(!auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, wire(unknownSource), afterBytes).valid,
        "UNKNOWN_SOURCE_POLICY_ACCEPTED"); negatives.emplace_back("unknown_source_policy");
    require(!auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, sourceBytes, afterBytes + "x").valid,
        "TRAILING_BYTES_ACCEPTED"); negatives.emplace_back("trailing_bytes");
    require(!auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, sourceBytes.substr(0, sourceBytes.size()-1), afterBytes).valid,
        "TRUNCATED_SOURCE_ACCEPTED"); negatives.emplace_back("truncated_source");
    require(wire(source) == sourceBytes && wire(after) == afterBytes, "IMMUTABLE_TEST_INPUTS_CHANGED");
    bool actualMatched = false;
    if (argc == 3)
    {
        const auto actualBytes = read(argv[2]); decode(actualBytes);
        const auto actual = auditZhangReceiverCheckpointBoundaryResetEquivalence(source.runtimeId, sourceBytes, actualBytes);
        require(actual.valid && actual.normalizedSourcePayloadSha256 == accepted.normalizedSourcePayloadSha256,
            "ACTUAL_RAW_AFTER_HAS_ADDITIONAL_DIFFERENCES:" + actual.failureReason);
        actualMatched = true;
    }
    std::cout << "RECEIVER_BOUNDARY_RESET_TEST_RESULT " << json::serialize(json::object{
        {"status", "PASS"}, {"receiver_count", source.receivers.size()}, {"ready_reset_count", ready},
        {"observation_workset_reset_receiver_count", worksets}, {"discarded_workset_observations", observations},
        {"source_sha256", accepted.sourcePayloadSha256}, {"normalized_source_sha256", accepted.normalizedSourcePayloadSha256},
        {"actual_raw_after_matched", actualMatched}, {"negative_cases", std::move(negatives)},
        {"immutable_inputs_unchanged", true}, {"PEA_epoch_executed", false}}) << '\n';
    return 0;
}
catch (const std::exception& error)
{ std::cerr << "RECEIVER_BOUNDARY_RESET_TEST_FAILED " << error.what() << '\n'; return 1; }
