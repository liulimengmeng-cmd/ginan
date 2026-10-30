#pragma once

#include "zhangUserPredictionChart.hpp"
#include <map>
#include <set>

inline const std::string& zhangPredictionMetadataPrefix()
{
    static const std::string prefix = "USER_PREDICTION_PRODUCT_V1:";
    return prefix;
}

inline std::string zhangPredictionMetadataKey(const std::string& receiver,
    int system, int signal, const std::string& satellite)
{
    return zhangPredictionMetadataPrefix() + std::to_string(receiver.size()) + ":" + receiver
        + ":" + std::to_string(system) + ":" + std::to_string(signal) + ":" + satellite;
}

inline std::string zhangPredictionEpochText(double epoch)
{
    std::ostringstream out;
    out << std::setprecision(17) << epoch;
    return out.str();
}

struct ZhangUserPredictionChartRecord
{
    ZhangUserPredictionChartEndpoint endpoint;
    std::string observedProof;
    std::string mode;
    std::string reason;
};

inline void zhangPredictionValidateRecord(const ZhangUserPredictionChartRecord& record)
{
    ZhangProductDatumTransition observed;
    if (!zhangPredictionChartEndpointValid(record.endpoint)
        || !zhangReadDatumWireProof(record.observedProof, observed)
        || zhangDatumEndpointIdentity(observed.target) != zhangDatumEndpointIdentity(record.endpoint.applied)
        || (record.mode != "INTEGER_NUMERIC" && record.mode != "RAW_UNSUPPORTED")
        || record.reason.empty() || record.reason.find_first_of("\r\n") != std::string::npos)
        throw std::runtime_error("PREDICTION_PRODUCT_RECORD_INVALID");
}

inline bool zhangPredictionLoadRecord(const std::map<std::string, std::string>& metadata,
    const std::string& base, ZhangUserPredictionChartRecord& destination)
{
    const std::string suffixes[] = {":LEDGER", ":BOUND_ENDPOINT", ":EPOCH", ":OBSERVED_PROOF", ":MODE", ":REASON"};
    std::size_t found = 0;
    for (const auto& suffix : suffixes) found += metadata.count(base + suffix);
    if (found == 0) return false;
    if (found != 6) throw std::runtime_error("PREDICTION_PRODUCT_RECORD_INCOMPLETE");
    ZhangUserPredictionChartRecord candidate;
    if (!zhangPredictionChartDecode(metadata.at(base + ":LEDGER"), candidate.endpoint))
        throw std::runtime_error("PREDICTION_PRODUCT_LEDGER_INVALID");
    candidate.observedProof = metadata.at(base + ":OBSERVED_PROOF");
    candidate.mode = metadata.at(base + ":MODE");
    candidate.reason = metadata.at(base + ":REASON");
    zhangPredictionValidateRecord(candidate);
    if (metadata.at(base + ":BOUND_ENDPOINT") != zhangDatumEndpointIdentity(candidate.endpoint.applied)
        || metadata.at(base + ":EPOCH") != zhangPredictionEpochText(candidate.endpoint.applied.epoch))
        throw std::runtime_error("PREDICTION_PRODUCT_MIRROR_MISMATCH");
    destination = std::move(candidate);
    return true;
}

// Destination changes only after validation/serialization. Production calls
// this within the existing private KF reference/product transaction.
inline void zhangPredictionStoreRecord(std::map<std::string, std::string>& metadata,
    const std::string& base, const ZhangUserPredictionChartRecord& record)
{
    zhangPredictionValidateRecord(record);
    const auto ledger = zhangPredictionChartEncode(record.endpoint);
    metadata[base + ":LEDGER"] = ledger;
    metadata[base + ":BOUND_ENDPOINT"] = zhangDatumEndpointIdentity(record.endpoint.applied);
    metadata[base + ":EPOCH"] = zhangPredictionEpochText(record.endpoint.applied.epoch);
    metadata[base + ":OBSERVED_PROOF"] = record.observedProof;
    metadata[base + ":MODE"] = record.mode;
    metadata[base + ":REASON"] = record.reason;
}

inline std::string zhangPredictionCaptureMetadata(
    const std::map<std::string, std::string>& metadata, double completedEpoch)
{
    if (!std::isfinite(completedEpoch)) throw std::runtime_error("PREDICTION_CHECKPOINT_TIME_INVALID");
    std::map<std::string, std::string> selected;
    std::set<std::string> bases;
    for (const auto& [key, value] : metadata)
    {
        if (!key.starts_with(zhangPredictionMetadataPrefix())) continue;
        const auto pos = key.rfind(':');
        if (pos == std::string::npos) throw std::runtime_error("PREDICTION_CHECKPOINT_KEY_INVALID");
        const auto suffix = key.substr(pos);
        if (suffix != ":LEDGER" && suffix != ":BOUND_ENDPOINT" && suffix != ":EPOCH"
            && suffix != ":OBSERVED_PROOF" && suffix != ":MODE" && suffix != ":REASON")
            throw std::runtime_error("PREDICTION_CHECKPOINT_UNKNOWN_FIELD");
        bases.insert(key.substr(0, pos)); selected.emplace(key, value);
    }
    if (selected.size() > 10000) throw std::runtime_error("PREDICTION_CHECKPOINT_SIZE_EXCEEDED");
    for (const auto& base : bases)
    {
        ZhangUserPredictionChartRecord record;
        if (!zhangPredictionLoadRecord(metadata, base, record)
            || record.endpoint.applied.epoch > completedEpoch + 1e-6)
            throw std::runtime_error("PREDICTION_CHECKPOINT_FUTURE_OR_MISSING_RECORD");
    }
    std::ostringstream out;
    out << "US20_PREDICTION_METADATA_V1 " << zhangPredictionEpochText(completedEpoch) << ' ' << selected.size();
    for (const auto& [key, value] : selected) out << ' ' << std::quoted(key) << ' ' << std::quoted(value);
    return out.str();
}

// The core restores the KF metadata itself. This mandatory independent mirror
// validates its exact causal contents before the core or globals are committed.
inline void zhangPredictionValidateMetadata(const std::string& payload,
    const std::map<std::string, std::string>& metadata, double completedEpoch)
{
    if (payload != zhangPredictionCaptureMetadata(metadata, completedEpoch))
        throw std::runtime_error("PREDICTION_CHECKPOINT_METADATA_MISMATCH");
}
