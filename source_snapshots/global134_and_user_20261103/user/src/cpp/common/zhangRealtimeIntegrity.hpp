#pragma once
#include "common/zhangCheckpoint.hpp"
#include <boost/json.hpp>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>

inline void zhangVerifyRealtimeNavigationManifest(
    const boost::json::object& manifest, const std::filesystem::path& directory)
{
    const auto* extension = manifest.if_contains("realtime_navigation");
    const char* mode = std::getenv("ZHANG_REALTIME_NAV_MODE");
    const bool required = mode && std::string(mode) != "OFF";
    if (!extension)
    {
        if (required) throw std::runtime_error("REALTIME_NAV_MANIFEST_MISSING");
        return;
    }
    const auto& entry = extension->as_object();
    const auto bytes = entry.at("bytes").to_number<std::uint64_t>();
    if (entry.at("schema").as_string() != "REALTIME_NAV_LEGACY_GEC_V1" ||
        bytes == 0 || bytes > 16 * 1024 * 1024)
        throw std::runtime_error("REALTIME_NAV_MANIFEST_INVALID");
    const auto path = directory / "navigation.json";
    std::string failure;
    if (!std::filesystem::is_regular_file(path) || std::filesystem::file_size(path) != bytes ||
        zhangCheckpointFileSha256(path.string(), &failure) != entry.at("sha256").as_string())
        throw std::runtime_error("REALTIME_NAV_SNAPSHOT_HASH_MISMATCH");
    if (required)
    {
        const char* generation = std::getenv("ZHANG_REALTIME_GENERATION");
        const char* model = std::getenv("ZHANG_REALTIME_MODEL_SHA256");
        if (!generation || !model || entry.at("generation").as_string() != generation ||
            entry.at("model_sha256").as_string() != model || entry.at("mode").as_string() != mode)
            throw std::runtime_error("REALTIME_NAV_GENERATION_OR_MODEL_MISMATCH");
    }
}
