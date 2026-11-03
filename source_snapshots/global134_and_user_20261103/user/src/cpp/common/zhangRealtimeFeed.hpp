#pragma once

#include "common/gTime.hpp"
#include "common/zhangR51CommittedBundle.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>

// One process pins one immutable commit before constructing any observation
// equations. Per-measurement lookups never wait, rescan or select another epoch.
namespace zhangRealtime
{
namespace fs = std::filesystem;

inline std::string mode()
{
    const char* value = std::getenv("ZHANG_REALTIME_MODE");
    const std::string result = value ? value : "OFF";
    if (result != "OFF" && result != "ON" && result != "BUNDLE_REPLAY")
        throw std::runtime_error("REALTIME_MODE_INVALID");
    return result;
}

inline bool enabled() { return mode() != "OFF"; }

inline long long epochKey(const GTime& time)
{
    const long double value = time.bigTime;
    if (!std::isfinite(value) || std::abs(value - std::round(value)) > 1e-6L)
        throw std::runtime_error("REALTIME_NON_INTEGER_EPOCH");
    return std::llround(value);
}

inline fs::path root()
{
    const char* value = std::getenv("ZHANG_REALTIME_PRODUCT_ROOT");
    if (!value || !*value) throw std::runtime_error("REALTIME_PRODUCT_ROOT_REQUIRED");
    const auto result = fs::canonical(value);
    if (!fs::is_directory(result)) throw std::runtime_error("REALTIME_PRODUCT_ROOT_INVALID");
    return result;
}

inline double waitSeconds()
{
    const char* value = std::getenv("ZHANG_REALTIME_PRODUCT_WAIT_SECONDS");
    if (!value) return 20;
    std::size_t consumed = 0;
    const double seconds = std::stod(value, &consumed);
    if (value[consumed] != '\0' || !std::isfinite(seconds) || seconds < 0 || seconds > 60)
        throw std::runtime_error("REALTIME_WAIT_SECONDS_INVALID");
    return seconds;
}

struct Pin
{
    long long epoch = std::numeric_limits<long long>::min();
    ZhangR51VerifiedBundle bundle;
    std::string commitHash;
    std::string solution;
    std::string reason = "NOT_PREPARED";
    double waitMilliseconds = 0;
    bool prepared = false;
};

inline Pin pinned;
inline std::mutex pinMutex;

inline Pin discover(const GTime& time)
{
    Pin result;
    result.epoch = epochKey(time);
    result.prepared = true;
    result.reason = "COMMITTED_EPOCH_UNAVAILABLE";
    const auto base = root();
    std::string prefix = time.to_string(0);
    for (auto& c : prefix) if (c == ' ' || c == ':') c = '_';
    prefix += '_';
    fs::path selected;
    std::size_t entries = 0;
    for (const auto& entry : fs::directory_iterator(base))
    {
        if (++entries > 8192) throw std::runtime_error("REALTIME_ACTIVE_DIRECTORY_LIMIT");
        const auto name = entry.path().filename().string();
        if (!name.starts_with(prefix) || name.ends_with(".pending") || !entry.is_directory()) continue;
        if (!fs::is_regular_file(entry.path() / "COMMIT.json")) continue;
        if (!selected.empty()) throw std::runtime_error("REALTIME_DUPLICATE_COMMITTED_EPOCH");
        selected = fs::canonical(entry.path());
        if (selected.parent_path() != base) throw std::runtime_error("REALTIME_BUNDLE_ESCAPES_ROOT");
    }
    if (selected.empty()) return result;
    const std::string before = zhangR51BundleFileHash(selected / "COMMIT.json");
    result.bundle = zhangR51VerifyCommittedBundle(selected, false);
    if (!result.bundle.valid)
        throw std::runtime_error("REALTIME_COMMITTED_BUNDLE_INVALID:" + result.bundle.failureReason);
    if (result.bundle.epoch != time.to_string(0)) throw std::runtime_error("REALTIME_COMMIT_EPOCH_MISMATCH");
    if (before != zhangR51BundleFileHash(selected / "COMMIT.json"))
        throw std::runtime_error("REALTIME_COMMIT_CHANGED_DURING_VERIFY");
    result.commitHash = before;
    std::ifstream manifestInput(selected / "COMMIT.json");
    std::ostringstream manifestText; manifestText << manifestInput.rdbuf();
    const auto manifest = boost::json::parse(manifestText.str()).as_object();
    const auto* regional = manifest.if_contains("regional_joint");
    if (!regional || !regional->is_object()) throw std::runtime_error("REALTIME_REGIONAL_MANIFEST_REQUIRED");
    const auto& files = regional->as_object().at("files").as_object();
    for (const auto* name : {"regional_parameters.csv", "regional_covariance.csv", "regional_linear_map.csv",
                            "regional_quality.csv", "regional_contract.txt"})
    {
        const auto& descriptor = files.at(name).as_object();
        const auto path = selected / name;
        if (fs::file_size(path) != descriptor.at("bytes").to_number<std::uint64_t>() ||
            zhangR51BundleFileHash(path) != descriptor.at("sha256").as_string())
            throw std::runtime_error(std::string("REALTIME_REGIONAL_HASH_MISMATCH:") + name);
    }
    std::ifstream contract(selected / "regional_contract.txt");
    std::string line;
    while (std::getline(contract, line))
        if (line.starts_with("solution="))
        {
            if (!result.solution.empty()) throw std::runtime_error("REALTIME_DUPLICATE_SOLUTION");
            result.solution = line.substr(9);
        }
    if ((result.solution != "FLOAT" && result.solution != "PRODUCT_FIXED") ||
        (result.solution == "PRODUCT_FIXED" && !result.bundle.ar))
        throw std::runtime_error("REALTIME_PUBLISHED_SOLUTION_AUTHORITY_MISMATCH");
    result.reason = "NONE";
    return result;
}

inline Pin begin(const GTime& time)
{
    std::lock_guard<std::mutex> lock(pinMutex);
    const long long epoch = epochKey(time);
    if (pinned.prepared && pinned.epoch == epoch) return pinned;
    if (pinned.prepared && epoch < pinned.epoch)
        throw std::runtime_error("REALTIME_EPOCH_REGRESSION");
    const auto started = std::chrono::steady_clock::now();
    double allowance = waitSeconds();
    if (mode() == "ON")
    {
        const double age = (timeGet() - time).to_double();
        if (age < -2) throw std::runtime_error("REALTIME_FUTURE_OBSERVATION_EPOCH");
        if (age >= 90)
        {
            pinned = {};
            pinned.epoch = epoch; pinned.prepared = true;
            pinned.reason = "END_TO_END_DEADLINE_EXCEEDED";
            return pinned;
        }
        allowance = std::min(allowance, 90 - age);
    }
    const auto deadline = started + std::chrono::duration<double>(allowance);
    do
    {
        pinned = discover(time);
        if (pinned.bundle.valid) break;
        if (std::chrono::steady_clock::now() >= deadline) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    } while (true);
    pinned.waitMilliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return pinned;
}

inline Pin current(const GTime& time)
{
    std::lock_guard<std::mutex> lock(pinMutex);
    if (!pinned.prepared || pinned.epoch != epochKey(time))
        throw std::runtime_error("REALTIME_EPOCH_NOT_PREPARED");
    return pinned;
}

inline void verifyPinned(const GTime& time)
{
    const auto value = current(time);
    if (!value.bundle.valid) return;
    const auto again = discover(time);
    if (!again.bundle.valid || again.bundle.directory != value.bundle.directory ||
        again.commitHash != value.commitHash)
        throw std::runtime_error("REALTIME_PIN_IDENTITY_CHANGED");
}
} // namespace zhangRealtime

// Implemented in the existing user translation units; forces all file reads
// and covariance factorization to finish before equation construction begins.
bool prepareZhangRealtimeUserProducts(const GTime& time);
bool prepareZhangRealtimeRegional(const GTime& time);
