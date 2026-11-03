#pragma once

#include <chrono>
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#endif

std::string zhangCheckpointSha256(const std::string& bytes);

/** Cold checkpoint evidence only. There is intentionally no KF insertion,
 * conditioning, Pzz restore, CURRENT AR registry, or permission conversion.
 * The parent is the complete already-written checkpoint payload, which keeps
 * the actual Gaussian joint state including every cross covariance. */
namespace zhang_history_archive
{
inline constexpr const char* directoryName = "history_evidence";
inline constexpr std::size_t maximumObjectBytes = 64ULL * 1024 * 1024;
inline constexpr std::size_t maximumManifestBytes = 4ULL * 1024 * 1024;
inline constexpr std::size_t maximumObjects = 32768;

struct Scope
{
    std::string owner, epoch, checkpointId, parentCheckpointId;
    std::string snapshotPayloadSha256;
    bool operator==(const Scope&) const = default;
};

struct Evidence
{
    std::string kind, metadataKey, payload;
    std::vector<std::string> parents;
    bool operator==(const Evidence&) const = default;
};

using VerifySnapshot = std::function<bool(const std::string&, std::string&)>;

inline bool shaValid(const std::string& value)
{
    return value.size() == 64 && value.find_first_not_of("0123456789abcdef") == std::string::npos;
}

inline bool scopeValid(const Scope& scope)
{
    return !scope.owner.empty() && !scope.epoch.empty() && !scope.checkpointId.empty() &&
        scope.parentCheckpointId != scope.checkpointId && shaValid(scope.snapshotPayloadSha256);
}

inline std::string encode(const Evidence& evidence)
{
    std::ostringstream out;
    out << "HISTORICAL_EVIDENCE_ONLY_V1 " << std::quoted(evidence.kind) << ' '
        << std::quoted(evidence.metadataKey) << ' ' << evidence.parents.size();
    for (const auto& parent : evidence.parents) out << ' ' << parent;
    out << ' ' << std::quoted(evidence.payload) << '\n';
    return out.str();
}

inline bool decode(const std::string& bytes, const Scope& scope, Evidence& evidence)
{
    std::istringstream in(bytes);
    std::string magic;
    std::size_t parents = 0;
    Evidence candidate;
    if (!(in >> magic >> std::quoted(candidate.kind) >> std::quoted(candidate.metadataKey) >> parents) ||
        magic != "HISTORICAL_EVIDENCE_ONLY_V1" ||
        (candidate.kind != "TRANSPORT_ONLY_CERTIFICATE" && candidate.kind != "FLOAT_TARGET_DEFINITION" &&
         candidate.kind != "COMPLETE_PARENT_METADATA_REFERENCE") ||
        candidate.metadataKey.empty() || parents != 1) return false;
    candidate.parents.resize(parents);
    for (auto& parent : candidate.parents)
        if (!(in >> parent) || parent != scope.snapshotPayloadSha256) return false;
    if (!(in >> std::quoted(candidate.payload))) return false;
    in >> std::ws;
    if (!in.eof() || encode(candidate) != bytes) return false;
    evidence = std::move(candidate);
    return true;
}

inline bool readBytes(const std::filesystem::path& path, std::size_t limit, std::string& bytes)
{
    std::error_code error;
    if (!std::filesystem::is_regular_file(std::filesystem::symlink_status(path, error)) || error)
        return false;
    const auto size = std::filesystem::file_size(path, error);
    if (error || size > limit) return false;
    std::ifstream input(path, std::ios::binary);
    if (!input) return false;
    // Allocation is bounded by the checked size. A file that grows after
    // stat is rejected by the final EOF check, never read without a limit.
    std::string candidate(static_cast<std::size_t>(size), '\0');
    if (size) input.read(candidate.data(), static_cast<std::streamsize>(size));
    if ((size && static_cast<std::uintmax_t>(input.gcount()) != size) || input.bad() ||
        input.peek() != std::char_traits<char>::eof()) return false;
    bytes = std::move(candidate);
    return true;
}

inline bool syncFile(const std::filesystem::path& path)
{
#ifdef _WIN32
    const int fd = _wopen(path.c_str(), _O_RDWR | _O_BINARY);
    if (fd < 0) return false;
    const bool good = _commit(fd) == 0;
    _close(fd);
#else
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;
    const bool good = ::fsync(fd) == 0;
    ::close(fd);
#endif
    return good;
}

inline bool syncDirectory(const std::filesystem::path& path)
{
#ifdef _WIN32
    // File commits and the subsequent MOVEFILE_WRITE_THROUGH rename carry
    // the persistence boundary on Windows; POSIX also syncs each directory.
    (void)path;
    return true;
#else
    const int fd = ::open(path.c_str(), O_RDONLY | O_DIRECTORY);
    if (fd < 0) return false;
    const bool good = ::fsync(fd) == 0;
    ::close(fd);
    return good;
#endif
}

inline bool renameDurably(const std::filesystem::path& source, const std::filesystem::path& target)
{
#ifdef _WIN32
    if (!MoveFileExW(source.c_str(), target.c_str(), MOVEFILE_WRITE_THROUGH)) return false;
#else
    std::error_code error;
    std::filesystem::rename(source, target, error);
    if (error) return false;
#endif
    return syncDirectory(target.parent_path());
}

inline bool writeImmutableWithTemporary(const std::filesystem::path& path,
    const std::string& bytes, const std::filesystem::path& temporary)
{
    // The temp must be in the same directory for an atomic no-replace
    // installation. O_EXCL means a collision cannot truncate another file.
    if (temporary.parent_path() != path.parent_path() || temporary == path) return false;
#ifdef _WIN32
    const int fd = _wopen(temporary.c_str(), _O_WRONLY | _O_CREAT | _O_EXCL | _O_BINARY,
        _S_IREAD | _S_IWRITE);
#else
    const int fd = ::open(temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0600);
#endif
    if (fd < 0) return false; // do not remove a temp owned by somebody else
    bool good = true;
    std::size_t offset = 0;
    while (offset < bytes.size())
    {
        const auto count = std::min<std::size_t>(bytes.size() - offset, 1024 * 1024);
#ifdef _WIN32
        const auto written = _write(fd, bytes.data() + offset, static_cast<unsigned int>(count));
#else
        const auto written = ::write(fd, bytes.data() + offset, count);
#endif
        if (written < 0 && errno == EINTR) continue;
        if (written <= 0) { good = false; break; }
        offset += static_cast<std::size_t>(written);
    }
#ifdef _WIN32
    good = good && _commit(fd) == 0;
    good = _close(fd) == 0 && good;
#else
    good = good && ::fsync(fd) == 0;
    good = ::close(fd) == 0 && good;
#endif
    if (!good) return false; // preserve this failed temp for investigation
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_WRITE_THROUGH)) return false;
#else
    // link is atomic and never replaces an existing destination, even if
    // another writer installed it AFTER this function created its temp.
    if (::link(temporary.c_str(), path.c_str()) != 0) return false;
    if (::unlink(temporary.c_str()) != 0) return false;
#endif
    return syncDirectory(path.parent_path());
}

inline bool writeImmutable(const std::filesystem::path& path, const std::string& bytes)
{
    const std::filesystem::path temporary = path.string() + ".tmp." +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    return writeImmutableWithTemporary(path, bytes, temporary);
}

inline std::string manifestBody(const Scope& scope, const std::vector<std::string>& objects)
{
    std::ostringstream out;
    out << "HISTORICAL_EVIDENCE_MANIFEST_V1 " << std::quoted(scope.owner) << ' '
        << std::quoted(scope.epoch) << ' ' << std::quoted(scope.checkpointId) << ' '
        << std::quoted(scope.parentCheckpointId) << ' ' << scope.snapshotPayloadSha256 << ' '
        << objects.size() << '\n';
    for (const auto& object : objects) out << object << '\n';
    return out.str();
}

inline bool read(const std::filesystem::path& directory, const Scope& expected,
    const VerifySnapshot& verifySnapshot, std::vector<Evidence>& result, std::string& reason)
{
    if (!scopeValid(expected) || !verifySnapshot)
    { reason = "HISTORY_SCOPE_OR_PARENT_VERIFIER_INVALID"; return false; }
    std::string manifest;
    if (!readBytes(directory / "manifest", maximumManifestBytes, manifest) || manifest.size() < 65 ||
        manifest[64] != '\n' || !shaValid(manifest.substr(0, 64)) ||
        zhangCheckpointSha256(manifest.substr(65)) != manifest.substr(0, 64))
    { reason = "HISTORY_MANIFEST_MISSING_PARTIAL_OR_CORRUPT"; return false; }
    std::istringstream in(manifest.substr(65));
    Scope scope;
    std::string magic;
    std::size_t count = 0;
    if (!(in >> magic >> std::quoted(scope.owner) >> std::quoted(scope.epoch) >>
        std::quoted(scope.checkpointId) >> std::quoted(scope.parentCheckpointId) >>
        scope.snapshotPayloadSha256 >> count) || magic != "HISTORICAL_EVIDENCE_MANIFEST_V1" ||
        scope != expected || count > maximumObjects)
    { reason = "HISTORY_MANIFEST_SCOPE_OR_COUNT_MISMATCH"; return false; }
    if (!verifySnapshot(scope.snapshotPayloadSha256, reason))
    { reason = "HISTORY_COMPLETE_SNAPSHOT_PARENT_UNAVAILABLE:" + reason; return false; }
    std::vector<Evidence> candidate;
    std::vector<std::string> ids;
    std::set<std::string> unique, metadataKeys;
    for (std::size_t index = 0; index < count; ++index)
    {
        std::string id, bytes;
        Evidence evidence;
        if (!(in >> id) || !shaValid(id) || !unique.insert(id).second ||
            !readBytes(directory / "objects" / (id + ".payload"), maximumObjectBytes, bytes) ||
            zhangCheckpointSha256(bytes) != id || !decode(bytes, scope, evidence) ||
            !metadataKeys.insert(evidence.metadataKey).second)
        { reason = "HISTORY_OBJECT_MISSING_PARTIAL_CORRUPT_OR_PARENT_MISMATCH"; return false; }
        ids.push_back(id);
        candidate.push_back(std::move(evidence));
    }
    in >> std::ws;
    if (!in.eof() || manifestBody(scope, ids) != manifest.substr(65))
    { reason = "HISTORY_MANIFEST_TRAILING_OR_NONCANONICAL_DATA"; return false; }
    if (!verifySnapshot(scope.snapshotPayloadSha256, reason))
    { reason = "HISTORY_PARENT_CHANGED_DURING_READ:" + reason; return false; }
    result = std::move(candidate); // output changes only after full verification
    reason = "HISTORY_VERIFIED_TRANSPORT_ONLY_NO_CURRENT_AR_AUTHORITY";
    return true;
}

inline bool write(const std::filesystem::path& directory, const Scope& scope,
    const std::vector<Evidence>& evidence, const VerifySnapshot& verifySnapshot, std::string& reason)
{
    if (!scopeValid(scope) || evidence.size() > maximumObjects || !verifySnapshot ||
        !verifySnapshot(scope.snapshotPayloadSha256, reason))
    { reason = "HISTORY_WRITE_SCOPE_OR_COMPLETE_PARENT_INVALID:" + reason; return false; }
    std::error_code error;
    if (std::filesystem::exists(directory, error) || error ||
        !std::filesystem::create_directory(directory, error) || error ||
        !std::filesystem::create_directory(directory / "objects", error) || error)
    { reason = "HISTORY_IMMUTABLE_DIRECTORY_CREATE_FAILED"; return false; }
    std::vector<std::string> ids;
    std::set<std::string> unique, metadataKeys;
    for (const auto& item : evidence)
    {
        const auto bytes = encode(item);
        Evidence decoded;
        const auto id = zhangCheckpointSha256(bytes);
        if (bytes.size() > maximumObjectBytes || !shaValid(id) || !decode(bytes, scope, decoded) ||
            !unique.insert(id).second || !metadataKeys.insert(item.metadataKey).second ||
            !writeImmutable(directory / "objects" / (id + ".payload"), bytes))
        { reason = "HISTORY_OBJECT_DURABLE_WRITE_FAILED"; return false; }
        ids.push_back(id);
    }
    const auto body = manifestBody(scope, ids);
    if (body.size() + 65 > maximumManifestBytes ||
        !verifySnapshot(scope.snapshotPayloadSha256, reason) ||
        !syncDirectory(directory / "objects") ||
        !writeImmutable(directory / "manifest", zhangCheckpointSha256(body) + '\n' + body) ||
        !syncDirectory(directory))
    { reason = "HISTORY_MANIFEST_DURABLE_COMMIT_FAILED:" + reason; return false; }
    std::vector<Evidence> verified;
    return read(directory, scope, verifySnapshot, verified, reason) && verified == evidence;
}
}
