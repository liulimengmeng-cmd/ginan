#pragma once

#include <charconv>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

// Counters are metadata, not additional stochastic states. Pure reference
// chart changes must not increment the physical counter. They survive state
// removal/re-add and are copied into the same-epoch disposable AR branch.
inline std::string zhangUserArcGenerationKey(
    bool internalPhysical,
    const std::string& receiver,
    const std::string& satellite,
    int code)
{
    return std::string(internalPhysical ? "USER_AR_PHYSICAL_GENERATION:"
                                        : "USER_AR_STATE_GENERATION:") +
        receiver + ":" + satellite + ":" + std::to_string(code);
}

// Integer authorization is neither physical tracking nor a stochastic FLOAT
// coordinate lifetime. A real-only certified chart map may preserve x/P while
// revoking any earlier integer authorization.
inline std::string zhangUserIntegerAuthorizationGenerationKey(
    const std::string& receiver, const std::string& satellite, int code)
{
    return "USER_AR_INTEGER_AUTHORIZATION_GENERATION:" + receiver + ":" +
        satellite + ":" + std::to_string(code);
}

inline std::uint64_t zhangUserArcGeneration(
    const std::map<std::string, std::string>& metadata,
    const std::string& key)
{
    const auto found = metadata.find(key);
    if (found == metadata.end()) return 0;
    std::uint64_t value = 0;
    const auto& text = found->second;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        throw std::runtime_error("INVALID_USER_AR_GENERATION_METADATA");
    return value;
}

inline void zhangIncrementUserArcGeneration(
    std::map<std::string, std::string>& metadata,
    const std::string& key)
{
    const auto value = zhangUserArcGeneration(metadata, key);
    if (value == std::numeric_limits<std::uint64_t>::max())
        throw std::runtime_error("USER_AR_GENERATION_OVERFLOW");
    metadata[key] = std::to_string(value + 1);
}
