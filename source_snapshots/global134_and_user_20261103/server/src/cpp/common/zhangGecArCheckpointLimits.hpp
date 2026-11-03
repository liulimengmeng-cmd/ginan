#pragma once

#include <cstddef>
#include <cstdint>

// Wire-byte contract, not a memory allowance. V2 factors seven identical
// binding scope fields into their enclosing evidence; no equation, numerical
// row, full state-order identity or proof is removed. The 1 GiB wire cap is
// unchanged and applies identically to legacy and scoped producers/readers.
inline constexpr std::size_t ZHANG_GEC_AR_RUNTIME_MAX_PAYLOAD_BYTES =
    std::size_t{1024} * 1024 * 1024;

inline constexpr bool zhangGecArRuntimePayloadSizeValid(std::size_t bytes) noexcept
{
    return bytes != 0 && bytes <= ZHANG_GEC_AR_RUNTIME_MAX_PAYLOAD_BYTES;
}

// A small scoped wire must not request unbounded reconstruction of long
// binding identities. This is a separate decoded-string allocation bound,
// not permission to increase the wire cap or alter native resource limits.
inline constexpr std::uint64_t ZHANG_GEC_AR_MAX_EXPANDED_BINDING_SCOPE_BYTES =
    std::uint64_t{4} * 1024 * 1024 * 1024;

inline bool zhangGecArAccumulateExpandedScopeBytes(std::size_t count,
    std::size_t bytesPerBinding, std::uint64_t& total) noexcept
{
    constexpr auto cap = ZHANG_GEC_AR_MAX_EXPANDED_BINDING_SCOPE_BYTES;
    if (total > cap || (bytesPerBinding != 0 && count > (cap - total) / bytesPerBinding))
        return false;
    total += static_cast<std::uint64_t>(count) * bytesPerBinding;
    return true;
}
