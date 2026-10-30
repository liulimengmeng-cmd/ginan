#pragma once

#include <cstddef>

// Wire-byte contract, not a memory allowance: the native process still runs
// under its independent address-space/RSS limits. A real 134-station G/E/C
// epoch contains about 468 MiB of complete current-epoch evidence. Keep the
// existing format and every physical/numeric/proof term; bound both producers
// and consumers identically instead of writing a section we cannot read.
inline constexpr std::size_t ZHANG_GEC_AR_RUNTIME_MAX_PAYLOAD_BYTES =
    std::size_t{1024} * 1024 * 1024;

inline constexpr bool zhangGecArRuntimePayloadSizeValid(std::size_t bytes) noexcept
{
    return bytes != 0 && bytes <= ZHANG_GEC_AR_RUNTIME_MAX_PAYLOAD_BYTES;
}
