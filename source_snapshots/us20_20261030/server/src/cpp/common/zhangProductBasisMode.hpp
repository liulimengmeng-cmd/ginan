#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

/** The product basis is an explicit model choice, not an AR certificate.
 * SINGLE uses the actual current estimator S-basis as the only authority;
 * LEGACY alone permits a separate product tree.  Never interpret an unknown
 * spelling as legacy, false, or an automatic fallback.  This header is pure
 * and deliberately independent of acsConfig, KFState and graph types. */
enum class ZhangProductBasisMode
{
    SINGLE_S_BASIS_NATIVE,
    LEGACY_PRODUCT_TREE
};

// This default changes ONLY in the isolated P0/P1 candidate source.  An old
// experiment/config/checkpoint is not thereby migrated or validated.
inline constexpr const char* ZhangDefaultProductBasisModeName="SINGLE_S_BASIS_NATIVE";

inline ZhangProductBasisMode zhangParseProductBasisMode(std::string_view setting)
{
    if(setting=="SINGLE_S_BASIS_NATIVE") return ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE;
    if(setting=="LEGACY_PRODUCT_TREE") return ZhangProductBasisMode::LEGACY_PRODUCT_TREE;
    throw std::invalid_argument(
        "ZHANG_PRODUCT_BASIS_MODE_INVALID: '"+std::string(setting)+
        "'; expected SINGLE_S_BASIS_NATIVE or explicit LEGACY_PRODUCT_TREE");
}

inline const char* zhangProductBasisModeName(ZhangProductBasisMode mode)
{
    switch(mode)
    {
        case ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE: return "SINGLE_S_BASIS_NATIVE";
        case ZhangProductBasisMode::LEGACY_PRODUCT_TREE: return "LEGACY_PRODUCT_TREE";
    }
    throw std::invalid_argument("ZHANG_PRODUCT_BASIS_MODE_INVALID_ENUM");
}

inline bool zhangSingleSBasisNativeEnabled(std::string_view setting)
{
    return zhangParseProductBasisMode(setting)==ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE;
}

inline bool zhangSingleSBasisNativeEnabled(ZhangProductBasisMode mode)
{
    // Validate even enum inputs: an invalid cast must not silently select the
    // legacy path in a caller's if/else branch.
    (void)zhangProductBasisModeName(mode);
    return mode==ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE;
}
