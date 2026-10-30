// Pending execution: no build or test run has been authorised for P0/P1.
// This tests the pure model selector, not graph correctness, checkpoint
// migration, integer authority, publication, or scientific product accuracy.
#include "common/zhangProductBasisMode.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
int checks=0;
void require(bool condition,const char* message)
{
    ++checks;
    if(!condition) throw std::runtime_error(message);
}
template<class F> bool rejected(F&& action)
{
    try {action();}
    catch(const std::invalid_argument&) {return true;}
    return false;
}
} // namespace

int main()
{
    try
    {
        require(std::string(ZhangDefaultProductBasisModeName)=="SINGLE_S_BASIS_NATIVE",
            "isolated candidate default silently restored the legacy product tree");
        require(zhangParseProductBasisMode("SINGLE_S_BASIS_NATIVE")==
            ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE,
            "literal native mode did not select native");
        require(zhangParseProductBasisMode("LEGACY_PRODUCT_TREE")==
            ZhangProductBasisMode::LEGACY_PRODUCT_TREE,
            "explicit literal legacy mode was not retained");
        require(zhangSingleSBasisNativeEnabled("SINGLE_S_BASIS_NATIVE") &&
            !zhangSingleSBasisNativeEnabled("LEGACY_PRODUCT_TREE"),
            "mode flag disagrees with the two explicit model definitions");
        require(zhangSingleSBasisNativeEnabled(ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE) &&
            !zhangSingleSBasisNativeEnabled(ZhangProductBasisMode::LEGACY_PRODUCT_TREE),
            "enum selector disagrees with the canonical model definitions");
        require(std::string(zhangProductBasisModeName(
            ZhangProductBasisMode::SINGLE_S_BASIS_NATIVE))=="SINGLE_S_BASIS_NATIVE" &&
            std::string(zhangProductBasisModeName(
            ZhangProductBasisMode::LEGACY_PRODUCT_TREE))=="LEGACY_PRODUCT_TREE",
            "canonical mode names cannot identify distinct checkpoint schemas");
        const std::vector<std::string> invalid={
            "", "AUTO", "NATIVE", "SINGLE", "PRODUCT_TREE", "HYBRID_STABLE",
            "single_s_basis_native", "legacy_product_tree", "SINGLE_S_BASIS_NATIVE ",
            " LEGACY_PRODUCT_TREE", "SINGLE_S_BASIS_NATIVE\n", "0", "1", "true", "false"
        };
        for(const auto& setting:invalid)
        {
            require(rejected([&] {(void)zhangParseProductBasisMode(setting);}),
                "unknown or aliased model value was accepted");
            require(rejected([&] {(void)zhangSingleSBasisNativeEnabled(setting);}),
                "unknown value silently selected the legacy false branch");
        }
        std::string embedded="SINGLE_S_BASIS_NATIVE";
        embedded.push_back('\0');embedded+="LEGACY_PRODUCT_TREE";
        require(rejected([&] {(void)zhangParseProductBasisMode(embedded);}),
            "embedded-null model input was silently truncated");
        const auto invalidEnum=static_cast<ZhangProductBasisMode>(-1);
        require(rejected([&] {(void)zhangProductBasisModeName(invalidEnum);}) &&
            rejected([&] {(void)zhangSingleSBasisNativeEnabled(invalidEnum);}),
            "invalid enum silently selected a legacy path");
        const std::string preserved="LEGACY_PRODUCT_TREE";
        (void)zhangParseProductBasisMode(preserved);
        require(preserved=="LEGACY_PRODUCT_TREE","pure selector mutated its input");
        std::cout<<"single_s_basis_mode: "<<checks<<" checks passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"single_s_basis_mode failed: "<<error.what()<<"\n";
        return 1;
    }
}
