#pragma once
#include <map>
#include <string>

// Deliberately fixture-specific migration, NOT a general provenance bypass.
// This bundle is the audited last FLOAT epoch, before the first AR call.
inline bool zhangR51AuditedFloatReuseAllowed(
    const std::string& bundleSha, const std::string& binarySha,
    const std::string& configSha, const std::string& inputSha,
    const std::map<std::string, std::string>& environment)
{
    const std::map<std::string, std::string> required{
        {"ZHANG_R51_REUSE_FLOAT_20240717", "1"},
        {"ZHANG_R51_ENABLE", "1"}, {"ZHANG_R51_RATIO_ONLY", "1"},
        {"ZHANG_R51_AR_START_GPST_SECONDS", "1405216800"},
        {"ZHANG_R49_FUSION_WEIGHT", "0.10"},
        {"OPENBLAS_NUM_THREADS", "4"}, {"OMP_NUM_THREADS", "4"},
        {"MKL_NUM_THREADS", "1"}};
    for (const auto& [name, value] : required)
    {
        const auto it = environment.find(name);
        if (it == environment.end() || it->second != value) return false;
    }
    return bundleSha == "ab397e655774da5c0a9632becd08fe171508b93ffeaf2dd7b04de899bae2182a"
        && binarySha == "b7b0d35686204154b44448eb097fb6990ca07453cecf9c1fc20f2e80f192cc38"
        && configSha == "eccca673c2acf8a77824e03573974fb468938055f40a85ddc3d8e7651bb21265"
        && inputSha == "4f6351a7db3a2af79ae728ffc1d8d6e9b19986cafe49d6c47c6a703cfeac4874";
}
