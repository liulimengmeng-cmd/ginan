#include "common/zhangR51FloatReuse.hpp"

#include <iostream>
#include <stdexcept>

int main()
{
    std::map<std::string, std::string> environment{
        {"ZHANG_R51_REUSE_FLOAT_20240717", "1"},
        {"ZHANG_R51_ENABLE", "1"},
        {"ZHANG_R51_RATIO_ONLY", "1"},
        {"ZHANG_R51_AR_START_GPST_SECONDS", "1405216800"},
        {"ZHANG_R49_FUSION_WEIGHT", "0.10"},
        {"OPENBLAS_NUM_THREADS", "4"},
        {"OMP_NUM_THREADS", "4"},
        {"MKL_NUM_THREADS", "1"}};
    std::string hashes[]{
        "ab397e655774da5c0a9632becd08fe171508b93ffeaf2dd7b04de899bae2182a",
        "b7b0d35686204154b44448eb097fb6990ca07453cecf9c1fc20f2e80f192cc38",
        "02bbaabbe2a836ec3dfc08d2ad206326c9dbddf41237692e0d3806b0043b764f",
        "4f6351a7db3a2af79ae728ffc1d8d6e9b19986cafe49d6c47c6a703cfeac4874"};
    const auto allowed = [&] {
        return zhangR51AuditedFloatReuseAllowed(
            hashes[0], hashes[1], hashes[2], hashes[3], environment);
    };
    const auto require = [](bool condition, const char* reason) {
        if (!condition) throw std::runtime_error(reason);
    };

    require(allowed(), "audited 4/4 source rejected");
    for (auto& hash : hashes)
    {
        const auto original = hash;
        hash = "WRONG";
        require(!allowed(), "wrong source hash accepted under 4/4");
        hash = original;
    }
    const auto originalEnvironment = environment;
    for (const auto& [name, value] : originalEnvironment)
    {
        environment[name] = "WRONG";
        require(!allowed(), "wrong source environment accepted");
        environment = originalEnvironment;
        environment.erase(name);
        require(!allowed(), "missing source environment accepted");
        environment = originalEnvironment;
    }

    environment["OMP_NUM_THREADS"] = "8";
    environment["OPENBLAS_NUM_THREADS"] = "1";
    require(!allowed(), "8/1 without explicit experiment marker accepted");
    environment["ZHANG_R51_REUSE_FLOAT_THREADS_8_1"] = "1";
    require(allowed(), "controlled 8/1 source rejected");
    for (auto& hash : hashes)
    {
        const auto original = hash;
        hash = "WRONG";
        require(!allowed(), "wrong source hash accepted under 8/1");
        hash = original;
    }
    environment["OMP_NUM_THREADS"] = "7";
    require(!allowed(), "7/1 accepted");
    environment["OMP_NUM_THREADS"] = "8";
    environment["OPENBLAS_NUM_THREADS"] = "4";
    require(!allowed(), "8/4 accepted");
    environment["OPENBLAS_NUM_THREADS"] = "1";
    environment["ZHANG_R51_REUSE_FLOAT_THREADS_8_1"] = "0";
    require(!allowed(), "8/1 with invalid marker accepted");
    environment["ZHANG_R51_REUSE_FLOAT_THREADS_8_1"] = "1";
    environment["OMP_NUM_THREADS"] = "4";
    environment["OPENBLAS_NUM_THREADS"] = "4";
    require(!allowed(), "4/4 with 8/1 marker accepted");

    std::cout << "PASS pinned checkpoint, 4/4 baseline and explicit 8/1 experiment\n";
}
