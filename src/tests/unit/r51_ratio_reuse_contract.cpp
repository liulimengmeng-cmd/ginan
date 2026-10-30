#include "common/zhangR51FloatReuse.hpp"
#include <iostream>
#include <stdexcept>

int main()
{
    const auto require = [](bool condition, const char* reason) {
        if (!condition) throw std::runtime_error(reason);
    };
    std::map<std::string, std::string> environment{
        {"ZHANG_R51_REUSE_RATIO_20240717", "1"}, {"ZHANG_R51_ENABLE", "1"},
        {"ZHANG_R51_RATIO_ONLY", "1"}, {"ZHANG_R51_AR_START_GPST_SECONDS", "1405216800"},
        {"ZHANG_R49_FUSION_WEIGHT", "0.10"}, {"OPENBLAS_NUM_THREADS", "4"},
        {"OMP_NUM_THREADS", "4"}, {"MKL_NUM_THREADS", "1"}};
    std::string hashes[]{
        "9e5ea901db26708ecb2f772b7ad5efde886d919491811a536302a915ec776c3d",
        "4463f3a4c97f0e213c15b995d9e215a873178428a0cf7bcab7e436c98f96b80a",
        "6ee5be509848128c71ff23e267b9e10df839c746ec06cd363fe3462ccbbbe841",
        "4f6351a7db3a2af79ae728ffc1d8d6e9b19986cafe49d6c47c6a703cfeac4874"};
    const auto allowed = [&]() {
        return zhangR51AuditedRatioReuseAllowed(hashes[0], hashes[1], hashes[2], hashes[3], environment);
    };
    require(allowed(), "pinned ratio-only fixture rejected");
    for (auto& hash : hashes)
    {
        const auto original = hash;
        hash = "WRONG";
        require(!allowed(), "wrong ratio-only source hash accepted");
        hash = original;
    }
    const auto originalEnvironment = environment;
    for (const auto& [name, value] : originalEnvironment)
    {
        environment[name] = "WRONG";
        require(!allowed(), "wrong ratio-only environment accepted");
        environment = originalEnvironment;
        environment.erase(name);
        require(!allowed(), "missing ratio-only environment accepted");
        environment = originalEnvironment;
    }
    environment["ZHANG_R51_RATIO_THRESHOLD"] = "2";
    environment["ZHANG_R51_EFFECTIVE_RATIO_THRESHOLD"] = "2";
    require(!allowed(), "fixed ratio=3 checkpoint reused under ratio=2 as same-policy migration");
    environment["ZHANG_R51_RATIO_THRESHOLD"] = "3";
    require(!allowed(), "recorded effective ratio mismatch accepted");
    environment["ZHANG_R51_EFFECTIVE_RATIO_THRESHOLD"] = "3";
    require(allowed(), "explicit unchanged ratio=3 policy rejected");
    environment["ZHANG_R51_RATIO_THRESHOLD"] = "nan";
    require(!allowed(), "invalid ratio policy accepted");
    std::cout << "PASS pinned ratio-only checkpoint, source mismatch rejection and ratio-policy isolation\n";
}
