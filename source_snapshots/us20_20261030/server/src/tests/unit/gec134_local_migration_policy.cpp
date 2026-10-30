// Pure migration-policy test: no PEA, KF, checkpoint importer or target
// canonical builder is linked. Fixtures come independently from the real
// immutable source canonical bytes and explicit target changes.
#include "zhangGec134Migration.hpp"
#include <openssl/evp.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

std::string zhangCheckpointSha256(const std::string& input)
{
    unsigned char output[EVP_MAX_MD_SIZE]; unsigned int size = 0;
    if (!EVP_Digest(input.data(), input.size(), output, &size, EVP_sha256(), nullptr) || size != 32)
        throw std::runtime_error("Independent standard SHA-256 failed");
    std::ostringstream text;
    for (unsigned int i = 0; i < size; ++i)
        text << std::hex << std::setw(2) << std::setfill('0') << unsigned(output[i]);
    return text.str();
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2) throw std::runtime_error("Expected the independent fixture JSON path");
        const auto fixture = YAML::LoadFile(argv[1]);
        std::map<std::string, std::string> base;
        for (const auto& pair : fixture["environment"])
            base.emplace(pair.first.as<std::string>(), pair.second.as<std::string>());
        const auto initial = fixture["cases"][0]["target"].as<std::string>();
        std::size_t checks = 0, localChecks = 0, serverChecks = 0, positives = 0;
        for (const auto& item : fixture["cases"])
        {
            const auto name = item["name"].as<std::string>();
            auto read = [&](const char* key) {
                return item[key] ? item[key].as<std::string>() : fixture[key].as<std::string>();
            };
            auto env = base;
            if (item["env_overrides"]) for (const auto& pair : item["env_overrides"])
                env[pair.first.as<std::string>()] = pair.second.as<std::string>();
            if (item["remove_env"]) for (const auto& key : item["remove_env"])
                env.erase(key.as<std::string>());
            // Warm with an accepted SERVER token first. Every rejection must
            // invalidate it and clear its target text, not leave stale authority.
            ZhangGec134MigrationAuthority authority; std::string reason;
            if (!validateZhangGec134Migration(fixture["source_bundle_sha"].as<std::string>(),
                    fixture["source_binary_sha"].as<std::string>(), fixture["source_config_sha"].as<std::string>(),
                    fixture["source_input_sha"].as<std::string>(), fixture["source_config_text"].as<std::string>(),
                    initial, fixture["target_input_sha"].as<std::string>(), base, authority, reason))
                throw std::runtime_error("Unchanged server positive fixture rejected: " + reason);
            const auto target = item["target"].as<std::string>();
            const bool got = validateZhangGec134Migration(read("source_bundle_sha"), read("source_binary_sha"),
                read("source_config_sha"), read("source_input_sha"), read("source_config_text"), target,
                read("target_input_sha"), env, authority, reason);
            const bool expected = item["expected"].as<bool>();
            if (got != expected || authority.valid() != expected ||
                (expected && (authority.targetConfigText() != target || !reason.empty())) ||
                (!expected && (reason.empty() || !authority.targetConfigText().empty())))
                throw std::runtime_error("Policy case failed: " + name + " reason=" + reason);
            if (expected && !authority.matchesSource(ZHANG_GEC134_SOURCE_MODE,
                    ZHANG_GEC134_SOURCE_BINARY_SHA256, ZHANG_GEC134_SOURCE_CONFIG_SHA256,
                    ZHANG_GEC134_SOURCE_INPUT_SHA256))
                throw std::runtime_error("Accepted token lost immutable source binding");
            ++checks;
            if (item["profile"].as<std::string>() == "local") ++localChecks; else ++serverChecks;
            if (expected) ++positives;
            std::cout << "PASS " << name << " accepted=" << got << " reason=" << reason << '\n';
        }
        std::cout << "GEC134_LOCAL_POLICY_PASSED checks=" << checks << " positive=" << positives
                  << " server_regression=" << serverChecks << " local=" << localChecks
                  << " independent_sha256=1 exact_target_token=1 no_PEA_or_KF=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GEC134_LOCAL_POLICY_FAIL " << error.what() << '\n';
        return 1;
    }
}
