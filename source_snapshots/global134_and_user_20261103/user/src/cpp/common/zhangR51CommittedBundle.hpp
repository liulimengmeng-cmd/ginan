#include "common/zhangRealtimeIntegrity.hpp"
#include "common/zhangRealtimeHandoffProof.hpp"
#pragma once

#include "common/zhangCheckpoint.hpp"
#include "common/zhangIntegerDecisionProof.hpp"
#include "common/zhangMultiGnssProductManifest.hpp"

#include <boost/json.hpp>

#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

/** A verified epoch directory is the only atomic R51 publication unit.
 * Products and covariance must be consumed from these paired paths. */
struct ZhangR51VerifiedBundle
{
    bool valid = false;
    bool ar = false;
    std::filesystem::path directory;
    std::filesystem::path products;
    std::filesystem::path covariance;
    std::filesystem::path certificate;
    std::string epoch;
    std::string route;
    std::string root;
    std::string failureReason = "NOT_EVALUATED";
};

inline std::string zhangR51BundleFileHash(const std::filesystem::path& path)
{
    std::string reason;
    const std::string hash =
        zhangCheckpointFileSha256(path.string(), &reason);
    if (hash.empty())
        throw std::runtime_error("R51_BUNDLE_HASH_FAILED:" + reason);
    return hash;
}

inline ZhangR51VerifiedBundle zhangR51VerifyBundleContents(
    const std::filesystem::path& directory, bool requireAr,
    bool allowPending)
{
    ZhangR51VerifiedBundle result;
    result.directory = directory;
    result.products = directory / "products.csv";
    result.covariance = directory / "covariance.csv";
    result.certificate = directory / "certificate.txt";
    auto fail = [&](const char* reason)
    {
        result.failureReason = reason;
        return result;
    };
    try
    {
        if (!std::filesystem::is_directory(directory))
            return fail("R51_BUNDLE_DIRECTORY_MISSING");
        if (!allowPending &&
            directory.filename().string().ends_with(".pending"))
            return fail("R51_BUNDLE_NOT_COMMITTED");
        const auto manifestPath = directory / "COMMIT.json";
        const auto sumsPath = directory / "SHA256SUMS";
        for (const auto& path : {result.products, result.certificate,
                 manifestPath, sumsPath})
            if (!std::filesystem::is_regular_file(path))
                return fail("R51_BUNDLE_REQUIRED_FILE_MISSING");

        std::ifstream manifestInput(manifestPath, std::ios::binary);
        std::ostringstream manifestText;
        manifestText << manifestInput.rdbuf();
        if (!manifestInput || manifestInput.bad())
            return fail("R51_BUNDLE_MANIFEST_READ_FAILED");
        const auto parsed = boost::json::parse(manifestText.str());
        if (!parsed.is_object())
            return fail("R51_BUNDLE_MANIFEST_INVALID");
        const auto& manifest = parsed.as_object();
        const auto& epoch = manifest.at("epoch");
        const auto& route = manifest.at("route");
        const auto& root = manifest.at("root");
        const auto& ar = manifest.at("ar");
        const auto& productBytes = manifest.at("products_bytes");
        const auto& covarianceBytes = manifest.at("covariance_bytes");
        if (!epoch.is_string() || !route.is_string() || !root.is_string() ||
            !ar.is_bool() || !productBytes.is_int64() ||
            !covarianceBytes.is_int64())
            return fail("R51_BUNDLE_MANIFEST_INVALID");
        result.epoch = epoch.as_string().c_str();
        result.route = route.as_string().c_str();
        result.root = root.as_string().c_str();
        result.ar = ar.as_bool();
        if (result.epoch.empty() || result.route.empty() ||
            result.root.empty() || productBytes.as_int64() <= 0 ||
            covarianceBytes.as_int64() < 0 ||
            (requireAr && !result.ar))
            return fail("R51_BUNDLE_MANIFEST_INVALID");
        if (result.ar && covarianceBytes.as_int64() <= 0)
            return fail("R51_BUNDLE_COVARIANCE_MISSING");
        if (std::filesystem::file_size(result.products) !=
            static_cast<std::uintmax_t>(productBytes.as_int64()))
            return fail("R51_BUNDLE_PRODUCT_SIZE_MISMATCH");
        const bool covarianceExists =
            std::filesystem::is_regular_file(result.covariance);
        if ((result.ar || requireAr) && !covarianceExists)
            return fail("R51_BUNDLE_COVARIANCE_MISSING");
        if (covarianceExists &&
            std::filesystem::file_size(result.covariance) !=
                static_cast<std::uintmax_t>(covarianceBytes.as_int64()))
            return fail("R51_BUNDLE_COVARIANCE_SIZE_MISMATCH");
        if (!covarianceExists && covarianceBytes.as_int64() != 0)
            return fail("R51_BUNDLE_COVARIANCE_SIZE_MISMATCH");

        std::ifstream sumsInput(sumsPath);
        if (!sumsInput) return fail("R51_BUNDLE_HASH_MANIFEST_READ_FAILED");
        std::map<std::string, std::string> declared;
        std::string line;
        while (std::getline(sumsInput, line))
        {
            if (line.size() < 67 || line.substr(64, 2) != "  ")
                return fail("R51_BUNDLE_HASH_MANIFEST_INVALID");
            const std::string name = line.substr(66);
            if (name != "products.csv" && name != "covariance.csv" &&
                name != "certificate.txt" && name != "COMMIT.json")
                return fail("R51_BUNDLE_HASH_MANIFEST_INVALID");
            if (!declared.emplace(name, line.substr(0, 64)).second)
                return fail("R51_BUNDLE_HASH_MANIFEST_INVALID");
        }
        if (sumsInput.bad())
            return fail("R51_BUNDLE_HASH_MANIFEST_READ_FAILED");
        if (declared.size() != (covarianceExists ? 4U : 3U))
            return fail("R51_BUNDLE_HASH_MANIFEST_INCOMPLETE");
        for (const auto& name : {"products.csv", "certificate.txt",
                 "COMMIT.json"})
            if (!declared.contains(name) ||
                declared.at(name) != zhangR51BundleFileHash(directory / name))
                return fail("R51_BUNDLE_HASH_MISMATCH");
        if (covarianceExists &&
            (!declared.contains("covariance.csv") ||
             declared.at("covariance.csv") !=
                zhangR51BundleFileHash(result.covariance)))
            return fail("R51_BUNDLE_HASH_MISMATCH");
        // Integrity hashes prove file identity, not that every dictionary
        // reference resolves to the same original proof statement. Validate
        // the optional INTERNED_V1 wire semantics after the existing checks;
        // old inline proof_detail certificates remain supported.
        std::ifstream certificateInput(result.certificate, std::ios::binary);
        ZhangHandoffProofCertificate handoffProofs;
        std::string proofReason;
        if (!certificateInput || !zhangReadHandoffProofCertificate(
                certificateInput, handoffProofs, proofReason))
        {
            result.failureReason = "R51_BUNDLE_HANDOFF_PROOF_INVALID:" + proofReason;
            return result;
        }
        if (const auto schema = manifest.if_contains("schema_version"))
        {
            if (!schema->is_int64() || schema->as_int64() != 2 ||
                !manifest.contains("systems") || !manifest.at("systems").is_object())
                return fail("R51_BUNDLE_SYSTEM_SCHEMA_INVALID");
            std::ifstream productRows(result.products);
            std::ifstream certificateRows(result.certificate, std::ios::binary);
            std::ostringstream certificateBytes;
            certificateBytes << certificateRows.rdbuf();
            std::string systemReason;
            if (!productRows || !certificateRows || !zhangMultiGnssValidateSystemManifest(
                    productRows, certificateBytes.str(), manifest.at("systems").as_object(),
                    result.ar, systemReason))
            {
                result.failureReason = "R51_BUNDLE_SYSTEM_AUTHORITY_INVALID:" + systemReason;
                return result;
            }
        }
        else
        {
            std::ifstream productRows(result.products);
            if (!productRows || manifest.contains("systems") ||
                zhangMultiGnssProductSystems(productRows).size() > 1)
                return fail("R51_BUNDLE_MULTI_SYSTEM_SCHEMA_REQUIRED");
        }
        zhangVerifyRealtimeNavigationManifest(manifest, directory);
        result.valid = true;
        result.failureReason = "NONE";
        return result;
    }
    catch (...)
    {
        return fail("R51_BUNDLE_VERIFICATION_EXCEPTION");
    }
}

/** Consumer API. A `.pending` directory never authorizes user products. */
inline ZhangR51VerifiedBundle zhangR51VerifyCommittedBundle(
    const std::filesystem::path& directory, bool requireAr)
{
    return zhangR51VerifyBundleContents(directory, requireAr, false);
}
