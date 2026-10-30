#include "common/zhangR51CommittedBundle.hpp"

#include <chrono>
#include <iostream>

static void require(bool condition, const char* why)
{
    if (!condition) throw std::runtime_error(why);
}

int main()
{
    namespace fs = std::filesystem;
    const auto directory = fs::temp_directory_path() / ("r51-handoff-proof-bundle-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    require(fs::create_directory(directory), "unique test directory could not be created");
    auto write = [](const fs::path& path, const std::string& bytes)
    {
        std::ofstream output(path, std::ios::binary);
        output.write(bytes.data(), bytes.size());
        require(bool(output), "fixture write failed");
    };
    write(directory / "products.csv", "solution,satellite\nPRODUCT_FIXED,G02\n");
    write(directory / "covariance.csv", "covariance\n1\n");
    boost::json::object manifest;
    manifest["epoch"] = "2019-01-01 02:00:00";
    manifest["route"] = "CURRENT_EPOCH_SINGLE_BASIS";
    manifest["root"] = "source-root";
    manifest["ar"] = true;
    manifest["products_bytes"] = std::int64_t(fs::file_size(directory / "products.csv"));
    manifest["covariance_bytes"] = std::int64_t(fs::file_size(directory / "covariance.csv"));
    write(directory / "COMMIT.json", boost::json::serialize(manifest));
    auto checksums = [&]()
    {
        std::ostringstream sums;
        for (const auto& name : {"products.csv", "covariance.csv", "certificate.txt", "COMMIT.json"})
            sums << zhangR51BundleFileHash(directory / name) << "  " << name << '\n';
        write(directory / "SHA256SUMS", sums.str());
    };
    std::string reason;
    ZhangDecisionProofDictionary dictionary;
    const std::string payload = "id=\"accepted-1\" source=\"epoch\" statement=\"a=1\nfull witness\"";
    require(dictionary.intern("accepted-1", payload, reason), "fixture intern failed");
    auto certificate = [&](const std::string& id)
    {
        std::ostringstream text;
        text << "handoff_proof_format=INTERNED_V1\n"
            << "handoff_satellite=G02 proof_count=1 proof=" << std::quoted(id)
            << " proof_detail_ref=" << std::quoted(id) << '\n'
            << "handoff_proof_dictionary_begin\n";
        dictionary.write(text);
        text << "handoff_proof_dictionary_end\n";
        return text.str();
    };
    write(directory / "certificate.txt", certificate("accepted-1"));
    checksums();
    require(zhangR51VerifyCommittedBundle(directory, true).valid,
        "valid interned bundle failed production committed reader");
    write(directory / "certificate.txt", certificate("unknown"));
    checksums(); // Hashes are correct: the rejection must be semantic.
    const auto missing = zhangR51VerifyCommittedBundle(directory, true);
    require(!missing.valid && missing.failureReason ==
        "R51_BUNDLE_HANDOFF_PROOF_INVALID:INTERNED_PROOF_REFERENCE_MISSING",
        "correctly hashed unknown proof reference accepted");

    // A conflicting record cannot be hidden behind an otherwise valid SHA.
    std::ostringstream conflict;
    conflict << "handoff_proof_format=INTERNED_V1\n"
        << "handoff_satellite=G02 proof_count=1 proof=\"accepted-1\" "
        << "proof_detail_ref=\"accepted-1\"\n"
        << "handoff_proof_dictionary_begin\nINTERNED_V1\n2\n";
    auto field = [&](const std::string& bytes)
    {
        conflict << bytes.size() << '\n';
        conflict.write(bytes.data(), bytes.size());
        conflict.put('\n');
    };
    field("accepted-1"); field(payload);
    field("accepted-1"); field(payload + "changed content");
    conflict << "handoff_proof_dictionary_end\n";
    write(directory / "certificate.txt", conflict.str());
    checksums();
    const auto collision = zhangR51VerifyCommittedBundle(directory, true);
    require(!collision.valid && collision.failureReason ==
        "R51_BUNDLE_HANDOFF_PROOF_INVALID:INTERNED_PROOF_ID_CONTENT_CONFLICT",
        "correctly hashed same-ID different-content dictionary accepted");

    std::ostringstream legacy;
    legacy << "handoff_satellite=G02 proof_count=1 proof=\"accepted-1\" proof_detail="
        << std::quoted(payload) << '\n';
    write(directory / "certificate.txt", legacy.str());
    checksums();
    require(zhangR51VerifyCommittedBundle(directory, true).valid,
        "legacy inline proof bundle compatibility failed");
    write(directory / "products.csv", "solution,satellite\nPRODUCT_FIXED,G02\nPRODUCT_FIXED,E11\n");
    manifest["products_bytes"] = std::int64_t(fs::file_size(directory / "products.csv"));
    write(directory / "COMMIT.json", boost::json::serialize(manifest));
    checksums();
    const auto downgraded = zhangR51VerifyCommittedBundle(directory, true);
    require(!downgraded.valid && downgraded.failureReason ==
        "R51_BUNDLE_MULTI_SYSTEM_SCHEMA_REQUIRED",
        "multi-system bundle bypassed system authority via legacy schema");
    write(directory / "products.csv", "solution,satellite\nPRODUCT_FIXED,G02\n");
    manifest["products_bytes"] = std::int64_t(fs::file_size(directory / "products.csv"));
    write(directory / "COMMIT.json", boost::json::serialize(manifest));
    checksums();
    write(directory / "certificate.txt", legacy.str() + "same hash must fail\n");
    const auto corrupted = zhangR51VerifyCommittedBundle(directory, true);
    require(!corrupted.valid && corrupted.failureReason == "R51_BUNDLE_HASH_MISMATCH",
        "proof parser weakened pre-existing SHA identity validation");
    for (const auto& name : {"products.csv", "covariance.csv", "certificate.txt", "COMMIT.json", "SHA256SUMS"})
        fs::remove(directory / name);
    fs::remove(directory);
    std::cout << "PASS handoff_bundle_reader real_sha_legacy_new_format_and_unknown_rejection\n";
}
