#include "common/zhangMultiGnssProductManifest.hpp"
#include "common/zhangMultiGnssSelection.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool pass, const char* reason)
{ if (!pass) throw std::runtime_error(reason); }

int main()
{
    try
    {
        ACSConfig config;
        require(zhangMultiGnssSatelliteAllowed(config, SatSys("C01")), "legacy selection changed");
        config.zhangFullRank.bds3_only = true;
        require(!zhangMultiGnssSatelliteAllowed(config, SatSys("C19")), "empty whitelist not fail-closed");
        config.zhangFullRank.bds3_satellite_whitelist = {"C19", "C45"};
        require(zhangMultiGnssSatelliteAllowed(config, SatSys("C19")), "approved satellite rejected");
        require(!zhangMultiGnssSatelliteAllowed(config, SatSys("C20")), "PRN threshold substituted for whitelist");
        require(!zhangMultiGnssSatelliteAllowed(config, SatSys("C01")), "unlisted BDS accepted");
        require(zhangMultiGnssSatelliteAllowed(config, SatSys("E01")), "BDS selection crossed system boundary");

        const ZhangMultiGnssManifestSignals signals{
            {"GPS", {"L1C", "L2W"}}, {"GAL", {"L1C", "L5Q"}},
            {"BDS", {"L1P", "L5P"}}};
        const std::string header = "solution,satellite,observable,pppar_usable,integer_valid,dual_frequency_ar_valid,numeric_valid,branch_valid,datum_version,integer_component_version,public_phase_segment_id,invalid_reason\n";
        const std::string body =
            "PRODUCT_FIXED,G01,L1C,1,1,1,1,1,0,1,G01-L1C,NONE\n"
            "PRODUCT_FIXED,G01,L2W,1,1,1,1,1,0,1,G01-L2W,NONE\n"
            "PRODUCT_FIXED,E01,L1C,0,0,0,1,1,2,0,E01-L1C,NO_CURRENT_CERTIFICATE\n"
            "PRODUCT_FIXED,E01,L5Q,0,0,0,1,1,2,0,E01-L5Q,NO_CURRENT_CERTIFICATE\n"
            "PRODUCT_FIXED,C19,L1P,1,1,1,1,1,3,4,C19-L1P,NONE\n"
            "PRODUCT_FIXED,C19,L5P,1,1,1,1,1,3,4,C19-L5P,NONE\n";
        const std::string certificate =
            "multi_system_certificate_begin system=GPS\nGPS independent proof\n"
            "multi_system_certificate_end system=GPS\n"
            "multi_system_certificate_begin system=BDS\nBDS independent proof\n"
            "multi_system_certificate_end system=BDS\n";
        std::istringstream source(header + body);
        auto manifest = zhangMultiGnssBuildSystemManifest(source, certificate, signals);
        // Independent expected counts/statuses: E01 shares L1C and PRN 01
        // with G01, but may not acquire any of GPS's authority.
        require(manifest.size() == 3, "system disappeared from manifest");
        for (const char* system : {"GPS", "GAL", "BDS"})
            require(manifest.at(system).as_object().at("fixed_phase_rows").as_int64() == 2,
                "independent per-system row count wrong");
        require(manifest.at("GPS").as_object().at("pppar_usable").as_bool(), "GPS certificate lost");
        require(!manifest.at("GAL").as_object().at("any_pppar_usable").as_bool(), "GAL inherited GPS L1C authority");
        require(manifest.at("GAL").as_object().at("system_certificate_sha256").is_null(), "missing GAL proof fabricated");
        require(manifest.at("BDS").as_object().at("pppar_usable").as_bool(), "unrelated BDS lost authority");
        require(manifest.at("GPS").as_object().at("system_certificate_sha256") !=
            manifest.at("BDS").as_object().at("system_certificate_sha256"), "certificate identity mixed");
        std::string reason;
        std::istringstream valid(header + body);
        require(zhangMultiGnssValidateSystemManifest(valid, certificate, manifest, true, reason), reason.c_str());
        auto forged = manifest;
        forged.at("GAL").as_object()["any_pppar_usable"] = true;
        std::istringstream forgedRows(header + body);
        require(!zhangMultiGnssValidateSystemManifest(forgedRows, certificate, forged, true, reason), "forged GAL status accepted");
        std::istringstream wrongTopFlag(header + body);
        require(!zhangMultiGnssValidateSystemManifest(wrongTopFlag, certificate, manifest, false, reason), "top-level AR mismatch accepted");
        auto mustReject = [&](std::string csv, std::string proof)
        {
            bool rejected = false;
            try { std::istringstream stream(csv); zhangMultiGnssBuildSystemManifest(stream, proof, signals); }
            catch (const std::exception&) { rejected = true; }
            require(rejected, "malformed multi-system product/proof accepted");
        };
        // The writer may retain a single-frequency PPP-only degradation
        // while unrelated satellites/systems still publish dual-frequency AR.
        auto partialBody = body;
        const std::string absent = "PRODUCT_FIXED,E01,L5Q,0,0,0,1,1,2,0,E01-L5Q,NO_CURRENT_CERTIFICATE\n";
        partialBody.erase(partialBody.find(absent), absent.size());
        std::istringstream partialSource(header + partialBody);
        auto partialManifest = zhangMultiGnssBuildSystemManifest(partialSource, certificate, signals);
        require(partialManifest.at("GAL").as_object().at("fixed_phase_rows").as_int64() == 1,
            "missing frequency was synthesized");
        require(!partialManifest.at("GAL").as_object().at("any_pppar_usable").as_bool() &&
            partialManifest.at("GPS").as_object().at("pppar_usable").as_bool() &&
            partialManifest.at("BDS").as_object().at("pppar_usable").as_bool(),
            "PPP-only missing band discarded or authorized unrelated AR");
        bool recordedMissing = false;
        for (const auto& value : partialManifest.at("GAL").as_object().at("failure_reasons").as_array())
            recordedMissing |= value.as_string() == "FREQUENCY_PARTNER_MISSING";
        require(recordedMissing, "single-frequency degradation was not recorded");
        std::istringstream partialValidation(header + partialBody);
        require(zhangMultiGnssValidateSystemManifest(partialValidation, certificate, partialManifest, true, reason),
            "valid partial PPP-only manifest rejected by consumer");
        const std::string normal = "PRODUCT_FIXED,E01,L1C,0,0,0,1,1,2,0,E01-L1C,NO_CURRENT_CERTIFICATE\n";
        auto assertedDual = partialBody;
        assertedDual.replace(assertedDual.find(normal), normal.size(),
            "PRODUCT_FIXED,E01,L1C,0,0,1,1,1,2,0,E01-L1C,NO_CURRENT_CERTIFICATE\n");
        mustReject(header + assertedDual, certificate);
        auto assertedAr = partialBody;
        assertedAr.replace(assertedAr.find(normal), normal.size(),
            "PRODUCT_FIXED,E01,L1C,1,1,1,1,1,2,0,E01-L1C,NONE\n");
        mustReject(header + assertedAr, certificate);
        auto asymmetricPair = body;
        asymmetricPair.replace(asymmetricPair.find(normal), normal.size(),
            "PRODUCT_FIXED,E01,L1C,1,1,1,1,1,2,0,E01-L1C,NONE\n");
        mustReject(header + asymmetricPair, certificate);
        auto withArColumn = [](const std::string& csv, bool forgeGalileo)
        {
            std::istringstream input(csv);
            std::ostringstream output;
            std::string line;
            std::getline(input, line); output << line << ",ar_valid\n";
            while (std::getline(input, line))
            {
                const auto fields = zhangMultiGnssCsvFields(line);
                output << line << ',' << (forgeGalileo && fields[1] == "E01" ? "1" : fields[3]) << '\n';
            }
            return output.str();
        };
        std::istringstream explicitArSource(withArColumn(header + partialBody, false));
        zhangMultiGnssBuildSystemManifest(explicitArSource, certificate, signals);
        mustReject(withArColumn(header + partialBody, true), certificate);
        mustReject(header + body + "PRODUCT_FIXED,G01,L1C,1,1,1,1,1,0,1,G01-L1C,NONE\n", certificate);
        auto wrongTrackingCode = header + body;
        wrongTrackingCode.replace(wrongTrackingCode.find("L5Q"), 3, "L5X");
        mustReject(wrongTrackingCode, certificate);
        mustReject(header + body, "no per-system certificate sections\n");
        mustReject(header + body, certificate + "multi_system_certificate_begin system=GPS\n");
        std::cout << "multi_gnss_product_manifest: PASS (independent G/E/C flags, hashes, explicit selection, failure isolation)\n";
        return 0;
    }
    catch (const std::exception& error)
    { std::cerr << error.what() << '\n'; return 1; }
}
