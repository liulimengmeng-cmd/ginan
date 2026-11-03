#pragma once

#include "common/zhangCheckpoint.hpp"
#include "common/satSys.hpp"
#include <boost/json.hpp>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using ZhangMultiGnssManifestSignals =
    std::map<std::string, std::vector<std::string>>;

inline std::vector<std::string> zhangMultiGnssCsvFields(const std::string& line)
{
    std::vector<std::string> fields;
    std::string value;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];
        if (c == '"')
        {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"')
            { value += '"'; ++i; }
            else quoted = !quoted;
        }
        else if (c == ',' && !quoted) { fields.push_back(value); value.clear(); }
        else if (c != '\r') value += c;
    }
    if (quoted) throw std::runtime_error("MULTI_GNSS_CSV_UNTERMINATED_QUOTE");
    fields.push_back(value);
    return fields;
}

inline std::string zhangMultiGnssCertificateSection(
    const std::string& certificate, const std::string& system)
{
    const std::string begin = "multi_system_certificate_begin system=" + system + "\n";
    const std::string end = "multi_system_certificate_end system=" + system + "\n";
    const auto first = certificate.find(begin);
    if (first == std::string::npos) return {};
    if ((first != 0 && certificate[first - 1] != '\n') ||
        certificate.find(begin, first + begin.size()) != std::string::npos)
        throw std::runtime_error("MULTI_GNSS_CERTIFICATE_SECTION_AMBIGUOUS");
    const auto last = certificate.find(end, first + begin.size());
    if (last == std::string::npos ||
        certificate.find(end, last + end.size()) != std::string::npos)
        throw std::runtime_error("MULTI_GNSS_CERTIFICATE_SECTION_UNCLOSED");
    return certificate.substr(first, last + end.size() - first);
}

// A missing schema is a legacy single-system exemption, not permission to
// silently drop the per-system authority of a multi-system product bundle.
inline std::set<E_Sys> zhangMultiGnssProductSystems(std::istream& products)
{
    std::string line;
    if (!std::getline(products, line))
        throw std::runtime_error("MULTI_GNSS_PRODUCT_HEADER_MISSING");
    const auto header = zhangMultiGnssCsvFields(line);
    const auto column = std::find(header.begin(), header.end(), "satellite");
    if (column == header.end() || std::count(header.begin(), header.end(), "satellite") != 1)
        throw std::runtime_error("MULTI_GNSS_PRODUCT_SATELLITE_COLUMN_INVALID");
    const auto index = static_cast<std::size_t>(column - header.begin());
    std::set<E_Sys> systems;
    while (std::getline(products, line))
    {
        if (line.empty()) continue;
        const auto fields = zhangMultiGnssCsvFields(line);
        if (fields.size() != header.size())
            throw std::runtime_error("MULTI_GNSS_PRODUCT_ROW_WIDTH_INVALID");
        const SatSys satellite(fields[index].c_str());
        if (satellite.prn <= 0 || satellite.id() != fields[index])
            throw std::runtime_error("MULTI_GNSS_PRODUCT_SATELLITE_INVALID");
        systems.insert(satellite.sys);
    }
    if (products.bad()) throw std::runtime_error("MULTI_GNSS_PRODUCT_READ_FAILED");
    return systems;
}

/** Build system metadata from the actual final PRODUCT_FIXED rows, not the
 * presence of a network-wide AR flag. The source certificate sections remain
 * independently bound, while covariance is still one complete joint matrix. */
inline boost::json::object zhangMultiGnssBuildSystemManifest(
    std::istream& products, const std::string& certificate,
    const ZhangMultiGnssManifestSignals& configuredSignals)
{
    struct Summary
    {
        std::int64_t rows = 0, usable = 0;
        std::map<std::string, std::set<std::string>> satellites;
        std::map<std::string, std::set<std::string>> usableSignals;
        std::set<std::string> satellitesClaimingAr;
        std::set<std::string> reasons, datumVersions, integerVersions, segments;
    };
    std::map<std::string, Summary> summaries;
    for (const auto& [system, signals] : configuredSignals)
    {
        if (signals.size() != 2 || signals[0] == signals[1])
            throw std::runtime_error("MULTI_GNSS_MANIFEST_BASELINE_INVALID");
        summaries[system];
    }
    std::string line;
    if (!std::getline(products, line))
        throw std::runtime_error("MULTI_GNSS_PRODUCT_HEADER_MISSING");
    const auto header = zhangMultiGnssCsvFields(line);
    std::map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < header.size(); ++i)
        if (!columns.emplace(header[i], i).second)
            throw std::runtime_error("MULTI_GNSS_PRODUCT_DUPLICATE_COLUMN");
    for (const char* name : {"solution", "satellite", "observable", "pppar_usable",
            "integer_valid", "dual_frequency_ar_valid", "numeric_valid", "branch_valid",
            "datum_version", "integer_component_version", "public_phase_segment_id", "invalid_reason"})
        if (!columns.contains(name))
            throw std::runtime_error(std::string("MULTI_GNSS_PRODUCT_COLUMN_MISSING:") + name);
    while (std::getline(products, line))
    {
        if (line.empty()) continue;
        const auto fields = zhangMultiGnssCsvFields(line);
        if (fields.size() != header.size())
            throw std::runtime_error("MULTI_GNSS_PRODUCT_ROW_WIDTH_INVALID");
        auto field = [&](const char* name) -> const std::string& { return fields.at(columns.at(name)); };
        if (field("solution") != "PRODUCT_FIXED") continue;
        const SatSys satellite(field("satellite").c_str());
        if (satellite.prn <= 0 || satellite.id() != field("satellite"))
            throw std::runtime_error("MULTI_GNSS_PRODUCT_SATELLITE_INVALID");
        const std::string system = enum_to_string(satellite.sys);
        const auto configured = configuredSignals.find(system);
        if (configured == configuredSignals.end())
            throw std::runtime_error("MULTI_GNSS_PRODUCT_UNDECLARED_SYSTEM");
        const auto& code = field("observable");
        if (code != configured->second[0] && code != configured->second[1])
            throw std::runtime_error("MULTI_GNSS_PRODUCT_UNDECLARED_SIGNAL");
        auto& summary = summaries[system];
        if (!summary.satellites[field("satellite")].insert(code).second)
            throw std::runtime_error("MULTI_GNSS_PRODUCT_DUPLICATE_SIGNAL");
        ++summary.rows;
        auto flag = [&](const char* name)
        {
            if (field(name) != "0" && field(name) != "1")
                throw std::runtime_error("MULTI_GNSS_PRODUCT_FLAG_INVALID");
            return field(name) == "1";
        };
        const bool usable = flag("pppar_usable");
        const bool dualAr = flag("dual_frequency_ar_valid");
        // ar_valid is present in production CSV65, but absent from supported
        // minimal legacy test/metadata fixtures. It cannot override pppar.
        const bool ar = columns.contains("ar_valid") ? flag("ar_valid") : usable;
        if (ar != usable)
            throw std::runtime_error("MULTI_GNSS_PRODUCT_AR_AUTHORITY_INCONSISTENT");
        if (usable || dualAr || ar)
            summary.satellitesClaimingAr.insert(field("satellite"));
        if (usable)
        {
            for (const char* name : {"integer_valid", "dual_frequency_ar_valid", "numeric_valid", "branch_valid"})
                if (field(name) != "1")
                    throw std::runtime_error("MULTI_GNSS_PRODUCT_AR_AUTHORITY_INCONSISTENT");
            ++summary.usable;
            summary.usableSignals[field("satellite")].insert(code);
        }
        else summary.reasons.insert(field("invalid_reason"));
        summary.datumVersions.insert(field("datum_version"));
        summary.integerVersions.insert(field("integer_component_version"));
        summary.segments.insert(field("public_phase_segment_id"));
    }
    boost::json::object systems;
    auto array = [](const auto& values)
    {
        boost::json::array result;
        for (const auto& value : values) result.emplace_back(value);
        return result;
    };
    const auto epochHash = zhangCheckpointSha256(certificate);
    for (const auto& [system, summary] : summaries)
    {
        bool missingPartner = false;
        for (const auto& [satellite, codes] : summary.satellites)
        {
            if (codes.size() != 2 && summary.satellitesClaimingAr.contains(satellite))
                throw std::runtime_error("MULTI_GNSS_PRODUCT_FREQUENCY_PARTNER_MISSING");
            // A real single-frequency PPP-only row is a legal degradation,
            // not a reason to discard unrelated G/E/C products. No missing
            // phase row/value/covariance is synthesized here.
            missingPartner |= codes.size() != 2;
            const auto usable = summary.usableSignals.find(satellite);
            if (usable != summary.usableSignals.end() && usable->second.size() != 2)
                throw std::runtime_error("MULTI_GNSS_PRODUCT_ASYMMETRIC_DUAL_AR_AUTHORITY");
        }
        const auto section = zhangMultiGnssCertificateSection(certificate, system);
        if (summary.usable > 0 && section.empty())
            throw std::runtime_error("MULTI_GNSS_PRODUCT_SYSTEM_CERTIFICATE_MISSING");
        boost::json::object item;
        item["signals"] = array(configuredSignals.at(system));
        item["fixed_phase_rows"] = summary.rows;
        item["pppar_usable_phase_rows"] = summary.usable;
        item["any_pppar_usable"] = summary.usable > 0;
        item["pppar_usable"] = summary.rows > 0 && summary.usable == summary.rows;
        item["datum_versions"] = array(summary.datumVersions);
        item["integer_component_versions"] = array(summary.integerVersions);
        item["public_phase_segments"] = array(summary.segments);
        auto reasons = summary.reasons;
        if (missingPartner) reasons.insert("FREQUENCY_PARTNER_MISSING");
        if (summary.rows == 0) reasons.insert("NO_FIXED_PRODUCT_ROWS");
        item["failure_reasons"] = array(reasons);
        item["epoch_certificate_sha256"] = epochHash;
        if (section.empty()) item["system_certificate_sha256"] = nullptr;
        else item["system_certificate_sha256"] = zhangCheckpointSha256(section);
        systems[system] = std::move(item);
    }
    return systems;
}

inline bool zhangMultiGnssValidateSystemManifest(
    std::istream& products, const std::string& certificate,
    const boost::json::object& declared, bool epochAr, std::string& reason)
{
    try
    {
        ZhangMultiGnssManifestSignals configured;
        for (const auto& entry : declared)
        {
            if (!entry.value().is_object()) throw std::runtime_error("MULTI_GNSS_SYSTEM_STATUS_INVALID");
            const auto& item = entry.value().as_object();
            const auto& signals = item.at("signals").as_array();
            for (const auto& signal : signals)
                configured[std::string(entry.key())].emplace_back(signal.as_string().c_str());
        }
        if (configured.size() < 2) throw std::runtime_error("MULTI_GNSS_SYSTEM_STATUS_INCOMPLETE");
        const auto actual = zhangMultiGnssBuildSystemManifest(products, certificate, configured);
        if (actual != declared) throw std::runtime_error("MULTI_GNSS_SYSTEM_STATUS_PRODUCT_MISMATCH");
        bool any = false;
        for (const auto& entry : actual) any |= entry.value().as_object().at("any_pppar_usable").as_bool();
        if (any != epochAr) throw std::runtime_error("MULTI_GNSS_GLOBAL_AR_STATUS_MISMATCH");
        reason = "NONE";
        return true;
    }
    catch (const std::exception& error) { reason = error.what(); return false; }
}
