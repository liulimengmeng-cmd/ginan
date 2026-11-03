#pragma once

#include <map>
#include <string>
#include <vector>

// An endpoint's FLOAT coordinate need not be integer relative to the KF
// chart reference. Only a PAIR inside a currently certified product component
// is authorised as an integer satellite difference. These certificates contain
// no integer value, historical constraint or FLOAT-state mutation.
struct ZhangUserGonsisComponentCertificate
{
    bool pppArUsable = false;
    bool dualFrequencyArValid = false;
    std::string component;
    int version = -1;
    int alignmentGeneration = -1;
};

inline bool zhangUserGonsisComponentCertificateValid(
    const ZhangUserGonsisComponentCertificate& certificate)
{
    return certificate.pppArUsable && certificate.dualFrequencyArValid &&
        !certificate.component.empty() && certificate.component != "NONE" &&
        certificate.component != "UNRESOLVED" && certificate.component != "UNAVAILABLE" &&
        certificate.version >= 0 && certificate.alignmentGeneration >= 0;
}

inline bool zhangUserGonsisComponentPairCertified(
    const ZhangUserGonsisComponentCertificate& first,
    const ZhangUserGonsisComponentCertificate& second)
{
    return zhangUserGonsisComponentCertificateValid(first) &&
        zhangUserGonsisComponentCertificateValid(second) &&
        first.component == second.component && first.version == second.version &&
        first.alignmentGeneration == second.alignmentGeneration;
}

inline std::string zhangUserGonsisComponentToken(
    const ZhangUserGonsisComponentCertificate& certificate)
{
    if (!zhangUserGonsisComponentCertificateValid(certificate)) return {};
    // Hex encoding is unambiguous and never introduces whitespace/trace
    // delimiters, even for a provider-supplied component name.
    static constexpr char hex[] = "0123456789ABCDEF";
    std::string name;
    for (unsigned char byte : certificate.component)
    {
        name += hex[byte >> 4];
        name += hex[byte & 15];
    }
    return "C" + std::to_string(certificate.component.size()) + "_" + name +
        "_V" + std::to_string(certificate.version) + "_A" +
        std::to_string(certificate.alignmentGeneration);
}

inline std::string zhangUserGonsisDualComponentDomain(
    const ZhangUserGonsisComponentCertificate& firstBand,
    const ZhangUserGonsisComponentCertificate& secondBand)
{
    const auto first = zhangUserGonsisComponentToken(firstBand);
    const auto second = zhangUserGonsisComponentToken(secondBand);
    if (first.empty() || second.empty()) return {};
    // The two bands may have differently named certified lattices. A PAIR
    // must match BOTH band tuples; equality between the bands is not required.
    return "F1_" + first + "_F2_" + second;
}

inline std::map<std::string, std::vector<int>> zhangUserGonsisPartitionDomains(
    const std::vector<std::string>& domains)
{
    std::map<std::string, std::vector<int>> groups;
    for (int index = 0; index < static_cast<int>(domains.size()); ++index)
        if (!domains[index].empty()) groups[domains[index]].push_back(index);
    return groups;
}
