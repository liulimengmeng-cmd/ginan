#include "common/zhangCurrentEvidenceAuthoritySnapshot.hpp"

#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace
{
struct Authority
{
    std::string posterior;
    static inline bool failCopy = false;

    Authority() = default;
    explicit Authority(std::string value) : posterior(std::move(value)) {}
    Authority(const Authority& other) : posterior(other.posterior)
    {
        if (failCopy) throw std::runtime_error("injected authority copy failure");
    }
    Authority& operator=(const Authority&) = default;
};

using Evidence = std::map<std::pair<std::string, int>, std::string>;
using Authorities = std::map<std::string, Authority>;
}

int main()
{
    Evidence evidence{{{"network", 1}, "accepted-P0"}, {{"foreign", 2}, "foreign-P0"}};
    Authorities authorities{{"network", Authority("accepted-P0")}, {"foreign", Authority("foreign-P0")}};
    const ZhangCurrentEvidenceAuthoritySnapshot baseline(evidence, authorities);

    // A rejected candidate must not leave its new authority pointing at the
    // restored previous evidence, including repeated rollback paths.
    evidence.at({"network", 1}) = "candidate-P1";
    authorities.at("network").posterior = "candidate-P1";
    evidence.erase({"foreign", 2});
    authorities.erase("foreign");
    baseline.restore(evidence, authorities);
    assert(evidence.at({"network", 1}) == "accepted-P0");
    assert(authorities.at("network").posterior == "accepted-P0");
    assert(evidence.at({"foreign", 2}) == "foreign-P0");
    assert(authorities.at("foreign").posterior == "foreign-P0");
    baseline.restore(evidence, authorities);
    assert(evidence.size() == 2 && authorities.size() == 2);

    // If preparation of the second registry fails, neither live registry is
    // changed. In particular evidence is not rolled back without authority.
    evidence.at({"network", 1}) = "candidate-P2";
    authorities.at("network").posterior = "candidate-P2";
    Authority::failCopy = true;
    bool threw = false;
    try { baseline.restore(evidence, authorities); }
    catch (const std::runtime_error&) { threw = true; }
    Authority::failCopy = false;
    assert(threw);
    assert(evidence.at({"network", 1}) == "candidate-P2");
    assert(authorities.at("network").posterior == "candidate-P2");
    baseline.restore(evidence, authorities);
    assert(evidence.at({"network", 1}) == "accepted-P0");
    assert(authorities.at("network").posterior == "accepted-P0");

    // The epoch-only presearch baseline has no current evidence or authority.
    Evidence emptyEvidence;
    Authorities emptyAuthorities;
    const ZhangCurrentEvidenceAuthoritySnapshot epochStart(emptyEvidence, emptyAuthorities);
    emptyEvidence.emplace(std::make_pair("network", 1), "rejected-current-P");
    emptyAuthorities.emplace("network", Authority("rejected-current-P"));
    epochStart.restore(emptyEvidence, emptyAuthorities);
    assert(emptyEvidence.empty() && emptyAuthorities.empty());

    // Pending diagnostics are preserved without manufacturing authority.
    Evidence diagnostics{{{"network", 1}, "pending-invalid-diagnostic"}};
    Authorities noAuthority;
    const ZhangCurrentEvidenceAuthoritySnapshot pending(diagnostics, noAuthority);
    diagnostics.clear();
    noAuthority.emplace("network", Authority("uncommitted"));
    pending.restore(diagnostics, noAuthority);
    assert(diagnostics.at({"network", 1}) == "pending-invalid-diagnostic");
    assert(noAuthority.empty());

    // A successful transaction does not call restore. Snapshot destruction
    // has no implicit action and cannot undo already committed publication.
    {
        const ZhangCurrentEvidenceAuthoritySnapshot unused(evidence, authorities);
        evidence.at({"network", 1}) = "committed-P3";
        authorities.at("network").posterior = "committed-P3";
    }
    assert(evidence.at({"network", 1}) == "committed-P3");
    assert(authorities.at("network").posterior == "committed-P3");
    std::cout << "CURRENT_EVIDENCE_AUTHORITY_SNAPSHOT_PASS" << std::endl;
}
