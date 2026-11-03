#pragma once

#include <type_traits>
#include <utility>

// Current evidence and the joint-delivery authority referencing it are one
// transactional state. This snapshot preserves that state; it never validates
// a certificate, creates authority, or makes old evidence current. Consumers
// and the checkpoint codec still enforce epoch/root/order/posterior/chart.
template<class EvidenceRegistry, class AuthorityRegistry>
class ZhangCurrentEvidenceAuthoritySnapshot
{
    EvidenceRegistry evidence;
    AuthorityRegistry authorities;

public:
    ZhangCurrentEvidenceAuthoritySnapshot(
        const EvidenceRegistry& currentEvidence,
        const AuthorityRegistry& currentAuthorities)
        : evidence(currentEvidence), authorities(currentAuthorities)
    {}

    void restore(EvidenceRegistry& currentEvidence,
                 AuthorityRegistry& currentAuthorities) const
    {
        static_assert(std::is_nothrow_swappable_v<EvidenceRegistry> &&
                      std::is_nothrow_swappable_v<AuthorityRegistry>,
                      "Evidence/authority rollback requires nonthrowing swaps");
        // Allocate both replacement registries before changing either live
        // registry. A failed allocation/copy cannot leave a half rollback.
        auto restoredEvidence = evidence;
        auto restoredAuthorities = authorities;
        using std::swap;
        swap(currentEvidence, restoredEvidence);
        swap(currentAuthorities, restoredAuthorities);
    }
};
