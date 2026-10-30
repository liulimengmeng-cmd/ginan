#include "common/zhangFullJointHandoff.hpp"
#include "common/zhangProductPublicationGate.hpp"
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message)
{ if (!condition) throw std::runtime_error(message); }

static ZhangProductIntegerTransition increment(
    int oldVersion, int newVersion, E_ObsCode signal)
{
    ZhangProductIntegerTransition out;
    out.observable = signal;
    out.valid = true;
    out.physicalEdges = {{"A", SatSys(E_Sys::GPS, 15)}, {"A", SatSys(E_Sys::GPS, 15)}};
    out.physicalArcVersions = {oldVersion, newVersion};
    out.coefficients = {-1, 1};
    return out;
}

int main()
{
    const auto signal = E_ObsCode::L1C;
    ZhangProductIntegerTransition zero;
    zero.valid = true;
    zero.observable = signal;
    auto chain = zhangAnchorFullJointStructuralChain("native322", 322, "public17", zero);
    require(chain.anchoredTo("native322", 322, "public17"), "initial public anchor missing");
    std::string reason;
    require(zhangAdvanceFullJointStructuralChain(chain, "native322", 322,
        "native323", 323, 323, increment(2, 5, signal), true, reason),
        "first unpublished chart link rejected");
    require(chain.publishedBackend == 322 && chain.tailBackend == 323 && chain.linkCount == 1,
        "structural advancement mutated last published endpoint");
    // This is the actual gap failure: publication stayed at 322, receipt begins
    // at 323. The public endpoint must NOT be used as the structural chain tail.
    require(!zhangFullJointTransportTailMatches("native322", "native323", 322, 323, 324, 324),
        "old single-step matcher unexpectedly bridges a publication gap");
    require(zhangAdvanceFullJointStructuralChain(chain, "native323", 323,
        "native324", 324, 324, increment(5, 9, signal), true, reason),
        "322 to unpublished 323 to recovered 324 chain rejected");
    require(chain.publishedBackend == 322 && chain.tailBackend == 324 && chain.linkCount == 2,
        "recovery moved the public endpoint before publication commit");
    require(chain.nativeToSegment.coefficients.size() == 2,
        "intermediate unavailable physical arc did not cancel exactly");
    // Independent substitution oracle: evaluate the symbol row directly, and
    // compare to (-N2+N5)+(-N5+N9), without reusing composition as the reference.
    std::size_t cases = 0;
    for (int n2 = -4; n2 <= 4; ++n2)
    for (int n5 = -3; n5 <= 3; ++n5)
    for (int n9 = -4; n9 <= 4; ++n9)
    {
        ZhangExactInteger actual = chain.nativeToSegment.affineOffsetCycles;
        for (std::size_t i = 0; i < chain.nativeToSegment.coefficients.size(); ++i)
        {
            const auto version = chain.nativeToSegment.physicalArcVersions[i];
            require(version == 2 || version == 9, "intermediate version survived reduction");
            actual += chain.nativeToSegment.coefficients[i] * (version == 2 ? n2 : n9);
        }
        require(actual == (-n2 + n5) + (-n5 + n9), "chain disagrees with independent substitution");
        ++cases;
    }
    auto unchanged = chain;
    require(!zhangAdvanceFullJointStructuralChain(chain, "native324", 324,
        "native325", 325, 325, increment(9, 10, signal), false, reason),
        "unproven same-time moment map advanced the chain");
    require(chain.tailBackend == unchanged.tailBackend &&
        chain.nativeToSegment.coefficients == unchanged.nativeToSegment.coefficients,
        "failed structural link partially committed");
    require(!zhangAdvanceFullJointStructuralChain(chain, "wrong-native", 324,
        "native325", 325, 325, increment(9, 10, signal), true, reason),
        "changed physical fingerprint accepted with matching backend number");
    require(!zhangAdvanceFullJointStructuralChain(chain, "native324", 324,
        "native325", 325, 326, increment(9, 10, signal), true, reason),
        "stale event receipt accepted on a later current chart");
    require(!zhangAdvanceFullJointStructuralChain(chain, "native324", 324,
        "native325", 325, 325, increment(9, 10, E_ObsCode::L2W), true, reason),
        "L2 structural increment was applied to L1 chain");
    require(!chain.anchoredTo("native322", 322, "public18"),
        "chain reattached to an unrelated public segment");
    auto newPublic = zhangAnchorFullJointStructuralChain("native324", 324, "public18", zero);
    require(newPublic.linkCount == 0 && !newPublic.anchoredTo("native322", 322, "public17"),
        "new public segment can silently jump back to old segment");
    auto limit = chain;
    limit.linkCount = ZhangFullJointStructuralChain::maximumLinks;
    require(!zhangAdvanceFullJointStructuralChain(limit, "native324", 324,
        "native325", 325, 325, increment(9, 10, signal), true, reason),
        "structural history grew beyond its configured hard bound");
    // Geometry alone is still not a certificate. Current-epoch integer proof
    // must be supplied separately even after the endpoint chain is repaired.
    ZhangProductPublicationEvidence gate;
    gate.requireFullJointHandoff = true;
    gate.priorPublishedIdentity = true;
    gate.nativeProductValid = gate.currentSpatialCertificateValid = gate.numericalStepValid = true;
    gate.conditionalProductMapProven = true;
    require(zhangProductPublicationGate(gate).route == ZhangProductPublicationRoute::NEW_SEGMENT,
        "structural chain reused a historical integer decision");
    std::cout << "STRUCTURAL_PUBLICATION_GAP_PASS substitutions=" << cases
        << " public_endpoint=322 suspended_tail=323 recovered_tail=324"
        << " history_integer_authority=0 signal_separation=1 fail_closed=1\n";
}
