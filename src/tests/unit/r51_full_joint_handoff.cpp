#include "common/zhangFullJointHandoff.hpp"
#include "common/zhangProductPublicationGate.hpp"
#include "common/zhangProductConditionalMap.hpp"
#include <iostream>
#include <stdexcept>

static void require(bool value, const char* reason)
{ if (!value) throw std::runtime_error(reason); }

int main()
{
    require(zhangFullJointTransportTailMatches("physical-A", "physical-A", 7, 7, 8, 8),
        "complete consecutive native chart chain rejected");
    require(!zhangFullJointTransportTailMatches("physical-A", "physical-A", 7, 8, 9, 9),
        "suspended publication gap silently attached to wrong chart tail");
    require(!zhangFullJointTransportTailMatches("physical-A", "physical-B", 7, 7, 8, 8),
        "same backend generation bypassed changed physical target tail");
    require(!zhangFullJointTransportTailMatches("physical-A", "physical-A", 7, 7, 8, 9),
        "stale event receipt accepted on a later chart");
    const ZhangGraphEdge oldEdge{"A", SatSys(E_Sys::GPS, 13)};
    const ZhangGraphEdge newEdge{"B", SatSys(E_Sys::GPS, 13)};
    ZhangProductIntegerTransition first;
    first.observable = E_ObsCode::L2W;
    first.valid = true;
    first.physicalEdges = {oldEdge, newEdge};
    first.physicalArcVersions = {2, 5};
    first.coefficients = {-1, 1};
    auto second = first;
    second.physicalEdges = {newEdge, oldEdge};
    second.physicalArcVersions = {5, 7};
    second.coefficients = {-1, 1};
    const auto total = zhangComposeFullJointDifference(first, second);
    require(total.valid && total.coefficients.size() == 2,
        "intermediate physical arc was not cancelled exactly");
    // Independent substitution oracle, not a second call to composition.
    for (int n2 = -3; n2 <= 3; n2++)
    for (int n5 = -2; n5 <= 2; n5++)
    for (int n7 = -3; n7 <= 3; n7++)
    {
        int actual = 0;
        for (std::size_t i = 0; i < total.coefficients.size(); i++)
            actual += total.coefficients[i].convert_to<int>() *
                (total.physicalArcVersions[i] == 2 ? n2 : n7);
        require(actual == (-n2 + n5) + (-n5 + n7),
            "direct and chained complete segment target disagree");
    }
    auto wrongSignal = second;
    wrongSignal.observable = E_ObsCode::L1C;
    require(!zhangComposeFullJointDifference(first, wrongSignal).valid,
        "different physical signal identities were combined");
    ZhangProductPublicationEvidence e;
    e.nativeProductValid = e.currentSpatialCertificateValid = e.numericalStepValid = true;
    e.priorPublishedIdentity = e.samePhysicalSegment = e.sameProductFunctional = e.sameBackendChart = true;
    e.requireFullJointHandoff = true;
    require(zhangProductPublicationGate(e).route == ZhangProductPublicationRoute::NEW_SEGMENT,
        "unchanged relative functional bypassed complete joint handoff");
    e.conditionalProductMapProven = e.exactSameTimeMapProven = true;
    require(zhangProductPublicationGate(e).route == ZhangProductPublicationRoute::NEW_SEGMENT,
        "geometry-only transport was mistaken for current integer authorization");
    e.fullJointHandoffProven = true;
    e.structuralProductMapProven = true;
    e.samePhysicalSegment = false;
    require(zhangProductPublicationGate(e).route == ZhangProductPublicationRoute::SAME_SEGMENT &&
        !e.currentEpochBridgeProven,
        "structural map identity was converted to a fake integer decision");
    std::array<ZhangCoordinateRow, 2> rows;
    require(zhangCheckProductOverlapNumericRows(rows, rows, rows,
        {false, false}, false, false).valid,
        "two structural-zero signals require invented retired auxiliaries");
    require(!zhangCheckProductOverlapNumericRows(rows, rows, rows,
        {true, true}, false, true).valid,
        "real retired-arc requirement was ignored");
    std::cout << "FULL_JOINT_HANDOFF_PASS substitutions=245 signal_separation=1 structural_zero=1\n";
}
