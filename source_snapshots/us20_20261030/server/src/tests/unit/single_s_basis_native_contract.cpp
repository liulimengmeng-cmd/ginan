// SOURCE ONLY: these contract assertions have NOT been compiled or executed.
// The graph and physical datum path below are literal independent fixtures;
// they are not reconstructed by calling the production tree selector twice.
#include "common/zhangSingleSBasisNativeProduct.hpp"
#include "common/zhangProductConditioningReceipt.hpp"
#include "common/zhangR47Candidate.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "pea/zhangReference.hpp"

#include <stdexcept>

int main()
{
    auto require = [](bool value, const char* reason)
    { if (!value) throw std::runtime_error(reason); };
    const SatSys g1(E_Sys::GPS, 1), g2(E_Sys::GPS, 2);
    const ZhangGraphEdge a1{"AAAA", g1}, a2{"AAAA", g2};
    const ZhangGraphEdge b1{"BBBB", g1}, b2{"BBBB", g2};
    ZhangGraphBasis s;
    s.rootReceiver = "AAAA";
    s.edges = {a1, a2, b1, b2};
    s.treeEdges = {a1, a2, b1};
    s.receivers = {"AAAA", "BBBB"};
    s.satellites = {g1, g2};
    s.componentCount = 1;
    s.connected = true;
    const auto native = ProductRelationBasisBuilder::buildSingleSBasisNative(
        s, "AUTHORITATIVE-S-SNAPSHOT", g1, E_Sys::GPS, E_ObsCode::L1C);
    require(native.valid && native.singleSBasisNative &&
        native.committedSBasisIdentity == "AUTHORITATIVE-S-SNAPSHOT" &&
        native.zeroDatumCorrection && native.exactRank == 0 &&
        native.datumCorrectionRank == 0,
        "Same committed S did not produce an explicitly native zero-G map");
    require(native.satellites == std::vector<SatSys>({g1, g2}) &&
        native.namedRelations.size() == 1 && native.exactRowBasis.empty() &&
        native.affineRecoveryParents.empty() && native.independentNamedIndices.empty(),
        "Zero G lost its catalogue or manufactured accepted integer rows");
    require(!ProductRelationBasisBuilder::buildSingleSBasisNative(
        s, "", g1, E_Sys::GPS, E_ObsCode::L1C).valid,
        "Native builder accepted an unbound S snapshot");

    ZhangZeroDatumNativeRequest request;
    request.requested = true;
    request.satellites = {g1, g2};
    request.conditioningRank = 1;
    request.singleSBasisNative = true;
    request.committedSBasisIdentity = "AUTHORITATIVE-S-SNAPSHOT";
    ZhangR47Candidate noDecision;
    noDecision.productConsequences.zeroDatumNative = request;
    require(!zhangR47CandidateContractValid(noDecision),
        "A native zero-G request alone became an integer transaction");
    ZhangProductConditioningIdentity identity{
        "FLOAT-OWNER", 1000, "FLOAT-POSTERIOR", "CURRENT-CHART",
        "ACTUAL-FLOAT-XP", "FULL-STATE-ORDER", "ACTUAL-FIXED-XP", "FULL-STATE-ORDER"};
    ZhangProductConditioningReceipt noReceipt;
    ZhangCurrentEpochNumericConstraint zeroIdentity;
    zeroIdentity.value = 0;
    std::vector<ZhangCurrentEpochNumericConstraint> onlyZeroIdentity{zeroIdentity};
    std::string reason;
    require(!noReceipt.bind(identity, onlyZeroIdentity, reason),
        "A structural 0=0 issued a numerical conditioning receipt");
    require(!zhangZeroDatumNativePosteriorVerified(request,
        noReceipt, identity, onlyZeroIdentity, reason),
        "A fabricated zero conditioning recipe authorized native AR");

    ZhangGraphIntegerContext source;
    source.basis = source.productBasis = s;
    source.initialized = true;
    source.productBasisMode = "SINGLE_S_BASIS_NATIVE";
    source.committedSBasisIdentity = request.committedSBasisIdentity;
    source.runtimeIdentity = "AUTHORITATIVE-FLOAT-OWNER";
    source.snapshotEpoch = "2024-07-17T02:51:00";
    source.representationVersion = source.productDatumVersion = 3;
    source.committedSGaugeVersion = 3;
    source.floatGaugeVersion = 2;
    source.integerComponentVersion = 1; // Not a product-datum ordering bound.
    source.signalObservationFinalized = true;
    for (const auto& edge : s.edges)
        for (const auto code : {E_ObsCode::L1C, E_ObsCode::L2W})
            source.signalArcVersions[edge][static_cast<int>(code)] = 7;
    source.signalObservationEdges[static_cast<int>(E_ObsCode::L1C)] = s.edges;
    source.signalObservationEdges[static_cast<int>(E_ObsCode::L2W)] = s.edges;
    auto branch = source;
    branch.runtimeIdentity = "DISPOSABLE-CONDITIONED-BRANCH";
    branch.committedSBasisIdentity = "BRANCH-LEASE-S-SNAPSHOT";
    require(zhangSingleSBasisNativeSnapshotMatches(request, source, reason),
        "Actual authoritative S snapshot was rejected");
    require(zhangSingleSBasisNativeGraphCoordinatesMatch(source, branch, reason),
        "Equal physical S coordinates were rejected solely for another runtime lease");
    auto moreRepresentationEvents = source;
    moreRepresentationEvents.representationVersion = 4;
    require(zhangSingleSBasisNativeViewValid(moreRepresentationEvents, reason),
        "Full chart generation was incorrectly required to equal the native S-gauge datum");
    auto wrongIdentity = source;
    wrongIdentity.committedSBasisIdentity = "OTHER-EPOCH-OR-S-SNAPSHOT";
    require(!zhangSingleSBasisNativeSnapshotMatches(request, wrongIdentity, reason),
        "A request was rebound to another authoritative S identity");

    // Literal S datum for G02-G01 is N(AAAA,G02)-N(AAAA,G01).
    ZhangProductIntegerFunctional l1;
    l1.valid = true;
    l1.satellite = g2;
    l1.referenceSatellite = g1;
    l1.observable = E_ObsCode::L1C;
    l1.physicalEdges = {a1, a2, b1, b2};
    l1.networkCoefficients = {-1, 1, 0, 0};
    l1.physicalArcVersions = {7, 7, 7, 7};
    auto l2 = l1;
    l2.observable = E_ObsCode::L2W;
    require(zhangSingleSBasisNativeSignalSupport(source, g2, E_ObsCode::L1C, l1, reason) &&
        zhangSingleSBasisNativeSignalSupport(source, g2, E_ObsCode::L2W, l2, reason),
        "Complete current two-signal physical paths were rejected");

    auto retainedButNotCurrent = source;
    retainedButNotCurrent.signalObservationEdges[static_cast<int>(E_ObsCode::L1C)].erase(a2);
    require(retainedButNotCurrent.signalArcVersions.at(a2).at(static_cast<int>(E_ObsCode::L1C)) == 7 &&
        retainedButNotCurrent.signalObservationEdges.at(static_cast<int>(E_ObsCode::L1C)).contains(b2),
        "Fixture failed to retain the physical arc and an alternate current G02 observation");
    require(!zhangSingleSBasisNativeSignalSupport(
        retainedButNotCurrent, g2, E_ObsCode::L1C, l1, reason) &&
        reason == "SINGLE_S_BASIS_NATIVE_DATUM_PATH_SIGNAL_NOT_CURRENTLY_ADMITTED",
        "A retained state or alternate observation silently replaced the committed S path");
    auto noSignal = source;
    noSignal.signalObservationEdges.erase(static_cast<int>(E_ObsCode::L2W));
    require(!zhangSingleSBasisNativeSignalSupport(noSignal, g2, E_ObsCode::L2W, l2, reason),
        "Missing current signal support authorized a dual-frequency native row");
    auto wrongArc = source;
    wrongArc.signalArcVersions[a2][static_cast<int>(E_ObsCode::L1C)] = 8;
    require(!zhangSingleSBasisNativeSignalSupport(wrongArc, g2, E_ObsCode::L1C, l1, reason) &&
        reason == "SINGLE_S_BASIS_NATIVE_DATUM_PATH_PHYSICAL_ARC_VERSION_MISMATCH",
        "A physical arc version mismatch inherited old datum support");
    auto unfinalized = source;
    unfinalized.signalObservationFinalized = false;
    require(!zhangSingleSBasisNativeSignalSupport(unfinalized, g2, E_ObsCode::L1C, l1, reason),
        "Unfinalized candidate observations were treated as final factor admission");

    auto independentProductTree = source;
    independentProductTree.productBasis.treeEdges = {a1, b1, b2};
    require(!zhangSingleSBasisNativeViewValid(independentProductTree, reason) &&
        !zhangSingleSBasisNativeSnapshotMatches(request, independentProductTree, reason),
        "An independent P tree was accepted as the native committed S view");
    auto nonzeroGamma = l1;
    nonzeroGamma.affineOffsetCycles = 1;
    require(!zhangSingleSBasisNativeSignalSupport(source, g2, E_ObsCode::L1C, nonzeroGamma, reason),
        "Native g=0 silently accepted an extra affine integer alignment");
}
