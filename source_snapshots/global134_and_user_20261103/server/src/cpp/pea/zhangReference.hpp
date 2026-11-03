#pragma once

#include "common/zhangP0InvocationCache.hpp"

#include <cstdint>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include "common/eigenIncluder.hpp"
#include "common/enums.h"
#include "common/satSys.hpp"
#include "common/zhangFullRank.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangIntegerAudit.hpp"

struct KFState;
struct ZhangProductIntegerTransition;
struct ReceiverMap;
using Trace = std::ostream;

struct ZhangPreparedCertifiedTreeTransportHistory
{
    std::string authoritativeOwner, sourceRoot, sourceOrder;
    long int sourceEpoch = 0;
    std::map<std::string, std::string> sealedDirectories;
    bool sourceCertificateVerified = false;
};

/** Stage from the unmodified authoritative FLOAT before AR_EXIT retirement.
 * This call does NOT write a historical directory or change the KF state. */
bool zhangPrepareCertifiedTreeTransportHistory(Trace& trace, const KFState& authoritativeFloat,
    const std::map<E_Sys, std::pair<std::string, std::string>>& committedScopes,
    ZhangPreparedCertifiedTreeTransportHistory& prepared, std::string* failureReason = nullptr);

/** Call ONLY after the atomic fixed-product bundle COMMIT. The destination
 * must be the authoritative FLOAT state used by the next graph update.
 * Each explicit source scope must belong to that same owner; private branch
 * promotion without a separate owner-transfer proof is rejected. Archived
 * evidence is transport-only and never restores held/current integer rows. */
bool zhangCommitCertifiedTreeTransportHistory(Trace& trace, KFState& authoritativeFloat,
    const ZhangPreparedCertifiedTreeTransportHistory& prepared,
    std::string* failureReason = nullptr);

/** Validate the same post-publication owner/epoch/seals and draft only the
 * metadata map. The caller must commit it in its explicit metadata-only
 * authority interval; no state, integer registry, or graph is changed here. */
bool zhangStageCertifiedTreeTransportHistoryCommit(const KFState& authoritativeFloat,
    const ZhangPreparedCertifiedTreeTransportHistory& prepared,
    std::map<std::string,std::string>& proposedMetadata,
    std::string* failureReason = nullptr);

struct ZhangHistoricalTransportBridgeProof
{
    bool proven = false;
    bool currentArAuthority = false;
    ZhangExactInteger value = 0;
    std::string owner, reason = "NO_COMMITTED_TRANSPORT_CERTIFICATES";
    long int lastSourceEpoch = 0;
    std::vector<std::string> proofIds;
    std::shared_ptr<const std::vector<std::string>> proofStatements;
};

/** Exact prior-COMMIT consequence for a VERSIONED product-datum difference.
 * The caller still needs a full FLOAT map, actual writer binding, CURRENT
 * native AR certificate and full-Q publication. This object cannot admit a
 * present-epoch integer equation or modify a FLOAT state. */
ZhangHistoricalTransportBridgeProof zhangQueryCertifiedHistoricalTransport(
    const KFState& state, E_Sys system, const ZhangProductIntegerTransition& difference);

inline constexpr std::uint32_t ZHANG_GRAPH_CHECKPOINT_SCHEMA_VERSION = 5;
inline constexpr const char* ZHANG_GRAPH_CHECKPOINT_SECTION_NAME =
    "zhang_graph_runtime";

/** Snapshot needed to express held integer rows in invariant physical arcs. */
struct ZhangGraphIntegerContext
{
    ZhangGraphBasis                  basis;
    // SINGLE_S_BASIS_NATIVE: exact read-only view of basis, never a second
    // selector. LEGACY_PRODUCT_TREE explicitly retains the old selector.
    ZhangGraphBasis                  productBasis;
    std::string                      productBasisMode;
    // Lossless, length-framed identity of the complete immutable S snapshot:
    // runtime/epoch, versions, graph catalogue, signal arcs and current rows.
    // A branch may carry its authoritative owner's identity; it cannot infer
    // equality from treeEdges alone or from a runtime alias.
    std::string                      committedSBasisIdentity;
    // Legacy edge topology generations. Never use these as a physical
    // ambiguity-arc identity in the signal-versioned mode.
    std::map<ZhangGraphEdge, int>    arcVersions;
    // Authoritative physical ambiguity identities: (receiver, satellite,
    // observable, signalArcVersions[edge][observable]).
    std::map<ZhangGraphEdge, std::map<int, int>> signalArcVersions;
    // Currently admitted, QC-passing observations for each signal. This can
    // differ from the common dual-frequency graph after a single-band slip.
    std::map<int, std::set<ZhangGraphEdge>> signalObservationEdges;
    bool                             signalObservationFinalized = false;
    int                              eventId = 0;
    int                              productDatumVersion = 0;
    // Native S-gauge generation changes only with root/tree or physical
    // signal identities ON that tree, not a pure chord-catalogue insertion.
    int                              committedSGaugeVersion = 0;
    int                              representationVersion = 0;
    int                              floatGaugeVersion = 0;
    int                              integerComponentVersion = 0;
    // Invocation-local component selector, not a second FLOAT/KF chart.
    std::string                      componentId;
    std::string                      runtimeIdentity;
    std::string                      snapshotEpoch;
    std::string                      productBasisReason;
    int                              sourceComponentCount = 1;
    int                              cycleRank = 0;
    bool                             containsSourceRoot = true;
    bool                             productBasisRebuilt = false;
    bool                             initialized = false;
};

/** Enumerate real connected subgraphs of the committed FLOAT forest.  This
 * copies structural/physical identity only; it does not change x, P, runtime
 * graphs, or certify any relation between distinct components. */
std::vector<ZhangGraphIntegerContext> zhangGraphIntegerComponentContexts(
    const KFState& state, E_Sys system, std::string* failureReason = nullptr);

/** Explicit read-only selection for one AR invocation.  Existing graph
 * consumers, including cloned disposable branches, see this same immutable
 * snapshot until scope exit.  No selection survives to another epoch. */
class ZhangGraphIntegerContextScope
{
    std::uint64_t token_ = 0;
public:
    ZhangGraphIntegerContextScope(const KFState& state, E_Sys system,
        const ZhangGraphIntegerContext& context);
    ~ZhangGraphIntegerContextScope();
    ZhangGraphIntegerContextScope(const ZhangGraphIntegerContextScope&) = delete;
    ZhangGraphIntegerContextScope& operator=(const ZhangGraphIntegerContextScope&) = delete;
};

/** Read-only observation support attached to one physical ambiguity-arc
 * version.  These counters are diagnostic only and are never consulted by the
 * graph controller, Kalman filter, or ambiguity resolver. */
struct ZhangPhysicalArcDiagnostic
{
    ZhangGraphEdge edge;
    E_ObsCode      observable = E_ObsCode::NONE;
    int            arcVersion = 0;
    int            ageEpochs = 0;
    int            observationEpochs = 0;
};

/** Event-side support summary for one changed satellite product functional. */
struct ZhangProductFunctionalEventDiagnostic
{
    std::vector<ZhangPhysicalArcDiagnostic> oldSupport;
    std::vector<ZhangPhysicalArcDiagnostic> newSupport;
    int                                      commonObservationEpochs = 0;
    std::string                              eventCause = "UNCLASSIFIED";
};

/** Ephemeral receipt for one physical graph event in the current FLOAT epoch.
 * It is indexed by runtime, epoch, constellation, signal and satellite.  The
 * old row/moments are captured before state retirement.  The actual writer
 * projection is not known until the product target chart is committed; a
 * state-coordinate pullback alone must not set sameTimeMapProven.
 * Offsets and moments use metres.  No receipt from a previous epoch is valid. */
struct ZhangGraphProductEventReceipt
{
    GTime time;
    E_Sys system = E_Sys::NONE;
    E_ObsCode observable = E_ObsCode::NONE;
    SatSys satellite;
    std::string eventKind = "UNCLASSIFIED";
    std::string failureReason = "NOT_EVALUATED";
    bool physicalArcChanged = false;
    bool representableSameTimeMapProven = false;
    bool oldProductProjectionCaptured = false;
    bool actualWriterProjectionBound = false;
    bool actualWriterRowMatches = false;
    bool sameTimeMapProven = false;
    bool conditionalProductMapProven = false;
    std::string conditionalMapReason = "NOT_EVALUATED";
    bool fullJointMomentTransportProven = false;
    std::string jointMomentReason = "NOT_EVALUATED";
    int oldChartGeneration = -1;
    int newChartGeneration = -1;
    // Actual graph event counter before and after this transaction. These
    // are not inferred from one another: product-only chart changes may keep
    // the same backend event generation.
    int oldBackendGeneration = -1;
    int newBackendGeneration = -1;
    std::string fullJointTransportId;
    ZhangCoordinateRow oldPhaseRow;
    ZhangCoordinateRow newPhaseRow;
    ZhangCoordinateRow actualWriterPhaseRow;
    ZhangCoordinateRow bridgeStateRowCycles;
    // The old product target is a versioned physical integer functional,
    // independent of the old S-tree path used by the KF transform.
    SatSys productReferenceSatellite;
    std::vector<ZhangGraphEdge> oldProductFunctionalEdges;
    std::vector<int> oldProductFunctionalArcVersions;
    ZhangExactVector oldProductFunctionalCoefficients;
    ZhangExactInteger oldProductFunctionalAffineOffsetCycles = 0;
    double oldOffsetMetres = 0;
    double newOffsetMetres = 0;
    double oldPhaseMeanMetres = 0;
    double oldPhaseVarianceMetres2 = 0;
    double oldClockMeanMetres = 0;
    double oldClockVarianceMetres2 = 0;
    double oldClockPhaseCovarianceMetres2 = 0;
    // Cross covariance to every pre-event satellite phase row (both signals).
    std::map<std::pair<SatSys, E_ObsCode>, double>
        oldPhaseCrossCovarianceMetres2;
    std::map<SatSys, double> oldPhaseClockCrossCovarianceMetres2;
    std::map<SatSys, double> oldClockClockCrossCovarianceMetres2;
    std::vector<ZhangGraphEdge> physicalPathEdges;
    ZhangExactVector physicalPathCoefficients;
    std::vector<int> oldPhysicalArcVersions;
    std::vector<int> newPhysicalArcVersions;
    std::vector<ZhangGraphEdge> bridgePhysicalEdges;
    std::vector<int> bridgePhysicalArcVersions;
    std::vector<int> bridgePhysicalCoefficients;
    // A retired old edge can remain as a bounded correlated FLOAT auxiliary.
    // These keys are not current observed physical arcs or integer facts.
    std::map<ZhangGraphEdge, KFKey> retainedOldArcAuxiliaryKeys;
};

bool zhangGraphProductEventReceipt(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    const SatSys& satellite,
    ZhangGraphProductEventReceipt& receipt);

/** One versioned physical phase-edge row in the current FLOAT posterior.
 * Current edges use the committed tree/chord chart; an old retired version is
 * readable only from a registered correlated PRODUCT_HANDOFF auxiliary whose
 * active dependency pins it, or whose inactive cache TTL has not expired.
 * This numeric lookup confers no current-observation or integer authority. */
bool zhangGraphPhysicalArcFloatRow(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    const ZhangGraphEdge& edge,
    int arcVersion,
    ZhangCoordinateRow& rowMetres,
    std::string& failureReason);

/** Match one proposed integer bridge to the graph event's exact signed
 * versioned-arc product-path difference. This checks structure only; it does
 * not fix or certify the integer value. */
bool zhangGraphProductEventTransitionWitness(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    const SatSys& satellite,
    const SatSys& referenceSatellite,
    const ZhangProductIntegerTransition& transition,
    std::string& failureReason);

/** Atomically marginalise expired/unreferenced old-arc auxiliaries across
 * every system, independently of current observation availability. Live rows
 * remain correlated FLOAT coordinates; no integer certificate is created. */
bool retireZhangGraphProductEventAuxiliaries(
    Trace& trace,
    KFState& state,
    const std::string& stage = "UNSPECIFIED");

/** Compatibility entry: the system argument no longer narrows lifecycle
 * maintenance. A global registry must be retired and validated globally. */
bool retireZhangGraphProductEventAuxiliaries(
    Trace& trace,
    KFState& state,
    E_Sys system);

/** Validate the bounded FLOAT identity registry against actual KF states.
 * Used at AR exit and checkpoint boundaries instead of requiring zero rows. */
bool zhangGraphProductHandoffStateValid(
    const KFState& state, std::string& failureReason);

/** Refresh complete active public-product support before destructive graph
 * transforms. These APIs retain physical identities only, never integer RHS. */
void zhangClearGraphProductHandoffDependencies(KFState& state);
bool zhangRegisterGraphProductHandoffDependencies(
    KFState& state, E_Sys system, E_ObsCode observable,
    const std::vector<ZhangGraphEdge>& edges,
    const std::vector<int>& arcVersions);

/** Read-only scalar product and covariance diagnostics at a named FLOAT
 * processing stage; never modifies graph, KF, or integer evidence. */
void traceZhangCurrentGraphProductProjection(
    Trace& trace,
    const KFState& state,
    const char* stage);

/** Select usable Zhang references and apply an exact state/covariance S-transform when they change.
 *
 * This function must be called after the epoch prediction and before constructing measurements for
 * the same epoch.
 */
void updateZhangFullRankReferences(
    Trace&       trace,
    ReceiverMap& receiverMap,
    KFState&     kfState
);

/** Whether the active general-tree controller models this baseline edge.
 *
 * Fixed-star mode always returns true. General-tree mode uses the committed
 * phase observation set, which may span several uncertain FLOAT components.
 * This does not authorize integer products between those components.
 */
bool zhangGraphModelsObservation(
    const KFState&    kfState,
    const std::string& receiver,
    const SatSys&      satellite,
    E_ObsCode          code
);

// The Zhang graph transaction owns baseline physical arc retirement, including
// tree arcs which have no independent AMBIGUITY state for generic cleanup.
bool zhangGraphOwnsAmbiguityLifecycle(E_Sys system, E_ObsCode code);

// Disposable-state regression entry point for the production pivot/retirement audit.
bool applyZhangGraphRetirementForAudit(
    Trace& trace, KFState& state, E_Sys system,
    const std::vector<E_ObsCode>& codes, const ZhangGraphBasis& oldBasis,
    const ZhangGraphBasis& newBasis, const std::set<ZhangGraphEdge>& retired,
    std::string& failureReason);

/** Whether this edge owns an integer fundamental-cycle ambiguity state. */
bool zhangGraphRetainsAmbiguity(
    const KFState&    kfState,
    const std::string& receiver,
    const SatSys&      satellite,
    E_ObsCode          code
);

/** Read-only first rejection reason for the direct current-cycle admission.
 * Diagnostic only: retained stochastic support is evaluated separately. */
std::string zhangGraphAmbiguityAdmissionReason(
    const KFState& state, const std::string& receiver,
    const SatSys& satellite, E_ObsCode code);

/** Whether a satellite belongs to the retained datum-state component and may
 * therefore be emitted as a current network product. */
bool zhangGraphProductSatelliteActive(
    const KFState& kfState,
    const SatSys&  satellite
);

/** Current graph and physical-arc versions for exact held-lattice transport. */
bool zhangGraphIntegerContext(
    const KFState&             kfState,
    E_Sys                      system,
    ZhangGraphIntegerContext& context
);

/** Read one exact signal-arc version. Missing identity never falls back to
 * the edge topology generation or to version zero. */
std::optional<int> zhangGraphSignalArcVersion(
    const ZhangGraphIntegerContext& context,
    const ZhangGraphEdge& edge,
    E_ObsCode observable);

/** Known signal-specific physical versions over the FLOAT/product basis union.
 * This is an identity catalogue, NOT current observation/AR admission. Missing
 * versions remain unknown; final QC signal support is checked separately. */
bool zhangGraphSignalArcVersions(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    std::map<ZhangGraphEdge, int>& versions,
    std::string& failureReason);

/** Commit the signal rows actually retained by the final /PPP KF factor.
 * A robustly deweighted row must carry the `zhangPppDeweighted` metadata
 * marker; a merely present pre-QC observation is not current AR evidence. */
bool zhangFinalizeGraphSignalObservationAdmission(
    const KFState& state,
    E_Sys system,
    const KFMeas& finalMeasurement,
    std::string& failureReason);

/** Old graph checkpoints cannot recover pre-migration signal identities.
 * The FLOAT state is retained, but prior integer/product continuity must be
 * discarded by the other checkpoint owners before publication. */
bool zhangGraphRequiresFreshProductIdentity(const KFState& state);

/** Diagnose the physical support that was removed/introduced by a product
 * functional change.  Inputs contain only non-zero support arcs. */
bool zhangProductFunctionalEventDiagnostic(
    const KFState&                       kfState,
    E_Sys                                system,
    const std::vector<ZhangGraphEdge>&   oldEdges,
    const std::vector<int>&              oldArcVersions,
    const std::vector<ZhangGraphEdge>&   newEdges,
    const std::vector<int>&              newArcVersions,
    ZhangProductFunctionalEventDiagnostic& diagnostic
);

bool zhangProductFunctionalEventDiagnostic(
    const KFState&                       kfState,
    E_Sys                                system,
    E_ObsCode                            observable,
    const std::vector<ZhangGraphEdge>&   oldEdges,
    const std::vector<int>&              oldArcVersions,
    const std::vector<ZhangGraphEdge>&   newEdges,
    const std::vector<int>&              newArcVersions,
    ZhangProductFunctionalEventDiagnostic& diagnostic
);

/** Copy graph coordinates to a distinct, explicitly identified branch.
 * A rejected clone leaves both source and destination graph runtimes intact. */
[[nodiscard]] bool cloneZhangGraphRuntime(
    const KFState& source,
    const KFState& destination,
    std::string* failureReason = nullptr
);
/** Rebind V2 whole-FLOAT metadata only on an exact, unconditioned copy.
 * Legacy empty metadata uses the unchanged graph-clone contract. */
[[nodiscard]] bool cloneZhangGraphRuntime(
    const KFState& source,
    KFState& destination,
    std::string* failureReason = nullptr
);

/** Remove graph coordinates belonging to a disposable state copy. */
void eraseZhangGraphRuntime(const KFState& state);

/** Exact, small snapshot of FLOAT S coordinates and physical arc lifecycle.
 * Product target trees, integer certificates and product authorization are
 * deliberately excluded: updating them cannot change the FLOAT chart. */
struct ZhangFloatGraphAuthorityToken
{
    bool valid = false;
    std::string runtimeIdentity;
    std::string physicalState;
};

bool captureZhangFloatGraphAuthorityToken(
    const KFState& state,
    ZhangFloatGraphAuthorityToken& token,
    std::string& failureReason);

bool zhangFloatGraphAuthorityUnchanged(
    const KFState& state,
    const ZhangFloatGraphAuthorityToken& token,
    std::string& failureReason);

/** Export every Zhang graph/outage runtime owned by one filter into a
 * versioned, pointer-free checkpoint section.
 *
 * runtimeId is the stable identity supplied by the checkpoint orchestrator;
 * an empty identity is rejected.  The returned payload never contains the
 * address of state.
 */
bool exportZhangGraphCheckpointSection(
    const KFState&      state,
    const std::string& runtimeId,
    std::string&       payload,
    std::string&       failureReason
);

/** Decode and fully validate a graph checkpoint section without changing any
 * live runtime state.  Bundle restoration must preflight every section before
 * committing any one of them.
 */
bool validateZhangGraphCheckpointSection(
    const std::string& runtimeId,
    const std::string& payload,
    std::string&       failureReason
);

/** Restore a Zhang graph/outage checkpoint section onto one configured filter.
 *
 * The embedded runtime identity must exactly match runtimeId.  The complete
 * section is decoded and validated before the destination's runtime bindings
 * are replaced atomically; a false return leaves all live runtime maps
 * unchanged.
 */
bool importZhangGraphCheckpointSection(
    KFState&           state,
    const std::string& runtimeId,
    const std::string& payload,
    std::string&       failureReason
);

/** Apply one exact, dimension-preserving graph S-basis transform to a
 * disposable audit branch and return the numerical state transform.  This
 * entry point deliberately does not update product-continuity history or any
 * live graph controller state. */
bool applyZhangGraphBasisTransformForAudit(
    Trace&                         trace,
    KFState&                       branch,
    E_Sys                          system,
    const std::vector<E_ObsCode>&  baselineObservables,
    const ZhangGraphBasis&         oldBasis,
    const ZhangGraphBasis&         newBasis,
    SparseMatrix<double>&          transform,
    std::string&                   failureReason
);

// Stochastic support is distinct from direct, currently observed AR eligibility.
bool zhangGraphStochasticSupportValid(const KFState&, const std::string&,
 const SatSys&, E_ObsCode);

// Lightweight, read-only identity for invocation-local graph snapshots.
ZhangP0GraphStamp zhangGraphIntegerContextStamp(
	const KFState& state,
	E_Sys system);
