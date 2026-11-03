#include "common/zhangP0ResourceProbe.hpp"
#include "common/zhangProductFloatTargetBatch.hpp"
#include "common/zhangFloatTargetResolverBatch.hpp"
#include "pea/zhangReference.hpp"
#include "common/zhangR48ProductDatum.hpp"
#include "common/zhangActiveGraphBasis.hpp"
#include "common/zhangSystemClockDatum.hpp"
#include "common/zhangEpochComponentContext.hpp"
#include "common/zhangPreselectedProductTarget.hpp"
#include "common/zhangCommittedSGraphEvent.hpp"
#include "common/zhangProductBasisMode.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iterator>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/log/trivial.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/trace.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangFullRank.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductPublicationGate.hpp"
#include "common/zhangProductConditionalMap.hpp"
#include "common/zhangProductChartIdentity.hpp"
#include "common/zhangProductHandoffFloatRegistry.hpp"
#include "common/zhangIntegerSupportResidualAudit.hpp"
#include "common/zhangIntegerAudit.hpp"
#include "common/zhangCertifiedTreeTransport.hpp"
#include "common/zhangR51Integrity.hpp"
#include "ambres/GNSSambres.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "pea/zhangPppAr.hpp"

using std::map;
using std::set;
using std::string;
using std::vector;

namespace
{
struct ReferenceAvailability
{
    map<string, set<SatSys>> satellitesByReceiver;
    map<SatSys, double>      elevationScore;
    set<ZhangGraphEdge>      edges;
    set<ZhangGraphEdge>      codeEdges;
    set<ZhangGraphEdge>      rawEdges;
    set<ZhangGraphEdge>      discontinuousEdges;
    map<ZhangGraphEdge, set<E_ObsCode>> discontinuitySignals;
    map<ZhangGraphEdge, vector<SigStat*>> discontinuityStats;
    set<ZhangGraphEdge>      qcExcludedEdges;
    set<ZhangGraphEdge>      elevationExcludedEdges;
    set<ZhangGraphEdge>      signalUnavailableEdges;
    set<ZhangGraphEdge>      phaseQcEdges;
    set<ZhangGraphEdge>      rejectedPhaseEdges;
    set<ZhangGraphEdge>      continuityUnknownEdges;
    map<int, set<ZhangGraphEdge>> signalEdges;
    map<ZhangGraphEdge, double> edgeQuality;
};

struct ReferenceOutageState
{
    int receiverEpochs = 0;
    int satelliteEpochs = 0;
};

using ZhangGraphRuntimeKey = std::pair<string, E_Sys>;

map<ZhangGraphRuntimeKey, ReferenceOutageState> outageStateMap;

struct GraphRuntimeState
{
    // Immutable after initialization/import; a mode switch needs an explicit
    // migration, not reinterpretation of already estimated KF coordinates.
    string productBasisMode;
    struct EdgeHistory
    {
        int continuousEpochs = 0;
        int outageEpochs = 0;
        int arcVersion = 0;
    };

    struct ArcVersionObservationHistory
    {
        int         firstObservedEpoch = -1;
        int         lastObservedEpoch = -1;
        int         observationEpochs = 0;
        vector<int> observedEpochs;
    };

    ZhangGraphBasis          basis;
    ZhangGraphBasis          activeBasis;
    ZhangGraphBasis          productBasis;
    // Legacy topology-generation snapshot, never a physical arc identity.
    map<ZhangGraphEdge, int> productArcVersions;
    map<ZhangGraphEdge, map<int, int>> productSignalArcVersions;
    set<string>              productCoreReceivers;
    set<ZhangGraphEdge>      observationEdges;
    set<ZhangGraphEdge>      stateEdges;
    map<ZhangGraphEdge, EdgeHistory> edgeHistory;
    map<pair<ZhangGraphEdge, int>, ArcVersionObservationHistory>
                                arcVersionObservationHistory;
    map<std::tuple<ZhangGraphEdge, int, int>,
        ArcVersionObservationHistory> signalArcVersionObservationHistory;
    int                         epochIndex = -1;
    string                      lastProductEventCause = "INITIALISE";
    bool                     initialized = false;
    int                      deferredEpochs = 0;
    int                      datumVersion = 0;
    int                      representationVersion = 0;
    int                      floatGaugeVersion = 0;
    int                      integerComponentVersion = 0;
    int                      eventCounter = 0;
    int                      productDatumVersion = 0;
    int                      committedSGaugeVersion = 0;
    bool                     productInitialized = false;
    // Product authorization is bound to a component revision; the product
    // gauge and integer component are distinct histories, not equal counters.
    int productAuthorityComponentVersion = -1;
    map<ZhangGraphEdge, map<int, int>> phaseArcVersions;
    map<ZhangGraphEdge, map<int, int>> phaseAdmissions;
    map<int, set<ZhangGraphEdge>> signalObservationEdges;
    bool signalObservationFinalized = false;
    std::uint64_t signalObservationFactorSequence = 0;
    bool requiresFreshProductIdentity = false;
    string phaseAdmissionTime;
    int phaseAdmissionChart = -1;
    int                      treeSlipShadowEvents = 0;
};

map<ZhangGraphRuntimeKey, GraphRuntimeState> graphStateMap;
string zhangGraphRuntimeId(const KFState& state);

bool singleSBasisNative()
{
    return zhangSingleSBasisNativeEnabled(acsConfig.zhangFullRank.product_basis_mode);
}

string configuredProductBasisMode()
{
    return zhangProductBasisModeName(
        zhangParseProductBasisMode(acsConfig.zhangFullRank.product_basis_mode));
}

// This cache confers no integer certification.  It is only an exact mirror of
// the ONE committed S chart; current observation support remains separately
// signal/QC bound.  In particular a retained tree edge is not a current row.
void syncSingleSBasisProductView(GraphRuntimeState& runtime)
{
    runtime.activeBasis = runtime.basis;
    runtime.productBasis = runtime.basis;
    runtime.productArcVersions.clear();
    runtime.productSignalArcVersions.clear();
    runtime.productCoreReceivers.clear();
    for (const auto& edge : runtime.basis.treeEdges)
    {
        const auto history = runtime.edgeHistory.find(edge);
        if (history != runtime.edgeHistory.end())
            runtime.productArcVersions.emplace(edge, history->second.arcVersion);
        const auto signals = runtime.phaseArcVersions.find(edge);
        if (signals != runtime.phaseArcVersions.end())
            runtime.productSignalArcVersions.emplace(edge, signals->second);
    }
    runtime.productDatumVersion = runtime.committedSGaugeVersion;
    runtime.productInitialized = runtime.initialized;
    runtime.productAuthorityComponentVersion = runtime.initialized
        ? runtime.integerComponentVersion : -1;
    runtime.requiresFreshProductIdentity = false;
}

void appendIdentityToken(std::ostringstream& out, const string& token)
{
    out << token.size() << ':' << token;
}

string committedSBasisIdentity(const ZhangGraphIntegerContext& context)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    auto token = [&](const auto& value)
    {
        std::ostringstream field;
        field.imbue(std::locale::classic());
        field << value;
        appendIdentityToken(out, field.str());
    };
    token("COMMITTED_S_SNAPSHOT_V1");
    token(context.productBasisMode); token(context.runtimeIdentity);
    token(context.snapshotEpoch); token(context.componentId);
    token(context.eventId); token(context.representationVersion);
    token(context.floatGaugeVersion); token(context.integerComponentVersion);
    token(context.productDatumVersion); token(context.sourceComponentCount);
    token(context.committedSGaugeVersion);
    const auto& basis = context.basis;
    token(basis.rootReceiver); token(basis.componentCount); token(basis.connected);
    token(basis.receivers.size());
    for (const auto& receiver : basis.receivers) token(receiver);
    token(basis.satellites.size());
    for (const auto& satellite : basis.satellites)
    { token(int(satellite.sys)); token(int(satellite.prn)); }
    auto edgeToken = [&](const ZhangGraphEdge& edge)
    { token(edge.receiver); token(int(edge.satellite.sys)); token(int(edge.satellite.prn)); };
    token(basis.edges.size()); for (const auto& edge : basis.edges) edgeToken(edge);
    token(basis.treeEdges.size()); for (const auto& edge : basis.treeEdges) edgeToken(edge);
    token(context.arcVersions.size());
    for (const auto& [edge, serial] : context.arcVersions)
    { edgeToken(edge); token(serial); }
    token(context.signalArcVersions.size());
    for (const auto& [edge, signals] : context.signalArcVersions)
    {
        edgeToken(edge); token(signals.size());
        for (const auto& [code, serial] : signals) { token(code); token(serial); }
    }
    token(context.signalObservationFinalized); token(context.signalObservationEdges.size());
    for (const auto& [code, edges] : context.signalObservationEdges)
    { token(code); token(edges.size()); for (const auto& edge : edges) edgeToken(edge); }
    return out.str();
}

struct SelectedGraphContext
{
    std::uint64_t token = 0;
    std::shared_ptr<const ZhangGraphIntegerContext> context;
};
std::mutex selectedGraphContextMutex;
std::uint64_t nextSelectedGraphContextToken = 0;
map<ZhangGraphRuntimeKey, vector<SelectedGraphContext>> selectedGraphContexts;

std::shared_ptr<const ZhangGraphIntegerContext> selectedGraphContext(
    const KFState& state, E_Sys system)
{
    const string runtime = zhangGraphRuntimeId(state);
    std::lock_guard<std::mutex> lock(selectedGraphContextMutex);
    const auto found = selectedGraphContexts.find({runtime, system});
    if (found == selectedGraphContexts.end() || found->second.empty()) return {};
    const auto& context = found->second.back().context;
    // A branch may inherit a selected view, but no view crosses an epoch.
    if (context->snapshotEpoch != state.time.to_string(9) ||
        context->productBasisMode != configuredProductBasisMode()) return {};
    return context;
}

string zhangGraphRuntimeId(const KFState& state);

using GraphProductReceiptKey =
    std::tuple<string, string, E_Sys, E_ObsCode, SatSys>;
map<GraphProductReceiptKey, ZhangGraphProductEventReceipt>
    graphProductEventReceipts;
map<ZhangGraphRuntimeKey, map<ZhangGraphEdge, int>>
    graphEventPreArcVersions;
map<ZhangGraphRuntimeKey, map<ZhangGraphEdge, map<int, int>>>
    graphEventPreSignalArcVersions;

// All callers are pre-AR graph transactions on the authoritative FLOAT root.
// The metadata travels with that KFState (including checkpoint and disposable
// copies), not in a process-global registry that an AR branch could mutate.
bool limitGraphHandoffFloatCapture(
    Trace& trace, const KFState& state,
    map<KFKey, ZhangCoordinateRow>& transform,
    map<std::pair<E_ObsCode, ZhangGraphEdge>, KFKey>& auxiliaries,
    map<string, string>& projectedMetadata, string& failureReason)
{
    const auto dependencies = zhangHandoffFloatDependencies(state);
    for (const auto& [identity, key] : auxiliaries)
    {
        ZhangHandoffFloatIdentity physical;
        if (!transform.contains(key) ||
            !zhangHandoffFloatIdentityFromKey(key, physical) ||
            physical.edge != identity.second || physical.observable != int(identity.first))
        {
            failureReason = "HANDOFF_FLOAT_CAPTURE_DESTINATION_IDENTITY_MISMATCH";
            trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time=" << state.time.to_string(0)
                  << " status=REJECTED reason=" << failureReason
                  << " observable=" << int(identity.first) << " key_num=" << key.num
                  << " edge=" << identity.second.receiver << '/' << identity.second.satellite.id()
                  << " posterior_preserved=1 integer_certification=0";
            return false;
        }
    }
    const auto plan = zhangPlanHandoffFloatCapture(state, transform);
    if (!plan.valid)
    {
        failureReason = plan.failureReason;
        trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time=" << state.time.to_string(0)
              << " status=REJECTED reason=" << plan.failureReason
              << " posterior_preserved=1 integer_certification=0";
        return false;
    }
    for (const auto& key : plan.discarded)
    {
        ZhangHandoffFloatIdentity physical;
        zhangHandoffFloatIdentityFromKey(key, physical);
        transform.erase(key);
        trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time="
              << state.time.to_string(0) << " status=NOT_RETAINED key=" << key
              << " edge=" << physical.edge.receiver << '/' << physical.edge.satellite.id()
              << " observable=" << physical.observable << " key_num=" << key.num
              << " arc_version=" << physical.arcVersion
              << " active_product_dependency=" << dependencies.contains(physical)
              << " reason=" << plan.discardReasons.at(key)
              << " canonical_key=" << (plan.aliases.contains(key) ?
                  plan.aliases.at(key).commaString() : "NONE")
              << " physical_dependency_preserved=" << plan.aliases.contains(key)
              << " integer_certification=0";
    }
    for (auto it = auxiliaries.begin(); it != auxiliaries.end();)
        if (plan.aliases.contains(it->second))
        { it->second = plan.aliases.at(it->second); ++it; }
        else it = plan.discarded.contains(it->second) ? auxiliaries.erase(it) : std::next(it);
    projectedMetadata = plan.metadata;
    trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time=" << state.time.to_string(0)
          << " status=PLANNED active_retained=" << plan.active
          << " inactive_retained=" << plan.inactive
          << " inactive_cache_budget=" << ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES
          << " active_above_cache_budget=" << (plan.active > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES)
          << " active_identity_directory_bound=" << ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES
          << " integer_certification=0";
    return true;
}

void storeGraphProductEventReceipts(
    const KFState& state,
    E_Sys system,
    vector<ZhangGraphProductEventReceipt> receipts)
{
    const string runtimeId = zhangGraphRuntimeId(state);
    const string epoch = state.time.to_string(9);
    for (auto it = graphProductEventReceipts.begin();
         it != graphProductEventReceipts.end();)
    {
        if (std::get<0>(it->first) == runtimeId
            && std::get<2>(it->first) == system
            && std::get<1>(it->first) != epoch)
            it = graphProductEventReceipts.erase(it);
        else
            ++it;
    }
    for (auto& receipt : receipts)
    {
        const GraphProductReceiptKey identity = {
            runtimeId, epoch, system, receipt.observable,
            receipt.satellite};
        // A second state transform in the same epoch requires composition of
        // both affine maps.  Until such composition exists, fail closed rather
        // than letting the last receipt masquerade as the whole handoff.
        if (graphProductEventReceipts.contains(identity))
        {
            receipt.conditionalProductMapProven = false;
            receipt.sameTimeMapProven = false;
            receipt.actualWriterRowMatches = false;
            receipt.representableSameTimeMapProven = false;
            receipt.fullJointMomentTransportProven = false;
            receipt.failureReason = "MULTIPLE_GRAPH_EVENTS_UNCOMPOSED";
            receipt.conditionalMapReason =
                "MULTIPLE_GRAPH_EVENTS_UNCOMPOSED";
            receipt.jointMomentReason =
                "MULTIPLE_GRAPH_EVENTS_UNCOMPOSED";
        }
        graphProductEventReceipts[identity] = std::move(receipt);
    }
}

void invalidateZhangProduct(GraphRuntimeState& runtime, const string& reason)
{
    if (runtime.productInitialized || !runtime.productBasis.edges.empty() ||
        !runtime.productArcVersions.empty() ||
        !runtime.productSignalArcVersions.empty() ||
        !runtime.productCoreReceivers.empty())
        runtime.integerComponentVersion++;
    runtime.productInitialized = false;
    runtime.productBasis = {};
    runtime.productArcVersions.clear();
    runtime.productSignalArcVersions.clear();
    runtime.productCoreReceivers.clear();
    runtime.productAuthorityComponentVersion = -1;
    runtime.lastProductEventCause = reason;
}


string zhangGraphRuntimeId(const KFState& state)
{
    auto branch = state.metaDataMap.find(
        ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA);
    if (branch != state.metaDataMap.end())
    {
        return branch->second;
    }
    return zhangCheckpointRuntimeId(state);
}

bool validZhangGraphRuntimeId(const string& runtimeId)
{
    return !runtimeId.empty() && runtimeId.size() <= 512 &&
        std::all_of(runtimeId.begin(), runtimeId.end(), [](unsigned char ch)
        {
            return ch >= 0x21 && ch <= 0x7e;
        });
}

set<E_Sys> zhangConfiguredSpanningTreeSystems()
{
    set<E_Sys> systems;
    if (!acsConfig.zhangFullRank.enable) return systems;
    for (const auto& [system, options] : acsConfig.zhangFullRank.sysOpts)
    {
        const auto processed = acsConfig.process_sys.find(system);
        if (system != E_Sys::NONE && options.use_spanning_tree &&
            processed != acsConfig.process_sys.end() && processed->second)
            systems.insert(system);
    }
    return systems;
}

bool zhangHasEstimatedGraphPhaseState(
    const KFState& state, E_Sys system, const vector<E_ObsCode>& observables)
{
    for (const auto& [key, index] : state.kfIndexMap)
        if (key.Sat.sys == system &&
            ((key.type == KF::AMBIGUITY && key.Sat.prn > 0 && !key.str.empty()) ||
             (key.type == KF::PHASE_BIAS &&
              ((key.Sat.prn == 0 && !key.str.empty()) ||
               (key.Sat.prn > 0 && key.str.empty())))) &&
            zhangFullRankUsesObservable(static_cast<E_ObsCode>(key.num), observables))
            return true;
    return false;
}

constexpr char ZHANG_GRAPH_CHECKPOINT_MAGIC[] =
    "GINAN_ZHANG_GRAPH_RUNTIME";
constexpr char ZHANG_GRAPH_V4_PAYLOAD_PREFIX[] =
    "GINAN_ZHANG_GRAPH_V4\n";
constexpr char ZHANG_GRAPH_V5_PAYLOAD_PREFIX[] =
    "GINAN_ZHANG_GRAPH_V5\n";

struct ZhangGraphCheckpointSatellite
{
    int system = static_cast<int>(E_Sys::NONE);
    int prn = 0;

    bool operator<(const ZhangGraphCheckpointSatellite& other) const
    {
        return std::tie(system, prn) < std::tie(other.system, other.prn);
    }

    bool operator==(const ZhangGraphCheckpointSatellite& other) const
    {
        return system == other.system && prn == other.prn;
    }

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & system;
        archive & prn;
    }
};

struct ZhangGraphCheckpointEdge
{
    string receiver;
    ZhangGraphCheckpointSatellite satellite;

    bool operator<(const ZhangGraphCheckpointEdge& other) const
    {
        return std::tie(receiver, satellite) <
               std::tie(other.receiver, other.satellite);
    }

    bool operator==(const ZhangGraphCheckpointEdge& other) const
    {
        return receiver == other.receiver && satellite == other.satellite;
    }

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & receiver;
        archive & satellite;
    }
};

struct ZhangGraphCheckpointBasis
{
    string                          rootReceiver;
    set<ZhangGraphCheckpointEdge>   edges;
    set<ZhangGraphCheckpointEdge>   treeEdges;
    set<string>                     receivers;
    set<ZhangGraphCheckpointSatellite> satellites;
    int                             componentCount = 0;
    bool                            connected = false;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & rootReceiver;
        archive & edges;
        archive & treeEdges;
        archive & receivers;
        archive & satellites;
        archive & componentCount;
        archive & connected;
    }
};

struct ZhangGraphCheckpointEdgeHistory
{
    int continuousEpochs = 0;
    int outageEpochs = 0;
    int arcVersion = 0;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & continuousEpochs;
        archive & outageEpochs;
        archive & arcVersion;
    }
};

struct ZhangGraphCheckpointArcObservationHistory
{
    int         firstObservedEpoch = -1;
    int         lastObservedEpoch = -1;
    int         observationEpochs = 0;
    vector<int> observedEpochs;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & firstObservedEpoch;
        archive & lastObservedEpoch;
        archive & observationEpochs;
        archive & observedEpochs;
    }
};

struct ZhangGraphCheckpointOutageState
{
    int system = static_cast<int>(E_Sys::NONE);
    int receiverEpochs = 0;
    int satelliteEpochs = 0;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & system;
        archive & receiverEpochs;
        archive & satelliteEpochs;
    }
};

struct ZhangGraphCheckpointRuntimeState
{
    int system = static_cast<int>(E_Sys::NONE);
    ZhangGraphCheckpointBasis basis;
    ZhangGraphCheckpointBasis activeBasis;
    ZhangGraphCheckpointBasis productBasis;
    map<ZhangGraphCheckpointEdge, int> productArcVersions;
    set<string> productCoreReceivers;
    set<ZhangGraphCheckpointEdge> observationEdges;
    set<ZhangGraphCheckpointEdge> stateEdges;
    map<ZhangGraphCheckpointEdge, ZhangGraphCheckpointEdgeHistory> edgeHistory;
    map<pair<ZhangGraphCheckpointEdge, int>,
        ZhangGraphCheckpointArcObservationHistory>
        arcVersionObservationHistory;
    int    epochIndex = -1;
    string lastProductEventCause = "INITIALISE";
    bool   initialized = false;
    int    deferredEpochs = 0;
    int    datumVersion = 0;
    int    representationVersion = 0;
    int    floatGaugeVersion = 0;
    int    integerComponentVersion = 0;
    int    eventCounter = 0;
    int    productDatumVersion = 0;
    bool   productInitialized = false;
    int productAuthorityComponentVersion = -1;
    map<ZhangGraphCheckpointEdge, map<int, int>> phaseArcVersions;
    map<ZhangGraphCheckpointEdge, map<int, int>> phaseAdmissions;
    string phaseAdmissionTime;
    int phaseAdmissionChart = -1;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & system;
        archive & basis;
        archive & activeBasis;
        archive & productBasis;
        archive & productArcVersions;
        archive & productCoreReceivers;
        archive & observationEdges;
        archive & stateEdges;
        archive & edgeHistory;
        archive & arcVersionObservationHistory;
        archive & epochIndex;
        archive & lastProductEventCause;
        archive & initialized;
        archive & deferredEpochs;
        archive & datumVersion;
        archive & representationVersion;
        archive & floatGaugeVersion;
        archive & integerComponentVersion;
        archive & eventCounter;
        archive & productDatumVersion;
        archive & productInitialized;
        archive & productAuthorityComponentVersion;
        archive & phaseArcVersions;
        archive & phaseAdmissions;
        archive & phaseAdmissionTime;
        archive & phaseAdmissionChart;
    }
};

struct ZhangGraphCheckpointPayload
{
    string magic = ZHANG_GRAPH_CHECKPOINT_MAGIC;
    std::uint32_t schemaVersion = 3;
    string runtimeId;
    vector<ZhangGraphCheckpointOutageState> outageStates;
    vector<ZhangGraphCheckpointRuntimeState> graphStates;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & magic;
        archive & schemaVersion;
        archive & runtimeId;
        archive & outageStates;
        archive & graphStates;
    }
};

// The schema-3 Boost field layout above is intentionally frozen. Appending
// fields to it would make the old 01:59:30 archive ambiguous under no-header
// deserialisation. Schema 4 wraps the old layout and adds signal identities.
struct ZhangGraphCheckpointRuntimeV4
{
    ZhangGraphCheckpointRuntimeState common;
    map<ZhangGraphCheckpointEdge, map<int, int>>
        productSignalArcVersions;
    map<pair<ZhangGraphCheckpointEdge, pair<int, int>>,
        ZhangGraphCheckpointArcObservationHistory>
        signalArcVersionObservationHistory;
    map<int, set<ZhangGraphCheckpointEdge>> signalObservationEdges;
    bool requiresFreshProductIdentity = false;
    bool signalIdentityAuthoritative = false;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & common;
        archive & productSignalArcVersions;
        archive & signalArcVersionObservationHistory;
        archive & signalObservationEdges;
        archive & requiresFreshProductIdentity;
        archive & signalIdentityAuthoritative;
    }
};

struct ZhangGraphCheckpointPayloadV4
{
    string magic = ZHANG_GRAPH_CHECKPOINT_MAGIC;
    std::uint32_t schemaVersion = 4; // frozen V4 wire layout, NOT latest schema
    string runtimeId;
    vector<ZhangGraphCheckpointOutageState> outageStates;
    vector<ZhangGraphCheckpointRuntimeV4> graphStates;

    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int& version)
    {
        archive & magic;
        archive & schemaVersion;
        archive & runtimeId;
        archive & outageStates;
        archive & graphStates;
    }
};

// V5 adds an explicit algorithm contract around the frozen V4 layout.  A
// legacy two-tree archive is never relabelled as a single-S estimator.
struct ZhangGraphCheckpointPayloadV5
{
    string magic = ZHANG_GRAPH_CHECKPOINT_MAGIC;
    std::uint32_t schemaVersion = ZHANG_GRAPH_CHECKPOINT_SCHEMA_VERSION;
    string productBasisMode;
    map<int, int> committedSGaugeVersions;
    ZhangGraphCheckpointPayloadV4 graph;
    template <class ARCHIVE>
    void serialize(ARCHIVE& archive, const unsigned int&)
    {
        archive & magic;
        archive & schemaVersion;
        archive & productBasisMode;
        archive & committedSGaugeVersions;
        archive & graph;
    }
};

bool zhangGraphCheckpointSystemValid(int system)
{
    return system >= static_cast<int>(E_Sys::NONE) &&
           system <= static_cast<int>(E_Sys::COMB);
}

ZhangGraphCheckpointSatellite zhangGraphCheckpointSatellite(
    const SatSys& satellite)
{
    return {
        static_cast<int>(satellite.sys),
        static_cast<int>(satellite.prn)};
}

SatSys zhangGraphCheckpointSatellite(
    const ZhangGraphCheckpointSatellite& satellite)
{
    return SatSys(static_cast<E_Sys>(satellite.system), satellite.prn);
}

ZhangGraphCheckpointEdge zhangGraphCheckpointEdge(
    const ZhangGraphEdge& edge)
{
    return {edge.receiver, zhangGraphCheckpointSatellite(edge.satellite)};
}

ZhangGraphEdge zhangGraphCheckpointEdge(
    const ZhangGraphCheckpointEdge& edge)
{
    return {edge.receiver, zhangGraphCheckpointSatellite(edge.satellite)};
}

ZhangGraphCheckpointBasis zhangGraphCheckpointBasis(
    const ZhangGraphBasis& basis)
{
    ZhangGraphCheckpointBasis snapshot;
    snapshot.rootReceiver = basis.rootReceiver;
    for (const auto& edge : basis.edges)
    {
        snapshot.edges.insert(zhangGraphCheckpointEdge(edge));
    }
    for (const auto& edge : basis.treeEdges)
    {
        snapshot.treeEdges.insert(zhangGraphCheckpointEdge(edge));
    }
    snapshot.receivers = basis.receivers;
    for (const auto& satellite : basis.satellites)
    {
        snapshot.satellites.insert(
            zhangGraphCheckpointSatellite(satellite));
    }
    snapshot.componentCount = basis.componentCount;
    snapshot.connected = basis.connected;
    return snapshot;
}

ZhangGraphBasis zhangGraphCheckpointBasis(
    const ZhangGraphCheckpointBasis& snapshot)
{
    ZhangGraphBasis basis;
    basis.rootReceiver = snapshot.rootReceiver;
    for (const auto& edge : snapshot.edges)
    {
        basis.edges.insert(zhangGraphCheckpointEdge(edge));
    }
    for (const auto& edge : snapshot.treeEdges)
    {
        basis.treeEdges.insert(zhangGraphCheckpointEdge(edge));
    }
    basis.receivers = snapshot.receivers;
    for (const auto& satellite : snapshot.satellites)
    {
        basis.satellites.insert(
            zhangGraphCheckpointSatellite(satellite));
    }
    basis.componentCount = snapshot.componentCount;
    basis.connected = snapshot.connected;
    return basis;
}

bool zhangGraphCheckpointEdgeValid(
    const ZhangGraphCheckpointEdge& edge,
    int                              system,
    string&                          failureReason,
    const string&                    field)
{
    if (edge.receiver.empty())
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_EMPTY_RECEIVER:" + field;
        return false;
    }
    if (!zhangGraphCheckpointSystemValid(edge.satellite.system) ||
        edge.satellite.system != system || edge.satellite.prn <= 0)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_INVALID_SATELLITE:" + field;
        return false;
    }
    return true;
}

bool zhangGraphCheckpointBasisValid(
    const ZhangGraphCheckpointBasis& snapshot,
    int                               system,
    string&                           failureReason,
    const string&                     field)
{
    if (snapshot.componentCount < 0)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_NEGATIVE_COMPONENTS:" + field;
        return false;
    }

    set<string> derivedReceivers;
    set<ZhangGraphCheckpointSatellite> derivedSatellites;
    for (const auto& edge : snapshot.edges)
    {
        if (!zhangGraphCheckpointEdgeValid(
                edge, system, failureReason, field + ".edges"))
        {
            return false;
        }
        derivedReceivers.insert(edge.receiver);
        derivedSatellites.insert(edge.satellite);
    }
    for (const auto& edge : snapshot.treeEdges)
    {
        if (!zhangGraphCheckpointEdgeValid(
                edge, system, failureReason, field + ".treeEdges") ||
            snapshot.edges.find(edge) == snapshot.edges.end())
        {
            if (failureReason.empty())
            {
                failureReason =
                    "ZHANG_GRAPH_CHECKPOINT_TREE_EDGE_OUTSIDE_GRAPH:" + field +
                    ":" + edge.receiver + ":" +
                    zhangGraphCheckpointSatellite(edge.satellite).id();
            }
            return false;
        }
    }
    for (const auto& satellite : snapshot.satellites)
    {
        if (!zhangGraphCheckpointSystemValid(satellite.system) ||
            satellite.system != system || satellite.prn <= 0)
        {
            failureReason =
                "ZHANG_GRAPH_CHECKPOINT_INVALID_BASIS_SATELLITE:" + field;
            return false;
        }
    }
    if (derivedReceivers != snapshot.receivers ||
        derivedSatellites != snapshot.satellites)
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_BASIS_NODE_SET_MISMATCH:" + field;
        return false;
    }

    if (snapshot.edges.empty())
    {
        if (!snapshot.treeEdges.empty() || !snapshot.receivers.empty() ||
            !snapshot.satellites.empty() || snapshot.componentCount != 0 ||
            snapshot.connected)
        {
            failureReason =
                "ZHANG_GRAPH_CHECKPOINT_INVALID_EMPTY_BASIS:" + field;
            return false;
        }
        return true;
    }
    // A FLOAT forest may retain an external historical reference. Its local
    // roots are uncertain states, not new hard-zero datums. Product charts
    // still require their reference to be a node of the authorized graph.
    if (snapshot.rootReceiver.empty() ||
        (field == "productBasis" && !snapshot.receivers.contains(snapshot.rootReceiver)))
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_INVALID_BASIS_ROOT:" + field;
        return false;
    }

    const ZhangGraphBasis basis = zhangGraphCheckpointBasis(snapshot);
    const ZhangGraphBasis derived =
        zhangBuildSpanningTree(basis.edges, basis.rootReceiver);
    if (derived.componentCount != basis.componentCount ||
        derived.connected != basis.connected)
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_BASIS_CONNECTIVITY_MISMATCH:" + field;
        return false;
    }

    const ZhangGraphBasis derivedTree =
        zhangBuildSpanningTree(basis.treeEdges, basis.rootReceiver);
    const int nodeCount = static_cast<int>(
        basis.receivers.size() + basis.satellites.size());
    if (derivedTree.receivers != basis.receivers ||
        derivedTree.satellites != basis.satellites ||
        derivedTree.componentCount != basis.componentCount ||
        static_cast<int>(basis.treeEdges.size()) !=
            nodeCount - basis.componentCount)
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_INVALID_BASIS_FOREST:" + field;
        return false;
    }
    return true;
}

bool zhangGraphCheckpointRuntimeValid(
    const ZhangGraphCheckpointRuntimeState& snapshot,
    string&                                  failureReason,
    bool singleS = false)
{
    failureReason.clear();
    if (!zhangGraphCheckpointSystemValid(snapshot.system))
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_INVALID_SYSTEM";
        return false;
    }
    if (!zhangGraphCheckpointBasisValid(
            snapshot.basis, snapshot.system, failureReason, "basis") ||
        !zhangGraphCheckpointBasisValid(
            snapshot.activeBasis,
            snapshot.system,
            failureReason,
            "activeBasis") ||
        !zhangGraphCheckpointBasisValid(
            snapshot.productBasis,
            snapshot.system,
            failureReason,
            singleS ? "committedSProductView" : "productBasis"))
    {
        return false;
    }
    if (snapshot.epochIndex < -1 || snapshot.deferredEpochs < 0 ||
        snapshot.datumVersion < 0 || snapshot.eventCounter < 0 ||
        snapshot.productDatumVersion < 0 ||
        snapshot.representationVersion < 0 ||
        snapshot.floatGaugeVersion < 0 ||
        snapshot.integerComponentVersion < 0)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_NEGATIVE_COUNTER";
        return false;
    }
    if (snapshot.floatGaugeVersion != snapshot.datumVersion ||
        (!singleS && snapshot.integerComponentVersion < snapshot.productDatumVersion) ||
        (snapshot.productInitialized && snapshot.productAuthorityComponentVersion !=
            snapshot.integerComponentVersion) ||
        (!snapshot.productInitialized && snapshot.productAuthorityComponentVersion != -1))
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_VERSION_SEMANTICS_MISMATCH";
        return false;
    }
    if (singleS &&
        (!zhangSameCommittedGraphBasis(zhangGraphCheckpointBasis(snapshot.basis),
             zhangGraphCheckpointBasis(snapshot.activeBasis)) ||
         !zhangSameCommittedGraphBasis(zhangGraphCheckpointBasis(snapshot.basis),
             zhangGraphCheckpointBasis(snapshot.productBasis)) ||
         snapshot.productInitialized != snapshot.initialized ||
         !snapshot.productCoreReceivers.empty()))
    {
        failureReason = "ZHANG_SINGLE_S_CHECKPOINT_HAS_INDEPENDENT_GRAPH_VIEW";
        return false;
    }
    if (snapshot.lastProductEventCause.empty())
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_EMPTY_EVENT_CAUSE";
        return false;
    }

    for (const auto& [edge, signals] : snapshot.phaseArcVersions)
    {
        if (!zhangGraphCheckpointEdgeValid(edge, snapshot.system, failureReason, "phaseArcVersions"))
            return false;
        for (const auto& [code, serial] : signals)
            if (code <= 0 || serial < 0)
            { failureReason = "ZHANG_GRAPH_INVALID_SIGNAL_ARC"; return false; }
    }
    for (const auto& [edge, signals] : snapshot.phaseAdmissions)
    {
        const auto arcs = snapshot.phaseArcVersions.find(edge);
        if (!snapshot.observationEdges.contains(edge) || arcs == snapshot.phaseArcVersions.end() ||
            snapshot.phaseAdmissionTime.empty() || snapshot.phaseAdmissionChart != snapshot.representationVersion)
        { failureReason = "ZHANG_GRAPH_STALE_PHASE_ADMISSION"; return false; }
        for (const auto& [code, serial] : signals)
            if (!arcs->second.contains(code) || arcs->second.at(code) != serial)
            { failureReason = "ZHANG_GRAPH_PHASE_ARC_IDENTITY_MISMATCH"; return false; }
    }
    for (const auto& edge : snapshot.observationEdges)
    {
        if (!zhangGraphCheckpointEdgeValid(
                edge,
                snapshot.system,
                failureReason,
                "observationEdges"))
        {
            return false;
        }
    }
    for (const auto& edge : snapshot.stateEdges)
    {
        if (!zhangGraphCheckpointEdgeValid(
                edge, snapshot.system, failureReason, "stateEdges"))
        {
            return false;
        }
    }
    if (!std::includes(
            snapshot.stateEdges.begin(),
            snapshot.stateEdges.end(),
            snapshot.observationEdges.begin(),
            snapshot.observationEdges.end()))
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_OBSERVATION_OUTSIDE_STATE_EDGES";
        return false;
    }

    for (const auto& [edge, history] : snapshot.edgeHistory)
    {
        if (!zhangGraphCheckpointEdgeValid(
                edge, snapshot.system, failureReason, "edgeHistory") ||
            history.continuousEpochs < 0 || history.outageEpochs < 0 ||
            history.arcVersion < 0)
        {
            if (failureReason.empty())
            {
                failureReason =
                    "ZHANG_GRAPH_CHECKPOINT_INVALID_EDGE_HISTORY";
            }
            return false;
        }
    }

    for (const auto& [versionedEdge, history] :
         snapshot.arcVersionObservationHistory)
    {
        const auto& edge = versionedEdge.first;
        const int arcVersion = versionedEdge.second;
        if (!zhangGraphCheckpointEdgeValid(
                edge,
                snapshot.system,
                failureReason,
                "arcVersionObservationHistory"))
        {
            return false;
        }
        auto current = snapshot.edgeHistory.find(edge);
        if (arcVersion < 0 || current == snapshot.edgeHistory.end() ||
            arcVersion > current->second.arcVersion ||
            history.observationEpochs < 0 ||
            history.observationEpochs !=
                static_cast<int>(history.observedEpochs.size()))
        {
            failureReason =
                "ZHANG_GRAPH_CHECKPOINT_INVALID_ARC_OBSERVATION_HISTORY";
            return false;
        }
        if (history.observedEpochs.empty())
        {
            if (history.firstObservedEpoch != -1 ||
                history.lastObservedEpoch != -1)
            {
                failureReason =
                    "ZHANG_GRAPH_CHECKPOINT_EMPTY_ARC_HISTORY_BOUNDS";
                return false;
            }
            continue;
        }
        if (history.firstObservedEpoch != history.observedEpochs.front() ||
            history.lastObservedEpoch != history.observedEpochs.back() ||
            history.firstObservedEpoch < 0 ||
            history.lastObservedEpoch > snapshot.epochIndex ||
            !std::is_sorted(
                history.observedEpochs.begin(),
                history.observedEpochs.end()) ||
            std::adjacent_find(
                history.observedEpochs.begin(),
                history.observedEpochs.end()) !=
                history.observedEpochs.end())
        {
            failureReason =
                "ZHANG_GRAPH_CHECKPOINT_INCONSISTENT_ARC_EPOCHS";
            return false;
        }
    }

    for (const auto& [edge, arcVersion] : snapshot.productArcVersions)
    {
        auto current = snapshot.edgeHistory.find(edge);
        auto versionHistory = snapshot.arcVersionObservationHistory.find(
            {edge, arcVersion});
        if (!zhangGraphCheckpointEdgeValid(
                edge,
                snapshot.system,
                failureReason,
                "productArcVersions") ||
            snapshot.productBasis.treeEdges.find(edge) ==
                snapshot.productBasis.treeEdges.end() ||
            arcVersion < 0 || current == snapshot.edgeHistory.end() ||
            arcVersion > current->second.arcVersion ||
            (singleS && arcVersion != current->second.arcVersion) ||
            (!singleS && versionHistory == snapshot.arcVersionObservationHistory.end()))
        {
            if (failureReason.empty())
            {
                failureReason =
                    "ZHANG_GRAPH_CHECKPOINT_INVALID_PRODUCT_ARC_VERSION";
            }
            failureReason += ":edge=" + edge.receiver + ":" +
                zhangGraphCheckpointSatellite(edge.satellite).id() +
                ":product_arc=" + std::to_string(arcVersion) +
                ":current_arc=" + (current == snapshot.edgeHistory.end()
                    ? string("MISSING") : std::to_string(current->second.arcVersion)) +
                ":tree_present=" + std::to_string(snapshot.productBasis.treeEdges.contains(edge)) +
                ":version_history_present=" +
                    std::to_string(versionHistory != snapshot.arcVersionObservationHistory.end());
            return false;
        }
    }
    if (snapshot.productInitialized)
    {
        if ((!singleS && !snapshot.productBasis.connected) ||
            snapshot.productArcVersions.size() !=
                snapshot.productBasis.treeEdges.size())
        {
            failureReason =
                "ZHANG_GRAPH_CHECKPOINT_INCOMPLETE_PRODUCT_DATUM";
            return false;
        }
    }
    else if (!snapshot.productArcVersions.empty() ||
             !snapshot.productCoreReceivers.empty() || !snapshot.productBasis.edges.empty())
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_UNINITIALIZED_PRODUCT_HAS_STATE";
        return false;
    }
    if (!std::includes(
            snapshot.productBasis.receivers.begin(),
            snapshot.productBasis.receivers.end(),
            snapshot.productCoreReceivers.begin(),
            snapshot.productCoreReceivers.end()))
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_CORE_OUTSIDE_PRODUCT_RECEIVERS";
        return false;
    }
    return true;
}

bool zhangGraphCheckpointPayloadValid(
    const ZhangGraphCheckpointPayload& snapshot,
    const string&                       runtimeId,
    string&                             failureReason)
{
    if (snapshot.magic != ZHANG_GRAPH_CHECKPOINT_MAGIC)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_MAGIC_MISMATCH";
        return false;
    }
    if (snapshot.schemaVersion != 3)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_SCHEMA_MISMATCH";
        return false;
    }
    if (snapshot.runtimeId.empty() || snapshot.runtimeId != runtimeId)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_RUNTIME_ID_MISMATCH";
        return false;
    }

    set<int> outageSystems;
    for (const auto& outage : snapshot.outageStates)
    {
        if (!zhangGraphCheckpointSystemValid(outage.system) ||
            !outageSystems.insert(outage.system).second ||
            outage.receiverEpochs < 0 || outage.satelliteEpochs < 0)
        {
            failureReason = "ZHANG_GRAPH_CHECKPOINT_INVALID_OUTAGE_STATE";
            return false;
        }
    }

    set<int> graphSystems;
    for (const auto& runtime : snapshot.graphStates)
    {
        if (!graphSystems.insert(runtime.system).second)
        {
            failureReason = "ZHANG_GRAPH_CHECKPOINT_DUPLICATE_GRAPH_SYSTEM";
            return false;
        }
        if (!zhangGraphCheckpointRuntimeValid(runtime, failureReason))
        {
            return false;
        }
    }
    return true;
}

bool zhangGraphCheckpointPayloadValid(
    const ZhangGraphCheckpointPayloadV4& snapshot,
    const string& runtimeId,
    string& failureReason,
    bool singleS = false)
{
    if (snapshot.magic != ZHANG_GRAPH_CHECKPOINT_MAGIC ||
        snapshot.schemaVersion != 4 ||
        snapshot.runtimeId.empty() || snapshot.runtimeId != runtimeId)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_V4_HEADER_INVALID";
        return false;
    }
    set<int> outageSystems;
    for (const auto& outage : snapshot.outageStates)
        if (!zhangGraphCheckpointSystemValid(outage.system) ||
            !outageSystems.insert(outage.system).second ||
            outage.receiverEpochs < 0 || outage.satelliteEpochs < 0)
        {
            failureReason = "ZHANG_GRAPH_CHECKPOINT_INVALID_OUTAGE_STATE";
            return false;
        }
    set<int> graphSystems;
    for (const auto& state : snapshot.graphStates)
    {
        const auto& common = state.common;
        if (!graphSystems.insert(common.system).second ||
            !zhangGraphCheckpointRuntimeValid(common, failureReason, singleS))
        {
            if (failureReason.empty())
                failureReason = "ZHANG_GRAPH_CHECKPOINT_DUPLICATE_GRAPH_SYSTEM";
            return false;
        }
        if (state.requiresFreshProductIdentity && common.productInitialized)
        {
            failureReason = "ZHANG_GRAPH_LEGACY_PRODUCT_IDENTITY_STILL_ACTIVE";
            return false;
        }
        if (singleS && state.requiresFreshProductIdentity)
        {
            failureReason = "ZHANG_SINGLE_S_CHECKPOINT_REQUIRES_EXPLICIT_MIGRATION";
            return false;
        }
        if (singleS && !state.signalIdentityAuthoritative)
        {
            failureReason = "ZHANG_SINGLE_S_CHECKPOINT_SIGNAL_IDENTITY_NOT_AUTHORITATIVE";
            return false;
        }
        if (state.signalIdentityAuthoritative &&
            common.productInitialized &&
            state.productSignalArcVersions.size() !=
                common.productBasis.treeEdges.size())
        {
            failureReason = "ZHANG_GRAPH_PRODUCT_SIGNAL_ARCS_INCOMPLETE";
            return false;
        }
        if ((!common.productInitialized ||
             !state.signalIdentityAuthoritative) &&
            !state.productSignalArcVersions.empty())
        {
            failureReason = "ZHANG_GRAPH_INACTIVE_PRODUCT_HAS_SIGNAL_ARCS";
            return false;
        }
        for (const auto& [edge, signals] :
             state.productSignalArcVersions)
        {
            const auto current = common.phaseArcVersions.find(edge);
            if (!common.productBasis.treeEdges.contains(edge) ||
                current == common.phaseArcVersions.end() ||
                signals.empty())
            {
                failureReason = "ZHANG_GRAPH_PRODUCT_SIGNAL_ARC_INVALID";
                return false;
            }
            for (const auto& [code, version] : signals)
                if (code <= 0 || version < 0 ||
                    !current->second.contains(code) ||
                    version > current->second.at(code))
                {
                    failureReason = "ZHANG_GRAPH_PRODUCT_SIGNAL_VERSION_INVALID";
                    return false;
                }
            if (singleS && signals != current->second)
            {
                failureReason = "ZHANG_SINGLE_S_PRODUCT_VIEW_SIGNAL_IDENTITY_MISMATCH";
                return false;
            }
        }
        for (const auto& [identity, history] :
             state.signalArcVersionObservationHistory)
        {
            const auto& edge = identity.first;
            const int code = identity.second.first;
            const int version = identity.second.second;
            const auto current = common.phaseArcVersions.find(edge);
            if (!zhangGraphCheckpointEdgeValid(edge, common.system,
                    failureReason, "signalArcVersionObservationHistory") ||
                code <= 0 || version < 0 ||
                current == common.phaseArcVersions.end() ||
                !current->second.contains(code) ||
                version > current->second.at(code) ||
                history.observationEpochs < 0 ||
                history.observationEpochs !=
                    int(history.observedEpochs.size()))
            {
                failureReason = "ZHANG_GRAPH_SIGNAL_ARC_HISTORY_INVALID";
                return false;
            }
        }
        for (const auto& [code, edges] : state.signalObservationEdges)
            for (const auto& edge : edges)
                if (code <= 0 || !common.basis.edges.contains(edge) ||
                    !common.stateEdges.contains(edge) ||
                    !common.observationEdges.contains(edge) ||
                    !common.phaseAdmissions.contains(edge) ||
                    !common.phaseAdmissions.at(edge).contains(code))
                {
                    failureReason = "ZHANG_GRAPH_SIGNAL_OBSERVATION_INVALID";
                    return false;
                }
    }
    return true;
}

template <typename TYPE>
bool zhangGraphCheckpointSerialize(
    const TYPE& value,
    string&     payload,
    string&     failureReason)
{
    try
    {
        std::ostringstream output(std::ios::binary | std::ios::out);
        boost::archive::binary_oarchive archive(
            output, boost::archive::no_header);
        archive << value;
        payload = output.str();
        return true;
    }
    catch (const std::exception& exception)
    {
        payload.clear();
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_SERIALIZE_FAILED:" +
            string(exception.what());
        return false;
    }
}

template <typename TYPE>
bool zhangGraphCheckpointDeserialize(
    const string& payload,
    TYPE&         value,
    string&       failureReason)
{
    try
    {
        std::istringstream input(
            payload, std::ios::binary | std::ios::in);
        {
            boost::archive::binary_iarchive archive(
                input, boost::archive::no_header);
            archive >> value;
        }
        if (input.peek() != std::char_traits<char>::eof())
        {
            failureReason = "ZHANG_GRAPH_CHECKPOINT_TRAILING_BYTES";
            return false;
        }
        return true;
    }
    catch (const std::exception& exception)
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_DESERIALIZE_FAILED:" +
            string(exception.what());
        return false;
    }
}

bool zhangGraphCheckpointDecodeAndValidate(
    const string&                       runtimeId,
    const string&                       payload,
    ZhangGraphCheckpointPayloadV4& snapshot,
    bool& legacyIdentityReset,
    string&                             failureReason)
{
    failureReason.clear();
    legacyIdentityReset = false;
    if (runtimeId.empty())
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_EMPTY_RUNTIME_ID";
        return false;
    }
    if (payload.empty())
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_EMPTY_PAYLOAD";
        return false;
    }
    const string mode = configuredProductBasisMode();
    const string v5prefix = ZHANG_GRAPH_V5_PAYLOAD_PREFIX;
    if (payload.starts_with(v5prefix))
    {
        ZhangGraphCheckpointPayloadV5 envelope;
        if (!zhangGraphCheckpointDeserialize(payload.substr(v5prefix.size()),
                envelope, failureReason)) return false;
        if (envelope.magic != ZHANG_GRAPH_CHECKPOINT_MAGIC ||
            envelope.schemaVersion != ZHANG_GRAPH_CHECKPOINT_SCHEMA_VERSION ||
            envelope.productBasisMode != mode)
        {
            failureReason = "ZHANG_GRAPH_CHECKPOINT_PRODUCT_BASIS_MODE_OR_V5_HEADER_MISMATCH";
            return false;
        }
        snapshot = std::move(envelope.graph);
        if (singleSBasisNative())
        {
            if (envelope.committedSGaugeVersions.size() != snapshot.graphStates.size())
            {
                failureReason = "ZHANG_SINGLE_S_CHECKPOINT_GAUGE_GENERATIONS_INCOMPLETE";
                return false;
            }
            for (const auto& state : snapshot.graphStates)
            {
                const auto gauge = envelope.committedSGaugeVersions.find(state.common.system);
                if (gauge == envelope.committedSGaugeVersions.end() || gauge->second < 0 ||
                    gauge->second != state.common.productDatumVersion)
                {
                    failureReason = "ZHANG_SINGLE_S_CHECKPOINT_PRODUCT_NOT_BOUND_TO_S_GAUGE";
                    return false;
                }
            }
        }
        else if (!envelope.committedSGaugeVersions.empty())
        {
            failureReason = "ZHANG_LEGACY_CHECKPOINT_HAS_NATIVE_S_GAUGE_METADATA";
            return false;
        }
        if (!zhangGraphCheckpointPayloadValid(snapshot, runtimeId, failureReason,
                singleSBasisNative())) return false;
        return true;
    }
    if (singleSBasisNative())
    {
        failureReason = "ZHANG_SINGLE_S_CHECKPOINT_REJECTED_LEGACY_V3_V4_REQUIRES_EXPLICIT_MIGRATION";
        return false;
    }
    const string prefix = ZHANG_GRAPH_V4_PAYLOAD_PREFIX;
    if (payload.starts_with(prefix))
    {
        if (!zhangGraphCheckpointDeserialize(
                payload.substr(prefix.size()), snapshot, failureReason) ||
            !zhangGraphCheckpointPayloadValid(
                snapshot, runtimeId, failureReason))
            return false;
        legacyIdentityReset = zhangR51SingleBasisHandoffEnabled() &&
            std::any_of(snapshot.graphStates.begin(),
                snapshot.graphStates.end(),
                [](const auto& graph)
                { return !graph.signalIdentityAuthoritative; });
        return true;
    }
    // Unprefixed Boost archives are the frozen schema-3 layout. Its own
    // magic/schema checks below reject any unrelated or corrupt payload.
    ZhangGraphCheckpointPayload legacy;
    if (!zhangGraphCheckpointDeserialize(
            payload, legacy, failureReason) ||
        !zhangGraphCheckpointPayloadValid(
            legacy, runtimeId, failureReason))
        return false;
    snapshot = {};
    snapshot.runtimeId = legacy.runtimeId;
    snapshot.outageStates = std::move(legacy.outageStates);
    for (auto& state : legacy.graphStates)
    {
        ZhangGraphCheckpointRuntimeV4 migrated;
        migrated.common = std::move(state);
        migrated.requiresFreshProductIdentity = true;
        snapshot.graphStates.push_back(std::move(migrated));
    }
    legacyIdentityReset = true;
    return true;
}

bool slipIsExcluded(const SigStat::SlipStat& slip)
{
    if (!slip.any)
    {
        return false;
    }

    return
        (acsConfig.exclude.LLI         && slip.LLI)        ||
        (acsConfig.exclude.GF          && slip.GF)         ||
        (acsConfig.exclude.MW          && slip.MW)         ||
        (acsConfig.exclude.SCDIA       && slip.SCDIA)      ||
        (acsConfig.exclude.retrack     && slip.retrack)    ||
        (acsConfig.exclude.single_freq && slip.singleFreq);
}

bool signalIsUsable(const GObs& obs, E_ObsCode code)
{
    for (auto& [frequency, signal] : obs.sigs)
    {
        if (signal.code != code || signal.L == 0 || signal.invalid || (signal.LLI & 2))
        {
            continue;
        }

        if (obs.satStat_ptr)
        {
            auto slipIt = obs.satStat_ptr->sigStatMap.find(ft2string(frequency));
            if (slipIt != obs.satStat_ptr->sigStatMap.end() &&
                (slipIt->second.phaseQuarantined || slipIt->second.tracking.halfCycle ||
                 slipIsExcluded(slipIt->second.slip)))
            {
                continue;
            }
        }

        return true;
    }

    return false;
}

bool slipRequiresArcRetirement(const SigStat::SlipStat& slip)
{
    const auto& reset = acsConfig.ambErrors.resetOnSlip;
    return (reset.LLI && slip.LLI) || (reset.retrack && slip.retrack) ||
        (reset.GF && slip.GF) || (reset.MW && slip.MW) ||
        (reset.SCDIA && slip.SCDIA) || (reset.single_freq && slip.singleFreq);
}

bool signalRequiresArcRetirement(const GObs& obs, E_ObsCode code)
{
    if (!obs.satStat_ptr)
    {
        return false;
    }

    for (const auto& [frequency, signal] : obs.sigs)
    {
        if (signal.code != code)
        {
            continue;
        }

        auto slipIt = obs.satStat_ptr->sigStatMap.find(ft2string(frequency));
        if (slipIt != obs.satStat_ptr->sigStatMap.end() &&
            (slipIt->second.tracking.pendingBreak || slipRequiresArcRetirement(slipIt->second.slip) ||
             slipRequiresArcRetirement(slipIt->second.savedSlip)))
        {
            return true;
        }
    }

    return false;
}

ReferenceAvailability referenceAvailability(
    ReceiverMap&                         receiverMap,
    E_Sys                                sys,
    const vector<E_ObsCode>&             baselineObservables
)
{
    ReferenceAvailability availability;

    for (auto& [id, receiver] : receiverMap)
    {
        if (!receiver.ready || receiver.obsList.empty())
        {
            continue;
        }

        auto& receiverOptions = acsConfig.getRecOpts(id);
        for (auto& obs : only<GObs>(receiver.obsList))
        {
            if (obs.Sat.sys != sys)
            {
                continue;
            }

            ZhangGraphEdge edge{id, obs.Sat};
            availability.rawEdges.insert(edge);
            bool hardDiscontinuity = false;
            for (E_ObsCode code : baselineObservables)
            {
                const auto phase = obs.satStat_ptr
                    ? obs.satStat_ptr->sigStatMap.find(enum_to_string(code))
                    : decltype(obs.satStat_ptr->sigStatMap)::iterator{};
                const bool rejectedPhase = obs.satStat_ptr &&
                    phase != obs.satStat_ptr->sigStatMap.end() &&
                    phase->second.phaseRejectCount >= acsConfig.ambErrors.phase_reject_limit;
                if(rejectedPhase && zhangEvidenceSlipQcEnabled())
                    availability.rejectedPhaseEdges.insert(edge);
                if(zhangEvidenceSlipQcEnabled() && obs.satStat_ptr)
                {
                    const auto frequency=code2Freq[obs.Sat.sys][code];
                    const auto& status=obs.satStat_ptr->sigStatMap[ft2string(frequency)];
                    if(status.phaseContinuityUnproven || status.slip.retrack ||
                        (status.slip.LLI && status.tracking.halfCycle))
                        availability.continuityUnknownEdges.insert(edge);
                }
                // Repeated filter rejection can still require a conservative
                // model-state retirement, but is not proof of physical slip.
                if (signalRequiresArcRetirement(obs, code) || rejectedPhase)
                {
                    hardDiscontinuity = true;
                    availability.discontinuitySignals[edge].insert(code);
                    const auto frequency = code2Freq[obs.Sat.sys][code];
                    availability.discontinuityStats[edge].push_back(
                        &obs.satStat_ptr->sigStatMap[ft2string(frequency)]);
                    availability.discontinuityStats[edge].push_back(
                        &obs.satStat_ptr->sigStatMap[enum_to_string(code)]);
                }
            }
            if (hardDiscontinuity)
            {
                availability.discontinuousEdges.insert(edge);
            }

            if (obs.exclude)
            {
                availability.qcExcludedEdges.insert(edge);
                continue;
            }

            // The graph must be built from edges that can actually create PPP rows.  The
            // elevation exclusion in receiverUducGnss is evaluated after the generic
            // observation flags, so checking only obs.exclude may select an unmodelled tree
            // edge and leave the satellite phase state rank deficient.
            if (acsConfig.exclude.elevation &&
                obs.satStat_ptr &&
                obs.satStat_ptr->el < receiverOptions.elevation_mask_deg * D2R)
            {
                availability.elevationExcludedEdges.insert(edge);
                continue;
            }

            // Code support has its own signal/QC requirements. Absence of a
            // carrier or a phase chart must not erase a valid code row.
            for (const auto& [frequency, signal] : obs.sigs)
                if (signal.P != 0 && !signal.invalid &&
                    zhangFullRankUsesObservable(signal.code, baselineObservables))
                    availability.codeEdges.insert(edge);

            bool usable = true;
            for (E_ObsCode code : baselineObservables)
            {
                const bool signalUsable = signalIsUsable(obs, code);
                if (signalUsable)
                    availability.signalEdges[static_cast<int>(code)].insert(edge);
                usable &= signalUsable;
            }

            if (!usable)
            {
                if(zhangEvidenceSlipQcEnabled() && obs.satStat_ptr)
                    for(const auto& [frequency,signal]:obs.sigs)
                        if(signal.L!=0 && obs.satStat_ptr->sigStatMap[ft2string(frequency)].phaseQuarantined)
                            availability.phaseQcEdges.insert(edge);
                availability.signalUnavailableEdges.insert(edge);
                continue;
            }

            availability.satellitesByReceiver[id].insert(obs.Sat);
            availability.edges.insert(edge);
            if (obs.satStat_ptr)
            {
                availability.elevationScore[obs.Sat] += obs.satStat_ptr->el;
                availability.edgeQuality[edge] = obs.satStat_ptr->el;
            }
        }
    }

    return availability;
}

set<SatSys> commonSatellites(const ReferenceAvailability& availability)
{
    set<SatSys> common;
    bool first = true;

    for (auto& [receiver, satellites] : availability.satellitesByReceiver)
    {
        if (first)
        {
            common = satellites;
            first = false;
            continue;
        }

        set<SatSys> intersection;
        std::set_intersection(
            common.begin(),
            common.end(),
            satellites.begin(),
            satellites.end(),
            std::inserter(intersection, intersection.begin())
        );
        common = std::move(intersection);
    }

    return common;
}

vector<string> orderedReceivers(
    const ReferenceAvailability& availability,
    const vector<string>&        candidates
)
{
    vector<string> ordered;
    set<string>    inserted;

    for (auto& candidate : candidates)
    {
        auto it = availability.satellitesByReceiver.find(candidate);
        if (it != availability.satellitesByReceiver.end() &&
            !it->second.empty() &&
            inserted.insert(candidate).second)
        {
            ordered.push_back(candidate);
        }
    }

    vector<std::pair<string, size_t>> remaining;
    for (auto& [receiver, satellites] : availability.satellitesByReceiver)
    {
        if (inserted.find(receiver) == inserted.end())
        {
            remaining.emplace_back(receiver, satellites.size());
        }
    }

    std::sort(
        remaining.begin(),
        remaining.end(),
        [](const auto& left, const auto& right)
        {
            if (left.second != right.second)
            {
                return left.second > right.second;
            }
            return left.first < right.first;
        }
    );
    for (auto& [receiver, count] : remaining)
    {
        ordered.push_back(receiver);
    }

    return ordered;
}

vector<SatSys> orderedSatellites(
    const set<SatSys>&            common,
    const map<SatSys, double>&    elevationScore,
    const vector<string>&         candidates
)
{
    vector<SatSys> ordered;
    set<SatSys>    inserted;

    for (auto& candidateId : candidates)
    {
        SatSys candidate(candidateId.c_str());
        if (common.find(candidate) != common.end() && inserted.insert(candidate).second)
        {
            ordered.push_back(candidate);
        }
    }

    vector<std::pair<SatSys, double>> remaining;
    for (auto& satellite : common)
    {
        if (inserted.find(satellite) == inserted.end())
        {
            auto scoreIt = elevationScore.find(satellite);
            remaining.emplace_back(
                satellite,
                scoreIt == elevationScore.end() ? 0 : scoreIt->second
            );
        }
    }

    std::sort(
        remaining.begin(),
        remaining.end(),
        [](const auto& left, const auto& right)
        {
            if (left.second != right.second)
            {
                return left.second > right.second;
            }
            return left.first < right.first;
        }
    );
    for (auto& [satellite, score] : remaining)
    {
        ordered.push_back(satellite);
    }

    return ordered;
}

void addCoefficient(
    map<KFKey, double>& coefficients,
    const KFKey&        source,
    double              value
)
{
    if (value == 0)
    {
        return;
    }

    coefficients[source] += value;
    if (coefficients[source] == 0)
    {
        coefficients.erase(source);
    }
}

bool addRequiredCoefficient(
    const KFState&      kfState,
    map<KFKey, double>& coefficients,
    const KFKey&        source,
    double              value
)
{
    if (value == 0)
    {
        return true;
    }

    if (kfState.kfIndexMap.find(source) == kfState.kfIndexMap.end())
    {
        return false;
    }

    addCoefficient(coefficients, source, value);
    return true;
}

bool transformClockDatum(
    const KFState&                         kfState,
    map<KFKey, map<KFKey, double>>&        transform,
    E_Sys                                  system,
    const string&                          oldReceiver,
    const string&                          newReceiver,
    KF                                     receiverType,
    KF                                     satelliteType
)
{
    set<int> components;
    set<string> receivers = {oldReceiver};
    vector<KFKey> satellites;

    for (auto& [key, index] : kfState.kfIndexMap)
    {
        if (key.type == receiverType &&
            zhangReceiverClockMatchesSystem(acsConfig, key, system))
        {
            components.insert(key.num);
            receivers.insert(key.str);
        }
        else if (key.type == satelliteType && key.Sat.sys == system)
        {
            components.insert(key.num);
            satellites.push_back(key);
        }
    }

    if (components.empty())
    {
        return true;
    }

    receivers.insert(newReceiver);

    for (int component : components)
    {
        KFKey newReferenceKey = zhangReceiverClockKey(
            acsConfig, newReceiver, system, component, receiverType);

        if (newReceiver != oldReceiver &&
            kfState.kfIndexMap.find(newReferenceKey) == kfState.kfIndexMap.end())
        {
            return false;
        }

        for (auto& receiver : receivers)
        {
            if (receiver == newReceiver)
            {
                continue;
            }

            KFKey destination = zhangReceiverClockKey(
                acsConfig, receiver, system, component, receiverType);

            auto& coefficients = transform[destination];

            if (receiver != oldReceiver)
            {
                KFKey oldReceiverKey = destination;
                if (!addRequiredCoefficient(kfState, coefficients, oldReceiverKey, +1))
                {
                    return false;
                }
            }

            if (newReceiver != oldReceiver)
            {
                addCoefficient(coefficients, newReferenceKey, -1);
            }
        }

        for (auto& satellite : satellites)
        {
            if (satellite.num != component)
            {
                continue;
            }

            auto& coefficients = transform[satellite];
            addCoefficient(coefficients, satellite, +1);
            if (newReceiver != oldReceiver)
            {
                addCoefficient(coefficients, newReferenceKey, -1);
            }
        }
    }

    return true;
}

double zhangWavelength(E_Sys sys, E_ObsCode code)
{
    auto sysFrequencyIt = code2Freq.find(sys);
    if (sysFrequencyIt == code2Freq.end())
    {
        return 0;
    }

    auto frequencyIt = sysFrequencyIt->second.find(code);
    if (frequencyIt == sysFrequencyIt->second.end())
    {
        return 0;
    }

    auto wavelengthIt = genericWavelength.find(frequencyIt->second);
    if (wavelengthIt == genericWavelength.end())
    {
        return 0;
    }

    return wavelengthIt->second;
}

/** Lightweight staged product moments from the authoritative FLOAT x/P.
 * This reads only scalar entries and does not clone the (potentially large)
 * network covariance. Corrections follow the writer's clock-minus-phase sign.
 * Full graph-event cross covariance is preserved by the state transform and
 * checked separately; these rows are diagnostics, not integer certificates. */
void traceZhangGraphProductProjection(
    Trace& trace,
    const KFState& state,
    E_Sys system,
    const vector<E_ObsCode>& observables,
    const char* stage)
{
    if (!zhangR51SingleBasisHandoffEnabled()
        || !acsConfig.zhangPppAr.output_diagnostics) return;
    set<SatSys> satellites;
    for (const auto& [key, index] : state.kfIndexMap)
        if (key.type == KF::PHASE_BIAS
            && key.Sat.sys == system && key.Sat.prn > 0
            && key.str.empty()
            && zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num), observables))
            satellites.insert(key.Sat);
    auto stateIndex = [&](const KFKey& key) -> int
    {
        const auto found = state.kfIndexMap.find(key);
        if (found == state.kfIndexMap.end()
            || found->second < 0 || found->second >= state.x.size()
            || found->second >= state.P.rows()
            || found->second >= state.P.cols()) return -1;
        return found->second;
    };
    for (const auto& satellite : satellites)
    {
        KFKey clockKey;
        clockKey.type = KF::SAT_CLOCK;
        clockKey.Sat = satellite;
        const int clock = stateIndex(clockKey);
        map<E_ObsCode, int> phaseIndices;
        for (const auto code : observables)
        {
            const int phase = stateIndex(
                zhangTransportSatelliteKey(code, satellite));
            if (phase < 0) continue;
            phaseIndices[code] = phase;
            if (clock < 0)
            {
                trace << "\nZHANG_PRODUCT_PROJECTION_STAGE time="
                      << state.time.to_string(0)
                      << " stage=" << stage
                      << " system=" << enum_to_string(system)
                      << " satellite=" << satellite.id()
                      << " signal=" << enum_to_string(code)
                      << " valid=0 reason=SAT_CLOCK_STATE_MISSING";
                continue;
            }
            const double correctionVariance = state.P(clock, clock)
                + state.P(phase, phase)
                - 2 * state.P(clock, phase);
            const bool finiteProjection =
                std::isfinite(state.x(phase))
                && std::isfinite(state.x(clock))
                && std::isfinite(state.P(phase, phase))
                && std::isfinite(state.P(clock, clock))
                && std::isfinite(state.P(clock, phase))
                && std::isfinite(correctionVariance)
                && correctionVariance >= -1e-10;
            trace << "\nZHANG_PRODUCT_PROJECTION_STAGE time="
                  << state.time.to_string(0)
                  << " stage=" << stage
                  << " system=" << enum_to_string(system)
                  << " satellite=" << satellite.id()
                  << " signal=" << enum_to_string(code)
                  << " valid=" << finiteProjection
                  << " phase_mean_m=" << state.x(phase)
                  << " clock_mean_m=" << state.x(clock)
                  << " correction_mean_m="
                  << state.x(clock) - state.x(phase)
                  << " phase_var_m2=" << state.P(phase, phase)
                  << " clock_var_m2=" << state.P(clock, clock)
                  << " phase_clock_cov_m2=" << state.P(phase, clock)
                  << " correction_var_m2=" << correctionVariance
                  << " state_dim=" << state.x.size()
                  << " factor_commit_sequence="
                  << state.factorCommitSequence
                  << " integer_certification=0 publication_authority=0";
        }
        if (clock >= 0 && observables.size() >= 2)
        {
            const auto first = phaseIndices.find(observables[0]);
            const auto second = phaseIndices.find(observables[1]);
            if (first != phaseIndices.end() && second != phaseIndices.end())
            {
                const int a = first->second;
                const int b = second->second;
                const double correctionCovariance =
                    state.P(clock, clock) - state.P(clock, a)
                    - state.P(clock, b) + state.P(a, b);
                trace << "\nZHANG_PRODUCT_PROJECTION_PAIR_STAGE time="
                      << state.time.to_string(0)
                      << " stage=" << stage
                      << " system=" << enum_to_string(system)
                      << " satellite=" << satellite.id()
                      << " first=" << enum_to_string(observables[0])
                      << " second=" << enum_to_string(observables[1])
                      << " phase_cross_cov_m2=" << state.P(a, b)
                      << " correction_cross_cov_m2="
                      << correctionCovariance;
            }
        }
    }
}

/** One product target chart is a family of rows on the sole KF S posterior.
 * The reference satellite carries a zero G row; every other satellite has
 * phase_s + wavelength * G_s * k.  An unavailable chord column is not zero. */
bool zhangGraphProductProjectionRows(
    const KFState& state,
    E_Sys system,
    const vector<E_ObsCode>& observables,
    const ZhangGraphBasis& stateBasis,
    const ZhangGraphBasis& targetChart,
    bool requireCurrentKeys,
    map<pair<SatSys, E_ObsCode>, ZhangCoordinateRow>& rows,
    string& reason)
{
    rows.clear();
    const auto target = zhangBuildSatelliteProductTarget(
        stateBasis, targetChart);
    if (!target.valid || target.productTreeTransported)
    {
        reason = target.valid
            ? "TARGET_CHART_STALE_TREE_TRANSPORT_NOT_AUTHORISED"
            : target.failureReason;
        return false;
    }
    map<SatSys, ZhangExactVector> directions;
    directions[target.referenceSatellite] =
        ZhangExactVector(target.currentChords.size());
    for (size_t row = 0; row < target.targetSatellites.size(); row++)
        directions[target.targetSatellites[row]] = target.matrix[row];
    bool complete = true;
    for (const auto code : observables)
    {
        const double lambda = zhangWavelength(system, code);
        if (!(lambda > 0) || !std::isfinite(lambda))
        {
            reason = "PRODUCT_WAVELENGTH_UNAVAILABLE";
            rows.clear();
            return false;
        }
        for (const auto& satellite : stateBasis.satellites)
        {
            const auto direction = directions.find(satellite);
            if (direction == directions.end() ||
                direction->second.size() != target.currentChords.size())
            {
                reason = "TARGET_SATELLITE_DIRECTION_MISSING";
                rows.clear();
                return false;
            }
            ZhangCoordinateRow row;
            const KFKey phase = zhangTransportSatelliteKey(code, satellite);
            if (requireCurrentKeys && !state.kfIndexMap.contains(phase))
            {
                complete = false;
                continue;
            }
            row[phase] = 1;
            bool rowAvailable = true;
            for (size_t chord = 0; chord < target.currentChords.size(); chord++)
            {
                const auto& exact = direction->second[chord];
                if (exact == 0) continue;
                const KFKey ambiguity = zhangTransportAmbiguityKey(
                    code, target.currentChords[chord]);
                if (requireCurrentKeys &&
                    !state.kfIndexMap.contains(ambiguity))
                {
                    rowAvailable = false;
                    break;
                }
                const double coefficient = lambda * exact.convert_to<double>();
                if (!std::isfinite(coefficient))
                {
                    reason = "TARGET_EXACT_COEFFICIENT_OUT_OF_RANGE";
                    rows.clear();
                    return false;
                }
                row[ambiguity] += coefficient;
            }
            if (rowAvailable)
                rows[{satellite, code}] = std::move(row);
            else
                complete = false;
        }
    }
    reason = complete
        ? "EXACT_TARGET_CHART_ROWS"
        : "PARTIAL_TARGET_CHART_ROWS_MISSING_STATE";
    return !rows.empty();
}

/** Preserve the FLOAT event's old satellite phase moments before any tree
 * edge or ambiguity coordinate is retired.  An old tree path is an exact
 * signed sum of physical edge rows.  Re-expressing that path on the new chart
 * is only possible when every versioned arc survives; this is deliberately
 * stricter than a numeric fit to post-event observations. */
using ZhangSignedPhysicalPath = map<ZhangGraphEdge, int>;

map<SatSys, ZhangSignedPhysicalPath> zhangSatelliteTreePaths(
    const ZhangGraphBasis& basis)
{
    map<string, ZhangSignedPhysicalPath> receiverPaths;
    map<SatSys, ZhangSignedPhysicalPath> satellitePaths;
    receiverPaths[basis.rootReceiver] = {};
    bool progress = true;
    while (progress)
    {
        progress = false;
        for (const auto& edge : basis.treeEdges)
        {
            const auto receiver = receiverPaths.find(edge.receiver);
            const auto satellite = satellitePaths.find(edge.satellite);
            if ((receiver != receiverPaths.end()) ==
                (satellite != satellitePaths.end())) continue;
            ZhangSignedPhysicalPath expression = receiver != receiverPaths.end()
                ? receiver->second : satellite->second;
            for (auto& [pathEdge, coefficient] : expression)
                coefficient = -coefficient;
            expression[edge]++;
            if (receiver != receiverPaths.end())
                satellitePaths.emplace(edge.satellite, std::move(expression));
            else
                receiverPaths.emplace(edge.receiver, std::move(expression));
            progress = true;
        }
    }
    return satellitePaths;
}

vector<ZhangGraphProductEventReceipt> prepareGraphProductEventReceipts(
    const KFState& oldState,
    E_Sys system,
    const vector<E_ObsCode>& observables,
    const ZhangGraphBasis& oldBasis,
    const ZhangGraphBasis& oldProductTarget,
    const ZhangGraphBasis& newBasis,
    const map<ZhangGraphEdge, map<int, int>>& oldVersions,
    const map<ZhangGraphEdge, map<int, int>>& currentVersions,
    const set<ZhangGraphEdge>& retiredArcs,
    const ZhangGraphCoordinatePlan& plan,
    const string& eventKind)
{
    const auto oldSatellitePaths = zhangSatelliteTreePaths(oldBasis);
    using ProductKey = pair<SatSys, E_ObsCode>;
    map<ProductKey, ZhangCoordinateRow> oldProductRows;
    string oldProjectionReason;
    const bool oldProjectionAvailable = zhangGraphProductProjectionRows(
        oldState, system, observables, oldBasis, oldProductTarget, true,
        oldProductRows, oldProjectionReason);
    const SatSys productReference = oldProductTarget.satellites.empty()
        ? SatSys() : *oldProductTarget.satellites.begin();
    auto signalVersion = [](const auto& identities,
                            const ZhangGraphEdge& edge, E_ObsCode code)
        -> std::optional<int>
    {
        const auto found = identities.find(edge);
        if (found == identities.end()) return std::nullopt;
        const auto signal = found->second.find(int(code));
        if (signal == found->second.end()) return std::nullopt;
        return signal->second;
    };
    auto projectedMean = [&](const ZhangCoordinateRow& row)
        -> std::optional<double>
    {
        double value = 0;
        for (const auto& [key, coefficient] : row)
        {
            const auto state = oldState.kfIndexMap.find(key);
            if (state == oldState.kfIndexMap.end() ||
                state->second < 0 || state->second >= oldState.x.size())
                return std::nullopt;
            value += coefficient * oldState.x(state->second);
        }
        if (!std::isfinite(value)) return std::nullopt;
        return value;
    };
    auto projectedCovariance = [&](const ZhangCoordinateRow& left,
                                   const ZhangCoordinateRow& right)
        -> std::optional<double>
    {
        double value = 0;
        for (const auto& [first, a] : left)
        {
            const auto i = oldState.kfIndexMap.find(first);
            if (i == oldState.kfIndexMap.end() || i->second < 0 ||
                i->second >= oldState.P.rows()) return std::nullopt;
            for (const auto& [second, b] : right)
            {
                const auto j = oldState.kfIndexMap.find(second);
                if (j == oldState.kfIndexMap.end() || j->second < 0 ||
                    j->second >= oldState.P.cols()) return std::nullopt;
                value += a * b * oldState.P(i->second, j->second);
            }
        }
        if (!std::isfinite(value)) return std::nullopt;
        return value;
    };

    vector<ZhangGraphProductEventReceipt> receipts;
    for (const auto observable : observables)
    for (const auto& satellite : oldBasis.satellites)
    {
        const auto oldProductFunctionals = oldProductTarget.connected
            ? zhangBuildProductIntegerFunctionals(
                oldProductTarget, oldVersions, observable,
                productReference)
            : map<SatSys, ZhangProductIntegerFunctional>{};
        ZhangGraphProductEventReceipt receipt;
        receipt.time = oldState.time;
        receipt.system = system;
        receipt.observable = observable;
        receipt.satellite = satellite;
        receipt.eventKind = eventKind;
        const auto owner = graphStateMap.find(
            {zhangGraphRuntimeId(oldState), system});
        if (owner != graphStateMap.end())
        {
            receipt.oldChartGeneration = owner->second.representationVersion;
            receipt.oldBackendGeneration = owner->second.eventCounter;
            receipt.newBackendGeneration = owner->second.eventCounter;
        }
        receipt.fullJointTransportId = zhangGraphRuntimeId(oldState) + ":" +
            oldState.time.to_string(9) + ":" + enum_to_string(system) +
            ":" + eventKind;
        const auto productKey = ProductKey{satellite, observable};
        const auto oldFunctional = oldProductFunctionals.find(satellite);
        receipt.productReferenceSatellite = productReference;
        if (oldFunctional != oldProductFunctionals.end() &&
            oldFunctional->second.valid)
        {
            const auto& functional = oldFunctional->second;
            receipt.oldProductFunctionalEdges = functional.physicalEdges;
            receipt.oldProductFunctionalArcVersions =
                functional.physicalArcVersions;
            receipt.oldProductFunctionalCoefficients =
                functional.networkCoefficients;
            receipt.oldProductFunctionalAffineOffsetCycles =
                functional.affineOffsetCycles;
            for (size_t index = 0;
                 index < functional.physicalEdges.size(); index++)
                if (functional.networkCoefficients[index] != 0 &&
                    !signalVersion(oldVersions,
                        functional.physicalEdges[index], observable))
                    receipt.failureReason =
                        "OLD_PRODUCT_FUNCTIONAL_ARC_VERSION_MISSING";
        }
        const auto oldRow = oldProductRows.find(productKey);
        receipt.oldProductProjectionCaptured = oldProjectionAvailable &&
            oldRow != oldProductRows.end();
        if (receipt.oldProductProjectionCaptured)
            receipt.oldPhaseRow = oldRow->second;
        const auto referencePath = oldSatellitePaths.find(productReference);
        if (receipt.failureReason != "NOT_EVALUATED")
        {
            // Preserve the more specific missing-version failure above.
        }
        else if (!receipt.oldProductProjectionCaptured)
            receipt.failureReason = "OLD_PRODUCT_TARGET_PROJECTION_NOT_CAPTURED:" +
                oldProjectionReason;
        else if (oldFunctional == oldProductFunctionals.end() ||
            !oldFunctional->second.valid)
            receipt.failureReason = "OLD_VERSIONED_PRODUCT_FUNCTIONAL_MISSING";
        else if (referencePath == oldSatellitePaths.end())
            receipt.failureReason = "OLD_S_REFERENCE_PHYSICAL_PATH_MISSING";
        else if (!newBasis.satellites.contains(satellite))
            receipt.failureReason = "SATELLITE_ABSENT_FROM_NEW_CHART";
        else if (oldBasis.rootReceiver != newBasis.rootReceiver)
            receipt.failureReason = "PRODUCT_ROOT_CHANGED_WITHOUT_GAUGE_MAP";
        else if (!(zhangWavelength(system, observable) > 0))
            receipt.failureReason = "PRODUCT_WAVELENGTH_UNAVAILABLE";

        if (receipt.oldProductProjectionCaptured)
        {
            const auto mean = projectedMean(receipt.oldPhaseRow);
            const auto variance = projectedCovariance(
                receipt.oldPhaseRow, receipt.oldPhaseRow);
            if (!mean || !variance)
            {
                receipt.failureReason = "OLD_PRODUCT_MOMENTS_NOT_FINITE";
            }
            else
            {
                receipt.oldPhaseMeanMetres = *mean;
                receipt.oldPhaseVarianceMetres2 = *variance;
            }
            for (const auto& [otherKey, otherRow] : oldProductRows)
            {
                const auto covariance = projectedCovariance(
                    receipt.oldPhaseRow, otherRow);
                if (covariance)
                    receipt.oldPhaseCrossCovarianceMetres2[otherKey] =
                        *covariance;
            }
            for (const auto& [key, index] : oldState.kfIndexMap)
            {
                if (key.type != KF::SAT_CLOCK || key.Sat.sys != system ||
                    key.Sat.prn <= 0 || !key.str.empty() || key.num != 0)
                    continue;
                const ZhangCoordinateRow clockRow{{key, 1}};
                const auto covariance = projectedCovariance(
                    receipt.oldPhaseRow, clockRow);
                if (covariance)
                    receipt.oldPhaseClockCrossCovarianceMetres2[key.Sat] =
                        *covariance;
            }
        }
        KFKey clockKey;
        clockKey.type = KF::SAT_CLOCK;
        clockKey.Sat = satellite;
        const auto oldClock = oldState.kfIndexMap.find(clockKey);
        if (oldClock == oldState.kfIndexMap.end())
        {
            if (receipt.failureReason == "NOT_EVALUATED")
                receipt.failureReason = "OLD_SATELLITE_CLOCK_STATE_MISSING";
        }
        else
        {
            const int index = oldClock->second;
            receipt.oldClockMeanMetres = oldState.x(index);
            receipt.oldClockVarianceMetres2 = oldState.P(index, index);
            if (receipt.oldProductProjectionCaptured)
            {
                const ZhangCoordinateRow clockRow{{clockKey, 1}};
                const auto covariance = projectedCovariance(
                    clockRow, receipt.oldPhaseRow);
                if (covariance)
                    receipt.oldClockPhaseCovarianceMetres2 = *covariance;
            }
            for (const auto& [key, otherIndex] : oldState.kfIndexMap)
                if (key.type == KF::SAT_CLOCK && key.Sat.sys == system
                    && key.Sat.prn > 0 && key.str.empty() && key.num == 0)
                    receipt.oldClockClockCrossCovarianceMetres2[key.Sat] =
                        oldState.P(index, otherIndex);
        }
        if (receipt.failureReason != "NOT_EVALUATED")
        {
            receipts.push_back(std::move(receipt));
            continue;
        }

        // The old absolute product phase is the old S reference-node
        // potential plus the versioned old P satellite-minus-reference path.
        // It is not the old S path of this satellite.
        map<ZhangGraphEdge, ZhangExactInteger> oldPhysicalPath;
        for (const auto& [edge, coefficient] : referencePath->second)
            oldPhysicalPath[edge] += coefficient;
        for (size_t term = 0;
             term < receipt.oldProductFunctionalEdges.size(); term++)
            oldPhysicalPath[receipt.oldProductFunctionalEdges[term]] +=
                receipt.oldProductFunctionalCoefficients[term];
        ZhangCoordinateRow reconstructedOldRow;
        const double lambda = zhangWavelength(system, observable);
        for (const auto& [edge, exactCoefficient] : oldPhysicalPath)
        {
            if (exactCoefficient == 0) continue;
            const auto previous = signalVersion(oldVersions, edge, observable);
            const auto current = signalVersion(currentVersions, edge, observable);
            if (!previous ||
                !oldBasis.edges.contains(edge))
            {
                receipt.failureReason =
                    "OLD_PRODUCT_PHYSICAL_ARC_OR_VERSION_MISSING";
                break;
            }
            const double coefficient = exactCoefficient.convert_to<double>();
            if (!std::isfinite(coefficient))
            {
                receipt.failureReason =
                    "OLD_PRODUCT_PHYSICAL_COEFFICIENT_OUT_OF_RANGE";
                break;
            }
            receipt.physicalPathEdges.push_back(edge);
            receipt.physicalPathCoefficients.push_back(exactCoefficient);
            receipt.oldPhysicalArcVersions.push_back(*previous);
            receipt.newPhysicalArcVersions.push_back(
                current.value_or(-1));
            ZhangCoordinateRow oldArc;
            if (edge.receiver != oldBasis.rootReceiver)
                oldArc[zhangTransportReceiverKey(
                    system, observable, edge.receiver)] = 1;
            oldArc[zhangTransportSatelliteKey(
                observable, edge.satellite)] = 1;
            if (!oldBasis.treeEdges.contains(edge))
                oldArc[zhangTransportAmbiguityKey(
                    observable, edge)] = lambda;
            zhangAddCoordinateRow(reconstructedOldRow,
                oldArc, coefficient);

            const bool physicalIdentityChanged =
                !current || *current != *previous;
            const bool needsAuxiliary = physicalIdentityChanged ||
                !newBasis.edges.contains(edge);
            receipt.physicalArcChanged |= physicalIdentityChanged;
            if (needsAuxiliary)
            {
                const auto auxiliary =
                    plan.retainedOldArcAuxiliaryKeys.find(
                        {observable, edge});
                if (auxiliary == plan.retainedOldArcAuxiliaryKeys.end() ||
                    !plan.transform.contains(auxiliary->second))
                {
                    receipt.failureReason =
                        "OLD_PRODUCT_PATH_ARC_RETIRED_WITHOUT_CORRELATED_AUXILIARY";
                    break;
                }
                receipt.retainedOldArcAuxiliaryKeys[edge] =
                    auxiliary->second;
                receipt.newPhaseRow[auxiliary->second] += coefficient;
                continue;
            }
            ZhangCoordinateRow newArc;
            if (edge.receiver != newBasis.rootReceiver)
                newArc[zhangTransportReceiverKey(
                    system, observable, edge.receiver)] = 1;
            newArc[zhangTransportSatelliteKey(
                observable, edge.satellite)] = 1;
            if (!newBasis.treeEdges.contains(edge))
                newArc[zhangTransportAmbiguityKey(
                    observable, edge)] = lambda;
            zhangAddCoordinateRow(receipt.newPhaseRow,
                newArc, coefficient);
        }
        if (receipt.failureReason == "NOT_EVALUATED")
        {
            ZhangCoordinateRow oldResidual = reconstructedOldRow;
            zhangAddCoordinateRow(oldResidual, receipt.oldPhaseRow, -1);
            for (const auto& [key, coefficient] : oldResidual)
                if (!std::isfinite(coefficient) ||
                    std::abs(coefficient) > 1e-10)
                    receipt.failureReason =
                        "OLD_PRODUCT_PATH_AND_G_PROJECTION_MISMATCH";
        }
        if (receipt.failureReason == "NOT_EVALUATED")
        {
            const auto abstractProof = zhangVerifyGraphFunctionalTransport(
                plan, receipt.oldPhaseRow, receipt.newPhaseRow,
                receipt.oldOffsetMetres, receipt.newOffsetMetres);
            receipt.representableSameTimeMapProven =
                abstractProof.sameTimeMapProven;
            receipt.conditionalMapReason =
                "ACTUAL_PRODUCT_WRITER_PROJECTION_NOT_BOUND";
            receipt.failureReason = !abstractProof.sameTimeMapProven
                ? abstractProof.failureReason
                : "ACTUAL_PRODUCT_WRITER_PROJECTION_NOT_BOUND";
        }
        receipts.push_back(std::move(receipt));
    }
    return receipts;
}

/** Bind a graph-event receipt only after the product target chart has been
 * committed.  The event transform was applied earlier, when its eventual
 * writer row was unknown.  Equality here is row equality on the same KF
 * posterior; an exact integer decision for changed physical functionals is
 * still a separate prerequisite in the writer. */
void bindGraphProductEventWriterRows(
    Trace& trace,
    const KFState& state,
    E_Sys system,
    const vector<E_ObsCode>& observables,
    const ZhangGraphBasis& stateBasis,
    const GraphRuntimeState& runtime)
{
    map<pair<SatSys, E_ObsCode>, ZhangCoordinateRow> writerRows;
    string projectionReason;
    const bool projectionsAvailable = zhangGraphProductProjectionRows(
        state, system, observables, stateBasis, runtime.productBasis, true,
        writerRows, projectionReason);
    const SatSys reference = runtime.productBasis.satellites.empty()
        ? SatSys() : *runtime.productBasis.satellites.begin();
    const auto newSatellitePaths = zhangSatelliteTreePaths(stateBasis);
    const string runtimeId = zhangGraphRuntimeId(state);
    const string epoch = state.time.to_string(9);
    int bound = 0;
    int proven = 0;
    for (auto& [identity, receipt] : graphProductEventReceipts)
    {
        if (std::get<0>(identity) != runtimeId ||
            std::get<1>(identity) != epoch ||
            std::get<2>(identity) != system) continue;
        const auto row = writerRows.find({receipt.satellite,
            receipt.observable});
        if (!projectionsAvailable || row == writerRows.end())
        {
            receipt.actualWriterProjectionBound = false;
            receipt.actualWriterRowMatches = false;
            receipt.sameTimeMapProven = false;
            receipt.conditionalProductMapProven = false;
            receipt.failureReason = "ACTUAL_PRODUCT_WRITER_ROW_UNAVAILABLE:" +
                projectionReason;
            continue;
        }
        receipt.actualWriterPhaseRow = row->second;
        receipt.actualWriterProjectionBound = true;
        receipt.newChartGeneration = runtime.representationVersion;
        receipt.newBackendGeneration = runtime.eventCounter;
        if (!std::isfinite(receipt.oldOffsetMetres) ||
            !std::isfinite(receipt.newOffsetMetres) ||
            receipt.oldOffsetMetres != 0 ||
            receipt.newOffsetMetres != 0)
        {
            receipt.sameTimeMapProven = false;
            receipt.conditionalProductMapProven = false;
            receipt.failureReason = "NONZERO_OR_NONFINITE_PRODUCT_OFFSET_NOT_MODELED";
            continue;
        }
        bound++;
        ZhangCoordinateRow difference = receipt.actualWriterPhaseRow;
        zhangAddCoordinateRow(difference, receipt.newPhaseRow, -1);
        receipt.bridgePhysicalEdges.clear();
        receipt.bridgePhysicalArcVersions.clear();
        receipt.bridgePhysicalCoefficients.clear();
        receipt.bridgeStateRowCycles.clear();
        bool rowEqual = !receipt.newPhaseRow.empty();
        for (const auto& [key, coefficient] : difference)
            if (!std::isfinite(coefficient) ||
                std::abs(coefficient) > 1e-10) rowEqual = false;
        receipt.actualWriterRowMatches = rowEqual;
        receipt.sameTimeMapProven =
            receipt.oldProductProjectionCaptured &&
            receipt.representableSameTimeMapProven &&
            receipt.fullJointMomentTransportProven && rowEqual;
        receipt.conditionalProductMapProven = false;
        receipt.conditionalMapReason =
            "CURRENT_PRODUCT_FUNCTIONAL_BRIDGE_NOT_EVALUATED";
        if (!receipt.sameTimeMapProven &&
            receipt.oldProductProjectionCaptured &&
            receipt.representableSameTimeMapProven &&
            receipt.fullJointMomentTransportProven)
        {
            const auto currentFunctionals = runtime.productBasis.connected
                ? zhangBuildProductIntegerFunctionals(
                    runtime.productBasis, runtime.phaseArcVersions,
                    receipt.observable, reference)
                : map<SatSys, ZhangProductIntegerFunctional>{};
            ZhangProductIntegerFunctional previous;
            previous.satellite = receipt.satellite;
            previous.referenceSatellite = receipt.productReferenceSatellite;
            previous.observable = receipt.observable;
            previous.physicalEdges = receipt.physicalPathEdges;
            previous.physicalArcVersions = receipt.oldPhysicalArcVersions;
            previous.networkCoefficients =
                receipt.physicalPathCoefficients;
            previous.affineOffsetCycles =
                receipt.oldProductFunctionalAffineOffsetCycles;
            previous.valid =
                previous.physicalEdges.size() ==
                    previous.physicalArcVersions.size() &&
                previous.physicalEdges.size() ==
                    previous.networkCoefficients.size();
            const auto current = currentFunctionals.find(receipt.satellite);
            const auto referencePath = newSatellitePaths.find(reference);
            if (!previous.valid || current == currentFunctionals.end() ||
                referencePath == newSatellitePaths.end())
                receipt.conditionalMapReason =
                    "PRODUCT_FUNCTIONAL_ENDPOINT_MISSING";
            else
            {
                ZhangProductIntegerFunctional currentAbsolute;
                currentAbsolute.satellite = receipt.satellite;
                currentAbsolute.referenceSatellite = reference;
                currentAbsolute.observable = receipt.observable;
                currentAbsolute.affineOffsetCycles =
                    current->second.affineOffsetCycles;
                map<ZhangGraphEdge, ZhangExactInteger> absolutePath;
                for (const auto& [edge, coefficient] :
                     referencePath->second)
                    absolutePath[edge] += coefficient;
                for (size_t term = 0;
                     term < current->second.physicalEdges.size(); term++)
                    absolutePath[current->second.physicalEdges[term]] +=
                        current->second.networkCoefficients[term];
                currentAbsolute.valid = true;
                for (const auto& [edge, coefficient] : absolutePath)
                {
                    if (coefficient == 0) continue;
                    const auto identity = runtime.phaseArcVersions.find(edge);
                    const auto signal = identity ==
                            runtime.phaseArcVersions.end()
                        ? map<int, int>::const_iterator{}
                        : identity->second.find(int(receipt.observable));
                    if (identity == runtime.phaseArcVersions.end() ||
                        signal == identity->second.end() ||
                        signal->second < 0)
                    {
                        currentAbsolute.valid = false;
                        break;
                    }
                    currentAbsolute.physicalEdges.push_back(edge);
                    currentAbsolute.physicalArcVersions.push_back(
                        signal->second);
                    currentAbsolute.networkCoefficients.push_back(
                        coefficient);
                }
                if (!currentAbsolute.valid)
                {
                    receipt.conditionalMapReason =
                        "CURRENT_ABSOLUTE_SIGNAL_ARC_VERSION_MISSING";
                    continue;
                }
                const auto transition =
                    zhangProductIntegerFunctionalDifference(
                        previous, currentAbsolute);
                const double lambda = zhangWavelength(
                    system, receipt.observable);
                auto physicalRow = [&](const ZhangGraphEdge& edge,
                                       int version,
                                       ZhangCoordinateRow& row,
                                       string& reason)
                {
                    row.clear();
                    const auto identity = runtime.phaseArcVersions.find(edge);
                    const auto signal = identity ==
                            runtime.phaseArcVersions.end()
                        ? map<int, int>::const_iterator{}
                        : identity->second.find(int(receipt.observable));
                    if (identity != runtime.phaseArcVersions.end() &&
                        signal != identity->second.end() &&
                        signal->second == version &&
                        stateBasis.edges.contains(edge))
                    {
                        if (edge.receiver != stateBasis.rootReceiver)
                            row[zhangTransportReceiverKey(
                                system, receipt.observable,
                                edge.receiver)] = 1;
                        row[zhangTransportSatelliteKey(
                            receipt.observable, edge.satellite)] = 1;
                        if (!stateBasis.treeEdges.contains(edge))
                            row[zhangTransportAmbiguityKey(
                                receipt.observable, edge)] = lambda;
                    }
                    else
                    {
                        const auto auxiliary =
                            receipt.retainedOldArcAuxiliaryKeys.find(edge);
                        bool exactOldVersion = false;
                        for (size_t old = 0;
                             old < receipt.physicalPathEdges.size(); old++)
                            if (receipt.physicalPathEdges[old] == edge &&
                                old < receipt.oldPhysicalArcVersions.size() &&
                                receipt.oldPhysicalArcVersions[old] == version)
                                exactOldVersion = true;
                        if (auxiliary ==
                                receipt.retainedOldArcAuxiliaryKeys.end() ||
                            !exactOldVersion)
                        {
                            reason =
                                "RETIRED_PRODUCT_ARC_AUXILIARY_MISSING";
                            return false;
                        }
                        row[auxiliary->second] = 1;
                    }
                    for (const auto& [key, coefficient] : row)
                        if (!state.kfIndexMap.contains(key) ||
                            !std::isfinite(coefficient))
                        {
                            reason =
                                "PRODUCT_BRIDGE_STATE_ROW_MISSING";
                            return false;
                        }
                    reason = "VERSIONED_PHYSICAL_ROW_IN_POSTERIOR";
                    return true;
                };
                const auto proof = zhangProveProductConditionalMap(
                    receipt.newPhaseRow, receipt.actualWriterPhaseRow,
                    transition, lambda,
                    receipt.fullJointMomentTransportProven, physicalRow);
                receipt.conditionalMapReason = proof.reason;
                receipt.conditionalProductMapProven = proof.valid;
                if (proof.valid)
                {
                    receipt.bridgeStateRowCycles =
                        proof.bridgeStateRowCycles;
                    receipt.bridgePhysicalEdges = proof.physicalEdges;
                    receipt.bridgePhysicalArcVersions =
                        proof.physicalArcVersions;
                    receipt.bridgePhysicalCoefficients =
                        proof.physicalCoefficients;
                }
                else
                {
                    receipt.bridgePhysicalEdges.clear();
                    receipt.bridgePhysicalArcVersions.clear();
                    receipt.bridgePhysicalCoefficients.clear();
                    receipt.bridgeStateRowCycles.clear();
                }
            }
        }
        receipt.failureReason = receipt.sameTimeMapProven
            ? "NONE" : receipt.conditionalProductMapProven
            ? "CONDITIONAL_CURRENT_INTEGER_BRIDGE_REQUIRED"
            : "ACTUAL_WRITER_PRODUCT_ROW_OR_JOINT_MOMENT_MISMATCH";
        proven += receipt.sameTimeMapProven;
    }
    if (bound > 0 || !writerRows.empty())
        trace << "\nZHANG_GRAPH_PRODUCT_WRITER_ROW_BIND time="
              << state.time.to_string(0)
              << " system=" << enum_to_string(system)
              << " bound=" << bound
              << " same_time_joint_proven=" << proven
              << " projection_status=" << projectionReason;
}

/** Numerically audit the entire phase/clock product block after one atomic
 * graph transform. The pre-event moments are already frozen in receipts, so
 * this never copies the network-sized KFState or covariance. Every included
 * phase row is the exact old physical functional expressed on the new chart;
 * this is not a claim that an unbridged new writer row has equal covariance. */
void validateGraphProductEventJointMoments(
    Trace& trace,
    const KFState& postState,
    vector<ZhangGraphProductEventReceipt>& receipts)
{
    vector<size_t> eligible;
    for (size_t i = 0; i < receipts.size(); i++)
        if (receipts[i].representableSameTimeMapProven)
            eligible.push_back(i);
    if (eligible.empty()) return;
    auto reject = [&](const string& reason)
    {
        for (const auto index : eligible)
        {
            auto& receipt = receipts[index];
            receipt.fullJointMomentTransportProven = false;
            receipt.jointMomentReason = reason;
            receipt.representableSameTimeMapProven = false;
            receipt.sameTimeMapProven = false;
            receipt.actualWriterRowMatches = false;
            receipt.conditionalProductMapProven = false;
            receipt.failureReason = reason;
            receipt.conditionalMapReason = reason;
        }
        trace << "\nZHANG_GRAPH_PRODUCT_JOINT_MOMENTS time="
              << postState.time.to_string(0)
              << " system=" << enum_to_string(receipts.front().system)
              << " rows=" << eligible.size()
              << " full_mean_covariance_proven=0 reason=" << reason;
    };
    if (postState.P.rows() != postState.x.size()
        || postState.P.cols() != postState.x.size())
    {
        reject("JOINT_PRODUCT_POSTERIOR_COVARIANCE_SHAPE_INVALID");
        return;
    }
    struct JointRow
    {
        vector<pair<int, double>> terms;
        size_t receiptIndex = 0;
        bool phase = false;
        double oldMean = 0;
    };
    vector<JointRow> rows;
    auto indexRow = [&](const ZhangCoordinateRow& symbolic,
                        JointRow& row) -> bool
    {
        for (const auto& [key, coefficient] : symbolic)
        {
            if (coefficient == 0) continue;
            const auto found = postState.kfIndexMap.find(key);
            if (found == postState.kfIndexMap.end()
                || found->second < 0
                || found->second >= postState.x.size()) return false;
            row.terms.emplace_back(found->second, coefficient);
        }
        return !row.terms.empty();
    };
    map<SatSys, size_t> clockSource;
    for (const auto index : eligible)
    {
        const auto& receipt = receipts[index];
        JointRow phase;
        phase.receiptIndex = index;
        phase.phase = true;
        phase.oldMean = receipt.oldPhaseMeanMetres;
        if (!indexRow(receipt.newPhaseRow, phase))
        {
            reject("JOINT_PHASE_DESTINATION_MISSING");
            return;
        }
        rows.push_back(std::move(phase));
        clockSource.emplace(receipt.satellite, index);
    }
    for (const auto& [satellite, index] : clockSource)
    {
        KFKey clock;
        clock.type = KF::SAT_CLOCK;
        clock.Sat = satellite;
        JointRow row;
        row.receiptIndex = index;
        row.oldMean = receipts[index].oldClockMeanMetres;
        if (!indexRow({{clock, 1}}, row))
        {
            reject("JOINT_CLOCK_DESTINATION_MISSING");
            return;
        }
        rows.push_back(std::move(row));
    }
    auto projectMean = [&](const JointRow& row)
    {
        double value = 0;
        for (const auto& [index, coefficient] : row.terms)
            value += coefficient * postState.x(index);
        return value;
    };
    auto projectCovariance = [&](const JointRow& first,
                                 const JointRow& second)
    {
        double value = 0;
        for (const auto& [a, ac] : first.terms)
        for (const auto& [b, bc] : second.terms)
            value += ac * bc * postState.P(a, b);
        return value;
    };
    auto oldCovariance = [&](const JointRow& first,
                             const JointRow& second,
                             double& value) -> bool
    {
        const auto& a = receipts[first.receiptIndex];
        const auto& b = receipts[second.receiptIndex];
        if (first.phase && second.phase)
        {
            const auto found = a.oldPhaseCrossCovarianceMetres2.find(
                {b.satellite, b.observable});
            if (found == a.oldPhaseCrossCovarianceMetres2.end())
                return false;
            value = found->second;
        }
        else if (first.phase || second.phase)
        {
            const auto& phase = first.phase ? a : b;
            const auto& clock = first.phase ? b : a;
            const auto found =
                phase.oldPhaseClockCrossCovarianceMetres2.find(
                    clock.satellite);
            if (found ==
                phase.oldPhaseClockCrossCovarianceMetres2.end())
                return false;
            value = found->second;
        }
        else
        {
            const auto found =
                a.oldClockClockCrossCovarianceMetres2.find(b.satellite);
            if (found == a.oldClockClockCrossCovarianceMetres2.end())
                return false;
            value = found->second;
        }
        return std::isfinite(value);
    };
    const int dimension = static_cast<int>(rows.size());
    VectorXd beforeMean(dimension);
    VectorXd afterMean(dimension);
    MatrixXd beforeCovariance(dimension, dimension);
    MatrixXd afterCovariance(dimension, dimension);
    for (int i = 0; i < dimension; i++)
    {
        beforeMean(i) = rows[i].oldMean;
        afterMean(i) = projectMean(rows[i]);
    }
    for (size_t i = 0; i < rows.size(); i++)
    for (size_t j = i; j < rows.size(); j++)
    {
        double before = 0;
        if (!oldCovariance(rows[i], rows[j], before))
        {
            reject("JOINT_PRODUCT_OLD_CROSS_COVARIANCE_MISSING");
            return;
        }
        const double after = projectCovariance(rows[i], rows[j]);
        beforeCovariance(i, j) = beforeCovariance(j, i) = before;
        afterCovariance(i, j) = afterCovariance(j, i) = after;
    }
    const auto agreement = zhangVerifyGraphJointMoments(
        beforeMean, beforeCovariance, afterMean, afterCovariance);
    if (!agreement.proven)
    {
        reject(agreement.failureReason);
        return;
    }
    for (const auto index : eligible)
    {
        receipts[index].fullJointMomentTransportProven = true;
        receipts[index].jointMomentReason =
            "COMPLETE_PHASE_CLOCK_MEAN_COVARIANCE_MATCH";
    }
    trace << "\nZHANG_GRAPH_PRODUCT_JOINT_MOMENTS time="
          << postState.time.to_string(0)
          << " system=" << enum_to_string(receipts.front().system)
          << " rows=" << rows.size()
          << " max_mean_error_m=" << agreement.maximumMeanError
          << " max_covariance_error_m2="
          << agreement.maximumCovarianceError
          << " full_mean_covariance_proven=1"
          << " published_new_functional_equivalence=0"
          << " integer_bridge_proven=0";
}

KFKey zhangReceiverPhaseKey(
    E_Sys       sys,
    E_ObsCode   code,
    const string& receiver
)
{
    KFKey key;
    key.type = KF::PHASE_BIAS;
    key.str  = receiver;
    key.Sat  = SatSys(sys, 0);
    key.num  = static_cast<int>(code);
    return key;
}

KFKey zhangSatellitePhaseKey(E_ObsCode code, const SatSys& satellite)
{
    KFKey key;
    key.type = KF::PHASE_BIAS;
    key.Sat  = satellite;
    key.num  = static_cast<int>(code);
    return key;
}

KFKey zhangAmbiguityKey(E_ObsCode code, const ZhangGraphEdge& edge)
{
    KFKey key;
    key.type = KF::AMBIGUITY;
    key.str  = edge.receiver;
    key.Sat  = edge.satellite;
    key.num  = static_cast<int>(code);
    return key;
}

void addExpression(
    map<KFKey, double>&       destination,
    const map<KFKey, double>& source,
    double                    scale = 1
)
{
    for (const auto& [key, coefficient] : source)
    {
        addCoefficient(destination, key, scale * coefficient);
    }
}

struct ZhangStateTransformAudit
{
    bool   evaluated = false;
    int    observableRows = 0;
    double meanAbsoluteNorm = std::numeric_limits<double>::infinity();
    double meanMaximumAbsolute = std::numeric_limits<double>::infinity();
    double covarianceAbsoluteNorm = std::numeric_limits<double>::infinity();
    double covarianceRelativeNorm = std::numeric_limits<double>::infinity();
    string failureReason = "NONE";
    string observable = "NONE";
    string receiver = "NONE";
    string satellite = "NONE";
    string missingKfKey = "NONE";
    int oldNodeCount = 0;
    int newNodeCount = 0;
    int representedEdgeCount = 0;
    int stateReconstructibleEdgeCount = 0;
};

/** Return the pre-event represented edges that can actually be reconstructed
 * from every requested observable in the immutable KF state.  A union of the
 * per-frequency ambiguity graphs is insufficient for a dual-frequency pivot:
 * the same physical edge must have both L1 and L2 state expressions. */
set<ZhangGraphEdge> zhangStateReconstructibleEdges(
    const KFState&                 kfState,
    E_Sys                          sys,
    const vector<E_ObsCode>&       baselineObservables,
    const ZhangGraphBasis&         oldBasis,
    const set<ZhangGraphEdge>&     representedEdges)
{
    auto present = [&](const KFKey& key) {
        return kfState.kfIndexMap.contains(key) && kfState.stateTransitionMap.contains(key);
    };
    vector<set<ZhangGraphEdge>> perObservableEdges;
    for (E_ObsCode code : baselineObservables)
    {
        set<ZhangGraphEdge> observableEdges;
        for (const auto& edge : representedEdges)
        {
            if (zhangWavelength(sys, code) <= 0)
            {
                continue;
            }
            if (edge.receiver != oldBasis.rootReceiver &&
                !present(zhangReceiverPhaseKey(sys, code, edge.receiver)))
            {
                continue;
            }
            if (!present(zhangSatellitePhaseKey(code, edge.satellite)))
            {
                continue;
            }
            if (!oldBasis.isTreeEdge(edge.receiver, edge.satellite) &&
                !present(zhangAmbiguityKey(code, edge)))
            {
                continue;
            }
            observableEdges.insert(edge);
        }
        perObservableEdges.push_back(std::move(observableEdges));
    }
    return zhangJointObservableEdgeIntersection(perObservableEdges);
}

/** Re-express the existing Zhang phase state in a new spanning-tree basis. */
bool transformZhangGraphBasis(
    Trace&                   trace,
    KFState&                 kfState,
    E_Sys                    sys,
    const vector<E_ObsCode>& baselineObservables,
    const ZhangGraphBasis&   oldBasis,
    const ZhangGraphBasis&   newTree,
    bool                     recordProductContinuity = true,
    const set<ZhangGraphEdge>& retiredEdges = {},
    ZhangStateTransformAudit*  audit = nullptr,
    KFState*                   rollbackSnapshot = nullptr,
    const map<ZhangGraphEdge, set<E_ObsCode>>* retiredSignals = nullptr
)
{
    if (zhangR51SingleBasisHandoffEnabled())
    {
        const auto owner = graphStateMap.find(
            {zhangGraphRuntimeId(kfState), sys});
        if (owner == graphStateMap.end() ||
            !owner->second.productInitialized ||
            owner->second.requiresFreshProductIdentity)
        {
            // A migrated schema-3 checkpoint has valid FLOAT x/P and its
            // old S-coordinate basis, but no authoritative per-signal arc
            // versions or published product branch. Re-express FLOAT with
            // the exact state transform below; there is no old integer
            // product to preserve or versioned auxiliary to certify.
            recordProductContinuity = false;
        }
    }
    auto signalRetired = [&](const ZhangGraphEdge& edge, E_ObsCode code)
    {
        if (!retiredEdges.contains(edge)) return false;
        if (!retiredSignals) return true;
        const auto signals = retiredSignals->find(edge);
        return signals == retiredSignals->end() ||
            signals->second.contains(code);
    };
    std::optional<KFState> preTransformState;
    const KFState* auditPreTransformState = nullptr;
    if (rollbackSnapshot)
    {
        *rollbackSnapshot = kfState;
        auditPreTransformState = rollbackSnapshot;
    }
    else if (audit)
    {
        preTransformState.emplace(kfState);
        auditPreTransformState = &*preTransformState;
    }
    if (audit)
    {
        *audit = {};
    }
    auto fail = [&](const string& reason,
                    E_ObsCode code = E_ObsCode::NONE,
                    const ZhangGraphEdge* edge = nullptr,
                    const KFKey* missingKey = nullptr)
    {
        if (audit)
        {
            audit->failureReason = reason;
            audit->observable = enum_to_string(code);
            audit->oldNodeCount = static_cast<int>(
                oldBasis.receivers.size() + oldBasis.satellites.size());
            audit->newNodeCount = static_cast<int>(
                newTree.receivers.size() + newTree.satellites.size());
            audit->representedEdgeCount = static_cast<int>(newTree.edges.size());
            if (edge)
            {
                audit->receiver = edge->receiver;
                audit->satellite = edge->satellite.id();
            }
            if (missingKey)
            {
                audit->missingKfKey = missingKey->commaString();
            }
        }
        return false;
    };
    vector<map<KFKey, double>> preObservableRows;
    vector<map<KFKey, double>> postObservableRows;
    for (const auto& edge : newTree.treeEdges)
        if (retiredEdges.contains(edge))
            return fail("RETIRED_ARC_CANNOT_SUPPLY_NEW_TREE", E_ObsCode::NONE, &edge);
    map<std::pair<SatSys, E_ObsCode>, double> oldSatellitePhases;
    for (const auto& [key, index] : kfState.kfIndexMap)
    {
        if (key.type == KF::PHASE_BIAS &&
            key.Sat.sys == sys &&
            key.Sat.prn > 0 &&
            key.str.empty() &&
            zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                baselineObservables
            ))
        {
            oldSatellitePhases[
                {key.Sat, static_cast<E_ObsCode>(key.num)}
            ] = kfState.x(index);
        }
    }

    map<KFKey, map<KFKey, double>> transform;
    map<std::pair<E_ObsCode, ZhangGraphEdge>, KFKey>
        retainedOldArcAuxiliaryKeys;

    auto isTargetPhase = [&](const KFKey& key)
    {
        if (key.type != KF::PHASE_BIAS ||
            !zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                baselineObservables
            ))
        {
            return false;
        }

        return key.Sat.sys == sys;
    };

    auto isTargetAmbiguity = [&](const KFKey& key)
    {
        return key.type == KF::AMBIGUITY &&
               key.Sat.sys == sys &&
               zhangFullRankUsesObservable(
                   static_cast<E_ObsCode>(key.num),
                   baselineObservables
               );
    };

    for (const auto& [key, index] : kfState.kfIndexMap)
    {
        if (!isTargetPhase(key) && !isTargetAmbiguity(key))
        {
            transform[key][key] = 1;
        }
    }

    for (E_ObsCode code : baselineObservables)
    {
        const double wavelength = zhangWavelength(sys, code);
        if (wavelength <= 0)
        {
            return fail("OBSERVABLE_WAVELENGTH_UNAVAILABLE", code);
        }

        set<string> stateReceivers = {oldBasis.rootReceiver};
        set<SatSys> stateSatellites;
        set<ZhangGraphEdge> ambiguityEdges;

        for (const auto& [key, index] : kfState.kfIndexMap)
        {
            if (key.num != static_cast<int>(code))
            {
                continue;
            }

            if (key.type == KF::PHASE_BIAS && key.Sat.sys == sys)
            {
                if (!key.str.empty() && key.Sat.prn == 0)
                {
                    stateReceivers.insert(key.str);
                }
                else if (key.str.empty() && key.Sat.prn > 0)
                {
                    stateSatellites.insert(key.Sat);
                }
            }
            else if (key.type == KF::AMBIGUITY && key.Sat.sys == sys)
            {
                ambiguityEdges.insert({key.str, key.Sat});
            }
        }

        set<ZhangGraphEdge> modelledEdges = oldBasis.treeEdges;
        modelledEdges.insert(ambiguityEdges.begin(), ambiguityEdges.end());

        // Nodes with no current baseline edge may be retired exactly by omitting their phase and
        // incident cycle states from the destination coordinate system.  Adding a genuinely new
        // node is handled separately as a leaf extension, because no old-state expression exists
        // for its phase datum.
        if (!std::includes(
                stateReceivers.begin(),
                stateReceivers.end(),
                newTree.receivers.begin(),
                newTree.receivers.end()
            ) ||
            !std::includes(
                stateSatellites.begin(),
                stateSatellites.end(),
                newTree.satellites.begin(),
                newTree.satellites.end()
            ))
        {
            BOOST_LOG_TRIVIAL(warning)
                << "Zhang graph transform cannot introduce a new phase node for "
                << enum_to_string(sys) << " " << enum_to_string(code);
            return fail("NEW_TREE_NODE_SET_MISMATCH", code);
        }

        for (auto edgeIt = modelledEdges.begin(); edgeIt != modelledEdges.end();)
        {
            if (newTree.receivers.find(edgeIt->receiver) == newTree.receivers.end() ||
                newTree.satellites.find(edgeIt->satellite) == newTree.satellites.end())
            {
                edgeIt = modelledEdges.erase(edgeIt);
            }
            else
            {
                ++edgeIt;
            }
        }

        const set<string>& receivers = newTree.receivers;
        const set<SatSys>& satellites = newTree.satellites;
        set<string> modelledReceivers;
        set<SatSys> modelledSatellites;
        for (const auto& edge : modelledEdges)
        {
            modelledReceivers.insert(edge.receiver);
            modelledSatellites.insert(edge.satellite);
        }

        if (modelledReceivers != receivers ||
            modelledSatellites != satellites)
        {
            BOOST_LOG_TRIVIAL(warning)
                << "Zhang graph transform node-set mismatch for "
                << enum_to_string(sys) << " " << enum_to_string(code)
                << ": state receivers/satellites=" << stateReceivers.size() << "/"
                << stateSatellites.size()
                << ", modelled=" << modelledReceivers.size() << "/"
                << modelledSatellites.size()
                << ", new tree=" << newTree.receivers.size() << "/"
                << newTree.satellites.size();
            return fail("OLD_AND_NEW_NODE_SET_MISMATCH", code);
        }

        for (const auto& edge : newTree.treeEdges)
        {
            if (modelledEdges.find(edge) == modelledEdges.end())
            {
                BOOST_LOG_TRIVIAL(warning)
                    << "Zhang graph transform lacks new tree edge "
                    << edge.receiver << "/" << edge.satellite.id() << " "
                    << enum_to_string(code);
                return fail("NEW_TREE_EDGE_NOT_REPRESENTED", code, &edge);
            }
        }

        auto oldEdgeExpression = [&](const ZhangGraphEdge& edge,
                                     map<KFKey, double>& expression)
        {
            if (edge.receiver != oldBasis.rootReceiver)
            {
                const KFKey receiverKey =
                    zhangReceiverPhaseKey(sys, code, edge.receiver);
                if (!addRequiredCoefficient(
                        kfState,
                        expression,
                        receiverKey,
                        +1
                    ))
                {
                    if (audit)
                    {
                        audit->failureReason = "OLD_NODE_STATE_MISSING";
                        audit->observable = enum_to_string(code);
                        audit->receiver = edge.receiver;
                        audit->satellite = edge.satellite.id();
                        audit->missingKfKey = receiverKey.commaString();
                    }
                    return false;
                }
            }

            const KFKey satelliteKey =
                zhangSatellitePhaseKey(code, edge.satellite);
            if (!addRequiredCoefficient(
                    kfState,
                    expression,
                    satelliteKey,
                    +1
                ))
            {
                if (audit)
                {
                    audit->failureReason = "OLD_NODE_STATE_MISSING";
                    audit->observable = enum_to_string(code);
                    audit->receiver = edge.receiver;
                    audit->satellite = edge.satellite.id();
                    audit->missingKfKey = satelliteKey.commaString();
                }
                return false;
            }

            if (!oldBasis.isTreeEdge(edge.receiver, edge.satellite))
            {
                const KFKey ambiguityKey = zhangAmbiguityKey(code, edge);
                if (!addRequiredCoefficient(
                        kfState,
                        expression,
                        ambiguityKey,
                        wavelength
                    ))
                {
                    if (audit)
                    {
                        audit->failureReason = "OBSERVABLE_STATE_MISSING";
                        audit->observable = enum_to_string(code);
                        audit->receiver = edge.receiver;
                        audit->satellite = edge.satellite.id();
                        audit->missingKfKey = ambiguityKey.commaString();
                    }
                    return false;
                }
            }

            return true;
        };

        map<ZhangGraphEdge, map<KFKey, double>> edgeExpressions;
        for (const auto& edge : modelledEdges)
        {
            if (!oldEdgeExpression(edge, edgeExpressions[edge]))
            {
                BOOST_LOG_TRIVIAL(warning)
                    << "Zhang graph transform cannot reconstruct old edge "
                    << edge.receiver << "/" << edge.satellite.id() << " "
                    << enum_to_string(code);
                return fail("OLD_EDGE_EXPRESSION_UNAVAILABLE", code, &edge);
            }
        }

        if (audit)
        {
            for (const auto& edge : newTree.edges)
            {
                if (signalRetired(edge, code)) continue;
                auto preExpression = edgeExpressions.find(edge);
                if (preExpression == edgeExpressions.end())
                {
                    return fail("PRE_AUDIT_EXPRESSION_MISSING", code, &edge);
                }
                preObservableRows.push_back(preExpression->second);

                map<KFKey, double> postExpression;
                if (edge.receiver != newTree.rootReceiver)
                {
                    addCoefficient(
                        postExpression,
                        zhangReceiverPhaseKey(sys, code, edge.receiver),
                        +1);
                }
                addCoefficient(
                    postExpression,
                    zhangSatellitePhaseKey(code, edge.satellite),
                    +1);
                if (!newTree.isTreeEdge(edge.receiver, edge.satellite))
                {
                    addCoefficient(
                        postExpression,
                        zhangAmbiguityKey(code, edge),
                        wavelength);
                }
                postObservableRows.push_back(std::move(postExpression));
            }
        }

        map<string, map<KFKey, double>> receiverExpressions;
        map<SatSys, map<KFKey, double>> satelliteExpressions;
        set<string> knownReceivers = {newTree.rootReceiver};
        set<SatSys> knownSatellites;

        bool progress = true;
        while (progress)
        {
            progress = false;
            for (const auto& edge : newTree.treeEdges)
            {
                bool receiverKnown =
                    knownReceivers.find(edge.receiver) != knownReceivers.end();
                bool satelliteKnown =
                    knownSatellites.find(edge.satellite) != knownSatellites.end();

                if (receiverKnown && !satelliteKnown)
                {
                    auto expression = edgeExpressions.at(edge);
                    addExpression(expression, receiverExpressions[edge.receiver], -1);
                    satelliteExpressions[edge.satellite] = std::move(expression);
                    knownSatellites.insert(edge.satellite);
                    progress = true;
                }
                else if (!receiverKnown && satelliteKnown)
                {
                    auto expression = edgeExpressions.at(edge);
                    addExpression(expression, satelliteExpressions[edge.satellite], -1);
                    receiverExpressions[edge.receiver] = std::move(expression);
                    knownReceivers.insert(edge.receiver);
                    progress = true;
                }
            }
        }

        if (knownReceivers != receivers || knownSatellites != satellites)
        {
            BOOST_LOG_TRIVIAL(warning)
                << "Zhang graph transform new tree is not connected in state space for "
                << enum_to_string(sys) << " " << enum_to_string(code)
                << ": reached receivers/satellites=" << knownReceivers.size() << "/"
                << knownSatellites.size()
                << ", expected=" << receivers.size() << "/" << satellites.size();
            return fail("NEW_TREE_NOT_CONNECTED_IN_STATE_SPACE", code);
        }

        for (const auto& receiver : receivers)
        {
            if (receiver == newTree.rootReceiver)
            {
                continue;
            }
            transform[zhangReceiverPhaseKey(sys, code, receiver)] =
                receiverExpressions.at(receiver);
        }

        for (const auto& satellite : satellites)
        {
            transform[zhangSatellitePhaseKey(code, satellite)] =
                satelliteExpressions.at(satellite);
        }

        for (const auto& edge : modelledEdges)
        {
            if (newTree.isTreeEdge(edge.receiver, edge.satellite))
            {
                continue;
            }

            // The expression is deliberately constructed above while the
            // pre-event arc still exists.  Omitting only its destination here
            // makes pivot+retire one atomic T*x / T*P*T' operation: common
            // phase states and every surviving cross-covariance are retained,
            // while the invalid old ambiguity coordinate is not copied into
            // the post-event state.
            if (signalRetired(edge, code))
            {
                continue;
            }

            auto expression = edgeExpressions.at(edge);
            addExpression(expression, receiverExpressions[edge.receiver], -1);
            addExpression(expression, satelliteExpressions[edge.satellite], -1);

            map<KFKey, double> ambiguityExpression;
            addExpression(ambiguityExpression, expression, 1 / wavelength);
            transform[zhangAmbiguityKey(code, edge)] = std::move(ambiguityExpression);
        }
    }

    if (zhangR51SingleBasisHandoffEnabled() && recordProductContinuity)
    {
        const auto previous = graphEventPreSignalArcVersions.find(
            {zhangGraphRuntimeId(kfState), sys});
        const auto owner = graphStateMap.find(
            {zhangGraphRuntimeId(kfState), sys});
        set<ZhangGraphEdge> auxiliaryCandidates = oldBasis.treeEdges;
        if (owner != graphStateMap.end())
            auxiliaryCandidates.insert(
                owner->second.productBasis.edges.begin(),
                owner->second.productBasis.edges.end());
        for (const auto& dependency : zhangHandoffFloatDependencies(kfState))
            if (dependency.edge.satellite.sys == sys &&
                oldBasis.edges.contains(dependency.edge))
                auxiliaryCandidates.insert(dependency.edge);
        for (const auto code : baselineObservables)
        for (const auto& edge : auxiliaryCandidates)
        {
            // Older public support may already be represented by a retained
            // auxiliary rather than the current source graph. Do not invent a
            // current physical row for that historical version.
            if (!oldBasis.edges.contains(edge)) continue;
            if (!signalRetired(edge, code)
                && newTree.edges.contains(edge)) continue;
            const auto oldSignal = previous == graphEventPreSignalArcVersions.end()
                ? std::optional<int>{}
                : [&]() -> std::optional<int>
                {
                    const auto edgeVersion = previous->second.find(edge);
                    if (edgeVersion == previous->second.end())
                        return std::nullopt;
                    const auto codeVersion = edgeVersion->second.find(int(code));
                    return codeVersion == edgeVersion->second.end()
                        ? std::nullopt
                        : std::optional<int>(codeVersion->second);
                }();
            if (!oldSignal
                || !previous->second.contains(edge)
                || !oldBasis.edges.contains(edge))
            {
                if (oldBasis.treeEdges.contains(edge))
                    return fail("OLD_PRODUCT_OR_TREE_ARC_VERSION_MISSING", code, &edge);
                trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time=" << kfState.time.to_string(0)
                      << " edge=" << edge.receiver << '/' << edge.satellite.id()
                      << " observable=" << enum_to_string(code)
                      << " status=NOT_RETAINED reason=SOURCE_SIGNAL_VERSION_UNAVAILABLE";
                continue;
            }
            const KFKey auxiliary = zhangGraphOldArcAuxiliaryKey(
                kfState, code, edge, *oldSignal);
            KFKey alreadyRetained;
            if (zhangHandoffFloatFind(kfState,
                    {edge, int(code), *oldSignal}, alreadyRetained))
            {
                retainedOldArcAuxiliaryKeys[{code, edge}] = alreadyRetained;
                continue;
            }
            if (transform.contains(auxiliary)
                || kfState.kfIndexMap.contains(auxiliary))
                return fail("OLD_TREE_ARC_AUXILIARY_DUPLICATE", code, &edge);
            ZhangCoordinateRow oldPhysical;
            if (edge.receiver != oldBasis.rootReceiver)
            {
                const KFKey receiver = zhangTransportReceiverKey(
                    sys, code, edge.receiver);
                if (!kfState.kfIndexMap.contains(receiver))
                    return fail("OLD_TREE_ARC_RECEIVER_MISSING", code, &edge);
                oldPhysical[receiver] = 1;
            }
            const KFKey satellite = zhangTransportSatelliteKey(
                code, edge.satellite);
            if (!kfState.kfIndexMap.contains(satellite))
                return fail("OLD_TREE_ARC_SATELLITE_MISSING", code, &edge);
            oldPhysical[satellite] = 1;
            if (!oldBasis.treeEdges.contains(edge))
            {
                const KFKey ambiguity = zhangTransportAmbiguityKey(
                    code, edge);
                if (!kfState.kfIndexMap.contains(ambiguity))
                {
                    trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time=" << kfState.time.to_string(0)
                          << " edge=" << edge.receiver << '/' << edge.satellite.id()
                          << " observable=" << enum_to_string(code)
                          << " status=NOT_RETAINED reason=SOURCE_CHORD_SIGNAL_STATE_UNAVAILABLE";
                    continue;
                }
                oldPhysical[ambiguity] = zhangWavelength(sys, code);
            }
            transform[auxiliary] = std::move(oldPhysical);
            retainedOldArcAuxiliaryKeys[{code, edge}] = auxiliary;
        }
    }

    map<string, string> captureMetadata;
    string captureFailureReason;
    if (!limitGraphHandoffFloatCapture(
            trace, kfState, transform, retainedOldArcAuxiliaryKeys,
            captureMetadata, captureFailureReason))
        return fail(captureFailureReason.c_str());

    // Validate every audit destination and its pullback before committing any
    // KF/factor mutation. In particular, retired chords have no destination.
    // A failed post-commit lookup must never leave a new posterior with an old graph.
    for (std::size_t i = 0; i < postObservableRows.size(); i++)
    {
        map<KFKey, double> pulledBack;
        for (const auto& [key, coefficient] : postObservableRows[i])
        {
            const auto destination = transform.find(key);
            if (destination == transform.end())
                return fail("AUDIT_DESTINATION_MISSING_BEFORE_COMMIT", E_ObsCode::NONE, nullptr, &key);
            addExpression(pulledBack, destination->second, coefficient);
        }
        addExpression(pulledBack, preObservableRows[i], -1);
        for (const auto& [key, coefficient] : pulledBack)
            if (std::abs(coefficient) > 1e-10)
                return fail("AUDIT_PULLBACK_MISMATCH_BEFORE_COMMIT", E_ObsCode::NONE, nullptr, &key);
    }

    vector<ZhangGraphProductEventReceipt> eventReceipts;
    if (zhangR51SingleBasisHandoffEnabled() && recordProductContinuity)
    {
        const auto owner = graphStateMap.find(
            {zhangGraphRuntimeId(kfState), sys});
        map<ZhangGraphEdge, map<int, int>> oldVersions;
        map<ZhangGraphEdge, map<int, int>> currentVersions;
        const auto oldIdentity = graphEventPreSignalArcVersions.find(
            {zhangGraphRuntimeId(kfState), sys});
        if (oldIdentity != graphEventPreSignalArcVersions.end())
            oldVersions = oldIdentity->second;
        if (owner != graphStateMap.end())
            currentVersions = owner->second.phaseArcVersions;
        for (const auto& edge : retiredEdges)
            for (const auto code : baselineObservables)
                if (signalRetired(edge, code))
                {
                    const auto oldEdge = oldVersions.find(edge);
                    if (oldEdge == oldVersions.end() ||
                        !oldEdge->second.contains(int(code))) continue;
                    currentVersions[edge][int(code)] =
                        std::max(currentVersions[edge][int(code)],
                            oldEdge->second.at(int(code)) + 1);
                }
        ZhangGraphCoordinatePlan productPlan;
        productPlan.valid = true;
        productPlan.transform = transform;
        productPlan.retainedOldArcAuxiliaryKeys =
            retainedOldArcAuxiliaryKeys;
        eventReceipts = prepareGraphProductEventReceipts(
            kfState, sys, baselineObservables, oldBasis,
            owner == graphStateMap.end()
                ? ZhangGraphBasis{} : owner->second.productBasis,
            newTree,
            oldVersions, currentVersions, retiredEdges, productPlan,
            retiredEdges.empty() ? "TREE_REBASE" : "TREE_ARC_RETIREMENT");
    }

    string label =
        "Zhang graph " + enum_to_string(sys) +
        (retiredEdges.empty()
            ? " tree exchange"
            : " pivot before retired arc projection");
    bool applied = kfState.applyStateTransform(trace, transform, label, {}, &captureMetadata);
    if (!applied)
    {
        return fail("STATE_TRANSFORM_APPLICATION_FAILED");
    }

    if (audit)
    {
        using IndexedRow = vector<std::pair<int, double>>;
        auto indexRows = [](
            const KFState& state,
            const vector<map<KFKey, double>>& symbolic,
            vector<IndexedRow>& indexed)
        {
            indexed.clear();
            indexed.reserve(symbolic.size());
            for (const auto& row : symbolic)
            {
                IndexedRow result;
                for (const auto& [key, coefficient] : row)
                {
                    auto found = state.kfIndexMap.find(key);
                    if (found == state.kfIndexMap.end())
                    {
                        return false;
                    }
                    result.emplace_back(found->second, coefficient);
                }
                indexed.push_back(std::move(result));
            }
            return true;
        };

        vector<IndexedRow> preRows;
        vector<IndexedRow> postRows;
        if (!auditPreTransformState ||
            !indexRows(*auditPreTransformState, preObservableRows, preRows) ||
            !indexRows(kfState, postObservableRows, postRows) ||
            preRows.size() != postRows.size())
        {
            if (auditPreTransformState)
                static_cast<KFState_&>(kfState) = static_cast<const KFState_&>(*auditPreTransformState);
            return fail("POST_AUDIT_EXPRESSION_MISSING");
        }

        auto projectedMean = [](const KFState& state, const IndexedRow& row)
        {
            double value = 0;
            for (const auto& [index, coefficient] : row)
            {
                value += coefficient * state.x(index);
            }
            return value;
        };
        auto projectedCovariance = [](
            const KFState& state,
            const IndexedRow& first,
            const IndexedRow& second)
        {
            double value = 0;
            for (const auto& [firstIndex, firstCoefficient] : first)
            for (const auto& [secondIndex, secondCoefficient] : second)
            {
                value += firstCoefficient * secondCoefficient *
                    state.P(firstIndex, secondIndex);
            }
            return value;
        };

        long double meanSquared = 0;
        double meanMaximum = 0;
        for (size_t row = 0; row < preRows.size(); row++)
        {
            const double difference =
                projectedMean(kfState, postRows[row]) -
                projectedMean(*auditPreTransformState, preRows[row]);
            meanSquared += static_cast<long double>(difference) * difference;
            meanMaximum = std::max(meanMaximum, std::abs(difference));
        }

        long double covarianceDifferenceSquared = 0;
        long double covarianceReferenceSquared = 0;
        for (size_t row = 0; row < preRows.size(); row++)
        for (size_t column = 0; column < preRows.size(); column++)
        {
            const double before = projectedCovariance(
                *auditPreTransformState, preRows[row], preRows[column]);
            const double after = projectedCovariance(
                kfState, postRows[row], postRows[column]);
            const long double difference =
                static_cast<long double>(after) - before;
            covarianceDifferenceSquared += difference * difference;
            covarianceReferenceSquared +=
                static_cast<long double>(before) * before;
        }

        audit->evaluated = true;
        audit->observableRows = static_cast<int>(preRows.size());
        audit->meanAbsoluteNorm =
            std::sqrt(static_cast<double>(meanSquared));
        audit->meanMaximumAbsolute = meanMaximum;
        audit->covarianceAbsoluteNorm =
            std::sqrt(static_cast<double>(covarianceDifferenceSquared));
        audit->covarianceRelativeNorm =
            covarianceReferenceSquared > 0
                ? std::sqrt(static_cast<double>(
                    covarianceDifferenceSquared /
                    covarianceReferenceSquared))
                : audit->covarianceAbsoluteNorm;
    }

    // Raw capture and remapped whole-target metadata were committed together
    // inside applyStateTransform; never overwrite the latter with the old map.
    if (!eventReceipts.empty())
        validateGraphProductEventJointMoments(
            trace, kfState, eventReceipts);
    if (!eventReceipts.empty())
        storeGraphProductEventReceipts(
            kfState, sys, std::move(eventReceipts));

    map<E_ObsCode, map<SatSys, double>> correctionChanges;
    for (const auto& [key, oldPhase] : oldSatellitePhases)
    {
        const auto& [satellite, code] = key;
        KFKey phaseKey = zhangSatellitePhaseKey(code, satellite);
        auto phaseIt = kfState.kfIndexMap.find(phaseKey);
        if (phaseIt == kfState.kfIndexMap.end())
        {
            continue;
        }

        double newPhase = kfState.x(phaseIt->second);
        correctionChanges[code][satellite] = -(newPhase - oldPhase);
    }
    for (const auto& [code, changes] : correctionChanges)
    {
        if (recordProductContinuity)
        {
            recordZhangExactPhaseTransforms(
                kfState.time, sys, code, changes
            );
        }
    }

    return true;
}

bool transportZhangGraphComponents(
    Trace& trace, KFState& state, E_Sys system,
    const vector<E_ObsCode>& observables, const ZhangGraphBasis& oldBasis,
    const set<ZhangGraphEdge>& eligible,
    const set<ZhangGraphEdge>& commonEligible,
    const map<int, set<ZhangGraphEdge>>& signalEligible,
    const set<ZhangGraphEdge>& retired,
    const map<ZhangGraphEdge, set<E_ObsCode>>& retiredSignals,
    ZhangGraphBasis& transported, ZhangGraphCoordinatePlan& plan)
{
    const auto reconstructible = zhangStateReconstructibleEdges(
        state, system, observables, oldBasis, oldBasis.edges);
    auto edges = eligible;
    for (const auto& edge : reconstructible)
        if (!zhangR51SingleBasisHandoffEnabled()
            && !retired.contains(edge)) edges.insert(edge);
    auto preferred = oldBasis.treeEdges;
    for (const auto& edge : retired) preferred.erase(edge);
    if (zhangR51SingleBasisHandoffEnabled())
    {
        // The shared phase datum is determined only by current observations
        // valid on every baseline signal. A one-frequency edge can augment
        // an already connected signal cycle, never create a common tree edge
        // or connect a new node/component on behalf of the other frequency.
        transported = zhangBuildSpanningTree(
            commonEligible, oldBasis.rootReceiver, preferred, {},
            reconstructible);
        set<ZhangGraphEdge> signalChords;
        std::set_difference(edges.begin(), edges.end(),
            commonEligible.begin(), commonEligible.end(),
            std::inserter(signalChords, signalChords.begin()));
        transported = zhangAugmentCommonTreeWithSignalChords(
            std::move(transported), signalChords);
    }
    else
        transported = zhangBuildSpanningTree(
            edges, oldBasis.rootReceiver, preferred, {}, reconstructible);
    map<int, set<ZhangGraphEdge>> destinationSignalEdges;
    if (zhangR51SingleBasisHandoffEnabled())
        for (const auto code : observables)
        {
            auto& destination = destinationSignalEdges[int(code)];
            for (const auto& edge : transported.edges)
                if (commonEligible.contains(edge) ||
                    (signalEligible.contains(int(code)) &&
                     signalEligible.at(int(code)).contains(edge)))
                    destination.insert(edge);
        }
    map<E_ObsCode, double> wavelengths;
    for (const auto code : observables) wavelengths[code] = zhangWavelength(system, code);
    map<ZhangGraphEdge, int> oldArcVersions;
    map<ZhangGraphEdge, map<int, int>> oldSignalArcVersions;
    set<ZhangGraphEdge> completeProductSupport;
    const auto owner = graphStateMap.find(
        {zhangGraphRuntimeId(state), system});
    const bool retainEventAuxiliaries =
        zhangR51SingleBasisHandoffEnabled() &&
        owner != graphStateMap.end() &&
        owner->second.productInitialized &&
        !owner->second.requiresFreshProductIdentity;
    if (retainEventAuxiliaries)
    {
        completeProductSupport = owner->second.productBasis.edges;
        for (auto it = completeProductSupport.begin();
             it != completeProductSupport.end();)
            it = oldBasis.edges.contains(*it) ? std::next(it)
                : completeProductSupport.erase(it);
        for (const auto& dependency : zhangHandoffFloatDependencies(state))
            if (dependency.edge.satellite.sys == system &&
                oldBasis.edges.contains(dependency.edge))
                completeProductSupport.insert(dependency.edge);
        const auto previous = graphEventPreArcVersions.find(
            {zhangGraphRuntimeId(state), system});
        if (previous != graphEventPreArcVersions.end())
            oldArcVersions = previous->second;
        const auto signalPrevious = graphEventPreSignalArcVersions.find(
            {zhangGraphRuntimeId(state), system});
        if (signalPrevious != graphEventPreSignalArcVersions.end())
            oldSignalArcVersions = signalPrevious->second;
    }
    plan = zhangPlanGraphCoordinateTransport(
        state, system, oldBasis, transported, retired, wavelengths,
        [&](const ZhangGraphEdge& edge, E_ObsCode code)
        {
            auto& options = acsConfig.getRecOpts(
                edge.receiver, {edge.satellite.sysName(), enum_to_string(code)});
            return initialStateFromConfig(options.ambiguity).P;
        },
        [&](const string& receiver, E_ObsCode code)
        {
            auto& options = acsConfig.getRecOpts(
                receiver, {SatSys(system, 0).sysName(), enum_to_string(code)});
            return initialStateFromConfig(options.phase_bias).P;
        }, &retiredSignals, retainEventAuxiliaries,
        retainEventAuxiliaries ? &oldArcVersions : nullptr,
        retainEventAuxiliaries ? &completeProductSupport : nullptr,
        retainEventAuxiliaries ? &oldSignalArcVersions : nullptr,
        zhangR51SingleBasisHandoffEnabled()
            ? &destinationSignalEdges : nullptr);
    if (!plan.valid)
    {
        trace << "\nZHANG_COMPONENT_TRANSPORT_REJECT time=" << state.time.to_string(0)
              << " reason=" << plan.failureReason << " posterior_preserved=1";
        return false;
    }
    for (const auto code : observables)
    for (const auto& edge : completeProductSupport)
    {
        const bool retiredForSignal = retired.contains(edge) &&
            (!retiredSignals.contains(edge) || retiredSignals.at(edge).contains(code));
        if ((retiredForSignal || !transported.edges.contains(edge)) &&
            !plan.retainedOldArcAuxiliaryKeys.contains({code, edge}))
            trace << "\nZHANG_HANDOFF_FLOAT_CAPTURE time=" << state.time.to_string(0)
                  << " edge=" << edge.receiver << '/' << edge.satellite.id()
                  << " observable=" << enum_to_string(code)
                  << " status=NOT_RETAINED reason=SOURCE_PHYSICAL_SIGNAL_UNREPRESENTABLE";
    }
    map<string, string> captureMetadata;
    string captureFailureReason;
    if (!limitGraphHandoffFloatCapture(trace, state, plan.transform,
            plan.retainedOldArcAuxiliaryKeys, captureMetadata, captureFailureReason))
    {
        state.lastFactorTransactionFailureReason = captureFailureReason;
        return false;
    }
    vector<ZhangGraphProductEventReceipt> eventReceipts;
    if (retainEventAuxiliaries)
    {
        map<ZhangGraphEdge, map<int, int>> oldVersions =
            oldSignalArcVersions;
        map<ZhangGraphEdge, map<int, int>> currentVersions;
        if (owner != graphStateMap.end())
            currentVersions = owner->second.phaseArcVersions;
        for (const auto& edge : retired)
            for (const auto code : observables)
            {
                const auto signals = retiredSignals.find(edge);
                if (signals != retiredSignals.end() &&
                    !signals->second.contains(code)) continue;
                const auto old = oldVersions.find(edge);
                if (old == oldVersions.end() ||
                    !old->second.contains(int(code))) continue;
                currentVersions[edge][int(code)] =
                    std::max(currentVersions[edge][int(code)],
                        old->second.at(int(code)) + 1);
            }
        eventReceipts = prepareGraphProductEventReceipts(
            state, system, observables, oldBasis,
            owner == graphStateMap.end()
                ? ZhangGraphBasis{} : owner->second.productBasis,
            transported,
            oldVersions, currentVersions, retired, plan,
            "COMPONENT_REPARAMETERISATION");
    }
    if (!state.applyStateTransform(trace, plan.transform,
            "ZHANG_COMPONENT_REPARAMETERISATION", plan.independentSourceVariances,
            &captureMetadata))
        return false;
    if (!eventReceipts.empty())
        validateGraphProductEventJointMoments(
            trace, state, eventReceipts);
    if (!eventReceipts.empty())
        storeGraphProductEventReceipts(
            state, system, std::move(eventReceipts));
    trace << "\nZHANG_COMPONENT_TRANSPORT time=" << state.time.to_string(0)
          << " component_count=" << transported.componentCount
          << " retained_physical_rows=" << plan.survivingPhysicalRows.size()
          << " fresh_edges=" << plan.freshEdges.size()
          << " independent_prior_sources=" << plan.independentSourceVariances.size()
          << " component_root_count=" << plan.componentRoots.size()
          << " symbolic_prediction_identity=1 full_cross_covariance_preserved=1"
          << " legacy_float_rebase="
          << (zhangR51SingleBasisHandoffEnabled() &&
              !retainEventAuxiliaries)
          << " old_integer_auxiliary_authority="
          << retainEventAuxiliaries
          << " integer_certification=0 observation_used_for_initialisation=0";
    return true;
}

bool resetZhangGraphPhaseCoordinates(
    Trace&                   trace,
    KFState&                 kfState,
    E_Sys                    sys,
    const vector<E_ObsCode>& baselineObservables,
    const set<string>&        affectedReceivers,
    const set<SatSys>&        affectedSatellites
)
{
    map<KFKey, map<KFKey, double>> transform;
    int removedStates = 0;
    for (const auto& [key, index] : kfState.kfIndexMap)
    {
        bool targetReceiverPhase =
            key.type == KF::PHASE_BIAS &&
            key.Sat.sys == sys &&
            key.Sat.prn == 0 &&
            !key.str.empty() &&
            affectedReceivers.find(key.str) != affectedReceivers.end() &&
            zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                baselineObservables
            );
        bool targetSatellitePhase =
            key.type == KF::PHASE_BIAS &&
            key.Sat.sys == sys &&
            key.Sat.prn > 0 &&
            key.str.empty() &&
            affectedSatellites.find(key.Sat) != affectedSatellites.end() &&
            zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                baselineObservables
            );
        bool targetAmbiguity =
            key.type == KF::AMBIGUITY &&
            key.Sat.sys == sys &&
            (affectedReceivers.find(key.str) != affectedReceivers.end() ||
             affectedSatellites.find(key.Sat) != affectedSatellites.end()) &&
            zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                baselineObservables
            );

        if (!targetReceiverPhase && !targetSatellitePhase && !targetAmbiguity)
        {
            transform[key][key] = 1;
        }
        else
        {
            removedStates++;
        }
    }

    return removedStates > 0 &&
           !transform.empty() &&
           kfState.applyStateTransform(
               trace,
               transform,
               "Zhang graph local phase-coordinate reinitialisation"
           );
}

bool transformZhangDatum(
    Trace&                       trace,
    KFState&                     kfState,
    E_Sys                        sys,
    const vector<E_ObsCode>&     baselineObservables,
    const string&                oldReceiver,
    const SatSys&                oldSatellite,
    const string&                newReceiver,
    const SatSys&                newSatellite
)
{
    map<KFKey, map<KFKey, double>> transform;

    auto isZhangPhaseState = [&](const KFKey& key)
    {
        if (key.type != KF::PHASE_BIAS ||
            !zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                baselineObservables
            ))
        {
            return false;
        }

        if (!key.str.empty())
        {
            return key.Sat.sys == sys && key.Sat.prn == 0;
        }

        return key.Sat.sys == sys && key.Sat.prn > 0;
    };

    auto isZhangAmbiguity = [&](const KFKey& key)
    {
        return key.type == KF::AMBIGUITY &&
               key.Sat.sys == sys &&
               zhangFullRankUsesObservable(
                   static_cast<E_ObsCode>(key.num),
                   baselineObservables
               );
    };

    for (auto& [key, index] : kfState.kfIndexMap)
    {
        const bool clockState =
            ((key.type == KF::REC_CLOCK || key.type == KF::REC_CLOCK_RATE) &&
             zhangReceiverClockMatchesSystem(acsConfig, key, sys)) ||
            ((key.type == KF::SAT_CLOCK || key.type == KF::SAT_CLOCK_RATE) &&
             key.Sat.sys == sys);

        if (clockState || isZhangPhaseState(key) || isZhangAmbiguity(key))
        {
            continue;
        }

        transform[key][key] = 1;
    }

    if (!transformClockDatum(
            kfState,
            transform,
            sys,
            oldReceiver,
            newReceiver,
            KF::REC_CLOCK,
            KF::SAT_CLOCK
        ) ||
        !transformClockDatum(
            kfState,
            transform,
            sys,
            oldReceiver,
            newReceiver,
            KF::REC_CLOCK_RATE,
            KF::SAT_CLOCK_RATE
        ))
    {
        BOOST_LOG_TRIVIAL(error)
            << "Cannot change Zhang receiver reference from " << oldReceiver << " to "
            << newReceiver << ": required receiver clock state is absent";
        return false;
    }

    for (E_ObsCode code : baselineObservables)
    {
        auto sysFrequencyIt = code2Freq.find(sys);
        if (sysFrequencyIt == code2Freq.end())
        {
            return false;
        }
        auto frequencyIt = sysFrequencyIt->second.find(code);
        if (frequencyIt == sysFrequencyIt->second.end())
        {
            return false;
        }
        auto wavelengthIt = genericWavelength.find(frequencyIt->second);
        if (wavelengthIt == genericWavelength.end())
        {
            return false;
        }
        const double wavelength = wavelengthIt->second;

        set<string> receivers = {oldReceiver};
        set<SatSys> satellites;
        set<std::pair<string, SatSys>> ambiguityEdges;

        for (auto& [key, index] : kfState.kfIndexMap)
        {
            if (key.type == KF::PHASE_BIAS && key.num == static_cast<int>(code))
            {
                if (!key.str.empty() && key.Sat.sys == sys && key.Sat.prn == 0)
                {
                    receivers.insert(key.str);
                }
                else if (key.str.empty() && key.Sat.sys == sys && key.Sat.prn > 0)
                {
                    satellites.insert(key.Sat);
                }
            }
            else if (
                key.type == KF::AMBIGUITY &&
                key.num == static_cast<int>(code) &&
                key.Sat.sys == sys
            )
            {
                ambiguityEdges.emplace(key.str, key.Sat);
            }
        }

        auto receiverPhaseKey = [&](const string& receiver)
        {
            KFKey key;
            key.type = KF::PHASE_BIAS;
            key.str  = receiver;
            key.Sat  = SatSys(sys, 0);
            key.num  = static_cast<int>(code);
            return key;
        };

        auto satellitePhaseKey = [&](const SatSys& satellite)
        {
            KFKey key;
            key.type = KF::PHASE_BIAS;
            key.Sat  = satellite;
            key.num  = static_cast<int>(code);
            return key;
        };

        auto ambiguityKey = [&](const string& receiver, const SatSys& satellite)
        {
            KFKey key;
            key.type = KF::AMBIGUITY;
            key.str  = receiver;
            key.Sat  = satellite;
            key.num  = static_cast<int>(code);
            return key;
        };

        if (receivers.find(newReceiver) == receivers.end() ||
            satellites.find(newSatellite) == satellites.end())
        {
            BOOST_LOG_TRIVIAL(debug)
                << "Zhang S-transform candidate lacks a phase tree state for "
                << newReceiver << "/" << newSatellite.id() << " " << enum_to_string(code);
            return false;
        }

        // The current filter may be sparse after outages and arc resets.  Its phase states form
        // a graph: the old reference receiver row and satellite column are the tree edges, while
        // each retained DD ambiguity is a non-tree edge.  Reject internally inconsistent dormant
        // ambiguities rather than silently inventing a missing receiver/satellite phase state.
        for (auto& [receiver, satellite] : ambiguityEdges)
        {
            if (receivers.find(receiver) == receivers.end() ||
                satellites.find(satellite) == satellites.end())
            {
                BOOST_LOG_TRIVIAL(debug)
                    << "Zhang S-transform found an ambiguity without both phase tree states: "
                    << receiver << "/" << satellite.id() << " " << enum_to_string(code);
                return false;
            }
        }

        auto edgeExpression = [&](const string& receiver,
                                  const SatSys& satellite,
                                  map<KFKey, double>& coefficients)
        {
            if (receiver != oldReceiver &&
                !addRequiredCoefficient(
                    kfState,
                    coefficients,
                    receiverPhaseKey(receiver),
                    +1
                ))
            {
                return false;
            }

            if (!addRequiredCoefficient(
                    kfState,
                    coefficients,
                    satellitePhaseKey(satellite),
                    +1
                ))
            {
                return false;
            }

            if (receiver != oldReceiver && satellite != oldSatellite)
            {
                if (ambiguityEdges.find({receiver, satellite}) == ambiguityEdges.end() ||
                    !addRequiredCoefficient(
                        kfState,
                        coefficients,
                        ambiguityKey(receiver, satellite),
                        wavelength
                    ))
                {
                    return false;
                }
            }

            return true;
        };

        map<SatSys, map<KFKey, double>> newSatelliteExpressions;
        for (auto& satellite : satellites)
        {
            if (!edgeExpression(
                    newReceiver,
                    satellite,
                    newSatelliteExpressions[satellite]
                ))
            {
                BOOST_LOG_TRIVIAL(debug)
                    << "Zhang S-transform candidate lacks receiver tree edge "
                    << newReceiver << "/" << satellite.id() << " " << enum_to_string(code);
                return false;
            }

            transform[satellitePhaseKey(satellite)] = newSatelliteExpressions[satellite];
        }

        map<string, map<KFKey, double>> newReceiverExpressions;
        for (auto& receiver : receivers)
        {
            if (receiver == newReceiver)
            {
                continue;
            }

            map<KFKey, double> receiverAtNewSatellite;
            map<KFKey, double> newReferenceAtNewSatellite;
            if (!edgeExpression(receiver, newSatellite, receiverAtNewSatellite) ||
                !edgeExpression(newReceiver, newSatellite, newReferenceAtNewSatellite))
            {
                BOOST_LOG_TRIVIAL(debug)
                    << "Zhang S-transform candidate lacks satellite tree edge "
                    << receiver << "/" << newSatellite.id() << " " << enum_to_string(code);
                return false;
            }

            auto& expression = newReceiverExpressions[receiver];
            for (auto& [source, coefficient] : receiverAtNewSatellite)
            {
                addCoefficient(expression, source, coefficient);
            }
            for (auto& [source, coefficient] : newReferenceAtNewSatellite)
            {
                addCoefficient(expression, source, -coefficient);
            }

            transform[receiverPhaseKey(receiver)] = expression;
        }

        // Re-express exactly the graph edges represented by the old state.  Removing the new tree
        // edges and retaining all remaining graph edges gives the same number of independent DD
        // states even when the receiver/satellite network is not a complete Cartesian product.
        set<std::pair<string, SatSys>> modelledEdges = ambiguityEdges;
        for (auto& satellite : satellites)
        {
            modelledEdges.emplace(oldReceiver, satellite);
        }
        for (auto& receiver : receivers)
        {
            modelledEdges.emplace(receiver, oldSatellite);
        }

        for (auto& [receiver, satellite] : modelledEdges)
        {
            if (receiver == newReceiver || satellite == newSatellite)
            {
                continue;
            }

            map<KFKey, double> edge;
            if (!edgeExpression(receiver, satellite, edge))
            {
                return false;
            }

            auto& expression = transform[ambiguityKey(receiver, satellite)];
            for (auto& [source, coefficient] : edge)
            {
                addCoefficient(expression, source, coefficient / wavelength);
            }
            for (auto& [source, coefficient] : newReceiverExpressions[receiver])
            {
                addCoefficient(expression, source, -coefficient / wavelength);
            }
            for (auto& [source, coefficient] : newSatelliteExpressions[satellite])
            {
                addCoefficient(expression, source, -coefficient / wavelength);
            }
        }
    }

    string label =
        "Zhang " + enum_to_string(sys) + " " + oldReceiver + "/" + oldSatellite.id() + " -> " +
        newReceiver + "/" + newSatellite.id();

    if (transform.size() != kfState.kfIndexMap.size())
    {
        BOOST_LOG_TRIVIAL(debug)
            << "Zhang S-transform is not dimension preserving for " << label << ": "
            << kfState.kfIndexMap.size() << " -> " << transform.size();
        return false;
    }

    return kfState.applyStateTransform(trace, transform, label);
}

string certifiedTransportMetadataKey(E_Sys system)
{
    return string(zhang_transport_certificate::metadataKey) + "/" + enum_to_string(system);
}

bool loadCertifiedTransportDirectory(const KFState& state, E_Sys system,
    zhang_transport_certificate::Directory& directory, string& reason)
{
    const auto found = state.metaDataMap.find(certifiedTransportMetadataKey(system));
    if (found == state.metaDataMap.end())
    { reason = "NO_COMMITTED_TRANSPORT_DIRECTORY"; return false; }
    const auto& envelope = found->second;
    if (envelope.size() < 65 || envelope[64] != '\n' ||
        zhangCheckpointSha256(envelope.substr(65)) != envelope.substr(0, 64))
    { reason = "TRANSPORT_DIRECTORY_SEAL_MISMATCH"; return false; }
    if (!zhang_transport_certificate::Directory::decode(envelope.substr(65),
        zhangGraphRuntimeId(state), static_cast<int>(system), directory, reason)) return false;
    if (directory.lastCommittedEpoch >= static_cast<long int>(state.time.bigTime))
    { reason = "TRANSPORT_DIRECTORY_NOT_A_PRIOR_COMMIT"; return false; }
    return true;
}

struct CertifiedTransportQueryCache
{
    string seal;
    long int readEpoch = 0;
    std::size_t envelopeBytes = 0;
    zhang_transport_certificate::Directory directory;
    std::shared_ptr<const zhang_transport_certificate::Query> query;
    vector<string> proofIds;
    std::shared_ptr<const vector<string>> proofStatements;
};

std::shared_ptr<const CertifiedTransportQueryCache> certifiedTransportQuery(
    const KFState& state, E_Sys system, string& reason)
{
    static map<std::pair<string, E_Sys>, std::shared_ptr<const CertifiedTransportQueryCache>> cache;
    const auto metadata = state.metaDataMap.find(certifiedTransportMetadataKey(system));
    if (metadata == state.metaDataMap.end()) { reason = "NO_COMMITTED_TRANSPORT_DIRECTORY"; return {}; }
    const string owner = zhangGraphRuntimeId(state);
    const auto epoch = static_cast<long int>(state.time.bigTime);
    auto& prior = cache[{owner, system}];
    // The owned metadata is immutable during graph/AR execution. A changed
    // owner, epoch, seal or size must pass the full read-only integrity gate.
    if (prior && prior->readEpoch == epoch && prior->envelopeBytes == metadata->second.size() &&
        prior->seal == metadata->second.substr(0, 64))
    { reason = prior->query->reason; return prior; }
    auto built = std::make_shared<CertifiedTransportQueryCache>();
    if (!loadCertifiedTransportDirectory(state, system, built->directory, reason)) return {};
    built->seal = metadata->second.substr(0, 64); built->envelopeBytes = metadata->second.size(); built->readEpoch = epoch;
    built->query = std::make_shared<const zhang_transport_certificate::Query>(built->directory);
    auto statements = std::make_shared<vector<string>>();
    auto add = [&](const string& id, const string& source, const string& payload)
    {
        std::ostringstream statement;
        statement << "id=" << std::quoted(id) << " source=" << std::quoted(source)
            << " statement=" << std::quoted(payload);
        built->proofIds.push_back(id); statements->push_back(statement.str());
    };
    for (const auto& [id, scope] : built->directory.sourceScopes)
    {
        std::ostringstream payload;
        payload << "TRANSPORT_ONLY_SOURCE epoch=" << scope.committedEpoch << " posterior=" << std::quoted(scope.posterior)
            << " chart=" << std::quoted(scope.chart) << " float_root=" << std::quoted(scope.floatRoot)
            << " state_order=" << std::quoted(scope.stateOrder);
        add("TRANSPORT_SOURCE|" + id, owner + "|PRIOR_PRODUCT_COMMIT", payload.str());
    }
    for (const auto& [id, payload] : built->directory.proofAtoms)
        add("TRANSPORT_DECISION|" + id, owner + "|PRIOR_PRODUCT_COMMIT", payload);
    for (std::size_t i = 0; i < built->directory.facts.size(); ++i)
    {
        const auto& fact = built->directory.facts[i];
        std::ostringstream payload;
        payload << "TRANSPORT_ONLY_PHYSICAL_RELATION source_scope=" << std::quoted(fact.sourceScopeId) << " terms=";
        for (const auto& [id, coefficient] : fact.physical) payload << std::quoted(id) << ':' << coefficient << ',';
        payload << " rhs=" << fact.value << " accepted_parents=";
        for (const auto& id : fact.proofIds) payload << std::quoted("TRANSPORT_DECISION|" + id) << ',';
        add("TRANSPORT_FACT|" + built->seal + "|" + std::to_string(i), owner + "|PRIOR_PRODUCT_COMMIT", payload.str());
    }
    built->proofStatements = std::move(statements);
    reason = built->query->reason;
    prior = built;
    return built;
}

zhang_transport_certificate::TreePlan certifiedTransportTreePlan(
    Trace& trace, const KFState& state, E_Sys system, const GraphRuntimeState& runtime,
    const set<ZhangGraphEdge>& available, const set<ZhangGraphEdge>& broken,
    const set<ZhangGraphEdge>& nonContinuousCandidates,
    const string& root, const map<ZhangGraphEdge, double>& quality,
    const vector<E_ObsCode>& codes,
    const map<ZhangGraphEdge, map<int, int>>& oldVersions,
    const map<int, set<ZhangGraphEdge>>& signalEdges)
{
    using namespace zhang_transport_certificate;
    string reason;
    const auto cached = certifiedTransportQuery(state, system, reason);
    const auto query = cached ? cached->query : nullptr;
    map<int, DirectionContract> perSignal;
    map<int, bool> completeSignalTargets;
    auto certify = [&](const ZhangGraphEdge& edge)
    {
        // A fresh observation of an endpoint-identical edge is not its old
        // physical arc. Pending breaks precede version consumption, so an
        // equality check against the old version map would be insufficient.
        // Keep this candidate legal for FLOAT restart, but never rank it as
        // an old-COMMIT-certified replacement.
        if (nonContinuousCandidates.contains(edge)) return false;
        const auto cycle = zhangFundamentalCycle(runtime.basis, edge);
        bool allSignals = !codes.empty() && !cycle.empty() && query && query->feasible;
        for (const auto code : codes)
        {
            Row target;
            bool identitiesComplete = !cycle.empty();
            const auto currentSignal = signalEdges.find(static_cast<int>(code));
            // The candidate itself must still carry its pre-event physical
            // arc. Broken OLD tree terms remain valid historical variables.
            identitiesComplete &= currentSignal != signalEdges.end() && currentSignal->second.contains(edge);
            for (const auto& [term, coefficient] : cycle)
            {
                const auto versions = oldVersions.find(term);
                if (versions == oldVersions.end() || !versions->second.contains(static_cast<int>(code)))
                { identitiesComplete = false; break; }
                target[enum_to_string(code) + "|" + term.receiver + "|" + term.satellite.id() +
                    "|V" + std::to_string(versions->second.at(static_cast<int>(code)))] += coefficient;
            }
            auto& contract = perSignal[static_cast<int>(code)];
            completeSignalTargets.try_emplace(static_cast<int>(code), true);
            completeSignalTargets[static_cast<int>(code)] = completeSignalTargets[static_cast<int>(code)] && identitiesComplete;
            contract.owner = zhangGraphRuntimeId(state);
            contract.sourceEpoch = std::to_string(cached ? cached->directory.lastCommittedEpoch : 0);
            contract.targetEpoch = state.time.to_string(9);
            contract.requestedDirections.push_back(target);
            ZhangExactInteger value;
            const bool known = identitiesComplete && query && query->consequence(target, value);
            if (known) { contract.knownDirections.push_back(target); contract.knownValues.push_back(value); }
            else contract.unknownDirections.push_back(target);
            contract.reason = identitiesComplete ? (known ? "HISTORICAL_TRANSPORT_ONLY_CERTIFIED" :
                query && query->feasible ? "HISTORICAL_TRANSPORT_DIRECTION_UNCERTIFIED" : reason)
                                                : "NO_CONTINUOUS_SIGNAL_VERSION_OR_OLD_ENDPOINT";
            allSignals &= known;
        }
        return allSignals;
    };
    auto plan = select(runtime.basis, available, broken, root, quality, certify);
    if (query && query->feasible)
        for (auto& [signal, contract] : perSignal)
            if (completeSignalTargets[signal])
            {
                auto classified = query->classify(contract.requestedDirections);
                classified.owner = contract.owner;
                classified.sourceEpoch = contract.sourceEpoch;
                classified.targetEpoch = contract.targetEpoch;
                contract = std::move(classified);
            }
    vector<Row> jointTargets;
    bool completeJointTargets = !perSignal.empty();
    for (const auto& [signal, contract] : perSignal)
    {
        completeJointTargets &= completeSignalTargets[signal];
        jointTargets.insert(jointTargets.end(), contract.requestedDirections.begin(), contract.requestedDirections.end());
    }
    if (query && query->feasible && completeJointTargets)
    {
        plan.jointSignalDirections = query->classify(jointTargets);
        plan.jointSignalDirections.owner = zhangGraphRuntimeId(state);
        plan.jointSignalDirections.sourceEpoch = std::to_string(cached->directory.lastCommittedEpoch);
        plan.jointSignalDirections.targetEpoch = state.time.to_string(9);
    }
    plan.perSignal = std::move(perSignal);
    trace << "\nZHANG_CERTIFIED_TREE_TRANSPORT time=" << state.time.to_string(0)
        << " system=" << enum_to_string(system) << " owner=" << zhangGraphRuntimeId(state)
        << " source_role=PRIOR_COMMIT_TRANSPORT_ONLY current_ar_authority=0"
        << " certified_replacements=" << plan.certifiedReplacements.size()
        << " float_fallback_replacements=" << plan.fallbackReplacements.size()
        << " new_endpoint_extensions=" << plan.newEndpointExtensions.size()
        << " joint_ranks_evaluated=" << plan.jointSignalDirections.ranksEvaluated
        << " joint_known_rank=" << plan.jointSignalDirections.knownRank
        << " joint_unknown_rank=" << plan.jointSignalDirections.unknownRank
        << " all_old_endpoints_certified=" << plan.allOldEndpointsCertified << " reason=" << reason;
    for (const auto& [signal, contract] : plan.perSignal)
        trace << "\nZHANG_CERTIFIED_TREE_SIGNAL_DIRECTIONS time=" << state.time.to_string(0)
            << " system=" << enum_to_string(system) << " signal=" << enum_to_string(static_cast<E_ObsCode>(signal))
            << " known_directions=" << contract.knownDirections.size()
            << " unknown_directions=" << contract.unknownDirections.size()
            << " ranks_evaluated=" << contract.ranksEvaluated
            << " known_rank=" << contract.knownRank << " unknown_rank=" << contract.unknownRank
            << " current_ar_authority=0 reason=" << contract.reason;
    return plan;
}

void updateZhangGraphBasis(
    Trace&                             trace,
    KFState&                           kfState,
    E_Sys                              sys,
    const ZhangFullRankSystemOptions& options,
    const ReferenceAvailability&       availability
)
{
    const string runtimeId = zhangGraphRuntimeId(kfState);
    const bool hasEstimatedPhaseState = zhangHasEstimatedGraphPhaseState(
        kfState, sys, options.baseline_observables);
    if (!validZhangGraphRuntimeId(runtimeId))
    {
        BOOST_LOG_TRIVIAL(error)
            << "ZHANG_GRAPH_RUNTIME_ID_UNBOUND sys=" << enum_to_string(sys);
        if (hasEstimatedPhaseState)
            throw std::runtime_error("ZHANG_FLOAT_GRAPH_RUNTIME_ID_INVALID_WITH_ESTIMATED_PHASE_STATE");
        return;
    }
    const auto existingRuntime = graphStateMap.find({runtimeId, sys});
    if (hasEstimatedPhaseState &&
        (existingRuntime == graphStateMap.end() || !existingRuntime->second.initialized))
    {
        // A restored checkpoint imports its initialized S topology before
        // this entry point. Existing phase coordinates cannot be interpreted
        // in a newly selected tree without an explicit state transform.
        trace << "\nZHANG_FLOAT_GRAPH_AUTHORITY_LOST time=" << kfState.time.to_string(0)
              << " system=" << enum_to_string(sys)
              << " runtime=" << runtimeId
              << " estimated_phase_state=1 action=ABORT_BEFORE_GRAPH_MUTATION";
        throw std::runtime_error("ZHANG_FLOAT_GRAPH_RUNTIME_MISSING_WITH_ESTIMATED_PHASE_STATE");
    }
    auto& runtime = graphStateMap[{runtimeId, sys}];
    const string productMode = configuredProductBasisMode();
    const bool singleS = singleSBasisNative();
    if (runtime.initialized && runtime.productBasisMode != productMode)
        throw std::runtime_error("ZHANG_PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION");
    runtime.productBasisMode = productMode;
    const auto previousSGaugeSignalArcs = runtime.phaseArcVersions;
    for (const auto& [signal, edges] : availability.signalEdges)
        for (const auto& edge : edges)
            runtime.phaseArcVersions[edge].try_emplace(signal, 0);
    if (zhangR51SingleBasisHandoffEnabled())
    {
        auto& before = graphEventPreArcVersions[{runtimeId, sys}];
        before.clear();
        for (const auto& [edge, history] : runtime.edgeHistory)
            before[edge] = history.arcVersion;
        graphEventPreSignalArcVersions[{runtimeId, sys}] =
            runtime.phaseArcVersions;
    }
    runtime.epochIndex++;
    struct GraphTiming
    {
        Trace& trace;
        const KFState& state;
        const GraphRuntimeState& runtime;
        const ReferenceAvailability& availability;
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        ~GraphTiming()
        {
            if (!acsConfig.zhangPppAr.output_diagnostics) return;
            trace << "\nZHANG_GRAPH_TIMING time=" << state.time.to_string(0)
                  << " seconds=" << std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-start).count()
                  << " states=" << state.x.size()
                  << " raw_edges=" << availability.rawEdges.size()
                  << " eligible_phase_edges=" << availability.edges.size()
                  << " eligible_code_edges=" << availability.codeEdges.size()
                  << " admitted_phase_edges=" << runtime.observationEdges.size()
                  << " float_components=" << runtime.basis.componentCount;
        }
    } graphTiming{trace, kfState, runtime, availability};

    auto traceCanonicalAudit = [&](const ZhangGraphBasis& basis,
                                   const string&           action,
                                   bool                    exactTransition)
    {
        if (!acsConfig.zhangPppAr.output_diagnostics)
        {
            return;
        }

        ZhangCanonicalIntegerAudit audit = zhangCanonicalIntegerAudit(basis);
        if (!audit.valid)
        {
            trace << "\nZHANG_CANONICAL_INTEGER_AUDIT time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " action=" << action
                  << " valid=0 reason=canonical_graph_construction_failed"
                  << " detail=" << audit.failureReason;
            return;
        }

        const string componentId = zhangIntegerComponentId(basis);
        for (E_ObsCode code : options.baseline_observables)
        {
            trace << "\nZHANG_CANONICAL_INTEGER_AUDIT time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " signal=" << enum_to_string(code)
                  << " action=" << action
                  << " valid=1"
                  << " component_id=" << componentId
                  << " root_node=R:" << basis.rootReceiver
                  << " datum_version=" << runtime.datumVersion
                  << " representation_version="
                  << runtime.representationVersion
                  << " float_gauge_version="
                  << runtime.floatGaugeVersion
                  << " integer_component_version="
                  << runtime.integerComponentVersion
                  << " arcs=" << basis.edges.size()
                  << " tree_datum_integers=" << audit.treeEdges.size()
                  << " cycle_integers=" << audit.chordEdges.size()
                  << " satellite_datum_rows="
                  << audit.satelliteDatumSingleDifferences.size()
                  << " satellite_fix_quotient_rows="
                  << audit.satelliteFixQuotient.size()
                  << " satellite_fix_quotient_nonzero_rows=0"
                  << " canonical_to_arc_fingerprint="
                  << audit.canonicalToArcFingerprint
                  << " datum_mapping_fingerprint="
                  << audit.datumMappingFingerprint
                  << " fix_quotient_fingerprint="
                  << audit.fixQuotientFingerprint
                  << " dense_canonical_materialised="
                  << audit.denseCanonicalMaterialised
                  << " canonical_to_arc_unimodular=STRUCTURAL_UNIT_BLOCK"
                  << " exact_epoch_transition=" << exactTransition;

            for (const auto& edge : audit.treeEdges)
            {
                int arcVersion = runtime.edgeHistory[edge].arcVersion;
                trace << "\nZHANG_CANONICAL_INTEGER_COORDINATE time="
                      << kfState.time.to_string(0)
                      << " signal=" << enum_to_string(code)
                      << " component_id=" << componentId
                      << " type=TREE_DATUM"
                      << " integer_id=D:" << edge.receiver << ":"
                      << edge.satellite.id() << ":A" << arcVersion
                      << " arc_id=" << edge.receiver << ":"
                      << edge.satellite.id() << ":" << enum_to_string(code)
                      << ":A" << arcVersion;
            }
            for (const auto& edge : audit.chordEdges)
            {
                int arcVersion = runtime.edgeHistory[edge].arcVersion;
                trace << "\nZHANG_CANONICAL_INTEGER_COORDINATE time="
                      << kfState.time.to_string(0)
                      << " signal=" << enum_to_string(code)
                      << " component_id=" << componentId
                      << " type=CYCLE"
                      << " integer_id=K:" << edge.receiver << ":"
                      << edge.satellite.id() << ":A" << arcVersion
                      << " arc_id=" << edge.receiver << ":"
                      << edge.satellite.id() << ":" << enum_to_string(code)
                      << ":A" << arcVersion;
            }
        }
    };

    auto edgeList = [](const set<ZhangGraphEdge>& edges)
    {
        std::ostringstream stream;
        bool first = true;
        for (const auto& edge : edges)
        {
            stream << (first ? "" : ",") << edge.receiver << ":"
                   << edge.satellite.id();
            first = false;
        }
        return first ? string("NONE") : stream.str();
    };
    auto receiverList = [](const set<string>& receivers)
    {
        std::ostringstream stream;
        bool first = true;
        for (const auto& receiver : receivers)
        {
            stream << (first ? "" : ",") << receiver;
            first = false;
        }
        return first ? string("NONE") : stream.str();
    };
    auto satelliteList = [](const set<SatSys>& satellites)
    {
        std::ostringstream stream;
        bool first = true;
        for (const auto& satellite : satellites)
        {
            stream << (first ? "" : ",") << satellite.id();
            first = false;
        }
        return first ? string("NONE") : stream.str();
    };
    auto traceGraphEvent = [&](const string&              action,
                               const string&              reason,
                               const set<ZhangGraphEdge>& removedTreeEdges,
                               const set<ZhangGraphEdge>& replacementEdges,
                               const set<string>&         resetReceivers,
                               const set<SatSys>&         resetSatellites,
                               int                        removedIntegerColumns,
                               bool                       exactTransform)
    {
        runtime.eventCounter++;
        if (zhangR51SingleBasisHandoffEnabled())
        {
            const auto runtimeId = zhangGraphRuntimeId(kfState);
            const auto epoch = kfState.time.to_string(9);
            for (auto& [identity, receipt] : graphProductEventReceipts)
                if (std::get<0>(identity) == runtimeId &&
                    std::get<1>(identity) == epoch &&
                    std::get<2>(identity) == sys)
                {
                    receipt.newBackendGeneration = runtime.eventCounter;
                    receipt.newChartGeneration =
                        runtime.representationVersion;
                }
        }
        if (!acsConfig.zhangPppAr.output_diagnostics)
        {
            return;
        }
        trace << "\nZHANG_GRAPH_INTEGER_EVENT time="
              << kfState.time.to_string(0)
              << " system=" << enum_to_string(sys)
              << " event_id=" << runtime.eventCounter
              << " event_type=" << action
              << " reason=" << reason
              << " affected_tree_edges=" << edgeList(removedTreeEdges)
              << " replacement_edges=" << edgeList(replacementEdges)
              << " local_reset_nodes=" << receiverList(resetReceivers)
              << " local_reset_satellites=" << satelliteList(resetSatellites)
              << " removed_integer_columns=" << removedIntegerColumns
              << " representation_version=" << runtime.representationVersion
              << " float_gauge_version=" << runtime.floatGaugeVersion
              << " integer_component_version="
              << runtime.integerComponentVersion
              << " held_rows_touched=DEFERRED_TO_HELD_LATTICE_EVENT"
              << " held_rows_removed=DEFERRED_TO_HELD_LATTICE_EVENT"
              << " exact_unimodular_transform_available=" << exactTransform
              << " held_lattice_storage=PHYSICAL_ARC_VERSION_HNF";
    };

    // The filter state graph is deliberately allowed to outlive the instantaneous
    // observation graph.  Old tree edges and retained cycle ambiguities provide exact
    // integer coordinates across short data gaps, but an excluded slip is a hard arc break.
    set<ZhangGraphEdge> modelledEdges = runtime.basis.treeEdges;
    for (const auto& [key, index] : kfState.kfIndexMap)
    {
        if (key.type == KF::AMBIGUITY &&
            key.Sat.sys == sys &&
            zhangFullRankUsesObservable(
                static_cast<E_ObsCode>(key.num),
                options.baseline_observables
            ))
        {
            modelledEdges.insert({key.str, key.Sat});
        }
    }

    map<ZhangGraphEdge, int> preEventArcVersions;
    for (const auto& edge : availability.discontinuousEdges)
    {
        auto& history = runtime.edgeHistory[edge];
        preEventArcVersions[edge] = history.arcVersion;
        history.continuousEpochs = 0;
        history.outageEpochs = options.state_edge_grace_epochs + 1;
    }
    struct PhaseAdmissionGuard
    {
        GraphRuntimeState& runtime;
        const KFState& state;
        const ReferenceAvailability& availability;
        const map<ZhangGraphEdge, int>& previous;
        const vector<E_ObsCode>& codes;
        ~PhaseAdmissionGuard()
        {
            runtime.phaseAdmissions.clear();
            runtime.signalObservationEdges.clear();
            runtime.signalObservationFinalized = false;
            runtime.signalObservationFactorSequence = 0;
            runtime.phaseAdmissionTime = state.time.to_string(9);
            runtime.phaseAdmissionChart = runtime.representationVersion;
            set<ZhangGraphEdge> safe;
            const auto graphSupported = zhangGraphSupportedSignalEdges(
                availability.signalEdges, runtime.stateEdges,
                runtime.basis.edges);
            for (const auto code : codes)
            {
                const auto observed = graphSupported.find(int(code));
                if (observed == graphSupported.end()) continue;
                for (const auto& edge : observed->second)
                {
                    const auto before = previous.find(edge);
                    const auto after = runtime.edgeHistory.find(edge);
                    const bool retiredButUncommitted =
                        availability.discontinuousEdges.contains(edge) &&
                        (before == previous.end() || after == runtime.edgeHistory.end() ||
                         after->second.arcVersion <= before->second);
                    const auto changed = availability.discontinuitySignals.find(edge);
                    if (retiredButUncommitted && (changed == availability.discontinuitySignals.end() ||
                        changed->second.contains(code))) continue;
                    const auto physical = runtime.phaseArcVersions.find(edge);
                    if (physical == runtime.phaseArcVersions.end()) continue;
                    const auto signal = physical->second.find(int(code));
                    // An absent physical identity is unknown, not version 0.
                    if (signal == physical->second.end() || signal->second < 0) continue;
                    runtime.phaseAdmissions[edge][int(code)] = signal->second;
                    runtime.signalObservationEdges[int(code)].insert(edge);
                    const auto topology = runtime.edgeHistory.find(edge);
                    if (topology != runtime.edgeHistory.end())
                    {
                        auto& observation =
                            runtime.arcVersionObservationHistory[
                                {edge, topology->second.arcVersion}];
                        if (observation.firstObservedEpoch < 0)
                            observation.firstObservedEpoch =
                                runtime.epochIndex;
                        if (observation.lastObservedEpoch !=
                            runtime.epochIndex)
                        {
                            observation.lastObservedEpoch =
                                runtime.epochIndex;
                            observation.observedEpochs.push_back(
                                runtime.epochIndex);
                            observation.observationEpochs =
                                int(observation.observedEpochs.size());
                        }
                    }
                    auto& observation =
                        runtime.signalArcVersionObservationHistory[
                            {edge, int(code),
                             runtime.phaseArcVersions[edge][int(code)]}];
                    if (observation.firstObservedEpoch < 0)
                        observation.firstObservedEpoch = runtime.epochIndex;
                    if (observation.lastObservedEpoch != runtime.epochIndex)
                    {
                        observation.lastObservedEpoch = runtime.epochIndex;
                        observation.observedEpochs.push_back(
                            runtime.epochIndex);
                        observation.observationEpochs =
                            int(observation.observedEpochs.size());
                    }
                    safe.insert(edge);
                }
            }
            runtime.observationEdges = std::move(safe);
        }
    } admissionGuard{runtime, kfState, availability, preEventArcVersions, options.baseline_observables};
    // Constructed AFTER admissionGuard: on every exit the committed S view is
    // frozen first, then current per-signal admission uses that same catalogue
    // before PPP measurement construction.  No KF state/covariance writes are
    // performed by this guard.
    struct CommittedSViewGuard
    {
        Trace& trace;
        const KFState& state;
        E_Sys system;
        GraphRuntimeState& runtime;
        const ReferenceAvailability& availability;
        const vector<E_ObsCode>& codes;
        ZhangGraphBasis before;
        bool wasInitialized;
        int representation;
        int event;
        bool singleS;
        map<ZhangGraphEdge, map<int, int>> previousSignalArcs;
        int exceptions = std::uncaught_exceptions();
        ~CommittedSViewGuard()
        {
            if (std::uncaught_exceptions() != exceptions) return;
            const auto change = zhangClassifyCommittedSGraphEvent(
                before, runtime.basis, wasInitialized,
                !availability.discontinuousEdges.empty());
            const bool catalogueChanged = !zhangSameCommittedGraphBasis(before, runtime.basis);
            if (singleS)
            {
                bool gaugeChanged = before.rootReceiver != runtime.basis.rootReceiver ||
                    before.treeEdges != runtime.basis.treeEdges;
                if (!gaugeChanged)
                    for (const auto& edge : runtime.basis.treeEdges)
                    {
                        const auto oldArc = previousSignalArcs.find(edge);
                        const auto newArc = runtime.phaseArcVersions.find(edge);
                        if ((oldArc == previousSignalArcs.end()) !=
                                (newArc == runtime.phaseArcVersions.end()) ||
                            (oldArc != previousSignalArcs.end() &&
                             newArc != runtime.phaseArcVersions.end() && oldArc->second != newArc->second))
                            gaugeChanged = true;
                    }
                if (wasInitialized && gaugeChanged) ++runtime.committedSGaugeVersion;
                // Existing exact transports already increment their chart.
                // A chord insertion also changes the coordinate catalogue,
                // even though all old means/Q/cross blocks remain untouched.
                if (catalogueChanged && wasInitialized &&
                    runtime.representationVersion == representation)
                    ++runtime.representationVersion;
                // Product backend generation is NOT the catalogue revision.
                // A pure chord augmentation changes no existing S potential
                // or native product row: new-coordinate coefficients are
                // exactly zero. Keep its backend event/gauge unchanged. Real
                // transports/physical retirements already issue their event.
                if (wasInitialized && gaugeChanged && runtime.eventCounter == event)
                    ++runtime.eventCounter;
                syncSingleSBasisProductView(runtime);
                if (change.kind != "TOPOLOGY_UNCHANGED")
                    runtime.lastProductEventCause = change.kind;
                if (zhangR51SingleBasisHandoffEnabled())
                    bindGraphProductEventWriterRows(trace, state, system, codes,
                        runtime.basis, runtime);
            }
            trace << "\nZHANG_COMMITTED_S_GRAPH_EVENT time=" << state.time.to_string(0)
                  << " system=" << enum_to_string(system)
                  << " product_basis_mode=" << runtime.productBasisMode
                  << " kind=" << change.kind
                  << " physical_arc_break_requested=" << change.physicalArcBreak
                  << " chord_catalogue_changed=" << change.chordCatalogueChanged
                  << " node_augmentation=" << change.nodeAugmentation
                  << " tree_replacement=" << change.treeReplacement
                  << " component_split=" << change.componentSplit
                  << " component_reconnect=" << change.componentReconnect
                  << " added_edges=" << change.addedEdges.size()
                  << " removed_edges=" << change.removedEdges.size()
                  << " representation_version=" << runtime.representationVersion
                  << " committed_s_gauge_version=" << runtime.committedSGaugeVersion
                  << " old_coordinates_modified_by_catalogue_commit=0"
                     " current_signal_admission=SEPARATE_PRE_FACTOR";
        }
    } committedSGuard{trace, kfState, sys, runtime, availability,
        options.baseline_observables, runtime.basis, runtime.initialized,
        runtime.representationVersion, runtime.eventCounter, singleS, previousSGaugeSignalArcs};
    struct ArcEventConsumptionGuard
    {
        GraphRuntimeState& runtime;
        const ReferenceAvailability& availability;
        const map<ZhangGraphEdge, int>& previous;
        ~ArcEventConsumptionGuard()
        {
            for (const auto& [edge, version] : previous)
            {
                const auto current = runtime.edgeHistory.find(edge);
                if (current == runtime.edgeHistory.end() || current->second.arcVersion <= version)
                    continue;
                const auto stats = availability.discontinuityStats.find(edge);
                if (stats == availability.discontinuityStats.end()) continue;
                for (auto* stat : stats->second)
                {
                    stat->savedSlip.any = false;
                    stat->tracking.pendingBreak = false;
                    stat->phaseRejectCount = 0;
                }
            }
        }
    } consumeCommittedArcEvents{runtime, availability, preEventArcVersions};
    bool arcVersionsAdvanced = false;
    auto advanceDiscontinuousArcVersions = [&]()
    {
        if (arcVersionsAdvanced) return;
        for (const auto& edge : availability.discontinuousEdges)
            runtime.edgeHistory[edge].arcVersion++;
        for (const auto& [edge, signals] : availability.discontinuitySignals)
            if (availability.discontinuousEdges.contains(edge))
                for (const auto code : signals)
                    runtime.phaseArcVersions[edge][int(code)]++;
        arcVersionsAdvanced = true;
    };

    set<ZhangGraphEdge> trackedEdges = modelledEdges;
    trackedEdges.insert(availability.edges.begin(), availability.edges.end());
    for (const auto& [edge, history] : runtime.edgeHistory)
    {
        trackedEdges.insert(edge);
    }

    for (const auto& edge : trackedEdges)
    {
        if (availability.discontinuousEdges.find(edge) !=
            availability.discontinuousEdges.end())
        {
            continue;
        }

        auto& history = runtime.edgeHistory[edge];
        if (availability.edges.find(edge) != availability.edges.end())
        {
            history.continuousEpochs++;
            history.outageEpochs = 0;
            auto& observations = runtime.arcVersionObservationHistory[
                {edge, history.arcVersion}];
            if (observations.firstObservedEpoch < 0)
            {
                observations.firstObservedEpoch = runtime.epochIndex;
            }
            observations.lastObservedEpoch = runtime.epochIndex;
            observations.observationEpochs++;
            observations.observedEpochs.push_back(runtime.epochIndex);
        }
        else
        {
            history.outageEpochs++;
        }
    }

    for(const auto& code:options.baseline_observables) {
        ZhangGraphEdge edge{"KIRI",SatSys("G27")};
        auto hist=runtime.edgeHistory.find(edge);
        int lastObserved=-1;
        if(hist!=runtime.edgeHistory.end()) {
            auto obs=runtime.arcVersionObservationHistory.find({edge,hist->second.arcVersion});
            if(obs!=runtime.arcVersionObservationHistory.end())lastObserved=obs->second.lastObservedEpoch;
        }
        BOOST_LOG_TRIVIAL(info)<<"R48_STATE_AVAILABILITY time="<<kfState.time.to_string(0)
            <<" key=KIRI:G27:"<<enum_to_string(code)
            <<" phase=GRAPH_PRE_TRANSITION raw_observation_present="<<availability.rawEdges.count(edge)
            <<" qc_observation_valid="<<availability.edges.count(edge)
            <<" physical_arc_version="<<(hist==runtime.edgeHistory.end()?-1:hist->second.arcVersion)
            <<" last_valid_observation_epoch="<<lastObserved
            <<" last_valid_observation_time="<<(lastObserved<0?"UNOBSERVED":
                (kfState.time-(runtime.epochIndex-lastObserved)*acsConfig.epoch_interval).to_string(0))
            <<" discontinuity="<<availability.discontinuousEdges.count(edge)
            <<" qc_excluded="<<availability.qcExcludedEdges.count(edge)
            <<" signal_unavailable="<<availability.signalUnavailableEdges.count(edge)
            <<" graph_role="<<(runtime.basis.treeEdges.count(edge)?"TREE":runtime.basis.edges.count(edge)?"CHORD":"ABSENT")
            <<" cycle_chart_id="<<runtime.representationVersion
            <<" branch_id="<<zhangGraphRuntimeId(kfState);
    }
    set<ZhangGraphEdge> stateCandidates = availability.edges;
    const auto liveRepresentedEdges = zhangStateReconstructibleEdges(
        kfState, sys, options.baseline_observables, runtime.basis, modelledEdges);
    for (const auto& edge : modelledEdges)
    {
        auto historyIt = runtime.edgeHistory.find(edge);
        if (!zhangR51SingleBasisHandoffEnabled()
            && liveRepresentedEdges.contains(edge)
            && historyIt != runtime.edgeHistory.end() &&
            historyIt->second.outageEpochs <= options.state_edge_grace_epochs)
        {
            stateCandidates.insert(edge);
        }
    }
    // A post-slip observation is a new physical arc.  It may be introduced as
    // a fresh chord after the transaction, but it cannot participate in the
    // pre-event replacement tree or keep the old arc alive through grace.
    for (const auto& edge : availability.discontinuousEdges)
        stateCandidates.erase(edge);

    set<ZhangGraphEdge> stateEdges = stateCandidates;
    set<ZhangGraphEdge> observationEdges;
    std::set_intersection(
        availability.edges.begin(),
        availability.edges.end(),
        stateEdges.begin(),
        stateEdges.end(),
        std::inserter(observationEdges, observationEdges.begin())
    );

    auto retainOldTreeRootComponent =
        [&]()
        {
            set<ZhangGraphEdge> activeOldTreeEdges;
            std::set_intersection(
                stateEdges.begin(),
                stateEdges.end(),
                runtime.basis.treeEdges.begin(),
                runtime.basis.treeEdges.end(),
                std::inserter(activeOldTreeEdges, activeOldTreeEdges.begin())
            );

            set<ZhangGraphEdge> connectedTreeEdges =
                zhangRootComponentEdges(
                    activeOldTreeEdges,
                    runtime.basis.rootReceiver
                );

            set<string> connectedReceivers = {runtime.basis.rootReceiver};
            set<SatSys> connectedSatellites;
            for (const auto& edge : connectedTreeEdges)
            {
                connectedReceivers.insert(edge.receiver);
                connectedSatellites.insert(edge.satellite);
            }

            set<ZhangGraphEdge> safeEdges;
            for (const auto& edge : stateEdges)
            {
                if (connectedReceivers.find(edge.receiver) != connectedReceivers.end() &&
                    connectedSatellites.find(edge.satellite) != connectedSatellites.end())
                {
                    safeEdges.insert(edge);
                }
            }

            runtime.stateEdges = std::move(safeEdges);
            runtime.observationEdges.clear();
            std::set_intersection(
                observationEdges.begin(),
                observationEdges.end(),
                runtime.stateEdges.begin(),
                runtime.stateEdges.end(),
                std::inserter(
                    runtime.observationEdges,
                    runtime.observationEdges.begin()
                )
            );
            if (zhangR51SingleBasisHandoffEnabled() &&
                !availability.discontinuousEdges.empty())
            {
                // A failed physical-arc transaction has not consumed the
                // slip or advanced its version.  The retained FLOAT root
                // component is diagnostic only; the old product authority
                // cannot cross that unresolved event.
                invalidateZhangProduct(runtime,
                    "PHYSICAL_ARC_TRANSACTION_DEFERRED");
            }
        };

    const bool allCurrentArcsRetiredWithFreshObservations =
        zhangR51SingleBasisHandoffEnabled() && stateEdges.empty() &&
        !availability.edges.empty();
    if (stateEdges.empty() && !availability.edges.empty())
    {
        // All old arcs may have retired together. These are destination-only
        // fresh arcs: the pivot planner still excludes the complete retirement
        // set from old-coordinate transport. Holding here would leave saved
        // slip events pending forever and prevent the network from restarting.
        stateEdges = availability.edges;
        observationEdges = availability.edges;
    }
    if (stateEdges.empty())
    {
        // No admissible carrier rows does not justify destroying the previous
        // posterior. Keep pending physical events until a legal chart can be
        // formed; code observations still have their independent admission.
        runtime.observationEdges.clear();
        runtime.stateEdges.clear();
        invalidateZhangProduct(runtime, "NO_CURRENT_PHASE_SUPPORT");
        trace << "\nZHANG_GRAPH_EMPTY_SUPPORT time=" << kfState.time.to_string(0)
              << " pending_physical_events=" << availability.discontinuousEdges.size()
              << " action=HOLD_POSTERIOR_NO_PHASE_ROWS product_authorized=0";
        return;
    }

    map<ZhangGraphEdge, double> activeQuality;
    map<ZhangGraphEdge, int> persistence;
	map<ZhangGraphEdge, ZhangIntegerArcQuality> integerSupportQuality;
    for (const auto& edge : stateEdges)
    {
        auto qualityIt = availability.edgeQuality.find(edge);
        if (qualityIt != availability.edgeQuality.end())
        {
            activeQuality[edge] = qualityIt->second;
        }

        auto historyIt = runtime.edgeHistory.find(edge);
        if (historyIt != runtime.edgeHistory.end())
        {
            persistence[edge] = historyIt->second.continuousEpochs;
			ZhangIntegerArcQuality audit;
			audit.ageEpochs = historyIt->second.continuousEpochs;
			audit.outageCount = historyIt->second.outageEpochs;
			audit.slipCount = availability.discontinuousEdges.count(edge) ? 1 : 0;
			auto arcHistory = runtime.arcVersionObservationHistory.find(
				{edge, historyIt->second.arcVersion});
			if (arcHistory != runtime.arcVersionObservationHistory.end())
				audit.observations = arcHistory->second.observationEpochs;
			auto elevation = availability.edgeQuality.find(edge);
			if (elevation != availability.edgeQuality.end())
				audit.elevationScore = elevation->second;
			const auto residual = zhangIntegerSupportResidualSummary(
				&kfState, edge);
			audit.observations = std::min(
				audit.observations,
				std::min(residual.phaseSamples, residual.codeSamples));
			audit.phaseResidualRms = residual.phaseRms;
			audit.codeResidualRms = residual.codeRms;
			audit.phaseResidualMad = residual.phaseMad;
			audit.codeResidualMad = residual.codeMad;
			audit.whitenedResidualScore =
				residual.maximumWhitenedResidualScore;
			integerSupportQuality[edge] = std::move(audit);
        }
    }

    // The authoritative FLOAT datum graph must continue to retain every
    // usable observation.  Risk/quality selection below chooses only an
    // integer product target functional chart, never a second KF state chart.
    ZhangGraphBasis candidate = zhangBuildSpanningTree(
		stateEdges,
		options.reference_receiver,
		runtime.basis.treeEdges,
		activeQuality,
		options.prefer_historical_edges ? modelledEdges : set<ZhangGraphEdge>{},
		options.prefer_historical_edges ? persistence : map<ZhangGraphEdge, int>{}
	);
    set<ZhangGraphEdge> transportBrokenTree;
    for (const auto& edge : runtime.basis.treeEdges)
        if (availability.discontinuousEdges.contains(edge) || !availability.edges.contains(edge))
            transportBrokenTree.insert(edge);
    if (singleS && runtime.initialized && !transportBrokenTree.empty())
        candidate = certifiedTransportTreePlan(trace, kfState, sys, runtime,
            stateEdges, transportBrokenTree, availability.discontinuousEdges,
            options.reference_receiver, activeQuality,
            options.baseline_observables, previousSGaugeSignalArcs, availability.signalEdges).basis;

    if (zhangR51SingleBasisHandoffEnabled() && candidate.connected)
    {
        // Keep a currently admitted single-signal chord only when the FLOAT
        // chart already has that signal's ambiguity coordinate. A missing
        // L2 observation cannot erase a healthy L1 path, and an unmodelled
        // edge cannot be treated as an available integer relation.
        for (const auto& [signal, edges] : availability.signalEdges)
        for (const auto& edge : edges)
        {
            if (availability.edges.contains(edge) ||
                !runtime.basis.edges.contains(edge) ||
                !candidate.receivers.contains(edge.receiver) ||
                !candidate.satellites.contains(edge.satellite) ||
                candidate.treeEdges.contains(edge)) continue;
            const auto code = static_cast<E_ObsCode>(signal);
            const bool oldEndpointsRepresented =
                (edge.receiver == runtime.basis.rootReceiver ||
                 kfState.kfIndexMap.contains(
                     zhangTransportReceiverKey(
                         sys, code, edge.receiver))) &&
                kfState.kfIndexMap.contains(
                    zhangTransportSatelliteKey(code, edge.satellite));
            const bool oldTreePhysicalRow =
                oldEndpointsRepresented &&
                runtime.basis.treeEdges.contains(edge);
            const bool oldChordPhysicalRow =
                oldEndpointsRepresented &&
                !runtime.basis.treeEdges.contains(edge) &&
                kfState.kfIndexMap.contains(
                    zhangTransportAmbiguityKey(code, edge));
            if (!oldTreePhysicalRow && !oldChordPhysicalRow) continue;
            candidate.edges.insert(edge);
            stateEdges.insert(edge);
        }
    }

    // A tree transaction is triggered by every physical edge invalidation,
    // not only a confirmed dual-frequency slip.  The old implementation used
    // treeEdges intersect discontinuousEdges, which made station QC, signal
    // loss, receiver exit, and satellite set structurally unreachable in the
    // shadow/authoritative transaction path.
    set<ZhangGraphEdge> brokenTreeEdges;
    for (const auto& edge : runtime.basis.treeEdges)
    {
        if (availability.discontinuousEdges.contains(edge) ||
            !availability.edges.contains(edge))
        {
            brokenTreeEdges.insert(edge);
        }
    }
    if (zhangR51SingleBasisHandoffEnabled())
    {
        for (const auto& edge : availability.discontinuousEdges)
        {
            const auto history = runtime.edgeHistory.find(edge);
            const char* eventClass = runtime.basis.treeEdges.contains(edge)
                ? "TREE_PHYSICAL_ARC_BREAK"
                : runtime.basis.edges.contains(edge)
                ? "NON_TREE_PHYSICAL_ARC_BREAK"
                : "UNREPRESENTED_PHYSICAL_ARC_BREAK";
            trace << "\nZHANG_SINGLE_BASIS_PHYSICAL_EVENT time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " edge=" << edge.receiver << ":"
                  << edge.satellite.id()
                  << " class=" << eventClass
                  << " old_arc_version="
                  << (history == runtime.edgeHistory.end()
                      ? -1 : history->second.arcVersion)
                  << " current_epoch_integer_bridge=NOT_EVALUATED";
        }
    }

    auto updateProductDatum = [&](
        const string& reason,
        bool preserveIntegerComponent = false)
    {
        // The single-S mode has no autonomous P optimizer.  The read-only
        // compatibility view is synchronized AFTER the actual S commit, not
        // to a tentative candidate before a transform can still fail.
        if (singleS) return true;
        if (!candidate.connected || !candidate.receivers.contains(candidate.rootReceiver))
            return false;
        const ZhangGraphBasis oldProduct = runtime.productBasis;
        const map<ZhangGraphEdge, int> oldProductArcVersions =
            runtime.productArcVersions;
        set<ZhangGraphEdge> productEdges =
            zhangR51SingleBasisHandoffEnabled()
                ? availability.edges : candidate.edges;
        set<string> nextProductCoreReceivers;
        int productCoreMinimumSatelliteSupport = 0;
        if (options.product_core_min_satellite_support > 0)
        {
            ZhangProductReceiverCore core;
			if (options.product_integer_support_core)
			{
				const auto auditedCore = zhangBuildIntegerSupportCore(
					productEdges, candidate.rootReceiver,
					runtime.productCoreReceivers,
					options.product_core_min_satellite_support,
					integerSupportQuality, {}, activeQuality, persistence);
				BOOST_LOG_TRIVIAL(info)
					<< "ZHANG_INTEGER_SUPPORT_CORE sys=" << enum_to_string(sys)
					<< " qualified_edges=" << auditedCore.qualifiedEdges.size()
					<< " rejected_edges=" << auditedCore.rejectedEdges.size()
					<< " valid=" << auditedCore.valid
					<< " status=" << auditedCore.failureReason
					<< " float_graph_edges=" << candidate.edges.size()
					<< " current_dual_product_edges=" << productEdges.size();
				if (!auditedCore.valid)
				{
					BOOST_LOG_TRIVIAL(error)
						<< "ZHANG_PRODUCT_CORE sys=" << enum_to_string(sys)
						<< " status=REJECTED reason=" << auditedCore.failureReason;
					return false;
				}
				core = auditedCore.receiverCore;
			}
			else
			{
				core = zhangBuildProductReceiverCore(
					productEdges, candidate.rootReceiver,
					runtime.productCoreReceivers,
					options.product_core_min_satellite_support,
					activeQuality, persistence);
			}
            if (!core.connected)
            {
                BOOST_LOG_TRIVIAL(error)
                    << "ZHANG_PRODUCT_CORE sys=" << enum_to_string(sys)
                    << " status=REJECTED reason=CONNECTED_CORE_UNAVAILABLE";
                return false;
            }
            productEdges = core.edges;
            nextProductCoreReceivers = core.receivers;
            productCoreMinimumSatelliteSupport =
                core.minimumSatelliteSupport;
        }
        // A soft receiver-core choice must not hide a physically healthy
        // old tree before the datum proposal can compare it. Add only edges
        // passing the unchanged current integer-support and arc/QC gates.
        bool oldCoreTreeHealthy=runtime.productInitialized && oldProduct.connected &&
            oldProduct.rootReceiver==candidate.rootReceiver;
        for(const auto& edge:oldProduct.treeEdges) {
            const auto quality=integerSupportQuality.find(edge);
            const auto history=runtime.edgeHistory.find(edge);
            const auto previous=oldProductArcVersions.find(edge);
            oldCoreTreeHealthy &= candidate.edges.count(edge)>0 &&
                (!zhangR51SingleBasisHandoffEnabled() ||
                    availability.edges.contains(edge)) &&
                !availability.discontinuousEdges.count(edge) &&
                !availability.qcExcludedEdges.count(edge) &&
                !availability.elevationExcludedEdges.count(edge) &&
                history!=runtime.edgeHistory.end() && previous!=oldProductArcVersions.end();
            if(history!=runtime.edgeHistory.end() && previous!=oldProductArcVersions.end())
                oldCoreTreeHealthy &= history->second.arcVersion==previous->second;
            if(options.product_integer_support_core)
                oldCoreTreeHealthy &= quality!=integerSupportQuality.end() &&
                    zhangEvaluateIntegerSupportQuality(quality->second,{}).eligibleForIntegerSupport;
        }
        if(oldCoreTreeHealthy) {
            productEdges.insert(oldProduct.treeEdges.begin(),oldProduct.treeEdges.end());
            nextProductCoreReceivers.insert(oldProduct.receivers.begin(),oldProduct.receivers.end());
        }
        ZhangGraphBasis nextProduct = runtime.productInitialized
            ? (options.product_integer_support_core
				? zhangBuildRiskAwareSpanningTree(
					productEdges, candidate.rootReceiver,
					runtime.productBasis.treeEdges, integerSupportQuality)
				: zhangBuildSpanningTree(
                productEdges,
                candidate.rootReceiver,
                runtime.productBasis.treeEdges,
                activeQuality,
                modelledEdges,
                persistence))
            : (options.product_integer_support_core
				? zhangBuildRiskAwareSpanningTree(
					productEdges, candidate.rootReceiver, {}, integerSupportQuality)
				: zhangBuildRootedProductTree(
                productEdges,
                candidate.rootReceiver,
                {},
                activeQuality,
                modelledEdges,
                persistence));
        if (!nextProduct.connected)
        {
            return false;
        }
        if (zhangR51SingleBasisHandoffEnabled())
        {
            // A stochastic grace edge cannot define the current product
            // chart.  A broken arc with a fresh, QC-passing observation can:
            // it is a new physical arc, not a continuation of its old integer.
            // The hard-invalid/version checks below prevent old integer or
            // public-segment authority from crossing that break.
            for (const auto& edge : nextProduct.edges)
            {
                if (!zhangCurrentObservedProductEdge(edge,
                    availability.edges, availability.qcExcludedEdges,
                    availability.elevationExcludedEdges))
                {
                    trace << "\nZHANG_SINGLE_BASIS_PRODUCT_REJECT time="
                          << kfState.time.to_string(0)
                          << " edge=" << edge.receiver << ":"
                          << edge.satellite.id()
                          << " reason=NONCURRENT_PRODUCT_SUPPORT";
                    return false;
                }
            }
        }

        std::set<ZhangGraphEdge> hardInvalid=availability.discontinuousEdges;
        hardInvalid.insert(availability.qcExcludedEdges.begin(),availability.qcExcludedEdges.end());
        hardInvalid.insert(availability.elevationExcludedEdges.begin(),availability.elevationExcludedEdges.end());
        for(const auto& edge:oldProduct.treeEdges) {
            auto oldVersion=oldProductArcVersions.find(edge);
            auto now=runtime.edgeHistory.find(edge);
            if(oldVersion==oldProductArcVersions.end() || now==runtime.edgeHistory.end() ||
               oldVersion->second!=now->second.arcVersion) hardInvalid.insert(edge);
        }
        auto proposal=proposeProductDatum(oldProduct,nextProduct,hardInvalid,runtime.productInitialized);
        auditProductDatumTransport(proposal,oldProduct);
        BOOST_LOG_TRIVIAL(info)<<"R48_TREE_TRANSPORT time="<<kfState.time.to_string(0)
            <<" reason="<<reason<<" old_healthy="<<proposal.oldHealthy
            <<" proposed_changed="<<(oldProduct.treeEdges!=nextProduct.treeEdges)
            <<" hard_invalid_edges="<<hardInvalid.size()
            <<" integer_transport_proven="<<proposal.transportProven()
            <<" numeric_frontend_consistency_valid="<<proposal.transportProven()
            <<" status="<<proposal.status;
        // First mutation remains below, after the candidate has been audited.
        // The product tree selects physical integer functionals only.  It is
        // allowed to differ from the single committed estimation S-basis:
        // forcing the trees equal makes the exact S-to-product target G zero
        // and destroys the product-closure search lattice.  No KF coordinate
        // transform or integer authority is inherited by this selection.
        nextProduct = commitProductDatum(proposal);
        map<ZhangGraphEdge, int> nextProductArcVersions;
        map<ZhangGraphEdge, map<int, int>> nextProductSignalArcVersions;
        for (const auto& edge : nextProduct.treeEdges)
        {
            auto history = runtime.edgeHistory.find(edge);
            if (history != runtime.edgeHistory.end())
            {
                nextProductArcVersions[edge] = history->second.arcVersion;
            }
            if (zhangR51SingleBasisHandoffEnabled())
            {
                const auto physical = runtime.phaseArcVersions.find(edge);
                if (physical == runtime.phaseArcVersions.end())
                    return false;
                for (const auto code : options.baseline_observables)
                {
                    const auto signal = physical->second.find(int(code));
                    if (signal == physical->second.end() ||
                        signal->second < 0)
                        return false;
                    nextProductSignalArcVersions[edge][int(code)] =
                        signal->second;
                }
            }
        }

        set<ZhangGraphEdge> removedProductEdges;
        set<ZhangGraphEdge> addedProductEdges;
        set<ZhangGraphEdge> versionChangedEdges;
        if (runtime.productInitialized)
        {
            std::set_difference(
                oldProduct.treeEdges.begin(), oldProduct.treeEdges.end(),
                nextProduct.treeEdges.begin(), nextProduct.treeEdges.end(),
                std::inserter(removedProductEdges, removedProductEdges.begin())
            );
            std::set_difference(
                nextProduct.treeEdges.begin(), nextProduct.treeEdges.end(),
                oldProduct.treeEdges.begin(), oldProduct.treeEdges.end(),
                std::inserter(addedProductEdges, addedProductEdges.begin())
            );
            for (const auto& edge : nextProduct.treeEdges)
            {
                auto oldVersion = oldProductArcVersions.find(edge);
                auto newVersion = nextProductArcVersions.find(edge);
                if (oldProduct.treeEdges.find(edge) !=
                        oldProduct.treeEdges.end() &&
                    oldVersion != oldProductArcVersions.end() &&
                    newVersion != nextProductArcVersions.end() &&
                    oldVersion->second != newVersion->second)
                {
                    versionChangedEdges.insert(edge);
                }
            }
        }

        bool preserved = runtime.productInitialized;
        if (preserved)
        {
            set<string> commonReceivers;
            set<SatSys> commonSatellites;
            std::set_intersection(
                runtime.productBasis.receivers.begin(),
                runtime.productBasis.receivers.end(),
                nextProduct.receivers.begin(),
                nextProduct.receivers.end(),
                std::inserter(commonReceivers, commonReceivers.end())
            );
            std::set_intersection(
                runtime.productBasis.satellites.begin(),
                runtime.productBasis.satellites.end(),
                nextProduct.satellites.begin(),
                nextProduct.satellites.end(),
                std::inserter(commonSatellites, commonSatellites.end())
            );
            set<ZhangGraphEdge> survivingCommonGraph;
            for (const auto& edge : nextProduct.edges)
            {
                if (commonReceivers.find(edge.receiver) !=
                        commonReceivers.end() &&
                    commonSatellites.find(edge.satellite) !=
                        commonSatellites.end() &&
                    availability.discontinuousEdges.find(edge) ==
                        availability.discontinuousEdges.end())
                {
                    survivingCommonGraph.insert(edge);
                }
            }
            const auto commonStateBasis = zhangBuildSpanningTree(
                survivingCommonGraph, candidate.rootReceiver);
            set<ZhangGraphEdge> slippedOldProductEdges;
            std::set_intersection(
                oldProduct.treeEdges.begin(), oldProduct.treeEdges.end(),
                availability.discontinuousEdges.begin(),
                availability.discontinuousEdges.end(),
                std::inserter(slippedOldProductEdges,
                              slippedOldProductEdges.begin()));
            // Product continuity is a functional property, not old-tree
            // containment.  An exact same-arc S-basis re-expression preserves
            // the common product lattice whenever the surviving represented
            // state graph spans the same common nodes.  A physical arc-version
            // change remains float-only until an integer bridge is certified.
            preserved =
                runtime.productBasis.rootReceiver == candidate.rootReceiver &&
                commonStateBasis.connected &&
                commonStateBasis.receivers == commonReceivers &&
                commonStateBasis.satellites == commonSatellites &&
                versionChangedEdges.empty() &&
                slippedOldProductEdges.empty();
        }
        const bool graphContinuity=preserved;
        preserved=preserved && proposal.transportProven();
        BOOST_LOG_TRIVIAL(info)<<"R48_DATUM_CONTINUITY time="<<kfState.time.to_string(0)
            <<" graph_continuity_valid="<<graphContinuity
            <<" arc_segment_continuity_valid="<<hardInvalid.empty()
            <<" integer_transport_proven="<<proposal.transportProven()
            <<" preserve_integer_component_requested="<<preserveIntegerComponent;
        if (runtime.productInitialized && !preserved)
        {
            runtime.productDatumVersion++;
            runtime.integerComponentVersion++;
        }
        bool changed = !runtime.productInitialized ||
            runtime.productBasis.treeEdges != nextProduct.treeEdges ||
            runtime.productBasis.receivers != nextProduct.receivers ||
            runtime.productBasis.satellites != nextProduct.satellites ||
            !versionChangedEdges.empty();

        ZhangSatelliteSupportMetrics oldSatelliteMetrics;
        ZhangSatelliteSupportMetrics newSatelliteMetrics;
        if (changed && acsConfig.zhangPppAr.output_diagnostics)
        {
            oldSatelliteMetrics = zhangSatelliteSupportMetrics(oldProduct.edges);
            newSatelliteMetrics = zhangSatelliteSupportMetrics(nextProduct.edges);
        }

        struct ProductEdgeDiagnostic
        {
            std::optional<ZhangGraphEdge> oldEdge;
            std::optional<ZhangGraphEdge> newEdge;
            string                        eventReason;
            string                        signal = "ALL_BASELINE";
            int                           oldArcVersion = -1;
            int                           newArcVersion = -1;
            int                           oldAlternativePaths = 0;
            int                           newAlternativePaths = 0;
        };
        vector<ZhangGraphEdge> oldChanged(
            removedProductEdges.begin(), removedProductEdges.end()
        );
        oldChanged.insert(
            oldChanged.end(),
            versionChangedEdges.begin(), versionChangedEdges.end()
        );
        vector<ZhangGraphEdge> newChanged(
            addedProductEdges.begin(), addedProductEdges.end()
        );
        newChanged.insert(
            newChanged.end(),
            versionChangedEdges.begin(), versionChangedEdges.end()
        );

        auto rawContainsReceiver = [&](const string& receiver)
        {
            return std::any_of(
                availability.rawEdges.begin(), availability.rawEdges.end(),
                [&](const auto& edge) { return edge.receiver == receiver; }
            );
        };
        auto rawContainsSatellite = [&](const SatSys& satellite)
        {
            return std::any_of(
                availability.rawEdges.begin(), availability.rawEdges.end(),
                [&](const auto& edge) { return edge.satellite == satellite; }
            );
        };
        auto classifyOldEdge = [&](const ZhangGraphEdge& edge,
                                   int oldAlternativePaths)
        {
            if(availability.rejectedPhaseEdges.count(edge))
                return string("FILTER_PHASE_REJECTION_NOT_PROVEN_SLIP");
            if(availability.continuityUnknownEdges.count(edge))
                return string("PHASE_CONTINUITY_UNKNOWN_CONSERVATIVE_RETIREMENT");
            if (versionChangedEdges.find(edge) != versionChangedEdges.end() ||
                availability.discontinuousEdges.find(edge) !=
                    availability.discontinuousEdges.end())
            {
                return string("CONFIRMED_CYCLE_SLIP");
            }
            if (availability.qcExcludedEdges.find(edge) !=
                    availability.qcExcludedEdges.end() ||
                availability.elevationExcludedEdges.find(edge) !=
                    availability.elevationExcludedEdges.end())
            {
                return string("STATION_QC_REMOVAL");
            }
            if(availability.phaseQcEdges.count(edge))
                return string("PHASE_QC_QUARANTINE_NOT_RAW_OBSERVATION_LOSS");
            if (availability.signalUnavailableEdges.find(edge) !=
                    availability.signalUnavailableEdges.end())
            {
                return string("TEMPORARY_OBSERVATION_LOSS");
            }
            if (!rawContainsSatellite(edge.satellite))
            {
                return string("SATELLITE_RISE_SET");
            }
            if (!rawContainsReceiver(edge.receiver) ||
                availability.rawEdges.find(edge) == availability.rawEdges.end())
            {
                return string("TEMPORARY_OBSERVATION_LOSS");
            }
            if (oldAlternativePaths == 0)
            {
                return string("PRODUCT_EDGE_NO_ALTERNATIVE_SUPPORT");
            }
            if (nextProduct.satellites.size() < oldProduct.satellites.size() ||
                nextProduct.receivers.size() < oldProduct.receivers.size())
            {
                return string("COMPONENT_SPLIT");
            }
            return string("TREE_REOPTIMIZATION");
        };

        vector<ProductEdgeDiagnostic> edgeDiagnostics;
        const size_t diagnosticCount = std::max(
            oldChanged.size(), newChanged.size()
        );
        set<string> classifiedReasons;
        for (size_t index = 0; index < diagnosticCount; index++)
        {
            ProductEdgeDiagnostic diagnostic;
            if (index < oldChanged.size())
            {
                diagnostic.oldEdge = oldChanged[index];
                auto version = oldProductArcVersions.find(*diagnostic.oldEdge);
                if (version != oldProductArcVersions.end())
                {
                    diagnostic.oldArcVersion = version->second;
                }
                diagnostic.oldAlternativePaths =
                    zhangAlternativePhysicalPathCount(
                        oldProduct.edges, *diagnostic.oldEdge
                    );
                diagnostic.eventReason = classifyOldEdge(
                    *diagnostic.oldEdge,
                    diagnostic.oldAlternativePaths
                );
                auto signals = availability.discontinuitySignals.find(
                    *diagnostic.oldEdge
                );
                if (signals != availability.discontinuitySignals.end())
                {
                    std::ostringstream listed;
                    for (E_ObsCode code : signals->second)
                    {
                        listed << (listed.tellp() > 0 ? "," : "")
                               << enum_to_string(code);
                    }
                    diagnostic.signal = listed.str();
                }
            }
            if (index < newChanged.size())
            {
                diagnostic.newEdge = newChanged[index];
                auto version = nextProductArcVersions.find(*diagnostic.newEdge);
                if (version != nextProductArcVersions.end())
                {
                    diagnostic.newArcVersion = version->second;
                }
                diagnostic.newAlternativePaths =
                    zhangAlternativePhysicalPathCount(
                        nextProduct.edges, *diagnostic.newEdge
                    );
            }
            if (diagnostic.eventReason.empty())
            {
                diagnostic.eventReason =
                    nextProduct.satellites.size() > oldProduct.satellites.size() ||
                    nextProduct.receivers.size() > oldProduct.receivers.size()
                        ? "COMPONENT_MERGE"
                        : "TREE_REOPTIMIZATION";
            }
            classifiedReasons.insert(diagnostic.eventReason);
            edgeDiagnostics.push_back(std::move(diagnostic));
        }

        struct TreeSupportSummary
        {
            int minimum = 0;
            int maximum = 0;
            int bridgeCount = 0;
            double mean = 0;
        };
        auto treeSupportSummary = [&](const ZhangGraphBasis& basis)
        {
            TreeSupportSummary summary;
            if (basis.treeEdges.empty())
            {
                return summary;
            }
            summary.minimum = std::numeric_limits<int>::max();
            long long total = 0;
            for (const auto& edge : basis.treeEdges)
            {
                int support = 1 + zhangAlternativePhysicalPathCount(
                    basis.edges, edge
                );
                summary.minimum = std::min(summary.minimum, support);
                summary.maximum = std::max(summary.maximum, support);
                summary.bridgeCount += support == 1;
                total += support;
            }
            summary.mean = static_cast<double>(total) / basis.treeEdges.size();
            return summary;
        };
        TreeSupportSummary oldTreeSupport;
        TreeSupportSummary newTreeSupport;
        if (changed && acsConfig.zhangPppAr.output_diagnostics)
        {
            oldTreeSupport = treeSupportSummary(oldProduct);
            newTreeSupport = treeSupportSummary(nextProduct);
        }

        // A change of the read-only P target chart is an event even when the
        // sole KF S chart did not move.  Capture the old *actual* product row
        // before committing the new target, then let the ordinary writer-row
        // binder compare the new row with it.  The identity coordinate plan
        // proves full joint x/P preservation without claiming that different
        // P functionals themselves have equal means or covariances.  A broken
        // physical arc has no identity-plan auxiliary and therefore remains
        // fail-closed until a correlated event transform supplies one.
        if (zhangR51SingleBasisHandoffEnabled() && changed &&
            runtime.productInitialized)
        {
            const string runtimeId = zhangGraphRuntimeId(kfState);
            const string epoch = kfState.time.to_string(9);
            bool priorGraphReceipt = false;
            for (const auto& [identity, receipt] :
                 graphProductEventReceipts)
                if (std::get<0>(identity) == runtimeId &&
                    std::get<1>(identity) == epoch &&
                    std::get<2>(identity) == sys)
                {
                    priorGraphReceipt = true;
                    break;
                }
            if (!priorGraphReceipt)
            {
                map<ZhangGraphEdge, int> currentVersions;
                for (const auto& [edge, history] : runtime.edgeHistory)
                    currentVersions[edge] = history.arcVersion;
                const auto identityPrecondition =
                    zhangCheckProductChartIdentityReceiptPreconditions(
                        runtime.basis, candidate, oldProduct, nextProduct,
                        oldProductArcVersions, currentVersions,
                        modelledEdges, availability.discontinuousEdges);
                if (!identityPrecondition.valid)
                {
                    trace << "\nZHANG_PRODUCT_CHART_ONLY_RECEIPT time="
                          << epoch << " system=" << enum_to_string(sys)
                          << " status=REJECTED reason="
                          << identityPrecondition.reason
                          << " integer_bridge_proven=0";
                }
                else
                {
                    ZhangGraphBasis receiptBasis = runtime.basis;
                    receiptBasis.edges = modelledEdges;
                    ZhangGraphCoordinatePlan identity;
                    identity.valid = true;
                    for (const auto& [key, index] : kfState.kfIndexMap)
                        identity.transform[key][key] = 1;
                    auto receipts = prepareGraphProductEventReceipts(
                        kfState, sys, options.baseline_observables,
                        receiptBasis, oldProduct, receiptBasis,
                        runtime.phaseArcVersions,
                        runtime.phaseArcVersions,
                        availability.discontinuousEdges, identity,
                        "PRODUCT_CHART_REOPTIMISATION");
                    if (!receipts.empty())
                    {
                        validateGraphProductEventJointMoments(
                            trace, kfState, receipts);
                        storeGraphProductEventReceipts(
                            kfState, sys, std::move(receipts));
                        trace << "\nZHANG_PRODUCT_CHART_ONLY_RECEIPT time="
                              << epoch << " system=" << enum_to_string(sys)
                              << " old_tree_edges=" << oldProduct.treeEdges.size()
                              << " new_tree_edges=" << nextProduct.treeEdges.size()
                              << " s_tree_unchanged=1 full_joint_identity_plan=1"
                              << " integer_bridge_proven=0";
                    }
                }
            }
        }

        runtime.productBasis = std::move(nextProduct);
        runtime.productArcVersions = std::move(nextProductArcVersions);
        runtime.productSignalArcVersions =
            std::move(nextProductSignalArcVersions);
        runtime.productCoreReceivers = std::move(nextProductCoreReceivers);
        runtime.productInitialized = true;
        runtime.productAuthorityComponentVersion = runtime.integerComponentVersion;
        runtime.requiresFreshProductIdentity = false;
        if (zhangR51SingleBasisHandoffEnabled())
            bindGraphProductEventWriterRows(
                trace, kfState, sys, options.baseline_observables,
                candidate, runtime);
        if (changed)
        {
            std::ostringstream eventCause;
            for (const auto& classified : classifiedReasons)
            {
                eventCause << (eventCause.tellp() > 0 ? "," : "")
                           << classified;
            }
            runtime.lastProductEventCause = classifiedReasons.empty()
                ? "INITIALISE"
                : eventCause.str();
        }
        if (changed && acsConfig.zhangPppAr.output_diagnostics)
        {
            std::ostringstream reasons;
            for (const auto& classified : classifiedReasons)
            {
                reasons << (reasons.tellp() > 0 ? "," : "") << classified;
            }
            trace << "\nZHANG_PRODUCT_DATUM_EVENT time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " reason=" << reason
                  << " classified_reasons="
                  << (classifiedReasons.empty() ? "INITIALISE" : reasons.str())
                  << " datum_version=" << runtime.productDatumVersion
                  << " continuity_preserved=" << preserved
                  << " component_id="
                  << zhangIntegerComponentId(runtime.productBasis)
                  << " product_tree_edges="
                  << runtime.productBasis.treeEdges.size()
                  << " target_tree_equals_s_tree="
                  << (runtime.productBasis.treeEdges == candidate.treeEdges)
                  << " target_tree_current_chords="
                  << std::count_if(
                         runtime.productBasis.treeEdges.begin(),
                         runtime.productBasis.treeEdges.end(),
                         [&](const auto& edge)
                         { return !candidate.treeEdges.contains(edge); })
                  << " product_receivers="
                  << runtime.productBasis.receivers.size()
                  << " product_satellites="
                  << runtime.productBasis.satellites.size()
                  << " product_core_enabled="
                  << (options.product_core_min_satellite_support > 0)
                  << " product_core_receivers="
                  << runtime.productCoreReceivers.size()
                  << " product_core_min_satellite_support="
                  << productCoreMinimumSatelliteSupport
                  << " old_product_tree_edges=" << edgeList(removedProductEdges)
                  << " new_product_tree_edges=" << edgeList(addedProductEdges)
                  << " arc_version_changes=" << edgeList(versionChangedEdges)
                  << " old_tree_min_support=" << oldTreeSupport.minimum
                  << " old_tree_mean_support=" << oldTreeSupport.mean
                  << " old_tree_bridge_count=" << oldTreeSupport.bridgeCount
                  << " new_tree_min_support=" << newTreeSupport.minimum
                  << " new_tree_mean_support=" << newTreeSupport.mean
                  << " new_tree_bridge_count=" << newTreeSupport.bridgeCount
                  << " old_satellite_bridge_count="
                  << oldSatelliteMetrics.bridgeEdges.size()
                  << " new_satellite_bridge_count="
                  << newSatelliteMetrics.bridgeEdges.size()
                  << " old_satellite_edge_connectivity="
                  << oldSatelliteMetrics.edgeConnectivity
                  << " new_satellite_edge_connectivity="
                  << newSatelliteMetrics.edgeConnectivity;

            for (const auto& diagnostic : edgeDiagnostics)
            {
                auto edgeText = [](const std::optional<ZhangGraphEdge>& edge)
                {
                    return edge
                        ? edge->receiver + ":" + edge->satellite.id()
                        : string("NONE");
                };
                trace << "\nZHANG_PRODUCT_DATUM_EDGE_EVENT time="
                      << kfState.time.to_string(0)
                      << " system=" << enum_to_string(sys)
                      << " old_product_tree_edge="
                      << edgeText(diagnostic.oldEdge)
                      << " new_product_tree_edge="
                      << edgeText(diagnostic.newEdge)
                      << " event_reason=" << diagnostic.eventReason
                      << " receiver="
                      << (diagnostic.oldEdge
                              ? diagnostic.oldEdge->receiver
                              : diagnostic.newEdge
                                    ? diagnostic.newEdge->receiver
                                    : "NONE")
                      << " satellite="
                      << (diagnostic.oldEdge
                              ? diagnostic.oldEdge->satellite.id()
                              : diagnostic.newEdge
                                    ? diagnostic.newEdge->satellite.id()
                                    : "NONE")
                      << " signal=" << diagnostic.signal
                      << " old_arc_version=" << diagnostic.oldArcVersion
                      << " new_arc_version=" << diagnostic.newArcVersion
                      << " old_support_count="
                      << (diagnostic.oldEdge
                              ? 1 + diagnostic.oldAlternativePaths : 0)
                      << " new_support_count="
                      << (diagnostic.newEdge
                              ? 1 + diagnostic.newAlternativePaths : 0)
                      << " old_alternative_exact_paths="
                      << diagnostic.oldAlternativePaths
                      << " new_alternative_exact_paths="
                      << diagnostic.newAlternativePaths
                      << " bridge_before="
                      << (diagnostic.oldEdge &&
                          diagnostic.oldAlternativePaths == 0)
                      << " bridge_after="
                      << (diagnostic.newEdge &&
                          diagnostic.newAlternativePaths == 0)
                      << " component_size_before="
                      << oldSatelliteMetrics.largestComponent
                      << " component_size_after="
                      << newSatelliteMetrics.largestComponent
                  << " datum_version_changed=" << !preserved;
            }

            trace << "\nZHANG_PRODUCT_GRAPH_REDUNDANCY time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " product_satellites="
                  << newSatelliteMetrics.satellites.size()
                  << " product_relation_edges="
                  << newSatelliteMetrics.supportCounts.size()
                  << " mean_support_count="
                  << newSatelliteMetrics.meanSupport
                  << " min_support_count="
                  << newSatelliteMetrics.minimumSupport
                  << " max_support_count="
                  << newSatelliteMetrics.maximumSupport
                  << " bridge_count="
                  << newSatelliteMetrics.bridgeEdges.size()
                  << " edge_connectivity="
                  << newSatelliteMetrics.edgeConnectivity
                  << " component_count="
                  << newSatelliteMetrics.componentCount
                  << " largest_component="
                  << newSatelliteMetrics.largestComponent
                  << " product_tree_min_support="
                  << newTreeSupport.minimum
                  << " product_tree_mean_support="
                  << newTreeSupport.mean
                  << " product_tree_bridge_count="
                  << newTreeSupport.bridgeCount
                  << " datum_version=" << runtime.productDatumVersion;

            trace << "\nZHANG_PRODUCT_FUNCTIONAL_CONTINUITY time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " representation_changed=" << changed
                  << " same_arc_state_graph_connected=" << preserved
                  << " float_gauge_continuous=1"
                  << " integer_cross_component_valid=" << preserved
                  << " classification="
                  << (preserved
                        ? "EXACT_S_TRANSFORM_SAME_ARC"
                        : "FLOAT_ONLY_PENDING_INTEGER_TRANSPORT")
                  << " product_functional_suspend=" << !preserved;
        }
        return true;
    };

    auto auditRetiredProductArcs = [&](
        const ZhangGraphBasis& oldProductBasis,
        const map<ZhangGraphEdge, int>& oldProductArcVersions,
        const set<ZhangGraphEdge>& retiredEdges,
        const string&              reason)
    {
        if (oldProductBasis.treeEdges.empty() || !runtime.productInitialized)
        {
            return false;
        }

        set<ZhangGraphEdge> retiredProductTreeEdges;
        std::set_intersection(
            oldProductBasis.treeEdges.begin(),
            oldProductBasis.treeEdges.end(),
            retiredEdges.begin(),
            retiredEdges.end(),
            std::inserter(
                retiredProductTreeEdges,
                retiredProductTreeEdges.begin()));
        if (retiredProductTreeEdges.empty())
        {
            return false;
        }

        // updateProductDatum() has already rebuilt the current product chart
        // with current physical arc versions and audited the datum transport.
        // Endpoint overlap with an old tree edge is not, by itself, evidence
        // that every current product relation is invalid.  Keep the rebuilt
        // structural chart so exact physical-history projection and current
        // recertification can operate; formal writers still require a current
        // product certificate and therefore fail closed independently.
        if (acsConfig.zhangPppAr.output_diagnostics)
        {
            trace << "\nZHANG_INTEGER_CONTINUITY_EVENT time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " reason=" << reason
                  << " retired_product_tree_edges="
                  << edgeList(retiredProductTreeEdges)
                  << " retired_product_arc_versions=";
            for (const auto& edge : retiredProductTreeEdges)
            {
                const auto version = oldProductArcVersions.find(edge);
                trace << edge.receiver << ":" << edge.satellite.id()
                      << ":V" << (version == oldProductArcVersions.end()
                            ? -1 : version->second) << ";";
            }
            trace
                  << " integer_component_version="
                  << runtime.integerComponentVersion
                  << " float_gauge_version=" << runtime.floatGaugeVersion
                  << " surviving_product_edges="
                  << runtime.productBasis.edges.size()
                  << " surviving_product_tree_rank="
                  << runtime.productBasis.treeEdges.size()
                  << " surviving_component_count="
                  << runtime.productBasis.componentCount
                  << " component_split="
                  << (runtime.productBasis.componentCount > 1)
                  << " float_gauge_continuous=1"
                  << " structural_product_chart_preserved=1"
                  << " endpoint_only_global_invalidation=0"
                  << " ar_valid=DEFERRED_TO_CURRENT_PRODUCT_CERTIFICATE"
                  << " pending_component_gauge_bridge=1"
                  << " pending_besd_bridge=1"
                  << " float_posterior_preserved=1";
        }
        return true;
    };

    auto localReinitialise = [&](const string& reason)
    {
        const auto savedRuntime = runtime;
        const auto oldProductBasis = runtime.productBasis;
        const auto oldProductArcVersions = runtime.productArcVersions;
        int removedColumns = 0;
        const auto oldIndex = kfState.kfIndexMap;
        ZhangGraphBasis oldBasis = runtime.basis;
        oldBasis.edges = modelledEdges;
        auto eligible = stateCandidates;
        eligible.insert(availability.edges.begin(), availability.edges.end());
        if (zhangR51SingleBasisHandoffEnabled())
            eligible.insert(candidate.edges.begin(), candidate.edges.end());
        ZhangGraphBasis transported;
        ZhangGraphCoordinatePlan plan;
        if (!transportZhangGraphComponents(trace, kfState, sys,
                options.baseline_observables, oldBasis, eligible,
                availability.edges, availability.signalEdges,
                availability.discontinuousEdges, availability.discontinuitySignals, transported, plan))
        {
            runtime = savedRuntime;
            return false;
        }
        // The covariance transform and all independent fresh priors have one
        // accepted factor transaction. Only now consume physical arc events.
        advanceDiscontinuousArcVersions();
        auto retiredCoordinates = availability.discontinuousEdges;
        for (const auto& edge : plan.freshEdges)
            if (oldBasis.edges.contains(edge) && !retiredCoordinates.contains(edge))
            {
                // An expired/missing coordinate cannot inherit authorization
                // from an old integer simply because its endpoint names match.
                runtime.edgeHistory[edge].arcVersion++;
                for (const auto code : options.baseline_observables)
                    if (!plan.survivingPhysicalRows.contains({code,edge}))
                        runtime.phaseArcVersions[edge][int(code)]++;
                retiredCoordinates.insert(edge);
            }
        candidate = transported;
        stateEdges = transported.edges;
        observationEdges.clear();
        std::set_intersection(availability.edges.begin(), availability.edges.end(),
            transported.edges.begin(), transported.edges.end(),
            std::inserter(observationEdges, observationEdges.begin()));
        if (candidate.connected && candidate.receivers.contains(candidate.rootReceiver))
        {
            if (!updateProductDatum("component_reparameterise"))
                invalidateZhangProduct(runtime, "COMPONENT_PRODUCT_DATUM_UNAVAILABLE");
        }
        else
        {
            invalidateZhangProduct(runtime, "FLOAT_COMPONENT_SPLIT");
        }
        // FLOAT transport alone never certifies a new inter-component bridge.
        // Preserve the rebuilt structural chart and let exact physical
        // projection plus the formal product certificate decide which integer
        // relations, if any, remain publishable.
        auditRetiredProductArcs(oldProductBasis, oldProductArcVersions,
            retiredCoordinates,
            "PHYSICAL_ARC_RETIRED_DURING_COMPONENT_TRANSPORT");
        runtime.basis = transported;
        // The observation basis must never select a different chart from the
        // committed one: keep datum-only tree edges out of observationEdges.
        runtime.activeBasis = transported;
        runtime.observationEdges = observationEdges;
        runtime.stateEdges = stateEdges;
        runtime.deferredEpochs = 0;
        runtime.representationVersion++;
        for (const auto& [key, index] : oldIndex)
            removedColumns += key.type == KF::AMBIGUITY && !kfState.kfIndexMap.contains(key);
        traceGraphEvent("component_reparameterise", reason,
            brokenTreeEdges, plan.freshEdges, {}, {},
            removedColumns, false);
        trace << "\nZHANG_COMPONENT_TRANSPORT_COMMIT time=" << kfState.time.to_string(0)
              << " raw_phase_edges=" << availability.edges.size()
              << " raw_code_edges=" << availability.codeEdges.size()
              << " admitted_phase_edges=" << observationEdges.size()
              << " retained_physical_rows=" << plan.survivingPhysicalRows.size()
              << " retired_arc_coordinates=" << retiredCoordinates.size()
              << " fresh_edges=" << plan.freshEdges.size()
              << " node_rectangle_reset=0 integer_bridge_certified=0";
        return true;
    };

    if (allCurrentArcsRetiredWithFreshObservations &&
        runtime.initialized && hasEstimatedPhaseState)
    {
        // The pre-event carrier graph is empty, but the new physical arcs
        // have valid observations.  Rebuild them as destination-only fresh
        // priors in one state/factor transaction.  Holding here would leave
        // their old versions pending forever; directly reusing the old chart
        // would instead inherit retired integer coordinates.
        if (!localReinitialise("ALL_CURRENT_ARCS_RETIRED"))
        {
            runtime.observationEdges.clear();
            runtime.stateEdges.clear();
            invalidateZhangProduct(runtime,
                "ALL_CURRENT_ARCS_RETIREMENT_DEFERRED");
            trace << "\nZHANG_GRAPH_EMPTY_SUPPORT time="
                  << kfState.time.to_string(0)
                  << " pending_physical_events="
                  << availability.discontinuousEdges.size()
                  << " action=HOLD_POSTERIOR_NO_PHASE_ROWS"
                  << " product_authorized=0";
        }
        return;
    }

    bool partialSignalRetirement = false;
    for (const auto& [edge, codes] : availability.discontinuitySignals)
        partialSignalRetirement |= !codes.empty() && codes.size() < options.baseline_observables.size();
    if (runtime.initialized && hasEstimatedPhaseState && partialSignalRetirement)
    {
        if (!localReinitialise("SIGNAL_LEVEL_ARC_RETIREMENT"))
            retainOldTreeRootComponent();
        return;
    }
    if (!runtime.initialized || !hasEstimatedPhaseState)
    {
        advanceDiscontinuousArcVersions();
        if (!updateProductDatum("initialise"))
            invalidateZhangProduct(runtime, "INITIALISE_PRODUCT_UNAVAILABLE");
        runtime.basis            = candidate;
        runtime.activeBasis      = candidate;
        runtime.observationEdges = observationEdges;
        runtime.stateEdges       = stateEdges;
        runtime.initialized = true;
        runtime.deferredEpochs = 0;
        traceGraphEvent(
            "initialise",
            "initial_component",
            {},
            candidate.treeEdges,
            {},
            {},
            0,
            true
        );
        traceCanonicalAudit(runtime.basis, "initialise", true);

        BOOST_LOG_TRIVIAL(info)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=initialise"
            << " nodes=" << candidate.receivers.size() + candidate.satellites.size()
            << " edges=" << candidate.edges.size()
            << " tree_edges=" << candidate.treeEdges.size()
            << " cycles=" << candidate.edges.size() - candidate.treeEdges.size();
        return;
    }

    std::ostringstream preEventVersions;
    bool firstPreEventVersion = true;
    for (const auto& [edge, version] : preEventArcVersions)
    {
        preEventVersions
            << (firstPreEventVersion ? "" : ";")
            << edge.receiver << "|" << edge.satellite.id()
            << "|A" << version;
        firstPreEventVersion = false;
    }

    // A slipped chord does not require a tree pivot, but its old ambiguity
    // coordinate must still be retired atomically before the post-slip sample
    // is allowed to create the next arc version.
    if (brokenTreeEdges.empty() && !availability.discontinuousEdges.empty())
    {
        set<ZhangGraphEdge> retiredChordEdges;
        std::set_intersection(
            modelledEdges.begin(), modelledEdges.end(),
            availability.discontinuousEdges.begin(),
            availability.discontinuousEdges.end(),
            std::inserter(retiredChordEdges, retiredChordEdges.begin()));
        if (!retiredChordEdges.empty())
        {
            ZhangGraphBasis oldBasis = runtime.basis;
            oldBasis.edges = modelledEdges;
            oldBasis.receivers.clear();
            oldBasis.satellites.clear();
            for (const auto& edge : modelledEdges)
            {
                oldBasis.receivers.insert(edge.receiver);
                oldBasis.satellites.insert(edge.satellite);
            }
            ZhangGraphBasis retainedBasis = runtime.basis;
            retainedBasis.edges.clear();
            std::set_difference(
                modelledEdges.begin(), modelledEdges.end(),
                retiredChordEdges.begin(), retiredChordEdges.end(),
                std::inserter(retainedBasis.edges, retainedBasis.edges.begin()));
            auto postRetirementStateEdges = stateEdges;
            for (const auto& edge : availability.discontinuousEdges)
                if (availability.edges.contains(edge) &&
                    retainedBasis.receivers.contains(edge.receiver) &&
                    retainedBasis.satellites.contains(edge.satellite))
                    postRetirementStateEdges.insert(edge);
            const auto eligibilityParts = partitionZhangEligibility(
                retainedBasis, postRetirementStateEdges);
            const auto activeAfterRetirement = zhangActiveGraphBasisWithAugmentation(
                retainedBasis, postRetirementStateEdges);
            if (!activeAfterRetirement.valid)
            {
                trace << "\nR48_ACTIVE_BASIS_REJECT time=" << kfState.time.to_string(0)
                      << " phase=CHORD_RETIRE reason=" << activeAfterRetirement.failureReason;
                if (!localReinitialise("CHORD_RETIRE_CHART_AUGMENTATION"))
                    retainOldTreeRootComponent();
                return;
            }
            if (transformZhangGraphBasis(
                    trace,
                    kfState,
                    sys,
                    options.baseline_observables,
                    oldBasis,
                    retainedBasis,
                    true,
                    retiredChordEdges,
                    nullptr,
                    nullptr,
                    &availability.discontinuitySignals))
            {
                advanceDiscontinuousArcVersions();
                for (const auto& edge : availability.discontinuousEdges)
                {
                    if (availability.edges.find(edge) !=
                            availability.edges.end() &&
                        retainedBasis.receivers.find(edge.receiver) !=
                            retainedBasis.receivers.end() &&
                        retainedBasis.satellites.find(edge.satellite) !=
                            retainedBasis.satellites.end())
                    {
                        stateEdges.insert(edge);
                        observationEdges.insert(edge);
                    }
                }
                candidate = activeAfterRetirement.basis;
                trace << "\nR48_ACTIVE_BASIS_COMMIT time=" << kfState.time.to_string(0)
                      << " phase=CHORD_RETIRE retained_datum_only_edges="
                      << edgeList(activeAfterRetirement.retainedDatumOnlyEdges)
                      << " coordinate_tree_unchanged=1 observation_eligibility_unchanged=1";
                if (!updateProductDatum("chord_arc_retire"))
                {
                    invalidateZhangProduct(runtime, "PENDING_BRIDGE_AFTER_CHORD_ARC_RETIRE");
                }
                runtime.representationVersion++;
                runtime.basis = zhangRepresentedBasisAfterAugmentation(retainedBasis, candidate);
                trace << "\nZHANG_GRAPH_NODE_AUGMENTATION time="
                      << kfState.time.to_string(0) << " phase=CHORD_RETIRE"
                      << " pending_edges=" << edgeList(eligibilityParts.requiresAugmentation)
                      << " old_chart_transported_first=1 integer_certification=0";
                runtime.activeBasis = candidate;
                runtime.observationEdges = observationEdges;
                runtime.stateEdges = stateEdges;
                runtime.deferredEpochs = 0;
                traceGraphEvent(
                    "chord_arc_retire",
                    "non_tree_arc_break",
                    {},
                    {},
                    {},
                    {},
                    static_cast<int>(retiredChordEdges.size() *
                        options.baseline_observables.size()),
                    true);
                return;
            }
            BOOST_LOG_TRIVIAL(error)
                << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
                << " action=defer chord retirement: atomic projection failed";
            retainOldTreeRootComponent();
            return;
        }
        advanceDiscontinuousArcVersions();
    }
    else if (brokenTreeEdges.empty())
    {
        advanceDiscontinuousArcVersions();
    }

    if (brokenTreeEdges.empty() &&
        candidate.rootReceiver == runtime.basis.rootReceiver &&
        candidate.treeEdges == runtime.basis.treeEdges &&
        candidate.receivers == runtime.basis.receivers &&
        candidate.satellites == runtime.basis.satellites)
    {
        ZhangGraphBasis augmented;
        if (!zhangAugmentCommittedSChords(runtime.basis, candidate, modelledEdges, augmented))
            throw std::runtime_error("ZHANG_UNCHANGED_TREE_CHORD_AUGMENTATION_INVALID");
        // The missing commit here formerly discarded a new chord from
        // PhaseAdmissionGuard/zhangGraphRetainsAmbiguity.  Catalogue extension
        // lets the measurement builder initialize ONLY the fresh ambiguity;
        // no old state, full covariance or cross block is reinitialized.
        runtime.basis = std::move(augmented);
        if (!updateProductDatum("graph_update"))
            invalidateZhangProduct(runtime, "GRAPH_UPDATE_PRODUCT_UNAVAILABLE");
        runtime.activeBasis      = candidate;
        runtime.observationEdges = observationEdges;
        runtime.stateEdges       = stateEdges;
        runtime.deferredEpochs = 0;
        return;
    }

    // A new receiver or satellite may extend the existing tree before any state for that node
    // exists.  Such a leaf extension leaves every existing coordinate unchanged.
    bool oldTreeRetained = std::includes(
        candidate.treeEdges.begin(),
        candidate.treeEdges.end(),
        runtime.basis.treeEdges.begin(),
        runtime.basis.treeEdges.end()
    );
    bool leafExtension = brokenTreeEdges.empty() && oldTreeRetained &&
        candidate.rootReceiver == runtime.basis.rootReceiver;
    if (brokenTreeEdges.empty() && leafExtension)
    {
        for (const auto& edge : candidate.treeEdges)
        {
            if (runtime.basis.treeEdges.find(edge) != runtime.basis.treeEdges.end())
            {
                continue;
            }

            bool newReceiver =
                runtime.basis.receivers.find(edge.receiver) == runtime.basis.receivers.end();
            bool newSatellite =
                runtime.basis.satellites.find(edge.satellite) == runtime.basis.satellites.end();
            if (!newReceiver && !newSatellite)
            {
                leafExtension = false;
                break;
            }
        }
    }

    if (leafExtension)
    {
        set<ZhangGraphEdge> addedTreeEdges;
        std::set_difference(
            candidate.treeEdges.begin(),
            candidate.treeEdges.end(),
            runtime.basis.treeEdges.begin(),
            runtime.basis.treeEdges.end(),
            std::inserter(addedTreeEdges, addedTreeEdges.begin())
        );
        if (singleS && zhangR51SingleBasisHandoffEnabled() &&
            runtime.productInitialized)
        {
            const auto runtimeId = zhangGraphRuntimeId(kfState);
            const auto epoch = kfState.time.to_string(9);
            const bool alreadyCaptured = std::any_of(
                graphProductEventReceipts.begin(), graphProductEventReceipts.end(),
                [&](const auto& entry)
                {
                    return std::get<0>(entry.first) == runtimeId &&
                        std::get<1>(entry.first) == epoch &&
                        std::get<2>(entry.first) == sys;
                });
            // An earlier real graph transform may already have captured the
            // published source in this epoch.  A leaf catalogue extension
            // leaves its destination KF rows unchanged: keep that earliest
            // source and let the ordinary writer binder validate it again.
            // Do not replace it with a later source, or claim composition of
            // two arbitrary graph transforms (which remains fail-closed).
            if (!alreadyCaptured)
            {
                // The catalogue will acquire new nodes, but no existing KF
                // coordinate is transformed here. Capture the OLD actual
                // product rows and their entire phase/clock joint block before
                // committing either target catalogue. A new node has no old
                // product endpoint and is deliberately excluded by oldBasis.
                //
                // Publishing no receipt while incrementing the event counter
                // used to strand old products at the previous backend. This
                // receipt is only FLOAT evidence; the ordinary writer binder
                // validates it after the new catalogue commits. It does not
                // authorize an integer or replace any physical-arc check.
                ZhangGraphCoordinatePlan identity;
                identity.valid = true;
                for (const auto& [key, index] : kfState.kfIndexMap)
                    identity.transform[key][key] = 1;
                auto receipts = prepareGraphProductEventReceipts(
                    kfState, sys, options.baseline_observables,
                    runtime.basis, runtime.productBasis, candidate,
                    previousSGaugeSignalArcs, runtime.phaseArcVersions,
                    availability.discontinuousEdges, identity,
                    "LEAF_EXTENSION_EXISTING_PRODUCT_IDENTITY");
                if (!receipts.empty())
                {
                    validateGraphProductEventJointMoments(
                        trace, kfState, receipts);
                    storeGraphProductEventReceipts(
                        kfState, sys, std::move(receipts));
                    trace << "\nZHANG_PRODUCT_LEAF_IDENTITY_RECEIPT time="
                          << kfState.time.to_string(0)
                          << " system=" << enum_to_string(sys)
                          << " old_satellites=" << runtime.basis.satellites.size()
                          << " new_satellites=" << candidate.satellites.size()
                          << " old_state_dimension=" << kfState.x.size()
                          << " old_coordinates_modified=0 new_node_old_endpoints=0"
                             " actual_writer_bound=DEFERRED integer_bridge_proven=0";
                }
            }
            else
                trace << "\nZHANG_PRODUCT_LEAF_IDENTITY_RECEIPT time="
                      << epoch << " system=" << enum_to_string(sys)
                      << " status=EXISTING_SOURCE_PRESERVED"
                         " existing_destination_rows_modified=0"
                         " actual_writer_bound=DEFERRED integer_bridge_proven=0";
        }
        if (!updateProductDatum("leaf_extension"))
            invalidateZhangProduct(runtime, "LEAF_EXTENSION_PRODUCT_UNAVAILABLE");
        runtime.basis            = zhangRepresentedBasisAfterAugmentation(runtime.basis, candidate);
        runtime.activeBasis      = candidate;
        runtime.observationEdges = observationEdges;
        runtime.stateEdges       = stateEdges;
        runtime.deferredEpochs = 0;
        runtime.representationVersion++;
        traceGraphEvent(
            "leaf_extension",
            "new_leaf_node",
            {},
            addedTreeEdges,
            {},
            {},
            0,
            true
        );
        traceCanonicalAudit(runtime.basis, "leaf_extension", true);
        BOOST_LOG_TRIVIAL(info)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=leaf_extension"
            << " tree_edges=" << candidate.treeEdges.size();
        return;
    }

    auto detachedNodes = [&]()
    {
        set<ZhangGraphEdge> retainedOldTree;
        std::set_intersection(
            runtime.basis.treeEdges.begin(),
            runtime.basis.treeEdges.end(),
            stateEdges.begin(),
            stateEdges.end(),
            std::inserter(retainedOldTree, retainedOldTree.begin())
        );
        set<ZhangGraphEdge> rootTree =
            zhangRootComponentEdges(retainedOldTree, runtime.basis.rootReceiver);

        set<string> rootReceivers = {runtime.basis.rootReceiver};
        set<SatSys> rootSatellites;
        for (const auto& edge : rootTree)
        {
            rootReceivers.insert(edge.receiver);
            rootSatellites.insert(edge.satellite);
        }

        set<string> affectedReceivers;
        set<SatSys> affectedSatellites;
        std::set_difference(
            runtime.basis.receivers.begin(),
            runtime.basis.receivers.end(),
            rootReceivers.begin(),
            rootReceivers.end(),
            std::inserter(affectedReceivers, affectedReceivers.begin())
        );
        std::set_difference(
            runtime.basis.satellites.begin(),
            runtime.basis.satellites.end(),
            rootSatellites.begin(),
            rootSatellites.end(),
            std::inserter(affectedSatellites, affectedSatellites.begin())
        );
        return std::make_pair(affectedReceivers, affectedSatellites);
    };

    auto legacyLocalReinitialise = [&](const string& reason)
    {
        advanceDiscontinuousArcVersions();
        auto [affectedReceivers, affectedSatellites] = detachedNodes();
        if (affectedReceivers.empty() && affectedSatellites.empty())
        {
            return false;
        }

        set<ZhangGraphEdge> removedTreeEdges;
        set<ZhangGraphEdge> replacementEdges;
        std::set_difference(
            runtime.basis.treeEdges.begin(),
            runtime.basis.treeEdges.end(),
            candidate.treeEdges.begin(),
            candidate.treeEdges.end(),
            std::inserter(removedTreeEdges, removedTreeEdges.begin())
        );
        std::set_difference(
            candidate.treeEdges.begin(),
            candidate.treeEdges.end(),
            runtime.basis.treeEdges.begin(),
            runtime.basis.treeEdges.end(),
            std::inserter(replacementEdges, replacementEdges.begin())
        );
        int removedIntegerColumns = 0;
        for (const auto& [key, index] : kfState.kfIndexMap)
        {
            removedIntegerColumns +=
                key.type == KF::AMBIGUITY &&
                key.Sat.sys == sys &&
                zhangFullRankUsesObservable(
                    static_cast<E_ObsCode>(key.num),
                    options.baseline_observables
                ) &&
                (affectedReceivers.find(key.str) != affectedReceivers.end() ||
                 affectedSatellites.find(key.Sat) != affectedSatellites.end());
        }

		// The proposed product tree is known before the rectangular local
		// phase-coordinate reset.  Capture its immutable satellite functionals
		// now, while the old coordinate system still contains every historical
		// chord needed to define them.  Waiting until product output runs after
		// the reset made a subset of otherwise mature event targets unavailable
		// at t0 (notably the 2019-199 02:51 event group).
		set<ZhangGraphEdge> proposedProductEdges = candidate.edges;
		bool proposedProductAvailable = true;
		if (!singleS && options.product_core_min_satellite_support > 0)
		{
			const auto proposedCore = options.product_integer_support_core
				? zhangBuildIntegerSupportCore(
					candidate.edges, candidate.rootReceiver,
					runtime.productCoreReceivers,
					options.product_core_min_satellite_support,
					integerSupportQuality, {}, activeQuality, persistence).receiverCore
				: zhangBuildProductReceiverCore(
					candidate.edges, candidate.rootReceiver,
					runtime.productCoreReceivers,
					options.product_core_min_satellite_support,
					activeQuality, persistence);
			if (proposedCore.connected)
			{
				proposedProductEdges = proposedCore.edges;
			}
			else
			{
				proposedProductAvailable = false;
			}
		}
		ZhangGraphBasis proposedProduct = singleS ? candidate : proposedProductAvailable &&
			runtime.productInitialized
			? (options.product_integer_support_core
				? zhangBuildRiskAwareSpanningTree(
					proposedProductEdges, candidate.rootReceiver,
					runtime.productBasis.treeEdges, integerSupportQuality)
				: zhangBuildSpanningTree(
				proposedProductEdges,
				candidate.rootReceiver,
				runtime.productBasis.treeEdges,
				activeQuality,
				modelledEdges,
				persistence))
			: proposedProductAvailable
			? (options.product_integer_support_core
				? zhangBuildRiskAwareSpanningTree(
					proposedProductEdges, candidate.rootReceiver, {},
					integerSupportQuality)
				: zhangBuildRootedProductTree(
				proposedProductEdges,
				candidate.rootReceiver,
				{},
				activeQuality,
				modelledEdges,
				persistence))
			: ZhangGraphBasis{};
		if (proposedProductAvailable && proposedProduct.connected)
		{
			map<ZhangGraphEdge, map<int, int>> previousProductArcVersions;
			const auto frozenSignal = graphEventPreSignalArcVersions.find(
				{runtimeId, sys});
			if (frozenSignal != graphEventPreSignalArcVersions.end())
				previousProductArcVersions = frozenSignal->second;
			for (const auto& [edge, signals] :
			     runtime.productSignalArcVersions)
				previousProductArcVersions[edge] = signals;
			const auto& proposedArcVersions = runtime.phaseArcVersions;
			registerZhangCandidateProductSnapshotsBeforeCoordinateReset(
				trace,
				kfState,
				kfState,
				sys,
				options.baseline_observables,
				runtime.basis,
				runtime.productBasis,
				previousProductArcVersions,
				proposedProduct,
				proposedArcVersions,
				kfState.time);
		}

        if (!resetZhangGraphPhaseCoordinates(
                trace,
                kfState,
                sys,
                options.baseline_observables,
                affectedReceivers,
                affectedSatellites
            ))
        {
            return false;
        }

        const size_t preservedReceivers =
            runtime.basis.receivers.size() - affectedReceivers.size();
        const size_t preservedSatellites =
            runtime.basis.satellites.size() - affectedSatellites.size();

        if (!updateProductDatum("local_reinitialise"))
            invalidateZhangProduct(runtime, "LOCAL_REINITIALISE_PRODUCT_UNAVAILABLE");
        runtime.basis            = candidate;
        runtime.activeBasis      = candidate;
        runtime.observationEdges = observationEdges;
        runtime.stateEdges       = stateEdges;
        runtime.deferredEpochs   = 0;
        runtime.datumVersion++;
        runtime.representationVersion++;
        runtime.floatGaugeVersion++;
        traceGraphEvent(
            "local_reinitialise",
            reason,
            removedTreeEdges,
            replacementEdges,
            affectedReceivers,
            affectedSatellites,
            removedIntegerColumns,
            false
        );
        traceCanonicalAudit(runtime.basis, "local_reinitialise", false);
        recordZhangPhaseReinitialisation(
            kfState.time,
            sys,
            options.baseline_observables,
            reason,
            affectedSatellites
        );
        BOOST_LOG_TRIVIAL(warning)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=local_reinitialise"
            << " reason=" << reason
            << " affected_receivers=" << affectedReceivers.size()
            << " affected_satellites=" << affectedSatellites.size()
            << " preserved_receivers=" << preservedReceivers
            << " preserved_satellites=" << preservedSatellites
            << " phase_datum_discontinuity=local";
        return true;
    };

    // A disconnected FLOAT forest keeps each component's internal rows and
    // historical cross-covariances. It is not a globally certified AR chart.
    if (!candidate.connected || !runtime.basis.connected)
    {
        if (!localReinitialise("FLOAT_COMPONENT_SPLIT_OR_REJOIN"))
        {
            runtime.observationEdges.clear();
            invalidateZhangProduct(runtime, "FLOAT_COMPONENT_TRANSPORT_FAILED");
        }
        return;
    }

    if (!brokenTreeEdges.empty())
    {
        const set<ZhangGraphEdge> stateReconstructibleEdges =
            zhangStateReconstructibleEdges(
                kfState,
                sys,
                options.baseline_observables,
                runtime.basis,
                modelledEdges);
        set<ZhangGraphEdge> nonReconstructibleEdges;
        std::set_difference(
            modelledEdges.begin(), modelledEdges.end(),
            stateReconstructibleEdges.begin(), stateReconstructibleEdges.end(),
            std::inserter(nonReconstructibleEdges,
                          nonReconstructibleEdges.begin()));
        set<ZhangGraphEdge> pivotCandidates = stateReconstructibleEdges;
        if (zhangR51SingleBasisHandoffEnabled())
        {
            // Reconstructibility of a stale old FLOAT row is not current
            // physical eligibility. It may be retained as an event auxiliary
            // but cannot become a replacement tree edge or product support.
            pivotCandidates.clear();
            std::set_intersection(
                stateReconstructibleEdges.begin(),
                stateReconstructibleEdges.end(),
                stateCandidates.begin(), stateCandidates.end(),
                std::inserter(pivotCandidates, pivotCandidates.end()));
        }
        auto pivot = zhangPlanPivotBeforeRetire(
            runtime.basis,
            pivotCandidates,
            brokenTreeEdges);
        if (singleS && pivot.required)
        {
            const auto certified = certifiedTransportTreePlan(trace, kfState, sys, runtime,
                pivotCandidates, brokenTreeEdges, availability.discontinuousEdges,
                runtime.basis.rootReceiver, activeQuality,
                options.baseline_observables, previousSGaugeSignalArcs, availability.signalEdges);
            // FLOAT fallback remains legal and undergoes the existing full
            // mean/covariance transform. It is not an integer authorization.
            if (certified.basis.connected && certified.basis.receivers == runtime.basis.receivers &&
                certified.basis.satellites == runtime.basis.satellites)
            {
                pivot.replacementBasis = certified.basis;
                pivot.replacementEdges.clear();
                std::set_difference(certified.basis.treeEdges.begin(), certified.basis.treeEdges.end(),
                    runtime.basis.treeEdges.begin(), runtime.basis.treeEdges.end(),
                    std::inserter(pivot.replacementEdges, pivot.replacementEdges.begin()));
                pivot.connected = true;
                pivot.failureReason = "NONE";
            }
        }
        if (acsConfig.zhangPppAr.tree_slip_shadow_replay &&
            runtime.treeSlipShadowEvents <
                acsConfig.zhangPppAr.tree_slip_shadow_max_events)
        {
            runtime.treeSlipShadowEvents++;
            ZhangStateTransformAudit shadowAudit;
            bool exactTransformAvailable = false;
            if (pivot.connected)
            {
                KFState shadowState = kfState;
                ZhangGraphBasis shadowOldBasis = runtime.basis;
                shadowOldBasis.edges = modelledEdges;
                ZhangGraphBasis shadowNewBasis = pivot.replacementBasis;
                shadowNewBasis.edges = pivot.survivingRepresentedEdges;
                set<ZhangGraphEdge> retiredEdges;
                std::set_intersection(
                    modelledEdges.begin(), modelledEdges.end(),
                    availability.discontinuousEdges.begin(),
                    availability.discontinuousEdges.end(),
                    std::inserter(retiredEdges, retiredEdges.begin()));
                exactTransformAvailable = transformZhangGraphBasis(
                    trace,
                    shadowState,
                    sys,
                    options.baseline_observables,
                    shadowOldBasis,
                    shadowNewBasis,
                    false,
                    retiredEdges,
                    &shadowAudit);
            }
            shadowAudit.stateReconstructibleEdgeCount =
                static_cast<int>(stateReconstructibleEdges.size());

            const auto representedForest = zhangBuildSpanningTree(
                pivot.survivingRepresentedEdges,
                runtime.basis.rootReceiver);
            const int missingRepresentedNodes =
                static_cast<int>(
                    runtime.basis.receivers.size() +
                    runtime.basis.satellites.size() -
                    representedForest.receivers.size() -
                    representedForest.satellites.size());
            const int representedComponents =
                representedForest.componentCount + missingRepresentedNodes;
            const auto certifiedCore = zhangBuildIntegerSupportCore(
                pivot.survivingRepresentedEdges,
                runtime.basis.rootReceiver,
                runtime.productCoreReceivers,
                1,
                integerSupportQuality,
                {},
                activeQuality,
                persistence);
            const auto certifiedForest = zhangBuildSpanningTree(
                certifiedCore.qualifiedEdges,
                runtime.basis.rootReceiver);

            int slipEventEdges = 0;
            int l1OnlySlipEdges = 0;
            int l2OnlySlipEdges = 0;
            int dualFrequencySlipEdges = 0;
            int stationQcRemovalEdges = 0;
            int temporaryObservationLossEdges = 0;
            int receiverExitEdges = 0;
            int satelliteRiseSetEdges = 0;
            const auto& eventCodes =
                acsConfig.zhangFullRank.sysOpts.at(sys).baseline_observables;
            for (const auto& edge : brokenTreeEdges)
            {
                if (availability.discontinuousEdges.contains(edge))
                {
                    slipEventEdges++;
                    const auto signals = availability.discontinuitySignals.find(edge);
                    const bool l1 = signals != availability.discontinuitySignals.end() &&
                        eventCodes.size() == 2 && signals->second.contains(eventCodes[0]);
                    const bool l2 = signals != availability.discontinuitySignals.end() &&
                        eventCodes.size() == 2 && signals->second.contains(eventCodes[1]);
                    l1OnlySlipEdges += l1 && !l2;
                    l2OnlySlipEdges += l2 && !l1;
                    dualFrequencySlipEdges += l1 && l2;
                    continue;
                }
                if (availability.qcExcludedEdges.contains(edge) ||
                    availability.elevationExcludedEdges.contains(edge))
                {
                    stationQcRemovalEdges++;
                    continue;
                }
                if (availability.signalUnavailableEdges.contains(edge))
                {
                    temporaryObservationLossEdges++;
                    continue;
                }
                bool receiverStillObserved = false;
                bool satelliteStillObserved = false;
                for (const auto& rawEdge : availability.rawEdges)
                {
                    receiverStillObserved |= rawEdge.receiver == edge.receiver;
                    satelliteStillObserved |= rawEdge.satellite == edge.satellite;
                }
                if (!receiverStillObserved) receiverExitEdges++;
                else if (!satelliteStillObserved) satelliteRiseSetEdges++;
                else temporaryObservationLossEdges++;
            }

            trace << "\nZHANG_TREE_SLIP_SHADOW_REPLAY time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " event_index=" << runtime.treeSlipShadowEvents
                  << " broken_tree_edge_count=" << brokenTreeEdges.size()
                  << " slip_edge_count=" << slipEventEdges
                  << " l1_only_slip_edge_count=" << l1OnlySlipEdges
                  << " l2_only_slip_edge_count=" << l2OnlySlipEdges
                  << " dual_frequency_slip_edge_count="
                  << dualFrequencySlipEdges
                  << " station_qc_removal_edge_count="
                  << stationQcRemovalEdges
                  << " temporary_observation_loss_edge_count="
                  << temporaryObservationLossEdges
                  << " receiver_exit_edge_count=" << receiverExitEdges
                  << " satellite_rise_set_edge_count="
                  << satelliteRiseSetEdges
                  << " old_tree_components_after_deletion="
                  << 1 + brokenTreeEdges.size()
                  << " represented_same_arc_component_count="
                  << representedComponents
                  << " abstract_represented_edge_count=" << modelledEdges.size()
                  << " state_reconstructible_edge_count="
                  << stateReconstructibleEdges.size()
                  << " non_reconstructible_edges="
                  << edgeList(nonReconstructibleEdges)
                  << " replacement_rank="
                  << pivot.replacementBasis.treeEdges.size()
                  << " exact_state_transform_available="
                  << exactTransformAvailable
                  << " certified_transport_rank="
                  << certifiedForest.treeEdges.size()
                  << " float_only_transport_rank="
                  << pivot.replacementBasis.treeEdges.size()
                  << " observable_rows=" << shadowAudit.observableRows
                  << " mean_absolute_norm=" << shadowAudit.meanAbsoluteNorm
                  << " mean_maximum_absolute="
                  << shadowAudit.meanMaximumAbsolute
                  << " covariance_absolute_norm="
                  << shadowAudit.covarianceAbsoluteNorm
                  << " covariance_relative_norm="
                  << shadowAudit.covarianceRelativeNorm
                  << " transform_failure_reason="
                  << shadowAudit.failureReason
                  << " failure_observable=" << shadowAudit.observable
                  << " failure_receiver=" << shadowAudit.receiver
                  << " failure_satellite=" << shadowAudit.satellite
                  << " missing_kf_key=" << shadowAudit.missingKfKey
                  << " old_node_count=" << shadowAudit.oldNodeCount
                  << " new_node_count=" << shadowAudit.newNodeCount
                  << " authoritative_feedback=0";

            // Experiment 4 must not exercise the new authoritative branch.
            // Preserve the legacy component-reset behaviour on the real state.
            if (legacyLocalReinitialise("tree_edge_arc_break"))
            {
                return;
            }
            retainOldTreeRootComponent();
            return;
        }
        if (!acsConfig.zhangPppAr.tree_slip_pivot_before_retire)
        {
            if (acsConfig.zhangPppAr.output_diagnostics)
            {
                trace << "\nZHANG_TREE_SLIP_TRANSACTION time="
                      << kfState.time.to_string(0)
                      << " system=" << enum_to_string(sys)
                      << " exact_state_transform=0"
                      << " status=LEGACY_LOCAL_REINITIALISE_CONTROL"
                      << " float_gauge_continuous=0"
                      << " integer_cross_component_valid=0";
            }
            if (legacyLocalReinitialise("tree_edge_arc_break"))
            {
                return;
            }
            retainOldTreeRootComponent();
            return;
        }
        if (pivot.connected)
        {
            ZhangGraphBasis oldBasis = runtime.basis;
            oldBasis.edges = modelledEdges;
            oldBasis.receivers.clear();
            oldBasis.satellites.clear();
            for (const auto& edge : modelledEdges)
            {
                oldBasis.receivers.insert(edge.receiver);
                oldBasis.satellites.insert(edge.satellite);
            }
            ZhangGraphBasis transformedBasis = pivot.replacementBasis;
            transformedBasis.edges = pivot.survivingRepresentedEdges;
            set<ZhangGraphEdge> retiredRepresentedEdges;
            std::set_intersection(
                modelledEdges.begin(), modelledEdges.end(),
                availability.discontinuousEdges.begin(),
                availability.discontinuousEdges.end(),
                std::inserter(retiredRepresentedEdges,
                              retiredRepresentedEdges.begin()));
            auto postPivotStateEdges = stateEdges;
            for (const auto& edge : availability.discontinuousEdges)
                if (availability.edges.contains(edge) &&
                    transformedBasis.receivers.contains(edge.receiver) &&
                    transformedBasis.satellites.contains(edge.satellite))
                    postPivotStateEdges.insert(edge);
            const auto eligibilityParts = partitionZhangEligibility(
                transformedBasis, postPivotStateEdges);
            const auto activeAfterPivot = zhangActiveGraphBasisWithAugmentation(
                transformedBasis, postPivotStateEdges);
            if (!activeAfterPivot.valid)
            {
                trace << "\nR48_ACTIVE_BASIS_REJECT time=" << kfState.time.to_string(0)
                      << " phase=PIVOT reason=" << activeAfterPivot.failureReason;
                if (!localReinitialise("PIVOT_CHART_AUGMENTATION"))
                    retainOldTreeRootComponent();
                return;
            }
            const int representationBefore = runtime.representationVersion;
            const int floatGaugeBefore = runtime.floatGaugeVersion;
            const int integerComponentBefore =
                runtime.integerComponentVersion;

            KFState transactionKfState;
            ZhangStateTransformAudit authoritativeAudit;
            if (transformZhangGraphBasis(
                    trace,
                    kfState,
                    sys,
                    options.baseline_observables,
                    oldBasis,
                    transformedBasis,
                    true,
                    retiredRepresentedEdges,
                    &authoritativeAudit,
                    &transactionKfState,
                    &availability.discontinuitySignals))
            {
                // Only after the atomic S-transform/projection succeeds does
                // the physical arc identity advance.  The post-slip sample is
                // then admitted as a fresh chord observation, never as a
                // replacement-tree edge.
                advanceDiscontinuousArcVersions();
				candidate = transformedBasis;
				candidate.edges = pivot.survivingRepresentedEdges;
				candidate.componentCount = 1;
				candidate.connected = true;
				modelledEdges = pivot.survivingRepresentedEdges;
                if (!updateProductDatum("pivot_before_retire", true))
                    invalidateZhangProduct(runtime, "PIVOT_PRODUCT_DATUM_UNAVAILABLE");
                for (const auto& edge : availability.discontinuousEdges)
                {
                    if (availability.edges.find(edge) ==
                            availability.edges.end() ||
                        transformedBasis.receivers.find(edge.receiver) ==
                            transformedBasis.receivers.end() ||
                        transformedBasis.satellites.find(edge.satellite) ==
                            transformedBasis.satellites.end())
                    {
                        continue;
                    }
                    stateEdges.insert(edge);
                    observationEdges.insert(edge);
                }
                candidate = activeAfterPivot.basis;
                trace << "\nR48_ACTIVE_BASIS_COMMIT time=" << kfState.time.to_string(0)
                      << " phase=PIVOT retained_datum_only_edges="
                      << edgeList(activeAfterPivot.retainedDatumOnlyEdges)
                      << " coordinate_tree_unchanged=1 observation_eligibility_unchanged=1";
                runtime.representationVersion++;
                runtime.basis = zhangRepresentedBasisAfterAugmentation(transformedBasis, candidate);
                trace << "\nZHANG_GRAPH_NODE_AUGMENTATION time="
                      << kfState.time.to_string(0) << " phase=PIVOT"
                      << " pending_edges=" << edgeList(eligibilityParts.requiresAugmentation)
                      << " old_chart_transported_first=1 integer_certification=0";
                runtime.activeBasis = candidate;
                runtime.observationEdges = observationEdges;
                runtime.stateEdges = stateEdges;
                runtime.deferredEpochs = 0;
                traceGraphEvent(
                    "pivot_before_retire",
                    retiredRepresentedEdges.empty()
                        ? "tree_edge_observation_support_change"
                        : "tree_edge_physical_arc_break",
                    pivot.slippedTreeEdges,
                    pivot.replacementEdges,
                    {},
                    {},
                    static_cast<int>(retiredRepresentedEdges.size() *
                        options.baseline_observables.size()),
                    true);
                trace << "\nZHANG_TREE_SLIP_TRANSACTION time="
                      << kfState.time.to_string(0)
                      << " system=" << enum_to_string(sys)
                      << " pre_event_represented_edges=" << oldBasis.edges.size()
                      << " slipped_edges=" << retiredRepresentedEdges.size()
                      << " surviving_represented_edges="
                      << pivot.survivingRepresentedEdges.size()
                      << " replacement_rank="
                      << pivot.replacementBasis.treeEdges.size()
                      << " component_count="
                      << pivot.replacementBasis.componentCount
                      << " pre_event_arc_versions="
                      << (firstPreEventVersion
                            ? "NONE"
                            : preEventVersions.str())
                      << " exact_state_transform=1"
                      << " retired_old_arc_columns="
                      << retiredRepresentedEdges.size() *
                            options.baseline_observables.size()
                      << " fresh_chord_edges="
                      << availability.discontinuousEdges.size()
                      << " representation_version_before="
                      << representationBefore
                      << " representation_version_after="
                      << runtime.representationVersion
                      << " float_gauge_version_before=" << floatGaugeBefore
                      << " float_gauge_version_after="
                      << runtime.floatGaugeVersion
                      << " integer_component_version_before="
                      << integerComponentBefore
                      << " integer_component_version_after="
                      << runtime.integerComponentVersion
                      << " phase_posterior_preserved=1"
                      << " covariance_cross_terms_preserved=1"
                      << " action=PIVOT_THEN_RETIRE_THEN_FRESH_CHORD";
                traceCanonicalAudit(candidate, "pivot_before_retire", true);
                return;
            }
            trace << "\nZHANG_TREE_SLIP_TRANSACTION time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " exact_state_transform=0"
                  << " status=TRANSACTION_ABORTED_POSTERIOR_PRESERVED"
                  << " reason=" << authoritativeAudit.failureReason
                  << " observable=" << authoritativeAudit.observable
                  << " receiver=" << authoritativeAudit.receiver
                  << " satellite=" << authoritativeAudit.satellite
                  << " missing_kf_key=" << authoritativeAudit.missingKfKey
                  << " abstract_represented_edge_count=" << modelledEdges.size()
                  << " state_reconstructible_edge_count="
                  << stateReconstructibleEdges.size();
            if (localReinitialise("CONNECTED_CHART_EXACT_TRANSPORT_UNAVAILABLE")) return;
            retainOldTreeRootComponent();
            BOOST_LOG_TRIVIAL(error)
                << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
                << " action=abort pivot-before-retire transaction"
                << " reason=exact transform failed on connected represented graph"
                << ", old posterior and arc versions preserved";
            return;
        }
        else if (acsConfig.zhangPppAr.output_diagnostics)
        {
            trace << "\nZHANG_TREE_SLIP_TRANSACTION time="
                  << kfState.time.to_string(0)
                  << " system=" << enum_to_string(sys)
                  << " exact_state_transform=0"
                  << " surviving_represented_edges="
                  << pivot.survivingRepresentedEdges.size()
                  << " status=COMPONENT_SPLIT_ALLOWED"
                  << " reason=" << pivot.failureReason;
        }
        if (localReinitialise("tree_edge_arc_break"))
        {
            return;
        }

        retainOldTreeRootComponent();
        BOOST_LOG_TRIVIAL(error)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=defer tree exchange: local arc-break reset failed"
            << ", broken_tree_edges=" << brokenTreeEdges.size()
            << ", retained_safe_edges=" << runtime.observationEdges.size();
        return;
    }

    if (options.core_skeleton)
    {
        runtime.deferredEpochs = 0;
        if (localReinitialise("CORE_SKELETON_PHYSICAL_SUPPORT_CHANGE")) return;
        retainOldTreeRootComponent();
        BOOST_LOG_TRIVIAL(warning)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=core_skeleton_hold"
            << " reason=non_leaf_tree_change"
            << " retained_observation_edges=" << runtime.observationEdges.size();
        return;
    }

    bool newTreeRepresented = std::includes(
        modelledEdges.begin(),
        modelledEdges.end(),
        candidate.treeEdges.begin(),
        candidate.treeEdges.end()
    );
    if (!newTreeRepresented)
    {
        runtime.deferredEpochs++;
        if (localReinitialise("replacement_edge_without_prior_state"))
        {
            return;
        }

        retainOldTreeRootComponent();
        BOOST_LOG_TRIVIAL(warning)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=defer tree exchange: a replacement edge has no prior state"
            << ", retained_safe_edges=" << runtime.observationEdges.size();
        return;
    }

    ZhangGraphBasis oldBasis = runtime.basis;
    oldBasis.edges           = modelledEdges;
    oldBasis.receivers.clear();
    oldBasis.satellites.clear();
    for (const auto& edge : modelledEdges)
    {
        oldBasis.receivers.insert(edge.receiver);
        oldBasis.satellites.insert(edge.satellite);
    }

    ZhangGraphBasis transformedBasis = candidate;
    transformedBasis.edges           = modelledEdges;

    if (!transformZhangGraphBasis(
            trace,
            kfState,
            sys,
            options.baseline_observables,
            oldBasis,
            transformedBasis
        ))
    {
        runtime.deferredEpochs++;
        if (localReinitialise("exact_state_transform_unavailable"))
        {
            return;
        }

        retainOldTreeRootComponent();
        BOOST_LOG_TRIVIAL(warning)
            << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
            << " action=defer tree exchange: exact state transform failed"
            << ", retained_safe_edges=" << runtime.observationEdges.size();
        return;
    }

    set<ZhangGraphEdge> removedTreeEdges;
    set<ZhangGraphEdge> replacementEdges;
    std::set_difference(
        oldBasis.treeEdges.begin(),
        oldBasis.treeEdges.end(),
        candidate.treeEdges.begin(),
        candidate.treeEdges.end(),
        std::inserter(removedTreeEdges, removedTreeEdges.begin())
    );
    std::set_difference(
        candidate.treeEdges.begin(),
        candidate.treeEdges.end(),
        oldBasis.treeEdges.begin(),
        oldBasis.treeEdges.end(),
        std::inserter(replacementEdges, replacementEdges.begin())
    );

    if (!updateProductDatum("tree_exchange"))
        invalidateZhangProduct(runtime, "TREE_EXCHANGE_PRODUCT_UNAVAILABLE");
    runtime.basis            = zhangRepresentedBasisAfterAugmentation(transformedBasis, candidate);
    runtime.activeBasis      = candidate;
    runtime.observationEdges = observationEdges;
    runtime.stateEdges       = stateEdges;
    runtime.deferredEpochs = 0;
    runtime.representationVersion++;
    traceGraphEvent(
        "tree_exchange",
        "exact_state_transform",
        removedTreeEdges,
        replacementEdges,
        {},
        {},
        0,
        true
    );
    // transformedBasis retains stale ambiguity arcs solely to make the state
    // transform dimension preserving.  Canonical integer coordinates belong
    // to the active retained component represented by candidate/stateEdges.
    traceCanonicalAudit(candidate, "tree_exchange", true);

    BOOST_LOG_TRIVIAL(info)
        << "ZHANG_GRAPH_BASIS sys=" << enum_to_string(sys)
        << " action=tree_exchange"
        << " modelled_edges=" << modelledEdges.size()
        << " tree_edges=" << transformedBasis.treeEdges.size()
        << " cycles=" << modelledEdges.size() - transformedBasis.treeEdges.size();
}
}  // namespace

bool zhangPrepareCertifiedTreeTransportHistory(Trace& trace, const KFState& authoritativeFloat,
    const std::map<E_Sys, std::pair<std::string, std::string>>& committedScopes,
    ZhangPreparedCertifiedTreeTransportHistory& prepared,
    std::string* failureReason)
{
    using namespace zhang_transport_certificate;
    const string owner = zhangGraphRuntimeId(authoritativeFloat);
    const auto epoch = static_cast<long int>(authoritativeFloat.time.bigTime);
    auto fail = [&](const string& reason)
    {
        if (failureReason) *failureReason = reason;
        trace << "\nZHANG_TRANSPORT_CERTIFICATE_COMMIT time=" << authoritativeFloat.time.to_string(0)
            << " owner=" << owner << " status=NOT_IMPORTED reason=" << reason
            << " old_directory_retained=1 current_ar_authority=0";
        return false;
    };
    prepared = {};
    if (committedScopes.empty()) return true; // No empty tag/schema mutation on checkpoint recovery.
    if (!validZhangGraphRuntimeId(owner) ||
        authoritativeFloat.metaDataMap.contains(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA) ||
        !authoritativeFloat.lastPppTransactionCommitted || !authoritativeFloat.lastPppQcConverged)
        return fail("TRANSPORT_COMMIT_DESTINATION_NOT_AUTHORITATIVE_FLOAT");
    // kfState is the authoritative FLOAT state, never the conditioned AR copy.
    // The atomic product COMMIT caller pins the admitted source scopes below.
    const auto root = zhangR51NumericRoot(authoritativeFloat.x, authoritativeFloat.P);
    const auto order = zhangCurrentEvidenceStateOrderId(authoritativeFloat.kfIndexMap);
    map<string, string> proposedMetadata;
    for (const auto& [system, scope] : committedScopes)
    {
        const auto found = zhangCurrentEpochEvidenceRegistry().find({owner, system});
        if (found == zhangCurrentEpochEvidenceRegistry().end())
            return fail("TRANSPORT_COMMIT_SAME_OWNER_CURRENT_EVIDENCE_MISSING");
        const auto& evidence = found->second;
        string bindingReason;
        if (!evidence.valid || evidence.runtimeId != owner || evidence.system != system ||
            !evidence.matches(epoch, scope.first, scope.second) || evidence.floatStateRoot != root ||
            !evidence.floatPosteriorPsdValidated || evidence.floatPosteriorPsdRootId != root ||
            evidence.floatPosteriorPsdStateOrderId != order ||
            !evidence.validatePhysicalNumericBindings(bindingReason) ||
            evidence.rows.size() != evidence.values.size() || evidence.rows.size() != evidence.rowProofs.size() ||
            evidence.eventOverlapRows.size() != evidence.eventOverlapValues.size() ||
            evidence.eventOverlapRows.size() != evidence.eventOverlapProofs.size())
            return fail("TRANSPORT_COMMIT_CURRENT_CERTIFICATE_OR_FULL_ROOT_MISMATCH:" + bindingReason);
        Directory directory;
        std::shared_ptr<const CertifiedTransportQueryCache> priorDomain;
        const auto metadata = authoritativeFloat.metaDataMap.find(certifiedTransportMetadataKey(system));
        if (metadata != authoritativeFloat.metaDataMap.end())
        {
            string reason;
            priorDomain = certifiedTransportQuery(authoritativeFloat, system, reason);
            if (!priorDomain || !priorDomain->query->feasible)
                return fail("TRANSPORT_COMMIT_PRIOR_DOMAIN_INVALID:" + reason);
            directory = priorDomain->directory;
        }
        else { directory.authoritativeOwner = owner; directory.system = static_cast<int>(system); }
        const auto priorFactCount = directory.facts.size();
        std::size_t impliedRowsSkipped = 0;
        map<Row, ZhangExactInteger> existing;
        for (const auto& fact : directory.facts) existing.emplace(fact.physical, fact.value);
        const SourceScope sourceScope{epoch, scope.first, scope.second, root, order};
        const string sourceScopeId = zhangCheckpointSha256(owner + "|" + enum_to_string(system) + "|" +
            std::to_string(epoch) + "|" + scope.first + "|" + scope.second + "|" + root + "|" + order);
        auto append = [&](const auto& rows, const auto& values, const auto& groups)
        {
            for (std::size_t r = 0; r < rows.size(); ++r)
            {
                const auto old = existing.find(rows[r]);
                if (old != existing.end()) { if (old->second != values[r]) return false; else continue; }
                // LAMBDA may select a different equivalent integer basis at
                // every epoch. Exact-row deduplication alone would archive
                // all of those redundant bases and repeatedly copy complete
                // source orders/proof DAGs. This directory has TRANSPORT
                // consequence authority only: a feasible old domain already
                // proving this constant needs no new fact or current RHS.
                // Do not use this rule for CURRENT integer-combination proof.
                ZhangExactInteger impliedValue;
                if (priorDomain && priorDomain->query->consequence(rows[r], impliedValue))
                {
                    if (impliedValue != values[r]) return false;
                    ++impliedRowsSkipped;
                    continue;
                }
                const auto closure = zhangDecisionRiskClosure(groups[r]);
                if (!closure.valid || closure.atoms.empty() || rows[r].empty()) return false;
                Fact fact;
                fact.physical = rows[r]; fact.value = values[r]; fact.sourceScopeId = sourceScopeId;
                if (!directory.sourceScopes.contains(sourceScopeId)) directory.sourceScopes.emplace(sourceScopeId, sourceScope);
                const string provenance = owner + "|PPP_COMMITTED|E" + std::to_string(epoch);
                for (const auto& [id, atom] : closure.atoms)
                {
                    if (!atom || (atom->observationProvenance != provenance &&
                        atom->observationProvenance.rfind(provenance + "|", 0) != 0)) return false;
                    const string content = zhangDecisionProofAtomContent(*atom);
                    const auto [old, inserted] = directory.proofAtoms.emplace(id, content);
                    if (!inserted && old->second != content) return false;
                }
                // DAG roots reference the lossless shared atom dictionary;
                // do not copy the entire transitive closure for every row.
                set<string> uniqueRoots;
                for (const auto& proof : groups[r])
                    if (uniqueRoots.insert(proof->id).second) fact.proofIds.push_back(proof->id);
                existing.emplace(fact.physical, fact.value);
                directory.facts.push_back(std::move(fact));
                if (directory.facts.size() > maximumFacts) return false;
            }
            return true;
        };
        if (!append(evidence.rows, evidence.values, evidence.rowProofs) ||
            !append(evidence.eventOverlapRows, evidence.eventOverlapValues, evidence.eventOverlapProofs))
            return fail("TRANSPORT_COMMIT_CONFLICT_UNBOUND_PROOF_OR_CAPACITY");
        if (directory.facts.empty()) continue;
        if (directory.facts.size() == priorFactCount)
        {
            // No new information was imported. Preserve both the exact seal
            // and its genuine last source epoch; do not label old facts as a
            // fresh epoch's certificate merely because they were re-queried.
            trace << "\nZHANG_TRANSPORT_CERTIFICATE_REUSE time=" << authoritativeFloat.time.to_string(0)
                << " owner=" << owner << " system=" << enum_to_string(system)
                << " prior_facts=" << priorFactCount << " implied_rows_skipped=" << impliedRowsSkipped
                << " sealed_bytes_unchanged=1 source_epoch_unchanged=1 current_ar_authority=0";
            continue;
        }
        directory.lastCommittedEpoch = epoch;
        const Query proposedDomain(directory);
        if (!proposedDomain.feasible)
            return fail("TRANSPORT_COMMIT_HISTORY_INTEGER_FEASIBILITY_REJECTED:" + proposedDomain.reason);
        const string encoded = directory.encode();
        if (encoded.size() > maximumBytes) return fail("TRANSPORT_COMMIT_BYTE_CAPACITY_OLD_DIRECTORY_RETAINED");
        proposedMetadata.emplace(certifiedTransportMetadataKey(system), zhangCheckpointSha256(encoded) + "\n" + encoded);
        trace << "\nZHANG_TRANSPORT_CERTIFICATE_STAGED time=" << authoritativeFloat.time.to_string(0)
            << " source_owner=" << evidence.runtimeId << " destination_owner=" << owner
            << " identity_transfer=1 role=CURRENT_ACCEPTED_TO_TRANSPORT_ONLY source_epoch=" << epoch
            << " source_posterior=" << scope.first << " source_chart=" << scope.second
            << " source_root=" << root << " source_order_sha256=" << zhangCheckpointSha256(order)
            << " facts=" << directory.facts.size() << " current_ar_authority=0";
    }
    prepared.authoritativeOwner = owner; prepared.sourceEpoch = epoch;
    prepared.sourceRoot = root; prepared.sourceOrder = order;
    prepared.sealedDirectories = std::move(proposedMetadata);
    prepared.sourceCertificateVerified = true;
    if (failureReason) *failureReason = "TRANSPORT_ONLY_CURRENT_SOURCE_VERIFIED_NOT_YET_IMPORTED";
    return true;
}

bool zhangStageCertifiedTreeTransportHistoryCommit(const KFState& authoritativeFloat,
    const ZhangPreparedCertifiedTreeTransportHistory& prepared,
    std::map<std::string,std::string>& proposedMetadata, std::string* failureReason)
{
    if (!prepared.sourceCertificateVerified || prepared.sealedDirectories.empty())
    { if (failureReason) *failureReason = "NO_PREPARED_TRANSPORT_FACTS"; return false; }
    string reason;
    const string owner = zhangGraphRuntimeId(authoritativeFloat);
    if (owner != prepared.authoritativeOwner ||
        static_cast<long int>(authoritativeFloat.time.bigTime) != prepared.sourceEpoch ||
        authoritativeFloat.metaDataMap.contains(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA) ||
        !authoritativeFloat.lastPppTransactionCommitted || !authoritativeFloat.lastPppQcConverged ||
        !zhangValidateCertifiedTreeTransportMetadata(prepared.sealedDirectories, owner, reason))
    {
        if (failureReason) *failureReason = "TRANSPORT_COMMIT_OWNER_EPOCH_OR_STAGED_SEAL_MISMATCH:" + reason;
        return false;
    }
    // Only the source posterior was used for numerical certification. The
    // current FLOAT owner may now be a validated AR_EXIT marginal; historical
    // facts concern immutable physical arcs, not those removed coordinates.
    proposedMetadata = authoritativeFloat.metaDataMap;
    for (const auto& [key, value] : prepared.sealedDirectories) proposedMetadata[key] = value;
    if (failureReason) *failureReason = "TRANSPORT_ONLY_COMMIT_READY_NOT_YET_IMPORTED";
    return true;
}

bool zhangCommitCertifiedTreeTransportHistory(Trace& trace, KFState& authoritativeFloat,
    const ZhangPreparedCertifiedTreeTransportHistory& prepared, std::string* failureReason)
{
    std::map<std::string,std::string> proposedMetadata;
    if (!zhangStageCertifiedTreeTransportHistoryCommit(authoritativeFloat, prepared,
            proposedMetadata, failureReason)) return false;
    authoritativeFloat.metaDataMap.swap(proposedMetadata);
    const string owner = zhangGraphRuntimeId(authoritativeFloat);
    trace << "\nZHANG_TRANSPORT_CERTIFICATE_OWNER_TRANSFER time=" << authoritativeFloat.time.to_string(0)
        << " source_owner=" << prepared.authoritativeOwner << " destination_owner=" << owner
        << " source_root=" << prepared.sourceRoot << " source_epoch=" << prepared.sourceEpoch
        << " identity_owner_transfer=1 source_posterior_verified=1 atomic_product_commit_complete=1"
        << " target_role=TRANSPORT_ONLY current_ar_authority=0 scopes=" << prepared.sealedDirectories.size();
    if (failureReason) *failureReason = "TRANSPORT_ONLY_COMMIT_IMPORTED";
    return true;
}

ZhangHistoricalTransportBridgeProof zhangQueryCertifiedHistoricalTransport(
    const KFState& state, E_Sys system, const ZhangProductIntegerTransition& difference)
{
    ZhangHistoricalTransportBridgeProof result;
    result.owner = zhangGraphRuntimeId(state);
    if (!difference.valid || difference.observable == E_ObsCode::NONE ||
        difference.physicalEdges.size() != difference.physicalArcVersions.size() ||
        difference.physicalEdges.size() != difference.coefficients.size() ||
        state.metaDataMap.contains(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA))
    { result.reason = "HISTORICAL_TRANSPORT_TARGET_OWNER_OR_SHAPE_INVALID"; return result; }
    zhang_transport_certificate::Row target;
    for (std::size_t i = 0; i < difference.physicalEdges.size(); ++i)
    {
        if (difference.coefficients[i] == 0) continue;
        const auto& edge = difference.physicalEdges[i];
        if (edge.receiver.empty() || edge.satellite.sys != system || difference.physicalArcVersions[i] < 0)
        { result.reason = "HISTORICAL_TRANSPORT_PHYSICAL_VERSION_INVALID"; return result; }
        target[enum_to_string(difference.observable) + "|" + edge.receiver + "|" + edge.satellite.id() +
            "|V" + std::to_string(difference.physicalArcVersions[i])] += difference.coefficients[i];
    }
    for (auto it = target.begin(); it != target.end();)
        if (it->second == 0) it = target.erase(it); else ++it;
    if (target.empty()) { result.reason = "HISTORICAL_TRANSPORT_EMPTY_TARGET_NOT_AN_AR_CERTIFICATE"; return result; }
    const auto cached = certifiedTransportQuery(state, system, result.reason);
    if (!cached || !cached->query->consequence(target, result.value))
    {
        if (cached && cached->query->feasible) result.reason = "HISTORICAL_TRANSPORT_DIRECTION_UNCERTIFIED";
        return result;
    }
    result.value += difference.affineOffsetCycles;
    result.proven = true;
    result.lastSourceEpoch = cached->directory.lastCommittedEpoch;
    result.proofIds = cached->proofIds;
    result.proofStatements = cached->proofStatements;
    result.reason = "PRIOR_COMMIT_EXACT_DATUM_TRANSPORT_ONLY_CURRENT_NATIVE_AR_STILL_REQUIRED";
    return result;
}

std::vector<ZhangGraphIntegerContext> zhangGraphIntegerComponentContexts(
    const KFState& state, E_Sys system, std::string* failureReason)
{
    auto reject = [&](const string& reason)
    {
        if (failureReason) *failureReason = reason;
        return vector<ZhangGraphIntegerContext>{};
    };
    const string runtimeId = zhangGraphRuntimeId(state);
    const auto owner = graphStateMap.find({runtimeId, system});
    if (owner == graphStateMap.end() || !owner->second.initialized)
        return reject("CURRENT_GRAPH_RUNTIME_NOT_INITIALIZED");
    const auto& runtime = owner->second;
    if (runtime.productBasisMode != configuredProductBasisMode())
        return reject("PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION");
    const bool singleS = singleSBasisNative();
    const auto& source = singleS || runtime.activeBasis.edges.empty()
        ? runtime.basis : runtime.activeBasis;
    // A view is not a new S parameterisation. Refuse stale/imported graphs
    // whose tree/root do not describe the committed KF coordinate system.
    if (source.rootReceiver != runtime.basis.rootReceiver ||
        source.treeEdges != runtime.basis.treeEdges)
        return reject("ACTIVE_BASIS_NOT_COMMITTED_COORDINATE_TREE");
    auto options = acsConfig.zhangFullRank.sysOpts.find(system);
    if (options == acsConfig.zhangFullRank.sysOpts.end() ||
        options->second.baseline_observables.size() != 2)
        return reject("COMPONENT_BASELINE_OBSERVABLES_NOT_CONFIGURED");
    if (!runtime.signalObservationFinalized ||
        runtime.phaseAdmissionTime != state.time.to_string(9) ||
        runtime.signalObservationFactorSequence == 0)
        return reject("CURRENT_SIGNAL_FACTOR_ROWS_NOT_FINALIZED");
    const auto firstCode = int(options->second.baseline_observables[0]);
    const auto secondCode = int(options->second.baseline_observables[1]);
    const auto first = runtime.signalObservationEdges.find(firstCode);
    const auto second = runtime.signalObservationEdges.find(secondCode);
    set<ZhangGraphEdge> currentProductEdges;
    if (first != runtime.signalObservationEdges.end() &&
        second != runtime.signalObservationEdges.end())
        std::set_intersection(first->second.begin(), first->second.end(),
            second->second.begin(), second->second.end(),
            std::inserter(currentProductEdges, currentProductEdges.end()));
    for (auto edge = currentProductEdges.begin(); edge != currentProductEdges.end();)
    {
        const auto serial = runtime.phaseArcVersions.find(*edge);
        bool represented = source.edges.contains(*edge) && serial != runtime.phaseArcVersions.end();
        if (serial != runtime.phaseArcVersions.end())
            represented &= serial->second.contains(firstCode) && serial->second.contains(secondCode);
        for (const auto code : options->second.baseline_observables)
        {
            represented &= state.kfIndexMap.contains(zhangTransportSatelliteKey(code, edge->satellite));
            if (edge->receiver != source.rootReceiver)
                represented &= state.kfIndexMap.contains(zhangTransportReceiverKey(system, code, edge->receiver));
            if (!source.treeEdges.contains(*edge))
                represented &= state.kfIndexMap.contains(zhangTransportAmbiguityKey(code, *edge));
        }
        edge = represented ? std::next(edge) : currentProductEdges.erase(edge);
    }
    vector<ZhangGraphIntegerContext> result;
    for (const auto& component : zhangEpochComponentGraphs(
            source, runtime.productBasis, currentProductEdges, singleS))
    {
        if (!component.valid) continue;
        ZhangGraphIntegerContext context;
        context.basis = component.basis;
        context.productBasis = component.productBasis;
        context.productBasisMode = runtime.productBasisMode;
        context.componentId = component.componentId;
        context.runtimeIdentity = runtimeId;
        context.snapshotEpoch = state.time.to_string(9);
        context.productBasisReason = component.failureReason;
        context.sourceComponentCount = component.sourceComponentCount;
        context.cycleRank = component.cycleRank;
        context.containsSourceRoot = component.containsSourceRoot;
        context.productBasisRebuilt = component.productBasisRebuilt;
        context.eventId = runtime.eventCounter;
        context.productDatumVersion = runtime.productDatumVersion;
        context.committedSGaugeVersion = runtime.committedSGaugeVersion;
        context.representationVersion = runtime.representationVersion;
        context.floatGaugeVersion = runtime.floatGaugeVersion;
        context.integerComponentVersion = runtime.integerComponentVersion;
        context.signalObservationFinalized = runtime.signalObservationFinalized;
        for (const auto& edge : component.basis.edges)
        {
            const auto legacy = runtime.edgeHistory.find(edge);
            if (legacy != runtime.edgeHistory.end())
                context.arcVersions[edge] = legacy->second.arcVersion;
            const auto physical = runtime.phaseArcVersions.find(edge);
            if (physical != runtime.phaseArcVersions.end())
                context.signalArcVersions[edge] = physical->second;
            for (const auto& [code, admitted] : runtime.signalObservationEdges)
                if (admitted.contains(edge)) context.signalObservationEdges[code].insert(edge);
        }
        context.initialized = true;
        context.committedSBasisIdentity = committedSBasisIdentity(context);
        result.push_back(std::move(context));
    }
    if (failureReason) *failureReason = result.empty()
        ? "NO_VALID_COMMITTED_INTEGER_COMPONENT" : "NONE";
    return result;
}

ZhangGraphIntegerContextScope::ZhangGraphIntegerContextScope(
    const KFState& state, E_Sys system, const ZhangGraphIntegerContext& context)
{
    const string runtime = zhangGraphRuntimeId(state);
    if (!context.initialized || context.snapshotEpoch != state.time.to_string(9) ||
        context.runtimeIdentity != runtime || !context.basis.connected ||
        !context.basis.receivers.contains(context.basis.rootReceiver))
        throw std::invalid_argument("INVALID_EPOCH_COMPONENT_CONTEXT_SELECTION");
    const auto owner = graphStateMap.find({runtime, system});
    if (owner == graphStateMap.end() || !owner->second.initialized ||
        owner->second.productBasisMode != configuredProductBasisMode() ||
        context.productBasisMode != owner->second.productBasisMode ||
        context.eventId != owner->second.eventCounter ||
        context.representationVersion != owner->second.representationVersion ||
        context.productDatumVersion != owner->second.productDatumVersion ||
        context.committedSGaugeVersion != owner->second.committedSGaugeVersion ||
        context.integerComponentVersion != owner->second.integerComponentVersion)
        throw std::invalid_argument("STALE_EPOCH_COMPONENT_CONTEXT_SELECTION");
    if (context.committedSBasisIdentity.empty() ||
        context.committedSBasisIdentity != committedSBasisIdentity(context))
        throw std::invalid_argument("COMPONENT_CONTEXT_COMPLETE_IDENTITY_MISMATCH");
    if (singleSBasisNative() &&
        !zhangSameCommittedGraphBasis(context.basis, context.productBasis))
        throw std::invalid_argument("SINGLE_S_COMPONENT_HAS_INDEPENDENT_PRODUCT_TREE");
    const auto& committed = owner->second.basis;
    for (const auto& edge : context.basis.edges)
        if (edge.satellite.sys != system || !committed.edges.contains(edge) ||
            context.basis.treeEdges.contains(edge) != committed.treeEdges.contains(edge))
            throw std::invalid_argument("COMPONENT_CONTEXT_SYSTEM_MISMATCH");
    for (const auto& [edge, signals] : context.signalArcVersions)
    {
        const auto current = owner->second.phaseArcVersions.find(edge);
        if (current == owner->second.phaseArcVersions.end() || current->second != signals)
            throw std::invalid_argument("COMPONENT_CONTEXT_PHYSICAL_SIGNAL_IDENTITY_MISMATCH");
    }
    // Validate against a freshly constructed read-only component, not merely
    // against a self-consistent identity string supplied by the caller.
    const auto available = zhangGraphIntegerComponentContexts(state, system);
    const bool exactView = std::any_of(available.begin(), available.end(),
        [&](const auto& expected)
        { return expected.committedSBasisIdentity == context.committedSBasisIdentity; });
    if (!exactView)
        throw std::invalid_argument("COMPONENT_CONTEXT_NOT_A_COMPLETE_COMMITTED_COMPONENT_VIEW");
    auto snapshot = std::make_shared<const ZhangGraphIntegerContext>(context);
    std::lock_guard<std::mutex> lock(selectedGraphContextMutex);
    token_ = ++nextSelectedGraphContextToken;
    selectedGraphContexts[{runtime, system}].push_back({token_, std::move(snapshot)});
}

ZhangGraphIntegerContextScope::~ZhangGraphIntegerContextScope()
{
    std::lock_guard<std::mutex> lock(selectedGraphContextMutex);
    for (auto it = selectedGraphContexts.begin(); it != selectedGraphContexts.end();)
    {
        std::erase_if(it->second, [&](const auto& entry) { return entry.token == token_; });
        it = it->second.empty() ? selectedGraphContexts.erase(it) : std::next(it);
    }
}

void traceZhangCurrentGraphProductProjection(
    Trace& trace,
    const KFState& state,
    const char* stage)
{
    if (!zhangR51SingleBasisHandoffEnabled()
        || !acsConfig.zhangFullRank.enable) return;
    for (const auto& [system, options] : acsConfig.zhangFullRank.sysOpts)
        if (options.use_spanning_tree && acsConfig.process_sys[system])
            traceZhangGraphProductProjection(
                trace, state, system, options.baseline_observables, stage);
}

bool applyZhangGraphRetirementForAudit(
    Trace& trace, KFState& state, E_Sys system,
    const vector<E_ObsCode>& codes, const ZhangGraphBasis& oldBasis,
    const ZhangGraphBasis& newBasis, const set<ZhangGraphEdge>& retired,
    string& failureReason)
{
    ZhangStateTransformAudit audit;
    const bool accepted = transformZhangGraphBasis(
        trace, state, system, codes, oldBasis, newBasis, false, retired, &audit);
    failureReason = audit.failureReason;
    return accepted;
}

bool applyZhangGraphBasisTransformForAudit(
    Trace&                        trace,
    KFState&                      branch,
    E_Sys                         system,
    const vector<E_ObsCode>&      baselineObservables,
    const ZhangGraphBasis&        oldBasis,
    const ZhangGraphBasis&        newBasis,
    SparseMatrix<double>&         transform,
    string&                       failureReason)
{
    transform.resize(0, 0);
    int callbackCount = 0;
    branch.exactStateTransformCallback =
        [&](const KFState&,
            GTime,
            const map<KFKey, int>&,
            const map<KFKey, int>&,
            const SparseMatrix<double>& applied,
            const string&,
			const VectorXd&,
			const MatrixXd&,
			std::uint64_t,
			std::uint64_t)
        {
            transform = applied;
            callbackCount++;
            return true;
        };
    if (!transformZhangGraphBasis(
            trace,
            branch,
            system,
            baselineObservables,
            oldBasis,
            newBasis,
            false))
    {
        failureReason = "E29_A2_EXACT_S_TRANSFORM_FAILED";
        return false;
    }
    if (callbackCount != 1
        || transform.rows() != branch.x.size()
        || transform.cols() == 0
        || !transform.isCompressed())
    {
        failureReason = "E29_A2_TRANSFORM_CALLBACK_MISMATCH";
        return false;
    }
    failureReason = "NONE";
    return true;
}

void updateZhangFullRankReferences(
    Trace&       trace,
    ReceiverMap& receiverMap,
    KFState&     kfState
)
{
    if (!acsConfig.zhangFullRank.enable)
    {
        return;
    }

    const string runtimeId = zhangGraphRuntimeId(kfState);
    if (runtimeId.empty())
    {
        BOOST_LOG_TRIVIAL(error) << "ZHANG_GRAPH_RUNTIME_ID_UNBOUND";
        return;
    }

    // Reject live mode changes BEFORE the epoch transaction, dependency sweep
    // or auxiliary retirement can mutate authoritative FLOAT metadata/state.
    const string mode = configuredProductBasisMode();
    for (const auto& [identity, runtime] : graphStateMap)
        if (identity.first == runtimeId && runtime.initialized &&
            runtime.productBasisMode != mode)
            throw std::runtime_error("ZHANG_PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION");

    // EPOCH_ONLY_AR may not expose a previous epoch's accepted integers to
    // this graph event.  Graph update precedes the later AR search, so the
    // idempotent begin must run here rather than only in ppp_ambres.cpp.
    zhangBeginEpochOnlyIntegerTransaction(kfState);
    // Refresh identities from complete active public functionals before any
    // old chord/tree coordinate can disappear. No accepted integer RHS is
    // imported by this hook.
    {
        ZhangPhaseTimer directoryTimer(trace,"PRE_GRAPH_FLOAT_TARGET_DIRECTORY_TOTAL",kfState.time.to_string(0));
        zhangRegisterPublishedProductHandoffDependencies(kfState);
    }

    // All systems share one FLOAT registry. Sweep before availability can
    // skip a system, and before any system-local graph event is evaluated.
    if (!retireZhangGraphProductEventAuxiliaries(trace, kfState, "PRE_GRAPH"))
    {
        trace << "\nZHANG_HANDOFF_FLOAT_PRE_GRAPH_SWEEP_FAILED time="
              << kfState.time.to_string(0) << " reason="
              << kfState.lastFactorTransactionFailureReason;
        trace.flush();
        throw std::runtime_error("ZHANG_HANDOFF_FLOAT_PRE_GRAPH_SWEEP_FAILED:" +
            kfState.lastFactorTransactionFailureReason);
    }

    for (auto& [sys, options] : acsConfig.zhangFullRank.sysOpts)
    {
        if (!acsConfig.process_sys[sys])
        {
            continue;
        }

        ReferenceAvailability availability =
            referenceAvailability(receiverMap, sys, options.baseline_observables);
        if (availability.satellitesByReceiver.empty())
        {
            continue;
        }

        if (options.use_spanning_tree)
        {
            traceZhangGraphProductProjection(trace, kfState, sys,
                options.baseline_observables, "PRE_GRAPH");
            updateZhangGraphBasis(trace, kfState, sys, options, availability);
            traceZhangGraphProductProjection(trace, kfState, sys,
                options.baseline_observables, "POST_GRAPH");
            continue;
        }

        if (!options.auto_reference_switch)
        {
            continue;
        }

        set<SatSys> common = commonSatellites(availability);
        if (common.empty())
        {
            BOOST_LOG_TRIVIAL(warning)
                << "ZHANG_REFERENCE_SWITCH sys=" << enum_to_string(sys)
                << " skipped: no baseline satellite is common to all active receivers";
            continue;
        }

        auto& outage = outageStateMap[{runtimeId, sys}];

        bool receiverAvailable =
            availability.satellitesByReceiver.find(options.reference_receiver) !=
            availability.satellitesByReceiver.end();

        SatSys oldSatellite(options.reference_satellite.c_str());
        bool satelliteAvailable = common.find(oldSatellite) != common.end();

        outage.receiverEpochs = receiverAvailable ? 0 : outage.receiverEpochs + 1;
        outage.satelliteEpochs = satelliteAvailable ? 0 : outage.satelliteEpochs + 1;

        bool changeReceiver = outage.receiverEpochs >= options.reference_outage_epochs;
        bool changeSatellite = outage.satelliteEpochs >= options.reference_outage_epochs;
        if (!changeReceiver && !changeSatellite)
        {
            continue;
        }

        vector<string> receiverChoices = {options.reference_receiver};
        if (changeReceiver)
        {
            receiverChoices =
                orderedReceivers(availability, options.reference_receiver_candidates);
        }

        vector<SatSys> satelliteChoices = {oldSatellite};
        if (changeSatellite)
        {
            satelliteChoices =
                orderedSatellites(
                    common,
                    availability.elevationScore,
                    options.reference_satellite_candidates
                );
        }

        if (receiverChoices.empty() || satelliteChoices.empty())
        {
            BOOST_LOG_TRIVIAL(warning)
                << "ZHANG_REFERENCE_SWITCH sys=" << enum_to_string(sys)
                << " skipped: no valid replacement reference";
            continue;
        }

        string oldReceiver = options.reference_receiver;

        bool hasEstimatedZhangState = false;
        for (auto& [key, index] : kfState.kfIndexMap)
        {
            if (key.type == KF::AMBIGUITY && key.Sat.sys == sys)
            {
                hasEstimatedZhangState = true;
                break;
            }
        }

        string newReceiver;
        SatSys newSatellite;
        bool   transformed = !hasEstimatedZhangState;

        for (auto& receiverCandidate : receiverChoices)
        {
            for (auto& satelliteCandidate : satelliteChoices)
            {
                if (hasEstimatedZhangState &&
                    !transformZhangDatum(
                        trace,
                        kfState,
                        sys,
                        options.baseline_observables,
                        oldReceiver,
                        oldSatellite,
                        receiverCandidate,
                        satelliteCandidate
                    ))
                {
                    continue;
                }

                newReceiver  = receiverCandidate;
                newSatellite = satelliteCandidate;
                transformed  = hasEstimatedZhangState;
                break;
            }

            if (!newReceiver.empty())
            {
                break;
            }
        }

        if (newReceiver.empty() || newSatellite.prn <= 0)
        {
            BOOST_LOG_TRIVIAL(warning)
                << "ZHANG_REFERENCE_SWITCH sys=" << enum_to_string(sys)
                << " deferred: no candidate has the complete phase-state tree required for an "
                   "exact transform; retaining "
                << oldReceiver << "/" << oldSatellite.id();
            continue;
        }

        options.reference_receiver  = newReceiver;
        options.reference_satellite = newSatellite.id();
        outage = {};

        BOOST_LOG_TRIVIAL(info)
            << "ZHANG_REFERENCE_SWITCH sys=" << enum_to_string(sys)
            << " old_receiver=" << oldReceiver
            << " new_receiver=" << newReceiver
            << " old_satellite=" << oldSatellite.id()
            << " new_satellite=" << newSatellite.id()
            << " transformed=" << transformed;
    }
    // One GLOBAL target-space compression after all system-local graph
    // transactions/receipt bindings, before any new FLOAT factor is built.
    // This does not merge the preceding per-system graph transforms or erase
    // their transient capture peaks.
    string compressionReason;
    if (!zhangCompressGraphProductFloatTargets(trace, kfState, compressionReason))
        throw std::runtime_error("ZHANG_FLOAT_TARGET_COMPRESSION_METADATA_INVALID:"+compressionReason);
}

bool zhangGraphOwnsAmbiguityLifecycle(E_Sys system, E_ObsCode code)
{
    if (!acsConfig.zhangFullRank.enable) return false;
    const auto options = acsConfig.zhangFullRank.sysOpts.find(system);
    return options != acsConfig.zhangFullRank.sysOpts.end() &&
        options->second.use_spanning_tree &&
        zhangFullRankUsesObservable(code, options->second.baseline_observables);
}

bool zhangGraphModelsObservation(
    const KFState&     kfState,
    const std::string& receiver,
    const SatSys&      satellite,
    E_ObsCode          code
)
{
    if (!acsConfig.zhangFullRank.enable)
    {
        return true;
    }

    auto optionsIt = acsConfig.zhangFullRank.sysOpts.find(satellite.sys);
    if (optionsIt == acsConfig.zhangFullRank.sysOpts.end() ||
        !zhangFullRankUsesObservable(code, optionsIt->second.baseline_observables) ||
        !optionsIt->second.use_spanning_tree)
    {
        return true;
    }

    const string runtimeId = zhangGraphRuntimeId(kfState);
    auto stateIt = graphStateMap.find({runtimeId, satellite.sys});
    if (stateIt == graphStateMap.end() || !stateIt->second.initialized)
    {
        return false;
    }

    const auto& runtime = stateIt->second;
    if (runtime.productBasisMode != configuredProductBasisMode()) return false;
    const ZhangGraphEdge edge{receiver, satellite};
    if (runtime.phaseAdmissionTime != kfState.time.to_string(9) ||
        runtime.phaseAdmissionChart != runtime.representationVersion ||
        !runtime.observationEdges.contains(edge)) return false;
    const auto token = runtime.phaseAdmissions.find(edge);
    const auto arc = runtime.phaseArcVersions.find(edge);
    return token != runtime.phaseAdmissions.end() && arc != runtime.phaseArcVersions.end() &&
        token->second.contains(int(code)) && arc->second.contains(int(code)) &&
        token->second.at(int(code)) == arc->second.at(int(code));
}

bool zhangGraphRetainsAmbiguity(
    const KFState&     kfState,
    const std::string& receiver,
    const SatSys&      satellite,
    E_ObsCode          code
)
{
    auto optionsIt = acsConfig.zhangFullRank.sysOpts.find(satellite.sys);
    if (optionsIt == acsConfig.zhangFullRank.sysOpts.end() ||
        !optionsIt->second.use_spanning_tree)
    {
        return true;
    }

    const string runtimeId = zhangGraphRuntimeId(kfState);
    auto stateIt = graphStateMap.find({runtimeId, satellite.sys});
    if (stateIt == graphStateMap.end() || !stateIt->second.initialized)
    {
        return false;
    }

    ZhangGraphEdge edge{receiver, satellite};
    if (stateIt->second.productBasisMode != configuredProductBasisMode()) return false;
    const auto selected = selectedGraphContext(kfState, satellite.sys);
    const auto& basis = selected ? selected->basis : stateIt->second.basis;
    const auto admission = stateIt->second.phaseAdmissions.find(edge);
    return stateIt->second.observationEdges.contains(edge) &&
           admission != stateIt->second.phaseAdmissions.end() &&
           admission->second.contains(int(code)) &&
           basis.edges.contains(edge) &&
           !basis.treeEdges.contains(edge);
}

std::string zhangGraphAmbiguityAdmissionReason(
    const KFState& state, const std::string& receiver,
    const SatSys& satellite, E_ObsCode code)
{
    const auto owner = graphStateMap.find({zhangGraphRuntimeId(state), satellite.sys});
    if (owner == graphStateMap.end() || !owner->second.initialized)
        return "GRAPH_RUNTIME_NOT_INITIALIZED";
    const auto& runtime = owner->second;
    if (runtime.productBasisMode != configuredProductBasisMode())
        return "PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION";
    const auto selected = selectedGraphContext(state, satellite.sys);
    const auto& basis = selected ? selected->basis : runtime.basis;
    const ZhangGraphEdge edge{receiver, satellite};
    if (!basis.edges.contains(edge)) return selected
        ? "OUTSIDE_SELECTED_INTEGER_COMPONENT" : "NOT_IN_COMMITTED_PHYSICAL_BASIS";
    if (basis.treeEdges.contains(edge)) return "S_TREE_EDGE_HAS_NO_CYCLE_COORDINATE";
    if (!runtime.observationEdges.contains(edge)) return "NO_CURRENT_OBSERVATION_SUPPORT";
    const auto admission = runtime.phaseAdmissions.find(edge);
    if (admission == runtime.phaseAdmissions.end() || !admission->second.contains(int(code)))
        return "CURRENT_SIGNAL_QC_ADMISSION_MISSING";
    const auto physical = runtime.phaseArcVersions.find(edge);
    if (physical == runtime.phaseArcVersions.end() || !physical->second.contains(int(code)))
        return "CURRENT_SIGNAL_ARC_VERSION_UNKNOWN";
    if (admission->second.at(int(code)) != physical->second.at(int(code)))
        return "CURRENT_SIGNAL_ARC_VERSION_MISMATCH";
    return "DIRECT_CURRENT_CYCLE_ELIGIBLE";
}

bool zhangGraphProductSatelliteActive(
    const KFState& kfState,
    const SatSys&  satellite
)
{
    auto optionsIt = acsConfig.zhangFullRank.sysOpts.find(satellite.sys);
    if (optionsIt == acsConfig.zhangFullRank.sysOpts.end() ||
        !optionsIt->second.use_spanning_tree)
    {
        return true;
    }

    if (const auto selected = selectedGraphContext(kfState, satellite.sys))
        return selected->initialized && selected->productBasis.connected &&
            selected->productBasis.receivers.contains(selected->productBasis.rootReceiver) &&
            selected->productBasis.satellites.contains(satellite);
    const string runtimeId = zhangGraphRuntimeId(kfState);
    auto stateIt = graphStateMap.find({runtimeId, satellite.sys});
    if (stateIt == graphStateMap.end() || !stateIt->second.initialized || !stateIt->second.productInitialized ||
        !stateIt->second.productBasis.connected ||
        !stateIt->second.productBasis.receivers.contains(stateIt->second.productBasis.rootReceiver) ||
        stateIt->second.productBasisMode != configuredProductBasisMode() ||
        stateIt->second.productAuthorityComponentVersion != stateIt->second.integerComponentVersion)
    {
        return false;
    }

    for (const auto& edge : stateIt->second.productBasis.edges)
    {
        if (edge.satellite == satellite)
        {
            return true;
        }
    }
    return false;
}

ZhangP0GraphStamp zhangGraphIntegerContextStamp(
	const KFState& state,
	E_Sys system)
{
	ZhangP0GraphStamp stamp;
	stamp.runtimeIdentity = zhangGraphRuntimeId(state);
	stamp.system = static_cast<int>(system);
	std::ostringstream epoch;
	epoch << std::hexfloat << state.time.bigTime;
	stamp.epochIdentity = epoch.str();

    if (const auto selected = selectedGraphContext(state, system))
    {
        stamp.runtimeIdentity += "|INTEGER_COMPONENT=" + selected->componentId;
        stamp.event = selected->eventId;
        stamp.representation = selected->representationVersion;
        stamp.productDatum = selected->productDatumVersion;
        stamp.floatGauge = selected->floatGaugeVersion;
        stamp.integerComponent = selected->integerComponentVersion;
        stamp.available = selected->initialized && selected->basis.connected;
        return stamp;
    }

	auto found = graphStateMap.find({stamp.runtimeIdentity, system});
	if (found == graphStateMap.end())
	{
		return stamp;
	}

	const auto& runtime = found->second;
    if (runtime.productBasisMode != configuredProductBasisMode()) return stamp;
	stamp.event = runtime.eventCounter;
	stamp.representation = runtime.representationVersion;
	stamp.productDatum = runtime.productDatumVersion;
	stamp.floatGauge = runtime.floatGaugeVersion;
	stamp.integerComponent = runtime.integerComponentVersion;
	stamp.runtimeEpoch = runtime.epochIndex;
	const auto& basis = !singleSBasisNative() && runtime.activeBasis.connected
		? runtime.activeBasis
		: runtime.basis;
	stamp.available = runtime.initialized
		&& basis.connected
		&& basis.receivers.contains(basis.rootReceiver);
	return stamp;
}

bool zhangGraphIntegerContext(
    const KFState&             kfState,
    E_Sys                      system,
    ZhangGraphIntegerContext& context
)
{
    context = {};
    if (const auto selected = selectedGraphContext(kfState, system))
    {
        context = *selected;
        return context.initialized;
    }
    const string runtimeId = zhangGraphRuntimeId(kfState);
    auto stateIt = graphStateMap.find({runtimeId, system});
    if (stateIt == graphStateMap.end() || !stateIt->second.initialized)
    {
        return false;
    }

    if (stateIt->second.productBasisMode != configuredProductBasisMode()) return false;
    context.productBasisMode = stateIt->second.productBasisMode;
    context.runtimeIdentity = runtimeId;
    context.snapshotEpoch = kfState.time.to_string(9);
    context.basis   = !singleSBasisNative() && stateIt->second.activeBasis.connected
        ? stateIt->second.activeBasis
        : stateIt->second.basis;
    // FLOAT forests retain component-local information and uncertain offsets.
    // This API serves the existing single-datum integer contract; a forest is
    // not an authorization to search or publish cross-component integers.
    if (!context.basis.connected || !context.basis.receivers.contains(context.basis.rootReceiver)) return false;
    context.productBasis = singleSBasisNative() ? context.basis : stateIt->second.productBasis;
    context.eventId = stateIt->second.eventCounter;
    context.productDatumVersion = stateIt->second.productDatumVersion;
    context.committedSGaugeVersion = stateIt->second.committedSGaugeVersion;
    context.representationVersion =
        stateIt->second.representationVersion;
    context.floatGaugeVersion = stateIt->second.floatGaugeVersion;
    context.integerComponentVersion =
        stateIt->second.integerComponentVersion;
    for (const auto& [edge, history] : stateIt->second.edgeHistory)
    {
        context.arcVersions[edge] = history.arcVersion;
    }
    context.signalArcVersions = stateIt->second.phaseArcVersions;
    context.signalObservationEdges =
        stateIt->second.signalObservationEdges;
    context.signalObservationFinalized =
        stateIt->second.signalObservationFinalized;
    if (zhangR51ProductRepairEnabled() &&
        !context.signalObservationFinalized)
        context.signalObservationEdges.clear();
    context.initialized = true;
    context.committedSBasisIdentity = committedSBasisIdentity(context);
    return true;
}

std::optional<int> zhangGraphSignalArcVersion(
    const ZhangGraphIntegerContext& context,
    const ZhangGraphEdge& edge,
    E_ObsCode observable)
{
    if (!context.initialized || observable == E_ObsCode::NONE)
        return std::nullopt;
    const auto physical = context.signalArcVersions.find(edge);
    if (physical == context.signalArcVersions.end())
        return std::nullopt;
    const auto signal = physical->second.find(int(observable));
    if (signal == physical->second.end() || signal->second < 0)
        return std::nullopt;
    return signal->second;
}

bool zhangGraphSignalArcVersions(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    map<ZhangGraphEdge, int>& versions,
    string& failureReason)
{
    versions.clear();
    failureReason = "NOT_EVALUATED";
    ZhangGraphIntegerContext context;
    if (!zhangGraphIntegerContext(state, system, context))
    {
        failureReason = "SIGNAL_GRAPH_CONTEXT_UNAVAILABLE";
        return false;
    }
    // This is a known physical-identity catalogue, not a QC or current-row
    // admission query. FLOAT cycle columns expand through context.basis;
    // products use context.productBasis. Both must see every already-known
    // version in their union, while a genuinely unknown signal stays absent
    // (never V0). Final phase-row support is tracked separately in
    // signalObservationEdges and checked by current-evidence admission.
    set<ZhangGraphEdge> cataloguedEdges = context.basis.edges;
    cataloguedEdges.insert(context.productBasis.edges.begin(),
        context.productBasis.edges.end());
    for (const auto& edge : cataloguedEdges)
    {
        const auto version = zhangGraphSignalArcVersion(
            context, edge, observable);
        if (version) versions[edge] = *version;
    }
    if (versions.empty())
    {
        failureReason = "NO_KNOWN_SIGNAL_PHYSICAL_ARCS";
        return false;
    }
    failureReason = "KNOWN_SIGNAL_PHYSICAL_ARC_CATALOGUE";
    return true;
}

bool zhangFinalizeGraphSignalObservationAdmission(
    const KFState& state,
    E_Sys system,
    const KFMeas& finalMeasurement,
    string& failureReason)
{
    failureReason = "NOT_EVALUATED";
    const auto owner = graphStateMap.find(
        {zhangGraphRuntimeId(state), system});
    if (owner == graphStateMap.end() || !owner->second.initialized)
    {
        failureReason = "SIGNAL_GRAPH_CONTEXT_UNAVAILABLE";
        return false;
    }
    auto& runtime = owner->second;
    if (runtime.productBasisMode != configuredProductBasisMode())
    {
        failureReason = "PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION";
        return false;
    }
    runtime.signalObservationFinalized = false;
    runtime.signalObservationFactorSequence = 0;
    const bool factorValid = state.lastPppTransactionCommitted &&
        state.lastPppQcConverged &&
        finalMeasurement.time.to_string(9) == state.time.to_string(9) &&
        runtime.phaseAdmissionTime == state.time.to_string(9) &&
        runtime.phaseAdmissionChart == runtime.representationVersion &&
        finalMeasurement.obsKeys.size() ==
            size_t(finalMeasurement.H.rows()) &&
        finalMeasurement.metaDataMaps.size() ==
            finalMeasurement.obsKeys.size() &&
        finalMeasurement.R.rows() == finalMeasurement.H.rows() &&
        finalMeasurement.R.cols() == finalMeasurement.H.rows();
    map<int, set<ZhangGraphEdge>> retained = factorValid
        ? zhangFinalAdmittedSignalFactorRows(
            system, finalMeasurement, runtime.signalObservationEdges,
            runtime.phaseAdmissions)
        : map<int, set<ZhangGraphEdge>>{};
    // A final KF row can still reference a formerly modelled FLOAT edge that
    // was removed from the active graph component in this transaction. Keep
    // the raw factor in the KF; do not promote it into graph integer support.
    retained = zhangGraphSupportedSignalEdges(
        retained, runtime.stateEdges, runtime.basis.edges);

    auto eraseEpoch = [&](GraphRuntimeState::ArcVersionObservationHistory& history)
    {
        std::erase(history.observedEpochs, runtime.epochIndex);
        history.observationEpochs = int(history.observedEpochs.size());
        history.firstObservedEpoch = history.observedEpochs.empty()
            ? -1 : history.observedEpochs.front();
        history.lastObservedEpoch = history.observedEpochs.empty()
            ? -1 : history.observedEpochs.back();
    };
    for (const auto& [code, edges] : runtime.signalObservationEdges)
        for (const auto& edge : edges)
            if (!retained[code].contains(edge))
            {
                const auto version = runtime.phaseArcVersions.find(edge);
                if (version == runtime.phaseArcVersions.end() ||
                    !version->second.contains(code)) continue;
                const auto history =
                    runtime.signalArcVersionObservationHistory.find(
                        {edge, code, version->second.at(code)});
                if (history !=
                    runtime.signalArcVersionObservationHistory.end())
                    eraseEpoch(history->second);
            }
    set<ZhangGraphEdge> retainedUnion;
    for (const auto& [code, edges] : retained)
        retainedUnion.insert(edges.begin(), edges.end());
    for (const auto& edge : runtime.observationEdges)
        if (!retainedUnion.contains(edge))
        {
            const auto topology = runtime.edgeHistory.find(edge);
            if (topology == runtime.edgeHistory.end()) continue;
            const auto history = runtime.arcVersionObservationHistory.find(
                {edge, topology->second.arcVersion});
            if (history != runtime.arcVersionObservationHistory.end())
                eraseEpoch(history->second);
        }
    for (auto it = runtime.phaseAdmissions.begin();
         it != runtime.phaseAdmissions.end();)
    {
        for (auto signal = it->second.begin();
             signal != it->second.end();)
            signal = retained[signal->first].contains(it->first)
                ? std::next(signal) : it->second.erase(signal);
        it = it->second.empty()
            ? runtime.phaseAdmissions.erase(it) : std::next(it);
    }
    runtime.signalObservationEdges = std::move(retained);
    runtime.observationEdges = std::move(retainedUnion);
    if (!factorValid)
    {
        failureReason = "FINAL_PPP_FACTOR_OR_QC_NOT_PROVEN";
        return false;
    }
    runtime.signalObservationFactorSequence = state.factorCommitSequence;
    runtime.signalObservationFinalized = true;
    if (singleSBasisNative())
    {
        // Final QC changes current rows, not the only committed S chart.
        syncSingleSBasisProductView(runtime);
    }
    else if (runtime.productInitialized)
    {
        // updateProductDatum runs before measurement construction.  A
        // candidate edge admitted there can still fail to enter the actual
        // committed S chart/final factor.  Revoke only that now-invalid
        // preselection; the current AR scope may select a new legal target
        // from the finalized represented component without mutating FLOAT.
        const auto& committed = runtime.activeBasis.edges.empty()
            ? runtime.basis : runtime.activeBasis;
        set<ZhangGraphEdge> representedTreeEdges;
        map<ZhangGraphEdge, int> currentArcs;
        set<std::pair<ZhangGraphEdge, int>> knownArcHistories;
        vector<int> requiredSignals;
        const auto configured = acsConfig.zhangFullRank.sysOpts.find(system);
        if (configured != acsConfig.zhangFullRank.sysOpts.end())
            for (const auto code : configured->second.baseline_observables)
                requiredSignals.push_back(int(code));
        for (const auto& edge : runtime.productBasis.treeEdges)
        {
            bool represented = committed.edges.contains(edge);
            for (const auto signal : requiredSignals)
            {
                const auto code = static_cast<E_ObsCode>(signal);
                represented &= state.kfIndexMap.contains(
                    zhangTransportSatelliteKey(code, edge.satellite));
                // The actual forest's non-global root potential is a
                // stochastic state, not a zero imposed by a local AR view.
                if (edge.receiver != committed.rootReceiver)
                    represented &= state.kfIndexMap.contains(
                        zhangTransportReceiverKey(system, code, edge.receiver));
                if (!committed.treeEdges.contains(edge))
                    represented &= state.kfIndexMap.contains(
                        zhangTransportAmbiguityKey(code, edge));
            }
            if (represented) representedTreeEdges.insert(edge);
            const auto current = runtime.edgeHistory.find(edge);
            if (current != runtime.edgeHistory.end()) currentArcs[edge] = current->second.arcVersion;
            const auto selected = runtime.productArcVersions.find(edge);
            if (selected != runtime.productArcVersions.end() &&
                runtime.arcVersionObservationHistory.contains({edge, selected->second}))
                knownArcHistories.insert({edge, selected->second});
        }
        const auto productAudit = zhangAuditPreselectedProductTarget(
            committed, runtime.productBasis, representedTreeEdges,
            runtime.productArcVersions, currentArcs, knownArcHistories,
            runtime.productSignalArcVersions, runtime.phaseArcVersions,
            requiredSignals);
        if (!productAudit.valid)
        {
            const auto selected = runtime.productArcVersions.find(productAudit.edge);
            const auto current = currentArcs.find(productAudit.edge);
            const int selectedArc = selected == runtime.productArcVersions.end() ? -1 : selected->second;
            BOOST_LOG_TRIVIAL(warning)
                << "ZHANG_FINAL_FACTOR_PRODUCT_SELECTOR_REVOKED time=" << state.time.to_string(0)
                << " system=" << enum_to_string(system)
                << " reason=" << productAudit.reason
                << " edge=" << (productAudit.hasEdge
                    ? productAudit.edge.receiver + ":" + productAudit.edge.satellite.id() : "NONE")
                << " selected_arc=" << selectedArc
                << " current_arc=" << (current == currentArcs.end() ? -1 : current->second)
                << " selected_arc_history_present=" << knownArcHistories.contains({productAudit.edge, selectedArc})
                << " product_satellites=" << runtime.productBasis.satellites.size()
                << " committed_satellites=" << committed.satellites.size()
                << " float_state_mutated=0 covariance_mutated=0 physical_arcs_mutated=0"
                   " observation_history_mutated=0 current_ar_reselection_allowed=1";
            invalidateZhangProduct(runtime, "FINAL_FACTOR_PRESELECTED_PRODUCT_INVALID:" + productAudit.reason);
        }
    }
    failureReason = "FINAL_PPP_PHASE_FACTOR_ROWS_BOUND";
    return true;
}

bool zhangGraphRequiresFreshProductIdentity(const KFState& state)
{
    const string runtimeId = zhangGraphRuntimeId(state);
    for (const auto& [identity, runtime] : graphStateMap)
        if (identity.first == runtimeId &&
            runtime.requiresFreshProductIdentity)
            return true;
    return false;
}

bool zhangGraphProductEventReceipt(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    const SatSys& satellite,
    ZhangGraphProductEventReceipt& receipt)
{
    if (!zhangR51SingleBasisHandoffEnabled()) return false;
    const auto found = graphProductEventReceipts.find({
        zhangGraphRuntimeId(state), state.time.to_string(9),
        system, observable, satellite});
    if (found == graphProductEventReceipts.end()) return false;
    receipt = found->second;
    return true;
}

bool zhangGraphPhysicalArcFloatRow(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    const ZhangGraphEdge& edge,
    int arcVersion,
    ZhangCoordinateRow& rowMetres,
    string& failureReason)
{
    rowMetres.clear();
    failureReason = "NOT_EVALUATED";
    if (!zhangR51SingleBasisHandoffEnabled())
    {
        failureReason = "SINGLE_BASIS_HANDOFF_DISABLED";
        return false;
    }
    if (edge.satellite.sys != system || arcVersion < 0)
    {
        failureReason = "INVALID_VERSIONED_PHYSICAL_EDGE";
        return false;
    }
    const auto owner = graphStateMap.find(
        {zhangGraphRuntimeId(state), system});
    if (owner == graphStateMap.end() || !owner->second.initialized)
    {
        failureReason = "CURRENT_GRAPH_UNAVAILABLE";
        return false;
    }
    const auto& runtime = owner->second;
    if (runtime.productBasisMode != configuredProductBasisMode())
    {
        failureReason = "PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION";
        return false;
    }
    if (zhangR51ProductRepairEnabled() &&
        (!runtime.signalObservationFinalized ||
        runtime.phaseAdmissionTime != state.time.to_string(9) ||
        runtime.signalObservationFactorSequence == 0))
    {
        failureReason = "CURRENT_SIGNAL_FACTOR_ROWS_NOT_FINALIZED";
        return false;
    }
    const auto& basis = !singleSBasisNative() && runtime.activeBasis.connected
        ? runtime.activeBasis : runtime.basis;
    const auto identity = runtime.phaseArcVersions.find(edge);
    const auto signal = identity == runtime.phaseArcVersions.end()
        ? map<int, int>::const_iterator{}
        : identity->second.find(int(observable));
    if (identity != runtime.phaseArcVersions.end()
        && signal != identity->second.end()
        && signal->second == arcVersion
        && runtime.observationEdges.contains(edge)
        && basis.edges.contains(edge))
    {
        // An edge-level observation flag is insufficient for a dual-signal
        // integer decision.  The exact signal/version must have been admitted
        // to this epoch's phase model, not merely be present in the graph.
        const int signal = static_cast<int>(observable);
        const auto admission = runtime.phaseAdmissions.find(edge);
        const auto serial = runtime.phaseArcVersions.find(edge);
        if (runtime.phaseAdmissionTime != state.time.to_string(9) ||
            runtime.phaseAdmissionChart != runtime.representationVersion ||
            admission == runtime.phaseAdmissions.end() ||
            serial == runtime.phaseArcVersions.end() ||
            !admission->second.contains(signal) ||
            !serial->second.contains(signal) ||
            admission->second.at(signal) != serial->second.at(signal))
        {
            failureReason = "CURRENT_PHYSICAL_SIGNAL_NOT_ADMITTED";
            return false;
        }
        const double wavelength = zhangWavelength(system, observable);
        if (!(wavelength > 0))
        {
            failureReason = "PRODUCT_WAVELENGTH_UNAVAILABLE";
            return false;
        }
        if (edge.receiver != basis.rootReceiver)
            rowMetres[zhangTransportReceiverKey(
                system, observable, edge.receiver)] = 1;
        rowMetres[zhangTransportSatelliteKey(
            observable, edge.satellite)] = 1;
        if (!basis.treeEdges.contains(edge))
            rowMetres[zhangTransportAmbiguityKey(
                observable, edge)] = wavelength;
        for (const auto& [key, coefficient] : rowMetres)
            if (state.kfIndexMap.find(key) == state.kfIndexMap.end())
            {
                rowMetres.clear();
                failureReason = "CURRENT_PHYSICAL_ROW_STATE_MISSING";
                return false;
            }
        failureReason = "CURRENT_VERSIONED_ARC_ROW";
        return true;
    }

    KFKey retained;
    if (zhangHandoffFloatFind(state,
            {edge, int(observable), arcVersion}, retained))
    {
        // This is a stochastic coordinate only. The current-domain exact
        // integer proof and current observation gates remain at the caller.
        rowMetres[retained] = 1;
        failureReason = "CORRELATED_RETIRED_ARC_AUXILIARY_ROW";
        return true;
    }
    failureReason = "VERSIONED_PHYSICAL_ARC_NOT_IN_CURRENT_POSTERIOR";
    return false;
}

bool zhangGraphProductEventTransitionWitness(
    const KFState& state,
    E_Sys system,
    E_ObsCode observable,
    const SatSys& satellite,
    const SatSys& referenceSatellite,
    const ZhangProductIntegerTransition& transition,
    string& failureReason)
{
    failureReason = "NOT_EVALUATED";
    const auto owner = graphStateMap.find(
        {zhangGraphRuntimeId(state), system});
    if (owner == graphStateMap.end() ||
        (zhangR51ProductRepairEnabled() &&
         (!owner->second.signalObservationFinalized ||
          owner->second.phaseAdmissionTime != state.time.to_string(9) ||
          owner->second.signalObservationFactorSequence == 0)))
    {
        failureReason = "CURRENT_SIGNAL_FACTOR_ROWS_NOT_FINALIZED";
        return false;
    }
    if (satellite.sys != system
        || referenceSatellite.sys != system)
    {
        failureReason = "INVALID_RELATIVE_PRODUCT_ENDPOINTS";
        return false;
    }
    if (!transition.valid
        || transition.physicalEdges.size() != transition.coefficients.size()
        || transition.physicalEdges.size()
            != transition.physicalArcVersions.size())
    {
        failureReason = "INVALID_VERSIONED_INTEGER_TRANSITION";
        return false;
    }
    if (transition.affineOffsetCycles != 0)
    {
        failureReason = "NONZERO_AFFINE_OFFSET_NOT_MODELED";
        return false;
    }
    ZhangGraphProductEventReceipt subject;
    ZhangGraphProductEventReceipt reference;
    if (!zhangGraphProductEventReceipt(
            state, system, observable, satellite, subject)
        || (satellite != referenceSatellite &&
            !zhangGraphProductEventReceipt(
                state, system, observable,
                referenceSatellite, reference)))
    {
        failureReason = "SAME_EPOCH_GRAPH_RECEIPT_MISSING";
        return false;
    }
    auto rowHasMapProof = [](const ZhangGraphProductEventReceipt& row)
    {
        return row.conditionalProductMapProven
            || (row.sameTimeMapProven
                && row.bridgePhysicalEdges.empty());
    };
    if (!rowHasMapProof(subject) ||
        (satellite != referenceSatellite && !rowHasMapProof(reference)) ||
        !subject.fullJointMomentTransportProven ||
        (satellite != referenceSatellite &&
            !reference.fullJointMomentTransportProven))
    {
        failureReason = "RELATIVE_PRODUCT_CONDITIONAL_MAP_UNPROVEN";
        return false;
    }
    if (satellite != referenceSatellite &&
        (subject.time != reference.time ||
         subject.fullJointTransportId != reference.fullJointTransportId ||
         subject.oldChartGeneration != reference.oldChartGeneration ||
         subject.newChartGeneration != reference.newChartGeneration))
    {
        failureReason = "INCONSISTENT_PRODUCT_EVENT_RECEIPTS";
        return false;
    }
    using VersionedEdge = pair<ZhangGraphEdge, int>;
    map<VersionedEdge, ZhangExactInteger> expected;
    auto accumulate = [&](const ZhangGraphProductEventReceipt& receipt,
                          int sign)
    {
        for (size_t i = 0; i < receipt.bridgePhysicalEdges.size(); i++)
            expected[{receipt.bridgePhysicalEdges[i],
                receipt.bridgePhysicalArcVersions[i]}] +=
                sign * receipt.bridgePhysicalCoefficients[i];
    };
    accumulate(subject, +1);
    if (satellite != referenceSatellite)
        accumulate(reference, -1);
    map<VersionedEdge, ZhangExactInteger> submitted;
    for (size_t i = 0; i < transition.physicalEdges.size(); i++)
    {
        if (transition.physicalArcVersions[i] < 0)
        {
            failureReason = "INVALID_TRANSITION_ARC_VERSION";
            return false;
        }
        submitted[{transition.physicalEdges[i],
            transition.physicalArcVersions[i]}] +=
            transition.coefficients[i];
    }
    for (auto it = expected.begin(); it != expected.end();)
        it = it->second == 0 ? expected.erase(it) : std::next(it);
    for (auto it = submitted.begin(); it != submitted.end();)
        it = it->second == 0 ? submitted.erase(it) : std::next(it);
    if (expected.empty())
    {
        failureReason = submitted.empty()
            ? "STRUCTURAL_ZERO_BRIDGE_NOT_REQUIRED"
            : "NO_INTEGER_BRIDGE_REQUIRED";
        return submitted.empty();
    }
    if (expected != submitted)
    {
        failureReason = "VERSIONED_RELATIVE_BRIDGE_MISMATCH";
        return false;
    }
    failureReason = "EXACT_RELATIVE_GRAPH_BRIDGE_SUPPORT";
    return true;
}

bool retireZhangGraphProductEventAuxiliaries(
    Trace& trace,
    KFState& state,
    const std::string& stage)
{
    const auto reject = [&](const string& phase, const string& reason)
    {
        state.lastFactorTransactionFailureReason = reason.empty()
            ? "HANDOFF_FLOAT_RETIREMENT_REJECTED_WITHOUT_REASON" : reason;
        trace << "\nZHANG_PRODUCT_HANDOFF_AUXILIARY_RETIRE time="
              << state.time.to_string(0)
              << " scope=GLOBAL stage=" << stage
              << " status=REJECTED phase=" << phase
              << " reason=" << state.lastFactorTransactionFailureReason;
        trace.flush();
        return false;
    };
    string reason;
    if (!zhangR51SingleBasisHandoffEnabled())
    {
        // Preserve the former helper's feature-off no-op semantics. Callers
        // with an existing strict output/AR boundary keep that check there.
        return true;
    }
    // An expired but otherwise valid coordinate is a retirement candidate,
    // never an excuse to discard corrupt identity/model/moment metadata.
    if (!zhangProductFloatTargetMetadataValid(state, reason))
        return reject("ENTRY_WHOLE_TARGET_METADATA", reason);
    if (!zhangHandoffFloatRegistryViewValid(
            state, state.metaDataMap, {}, true, reason))
        return reject("ENTRY_REGISTRY", reason);
    unsigned retentionSeconds = 0;
    if (!zhangHandoffFloatRetentionSeconds(retentionSeconds, &reason))
        return reject("RETENTION_POLICY", reason);
    const auto dependencies = zhangHandoffFloatDependencies(state);
    map<KFKey, ZhangCoordinateRow> keep;
    set<string> removedRecords;
    set<KFKey> removedKeys;
    int retired = 0;
    int retained = 0;
    int retiredDependencies = 0;
    std::size_t retainedActive = 0;
    std::size_t retainedInactive = 0;
    for (const auto& [key, index] : state.kfIndexMap)
    {
        if (key.type == KF::PRODUCT_HANDOFF)
        {
            if (zhangIsProductFloatCarrierKey(key))
            {
                // Carrier lifetime belongs to complete target-space planning,
                // never to an individual physical-arc cache identity or TTL.
                if (!zhangProductFloatCarrierCoordinateValid(state, state.metaDataMap, key, reason))
                    return reject("ENTRY_CARRIER", reason);
                keep[key][key] = 1;
                continue;
            }
            ZhangHandoffFloatRecord record;
            if (!zhangHandoffFloatReadRecord(state, key, record))
                return reject("ENTRY_IDENTITY", "ZHANG_HANDOFF_FLOAT_REGISTRY_MISSING_OR_CORRUPT");
            // A declared active physical identity pins the actual KF
            // coordinate and its full cross-P, independent of cache TTL or
            // constellation order. Neither lookup nor publication refreshes
            // birth. Unregistered identities keep only the original same-
            // epoch inactive cache lifecycle and global lexical cache budget.
            const bool requested = dependencies.contains(record.identity);
            const bool live = retentionSeconds > 0 &&
                (requested || (record.birthTime == state.time.bigTime &&
                    retainedInactive < ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES));
            if (!live)
            {
                trace << "\nZHANG_HANDOFF_FLOAT_EXPIRY time="
                      << state.time.to_string(0)
                      << " scope=GLOBAL stage=" << stage << " status=PLANNED key=" << key
                      << " observable=" << record.identity.observable
                      << " key_num=" << key.num << " arc_version=" << record.identity.arcVersion
                      << " active_product_dependency=" << requested
                      << " age_seconds=" << double(state.time.bigTime - record.birthTime)
                      << " reason=" << (retentionSeconds == 0 ? "EXPLICIT_RETENTION_DISABLED" :
                          record.birthTime == state.time.bigTime ? "INACTIVE_CACHE_CAPACITY" :
                          "NO_ACTIVE_PRODUCT_DEPENDENCY")
                      << " integer_certification=0";
                removedRecords.insert(zhangHandoffFloatRecordKey(key));
                removedKeys.insert(key);
                ++retired;
                retiredDependencies += requested;
                continue;
            }
            ++retained;
            retainedActive += requested;
            retainedInactive += !requested;
        }
        keep[key][key] = 1;
    }
    if (retired == 0)
    {
        if (!zhangHandoffFloatRegistryValid(state, reason))
            return reject("FINAL_AGE_CONTRACT", reason);
        trace << "\nZHANG_PRODUCT_HANDOFF_AUXILIARY_RETIRE time="
              << state.time.to_string(0) << " scope=GLOBAL stage=" << stage
              << " status=UNCHANGED retired=0 retained=" << retained
              << " active_retained=" << retainedActive
              << " inactive_retained=" << retainedInactive
              << " active_above_cache_budget=" << (retainedActive > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES)
              << " retention_seconds=" << retentionSeconds
              << " strict_projected_registry_valid=1 integer_certification=0";
        trace.flush();
        return true;
    }
    if (keep.empty())
        return reject("PROJECTION_PLAN", "HANDOFF_FLOAT_RETIREMENT_EMPTY_FILTER");
    // Only the small metadata map is staged. The view validates the actual
    // surviving x/P principal submatrix and constant models without another
    // full KF/P copy. Retired dependency identities stay explicit; this sweep
    // does not certify continued publication of a removed direction.
    auto projectedMetadata = state.metaDataMap;
    for (const auto& record : removedRecords) projectedMetadata.erase(record);
    if (!zhangHandoffFloatRegistryViewValid(
            state, projectedMetadata, removedKeys, false, reason))
        return reject("PROJECTED_REGISTRY", reason);
    if (!state.applyStateTransform(
        trace, keep, "ZHANG_PRODUCT_HANDOFF_AUXILIARY_RETIREMENT", {}, &projectedMetadata))
        return reject("STATE_TRANSFORM", state.lastFactorTransactionFailureReason);
    // applyStateTransform already snapshots/rolls back the complete KF and
    // callback callable values. After acceptance, metadata commit cannot fail;
    // there is deliberately no new rejecting validation after this commit.
    trace << "\nZHANG_PRODUCT_HANDOFF_AUXILIARY_RETIRE time="
          << state.time.to_string(0)
          << " scope=GLOBAL stage=" << stage
          << " status=COMMITTED count=" << retired << " retired=" << retired
          << " retained_correlated_float=" << retained
          << " retained=" << retained
          << " active_retained=" << retainedActive
          << " inactive_retained=" << retainedInactive
          << " active_above_cache_budget=" << (retainedActive > ZHANG_HANDOFF_FLOAT_MAX_AUXILIARIES)
          << " retired_active_product_dependencies=" << retiredDependencies
          << " retention_seconds=" << retentionSeconds
          << " inactive_maximum_seconds=" << double(ZHANG_HANDOFF_FLOAT_MAX_AGE_SECONDS)
          << " active_retention=DECLARED_DEPENDENCY_LIFETIME"
          << " strict_projected_registry_valid=1"
          << " expired_dependency_requires_publication_revalidation=" << (retiredDependencies != 0)
          << " cross_covariance_marginalised=1 integer_certification=0";
    trace.flush();
    return true;
}

bool retireZhangGraphProductEventAuxiliaries(
    Trace& trace, KFState& state, E_Sys)
{
    return retireZhangGraphProductEventAuxiliaries(trace, state, "COMPAT_SYSTEM_API");
}

bool zhangGraphProductHandoffStateValid(
    const KFState& state, string& failureReason)
{
    return zhangHandoffFloatRegistryValid(state, failureReason) &&
        zhangProductFloatTargetMetadataValid(state, failureReason);
}

void zhangClearGraphProductHandoffDependencies(KFState& state)
{
    for (auto it = state.metaDataMap.begin(); it != state.metaDataMap.end();)
        it = it->first.starts_with(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX)
            ? state.metaDataMap.erase(it) : std::next(it);
}

bool zhangRegisterGraphProductHandoffDependencies(
    KFState& state, E_Sys system, E_ObsCode observable,
    const vector<ZhangGraphEdge>& edges, const vector<int>& versions)
{
    if (edges.size() != versions.size()) return false;
    map<string, string> additions;
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        if (edges[i].satellite.sys != system || versions[i] < 0 ||
            edges[i].receiver.empty()) return false;
        const auto identity = zhangHandoffFloatIdentityText(
            {edges[i], int(observable), versions[i]});
        additions[std::string(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX) + identity] = identity;
    }
    std::size_t total = additions.size();
    for (const auto& [name, value] : state.metaDataMap)
        if (name.starts_with(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX) &&
            !additions.contains(name)) ++total;
    if (total > ZHANG_HANDOFF_FLOAT_MAX_DEPENDENCIES) return false;
    state.metaDataMap.insert(additions.begin(), additions.end());
    return true;
}

namespace
{
bool resolveCapturePhysicalFloatArc(const KFState& state,
    const ZhangProductFloatVersionedTerm& term,ZhangCoordinateRow& row)
{
    row.clear();
    const double wavelength=zhangWavelength(term.edge.satellite.sys,term.observable);
    if(!(wavelength>0) || !std::isfinite(wavelength) || term.arcVersion<0)return false;
    const auto owner=graphStateMap.find({zhangGraphRuntimeId(state),term.edge.satellite.sys});
    if(owner!=graphStateMap.end() && owner->second.initialized)
    {
        const auto& runtime=owner->second;
        if(runtime.productBasisMode!=configuredProductBasisMode())return false;
        const auto edge=runtime.phaseArcVersions.find(term.edge);
        const auto& basis=!singleSBasisNative() && runtime.activeBasis.connected
            ?runtime.activeBasis:runtime.basis;
        if(edge!=runtime.phaseArcVersions.end() && edge->second.contains(int(term.observable)) &&
            edge->second.at(int(term.observable))==term.arcVersion && basis.edges.contains(term.edge))
        {
            if(term.edge.receiver!=basis.rootReceiver)
                row[zhangTransportReceiverKey(term.edge.satellite.sys,term.observable,term.edge.receiver)]=1;
            row[zhangTransportSatelliteKey(term.observable,term.edge.satellite)]=1;
            if(!basis.treeEdges.contains(term.edge))
                row[zhangTransportAmbiguityKey(term.observable,term.edge)]=
                    wavelength;
            bool valid=true;
            for(const auto& [key,value]:row)valid &= state.kfIndexMap.contains(key) && std::isfinite(value);
            if(valid)return true;
            row.clear();
        }
    }
    KFKey retained;
    if(zhangHandoffFloatFind(state,{term.edge,int(term.observable),term.arcVersion},retained))
    {row[retained]=1;return true;}
    return false;
}
bool sameWholeFloatMeaning(const ZhangProductFloatTargetRequest& a,
    const ZhangProductFloatTargetRequest& b)
{
    return a.targetId==b.targetId && a.consumer==b.consumer && a.lifecycleId==b.lifecycleId &&
        a.unit==b.unit && a.symbolicSpace==b.symbolicSpace &&
        a.physicalCoefficients==b.physicalCoefficients &&
        a.nuisanceCoefficients==b.nuisanceCoefficients && a.physicalScale==b.physicalScale && a.symbolicOffset==b.symbolicOffset;
}
}

bool zhangResolveGraphProductHandoffTarget(const KFState& state,
    const ZhangProductFloatTargetRequest& request,std::map<KFKey,double>& row,
    double& offset,std::string& reason)
{
    ZhangPhaseTimer timer("WHOLE_FLOAT_TARGET_RESOLVE");
    row.clear();offset=0;
    if(!zhangProductFloatTargetMetadataValid(state,reason))return false;
    std::vector<ZhangProductFloatTargetRequest> registered;
    if(!zhangReadFloatTargets(state.metaDataMap,registered,reason))return false;
    const auto owner=zhangGraphRuntimeId(state),order=zhangFloatCarrierStateOrder(state.kfIndexMap);
    auto usable=[&](const ZhangProductFloatTargetRequest& item) {
        if(!item.numericBindingComplete || item.scope.runtimeId!=owner || item.scope.stateOrderId!=order)return false;
        for(const auto& [key,value]:item.numericRow)
            if(!state.kfIndexMap.contains(key) || !std::isfinite(value))return false;
        return std::isfinite(item.numericOffset);
    };
    for(const auto& item:registered)
        if(sameWholeFloatMeaning(item,request) && usable(item))
        {row=item.numericRow;offset=item.numericOffset;reason="COMPLETE_REGISTERED_FLOAT_TARGET_BOUND";return true;}
    if(request.symbolicSpace!=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE ||
        request.unit!="m" || !request.nuisanceCoefficients.empty())
    {reason="WHOLE_FLOAT_TARGET_SYMBOLIC_BINDING_UNAVAILABLE";return false;}
    std::map<ZhangProductFloatVersionedTerm,std::size_t> columns;
    for(const auto& [term,value]:request.physicalCoefficients)if(value!=0)columns.try_emplace(term,columns.size());
    std::vector<const ZhangProductFloatTargetRequest*> candidates;
    for(const auto& item:registered)
        if(item.symbolicSpace==ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE &&
            item.nuisanceCoefficients.empty() && usable(item))
        {
            candidates.push_back(&item);
            for(const auto& [term,value]:item.physicalCoefficients)
                if(value!=0)columns.try_emplace(term,columns.size());
        }
    ZhangBridgeSpanMatrix symbols;
    std::vector<ZhangCoordinateRow> numeric;
    std::vector<double> constants;
    for(const auto& [term,column]:columns)
    {
        ZhangCoordinateRow physical;
        if(!resolveCapturePhysicalFloatArc(state,term,physical))continue;
        ZhangBridgeSpanRow symbol(columns.size());symbol[column]=1;
        symbols.push_back(std::move(symbol));numeric.push_back(std::move(physical));constants.push_back(0);
    }
    for(const auto* item:candidates)
    {
        ZhangBridgeSpanRow symbol(columns.size());
        for(const auto& [term,value]:item->physicalCoefficients)if(value!=0)
            symbol[columns.at(term)]=item->physicalScale*ZhangFloatCarrierRational(value);
        symbols.push_back(std::move(symbol));numeric.push_back(item->numericRow);
        constants.push_back(item->numericOffset-item->symbolicOffset.convert_to<double>());
    }
    ZhangBridgeSpanRow target(columns.size());
    for(const auto& [term,value]:request.physicalCoefficients)if(value!=0)
        target[columns.at(term)]=request.physicalScale*ZhangFloatCarrierRational(value);
    const auto analysis=analyzeZhangBridgeSupportOnActualSpan({},symbols,{target},columns.size());
    if(!analysis.valid || !analysis.representable)
    {reason="WHOLE_FLOAT_TARGET_OUTSIDE_ACTUAL_PHASE_EDGE_SPAN";return false;}
    offset=request.symbolicOffset.convert_to<double>();
    for(std::size_t i=0;i<numeric.size();++i)
    {
        const double coefficient=analysis.numericCombinations[0][i].convert_to<double>();
        if(!std::isfinite(coefficient)){reason="WHOLE_FLOAT_TARGET_COEFFICIENT_OUT_OF_RANGE";return false;}
        zhangAddCoordinateRow(row,numeric[i],coefficient);offset+=coefficient*constants[i];
    }
    if(!std::isfinite(offset)){reason="WHOLE_FLOAT_TARGET_OFFSET_NONFINITE";row.clear();return false;}
    reason="WHOLE_PHASE_EDGE_FLOAT_TARGET_RESOLVED_WITH_EXACT_SPAN_WITNESS";return true;
}

namespace
{
// The public one-request resolver above remains the byte-unchanged oracle.
// This private context exists only inside one synchronous const batch. There
// are no callbacks or writes to the source/graph in the query path below.
class ImmutableFloatTargetResolverContext
{
    const KFState& state;
    ZhangFloatTargetResolverBatchStats& stats;
    const std::map<std::string,std::string> sourceMetadata;
    const std::map<KFKey,int> sourceIndices;
    const std::uint64_t sourceSequence;
    const GTime sourceTime;
    std::string owner,order,physicalSnapshot,carrierSnapshot;
    std::vector<ZhangProductFloatTargetRequest> registered;
    std::vector<bool> usable;
    std::map<std::string,std::vector<std::size_t>> byId;
    std::vector<std::size_t> physicalCandidates;
    std::map<ZhangProductFloatVersionedTerm,std::pair<bool,ZhangCoordinateRow>> physicalRows;
    using Revision=std::tuple<bool,int,int,int,int,int,int>;
    std::map<E_Sys,Revision> revisions;

    std::string physicalIdentity() const
    {
        std::ostringstream out;
        // Snapshot the literal setting without introducing a mode parse on a
        // registered/SOURCE_KF hit (the old resolver never parses it there).
        out<<std::quoted(acsConfig.zhangFullRank.product_basis_mode);
        for(const auto& [key,runtime]:graphStateMap)if(key.first==owner)
        {
            out<<' '<<int(key.second)<<' '<<runtime.initialized<<' '<<std::quoted(runtime.productBasisMode);
            auto basis=[&](const ZhangGraphBasis& b) {
                out<<' '<<b.connected<<' '<<std::quoted(b.rootReceiver)<<' '<<b.edges.size();
                for(const auto& edge:b.edges)out<<' '<<std::quoted(edge.receiver)<<' '<<std::quoted(edge.satellite.id());
                out<<' '<<b.treeEdges.size();
                for(const auto& edge:b.treeEdges)out<<' '<<std::quoted(edge.receiver)<<' '<<std::quoted(edge.satellite.id());
            };
            basis(runtime.basis);basis(runtime.activeBasis);
            out<<' '<<runtime.phaseArcVersions.size();
            for(const auto& [edge,signals]:runtime.phaseArcVersions)
            {
                out<<' '<<std::quoted(edge.receiver)<<' '<<std::quoted(edge.satellite.id())<<' '<<signals.size();
                for(const auto& [code,version]:signals)out<<' '<<code<<' '<<version;
            }
        }
        return out.str();
    }
    std::string carrierModelIdentity() const
    {
        std::ostringstream out;
        // Raw bits preserve even irrelevant NaN/signed-zero representations;
        // no stricter model policy than the unchanged validator is introduced.
        auto number=[&](double value) {std::uint64_t bits;std::memcpy(&bits,&value,sizeof(bits));out<<' '<<bits;};
        for(const auto& [key,index]:state.kfIndexMap)if(zhangIsProductFloatCarrierKey(key))
        {
            out<<std::quoted(static_cast<std::string>(key));
            const auto model=state.stateTransitionMap.find(key);
            out<<' '<<(model!=state.stateTransitionMap.end());
            if(model!=state.stateTransitionMap.end())out<<' '<<model->second.size();
            if(model!=state.stateTransitionMap.end())for(const auto& [source,powers]:model->second)
            {
                out<<' '<<std::quoted(static_cast<std::string>(source))<<' '<<powers.size();
                for(const auto& [power,value]:powers){out<<' '<<power;number(value);}
            }
            for(const auto* values:{&state.procNoiseMap,&state.gaussMarkovTauMap})
            {
                const auto value=values->find(key);out<<' '<<(value!=values->end());
                if(value!=values->end())number(value->second);
            }
        }
        return out.str();
    }
    static Revision revision(const GraphRuntimeState& r)
    {return {r.initialized,r.epochIndex,r.representationVersion,r.floatGaugeVersion,
        r.integerComponentVersion,r.eventCounter,r.committedSGaugeVersion};}
public:
    ImmutableFloatTargetResolverContext(const KFState& input,ZhangFloatTargetResolverBatchStats& work)
        :state(input),stats(work),sourceMetadata(input.metaDataMap),sourceIndices(input.kfIndexMap),
         sourceSequence(input.factorCommitSequence),sourceTime(input.time){}

    bool prepare(std::string& reason)
    {
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_SOURCE_VALIDATE");
            ++stats.sourceValidatorCalls;
            if(!zhangProductFloatTargetMetadataValid(state,reason))return false;
        }
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_DIRECTORY_DECODE");
            ++stats.explicitDirectoryDecodeCalls;
            if(!zhangReadFloatTargets(state.metaDataMap,registered,reason))return false;
            stats.sourceTargets=registered.size();
        }
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_OWNER_ORDER");
            owner=zhangGraphRuntimeId(state);order=zhangFloatCarrierStateOrder(state.kfIndexMap);
            ++stats.stateOrderBuilds;
        }
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_INDEX");
            usable.resize(registered.size());
            for(std::size_t i=0;i<registered.size();++i)
            {
                const auto& item=registered[i];
                bool valid=item.numericBindingComplete && item.scope.runtimeId==owner &&
                    item.scope.stateOrderId==order && std::isfinite(item.numericOffset);
                for(const auto& [key,value]:item.numericRow)
                    valid &= state.kfIndexMap.contains(key) && std::isfinite(value);
                usable[i]=valid;byId[item.targetId].push_back(i);++stats.indexedTargets;
                if(valid && item.symbolicSpace==ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE &&
                    item.nuisanceCoefficients.empty())physicalCandidates.push_back(i);
            }
        }
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_SOURCE_SNAPSHOT");
            physicalSnapshot=physicalIdentity();carrierSnapshot=carrierModelIdentity();
            for(const auto& [key,runtime]:graphStateMap)if(key.first==owner)
                revisions.emplace(key.second,revision(runtime));
        }
        return unchanged(reason,true);
    }

    bool unchanged(std::string& reason,bool complete=false)
    {
        ++stats.sourceGuardChecks;
        if(state.factorCommitSequence!=sourceSequence || state.time!=sourceTime ||
            zhangGraphRuntimeId(state)!=owner || state.metaDataMap.size()!=sourceMetadata.size() ||
            state.kfIndexMap.size()!=sourceIndices.size())
        {reason="WHOLE_FLOAT_RESOLVER_SOURCE_CHANGED_DURING_BATCH";return false;}
        for(const auto& [system,expected]:revisions)
        {
            const auto found=graphStateMap.find({owner,system});
            if(found==graphStateMap.end() || revision(found->second)!=expected)
            {reason="WHOLE_FLOAT_RESOLVER_PHYSICAL_CHART_CHANGED_DURING_BATCH";return false;}
        }
        bool exactIndices=true;
        if(complete)
        {
            auto actual=state.kfIndexMap.begin(),expected=sourceIndices.begin();
            for(;actual!=state.kfIndexMap.end();++actual,++expected)
            {
                const auto& a=actual->first;const auto& b=expected->first;
                exactIndices &= actual->second==expected->second && a.type==b.type && a.Sat==b.Sat &&
                    a.str==b.str && a.num==b.num && a.comment==b.comment;
            }
        }
        if(complete && (state.metaDataMap!=sourceMetadata || !exactIndices ||
            physicalIdentity()!=physicalSnapshot || carrierModelIdentity()!=carrierSnapshot))
        {reason="WHOLE_FLOAT_RESOLVER_SOURCE_CHANGED_DURING_BATCH";return false;}
        return true;
    }

    ZhangFloatTargetResolution resolve(const ZhangProductFloatTargetRequest& request)
    {
        ZhangFloatTargetResolution result;
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_REGISTERED_LOOKUP");
            // The index only narrows the original map-order scan. Every match
            // still compares ALL sameWholeFloatMeaning fields; no id grants
            // an alias, and repeated requests remain repeated results.
            if(const auto found=byId.find(request.targetId);found!=byId.end())
                for(const auto index:found->second)
                    if(sameWholeFloatMeaning(registered[index],request) && usable[index])
                    {
                        result.resolved=true;result.row=registered[index].numericRow;
                        result.offset=registered[index].numericOffset;
                        result.reason="COMPLETE_REGISTERED_FLOAT_TARGET_BOUND";
                        ++stats.registeredMatches;return result;
                    }
        }
        if(request.symbolicSpace!=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE ||
            request.unit!="m" || !request.nuisanceCoefficients.empty())
        {result.reason="WHOLE_FLOAT_TARGET_SYMBOLIC_BINDING_UNAVAILABLE";return result;}
        ++stats.spanQueries;
        std::map<ZhangProductFloatVersionedTerm,std::size_t> columns;
        ZhangBridgeSpanMatrix symbols;std::vector<ZhangCoordinateRow> numeric;
        std::vector<double> constants;
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_SPAN_BUILD");
            for(const auto& [term,value]:request.physicalCoefficients)if(value!=0)
                columns.try_emplace(term,columns.size());
            for(const auto index:physicalCandidates)
                for(const auto& [term,value]:registered[index].physicalCoefficients)if(value!=0)
                    columns.try_emplace(term,columns.size());
            for(const auto& [term,column]:columns)
            {
                auto found=physicalRows.find(term);
                if(found==physicalRows.end())
                {
                    ZhangCoordinateRow physical;++stats.physicalLookups;
                    const bool valid=resolveCapturePhysicalFloatArc(state,term,physical);
                    found=physicalRows.emplace(term,std::make_pair(valid,std::move(physical))).first;
                }
                if(!found->second.first)continue;
                ZhangBridgeSpanRow symbol(columns.size());symbol[column]=1;
                symbols.push_back(std::move(symbol));numeric.push_back(found->second.second);constants.push_back(0);
            }
            for(const auto index:physicalCandidates)
            {
                const auto& item=registered[index];ZhangBridgeSpanRow symbol(columns.size());
                for(const auto& [term,value]:item.physicalCoefficients)if(value!=0)
                    symbol[columns.at(term)]=item.physicalScale*ZhangFloatCarrierRational(value);
                symbols.push_back(std::move(symbol));numeric.push_back(item.numericRow);
                constants.push_back(item.numericOffset-item.symbolicOffset.convert_to<double>());
            }
        }
        ZhangBridgeSpanRow target(columns.size());
        for(const auto& [term,value]:request.physicalCoefficients)if(value!=0)
            target[columns.at(term)]=request.physicalScale*ZhangFloatCarrierRational(value);
        const auto analysis=[&]() {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_EXACT_SPAN_QUERY");
            return analyzeZhangBridgeSupportOnActualSpan({},symbols,{target},columns.size());
        }();
        if(!analysis.valid || !analysis.representable)
        {result.reason="WHOLE_FLOAT_TARGET_OUTSIDE_ACTUAL_PHASE_EDGE_SPAN";return result;}
        {
            ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_NUMERIC_RECONSTRUCTION");
            result.offset=request.symbolicOffset.convert_to<double>();
            for(std::size_t i=0;i<numeric.size();++i)
            {
                const double coefficient=analysis.numericCombinations[0][i].convert_to<double>();
                if(!std::isfinite(coefficient))
                {result.reason="WHOLE_FLOAT_TARGET_COEFFICIENT_OUT_OF_RANGE";return result;}
                zhangAddCoordinateRow(result.row,numeric[i],coefficient);result.offset+=coefficient*constants[i];
            }
        }
        if(!std::isfinite(result.offset))
        {result.reason="WHOLE_FLOAT_TARGET_OFFSET_NONFINITE";result.row.clear();return result;}
        result.resolved=true;result.reason="WHOLE_PHASE_EDGE_FLOAT_TARGET_RESOLVED_WITH_EXACT_SPAN_WITNESS";
        return result;
    }
    std::vector<ZhangProductFloatTargetRequest> takeDirectory(){return std::move(registered);}
};
}

bool zhangResolveGraphProductHandoffTargetsBatch(const KFState& state,
    const std::vector<ZhangProductFloatTargetRequest>& requests,
    std::vector<ZhangFloatTargetResolution>& output,std::string& reason,
    ZhangFloatTargetResolverBatchStats* statistics,
    std::vector<ZhangProductFloatTargetRequest>* registeredDirectory)
{
    ZhangPhaseTimer total("WHOLE_FLOAT_TARGET_RESOLVE_BATCH");
    output.clear();if(registeredDirectory)registeredDirectory->clear();
    ZhangFloatTargetResolverBatchStats local;
    auto& stats=statistics?*statistics:local;stats={};stats.requestedTargets=requests.size();
    auto context=[&]() {
        ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_CREATE_CONTEXT");
        return std::make_unique<ImmutableFloatTargetResolverContext>(state,stats);
    }();
    if(!context->prepare(reason))return false;
    std::vector<ZhangFloatTargetResolution> staged;staged.reserve(requests.size());
    for(const auto& request:requests)
    {
        if(!context->unchanged(reason))return false;
        staged.push_back(context->resolve(request));
        if(staged.back().resolved)++stats.resolvedTargets;else ++stats.unboundTargets;
    }
    {
        ZhangPhaseTimer timer("WHOLE_FLOAT_BATCH_FINAL_SOURCE_GUARD");
        if(!context->unchanged(reason,true))return false;
    }
    if(registeredDirectory)*registeredDirectory=context->takeDirectory();
    output.swap(staged);reason="IMMUTABLE_FLOAT_TARGET_BATCH_RESOLVED";
    return true;
}

bool zhangRegisterGraphProductHandoffTarget(KFState& state,
    const ZhangProductFloatTargetRequest& requested,std::string& reason)
{
    if(!zhangProductFloatTargetMetadataValid(state,reason))return false;
    if(requested.targetId.empty() || requested.unit.empty())
    {reason="WHOLE_FLOAT_TARGET_REGISTRATION_ID_OR_UNIT_MISSING";return false;}
    auto target=requested;
    const auto owner=zhangGraphRuntimeId(state),order=zhangFloatCarrierStateOrder(state.kfIndexMap);
    std::uint64_t generation=state.factorCommitSequence;
    if(const auto previous=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);
        previous!=state.metaDataMap.end())
    {
        try {std::size_t used=0;generation=std::stoull(previous->second,&used);
            if(used!=previous->second.size() || generation>state.factorCommitSequence)
            {reason="FLOAT_TARGET_REGISTRATION_GENERATION_INVALID";return false;}}
        catch (...){reason="FLOAT_TARGET_REGISTRATION_GENERATION_INVALID";return false;}
    }
    if(owner.empty() || order.empty() || !std::isfinite(target.numericOffset) ||
        (!target.scope.runtimeId.empty() && target.scope.runtimeId!=owner) ||
        (target.numericBindingComplete && (target.scope.runtimeId!=owner ||
            target.scope.stateOrderId!=order || target.scope.chartId.empty() ||
            target.scope.coordinateGeneration!=generation)))
    {reason="WHOLE_FLOAT_TARGET_REGISTRATION_OWNER_INVALID";return false;}
    if(!target.numericBindingComplete)
    {
        double offset=0;ZhangCoordinateRow row;
        if(zhangResolveGraphProductHandoffTarget(state,target,row,offset,reason))
        {target.numericRow=std::move(row);target.numericOffset=offset;target.numericBindingComplete=true;}
    }
    if(target.numericBindingComplete)
        for(const auto& [key,value]:target.numericRow)
            if(!state.kfIndexMap.contains(key) || !std::isfinite(value))
            {reason="WHOLE_FLOAT_TARGET_REGISTRATION_NUMERIC_ROW_INVALID";return false;}
    target.scope.runtimeId=owner;target.scope.stateOrderId=order;
    target.scope.coordinateGeneration=generation;
    if(target.numericBindingComplete && target.scope.chartId.empty())
        target.scope.chartId="CAPTURED_WHOLE_PHASE_EDGE_FLOAT_SOURCE:"+std::to_string(generation);
    if(target.numericBindingComplete)target.bindingReason="COMPLETE_CURRENT_FLOAT_ROW_REGISTERED";
    else target.bindingReason="ACTIVE_WHOLE_FLOAT_TARGET_UNBOUND:"+target.bindingReason;
    auto metadata=state.metaDataMap;
    metadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
    metadata[ZHANG_FLOAT_TARGET_SCOPE_V2]=order;
    metadata[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(generation);
    metadata[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(state.factorCommitSequence);
    state.metaDataMap.swap(metadata);
    reason=target.numericBindingComplete?"BOUND_WHOLE_FLOAT_TARGET_REGISTERED":"UNBOUND_WHOLE_FLOAT_TARGET_RETAINED";
    return true;
}

namespace
{
// Kept separate from the old validator: the test can still ask the original
// MetadataValid to independently accept the actual post-batch KFState. This
// view never copies the network P merely to supply a staged metadata map.
bool validateFloatTargetBatchView(const KFState& state,
    const std::map<std::string,std::string>& metadata,
    const std::string& owner,const std::string& order,std::string& reason)
{
    if(hasZhangFloatTargetV2State(metadata))
    {
        std::uint64_t saved=0;
        if(!readZhangFloatTargetSavedFactorSequence(metadata,owner,order,saved,reason) ||
            saved!=state.factorCommitSequence)
        {if(reason=="FLOAT_TARGET_ACTUAL_FACTOR_SEQUENCE_VALIDATED")reason="FLOAT_TARGET_SAVED_SEQUENCE_NOT_ACTUAL";return false;}
    }
    std::vector<ZhangProductFloatTargetRequest> targets;
    if(!zhangReadFloatTargets(metadata,targets,reason))return false;
    std::optional<std::uint64_t> generation;
    if(const auto entry=metadata.find(ZHANG_FLOAT_TARGET_GENERATION_V2);entry!=metadata.end())
        try {std::size_t used=0;const auto value=std::stoull(entry->second,&used);
            if(used==entry->second.size())generation=value;}catch(...){}
    for(const auto& target:targets)
    {
        if(target.scope.runtimeId!=owner){reason="FLOAT_TARGET_OWNER_MISMATCH";return false;}
        if(!target.numericBindingComplete)continue;
        if(!generation || target.scope.coordinateGeneration!=*generation || *generation>state.factorCommitSequence)
        {reason="FLOAT_TARGET_COMPLETE_GENERATION_MISMATCH";return false;}
        if(target.scope.stateOrderId!=order){reason="FLOAT_TARGET_COMPLETE_ORDER_MISMATCH";return false;}
        for(const auto& [key,value]:target.numericRow)
            if(!std::isfinite(value) || !state.kfIndexMap.contains(key))
            {reason="FLOAT_TARGET_BOUND_KEY_OR_COEFFICIENT_INVALID";return false;}
    }
    std::set<std::string> records;
    for(const auto& [key,index]:state.kfIndexMap)if(zhangIsProductFloatCarrierKey(key))
    {
        if(!zhangProductFloatCarrierCoordinateValid(state,metadata,key,reason))return false;
        records.insert(zhangProductFloatCarrierRecordKey(key));
    }
    for(const auto& [name,value]:metadata)
        if(name.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX) && !records.contains(name))
        {reason="FLOAT_CARRIER_ORPHAN_DEFINITION";return false;}
    reason="WHOLE_FLOAT_TARGET_AND_CARRIER_METADATA_VALID";return true;
}

// The old public resolver is deliberately unchanged. This transaction-local
// equivalent reads the already validated directory, including earlier ordered
// upserts. Its map order is exactly the original metadata-key order.
bool resolveFloatTargetOnBatchDirectory(const KFState& state,
    const ZhangProductFloatTargetRequest& request,
    const std::map<std::string,ZhangProductFloatTargetRequest>& directory,
    const std::string& owner,const std::string& order,
    ZhangCoordinateRow& row,double& offset,std::string& reason)
{
    row.clear();offset=0;
    auto usable=[&](const ZhangProductFloatTargetRequest& item) {
        if(!item.numericBindingComplete || item.scope.runtimeId!=owner || item.scope.stateOrderId!=order)return false;
        for(const auto& [key,value]:item.numericRow)
            if(!state.kfIndexMap.contains(key) || !std::isfinite(value))return false;
        return std::isfinite(item.numericOffset);
    };
    for(const auto& [id,item]:directory)
        if(sameWholeFloatMeaning(item,request) && usable(item))
        {row=item.numericRow;offset=item.numericOffset;reason="COMPLETE_REGISTERED_FLOAT_TARGET_BOUND";return true;}
    if(request.symbolicSpace!=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE ||
        request.unit!="m" || !request.nuisanceCoefficients.empty())
    {reason="WHOLE_FLOAT_TARGET_SYMBOLIC_BINDING_UNAVAILABLE";return false;}
    std::map<ZhangProductFloatVersionedTerm,std::size_t> columns;
    for(const auto& [term,value]:request.physicalCoefficients)if(value!=0)columns.try_emplace(term,columns.size());
    std::vector<const ZhangProductFloatTargetRequest*> candidates;
    for(const auto& [id,item]:directory)
        if(item.symbolicSpace==ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE &&
            item.nuisanceCoefficients.empty() && usable(item))
        {
            candidates.push_back(&item);
            for(const auto& [term,value]:item.physicalCoefficients)
                if(value!=0)columns.try_emplace(term,columns.size());
        }
    ZhangBridgeSpanMatrix symbols;
    std::vector<ZhangCoordinateRow> numeric;
    std::vector<double> constants;
    for(const auto& [term,column]:columns)
    {
        ZhangCoordinateRow physical;
        if(!resolveCapturePhysicalFloatArc(state,term,physical))continue;
        ZhangBridgeSpanRow symbol(columns.size());symbol[column]=1;
        symbols.push_back(std::move(symbol));numeric.push_back(std::move(physical));constants.push_back(0);
    }
    for(const auto* item:candidates)
    {
        ZhangBridgeSpanRow symbol(columns.size());
        for(const auto& [term,value]:item->physicalCoefficients)if(value!=0)
            symbol[columns.at(term)]=item->physicalScale*ZhangFloatCarrierRational(value);
        symbols.push_back(std::move(symbol));numeric.push_back(item->numericRow);
        constants.push_back(item->numericOffset-item->symbolicOffset.convert_to<double>());
    }
    ZhangBridgeSpanRow target(columns.size());
    for(const auto& [term,value]:request.physicalCoefficients)if(value!=0)
        target[columns.at(term)]=request.physicalScale*ZhangFloatCarrierRational(value);
    const auto analysis=analyzeZhangBridgeSupportOnActualSpan({},symbols,{target},columns.size());
    if(!analysis.valid || !analysis.representable)
    {reason="WHOLE_FLOAT_TARGET_OUTSIDE_ACTUAL_PHASE_EDGE_SPAN";return false;}
    offset=request.symbolicOffset.convert_to<double>();
    for(std::size_t i=0;i<numeric.size();++i)
    {
        const double coefficient=analysis.numericCombinations[0][i].convert_to<double>();
        if(!std::isfinite(coefficient)){reason="WHOLE_FLOAT_TARGET_COEFFICIENT_OUT_OF_RANGE";return false;}
        zhangAddCoordinateRow(row,numeric[i],coefficient);offset+=coefficient*constants[i];
    }
    if(!std::isfinite(offset)){reason="WHOLE_FLOAT_TARGET_OFFSET_NONFINITE";row.clear();return false;}
    reason="WHOLE_PHASE_EDGE_FLOAT_TARGET_RESOLVED_WITH_EXACT_SPAN_WITNESS";return true;
}
}

bool zhangRefreshGraphProductHandoffTargetsBatch(KFState& state,
    const std::vector<ZhangProductFloatTargetRequest>& requested,
    std::string& reason,ZhangFloatTargetBatchStats* statistics)
{
    ZhangPhaseTimer total("PRE_GRAPH_FLOAT_TARGET_BATCH_REFRESH");
    ZhangFloatTargetBatchStats local;
    auto& stats=statistics?*statistics:local;stats={};
    stats.requestedTargets=requested.size();stats.sourceMetadataEntries=state.metaDataMap.size();
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_SOURCE_VALIDATE");
        ++stats.sourceValidatorCalls;
        if(!zhangProductFloatTargetMetadataValid(state,reason))return false;
    }
    const auto owner=zhangGraphRuntimeId(state),order=zhangFloatCarrierStateOrder(state.kfIndexMap);
    stats.stateOrderBytes=order.size();
    std::map<std::string,ZhangProductFloatTargetRequest> directory;
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_DIRECTORY_DECODE");
        std::vector<ZhangProductFloatTargetRequest> parsed;
        ++stats.sourceDirectoryDecodeCalls;
        if(!zhangReadFloatTargets(state.metaDataMap,parsed,reason))return false;
        stats.sourceTargets=parsed.size();
        for(auto& item:parsed)directory.emplace(item.targetId,std::move(item));
    }
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_CLEAR_ACTIVE");
        for(auto& [id,item]:directory)item.active=false;
    }
    std::uint64_t generation=state.factorCommitSequence;
    if(const auto previous=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);previous!=state.metaDataMap.end())
        generation=std::stoull(previous->second); // source validator has proved this decimal
    std::set<std::string> requestedIds;
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_REGISTER_TYPED");
        for(const auto& input:requested)
        {
            auto target=input;
            if(target.targetId.empty() || target.unit.empty())
            {reason="WHOLE_FLOAT_TARGET_REGISTRATION_ID_OR_UNIT_MISSING";return false;}
            if(owner.empty() || order.empty() || !std::isfinite(target.numericOffset) ||
                (!target.scope.runtimeId.empty() && target.scope.runtimeId!=owner) ||
                (target.numericBindingComplete && (target.scope.runtimeId!=owner ||
                    target.scope.stateOrderId!=order || target.scope.chartId.empty() ||
                    target.scope.coordinateGeneration!=generation)))
            {reason="WHOLE_FLOAT_TARGET_REGISTRATION_OWNER_INVALID";return false;}
            if(!target.numericBindingComplete)
            {
                ZhangPhaseTimer resolver("PRE_GRAPH_FLOAT_TARGET_UNBOUND_RESOLVE");
                ++stats.resolverCalls;
                double offset=0;ZhangCoordinateRow row;
                if(resolveFloatTargetOnBatchDirectory(state,target,directory,owner,order,row,offset,reason))
                {target.numericRow=std::move(row);target.numericOffset=offset;target.numericBindingComplete=true;}
            }
            if(target.numericBindingComplete)
                for(const auto& [key,value]:target.numericRow)
                    if(!state.kfIndexMap.contains(key) || !std::isfinite(value))
                    {reason="WHOLE_FLOAT_TARGET_REGISTRATION_NUMERIC_ROW_INVALID";return false;}
            target.scope.runtimeId=owner;target.scope.stateOrderId=order;target.scope.coordinateGeneration=generation;
            if(target.numericBindingComplete && target.scope.chartId.empty())
                target.scope.chartId="CAPTURED_WHOLE_PHASE_EDGE_FLOAT_SOURCE:"+std::to_string(generation);
            if(target.numericBindingComplete)target.bindingReason="COMPLETE_CURRENT_FLOAT_ROW_REGISTERED";
            else target.bindingReason="ACTIVE_WHOLE_FLOAT_TARGET_UNBOUND:"+target.bindingReason;
            // The old loop would decode this entry at the next iteration.
            // Do not let a later duplicate repair an invalid intermediate
            // wire entry and thereby widen that old acceptance domain.
            // This is only the existing codec gate: unbound numeric rows are
            // still NOT required to name an existing KF column.
            ++stats.requestWireValidationCalls;
            ZhangProductFloatTargetRequest wireCheck;
            if(!zhangDecodeFloatTarget(zhangEncodeFloatTarget(target),wireCheck))
            {reason="FLOAT_TARGET_BATCH_REQUEST_WIRE_INVALID";return false;}
            if(!requestedIds.insert(target.targetId).second)++stats.duplicateRequests;
            if(target.numericBindingComplete)++stats.registeredBound;else ++stats.registeredUnbound;
            directory.insert_or_assign(target.targetId,std::move(target));
        }
    }
    std::map<std::string,std::string> metadata;
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_DIRECTORY_SERIALIZE");
        metadata=state.metaDataMap; // the one complete metadata copy, not one per target
        for(const auto& [id,target]:directory)
        {
            auto encoded=zhangEncodeFloatTarget(target);
            ++stats.serializedTargets;stats.serializedTargetBytes+=encoded.size();
            metadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+id]=std::move(encoded);
        }
        if(!requested.empty())
        {
            metadata[ZHANG_FLOAT_TARGET_SCOPE_V2]=order;
            metadata[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(generation);
            metadata[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(state.factorCommitSequence);
        }
    }
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_FINAL_VIEW_VALIDATE");
        ++stats.finalValidatorCalls;
        if(!validateFloatTargetBatchView(state,metadata,owner,order,reason))return false;
    }
    stats.finalTargets=directory.size();stats.finalMetadataEntries=metadata.size();
    {
        ZhangPhaseTimer timer("PRE_GRAPH_FLOAT_TARGET_ATOMIC_SWAP");
        state.metaDataMap.swap(metadata);++stats.atomicSwaps;
    }
    if(auto* trace=ZhangPhaseTimer::currentTrace())
        *trace<<"\nZHANG_FLOAT_TARGET_BATCH epoch="<<std::quoted(state.time.to_string(0))
            <<" status=COMMITTED requested="<<stats.requestedTargets<<" source_targets="<<stats.sourceTargets
            <<" final_targets="<<stats.finalTargets<<" source_metadata_entries="<<stats.sourceMetadataEntries
            <<" final_metadata_entries="<<stats.finalMetadataEntries<<" state_order_bytes="<<stats.stateOrderBytes
            <<" duplicate_requests="<<stats.duplicateRequests<<" bound_upserts="<<stats.registeredBound
            <<" unbound_upserts="<<stats.registeredUnbound<<" unbound_resolver_calls="<<stats.resolverCalls
            <<" source_validator_calls="<<stats.sourceValidatorCalls<<" source_explicit_decode_calls="<<stats.sourceDirectoryDecodeCalls
            <<" final_validator_calls="<<stats.finalValidatorCalls<<" serialized_targets="<<stats.serializedTargets
            <<" request_wire_validation_calls="<<stats.requestWireValidationCalls
            <<" serialized_target_bytes="<<stats.serializedTargetBytes<<" atomic_swaps="<<stats.atomicSwaps
            <<" x_P_dx_changed=0 integer_domain_changed=0";
    reason="WHOLE_FLOAT_TARGET_DIRECTORY_BATCH_REFRESHED";return true;
}


bool zhangCompressGraphProductFloatTargets(std::ostream& trace,KFState& state,std::string& reason)
{
    ZhangPhaseTimer timer(trace,"GLOBAL_FLOAT_TARGET_CARRIER_COMPRESSION");
    if(!zhangR51SingleBasisHandoffEnabled()) {reason="FEATURE_OFF";return true;}
    if(!zhangProductFloatTargetMetadataValid(state,reason))return false;
    const auto owner=zhangGraphRuntimeId(state),order=zhangFloatCarrierStateOrder(state.kfIndexMap);
    std::vector<ZhangProductFloatTargetRequest> consumers,registered;
    auto skip=[&](const std::string& why) {
        reason=why;trace<<"\nZHANG_FLOAT_TARGET_COMPRESSION status=RAW_PRESERVED reason="<<std::quoted(why)
            <<" integer_domain_unchanged=1 full_covariance_unchanged=1";return true;};
    if(!zhangCollectActiveProductFloatConsumerTargets(state,consumers,reason))return skip("CONSUMER_DIRECTORY_INCOMPLETE:"+reason);
    if(!zhangReadFloatTargets(state.metaDataMap,registered,reason))return false;
    std::uint64_t generation=state.factorCommitSequence;
    if(const auto item=state.metaDataMap.find(ZHANG_FLOAT_TARGET_GENERATION_V2);item!=state.metaDataMap.end())
        generation=std::stoull(item->second); // strict schema already checked above
    std::map<std::string,ZhangProductFloatTargetRequest> all;
    auto add=[&](ZhangProductFloatTargetRequest target) {
        if(!target.active)return true;
        if(!target.numericBindingComplete || !std::isfinite(target.numericOffset) || target.scope.runtimeId!=owner || target.scope.stateOrderId!=order ||
            target.scope.coordinateGeneration!=generation)return false;
        for(const auto& [key,value]:target.numericRow)
            if(!state.kfIndexMap.contains(key) || !std::isfinite(value))return false;
        const auto existing=all.find(target.targetId);
        if(existing!=all.end())return sameWholeFloatMeaning(existing->second,target) &&
            existing->second.numericRow==target.numericRow && existing->second.numericOffset==target.numericOffset;
        all.emplace(target.targetId,std::move(target));return true;};
    for(const auto& target:registered)
        if(target.active && !target.consumer.starts_with("GRAPH_RECEIPT_") && !add(target))
            return skip("ACTIVE_REGISTERED_TARGET_UNBOUND:"+target.targetId);
    for(auto& target:consumers)if(!add(std::move(target)))return skip("ACTIVE_CONSUMER_TARGET_UNBOUND_OR_CONFLICT");
    auto stagedReceipts=graphProductEventReceipts;
    std::map<std::string,ZhangCoordinateRow*> receiptRows;
    auto receiptTarget=[&](const GraphProductReceiptKey& key,const std::string& role,
        ZhangCoordinateRow& row,bool available) {
        const auto& [runtime,epoch,system,code,satellite]=key;
        if(!available || row.empty())return false;
        ZhangProductFloatTargetRequest target;
        target.targetId="GRAPH_RECEIPT_"+role+"|"+epoch+"|"+satellite.id()+"|"+enum_to_string(code);
        target.consumer="GRAPH_RECEIPT_"+role;target.lifecycleId=epoch;target.unit=role=="BRIDGE"?"cycles":"m";
        target.scope={owner,"POST_GRAPH_RECEIPT",order,generation};
        target.numericRow=row;target.numericBindingComplete=true;
        const auto id=target.targetId;
        if(!add(std::move(target)))return false;
        receiptRows.emplace(id,&row);return true;};
    for(auto& [key,receipt]:stagedReceipts)
    {
        if(std::get<0>(key)!=owner || receipt.time!=state.time)continue;
        if(receipt.oldProductProjectionCaptured &&
            !receiptTarget(key,"OLD_MAPPED",receipt.newPhaseRow,receipt.representableSameTimeMapProven))
            return skip("ACTIVE_EVENT_OLD_FLOAT_TARGET_NOT_COMPLETELY_MAPPED");
        if(receipt.actualWriterProjectionBound &&
            !receiptTarget(key,"WRITER",receipt.actualWriterPhaseRow,true))
            return skip("ACTIVE_EVENT_WRITER_FLOAT_TARGET_NOT_COMPLETELY_MAPPED");
        if(!receipt.bridgeStateRowCycles.empty() && !receiptTarget(key,"BRIDGE",receipt.bridgeStateRowCycles,true))
            return skip("ACTIVE_EVENT_BRIDGE_FLOAT_TARGET_NOT_COMPLETELY_MAPPED");
    }
    std::vector<KFKey> sourceKeys(state.kfIndexMap.size());
    for(const auto& [key,index]:state.kfIndexMap)sourceKeys.at(index)=key;
    const auto dependencies=zhangHandoffFloatDependencies(state);
    std::vector<std::size_t> protectedIndices;
    for(const auto& [key,index]:state.kfIndexMap)
    {
        if(key.type!=KF::PRODUCT_HANDOFF){protectedIndices.push_back(index);continue;}
        if(zhangIsProductFloatCarrierKey(key))continue;
        ZhangHandoffFloatRecord record;
        if(!zhangHandoffFloatReadRecord(state,key,record)){reason="RAW_AUXILIARY_REGISTRY_INVALID";return false;}
        if(!dependencies.contains(record.identity))protectedIndices.push_back(index); // bounded inactive cache
    }
    std::vector<ZhangFloatCarrierIndexedTarget> indexed;
    for(const auto& [id,target]:all)
    {
        ZhangFloatCarrierIndexedTarget row;row.targetId=id;row.unit=target.unit;row.complete=true;
        row.offset=*zhangFloatCarrierExactDouble(target.numericOffset);
        for(const auto& [key,value]:target.numericRow)
        {
            const auto exact=zhangFloatCarrierExactDouble(value);
            if(!exact)return skip("NONFINITE_COMPLETE_TARGET_COEFFICIENT");
            if(*exact!=0)row.sourceRow.emplace(state.kfIndexMap.at(key),*exact);
        }
        indexed.push_back(std::move(row));
    }
    const auto plan=planZhangProductFloatCarrierRows(sourceKeys.size(),protectedIndices,indexed);
    if(!plan.valid)return skip("EXACT_ROW_PLAN_DECLINED:"+plan.reason);
    if(!plan.reducesDimension)return skip("NO_TARGET_SPACE_RANK_REDUCTION");
    std::map<KFKey,ZhangCoordinateRow> transform;
    std::vector<KFKey> destinationKeys;
    std::vector<std::pair<KFKey,ZhangProductFloatTargetRequest>> symbolicBindings;
    auto workingMetadata=state.metaDataMap;
    for(auto item=workingMetadata.begin();item!=workingMetadata.end();)
        if(item->first.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX))item=workingMetadata.erase(item);else ++item;
    for(std::size_t row=0;row<plan.transformRows.size();++row)
    {
        KFKey key;
        if(row<plan.protectedCount)key=sourceKeys.at(protectedIndices.at(row));
        else
        {
            key.type=KF::PRODUCT_HANDOFF;
            key.str="ZHANG_FLOAT_CARRIER:"+owner+":"+std::to_string(state.factorCommitSequence+1)+":"+std::to_string(row-plan.protectedCount);
        }
        destinationKeys.push_back(key);
        for(const auto& [source,coefficient]:plan.transformRows[row])
        {
            const double value=coefficient.convert_to<double>();
            if(!std::isfinite(value) || zhangFloatCarrierExactDouble(value)!=std::optional<ZhangFloatCarrierRational>(coefficient))
                return skip("CARRIER_TRANSFORM_NOT_EXACTLY_DOUBLE_REPRESENTABLE");
            transform[key].emplace(sourceKeys.at(source),value);
        }
        if(row>=plan.protectedCount)
        {
            ZhangProductFloatTargetRequest definition;
            definition.targetId=key.str;definition.consumer="FLOAT_CARRIER_SOURCE_DEFINITION";
            definition.lifecycleId=std::to_string(state.factorCommitSequence+1);
            definition.unit=plan.carrierUnits.at(row-plan.protectedCount);
            definition.scope={owner,"EXACT_SOURCE_TARGET_RESIDUAL",order,generation};
            definition.numericRow=transform.at(key);definition.active=false;
            definition.bindingReason="ARCHIVED_SOURCE_ROW_DEFINITION_NOT_A_PHYSICAL_INTEGER";
            std::map<ZhangProductFloatVersionedTerm,ZhangFloatCarrierRational> phase;
            for(const auto& [source,coefficient]:plan.transformRows[row])
            {
                const auto& oldKey=sourceKeys.at(source);
                if(zhangIsProductFloatCarrierKey(oldKey))
                {
                    const auto old=state.metaDataMap.find(zhangProductFloatCarrierRecordKey(oldKey));
                    ZhangProductFloatTargetRequest oldDefinition;
                    if(old==state.metaDataMap.end() || !zhangDecodeFloatTarget(old->second,oldDefinition) ||
                        oldDefinition.symbolicSpace!=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE ||
                        !oldDefinition.nuisanceCoefficients.empty())
                        return skip("OLD_CARRIER_WHOLE_PHASE_DEFINITION_UNAVAILABLE");
                    for(const auto& [term,value]:oldDefinition.physicalCoefficients)
                        phase[term]+=coefficient*oldDefinition.physicalScale*ZhangFloatCarrierRational(value);
                }
                else
                {
                    ZhangHandoffFloatRecord old;
                    if(!zhangHandoffFloatReadRecord(state,oldKey,old))return skip("RAW_PHASE_DEFINITION_UNAVAILABLE");
                    phase[{old.identity.edge,static_cast<E_ObsCode>(old.identity.observable),old.identity.arcVersion}]+=coefficient;
                }
            }
            ZhangExactInteger denominatorScale=1;
            for(auto item=phase.begin();item!=phase.end();)
                if(item->second==0)item=phase.erase(item);else
                {
                    ZhangExactInteger a=denominatorScale,b=denominator(item->second);
                    while(b!=0){ZhangExactInteger r=a%b;a=b;b=r;}
                    denominatorScale=denominatorScale/a*denominator(item->second);++item;
                }
            definition.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
            definition.physicalScale=ZhangFloatCarrierRational(1)/ZhangFloatCarrierRational(denominatorScale);
            for(const auto& [term,value]:phase)
            {
                const ZhangFloatCarrierRational integer=value*denominatorScale;
                definition.physicalCoefficients.emplace(term,numerator(integer));
            }
            workingMetadata[zhangProductFloatCarrierRecordKey(key)]=zhangEncodeFloatTarget(definition);
            auto binding=definition;
            binding.targetId="FLOAT_CARRIER_BINDING|"+key.str;
            binding.consumer="FLOAT_CARRIER_SYMBOLIC_BINDING";
            binding.numericBindingComplete=true;
            binding.bindingReason="EXACT_FLOAT_PHASE_COMBINATION_WITHOUT_INTEGER_AUTHORITY";
            symbolicBindings.emplace_back(key,std::move(binding));
        }
    }
    std::map<std::string,ZhangCoordinateRow> hints;
    for(const auto& resolved:plan.targets)
    {
        auto& row=hints[resolved.targetId];
        for(const auto& [index,coefficient]:resolved.destinationRow)
        {
            const double value=coefficient.convert_to<double>();
            if(!std::isfinite(value) || zhangFloatCarrierExactDouble(value)!=std::optional<ZhangFloatCarrierRational>(coefficient))
                return skip("TARGET_RECONSTRUCTION_NOT_EXACTLY_DOUBLE_REPRESENTABLE");
            row.emplace(destinationKeys.at(index),value);
        }
        if(const auto receipt=receiptRows.find(resolved.targetId);receipt!=receiptRows.end())*receipt->second=row;
    }
    for(const auto& [key,binding]:symbolicBindings)hints[binding.targetId]={{key,1}};
    for(const auto& [key,index]:state.kfIndexMap)
        if(key.type==KF::PRODUCT_HANDOFF && !zhangIsProductFloatCarrierKey(key) && !transform.contains(key))
            workingMetadata.erase(zhangHandoffFloatRecordKey(key));
    // Add the entire actual consumer set as SOURCE rows BEFORE the central
    // transaction; the override may contain only raw/carrier directory changes.
    auto sourceMetadata=state.metaDataMap;
    for(auto item=sourceMetadata.begin();item!=sourceMetadata.end();)
        if(item->first.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX))
        {
            ZhangProductFloatTargetRequest previous;
            if(!zhangDecodeFloatTarget(item->second,previous)){reason="SOURCE_TARGET_METADATA_INVALID";return false;}
            if(!previous.active || previous.consumer.starts_with("GRAPH_RECEIPT_"))item=sourceMetadata.erase(item);else ++item;
        }
        else ++item;
    for(const auto& [id,target]:all)sourceMetadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+id]=zhangEncodeFloatTarget(target);
    for(const auto& [key,binding]:symbolicBindings)
        sourceMetadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+binding.targetId]=zhangEncodeFloatTarget(binding);
    sourceMetadata[ZHANG_FLOAT_TARGET_SCOPE_V2]=order;
    sourceMetadata[ZHANG_FLOAT_TARGET_GENERATION_V2]=std::to_string(generation);
    sourceMetadata[ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2]=std::to_string(state.factorCommitSequence);
    // Copy the NEW source target table into the override unchanged. Central
    // owner/core equality checks then see one and the same source snapshot.
    for(auto item=workingMetadata.begin();item!=workingMetadata.end();)
        if(item->first.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX))item=workingMetadata.erase(item);else ++item;
    for(const auto& [name,value]:sourceMetadata)
        if(name.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || name==ZHANG_FLOAT_TARGET_SCOPE_V2 ||
            name==ZHANG_FLOAT_TARGET_GENERATION_V2 || name==ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2)
            workingMetadata[name]=value;
    auto callback=state.exactStateTransformCallback;
    const auto oldFailure=state.lastFactorTransactionFailureReason;
    state.metaDataMap.swap(sourceMetadata); // sourceMetadata now owns exact rollback metadata
    struct Restore {
        KFState& state;std::map<std::string,std::string>& metadata;
        decltype(KFState::exactStateTransformCallback)& callback;bool committed=false;
        ~Restore(){state.exactStateTransformCallback.swap(callback);if(!committed)state.metaDataMap.swap(metadata);}
    } restore{state,sourceMetadata,callback};
    state.exactStateTransformCallback=[&,original=callback](const KFState& source,GTime time,
        const std::map<KFKey,int>& sourceOrder,const std::map<KFKey,int>& destinationOrder,
        const SparseMatrix<double>& t,const std::string& label,const VectorXd& mean,
        const MatrixXd& covariance,std::uint64_t before,std::uint64_t after)
    {
        MatrixXd sourceRows=MatrixXd::Zero(all.size(),source.x.size());
        MatrixXd destinationRows=MatrixXd::Zero(all.size(),mean.size());
        int targetIndex=0;
        for(const auto& [id,target]:all)
        {
            for(const auto& [key,value]:target.numericRow)sourceRows(targetIndex,sourceOrder.at(key))=value;
            for(const auto& [key,value]:hints.at(id))destinationRows(targetIndex,destinationOrder.at(key))=value;
            ++targetIndex;
        }
        const VectorXd expectedMean=sourceRows*source.x,actualMean=destinationRows*mean;
        const MatrixXd expectedCross=(sourceRows*source.P)*t.transpose();
        const MatrixXd actualCross=destinationRows*covariance;
        const double meanError=(expectedMean-actualMean).norm()/std::max(1.0,expectedMean.norm());
        const double covarianceError=(expectedCross-actualCross).norm()/std::max(1.0,expectedCross.norm());
        if(!std::isfinite(meanError) || !std::isfinite(covarianceError) || meanError>1e-10 || covarianceError>1e-10)
        {state.lastFactorTransactionFailureReason="FLOAT_TARGET_COMPRESSION_COMPLETE_JOINT_MOMENT_AUDIT_FAILED";return false;}
        trace<<"\nZHANG_FLOAT_TARGET_COMPRESSION_JOINT_AUDIT targets="<<all.size()
            <<" destination_states="<<mean.size()<<" mean_relative_error="<<meanError
            <<" full_cross_covariance_relative_error="<<covarianceError<<" diagonal_only=0";
        return !original || original(source,time,sourceOrder,destinationOrder,t,label,mean,covariance,before,after);
    };
    if(!state.applyStateTransform(trace,transform,ZHANG_FLOAT_TARGET_COMPRESSION_LABEL,{},&workingMetadata,&hints))
    {
        const auto failure=state.lastFactorTransactionFailureReason;
        state.lastFactorTransactionFailureReason=oldFailure;
        return skip("OPTIONAL_STRICT_TRANSFORM_DECLINED:"+failure);
    }
    graphProductEventReceipts.swap(stagedReceipts);
    restore.committed=true;
    reason="GLOBAL_WHOLE_FLOAT_TARGET_SPACE_COMPRESSED";
    trace<<"\nZHANG_FLOAT_TARGET_COMPRESSION status=COMMITTED source_states="<<plan.sourceDimension
        <<" destination_states="<<state.x.size()<<" carriers="<<plan.carrierRank
        <<" active_targets="<<all.size()<<" graph_system_transforms_merged=0 transient_capture_peak_removed=0"
        <<" full_cross_covariance_preserved=1 integer_domain_unchanged=1";
    return true;
}

void zhangClearGraphProductHandoffTargets(KFState& state)
{
    std::vector<ZhangProductFloatTargetRequest> targets;std::string reason;
    if(!zhangProductFloatTargetMetadataValid(state,reason) ||
        !zhangReadFloatTargets(state.metaDataMap,targets,reason))
        throw std::runtime_error("FLOAT_TARGET_DIRECTORY_INVALID_BEFORE_REFRESH:"+reason);
    auto metadata=state.metaDataMap;
    for(auto& target:targets)
    {
        target.active=false;
        metadata[std::string(ZHANG_FLOAT_TARGET_V2_PREFIX)+target.targetId]=zhangEncodeFloatTarget(target);
    }
    state.metaDataMap.swap(metadata);
}

bool zhangProductFunctionalEventDiagnostic(
    const KFState&                         kfState,
    E_Sys                                  system,
    E_ObsCode                              observable,
    const vector<ZhangGraphEdge>&          oldEdges,
    const vector<int>&                     oldArcVersions,
    const vector<ZhangGraphEdge>&          newEdges,
    const vector<int>&                     newArcVersions,
    ZhangProductFunctionalEventDiagnostic& diagnostic
)
{
    diagnostic = {};
    const string runtimeId = zhangGraphRuntimeId(kfState);
    auto stateIt = graphStateMap.find({runtimeId, system});
    if (stateIt == graphStateMap.end()
        || oldEdges.size() != oldArcVersions.size()
        || newEdges.size() != newArcVersions.size())
    {
        return false;
    }

    const auto& runtime = stateIt->second;
    auto observationHistory = [&](const pair<ZhangGraphEdge, int>& arc)
        -> const GraphRuntimeState::ArcVersionObservationHistory*
    {
        if (observable != E_ObsCode::NONE)
        {
            const auto signal =
                runtime.signalArcVersionObservationHistory.find(
                    {arc.first, int(observable), arc.second});
            return signal == runtime.signalArcVersionObservationHistory.end()
                ? nullptr : &signal->second;
        }
        const auto legacy = runtime.arcVersionObservationHistory.find(arc);
        return legacy == runtime.arcVersionObservationHistory.end()
            ? nullptr : &legacy->second;
    };
    diagnostic.eventCause = runtime.lastProductEventCause;
    using VersionedEdge = pair<ZhangGraphEdge, int>;
    set<VersionedEdge> oldSet;
    set<VersionedEdge> newSet;
    for (size_t index = 0; index < oldEdges.size(); index++)
    {
        oldSet.emplace(oldEdges[index], oldArcVersions[index]);
    }
    for (size_t index = 0; index < newEdges.size(); index++)
    {
        newSet.emplace(newEdges[index], newArcVersions[index]);
    }

    vector<VersionedEdge> removed;
    vector<VersionedEdge> introduced;
    std::set_difference(
        oldSet.begin(), oldSet.end(), newSet.begin(), newSet.end(),
        back_inserter(removed));
    std::set_difference(
        newSet.begin(), newSet.end(), oldSet.begin(), oldSet.end(),
        back_inserter(introduced));

    auto appendDiagnostic = [&](const VersionedEdge& versioned,
                                vector<ZhangPhysicalArcDiagnostic>& output)
    {
        ZhangPhysicalArcDiagnostic item;
        item.edge = versioned.first;
        item.observable = observable;
        item.arcVersion = versioned.second;
        const auto* history = observationHistory(versioned);
        if (history)
        {
            item.observationEpochs = history->observationEpochs;
            if (history->firstObservedEpoch >= 0
                && history->lastObservedEpoch >=
                    history->firstObservedEpoch)
            {
                item.ageEpochs = history->lastObservedEpoch
                    - history->firstObservedEpoch + 1;
            }
        }
        output.push_back(std::move(item));
    };
    for (const auto& versioned : removed)
    {
        appendDiagnostic(versioned, diagnostic.oldSupport);
    }
    for (const auto& versioned : introduced)
    {
        appendDiagnostic(versioned, diagnostic.newSupport);
    }

    vector<int> commonEpochs;
    bool initialised = false;
    for (const auto& versioned : removed)
    {
        const auto* history = observationHistory(versioned);
        if (!history)
        {
            commonEpochs.clear();
            initialised = true;
            break;
        }
        if (!initialised)
        {
            commonEpochs = history->observedEpochs;
            initialised = true;
        }
        else
        {
            vector<int> intersection;
            std::set_intersection(
                commonEpochs.begin(), commonEpochs.end(),
                history->observedEpochs.begin(),
                history->observedEpochs.end(),
                back_inserter(intersection));
            commonEpochs = std::move(intersection);
        }
    }
    for (const auto& versioned : introduced)
    {
        const auto* history = observationHistory(versioned);
        if (!history)
        {
            commonEpochs.clear();
            initialised = true;
            break;
        }
        if (!initialised)
        {
            commonEpochs = history->observedEpochs;
            initialised = true;
        }
        else
        {
            vector<int> intersection;
            std::set_intersection(
                commonEpochs.begin(), commonEpochs.end(),
                history->observedEpochs.begin(),
                history->observedEpochs.end(),
                back_inserter(intersection));
            commonEpochs = std::move(intersection);
        }
    }
    diagnostic.commonObservationEpochs = initialised
        ? static_cast<int>(commonEpochs.size())
        : 0;
    return true;
}

bool zhangProductFunctionalEventDiagnostic(
    const KFState& kfState,
    E_Sys system,
    const vector<ZhangGraphEdge>& oldEdges,
    const vector<int>& oldArcVersions,
    const vector<ZhangGraphEdge>& newEdges,
    const vector<int>& newArcVersions,
    ZhangProductFunctionalEventDiagnostic& diagnostic)
{
    return zhangProductFunctionalEventDiagnostic(
        kfState, system, E_ObsCode::NONE, oldEdges, oldArcVersions,
        newEdges, newArcVersions, diagnostic);
}

static bool cloneZhangGraphRuntimeBase(
    const KFState& source,
    const KFState& destination,
    std::string* failureReason
)
{
    const string sourceRuntimeId = zhangGraphRuntimeId(source);
    const string destinationRuntimeId = zhangGraphRuntimeId(destination);
    auto reject = [&](const string& reason)
    {
        if (failureReason) *failureReason = reason;
        BOOST_LOG_TRIVIAL(error)
            << "ZHANG_GRAPH_CLONE_REJECT reason=" << reason
            << " source=" << sourceRuntimeId << " destination=" << destinationRuntimeId;
        return false;
    };
    if (!validZhangGraphRuntimeId(sourceRuntimeId) ||
        !validZhangGraphRuntimeId(destinationRuntimeId))
        return reject("GRAPH_RUNTIME_ID_INVALID");
    if (sourceRuntimeId == destinationRuntimeId)
        return reject("GRAPH_SOURCE_DESTINATION_RUNTIME_ALIAS");

    const auto configuredSystems = zhangConfiguredSpanningTreeSystems();
    bool estimatedGraphState = false;
    for (const auto system : configuredSystems)
    {
        if (!zhangHasEstimatedGraphPhaseState(source, system,
                acsConfig.zhangFullRank.sysOpts.at(system).baseline_observables)) continue;
        estimatedGraphState = true;
        const auto graph = graphStateMap.find({sourceRuntimeId, system});
        if (graph == graphStateMap.end() || !graph->second.initialized)
            return reject("GRAPH_SOURCE_RUNTIME_UNAVAILABLE");
        if (graph->second.productBasisMode != configuredProductBasisMode())
            return reject("GRAPH_SOURCE_PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION");
    }
    if (configuredSystems.empty())
    {
        if (failureReason) *failureReason = "NO_GRAPH_CONFIGURED";
        return true;
    }

    // Allocate all replacement nodes before removing any destination data.
    // merge() then transfers these preallocated map nodes without copying.
    map<ZhangGraphRuntimeKey, GraphRuntimeState> copies;
    bool initializedSource = false;
    try
    {
        for (const auto& [identity, runtime] : graphStateMap)
        {
            if (identity.first != sourceRuntimeId) continue;
            if (runtime.initialized && runtime.productBasisMode != configuredProductBasisMode())
                return reject("GRAPH_SOURCE_PRODUCT_BASIS_MODE_MISMATCH");
            initializedSource |= runtime.initialized;
            copies.emplace(ZhangGraphRuntimeKey{destinationRuntimeId, identity.second}, runtime);
        }
    }
    catch (const std::exception&)
    {
        return reject("GRAPH_RUNTIME_CLONE_ALLOCATION_FAILED");
    }
    if (copies.empty() || !initializedSource)
    {
        if (estimatedGraphState) return reject("GRAPH_SOURCE_RUNTIME_UNAVAILABLE");
        if (failureReason) *failureReason = "NO_ESTIMATED_GRAPH_STATE";
        return true;
    }
    // The selected integer view is not part of the mutable graph runtime or
    // checkpoint.  A disposable branch must nevertheless query the same
    // invocation snapshot as its authoritative capture owner.
    map<ZhangGraphRuntimeKey, vector<SelectedGraphContext>> selectedCopies;
    try
    {
        std::lock_guard<std::mutex> lock(selectedGraphContextMutex);
        for (const auto& [identity, selections] : selectedGraphContexts)
            if (identity.first == sourceRuntimeId)
                selectedCopies.emplace(ZhangGraphRuntimeKey{destinationRuntimeId, identity.second}, selections);
    }
    catch (const std::exception&)
    { return reject("GRAPH_SELECTED_CONTEXT_CLONE_ALLOCATION_FAILED"); }
    eraseZhangGraphRuntime(destination);
    graphStateMap.merge(copies);
    {
        std::lock_guard<std::mutex> lock(selectedGraphContextMutex);
        selectedGraphContexts.merge(selectedCopies);
    }
    if (failureReason) *failureReason = "NONE";
    return true;
}

bool cloneZhangGraphRuntime(const KFState& source, const KFState& destination,
    std::string* failureReason)
{
    if (hasZhangFloatTargetV2State(source.metaDataMap) || hasZhangFloatTargetV2State(destination.metaDataMap))
    {
        if (failureReason) *failureReason = "GRAPH_FLOAT_TARGET_CLONE_REQUIRES_MUTABLE_DESTINATION";
        return false;
    }
    return cloneZhangGraphRuntimeBase(source, destination, failureReason);
}

bool cloneZhangGraphRuntime(const KFState& source, KFState& destination,
    std::string* failureReason)
{
    auto reject = [&](const string& reason)
    { if (failureReason) *failureReason = reason; return false; };
    const bool sourceV2 = hasZhangFloatTargetV2State(source.metaDataMap);
    const bool destinationV2 = hasZhangFloatTargetV2State(destination.metaDataMap);
    if (!sourceV2 && !destinationV2)
        return cloneZhangGraphRuntimeBase(source, destination, failureReason);
    if (!sourceV2 || !destinationV2)
        return reject("GRAPH_FLOAT_TARGET_CLONE_DIRECTORY_NOT_IDENTICAL_COPY");
    string reason;
    if (!zhangProductFloatTargetMetadataValid(source, reason))
        return reject("GRAPH_FLOAT_TARGET_CLONE_SOURCE_INVALID:" + reason);
    const auto sourceOwner = zhangGraphRuntimeId(source), destinationOwner = zhangGraphRuntimeId(destination);
    if (!validZhangGraphRuntimeId(destinationOwner) || sourceOwner == destinationOwner)
        return reject("GRAPH_FLOAT_TARGET_CLONE_DESTINATION_IDENTITY_INVALID");
    auto rawEqual = [](const auto& a, const auto& b)
    { return a.rows() == b.rows() && a.cols() == b.cols() &&
        (a.size() == 0 || std::memcmp(a.data(), b.data(), std::size_t(a.size()) * sizeof(double)) == 0); };
    if (!rawEqual(source.x, destination.x) || !rawEqual(source.dx, destination.dx) ||
        !rawEqual(source.P, destination.P) || source.factorCommitSequence != destination.factorCommitSequence ||
        captureZhangCheckpointTime(source.time).bigTimeBytes != captureZhangCheckpointTime(destination.time).bigTimeBytes ||
        source.kfIndexMap.size() != destination.kfIndexMap.size())
        return reject("GRAPH_FLOAT_TARGET_CLONE_NOT_IDENTICAL_UNCONDITIONED_POSTERIOR");
    auto left = source.kfIndexMap.cbegin(), right = destination.kfIndexMap.cbegin();
    for (; left != source.kfIndexMap.cend(); ++left, ++right)
    {
        const auto& a = left->first; const auto& b = right->first;
        if (!(a == b) || left->second != right->second || a.comment != b.comment || a.rec_ptr != b.rec_ptr ||
            captureZhangCheckpointTime(a.estimatedTime).bigTimeBytes != captureZhangCheckpointTime(b.estimatedTime).bigTimeBytes)
            return reject("GRAPH_FLOAT_TARGET_CLONE_FULL_KEY_ORDER_MISMATCH");
    }
    auto v2Key = [](const string& key)
    { return key.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || key.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX) ||
        key == ZHANG_FLOAT_TARGET_SCOPE_V2 || key == ZHANG_FLOAT_TARGET_GENERATION_V2 ||
        key == ZHANG_FLOAT_TARGET_SAVED_FACTOR_SEQUENCE_V2; };
    for (const auto& [key, value] : source.metaDataMap) if (v2Key(key))
    {
        const auto found = destination.metaDataMap.find(key);
        if (found == destination.metaDataMap.end() || found->second != value)
            return reject("GRAPH_FLOAT_TARGET_CLONE_DIRECTORY_NOT_IDENTICAL_COPY");
    }
    for (const auto& [key, value] : destination.metaDataMap)
        if (v2Key(key) && !source.metaDataMap.contains(key))
            return reject("GRAPH_FLOAT_TARGET_CLONE_UNRELATED_DESTINATION_DIRECTORY");
    map<string, string> staged;
    string successReason;
    try
    {
        staged = destination.metaDataMap;
        if (failureReason) successReason = "GRAPH_AND_WHOLE_FLOAT_TARGET_CLONE_COMMITTED";
        for (auto& [key, value] : staged)
            if (key.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) || key.starts_with(ZHANG_FLOAT_CARRIER_V2_PREFIX))
            {
                ZhangProductFloatTargetRequest target;
                if (!zhangDecodeFloatTarget(value, target) || target.scope.runtimeId != sourceOwner)
                    return reject("GRAPH_FLOAT_TARGET_CLONE_SOURCE_RECORD_INVALID");
                target.scope.runtimeId = destinationOwner;
                value = zhangEncodeFloatTarget(target);
            }
        std::uint64_t saved = 0;
        if (!readZhangFloatTargetSavedFactorSequence(staged, destinationOwner,
                zhangFloatCarrierStateOrder(destination.kfIndexMap), saved, reason) || saved != destination.factorCommitSequence)
            return reject("GRAPH_FLOAT_TARGET_CLONE_STAGED_BINDING_INVALID:" + reason);
        for (const auto& [key, index] : destination.kfIndexMap)
            if (zhangIsProductFloatCarrierKey(key) && !zhangProductFloatCarrierCoordinateValid(destination, staged, key, reason))
                return reject("GRAPH_FLOAT_TARGET_CLONE_DESTINATION_CARRIER_INVALID:" + reason);
    }
    catch (const std::exception&)
    { return reject("GRAPH_FLOAT_TARGET_CLONE_METADATA_ALLOCATION_FAILED"); }
    if (!cloneZhangGraphRuntimeBase(source, destination, failureReason)) return false;
    // History remains owned by the authoritative FLOAT source. Only V2
    // non-integer numeric target/carrier records are rebound to this lease.
    static_assert(noexcept(destination.metaDataMap.swap(staged)));
    destination.metaDataMap.swap(staged);
    if (failureReason) failureReason->swap(successReason);
    return true;
}

bool captureZhangFloatGraphAuthorityToken(
    const KFState& state,
    ZhangFloatGraphAuthorityToken& token,
    std::string& failureReason)
{
    token = {};
    const string runtimeId = zhangGraphRuntimeId(state);
    if (!validZhangGraphRuntimeId(runtimeId))
    {
        failureReason = "FLOAT_GRAPH_RUNTIME_ID_INVALID";
        return false;
    }
    const auto configuredSystems = zhangConfiguredSpanningTreeSystems();
    if (configuredSystems.empty())
    {
        token.valid = true;
        token.runtimeIdentity = runtimeId;
        token.physicalState = "NO_GRAPH_CONFIGURED";
        failureReason = "NO_GRAPH_CONFIGURED";
        return true;
    }
    bool estimatedGraphState = false;
    for (const auto system : configuredSystems)
    {
        if (!zhangHasEstimatedGraphPhaseState(state, system,
                acsConfig.zhangFullRank.sysOpts.at(system).baseline_observables)) continue;
        estimatedGraphState = true;
        const auto graph = graphStateMap.find({runtimeId, system});
        if (graph == graphStateMap.end() || !graph->second.initialized)
        {
            failureReason = "FLOAT_GRAPH_RUNTIME_UNAVAILABLE";
            return false;
        }
    }
    std::ostringstream output;
    output.imbue(std::locale::classic());
    auto text = [&](const string& value) { output << value.size() << ':' << value << ';'; };
    auto satellite = [&](const SatSys& value) { output << int(value.sys) << ',' << value.prn << ';'; };
    auto edge = [&](const ZhangGraphEdge& value) { text(value.receiver); satellite(value.satellite); };
    auto edges = [&](const set<ZhangGraphEdge>& values)
    {
        output << values.size() << ';';
        for (const auto& value : values) edge(value);
    };
    auto basis = [&](const ZhangGraphBasis& value)
    {
        text(value.rootReceiver); edges(value.edges); edges(value.treeEdges);
        output << value.receivers.size() << ';';
        for (const auto& receiver : value.receivers) text(receiver);
        output << value.satellites.size() << ';';
        for (const auto& sat : value.satellites) satellite(sat);
        output << value.componentCount << ';' << value.connected << ';';
    };
    auto signalRows = [&](const map<ZhangGraphEdge, map<int, int>>& values)
    {
        output << values.size() << ';';
        for (const auto& [physicalEdge, signals] : values)
        {
            edge(physicalEdge); output << signals.size() << ';';
            for (const auto& [signal, version] : signals) output << signal << ',' << version << ';';
        }
    };
    text(runtimeId);
    output << std::hexfloat << std::setprecision(std::numeric_limits<long double>::max_digits10)
           << state.time.bigTime << ';';
    bool initializedRuntime = false;
    for (const auto& [identity, runtime] : graphStateMap)
    {
        if (identity.first != runtimeId) continue;
        if (runtime.initialized && runtime.productBasisMode != configuredProductBasisMode())
        {
            failureReason = "FLOAT_GRAPH_PRODUCT_BASIS_MODE_CHANGE_REQUIRES_EXPLICIT_MIGRATION";
            return false;
        }
        initializedRuntime |= runtime.initialized;
        output << int(identity.second) << ';' << runtime.initialized << ';'
               << runtime.datumVersion << ';' << runtime.representationVersion << ';'
               << runtime.floatGaugeVersion << ';' << runtime.eventCounter << ';'
               << runtime.epochIndex << ';' << runtime.deferredEpochs << ';';
        text(runtime.productBasisMode);
        basis(runtime.basis); basis(runtime.activeBasis);
        // A native product view is part of the S authority contract, not an
        // autonomous selector an AR branch may mutate without detection.
        if (runtime.productBasisMode == "SINGLE_S_BASIS_NATIVE")
        {
            basis(runtime.productBasis);
            output << runtime.productDatumVersion << ';'
                   << runtime.committedSGaugeVersion << ';'
                   << runtime.productInitialized << ';'
                   << runtime.productAuthorityComponentVersion << ';';
        }
        edges(runtime.observationEdges); edges(runtime.stateEdges);
        output << runtime.edgeHistory.size() << ';';
        for (const auto& [physicalEdge, history] : runtime.edgeHistory)
        {
            edge(physicalEdge);
            output << history.arcVersion << ';' << history.continuousEpochs << ';'
                   << history.outageEpochs << ';';
        }
        signalRows(runtime.phaseArcVersions); signalRows(runtime.phaseAdmissions);
        output << runtime.signalObservationEdges.size() << ';';
        for (const auto& [signal, physicalEdges] : runtime.signalObservationEdges)
        {
            output << signal << ';'; edges(physicalEdges);
        }
        output << runtime.signalObservationFinalized << ';'
               << runtime.signalObservationFactorSequence << ';';
        text(runtime.phaseAdmissionTime);
        output << runtime.phaseAdmissionChart << ';';
    }
    if (!initializedRuntime && estimatedGraphState)
    {
        failureReason = "FLOAT_GRAPH_RUNTIME_UNAVAILABLE";
        return false;
    }
    token.runtimeIdentity = runtimeId;
    token.physicalState = output.str();
    token.valid = true;
    failureReason = "NONE";
    return true;
}

bool zhangFloatGraphAuthorityUnchanged(
    const KFState& state,
    const ZhangFloatGraphAuthorityToken& token,
    std::string& failureReason)
{
    if (!token.valid)
    {
        failureReason = "FLOAT_GRAPH_AUTHORITY_TOKEN_INVALID";
        return false;
    }
    ZhangFloatGraphAuthorityToken current;
    if (!captureZhangFloatGraphAuthorityToken(state, current, failureReason)) return false;
    if (current.runtimeIdentity != token.runtimeIdentity || current.physicalState != token.physicalState)
    {
        failureReason = "FLOAT_GRAPH_TOPOLOGY_OR_PHYSICAL_ARC_CHANGED";
        return false;
    }
    failureReason = "NONE";
    return true;
}

void eraseZhangGraphRuntime(const KFState& state)
{
    const string runtimeId = zhangGraphRuntimeId(state);
    if (runtimeId.empty())
    {
        return;
    }
    {
        std::lock_guard<std::mutex> lock(selectedGraphContextMutex);
        for (auto it = selectedGraphContexts.begin(); it != selectedGraphContexts.end();)
            it = it->first.first == runtimeId
                ? selectedGraphContexts.erase(it) : std::next(it);
    }
    for (auto it = graphStateMap.begin(); it != graphStateMap.end();)
    {
        if (it->first.first == runtimeId)
        {
            it = graphStateMap.erase(it);
        }
        else
        {
            ++it;
        }
    }
    for (auto it = graphProductEventReceipts.begin();
         it != graphProductEventReceipts.end();)
        it = std::get<0>(it->first) == runtimeId
            ? graphProductEventReceipts.erase(it) : std::next(it);
    for (auto it = graphEventPreArcVersions.begin();
         it != graphEventPreArcVersions.end();)
        it = it->first.first == runtimeId
            ? graphEventPreArcVersions.erase(it) : std::next(it);
    for (auto it = graphEventPreSignalArcVersions.begin();
         it != graphEventPreSignalArcVersions.end();)
        it = it->first.first == runtimeId
            ? graphEventPreSignalArcVersions.erase(it) : std::next(it);
}

bool exportZhangGraphCheckpointSection(
    const KFState&      state,
    const std::string& runtimeId,
    std::string&       payload,
    std::string&       failureReason
)
{
    payload.clear();
    failureReason.clear();
    if (runtimeId.empty() || zhangGraphRuntimeId(state) != runtimeId)
    {
        failureReason = runtimeId.empty()
            ? "ZHANG_GRAPH_CHECKPOINT_EMPTY_RUNTIME_ID"
            : "ZHANG_GRAPH_CHECKPOINT_RUNTIME_ID_NOT_BOUND";
        return false;
    }
    // Auxiliary identities live in the core KF metadata, whose archive already
    // includes metaDataMap. Keep schema-3/4 graph bytes backward compatible;
    // validate the complete core/registry bijection instead of rejecting every
    // correlated retained FLOAT row.
    if (!zhangHandoffFloatRegistryValid(state, failureReason)) return false;
    // A validator success label is not the failure of the next independent
    // graph check.  In particular it must not mask a broken product arc key.
    failureReason.clear();

    ZhangGraphCheckpointPayloadV4 snapshot;
    snapshot.runtimeId = runtimeId;

    for (const auto& [identity, outage] : outageStateMap)
    {
        if (identity.first != runtimeId)
        {
            continue;
        }
        snapshot.outageStates.push_back({
            static_cast<int>(identity.second),
            outage.receiverEpochs,
            outage.satelliteEpochs});
    }

    for (const auto& [identity, runtime] : graphStateMap)
    {
        if (identity.first != runtimeId)
        {
            continue;
        }

        ZhangGraphCheckpointRuntimeV4 versionedOutput;
        if (runtime.productBasisMode != configuredProductBasisMode())
        {
            failureReason = "ZHANG_GRAPH_CHECKPOINT_RUNTIME_PRODUCT_BASIS_MODE_MISMATCH";
            return false;
        }
        if (singleSBasisNative() && runtime.productDatumVersion != runtime.committedSGaugeVersion)
        {
            failureReason = "ZHANG_SINGLE_S_CHECKPOINT_PRODUCT_NOT_BOUND_TO_S_GAUGE";
            return false;
        }
        auto& output = versionedOutput.common;
        output.system = static_cast<int>(identity.second);
        output.basis = zhangGraphCheckpointBasis(runtime.basis);
        output.activeBasis = zhangGraphCheckpointBasis(runtime.activeBasis);
        output.productBasis = zhangGraphCheckpointBasis(runtime.productBasis);
        for (const auto& [edge, arcVersion] : runtime.productArcVersions)
        {
            output.productArcVersions[
                zhangGraphCheckpointEdge(edge)] = arcVersion;
        }
        output.productCoreReceivers = runtime.productCoreReceivers;
        for (const auto& edge : runtime.observationEdges)
        {
            output.observationEdges.insert(
                zhangGraphCheckpointEdge(edge));
        }
        for (const auto& edge : runtime.stateEdges)
        {
            output.stateEdges.insert(zhangGraphCheckpointEdge(edge));
        }
        for (const auto& [edge, history] : runtime.edgeHistory)
        {
            output.edgeHistory[zhangGraphCheckpointEdge(edge)] = {
                history.continuousEpochs,
                history.outageEpochs,
                history.arcVersion};
        }
        for (const auto& [versionedEdge, history] :
             runtime.arcVersionObservationHistory)
        {
            output.arcVersionObservationHistory[
                {zhangGraphCheckpointEdge(versionedEdge.first),
                 versionedEdge.second}] = {
                    history.firstObservedEpoch,
                    history.lastObservedEpoch,
                    history.observationEpochs,
                    history.observedEpochs};
        }
        output.epochIndex = runtime.epochIndex;
        output.lastProductEventCause = runtime.lastProductEventCause;
        output.initialized = runtime.initialized;
        output.deferredEpochs = runtime.deferredEpochs;
        output.datumVersion = runtime.datumVersion;
        output.representationVersion = runtime.representationVersion;
        output.floatGaugeVersion = runtime.floatGaugeVersion;
        output.integerComponentVersion = runtime.integerComponentVersion;
        output.eventCounter = runtime.eventCounter;
        output.productDatumVersion = runtime.productDatumVersion;
        output.productInitialized = runtime.productInitialized;
        output.productAuthorityComponentVersion = runtime.productAuthorityComponentVersion;
        output.phaseAdmissionTime = runtime.phaseAdmissionTime;
        output.phaseAdmissionChart = runtime.phaseAdmissionChart;
        for (const auto& [edge, signals] : runtime.phaseArcVersions)
            output.phaseArcVersions[zhangGraphCheckpointEdge(edge)] = signals;
        for (const auto& [edge, signals] : runtime.phaseAdmissions)
            output.phaseAdmissions[zhangGraphCheckpointEdge(edge)] = signals;
        for (const auto& [edge, signals] :
             runtime.productSignalArcVersions)
            versionedOutput.productSignalArcVersions[
                zhangGraphCheckpointEdge(edge)] = signals;
        for (const auto& [identity, history] :
             runtime.signalArcVersionObservationHistory)
            versionedOutput.signalArcVersionObservationHistory[
                {zhangGraphCheckpointEdge(std::get<0>(identity)),
                 {std::get<1>(identity), std::get<2>(identity)}}] = {
                    history.firstObservedEpoch,
                    history.lastObservedEpoch,
                    history.observationEpochs,
                    history.observedEpochs};
        for (const auto& [code, edges] :
             runtime.signalObservationEdges)
            for (const auto& edge : edges)
                versionedOutput.signalObservationEdges[code].insert(
                    zhangGraphCheckpointEdge(edge));
        versionedOutput.requiresFreshProductIdentity =
            runtime.requiresFreshProductIdentity;
        versionedOutput.signalIdentityAuthoritative =
            singleSBasisNative() || zhangR51SingleBasisHandoffEnabled();
        snapshot.graphStates.push_back(std::move(versionedOutput));
    }

    auto bySystem = [](const auto& left, const auto& right)
    {
        return left.system < right.system;
    };
    std::sort(
        snapshot.outageStates.begin(),
        snapshot.outageStates.end(),
        bySystem);
    std::sort(
        snapshot.graphStates.begin(),
        snapshot.graphStates.end(),
        [](const auto& left, const auto& right)
        {
            return left.common.system < right.common.system;
        });

    if (!zhangGraphCheckpointPayloadValid(
            snapshot, runtimeId, failureReason, singleSBasisNative()))
    {
        return false;
    }
    string archive;
    ZhangGraphCheckpointPayloadV5 envelope;
    envelope.productBasisMode = configuredProductBasisMode();
    if (singleSBasisNative())
        for (const auto& [identity, runtime] : graphStateMap)
            if (identity.first == runtimeId)
                envelope.committedSGaugeVersions.emplace(int(identity.second), runtime.committedSGaugeVersion);
    envelope.graph = std::move(snapshot);
    if (!zhangGraphCheckpointSerialize(
            envelope, archive, failureReason))
        return false;
    payload = string(ZHANG_GRAPH_V5_PAYLOAD_PREFIX) + archive;
    return true;
}

bool validateZhangGraphCheckpointSection(
    const std::string& runtimeId,
    const std::string& payload,
    std::string&       failureReason
)
{
    ZhangGraphCheckpointPayloadV4 snapshot;
    bool legacyIdentityReset = false;
    return zhangGraphCheckpointDecodeAndValidate(
        runtimeId, payload, snapshot, legacyIdentityReset,
        failureReason);
}

bool importZhangGraphCheckpointSection(
    KFState&           state,
    const std::string& runtimeId,
    const std::string& payload,
    std::string&       failureReason
)
{
    if (zhangGraphRuntimeId(state) != runtimeId)
    {
        failureReason = "ZHANG_GRAPH_CHECKPOINT_RUNTIME_ID_NOT_BOUND";
        return false;
    }
    if (!zhangHandoffFloatRegistryValid(state, failureReason)) return false;
    ZhangGraphCheckpointPayloadV4 snapshot;
    bool legacyIdentityReset = false;
    if (!zhangGraphCheckpointDecodeAndValidate(
            runtimeId, payload, snapshot, legacyIdentityReset,
            failureReason))
    {
        return false;
    }

    vector<pair<E_Sys, ReferenceOutageState>> restoredOutages;
    vector<pair<E_Sys, GraphRuntimeState>> restoredGraphs;
    try
    {
        restoredOutages.reserve(snapshot.outageStates.size());
        for (const auto& input : snapshot.outageStates)
        {
            ReferenceOutageState output;
            output.receiverEpochs = input.receiverEpochs;
            output.satelliteEpochs = input.satelliteEpochs;
            restoredOutages.emplace_back(
                static_cast<E_Sys>(input.system), output);
        }

        restoredGraphs.reserve(snapshot.graphStates.size());
        for (const auto& versionedInput : snapshot.graphStates)
        {
            const auto& input = versionedInput.common;
            GraphRuntimeState output;
            output.productBasisMode = configuredProductBasisMode();
            output.basis = zhangGraphCheckpointBasis(input.basis);
            output.activeBasis =
                zhangGraphCheckpointBasis(input.activeBasis);
            output.productBasis =
                zhangGraphCheckpointBasis(input.productBasis);
            for (const auto& [edge, arcVersion] : input.productArcVersions)
            {
                output.productArcVersions[
                    zhangGraphCheckpointEdge(edge)] = arcVersion;
            }
            output.productCoreReceivers = input.productCoreReceivers;
            for (const auto& edge : input.observationEdges)
            {
                output.observationEdges.insert(
                    zhangGraphCheckpointEdge(edge));
            }
            for (const auto& edge : input.stateEdges)
            {
                output.stateEdges.insert(zhangGraphCheckpointEdge(edge));
            }
            for (const auto& [edge, history] : input.edgeHistory)
            {
                auto& target =
                    output.edgeHistory[zhangGraphCheckpointEdge(edge)];
                target.continuousEpochs = history.continuousEpochs;
                target.outageEpochs = history.outageEpochs;
                target.arcVersion = history.arcVersion;
            }
            for (const auto& [versionedEdge, history] :
                 input.arcVersionObservationHistory)
            {
                auto& target = output.arcVersionObservationHistory[
                    {zhangGraphCheckpointEdge(versionedEdge.first),
                     versionedEdge.second}];
                target.firstObservedEpoch = history.firstObservedEpoch;
                target.lastObservedEpoch = history.lastObservedEpoch;
                target.observationEpochs = history.observationEpochs;
                target.observedEpochs = history.observedEpochs;
            }
            output.epochIndex = input.epochIndex;
            output.lastProductEventCause = input.lastProductEventCause;
            output.initialized = input.initialized;
            output.deferredEpochs = input.deferredEpochs;
            output.datumVersion = input.datumVersion;
            output.representationVersion = input.representationVersion;
            output.floatGaugeVersion = input.floatGaugeVersion;
            output.integerComponentVersion =
                input.integerComponentVersion;
            output.eventCounter = input.eventCounter;
            output.productDatumVersion = input.productDatumVersion;
            output.committedSGaugeVersion = singleSBasisNative() ? input.productDatumVersion : 0;
            output.productInitialized = input.productInitialized;
            output.productAuthorityComponentVersion = input.productAuthorityComponentVersion;
            output.phaseAdmissionTime = input.phaseAdmissionTime;
            output.phaseAdmissionChart = input.phaseAdmissionChart;
            for (const auto& [edge, signals] : input.phaseArcVersions)
                output.phaseArcVersions[zhangGraphCheckpointEdge(edge)] = signals;
            for (const auto& [edge, signals] : input.phaseAdmissions)
                output.phaseAdmissions[zhangGraphCheckpointEdge(edge)] = signals;
            for (const auto& [edge, signals] :
                 versionedInput.productSignalArcVersions)
                output.productSignalArcVersions[
                    zhangGraphCheckpointEdge(edge)] = signals;
            for (const auto& [identity, history] :
                 versionedInput.signalArcVersionObservationHistory)
            {
                auto& target = output.signalArcVersionObservationHistory[
                    {zhangGraphCheckpointEdge(identity.first),
                     identity.second.first, identity.second.second}];
                target.firstObservedEpoch = history.firstObservedEpoch;
                target.lastObservedEpoch = history.lastObservedEpoch;
                target.observationEpochs = history.observationEpochs;
                target.observedEpochs = history.observedEpochs;
            }
            for (const auto& [code, edges] :
                 versionedInput.signalObservationEdges)
                for (const auto& edge : edges)
                    output.signalObservationEdges[code].insert(
                        zhangGraphCheckpointEdge(edge));
            output.requiresFreshProductIdentity =
                versionedInput.requiresFreshProductIdentity;
            if (legacyIdentityReset)
            {
                // Schema 3 did not bind integer/product rows to independent
                // signal arcs. Preserve FLOAT x/P and the S topology, but
                // discard only authority that cannot be reconstructed.
                output.productBasis = {};
                output.productArcVersions.clear();
                output.productSignalArcVersions.clear();
                output.productCoreReceivers.clear();
                output.productInitialized = false;
                output.productAuthorityComponentVersion = -1;
                output.productDatumVersion++;
                output.integerComponentVersion++;
                output.phaseArcVersions.clear();
                output.phaseAdmissions.clear();
                output.phaseAdmissionTime.clear();
                output.phaseAdmissionChart = -1;
                output.signalObservationEdges.clear();
                output.signalArcVersionObservationHistory.clear();
                output.requiresFreshProductIdentity = true;
                output.lastProductEventCause =
                    "LEGACY_SIGNAL_IDENTITY_REBUILT_FROM_FLOAT";
            }
            restoredGraphs.emplace_back(
                static_cast<E_Sys>(input.system), std::move(output));
        }

        // Build complete replacement maps first.  No live state is modified
        // until every allocation/conversion above and below has succeeded.
        auto replacementOutages = outageStateMap;
        for (auto it = replacementOutages.begin();
             it != replacementOutages.end();)
        {
            if (it->first.first == runtimeId)
            {
                it = replacementOutages.erase(it);
            }
            else
            {
                ++it;
            }
        }
        for (const auto& [system, outage] : restoredOutages)
        {
            replacementOutages[{runtimeId, system}] = outage;
        }

        auto replacementGraphs = graphStateMap;
        for (auto it = replacementGraphs.begin();
             it != replacementGraphs.end();)
        {
            if (it->first.first == runtimeId)
            {
                it = replacementGraphs.erase(it);
            }
            else
            {
                ++it;
            }
        }
        for (auto& [system, runtime] : restoredGraphs)
        {
            replacementGraphs[{runtimeId, system}] = std::move(runtime);
        }

        // std::map::swap with the standard allocator is noexcept.  These two
        // swaps are therefore the single commit point for the validated pair
        // of runtime maps.
        outageStateMap.swap(replacementOutages);
        graphStateMap.swap(replacementGraphs);
        return true;
    }
    catch (const std::exception& exception)
    {
        failureReason =
            "ZHANG_GRAPH_CHECKPOINT_RESTORE_PREPARE_FAILED:" +
            string(exception.what());
        return false;
    }
}

bool zhangGraphStochasticSupportValid(const KFState& state,const std::string& receiver,
 const SatSys& satellite,E_ObsCode code)
{
 auto o=acsConfig.zhangFullRank.sysOpts.find(satellite.sys);
 if(o==acsConfig.zhangFullRank.sysOpts.end() || !o->second.use_spanning_tree ||
    !zhangFullRankUsesObservable(code,o->second.baseline_observables))return false;
 auto it=graphStateMap.find({zhangGraphRuntimeId(state),satellite.sys});
 if(it==graphStateMap.end() || !it->second.initialized)return false;
 const auto& runtime=it->second;
 if(runtime.productBasisMode!=configuredProductBasisMode())return false;
 // A retained FLOAT coordinate remains available to numerical transport,
 // but it is not a current integer candidate in the single-S product mode.
 if(singleSBasisNative() && !zhangGraphModelsObservation(state,receiver,satellite,code))return false;
 const auto selected=selectedGraphContext(state,satellite.sys);
 const auto& basis=selected?selected->basis:
    (runtime.activeBasis.connected?runtime.activeBasis:runtime.basis);
 ZhangGraphEdge edge{receiver,satellite};
 auto history=runtime.edgeHistory.find(edge);
 // The controller removes discontinuous arcs from stateEdges before commit.
 // Requiring a propagated KF state at the caller forbids synthesizing columns.
 return runtime.stateEdges.count(edge)>0 && basis.edges.count(edge)>0 &&
    !basis.treeEdges.count(edge) && history!=runtime.edgeHistory.end() &&
    history->second.outageEpochs<=o->second.state_edge_grace_epochs;
}
