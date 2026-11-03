#pragma once
#include "pea/zhangReference.hpp"
namespace checkpoint_test_graph_wire {
constexpr char ZHANG_GRAPH_CHECKPOINT_MAGIC[] =
    "GINAN_ZHANG_GRAPH_RUNTIME";
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
}
