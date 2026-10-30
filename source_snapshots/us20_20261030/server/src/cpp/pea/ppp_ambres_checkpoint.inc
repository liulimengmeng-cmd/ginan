/* E29 Infra-0 checkpoint adapter for ppp_ambres.cpp private runtime.
 *
 * This file is included once, at the end of ppp_ambres.cpp.  It is not a
 * separately compiled translation unit: keeping the adapter here avoids
 * making the AR implementation's physical-lattice/window types public.
 */

#include "common/zhangMultiGnssFloatCheckpoint.hpp"
#include "common/zhangGecArCheckpointLimits.hpp"
#include <boost/json.hpp>

namespace
{
constexpr std::size_t ZHANG_AMBRES_MAX_CONTAINER = 200000;
constexpr std::size_t ZHANG_AMBRES_MAX_STRING = 1 << 20;
constexpr int ZHANG_E27_MAX_SATELLITE = 64;
constexpr int ZHANG_E27_MAX_FACTORS = 360;
constexpr double ZHANG_E27_LAG_SECONDS = 3600;

struct AmbresArcCheckpoint
{
    E_ObsCode code = E_ObsCode::NONE;
    string receiver;
    SatSys satellite;
    int version = 0;
    string coefficient;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & code;
        ar & receiver;
        ar & satellite;
        ar & version;
        ar & coefficient;
    }
};

struct AmbresHeldRowCheckpoint
{
    vector<AmbresArcCheckpoint> coefficients;
    string value;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & coefficients;
        ar & value;
    }
};

struct AmbresHeldLatticeCheckpoint
{
    string ownerRuntimeId;
    E_Sys system = E_Sys::NONE;
    vector<AmbresHeldRowCheckpoint> rows;
    int lastEventId = 0;
    bool consistent = true;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & system;
        ar & rows;
        ar & lastEventId;
        ar & consistent;
    }
};

struct AmbresHeldEvidenceCheckpoint
{
    string ownerRuntimeId;
    E_Sys system = E_Sys::NONE;
    AmbresHeldRowCheckpoint row;
    int confirmations = 0;
    long int lastEpoch = 0;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & system;
        ar & row;
        ar & confirmations;
        ar & lastEpoch;
    }
};

struct AmbresExactTermCheckpoint
{
    string column;
    string coefficient;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & column;
        ar & coefficient;
    }
};

struct AmbresSparseExactRowCheckpoint
{
    vector<AmbresExactTermCheckpoint> terms;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & terms;
    }
};

struct AmbresAdmissionCandidateCheckpoint
{
    int certificateKind = 0;
    string relationId;
    string satellite;
    string observable;
    AmbresSparseExactRowCheckpoint physicalCoefficients;
    string integerValue;
    bool exactIntegerEstimable = false;
    bool phaseSegmentCompatible = false;
    bool scalarReliabilityPassed = false;
    bool jointNisPassed = false;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & certificateKind;
        ar & relationId;
        ar & satellite;
        ar & observable;
        ar & physicalCoefficients;
        ar & integerValue;
        ar & exactIntegerEstimable;
        ar & phaseSegmentCompatible;
        ar & scalarReliabilityPassed;
        ar & jointNisPassed;
    }
};

struct AmbresTemporalConfirmationCheckpoint
{
    string relationKey;
    string integerValue;
    int consistentEpochs = 0;
    long int firstEpoch = 0;
    long int lastEpoch = 0;
    vector<string> supportSignatures;
    bool redundancyConfirmed = false;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & relationKey;
        ar & integerValue;
        ar & consistentEpochs;
        ar & firstEpoch;
        ar & lastEpoch;
        ar & supportSignatures;
        ar & redundancyConfirmed;
    }
};

struct AmbresAdmissionStateCheckpoint
{
    string ownerRuntimeId;
    E_Sys system = E_Sys::NONE;
    vector<AmbresSparseExactRowCheckpoint> certifiedRows;
    vector<string> certifiedValues;
    vector<string> certifiedSatellites;
    vector<string> certifiedRelationIds;
    vector<AmbresAdmissionCandidateCheckpoint> pendingCandidates;
    vector<AmbresTemporalConfirmationCheckpoint> temporalConfirmations;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & system;
        ar & certifiedRows;
        ar & certifiedValues;
        ar & certifiedSatellites;
        ar & certifiedRelationIds;
        ar & pendingCandidates;
        ar & temporalConfirmations;
    }
};

struct AmbresHeldUserWideLaneIntegerCheckpoint
{
    SatSys satellite;
    string value;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & satellite;
        ar & value;
    }
};

struct AmbresHeldUserWideLaneCheckpoint
{
    string ownerRuntimeId;
    string receiver;
    E_Sys system = E_Sys::NONE;
    bool ifAcceptance = false;
    SatSys reference;
    vector<AmbresHeldUserWideLaneIntegerCheckpoint> integers;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & receiver;
        ar & system;
        ar & ifAcceptance;
        ar & reference;
        ar & integers;
    }
};

struct AmbresGraphEdgeCheckpoint
{
    string receiver;
    SatSys satellite;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & receiver;
        ar & satellite;
    }
};

struct AmbresFunctionalCheckpoint
{
    SatSys satellite;
    SatSys referenceSatellite;
    vector<AmbresGraphEdgeCheckpoint> physicalEdges;
    vector<string> networkCoefficients;
    vector<int> physicalArcVersions;
    string affineOffsetCycles;
    int temporalBasisVersion = 0;
    bool valid = false;
    string failureReason;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & satellite;
        ar & referenceSatellite;
        ar & physicalEdges;
        ar & networkCoefficients;
        ar & physicalArcVersions;
        ar & affineOffsetCycles;
        ar & temporalBasisVersion;
        ar & valid;
        ar & failureReason;
    }
};

struct AmbresTransitionCheckpoint
{
    vector<AmbresGraphEdgeCheckpoint> physicalEdges;
    vector<int> physicalArcVersions;
    vector<string> coefficients;
    string affineOffsetCycles;
    bool valid = false;
    string failureReason;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & physicalEdges;
        ar & physicalArcVersions;
        ar & coefficients;
        ar & affineOffsetCycles;
        ar & valid;
        ar & failureReason;
    }
};

struct AmbresPendingTransitionCheckpoint
{
    ZhangCheckpointTime eventTime;
    E_Sys system = E_Sys::NONE;
    SatSys satellite;
    E_ObsCode observable = E_ObsCode::NONE;
    string eventId;
    AmbresFunctionalCheckpoint oldFunctional;
    AmbresFunctionalCheckpoint newFunctional;
    AmbresTransitionCheckpoint transition;
    string oldIdentity;
    string newIdentity;
    string oldSBasisFingerprint;
    string newSBasisFingerprint;
    string oldPhaseSegmentIdentity;
    string newPhaseSegmentIdentity;
    bool phaseSegmentChanged = false;
    string eventCause;
    int oldProductSegment = 0;
    int newProductSegment = 0;
    string oldSnapshotIdentity;
    string newSnapshotIdentity;
    string exactTransformChainId;
    int oldSnapshotReferenceCount = 0;
    int newSnapshotReferenceCount = 0;
    ZhangCheckpointTime expiryTime;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & eventTime;
        ar & system;
        ar & satellite;
        ar & observable;
        ar & eventId;
        ar & oldFunctional;
        ar & newFunctional;
        ar & transition;
        ar & oldIdentity;
        ar & newIdentity;
        ar & oldSBasisFingerprint;
        ar & newSBasisFingerprint;
        ar & oldPhaseSegmentIdentity;
        ar & newPhaseSegmentIdentity;
        ar & phaseSegmentChanged;
        ar & eventCause;
        ar & oldProductSegment;
        ar & newProductSegment;
        ar & oldSnapshotIdentity;
        ar & newSnapshotIdentity;
        ar & exactTransformChainId;
        ar & oldSnapshotReferenceCount;
        ar & newSnapshotReferenceCount;
        ar & expiryTime;
    }
};

struct AmbresActiveTransitionsCheckpoint
{
    string ownerRuntimeId;
    vector<AmbresPendingTransitionCheckpoint> transitions;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & transitions;
    }
};

struct AmbresE27ProductDatumCheckpoint
{
    SatSys satellite;
    int firstDiscontinuity = 0;
    int firstDatumVersion = 0;
    int secondDiscontinuity = 0;
    int secondDatumVersion = 0;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & satellite;
        ar & firstDiscontinuity;
        ar & firstDatumVersion;
        ar & secondDiscontinuity;
        ar & secondDatumVersion;
    }
};

struct AmbresE27LastValidCheckpoint
{
    SatSys satellite;
    ZhangCheckpointTime time;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & satellite;
        ar & time;
    }
};

struct AmbresE27RawFactorCheckpoint
{
    ZhangCheckpointTime time;
    vector<int> satellites;
    VectorXd values;
    MatrixXd covariance;
    vector<string> noiseKeys;
    VectorXd noiseVariances;
    MatrixXd satelliteNoiseDesign;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & time;
        ar & satellites;
        ar & values;
        ar & covariance;
        ar & noiseKeys;
        ar & noiseVariances;
        ar & satelliteNoiseDesign;
    }
};

struct AmbresE27RuntimeCheckpoint
{
    string ownerRuntimeId;
    string receiver;
    E_Sys system = E_Sys::NONE;
    string solution;
    vector<AmbresE27LastValidCheckpoint> lastValid;
    map<SatSys, int> arcVersions;
    vector<AmbresE27ProductDatumCheckpoint> productDatum;
    SatSys reference;
    map<int, int> accumulatorArcVersions;
    vector<AmbresE27RawFactorCheckpoint> factors;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & receiver;
        ar & system;
        ar & solution;
        ar & lastValid;
        ar & arcVersions;
        ar & productDatum;
        ar & reference;
        ar & accumulatorArcVersions;
        ar & factors;
    }
};

struct AmbresRelinkPriorCheckpoint
{
    string ownerRuntimeId;
    E_Sys system = E_Sys::NONE;
    SatSys anchor;
    SatSys satellite;
    ZhangCheckpointTime time;
    E_ObsCode firstCode = E_ObsCode::NONE;
    E_ObsCode secondCode = E_ObsCode::NONE;
    long int productDatumVersion = 0;
    double firstMean = 0;
    double firstVariance = 0;
    double wideLaneMean = 0;
    double wideLaneVariance = 0;
    double firstWideLaneCovariance = 0;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & system;
        ar & anchor;
        ar & satellite;
        ar & time;
        ar & firstCode;
        ar & secondCode;
        ar & productDatumVersion;
        ar & firstMean;
        ar & firstVariance;
        ar & wideLaneMean;
        ar & wideLaneVariance;
        ar & firstWideLaneCovariance;
    }
};

struct AmbresScalarIncrementCheckpoint
{
    ZhangCheckpointTime time;
    double information = 0;
    double natural = 0;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & time;
        ar & information;
        ar & natural;
    }
};

struct AmbresScalarAccumulatorCheckpoint
{
    string ownerRuntimeId;
    string key;
    vector<AmbresScalarIncrementCheckpoint> increments;
    long long integerHypothesis = 0;
    bool initialized = false;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & key;
        ar & increments;
        ar & integerHypothesis;
        ar & initialized;
    }
};

struct AmbresJointIncrementCheckpoint
{
    ZhangCheckpointTime time;
    std::array<double, 4> information{};
    std::array<double, 2> natural{};

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & time;
        ar & information;
        ar & natural;
    }
};

struct AmbresJointAccumulatorCheckpoint
{
    string ownerRuntimeId;
    string key;
    vector<AmbresJointIncrementCheckpoint> increments;
    long long wideLaneHypothesis = 0;
    long long firstHypothesis = 0;
    bool initialized = false;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & key;
        ar & increments;
        ar & wideLaneHypothesis;
        ar & firstHypothesis;
        ar & initialized;
    }
};

struct AmbresWhitenedObservationCheckpoint
{
    ZhangCheckpointTime time;
    double value = 0;
    double variance = 0;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & time;
        ar & value;
        ar & variance;
    }
};

struct AmbresWhitenedAccumulatorCheckpoint
{
    string ownerRuntimeId;
    string key;
    vector<AmbresWhitenedObservationCheckpoint> observations;
    vector<AmbresArcCheckpoint> physicalRow;
    string physicalSegmentIdentity;
    int lastGraphEventId = 0;
    int lastProductDatumVersion = 0;
    bool initialized = false;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & ownerRuntimeId;
        ar & key;
        ar & observations;
        ar & physicalRow;
        ar & physicalSegmentIdentity;
        ar & lastGraphEventId;
        ar & lastProductDatumVersion;
        ar & initialized;
    }
};

struct R51ProofAtomCheckpoint {
    string id,statement,provenance;double risk=1;vector<string> parents;
    template<class Archive>void serialize(Archive& ar,const unsigned int){ar&id&statement&provenance&risk&parents;}
};
struct R51PhysicalArchiveCheckpoint {
    string owner;E_Sys system=E_Sys::NONE;bool valid=false;
    vector<map<string,string>> rows;vector<string> values,roots;
    vector<R51ProofAtomCheckpoint> proofs;
    template<class Archive>void serialize(Archive& ar,const unsigned int){ar&owner&system&valid&rows&values&roots&proofs;}
};

struct AmbresCheckpointEnvelope
{
    std::uint32_t schemaVersion = ZHANG_AMBRES_CHECKPOINT_SCHEMA_VERSION;
    string runtimeId;
    bool filterError = false;
    bool transactionalConditioningFailed = false;
    string transactionalConditioningReason;
    bool arVerbose = false;
    double fixedAmbiguityVariance = 0;
    vector<R51PhysicalArchiveCheckpoint> r51Archives;
    vector<AmbresE27RuntimeCheckpoint> e27Runtimes;
    vector<AmbresHeldLatticeCheckpoint> heldLattices;
    vector<AmbresHeldEvidenceCheckpoint> heldEvidence;
    vector<AmbresAdmissionStateCheckpoint> admissionStates;
    vector<AmbresHeldUserWideLaneCheckpoint> heldUserWideLanes;
    vector<AmbresActiveTransitionsCheckpoint> activeTransitions;
    vector<AmbresRelinkPriorCheckpoint> relinkPriors;
    vector<AmbresScalarAccumulatorCheckpoint> scalarAccumulators;
    vector<AmbresJointAccumulatorCheckpoint> jointAccumulators;
    vector<AmbresWhitenedAccumulatorCheckpoint> whitenedAccumulators;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & schemaVersion;
        ar & runtimeId;
        ar & filterError;
        ar & transactionalConditioningFailed;
        ar & transactionalConditioningReason;
        ar & arVerbose;
        ar & fixedAmbiguityVariance;
        ar & r51Archives;
        ar & e27Runtimes;
        ar & heldLattices;
        ar & heldEvidence;
        ar & admissionStates;
        ar & heldUserWideLanes;
        ar & activeTransitions;
        ar & relinkPriors;
        ar & scalarAccumulators;
        ar & jointAccumulators;
        ar & whitenedAccumulators;
    }
};

/** Multi-GNSS FLOAT checkpoints may carry configuration, never private AR
 * runtime. In particular, an empty registered object is still a private
 * runtime entry and must not be silently dropped to pass this gate. */
bool multiGnssFloatAmbresCheckpointEmpty(
    const AmbresCheckpointEnvelope& snapshot, string& failureReason)
{
    if (!zhangMultiGnssFloatCheckpointActive(acsConfig)) return true;
    const auto reject = [&](const char* field)
    {
        failureReason = string("AMBRES_MULTI_GNSS_FLOAT_PRIVATE_STATE_NONEMPTY:") + field;
        return false;
    };
    // Existing writers skip empty held-user WL and empty nested registries.
    // Their registration itself is out of scope, on export and import alike.
#define ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(field) \
    if (!field.empty()) return reject("live_" #field)
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(r51PhysicalArchives);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangE27WideLaneRuntimes);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangPersistentHeldLattices);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangPersistentHeldEvidence);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangProductRelationAdmissionStateRegistry());
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangHeldUserWideLaneRegistry);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangActiveTemporalBesdTransitions);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangRelinkPriorMoments);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangRelinkShadowAccumulators);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangRelinkJointAccumulators);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangWhitenedWlAccumulators);
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangCurrentEpochEvidenceRegistry());
    ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY(zhangMultiGnssCurrentEvidenceAuthorities());
#undef ZHANG_FLOAT_AMBRES_REQUIRE_LIVE_EMPTY
    if (::filterError) return reject("live_filterError");
    if (zhangTransactionalConditioningFailed)
        return reject("live_transactionalConditioningFailed");
    if (!zhangTransactionalConditioningReason.empty())
        return reject("live_transactionalConditioningReason");
    if (snapshot.filterError) return reject("filterError");
    if (snapshot.transactionalConditioningFailed)
        return reject("transactionalConditioningFailed");
    if (!snapshot.transactionalConditioningReason.empty())
        return reject("transactionalConditioningReason");
#define ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(field) \
    if (!snapshot.field.empty()) return reject(#field)
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(r51Archives);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(e27Runtimes);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(heldLattices);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(heldEvidence);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(admissionStates);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(heldUserWideLanes);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(activeTransitions);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(relinkPriors);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(scalarAccumulators);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(jointAccumulators);
    ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY(whitenedAccumulators);
#undef ZHANG_FLOAT_AMBRES_REQUIRE_EMPTY
    // Preserve arVerbose/fixedAmbiguityVariance as configuration. The full
    // wire comparison fails closed if a future member adds private history.
    AmbresCheckpointEnvelope allowed;
    allowed.schemaVersion = snapshot.schemaVersion;
    allowed.runtimeId = snapshot.runtimeId;
    allowed.arVerbose = snapshot.arVerbose;
    allowed.fixedAmbiguityVariance = snapshot.fixedAmbiguityVariance;
    const auto wire = [](const AmbresCheckpointEnvelope& value)
    {
        std::ostringstream output(std::ios::binary | std::ios::out);
        boost::archive::binary_oarchive archive(output, boost::archive::no_header);
        archive << value;
        return output.str();
    };
    try
    {
        if (wire(snapshot) != wire(allowed)) return reject("unrecognized_private_history");
    }
    catch (const std::exception&)
    {
        failureReason = "AMBRES_MULTI_GNSS_FLOAT_PRIVATE_STATE_AUDIT_FAILED";
        return false;
    }
    return true;
}

std::mutex zhangAmbresCheckpointMutex;

bool ambresValidString(const string& value, bool allowEmpty = true)
{
    return value.size() <= ZHANG_AMBRES_MAX_STRING &&
        (allowEmpty || !value.empty());
}

bool ambresValidSystem(E_Sys system)
{
    return system == E_Sys::GPS ||
        (zhangMultiGnssArCheckpointActive(acsConfig) &&
         (system == E_Sys::GAL || system == E_Sys::BDS));
}

bool ambresValidGpsSatellite(const SatSys& satellite, bool allowNone = false)
{
    if (allowNone && satellite.prn == 0)
    {
        return satellite.sys == E_Sys::NONE || ambresValidSystem(satellite.sys);
    }
    return ambresValidSystem(satellite.sys) && satellite.prn > 0 &&
        satellite.prn <= ZHANG_E27_MAX_SATELLITE;
}

bool ambresValidObservable(E_ObsCode observable, E_Sys system = E_Sys::GPS)
{
    if (system == E_Sys::GPS)
        return observable == E_ObsCode::L1C || observable == E_ObsCode::L2W;
    if (!zhangMultiGnssArCheckpointActive(acsConfig)) return false;
    if (system == E_Sys::GAL)
        return observable == E_ObsCode::L1C || observable == E_ObsCode::L5Q;
    if (system == E_Sys::BDS)
        return observable == E_ObsCode::L1P || observable == E_ObsCode::L5P;
    return false;
}

bool ambresSatelliteUsesSystem(const SatSys& satellite, E_Sys system, bool allowNone = false)
{
    return ambresValidGpsSatellite(satellite, allowNone) &&
        ((allowNone && satellite.prn == 0 && satellite.sys == E_Sys::NONE) ||
         satellite.sys == system);
}

bool ambresRestoreFiniteTime(
    const ZhangCheckpointTime& snapshot,
    GTime& value)
{
    value = restoreZhangCheckpointTime(snapshot);
    return std::isfinite(value.bigTime);
}

string ambresExactText(const ZhangExactInteger& value)
{
    return value.str(0, std::ios_base::fmtflags(0));
}

bool ambresParseExact(
    const string& text,
    ZhangExactInteger& value,
    string& reason)
{
    if (text.empty() || text.size() > 4096)
    {
        reason = "AMBRES_CHECKPOINT_EXACT_INTEGER_LENGTH_INVALID";
        return false;
    }
    std::size_t first = text.front() == '-' ? 1 : 0;
    if (first == text.size() ||
        !std::all_of(text.begin() + first, text.end(), [](unsigned char ch)
        {
            return std::isdigit(ch) != 0;
        }))
    {
        reason = "AMBRES_CHECKPOINT_EXACT_INTEGER_TEXT_INVALID";
        return false;
    }
    try
    {
        value = ZhangExactInteger(text);
    }
    catch (const std::exception&)
    {
        reason = "AMBRES_CHECKPOINT_EXACT_INTEGER_PARSE_FAILED";
        return false;
    }
    return true;
}

AmbresArcCheckpoint ambresCheckpointArc(
    const ZhangPhysicalIntegerArc& arc,
    const ZhangExactInteger& coefficient)
{
    return {arc.code, arc.edge.receiver, arc.edge.satellite, arc.version,
        ambresExactText(coefficient)};
}

bool ambresRestoreArc(
    const AmbresArcCheckpoint& snapshot,
    ZhangPhysicalIntegerArc& arc,
    ZhangExactInteger& coefficient,
    string& reason)
{
    if (!ambresValidObservable(snapshot.code, snapshot.satellite.sys) ||
        !ambresValidString(snapshot.receiver, false) ||
        !ambresValidGpsSatellite(snapshot.satellite) || snapshot.version < 0 ||
        !ambresParseExact(snapshot.coefficient, coefficient, reason) ||
        coefficient == 0)
    {
        if (reason.empty())
        {
            reason = "AMBRES_CHECKPOINT_PHYSICAL_ARC_INVALID";
        }
        return false;
    }
    arc = {snapshot.code, {snapshot.receiver, snapshot.satellite},
        snapshot.version};
    return true;
}

AmbresHeldRowCheckpoint ambresCheckpointHeldRow(
    const ZhangPersistentHeldRow& row)
{
    AmbresHeldRowCheckpoint result;
    for (const auto& [arc, coefficient] : row.coefficients)
    {
        result.coefficients.push_back(ambresCheckpointArc(arc, coefficient));
    }
    result.value = ambresExactText(row.value);
    return result;
}

bool ambresRestoreHeldRow(
    const AmbresHeldRowCheckpoint& snapshot,
    ZhangPersistentHeldRow& row,
    string& reason)
{
    if (snapshot.coefficients.empty() ||
        snapshot.coefficients.size() > ZHANG_AMBRES_MAX_CONTAINER ||
        !ambresParseExact(snapshot.value, row.value, reason))
    {
        if (reason.empty())
        {
            reason = "AMBRES_CHECKPOINT_HELD_ROW_INVALID";
        }
        return false;
    }
    for (const auto& item : snapshot.coefficients)
    {
        ZhangPhysicalIntegerArc arc;
        ZhangExactInteger coefficient;
        if (!ambresRestoreArc(item, arc, coefficient, reason) ||
            !row.coefficients.emplace(arc, coefficient).second)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_HELD_ROW_DUPLICATE_ARC";
            }
            return false;
        }
    }
    return true;
}

bool ambresHeldRowUsesSystem(const ZhangPersistentHeldRow& row, E_Sys system)
{
    return std::all_of(row.coefficients.begin(), row.coefficients.end(),
        [&](const auto& item) { return item.first.edge.satellite.sys == system &&
            ambresValidObservable(item.first.code, system); });
}

bool ambresGecArchiveScope(const AmbresCheckpointEnvelope& snapshot, string& reason)
{
    if (!zhangMultiGnssArCheckpointActive(acsConfig)) return true;
    for (const auto& item : snapshot.r51Archives)
    {
        const auto codes = zhangCurrentNetworkBaselineCodes(item.system);
        const auto rejectScope = [&](const char* detail)
        {
            reason = string("AMBRES_GEC_R51_ARCHIVE_SCOPE_INVALID:") + detail +
                ":system=" + enum_to_string(item.system) + ":owner=" + item.owner +
                ":rows=" + std::to_string(item.rows.size()) +
                ":values=" + std::to_string(item.values.size());
            return false;
        };
        if (!ambresValidSystem(item.system)) return rejectScope("UNSUPPORTED_SYSTEM");
        if (!codes) return rejectScope("BASELINE_NOT_ACTIVE_OR_CONFIGURED");
        if (item.rows.size() != item.values.size()) return rejectScope("ROW_VALUE_COUNT_MISMATCH");
        if (item.rows.size() > ZHANG_AMBRES_MAX_CONTAINER) return rejectScope("ROW_COUNT_LIMIT");
        for (const auto& source : item.rows)
        {
            map<string,ZhangExactInteger> row;
            for (const auto& [id,text] : source)
            {
                ZhangExactInteger x;
                if (!ambresParseExact(text,x,reason) || x == 0)
                { reason = "AMBRES_GEC_R51_ARCHIVE_INTEGER_INVALID"; return false; }
                row.emplace(id,std::move(x));
            }
            if (!zhangMultiGnssPhysicalRowsUseSystem({row},item.system,*codes))
            { reason = "AMBRES_GEC_R51_ARCHIVE_SYSTEM_SIGNAL_INVALID"; return false; }
        }
    }
    return true;
}

AmbresSparseExactRowCheckpoint ambresCheckpointSparseExactRow(
    const map<string, ZhangExactInteger>& row)
{
    AmbresSparseExactRowCheckpoint result;
    for (const auto& [column, coefficient] : row)
    {
        result.terms.push_back({column, ambresExactText(coefficient)});
    }
    return result;
}

bool ambresRestoreSparseExactRow(
    const AmbresSparseExactRowCheckpoint& snapshot,
    map<string, ZhangExactInteger>& row,
    string& reason,
    bool allowEmpty = false)
{
    if ((!allowEmpty && snapshot.terms.empty())
        || snapshot.terms.size() > ZHANG_AMBRES_MAX_CONTAINER)
    {
        reason = "AMBRES_CHECKPOINT_SPARSE_EXACT_ROW_INVALID";
        return false;
    }
    for (const auto& term : snapshot.terms)
    {
        ZhangExactInteger coefficient;
        if (!ambresValidString(term.column, false)
            || !ambresParseExact(term.coefficient, coefficient, reason)
            || coefficient == 0
            || !row.emplace(term.column, coefficient).second)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_SPARSE_EXACT_TERM_INVALID";
            }
            return false;
        }
    }
    return true;
}

AmbresAdmissionCandidateCheckpoint ambresCheckpointAdmissionCandidate(
    const ZhangProductRelationAdmissionCandidate& candidate)
{
    AmbresAdmissionCandidateCheckpoint result;
    result.certificateKind = static_cast<int>(candidate.certificateKind);
    result.relationId = candidate.relationId;
    result.satellite = candidate.satellite;
    result.observable = candidate.observable;
    result.physicalCoefficients =
        ambresCheckpointSparseExactRow(candidate.physicalCoefficients);
    result.integerValue = ambresExactText(candidate.integerValue);
    result.exactIntegerEstimable = candidate.exactIntegerEstimable;
    result.phaseSegmentCompatible = candidate.phaseSegmentCompatible;
    result.scalarReliabilityPassed = candidate.scalarReliabilityPassed;
    result.jointNisPassed = candidate.jointNisPassed;
    return result;
}

bool ambresRestoreAdmissionCandidate(
    const AmbresAdmissionCandidateCheckpoint& snapshot,
    ZhangProductRelationAdmissionCandidate& candidate,
    string& reason, E_Sys system = E_Sys::GPS)
{
    if ((snapshot.certificateKind != static_cast<int>(
            TemporalCertificateKind::SELF_GAUGE_SHIFT)
            && snapshot.certificateKind != static_cast<int>(
                TemporalCertificateKind::INTER_SATELLITE_BRIDGE))
        || !ambresValidString(snapshot.relationId, false)
        || !ambresValidString(snapshot.satellite, false)
        || !ambresValidString(snapshot.observable, false)
        || !ambresParseExact(snapshot.integerValue,
            candidate.integerValue, reason)
        || !ambresRestoreSparseExactRow(snapshot.physicalCoefficients,
            candidate.physicalCoefficients, reason))
    {
        if (reason.empty())
        {
            reason = "AMBRES_CHECKPOINT_ADMISSION_CANDIDATE_INVALID";
        }
        return false;
    }
    candidate.certificateKind = static_cast<TemporalCertificateKind>(
        snapshot.certificateKind);
    candidate.relationId = snapshot.relationId;
    candidate.satellite = snapshot.satellite;
    candidate.observable = snapshot.observable;
    candidate.exactIntegerEstimable = snapshot.exactIntegerEstimable;
    candidate.phaseSegmentCompatible = snapshot.phaseSegmentCompatible;
    candidate.scalarReliabilityPassed = snapshot.scalarReliabilityPassed;
    candidate.jointNisPassed = snapshot.jointNisPassed;
    if (zhangMultiGnssArCheckpointActive(acsConfig))
    {
        const SatSys satellite(snapshot.satellite.c_str());
        const auto codes = zhangCurrentNetworkBaselineCodes(system);
        if (!ambresSatelliteUsesSystem(satellite, system) || !codes ||
            (snapshot.observable != enum_to_string((*codes)[0]) &&
             snapshot.observable != enum_to_string((*codes)[1])) ||
            !zhangMultiGnssPhysicalRowsUseSystem({candidate.physicalCoefficients}, system, *codes))
        {
            reason = "AMBRES_CHECKPOINT_ADMISSION_SYSTEM_SIGNAL_INVALID";
            return false;
        }
    }
    return true;
}

AmbresFunctionalCheckpoint ambresCheckpointFunctional(
    const ZhangProductIntegerFunctional& value)
{
    AmbresFunctionalCheckpoint result;
    result.satellite = value.satellite;
    result.referenceSatellite = value.referenceSatellite;
    for (const auto& edge : value.physicalEdges)
    {
        result.physicalEdges.push_back({edge.receiver, edge.satellite});
    }
    for (const auto& coefficient : value.networkCoefficients)
    {
        result.networkCoefficients.push_back(ambresExactText(coefficient));
    }
    result.physicalArcVersions = value.physicalArcVersions;
    result.affineOffsetCycles = ambresExactText(value.affineOffsetCycles);
    result.temporalBasisVersion = value.temporalBasisVersion;
    result.valid = value.valid;
    result.failureReason = value.failureReason;
    return result;
}

bool ambresRestoreFunctional(
    const AmbresFunctionalCheckpoint& snapshot,
    ZhangProductIntegerFunctional& value,
    string& reason)
{
    if (snapshot.physicalEdges.size() != snapshot.networkCoefficients.size() ||
        snapshot.physicalEdges.size() != snapshot.physicalArcVersions.size() ||
        snapshot.physicalEdges.size() > ZHANG_AMBRES_MAX_CONTAINER ||
        !ambresValidString(snapshot.failureReason))
    {
        reason = "AMBRES_CHECKPOINT_FUNCTIONAL_DIMENSION_INVALID";
        return false;
    }
    ZhangProductIntegerFunctional candidate;
    candidate.satellite = snapshot.satellite;
    candidate.referenceSatellite = snapshot.referenceSatellite;
    if (snapshot.valid &&
        (!ambresValidGpsSatellite(candidate.satellite) ||
         !ambresValidGpsSatellite(candidate.referenceSatellite) ||
         candidate.satellite.sys != candidate.referenceSatellite.sys))
    {
        reason = "AMBRES_CHECKPOINT_FUNCTIONAL_SATELLITE_INVALID";
        return false;
    }
    for (std::size_t index = 0; index < snapshot.physicalEdges.size(); index++)
    {
        const auto& edge = snapshot.physicalEdges[index];
        ZhangExactInteger coefficient;
        if (!ambresValidString(edge.receiver, false) ||
            !ambresValidGpsSatellite(edge.satellite) ||
            (snapshot.valid && edge.satellite.sys != candidate.satellite.sys) ||
            // Version zero is the established representation for a physical
            // product-basis edge that has no registered arc version yet.
            // It must round-trip exactly; only negative versions are invalid.
            snapshot.physicalArcVersions[index] < 0 ||
            !ambresParseExact(
                snapshot.networkCoefficients[index], coefficient, reason))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_FUNCTIONAL_EDGE_INVALID";
            }
            return false;
        }
        candidate.physicalEdges.push_back({edge.receiver, edge.satellite});
        candidate.networkCoefficients.push_back(coefficient);
    }
    candidate.physicalArcVersions = snapshot.physicalArcVersions;
    if (!ambresParseExact(
            snapshot.affineOffsetCycles, candidate.affineOffsetCycles, reason))
    {
        return false;
    }
    candidate.temporalBasisVersion = snapshot.temporalBasisVersion;
    candidate.valid = snapshot.valid;
    candidate.failureReason = snapshot.failureReason;
    value = std::move(candidate);
    return true;
}

AmbresTransitionCheckpoint ambresCheckpointTransition(
    const ZhangProductIntegerTransition& value)
{
    AmbresTransitionCheckpoint result;
    for (const auto& edge : value.physicalEdges)
    {
        result.physicalEdges.push_back({edge.receiver, edge.satellite});
    }
    result.physicalArcVersions = value.physicalArcVersions;
    for (const auto& coefficient : value.coefficients)
    {
        result.coefficients.push_back(ambresExactText(coefficient));
    }
    result.affineOffsetCycles = ambresExactText(value.affineOffsetCycles);
    result.valid = value.valid;
    result.failureReason = value.failureReason;
    return result;
}

bool ambresRestoreTransition(
    const AmbresTransitionCheckpoint& snapshot,
    ZhangProductIntegerTransition& value,
    string& reason)
{
    if (snapshot.physicalEdges.size() != snapshot.coefficients.size() ||
        snapshot.physicalEdges.size() != snapshot.physicalArcVersions.size() ||
        snapshot.physicalEdges.size() > ZHANG_AMBRES_MAX_CONTAINER ||
        !ambresValidString(snapshot.failureReason))
    {
        reason = "AMBRES_CHECKPOINT_TRANSITION_DIMENSION_INVALID";
        return false;
    }
    ZhangProductIntegerTransition candidate;
    E_Sys transitionSystem = E_Sys::NONE;
    for (std::size_t index = 0; index < snapshot.physicalEdges.size(); index++)
    {
        const auto& edge = snapshot.physicalEdges[index];
        ZhangExactInteger coefficient;
        if (!ambresValidString(edge.receiver, false) ||
            !ambresValidGpsSatellite(edge.satellite) ||
            (transitionSystem != E_Sys::NONE && edge.satellite.sys != transitionSystem) ||
            snapshot.physicalArcVersions[index] < 0 ||
            !ambresParseExact(snapshot.coefficients[index], coefficient, reason) ||
            coefficient == 0)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_TRANSITION_EDGE_INVALID";
            }
            return false;
        }
        candidate.physicalEdges.push_back({edge.receiver, edge.satellite});
        transitionSystem = edge.satellite.sys;
        candidate.coefficients.push_back(coefficient);
    }
    candidate.physicalArcVersions = snapshot.physicalArcVersions;
    if (!ambresParseExact(
            snapshot.affineOffsetCycles, candidate.affineOffsetCycles, reason))
    {
        return false;
    }
    candidate.valid = snapshot.valid;
    candidate.failureReason = snapshot.failureReason;
    value = std::move(candidate);
    return true;
}

AmbresPendingTransitionCheckpoint ambresCheckpointPending(
    const ZhangPendingProductTransition& value)
{
    AmbresPendingTransitionCheckpoint result;
    result.eventTime = captureZhangCheckpointTime(value.eventTime);
    result.system = value.system;
    result.satellite = value.satellite;
    result.observable = value.observable;
    result.eventId = value.eventId;
    result.oldFunctional = ambresCheckpointFunctional(value.oldFunctional);
    result.newFunctional = ambresCheckpointFunctional(value.newFunctional);
    result.transition = ambresCheckpointTransition(value.transition);
    result.oldIdentity = value.oldIdentity;
    result.newIdentity = value.newIdentity;
    result.oldSBasisFingerprint = value.oldSBasisFingerprint;
    result.newSBasisFingerprint = value.newSBasisFingerprint;
    result.oldPhaseSegmentIdentity = value.oldPhaseSegmentIdentity;
    result.newPhaseSegmentIdentity = value.newPhaseSegmentIdentity;
    result.phaseSegmentChanged = value.phaseSegmentChanged;
    result.eventCause = value.eventCause;
    result.oldProductSegment = value.oldProductSegment;
    result.newProductSegment = value.newProductSegment;
    result.oldSnapshotIdentity = value.oldSnapshotIdentity;
    result.newSnapshotIdentity = value.newSnapshotIdentity;
    result.exactTransformChainId = value.exactTransformChainId;
    result.oldSnapshotReferenceCount = value.oldSnapshotReferenceCount;
    result.newSnapshotReferenceCount = value.newSnapshotReferenceCount;
    result.expiryTime = captureZhangCheckpointTime(value.expiryTime);
    return result;
}

bool ambresRestorePending(
    const AmbresPendingTransitionCheckpoint& snapshot,
    ZhangPendingProductTransition& value,
    string& reason)
{
    ZhangPendingProductTransition candidate;
    if (!ambresValidSystem(snapshot.system) ||
        !ambresSatelliteUsesSystem(snapshot.satellite, snapshot.system) ||
        !ambresValidObservable(snapshot.observable, snapshot.system) ||
        !ambresValidString(snapshot.eventId, false) ||
        snapshot.oldSnapshotReferenceCount < 0 ||
        snapshot.newSnapshotReferenceCount < 0 ||
        !ambresRestoreFiniteTime(snapshot.eventTime, candidate.eventTime) ||
        !ambresRestoreFiniteTime(snapshot.expiryTime, candidate.expiryTime) ||
        candidate.expiryTime < candidate.eventTime)
    {
        reason = "AMBRES_CHECKPOINT_PENDING_TRANSITION_INVALID";
        return false;
    }
    candidate.system = snapshot.system;
    candidate.satellite = snapshot.satellite;
    candidate.observable = snapshot.observable;
    candidate.eventId = snapshot.eventId;
    if (!ambresRestoreFunctional(
            snapshot.oldFunctional, candidate.oldFunctional, reason) ||
        !ambresRestoreFunctional(
            snapshot.newFunctional, candidate.newFunctional, reason) ||
        !ambresRestoreTransition(
            snapshot.transition, candidate.transition, reason))
    {
        return false;
    }
    if (!candidate.oldFunctional.valid || !candidate.newFunctional.valid ||
        !candidate.transition.valid ||
        candidate.oldFunctional.satellite != candidate.satellite ||
        candidate.newFunctional.satellite != candidate.satellite ||
        std::any_of(candidate.transition.physicalEdges.begin(),candidate.transition.physicalEdges.end(),
            [&](const auto& edge) { return edge.satellite.sys != candidate.system; }))
    {
        reason = "AMBRES_CHECKPOINT_PENDING_FUNCTIONAL_REFERENCE_INVALID";
        return false;
    }
    const std::array<const string*, 11> texts = {
        &snapshot.oldIdentity, &snapshot.newIdentity,
        &snapshot.oldSBasisFingerprint, &snapshot.newSBasisFingerprint,
        &snapshot.oldPhaseSegmentIdentity, &snapshot.newPhaseSegmentIdentity,
        &snapshot.eventCause, &snapshot.oldSnapshotIdentity,
        &snapshot.newSnapshotIdentity, &snapshot.exactTransformChainId,
        &snapshot.transition.failureReason};
    if (!std::all_of(texts.begin(), texts.end(), [](const string* text)
        {
            return ambresValidString(*text);
        }))
    {
        reason = "AMBRES_CHECKPOINT_PENDING_TRANSITION_STRING_INVALID";
        return false;
    }
    candidate.oldIdentity = snapshot.oldIdentity;
    candidate.newIdentity = snapshot.newIdentity;
    candidate.oldSBasisFingerprint = snapshot.oldSBasisFingerprint;
    candidate.newSBasisFingerprint = snapshot.newSBasisFingerprint;
    candidate.oldPhaseSegmentIdentity = snapshot.oldPhaseSegmentIdentity;
    candidate.newPhaseSegmentIdentity = snapshot.newPhaseSegmentIdentity;
    candidate.phaseSegmentChanged = snapshot.phaseSegmentChanged;
    candidate.eventCause = snapshot.eventCause;
    candidate.oldProductSegment = snapshot.oldProductSegment;
    candidate.newProductSegment = snapshot.newProductSegment;
    candidate.oldSnapshotIdentity = snapshot.oldSnapshotIdentity;
    candidate.newSnapshotIdentity = snapshot.newSnapshotIdentity;
    candidate.exactTransformChainId = snapshot.exactTransformChainId;
    candidate.oldSnapshotReferenceCount = snapshot.oldSnapshotReferenceCount;
    candidate.newSnapshotReferenceCount = snapshot.newSnapshotReferenceCount;
    value = std::move(candidate);
    return true;
}

struct AmbresRestoredRuntime
{
    map<ZhangE27WideLaneKey, ZhangE27WideLaneRuntime> e27Runtimes;
    map<pair<string, E_Sys>, ZhangPersistentHeldLattice> heldLattices;
    map<
        pair<string, E_Sys>,
        map<ZhangPersistentHeldRow, ZhangPersistentHeldEvidence>> heldEvidence;
    map<pair<string, E_Sys>, ZhangProductRelationAdmissionState>
        admissionStates;
    map<ZhangHeldUserWideLaneKey, ZhangHeldUserWideLane>
        heldUserWideLanes;
    map<string, vector<ZhangPendingProductTransition>> activeTransitions;
    map<string, map<ZhangRelinkMomentKey, ZhangRelinkPriorMoment>> relinkPriors;
    map<string, map<string, ZhangRelinkShadowAccumulator>> scalarAccumulators;
    map<string, map<string, ZhangRelinkJointAccumulator>> jointAccumulators;
    map<string, map<string, ZhangWhitenedWlAccumulator>> whitenedAccumulators;
};

bool ambresRestoreEnvelope(
    const AmbresCheckpointEnvelope& snapshot,
    const string& rootRuntimeId,
    AmbresRestoredRuntime& restored,
    ZhangAmbresCheckpointResult& result);

void ambresPrepareCommittedMaps(
    const string& runtimeId,
    const AmbresRestoredRuntime& restored,
    AmbresRestoredRuntime& committed);

bool ambresHasForeignRuntime(const string& rootRuntimeId)
{
    for (const auto& [key, ignored] : zhangE27WideLaneRuntimes)
    {
        if (!zhangAmbresRuntimeBelongsTo(key.runtimeId, rootRuntimeId))
        {
            return true;
        }
    }
    const auto nestedForeign = [&](const auto& values, auto runtimeIdOf)
    {
        return std::any_of(values.begin(), values.end(), [&](const auto& item)
        {
            return !zhangAmbresRuntimeBelongsTo(
                runtimeIdOf(item.first), rootRuntimeId);
        });
    };
    return nestedForeign(zhangPersistentHeldLattices,
               [](const auto& key) -> const string& { return key.first; }) ||
        nestedForeign(zhangPersistentHeldEvidence,
               [](const auto& key) -> const string& { return key.first; }) ||
        nestedForeign(zhangProductRelationAdmissionStateRegistry(),
               [](const auto& key) -> const string& { return key.first; }) ||
        nestedForeign(zhangHeldUserWideLaneRegistry,
               [](const auto& key) -> const string&
               { return std::get<0>(key); }) ||
        nestedForeign(zhangActiveTemporalBesdTransitions,
               [](const auto& key) -> const string& { return key; }) ||
        nestedForeign(zhangRelinkPriorMoments,
               [](const auto& key) -> const string& { return key; }) ||
        nestedForeign(zhangRelinkShadowAccumulators,
               [](const auto& key) -> const string& { return key; }) ||
        nestedForeign(zhangRelinkJointAccumulators,
               [](const auto& key) -> const string& { return key; }) ||
        nestedForeign(zhangWhitenedWlAccumulators,
               [](const auto& key) -> const string& { return key; });
}

} // namespace

// This additional section does not alter either legacy archive. In particular,
// the old multi-GNSS FLOAT empty-state classifier still sees the original wire.
namespace
{
namespace gecjson = boost::json;

void gecArRequire(bool value, const char* reason)
{
    if (!value) throw std::runtime_error(reason);
}

const gecjson::object& gecArObject(const gecjson::value& value,
    std::initializer_list<const char*> fields)
{
    const auto& object = value.as_object();
    gecArRequire(object.size() == fields.size(), "GEC_AR_OBJECT_FIELDS_INVALID");
    for (const auto* field : fields)
        gecArRequire(object.contains(field), "GEC_AR_OBJECT_FIELD_MISSING");
    return object;
}

string gecArString(const gecjson::value& value, bool allowEmpty = true)
{
    const auto& text = value.as_string();
    string result(text.data(), text.size());
    gecArRequire(ambresValidString(result, allowEmpty), "GEC_AR_STRING_INVALID");
    return result;
}

const gecjson::array& gecArArray(const gecjson::value& value)
{
    const auto& array = value.as_array();
    gecArRequire(array.size() <= ZHANG_AMBRES_MAX_CONTAINER, "GEC_AR_CONTAINER_LIMIT");
    return array;
}

E_Sys gecArSystem(const gecjson::value& value, bool allowNone = false)
{
    const auto text = gecArString(value, false);
    if (allowNone && text == enum_to_string(E_Sys::NONE)) return E_Sys::NONE;
    for (auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
        if (text == enum_to_string(system)) return system;
    throw std::runtime_error("GEC_AR_SYSTEM_INVALID");
}

gecjson::value gecArPhysical(const ZhangCurrentEpochPhysicalRow& row)
{
    gecjson::object result;
    for (const auto& [id, coefficient] : row) result[id] = ambresExactText(coefficient);
    return result;
}

ZhangCurrentEpochPhysicalRow gecArPhysical(const gecjson::value& value, E_Sys system,
    bool allowEmpty = false)
{
    const auto& object = value.as_object();
    gecArRequire((allowEmpty || !object.empty()) && object.size() <= ZHANG_AMBRES_MAX_CONTAINER,
        "GEC_AR_PHYSICAL_ROW_SIZE_INVALID");
    ZhangCurrentEpochPhysicalRow row;
    for (const auto& entry : object)
    {
        ZhangExactInteger coefficient; string reason;
        const string id(entry.key().data(), entry.key().size());
        gecArRequire(ambresValidString(id, false) &&
            ambresParseExact(gecArString(entry.value(), false), coefficient, reason) && coefficient != 0,
            "GEC_AR_PHYSICAL_TERM_INVALID");
        row.emplace(id, std::move(coefficient));
    }
    const auto codes = zhangCurrentNetworkBaselineCodes(system);
    gecArRequire(codes && zhangMultiGnssPhysicalRowsUseSystem({row}, system, *codes),
        "GEC_AR_PHYSICAL_SYSTEM_SIGNAL_INVALID");
    return row;
}

gecjson::value gecArNumeric(const ZhangCoordinateRow& row)
{
    gecjson::array result;
    for (const auto& [key, coefficient] : row)
    {
        // The receiver address is replaced by its stable receiver ID. Identity
        // metadata and the full finite long-double timestamp remain lossless.
        const auto time = captureZhangCheckpointTime(key.estimatedTime);
        string timeHex; static constexpr char digits[] = "0123456789abcdef";
        for (auto byte : time.bigTimeBytes)
        { timeHex += digits[byte >> 4]; timeHex += digits[byte & 15]; }
        result.emplace_back(gecjson::object{{"type", static_cast<int>(key.type)},
            {"satellite", key.Sat.id()}, {"receiver", key.str}, {"number", key.num},
            {"comment", key.comment}, {"estimated_time", timeHex},
            {"receiver_pointer", key.rec_ptr ? key.rec_ptr->id : string{}}, {"coefficient", coefficient}});
    }
    return result;
}

ZhangCoordinateRow gecArNumeric(const gecjson::value& value, E_Sys system,
    const KFState& owner, int dimension)
{
    const auto& array = gecArArray(value);
    gecArRequire(!array.empty() && array.size() <= static_cast<std::size_t>(dimension),
        "GEC_AR_NUMERIC_ROW_SIZE_INVALID");
    ZhangCoordinateRow result;
    for (const auto& term : array)
    {
        const auto& item = gecArObject(term, {"type","satellite","receiver","number","comment",
            "estimated_time","receiver_pointer","coefficient"});
        KFKey key;
        const auto type = item.at("type").as_int64();
        gecArRequire(type == int(KF::AMBIGUITY) || type == int(KF::PHASE_BIAS) ||
            type == int(KF::PRODUCT_HANDOFF), "GEC_AR_NUMERIC_KEY_TYPE_INVALID");
        key.type = static_cast<KF>(type);
        key.Sat = SatSys(gecArString(item.at("satellite"), false).c_str());
        key.str = gecArString(item.at("receiver"));
        const auto number = item.at("number").as_int64();
        gecArRequire(number >= 0 && number <= std::numeric_limits<int>::max(),
            "GEC_AR_NUMERIC_KEY_NUMBER_INVALID");
        key.num = static_cast<int>(number);
        gecArRequire(key.Sat.sys == system && (key.Sat.prn == 0 || ambresValidGpsSatellite(key.Sat)) &&
            ambresValidObservable(static_cast<E_ObsCode>(key.num), system),
            "GEC_AR_NUMERIC_SYSTEM_SIGNAL_INVALID");
        key.comment = gecArString(item.at("comment"));
        const auto timeHex = gecArString(item.at("estimated_time"), false);
        ZhangCheckpointTime time;
        gecArRequire(timeHex.size() == time.bigTimeBytes.size()*2, "GEC_AR_KEY_TIME_INVALID");
        auto nibble = [](char c) -> unsigned {
            if (c >= '0' && c <= '9') return c-'0';
            if (c >= 'a' && c <= 'f') return c-'a'+10;
            throw std::runtime_error("GEC_AR_KEY_TIME_INVALID");
        };
        for (std::size_t i=0; i<time.bigTimeBytes.size(); ++i)
            time.bigTimeBytes[i] = static_cast<unsigned char>((nibble(timeHex[2*i])<<4)|nibble(timeHex[2*i+1]));
        key.estimatedTime = restoreZhangCheckpointTime(time);
        gecArRequire(std::isfinite(key.estimatedTime.bigTime), "GEC_AR_KEY_TIME_INVALID");
        const auto found = owner.kfIndexMap.find(key);
        gecArRequire(found != owner.kfIndexMap.end() && found->second >= 0 && found->second < dimension,
            "GEC_AR_NUMERIC_KEY_NOT_IN_FLOAT_CORE");
        const auto receiver = gecArString(item.at("receiver_pointer"));
        if (!receiver.empty())
        {
            gecArRequire(found->first.rec_ptr && found->first.rec_ptr->id == receiver,
                "GEC_AR_RECEIVER_POINTER_BINDING_INVALID");
            key.rec_ptr = found->first.rec_ptr;
        }
        const auto& numberValue = item.at("coefficient");
        const double coefficient = numberValue.is_double() ? numberValue.as_double() :
            numberValue.is_int64() ? static_cast<double>(numberValue.as_int64()) :
            numberValue.is_uint64() ? static_cast<double>(numberValue.as_uint64()) :
            std::numeric_limits<double>::quiet_NaN();
        gecArRequire(std::isfinite(coefficient) && coefficient != 0 && result.emplace(key, coefficient).second,
            "GEC_AR_NUMERIC_TERM_INVALID");
    }
    return result;
}

gecjson::value gecArProofGroups(const vector<ZhangDecisionProofs>& groups,
    map<string, ZhangDecisionProofPtr>& atoms)
{
    gecjson::array result;
    for (const auto& roots : groups)
    {
        const auto closure = zhangDecisionRiskClosure(roots);
        gecArRequire(closure.valid, "GEC_AR_PROOF_CLOSURE_INVALID");
        gecjson::array group;
        for (const auto& root : roots) group.emplace_back(root->id);
        for (const auto& [id, atom] : closure.atoms)
        {
            const auto prior = atoms.find(id);
            if (prior != atoms.end())
            {
                vector<string> oldParents, newParents;
                for (const auto& parent : prior->second->parents) oldParents.push_back(parent->id);
                for (const auto& parent : atom->parents) newParents.push_back(parent->id);
                gecArRequire(prior->second->originalStatement == atom->originalStatement &&
                    prior->second->observationProvenance == atom->observationProvenance &&
                    prior->second->conditionalFailureBound == atom->conditionalFailureBound && oldParents == newParents &&
                    prior->second.get() == atom.get(),
                    "GEC_AR_PROOF_ID_CONFLICT");
            }
            atoms[id] = atom;
        }
        result.emplace_back(std::move(group));
    }
    return result;
}

gecjson::object gecArEvidence(const ZhangCurrentEpochEvidenceKey& key,
    const ZhangCurrentEpochEvidence& value,
    map<string,ZhangDecisionProofPtr>& bundleProofAtoms)
{
    gecjson::object result{{"owner",key.first},{"key_system",enum_to_string(key.second)},
        {"runtime",value.runtimeId},{"system",enum_to_string(value.system)},
        {"epoch",value.epochGpstSeconds},{"posterior",value.posteriorId},{"chart_id",value.chartId},
        {"float_root",value.floatStateRoot},{"psd_validated",value.floatPosteriorPsdValidated},
        {"psd_root",value.floatPosteriorPsdRootId},{"state_order",value.floatPosteriorPsdStateOrderId},
        {"generation",value.constraintGeneration},{"valid",value.valid}};
    gecjson::array chords, expansions, unknown, active, signals;
    for (const auto& [column,id] : value.chart.chordIdentities)
        chords.emplace_back(gecjson::object{{"column",column},{"identity",id}});
    for (const auto& [column,row] : value.chart.expansions)
        expansions.emplace_back(gecjson::object{{"column",column},{"row",gecArPhysical(row)}});
    for (const auto& [column,reason] : value.chart.unexpandedColumns)
        unknown.emplace_back(gecjson::object{{"column",column},{"reason",reason}});
    for (const auto& id : value.chart.activePhysicalIdentities) active.emplace_back(id);
    for (const auto& id : value.chart.cataloguedSignals) signals.emplace_back(id);
    result["chart"] = gecjson::object{{"columns",value.chart.columns},{"chords",std::move(chords)},
        {"expansions",std::move(expansions)},{"unknown",std::move(unknown)},
        {"active",std::move(active)},{"signals",std::move(signals)},{"failure",value.chart.failureReason}};
    gecjson::array rows, values, eventRows, eventValues, witnesses, numeric, bindings;
    for (const auto& row : value.rows) rows.emplace_back(gecArPhysical(row));
    for (const auto& x : value.values) values.emplace_back(ambresExactText(x));
    for (const auto& row : value.eventOverlapRows) eventRows.emplace_back(gecArPhysical(row));
    for (const auto& x : value.eventOverlapValues) eventValues.emplace_back(ambresExactText(x));
    for (const auto& witness : value.eventOverlapWitnesses) witnesses.emplace_back(witness);
    for (const auto& row : value.numericConstraints)
        numeric.emplace_back(gecjson::object{{"row",gecArNumeric(row.rowCycles)},{"value",ambresExactText(row.value)}});
    for (const auto& row : value.physicalNumericBindings)
        bindings.emplace_back(gecjson::object{{"physical",gecArPhysical(row.physical)},
            {"numeric",gecArNumeric(row.rowCycles)},{"value",ambresExactText(row.value)},
            {"epoch",row.epochGpstSeconds},{"system",enum_to_string(row.system)},
            {"runtime",row.runtimeId},{"posterior",row.posteriorId},{"root",row.rootId},
            {"chart",row.chartId},{"order",row.stateOrderId}});
    result["rows"] = std::move(rows); result["values"] = std::move(values);
    result["event_rows"] = std::move(eventRows); result["event_values"] = std::move(eventValues);
    result["witnesses"] = std::move(witnesses); result["numeric"] = std::move(numeric);
    result["bindings"] = std::move(bindings);
    map<string,ZhangDecisionProofPtr> proofAtoms;
    result["proof_groups"] = gecArProofGroups(value.rowProofs, proofAtoms);
    result["event_proof_groups"] = gecArProofGroups(value.eventOverlapProofs, proofAtoms);
    gecjson::array proofs;
    for (const auto& [id,atom] : proofAtoms)
    {
        // Proof identity is bundle-wide, not constellation-local. A genuine
        // shared ancestor remains one immutable decision; distinct objects
        // claiming its ID are still rejected under the existing risk contract.
        const auto [prior,inserted] = bundleProofAtoms.emplace(id,atom);
        gecArRequire(inserted || prior->second.get() == atom.get(),
            "GEC_AR_CROSS_PARTITION_PROOF_ID_CONFLICT");
        gecjson::array parents;
        for (const auto& parent : atom->parents) parents.emplace_back(parent->id);
        proofs.emplace_back(gecjson::object{{"id",id},{"statement",atom->originalStatement},
            {"provenance",atom->observationProvenance},{"risk",atom->conditionalFailureBound},
            {"parents",std::move(parents)}});
    }
    result["proof_atoms"] = std::move(proofs);
    return result;
}

ZhangCurrentEpochEvidence gecArRestoreEvidence(const gecjson::value& item,
    const string& runtime, const KFState& owner, int dimension,
    const string& expectedRoot, const string& expectedOrder, ZhangCurrentEpochEvidenceKey& key,
    map<string,std::shared_ptr<ZhangIntegerDecisionProof>>& bundleProofAtoms)
{
    const auto& row = gecArObject(item, {"owner","key_system","runtime","system","epoch","posterior",
        "chart_id","float_root","psd_validated","psd_root","state_order","generation","valid",
        "chart","rows","values","event_rows","event_values","witnesses","numeric","bindings",
        "proof_groups","event_proof_groups","proof_atoms"});
    key = {gecArString(row.at("owner"),false),gecArSystem(row.at("key_system"))};
    gecArRequire(zhangAmbresRuntimeBelongsTo(key.first,runtime), "GEC_AR_FOREIGN_EVIDENCE_OWNER");
    ZhangCurrentEpochEvidence value;
    value.valid = row.at("valid").as_bool();
    value.runtimeId = gecArString(row.at("runtime")); value.system = gecArSystem(row.at("system"), !value.valid);
    value.epochGpstSeconds = row.at("epoch").as_int64();
    value.posteriorId = gecArString(row.at("posterior")); value.chartId = gecArString(row.at("chart_id"));
    value.floatStateRoot = gecArString(row.at("float_root"));
    value.floatPosteriorPsdValidated = row.at("psd_validated").as_bool();
    value.floatPosteriorPsdRootId = gecArString(row.at("psd_root"));
    value.floatPosteriorPsdStateOrderId = gecArString(row.at("state_order"));
    gecArRequire(row.at("generation").is_uint64() ||
        (row.at("generation").is_int64() && row.at("generation").as_int64() >= 0), "GEC_AR_GENERATION_INVALID");
    value.constraintGeneration = row.at("generation").is_uint64() ? row.at("generation").as_uint64() :
        static_cast<std::uint64_t>(row.at("generation").as_int64());
    gecArRequire(!value.valid || (value.runtimeId == key.first && value.system == key.second &&
        value.epochGpstSeconds == static_cast<long int>(owner.time.bigTime)), "GEC_AR_EVIDENCE_SCOPE_INVALID");
    if (value.valid)
        gecArRequire(!value.posteriorId.empty() && !value.chartId.empty() &&
            value.floatStateRoot == expectedRoot &&
            value.floatPosteriorPsdStateOrderId == expectedOrder &&
            (!value.floatPosteriorPsdValidated || value.floatPosteriorPsdRootId == value.floatStateRoot),
            "GEC_AR_EVIDENCE_FLOAT_ROOT_INVALID");
    const auto& chart = gecArObject(row.at("chart"), {"columns","chords","expansions","unknown","active","signals","failure"});
    const auto columns = chart.at("columns").as_int64();
    gecArRequire(columns >= 0 && columns <= dimension, "GEC_AR_CHART_DIMENSION_INVALID");
    value.chart.columns = static_cast<int>(columns); value.chart.failureReason = gecArString(chart.at("failure"));
    auto column = [&](const gecjson::object& object) {
        const auto c = object.at("column").as_int64();
        gecArRequire(c >= 0 && c < columns, "GEC_AR_CHART_COLUMN_INVALID"); return static_cast<int>(c);
    };
    for (const auto& term : gecArArray(chart.at("chords")))
    {
        const auto& object = gecArObject(term,{"column","identity"});
        const auto id = gecArString(object.at("identity"),false);
        gecArPhysical(gecjson::object{{id,"1"}},value.system);
        gecArRequire(value.chart.chordIdentities.emplace(column(object),id).second,"GEC_AR_DUPLICATE_CHORD");
    }
    for (const auto& term : gecArArray(chart.at("expansions")))
    {
        const auto& object = gecArObject(term,{"column","row"});
        gecArRequire(value.chart.expansions.emplace(column(object),gecArPhysical(object.at("row"),value.system)).second,
            "GEC_AR_DUPLICATE_EXPANSION");
    }
    for (const auto& term : gecArArray(chart.at("unknown")))
    {
        const auto& object = gecArObject(term,{"column","reason"}); const auto c = column(object);
        gecArRequire(!value.chart.expansions.contains(c) && !value.chart.chordIdentities.contains(c) &&
            value.chart.unexpandedColumns.emplace(c,gecArString(object.at("reason"),false)).second,
            "GEC_AR_UNKNOWN_COLUMN_CONFLICT");
    }
    for (const auto& term : gecArArray(chart.at("active")))
    {
        const auto id = gecArString(term,false); gecArPhysical(gecjson::object{{id,"1"}},value.system);
        gecArRequire(value.chart.activePhysicalIdentities.insert(id).second,"GEC_AR_DUPLICATE_ACTIVE_ARC");
    }
    for (const auto& term : gecArArray(chart.at("signals")))
    {
        const auto signal = gecArString(term,false); const auto codes = zhangCurrentNetworkBaselineCodes(value.system);
        gecArRequire(codes && (signal == enum_to_string((*codes)[0]) || signal == enum_to_string((*codes)[1])) &&
            value.chart.cataloguedSignals.insert(signal).second,"GEC_AR_CHART_SIGNAL_INVALID");
    }
    gecArRequire(value.chart.chordIdentities.size() == value.chart.expansions.size(), "GEC_AR_CHART_INCOMPLETE");
    for (const auto& [c,id] : value.chart.chordIdentities)
        gecArRequire(value.chart.expansions.contains(c) && value.chart.expansions.at(c).contains(id) &&
            value.chart.expansions.at(c).at(id) == 1, "GEC_AR_CHART_CHORD_COEFFICIENT_INVALID");
    auto restoreRows = [&](const char* name, auto& destination) {
        for (const auto& term : gecArArray(row.at(name))) destination.push_back(gecArPhysical(term,value.system));
    };
    auto restoreValues = [&](const char* name, auto& destination) {
        for (const auto& term : gecArArray(row.at(name)))
        { ZhangExactInteger x; string reason; gecArRequire(ambresParseExact(gecArString(term,false),x,reason),"GEC_AR_INTEGER_INVALID"); destination.push_back(std::move(x)); }
    };
    restoreRows("rows",value.rows); restoreRows("event_rows",value.eventOverlapRows);
    restoreValues("values",value.values); restoreValues("event_values",value.eventOverlapValues);
    for (const auto& physical : value.rows)
    {
        ZhangExactVector projected;
        gecArRequire(value.chart.project(physical,projected),"GEC_AR_CURRENT_ROW_CHART_INVALID");
    }
    for (const auto& term : gecArArray(row.at("witnesses"))) value.eventOverlapWitnesses.push_back(gecArString(term,false));
    map<string,std::shared_ptr<ZhangIntegerDecisionProof>> atoms;
    map<string,vector<string>> parents; map<string,set<string>> children; map<string,std::size_t> degree;
    for (const auto& term : gecArArray(row.at("proof_atoms")))
    {
        const auto& object = gecArObject(term,{"id","statement","provenance","risk","parents"});
        auto atom = std::make_shared<ZhangIntegerDecisionProof>(); atom->id = gecArString(object.at("id"),false);
        atom->originalStatement = gecArString(object.at("statement")); atom->observationProvenance = gecArString(object.at("provenance"));
        atom->conditionalFailureBound = object.at("risk").as_double();
        gecArRequire(std::isfinite(atom->conditionalFailureBound) && atom->conditionalFailureBound >= 0 &&
            atom->conditionalFailureBound <= 1 && atoms.emplace(atom->id,atom).second,"GEC_AR_PROOF_ATOM_INVALID");
        for (const auto& p : gecArArray(object.at("parents"))) parents[atom->id].push_back(gecArString(p,false));
        degree[atom->id] = 0;
    }
    for (const auto& [id,ps] : parents) for (const auto& p : ps)
    {
        gecArRequire(atoms.contains(p),"GEC_AR_PROOF_PARENT_MISSING");
        if (children[p].insert(id).second) ++degree[id];
    }
    vector<string> ready;
    for (const auto& [id,d] : degree) if (d == 0) ready.push_back(id);
    std::size_t visited = 0;
    for (std::size_t i=0; i<ready.size(); ++i)
    { ++visited; for (const auto& child : children[ready[i]]) if (--degree[child] == 0) ready.push_back(child); }
    gecArRequire(visited == atoms.size(),"GEC_AR_PROOF_CYCLE");
    // Parents precede their children in ready. Intern across the COMPLETE
    // section before exposing any row roots, retaining shared-ancestor pointer
    // identity required by zhangDecisionRiskClosure after a G/E/C joint union.
    for (const auto& id : ready)
    {
        auto& atom = atoms.at(id);
        for (const auto& p : parents[id]) atom->parents.push_back(atoms.at(p));
        const auto [prior,inserted] = bundleProofAtoms.emplace(id,atom);
        if (!inserted)
        {
            const auto& previous = prior->second;
            gecArRequire(previous->originalStatement == atom->originalStatement &&
                previous->observationProvenance == atom->observationProvenance &&
                previous->conditionalFailureBound == atom->conditionalFailureBound &&
                previous->parents == atom->parents,
                "GEC_AR_CROSS_PARTITION_PROOF_CONTENT_CONFLICT");
            atom = previous;
        }
    }
    set<string> usedAtoms;
    auto restoreGroups = [&](const char* name, auto& destination) {
        for (const auto& term : gecArArray(row.at(name)))
        {
            ZhangDecisionProofs group;
            for (const auto& idValue : gecArArray(term))
            { const auto id = gecArString(idValue,false); gecArRequire(atoms.contains(id),"GEC_AR_PROOF_ROOT_MISSING"); group.push_back(atoms.at(id)); }
            gecArRequire(!group.empty(),"GEC_AR_PHYSICAL_ROW_PROOF_MISSING");
            const auto closure = zhangDecisionRiskClosure(group); gecArRequire(closure.valid,"GEC_AR_PROOF_CLOSURE_INVALID");
            const auto currentSource = value.runtimeId + "|PPP_COMMITTED|E" + std::to_string(value.epochGpstSeconds);
            for (const auto& [id,atom] : closure.atoms)
            {
                gecArRequire(atom->observationProvenance == currentSource ||
                    atom->observationProvenance.rfind(currentSource+"|",0) == 0,"GEC_AR_PROOF_PROVENANCE_INVALID");
                usedAtoms.insert(id);
            }
            destination.push_back(std::move(group));
        }
    };
    restoreGroups("proof_groups",value.rowProofs); restoreGroups("event_proof_groups",value.eventOverlapProofs);
    gecArRequire(usedAtoms.size() == atoms.size() && value.rows.size() == value.values.size() &&
        value.rows.size() == value.rowProofs.size() && value.eventOverlapRows.size() == value.eventOverlapValues.size() &&
        value.eventOverlapRows.size() == value.eventOverlapProofs.size() &&
        value.eventOverlapRows.size() == value.eventOverlapWitnesses.size(), "GEC_AR_EVIDENCE_ROW_COUNTS_INVALID");
    for (const auto& term : gecArArray(row.at("numeric")))
    {
        const auto& object = gecArObject(term,{"row","value"}); ZhangExactInteger x; string reason;
        gecArRequire(ambresParseExact(gecArString(object.at("value"),false),x,reason),"GEC_AR_INTEGER_INVALID");
        value.numericConstraints.push_back({gecArNumeric(object.at("row"),value.system,owner,dimension),std::move(x)});
    }
    for (const auto& term : gecArArray(row.at("bindings")))
    {
        const auto& object = gecArObject(term,{"physical","numeric","value","epoch","system","runtime","posterior","root","chart","order"});
        ZhangCurrentEpochPhysicalNumericBinding binding; string reason;
        binding.physical = gecArPhysical(object.at("physical"),value.system);
        binding.rowCycles = gecArNumeric(object.at("numeric"),value.system,owner,dimension);
        gecArRequire(ambresParseExact(gecArString(object.at("value"),false),binding.value,reason),"GEC_AR_INTEGER_INVALID");
        binding.epochGpstSeconds = object.at("epoch").as_int64(); binding.system = gecArSystem(object.at("system"));
        binding.runtimeId = gecArString(object.at("runtime")); binding.posteriorId = gecArString(object.at("posterior"));
        binding.rootId = gecArString(object.at("root")); binding.chartId = gecArString(object.at("chart"));
        binding.stateOrderId = gecArString(object.at("order")); value.physicalNumericBindings.push_back(std::move(binding));
    }
    if (!value.rows.empty() || !value.eventOverlapRows.empty() || !value.physicalNumericBindings.empty())
    { string reason; gecArRequire(value.validatePhysicalNumericBindings(reason),"GEC_AR_PHYSICAL_NUMERIC_BINDING_INVALID"); }
    // Compiled HNF/domain objects contain no independent evidence. Rebuild on
    // first demand from the exact rows above; never serialize addresses/locks.
    value.domainCache = std::make_shared<ZhangCurrentEpochDomainCache>();
    return value;
}

using GecEvidenceMap = std::map<ZhangCurrentEpochEvidenceKey,ZhangCurrentEpochEvidence>;
using GecAuthorityMap = std::map<string,ZhangMultiGnssCurrentEvidenceAuthority>;

gecjson::object gecArSnapshot(const KFState& owner, const string& runtime)
{
    gecjson::array evidence, authorities;
    map<string,ZhangDecisionProofPtr> bundleProofAtoms;
    for (const auto& [key,value] : zhangCurrentEpochEvidenceRegistry())
    {
        gecArRequire(zhangAmbresRuntimeBelongsTo(key.first,runtime),"GEC_AR_FOREIGN_EVIDENCE_OWNER");
        evidence.emplace_back(gecArEvidence(key,value,bundleProofAtoms));
    }
    for (const auto& [key,value] : zhangMultiGnssCurrentEvidenceAuthorities())
    {
        gecArRequire(zhangAmbresRuntimeBelongsTo(key,runtime),"GEC_AR_FOREIGN_AUTHORITY_OWNER");
        gecjson::object systems;
        for (const auto& [system,ids] : value.systemIdentities)
            systems[enum_to_string(system)] = gecjson::array{ids.first,ids.second};
        authorities.emplace_back(gecjson::object{{"owner",key},{"epoch",value.epochGpstSeconds},
            {"float_root",value.floatRoot},{"state_order",value.floatStateOrderId},{"systems",std::move(systems)}});
    }
    return {{"schema",static_cast<int>(ZHANG_GEC_AR_RUNTIME_CHECKPOINT_SCHEMA_VERSION)},{"runtime",runtime},
        {"epoch",static_cast<long int>(owner.time.bigTime)},{"dimension",owner.x.size()},
        {"float_root",zhangR51NumericRoot(owner.x,owner.P)},
        {"state_order",zhangCurrentEvidenceStateOrderId(owner.kfIndexMap)},
        {"evidence",std::move(evidence)},{"authorities",std::move(authorities)}};
}

void gecArRestoreSnapshot(const gecjson::value& document, const KFState& owner,
    const string& runtime, int dimension, GecEvidenceMap& evidence, GecAuthorityMap& authorities)
{
    const auto& object = gecArObject(document,{"schema","runtime","epoch","dimension","float_root","state_order","evidence","authorities"});
    const auto expectedRoot = zhangR51NumericRoot(owner.x,owner.P);
    const auto expectedOrder = zhangCurrentEvidenceStateOrderId(owner.kfIndexMap);
    gecArRequire(object.at("schema").as_int64() == ZHANG_GEC_AR_RUNTIME_CHECKPOINT_SCHEMA_VERSION &&
        gecArString(object.at("runtime"),false) == runtime && dimension >= 0 && owner.x.size() == dimension &&
        owner.P.rows() == dimension && owner.P.cols() == dimension && owner.x.allFinite() && owner.P.allFinite() &&
        object.at("dimension").as_int64() == dimension && object.at("epoch").as_int64() == static_cast<long int>(owner.time.bigTime) &&
        gecArString(object.at("float_root"),false) == expectedRoot &&
        gecArString(object.at("state_order"),false) == expectedOrder,
        "GEC_AR_SECTION_FLOAT_CORE_IDENTITY_INVALID");
    map<string,std::shared_ptr<ZhangIntegerDecisionProof>> bundleProofAtoms;
    for (const auto& item : gecArArray(object.at("evidence")))
    {
        ZhangCurrentEpochEvidenceKey key;
        auto value = gecArRestoreEvidence(item,runtime,owner,dimension,expectedRoot,expectedOrder,key,bundleProofAtoms);
        gecArRequire(evidence.emplace(std::move(key),std::move(value)).second,"GEC_AR_DUPLICATE_EVIDENCE");
    }
    for (const auto& item : gecArArray(object.at("authorities")))
    {
        const auto& row = gecArObject(item,{"owner","epoch","float_root","state_order","systems"});
        const auto key = gecArString(row.at("owner"),false);
        gecArRequire(zhangAmbresRuntimeBelongsTo(key,runtime),"GEC_AR_FOREIGN_AUTHORITY_OWNER");
        ZhangMultiGnssCurrentEvidenceAuthority authority;
        authority.epochGpstSeconds = row.at("epoch").as_int64();
        authority.floatRoot = gecArString(row.at("float_root"),false);
        authority.floatStateOrderId = gecArString(row.at("state_order"),false);
        gecArRequire(authority.epochGpstSeconds == static_cast<long int>(owner.time.bigTime) &&
            authority.floatRoot == gecArString(object.at("float_root")) &&
            authority.floatStateOrderId == gecArString(object.at("state_order")),"GEC_AR_AUTHORITY_SCOPE_INVALID");
        const auto& systems = row.at("systems").as_object();
        gecArRequire(!systems.empty() && systems.size() <= 3,"GEC_AR_AUTHORITY_SYSTEM_COUNT_INVALID");
        for (const auto& entry : systems)
        {
            const auto system = gecArSystem(gecjson::string(entry.key())); const auto& ids = gecArArray(entry.value());
            gecArRequire(ids.size() == 2,"GEC_AR_AUTHORITY_IDENTITIES_INVALID");
            const auto found = evidence.find({key,system});
            gecArRequire(found != evidence.end() && found->second.valid &&
                found->second.posteriorId == gecArString(ids[0],false) &&
                found->second.chartId == gecArString(ids[1],false),"GEC_AR_AUTHORITY_EVIDENCE_MISSING");
            authority.systemIdentities.emplace(system,std::make_pair(gecArString(ids[0]),gecArString(ids[1])));
        }
        gecArRequire(authorities.emplace(key,std::move(authority)).second,"GEC_AR_DUPLICATE_AUTHORITY");
    }
}
} // namespace

bool exportZhangGecArRuntimeSection(const KFState& owner, const string& runtimeId,
    string& payload, string& failureReason)
{
    std::lock_guard<std::mutex> lock(zhangAmbresCheckpointMutex);
    payload.clear(); failureReason.clear();
    try
    {
        gecArRequire(zhangMultiGnssArCheckpointActive(acsConfig),"GEC_AR_SECTION_REQUIRES_EXPLICIT_AR_MODE");
        gecArRequire(validZhangAmbresRuntimeId(runtimeId) && zhangAmbresRuntimeId(owner) == runtimeId,
            "GEC_AR_SECTION_RUNTIME_ID_INVALID");
        auto document = gecArSnapshot(owner,runtimeId);
        GecEvidenceMap evidence; GecAuthorityMap authorities;
        gecArRestoreSnapshot(document,owner,runtimeId,owner.x.size(),evidence,authorities);
        payload = gecjson::serialize(document);
        gecArRequire(zhangGecArRuntimePayloadSizeValid(payload.size()),"GEC_AR_SECTION_SIZE_INVALID");
        return true;
    }
    catch (const std::exception& error)
    { failureReason = string("GEC_AR_SECTION_EXPORT_FAILED:")+error.what(); payload.clear(); return false; }
}

bool importZhangGecArRuntimeSection(KFState& owner, const string& runtimeId,
    const string& payload, bool preflightOnly, int stateDimension, string& failureReason)
{
    std::lock_guard<std::mutex> lock(zhangAmbresCheckpointMutex);
    failureReason.clear();
    try
    {
        gecArRequire(zhangMultiGnssArCheckpointActive(acsConfig),"GEC_AR_SECTION_REQUIRES_EXPLICIT_AR_MODE");
        gecArRequire(validZhangAmbresRuntimeId(runtimeId) && zhangAmbresRuntimeId(owner) == runtimeId,
            "GEC_AR_SECTION_RUNTIME_ID_INVALID");
        gecArRequire(zhangGecArRuntimePayloadSizeValid(payload.size()),"GEC_AR_SECTION_SIZE_INVALID");
        for (const auto& [key,value] : zhangCurrentEpochEvidenceRegistry())
            gecArRequire(zhangAmbresRuntimeBelongsTo(key.first,runtimeId),"GEC_AR_FOREIGN_LIVE_EVIDENCE_OWNER");
        for (const auto& [key,value] : zhangMultiGnssCurrentEvidenceAuthorities())
            gecArRequire(zhangAmbresRuntimeBelongsTo(key,runtimeId),"GEC_AR_FOREIGN_LIVE_AUTHORITY_OWNER");
        GecEvidenceMap evidence; GecAuthorityMap authorities;
        gecArRestoreSnapshot(gecjson::parse(payload),owner,runtimeId,stateDimension,evidence,authorities);
        if (preflightOnly) return true;
        // All parsing, scope/proof checks, and allocations finish before the
        // noexcept swaps. Failure never partially edits either live registry.
        zhangCurrentEpochEvidenceRegistry().swap(evidence);
        zhangMultiGnssCurrentEvidenceAuthorities().swap(authorities);
        return true;
    }
    catch (const std::exception& error)
    { failureReason = string("GEC_AR_SECTION_IMPORT_FAILED:")+error.what(); return false; }
}

bool zhangGecArRuntimeRegistriesEmpty()
{
    std::lock_guard<std::mutex> lock(zhangAmbresCheckpointMutex);
    return zhangCurrentEpochEvidenceRegistry().empty() &&
        zhangMultiGnssCurrentEvidenceAuthorities().empty();
}

ZhangAmbresCheckpointResult exportZhangAmbresCheckpointSection(
    const KFState& owner,
    const string& runtimeId,
    string& payload)
{
    std::lock_guard<std::mutex> lock(zhangAmbresCheckpointMutex);
    ZhangAmbresCheckpointResult result;
    payload.clear();
    if (!validZhangAmbresRuntimeId(runtimeId))
    {
        result.failureReason = "AMBRES_CHECKPOINT_RUNTIME_ID_INVALID";
        return result;
    }
    if (zhangAmbresRuntimeId(owner) != runtimeId)
    {
        result.failureReason = "AMBRES_CHECKPOINT_OWNER_RUNTIME_ID_MISMATCH";
        return result;
    }
    if (ambresHasForeignRuntime(runtimeId))
    {
        result.failureReason = "AMBRES_CHECKPOINT_RUNTIME_OWNER_CONFLICT";
        return result;
    }

    AmbresCheckpointEnvelope snapshot;
    snapshot.runtimeId = runtimeId;
    snapshot.filterError = ::filterError;
    snapshot.transactionalConditioningFailed =
        zhangTransactionalConditioningFailed;
    snapshot.transactionalConditioningReason =
        zhangTransactionalConditioningReason;
    snapshot.arVerbose = AR_VERBO;
    snapshot.fixedAmbiguityVariance = FIXED_AMB_VAR;

    if(zhangR51Enabled())for(const auto& [key,store]:r51PhysicalArchives) {
        if(!zhangAmbresRuntimeBelongsTo(key.first,runtimeId))continue;
        R51PhysicalArchiveCheckpoint item;item.owner=key.first;item.system=key.second;item.valid=store.valid;
        for(const auto& row:store.rows){map<string,string> terms;for(const auto& [id,x]:row)terms[id]=x.convert_to<string>();item.rows.push_back(terms);}
        for(const auto& x:store.values)item.values.push_back(x.convert_to<string>());
        const auto closure=zhangDecisionRiskClosure(store.parents);
        if(!closure.valid){result.failureReason="R51_ARCHIVE_PROOF_INVALID";return result;}
        for(const auto& root:store.parents)item.roots.push_back(root->id);
        for(const auto& [id,proof]:closure.atoms){R51ProofAtomCheckpoint atom;
            atom.id=id;atom.statement=proof->originalStatement;atom.provenance=proof->observationProvenance;atom.risk=proof->conditionalFailureBound;
            for(const auto& parent:proof->parents)atom.parents.push_back(parent->id);item.proofs.push_back(atom);}
        snapshot.r51Archives.push_back(std::move(item));
    }

    for (const auto& [key, runtime] : zhangE27WideLaneRuntimes)
    {
        if (!zhangAmbresRuntimeBelongsTo(key.runtimeId, runtimeId))
        {
            continue;
        }
        if (runtime.accumulator.factorCount() != runtime.rawFactors.size())
        {
            result.failureReason =
                "AMBRES_CHECKPOINT_E27_LIVE_FACTOR_LEDGER_MISMATCH";
            return result;
        }
        AmbresE27RuntimeCheckpoint item;
        item.ownerRuntimeId = key.runtimeId;
        item.receiver = key.receiver;
        item.system = key.system;
        item.solution = key.solution;
        for (const auto& [satellite, time] : runtime.lastValid)
        {
            item.lastValid.push_back({
                satellite, captureZhangCheckpointTime(time)});
        }
        item.arcVersions = runtime.arcVersion;
        for (const auto& [satellite, datum] : runtime.productDatum)
        {
            item.productDatum.push_back({
                satellite, std::get<0>(datum), std::get<1>(datum),
                std::get<2>(datum), std::get<3>(datum)});
        }
        item.reference = runtime.reference;
        item.accumulatorArcVersions = runtime.accumulatorArcVersions;
        for (const auto& factor : runtime.rawFactors)
        {
            item.factors.push_back({
                captureZhangCheckpointTime(factor.time), factor.satellites,
                factor.values, factor.covariance, factor.noiseKeys,
                factor.noiseVariances, factor.satelliteNoiseDesign});
        }
        snapshot.e27Runtimes.push_back(std::move(item));
    }

    for (const auto& [identity, lattice] : zhangPersistentHeldLattices)
    {
        if (!zhangAmbresRuntimeBelongsTo(identity.first, runtimeId))
        {
            continue;
        }
        AmbresHeldLatticeCheckpoint item;
        item.ownerRuntimeId = identity.first;
        item.system = identity.second;
        item.lastEventId = lattice.lastEventId;
        item.consistent = lattice.consistent;
        for (const auto& row : lattice.rows)
        {
            item.rows.push_back(ambresCheckpointHeldRow(row));
        }
        snapshot.heldLattices.push_back(std::move(item));
    }

    for (const auto& [identity, rows] : zhangPersistentHeldEvidence)
    {
        if (!zhangAmbresRuntimeBelongsTo(identity.first, runtimeId))
        {
            continue;
        }
        for (const auto& [row, evidence] : rows)
        {
            snapshot.heldEvidence.push_back({
                identity.first, identity.second,
                ambresCheckpointHeldRow(row), evidence.confirmations,
                evidence.lastEpoch});
        }
    }

    for (const auto& [identity, state] :
         zhangProductRelationAdmissionStateRegistry())
    {
        if (!zhangAmbresRuntimeBelongsTo(identity.first, runtimeId))
        {
            continue;
        }
        string stateReason;
        if (!ProductRelationAdmission::validateState(state, &stateReason))
        {
            result.failureReason =
                "AMBRES_CHECKPOINT_ADMISSION_STATE_INVALID:" + stateReason;
            return result;
        }
        AmbresAdmissionStateCheckpoint item;
        item.ownerRuntimeId = identity.first;
        item.system = identity.second;
        for (const auto& row : state.certifiedRows)
        {
            item.certifiedRows.push_back(
                ambresCheckpointSparseExactRow(row));
        }
        for (const auto& value : state.certifiedValues)
        {
            item.certifiedValues.push_back(ambresExactText(value));
        }
        item.certifiedSatellites.assign(
            state.certifiedSatellites.begin(),
            state.certifiedSatellites.end());
        item.certifiedRelationIds.assign(
            state.certifiedRelationIds.begin(),
            state.certifiedRelationIds.end());
        for (const auto& [relationId, candidate] : state.pendingCandidates)
        {
            item.pendingCandidates.push_back(
                ambresCheckpointAdmissionCandidate(candidate));
        }
        for (const auto& [relationKey, confirmation] :
             state.temporalCertificateConfirmations)
        {
            AmbresTemporalConfirmationCheckpoint saved;
            saved.relationKey = relationKey;
            saved.integerValue = ambresExactText(confirmation.integerValue);
            saved.consistentEpochs = confirmation.consistentEpochs;
            saved.firstEpoch = confirmation.firstEpoch;
            saved.lastEpoch = confirmation.lastEpoch;
            saved.supportSignatures.assign(
                confirmation.supportSignatures.begin(),
                confirmation.supportSignatures.end());
            saved.redundancyConfirmed = confirmation.redundancyConfirmed;
            item.temporalConfirmations.push_back(std::move(saved));
        }
        snapshot.admissionStates.push_back(std::move(item));
    }

    for (const auto& [identity, held] : zhangHeldUserWideLaneRegistry)
    {
        const auto& [ownerRuntimeId, receiver, system, ifAcceptance] =
            identity;
        if (!zhangAmbresRuntimeBelongsTo(ownerRuntimeId, runtimeId))
        {
            continue;
        }
        if (held.integers.empty())
        {
            continue;
        }
        AmbresHeldUserWideLaneCheckpoint item;
        item.ownerRuntimeId = ownerRuntimeId;
        item.receiver = receiver;
        item.system = system;
        item.ifAcceptance = ifAcceptance;
        item.reference = held.reference;
        for (const auto& [satellite, value] : held.integers)
        {
            item.integers.push_back({satellite, ambresExactText(value)});
        }
        snapshot.heldUserWideLanes.push_back(std::move(item));
    }

    for (const auto& [ownerRuntimeId, transitions] :
         zhangActiveTemporalBesdTransitions)
    {
        if (!zhangAmbresRuntimeBelongsTo(ownerRuntimeId, runtimeId))
        {
            continue;
        }
        AmbresActiveTransitionsCheckpoint item;
        item.ownerRuntimeId = ownerRuntimeId;
        for (const auto& transition : transitions)
        {
            item.transitions.push_back(ambresCheckpointPending(transition));
        }
        snapshot.activeTransitions.push_back(std::move(item));
    }

    for (const auto& [ownerRuntimeId, moments] : zhangRelinkPriorMoments)
    {
        if (!zhangAmbresRuntimeBelongsTo(ownerRuntimeId, runtimeId))
        {
            continue;
        }
        for (const auto& [key, moment] : moments)
        {
            snapshot.relinkPriors.push_back({
                ownerRuntimeId, key.system, key.anchor, key.satellite,
                captureZhangCheckpointTime(moment.time), moment.firstCode,
                moment.secondCode, moment.productDatumVersion,
                moment.firstMean, moment.firstVariance, moment.wideLaneMean,
                moment.wideLaneVariance, moment.firstWideLaneCovariance});
        }
    }

    for (const auto& [ownerRuntimeId, accumulators] :
         zhangRelinkShadowAccumulators)
    {
        if (!zhangAmbresRuntimeBelongsTo(ownerRuntimeId, runtimeId))
        {
            continue;
        }
        for (const auto& [key, accumulator] : accumulators)
        {
            AmbresScalarAccumulatorCheckpoint item;
            item.ownerRuntimeId = ownerRuntimeId;
            item.key = key;
            item.integerHypothesis = accumulator.integerHypothesis;
            item.initialized = accumulator.initialized;
            for (const auto& increment : accumulator.increments)
            {
                item.increments.push_back({
                    captureZhangCheckpointTime(increment.time),
                    increment.information, increment.natural});
            }
            snapshot.scalarAccumulators.push_back(std::move(item));
        }
    }

    for (const auto& [ownerRuntimeId, accumulators] :
         zhangRelinkJointAccumulators)
    {
        if (!zhangAmbresRuntimeBelongsTo(ownerRuntimeId, runtimeId))
        {
            continue;
        }
        for (const auto& [key, accumulator] : accumulators)
        {
            AmbresJointAccumulatorCheckpoint item;
            item.ownerRuntimeId = ownerRuntimeId;
            item.key = key;
            item.wideLaneHypothesis = accumulator.wideLaneHypothesis;
            item.firstHypothesis = accumulator.firstHypothesis;
            item.initialized = accumulator.initialized;
            for (const auto& increment : accumulator.increments)
            {
                AmbresJointIncrementCheckpoint saved;
                saved.time = captureZhangCheckpointTime(increment.time);
                saved.information = {
                    increment.information(0, 0), increment.information(0, 1),
                    increment.information(1, 0), increment.information(1, 1)};
                saved.natural = {
                    increment.natural(0), increment.natural(1)};
                item.increments.push_back(std::move(saved));
            }
            snapshot.jointAccumulators.push_back(std::move(item));
        }
    }

    for (const auto& [ownerRuntimeId, accumulators] :
         zhangWhitenedWlAccumulators)
    {
        if (!zhangAmbresRuntimeBelongsTo(ownerRuntimeId, runtimeId))
        {
            continue;
        }
        for (const auto& [key, accumulator] : accumulators)
        {
            AmbresWhitenedAccumulatorCheckpoint item;
            item.ownerRuntimeId = ownerRuntimeId;
            item.key = key;
            item.physicalSegmentIdentity =
                accumulator.physicalSegmentIdentity;
            item.lastGraphEventId = accumulator.lastGraphEventId;
            item.lastProductDatumVersion =
                accumulator.lastProductDatumVersion;
            item.initialized = accumulator.initialized;
            for (const auto& observation : accumulator.observations)
            {
                item.observations.push_back({
                    captureZhangCheckpointTime(observation.time),
                    observation.value, observation.variance});
            }
            for (const auto& [arc, coefficient] : accumulator.physicalRow)
            {
                item.physicalRow.push_back(
                    ambresCheckpointArc(arc, coefficient));
            }
            snapshot.whitenedAccumulators.push_back(std::move(item));
        }
    }

    if (!multiGnssFloatAmbresCheckpointEmpty(snapshot, result.failureReason))
    {
        return result;
    }
    if (!ambresGecArchiveScope(snapshot, result.failureReason)) return result;
    AmbresRestoredRuntime preflight;
    ZhangAmbresCheckpointResult preflightResult;
    if (!ambresRestoreEnvelope(
            snapshot, runtimeId, preflight, preflightResult))
    {
        return preflightResult;
    }
    try
    {
        std::ostringstream output(std::ios::binary | std::ios::out);
        boost::archive::binary_oarchive archive(
            output, boost::archive::no_header);
        archive << snapshot;
        payload = output.str();
    }
    catch (const std::exception& exception)
    {
        result.failureReason =
            string("AMBRES_CHECKPOINT_SERIALIZE_FAILED:") + exception.what();
        payload.clear();
        return result;
    }
    if (payload.empty())
    {
        result.failureReason = "AMBRES_CHECKPOINT_EMPTY_PAYLOAD";
        return result;
    }
    preflightResult.valid = true;
    return preflightResult;
}

ZhangAmbresCheckpointResult importZhangAmbresCheckpointSection(
    KFState& owner,
    const string& runtimeId,
    const string& payload,
    bool validateOnly)
{
    std::lock_guard<std::mutex> lock(zhangAmbresCheckpointMutex);
    ZhangAmbresCheckpointResult result;
    if (!validZhangAmbresRuntimeId(runtimeId))
    {
        result.failureReason = "AMBRES_CHECKPOINT_RUNTIME_ID_INVALID";
        return result;
    }
    if (payload.empty())
    {
        result.failureReason = "AMBRES_CHECKPOINT_EMPTY_PAYLOAD";
        return result;
    }
    if (zhangAmbresRuntimeId(owner) != runtimeId)
    {
        result.failureReason = "AMBRES_CHECKPOINT_OWNER_RUNTIME_ID_MISMATCH";
        return result;
    }
    if (ambresHasForeignRuntime(runtimeId))
    {
        result.failureReason = "AMBRES_CHECKPOINT_RUNTIME_OWNER_CONFLICT";
        return result;
    }

    AmbresCheckpointEnvelope snapshot;
    try
    {
        std::istringstream input(payload, std::ios::binary | std::ios::in);
        {
            boost::archive::binary_iarchive archive(
                input, boost::archive::no_header);
            archive >> snapshot;
        }
        if (input.peek() != std::char_traits<char>::eof())
        {
            result.failureReason = "AMBRES_CHECKPOINT_TRAILING_BYTES";
            return result;
        }
    }
    catch (const std::exception& exception)
    {
        result.failureReason =
            string("AMBRES_CHECKPOINT_DESERIALIZE_FAILED:") + exception.what();
        return result;
    }

    if (!multiGnssFloatAmbresCheckpointEmpty(snapshot, result.failureReason))
    {
        return result;
    }
    if (!ambresGecArchiveScope(snapshot, result.failureReason)) return result;
    decltype(r51PhysicalArchives) restoredR51;
    try {
        for(const auto& item:snapshot.r51Archives) {
            if(!zhangAmbresRuntimeBelongsTo(item.owner,runtimeId) || item.rows.size()!=item.values.size() ||
                item.rows.size()>ZHANG_AMBRES_MAX_CONTAINER || item.proofs.size()>ZHANG_AMBRES_MAX_CONTAINER)
                throw std::runtime_error("R51_ARCHIVE_INVALID");
            std::map<string,std::shared_ptr<ZhangIntegerDecisionProof>> atoms;
            for(const auto& dto:item.proofs){auto atom=std::make_shared<ZhangIntegerDecisionProof>();
                atom->id=dto.id;atom->originalStatement=dto.statement;atom->observationProvenance=dto.provenance;atom->conditionalFailureBound=dto.risk;
                if(!atoms.emplace(dto.id,atom).second)throw std::runtime_error("R51_PROOF_DUPLICATE");}
            for(const auto& dto:item.proofs)for(const auto& parent:dto.parents)atoms.at(dto.id)->parents.push_back(atoms.at(parent));
            ZhangR51PhysicalStore store;store.valid=item.valid;
            for(const auto& id:item.roots)store.parents.push_back(atoms.at(id));
            if(!zhangDecisionRiskClosure(store.parents).valid)throw std::runtime_error("R51_PROOF_CLOSURE_INVALID");
            for(int i=0;i<item.rows.size();++i){map<string,ZhangExactInteger> row;
                for(const auto& [id,x]:item.rows[i])row[id]=ZhangExactInteger(x);
                store.rows.push_back(row);store.values.push_back(ZhangExactInteger(item.values[i]));}
            if(!restoredR51.emplace(std::make_pair(item.owner,item.system),std::move(store)).second)
                throw std::runtime_error("R51_ARCHIVE_DUPLICATE");
        }
    } catch(const std::exception& error){result.failureReason=error.what();return result;}

    AmbresRestoredRuntime restored;
    if (!ambresRestoreEnvelope(snapshot, runtimeId, restored, result))
    {
        return result;
    }
    if (validateOnly)
    {
        result.valid = true;
        return result;
    }

    AmbresRestoredRuntime committed;
    string committedConditioningReason;
    try
    {
        ambresPrepareCommittedMaps(runtimeId, restored, committed);
        committedConditioningReason =
            snapshot.transactionalConditioningReason;
    }
    catch (const std::exception& exception)
    {
        result.failureReason =
            string("AMBRES_CHECKPOINT_COMMIT_PREPARE_FAILED:") +
            exception.what();
        return result;
    }

    if(zhangR51Enabled())r51PhysicalArchives.swap(restoredR51);
    zhangE27WideLaneRuntimes.swap(committed.e27Runtimes);
    zhangPersistentHeldLattices.swap(committed.heldLattices);
    zhangPersistentHeldEvidence.swap(committed.heldEvidence);
    zhangProductRelationAdmissionStateRegistry().swap(
        committed.admissionStates);
    zhangHeldUserWideLaneRegistry.swap(committed.heldUserWideLanes);
    zhangActiveTemporalBesdTransitions.swap(committed.activeTransitions);
    zhangRelinkPriorMoments.swap(committed.relinkPriors);
    zhangRelinkShadowAccumulators.swap(committed.scalarAccumulators);
    zhangRelinkJointAccumulators.swap(committed.jointAccumulators);
    zhangWhitenedWlAccumulators.swap(committed.whitenedAccumulators);
    ::filterError = snapshot.filterError;
    zhangTransactionalConditioningFailed =
        snapshot.transactionalConditioningFailed;
    zhangTransactionalConditioningReason.swap(
        committedConditioningReason);
    AR_VERBO = snapshot.arVerbose;
    FIXED_AMB_VAR = snapshot.fixedAmbiguityVariance;
    result.valid = true;
    return result;
}

ZhangAmbresCheckpointResult inspectZhangAmbresCheckpointSnapshotReferences(
    const string& runtimeId,
    const string& payload,
    ZhangCheckpointSnapshotReferenceSummary& summary)
{
    std::lock_guard<std::mutex> lock(zhangAmbresCheckpointMutex);
    summary = {};
    summary.sectionName = ZHANG_AMBRES_CHECKPOINT_SECTION_NAME;
    summary.runtimeId = runtimeId;
    ZhangAmbresCheckpointResult result;
    auto fail = [&](const string& reason)
    {
        summary.valid = false;
        summary.failureReason = reason;
        result.valid = false;
        result.failureReason = reason;
        return result;
    };
    if (!validZhangAmbresRuntimeId(runtimeId))
    {
        return fail("AMBRES_CHECKPOINT_RUNTIME_ID_INVALID");
    }
    if (payload.empty())
    {
        return fail("AMBRES_CHECKPOINT_EMPTY_PAYLOAD");
    }

    AmbresCheckpointEnvelope snapshot;
    try
    {
        std::istringstream input(payload, std::ios::binary | std::ios::in);
        {
            boost::archive::binary_iarchive archive(
                input, boost::archive::no_header);
            archive >> snapshot;
        }
        if (input.peek() != std::char_traits<char>::eof())
        {
            return fail("AMBRES_CHECKPOINT_TRAILING_BYTES");
        }
    }
    catch (const std::exception& exception)
    {
        return fail(
            string("AMBRES_CHECKPOINT_DESERIALIZE_FAILED:") +
            exception.what());
    }

    if (snapshot.schemaVersion != ZHANG_AMBRES_CHECKPOINT_SCHEMA_VERSION)
    {
        return fail("AMBRES_CHECKPOINT_SCHEMA_MISMATCH");
    }
    if (snapshot.runtimeId != runtimeId)
    {
        return fail("AMBRES_CHECKPOINT_RUNTIME_ID_MISMATCH");
    }
    if (snapshot.activeTransitions.size() > ZHANG_AMBRES_MAX_CONTAINER)
    {
        return fail("AMBRES_CHECKPOINT_ACTIVE_TRANSITION_CONTAINER_INVALID");
    }
    auto accountReference = [&](const string& identity, int declaredCount)
    {
        if (identity.empty())
        {
            return declaredCount == 0;
        }
        if (declaredCount <= 0)
        {
            return false;
        }
        auto& actual = summary.actualReferenceCounts[identity];
        if (actual == std::numeric_limits<std::size_t>::max())
        {
            return false;
        }
        actual++;
        // AMBRES computes this diagnostic within each branch.  Different
        // branch owners may legitimately report different positive counts
        // for the same physical snapshot, so bundle integrity is based on
        // the recomputed endpoint inventory, not this local diagnostic.
        if (declaredCount > 0)
        {
            summary.declaredReferenceCounts.emplace(
                identity, static_cast<std::size_t>(declaredCount));
        }
        return true;
    };
    set<string> ownerRuntimeIds;
    for (const auto& active : snapshot.activeTransitions)
    {
        if (!validZhangAmbresRuntimeId(active.ownerRuntimeId) ||
            !zhangAmbresRuntimeBelongsTo(active.ownerRuntimeId, runtimeId) ||
            active.transitions.size() > ZHANG_AMBRES_MAX_CONTAINER ||
            !ownerRuntimeIds.insert(active.ownerRuntimeId).second)
        {
            return fail("AMBRES_CHECKPOINT_ACTIVE_TRANSITION_HEADER_INVALID");
        }
        for (const auto& transition : active.transitions)
        {
            if (!accountReference(
                    transition.oldSnapshotIdentity,
                    transition.oldSnapshotReferenceCount) ||
                !accountReference(
                    transition.newSnapshotIdentity,
                    transition.newSnapshotReferenceCount))
            {
                return fail(
                    "AMBRES_CHECKPOINT_ACTIVE_SNAPSHOT_REFERENCE_INVALID");
            }
            summary.transitionCount++;
        }
    }
    summary.valid = true;
    summary.failureReason.clear();
    result.valid = true;
    return result;
}

namespace
{
bool ambresValidOwner(
    const string& ownerRuntimeId,
    const string& rootRuntimeId,
    string& reason)
{
    if (!validZhangAmbresRuntimeId(ownerRuntimeId) ||
        !zhangAmbresRuntimeBelongsTo(ownerRuntimeId, rootRuntimeId))
    {
        reason = "AMBRES_CHECKPOINT_OWNER_RUNTIME_ID_INVALID";
        return false;
    }
    return true;
}

bool ambresValidateFactor(
    const AmbresE27RawFactorCheckpoint& snapshot,
    ZhangE27WideLaneRawFactor& factor,
    string& reason)
{
    if (!ambresRestoreFiniteTime(snapshot.time, factor.time) ||
        snapshot.satellites.size() < 2 ||
        snapshot.satellites.size() > ZHANG_E27_MAX_SATELLITE ||
        snapshot.values.size() != static_cast<int>(snapshot.satellites.size()) ||
        snapshot.covariance.rows() != snapshot.values.size() ||
        snapshot.covariance.cols() != snapshot.values.size() ||
        snapshot.noiseKeys.size() !=
            static_cast<std::size_t>(snapshot.noiseVariances.size()) ||
        snapshot.satelliteNoiseDesign.rows() != snapshot.values.size() ||
        snapshot.satelliteNoiseDesign.cols() != snapshot.noiseVariances.size() ||
        !snapshot.values.allFinite() || !snapshot.covariance.allFinite() ||
        !snapshot.noiseVariances.allFinite() ||
        !snapshot.satelliteNoiseDesign.allFinite() ||
        (snapshot.noiseVariances.array() <= 0).any())
    {
        reason = "AMBRES_CHECKPOINT_E27_FACTOR_DIMENSION_OR_VALUE_INVALID";
        return false;
    }
    set<int> satellites;
    for (int satellite : snapshot.satellites)
    {
        if (satellite <= 0 || satellite > ZHANG_E27_MAX_SATELLITE ||
            !satellites.insert(satellite).second)
        {
            reason = "AMBRES_CHECKPOINT_E27_FACTOR_SATELLITE_INVALID";
            return false;
        }
    }
    set<string> noiseKeys;
    for (const string& key : snapshot.noiseKeys)
    {
        if (!ambresValidString(key, false) || !noiseKeys.insert(key).second)
        {
            reason = "AMBRES_CHECKPOINT_E27_FACTOR_NOISE_KEY_INVALID";
            return false;
        }
    }
    const MatrixXd symmetric =
        0.5 * (snapshot.covariance + snapshot.covariance.transpose());
    const double scale = std::max(1.0, symmetric.cwiseAbs().maxCoeff());
    if ((snapshot.covariance - snapshot.covariance.transpose()).norm() >
        1e-10 * scale)
    {
        reason = "AMBRES_CHECKPOINT_E27_FACTOR_COVARIANCE_ASYMMETRIC";
        return false;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(symmetric);
    if (solver.info() != Eigen::Success ||
        !solver.eigenvalues().allFinite() ||
        solver.eigenvalues().minCoeff() < -1e-10 * scale)
    {
        reason = "AMBRES_CHECKPOINT_E27_FACTOR_COVARIANCE_NOT_PSD";
        return false;
    }
    factor.time = restoreZhangCheckpointTime(snapshot.time);
    factor.satellites = snapshot.satellites;
    factor.values = snapshot.values;
    factor.covariance = snapshot.covariance;
    factor.noiseKeys = snapshot.noiseKeys;
    factor.noiseVariances = snapshot.noiseVariances;
    factor.satelliteNoiseDesign = snapshot.satelliteNoiseDesign;
    return true;
}

bool ambresRestoreE27Runtime(
    const AmbresE27RuntimeCheckpoint& snapshot,
    const string& rootRuntimeId,
    ZhangE27WideLaneKey& key,
    ZhangE27WideLaneRuntime& runtime,
    string& reason)
{
    if (!ambresValidOwner(snapshot.ownerRuntimeId, rootRuntimeId, reason) ||
        !ambresValidString(snapshot.receiver, false) ||
        !ambresValidString(snapshot.solution, false) ||
        !ambresValidSystem(snapshot.system) ||
        snapshot.lastValid.size() > ZHANG_AMBRES_MAX_CONTAINER ||
        snapshot.productDatum.size() > ZHANG_AMBRES_MAX_CONTAINER ||
        snapshot.factors.size() > ZHANG_E27_MAX_FACTORS ||
        !ambresSatelliteUsesSystem(snapshot.reference, snapshot.system, true))
    {
        if (reason.empty())
        {
            reason = "AMBRES_CHECKPOINT_E27_RUNTIME_HEADER_INVALID";
        }
        return false;
    }
    key = {snapshot.ownerRuntimeId, snapshot.receiver, snapshot.system,
        snapshot.solution};
    for (const auto& item : snapshot.lastValid)
    {
        GTime time;
        if (!ambresSatelliteUsesSystem(item.satellite, snapshot.system) ||
            !ambresRestoreFiniteTime(item.time, time) ||
            !runtime.lastValid.emplace(item.satellite, time).second)
        {
            reason = "AMBRES_CHECKPOINT_E27_LAST_VALID_INVALID";
            return false;
        }
    }
    for (const auto& [satellite, version] : snapshot.arcVersions)
    {
        if (!ambresSatelliteUsesSystem(satellite, snapshot.system) || version <= 0)
        {
            reason = "AMBRES_CHECKPOINT_E27_ARC_VERSION_INVALID";
            return false;
        }
        runtime.arcVersion[satellite] = version;
    }
    for (const auto& item : snapshot.productDatum)
    {
        if (!ambresSatelliteUsesSystem(item.satellite, snapshot.system) ||
            item.firstDiscontinuity < 0 || item.firstDatumVersion < 0 ||
            item.secondDiscontinuity < 0 || item.secondDatumVersion < 0 ||
            !runtime.productDatum.emplace(item.satellite, std::make_tuple(
                item.firstDiscontinuity, item.firstDatumVersion,
                item.secondDiscontinuity, item.secondDatumVersion)).second)
        {
            reason = "AMBRES_CHECKPOINT_E27_PRODUCT_DATUM_INVALID";
            return false;
        }
    }
    if (runtime.lastValid.size() != runtime.arcVersion.size() ||
        runtime.lastValid.size() != runtime.productDatum.size())
    {
        reason = "AMBRES_CHECKPOINT_E27_SATELLITE_REGISTRY_INCOMPLETE";
        return false;
    }
    runtime.reference = snapshot.reference;
    runtime.accumulatorArcVersions = snapshot.accumulatorArcVersions;
    if (runtime.accumulatorArcVersions.size() != runtime.arcVersion.size())
    {
        reason = "AMBRES_CHECKPOINT_E27_ACCUMULATOR_ARC_REGISTRY_INCOMPLETE";
        return false;
    }
    for (const auto& [satellite, version] : runtime.accumulatorArcVersions)
    {
        if (satellite <= 0 || satellite > ZHANG_E27_MAX_SATELLITE ||
            version <= 0)
        {
            reason = "AMBRES_CHECKPOINT_E27_ACCUMULATOR_ARC_INVALID";
            return false;
        }
        SatSys sat(E_Sys::GPS, satellite);
        auto found = runtime.arcVersion.find(sat);
        if (found == runtime.arcVersion.end() || found->second != version)
        {
            reason = "AMBRES_CHECKPOINT_E27_ACCUMULATOR_ARC_CONFLICT";
            return false;
        }
        runtime.accumulator.setArcVersion(satellite, version);
    }
    GTime previousTime;
    bool havePrevious = false;
    for (const auto& item : snapshot.factors)
    {
        ZhangE27WideLaneRawFactor factor;
        if (!ambresValidateFactor(item, factor, reason) ||
            (havePrevious && factor.time < previousTime))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_E27_FACTOR_CHRONOLOGY_INVALID";
            }
            return false;
        }
        havePrevious = true;
        previousTime = factor.time;
        const double replayTime = static_cast<double>(factor.time.bigTime);
        if (!std::isfinite(replayTime))
        {
            reason = "AMBRES_CHECKPOINT_E27_FACTOR_TIME_OUT_OF_RANGE";
            return false;
        }
        runtime.accumulator.addEpoch(
            replayTime, factor.satellites,
            factor.values, factor.covariance, factor.noiseKeys,
            factor.noiseVariances, factor.satelliteNoiseDesign);
        runtime.rawFactors.push_back(std::move(factor));
    }
    if (runtime.accumulator.factorCount() != runtime.rawFactors.size())
    {
        reason = "AMBRES_CHECKPOINT_E27_FACTOR_REPLAY_MISMATCH";
        return false;
    }
    if (runtime.rawFactors.size() > 1 &&
        (runtime.rawFactors.back().time - runtime.rawFactors.front().time)
                .to_double() > ZHANG_E27_LAG_SECONDS + 1e-6)
    {
        reason = "AMBRES_CHECKPOINT_E27_FACTOR_WINDOW_NOT_PRUNED";
        return false;
    }
    return true;
}

bool ambresValidateChronology(
    const vector<AmbresScalarIncrementCheckpoint>& increments,
    std::deque<ZhangRelinkShadowIncrement>& restored,
    string& reason)
{
    if (increments.size() > ZHANG_AMBRES_MAX_CONTAINER)
    {
        reason = "AMBRES_CHECKPOINT_SCALAR_WINDOW_TOO_LARGE";
        return false;
    }
    GTime previous;
    bool havePrevious = false;
    for (const auto& item : increments)
    {
        GTime time;
        if (!ambresRestoreFiniteTime(item.time, time) ||
            !std::isfinite(item.information) || item.information <= 0 ||
            !std::isfinite(item.natural) ||
            (havePrevious && !(previous < time)) ||
            (havePrevious &&
             acsConfig.zhangPppAr.multi_epoch_relink_shadow_max_gap_seconds > 0 &&
             (time - previous).to_double() >
                acsConfig.zhangPppAr.multi_epoch_relink_shadow_max_gap_seconds))
        {
            reason = "AMBRES_CHECKPOINT_SCALAR_INCREMENT_INVALID";
            return false;
        }
        restored.push_back({time, item.information, item.natural});
        previous = time;
        havePrevious = true;
    }
    return true;
}

bool ambresValidateJointChronology(
    const vector<AmbresJointIncrementCheckpoint>& increments,
    std::deque<ZhangRelinkJointIncrement>& restored,
    string& reason)
{
    if (increments.size() > ZHANG_AMBRES_MAX_CONTAINER)
    {
        reason = "AMBRES_CHECKPOINT_JOINT_WINDOW_TOO_LARGE";
        return false;
    }
    GTime previous;
    bool havePrevious = false;
    for (const auto& item : increments)
    {
        GTime time;
        Matrix2d information;
        Vector2d natural;
        information << item.information[0], item.information[1],
            item.information[2], item.information[3];
        natural << item.natural[0], item.natural[1];
        const double scale = std::max(1.0, information.cwiseAbs().maxCoeff());
        Eigen::SelfAdjointEigenSolver<Matrix2d> solver(
            0.5 * (information + information.transpose()));
        if (!ambresRestoreFiniteTime(item.time, time) ||
            !information.allFinite() || !natural.allFinite() ||
            std::abs(information(0, 1) - information(1, 0)) > 1e-10 * scale ||
            solver.info() != Eigen::Success ||
            solver.eigenvalues().minCoeff() < -1e-10 * scale ||
            solver.eigenvalues().maxCoeff() <= 0 ||
            (havePrevious && !(previous < time)) ||
            (havePrevious &&
             acsConfig.zhangPppAr.multi_epoch_relink_shadow_max_gap_seconds > 0 &&
             (time - previous).to_double() >
                acsConfig.zhangPppAr.multi_epoch_relink_shadow_max_gap_seconds))
        {
            reason = "AMBRES_CHECKPOINT_JOINT_INCREMENT_INVALID";
            return false;
        }
        restored.push_back({time, information, natural});
        previous = time;
        havePrevious = true;
    }
    return true;
}

bool ambresRestoreEnvelope(
    const AmbresCheckpointEnvelope& snapshot,
    const string& rootRuntimeId,
    AmbresRestoredRuntime& restored,
    ZhangAmbresCheckpointResult& result)
{
    string& reason = result.failureReason;
    const auto containerValid = [](std::size_t size)
    {
        return size <= ZHANG_AMBRES_MAX_CONTAINER;
    };
    if (snapshot.schemaVersion != ZHANG_AMBRES_CHECKPOINT_SCHEMA_VERSION)
    {
        reason = "AMBRES_CHECKPOINT_SCHEMA_MISMATCH";
        return false;
    }
    if (snapshot.runtimeId != rootRuntimeId)
    {
        reason = "AMBRES_CHECKPOINT_RUNTIME_ID_MISMATCH";
        return false;
    }
    if (snapshot.filterError ||
        !ambresValidString(snapshot.transactionalConditioningReason) ||
        (snapshot.transactionalConditioningFailed !=
            !snapshot.transactionalConditioningReason.empty()) ||
        !std::isfinite(snapshot.fixedAmbiguityVariance) ||
        snapshot.fixedAmbiguityVariance <= 0 ||
        !containerValid(snapshot.e27Runtimes.size()) ||
        !containerValid(snapshot.heldLattices.size()) ||
        !containerValid(snapshot.heldEvidence.size()) ||
        !containerValid(snapshot.admissionStates.size()) ||
        !containerValid(snapshot.heldUserWideLanes.size()) ||
        !containerValid(snapshot.activeTransitions.size()) ||
        !containerValid(snapshot.relinkPriors.size()) ||
        !containerValid(snapshot.scalarAccumulators.size()) ||
        !containerValid(snapshot.jointAccumulators.size()) ||
        !containerValid(snapshot.whitenedAccumulators.size()))
    {
        reason = "AMBRES_CHECKPOINT_ENVELOPE_VALUE_INVALID";
        return false;
    }

    for (const auto& item : snapshot.e27Runtimes)
    {
        ZhangE27WideLaneKey key;
        ZhangE27WideLaneRuntime runtime;
        if (!ambresRestoreE27Runtime(
                item, rootRuntimeId, key, runtime, reason) ||
            !restored.e27Runtimes.emplace(key, std::move(runtime)).second)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_E27_RUNTIME_DUPLICATE";
            }
            return false;
        }
        result.e27Factors += item.factors.size();
    }
    result.e27Runtimes = restored.e27Runtimes.size();

    for (const auto& item : snapshot.heldLattices)
    {
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            !ambresValidSystem(item.system) || item.lastEventId < 0 ||
            !item.consistent || item.rows.size() > ZHANG_AMBRES_MAX_CONTAINER)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_HELD_LATTICE_HEADER_INVALID";
            }
            return false;
        }
        ZhangPersistentHeldLattice lattice;
        lattice.lastEventId = item.lastEventId;
        lattice.consistent = item.consistent;
        for (const auto& row : item.rows)
        {
            ZhangPersistentHeldRow restoredRow;
            if (!ambresRestoreHeldRow(row, restoredRow, reason) ||
                !ambresHeldRowUsesSystem(restoredRow, item.system))
            {
                return false;
            }
            lattice.rows.push_back(std::move(restoredRow));
        }
        normalisePersistentHeldLattice(lattice);
        if (!lattice.consistent ||
            !restored.heldLattices.emplace(
                std::make_pair(item.ownerRuntimeId, item.system),
                std::move(lattice)).second)
        {
            reason = "AMBRES_CHECKPOINT_HELD_LATTICE_INCONSISTENT_OR_DUPLICATE";
            return false;
        }
        result.heldRows += item.rows.size();
    }
    result.heldLattices = restored.heldLattices.size();

    for (const auto& item : snapshot.heldEvidence)
    {
        ZhangPersistentHeldRow row;
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            !ambresValidSystem(item.system) || item.confirmations <= 0 ||
            item.confirmations >= std::max(
                1, acsConfig.zhangPppAr.promotion_confirmation_epochs) ||
            item.lastEpoch <= 0 ||
            !ambresRestoreHeldRow(item.row, row, reason) ||
            !ambresHeldRowUsesSystem(row, item.system))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_HELD_EVIDENCE_INVALID";
            }
            return false;
        }
        auto& evidence = restored.heldEvidence[
            {item.ownerRuntimeId, item.system}];
        if (!evidence.emplace(
                std::move(row), ZhangPersistentHeldEvidence{
                    item.confirmations, item.lastEpoch}).second)
        {
            reason = "AMBRES_CHECKPOINT_HELD_EVIDENCE_DUPLICATE";
            return false;
        }
        result.heldEvidenceRows++;
    }

    for (const auto& item : snapshot.admissionStates)
    {
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason)
            || !ambresValidSystem(item.system)
            || !containerValid(item.certifiedRows.size())
            || !containerValid(item.certifiedValues.size())
            || !containerValid(item.certifiedSatellites.size())
            || !containerValid(item.certifiedRelationIds.size())
            || !containerValid(item.pendingCandidates.size())
            || !containerValid(item.temporalConfirmations.size())
            || item.certifiedRows.size() != item.certifiedValues.size())
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_ADMISSION_HEADER_INVALID";
            }
            return false;
        }
        ZhangProductRelationAdmissionState state;
        for (const auto& row : item.certifiedRows)
        {
            map<string, ZhangExactInteger> restoredRow;
            if (!ambresRestoreSparseExactRow(row, restoredRow, reason))
            {
                return false;
            }
            state.certifiedRows.push_back(std::move(restoredRow));
        }
        for (const auto& text : item.certifiedValues)
        {
            ZhangExactInteger value;
            if (!ambresParseExact(text, value, reason))
            {
                return false;
            }
            state.certifiedValues.push_back(std::move(value));
        }
        for (const auto& satellite : item.certifiedSatellites)
        {
            if (!ambresValidString(satellite, false)
                || (zhangMultiGnssArCheckpointActive(acsConfig) &&
                    !ambresSatelliteUsesSystem(SatSys(satellite.c_str()), item.system))
                || !state.certifiedSatellites.insert(satellite).second)
            {
                reason = "AMBRES_CHECKPOINT_ADMISSION_SATELLITE_INVALID";
                return false;
            }
        }
        for (const auto& relationId : item.certifiedRelationIds)
        {
            if (!ambresValidString(relationId, false)
                || !state.certifiedRelationIds.insert(relationId).second)
            {
                reason = "AMBRES_CHECKPOINT_ADMISSION_RELATION_ID_INVALID";
                return false;
            }
        }
        for (const auto& saved : item.pendingCandidates)
        {
            ZhangProductRelationAdmissionCandidate candidate;
            if (!ambresRestoreAdmissionCandidate(saved, candidate, reason, item.system)
                || !state.pendingCandidates.emplace(
                    candidate.relationId, std::move(candidate)).second)
            {
                if (reason.empty())
                {
                    reason =
                        "AMBRES_CHECKPOINT_ADMISSION_CANDIDATE_DUPLICATE";
                }
                return false;
            }
        }
        for (const auto& saved : item.temporalConfirmations)
        {
            TemporalCertificateConfirmationState confirmation;
            if (!ambresValidString(saved.relationKey, false)
                || saved.consistentEpochs <= 0
                || saved.firstEpoch <= 0
                || saved.lastEpoch < saved.firstEpoch
                || !containerValid(saved.supportSignatures.size())
                || !ambresParseExact(saved.integerValue,
                    confirmation.integerValue, reason))
            {
                if (reason.empty())
                {
                    reason =
                        "AMBRES_CHECKPOINT_TEMPORAL_CONFIRMATION_INVALID";
                }
                return false;
            }
            confirmation.canonicalRelationKey = saved.relationKey;
            confirmation.consistentEpochs = saved.consistentEpochs;
            confirmation.firstEpoch = saved.firstEpoch;
            confirmation.lastEpoch = saved.lastEpoch;
            confirmation.redundancyConfirmed = saved.redundancyConfirmed;
            for (const auto& signature : saved.supportSignatures)
            {
                if (!ambresValidString(signature, false)
                    || !confirmation.supportSignatures.insert(
                        signature).second)
                {
                    reason =
                        "AMBRES_CHECKPOINT_TEMPORAL_SUPPORT_INVALID";
                    return false;
                }
            }
            if (!state.temporalCertificateConfirmations.emplace(
                    saved.relationKey, std::move(confirmation)).second)
            {
                reason =
                    "AMBRES_CHECKPOINT_TEMPORAL_CONFIRMATION_DUPLICATE";
                return false;
            }
        }
        string stateReason;
        if (zhangMultiGnssArCheckpointActive(acsConfig))
        {
            const auto codes = zhangCurrentNetworkBaselineCodes(item.system);
            if (!codes || !zhangMultiGnssPhysicalRowsUseSystem(state.certifiedRows, item.system, *codes))
            { reason = "AMBRES_CHECKPOINT_ADMISSION_PHYSICAL_SCOPE_INVALID"; return false; }
        }
        if (!ProductRelationAdmission::validateState(state, &stateReason)
            || !restored.admissionStates.emplace(
                make_pair(item.ownerRuntimeId, item.system),
                std::move(state)).second)
        {
            reason = stateReason.empty()
                ? "AMBRES_CHECKPOINT_ADMISSION_STATE_DUPLICATE"
                : "AMBRES_CHECKPOINT_ADMISSION_STATE_INVALID:" + stateReason;
            return false;
        }
        result.productRelationPendingCandidates +=
            item.pendingCandidates.size();
        result.temporalCertificateConfirmations +=
            item.temporalConfirmations.size();
    }
    result.productRelationAdmissionStates = restored.admissionStates.size();

    for (const auto& item : snapshot.heldUserWideLanes)
    {
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason)
            || !ambresValidString(item.receiver, false)
            || !ambresValidSystem(item.system)
            || !ambresSatelliteUsesSystem(item.reference, item.system)
            || item.integers.empty()
            || !containerValid(item.integers.size()))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_HELD_USER_WL_HEADER_INVALID";
            }
            return false;
        }
        ZhangHeldUserWideLane held;
        held.reference = item.reference;
        for (const auto& saved : item.integers)
        {
            ZhangExactInteger value;
            if (!ambresSatelliteUsesSystem(saved.satellite, item.system)
                || saved.satellite == held.reference
                || !ambresParseExact(saved.value, value, reason)
                || !held.integers.emplace(saved.satellite, value).second)
            {
                if (reason.empty())
                {
                    reason = "AMBRES_CHECKPOINT_HELD_USER_WL_INTEGER_INVALID";
                }
                return false;
            }
        }
        if (!restored.heldUserWideLanes.emplace(
                ZhangHeldUserWideLaneKey{item.ownerRuntimeId, item.receiver,
                    item.system, item.ifAcceptance},
                std::move(held)).second)
        {
            reason = "AMBRES_CHECKPOINT_HELD_USER_WL_DUPLICATE";
            return false;
        }
        result.heldUserWideLaneIntegers += item.integers.size();
    }
    result.heldUserWideLaneStates = restored.heldUserWideLanes.size();

    for (const auto& item : snapshot.activeTransitions)
    {
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            item.transitions.size() > ZHANG_AMBRES_MAX_CONTAINER)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_ACTIVE_TRANSITION_HEADER_INVALID";
            }
            return false;
        }
        vector<ZhangPendingProductTransition> transitions;
        set<string> identities;
        // Reference-count diagnostics are computed per AMBRES branch at
        // activation time.  Validate consistency within this owner only;
        // comparing the same snapshot identity across branches is invalid.
        map<string, int> declaredSnapshotReferences;
        auto validateDeclaredSnapshotReference = [&](
            const string& identity,
            int declaredCount)
        {
            if (identity.empty())
            {
                return declaredCount == 0;
            }
            if (declaredCount <= 0)
            {
                return false;
            }
            const auto [found, inserted] =
                declaredSnapshotReferences.emplace(identity, declaredCount);
            return inserted || found->second == declaredCount;
        };
        for (const auto& pending : item.transitions)
        {
            ZhangPendingProductTransition transition;
            if (!ambresRestorePending(pending, transition, reason))
            {
                return false;
            }
            const string identity = transition.eventId + "|" +
                transition.satellite.id() + "|" +
                enum_to_string(transition.observable);
            if (!identities.insert(identity).second)
            {
                reason = "AMBRES_CHECKPOINT_ACTIVE_TRANSITION_DUPLICATE";
                return false;
            }
            if (!validateDeclaredSnapshotReference(
                    transition.oldSnapshotIdentity,
                    transition.oldSnapshotReferenceCount) ||
                !validateDeclaredSnapshotReference(
                    transition.newSnapshotIdentity,
                    transition.newSnapshotReferenceCount))
            {
                reason =
                    "AMBRES_CHECKPOINT_ACTIVE_TRANSITION_REFERENCE_"
                    "DECLARATION_CONFLICT";
                return false;
            }
            transitions.push_back(std::move(transition));
        }
        if (!restored.activeTransitions.emplace(
                item.ownerRuntimeId, std::move(transitions)).second)
        {
            reason = "AMBRES_CHECKPOINT_ACTIVE_TRANSITION_OWNER_DUPLICATE";
            return false;
        }
        result.activeTemporalTransitions += item.transitions.size();
    }

    for (const auto& item : snapshot.relinkPriors)
    {
        GTime time;
        const bool valuesValid = std::isfinite(item.firstMean) &&
            std::isfinite(item.firstVariance) && item.firstVariance >= 0 &&
            std::isfinite(item.wideLaneMean) &&
            std::isfinite(item.wideLaneVariance) && item.wideLaneVariance >= 0 &&
            std::isfinite(item.firstWideLaneCovariance) &&
            std::isfinite(item.firstVariance * item.wideLaneVariance);
        const double covarianceScale = std::max(
            1.0, item.firstVariance * item.wideLaneVariance);
        const bool covarianceValid =
            item.firstWideLaneCovariance * item.firstWideLaneCovariance <=
                item.firstVariance * item.wideLaneVariance +
                    1e-10 * covarianceScale;
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            !ambresValidSystem(item.system) ||
            !ambresSatelliteUsesSystem(item.anchor, item.system) ||
            !ambresSatelliteUsesSystem(item.satellite, item.system) ||
            item.anchor == item.satellite ||
            !ambresValidObservable(item.firstCode, item.system) ||
            !ambresValidObservable(item.secondCode, item.system) ||
            item.firstCode == item.secondCode ||
            item.productDatumVersion < 0 || !valuesValid || !covarianceValid ||
            !ambresRestoreFiniteTime(item.time, time))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_RELINK_PRIOR_INVALID";
            }
            return false;
        }
        ZhangRelinkPriorMoment moment;
        moment.time = time;
        moment.firstCode = item.firstCode;
        moment.secondCode = item.secondCode;
        moment.productDatumVersion = item.productDatumVersion;
        moment.firstMean = item.firstMean;
        moment.firstVariance = item.firstVariance;
        moment.wideLaneMean = item.wideLaneMean;
        moment.wideLaneVariance = item.wideLaneVariance;
        moment.firstWideLaneCovariance = item.firstWideLaneCovariance;
        if (!restored.relinkPriors[item.ownerRuntimeId].emplace(
                ZhangRelinkMomentKey{
                    item.system, item.anchor, item.satellite}, moment).second)
        {
            reason = "AMBRES_CHECKPOINT_RELINK_PRIOR_DUPLICATE";
            return false;
        }
        result.relinkPriorMoments++;
    }

    for (const auto& item : snapshot.scalarAccumulators)
    {
        ZhangRelinkShadowAccumulator accumulator;
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            !ambresValidString(item.key, false) || !item.initialized ||
            item.increments.empty() ||
            static_cast<int>(item.increments.size()) >
                acsConfig.zhangPppAr.multi_epoch_relink_shadow_max_epochs ||
            !ambresValidateChronology(
                item.increments, accumulator.increments, reason))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_SCALAR_ACCUMULATOR_INVALID";
            }
            return false;
        }
        double accumulatedInformation = 0;
        double accumulatedNatural = 0;
        for (const auto& increment : accumulator.increments)
        {
            accumulatedInformation += increment.information;
            accumulatedNatural += increment.natural;
        }
        if (!(accumulatedInformation > 0) ||
            !std::isfinite(accumulatedInformation) ||
            !std::isfinite(accumulatedNatural))
        {
            reason = "AMBRES_CHECKPOINT_SCALAR_ACCUMULATED_INFORMATION_INVALID";
            return false;
        }
        accumulator.integerHypothesis = item.integerHypothesis;
        accumulator.initialized = item.initialized;
        if (!restored.scalarAccumulators[item.ownerRuntimeId].emplace(
                item.key, std::move(accumulator)).second)
        {
            reason = "AMBRES_CHECKPOINT_SCALAR_ACCUMULATOR_DUPLICATE";
            return false;
        }
        result.relinkScalarIncrements += item.increments.size();
    }

    for (const auto& item : snapshot.jointAccumulators)
    {
        ZhangRelinkJointAccumulator accumulator;
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            !ambresValidString(item.key, false) || !item.initialized ||
            item.increments.empty() ||
            static_cast<int>(item.increments.size()) >
                acsConfig.zhangPppAr.multi_epoch_relink_shadow_max_epochs ||
            !ambresValidateJointChronology(
                item.increments, accumulator.increments, reason))
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_JOINT_ACCUMULATOR_INVALID";
            }
            return false;
        }
        Matrix2d accumulatedInformation = Matrix2d::Zero();
        Vector2d accumulatedNatural = Vector2d::Zero();
        for (const auto& increment : accumulator.increments)
        {
            accumulatedInformation += increment.information;
            accumulatedNatural += increment.natural;
        }
        Eigen::SelfAdjointEigenSolver<Matrix2d> accumulatedSolver(
            0.5 * (accumulatedInformation +
                   accumulatedInformation.transpose()));
        if (!accumulatedInformation.allFinite() ||
            !accumulatedNatural.allFinite() ||
            accumulatedSolver.info() != Eigen::Success ||
            accumulatedSolver.eigenvalues().minCoeff() <= 0)
        {
            reason = "AMBRES_CHECKPOINT_JOINT_ACCUMULATED_INFORMATION_INVALID";
            return false;
        }
        accumulator.wideLaneHypothesis = item.wideLaneHypothesis;
        accumulator.firstHypothesis = item.firstHypothesis;
        accumulator.initialized = item.initialized;
        if (!restored.jointAccumulators[item.ownerRuntimeId].emplace(
                item.key, std::move(accumulator)).second)
        {
            reason = "AMBRES_CHECKPOINT_JOINT_ACCUMULATOR_DUPLICATE";
            return false;
        }
        result.relinkJointIncrements += item.increments.size();
    }

    for (const auto& item : snapshot.whitenedAccumulators)
    {
        if (!ambresValidOwner(item.ownerRuntimeId, rootRuntimeId, reason) ||
            !ambresValidString(item.key, false) ||
            !ambresValidString(item.physicalSegmentIdentity) ||
            item.lastGraphEventId < 0 || item.lastProductDatumVersion < 0 ||
            item.observations.size() > ZHANG_AMBRES_MAX_CONTAINER ||
            static_cast<int>(item.observations.size()) >
                acsConfig.zhangPppAr.whitened_wl_fixed_lag_max_observations ||
            item.physicalRow.size() > ZHANG_AMBRES_MAX_CONTAINER)
        {
            if (reason.empty())
            {
                reason = "AMBRES_CHECKPOINT_WHITENED_ACCUMULATOR_HEADER_INVALID";
            }
            return false;
        }
        ZhangWhitenedWlAccumulator accumulator;
        GTime previous;
        bool havePrevious = false;
        for (const auto& observation : item.observations)
        {
            GTime time;
            if (!ambresRestoreFiniteTime(observation.time, time) ||
                !std::isfinite(observation.value) ||
                !std::isfinite(observation.variance) ||
                observation.variance <= 0 ||
                (havePrevious && !(previous < time)))
            {
                reason = "AMBRES_CHECKPOINT_WHITENED_OBSERVATION_INVALID";
                return false;
            }
            accumulator.observations.push_back({
                time, observation.value, observation.variance});
            previous = time;
            havePrevious = true;
        }
        if (accumulator.observations.size() > 1 &&
            (accumulator.observations.back().time -
             accumulator.observations.front().time).to_double() >
                acsConfig.zhangPppAr.whitened_wl_fixed_lag_seconds + 1e-3)
        {
            reason = "AMBRES_CHECKPOINT_WHITENED_WINDOW_NOT_PRUNED";
            return false;
        }
        for (const auto& arcSnapshot : item.physicalRow)
        {
            ZhangPhysicalIntegerArc arc;
            ZhangExactInteger coefficient;
            if (!ambresRestoreArc(
                    arcSnapshot, arc, coefficient, reason) ||
                !accumulator.physicalRow.emplace(arc, coefficient).second)
            {
                if (reason.empty())
                {
                    reason = "AMBRES_CHECKPOINT_WHITENED_PHYSICAL_ROW_DUPLICATE";
                }
                return false;
            }
        }
        double accumulatedInformation = 0;
        double accumulatedNatural = 0;
        for (const auto& observation : accumulator.observations)
        {
            accumulatedInformation += 1 / observation.variance;
            accumulatedNatural += observation.value / observation.variance;
        }
        if ((!accumulator.observations.empty() &&
             !(accumulatedInformation > 0)) ||
            !std::isfinite(accumulatedInformation) ||
            !std::isfinite(accumulatedNatural))
        {
            reason = "AMBRES_CHECKPOINT_WHITENED_INFORMATION_INVALID";
            return false;
        }
        accumulator.physicalSegmentIdentity = item.physicalSegmentIdentity;
        accumulator.lastGraphEventId = item.lastGraphEventId;
        accumulator.lastProductDatumVersion = item.lastProductDatumVersion;
        accumulator.initialized = item.initialized;
        if (accumulator.initialized &&
            (accumulator.physicalRow.empty() ||
             accumulator.physicalSegmentIdentity.empty()))
        {
            reason = "AMBRES_CHECKPOINT_WHITENED_INITIALIZED_ROW_EMPTY";
            return false;
        }
        if (!restored.whitenedAccumulators[item.ownerRuntimeId].emplace(
                item.key, std::move(accumulator)).second)
        {
            reason = "AMBRES_CHECKPOINT_WHITENED_ACCUMULATOR_DUPLICATE";
            return false;
        }
        result.whitenedWlObservations += item.observations.size();
    }
    return true;
}

template<class Map, class Predicate>
void ambresEraseIf(Map& values, Predicate predicate)
{
    for (auto iterator = values.begin(); iterator != values.end();)
    {
        if (predicate(*iterator))
        {
            iterator = values.erase(iterator);
        }
        else
        {
            ++iterator;
        }
    }
}

void ambresPrepareCommittedMaps(
    const string& runtimeId,
    const AmbresRestoredRuntime& restored,
    AmbresRestoredRuntime& committed)
{
    committed.e27Runtimes = zhangE27WideLaneRuntimes;
    committed.heldLattices = zhangPersistentHeldLattices;
    committed.heldEvidence = zhangPersistentHeldEvidence;
    committed.admissionStates =
        zhangProductRelationAdmissionStateRegistry();
    committed.heldUserWideLanes = zhangHeldUserWideLaneRegistry;
    committed.activeTransitions = zhangActiveTemporalBesdTransitions;
    committed.relinkPriors = zhangRelinkPriorMoments;
    committed.scalarAccumulators = zhangRelinkShadowAccumulators;
    committed.jointAccumulators = zhangRelinkJointAccumulators;
    committed.whitenedAccumulators = zhangWhitenedWlAccumulators;

    ambresEraseIf(committed.e27Runtimes, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first.runtimeId, runtimeId);
    });
    ambresEraseIf(committed.heldLattices, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first.first, runtimeId);
    });
    ambresEraseIf(committed.heldEvidence, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first.first, runtimeId);
    });
    ambresEraseIf(committed.admissionStates, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first.first, runtimeId);
    });
    ambresEraseIf(committed.heldUserWideLanes, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(
            std::get<0>(item.first), runtimeId);
    });
    ambresEraseIf(committed.activeTransitions, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first, runtimeId);
    });
    ambresEraseIf(committed.relinkPriors, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first, runtimeId);
    });
    ambresEraseIf(committed.scalarAccumulators, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first, runtimeId);
    });
    ambresEraseIf(committed.jointAccumulators, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first, runtimeId);
    });
    ambresEraseIf(committed.whitenedAccumulators, [&](const auto& item)
    {
        return zhangAmbresRuntimeBelongsTo(item.first, runtimeId);
    });

    committed.e27Runtimes.insert(
        restored.e27Runtimes.begin(), restored.e27Runtimes.end());
    committed.heldLattices.insert(
        restored.heldLattices.begin(), restored.heldLattices.end());
    committed.heldEvidence.insert(
        restored.heldEvidence.begin(), restored.heldEvidence.end());
    committed.admissionStates.insert(
        restored.admissionStates.begin(), restored.admissionStates.end());
    committed.heldUserWideLanes.insert(
        restored.heldUserWideLanes.begin(),
        restored.heldUserWideLanes.end());
    committed.activeTransitions.insert(
        restored.activeTransitions.begin(), restored.activeTransitions.end());
    committed.relinkPriors.insert(
        restored.relinkPriors.begin(), restored.relinkPriors.end());
    committed.scalarAccumulators.insert(
        restored.scalarAccumulators.begin(),
        restored.scalarAccumulators.end());
    committed.jointAccumulators.insert(
        restored.jointAccumulators.begin(), restored.jointAccumulators.end());
    committed.whitenedAccumulators.insert(
        restored.whitenedAccumulators.begin(),
        restored.whitenedAccumulators.end());
}

} // namespace
