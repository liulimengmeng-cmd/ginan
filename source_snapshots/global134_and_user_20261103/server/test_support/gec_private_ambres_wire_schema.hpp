#pragma once
namespace gec_test_ambres_wire {
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
}
