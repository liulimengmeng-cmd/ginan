// Generated schema declarations only; not a reference gate.
#pragma once
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include "common/zhangFactorCapture.hpp"
#include "common/zhangPhaseContinuity.hpp"
#include "common/zhangHybridService.hpp"
#include "common/zhangProductNumericalContinuity.hpp"
using std::map; using std::set; using std::string; using std::vector;
namespace boost::serialization
{
template<class Archive>
void serialize(Archive& ar, ZhangCapturedStateKey& value, const unsigned int)
{
	ar & value.type;
	ar & value.satellite;
	ar & value.receiver;
	ar & value.number;
}

template<class Archive>
void serialize(Archive& ar, SparseMatrix<double>& value, const unsigned int)
{
	int rows = value.rows();
	int columns = value.cols();
	std::vector<int> rowIndices;
	std::vector<int> columnIndices;
	std::vector<double> coefficients;
	if constexpr (Archive::is_saving::value)
	{
		rowIndices.reserve(value.nonZeros());
		columnIndices.reserve(value.nonZeros());
		coefficients.reserve(value.nonZeros());
		for (int outer = 0; outer < value.outerSize(); outer++)
		for (SparseMatrix<double>::InnerIterator entry(value, outer);
			 entry; ++entry)
		{
			rowIndices.push_back(entry.row());
			columnIndices.push_back(entry.col());
			coefficients.push_back(entry.value());
		}
	}
	ar & rows;
	ar & columns;
	ar & rowIndices;
	ar & columnIndices;
	ar & coefficients;
	if constexpr (Archive::is_loading::value)
	{
		if (rows < 0 || columns < 0 || rows > 1000000 || columns > 1000000
		 || rowIndices.size() != columnIndices.size()
		 || rowIndices.size() != coefficients.size()
		 || rowIndices.size() > 100000000)
		{
			throw std::runtime_error("INVALID_ZHANG_SPARSE_MATRIX_ARCHIVE");
		}
		std::vector<Triplet<double>> triplets;
		triplets.reserve(coefficients.size());
		for (std::size_t index = 0; index < coefficients.size(); index++)
		{
			if (rowIndices[index] < 0 || rowIndices[index] >= rows
			 || columnIndices[index] < 0 || columnIndices[index] >= columns
			 || !std::isfinite(coefficients[index]))
			{
				throw std::runtime_error("INVALID_ZHANG_SPARSE_MATRIX_ENTRY");
			}
			triplets.emplace_back(
				rowIndices[index], columnIndices[index], coefficients[index]);
		}
		value.resize(rows, columns);
		value.setFromTriplets(triplets.begin(), triplets.end());
	}
}

template<class Archive>
void serialize(
	Archive& ar, ZhangCapturedPhysicalArcVersion& value, const unsigned int)
{
	ar & value.arc;
	ar & value.version;
}

template<class Archive>
void serialize(Archive& ar, ZhangInnovationScaleGroup& value,
	const unsigned int)
{
	ar & value.identity;
	ar & value.blocks;
	ar & value.samples;
	ar & value.marginalStandardisedSquaredSum;
	ar & value.maximumAbsoluteRatio;
}

template<class Archive>
void serialize(Archive& ar, ZhangCapturedSnapshotOperation& value,
	const unsigned int)
{
	ar & value.kind;
	ar & value.afterEventSequence;
	ar & value.operationSequence;
	ar & value.identities;
	ar & value.physicalVersions;
	ar & value.rows;
	ar & value.offsets;
}

template<class Archive>
void serialize(Archive& ar, ZhangCanonicalSatelliteRelation& value,
	const unsigned int)
{
	ar & value.anchor;
	ar & value.satellite;
}

template<class Archive>
void serialize(Archive& ar, ZhangPersistentProductDatumCheckpointState& value,
	const unsigned int)
{
	ar & value.system;
	ar & value.observable;
	ar & value.relation;
	ar & value.initialised;
	ar & value.version;
	ar & value.anchorPhaseSegment;
	ar & value.satellitePhaseSegment;
	ar & value.anchorDatumVersion;
	ar & value.satelliteDatumVersion;
}

template<class Archive>
void serialize(Archive& ar, ZhangPersistentProductDatumCheckpoint& value,
	const unsigned int)
{
	ar & value.canonicalRelations;
	ar & value.datumStates;
}

template<class Archive>
void serialize(Archive& ar, ZhangSatellitePhaseSegment& value,
	const unsigned int)
{
	ar & value.satellite;
	ar & value.segment;
}

template<class Archive>
void serialize(Archive& ar, ZhangSatelliteDatumRelation& value,
	const unsigned int)
{
	ar & value.a;
	ar & value.b;
	ar & value.difference;
	ar & value.promoted;
	ar & value.provenance;
}

template<class Archive>
void serialize(Archive& ar, ZhangSatelliteDatumManagerCheckpoint& value,
	const unsigned int)
{
	ar & value.system;
	ar & value.observable;
	ar & value.currentSegments;
	ar & value.discontinuityCounters;
	ar & value.datumVersions;
	ar & value.alignmentCycles;
	ar & value.alignmentKnown;
	ar & value.precisionValid;
	ar & value.relations;
	ar & value.redundantRelations;
	ar & value.eventCounts;
	ar & value.conflictCount;
	ar & value.topologyVersion;
	ar & value.alignmentGeneration;
}

template<class Archive>
void serialize(Archive& ar, ZhangHybridRealGaugeCheckpoint& value,
	const unsigned int)
{
	ar & value.initialized;
	ar & value.generation;
	ar & value.previousValues;
	ar & value.previousSegments;
}

template<class Archive>
void serialize(Archive& ar, ZhangProductDatumVersionTracker& value,
	const unsigned int)
{
	ar & value.initialized;
	ar & value.version;
}

} // namespace boost::serialization
namespace checkpoint_test_wire {
struct ProductKey
{
    SatSys    satellite;
    E_ObsCode observable = E_ObsCode::NONE;

    bool operator<(const ProductKey& other) const
    {
        return std::tie(satellite, observable) <
               std::tie(other.satellite, other.observable);
    }

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & satellite;
		ar & observable;
	}
};

struct E27RawNoiseRow
{
    vector<string> stampedKeys;
    VectorXd coefficients;
    VectorXd variances;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & stampedKeys;
		ar & coefficients;
		ar & variances;
	}
};

struct PromotionEvidenceKey
{
    E_Sys      system = E_Sys::NONE;
    E_ObsCode  observable = E_ObsCode::NONE;
    SatSys     a;
    int        segmentA = 0;
    SatSys     b;
    int        segmentB = 0;

    bool operator<(const PromotionEvidenceKey& other) const
    {
        return std::tie(system, observable, a, segmentA, b, segmentB) <
               std::tie(
                   other.system, other.observable,
                   other.a, other.segmentA, other.b, other.segmentB
               );
    }

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & system;
		ar & observable;
		ar & a;
		ar & segmentA;
		ar & b;
		ar & segmentB;
	}
};

struct PromotionEvidence
{
    long long difference = 0;
    long int  lastEpoch = 0;
    int       confirmations = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & difference;
		ar & lastEpoch;
		ar & confirmations;
	}
};

struct ProductHistoryKey
{
    string     solution;
    SatSys     satellite;
    E_ObsCode  observable = E_ObsCode::NONE;

    bool operator<(const ProductHistoryKey& other) const
    {
        return std::tie(solution, satellite, observable) <
               std::tie(other.solution, other.satellite, other.observable);
    }

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & solution;
		ar & satellite;
		ar & observable;
	}
};

struct UserReferenceKey
{
    string         runtimeId;
    string         receiver;
    E_Sys          sys = E_Sys::NONE;
    E_ObsCode      observable = E_ObsCode::NONE;

    bool operator<(const UserReferenceKey& other) const
    {
        return std::tie(runtimeId, receiver, sys, observable) <
               std::tie(
                   other.runtimeId, other.receiver,
                   other.sys, other.observable);
    }

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & runtimeId;
		ar & receiver;
		ar & sys;
		ar & observable;
	}
};

struct UserReferenceState
{
    SatSys reference;
    int    productCounter = -1;
    int    datumVersion = -1;
    map<SatSys, std::pair<int, int>> satelliteDatum;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & reference;
		ar & productCounter;
		ar & datumVersion;
		ar & satelliteDatum;
	}
};

struct UserDualReferenceKey
{
    string runtimeId;
    string receiver;
    E_Sys system = E_Sys::NONE;

    bool operator<(const UserDualReferenceKey& other) const
    {
        return std::tie(runtimeId, receiver, system) <
            std::tie(other.runtimeId, other.receiver, other.system);
    }

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & runtimeId;
		ar & receiver;
		ar & system;
	}
};

struct ZhangCapturedFactorEventCheckpoint
{
	ZhangCapturedFactorKind kind = ZhangCapturedFactorKind::MEASUREMENT;
	ZhangCheckpointTime time;
	std::size_t sequence = 0;
	std::uint64_t beforeCommitSequence = 0;
	std::uint64_t afterCommitSequence = 0;
	string factorId;
	string label;
	vector<ZhangCapturedStateKey> sourceKeys;
	vector<ZhangCapturedStateKey> destinationKeys;
	SparseMatrix<double> design;
	SparseMatrix<double> covariance;
	VectorXd rightHandSide;
	VectorXd prefitRatios;
	vector<ZhangCapturedStateKey> observationKeys;
	bool dimensionPreserving = false;
	bool nonsingularCoordinateTransform = false;
	bool preserveUnrepresentablePersistentTargets = false;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & kind;
		ar & time;
		ar & sequence;
		ar & beforeCommitSequence;
		ar & afterCommitSequence;
		ar & factorId;
		ar & label;
		ar & sourceKeys;
		ar & destinationKeys;
		ar & design;
		ar & covariance;
		ar & rightHandSide;
		ar & prefitRatios;
		ar & observationKeys;
		ar & dimensionPreserving;
		ar & nonsingularCoordinateTransform;
		ar & preserveUnrepresentablePersistentTargets;
	}
};

struct ZhangCapturedPhysicalTargetCheckpoint
{
	ZhangCheckpointTime time;
	std::size_t afterEventSequence = 0;
	string identity;
	string physicalArcSignature;
	string phaseSegmentIdentity;
	vector<ZhangCapturedPhysicalArcVersion> physicalArcVersions;
	bool resetPhysicalIdentity = false;
	bool continuedAcrossCoordinateChange = false;
	vector<ZhangCapturedStateKey> keys;
	SparseMatrix<double> row;
	double offset = 0;
	double mean = 0;
	double variance = 0;
	int unresolvedIntegerGaugeRank = 0;
	string integerGaugeIdentity;
	string separatorIdentity;
	string canonicalCoordinateIdentity;
	string productDatumIdentity;
	int productDatumVersion = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & time;
		ar & afterEventSequence;
		ar & identity;
		ar & physicalArcSignature;
		ar & phaseSegmentIdentity;
		ar & physicalArcVersions;
		ar & resetPhysicalIdentity;
		ar & continuedAcrossCoordinateChange;
		ar & keys;
		ar & row;
		ar & offset;
		ar & mean;
		ar & variance;
		ar & unresolvedIntegerGaugeRank;
		ar & integerGaugeIdentity;
		ar & separatorIdentity;
		ar & canonicalCoordinateIdentity;
		ar & productDatumIdentity;
		ar & productDatumVersion;
	}
};

struct ZhangCapturedUnresolvedIntegerDatumCheckpoint
{
	ZhangCheckpointTime time;
	std::size_t afterEventSequence = 0;
	string identity;
	int missingGaugeRank = 1;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & time;
		ar & afterEventSequence;
		ar & identity;
		ar & missingGaugeRank;
	}
};

struct ZhangCapturedRetainedTargetBlockCheckpoint
{
	ZhangCheckpointTime time;
	std::size_t afterEventSequence = 0;
	int targetCount = 0;
	int informationRank = 0;
	int residualDof = 0;
	int projectedGaugeRank = 0;
	bool likelihoodValid = false;
	bool valid = false;
	double whitenedSquaredNorm = std::numeric_limits<double>::quiet_NaN();
	VectorXd whitenedResidual;
	vector<string> separatorIdentities;
	vector<string> gaugeIdentities;
	vector<bool> absoluteValid;
	vector<double> coordinateOffsets;
	MatrixXd likelihoodDesign;
	VectorXd likelihoodObservation;
	MatrixXd likelihoodCovariance;
	string failureReason;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & time;
		ar & afterEventSequence;
		ar & targetCount;
		ar & informationRank;
		ar & residualDof;
		ar & projectedGaugeRank;
		ar & likelihoodValid;
		ar & valid;
		ar & whitenedSquaredNorm;
		ar & whitenedResidual;
		ar & separatorIdentities;
		ar & gaugeIdentities;
		ar & absoluteValid;
		ar & coordinateOffsets;
		ar & likelihoodDesign;
		ar & likelihoodObservation;
		ar & likelihoodCovariance;
		ar & failureReason;
	}
};

struct ZhangFactorCaptureCheckpoint
{
	std::size_t maximumEvents = 0;
	vector<ZhangCapturedStateKey> initialKeys;
	VectorXd initialMean;
	MatrixXd initialCovariance;
	vector<ZhangCapturedStateKey> currentKeys;
	VectorXd replayMean;
	MatrixXd replayCovariance;
	std::deque<ZhangCapturedFactorEventCheckpoint> events;
	std::deque<ZhangCapturedSnapshotOperation> snapshotOperations;
	std::deque<ZhangCapturedPhysicalTargetCheckpoint> physicalTargets;
	std::deque<ZhangCapturedUnresolvedIntegerDatumCheckpoint>
		unresolvedIntegerDatums;
	std::deque<ZhangCapturedRetainedTargetBlockCheckpoint>
		retainedTargetBlocks;
	ZhangCapturedRetainedTargetBlockCheckpoint currentRetainedTargetBlock;
	vector<ZhangInnovationScaleGroup> innovationScaleGroups;
	VectorXd lastMeasurementPriorMean;
	MatrixXd lastMeasurementPriorCovariance;
	ZhangCheckpointTime lastMeasurementTime;
	std::size_t lastMeasurementTargetStart = 0;
	string lastFailure;
	string lastTargetDispositionReason;
	double maximumReplayPriorMeanRelativeError = 0;
	double maximumReplayPriorCovarianceRelativeError = 0;
	double maximumCommitMeanRelativeError = 0;
	double maximumCommitCovarianceRelativeError = 0;
	std::uint64_t lastCommitSequence = 0;
	bool commitSequenceAnchored = false;
	double maximumTargetMeanRelativeError = 0;
	double maximumTargetVarianceRelativeError = 0;
	double maximumRawSquareRootMeanRelativeError = 0;
	double maximumRawSquareRootCovarianceRelativeError = 0;
	double maximumPersistentTransformMeanRelativeError = 0;
	double maximumPersistentTransformCovarianceRelativeError = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & maximumEvents;
		ar & initialKeys;
		ar & initialMean;
		ar & initialCovariance;
		ar & currentKeys;
		ar & replayMean;
		ar & replayCovariance;
		ar & events;
		ar & snapshotOperations;
		ar & physicalTargets;
		ar & unresolvedIntegerDatums;
		ar & retainedTargetBlocks;
		ar & currentRetainedTargetBlock;
		ar & innovationScaleGroups;
		ar & lastMeasurementPriorMean;
		ar & lastMeasurementPriorCovariance;
		ar & lastMeasurementTime;
		ar & lastMeasurementTargetStart;
		ar & lastFailure;
		ar & lastTargetDispositionReason;
		ar & maximumReplayPriorMeanRelativeError;
		ar & maximumReplayPriorCovarianceRelativeError;
		ar & maximumCommitMeanRelativeError;
		ar & maximumCommitCovarianceRelativeError;
		ar & lastCommitSequence;
		ar & commitSequenceAnchored;
		ar & maximumTargetMeanRelativeError;
		ar & maximumTargetVarianceRelativeError;
		ar & maximumRawSquareRootMeanRelativeError;
		ar & maximumRawSquareRootCovarianceRelativeError;
		ar & maximumPersistentTransformMeanRelativeError;
		ar & maximumPersistentTransformCovarianceRelativeError;
	}
};

struct ZhangPhaseContinuityCheckpoint
{
	int counter = 0;
	long long integerShiftCycles = 0;
	double fractionalShiftCycles = 0;
	int datumVersion = 0;
	ZhangCheckpointTime validFrom;
	int iod = 0;
	string resetReason;
	int stabilizationRemaining = 0;
	bool hasFixedDatum = false;
	ZhangCheckpointTime lastEpoch;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & counter;
		ar & integerShiftCycles;
		ar & fractionalShiftCycles;
		ar & datumVersion;
		ar & validFrom;
		ar & iod;
		ar & resetReason;
		ar & stabilizationRemaining;
		ar & hasFixedDatum;
		ar & lastEpoch;
	}
};

struct GlobalContinuityCheckpoint
{
	int counter = 0;
	int datumVersion = 0;
	ZhangCheckpointTime validFrom;
	int iod = 0;
	string resetReason;
	int stabilizationRemaining = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & counter;
		ar & datumVersion;
		ar & validFrom;
		ar & iod;
		ar & resetReason;
		ar & stabilizationRemaining;
	}
};

struct ZhangPendingSnapshotPinsCheckpoint
{
	ZhangCheckpointTime eventTime;
	set<string> identities;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & eventTime;
		ar & identities;
	}
};

struct E27NoiseSensitivityCheckpoint
{
	ZhangCheckpointTime time;
	VectorXd stateDerivative;
	double variance = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & time;
		ar & stateDerivative;
		ar & variance;
	}
};

struct E27JointNoiseRuntimeCheckpoint
{
	map<string, E27NoiseSensitivityCheckpoint> sensitivities;
	string failureReason;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & sensitivities;
		ar & failureReason;
	}
};

struct ProductHistoryCheckpoint
{
	ZhangCheckpointTime time;
	double correction = 0;
	int discontinuityCounter = 0;
	int datumVersion = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & time;
		ar & correction;
		ar & discontinuityCounter;
		ar & datumVersion;
	}
};

struct ProductPublicationIdentityCheckpoint
{
	string phaseSegmentId;
	string productFunctional;
	int backendGeneration = -1;
	double correctionSigma = 0;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & phaseSegmentId;
		ar & productFunctional;
		ar & backendGeneration;
		ar & correctionSigma;
	}
};

struct ProductNumericalHistoryCheckpoint
{
    ZhangProductNumericalDynamics dynamics;
    string publicationProof;
    template<class Archive>
    void serialize(Archive& archive, const unsigned int)
    { archive & dynamics; archive & publicationProof; }
};

struct ZhangGraphEdgeCheckpoint
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

struct ZhangProductIntegerFunctionalCheckpoint
{
	SatSys satellite;
	SatSys referenceSatellite;
	vector<ZhangGraphEdgeCheckpoint> physicalEdges;
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

struct ZhangProductIntegerTransitionCheckpoint
{
	vector<ZhangGraphEdgeCheckpoint> physicalEdges;
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

struct ZhangFullJointStructuralChainCheckpoint
{
	string publishedNativeFingerprint, publishedSegment, tailNativeFingerprint;
	int publishedBackend = -1, tailBackend = -1;
	std::uint64_t linkCount = 0;
	ZhangProductIntegerTransitionCheckpoint nativeToSegment;
	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & publishedNativeFingerprint & publishedSegment & tailNativeFingerprint;
		ar & publishedBackend & tailBackend & linkCount & nativeToSegment;
	}
};

struct ZhangPendingProductTransitionCheckpoint
{
	ZhangCheckpointTime eventTime;
	E_Sys system = E_Sys::NONE;
	SatSys satellite;
	E_ObsCode observable = E_ObsCode::NONE;
	string eventId;
	ZhangProductIntegerFunctionalCheckpoint oldFunctional;
	ZhangProductIntegerFunctionalCheckpoint newFunctional;
	ZhangProductIntegerTransitionCheckpoint transition;
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

struct E27RawNoiseRowCheckpoint
{
	string epoch;
	string receiver;
	E_Sys system = E_Sys::NONE;
	SatSys satellite;
	E27RawNoiseRow row;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & epoch;
		ar & receiver;
		ar & system;
		ar & satellite;
		ar & row;
	}
};

struct ZhangHybridRealGaugeCheckpointEntry
{
	string solution;
	E_Sys system = E_Sys::NONE;
	string parameter;
	E_ObsCode observable = E_ObsCode::NONE;
	ZhangHybridRealGaugeCheckpoint gauge;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & solution;
		ar & system;
		ar & parameter;
		ar & observable;
		ar & gauge;
	}
};

struct ZhangPppArCheckpointEnvelope
{
	std::uint32_t schemaVersion = ZHANG_PPP_AR_CHECKPOINT_SCHEMA_VERSION;
	string runtimeId;
	vector<ZhangPendingProductTransitionCheckpoint> pendingTransitions;
	bool hasPendingSnapshotPins = false;
	ZhangPendingSnapshotPinsCheckpoint pendingPins;
	bool factorCaptureConfigured = false;
	bool hasFactorCapture = false;
	ZhangFactorCaptureCheckpoint factorCapture;
	bool hasPersistentProductDatumRegistry = false;
	ZhangPersistentProductDatumCheckpoint persistentProductDatumRegistry;
	vector<E27RawNoiseRowCheckpoint> e27RawRows;
	bool hasE27JointNoiseRuntime = false;
	E27JointNoiseRuntimeCheckpoint e27JointNoiseRuntime;

	map<ProductKey, ZhangPhaseContinuityCheckpoint> continuity;
	map<E_Sys, ZhangProductDatumVersionTracker> productDatumVersionTrackers;
	map<ProductKey, string> physicalFunctionalIdentities;
	map<ProductKey, ZhangProductIntegerFunctionalCheckpoint>
		physicalFunctionals;
	map<ProductKey, string> sBasisFingerprints;
	map<ProductKey, string> phaseSegmentIdentities;
	map<ProductKey, long long> treeAlignmentCycles;
	map<ProductKey, string> snapshotIdentities;
	map<std::pair<E_Sys, E_ObsCode>, GlobalContinuityCheckpoint>
		globalContinuity;
	vector<ZhangSatelliteDatumManagerCheckpoint> satelliteDatumManagers;
	vector<ZhangHybridRealGaugeCheckpointEntry> hybridRealGauges;
	map<PromotionEvidenceKey, PromotionEvidence> promotionEvidence;
	map<PromotionEvidenceKey, PromotionEvidence> relinkEvidence;
	map<ProductHistoryKey, ProductHistoryCheckpoint> productHistory;
	map<UserReferenceKey, UserReferenceState> userReferences;
	map<UserDualReferenceKey, SatSys> userDualReferences;
	map<ProductHistoryKey, ProductPublicationIdentityCheckpoint>
		productPublicationIdentity;
	map<ProductKey, ZhangProductIntegerFunctionalCheckpoint>
		publishedPhysicalFunctionals;
	map<ProductKey, ZhangProductIntegerFunctionalCheckpoint> segmentPhysicalFunctionals;
	map<ProductKey, ZhangProductIntegerTransitionCheckpoint> fullJointDifferences;
	map<ProductKey, ZhangFullJointStructuralChainCheckpoint> structuralTransportChains;
    map<ProductHistoryKey, ProductNumericalHistoryCheckpoint> productNumericalHistory;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int)
	{
		ar & schemaVersion;
		ar & runtimeId;
		ar & pendingTransitions;
		ar & hasPendingSnapshotPins;
		ar & pendingPins;
		ar & factorCaptureConfigured;
		ar & hasFactorCapture;
		ar & factorCapture;
		ar & hasPersistentProductDatumRegistry;
		ar & persistentProductDatumRegistry;
		ar & e27RawRows;
		ar & hasE27JointNoiseRuntime;
		ar & e27JointNoiseRuntime;
		ar & continuity;
		ar & productDatumVersionTrackers;
		ar & physicalFunctionalIdentities;
		ar & physicalFunctionals;
		ar & sBasisFingerprints;
		ar & phaseSegmentIdentities;
		ar & treeAlignmentCycles;
		ar & snapshotIdentities;
		ar & globalContinuity;
		ar & satelliteDatumManagers;
		ar & hybridRealGauges;
		ar & promotionEvidence;
		ar & relinkEvidence;
		ar & productHistory;
		ar & userReferences;
		ar & userDualReferences;
		if (schemaVersion >= 4)
		{
			ar & productPublicationIdentity;
			ar & publishedPhysicalFunctionals;
		}
		if (schemaVersion >= 5)
		{
			ar & segmentPhysicalFunctionals;
			ar & fullJointDifferences;
		}
		if (schemaVersion >= 6) ar & structuralTransportChains;
        if (schemaVersion >= 7) ar & productNumericalHistory;
	}
};
}
