#pragma once
// Test-only independent schema mirror, frozen from zhangReceiverCheckpoint.cpp
// through ZhangReceiverRuntimeEnvelope. No live capture/import implementation.
#include "pea/zhangReceiverCheckpoint.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include "common/acsConfig.hpp"
#include "common/navigation.hpp"
#include "common/receiver.hpp"
#include "common/sinex.hpp"
#include "common/zhangCheckpoint.hpp"

namespace ReceiverAuditSchema
{
using std::array;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

constexpr std::size_t MAX_RECEIVER_RUNTIME_PAYLOAD_BYTES =
	std::size_t{2} * 1024 * 1024 * 1024;
constexpr const char* POST_EPOCH_BOUNDARY = "POST_EPOCH_COMMITTED";
constexpr const char* OBSERVATION_POLICY = "RESET_POST_EPOCH_WORKSET";
constexpr const char* CACHE_POLICY = "RESET_DERIVED_CACHES";
constexpr const char* READY_POLICY = "RESET_NOT_READY";
constexpr const char* CONFIG_POINTER_POLICY =
	"REBIND_FROM_CONFIGURED_DESTINATION";
constexpr const char* SATELLITE_SCOPE = "PRECISE_PRODUCT_RUNTIME_ONLY";

struct ZhangRuntimeAttitude
{
	ZhangCheckpointTime startTime;
	double startSign = 0;
	double startYaw = 0;
	double startYawRate = 0;
	ZhangCheckpointTime excludeTime;
	double nominalYaw = 0;
	double modelYaw = 0;
	ZhangCheckpointTime modelYawTime;
	bool modelYawValid = false;
	VectorEcef eXBody;
	VectorEcef eYBody;
	VectorEcef eZBody;
	VectorEcef eXAnt;
	VectorEcef eYAnt;
	VectorEcef eZAnt;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & startTime;
		ar & startSign;
		ar & startYaw;
		ar & startYawRate;
		ar & excludeTime;
		ar & nominalYaw;
		ar & modelYaw;
		ar & modelYawTime;
		ar & modelYawValid;
		ar & eXBody;
		ar & eYBody;
		ar & eZBody;
		ar & eXAnt;
		ar & eYAnt;
		ar & eZAnt;
	}
};

struct ZhangRuntimeLinearCombination
{
	bool valid = false;
	double gfPhaseMetres = 0;
	double gfCodeMetres = 0;
	double wlPhaseMetres = 0;
	double wlPhaseCycles = 0;
	double wlCodeMetres = 0;
	double wlCodeCycles = 0;
	double nlPhaseMetres = 0;
	double nlPhaseCycles = 0;
	double nlCodeMetres = 0;
	double nlCodeCycles = 0;
	double ifPhaseMetres = 0;
	double ifCodeMetres = 0;
	double mwMetres = 0;
	double mwCycles = 0;
	double wavelengthA = 0;
	double wavelengthB = 0;
	double wavelengthWl = 0;
	double wavelengthNl = 0;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & valid;
		ar & gfPhaseMetres;
		ar & gfCodeMetres;
		ar & wlPhaseMetres;
		ar & wlPhaseCycles;
		ar & wlCodeMetres;
		ar & wlCodeCycles;
		ar & nlPhaseMetres;
		ar & nlPhaseCycles;
		ar & nlCodeMetres;
		ar & nlCodeCycles;
		ar & ifPhaseMetres;
		ar & ifCodeMetres;
		ar & mwMetres;
		ar & mwCycles;
		ar & wavelengthA;
		ar & wavelengthB;
		ar & wavelengthWl;
		ar & wavelengthNl;
	}
};

struct ZhangRuntimeLcState
{
	ZhangCheckpointTime time;
	SatSys satellite;
	map<int, double> phaseMetres;
	map<int, double> codeMetres;
	map<int, double> multipathMetres;
	map<pair<int, int>, ZhangRuntimeLinearCombination> combinations;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & time;
		ar & satellite;
		ar & phaseMetres;
		ar & codeMetres;
		ar & multipathMetres;
		ar & combinations;
	}
};

struct ZhangRuntimeSlipFilter
{
	array<double, 3> state{};
	array<array<double, 3>, 3> covariance{};
	int slip = 0;
	array<int, 3> ambiguity{};
	int epochCount = 0;
	ZhangRuntimeLcState previousCombination;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & state;
		ar & covariance;
		ar & slip;
		ar & ambiguity;
		ar & epochCount;
		ar & previousCombination;
	}
};

struct ZhangRuntimeSignalStatus
{
	unsigned int savedSlipFlags = 0;
	unsigned int slipFlags = 0;
	unsigned int phaseRejectCount = 0;
    ZhangSignalTracking tracking;
    bool phaseQuarantined=false, codeQuarantined=false, previousPhaseQuarantined=false;
    bool phaseContinuityUnproven=false;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & savedSlipFlags;
		ar & slipFlags;
		ar & phaseRejectCount;
        ar & tracking & phaseQuarantined & codeQuarantined & previousPhaseQuarantined;
        if(version>=1)ar & phaseContinuityUnproven;
	}
};

struct ZhangRuntimeSatelliteStatus
{
	double ambiguityVariance = 0;
	double geometryFreeAmbiguity = 0;
	ZhangCheckpointTime lastObservationTime;
	double externalIonosphere = 0;
	double externalIonosphereVariance = 0;

	double mwSlipMean = 0;
	double mwSlipVariance = 0;
	double emwSlipMean = 0;
	double emwSlipVariance = 0;
	array<int, 3> repairedAmbiguity{};
	double melbourneWubbena = 0;
	double geometryFree = 0;
    ZhangCombinationDetector acceptedGf, acceptedMw;
    ZhangSlipEvidenceState slipEvidence;
	ZhangRuntimeSlipFilter slipFilter;
	ZhangRuntimeLcState previousCombination;
	ZhangRuntimeLcState currentCombination;

	double azimuth = 0;
	double elevation = 0;
	double phaseWindup = 0;
	double wetMapping = 0;
	array<double, 2> wetGradientMapping{};
	VectorEcef lineOfSight;
	ZhangCheckpointTime lastIonosphereTime;
	double deltaIonosphere = 0;
	double sigmaIonosphere = 0;
	double previousStec = 0;
	double nadir = 0;
	bool slip = false;
	map<string, ZhangRuntimeSignalStatus> signalStatuses;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & ambiguityVariance;
		ar & geometryFreeAmbiguity;
		ar & lastObservationTime;
		ar & externalIonosphere;
		ar & externalIonosphereVariance;
		ar & mwSlipMean;
		ar & mwSlipVariance;
		ar & emwSlipMean;
		ar & emwSlipVariance;
		ar & repairedAmbiguity;
		ar & melbourneWubbena;
		ar & geometryFree;
        ar & acceptedGf & acceptedMw;
		ar & slipFilter;
		ar & previousCombination;
		ar & currentCombination;
		ar & azimuth;
		ar & elevation;
		ar & phaseWindup;
		ar & wetMapping;
		ar & wetGradientMapping;
		ar & lineOfSight;
		ar & lastIonosphereTime;
		ar & deltaIonosphere;
		ar & sigmaIonosphere;
		ar & previousStec;
		ar & nadir;
		ar & slip;
		ar & signalStatuses;
        if(version>=1) ar & slipEvidence;
	}
};

} // namespace (class version 0 remains readable from the old FLOAT checkpoint)
BOOST_CLASS_VERSION(ReceiverAuditSchema::ZhangRuntimeSatelliteStatus, 1)
BOOST_CLASS_VERSION(ReceiverAuditSchema::ZhangRuntimeSignalStatus, 1)
namespace ReceiverAuditSchema
{

template <typename TYPE>
struct ZhangRuntimeMetadataField
{
	TYPE value{};
	bool valid = false;
	int winningSource = 0;
	std::uint32_t sourceMask = 0;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & value;
		ar & valid;
		ar & winningSource;
		ar & sourceMask;
	}
};

struct ZhangRuntimeReceiverMetadata
{
	vector<int> sourcePriority;
	ZhangRuntimeMetadataField<string> receiverType;
	ZhangRuntimeMetadataField<string> receiverFirmware;
	ZhangRuntimeMetadataField<string> receiverSerial;
	ZhangRuntimeMetadataField<string> antennaDescriptor;
	ZhangRuntimeMetadataField<string> antennaSerial;
	ZhangRuntimeMetadataField<string> markerName;
	ZhangRuntimeMetadataField<string> markerNumber;
	ZhangRuntimeMetadataField<Vector3d> antennaDelta;
	ZhangRuntimeMetadataField<Vector3d> stationPosition;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & sourcePriority;
		ar & receiverType;
		ar & receiverFirmware;
		ar & receiverSerial;
		ar & antennaDescriptor;
		ar & antennaSerial;
		ar & markerName;
		ar & markerNumber;
		ar & antennaDelta;
		ar & stationPosition;
	}
};

struct ZhangRuntimeKfBindingShape
{
	map<ZhangCheckpointKfKey, bool> receiverBindings;
	map<string, bool> filterChunkTraceBindings;
	std::size_t stateRejectCallbackCount = 0;
	std::size_t measurementRejectCallbackCount = 0;
	bool acceptedMeasurementCallback = false;
	bool stateTransitionCallback = false;
	bool exactStateTransformCallback = false;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & receiverBindings;
		ar & filterChunkTraceBindings;
		ar & stateRejectCallbackCount;
		ar & measurementRejectCallbackCount;
		ar & acceptedMeasurementCallback;
		ar & stateTransitionCallback;
		ar & exactStateTransformCallback;
	}
};

struct ZhangRuntimeSolution
{
	ZhangCheckpointTime sppTime;
	ZhangCheckpointKfCore sppState;
	ZhangRuntimeKfBindingShape sppBindingShape;
	VectorEcef sppPosition;
	double sppClock = 0;
	int clockReferenceSystem = 0;
	double sppPppClockOffset = 0;
	bool clockAdjustmentReady = false;
	int status = 0;
	int measurementCount = 0;
	double gdop = 0;
	double pdop = 0;
	double hdop = 0;
	double vdop = 0;
	double horizontalProtectionLevel = -1;
	double verticalProtectionLevel = -1;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & sppTime;
		ar & sppState;
		ar & sppBindingShape;
		ar & sppPosition;
		ar & sppClock;
		ar & clockReferenceSystem;
		ar & sppPppClockOffset;
		ar & clockAdjustmentReady;
		ar & status;
		ar & measurementCount;
		ar & gdop;
		ar & pdop;
		ar & hdop;
		ar & vdop;
		ar & horizontalProtectionLevel;
		ar & verticalProtectionLevel;
	}
};

struct ZhangRuntimeSinexState
{
	string siteBindingFingerprint;
	string receiverBindingFingerprint;
	string antennaBindingFingerprint;
	string eccentricityBindingFingerprint;
	array<double, 3> start{};
	array<double, 3> stop{};
	bool primary = false;
	VectorEcef position;
	VectorEcef variance;
	VectorEcef velocity;
	ZhangCheckpointTime referenceEpoch;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & siteBindingFingerprint;
		ar & receiverBindingFingerprint;
		ar & antennaBindingFingerprint;
		ar & eccentricityBindingFingerprint;
		ar & start;
		ar & stop;
		ar & primary;
		ar & position;
		ar & variance;
		ar & velocity;
		ar & referenceEpoch;
	}
};

struct ZhangRuntimeReceiver
{
	string mapKey;
	string id;
	bool pseudoReceiver = false;
	bool invalid = false;
	ZhangRuntimeSinexState sinex;
	ZhangRuntimeReceiverMetadata metadata;
	map<string, string> metadataMap;
	string source;

	ZhangCheckpointTime firstEpoch;
	ZhangCheckpointTime lastEpoch;
	int epochCount = 0;
	int observationCount = 0;
	int slipCount = 0;
	map<int, int> codeCount;
	map<string, int> satelliteCount;
	int receiverErrorEpochs = 0;
	int receiverErrorCount = 0;

	ZhangRuntimeSolution solution;
	string antennaType;
	string receiverType;
	string antennaId;
	map<SatSys, ZhangRuntimeSatelliteStatus> satelliteStatuses;
	VectorEnu antennaDelta;
	ZhangRuntimeAttitude attitude;

	bool primaryApriori = false;
	array<double, 3> aprioriTime{};
	double aprioriClock = 0;
	double aprioriClockVariance = 0;
	Vector3d aprioriPosition = Vector3d::Zero();
	Matrix3d aprioriPositionVariance = Matrix3d::Zero();
	Vector3d minimumConstraintApriori = Vector3d::Zero();
	Vector3d geodeticPosition = Vector3d::Zero();
	bool ready = false;
	Vector3d antennaBoresight = Vector3d::Zero();
	Vector3d antennaAzimuth = Vector3d::Zero();
	string traceFilename;
	string jsonTraceFilename;
	string sppOutputFile;
	map<SatSys, ZhangCheckpointTime> savedSlips;
	map<int, vector<int>> trackedSignals;
	unsigned int failureFlags = 0;
	std::size_t resetObservationCount = 0;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & mapKey;
		ar & id;
		ar & pseudoReceiver;
		ar & invalid;
		ar & sinex;
		ar & metadata;
		ar & metadataMap;
		ar & source;
		ar & firstEpoch;
		ar & lastEpoch;
		ar & epochCount;
		ar & observationCount;
		ar & slipCount;
		ar & codeCount;
		ar & satelliteCount;
		ar & receiverErrorEpochs;
		ar & receiverErrorCount;
		ar & solution;
		ar & antennaType;
		ar & receiverType;
		ar & antennaId;
		ar & satelliteStatuses;
		ar & antennaDelta;
		ar & attitude;
		ar & primaryApriori;
		ar & aprioriTime;
		ar & aprioriClock;
		ar & aprioriClockVariance;
		ar & aprioriPosition;
		ar & aprioriPositionVariance;
		ar & minimumConstraintApriori;
		ar & geodeticPosition;
		ar & ready;
		ar & antennaBoresight;
		ar & antennaAzimuth;
		ar & traceFilename;
		ar & jsonTraceFilename;
		ar & sppOutputFile;
		ar & savedSlips;
		ar & trackedSignals;
		ar & failureFlags;
		ar & resetObservationCount;
	}
};

struct ZhangReceiverRuntimeEnvelope
{
	std::uint32_t schemaVersion =
		ZHANG_RECEIVER_RUNTIME_CHECKPOINT_SCHEMA_VERSION;
	string sectionName = ZHANG_RECEIVER_RUNTIME_CHECKPOINT_SECTION_NAME;
	string runtimeId;
	string boundary = POST_EPOCH_BOUNDARY;
	string observationPolicy = OBSERVATION_POLICY;
	string cachePolicy = CACHE_POLICY;
	string readyPolicy = READY_POLICY;
	string configurationPointerPolicy = CONFIG_POINTER_POLICY;
	map<string, set<string>> customAliasesMap;
	vector<ZhangRuntimeReceiver> receivers;

	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int& version)
	{
		ar & schemaVersion;
		ar & sectionName;
		ar & runtimeId;
		ar & boundary;
		ar & observationPolicy;
		ar & cachePolicy;
		ar & readyPolicy;
		ar & configurationPointerPolicy;
		ar & customAliasesMap;
		ar & receivers;
	}
};

} // namespace ReceiverAuditSchema

