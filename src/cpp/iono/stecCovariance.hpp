#pragma once

#include <string>
#include <vector>

struct StecCovarianceStateRecord
{
    int         filterIndex = -1;
    std::string site;
    std::string satellite;
    int         stateNumber = 0;
    double      estimateTecu = 0;
};

struct StecCovarianceArContext
{
    bool        routineInvoked = false;
    int         eligibleAmbiguityCount = 0;
    int         integerAmbiguityCoordinateCount = 0;
    bool        receiverSingleDifferenceApplied = false;
    int         receiverDatumGroupCount = 0;
    int         droppedSingletonGroupCount = 0;
    int         resolvedCombinationCount = 0;
    bool        pseudoObservationsSubmitted = false;
    std::string mode = "OFF";
    double      configuredSuccessRateThreshold = 0;
    double      configuredSolutionRatioThreshold = 0;
    std::string diagnosticStatus = "NOT_RUN";
    int         selectedDecorrelatedAmbiguityCount = 0;
    int         integerCandidateCount = 0;
    double      bootstrappedSuccessRate = -1;
    double      bestSquaredNorm = -1;
    double      secondSquaredNorm = -1;
    double      solutionRatio = -1;
};

enum class E_StecCovarianceCsvStatus
{
    OK,
    NO_STATES,
    STATE_LIMIT_EXCEEDED,
    DIMENSION_MISMATCH,
    NONFINITE_VALUE,
    ASYMMETRY_EXCEEDED
};

struct StecCovarianceCsvEpoch
{
    int                                    gpsWeek = 0;
    double                                 gpsTow = 0;
    std::string                            posteriorStage;
    bool                                   arRoutineInvoked = false;
    int                                    arEligibleAmbiguityCount = 0;
    int                                    arIntegerAmbiguityCoordinateCount = 0;
    bool                                   arReceiverSingleDifferenceApplied = false;
    int                                    arReceiverDatumGroupCount = 0;
    int                                    arDroppedSingletonGroupCount = 0;
    int                                    arResolvedCombinationCount = 0;
    bool                                   arPseudoObservationsSubmitted = false;
    std::string                            arMode = "OFF";
    double                                 arConfiguredSuccessRateThreshold = 0;
    double                                 arConfiguredSolutionRatioThreshold = 0;
    std::string                            arDiagnosticStatus = "NOT_RUN";
    int                                    arSelectedDecorrelatedAmbiguityCount = 0;
    int                                    arIntegerCandidateCount = 0;
    double                                 arBootstrappedSuccessRate = -1;
    double                                 arBestSquaredNorm = -1;
    double                                 arSecondSquaredNorm = -1;
    double                                 arSolutionRatio = -1;
    std::vector<StecCovarianceStateRecord> states;
    std::vector<double>                    covarianceTecu2;
};

struct StecCovarianceCsvResult
{
    E_StecCovarianceCsvStatus status = E_StecCovarianceCsvStatus::DIMENSION_MISMATCH;
    std::string               payload;
    double                    maxAbsAsymmetryTecu2 = 0;
};

struct StecSatelliteDifferenceDatumRecord
{
    int         datumIndex = -1;
    std::string site;
    std::string constellation;
    int         stateNumber = 0;
    std::string referenceSatellite;
    int         memberCount = 0;
};

struct StecSatelliteDifferenceTransformRecord
{
    int    differenceLocalIndex = -1;
    int    sourceLocalIndex = -1;
    double coefficient = 0;
};

struct StecSatelliteDifferenceStateRecord
{
    int         localIndex = -1;
    int         datumIndex = -1;
    std::string site;
    std::string constellation;
    std::string targetSatellite;
    std::string referenceSatellite;
    int         stateNumber = 0;
    int         targetSourceLocalIndex = -1;
    int         referenceSourceLocalIndex = -1;
    double      estimateTecu = 0;
};

struct StecSatelliteDifferenceCovarianceEpoch
{
    int                                                 gpsWeek = 0;
    double                                              gpsTow = 0;
    std::string                                         posteriorStage;
    int                                                 sourceStateCount = 0;
    int                                                 singletonDatumCount = 0;
    std::vector<StecSatelliteDifferenceDatumRecord>     datums;
    std::vector<StecSatelliteDifferenceTransformRecord> transform;
    std::vector<StecSatelliteDifferenceStateRecord>     states;
    std::vector<double>                                 covarianceTecu2;
};

const char* stecCovarianceCsvStatusName(E_StecCovarianceCsvStatus status);

std::string stecCovarianceCsvSchema();

StecCovarianceCsvResult serializeStecCovarianceCsvEpoch(
    const StecCovarianceCsvEpoch& epoch,
    int                           maxStates,
    double                        relativeSymmetryTolerance = 1e-10
);

StecSatelliteDifferenceCovarianceEpoch buildStecSatelliteDifferenceCovarianceEpoch(
    const StecCovarianceCsvEpoch& sourceEpoch
);

std::string stecSatelliteDifferenceCovarianceCsvSchema();

StecCovarianceCsvResult serializeStecSatelliteDifferenceCovarianceCsvEpoch(
    const StecSatelliteDifferenceCovarianceEpoch& epoch,
    int                                           maxSourceStates,
    double                                        relativeSymmetryTolerance = 1e-10
);
