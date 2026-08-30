#include <cmath>
#include <iostream>
#include <limits>

#include "iono/stecCovariance.hpp"

namespace
{
StecCovarianceCsvEpoch twoStateEpoch()
{
    StecCovarianceCsvEpoch epoch;
    epoch.gpsWeek = 2323;
    epoch.gpsTow = 12345.5;
    epoch.posteriorStage = "FILTER_POSTERIOR_AFTER_AR_PSEUDOOBS_SUBMITTED_UNVERIFIED";
    epoch.arRoutineInvoked = true;
    epoch.arEligibleAmbiguityCount = 5;
    epoch.arIntegerAmbiguityCoordinateCount = 3;
    epoch.arReceiverSingleDifferenceApplied = true;
    epoch.arReceiverDatumGroupCount = 2;
    epoch.arDroppedSingletonGroupCount = 0;
    epoch.arResolvedCombinationCount = 2;
    epoch.arPseudoObservationsSubmitted = true;
    epoch.arMode = "LAMBDA_ALT";
    epoch.arConfiguredSuccessRateThreshold = 0.9999;
    epoch.arConfiguredSolutionRatioThreshold = 3;
    epoch.arDiagnosticStatus = "RESOLVED_RATIO_ACCEPTED";
    epoch.arSelectedDecorrelatedAmbiguityCount = 2;
    epoch.arIntegerCandidateCount = 2;
    epoch.arBootstrappedSuccessRate = 0.99995;
    epoch.arBestSquaredNorm = 1.25;
    epoch.arSecondSquaredNorm = 5;
    epoch.arSolutionRatio = 4;
    epoch.states = {
        {7, "ALIC", "G01", 0, 12.5},
        {9, "HOB2", "G03", 0, 8.25}
    };
    epoch.covarianceTecu2 = {4.0, -0.75, -0.75, 2.25};
    return epoch;
}

bool check(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
    }
    return condition;
}

bool near(double actual, double expected, double tolerance = 1e-12)
{
    return std::abs(actual - expected) <= tolerance;
}
}  // namespace

int main()
{
    bool passed = true;

    auto result = serializeStecCovarianceCsvEpoch(twoStateEpoch(), 8);
    passed &= check(result.status == E_StecCovarianceCsvStatus::OK, "valid covariance status");
    passed &= check(result.maxAbsAsymmetryTecu2 == 0, "valid covariance symmetry");
    passed &= check(
        result.payload.find("META,2323,12345.5,OK,2,3,0,") != std::string::npos,
        "metadata row"
    );
    passed &= check(
        result.payload.find(
            "FILTER_POSTERIOR_AFTER_AR_PSEUDOOBS_SUBMITTED_UNVERIFIED,1,5,3,1,2,0,2,1,"
            "LAMBDA_ALT,0.99990000000000001,3,RESOLVED_RATIO_ACCEPTED,2,2,"
            "0.99995000000000001,1.25,5,4"
        ) != std::string::npos,
        "ambiguity-resolution capture context"
    );
    passed &= check(
        result.payload.find("STATE,2323,12345.5,0,ALIC,G01,0,7,12.5,4") != std::string::npos,
        "state row"
    );
    passed &= check(
        result.payload.find("COV,2323,12345.5,0,1,-0.75") != std::string::npos,
        "covariance row"
    );

    auto limited = serializeStecCovarianceCsvEpoch(twoStateEpoch(), 1);
    passed &= check(
        limited.status == E_StecCovarianceCsvStatus::STATE_LIMIT_EXCEEDED,
        "state limit status"
    );
    passed &= check(
        limited.payload.find("STATE_LIMIT_EXCEEDED,2,0") != std::string::npos,
        "state limit metadata"
    );
    passed &= check(limited.payload.find("STATE,") == std::string::npos, "no partial state rows");
    passed &= check(limited.payload.find("COV,") == std::string::npos, "no partial covariance rows");

    auto asymmetric = twoStateEpoch();
    asymmetric.covarianceTecu2[2] = -0.5;
    auto asymmetricResult = serializeStecCovarianceCsvEpoch(asymmetric, 8, 1e-12);
    passed &= check(
        asymmetricResult.status == E_StecCovarianceCsvStatus::ASYMMETRY_EXCEEDED,
        "asymmetry status"
    );
    passed &= check(asymmetricResult.maxAbsAsymmetryTecu2 == 0.25, "asymmetry magnitude");

    auto nonfinite = twoStateEpoch();
    nonfinite.covarianceTecu2[1] = std::numeric_limits<double>::quiet_NaN();
    nonfinite.covarianceTecu2[2] = std::numeric_limits<double>::quiet_NaN();
    auto nonfiniteResult = serializeStecCovarianceCsvEpoch(nonfinite, 8);
    passed &= check(
        nonfiniteResult.status == E_StecCovarianceCsvStatus::NONFINITE_VALUE,
        "nonfinite status"
    );

    const auto schema = stecCovarianceCsvSchema();
    passed &= check(schema.find("GINAN_STEC_COVARIANCE_V3") != std::string::npos, "schema version");
    passed &= check(schema.find("estimate_tecu") != std::string::npos, "estimate units");
    passed &= check(schema.find("covariance_tecu2") != std::string::npos, "covariance units");
    passed &= check(schema.find("posterior_stage") != std::string::npos, "posterior semantics");
    passed &= check(
        schema.find("ar_resolved_combination_count") != std::string::npos,
        "ambiguity-resolution evidence fields"
    );
    passed &= check(
        schema.find("ar_integer_ambiguity_coordinate_count") != std::string::npos,
        "integer-estimable coordinate evidence"
    );

    StecCovarianceCsvEpoch oneReceiver;
    oneReceiver.gpsWeek = 2323;
    oneReceiver.gpsTow = 12345.5;
    oneReceiver.posteriorStage = "FILTER_POSTERIOR";
    oneReceiver.states = {
        {0, "A", "G01", 0, 10},
        {1, "A", "G02", 0, 12},
        {2, "A", "G03", 0, 15}
    };
    oneReceiver.covarianceTecu2 = {
        4,   1,   0.5,
        1,   9,   2,
        0.5, 2,   16
    };
    const auto oneReceiverDifference =
        buildStecSatelliteDifferenceCovarianceEpoch(oneReceiver);
    passed &= check(oneReceiverDifference.datums.size() == 1, "single STEC datum group");
    passed &= check(
        oneReceiverDifference.datums[0].referenceSatellite == "G01",
        "deterministic lowest-numbered STEC reference"
    );
    passed &= check(oneReceiverDifference.states.size() == 2, "three states make two differences");
    passed &= check(
        oneReceiverDifference.states[0].targetSatellite == "G02" &&
            oneReceiverDifference.states[1].targetSatellite == "G03",
        "deterministic STEC difference order"
    );
    passed &= check(
        near(oneReceiverDifference.states[0].estimateTecu, 2) &&
            near(oneReceiverDifference.states[1].estimateTecu, 5),
        "STEC satellite-difference estimates"
    );
    passed &= check(
        oneReceiverDifference.covarianceTecu2.size() == 4 &&
            near(oneReceiverDifference.covarianceTecu2[0], 11) &&
            near(oneReceiverDifference.covarianceTecu2[1], 4.5) &&
            near(oneReceiverDifference.covarianceTecu2[2], 4.5) &&
            near(oneReceiverDifference.covarianceTecu2[3], 19),
        "global D C D-transpose non-diagonal formula"
    );
    passed &= check(
        oneReceiverDifference.transform.size() == 4,
        "two sparse transform terms per STEC difference"
    );

    auto differenceSerialized = serializeStecSatelliteDifferenceCovarianceCsvEpoch(
        oneReceiverDifference,
        8
    );
    passed &= check(
        differenceSerialized.status == E_StecCovarianceCsvStatus::OK,
        "valid STEC satellite-difference covariance status"
    );
    passed &= check(
        differenceSerialized.payload.find("DATUM,2323,12345.5,0,A,G,0,G01,3,") !=
            std::string::npos,
        "STEC datum record"
    );
    passed &= check(
        differenceSerialized.payload.find("TRANSFORM,2323,12345.5,0,1,1") !=
            std::string::npos &&
            differenceSerialized.payload.find("TRANSFORM,2323,12345.5,0,0,-1") !=
                std::string::npos,
        "STEC transform records"
    );
    passed &= check(
        differenceSerialized.payload.find("SD_STATE,2323,12345.5,0,0,A,G,G02,G01") !=
            std::string::npos &&
            differenceSerialized.payload.find("SD_COV,2323,12345.5,0,1,4.5") !=
                std::string::npos,
        "STEC satellite-difference state and covariance records"
    );

    StecCovarianceCsvEpoch crossReceiver;
    crossReceiver.states = {
        {0, "A", "G01", 0, 10},
        {1, "A", "G02", 0, 12},
        {2, "B", "G03", 0, 20},
        {3, "B", "G04", 0, 23}
    };
    crossReceiver.covarianceTecu2 = {
        10, 0,   0.4, 0.1,
        0,  10,  0.2, 0.7,
        0.4,0.2, 10,  0,
        0.1,0.7, 0,   10
    };
    const auto crossReceiverDifference =
        buildStecSatelliteDifferenceCovarianceEpoch(crossReceiver);
    passed &= check(
        crossReceiverDifference.states.size() == 2,
        "two receiver groups produce two STEC differences"
    );
    passed &= check(
        near(crossReceiverDifference.covarianceTecu2[1], 0.8) &&
            near(crossReceiverDifference.covarianceTecu2[2], 0.8),
        "global transform retains cross-receiver STEC covariance"
    );

    StecCovarianceCsvEpoch grouped;
    grouped.states = {
        {0, "A", "G03", 0, 3},
        {1, "A", "G01", 0, 1},
        {2, "A", "G05", 1, 5},
        {3, "A", "G02", 1, 2},
        {4, "A", "E11", 0, 11},
        {5, "B", "G08", 0, 8}
    };
    grouped.covarianceTecu2.assign(36, 0);
    for (int index = 0; index < 6; index++)
    {
        grouped.covarianceTecu2[index * 6 + index] = 1;
    }
    const auto groupedDifference = buildStecSatelliteDifferenceCovarianceEpoch(grouped);
    passed &= check(groupedDifference.datums.size() == 4, "site-system-state STEC grouping");
    passed &= check(groupedDifference.singletonDatumCount == 2, "STEC singleton accounting");
    passed &= check(groupedDifference.states.size() == 2, "singletons create no difference rows");
    passed &= check(
        groupedDifference.transform.size() == 4,
        "multi-group transform contains only non-singleton rows"
    );
    bool foundGpsState0Reference = false;
    bool foundGpsState1Reference = false;
    for (const auto& datum : groupedDifference.datums)
    {
        if (datum.site == "A" && datum.constellation == "G" && datum.stateNumber == 0)
        {
            foundGpsState0Reference = datum.referenceSatellite == "G01";
        }
        if (datum.site == "A" && datum.constellation == "G" && datum.stateNumber == 1)
        {
            foundGpsState1Reference = datum.referenceSatellite == "G02";
        }
    }
    passed &= check(
        foundGpsState0Reference && foundGpsState1Reference,
        "state number separates deterministic STEC datum groups"
    );
    const auto groupedSerialized =
        serializeStecSatelliteDifferenceCovarianceCsvEpoch(groupedDifference, 8);
    passed &= check(
        groupedSerialized.status == E_StecCovarianceCsvStatus::OK &&
            groupedSerialized.payload.find("SINGLETON_DROPPED") != std::string::npos,
        "singleton datum is retained in a valid multi-group sidecar"
    );

    const auto differenceSchema = stecSatelliteDifferenceCovarianceCsvSchema();
    passed &= check(
        differenceSchema.find("GINAN_STEC_SATELLITE_DIFFERENCE_COVARIANCE_V1") !=
            std::string::npos,
        "STEC satellite-difference schema version"
    );
    passed &= check(
        differenceSchema.find("# DATUM") != std::string::npos &&
            differenceSchema.find("# TRANSFORM") != std::string::npos &&
            differenceSchema.find("# SD_STATE") != std::string::npos &&
            differenceSchema.find("# SD_COV") != std::string::npos,
        "STEC satellite-difference schema record types"
    );

    if (passed)
    {
        std::cout << "PASS: stec_covariance_tests\n";
        return 0;
    }
    return 1;
}
