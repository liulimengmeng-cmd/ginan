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
            "FILTER_POSTERIOR_AFTER_AR_PSEUDOOBS_SUBMITTED_UNVERIFIED,1,5,2,1,"
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
    passed &= check(schema.find("GINAN_STEC_COVARIANCE_V2") != std::string::npos, "schema version");
    passed &= check(schema.find("estimate_tecu") != std::string::npos, "estimate units");
    passed &= check(schema.find("covariance_tecu2") != std::string::npos, "covariance units");
    passed &= check(schema.find("posterior_stage") != std::string::npos, "posterior semantics");
    passed &= check(
        schema.find("ar_resolved_combination_count") != std::string::npos,
        "ambiguity-resolution evidence fields"
    );

    if (passed)
    {
        std::cout << "PASS: stec_covariance_tests\n";
        return 0;
    }
    return 1;
}
