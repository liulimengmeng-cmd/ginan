#include "iono/stecCovariance.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

using std::string;

namespace
{
string csvEscape(const string& value)
{
    if (value.find_first_of(",\"\r\n") == string::npos)
    {
        return value;
    }

    string escaped = "\"";
    for (char character : value)
    {
        if (character == '"')
        {
            escaped += '"';
        }
        escaped += character;
    }
    escaped += '"';
    return escaped;
}

void appendMeta(
    std::ostringstream&           output,
    const StecCovarianceCsvEpoch& epoch,
    E_StecCovarianceCsvStatus     status,
    long long                     upperTriangleCount,
    double                        maxAbsAsymmetry
)
{
    output << "META,"
           << epoch.gpsWeek << ','
           << std::setprecision(17) << epoch.gpsTow << ','
           << stecCovarianceCsvStatusName(status) << ','
           << epoch.states.size() << ','
           << upperTriangleCount << ','
           << maxAbsAsymmetry << ','
           << csvEscape(epoch.posteriorStage) << ','
           << epoch.arRoutineInvoked << ','
           << epoch.arEligibleAmbiguityCount << ','
           << epoch.arResolvedCombinationCount << ','
           << epoch.arPseudoObservationsSubmitted << ','
           << csvEscape(epoch.arMode) << ','
           << epoch.arConfiguredSuccessRateThreshold << ','
           << epoch.arConfiguredSolutionRatioThreshold << ','
           << csvEscape(epoch.arDiagnosticStatus) << ','
           << epoch.arSelectedDecorrelatedAmbiguityCount << ','
           << epoch.arIntegerCandidateCount << ','
           << epoch.arBootstrappedSuccessRate << ','
           << epoch.arBestSquaredNorm << ','
           << epoch.arSecondSquaredNorm << ','
           << epoch.arSolutionRatio << '\n';
}
}  // namespace

const char* stecCovarianceCsvStatusName(E_StecCovarianceCsvStatus status)
{
    switch (status)
    {
        case E_StecCovarianceCsvStatus::OK:
            return "OK";
        case E_StecCovarianceCsvStatus::NO_STATES:
            return "NO_STATES";
        case E_StecCovarianceCsvStatus::STATE_LIMIT_EXCEEDED:
            return "STATE_LIMIT_EXCEEDED";
        case E_StecCovarianceCsvStatus::DIMENSION_MISMATCH:
            return "DIMENSION_MISMATCH";
        case E_StecCovarianceCsvStatus::NONFINITE_VALUE:
            return "NONFINITE_VALUE";
        case E_StecCovarianceCsvStatus::ASYMMETRY_EXCEEDED:
            return "ASYMMETRY_EXCEEDED";
    }
    return "UNKNOWN";
}

string stecCovarianceCsvSchema()
{
    return
        "# GINAN_STEC_COVARIANCE_V2\n"
        "# META,gps_week,gps_tow,status,state_count,upper_triangle_count,"
        "max_abs_asymmetry_tecu2,posterior_stage,ar_routine_invoked,"
        "ar_eligible_ambiguity_count,ar_resolved_combination_count,"
        "ar_pseudoobservations_submitted,ar_mode,ar_configured_success_rate_threshold,"
        "ar_configured_solution_ratio_threshold,ar_diagnostic_status,"
        "ar_selected_decorrelated_ambiguity_count,ar_integer_candidate_count,"
        "ar_bootstrapped_success_rate,ar_best_squared_norm,ar_second_squared_norm,"
        "ar_solution_ratio\n"
        "# STATE,gps_week,gps_tow,local_index,site,satellite,state_number,"
        "filter_index,estimate_tecu,variance_tecu2\n"
        "# COV,gps_week,gps_tow,row_local_index,column_local_index,covariance_tecu2\n";
}

StecCovarianceCsvResult serializeStecCovarianceCsvEpoch(
    const StecCovarianceCsvEpoch& epoch,
    int                           maxStates,
    double                        relativeSymmetryTolerance
)
{
    StecCovarianceCsvResult result;
    std::ostringstream      output;
    const long long         stateCount = static_cast<long long>(epoch.states.size());

    if (stateCount == 0)
    {
        result.status = E_StecCovarianceCsvStatus::NO_STATES;
        appendMeta(output, epoch, result.status, 0, 0);
        result.payload = output.str();
        return result;
    }

    if (maxStates <= 0 || stateCount > maxStates)
    {
        result.status = E_StecCovarianceCsvStatus::STATE_LIMIT_EXCEEDED;
        appendMeta(output, epoch, result.status, 0, 0);
        result.payload = output.str();
        return result;
    }

    if (epoch.covarianceTecu2.size() != static_cast<size_t>(stateCount * stateCount))
    {
        result.status = E_StecCovarianceCsvStatus::DIMENSION_MISMATCH;
        appendMeta(output, epoch, result.status, 0, 0);
        result.payload = output.str();
        return result;
    }

    bool   finite = true;
    double maximumAbsoluteCovariance = 0;
    double maximumAbsoluteAsymmetry = 0;
    for (int row = 0; row < stateCount; row++)
    {
        finite = finite && std::isfinite(epoch.states[row].estimateTecu);
        for (int column = 0; column < stateCount; column++)
        {
            const double covariance = epoch.covarianceTecu2[row * stateCount + column];
            finite = finite && std::isfinite(covariance);
            maximumAbsoluteCovariance = std::max(maximumAbsoluteCovariance, std::abs(covariance));
            maximumAbsoluteAsymmetry = std::max(
                maximumAbsoluteAsymmetry,
                std::abs(covariance - epoch.covarianceTecu2[column * stateCount + row])
            );
        }
    }

    result.maxAbsAsymmetryTecu2 = maximumAbsoluteAsymmetry;
    const double symmetryLimit = std::max(0.0, relativeSymmetryTolerance) *
                                 std::max(1.0, maximumAbsoluteCovariance);

    if (!finite)
    {
        result.status = E_StecCovarianceCsvStatus::NONFINITE_VALUE;
    }
    else if (maximumAbsoluteAsymmetry > symmetryLimit)
    {
        result.status = E_StecCovarianceCsvStatus::ASYMMETRY_EXCEEDED;
    }
    else
    {
        result.status = E_StecCovarianceCsvStatus::OK;
    }

    const long long upperTriangleCount = stateCount * (stateCount + 1) / 2;
    appendMeta(output, epoch, result.status, upperTriangleCount, maximumAbsoluteAsymmetry);

    output << std::setprecision(17);
    for (int localIndex = 0; localIndex < stateCount; localIndex++)
    {
        const auto& state = epoch.states[localIndex];
        output << "STATE,"
               << epoch.gpsWeek << ','
               << epoch.gpsTow << ','
               << localIndex << ','
               << csvEscape(state.site) << ','
               << csvEscape(state.satellite) << ','
               << state.stateNumber << ','
               << state.filterIndex << ','
               << state.estimateTecu << ','
               << epoch.covarianceTecu2[localIndex * stateCount + localIndex] << '\n';
    }

    for (int row = 0; row < stateCount; row++)
    {
        for (int column = row; column < stateCount; column++)
        {
            output << "COV,"
                   << epoch.gpsWeek << ','
                   << epoch.gpsTow << ','
                   << row << ','
                   << column << ','
                   << epoch.covarianceTecu2[row * stateCount + column] << '\n';
        }
    }

    result.payload = output.str();
    return result;
}
