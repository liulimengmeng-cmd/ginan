#include "iono/stecCovariance.hpp"

#include "common/eigenIncluder.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>

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
           << epoch.arIntegerAmbiguityCoordinateCount << ','
           << epoch.arReceiverSingleDifferenceApplied << ','
           << epoch.arReceiverDatumGroupCount << ','
           << epoch.arDroppedSingletonGroupCount << ','
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

int satelliteNumber(const string& satellite)
{
    if (satellite.size() < 2)
    {
        return std::numeric_limits<int>::max();
    }

    try
    {
        size_t consumed = 0;
        const int number = std::stoi(satellite.substr(1), &consumed);
        if (consumed == satellite.size() - 1)
        {
            return number;
        }
    }
    catch (...)
    {
    }
    return std::numeric_limits<int>::max();
}

bool satelliteLess(const StecCovarianceStateRecord& left, const StecCovarianceStateRecord& right)
{
    const int leftNumber = satelliteNumber(left.satellite);
    const int rightNumber = satelliteNumber(right.satellite);
    if (leftNumber != rightNumber)
    {
        return leftNumber < rightNumber;
    }
    return left.satellite < right.satellite;
}

void appendSatelliteDifferenceMeta(
    std::ostringstream&                           output,
    const StecSatelliteDifferenceCovarianceEpoch& epoch,
    E_StecCovarianceCsvStatus                     status,
    long long                                     upperTriangleCount,
    double                                        maxAbsAsymmetry
)
{
    output << "META,"
           << epoch.gpsWeek << ','
           << std::setprecision(17) << epoch.gpsTow << ','
           << stecCovarianceCsvStatusName(status) << ','
           << epoch.sourceStateCount << ','
           << epoch.states.size() << ','
           << epoch.datums.size() << ','
           << epoch.singletonDatumCount << ','
           << upperTriangleCount << ','
           << maxAbsAsymmetry << ','
           << csvEscape(epoch.posteriorStage) << '\n';
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
        "# GINAN_STEC_COVARIANCE_V3\n"
        "# META,gps_week,gps_tow,status,state_count,upper_triangle_count,"
        "max_abs_asymmetry_tecu2,posterior_stage,ar_routine_invoked,"
        "ar_eligible_ambiguity_count,ar_integer_ambiguity_coordinate_count,"
        "ar_receiver_single_difference_applied,ar_receiver_datum_group_count,"
        "ar_dropped_singleton_group_count,ar_resolved_combination_count,"
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

StecSatelliteDifferenceCovarianceEpoch buildStecSatelliteDifferenceCovarianceEpoch(
    const StecCovarianceCsvEpoch& sourceEpoch
)
{
    using GroupKey = std::tuple<string, string, int>;

    StecSatelliteDifferenceCovarianceEpoch result;
    result.gpsWeek = sourceEpoch.gpsWeek;
    result.gpsTow = sourceEpoch.gpsTow;
    result.posteriorStage = sourceEpoch.posteriorStage;
    result.sourceStateCount = static_cast<int>(sourceEpoch.states.size());

    std::map<GroupKey, std::vector<int>> groups;
    for (int sourceIndex = 0; sourceIndex < result.sourceStateCount; sourceIndex++)
    {
        const auto& sourceState = sourceEpoch.states[sourceIndex];
        const string constellation = sourceState.satellite.empty()
            ? string()
            : sourceState.satellite.substr(0, 1);
        groups[{sourceState.site, constellation, sourceState.stateNumber}].push_back(sourceIndex);
    }

    for (auto& [groupKey, members] : groups)
    {
        std::sort(
            members.begin(),
            members.end(),
            [&](int left, int right)
            {
                return satelliteLess(sourceEpoch.states[left], sourceEpoch.states[right]);
            }
        );

        const auto& [site, constellation, stateNumber] = groupKey;
        const int datumIndex = static_cast<int>(result.datums.size());
        const int referenceSourceIndex = members.front();
        const auto& referenceState = sourceEpoch.states[referenceSourceIndex];

        result.datums.push_back(
            {
                datumIndex,
                site,
                constellation,
                stateNumber,
                referenceState.satellite,
                static_cast<int>(members.size())
            }
        );

        if (members.size() < 2)
        {
            result.singletonDatumCount++;
            continue;
        }

        for (size_t memberIndex = 1; memberIndex < members.size(); memberIndex++)
        {
            const int targetSourceIndex = members[memberIndex];
            const auto& targetState = sourceEpoch.states[targetSourceIndex];
            const int differenceLocalIndex = static_cast<int>(result.states.size());

            result.states.push_back(
                {
                    differenceLocalIndex,
                    datumIndex,
                    site,
                    constellation,
                    targetState.satellite,
                    referenceState.satellite,
                    stateNumber,
                    targetSourceIndex,
                    referenceSourceIndex,
                    targetState.estimateTecu - referenceState.estimateTecu
                }
            );
            result.transform.push_back({differenceLocalIndex, targetSourceIndex, +1});
            result.transform.push_back({differenceLocalIndex, referenceSourceIndex, -1});
        }
    }

    const int differenceStateCount = static_cast<int>(result.states.size());
    const long long expectedSourceCovarianceCount =
        static_cast<long long>(result.sourceStateCount) * result.sourceStateCount;
    if (differenceStateCount == 0 ||
        sourceEpoch.covarianceTecu2.size() !=
            static_cast<size_t>(expectedSourceCovarianceCount))
    {
        return result;
    }

    MatrixXd transform = MatrixXd::Zero(differenceStateCount, result.sourceStateCount);
    for (const auto& term : result.transform)
    {
        transform(term.differenceLocalIndex, term.sourceLocalIndex) = term.coefficient;
    }

    VectorXd sourceEstimate(result.sourceStateCount);
    MatrixXd sourceCovariance(result.sourceStateCount, result.sourceStateCount);
    for (int row = 0; row < result.sourceStateCount; row++)
    {
        sourceEstimate(row) = sourceEpoch.states[row].estimateTecu;
        for (int column = 0; column < result.sourceStateCount; column++)
        {
            sourceCovariance(row, column) =
                sourceEpoch.covarianceTecu2[row * result.sourceStateCount + column];
        }
    }

    const VectorXd differenceEstimate = transform * sourceEstimate;
    const MatrixXd differenceCovariance =
        transform * sourceCovariance * transform.transpose();

    result.covarianceTecu2.resize(differenceStateCount * differenceStateCount);
    for (int row = 0; row < differenceStateCount; row++)
    {
        result.states[row].estimateTecu = differenceEstimate(row);
        for (int column = 0; column < differenceStateCount; column++)
        {
            result.covarianceTecu2[row * differenceStateCount + column] =
                differenceCovariance(row, column);
        }
    }

    return result;
}

string stecSatelliteDifferenceCovarianceCsvSchema()
{
    return
        "# GINAN_STEC_SATELLITE_DIFFERENCE_COVARIANCE_V1\n"
        "# META,gps_week,gps_tow,status,source_state_count,difference_state_count,"
        "datum_count,singleton_datum_count,upper_triangle_count,"
        "max_abs_asymmetry_tecu2,posterior_stage\n"
        "# DATUM,gps_week,gps_tow,datum_index,site,constellation,state_number,"
        "reference_satellite,member_count,status\n"
        "# TRANSFORM,gps_week,gps_tow,difference_local_index,source_local_index,"
        "coefficient\n"
        "# SD_STATE,gps_week,gps_tow,local_index,datum_index,site,constellation,"
        "target_satellite,reference_satellite,state_number,target_source_local_index,"
        "reference_source_local_index,estimate_tecu,variance_tecu2\n"
        "# SD_COV,gps_week,gps_tow,row_local_index,column_local_index,covariance_tecu2\n";
}

StecCovarianceCsvResult serializeStecSatelliteDifferenceCovarianceCsvEpoch(
    const StecSatelliteDifferenceCovarianceEpoch& epoch,
    int                                           maxSourceStates,
    double                                        relativeSymmetryTolerance
)
{
    StecCovarianceCsvResult result;
    std::ostringstream output;
    const long long differenceStateCount = static_cast<long long>(epoch.states.size());

    if (epoch.sourceStateCount <= 0 || differenceStateCount == 0)
    {
        result.status = E_StecCovarianceCsvStatus::NO_STATES;
    }
    else if (maxSourceStates <= 0 || epoch.sourceStateCount > maxSourceStates)
    {
        result.status = E_StecCovarianceCsvStatus::STATE_LIMIT_EXCEEDED;
    }
    else if (epoch.covarianceTecu2.size() !=
             static_cast<size_t>(differenceStateCount * differenceStateCount))
    {
        result.status = E_StecCovarianceCsvStatus::DIMENSION_MISMATCH;
    }
    else
    {
        bool finite = true;
        double maximumAbsoluteCovariance = 0;
        double maximumAbsoluteAsymmetry = 0;
        for (int row = 0; row < differenceStateCount; row++)
        {
            finite = finite && std::isfinite(epoch.states[row].estimateTecu);
            for (int column = 0; column < differenceStateCount; column++)
            {
                const double covariance =
                    epoch.covarianceTecu2[row * differenceStateCount + column];
                finite = finite && std::isfinite(covariance);
                maximumAbsoluteCovariance =
                    std::max(maximumAbsoluteCovariance, std::abs(covariance));
                maximumAbsoluteAsymmetry = std::max(
                    maximumAbsoluteAsymmetry,
                    std::abs(
                        covariance - epoch.covarianceTecu2[
                            column * differenceStateCount + row
                        ]
                    )
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
    }

    const bool hasCompleteCovariance =
        epoch.covarianceTecu2.size() ==
        static_cast<size_t>(differenceStateCount * differenceStateCount);
    const long long upperTriangleCount = hasCompleteCovariance
        ? differenceStateCount * (differenceStateCount + 1) / 2
        : 0;
    appendSatelliteDifferenceMeta(
        output,
        epoch,
        result.status,
        upperTriangleCount,
        result.maxAbsAsymmetryTecu2
    );

    output << std::setprecision(17);
    for (const auto& datum : epoch.datums)
    {
        output << "DATUM,"
               << epoch.gpsWeek << ','
               << epoch.gpsTow << ','
               << datum.datumIndex << ','
               << csvEscape(datum.site) << ','
               << csvEscape(datum.constellation) << ','
               << datum.stateNumber << ','
               << csvEscape(datum.referenceSatellite) << ','
               << datum.memberCount << ','
               << (datum.memberCount >= 2 ? "DIFFERENCES_CREATED" : "SINGLETON_DROPPED")
               << '\n';
    }

    if (result.status != E_StecCovarianceCsvStatus::OK)
    {
        result.payload = output.str();
        return result;
    }

    for (const auto& term : epoch.transform)
    {
        output << "TRANSFORM,"
               << epoch.gpsWeek << ','
               << epoch.gpsTow << ','
               << term.differenceLocalIndex << ','
               << term.sourceLocalIndex << ','
               << term.coefficient << '\n';
    }

    for (const auto& state : epoch.states)
    {
        output << "SD_STATE,"
               << epoch.gpsWeek << ','
               << epoch.gpsTow << ','
               << state.localIndex << ','
               << state.datumIndex << ','
               << csvEscape(state.site) << ','
               << csvEscape(state.constellation) << ','
               << csvEscape(state.targetSatellite) << ','
               << csvEscape(state.referenceSatellite) << ','
               << state.stateNumber << ','
               << state.targetSourceLocalIndex << ','
               << state.referenceSourceLocalIndex << ','
               << state.estimateTecu << ','
               << epoch.covarianceTecu2[
                      state.localIndex * differenceStateCount + state.localIndex
                  ]
               << '\n';
    }

    for (int row = 0; row < differenceStateCount; row++)
    {
        for (int column = row; column < differenceStateCount; column++)
        {
            output << "SD_COV,"
                   << epoch.gpsWeek << ','
                   << epoch.gpsTow << ','
                   << row << ','
                   << column << ','
                   << epoch.covarianceTecu2[row * differenceStateCount + column]
                   << '\n';
        }
    }

    result.payload = output.str();
    return result;
}
