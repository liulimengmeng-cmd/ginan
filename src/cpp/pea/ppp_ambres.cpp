// #pragma GCC optimize ("O0")
/**------------------------------------------------------------------------------
 * reference :
 *     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
 *         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
 *         1995
 *     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
 *         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
 *-----------------------------------------------------------------------------*/

#include <algorithm>
#include <iostream>
#include <math.h>
#include <sstream>
#include "ambres/GNSSambres.hpp"
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/biases.hpp"
#include "common/common.hpp"
#include "common/eigenIncluder.hpp"
#include "common/trace.hpp"
#include "pea/ppp.hpp"

static bool filterError = false;

static void traceIntegerComplementRows(
    Trace&                 trace,
    const MatrixXd&        rowsInOriginalAmbiguities,
    const VectorXd&        fixedIntegers,
    const map<int, KFKey>& originalAmbiguityMap
)
{
    const int rowCount = rowsInOriginalAmbiguities.rows();
    const int columnCount = rowsInOriginalAmbiguities.cols();
    if (fixedIntegers.size() != rowCount ||
        static_cast<int>(originalAmbiguityMap.size()) != columnCount)
    {
        tracepdeex(
            1,
            trace,
            "\nPPP_AR INTEGER_COMPLEMENT_ROW stage=2 status=INVALID_DIMENSIONS "
            "rows=%d columns=%d fixed_integers=%d ambiguity_map=%d "
            "action=PROBE_ONLY_NOT_SUBMITTED",
            rowCount,
            columnCount,
            static_cast<int>(fixedIntegers.size()),
            static_cast<int>(originalAmbiguityMap.size())
        );
        return;
    }

    for (int row = 0; row < rowCount; row++)
    {
        int support = 0;
        bool integerCoefficients = true;
        std::ostringstream terms;
        for (int column = 0; column < columnCount; column++)
        {
            const double coefficient = rowsInOriginalAmbiguities(row, column);
            const long long roundedCoefficient = std::llround(coefficient);
            if (std::abs(coefficient - roundedCoefficient) > 1e-9)
            {
                integerCoefficients = false;
            }
            if (roundedCoefficient == 0)
            {
                continue;
            }
            const auto ambiguity = originalAmbiguityMap.find(column);
            if (ambiguity == originalAmbiguityMap.end())
            {
                integerCoefficients = false;
                continue;
            }
            const KFKey& key = ambiguity->second;
            terms << (roundedCoefficient >= 0 ? "+" : "")
                  << roundedCoefficient << " A(" << key.str << ","
                  << key.Sat.id() << "," << key.code() << ") ";
            support++;
        }

        const int traceOutputLevel = integerCoefficients && support > 0 ? 2 : 1;
        tracepdeex(
            traceOutputLevel,
            trace,
            "\nPPP_AR INTEGER_COMPLEMENT_ROW stage=2 row=%d rhs=%.17g "
            "support=%d status=%s terms=%saction=PROBE_ONLY_NOT_SUBMITTED",
            row,
            fixedIntegers(row),
            support,
            integerCoefficients && support > 0 ? "INTEGER_MAPPED" : "INVALID_MAPPING",
            terms.str().c_str()
        );
    }
}

static void traceDualFrequencyCandidateRows(
    Trace&                              trace,
    const DualFrequencyAmbiguityDatum&  datum,
    const MatrixXd&                     rowsInOriginalAmbiguities,
    const VectorXd&                     fixedIntegers,
    int                                 wideLaneCandidateCount,
    const char*                         candidateScope,
    const map<int, KFKey>&              originalAmbiguityMap
)
{
    if (rowsInOriginalAmbiguities.rows() != fixedIntegers.size() ||
        rowsInOriginalAmbiguities.cols() !=
            static_cast<int>(originalAmbiguityMap.size()))
    {
        tracepdeex(
            1,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=%s system=%s "
            "candidate_scope=%s status=INVALID_DIMENSIONS "
            "action=PROBE_ONLY_NOT_SUBMITTED",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            candidateScope
        );
        return;
    }

    for (int row = 0; row < rowsInOriginalAmbiguities.rows(); row++)
    {
        int support = 0;
        bool integerCoefficients = true;
        std::ostringstream terms;
        for (int column = 0; column < rowsInOriginalAmbiguities.cols(); column++)
        {
            const double coefficient = rowsInOriginalAmbiguities(row, column);
            const long long roundedCoefficient = std::llround(coefficient);
            if (std::abs(coefficient - roundedCoefficient) > 1e-9)
            {
                integerCoefficients = false;
            }
            if (roundedCoefficient == 0)
            {
                continue;
            }
            const auto ambiguity = originalAmbiguityMap.find(column);
            if (ambiguity == originalAmbiguityMap.end())
            {
                integerCoefficients = false;
                continue;
            }
            const KFKey& key = ambiguity->second;
            terms << (roundedCoefficient >= 0 ? "+" : "")
                  << roundedCoefficient << " A(" << key.str << ","
                  << key.Sat.id() << "," << key.code() << ") ";
            support++;
        }

        const int traceOutputLevel = integerCoefficients && support > 0 ? 2 : 1;
        tracepdeex(
            traceOutputLevel,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=%s system=%s "
            "reference=%s candidate_scope=%s family=%s row=%d rhs=%.17g "
            "support=%d status=%s "
            "terms=%saction=PROBE_ONLY_NOT_SUBMITTED",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            datum.pivot.id().c_str(),
            candidateScope,
            row < wideLaneCandidateCount ? "WIDE_LANE" : "SECOND_D2",
            row,
            fixedIntegers(row),
            support,
            integerCoefficients && support > 0 ? "INTEGER_MAPPED" : "INVALID_MAPPING",
            terms.str().c_str()
        );
    }
}

struct DualFrequencySubsetProbe
{
    vector<int> selectedEdgeIndices;
    MatrixXd groupTransform;
    MatrixXd combinedTransform;
    GinAR_mtx wideLaneAttempt;
    int wideLaneFixedCount = 0;
    ConditionedIntegerFamily conditionedSecondFamily;
    IterativeIntegerFamilyResolution secondFamilyResolution;
    int combinedRank = 0;
    bool combinedIntegerValued = false;
    bool secondFamilyUnimodular = false;
    bool fullCandidate = false;
    string diagnosticStatus = "NOT_RUN";
};

static DualFrequencySubsetProbe probeDualFrequencySubset(
    Trace&             trace,
    const GinAR_mtx&   originalAmbiguities,
    const MatrixXd&    fullGroupTransform,
    const vector<int>& selectedEdgeIndices,
    const GinAR_opt&   options
)
{
    DualFrequencySubsetProbe result;
    result.selectedEdgeIndices = selectedEdgeIndices;
    const int fullFamilyCount = fullGroupTransform.rows() / 2;
    const int familyCount = selectedEdgeIndices.size();
    if (familyCount <= 0 || fullGroupTransform.rows() != 2 * fullFamilyCount ||
        fullGroupTransform.cols() != originalAmbiguities.aflt.size())
    {
        result.diagnosticStatus = "INVALID_SUBSET_DIMENSIONS";
        return result;
    }

    result.groupTransform = MatrixXd::Zero(
        2 * familyCount,
        originalAmbiguities.aflt.size()
    );
    for (int row = 0; row < familyCount; row++)
    {
        const int edge = selectedEdgeIndices[row];
        if (edge < 0 || edge >= fullFamilyCount)
        {
            result = {};
            result.diagnosticStatus = "INVALID_SUBSET_EDGE_INDEX";
            return result;
        }
        result.groupTransform.row(row) = fullGroupTransform.row(edge);
        result.groupTransform.row(familyCount + row) =
            fullGroupTransform.row(fullFamilyCount + edge);
    }

    GinAR_mtx jointFamilies;
    jointFamilies.aflt = result.groupTransform * originalAmbiguities.aflt;
    jointFamilies.Paflt =
        result.groupTransform * originalAmbiguities.Paflt *
        result.groupTransform.transpose();
    result.wideLaneAttempt.aflt = jointFamilies.aflt.head(familyCount);
    result.wideLaneAttempt.Paflt = jointFamilies.Paflt.topLeftCorner(
        familyCount,
        familyCount
    );
    result.wideLaneFixedCount = GNSS_AR(
        trace,
        result.wideLaneAttempt,
        options
    );
    if (result.wideLaneFixedCount != familyCount)
    {
        result.diagnosticStatus = "WIDE_LANE_FAMILY_NOT_FULL";
        return result;
    }

    result.conditionedSecondFamily = conditionSecondIntegerFamilyOnFixedFirst(
        jointFamilies,
        familyCount,
        result.wideLaneAttempt
    );
    if (result.conditionedSecondFamily.diagnosticStatus !=
        "SECOND_INTEGER_FAMILY_CONDITIONED")
    {
        result.diagnosticStatus = "SECOND_FAMILY_CONDITIONING_REJECTED";
        return result;
    }

    GinAR_opt secondFamilyOptions = options;
    secondFamilyOptions.minimumDecorrelatedAmbiguityCount = 1;
    result.secondFamilyResolution = resolveIntegerFamilyIteratively(
        trace,
        result.conditionedSecondFamily.ambiguityResolution,
        secondFamilyOptions
    );
    const int secondFamilyFixedCount =
        result.secondFamilyResolution.fixedIntegers.size();
    result.combinedTransform = MatrixXd::Zero(
        result.wideLaneFixedCount + secondFamilyFixedCount,
        2 * familyCount
    );
    result.combinedTransform.topRows(result.wideLaneFixedCount) =
        result.conditionedSecondFamily.fixedFirstFamilyTransform;
    if (secondFamilyFixedCount > 0)
    {
        result.combinedTransform.bottomRows(secondFamilyFixedCount) =
            result.secondFamilyResolution.transformToInputCoordinates *
            result.conditionedSecondFamily.secondFamilyTransform;
    }
    result.combinedRank = result.combinedTransform.rows() > 0
        ? result.combinedTransform.fullPivLu().rank()
        : 0;
    result.combinedIntegerValued = result.combinedTransform.rows() > 0 &&
        (result.combinedTransform.array() -
         result.combinedTransform.array().round()).abs().maxCoeff() <= 1e-12;

    if (secondFamilyFixedCount == familyCount)
    {
        const MatrixXd& secondTransform =
            result.secondFamilyResolution.transformToInputCoordinates;
        const double determinant = secondTransform.determinant();
        result.secondFamilyUnimodular =
            (secondTransform.array() - secondTransform.array().round())
                    .abs()
                    .maxCoeff() <= 1e-12 &&
            std::isfinite(determinant) &&
            std::abs(std::abs(determinant) - 1) <= 1e-8;
    }
    result.fullCandidate =
        secondFamilyFixedCount == familyCount &&
        result.secondFamilyUnimodular &&
        result.combinedIntegerValued &&
        result.combinedRank == 2 * familyCount;
    result.diagnosticStatus = result.fullCandidate
        ? "FULL_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED"
        : "SECOND_D2_FAMILY_NOT_FULL";
    return result;
}

static void traceDualFrequencyDatumDiagnostic(
    Trace&                                  trace,
    const GinAR_mtx&                        originalAmbiguities,
    const DualFrequencyAmbiguityTransform&  transform,
    const GinAR_opt&                        options
)
{
    tracepdeex(
        2,
        trace,
        "\nPPP_AR DUAL_FREQUENCY_BASIS original=%d rows=%d expected_rank=%d "
        "actual_rank=%d complete_groups=%d incomplete_groups=%d paired_ambiguities=%d "
        "unmatched_ambiguities=%d integer_valued=%d full_row_rank=%d covers_all=%d "
        "status=%s action=PROBE_ONLY_NOT_SUBMITTED",
        static_cast<int>(originalAmbiguities.aflt.size()),
        static_cast<int>(transform.matrix.rows()),
        transform.expectedIntegerRank,
        transform.actualIntegerRank,
        transform.completeGroupCount,
        transform.incompleteGroupCount,
        transform.pairedAmbiguityCount,
        transform.unmatchedAmbiguityCount,
        transform.integerValued,
        transform.fullRowRank,
        transform.coversEligibleAmbiguities,
        transform.diagnosticStatus.c_str()
    );

    int fullVisibleCandidateGroupCount = 0;
    int fullVisibleCandidateRowCount = 0;
    int selectedCandidateGroupCount = 0;
    int selectedCandidateRowCount = 0;
    for (const auto& datum : transform.groups)
    {
        tracepdeex(
            2,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_GROUP receiver=%s system=%s first_signal=%s "
            "second_signal=%s reference=%s common_satellites=%d unmatched_ambiguities=%d "
            "wide_lane_rows=%d second_d2_rows=%d complete_graph=%d "
            "action=%s",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            enum_to_string(int_to_enum<E_ObsCode>(datum.firstObservation)).c_str(),
            enum_to_string(int_to_enum<E_ObsCode>(datum.secondObservation)).c_str(),
            datum.pivot.id().c_str(),
            datum.commonSatelliteCount,
            datum.unmatchedAmbiguityCount,
            datum.wideLaneRowCount,
            datum.complementRowCount,
            datum.completeSignalGraph,
            datum.completeSignalGraph ?
                "PROBE_ONLY_NOT_SUBMITTED" : "INCOMPLETE_GRAPH_NOT_PROBED"
        );
        if (!datum.completeSignalGraph || datum.wideLaneRowCount <= 0)
        {
            continue;
        }

        const int visibleFamilyCount = datum.wideLaneRowCount;
        const MatrixXd fullGroupTransform = transform.matrix.middleRows(
            datum.wideLaneRowOffset,
            visibleFamilyCount + datum.complementRowCount
        );
        vector<int> allVisibleEdges(visibleFamilyCount);
        for (int edge = 0; edge < visibleFamilyCount; edge++)
        {
            allVisibleEdges[edge] = edge;
        }
        DualFrequencySubsetProbe visibleProbe = probeDualFrequencySubset(
            trace,
            originalAmbiguities,
            fullGroupTransform,
            allVisibleEdges,
            options
        );
        DualFrequencySubsetProbe selectedProbe = visibleProbe;

        // Conventional partial ambiguity resolution is applied to an explicit
        // common-satellite subset.  Rows are ordered by the joint WL+d2 float
        // variance and the largest nested subset passing every integer/rank
        // gate is retained.  Arbitrary isolated rows are never called a datum.
        const int minimumSelectedFamilyCount =
            options.minimumDecorrelatedAmbiguityCount;
        if (!visibleProbe.fullCandidate &&
            visibleFamilyCount > minimumSelectedFamilyCount)
        {
            const MatrixXd jointCovariance =
                fullGroupTransform * originalAmbiguities.Paflt *
                fullGroupTransform.transpose();
            vector<pair<double, int>> rankedEdges;
            rankedEdges.reserve(visibleFamilyCount);
            for (int edge = 0; edge < visibleFamilyCount; edge++)
            {
                double quality =
                    jointCovariance(edge, edge) +
                    jointCovariance(visibleFamilyCount + edge,
                                    visibleFamilyCount + edge);
                if (!std::isfinite(quality))
                {
                    quality = HUGE_VAL;
                }
                rankedEdges.push_back({quality, edge});
            }
            std::stable_sort(
                rankedEdges.begin(),
                rankedEdges.end(),
                [](const auto& left, const auto& right)
                {
                    if (left.first != right.first)
                    {
                        return left.first < right.first;
                    }
                    return left.second < right.second;
                }
            );

            for (int subsetFamilyCount = visibleFamilyCount - 1;
                 subsetFamilyCount >= minimumSelectedFamilyCount;
                 subsetFamilyCount--)
            {
                vector<int> subsetEdges;
                subsetEdges.reserve(subsetFamilyCount);
                for (int rank = 0; rank < subsetFamilyCount; rank++)
                {
                    subsetEdges.push_back(rankedEdges[rank].second);
                }
                std::sort(subsetEdges.begin(), subsetEdges.end());
                DualFrequencySubsetProbe subsetProbe = probeDualFrequencySubset(
                    trace,
                    originalAmbiguities,
                    fullGroupTransform,
                    subsetEdges,
                    options
                );
                if (subsetProbe.fullCandidate)
                {
                    selectedProbe = std::move(subsetProbe);
                    break;
                }
            }
        }

        const int selectedFamilyCount = selectedProbe.selectedEdgeIndices.size();
        const int secondFamilyFixedCount =
            selectedProbe.secondFamilyResolution.fixedIntegers.size();
        const bool secondFamilyRowsIntegerValued = secondFamilyFixedCount == 0 ||
            (selectedProbe.secondFamilyResolution.transformToInputCoordinates.array() -
             selectedProbe.secondFamilyResolution.transformToInputCoordinates.array().round())
                    .abs()
                    .maxCoeff() <= 1e-12;
        const int secondFamilyRowRank = secondFamilyFixedCount > 0
            ? selectedProbe.secondFamilyResolution.transformToInputCoordinates
                  .fullPivLu().rank()
            : 0;
        const bool fullVisibleCandidate =
            visibleProbe.fullCandidate &&
            selectedFamilyCount == visibleFamilyCount;
        if (selectedProbe.fullCandidate)
        {
            selectedCandidateGroupCount++;
            selectedCandidateRowCount += selectedProbe.combinedRank;
            if (fullVisibleCandidate)
            {
                fullVisibleCandidateGroupCount++;
                fullVisibleCandidateRowCount += selectedProbe.combinedRank;
            }
            VectorXd combinedFixedIntegers(2 * selectedFamilyCount);
            combinedFixedIntegers <<
                selectedProbe.wideLaneAttempt.zfix,
                selectedProbe.secondFamilyResolution.fixedIntegers;
            traceDualFrequencyCandidateRows(
                trace,
                datum,
                selectedProbe.combinedTransform * selectedProbe.groupTransform,
                combinedFixedIntegers,
                selectedFamilyCount,
                fullVisibleCandidate ? "FULL_VISIBLE_GROUP" : "SELECTED_SUBSET",
                originalAmbiguities.ambmap
            );
        }
        else if (visibleProbe.secondFamilyResolution.fixedIntegers.size() > 0)
        {
            const MatrixXd acceptedSecondRows =
                visibleProbe.secondFamilyResolution.transformToInputCoordinates *
                visibleProbe.conditionedSecondFamily.secondFamilyTransform *
                visibleProbe.groupTransform;
            traceDualFrequencyCandidateRows(
                trace,
                datum,
                acceptedSecondRows,
                visibleProbe.secondFamilyResolution.fixedIntegers,
                0,
                "SECOND_D2_PARTIAL",
                originalAmbiguities.ambmap
            );
        }

        const GinAR_mtx& lastSecondAttempt =
            selectedProbe.secondFamilyResolution.lastAttempt;

        tracepdeex(
            2,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_STAGE receiver=%s system=%s reference=%s "
            "visible_family_count=%d selected_family_count=%d excluded_satellites=%d "
            "wide_lane_target=%d wide_lane_fixed=%d wide_lane_status=%s "
            "wide_lane_success_rate=%.17g wide_lane_ratio=%.17g "
            "second_d2_target=%d second_d2_fixed=%d second_d2_status=%s "
            "second_d2_minimum_decorrelated=%d "
            "second_d2_attempted_stages=%d second_d2_accepted_stages=%d "
            "second_d2_last_attempt_status=%s second_d2_success_rate=%.17g "
            "second_d2_ratio=%.17g second_d2_integer_valued=%d "
            "second_d2_row_rank=%d second_d2_unimodular=%d "
            "combined_integer_valued=%d "
            "combined_rank=%d target_rank=%d status=%s "
            "action=PROBE_ONLY_NOT_SUBMITTED",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            datum.pivot.id().c_str(),
            visibleFamilyCount,
            selectedFamilyCount,
            visibleFamilyCount - selectedFamilyCount,
            selectedFamilyCount,
            selectedProbe.wideLaneFixedCount,
            selectedProbe.wideLaneAttempt.diagnosticStatus.c_str(),
            selectedProbe.wideLaneAttempt.bootstrappedSuccessRate,
            selectedProbe.wideLaneAttempt.solutionRatio,
            selectedFamilyCount,
            secondFamilyFixedCount,
            selectedProbe.secondFamilyResolution.diagnosticStatus.c_str(),
            1,
            selectedProbe.secondFamilyResolution.attemptedStageCount,
            selectedProbe.secondFamilyResolution.acceptedStageCount,
            lastSecondAttempt.diagnosticStatus.c_str(),
            lastSecondAttempt.bootstrappedSuccessRate,
            lastSecondAttempt.solutionRatio,
            secondFamilyRowsIntegerValued,
            secondFamilyRowRank,
            selectedProbe.secondFamilyUnimodular,
            selectedProbe.combinedIntegerValued,
            selectedProbe.combinedRank,
            2 * selectedFamilyCount,
            selectedProbe.fullCandidate ?
                (fullVisibleCandidate ?
                    "FULL_VISIBLE_GROUP_INTEGER_DATUM_CANDIDATE_UNVERIFIED" :
                    "FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED") :
                selectedProbe.diagnosticStatus.c_str()
        );
    }

    const bool fullVisibleEpochCandidate =
        transform.diagnosticStatus == "FULL_DUAL_FREQUENCY_INTEGER_BASIS" &&
        fullVisibleCandidateGroupCount == transform.completeGroupCount &&
        fullVisibleCandidateRowCount == transform.expectedIntegerRank;
    const bool fullSelectedEpochCandidate =
        transform.incompleteGroupCount == 0 &&
        selectedCandidateGroupCount == transform.completeGroupCount;
    tracepdeex(
        2,
        trace,
        "\nPPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=%d "
        "target_groups=%d full_candidate_rows=%d target_rank=%d "
        "selected_candidate_groups=%d selected_candidate_rows=%d status=%s "
        "wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED",
        fullVisibleCandidateGroupCount,
        transform.completeGroupCount,
        fullVisibleCandidateRowCount,
        transform.expectedIntegerRank,
        selectedCandidateGroupCount,
        selectedCandidateRowCount,
        fullVisibleEpochCandidate ?
            "FULL_VISIBLE_INTEGER_DATUM_CANDIDATE_UNVERIFIED" :
            (fullSelectedEpochCandidate ?
                "FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED" :
                "NO_FULL_INTEGER_DATUM_CANDIDATE")
    );
}

bool recordFilterError(RejectCallbackDetails rejectDetails)
{
    filterError = true;

    return true;
}

bool applyBestIntegerAmbiguity(
    Trace&   trace,   ///< Debug trace
    KFState& kfState  ///< Reference to Kalman filter containing float solutions
)
{
    KFKey  bestKey;
    double smallestVar = 1e10;

    for (auto& [key, index] : kfState.kfIndexMap)
    {
        if (key.type != KF::AMBIGUITY)
        {
            continue;
        }

        double var = kfState.P(index, index);

        if (var > smallestVar || var < FIXED_AMB_VAR * 5)
        {
            continue;
        }

        smallestVar = var;
        bestKey     = key;
    }

    if (bestKey.type == KF::NONE)
    {
        return false;
    }

    KFMeasEntryList kfMeasEntryList;

    int index = kfState.kfIndexMap[bestKey];

    double closest = round(kfState.x(index));

    KFMeasEntry measEntry(&kfState);

    measEntry.obsKey = bestKey;

    measEntry.addDsgnEntry(bestKey, 1);

    measEntry.setValue(closest);
    measEntry.setNoise(FIXED_AMB_VAR);

    kfMeasEntryList.push_back(measEntry);

    KFMeas kfMeas(kfState, kfMeasEntryList, kfState.time);

    filterError = false;
    kfState.measRejectCallbacks.push_back(recordFilterError);
    {
        kfState.filterKalman(trace, kfMeas);
    }
    kfState.measRejectCallbacks.pop_back();

    if (filterError)
    {
        return false;
    }

    kfState.outputStates(trace, "/AR1");

    return true;
}

bool applyUCAmbiguities(
    Trace&     trace,    ///< Debug trace
    KFState&   kfState,  ///< Reference to Kalman filter containing float solutions
    GinAR_mtx& mtrx  ///< Reference to structure containing fixed ambiguities and Z transformations
)
{
    int nz = mtrx.zfix.size();
    int nx = mtrx.ambmap.size();

    tracepdeex(1, trace, "   %d out of %d ambiguities resolved, applying...\n", nz, nx);

    MatrixXd Z    = mtrx.Ztrs;
    VectorXd zfix = mtrx.zfix;

    if (Z.rows() != nz || Z.cols() != nx)
    {
        BOOST_LOG_TRIVIAL(error)
            << "PPP-AR integer pseudo-observation design dimensions are inconsistent";
        tracepdeex(
            1,
            trace,
            "\nPPP_AR PSEUDOOBS_DESIGN rows=%d original_ambiguities=%d "
            "z_rows=%d z_columns=%d status=INVALID_DIMENSIONS",
            nz,
            nx,
            static_cast<int>(Z.rows()),
            static_cast<int>(Z.cols())
        );
        return false;
    }

    vector<KFKey> ambiguityKeys;
    vector<int>   ambiguityStateIndices;
    ambiguityKeys.reserve(nx);
    ambiguityStateIndices.reserve(nx);
    for (int column = 0; column < nx; column++)
    {
        const auto ambiguity = mtrx.ambmap.find(column);
        if (ambiguity == mtrx.ambmap.end())
        {
            BOOST_LOG_TRIVIAL(error)
                << "PPP-AR integer pseudo-observation ambiguity map is not contiguous";
            tracepdeex(
                1,
                trace,
                "\nPPP_AR PSEUDOOBS_DESIGN rows=%d original_ambiguities=%d "
                "status=INVALID_AMBIGUITY_MAP",
                nz,
                nx
            );
            return false;
        }

        const auto state = kfState.kfIndexMap.find(ambiguity->second);
        if (state == kfState.kfIndexMap.end() || state->second < 0 ||
            state->second >= kfState.x.rows())
        {
            BOOST_LOG_TRIVIAL(error)
                << "PPP-AR integer pseudo-observation references a missing ambiguity state";
            tracepdeex(
                1,
                trace,
                "\nPPP_AR PSEUDOOBS_DESIGN rows=%d original_ambiguities=%d "
                "status=MISSING_AMBIGUITY_STATE",
                nz,
                nx
            );
            return false;
        }

        ambiguityKeys.push_back(ambiguity->second);
        ambiguityStateIndices.push_back(state->second);
    }

    if (AR_VERBO)
    {
        trace << "\n"
              << "zfix =" << "\n"
              << zfix.transpose() << "\n";
        trace << "\n"
              << "Ztrs =" << "\n"
              << Z << "\n";
    }

    KFMeasEntryList kfMeasEntryList;

    for (int i = 0; i < nz; i++)
    {
        double residual = zfix(i);

        KFMeasEntry measEntry(&kfState);

        measEntry.obsKey.type    = KF::Z_AMB;
        measEntry.obsKey.num     = i;
        measEntry.obsKey.comment = "Ambiguity Pseudoobs";

        // Each resolved integer combination is an independent pseudo-
        // observation.  Reusing one noise key here makes every row share the
        // same scalar noise source, producing FIXED_AMB_VAR * 11^T (rank one)
        // instead of the intended FIXED_AMB_VAR * I.
        measEntry.setNoise(FIXED_AMB_VAR);

        tracepdeex(4, trace, "      Applying:  ");

        for (int j = 0; j < nx; j++)
        {
            if (Z(i, j) == 0)
            {
                continue;
            }

            double ambiguity = 0;

            const KFKey& key = ambiguityKeys[j];
            kfState.getKFValue(key, ambiguity);

            residual -= Z(i, j) * ambiguity;

            tracepdeex(
                4,
                trace,
                "%+3.0f A(%s,%s,%3s) ",
                Z(i, j),
                key.str.c_str(),
                key.Sat.id().c_str(),
                key.code().c_str()
            );

            InitialState init;
            init.x = ambiguity;
            init.P = 3600;

            measEntry.addDsgnEntry(key, Z(i, j), init);
        }

        tracepdeex(4, trace, "= %+10.5f\n", zfix(i));

        measEntry.setInnov(residual);

        kfMeasEntryList.push_back(measEntry);
    }

    KFMeas kfMeas(kfState, kfMeasEntryList, kfState.time);

    MatrixXd expectedDesign = MatrixXd::Zero(nz, kfState.x.rows());
    for (int column = 0; column < nx; column++)
    {
        // Mirror KFMeasEntry::addDsgnEntry(), which accumulates coefficients
        // when more than one input column refers to the same state key.
        expectedDesign.col(ambiguityStateIndices[column]) += Z.col(column);
    }
    if (kfMeas.H.rows() != expectedDesign.rows() ||
        kfMeas.H.cols() != expectedDesign.cols() ||
        !kfMeas.H.isApprox(expectedDesign, 1e-12))
    {
        BOOST_LOG_TRIVIAL(error)
            << "PPP-AR integer pseudo-observation design does not match mapped Z*D constraints";
        tracepdeex(
            1,
            trace,
            "\nPPP_AR PSEUDOOBS_DESIGN rows=%d original_ambiguities=%d "
            "status=INVALID_Z_TIMES_D_DESIGN",
            nz,
            nx
        );
        return false;
    }

    tracepdeex(
        2,
        trace,
        "\nPPP_AR PSEUDOOBS_DESIGN rows=%d original_ambiguities=%d "
        "status=MATCHES_Z_TIMES_D",
        nz,
        nx
    );

    const MatrixXd expectedNoise = MatrixXd::Identity(nz, nz) * FIXED_AMB_VAR;
    if (kfMeas.R.rows() != expectedNoise.rows() ||
        kfMeas.R.cols() != expectedNoise.cols() ||
        !kfMeas.R.isApprox(expectedNoise, FIXED_AMB_VAR * 1e-12))
    {
        BOOST_LOG_TRIVIAL(error)
            << "PPP-AR integer pseudo-observation noise is not independent diagonal noise";
        tracepdeex(
            1,
            trace,
            "\nPPP_AR PSEUDOOBS_NOISE rows=%d status=INVALID_NON_DIAGONAL",
            nz
        );
        return false;
    }

    tracepdeex(
        2,
        trace,
        "\nPPP_AR PSEUDOOBS_NOISE rows=%d variance=%.3e status=INDEPENDENT_DIAGONAL",
        nz,
        FIXED_AMB_VAR
    );

    kfState.filterKalman(trace, kfMeas, "/AR", true);

    // filterKalman() has no acceptance return value.  This event therefore
    // records only that the checked pseudo-observation block was submitted and
    // the filter call returned; it is deliberately not a correct-fix claim.
    tracepdeex(
        2,
        trace,
        "\nPPP_AR PSEUDOOBS_SUBMISSION rows=%d "
        "status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED",
        nz
    );

    return true;
}

AmbiguityResolutionAttempt fixAndHoldAmbiguities(
    Trace&   trace,   ///< Debug trace
    KFState& kfState  ///< Filter state
)
{
    AmbiguityResolutionAttempt result;

    tracepdeex(3, trace, "%s: %s\n", __FUNCTION__, kfState.time.to_string().c_str());

    if (acsConfig.ambrOpts.mode == E_ARmode::OFF)
    {
        return result;
    }

    result.routineInvoked = true;

    GinAR_mtx        ARmtx;
    map<string, int> nsat;  // number of satellites visible by station
    map<SatSys, int> nsta;  // number of stations visible by satellite

    int         ind = 0;
    vector<int> indices;
    for (auto& [key, index] : kfState.kfIndexMap)
    {
        if (key.type != KF::AMBIGUITY)
        {
            continue;
        }

        if (acsConfig.solve_amb_for[key.Sat.sys] == false)
        {
            continue;
        }

        indices.push_back(index);

        ARmtx.ambmap[ind] = key;
        ind++;
    }

    result.eligibleAmbiguityCount = ind;
    if (ind == 0)
    {
        result.diagnosticStatus = "NO_ELIGIBLE_AMBIGUITIES";
        return result;
    }

    ARmtx.aflt  = kfState.x(indices);
    ARmtx.Paflt = kfState.P(indices, indices);

    const ReceiverAmbiguityTransform integerTransform =
        buildReceiverAmbiguityIntegerTransform(ARmtx, acsConfig.receiver_amb_pivot);
    result.integerAmbiguityCoordinateCount = integerTransform.matrix.rows();
    result.receiverSingleDifferenceApplied =
        integerTransform.singleDifferencedGroupCount > 0;
    result.receiverDatumGroupCount = integerTransform.singleDifferencedGroupCount;
    result.droppedSingletonGroupCount = integerTransform.droppedSingletonGroupCount;

    tracepdeex(
        2,
        trace,
        "\nPPP_AR INTEGER_COORDINATES original=%d integer=%d receiver_sd_groups=%d "
        "identity_groups=%d dropped_singletons=%d",
        ind,
        result.integerAmbiguityCoordinateCount,
        integerTransform.singleDifferencedGroupCount,
        integerTransform.identityGroupCount,
        integerTransform.droppedSingletonGroupCount
    );
    for (const auto& datum : integerTransform.groups)
    {
        if (!datum.singleDifferenced)
        {
            continue;
        }
        tracepdeex(
            2,
            trace,
            "\nPPP_AR RECEIVER_SD receiver=%s system=%s signal=%s reference=%s "
            "members=%d status=%s",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            enum_to_string(int_to_enum<E_ObsCode>(datum.observation)).c_str(),
            datum.pivot.id().c_str(),
            datum.memberCount,
            datum.memberCount >= 2 ? "INTEGER_COORDINATES_CREATED" : "SINGLETON_DROPPED"
        );
    }

    if (integerTransform.matrix.rows() == 0)
    {
        result.diagnosticStatus = "NO_INTEGER_ESTIMABLE_AMBIGUITIES";
        return result;
    }

    GinAR_mtx integerARmtx;
    integerARmtx.aflt = integerTransform.matrix * ARmtx.aflt;
    integerARmtx.Paflt =
        integerTransform.matrix * ARmtx.Paflt * integerTransform.matrix.transpose();
    const GinAR_mtx integerARInput = integerARmtx;

    GinAR_opt ARopt;
    ARopt.mode   = acsConfig.ambrOpts.mode;
    ARopt.sucthr = acsConfig.ambrOpts.succsThres;
    ARopt.ratthr = acsConfig.ambrOpts.ratioThres;
    ARopt.nset   = acsConfig.ambrOpts.lambda_set;
    ARopt.nitr   = acsConfig.ambrOpts.AR_max_itr;

    if (traceLevel > 4)
        AR_VERBO = true;

    if (acsConfig.ambrOpts.integer_complement_diagnostics)
    {
        const DualFrequencyAmbiguityTransform dualFrequencyTransform =
            buildDualFrequencyAmbiguityIntegerTransform(
                ARmtx,
                acsConfig.receiver_amb_pivot
            );
        traceDualFrequencyDatumDiagnostic(
            trace,
            ARmtx,
            dualFrequencyTransform,
            ARopt
        );
    }

    // Resolve and apply ambiguities
    int nfix = GNSS_AR(trace, integerARmtx, ARopt);
    result.resolvedCombinationCount = nfix;
    result.diagnosticStatus = integerARmtx.diagnosticStatus;
    result.selectedDecorrelatedAmbiguityCount =
        integerARmtx.selectedDecorrelatedAmbiguityCount;
    result.integerCandidateCount = integerARmtx.integerCandidateCount;
    result.bootstrappedSuccessRate = integerARmtx.bootstrappedSuccessRate;
    result.bestSquaredNorm = integerARmtx.bestSquaredNorm;
    result.secondSquaredNorm = integerARmtx.secondSquaredNorm;
    result.solutionRatio = integerARmtx.solutionRatio;
    if (acsConfig.ambrOpts.integer_complement_diagnostics && nfix > 0 &&
        nfix < result.integerAmbiguityCoordinateCount)
    {
        const ConditionalIntegerComplement complement =
            buildConditionalIntegerComplement(integerARInput, integerARmtx);
        if (complement.diagnosticStatus == "CONDITIONAL_INTEGER_COMPLEMENT_READY")
        {
            GinAR_mtx complementProbe = complement.ambiguityResolution;
            GinAR_opt complementProbeOptions = ARopt;
            // Probe the one- and two-dimensional remainder without weakening
            // the production stage-one minimum or submitting the result.
            complementProbeOptions.minimumDecorrelatedAmbiguityCount = 1;
            const int complementFixCount =
                GNSS_AR(trace, complementProbe, complementProbeOptions);
            tracepdeex(
                2,
                trace,
                "\nPPP_AR INTEGER_COMPLEMENT_DIAGNOSTIC stage1_rows=%d "
                "remaining_coordinates=%d stage2_probe_rows=%d combined_independent_rows=%d "
                "target_integer_rank=%d stage2_status=%s stage2_success_rate=%.17g "
                "stage2_minimum_decorrelated=%d stage2_selected_decorrelated=%d "
                "stage2_integer_candidates=%d "
                "stage2_best_squared_norm=%.17g stage2_second_squared_norm=%.17g "
                "stage2_ratio=%.17g action=PROBE_ONLY_NOT_SUBMITTED",
                nfix,
                static_cast<int>(complement.ambiguityResolution.aflt.size()),
                complementFixCount,
                nfix + complementFixCount,
                result.integerAmbiguityCoordinateCount,
                complementProbe.diagnosticStatus.c_str(),
                complementProbe.bootstrappedSuccessRate,
                complementProbeOptions.minimumDecorrelatedAmbiguityCount,
                complementProbe.selectedDecorrelatedAmbiguityCount,
                complementProbe.integerCandidateCount,
                complementProbe.bestSquaredNorm,
                complementProbe.secondSquaredNorm,
                complementProbe.solutionRatio
            );
            if (complementFixCount > 0)
            {
                const MatrixXd rowsInOriginalAmbiguities =
                    complementProbe.Ztrs *
                    complement.transformToInputCoordinates *
                    integerTransform.matrix;
                traceIntegerComplementRows(
                    trace,
                    rowsInOriginalAmbiguities,
                    complementProbe.zfix,
                    ARmtx.ambmap
                );
            }
        }
        else
        {
            tracepdeex(
                1,
                trace,
                "\nPPP_AR INTEGER_COMPLEMENT_DIAGNOSTIC stage1_rows=%d "
                "remaining_coordinates=%d stage2_probe_rows=0 combined_independent_rows=%d "
                "target_integer_rank=%d stage2_status=%s stage2_success_rate=-1 "
                "stage2_minimum_decorrelated=1 stage2_selected_decorrelated=0 "
                "stage2_integer_candidates=0 "
                "stage2_best_squared_norm=-1 stage2_second_squared_norm=-1 "
                "stage2_ratio=-1 action=PROBE_ONLY_NOT_SUBMITTED",
                nfix,
                result.integerAmbiguityCoordinateCount - nfix,
                nfix,
                result.integerAmbiguityCoordinateCount,
                complement.diagnosticStatus.c_str()
            );
        }
    }
    if (nfix > 0)
    {
        const int integerCoordinateCount = integerARmtx.Ztrs.cols();
        const int originalAmbiguityCount = integerTransform.matrix.cols();
        if (!mapIntegerAmbiguityConstraintsToOriginalState(
                integerARmtx,
                integerTransform,
                ARmtx.ambmap
            ))
        {
            tracepdeex(
                1,
                trace,
                "\nPPP_AR INTEGER_FEEDBACK rows=%d integer_coordinates=%d "
                "original_ambiguities=%d status=INVALID_DIMENSIONS",
                nfix,
                integerCoordinateCount,
                originalAmbiguityCount
            );
            result.diagnosticStatus = "INTEGER_FEEDBACK_MAPPING_INVALID";
            return result;
        }
        tracepdeex(
            2,
            trace,
            "\nPPP_AR INTEGER_FEEDBACK rows=%d integer_coordinates=%d "
            "original_ambiguities=%d status=Z_TIMES_D_MAPPED",
            nfix,
            integerCoordinateCount,
            originalAmbiguityCount
        );
        result.pseudoObservationsSubmitted =
            applyUCAmbiguities(trace, kfState, integerARmtx);
        if (!result.pseudoObservationsSubmitted)
        {
            result.diagnosticStatus = "PSEUDOOBS_MODEL_INVALID";
        }
    }

    while (0)
    {
        bool applied = applyBestIntegerAmbiguity(trace, kfState);

        if (applied == false)
        {
            break;
        }
    }

    return result;
}

bool queryBiasUC(
    Trace&   trace,    ///< debug stream
    GTime    time,     ///< time of biases
    KFState& kfState,  ///< filter state to take biases from
    SatSys   Sat,    ///< satellite (for receiver biases, sat.sys needs to be set to the appropriate
                     ///< system, and sat.prn must be 0)
    string     rec,  ///< receiver  (for satellite biases nees to be "")
    E_ObsCode  code,  ///< signal code
    double&    bias,  ///< bias value
    double&    var,   ///< bias variance
    E_MeasType type   ///< measurement type
)
{
    KFKey kfKey;
    kfKey.str = rec;
    kfKey.Sat = Sat;
    kfKey.num = static_cast<int>(code);

    if (Sat.prn == 0)  // todo? check if needed and reverse logic
    {
        auto& recOpts = acsConfig.getRecOpts(rec, {Sat.sysName(), enum_to_string(code)});

        if (type == CODE)
        {
            if (recOpts.codeBiasModel.enable == false)
                return true;

            InitialState init = initialStateFromConfig(recOpts.code_bias);
            if (init.estimate == false)
            {
                getBias(trace, time, rec, Sat, code, CODE, bias, var);
                return true;
            }

            kfKey.type = KF::CODE_BIAS;

            return kfState.getKFValue(kfKey, bias, &var) != E_Source::NONE;
        }

        if (type == PHAS)
        {
            if (recOpts.phaseBiasModel.enable == false)
                return true;

            InitialState init = initialStateFromConfig(recOpts.phase_bias);
            if (init.estimate == false)
            {
                getBias(trace, time, rec, Sat, code, PHAS, bias, var);

                return true;
            }

            kfKey.type = KF::PHASE_BIAS;

            return kfState.getKFValue(kfKey, bias, &var) != E_Source::NONE;
        }
    }
    else if (rec.empty())
    {
        auto& satOpts = acsConfig.getSatOpts(Sat);

        if (type == CODE)
        {
            if (!satOpts.codeBiasModel.enable)
                return true;

            InitialState init = initialStateFromConfig(satOpts.code_bias);
            if (init.estimate == false)
            {
                getBias(trace, time, Sat.id(), Sat, code, CODE, bias, var);
                return true;
            }

            kfKey.type       = KF::CODE_BIAS;
            E_Source passSrc = kfState.getKFValue(kfKey, bias, &var);
            bool     pass    = passSrc != E_Source::NONE;

            tracepdeex(
                5,
                trace,
                "\n Searching UC %s - %s",
                ((string)kfKey).c_str(),
                pass ? "found" : "not found"
            );

            return pass;
        }

        if (type == PHAS)
        {
            if (satOpts.phaseBiasModel.enable == false)
                return true;

            InitialState init = initialStateFromConfig(satOpts.phase_bias);
            if (init.estimate == false)
            {
                getBias(trace, time, Sat.id(), Sat, code, PHAS, bias, var);
                return true;
            }

            kfKey.type       = KF::PHASE_BIAS;
            E_Source passSrc = kfState.getKFValue(kfKey, bias, &var);
            bool     pass    = passSrc != E_Source::NONE;

            tracepdeex(
                5,
                trace,
                "\n Searching UC %s - %s",
                ((string)kfKey).c_str(),
                pass ? "found" : "not found"
            );

            return pass;
        }
    }

    return false;
}
