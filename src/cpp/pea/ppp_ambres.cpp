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
#include <cstdlib>
#include <cstdio>
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
#include "pppArJointGate.hpp"

static bool filterError = false;

// Opt-in, read-only event snapshots for exact same-prior constraint diagnostics.
// The directory must already exist. fopen("wx") refuses to overwrite evidence.
static void dumpIntegerConstraintSnapshot(
    Trace& trace, const KFState& state, const MatrixXd& H,
    const VectorXd& v, const MatrixXd& R, const VectorXd& rhs)
{
    const char* directory = std::getenv("GINAN_AR_SNAPSHOT_DIRECTORY");
    const char* times = std::getenv("GINAN_AR_SNAPSHOT_TIMES");
    if (!directory || !times)
        return;
    const string epoch = state.time.to_string(0);
    if (string(times).find("|" + epoch + "|") == string::npos)
        return;
    string name = epoch;
    std::replace(name.begin(), name.end(), ':', '-');
    std::replace(name.begin(), name.end(), ' ', '_');
    const string path = string(directory) + "/" + name + ".snapshot";
    FILE* file = std::fopen(path.c_str(), "wx");
    if (!file)
    {
        tracepdeex(1, trace, "\nPPP_AR SNAPSHOT status=OPEN_FAILED path=%s", path.c_str());
        return;
    }
    std::fprintf(file, "GINAN_AR_SNAPSHOT_V1\nGPST %s\n", epoch.c_str());
    auto matrix = [&](const char* label, const MatrixXd& values)
    {
        std::fprintf(file, "%s %ld %ld\n", label, long(values.rows()), long(values.cols()));
        for (int i = 0; i < values.rows(); i++)
        {
            for (int j = 0; j < values.cols(); j++)
                std::fprintf(file, j ? " %.17g" : "%.17g", values(i,j));
            std::fprintf(file, "\n");
        }
    };
    matrix("x", state.x);
    matrix("P", state.P);
    matrix("H", H);
    matrix("v", v);
    matrix("R", R);
    matrix("rhs", rhs);
    std::fprintf(file, "KEYS %ld\n", long(state.kfIndexMap.size()));
    for (const auto& [key, index] : state.kfIndexMap)
        std::fprintf(file, "%d %s %s %s %d\n", index,
            enum_to_string(key.type).c_str(), key.str.empty() ? "-" : key.str.c_str(),
            key.Sat.id().empty() ? "-" : key.Sat.id().c_str(), key.num);
    const bool success = std::fclose(file) == 0;
    tracepdeex(2, trace, "\nPPP_AR SNAPSHOT status=%s path=%s",
        success ? "WRITTEN_READ_ONLY" : "WRITE_FAILED", path.c_str());
}

static const char* ambiguityFeedbackStateBlock(KF type)
{
    if (type == KF::REC_POS || type == KF::REC_POS_RATE)
    {
        return "COORDINATE";
    }
    if (type == KF::IONO_STEC)
    {
        return "STEC";
    }
    if (type == KF::AMBIGUITY)
    {
        return "AMBIGUITY";
    }
    if (type == KF::TROP || type == KF::TROP_GRAD || type == KF::TROP_MODEL)
    {
        return "TROPOSPHERE";
    }
    if (type > KF::BEGIN_CLOCK_STATES && type < KF::END_CLOCK_STATES)
    {
        return "CLOCK";
    }
    return "OTHER";
}

static map<string, vector<int>> ambiguityFeedbackStateBlocks(const KFState& kfState)
{
    map<string, vector<int>> blocks;
    vector<bool> assigned(kfState.x.size(), false);
    for (const auto& [key, index] : kfState.kfIndexMap)
    {
        if (index < 0 || index >= kfState.x.size() || assigned[index])
        {
            continue;
        }
        assigned[index] = true;
        blocks[ambiguityFeedbackStateBlock(key.type)].push_back(index);
    }
    return blocks;
}

static double indexedTrace(const MatrixXd& matrix, const vector<int>& indices)
{
    double result = 0;
    for (int index : indices)
    {
        result += matrix(index, index);
    }
    return result;
}

static double indexedSquaredNorm(const VectorXd& values, const vector<int>& indices)
{
    double result = 0;
    for (int index : indices)
    {
        result += values(index) * values(index);
    }
    return result;
}

struct PhaseWindupIntegerCanonicalization
{
    VectorXd integerOffsets;
    int nonzeroOffsetCount = 0;
    int missingReceiverCount = 0;
    int missingSatelliteStateCount = 0;
    int nonfinitePhaseWindupCount = 0;
    bool complete = false;
};

static PhaseWindupIntegerCanonicalization buildPhaseWindupIntegerCanonicalization(
    Trace&            trace,
    const GinAR_mtx&  ambiguityResolution
)
{
    PhaseWindupIntegerCanonicalization result;
    const int ambiguityCount = ambiguityResolution.aflt.size();
    result.integerOffsets = VectorXd::Zero(ambiguityCount);

    for (const auto& [column, key] : ambiguityResolution.ambmap)
    {
        if (column < 0 || column >= ambiguityCount)
        {
            result.missingSatelliteStateCount++;
            continue;
        }
        if (key.rec_ptr == nullptr)
        {
            result.missingReceiverCount++;
            continue;
        }
        const auto satelliteState = key.rec_ptr->satStatMap.find(key.Sat);
        if (satelliteState == key.rec_ptr->satStatMap.end())
        {
            result.missingSatelliteStateCount++;
            continue;
        }

        const double continuousPhaseWindup = satelliteState->second.phw;
        if (!std::isfinite(continuousPhaseWindup))
        {
            result.nonfinitePhaseWindupCount++;
            continue;
        }

        // phaseWindup() chooses the continuous branch with
        //   phw = fraction + floor(previous - fraction + 0.5).
        // Removing the nearest integer winding maps every cold start to the
        // same [-0.5, 0.5) representative.  Since the measurement model uses
        // lambda*(phw + N), the corresponding ambiguity coordinate is
        // Nc = Nraw + winding.
        const double integerWinding = std::floor(continuousPhaseWindup + 0.5);
        result.integerOffsets(column) = integerWinding;
        if (integerWinding != 0)
        {
            result.nonzeroOffsetCount++;
        }
        tracepdeex(
            4,
            trace,
            "\nPPP_AR PHASE_WINDUP_CANONICAL_ROW receiver=%s satellite=%s signal=%s "
            "continuous_phw=%.17g integer_winding=%.0f raw_float=%.17g "
            "canonical_float=%.17g status=INTEGER_GAUGE_MAPPED",
            key.str.c_str(),
            key.Sat.id().c_str(),
            key.code().c_str(),
            continuousPhaseWindup,
            integerWinding,
            ambiguityResolution.aflt(column),
            ambiguityResolution.aflt(column) + integerWinding
        );
    }

    result.complete =
        static_cast<int>(ambiguityResolution.ambmap.size()) == ambiguityCount &&
        result.missingReceiverCount == 0 &&
        result.missingSatelliteStateCount == 0 &&
        result.nonfinitePhaseWindupCount == 0;
    const int summaryTraceLevel = result.complete ? 2 : 1;
    tracepdeex(
        summaryTraceLevel,
        trace,
        "\nPPP_AR PHASE_WINDUP_CANONICALIZATION ambiguities=%d nonzero_offsets=%d "
        "missing_receivers=%d missing_satellite_states=%d nonfinite_phase_windup=%d "
        "status=%s",
        ambiguityCount,
        result.nonzeroOffsetCount,
        result.missingReceiverCount,
        result.missingSatelliteStateCount,
        result.nonfinitePhaseWindupCount,
        result.complete ? "COMPLETE_INTEGER_GAUGE" : "INCOMPLETE_NOT_APPLIED"
    );
    return result;
}

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
    const GinAR_mtx&                    originalAmbiguities,
    const char*                         candidateAction
)
{
    if (rowsInOriginalAmbiguities.rows() != fixedIntegers.size() ||
        rowsInOriginalAmbiguities.cols() != originalAmbiguities.aflt.size() ||
        originalAmbiguities.Paflt.rows() != originalAmbiguities.aflt.size() ||
        originalAmbiguities.Paflt.cols() != originalAmbiguities.aflt.size() ||
        rowsInOriginalAmbiguities.cols() !=
            static_cast<int>(originalAmbiguities.ambmap.size()))
    {
        tracepdeex(
            1,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=%s system=%s "
            "candidate_scope=%s status=INVALID_DIMENSIONS action=%s",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            candidateScope,
            candidateAction
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
            const auto ambiguity = originalAmbiguities.ambmap.find(column);
            if (ambiguity == originalAmbiguities.ambmap.end())
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
        const VectorXd candidateRow = rowsInOriginalAmbiguities.row(row);
        const double floatValue = candidateRow.dot(originalAmbiguities.aflt);
        const double variance = (
            candidateRow.transpose() * originalAmbiguities.Paflt * candidateRow
        )(0, 0);
        const double formalSigma = variance >= 0 && std::isfinite(variance)
            ? std::sqrt(variance)
            : -1;
        const double floatMinusInteger = floatValue - fixedIntegers(row);
        tracepdeex(
            traceOutputLevel,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=%s system=%s "
            "reference=%s candidate_scope=%s family=%s row=%d rhs=%.17g "
            "float_value=%.17g float_minus_integer=%.17g formal_sigma=%.17g "
            "support=%d status=%s terms=%saction=%s",
            datum.receiver.c_str(),
            enum_to_string(datum.system).c_str(),
            datum.pivot.id().c_str(),
            candidateScope,
            row < wideLaneCandidateCount ? "WIDE_LANE" : "SECOND_D2",
            row,
            fixedIntegers(row),
            floatValue,
            floatMinusInteger,
            formalSigma,
            support,
            integerCoefficients && support > 0 ? "INTEGER_MAPPED" : "INVALID_MAPPING",
            terms.str().c_str(),
            candidateAction
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

static GinAR_mtx buildAndTraceDualFrequencyDatumCandidate(
    Trace&                                  trace,
    const GinAR_mtx&                        originalAmbiguities,
    const DualFrequencyAmbiguityTransform&  transform,
    const GinAR_opt&                        options
)
{
    GinAR_mtx selectedConstraints;
    selectedConstraints.Ztrs.resize(0, originalAmbiguities.aflt.size());
    selectedConstraints.zfix.resize(0);
    selectedConstraints.ambmap = originalAmbiguities.ambmap;
    const bool feedbackRequested =
        acsConfig.ambrOpts.dual_frequency_subset_feedback;
    const char* candidateAction = feedbackRequested
        ? "CANDIDATE_PENDING_SAFETY_GATES"
        : "PROBE_ONLY_NOT_SUBMITTED";

    tracepdeex(
        2,
        trace,
        "\nPPP_AR DUAL_FREQUENCY_BASIS original=%d rows=%d expected_rank=%d "
        "actual_rank=%d complete_groups=%d incomplete_groups=%d paired_ambiguities=%d "
        "unmatched_ambiguities=%d integer_valued=%d full_row_rank=%d covers_all=%d "
        "status=%s action=%s",
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
        transform.diagnosticStatus.c_str(),
        candidateAction
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
            datum.completeSignalGraph ? candidateAction : "INCOMPLETE_GRAPH_NOT_PROBED"
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
            const MatrixXd candidateRows =
                selectedProbe.combinedTransform * selectedProbe.groupTransform;
            const int previousConstraintCount = selectedConstraints.Ztrs.rows();
            const int addedConstraintCount = candidateRows.rows();
            selectedConstraints.Ztrs.conservativeResize(
                previousConstraintCount + addedConstraintCount,
                originalAmbiguities.aflt.size()
            );
            selectedConstraints.Ztrs.bottomRows(addedConstraintCount) = candidateRows;
            selectedConstraints.zfix.conservativeResize(
                previousConstraintCount + addedConstraintCount
            );
            selectedConstraints.zfix.tail(addedConstraintCount) = combinedFixedIntegers;
            traceDualFrequencyCandidateRows(
                trace,
                datum,
                candidateRows,
                combinedFixedIntegers,
                selectedFamilyCount,
                fullVisibleCandidate ? "FULL_VISIBLE_GROUP" : "SELECTED_SUBSET",
                originalAmbiguities,
                candidateAction
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
                originalAmbiguities,
                candidateAction
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
            "combined_rank=%d target_rank=%d status=%s action=%s",
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
                selectedProbe.diagnosticStatus.c_str(),
            candidateAction
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
        "feedback_requested=%d wrong_fix_certified=0 filter_feedback=0 action=%s",
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
                "NO_FULL_INTEGER_DATUM_CANDIDATE"),
        feedbackRequested,
        candidateAction
    );
    return selectedConstraints;
}

static bool validateDualFrequencyFeedbackPhaseBiasCoverage(
    Trace&          trace,
    const KFState&  kfState,
    const GinAR_mtx& constraints
)
{
    const int columnCount = constraints.Ztrs.cols();
    int usedAmbiguityCount = 0;
    int phaseBiasModelDisabledCount = 0;
    int missingPhaseBiasCount = 0;
    int invalidMapCount = 0;

    for (int column = 0; column < columnCount; column++)
    {
        if (constraints.Ztrs.rows() == 0 ||
            constraints.Ztrs.col(column).cwiseAbs().maxCoeff() <= 1e-12)
        {
            continue;
        }
        usedAmbiguityCount++;
        const auto ambiguity = constraints.ambmap.find(column);
        if (ambiguity == constraints.ambmap.end())
        {
            invalidMapCount++;
            continue;
        }
        const KFKey& key = ambiguity->second;
        auto& satOpts = acsConfig.getSatOpts(key.Sat, {key.code()});
        if (!satOpts.phaseBiasModel.enable)
        {
            phaseBiasModelDisabledCount++;
            continue;
        }

        double bias = 0;
        double variance = 0;
        const E_ObsCode code = int_to_enum<E_ObsCode>(key.num);
        const bool found = getBias(
            trace,
            kfState.time,
            key.Sat.id(),
            key.Sat,
            code,
            PHAS,
            bias,
            variance
        );
        if (!found || !std::isfinite(bias) || !std::isfinite(variance))
        {
            missingPhaseBiasCount++;
        }
    }

    const bool valid =
        usedAmbiguityCount > 0 &&
        invalidMapCount == 0 &&
        phaseBiasModelDisabledCount == 0 &&
        missingPhaseBiasCount == 0;
    const int gateTraceLevel = valid ? 2 : 1;
    tracepdeex(
        gateTraceLevel,
        trace,
        "\nPPP_AR DUAL_FREQUENCY_PHASE_BIAS_GATE used_ambiguities=%d "
        "phase_bias_model_disabled=%d missing_phase_bias=%d invalid_map=%d status=%s",
        usedAmbiguityCount,
        phaseBiasModelDisabledCount,
        missingPhaseBiasCount,
        invalidMapCount,
        valid ? "COMPLETE_PRODUCT_COVERAGE" : "REJECTED_INCOMPLETE_PRODUCT_COVERAGE"
    );
    return valid;
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
    GinAR_mtx& mtrx, ///< Reference to structure containing fixed ambiguities and Z transformations
    const char** rejectionReason = nullptr
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

    const bool traceFeedbackDiagnostics =
        acsConfig.ambrOpts.dual_frequency_feedback_diagnostics;
    dumpIntegerConstraintSnapshot(trace, kfState, kfMeas.H, kfMeas.V, kfMeas.R, zfix);
    // Experimental opt-in: preserve the baseline path unless explicitly enabled.
    // Reuse the filter's sigma convention, but screen the complete AR block before
    // any update or posterior row deweighting. Passing does not certify an integer.
    const char* jointGateOption = std::getenv("GINAN_AR_JOINT_NIS_GATE");
    const bool jointGateEnabled = jointGateOption && string(jointGateOption) == "1";
    if (jointGateEnabled)
    {
        const MatrixXd jointCovariance =
            kfMeas.H * kfState.P * kfMeas.H.transpose() + kfMeas.R;
        const auto gate = pppArJointGate(
            jointCovariance, kfMeas.V, kfState.chiSquareTest.sigma_threshold);
        tracepdeex(
            1, trace,
            "\nPPP_AR JOINT_GATE rows=%d nis=%.17g threshold=%.17g "
            "sigma_threshold=%.17g nominal_alpha=%.17g status=%s action=%s",
            nz, gate.nis, gate.threshold, kfState.chiSquareTest.sigma_threshold,
            gate.alpha, gate.status,
            gate.passed ? "CONTINUE_TO_FILTER_UNVERIFIED" : "KEEP_FLOAT_NOT_SUBMITTED"
        );
        if (!gate.passed)
        {
            if (rejectionReason)
                *rejectionReason = "JOINT_CONSISTENCY_GATE_REJECTED";
            return false;
        }
    }
    const VectorXd stateBefore = traceFeedbackDiagnostics
        ? kfState.x
        : VectorXd();
    const MatrixXd covarianceBefore = traceFeedbackDiagnostics
        ? kfState.P
        : MatrixXd();
    VectorXd shadowStateUpdate;
    MatrixXd shadowCovarianceReduction;
    map<string, vector<int>> stateBlocks;
    map<string, vector<int>> coupledStateBlocks;
    bool shadowValid = false;
    double jointNis = -1;

    if (traceFeedbackDiagnostics)
    {
        const MatrixXd innovationCovariance =
            kfMeas.H * covarianceBefore * kfMeas.H.transpose() + kfMeas.R;
        const LDLT<MatrixXd> innovationSolver(innovationCovariance);
        shadowValid =
            innovationSolver.info() == Eigen::Success &&
            innovationSolver.isPositive() &&
            innovationCovariance.allFinite() &&
            kfMeas.V.allFinite();
        if (shadowValid)
        {
            const MatrixXd stateMeasurementCovariance =
                covarianceBefore * kfMeas.H.transpose();
            const VectorXd solvedInnovation = innovationSolver.solve(kfMeas.V);
            const MatrixXd solvedCrossCovariance = innovationSolver.solve(
                stateMeasurementCovariance.transpose()
            );
            shadowStateUpdate = stateMeasurementCovariance * solvedInnovation;
            shadowCovarianceReduction =
                stateMeasurementCovariance * solvedCrossCovariance;
            shadowValid =
                shadowStateUpdate.allFinite() &&
                shadowCovarianceReduction.allFinite();
            jointNis = shadowValid ? kfMeas.V.dot(solvedInnovation) : -1;
        }

        const int shadowTraceLevel = shadowValid ? 2 : 1;
        tracepdeex(
            shadowTraceLevel,
            trace,
            "\nPPP_AR FEEDBACK_SHADOW_SUMMARY rows=%d innovation_norm=%.17g "
            "joint_nis=%.17g nis_per_row=%.17g nis_gate_configured=%d "
            "status=%s action=ANALYTIC_ONE_STEP_DIAGNOSTIC_ONLY",
            nz,
            kfMeas.V.norm(),
            jointNis,
            shadowValid && nz > 0 ? jointNis / nz : -1,
            jointGateEnabled,
            shadowValid ? "VALID_LINEAR_SHADOW" : "INVALID_INNOVATION_COVARIANCE"
        );

        if (shadowValid)
        {
            stateBlocks = ambiguityFeedbackStateBlocks(kfState);
            for (const auto& [block, indices] : stateBlocks)
            {
                const double priorTrace = indexedTrace(covarianceBefore, indices);
                const double predictedReduction =
                    indexedTrace(shadowCovarianceReduction, indices);
                auto& coupledIndices = coupledStateBlocks[block];
                for (int index : indices)
                {
                    const double scale = std::max(
                        1e-18,
                        std::abs(covarianceBefore(index, index)) * 1e-12
                    );
                    if (std::abs(shadowCovarianceReduction(index, index)) > scale)
                    {
                        coupledIndices.push_back(index);
                    }
                }
                const double coupledPriorTrace =
                    indexedTrace(covarianceBefore, coupledIndices);
                const double coupledPredictedReduction =
                    indexedTrace(shadowCovarianceReduction, coupledIndices);
                tracepdeex(
                    2,
                    trace,
                    "\nPPP_AR FEEDBACK_SHADOW_BLOCK block=%s states=%d "
                    "prior_trace=%.17g predicted_trace_reduction=%.17g "
                    "relative_trace_gain=%.17g predicted_shift_norm=%.17g "
                    "coupled_states=%d coupled_prior_trace=%.17g "
                    "coupled_predicted_trace_reduction=%.17g "
                    "coupled_relative_trace_gain=%.17g "
                    "status=VALID_LINEAR_SHADOW action=DIAGNOSTIC_ONLY",
                    block.c_str(),
                    static_cast<int>(indices.size()),
                    priorTrace,
                    predictedReduction,
                    priorTrace > 0 ? predictedReduction / priorTrace : -1,
                    std::sqrt(indexedSquaredNorm(shadowStateUpdate, indices)),
                    static_cast<int>(coupledIndices.size()),
                    coupledPriorTrace,
                    coupledPredictedReduction,
                    coupledPriorTrace > 0
                        ? coupledPredictedReduction / coupledPriorTrace
                        : -1
                );
            }
        }
    }

    kfState.filterKalman(trace, kfMeas, "/AR", true);

    if (traceFeedbackDiagnostics && shadowValid &&
        kfState.x.size() == stateBefore.size() &&
        kfState.P.rows() == covarianceBefore.rows() &&
        kfState.P.cols() == covarianceBefore.cols())
    {
        const VectorXd realisedStateUpdate = kfState.x - stateBefore;
        const MatrixXd realisedCovarianceReduction = covarianceBefore - kfState.P;
        for (const auto& [block, indices] : stateBlocks)
        {
            const double priorTrace = indexedTrace(covarianceBefore, indices);
            const double realisedReduction =
                indexedTrace(realisedCovarianceReduction, indices);
            const auto& coupledIndices = coupledStateBlocks[block];
            const double coupledPriorTrace =
                indexedTrace(covarianceBefore, coupledIndices);
            const double coupledRealisedReduction =
                indexedTrace(realisedCovarianceReduction, coupledIndices);
            tracepdeex(
                2,
                trace,
                "\nPPP_AR FEEDBACK_REALISED_BLOCK block=%s states=%d "
                "realised_trace_reduction=%.17g relative_trace_gain=%.17g "
                "realised_shift_norm=%.17g coupled_states=%d "
                "coupled_realised_trace_reduction=%.17g "
                "coupled_relative_trace_gain=%.17g status=FILTER_CALL_RETURNED "
                "action=DIAGNOSTIC_ONLY",
                block.c_str(),
                static_cast<int>(indices.size()),
                realisedReduction,
                priorTrace > 0 ? realisedReduction / priorTrace : -1,
                std::sqrt(indexedSquaredNorm(realisedStateUpdate, indices)),
                static_cast<int>(coupledIndices.size()),
                coupledRealisedReduction,
                coupledPriorTrace > 0
                    ? coupledRealisedReduction / coupledPriorTrace
                    : -1
            );
        }
        const MatrixXd shadowPosteriorCovariance =
            covarianceBefore - shadowCovarianceReduction;
        tracepdeex(
            2,
            trace,
            "\nPPP_AR FEEDBACK_SHADOW_REALISATION state_update_error_norm=%.17g "
            "covariance_error_norm=%.17g status=COMPARED_WITH_FILTER_RESULT "
            "action=DIAGNOSTIC_ONLY",
            (realisedStateUpdate - shadowStateUpdate).norm(),
            (kfState.P - shadowPosteriorCovariance).norm()
        );
    }
    else if (traceFeedbackDiagnostics)
    {
        tracepdeex(
            1,
            trace,
            "\nPPP_AR FEEDBACK_SHADOW_REALISATION status=NOT_COMPARABLE "
            "action=DIAGNOSTIC_ONLY"
        );
    }

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
    KFState& kfState, ///< Filter state
    bool     inRts    ///< True when called for a smoothed state
)
{
    AmbiguityResolutionAttempt result;

    tracepdeex(3, trace, "%s: %s\n", __FUNCTION__, kfState.time.to_string().c_str());

    if (acsConfig.ambrOpts.mode == E_ARmode::OFF)
    {
        return result;
    }

    result.routineInvoked = true;

    if (inRts && acsConfig.ambrOpts.canonicalize_phase_windup_integer)
    {
        tracepdeex(
            1,
            trace,
            "\nPPP_AR PHASE_WINDUP_CANONICALIZATION status="
            "RTS_EPOCH_WINDING_UNAVAILABLE action=AR_SKIPPED"
        );
        result.diagnosticStatus = "RTS_PHASE_WINDUP_CANONICALIZATION_UNAVAILABLE";
        return result;
    }

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

    VectorXd phaseWindupIntegerOffsets;
    if (acsConfig.ambrOpts.canonicalize_phase_windup_integer)
    {
        const PhaseWindupIntegerCanonicalization canonicalization =
            buildPhaseWindupIntegerCanonicalization(trace, ARmtx);
        if (!canonicalization.complete ||
            !canonicalizeIntegerAmbiguityFloats(
                ARmtx,
                canonicalization.integerOffsets
            ))
        {
            result.diagnosticStatus = "PHASE_WINDUP_CANONICALIZATION_INVALID";
            return result;
        }
        phaseWindupIntegerOffsets = canonicalization.integerOffsets;
    }

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

    if (acsConfig.ambrOpts.dual_frequency_subset_probe_only &&
        acsConfig.ambrOpts.dual_frequency_subset_feedback)
    {
        tracepdeex(
            1,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CONTROL status=INVALID_PROBE_AND_FEEDBACK_BOTH_ENABLED "
            "legacy_feedback=0 new_subset_feedback=0 action=NOT_SUBMITTED"
        );
        result.diagnosticStatus = "DUAL_FREQUENCY_PROBE_FEEDBACK_CONFIG_CONFLICT";
        return result;
    }

    GinAR_mtx dualFrequencySubsetCandidate;
    if (acsConfig.ambrOpts.integer_complement_diagnostics ||
        acsConfig.ambrOpts.dual_frequency_subset_probe_only ||
        acsConfig.ambrOpts.dual_frequency_subset_feedback)
    {
        const DualFrequencyAmbiguityTransform dualFrequencyTransform =
            buildDualFrequencyAmbiguityIntegerTransform(
                ARmtx,
                acsConfig.receiver_amb_pivot
            );
        dualFrequencySubsetCandidate = buildAndTraceDualFrequencyDatumCandidate(
            trace,
            ARmtx,
            dualFrequencyTransform,
            ARopt
        );
    }

    if (acsConfig.ambrOpts.dual_frequency_subset_probe_only)
    {
        tracepdeex(
            2,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CONTROL legacy_feedback=0 "
            "new_subset_feedback=0 status=FLOAT_STATE_PROBE_ONLY "
            "action=PROBE_ONLY_NOT_SUBMITTED"
        );
        result.diagnosticStatus = "DUAL_FREQUENCY_SUBSET_FLOAT_PROBE_ONLY";
        return result;
    }

    if (acsConfig.ambrOpts.dual_frequency_subset_feedback)
    {
        const int candidateRowCount = dualFrequencySubsetCandidate.Ztrs.rows();
        const int ambiguityCount = ARmtx.aflt.size();
        const bool validDimensions =
            candidateRowCount > 0 &&
            dualFrequencySubsetCandidate.Ztrs.cols() == ambiguityCount &&
            dualFrequencySubsetCandidate.zfix.size() == candidateRowCount &&
            static_cast<int>(dualFrequencySubsetCandidate.ambmap.size()) == ambiguityCount;
        const bool integerDesign = validDimensions &&
            (dualFrequencySubsetCandidate.Ztrs.array() -
             dualFrequencySubsetCandidate.Ztrs.array().round()).abs().maxCoeff() <= 1e-12;
        const bool fullRowRank = validDimensions &&
            dualFrequencySubsetCandidate.Ztrs.fullPivLu().rank() == candidateRowCount;

        if (!acsConfig.ambrOpts.canonicalize_phase_windup_integer ||
            phaseWindupIntegerOffsets.size() != ambiguityCount)
        {
            tracepdeex(
                1,
                trace,
                "\nPPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=%d legacy_feedback=0 "
                "new_subset_feedback=0 status=REJECTED_CANONICAL_WINDUP_REQUIRED "
                "action=NOT_SUBMITTED",
                candidateRowCount
            );
            result.diagnosticStatus = "DUAL_FREQUENCY_CANONICAL_WINDUP_REQUIRED";
            return result;
        }
        if (!validDimensions || !integerDesign || !fullRowRank)
        {
            tracepdeex(
                1,
                trace,
                "\nPPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=%d ambiguities=%d "
                "valid_dimensions=%d integer_design=%d full_row_rank=%d "
                "legacy_feedback=0 new_subset_feedback=0 status=REJECTED_INVALID_CANDIDATE_BLOCK "
                "action=NOT_SUBMITTED",
                candidateRowCount,
                ambiguityCount,
                validDimensions,
                integerDesign,
                fullRowRank
            );
            result.diagnosticStatus = candidateRowCount == 0
                ? "NO_DUAL_FREQUENCY_SUBSET_CANDIDATE"
                : "DUAL_FREQUENCY_SUBSET_CANDIDATE_INVALID";
            return result;
        }
        if (!validateDualFrequencyFeedbackPhaseBiasCoverage(
                trace,
                kfState,
                dualFrequencySubsetCandidate
            ))
        {
            tracepdeex(
                1,
                trace,
                "\nPPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=%d legacy_feedback=0 "
                "new_subset_feedback=0 status=REJECTED_PHASE_BIAS_GATE "
                "action=NOT_SUBMITTED",
                candidateRowCount
            );
            result.diagnosticStatus = "DUAL_FREQUENCY_PHASE_BIAS_GATE_REJECTED";
            return result;
        }

        // Candidate integers were searched in Nc = Nraw + k.  Convert the
        // right-hand sides back to the filter's continuous-windup state gauge.
        dualFrequencySubsetCandidate.zfix -=
            dualFrequencySubsetCandidate.Ztrs * phaseWindupIntegerOffsets;
        tracepdeex(
            2,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_FEEDBACK_MAP rows=%d ambiguities=%d "
            "status=CANONICAL_TO_FILTER_GAUGE",
            candidateRowCount,
            ambiguityCount
        );
        result.resolvedCombinationCount = candidateRowCount;
        result.diagnosticStatus = "DUAL_FREQUENCY_SUBSET_FEEDBACK_READY";
        const char* rejectionReason = "PSEUDOOBS_MODEL_REJECTED";
        result.pseudoObservationsSubmitted =
            applyUCAmbiguities(trace, kfState, dualFrequencySubsetCandidate, &rejectionReason);
        const int feedbackTraceLevel = result.pseudoObservationsSubmitted ? 2 : 1;
        tracepdeex(
            feedbackTraceLevel,
            trace,
            "\nPPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=%d legacy_feedback=0 "
            "new_subset_feedback=%d status=%s action=%s",
            candidateRowCount,
            result.pseudoObservationsSubmitted,
            result.pseudoObservationsSubmitted
                ? "SUBMITTED_UNVERIFIED"
                : rejectionReason,
            result.pseudoObservationsSubmitted
                ? "FILTER_CALL_RETURNED"
                : "NOT_SUBMITTED"
        );
        if (!result.pseudoObservationsSubmitted)
        {
            result.diagnosticStatus = rejectionReason;
        }
        return result;
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
                ARmtx.ambmap,
                phaseWindupIntegerOffsets
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
