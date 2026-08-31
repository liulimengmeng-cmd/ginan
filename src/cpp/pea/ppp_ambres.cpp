// #pragma GCC optimize ("O0")
/**------------------------------------------------------------------------------
 * reference :
 *     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
 *         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
 *         1995
 *     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
 *         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
 *-----------------------------------------------------------------------------*/

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
            const int complementFixCount = GNSS_AR(trace, complementProbe, ARopt);
            tracepdeex(
                2,
                trace,
                "\nPPP_AR INTEGER_COMPLEMENT_DIAGNOSTIC stage1_rows=%d "
                "remaining_coordinates=%d stage2_probe_rows=%d combined_independent_rows=%d "
                "target_integer_rank=%d stage2_status=%s stage2_success_rate=%.17g "
                "stage2_selected_decorrelated=%d stage2_integer_candidates=%d "
                "stage2_best_squared_norm=%.17g stage2_second_squared_norm=%.17g "
                "stage2_ratio=%.17g action=PROBE_ONLY_NOT_SUBMITTED",
                nfix,
                static_cast<int>(complement.ambiguityResolution.aflt.size()),
                complementFixCount,
                nfix + complementFixCount,
                result.integerAmbiguityCoordinateCount,
                complementProbe.diagnosticStatus.c_str(),
                complementProbe.bootstrappedSuccessRate,
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
                "stage2_selected_decorrelated=0 stage2_integer_candidates=0 "
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
