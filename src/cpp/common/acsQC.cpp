// #pragma GCC optimize ("O0")

#include "common/acsQC.hpp"
#include "common/zhangCycleSlipEvidence.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/enums.h"
#include "common/ephemeris.hpp"
#include "common/navigation.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/rinex.hpp"
#include "common/satStat.hpp"
#include "common/trace.hpp"
#include "rtklib/lambda.h"
#include "orbprop/coordinates.hpp"

#define THRES_MW_JUMP 10.0
#define PDEGAP 60.0
#define PDESLIPTHRESHOLD 0.5
#define PROC_NOISE_IONO 0.001

enum class E_SlipDiagReason
{
    NONE,
    RECEIVER_LLI,
    NO_FREQUENCY_PAIR,
    NO_FREQUENCIES,
    INVALID_LC,
    NO_PREVIOUS_GF,
    NO_PREVIOUS_MW,
    PHASE_JUMP,
    MW_JUMP,
    WITHIN_THRESHOLD,
    FIRST_EPOCH,
    TIME_GAP,
    LOW_ELEVATION,
    DUAL_FREQUENCY,
    TRIPLE_FREQUENCY,
    THIRD_FREQUENCY_REACQUIRED,
    FREQUENCY_REACQUIRED,
    SINGLE_FREQUENCY,
    MISSING_FREQUENCY_VARIANCE,
    INVALID_CODE_VARIANCE,
    INVALID_PHASE_VARIANCE,
    VALID_VARIANCE,
    LOM_WITHIN_THRESHOLD,
    LOM_OUTLIER,
    LAMBDA_FIXED,
    LAMBDA_FLOAT,
    KALMAN_WAITING
};

/** Emit one structured cycle-slip diagnostic trace line.
 *
 * These lines are intended to make detector decisions auditable without changing
 * the detector result.  They are emitted at trace level 2 with stable
 * key=value fields so trace files can be summarised with simple text tools.
 */
void traceSlipEvent(
    Trace&           trace,
    const char*      detector,
    const GObs&      obs,
    const char*      action,
    E_FType          frq1,
    E_FType          frq2,
    E_FType          frq3,
    double           value1,
    double           value2,
    double           threshold,
    E_SlipDiagReason reason
)
{
    tracepdeex(
        2,
        trace,
        "\nPDE-CS-DIAG detector=%s action=%s epoch=%s rec=%s sat=%s f1=%s f2=%s f3=%s "
        "value1=%.6f value2=%.6f threshold=%.6f reason=%s",
        detector,
        action,
        obs.time.to_string(2).c_str(),
        obs.mount.c_str(),
        obs.Sat.id().c_str(),
        enum_to_string(frq1).c_str(),
        enum_to_string(frq2).c_str(),
        enum_to_string(frq3).c_str(),
        value1,
        value2,
        threshold,
        enum_to_lowerstring(reason).c_str()
    );
}

/** Emit a structured diagnostic line for SCDIA internals.
 *
 * `PDE-SCDIA-DIAG` is intentionally separate from `PDE-CS-DIAG` so existing
 * detector-routing summaries remain stable while SCDIA outcomes become
 * machine-readable.
 */
void traceScdiaEvent(
    Trace&           trace,
    const GObs&      obs,
    const char*      action,
    E_FType          frq1,
    E_FType          frq2,
    E_FType          frq3,
    int              nf,
    E_FilterMode     filterMode,
    double           value1,
    double           value2,
    double           threshold,
    double           amb1,
    double           amb2,
    double           amb3,
    E_SlipDiagReason reason
)
{
    tracepdeex(
        2,
        trace,
        "\nPDE-SCDIA-DIAG action=%s epoch=%s rec=%s sat=%s nf=%d mode=%s f1=%s f2=%s f3=%s "
        "value1=%.6f value2=%.6f threshold=%.6f amb1=%.6f amb2=%.6f amb3=%.6f reason=%s",
        action,
        obs.time.to_string(2).c_str(),
        obs.mount.c_str(),
        obs.Sat.id().c_str(),
        nf,
        enum_to_string(filterMode).c_str(),
        enum_to_string(frq1).c_str(),
        enum_to_string(frq2).c_str(),
        enum_to_string(frq3).c_str(),
        value1,
        value2,
        threshold,
        amb1,
        amb2,
        amb3,
        enum_to_lowerstring(reason).c_str()
    );
}

double lomThreshold(int dof)
{
    const double chisqr_arr[100] = {
        10.8, 13.8, 16.3, 18.5, 20.5, 22.5, 24.3, 26.1, 27.9, 29.6, 31.3, 32.9, 34.5, 36.1, 37.7,
        39.3, 40.8, 42.3, 43.8, 45.3, 46.8, 48.3, 49.7, 51.2, 52.6, 54.1, 55.5, 56.9, 58.3, 59.7,
        61.1, 62.5, 63.9, 65.2, 66.6, 68.0, 69.3, 70.7, 72.1, 73.4, 74.7, 76.0, 77.3, 78.6, 80.0,
        81.3, 82.6, 84.0, 85.4, 86.7, 88.0, 89.3, 90.6, 91.9, 93.3, 94.7, 96.0, 97.4, 98.7, 100,
        101,  102,  103,  104,  105,  107,  108,  109,  110,  112,  113,  114,  115,  116,  118,
        119,  120,  122,  123,  125,  126,  127,  128,  129,  131,  132,  133,  134,  135,  137,
        138,  139,  140,  142,  143,  144,  145,  147,  148,  149
    };

    if (dof <= 0 || dof > 100)
        return 0;

    return chisqr_arr[dof - 1] / dof;
}

/** Select up to three configured frequency bands for a satellite system.
 *
 * This is the legacy frequency selector used by several modelling paths outside
 * slip detection.  It follows `acsConfig.code_priorities[sys]` and converts the
 * first unique codes into frequency bands, but it does not inspect a particular
 * observation to confirm that measurements are present.
 *
 * The outputs retain historical defaults (`F1`, `F2`, `F5`) when priorities are
 * incomplete.  Do not use this helper when the algorithm must know which
 * frequencies are actually observed at the current epoch; use `obsFreqs()` for
 * slip detection.
 *
 * @param[in]  sys Satellite system.
 * @param[out] ft1 First configured frequency band.
 * @param[out] ft2 Second configured frequency band.
 * @param[out] ft3 Third configured frequency band.
 *
 * @return false only when no code priorities exist for `sys`; true otherwise.
 */
bool satFreqs(E_Sys sys, E_FType& ft1, E_FType& ft2, E_FType& ft3)
{
    bool ft1Ready = false;
    bool ft2Ready = false;

    ft1 = F1;
    ft2 = F2;
    ft3 = F5;

    if (acsConfig.code_priorities.find(sys) == acsConfig.code_priorities.end())
        return false;

    for (auto& code : acsConfig.code_priorities[sys])
    {
        E_FType ft = code2Freq[sys][code];

        if (ft1Ready == false)
        {
            ft1      = ft;
            ft1Ready = true;
            continue;
        }

        if (ft == ft1)
            continue;

        if (ft2Ready == false)
        {
            ft2      = ft;
            ft2Ready = true;
            continue;
        }

        if (ft == ft2)
            continue;

        {
            ft3 = ft;
            break;
        }
    }

    return true;
}

/** Select observed frequency bands that are usable for slip detection.
 *
 * This helper is intentionally stricter than `satFreqs()`: it still honours
 * configured code priorities, but only returns frequency bands that are present
 * in the current observation, have a non-zero representative phase measurement,
 * and have a non-zero wavelength available in the satellite navigation data.
 *
 * This prevents GF/MW/PDE slip checks from evaluating assumed frequency bands
 * when an epoch is missing one of the configured signals.
 *
 * @param[in]  obs Observation to inspect.
 * @param[out] ft1 First observed usable frequency, or `NONE`.
 * @param[out] ft2 Second observed usable frequency, or `NONE`.
 * @param[out] ft3 Third observed usable frequency, or `NONE`.
 *
 * @return Number of usable observed frequencies written to the output
 *         parameters, from 0 to 3.
 */
int obsFreqs(const GObs& obs, E_FType& ft1, E_FType& ft2, E_FType& ft3)
{
    ft1 = NONE;
    ft2 = NONE;
    ft3 = NONE;

    E_Sys sys = obs.Sat.sys;
    if (acsConfig.code_priorities.find(sys) == acsConfig.code_priorities.end())
        return 0;

    if (obs.satNav_ptr == nullptr)
        return 0;

    auto sysCodeIt = code2Freq.find(sys);
    if (sysCodeIt == code2Freq.end())
        return 0;

    int count = 0;

    for (auto& code : acsConfig.code_priorities[sys])
    {
        auto codeIt = sysCodeIt->second.find(code);
        if (codeIt == sysCodeIt->second.end())
            continue;

        E_FType ft = codeIt->second;
        if (ft == NONE || ft == ft1 || ft == ft2 || ft == ft3)
            continue;

        auto sigIt = obs.sigs.find(ft);
        if (sigIt == obs.sigs.end() || sigIt->second.L == 0)
            continue;

        auto lamIt = obs.satNav_ptr->lamMap.find(ft);
        if (lamIt == obs.satNav_ptr->lamMap.end() || lamIt->second == 0)
            continue;

        if (count == 0)
        {
            ft1 = ft;
            count++;
            continue;
        }
        if (count == 1)
        {
            ft2 = ft;
            count++;
            continue;
        }
        {
            ft3 = ft;
            count++;
            break;
        }
    }

    return count;
}

struct SlipNoise
{
    double           sigmaCode  = 0;
    double           sigmaPhase = 0;
    E_SlipDiagReason reason     = E_SlipDiagReason::NONE;
};

/** Derive conservative slip-detector noise from the selected observation bands.
 *
 * Slip detection chooses concrete observed frequencies with `obsFreqs()`.  The
 * corresponding noise values must come from those same bands; using an
 * arbitrary first signal can pass zero or invalid variances into PDE/SCDIA and
 * produce NaNs in trace output.
 *
 * This helper requires each selected band to have finite, positive code and
 * phase variance.  It then uses the largest selected variance for each
 * observable type as a conservative single-noise input to the current PDE/SCDIA
 * equations.
 *
 * @param[in]  obs   Observation containing selected signals.
 * @param[in]  freqs Selected frequency bands.
 * @param[in]  nf    Number of selected frequencies to validate.
 * @param[out] noise Derived standard deviations, or failure reason.
 *
 * @return true when all selected bands have usable variances.
 */
bool slipNoise(const GObs& obs, const E_FType freqs[], int nf, SlipNoise& noise)
{
    double maxCodeVar  = 0;
    double maxPhaseVar = 0;

    for (int i = 0; i < nf; i++)
    {
        auto sigIt = obs.sigs.find(freqs[i]);
        if (sigIt == obs.sigs.end())
        {
            noise.reason = E_SlipDiagReason::MISSING_FREQUENCY_VARIANCE;
            return false;
        }

        double codeVar = sigIt->second.codeVar;
        double phasVar = sigIt->second.phasVar;

        if (!std::isfinite(codeVar) || codeVar <= 0)
        {
            noise.reason = E_SlipDiagReason::INVALID_CODE_VARIANCE;
            return false;
        }

        if (!std::isfinite(phasVar) || phasVar <= 0)
        {
            noise.reason = E_SlipDiagReason::INVALID_PHASE_VARIANCE;
            return false;
        }

        if (codeVar > maxCodeVar)
            maxCodeVar = codeVar;

        if (phasVar > maxPhaseVar)
            maxPhaseVar = phasVar;
    }

    noise.sigmaCode  = sqrt(maxCodeVar);
    noise.sigmaPhase = sqrt(maxPhaseVar);
    noise.reason     = E_SlipDiagReason::VALID_VARIANCE;

    return true;
}

/** Detect cycle slip by reported loss of lock
 */
void detslp_ll(
    Trace&   trace,   ///< Trace to output to
    ObsList& obsList  ///< List of observations to detect slips within
)
{
    if (obsList.empty())
    {
        tracepdeex(3, trace, "\n%s: epoch=? n=%zu (empty obsList)", __FUNCTION__, obsList.size());
        return;
    }

    // Find first non-null element for the timestamp
    std::string epoch = "?";
    for (const auto& sp : obsList)
    {
        if (sp)
        {
            epoch = sp->time.to_string(2);
            break;
        }
    }

    tracepdeex(3, trace, "\n%s: epoch=%s n=%zu", __FUNCTION__, epoch.c_str(), obsList.size());

    // 	auto begin_iter = boost::make_filter_iterator([]

    for (auto& obs : only<GObs>(obsList))
        for (auto& [ft, sig] : obs.sigs)
        {
            if (obs.exclude)
            {
                continue;
            }

            if (acsConfig.zhangFullRank.enable)
            {
                auto& status=obs.satStat_ptr->sigStatMap[ft2string(ft)];
                const auto decision=zhangDecidePhaseTracking(status.tracking,
                    sig.LLI, sig.L!=0, int(sig.code), double(obs.time.bigTime),
                    1.5*acsConfig.epoch_interval);
                status.phaseQuarantined = !decision.useCurrentPhase;
                if (decision.breakArc)
                {
                    if(zhangEvidenceSlipQcEnabled())
                        status.phaseContinuityUnproven=decision.gap || (sig.LLI&2) || !(sig.LLI&1);
                    if (decision.gap) {status.slip.retrack=true;status.savedSlip.retrack=true;}
                    else {status.slip.LLI=true;status.savedSlip.LLI=true;}
                }
                if (decision.breakArc || !decision.halfCycleResolved)
                    trace << "\nZHANG_PHASE_TRACKING time=" << obs.time.to_string(0)
                          << " receiver=" << obs.mount << " satellite=" << obs.Sat.id()
                          << " code=" << enum_to_string(sig.code) << " lli=" << int(sig.LLI)
                          << " break_arc=" << decision.breakArc
                          << " use_phase=" << decision.useCurrentPhase << " use_code=1"
                          << " half_cycle_resolved=" << decision.halfCycleResolved
                          << " gap_continuity_unproven=" << decision.gap
                          << " event_sequence=" << status.tracking.eventSequence;
                continue;
            }
            if (sig.L == 0 || (sig.LLI & 0x01) == 0)
            {
                continue;
            }

            traceSlipEvent(
                trace,
                "LLI",
                obs,
                "detected",
                ft,
                NONE,
                NONE,
                1,
                0,
                0,
                E_SlipDiagReason::RECEIVER_LLI
            );

            obs.satStat_ptr->sigStatMap[ft2string(ft)].slip.LLI      = true;
            obs.satStat_ptr->sigStatMap[ft2string(ft)].savedSlip.LLI = true;
        }
}

/** Detect cycle slip by geometry free phase jump
 */
void detslp_gf(
    Trace&   trace,   ///< Trace to output to
    ObsList& obsList  ///< List of observations to detect slips within
)
{
    if (obsList.empty())
    {
        tracepdeex(3, trace, "\n%s: epoch=? n=%zu (empty obsList)", __FUNCTION__, obsList.size());
        return;
    }

    // Find first non-null element for the timestamp
    std::string epoch = "?";
    for (const auto& sp : obsList)
    {
        if (sp)
        {
            epoch = sp->time.to_string(2);
            break;
        }
    }

    tracepdeex(3, trace, "\n%s: epoch=%s n=%zu", __FUNCTION__, epoch.c_str(), obsList.size());

    for (auto& obs : only<GObs>(obsList))
    {
        if (obs.exclude)
        {
            continue;
        }

        E_FType frq1;
        E_FType frq2;
        E_FType frq3;
        int     nf = obsFreqs(obs, frq1, frq2, frq3);
        if (nf < 2)
        {
            traceSlipEvent(
                trace,
                "GF",
                obs,
                "skipped",
                NONE,
                NONE,
                NONE,
                0,
                0,
                0,
                E_SlipDiagReason::NO_FREQUENCY_PAIR
            );
            continue;
        }

        S_LC& lc = getLC(obs.satStat_ptr->lc_new, frq1, frq2);

        double gf1 = lc.GF_Phas_m;
        if (lc.valid == false || gf1 == 0)
        {
            traceSlipEvent(
                trace,
                "GF",
                obs,
                "skipped",
                frq1,
                frq2,
                NONE,
                gf1,
                0,
                0,
                E_SlipDiagReason::INVALID_LC
            );
            continue;
        }

        if (acsConfig.zhangFullRank.enable)
        {
            auto& status=*obs.satStat_ptr;
            auto& first=status.sigStatMap[ft2string(frq1)];
            auto& second=status.sigStatMap[ft2string(frq2)];
            if (first.tracking.halfCycle || second.tracking.halfCycle) continue;
            const auto& opts=acsConfig.preprocOpts;
            const double variance=obs.sigs.at(frq1).phasVar+obs.sigs.at(frq2).phasVar;
            const auto decision=zhangInspectCombination(status.acceptedGf,
                double(obs.time.bigTime),gf1,variance,opts.zhang_gf_variance_rate,
                opts.zhang_gf_sigma_gate,opts.zhang_detector_max_gap,
                first.slip.LLI || second.slip.LLI || first.slip.retrack || second.slip.retrack);
            first.phaseQuarantined |= !decision.useSample;
            second.phaseQuarantined |= !decision.useSample;
            if (decision.breakArc)
            {
                first.slip.GF=second.slip.GF=true;
                first.savedSlip.GF=second.savedSlip.GF=true;
                first.tracking.pendingBreak=second.tracking.pendingBreak=true;
            }
            if (decision.useSample) status.gf=gf1;
            trace << "\nZHANG_GF_PREDICTION time=" << obs.time.to_string(0)
                  << " receiver=" << obs.mount << " satellite=" << obs.Sat.id()
                  << " value_m=" << gf1 << " prediction_m=" << decision.prediction
                  << " innovation_m=" << decision.innovation << " variance_m2=" << decision.variance
                  << " statistic=" << decision.statistic << " sigma_gate=" << opts.zhang_gf_sigma_gate
                  << " suspect=" << decision.suspect << " break_arc=" << decision.breakArc
                  << " accepted=" << decision.useSample << " cross_frequency_noise_assumed_independent=1";
            continue;
        }
        double gf0          = obs.satStat_ptr->gf;
        obs.satStat_ptr->gf = gf1;

        if (gf0 == 0)
        {
            traceSlipEvent(
                trace,
                "GF",
                obs,
                "initialised",
                frq1,
                frq2,
                NONE,
                gf1,
                gf0,
                0,
                E_SlipDiagReason::NO_PREVIOUS_GF
            );
            continue;
        }

        if (fabs(gf1 - gf0) > acsConfig.preprocOpts.slip_threshold)
        {
            obs.satStat_ptr->sigStatMap[ft2string(frq1)].slip.GF      = true;
            obs.satStat_ptr->sigStatMap[ft2string(frq2)].slip.GF      = true;
            obs.satStat_ptr->sigStatMap[ft2string(frq1)].savedSlip.GF = true;
            obs.satStat_ptr->sigStatMap[ft2string(frq2)].savedSlip.GF = true;
            traceSlipEvent(
                trace,
                "GF",
                obs,
                "detected",
                frq1,
                frq2,
                NONE,
                gf1,
                gf0,
                acsConfig.preprocOpts.slip_threshold,
                E_SlipDiagReason::PHASE_JUMP
            );
        }
        else
        {
            traceSlipEvent(
                trace,
                "GF",
                obs,
                "accepted",
                frq1,
                frq2,
                NONE,
                gf1,
                gf0,
                acsConfig.preprocOpts.slip_threshold,
                E_SlipDiagReason::WITHIN_THRESHOLD
            );
        }
    }
}

/** Detect slip by Melbourne-Wubbena linear combination jump
 */
void detslp_mw(
    Trace&   trace,   ///< Trace to output to
    ObsList& obsList  ///< List of observations to detect slips within
)
{
    if (obsList.empty())
    {
        tracepdeex(3, trace, "\n%s: epoch=? n=%zu (empty obsList)", __FUNCTION__, obsList.size());
        return;
    }

    // Find first non-null element for the timestamp
    std::string epoch = "?";
    for (const auto& sp : obsList)
    {
        if (sp)
        {
            epoch = sp->time.to_string(2);
            break;
        }
    }

    tracepdeex(3, trace, "\n%s: epoch=%s n=%zu", __FUNCTION__, epoch.c_str(), obsList.size());

    for (auto& obs : only<GObs>(obsList))
    {
        if (obs.exclude)
        {
            continue;
        }

        E_FType frq1;
        E_FType frq2;
        E_FType frq3;
        int     nf = obsFreqs(obs, frq1, frq2, frq3);
        if (nf < 2)
        {
            traceSlipEvent(
                trace,
                "MW",
                obs,
                "skipped",
                NONE,
                NONE,
                NONE,
                0,
                0,
                0,
                E_SlipDiagReason::NO_FREQUENCY_PAIR
            );
            continue;
        }

        S_LC& lc = getLC(obs.satStat_ptr->lc_new, frq1, frq2);

        double mw1 = lc.MW_c;
        if (lc.valid == false || mw1 == 0)
        {
            traceSlipEvent(
                trace,
                "MW",
                obs,
                "skipped",
                frq1,
                frq2,
                NONE,
                mw1,
                0,
                0,
                E_SlipDiagReason::INVALID_LC
            );
            continue;
        }

        if (acsConfig.zhangFullRank.enable)
        {
            auto& status=*obs.satStat_ptr;
            auto& first=status.sigStatMap[ft2string(frq1)];
            auto& second=status.sigStatMap[ft2string(frq2)];
            if (first.tracking.halfCycle || second.tracking.halfCycle) continue;
            const auto& a=obs.sigs.at(frq1);const auto& b=obs.sigs.at(frq2);
            const double sum=lc.lam_A+lc.lam_B;
            if (lc.lam_A<=0 || lc.lam_B<=0 || lc.lam_WL==0 || sum==0) continue;
            const double variance=a.phasVar/SQR(lc.lam_A)+b.phasVar/SQR(lc.lam_B)+
                (SQR(lc.lam_B/sum)*a.codeVar+SQR(lc.lam_A/sum)*b.codeVar)/SQR(lc.lam_WL);
            const auto& opts=acsConfig.preprocOpts;
            const auto decision=zhangInspectCombination(status.acceptedMw,
                double(obs.time.bigTime),mw1,variance,opts.zhang_mw_variance_rate,
                opts.zhang_mw_sigma_gate,opts.zhang_detector_max_gap,
                first.slip.LLI || second.slip.LLI || first.slip.retrack || second.slip.retrack);
            // MW alone cannot identify whether code or phase caused a spike.
            // Quarantine the suspect samples, without retiring a continuous arc.
            first.phaseQuarantined |= !decision.useSample;
            second.phaseQuarantined |= !decision.useSample;
            first.codeQuarantined |= !decision.useSample;
            second.codeQuarantined |= !decision.useSample;
            if (decision.breakArc)
            {
                first.slip.MW=second.slip.MW=true;
                first.savedSlip.MW=second.savedSlip.MW=true;
                first.tracking.pendingBreak=second.tracking.pendingBreak=true;
            }
            if (decision.useSample) status.mw=mw1;
            trace << "\nZHANG_MW_PREDICTION time=" << obs.time.to_string(0)
                  << " receiver=" << obs.mount << " satellite=" << obs.Sat.id()
                  << " innovation_cycles=" << decision.innovation << " variance_cycles2=" << decision.variance
                  << " statistic=" << decision.statistic << " suspect=" << decision.suspect
                  << " break_arc=" << decision.breakArc << " accepted=" << decision.useSample
                  << " source_identified=0 gf_and_mw_and_gate=0";
            continue;
        }
        double mw0          = obs.satStat_ptr->mw;
        obs.satStat_ptr->mw = mw1;

        if (mw0 == 0)
        {
            traceSlipEvent(
                trace,
                "MW",
                obs,
                "initialised",
                frq1,
                frq2,
                NONE,
                mw1,
                mw0,
                0,
                E_SlipDiagReason::NO_PREVIOUS_MW
            );
            continue;
        }

        if (fabs(mw1 - mw0) > THRES_MW_JUMP)
        {
            obs.satStat_ptr->sigStatMap[ft2string(frq1)].slip.MW      = true;
            obs.satStat_ptr->sigStatMap[ft2string(frq2)].slip.MW      = true;
            obs.satStat_ptr->sigStatMap[ft2string(frq1)].savedSlip.MW = true;
            obs.satStat_ptr->sigStatMap[ft2string(frq2)].savedSlip.MW = true;
            traceSlipEvent(
                trace,
                "MW",
                obs,
                "detected",
                frq1,
                frq2,
                NONE,
                mw1,
                mw0,
                THRES_MW_JUMP,
                E_SlipDiagReason::MW_JUMP
            );
        }
        else
        {
            traceSlipEvent(
                trace,
                "MW",
                obs,
                "accepted",
                frq1,
                frq2,
                NONE,
                mw1,
                mw0,
                THRES_MW_JUMP,
                E_SlipDiagReason::WITHIN_THRESHOLD
            );
        }
    }
}

/** Melbourne-Wenbunna (MW) measurement noise (m)
 */
double mwnoise(
    double sigcode,   ///< Code noise
    double sigphase,  ///< Phase noise
    double lam1,      ///< L1 wavelength
    double lam2       ///< L2 wavelength
)
{
    double a =
        lam2 * lam2 / (lam2 + lam1) / (lam2 + lam1) + lam1 * lam1 / (lam2 + lam1) / (lam2 + lam1);
    double b =
        lam2 * lam2 / (lam2 - lam1) / (lam2 - lam1) + lam1 * lam1 / (lam2 - lam1) / (lam2 - lam1);
    return SQRT(a * SQR(sigcode) + b * SQR(sigphase));
}

/** Single channel detection–identification–adaptation (DIA) for integer cycle slips
 */
void scdia(
    Trace&            trace,       ///< Trace to output to
    SatStat&          satStat,     ///< Persistant satellite status parameters
    lc_t&             lc,          ///< Linear combinations
    const GObs&       obs,         ///< Observation context for diagnostics
    map<int, double>& lam,         ///< Signal wavelength map
    double            sigmaPhase,  ///< Phase noise
    double            sigmaCode,   ///< Code noise
    int               nf,          ///< Number of frequencies
    E_FilterMode      filterMode,  ///< LSQ/Kalman filter flag
    E_FType           frq1,
    E_FType           frq2,
    E_FType           frq3
)
{
    if (nf == 0)
        return;

    lc_t* lc_pre_ptr;

    if (filterMode == E_FilterMode::LSQ)
        lc_pre_ptr = &satStat.lc_pre;
    else
        lc_pre_ptr = &satStat.flt.lc_pre;
    if (nf == 1)
        lc_pre_ptr = &satStat.flt.lc_pre;

    auto& lc_pre = *lc_pre_ptr;

    /* single frequency not supported in current PDE */
    if (nf == 1)
    {
        return;
    }

    E_FType freqs[3] = {frq1, frq2, frq3};

    /* m-rows measurements, n-cols unknowns */
    int      m = 2 * nf + 1;
    int      n = 2 + nf;
    VectorXd Z = VectorXd::Zero(m);
    MatrixXd R = MatrixXd::Identity(m, m);
    MatrixXd H = MatrixXd::Zero(m, n);

    double lam1 = lam[frq1];
    int    i    = 0;

    // phase and code
    for (int f = 0; f < nf; f++)
    {
        E_FType frqX = freqs[f];
        double  lamX = lam[frqX];

        Z[i]        = lc.L_m[frqX] - lc_pre.L_m[frqX];
        R(i, i)     = 1 / (2 * SQR(sigmaPhase));
        H(i, 0)     = 1;
        H(i, 1)     = -SQR(lamX) / SQR(lam1);
        H(i, 2 + f) = lamX;
        i++;

        Z[i]    = lc.P[frqX] - lc_pre.P[frqX];
        R(i, i) = 1 / (2 * SQR(sigmaCode));
        H(i, 0) = 1;
        H(i, 1) = +SQR(lamX) / SQR(lam1);
        i++;
    }

    // ionosphere
    {
        Z[i]    = satStat.dIono;
        R(i, i) = 1 / SQR(satStat.sigmaIono);
        H(i, 1) = 1;
        i++;
    }

    /* perform LOM test for outlier detection */
    /* design matrix for LOM test */
    MatrixXd Hlom  = H.leftCols(2);
    VectorXd v     = VectorXd::Zero(m);
    int      ind   = lsqqc(trace, Hlom.data(), R.data(), Z.data(), v.data(), m, 2, 0, 0);
    double   vtpv  = v.dot(R * v);
    int      dof   = m - 2;
    double   val   = dof > 0 ? vtpv / dof : 0;
    double   thres = lomThreshold(dof);
    if (ind == 0)
    {
        traceScdiaEvent(
            trace,
            obs,
            "accepted",
            frq1,
            frq2,
            frq3,
            nf,
            filterMode,
            vtpv,
            val,
            thres,
            0,
            0,
            0,
            E_SlipDiagReason::LOM_WITHIN_THRESHOLD
        );
        return;
    }

    traceScdiaEvent(
        trace,
        obs,
        "detected",
        frq1,
        frq2,
        frq3,
        nf,
        filterMode,
        vtpv,
        val,
        thres,
        0,
        0,
        0,
        E_SlipDiagReason::LOM_OUTLIER
    );

    satStat.sigStatMap[ft2string(frq1)].slip.SCDIA = true;
    satStat.sigStatMap[ft2string(frq2)].slip.SCDIA = true;
    if (nf == 3)
        satStat.sigStatMap[ft2string(frq3)].slip.SCDIA = true;
    satStat.sigStatMap[ft2string(frq1)].savedSlip.SCDIA = true;
    satStat.sigStatMap[ft2string(frq2)].savedSlip.SCDIA = true;
    if (nf == 3)
        satStat.sigStatMap[ft2string(frq3)].savedSlip.SCDIA = true;

    VectorXd xp = VectorXd::Zero(n);
    MatrixXd Pp = MatrixXd::Zero(n, n);

    if (filterMode == E_FilterMode::LSQ)
    {
        MatrixXd N  = MatrixXd::Zero(n, m);
        VectorXd N1 = VectorXd::Zero(n);
        matmul("TN", n, m, m, 1, H.data(), R.data(), 0, N.data());  /* H'*R */
        matmul("NN", n, n, m, 1, N.data(), H.data(), 0, Pp.data()); /* H'*R*H */
        matmul("NN", n, 1, m, 1, N.data(), Z.data(), 0, N1.data()); /* Nl=H'*R*Z */
        if (!matinv(Pp.data(), n))
        {
            matmul("NN", n, 1, n, 1, Pp.data(), N1.data(), 0, xp.data());
        }
        /* store float solution and vc matrix */
        matcpy(satStat.flt.a, xp.data() + 2, 1, nf);

        for (int i = 0; i < nf; i++)
            for (int j = 0; j < nf; j++)
                satStat.flt.Qa[i][j] = Pp.data()[(i + 2) * n + j + 2];
    }
    else
    {
        satStat.flt.ne++;
        if (satStat.flt.ne < 2)
        {
            satStat.flt.slip = 0;
            satStat.flt.ne   = 0;

            traceScdiaEvent(
                trace,
                obs,
                "waiting",
                frq1,
                frq2,
                frq3,
                nf,
                filterMode,
                satStat.flt.ne,
                2,
                0,
                0,
                0,
                0,
                E_SlipDiagReason::KALMAN_WAITING
            );
            return;
        }

        VectorXd x = VectorXd::Zero(n);
        matcpy(x.data() + 2, satStat.flt.a, 1, nf);

        /* time update */
        MatrixXd Px = MatrixXd::Zero(n, n);
        for (int i = 0; i < nf; i++)
            for (int j = 0; j < nf; j++)
                Px.data()[(i + 2) * n + j + 2] = satStat.flt.Qa[i][j];

        Px.data()[0]     = 1E6;
        Px.data()[1 + n] = 1E6;

        /* measurement-prediction */
        matmul("NN", m, 1, n, -1, H.data(), x.data(), 1, Z.data());

        /* transpose of desgin matrix */
        MatrixXd I  = MatrixXd::Identity(m, m);
        MatrixXd H1 = MatrixXd::Zero(n, m);
        matmul("TN", n, m, m, +1, H.data(), I.data(), 0, H1.data());

        /* measurement update */
        if (!matinv(R.data(), m))
            filter_(x.data(), Px.data(), H1.data(), Z.data(), R.data(), n, m, xp.data(), Pp.data());

        matcpy(satStat.flt.a, xp.data() + 2, 1, nf);

        for (int i = 0; i < nf; i++)
            for (int j = 0; j < nf; j++)
                satStat.flt.Qa[i][j] = Pp.data()[(i + 2) * n + j + 2];
    }

    /* ambiguity vector and its variance */
    VectorXd a = VectorXd::Zero(nf);
    matcpy(a.data(), xp.data() + n - nf, nf, 1);

    MatrixXd Qa = MatrixXd::Zero(nf, nf);
    for (int i = 0; i < nf; i++)
        for (int j = 0; j < nf; j++)
        {
            Qa.data()[i * nf + j] = Pp.data()[(n - nf + i) * n + j + n - nf];
        }

    /* integer cycle slip estimation */
    MatrixXd F = MatrixXd::Zero(nf, 2);
    double   s[2];
    bool     pass = false;
    lambda(trace, nf, 2, a.data(), Qa.data(), F.data(), s, acsConfig.predefined_fail, pass);

    double ratio = s[1] != 0 ? s[0] / s[1] : 0;
    traceScdiaEvent(
        trace,
        obs,
        pass ? "fixed" : "float",
        frq1,
        frq2,
        frq3,
        nf,
        filterMode,
        s[0],
        s[1],
        ratio,
        F.data()[0],
        nf > 1 ? F.data()[1] : 0,
        nf > 2 ? F.data()[2] : 0,
        pass ? E_SlipDiagReason::LAMBDA_FIXED : E_SlipDiagReason::LAMBDA_FLOAT
    );

    if (filterMode == E_FilterMode::LSQ)
    {
        /* least-squares */
        satStat.amb[0] = 0;
        satStat.amb[1] = 0;
        satStat.amb[2] = 0;
        tracepdeex(2, trace, "(freq=%d) ", nf);
        if (pass)
        {
            tracepdeex(2, trace, "fixed ");
            for (int i = 0; i < nf; i++)
                satStat.amb[i] = ROUND(F.data()[i]);
        }
    }
    else
    {
        /* kalman filter */
        satStat.flt.amb[0] = 0;
        satStat.flt.amb[1] = 0;
        satStat.flt.amb[2] = 0;
        if (pass)
        {
            memset(satStat.flt.a, 0, 3 * sizeof(double));
            memset(satStat.flt.Qa, 0, 9);  // todo? looks sketchy
            satStat.flt.slip |= 2;
            tracepdeex(1, trace, "     ACC fixed ");
            for (int i = 0; i < nf; i++)
            {
                satStat.flt.amb[i] = ROUND(F.data()[i]);
            }
        }
        tracepdeex(1, trace, "ACC epoch used=%2d\n", satStat.flt.ne);
        if (pass)
            satStat.flt.ne = 0;
    }
}

/** Cycle slip detection for dual-frequency
 */
void cycleslip2(
    Trace&   trace,    ///< Trace to output to
    SatStat& satStat,  ///< Persistant satellite status parameters
    lc_t&    lcBase,   ///< Linear combinations
    GObs&    obs       ///< Navigation object for this satellite
)
{
    string timeStr = lcBase.time.to_string(2);

    auto& recOpts = acsConfig.getRecOpts(obs.mount);

    double dt = (lcBase.time - satStat.lc_pre.time).to_double();

    if (dt < 20 || dt > PDEGAP)
    {
        // small interval or reset

        satStat.dIono = 0;
        // approximation of ionosphere residual

        satStat.sigmaIono = PROC_NOISE_IONO * SQRT(dt);
    }
    else
    {
        // medium interval ~30s

        if (satStat.dIono == 0)
        {
            satStat.sigmaIono = PROC_NOISE_IONO * SQRT(dt);
        }
    }

    if (satStat.sigmaIono == 0)
    {
        satStat.sigmaIono = 0.001;
    }

    E_FType frq1;
    E_FType frq2;
    E_FType frq3;
    int     nf = obsFreqs(obs, frq1, frq2, frq3);
    if (nf < 2)
    {
        return;
    }

    auto& lam = obs.satNav_ptr->lamMap;

    double lam1 = lam[frq1];
    double lam2 = lam[frq2];

    double lamw = lam1 * lam2 / (lam2 - lam1);  // todo? rename

    /* ionosphere coefficient */
    double coef = SQR(lam2) / SQR(lam1) - 1;

    E_FType   freqs[] = {frq1, frq2};
    SlipNoise noise;
    if (!slipNoise(obs, freqs, 2, noise))
    {
        traceSlipEvent(trace, "PDE", obs, "skipped", frq1, frq2, NONE, 0, 0, 0, noise.reason);
        return;
    }

    double sigmaCode  = noise.sigmaCode;
    double sigmaPhase = noise.sigmaPhase;

    double sigmaGF = 2 * sigmaPhase;

    S_LC lcNew = getLC(lcBase, frq1, frq2);
    S_LC lcPre = getLC(satStat.lc_pre, frq1, frq2);

    double mwNoise = mwnoise(sigmaCode, sigmaPhase, lam1, lam2);

    /* averaged MW measurement and noise */
    double fNw;
    if (acsConfig.preprocOpts.mw_proc_noise)
    {
        fNw = lcNew.MW_c - satStat.mwSlip.mean;
    }
    else
    {
        fNw = lcNew.MW_c - lcPre.MW_c;
    } /* Eq (6) in TN */

    /* clock jump */
    if (fabs(fNw * lamw) > 10e-3 * CLIGHT)
    {
        tracepdeex(1, trace, "Potential clock jump rather than cycle slip -cs2\n");
    }

    double deltaGF = lcNew.GF_Phas_m - lcPre.GF_Phas_m; /* Eq (9) in TN */

    tracepdeex(
        2,
        trace,
        "\nPDE-CS GPST DUAL  %s %4s %5.2f %5.3f %8.4f %7.4f %8.4f                           "
        "     ",
        timeStr.c_str(),
        lcBase.Sat.id().c_str(),
        satStat.el * R2D,
        lamw,
        deltaGF,
        fNw,
        sigmaGF
    );

    /* cycle slip detection */
    if (satStat.el >= recOpts.elevation_mask_deg * D2R)
    {
        scdia(
            trace,
            satStat,
            lcBase,
            obs,
            lam,
            sigmaPhase,
            sigmaCode,
            2,
            E_FilterMode::LSQ,
            frq1,
            frq2,
            frq3
        );
    }

    /* update TD ionosphere residual */
    if (satStat.sigStatMap[ft2string(frq1)].slip.any == 0 &&
        satStat.sigStatMap[ft2string(frq2)].slip.any == 0)
    {
        satStat.dIono     = deltaGF / coef;
        satStat.sigmaIono = sigmaGF / coef;
    }
}

/** Cycle slip detection and repair for triple-frequency
 */
void cycleslip3(
    Trace&   trace,    ///< Trace to output to
    SatStat& satStat,  ///< Persistant satellite status parameters
    lc_t&    lc,       ///< Linear combinations
    GObs&    obs       ///< Navigation object for this satellite
)
{
    string timeStr = lc.time.to_string(2);

    auto& recOpts = acsConfig.getRecOpts(obs.mount);

    double dt = (lc.time - satStat.lc_pre.time).to_double();

    /* small interval */
    if (dt < 20)
    {
        satStat.dIono = 0;

        /* approximation of ionosphere residual */
        satStat.sigmaIono = PROC_NOISE_IONO * SQRT(dt);
    }
    else
    {
        /* large interval */
        if (satStat.sigmaIono == 0)
        {
            satStat.sigmaIono = PROC_NOISE_IONO * SQRT(dt);
        }
    }

    if (satStat.sigmaIono == 0)
    {
        satStat.sigmaIono = 0.001;
    }

    E_FType frq1;
    E_FType frq2;
    E_FType frq3;
    int     nf = obsFreqs(obs, frq1, frq2, frq3);
    if (nf < 3)
        return;

    auto&  lam  = obs.satNav_ptr->lamMap;
    double lam1 = lam[frq1];
    double lam2 = lam[frq2];
    double lam5 = lam[frq3];

    /* TD MW noise (m) */
    double lamew = lam2 * lam5 / (lam5 - lam2);
    if (lamew < 0)
        lamew *= -1;

    E_FType   freqs[] = {frq1, frq2, frq3};
    SlipNoise noise;
    if (!slipNoise(obs, freqs, 3, noise))
    {
        traceSlipEvent(trace, "PDE", obs, "skipped", frq1, frq2, frq3, 0, 0, 0, noise.reason);
        return;
    }

    double sigmaCode  = noise.sigmaCode;
    double sigmaPhase = noise.sigmaPhase;

    double mwNoise12 = mwnoise(sigmaCode, sigmaPhase, lam1, lam2);
    double mwNoise15 = mwnoise(sigmaCode, sigmaPhase, lam1, lam5);
    double mwNoise25 = mwnoise(sigmaCode, sigmaPhase, lam2, lam5);

    double sigmaGF = 2 * sigmaPhase; /* TD GF noise */

    S_LC lc25new = getLC(lc, frq2, frq3);
    S_LC lc25pre = getLC(satStat.lc_pre, frq2, frq3);

    /* averaged EMW measurement and noise */
    double fNew;
    // 	double sigmaEMW;
    if (acsConfig.preprocOpts.mw_proc_noise)
    {
        fNew = lc25new.MW_c - satStat.emwSlip.mean;
    }
    else
    {
        fNew = lc25new.MW_c - lc25pre.MW_c;
    } /* Eq (13) in TN */

    double deltaGF25 = lc25new.GF_Phas_m - lc25pre.GF_Phas_m;

    /* clock jump */
    if (fabs(fNew * lamew) > 10e-3 * CLIGHT)
    {
        fprintf(stdout, "Potential clock jump rather than cycle slip -cs3\n");
        return;
    }

    /* ionosphere coefficient for L2 & L5 */
    double coef1 =
        SQR(CLIGHT / lam1) / SQR(CLIGHT / lam5) - SQR(CLIGHT / lam1) / SQR(CLIGHT / lam2);
    if (coef1 < 0)
        coef1 = -coef1;

    S_LC lcNew = getLC(lc, frq1, frq2);
    S_LC lcPre = getLC(satStat.lc_pre, frq1, frq2);

    double lamw = lam1 * lam2 / (lam2 - lam1);

    double coef = SQR(lam2) / SQR(lam1) - 1;  // ionosphere coefficient for L1 & LX

    /* averaged MW measurement and noise */
    double fNw;
    if (acsConfig.preprocOpts.mw_proc_noise)
    {
        fNw = lcNew.MW_c - satStat.mwSlip.mean;
    }
    else
    {
        fNw = lcNew.MW_c - lcPre.MW_c;
    } /* Eq (6) in TN */

    double deltaGF = lcNew.GF_Phas_m - lcPre.GF_Phas_m;

    tracepdeex(
        2,
        trace,
        "\nPDE-CS GPST TRIP  %s %4s %5.2f %5.3f %8.4f %7.4f %8.4f        %6.2f %8.4f %7.4f ",
        timeStr.c_str(),
        lc.Sat.id().c_str(),
        satStat.el * R2D,
        lamw,
        deltaGF,
        fNw,
        sigmaGF,
        lamew,
        deltaGF25,
        fNew
    );

    if (satStat.el >= recOpts.elevation_mask_deg * D2R)
    {
        scdia(
            trace,
            satStat,
            lc,
            obs,
            lam,
            sigmaPhase,
            sigmaCode,
            3,
            E_FilterMode::LSQ,
            frq1,
            frq2,
            frq3
        );
    }

    /* update TD ionosphere residual */
    if (satStat.sigStatMap[ft2string(frq1)].slip.any == 0 &&
        satStat.sigStatMap[ft2string(frq2)].slip.any == 0 &&
        satStat.sigStatMap[ft2string(frq3)].slip.any == 0)
    {
        satStat.dIono     = deltaGF / coef;
        satStat.sigmaIono = sigmaGF / coef;
    }
}

/** Cycle slip detection and repair
 */
void detectslip(
    Trace&   trace,    ///< Trace to output to
    SatStat& satStat,  ///< Persistant satellite status parameters
    lc_t&    lc_new,   ///< Linear combination for this epoch
    lc_t&    lc_old,   ///< Linear combination from previous epoch
    GObs&    obs       ///< Navigation object for this satellite
)
{
    bool  dualFreq = false;
    E_Sys sys      = lc_new.Sat.sys;

    char id[32];
    lc_new.Sat.getId(id);

    string timeStr = lc_new.time.to_string(2);

    auto& recOpts = acsConfig.getRecOpts(obs.mount);

    if (acsConfig.process_sys[sys] == false)
        return;

    E_FType frq1;
    E_FType frq2;
    E_FType frq3;
    int     nf = obsFreqs(obs, frq1, frq2, frq3);
    if (nf < 2)
    {
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "skipped",
            NONE,
            NONE,
            NONE,
            0,
            0,
            0,
            E_SlipDiagReason::NO_FREQUENCIES
        );
        return;
    }

    /* first epoch or large gap or low elevation */  // todo? initialisation stuff, remove
    if (satStat.lc_pre.time.bigTime == 0 || satStat.el < recOpts.elevation_mask_deg * D2R ||
        lc_new.time > lc_old.time + PDEGAP)
    {
        satStat.mwSlip  = {};
        satStat.emwSlip = {};

        if (lc_new.time > lc_old.time + PDEGAP)
            tracepdeex(
                1,
                trace,
                "\nPDE-CS GPST       %s %4s %5.2f --time gap --",
                timeStr.c_str(),
                id,
                satStat.el * R2D
            );
        if (satStat.el < recOpts.elevation_mask_deg * D2R)
            tracepdeex(
                1,
                trace,
                "\nPDE-CS GPST       %s %4s %5.2f --low_elevation --",
                timeStr.c_str(),
                id,
                satStat.el * R2D
            );
        else
            tracepdeex(
                1,
                trace,
                "\nPDE-CS GPST       %s %4s %5.2f --satStat.lc_pre.time.time --",
                timeStr.c_str(),
                id,
                satStat.el * R2D
            );

        E_SlipDiagReason reason = E_SlipDiagReason::FIRST_EPOCH;
        if (lc_new.time > lc_old.time + PDEGAP)
        {
            reason = E_SlipDiagReason::TIME_GAP;
        }
        else if (satStat.el < recOpts.elevation_mask_deg * D2R)
        {
            reason = E_SlipDiagReason::LOW_ELEVATION;
        }
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "initialised",
            frq1,
            frq2,
            frq3,
            satStat.el * R2D,
            recOpts.elevation_mask_deg,
            PDEGAP,
            reason
        );

        return;
    }

    if (nf == 2 && lc_new.L_m[frq1] != 0 && lc_new.L_m[frq2] != 0)
    {
        dualFreq = true;
    }

    if (dualFreq && lc_old.L_m[frq1] != 0 && lc_old.L_m[frq2] != 0)
    {
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "evaluating",
            frq1,
            frq2,
            NONE,
            2,
            0,
            0,
            E_SlipDiagReason::DUAL_FREQUENCY
        );
        cycleslip2(trace, satStat, lc_new, obs);

        /* update averaged MW noise when no cycle slip */
        if (satStat.sigStatMap[ft2string(frq1)].slip.any == 0 &&
            satStat.sigStatMap[ft2string(frq2)].slip.any == 0)
        {
            S_LC& lc12 = getLC(lc_new, frq1, frq2);
            lowPassFilter(satStat.mwSlip, lc12.MW_c, acsConfig.preprocOpts.mw_proc_noise);
        }
        else
        {
            satStat.mwSlip = {};
        }
    }
    /* track L5 again */
    else if (
        nf >= 3 && lc_new.L_m[frq1] != 0 && lc_new.L_m[frq2] != 0 && lc_new.L_m[frq3] != 0 &&
        lc_old.L_m[frq1] != 0 && lc_old.L_m[frq2] != 0 && lc_old.L_m[frq3] == 0
    )  // was zero, now not.
    {
        /* set slip flag for L5 (introduce new ambiguity for L5) */
        if (!acsConfig.zhangFullRank.enable)
        {
            satStat.sigStatMap[ft2string(frq3)].slip.retrack = true;
            satStat.sigStatMap[ft2string(frq3)].savedSlip.retrack = true;
        }
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "retracking",
            frq1,
            frq2,
            frq3,
            3,
            2,
            0,
            E_SlipDiagReason::THIRD_FREQUENCY_REACQUIRED
        );
        cycleslip2(trace, satStat, lc_new, obs);

        /* update averaged MW noise when no cycle slip */
        if (satStat.sigStatMap[ft2string(frq1)].slip.any == 0 &&
            satStat.sigStatMap[ft2string(frq2)].slip.any == 0)
        {
            S_LC& lc12 = getLC(lc_new, frq1, frq2);
            lowPassFilter(satStat.mwSlip, lc12.MW_c, acsConfig.preprocOpts.mw_proc_noise);
        }
        else
        {
            satStat.mwSlip = {};
        }
    }
    /* Triple-frequency */
    else if (
        nf >= 3 && lc_new.L_m[frq1] != 0 && lc_new.L_m[frq2] != 0 && lc_new.L_m[frq3] != 0 &&
        lc_old.L_m[frq1] != 0 && lc_old.L_m[frq2] != 0 && lc_old.L_m[frq3] != 0
    )
    {
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "evaluating",
            frq1,
            frq2,
            frq3,
            3,
            0,
            0,
            E_SlipDiagReason::TRIPLE_FREQUENCY
        );
        cycleslip3(trace, satStat, lc_new, obs);

        if (satStat.el * R2D > 30)
        {
            if (satStat.sigStatMap[ft2string(frq1)].slip.any == 2  // todo? check the 2
                && satStat.amb[0] == 0 && satStat.amb[1] == 0 && satStat.amb[2] == 0)
            {
                satStat.sigStatMap[ft2string(frq1)].slip.any = 0;
                satStat.sigStatMap[ft2string(frq2)].slip.any = 0;
                satStat.sigStatMap[ft2string(frq3)].slip.any = 0;
            }
        }

        /*update averaged MW25 noise when no cycle slip */
        if (satStat.sigStatMap[ft2string(frq1)].slip.any == 0 &&
            satStat.sigStatMap[ft2string(frq2)].slip.any == 0 &&
            satStat.sigStatMap[ft2string(frq3)].slip.any == 0)
        {
            S_LC& lc25 = getLC(lc_new, frq2, frq3);
            lowPassFilter(satStat.emwSlip, lc25.MW_c, acsConfig.preprocOpts.mw_proc_noise);
        }
        else
        {
            satStat.emwSlip = {};
        }
    }
    /* track L1 or L2 again, new rising satellite */
    else if (dualFreq && (lc_old.L_m[frq1] == 0 || lc_old.L_m[frq2] == 0))
    {
        satStat.flt.slip = 0;
        satStat.flt.ne   = 0;
        for (auto& [key, sigStat] : satStat.sigStatMap)
        {
            if (!acsConfig.zhangFullRank.enable)
            { sigStat.slip.retrack=true; sigStat.savedSlip.retrack=true; }
        }

        tracepdeex(
            1,
            trace,
            "\nPDE-CS GPST       %s %4s %5.2f --  re-tracking   --\n",
            timeStr.c_str(),
            id,
            satStat.el * R2D
        );
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "retracking",
            frq1,
            frq2,
            NONE,
            2,
            0,
            0,
            E_SlipDiagReason::FREQUENCY_REACQUIRED
        );
    }
    else
    {
        satStat.flt.slip = 0;
        satStat.flt.ne   = 0;
        for (auto& [key, sigStat] : satStat.sigStatMap)
        {
            if (!acsConfig.zhangFullRank.enable)
            { sigStat.slip.singleFreq=true; sigStat.savedSlip.singleFreq=true; }
        }

        tracepdeex(
            1,
            trace,
            "\nPDE-CS GPST       %s %4s %5.2f --single frequency--\n",
            timeStr.c_str(),
            id,
            satStat.el * R2D
        );
        traceSlipEvent(
            trace,
            "PDE",
            obs,
            "flagged",
            frq1,
            frq2,
            frq3,
            1,
            0,
            0,
            E_SlipDiagReason::SINGLE_FREQUENCY
        );
    }
}

void clearSlips(ObsList& obsList)
{
    // clear non-persistent status values.
    for (auto& obs : only<GObs>(obsList))
    {
        if (acsConfig.process_sys[obs.Sat.sys] == false)
        {
            continue;
        }

        auto& satOpts = acsConfig.getSatOpts(obs.Sat);

        if (satOpts.exclude)
        {
            continue;
        }

        for (auto& [sigKey, sigStat] : obs.satStat_ptr->sigStatMap)
        {
            SatStat& satStat = *(obs.satStat_ptr);

            satStat.slip     = false;  // todo? is this used?
            sigStat.slip.any = 0;
            sigStat.previousPhaseQuarantined=sigStat.phaseQuarantined;
            if(!sigStat.tracking.pendingBreak)sigStat.phaseContinuityUnproven=false;
            sigStat.phaseQuarantined=false;
            sigStat.codeQuarantined=false;
        }
    }
}

void detectZhangEvidenceSlips(Trace& trace, ObsList& obsList);

// Optional OFFLINE migration of old schema-2 checkpoints. Only detector history
// is copied; the authoritative receiver tracking, graph and KF are untouched.
// The E29 input manifest has already authenticated this same local source file.
// A separate empty Navigation is used by the RINEX header reader to avoid writes
// to shared navigation during parallel receiver preprocessing.
void warmZhangSlipEvidenceFromRaw(Trace& trace, ObsList& current)
{
    static thread_local bool warming=false;
    const char* setting=std::getenv("ZHANG_R51_SLIP_RAW_WARMUP");
    if(warming || !setting || std::string(setting)!="1")return;
    Receiver* receiver=nullptr;double now=0;bool hasPrevious=false;
    for(auto& obs:only<GObs>(current))
    {
        if(!obs.rec_ptr || !obs.satStat_ptr)continue;
        if(!obs.satStat_ptr->slipEvidence.accepted.empty())return;
        if(receiver && receiver!=obs.rec_ptr)return;
        receiver=obs.rec_ptr;now=obs.time.bigTime;
        hasPrevious |= obs.satStat_ptr->lc_pre.time.bigTime>0;
    }
    if(!receiver || !hasPrevious)return;
    std::string source=receiver->source;
    if(source.rfind("file://",0)==0)source.erase(0,7);
    if(source.empty() || source.find("://")!=std::string::npos)return;
    std::ifstream input(source);
    if(!input)
    {
        for(auto& obs:only<GObs>(current))
            if(obs.satStat_ptr)
                for(const auto& [frequency,signal]:obs.sigs)
                    obs.satStat_ptr->sigStatMap[ft2string(frequency)].phaseQuarantined=true;
        trace<<"\nZHANG_SLIP_RAW_WARMUP rec="<<receiver->id<<" copied=0 reason=SOURCE_UNREADABLE";
        return;
    }
    struct Guard {bool& active;Guard(bool& a):active(a){active=true;}~Guard(){active=false;}} guard(warming);
    const auto started=std::chrono::steady_clock::now();
    trace<<"\nZHANG_SLIP_RAW_WARMUP_BEGIN rec="<<receiver->id<<" resume_time="<<std::to_string(now)
        <<" historical_records_follow=1 scope=DETECTOR_ONLY";
    map<SatSys,SatStat> scratch;
    Navigation headerOnly;RinexStation station;
    char type=0;double version=0;E_Sys system=E_Sys::NONE;E_TimeSys timeSystem=E_TimeSys::GPST;
    map<E_Sys,map<int,CodeType>> codes;
    int epochs=0,restored=0;
    while(input.good())
    {
        ObsList parsed;readRnx(input,type,parsed,headerOnly,station,version,system,timeSystem,codes);
        if(parsed.empty())continue;
        const double time=parsed.front()->time.bigTime;
        if(time>=now)break;
        if(time<now-1920)continue; // fill the bounded 64-sample retrospective window
        ObsList prepared;
        for(auto& obs:only<GObs>(parsed))
        {
            if(obs.Sat.sys!=E_Sys::GPS || !acsConfig.process_sys[obs.Sat.sys])continue;
            const auto navigation=nav.satNavMap.find(obs.Sat);
            if(navigation==nav.satNavMap.end())continue; // never insert a global map in this parallel loop
            if(!navigation->second.lamMap.count(F1) || !navigation->second.lamMap.count(F2))continue;
            obs.mount=receiver->id;obs.rec_ptr=receiver;
            obs.satStat_ptr=&scratch[obs.Sat];obs.satNav_ptr=&navigation->second;
            for(const auto code:acsConfig.code_priorities[obs.Sat.sys])
            {
                const auto frequency=code2Freq[obs.Sat.sys][code];
                if(obs.sigs.count(frequency))continue;
                for(const auto& signal:obs.sigsLists[frequency])
                    if(signal.code==code && signal.L!=0){obs.sigs[frequency]=signal;break;}
            }
            if(!satPosClk(trace,obs.time,obs,nav,{E_Source::PRECISE},{E_Source::PRECISE},
                nullptr,nullptr,E_OffsetType::COM))continue;
            VectorEcef los;geodist(obs.rSatCom,receiver->aprioriPos,los);
            satazel(ecef2pos(receiver->aprioriPos),los,*obs.satStat_ptr);
            if(obs.satStat_ptr->el<acsConfig.getRecOpts(receiver->id).elevation_mask_deg*D2R)continue;
            obsVariance(obs);
            prepared.push_back(static_cast<std::shared_ptr<GObs>>(obs));
        }
        clearSlips(prepared);
        for(auto& obs:only<GObs>(prepared))
        {obs.satStat_ptr->lc_pre=obs.satStat_ptr->lc_new;obs.satStat_ptr->lc_new={};}
        obs2lcs(trace,prepared);detslp_ll(trace,prepared);detectZhangEvidenceSlips(trace,prepared);
        epochs++;
    }
    for(auto& obs:only<GObs>(current))
    {
        if(!obs.satStat_ptr)continue;
        const auto found=scratch.find(obs.Sat);
        if(found==scratch.end())
        {
            if(obs.satStat_ptr->lc_pre.time.bigTime>0)
                for(const auto& [frequency,signal]:obs.sigs)
                    obs.satStat_ptr->sigStatMap[ft2string(frequency)].phaseQuarantined=true;
            continue;
        }
        const auto& warm=found->second;
        const auto& previous=obs.satStat_ptr->lc_pre;
        bool same=warm.lc_new.time==previous.time && !warm.slipEvidence.accepted.empty();
        for(const auto& [frequency,value]:previous.L_m)
            same &= warm.lc_new.L_m.count(frequency) && std::abs(warm.lc_new.L_m.at(frequency)-value)<1e-7;
        for(const auto& [frequency,value]:previous.P)
            same &= warm.lc_new.P.count(frequency) && std::abs(warm.lc_new.P.at(frequency)-value)<1e-7;
        if(!same)
        {
            for(const auto& [frequency,signal]:obs.sigs)
                obs.satStat_ptr->sigStatMap[ft2string(frequency)].phaseQuarantined=true;
            trace<<"\nZHANG_SLIP_RAW_WARMUP rec="<<receiver->id<<" sat="<<obs.Sat.id()
                <<" copied=0 reason=PREVIOUS_RAW_LC_MISMATCH";
            continue;
        }
        obs.satStat_ptr->slipEvidence=warm.slipEvidence;restored++;
    }
    trace<<"\nZHANG_SLIP_RAW_WARMUP rec="<<receiver->id<<" epochs="<<epochs<<" copied="<<restored
        <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()
        <<" kf_changed=0 graph_changed=0 authoritative_tracking_changed=0";
    trace<<"\nZHANG_SLIP_RAW_WARMUP_END rec="<<receiver->id<<" historical_records_follow=0";
}

// Capture all receiver-link inputs before committing any detector state. L5 is
// an auxiliary QC observation only: it does not enter the configured AR model.
void detectZhangEvidenceSlips(Trace& trace, ObsList& obsList)
{
    warmZhangSlipEvidenceFromRaw(trace,obsList);
    struct Input { GObs* obs; E_FType first, second, third; ZhangSlipSample sample; bool hard=false; };
    std::vector<Input> inputs;
    for (auto& obs:only<GObs>(obsList))
    {
        if(obs.exclude || !obs.satStat_ptr || !obs.satNav_ptr || !obs.rec_ptr) continue;
        E_FType first,second,third;
        if(obsFreqs(obs,first,second,third)<2)
        {
            for(const auto& [frequency,signal]:obs.sigs)
                if(signal.L!=0)obs.satStat_ptr->sigStatMap[ft2string(frequency)].phaseQuarantined=true;
            continue;
        }
        const auto& one=obs.sigs.at(first); const auto& two=obs.sigs.at(second);
        if(one.LLI&2 || two.LLI&2) continue;
        if(obs.satStat_ptr->sigStatMap[ft2string(first)].phaseQuarantined ||
            obs.satStat_ptr->sigStatMap[ft2string(second)].phaseQuarantined)continue;
        const Sig* auxiliary=nullptr;
        E_FType auxiliaryFrequency=NONE;
        for(const auto& [frequency,signals]:obs.sigsLists)
        {
            if(frequency==first || frequency==second || !obs.satNav_ptr->lamMap.count(frequency) ||
                !std::isfinite(obs.satNav_ptr->lamMap.at(frequency)) || obs.satNav_ptr->lamMap.at(frequency)<=0)continue;
            for(const auto& signal:signals)
                if(signal.L!=0 && !signal.invalid && !(signal.LLI&3))
                { auxiliary=&signal;auxiliaryFrequency=frequency;break; }
            if(auxiliary) break;
        }
        const double l1=obs.satNav_ptr->lamMap[first], l2=obs.satNav_ptr->lamMap[second];
        if(l1<=0 || l2<=0 || l1==l2 || one.phasVar<=0 || two.phasVar<=0)
        {
            obs.satStat_ptr->sigStatMap[ft2string(first)].phaseQuarantined=true;
            obs.satStat_ptr->sigStatMap[ft2string(second)].phaseQuarantined=true;
            continue;
        }
        const double l3=auxiliary ? obs.satNav_ptr->lamMap[auxiliaryFrequency] : 0;
        const bool code=one.P!=0 && two.P!=0 && one.codeVar>0 && two.codeVar>0;
        // The network configuration prefers KALMAN satellite clocks. Those
        // already contain this link's estimation history and cannot provide an
        // independent slip check. Use a COPY with external precise clocks and
        // positions only, leaving the estimation observation itself unchanged.
        GObs independent=obs;
        if(independent.clkSource!=E_Source::PRECISE || independent.posSource!=E_Source::PRECISE)
            satPosClk(trace,independent.time,independent,nav,{E_Source::PRECISE},{E_Source::PRECISE},
                nullptr,nullptr,E_OffsetType::COM);
        VectorEcef los;
        const double range=geodist(independent.rSatCom,obs.rec_ptr->aprioriPos,los);
        const bool geometry=independent.ephPosValid && independent.ephClkValid &&
            independent.clkSource==E_Source::PRECISE && independent.posSource==E_Source::PRECISE && range>0;
        const int signature=static_cast<int>(one.code)*1000+static_cast<int>(two.code);
        auto sample=zhangMakeSlipSample(obs.time.bigTime,signature,{l1,l2,l3},
            {one.L*l1,two.L*l2,auxiliary ? auxiliary->L*l3 : 0},
            {one.phasVar,two.phasVar,auxiliary ? std::max(auxiliary->phasVar,one.phasVar) : 0},
            {one.P,two.P},{one.codeVar,two.codeVar},auxiliary,code,
            range-CLIGHT*independent.satClk,geometry);
        sample.thirdSignature=auxiliary ? static_cast<int>(auxiliary->code) : 0;
        sample.gfRateValid=one.D!=0 && two.D!=0 && std::isfinite(one.D) && std::isfinite(two.D);
        if(sample.gfRateValid)sample.gfRate=-l1*one.D+l2*two.D;
        // Schema-2 FLOAT checkpoints have the previous raw LC but no new QC
        // history. Rebuild its independent geometry once; do not blindly accept
        // the first new sample against an already-existing physical ambiguity.
        auto& evidence=obs.satStat_ptr->slipEvidence;
        const auto& previous=obs.satStat_ptr->lc_pre;
        if(evidence.accepted.empty() && geometry && previous.time.bigTime>0 &&
            obs.time>previous.time && (obs.time-previous.time).to_double()<=90 &&
            previous.L_m.count(first) && previous.L_m.count(second) &&
            previous.P.count(first) && previous.P.count(second) &&
            previous.L_m.at(first)!=0 && previous.L_m.at(second)!=0 &&
            previous.P.at(first)!=0 && previous.P.at(second)!=0)
        {
            GObs old=obs;old.time=previous.time;
            for(auto& [frequency,signal]:old.sigs)
            {
                signal.P=previous.P.count(frequency) ? previous.P.at(frequency) : 0;
                signal.L=previous.L_m.count(frequency) ? previous.L_m.at(frequency)/obs.satNav_ptr->lamMap[frequency] : 0;
            }
            SatStat oldStatus;old.satStat_ptr=&oldStatus;
            if(satPosClk(trace,old.time,old,nav,{E_Source::PRECISE},{E_Source::PRECISE},
                nullptr,nullptr,E_OffsetType::COM))
            {
                VectorEcef oldLos;
                const double oldRange=geodist(old.rSatCom,obs.rec_ptr->aprioriPos,oldLos);
                satazel(ecef2pos(obs.rec_ptr->aprioriPos),oldLos,oldStatus);
                const double scale=std::pow(std::sin(obs.satStat_ptr->el)/
                    std::max(.05,std::sin(oldStatus.el)),2);
                auto seed=zhangMakeSlipSample(previous.time.bigTime,signature,{l1,l2,0},
                    {previous.L_m.at(first),previous.L_m.at(second),0},
                    {one.phasVar*scale,two.phasVar*scale,0},
                    {previous.P.at(first),previous.P.at(second)},
                    {one.codeVar*scale,two.codeVar*scale},false,true,oldRange-CLIGHT*old.satClk,true);
                evidence.accepted.push_back(seed);evidence.screening.push_back(seed);
                tracepdeex(2,trace,"\nZHANG_SLIP_CHECKPOINT_SEED rec=%s sat=%s time=%s auxiliary_reconstructed=0",
                    obs.mount.c_str(),obs.Sat.id().c_str(),old.time.to_string().c_str());
            }
        }
        const bool hard=obs.satStat_ptr->sigStatMap[ft2string(first)].slip.any ||
            obs.satStat_ptr->sigStatMap[ft2string(second)].slip.any;
        inputs.push_back({&obs,first,second,auxiliaryFrequency,sample,hard});
    }
    std::vector<std::vector<ZhangSlipClockPeer>> peers(inputs.size());
    for(size_t i=0;i<inputs.size();i++)
        for(size_t j=0;j<inputs.size();j++)
        {
            if(i==j || inputs[i].obs->rec_ptr!=inputs[j].obs->rec_ptr ||
                inputs[j].hard || !inputs[j].sample.geometry ||
                inputs[i].obs->Sat.sys!=inputs[j].obs->Sat.sys) continue;
            const auto* reference=zhangSlipClockReference(inputs[i].obs->satStat_ptr->slipEvidence);
            const auto& history=inputs[j].obs->satStat_ptr->slipEvidence.screening;
            if(!reference)continue;
            const auto previous=std::find_if(history.begin(),history.end(),[&](const auto& old)
                {return old.time==reference->time && old.geometry && old.signature==inputs[j].sample.signature;});
            if(previous==history.end())continue;
            peers[i].push_back({previous->time,inputs[j].sample.rawIf-previous->rawIf,
                inputs[j].sample.covariance[0]+previous->covariance[0],history,inputs[j].sample,true});
        }
    for(size_t i=0;i<inputs.size();i++)
    {
        auto& input=inputs[i]; auto& obs=*input.obs;
        ZhangSlipTrendDiagnostic trend;
        auto decision=zhangInspectSlipEvidence(obs.satStat_ptr->slipEvidence,input.sample,peers[i],input.hard,{},&trend);
        const std::array<E_FType,3> frequencies{input.first,input.second,input.third};
        for(int f=0;f<2;f++)
        {
            auto& status=obs.satStat_ptr->sigStatMap[ft2string(frequencies[f])];
            status.phaseQuarantined|=!decision.usePhase;
            status.codeQuarantined|=!decision.useCode;
            if(decision.breakArc && decision.retire[f])
            {
                status.slip.SCDIA=true;status.savedSlip.SCDIA=true;
                status.tracking.pendingBreak=true;
                status.phaseContinuityUnproven=decision.action==ZhangSlipAction::UNRESOLVED_STEP;
            }
        }
        tracepdeex(2,trace,
            "\nZHANG_SLIP_EVIDENCE epoch=%s rec=%s sat=%s action=%s phase=%d code=%d break=%d "
            "clock_model=MATCHED_PEER_WINDOW clock_partners=%d if=%.9f gf=%.9f mw=%.9f triple=%.9f "
            "sigma_if=%.9f sigma_gf=%.9f sigma_mw=%.9f sigma_triple=%.9f "
            "phase_score=%.9f change_score=%.9f change_time=%.0f best=%.9f ratio_lower_bound=%.9f n1=%d n2=%d n3=%d "
            "doppler=%d doppler_assisted=%d doppler_error=%.9f doppler_sigma=%.9f doppler_score=%.9f doppler_integer_fit=%.9f "
            "early_stage=%d trend_qualified=%d trend_assisted=%d trend_window=%d trend_consistency=%.9f trend_score=%.9f trend_phase_support=%.9f trend_fit=%.9f",
            obs.time.to_string().c_str(),obs.mount.c_str(),obs.Sat.id().c_str(),
            zhangSlipActionName(decision.action),decision.usePhase,decision.useCode,decision.breakArc,
            decision.clockPartners,decision.innovation[0],decision.innovation[1],decision.innovation[2],
            decision.innovation[3],decision.sigma[0],decision.sigma[1],decision.sigma[2],decision.sigma[3],
            decision.phaseScore,decision.changePointScore,decision.changePointTime,decision.bestScore,decision.ratio,
            decision.integer[0],decision.integer[1],decision.integer[2],decision.dopplerAvailable,
            decision.dopplerAssisted,decision.dopplerInnovation,decision.dopplerSigma,
            decision.dopplerScore,decision.dopplerIntegerFit,
            trend.stage,trend.qualified,trend.assisted,trend.window,trend.consistency,trend.score,trend.support,trend.fit);
    }
}

/** Detect slips for multiple observations
 */
void detectslips(
    Trace&   trace,   ///< Trace to output to
    ObsList& obsList  ///< List of observations to detect slips within
)
{
    tracepdeex(2, trace, "\n   *-------- PDE cycle slip detection & repair --------*\n");

    detslp_ll(trace, obsList);
    if (acsConfig.zhangFullRank.enable && !acsConfig.zhangPppAr.user_adapter && zhangEvidenceSlipQcEnabled())
    {
        detectZhangEvidenceSlips(trace,obsList);
        for(auto& obs:only<GObs>(obsList))
            if(obs.satStat_ptr)
                for(const auto& [ft,sig]:obs.sigs)
                    obs.satStat_ptr->slip|=obs.satStat_ptr->sigStatMap[ft2string(ft)].slip.any!=0;
        return; // legacy GF persistence and SCDIA do not silently overrule new QC
    }
    detslp_gf(trace, obsList);
    detslp_mw(trace, obsList);

    tracepdeex(
        2,
        trace,
        "\nPDE-CS GPST       epoch                   prn  el   lamw    gf12    mw12     siggf  "
        "sigmw  "
        "lamew     gf25    mw25   "
        "            LC                   N1   N2   N5\n"
    );

    for (auto& obs : only<GObs>(obsList))
    {
        if (obs.exclude)
        {
            continue;
        }

        SatStat& satStat = *(obs.satStat_ptr);

        bool cleanForLegacyPredictor=true;
        if (acsConfig.zhangFullRank.enable)
            for (const auto& [ft,sig]:obs.sigs)
            {
                const auto& status=satStat.sigStatMap[ft2string(ft)];
                cleanForLegacyPredictor &= !status.phaseQuarantined && !status.previousPhaseQuarantined;
            }
        if (cleanForLegacyPredictor)
            detectslip(trace, satStat, satStat.lc_new, satStat.lc_pre, obs);

        for (auto& [ft, sig] : obs.sigs)
        {
            auto& sigStat = obs.satStat_ptr->sigStatMap[ft2string(ft)];

            if (sigStat.slip.any)
            {
                satStat.slip = true;
            }
        }
    }
}
