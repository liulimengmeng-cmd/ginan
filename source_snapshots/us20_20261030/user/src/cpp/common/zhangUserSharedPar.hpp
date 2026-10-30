#pragma once

#include "ambres/GNSSambres.hpp"
#include "common/zhangRatioOnly.hpp"
#include "common/zhangUserArPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct ZhangUserSharedParAudit
{
    int inputDimension = 0;
    int acceptedRank = 0;
    std::vector<int> testedDimensions;
    std::vector<std::string> trialReasons;
    std::vector<double> trialRatios;
    std::vector<double> trialBootstrapSuccess;
    std::vector<int> trialSupportSatellites;
    bool independentSuccessGate = false;
    int minimumRank = 1;
    int maximumTrials = 0;
    std::string reason = "NOT_EVALUATED";
};

// Ratio-only PAR must be explicit: the legacy alpha=0 path does not shrink a
// rejected ratio candidate. Each trial below re-solves the requested reliable
// LAMBDA suffix from the SAME full float mean/covariance. We neither truncate a
// previously found candidate nor replace a mixed suffix with principal rows.
inline int zhangResolveUserSharedPar(
    Trace& trace, GinAR_mtx& matrix, GinAR_opt options,
    ZhangUserSharedParAudit* auditOutput = nullptr,
    const MatrixXd* physicalSupportTransform = nullptr)
{
    if (!zhangRatioOnly())
        throw std::runtime_error("USER_SHARED_PAR_REQUIRES_RATIO_ONLY_POLICY");
    ZhangUserSharedParAudit audit;
    const auto policy = zhangUserArReliabilityPolicy();
    audit.independentSuccessGate = policy.enforceBootstrapSuccess;
    audit.minimumRank = policy.minimumParRank;
    audit.maximumTrials = policy.maximumParTrials;
    const GinAR_mtx frozen = matrix;
    const int dimension = frozen.aflt.size();
    audit.inputDimension = dimension;
    auto reject = [&](const std::string& reason)
    {
        matrix = frozen;
        matrix.Ztrs.resize(0, dimension);
        matrix.zfix.resize(0);
        audit.reason = reason;
        if (auditOutput) *auditOutput = audit;
        trace << "\nZHANG_USER_SHARED_PAR input_dimension=" << dimension
              << " accepted_rank=0 status=" << reason
              << " same_full_covariance=1 nis_evaluations=0 jitter=0";
        return 0;
    };
    if (dimension == 0) return reject("NO_SEARCHABLE_INTEGER_DIRECTIONS");
    if (frozen.Paflt.rows() != dimension || frozen.Paflt.cols() != dimension ||
        !frozen.aflt.allFinite() || !frozen.Paflt.allFinite() ||
        !std::isfinite(options.ratthr) || options.ratthr <= 1)
        return reject("PAR_INPUT_INVALID");
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(
        (0.5 * (frozen.Paflt + frozen.Paflt.transpose())).eval());
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
        return reject("PAR_INPUT_EIGEN_FAILED");
    const double positiveTolerance = std::max(1e-14,
        std::max(1.0, eigen.eigenvalues().cwiseAbs().maxCoeff()) *
        std::numeric_limits<double>::epsilon() * dimension * 64);
    if (eigen.eigenvalues().minCoeff() <= positiveTolerance)
        return reject("PAR_INPUT_NOT_POSITIVE_DEFINITE");
    options.mode = E_ARmode::LAMBDA_ALT;
    options.min_lambda_fix_count = policy.minimumParRank;
    options.force_bootstrap_success = policy.enforceBootstrapSuccess;
    options.maximum_lambda_conditional_sigma = policy.maximumConditionalSigmaCycles;
    options.lambda_candidate_nis_alpha = 0;
    options.nset = std::max(2, options.nset);
    options.lambda_candidate_row_ablation = "NONE";
    options.lambda_candidate_ablation_target_columns.clear();
    if (policy.safe && (!physicalSupportTransform ||
        physicalSupportTransform->rows() != dimension ||
        physicalSupportTransform->cols() % 2 || !physicalSupportTransform->allFinite()))
        return reject("SAFE_PAR_PHYSICAL_SUPPORT_REQUIRED");
    for (int tested = dimension; tested >= policy.minimumParRank; --tested)
    {
        if (policy.maximumParTrials > 0 && audit.testedDimensions.size() >= policy.maximumParTrials)
            return reject("PAR_TRIAL_BUDGET_EXHAUSTED");
        auto trial = frozen;
        options.max_lambda_fix_count = tested;
        const int accepted = GNSS_AR(trace, trial, options);
        if (trial.aflt.size() != dimension || trial.Paflt.rows() != dimension ||
            trial.Paflt.cols() != dimension ||
            (trial.aflt.array() != frozen.aflt.array()).any() ||
            (trial.Paflt.array() != frozen.Paflt.array()).any())
            throw std::runtime_error("USER_SHARED_PAR_FLOAT_INPUT_MUTATED");
        if (trial.searchDiagnostic.localNisExecuted || trial.lambda_candidate_nis_valid)
            throw std::runtime_error("USER_SHARED_PAR_UNEXPECTED_NIS_EXECUTION");
        audit.testedDimensions.push_back(tested);
        audit.trialReasons.push_back(trial.searchDiagnostic.reason);
        audit.trialRatios.push_back(trial.searchDiagnostic.ratio);
        audit.trialBootstrapSuccess.push_back(trial.lambda_selected_bootstrap_success);
        int supportSatellites = 0;
        if (accepted > 0 && physicalSupportTransform)
        {
            const MatrixXd support = trial.Ztrs * *physicalSupportTransform;
            for (int c = 0; c < support.cols(); c += 2)
                if (!support.middleCols(c, 2).isZero(1e-10)) ++supportSatellites;
        }
        audit.trialSupportSatellites.push_back(supportSatellites);
        const bool supportAccepted = !policy.safe ||
            supportSatellites >= policy.minimumSupportSatellites;
        const std::string trialReason = accepted > 0 && !supportAccepted
            ? "INSUFFICIENT_PHYSICAL_SUPPORT_SATELLITES" : trial.searchDiagnostic.reason;
        audit.trialReasons.back() = trialReason;
        trace << "\nZHANG_USER_SHARED_PAR_TRIAL input_dimension=" << dimension
              << " tested_dimension=" << tested << " accepted_rows=" << accepted
              << " ratio=" << trial.searchDiagnostic.ratio
              << " ratio_threshold=" << options.ratthr
              << " complete=" << trial.searchDiagnostic.ilsComplete
              << " reason=" << trialReason
              << " bootstrap_success=" << trial.lambda_selected_bootstrap_success
              << " success_threshold=" << options.sucthr
              << " independent_success_gate=" << options.force_bootstrap_success
              << " conditional_sigma_max_cycles=" << options.maximum_lambda_conditional_sigma
              << " support_satellites=" << supportSatellites
              << " minimum_rank=" << policy.minimumParRank
              << " max_trials=" << policy.maximumParTrials
              << " risk_claim=UNCALIBRATED_MODEL_CONDITIONAL_NOT_FFRT"
              << " same_full_covariance=1 nis_evaluations=0 jitter=0";
        if (accepted <= 0 || !supportAccepted)
        {
            // If the full request was already reduced by conditional-D gates,
            // do not spend the next trial re-solving the identical suffix.
            if (policy.safe && trial.searchDiagnostic.bootstrapEligibleRank > 0)
                tested = std::min(tested, trial.searchDiagnostic.bootstrapEligibleRank);
            continue;
        }
        if (accepted > tested || accepted < policy.minimumParRank || trial.Ztrs.rows() != accepted ||
            trial.Ztrs.cols() != dimension || trial.zfix.size() != accepted ||
            !trial.searchDiagnostic.ilsComplete || !trial.searchDiagnostic.ratioExecuted)
            throw std::runtime_error("USER_SHARED_PAR_ACCEPTED_DECISION_CONTRACT_INVALID");
        matrix = std::move(trial);
        audit.acceptedRank = accepted;
        audit.reason = accepted == dimension ? "FULL_ILS_ACCEPTED" : "PARTIAL_ILS_ACCEPTED";
        if (auditOutput) *auditOutput = audit;
        trace << "\nZHANG_USER_SHARED_PAR input_dimension=" << dimension
              << " accepted_rank=" << accepted << " attempts=" << audit.testedDimensions.size()
              << " status=" << audit.reason
              << " same_full_covariance=1 nis_evaluations=0 jitter=0";
        return accepted;
    }
    return reject(policy.safe ? "ALL_INDEPENDENT_RELIABILITY_PAR_SUFFIXES_REJECTED"
        : "ALL_RATIO_ONLY_PAR_SUFFIXES_REJECTED");
}
