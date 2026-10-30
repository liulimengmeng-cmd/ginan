#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

// Independent reimplementation of the uploaded GONSIS research strategy,
// not copied RTKLIB/GONSIS code. Provenance in GONSIS-master:
// src/ppp_ar.cpp:433-543 (soft WL and IF constraints), 648-651 (FLOAT
// position gate); src/pppos.cpp:75,1373-1414,1680-1688 (4-sigma fixed
// residual gate and FIX position gate). These checks are research acceptance
// gates, not a calibrated integer correctness probability. No observation
// update, persistent FLOAT mutation, NIS, jitter or eigenvalue clipping occurs.

struct ZhangUserGonsisOutput
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    Eigen::VectorXd mean;
    Eigen::MatrixXd covariance;
    Eigen::VectorXd innovation;
    Eigen::VectorXd constraintResidual;
    double maximumSoftClosureError = 0;
};

enum class ZhangUserGonsisMeasurementKind { CODE, PHASE, REGIONAL_IONO_SD, REGIONAL_TROP_ZTD };

struct ZhangUserGonsisResidualSummary
{
    Eigen::Index rowCount = 0;
    double maximumNormalizedResidual = 0;
    double maximumAbsoluteResidualMetres = 0;
    double minimumFourSigmaMarginMetres = std::numeric_limits<double>::quiet_NaN();
    double sigmaAtMaximumNormalizedResidual = std::numeric_limits<double>::quiet_NaN();
    Eigen::Index rowOfMaximumNormalizedResidual = -1;
};

struct ZhangUserGonsisResidualGate
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    Eigen::VectorXd fixedResidual;
    double maximumNormalizedResidual = 0;
    Eigen::Index codeRows = 0;
    Eigen::Index phaseRows = 0;
    Eigen::Index regionalIonoRows = 0;
    Eigen::Index expectedRegionalIonoRows = -1;
    Eigen::Index regionalTropRows = 0, expectedRegionalTropRows = -1;
    bool fullMeasurementCovarianceValidated = false;
    Eigen::VectorXd measurementSigmas;
    Eigen::VectorXd normalizedResiduals;
    Eigen::VectorXd fourSigmaMarginsMetres;
    ZhangUserGonsisResidualSummary codeSummary, phaseSummary, regionalIonoSummary, regionalTropSummary;
    Eigen::Index firstRejectedRow = -1;
};

// A separate builder metadata marker and the final measurement key must agree.
// Missing/mistyped regional rows are not a legitimate zero-pseudo geometry
// event. This helper intentionally has no KF or process-global dependencies.
inline bool zhangUserGonsisRegionalResidualRowContract(
    bool regionalEnabled, bool ionosphereType, bool builderMarker,
    bool regionalMetreLabel, bool& regionalRow, std::string& reason)
{
    regionalRow = false;
    if (!(builderMarker || regionalMetreLabel || (regionalEnabled && ionosphereType)))
    { reason = "NOT_A_REGIONAL_IONO_ROW"; return true; }
    if (!regionalEnabled || !ionosphereType || !builderMarker || !regionalMetreLabel)
    { reason = "GONSIS_REGIONAL_IONO_ROW_IDENTITY_OR_UNIT_MISMATCH"; return false; }
    regionalRow = true;
    reason = "CURRENT_REGIONAL_SD_L1_METRE_ROW_BOUND";
    return true;
}

struct ZhangUserGonsisPositionGate
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    double floatAverageSigma = std::numeric_limits<double>::quiet_NaN();
    double fixedSigmaNorm = std::numeric_limits<double>::quiet_NaN();
};

namespace zhangUserGonsisDetail
{
// A scale-relative roundoff budget. Zero eigenvalues are allowed in a prior
// or posterior (for example a deterministic constant state). An asymmetric
// or indefinite covariance outside this budget is rejected, never repaired.
inline std::string covarianceFailure(const Eigen::MatrixXd& covariance)
{
    const Eigen::Index dimension = covariance.rows();
    if (dimension <= 0 || covariance.cols() != dimension)
        return "COVARIANCE_DIMENSION_INVALID";
    if (!covariance.allFinite()) return "COVARIANCE_NONFINITE";
    double scale = covariance.cwiseAbs().maxCoeff();
    if (scale == 0) scale = 1;
    const Eigen::MatrixXd scaled = covariance / scale;
    const double matrixNorm = std::max(1.0,
        scaled.cwiseAbs().rowwise().sum().maxCoeff());
    const double tolerance = 128 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(dimension) * matrixNorm;
    if ((scaled - scaled.transpose()).cwiseAbs().maxCoeff() > tolerance)
        return "COVARIANCE_ASYMMETRIC";
    // Symmetrisation is local to the eigenvalue audit, not the returned state.
    const Eigen::MatrixXd symmetric =
        (0.5 * scaled + 0.5 * scaled.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetric);
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
        return "COVARIANCE_EIGEN_FAILED";
    if (eigen.eigenvalues().minCoeff() < -tolerance)
        return "COVARIANCE_NOT_PSD";
    return {};
}
} // namespace zhangUserGonsisDetail

// All constraint innovations come from this ONE immutable snapshot. The caller
// must concatenate the selected rows/values/variances for a retry and call this
// function again with the same frozen FLOAT moments, rather than conditioning
// a rejected candidate. Columns must use the same ordered state keys as FLOAT.
inline ZhangUserGonsisOutput zhangConditionUserGonsisOutput(
    const Eigen::VectorXd& floatMean, const Eigen::MatrixXd& floatCovariance,
    const Eigen::MatrixXd& constraintRows, const Eigen::VectorXd& constraintValues,
    const Eigen::VectorXd& constraintVariances,
    const Eigen::MatrixXd* correlatedConstraintNoise = nullptr)
{
    ZhangUserGonsisOutput out;
    const Eigen::Index dimension = floatMean.size();
    const Eigen::Index constraints = constraintRows.rows();
    if (dimension <= 0 || constraints <= 0 ||
        floatCovariance.rows() != dimension || floatCovariance.cols() != dimension ||
        constraintRows.cols() != dimension || constraintValues.size() != constraints ||
        constraintVariances.size() != constraints)
    { out.reason = "GONSIS_OUTPUT_DIMENSION_MISMATCH"; return out; }
    if (!floatMean.allFinite() || !floatCovariance.allFinite() ||
        !constraintRows.allFinite() || !constraintValues.allFinite() ||
        !constraintVariances.allFinite())
    { out.reason = "GONSIS_OUTPUT_NONFINITE_INPUT"; return out; }
    if ((constraintVariances.array() <= 0).any())
    { out.reason = "GONSIS_CONSTRAINT_VARIANCE_NOT_POSITIVE"; return out; }
    if (correlatedConstraintNoise)
    {
        if (correlatedConstraintNoise->rows()!=constraints || correlatedConstraintNoise->cols()!=constraints ||
            !zhangUserGonsisDetail::covarianceFailure(*correlatedConstraintNoise).empty())
        { out.reason="GONSIS_CORRELATED_CONSTRAINT_NOISE_INVALID"; return out; }
        Eigen::LLT<Eigen::MatrixXd> check(*correlatedConstraintNoise);
        if (check.info()!=Eigen::Success)
        { out.reason="GONSIS_CORRELATED_CONSTRAINT_NOISE_NOT_PD"; return out; }
    }
    const auto priorFailure = zhangUserGonsisDetail::covarianceFailure(floatCovariance);
    if (!priorFailure.empty())
    { out.reason = "GONSIS_PRIOR_" + priorFailure; return out; }

    out.innovation = constraintValues - constraintRows * floatMean;
    const Eigen::MatrixXd cross = floatCovariance * constraintRows.transpose();
    Eigen::MatrixXd innovationCovariance = constraintRows * cross;
    if (correlatedConstraintNoise) innovationCovariance += *correlatedConstraintNoise;
    else innovationCovariance.diagonal() += constraintVariances;
    if (!out.innovation.allFinite() || !cross.allFinite() ||
        !innovationCovariance.allFinite())
    { out.reason = "GONSIS_CONDITIONING_NONFINITE_SYSTEM"; return out; }
    const auto systemFailure =
        zhangUserGonsisDetail::covarianceFailure(innovationCovariance);
    if (!systemFailure.empty())
    { out.reason = "GONSIS_INNOVATION_" + systemFailure; return out; }
    Eigen::LDLT<Eigen::MatrixXd> factor(innovationCovariance);
    if (factor.info() != Eigen::Success || !factor.vectorD().allFinite() ||
        !factor.isPositive() || (factor.vectorD().array() <= 0).any())
    { out.reason = "GONSIS_INNOVATION_NOT_POSITIVE_DEFINITE"; return out; }
    const Eigen::VectorXd solvedInnovation = factor.solve(out.innovation);
    const Eigen::MatrixXd solvedCross = factor.solve(cross.transpose());
    if (factor.info() != Eigen::Success || !solvedInnovation.allFinite() ||
        !solvedCross.allFinite())
    { out.reason = "GONSIS_CONDITIONING_SOLVE_FAILED"; return out; }
    const Eigen::MatrixXd gain = solvedCross.transpose();
    out.mean = floatMean + cross * solvedInnovation;
    const Eigen::MatrixXd complement =
        Eigen::MatrixXd::Identity(dimension, dimension) - gain * constraintRows;
    // Joseph form avoids subtraction of nearly equal prior/reduction matrices.
    if (correlatedConstraintNoise)
        out.covariance = complement * floatCovariance * complement.transpose() +
            gain * (*correlatedConstraintNoise) * gain.transpose();
    else out.covariance = complement * floatCovariance * complement.transpose() +
        gain * constraintVariances.asDiagonal() * gain.transpose();
    if (!out.mean.allFinite() || !out.covariance.allFinite())
    { out.reason = "GONSIS_CONDITIONED_MOMENTS_NONFINITE"; return out; }
    const auto posteriorFailure = zhangUserGonsisDetail::covarianceFailure(out.covariance);
    if (!posteriorFailure.empty())
    { out.reason = "GONSIS_POSTERIOR_" + posteriorFailure; return out; }

    // Soft constraints generally DO NOT close to zero. Their correct remaining
    // residual is R Q^-1 (z-H*FLOAT), not an exact-integer equality/nullspace.
    out.constraintResidual = constraintValues - constraintRows * out.mean;
    Eigen::VectorXd expectedResidual;
    if (correlatedConstraintNoise) expectedResidual=(*correlatedConstraintNoise)*solvedInnovation;
    else expectedResidual=constraintVariances.array()*solvedInnovation.array();
    if (!out.constraintResidual.allFinite() || !expectedResidual.allFinite())
    { out.reason = "GONSIS_SOFT_CLOSURE_NONFINITE"; return out; }
    const long double rounding = 256 * std::numeric_limits<double>::epsilon() *
        static_cast<long double>(std::max(dimension, constraints));
    for (Eigen::Index row = 0; row < constraints; ++row)
    {
        long double scale = std::abs(static_cast<long double>(constraintValues(row))) +
            std::abs(static_cast<long double>(expectedResidual(row)));
        for (Eigen::Index column = 0; column < dimension; ++column)
            scale += std::abs(static_cast<long double>(constraintRows(row, column))) *
                std::abs(static_cast<long double>(out.mean(column)));
        const long double tolerance = rounding * std::max(1.0L, scale);
        const long double error = std::abs(
            static_cast<long double>(out.constraintResidual(row)) - expectedResidual(row));
        if (!std::isfinite(tolerance) || error > std::numeric_limits<double>::max())
        { out.reason = "GONSIS_SOFT_CLOSURE_SCALE_OVERFLOW"; return out; }
        out.maximumSoftClosureError = std::max(out.maximumSoftClosureError,
            static_cast<double>(error));
        if (error > tolerance)
        { out.reason = "GONSIS_SOFT_CLOSURE_FAILED"; return out; }
    }
    out.valid = true;
    out.reason = "GONSIS_SOFT_OUTPUT_VALID";
    return out;
}

// The FLOAT mean, residual and design must describe the same frozen posterior
// and ordered state columns. Variances are the retained measurement R diagonal,
// including the caller's measurement/model noise, NEVER R+H*P*H'. Every retained
// CODE, PHASE and supplied REGIONAL_IONO_SD row is tested, not only rows
// belonging to fixed satellites. Full same-sample R is retained/validated when
// supplied, but the acceptance rule remains rowwise 4 sigma: no inverse-R
// quadratic form, NIS or statistical independence certification is computed.
inline ZhangUserGonsisResidualGate zhangValidateUserGonsisResiduals(
    const Eigen::VectorXd& floatMean, const Eigen::VectorXd& fixedMean,
    const Eigen::MatrixXd& design, const Eigen::VectorXd& floatPostfitResidual,
    const Eigen::VectorXd& measurementVariances,
    const std::vector<ZhangUserGonsisMeasurementKind>& kinds,
    const Eigen::MatrixXd* fullMeasurementCovariance = nullptr,
    Eigen::Index expectedRegionalIonoRows = -1,
    Eigen::Index expectedRegionalTropRows = -1)
{
    ZhangUserGonsisResidualGate out;
    out.expectedRegionalIonoRows = expectedRegionalIonoRows;
    out.expectedRegionalTropRows = expectedRegionalTropRows;
    const Eigen::Index dimension = floatMean.size();
    const Eigen::Index measurements = design.rows();
    if (dimension <= 0 || measurements <= 0 || fixedMean.size() != dimension ||
        design.cols() != dimension || floatPostfitResidual.size() != measurements ||
        measurementVariances.size() != measurements ||
        kinds.size() != static_cast<std::size_t>(measurements))
    { out.reason = "GONSIS_RESIDUAL_DIMENSION_MISMATCH"; return out; }
    if (!floatMean.allFinite() || !fixedMean.allFinite() || !design.allFinite() ||
        !floatPostfitResidual.allFinite() || !measurementVariances.allFinite())
    { out.reason = "GONSIS_RESIDUAL_NONFINITE_INPUT"; return out; }
    if ((measurementVariances.array() <= 0).any())
    { out.reason = "GONSIS_MEASUREMENT_VARIANCE_NOT_POSITIVE"; return out; }
    for (const auto kind : kinds)
    {
        if (kind == ZhangUserGonsisMeasurementKind::CODE) ++out.codeRows;
        else if (kind == ZhangUserGonsisMeasurementKind::PHASE) ++out.phaseRows;
        else if (kind == ZhangUserGonsisMeasurementKind::REGIONAL_IONO_SD) ++out.regionalIonoRows;
        else if (kind == ZhangUserGonsisMeasurementKind::REGIONAL_TROP_ZTD) ++out.regionalTropRows;
        else { out.reason = "GONSIS_MEASUREMENT_KIND_INVALID"; return out; }
    }
    if (out.codeRows == 0 || out.phaseRows == 0)
    { out.reason = "GONSIS_RETAINED_CODE_PHASE_REQUIRED"; return out; }
    if (expectedRegionalIonoRows < -1 ||
        (expectedRegionalIonoRows >= 0 && out.regionalIonoRows != expectedRegionalIonoRows))
    { out.reason = "GONSIS_REGIONAL_IONO_CAPTURE_COUNT_MISMATCH"; return out; }
    if (expectedRegionalTropRows < -1 || (expectedRegionalTropRows>=0 && out.regionalTropRows!=expectedRegionalTropRows))
    { out.reason="GONSIS_REGIONAL_TROP_CAPTURE_COUNT_MISMATCH";return out; }
    if (out.regionalTropRows>0 && !fullMeasurementCovariance)
    { out.reason="GONSIS_REGIONAL_TROP_FULL_MEASUREMENT_COVARIANCE_REQUIRED";return out; }
    if (out.regionalIonoRows > 0 && !fullMeasurementCovariance)
    { out.reason = "GONSIS_REGIONAL_IONO_FULL_MEASUREMENT_COVARIANCE_REQUIRED"; return out; }
    if (fullMeasurementCovariance)
    {
        if (fullMeasurementCovariance->rows() != measurements ||
            fullMeasurementCovariance->cols() != measurements ||
            !zhangUserGonsisDetail::covarianceFailure(*fullMeasurementCovariance).empty())
        { out.reason = "GONSIS_FULL_MEASUREMENT_COVARIANCE_INVALID"; return out; }
        for (Eigen::Index row = 0; row < measurements; ++row)
        {
            const double scale = std::max({1.0, std::abs(measurementVariances(row)),
                std::abs((*fullMeasurementCovariance)(row,row))});
            if (std::abs((*fullMeasurementCovariance)(row,row)-measurementVariances(row)) >
                32*std::numeric_limits<double>::epsilon()*scale)
            { out.reason = "GONSIS_MEASUREMENT_COVARIANCE_DIAGONAL_MISMATCH"; return out; }
        }
        out.fullMeasurementCovarianceValidated = true;
    }
    const Eigen::VectorXd displacement = fixedMean - floatMean;
    out.fixedResidual = floatPostfitResidual - design * displacement;
    if (!displacement.allFinite() || !out.fixedResidual.allFinite())
    { out.reason = "GONSIS_FIXED_RESIDUAL_NONFINITE"; return out; }
    bool failed = false;
    bool regionalFailed = false, tropFailed = false;
    out.measurementSigmas.resize(measurements);
    out.normalizedResiduals.resize(measurements);
    out.fourSigmaMarginsMetres.resize(measurements);
    for (Eigen::Index row = 0; row < measurements; ++row)
    {
        const double sigma = std::sqrt(measurementVariances(row));
        const double magnitude = std::abs(out.fixedResidual(row));
        const double normalized = magnitude / sigma;
        const double margin = 4*sigma-magnitude;
        out.measurementSigmas(row) = sigma;
        out.normalizedResiduals(row) = normalized;
        out.fourSigmaMarginsMetres(row) = margin;
        out.maximumNormalizedResidual = std::max(out.maximumNormalizedResidual, normalized);
        auto& summary = kinds[row] == ZhangUserGonsisMeasurementKind::CODE ? out.codeSummary :
            kinds[row] == ZhangUserGonsisMeasurementKind::PHASE ? out.phaseSummary :
            kinds[row] == ZhangUserGonsisMeasurementKind::REGIONAL_TROP_ZTD ? out.regionalTropSummary : out.regionalIonoSummary;
        if (summary.rowCount == 0 || normalized > summary.maximumNormalizedResidual)
        {
            summary.maximumNormalizedResidual = normalized;
            summary.sigmaAtMaximumNormalizedResidual = sigma;
            summary.rowOfMaximumNormalizedResidual = row;
        }
        summary.minimumFourSigmaMarginMetres = summary.rowCount == 0 ? margin :
            std::min(summary.minimumFourSigmaMarginMetres, margin);
        summary.maximumAbsoluteResidualMetres = std::max(summary.maximumAbsoluteResidualMetres, magnitude);
        ++summary.rowCount;
        const bool rowFailed = magnitude > 4 * sigma; // Exactly 4 sigma is accepted.
        if (rowFailed && out.firstRejectedRow < 0) out.firstRejectedRow = row;
        failed |= rowFailed;
        regionalFailed |= rowFailed && kinds[row] == ZhangUserGonsisMeasurementKind::REGIONAL_IONO_SD;
        tropFailed |= rowFailed && kinds[row] == ZhangUserGonsisMeasurementKind::REGIONAL_TROP_ZTD;
    }
    if (failed)
    { out.reason = tropFailed ? "GONSIS_FIXED_REGIONAL_TROP_4SIGMA_FAILED" :
        regionalFailed ? "GONSIS_FIXED_REGIONAL_IONO_4SIGMA_FAILED" :
        "GONSIS_FIXED_RESIDUAL_4SIGMA_FAILED"; return out; }
    out.valid = true;
    out.reason = "GONSIS_FIXED_RESIDUAL_VALID";
    return out;
}

inline void zhangTraceUserGonsisResidualGate(std::ostream& stream,
    const std::string& context, const ZhangUserGonsisResidualGate& gate)
{
    const auto precision = stream.precision();
    stream << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "\nZHANG_USER_GONSIS_OUTPUT_RESIDUAL " << context
        << " code_rows=" << gate.codeRows << " phase_rows=" << gate.phaseRows
        << " regional_iono_rows=" << gate.regionalIonoRows
        << " expected_regional_iono_rows=" << gate.expectedRegionalIonoRows
        << " regional_trop_rows=" << gate.regionalTropRows << " expected_regional_trop_rows=" << gate.expectedRegionalTropRows
        << " full_same_epoch_R_validated=" << gate.fullMeasurementCovarianceValidated
        << " maximum_normalized_residual=" << gate.maximumNormalizedResidual
        << " first_rejected_row=" << gate.firstRejectedRow
        << " threshold_sigma=4 status=" << gate.reason
        << " residual_model=FINAL_FLOAT_LINEARIZATION_FULL_R_RETAINED_ROWWISE_GATE"
           " units=m observation_update=0 nis_evaluations=0 quadratic_forms=0";
    auto summary = [&](const char* name, const ZhangUserGonsisResidualSummary& value)
    {
        stream << "\nZHANG_USER_GONSIS_OUTPUT_RESIDUAL_KIND " << context
            << " kind=" << name << " rows=" << value.rowCount
            << " maximum_normalized_residual=" << value.maximumNormalizedResidual
            << " maximum_absolute_residual_m=" << value.maximumAbsoluteResidualMetres
            << " minimum_4sigma_margin_m=" << value.minimumFourSigmaMarginMetres
            << " sigma_at_max_normalized_m=" << value.sigmaAtMaximumNormalizedResidual
            << " row_at_max_normalized=" << value.rowOfMaximumNormalizedResidual
            << " threshold_sigma=4 nis_evaluations=0";
    };
    summary("CODE",gate.codeSummary); summary("PHASE",gate.phaseSummary);
    summary("REGIONAL_IONO_SD",gate.regionalIonoSummary);
    summary("REGIONAL_TROP_ZTD",gate.regionalTropSummary);
    stream.precision(precision);
}

inline ZhangUserGonsisPositionGate zhangValidateUserGonsisPosition(
    const Eigen::Vector3d& floatVariances, const Eigen::Vector3d& fixedVariances)
{
    ZhangUserGonsisPositionGate out;
    if (!floatVariances.allFinite() || !fixedVariances.allFinite())
    { out.reason = "GONSIS_POSITION_VARIANCE_NONFINITE"; return out; }
    if ((floatVariances.array() < 0).any() || (fixedVariances.array() < 0).any())
    { out.reason = "GONSIS_POSITION_VARIANCE_NEGATIVE"; return out; }
    const Eigen::Vector3d floatSigma = floatVariances.array().sqrt();
    const Eigen::Vector3d fixedSigma = fixedVariances.array().sqrt();
    out.floatAverageSigma = floatSigma(0) / 3 + floatSigma(1) / 3 + floatSigma(2) / 3;
    out.fixedSigmaNorm = std::hypot(fixedSigma(0), fixedSigma(1), fixedSigma(2));
    if (!std::isfinite(out.floatAverageSigma) || !std::isfinite(out.fixedSigmaNorm))
    { out.reason = "GONSIS_POSITION_SIGMA_NONFINITE"; return out; }
    if (out.floatAverageSigma > 0.5)
    { out.reason = "GONSIS_FLOAT_POSITION_SIGMA_FAILED"; return out; }
    if (out.fixedSigmaNorm >= 0.15)
    { out.reason = "GONSIS_FIXED_POSITION_SIGMA_FAILED"; return out; }
    out.valid = true;
    out.reason = "GONSIS_POSITION_SIGMA_VALID";
    return out;
}

inline ZhangUserGonsisPositionGate zhangValidateUserGonsisPosition(
    const Eigen::MatrixXd& floatCovariance, const Eigen::MatrixXd& fixedCovariance,
    const std::array<Eigen::Index, 3>& positionIndices)
{
    ZhangUserGonsisPositionGate out;
    const Eigen::Index dimension = floatCovariance.rows();
    if (dimension <= 0 || floatCovariance.cols() != dimension ||
        fixedCovariance.rows() != dimension || fixedCovariance.cols() != dimension)
    { out.reason = "GONSIS_POSITION_DIMENSION_MISMATCH"; return out; }
    if (!floatCovariance.allFinite() || !fixedCovariance.allFinite())
    { out.reason = "GONSIS_POSITION_VARIANCE_NONFINITE"; return out; }
    Eigen::Vector3d floatVariances, fixedVariances;
    for (std::size_t component = 0; component < positionIndices.size(); ++component)
    {
        const auto index = positionIndices[component];
        if (index < 0 || index >= dimension)
        { out.reason = "GONSIS_POSITION_INDEX_INVALID"; return out; }
        for (std::size_t previous = 0; previous < component; ++previous)
            if (positionIndices[previous] == index)
            { out.reason = "GONSIS_POSITION_INDEX_DUPLICATED"; return out; }
        floatVariances(component) = floatCovariance(index, index);
        fixedVariances(component) = fixedCovariance(index, index);
    }
    return zhangValidateUserGonsisPosition(floatVariances, fixedVariances);
}
