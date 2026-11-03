#include "common/zhangUserGonsisOutput.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

// Pure numerical contracts. The oracle below is an independently assembled
// information-system inverse and a hand-derived singular-prior case, not a
// repeated call to the production conditioner. No observation filter is run.
extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool accepted, const char* label)
    {
        ++checks;
        if (!accepted) throw std::runtime_error(std::string("FAIL ") + label);
    };
    try
    {
        using Kind = ZhangUserGonsisMeasurementKind;
        Eigen::VectorXd mean(2); mean << 0.4, -0.3;
        Eigen::MatrixXd covariance(2, 2); covariance << 2, 0.6, 0.6, 1.5;
        Eigen::MatrixXd rows(2, 2); rows << 1, -1, 0, 1;
        Eigen::VectorXd values(2); values << 2, -1;
        Eigen::VectorXd variances(2); variances << 0.04, 0.09;
        const auto savedMean = mean;
        const auto savedCovariance = covariance;
        const auto savedRows = rows;
        const auto savedValues = values;
        const auto savedVariances = variances;
        const auto output = zhangConditionUserGonsisOutput(
            mean, covariance, rows, values, variances);
        const Eigen::MatrixXd inversePrior = covariance.inverse();
        const Eigen::MatrixXd inverseNoise = variances.cwiseInverse().asDiagonal();
        const Eigen::MatrixXd expectedCovariance =
            (inversePrior + rows.transpose() * inverseNoise * rows).inverse();
        const Eigen::VectorXd expectedMean = expectedCovariance *
            (inversePrior * mean + rows.transpose() * inverseNoise * values);
        check(output.valid, "soft posterior accepted");
        check((output.mean - expectedMean).norm() < 2e-12,
            "full mean matches independent information inverse");
        check((output.covariance - expectedCovariance).norm() < 2e-12,
            "full covariance matches independent information inverse");
        check(output.constraintResidual.norm() > 1e-4,
            "soft target retains a nonzero residual");
        check((output.innovation - (values - rows * savedMean)).norm() == 0,
            "all innovations use frozen mean");
        check((mean - savedMean).norm() == 0 &&
              (covariance - savedCovariance).norm() == 0 &&
              (rows - savedRows).norm() == 0 && (values - savedValues).norm() == 0 &&
              (variances - savedVariances).norm() == 0,
            "all caller moments and constraints remain unchanged");

        // P=diag(0,2), H=(0,1), R=0.5, z=3 -> x=(7,2.4), P=diag(0,0.4).
        Eigen::VectorXd singularMean(2); singularMean << 7, 0;
        Eigen::MatrixXd singularP = Eigen::MatrixXd::Zero(2, 2);
        singularP(1, 1) = 2;
        Eigen::MatrixXd singularRows(1, 2); singularRows << 0, 1;
        Eigen::VectorXd singularValues(1); singularValues << 3;
        Eigen::VectorXd singularNoise(1); singularNoise << 0.5;
        const auto singular = zhangConditionUserGonsisOutput(
            singularMean, singularP, singularRows, singularValues, singularNoise);
        check(singular.valid && singular.mean(0) == 7 &&
              std::abs(singular.mean(1) - 2.4) < 1e-12 &&
              singular.covariance(0, 0) == 0 &&
              std::abs(singular.covariance(1, 1) - 0.4) < 1e-12,
            "PSD prior and zero posterior direction allowed");

        auto badP = covariance;
        badP(0, 1) += 0.01;
        check(!zhangConditionUserGonsisOutput(mean, badP, rows, values, variances).valid,
            "asymmetric prior fails without repair");
        badP = covariance; badP(0, 0) = -1;
        check(!zhangConditionUserGonsisOutput(mean, badP, rows, values, variances).valid,
            "indefinite prior fails without jitter");
        auto badVariances = variances;
        badVariances(0) = 0;
        check(!zhangConditionUserGonsisOutput(mean, covariance, rows, values, badVariances).valid,
            "zero constraint variance is not a soft constraint");
        badVariances(0) = -0.01;
        check(!zhangConditionUserGonsisOutput(mean, covariance, rows, values, badVariances).valid,
            "negative constraint variance rejected");
        badVariances(0) = std::numeric_limits<double>::quiet_NaN();
        check(!zhangConditionUserGonsisOutput(mean, covariance, rows, values, badVariances).valid,
            "nonfinite constraint variance rejected");
        check(!zhangConditionUserGonsisOutput(mean, covariance,
            Eigen::MatrixXd(0, 2), Eigen::VectorXd(0), Eigen::VectorXd(0)).valid,
            "empty constraint selection rejected");
        check(!zhangConditionUserGonsisOutput(mean, covariance,
            Eigen::MatrixXd::Identity(3, 3), values, variances).valid,
            "constraint dimensions must use frozen column order");

        // Joseph form must preserve a very small positive soft posterior.
        Eigen::VectorXd cancellationMean(1); cancellationMean << 0;
        Eigen::MatrixXd cancellationP(1, 1); cancellationP << 1e20;
        Eigen::MatrixXd cancellationRows(1, 1); cancellationRows << 1;
        Eigen::VectorXd cancellationValues(1); cancellationValues << 1;
        Eigen::VectorXd cancellationNoise(1); cancellationNoise << 1e-8;
        const auto cancellation = zhangConditionUserGonsisOutput(cancellationMean,
            cancellationP, cancellationRows, cancellationValues, cancellationNoise);
        check(cancellation.valid && cancellation.covariance(0, 0) > 0 &&
              std::abs(cancellation.covariance(0, 0) / 1e-8 - 1) < 1e-12,
            "Joseph posterior avoids prior reduction cancellation");

        Eigen::VectorXd floatMean(2); floatMean << 10, 20;
        Eigen::VectorXd fixedMean(2); fixedMean << 11, 18;
        Eigen::MatrixXd design(2, 2); design << 1, 0, 0, 1;
        Eigen::VectorXd floatResidual(2); floatResidual << 1.2, -1.9;
        Eigen::VectorXd measurementVariance(2); measurementVariance << 0.01, 0.04;
        const std::vector<Kind> kinds{Kind::PHASE, Kind::CODE};
        const auto residual = zhangValidateUserGonsisResiduals(floatMean, fixedMean,
            design, floatResidual, measurementVariance, kinds);
        check(residual.valid && residual.codeRows == 1 && residual.phaseRows == 1 &&
              std::abs(residual.fixedResidual(0) - 0.2) < 1e-14 &&
              std::abs(residual.fixedResidual(1) - 0.1) < 1e-14,
            "fixed residual uses displacement from the matching FLOAT snapshot");
        auto failingResidual = floatResidual;
        failingResidual(0) = 1.5;
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean, design,
            failingResidual, measurementVariance, kinds).valid,
            "fixed phase residual above four measurement sigma rejected");
        failingResidual = floatResidual; failingResidual(1) = -1;
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean, design,
            failingResidual, measurementVariance, kinds).valid,
            "fixed code residual above four measurement sigma rejected");

        // At exactly four sigma the source accepts; the next representable
        // residual is rejected. No variance from H*P*H' may inflate this gate.
        Eigen::VectorXd boundary(2); boundary << 4, 4;
        Eigen::VectorXd unitVariance = Eigen::VectorXd::Ones(2);
        check(zhangValidateUserGonsisResiduals(floatMean, floatMean, design,
            boundary, unitVariance, kinds).valid, "four sigma boundary inclusive");
        boundary(0) = std::nextafter(4.0, 5.0);
        check(!zhangValidateUserGonsisResiduals(floatMean, floatMean, design,
            boundary, unitVariance, kinds).valid, "above four sigma boundary rejected");
        auto badMeasurementVariance = measurementVariance;
        badMeasurementVariance(0) = -1;
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean, design,
            floatResidual, badMeasurementVariance, kinds).valid,
            "negative measurement variance rejected");
        badMeasurementVariance(0) = 0;
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean, design,
            floatResidual, badMeasurementVariance, kinds).valid,
            "zero measurement variance rejected");
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean, design,
            floatResidual, measurementVariance, {Kind::CODE, Kind::CODE}).valid,
            "phase support cannot be absent");
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean, design,
            floatResidual, measurementVariance, {Kind::PHASE, Kind::PHASE}).valid,
            "code support cannot be absent");
        check(!zhangValidateUserGonsisResiduals(floatMean, fixedMean,
            Eigen::MatrixXd(0, 2), Eigen::VectorXd(0), Eigen::VectorXd(0), {}).valid,
            "empty residual set fails closed");
        auto nonfiniteFixed = fixedMean;
        nonfiniteFixed(0) = std::numeric_limits<double>::infinity();
        check(!zhangValidateUserGonsisResiduals(floatMean, nonfiniteFixed, design,
            floatResidual, measurementVariance, kinds).valid,
            "nonfinite candidate mean rejected");

        // A stale/premixed mean can falsely appear to pass. The production
        // caller must bind mean, H and residual to ONE snapshot. With FLOAT=0,
        // fixed=1, r_FLOAT=0, sigma=.1 the true residual -1 must fail. Supplying
        // the fixed mean as FLOAT incorrectly reports zero; this is deliberately
        // not used as an oracle and documents the required provenance guard.
        Eigen::VectorXd zeroMean = Eigen::VectorXd::Zero(2);
        Eigen::VectorXd shiftedMean = Eigen::VectorXd::Ones(2);
        Eigen::VectorXd zeroResidual = Eigen::VectorXd::Zero(2);
        Eigen::VectorXd smallVariance = Eigen::VectorXd::Constant(2, 0.01);
        check(!zhangValidateUserGonsisResiduals(zeroMean, shiftedMean, design,
            zeroResidual, smallVariance, kinds).valid,
            "true frozen snapshot rejects an incompatible fixed candidate");
        check(zhangValidateUserGonsisResiduals(shiftedMean, shiftedMean, design,
            zeroResidual, smallVariance, kinds).valid,
            "premixed snapshot requires caller provenance guard");

        Eigen::Vector3d floatPositionVariance = Eigen::Vector3d::Constant(0.25);
        Eigen::Vector3d fixedPositionVariance; fixedPositionVariance << 0.01, 0, 0;
        const auto position = zhangValidateUserGonsisPosition(
            floatPositionVariance, fixedPositionVariance);
        check(position.valid && position.floatAverageSigma == 0.5 &&
              position.fixedSigmaNorm == 0.1, "FLOAT half metre boundary inclusive");
        floatPositionVariance = Eigen::Vector3d::Constant(0.2501);
        check(!zhangValidateUserGonsisPosition(floatPositionVariance,
            fixedPositionVariance).valid, "FLOAT average sigma above threshold rejected");
        floatPositionVariance = Eigen::Vector3d::Constant(0.01);
        fixedPositionVariance << 0.0225, 0, 0;
        check(!zhangValidateUserGonsisPosition(floatPositionVariance,
            fixedPositionVariance).valid, "FIX fifteen centimetre boundary strict");
        fixedPositionVariance << 0.0224, 0, 0;
        check(zhangValidateUserGonsisPosition(floatPositionVariance,
            fixedPositionVariance).valid, "FIX below threshold accepted");
        fixedPositionVariance(0) = -1e-20;
        check(!zhangValidateUserGonsisPosition(floatPositionVariance,
            fixedPositionVariance).valid, "negative coordinate variance rejected");
        fixedPositionVariance(0) = std::numeric_limits<double>::quiet_NaN();
        check(!zhangValidateUserGonsisPosition(floatPositionVariance,
            fixedPositionVariance).valid, "missing coordinate variance rejected");
        const Eigen::MatrixXd floatPositionP = 0.01 * Eigen::MatrixXd::Identity(4, 4);
        const Eigen::MatrixXd fixedPositionP = 0.001 * Eigen::MatrixXd::Identity(4, 4);
        check(zhangValidateUserGonsisPosition(floatPositionP, fixedPositionP,
            std::array<Eigen::Index, 3>{0, 1, 2}).valid, "key-selected position indices");
        check(!zhangValidateUserGonsisPosition(floatPositionP, fixedPositionP,
            std::array<Eigen::Index, 3>{0, 1, -1}).valid, "missing position index rejected");
        check(!zhangValidateUserGonsisPosition(floatPositionP, fixedPositionP,
            std::array<Eigen::Index, 3>{0, 1, 4}).valid, "outside position index rejected");
        check(!zhangValidateUserGonsisPosition(floatPositionP, fixedPositionP,
            std::array<Eigen::Index, 3>{0, 1, 1}).valid, "duplicate position axis rejected");

        std::cout << "USER_GONSIS_OUTPUT_CONTRACT_PASS checks=" << checks << '\n';
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
