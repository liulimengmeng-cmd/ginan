#pragma once

#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>

// A nominal Gaussian consistency screen, not an integer correctness test.
// Use the unchanged pre-update covariance; never repair it or deweight rows here.
struct PppArJointGateResult
{
    bool passed = false;
    double nis = -1;
    double threshold = -1;
    double alpha = -1;
    const char* status = "INVALID_INPUT";
};

inline PppArJointGateResult pppArJointGate(
    const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& innovation,
    double sigmaThreshold)
{
    PppArJointGateResult result;
    const auto rows = innovation.size();
    if (rows <= 0 || covariance.rows() != rows || covariance.cols() != rows ||
        !covariance.allFinite() || !innovation.allFinite() ||
        !std::isfinite(sigmaThreshold) || sigmaThreshold <= 0)
        return result;

    result.alpha = std::erfc(sigmaThreshold / std::sqrt(2.0));
    if (!(result.alpha > 0 && result.alpha < 1))
        return result;
    if (!covariance.isApprox(covariance.transpose(), 1e-10))
    {
        result.status = "ASYMMETRIC_COVARIANCE";
        return result;
    }
    const Eigen::LLT<Eigen::MatrixXd> solver(covariance);
    if (solver.info() != Eigen::Success)
    {
        result.status = "NON_POSITIVE_DEFINITE_COVARIANCE";
        return result;
    }
    const Eigen::VectorXd solved = solver.solve(innovation);
    result.nis = innovation.dot(solved);
    if (!solved.allFinite() || !std::isfinite(result.nis) || result.nis < 0)
    {
        result.status = "INVALID_NIS";
        return result;
    }
    const boost::math::chi_squared distribution(static_cast<double>(rows));
    result.threshold = boost::math::quantile(
        boost::math::complement(distribution, result.alpha));
    result.passed = result.nis <= result.threshold;
    result.status = result.passed ? "PASSED_NOMINAL_SCREEN" : "REJECTED_JOINT_NIS";
    return result;
}
