#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <boost/math/distributions/chi_squared.hpp>
#include "common/eigenIncluder.hpp"

struct ZhangIntegerCandidateNis
{
	bool   valid     = false;
	bool   deterministic = false;
	double nis       = std::numeric_limits<double>::quiet_NaN();
	double threshold = std::numeric_limits<double>::quiet_NaN();
	int    rank      = 0;
 double minEigenvalue=0, maxEigenvalue=0, rankTolerance=0, nullResidual=0;
 std::string status="INVALID_INPUT";
};

inline ZhangIntegerCandidateNis assessZhangIntegerCandidateNis(
	const VectorXd& innovation,
	MatrixXd        covariance,
	double          alpha)
{
	ZhangIntegerCandidateNis assessment;
	if (innovation.size() == 0 || covariance.rows() != innovation.size()
	 || covariance.cols() != innovation.size() || !innovation.allFinite() || !covariance.allFinite()
     || !(alpha>0 && alpha<1))
	{
		return assessment;
	}

	covariance = 0.5 * (covariance + covariance.transpose());
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(covariance);
	if (eigenSolver.info() != Eigen::Success
	 || !eigenSolver.eigenvalues().allFinite())
	{
		return assessment;
	}

	const double largestEigenvalue = eigenSolver.eigenvalues().maxCoeff();
	const double rankTolerance = std::max(
		1e-14,
		1e-12 * std::max(0.0, largestEigenvalue));
	assessment.minEigenvalue=eigenSolver.eigenvalues().minCoeff();
 assessment.maxEigenvalue=largestEigenvalue;assessment.rankTolerance=rankTolerance;
 if (assessment.minEigenvalue < -rankTolerance)
	{
        assessment.status="NON_PSD_COVARIANCE";
		return assessment;
	}

	VectorXd coordinates = eigenSolver.eigenvectors().transpose() * innovation;
	double maximumNullInnovation = 0;
	assessment.nis = 0;
	for (int index = 0; index < coordinates.size(); index++)
	{
		const double eigenvalue = eigenSolver.eigenvalues()(index);
		if (eigenvalue > rankTolerance)
		{
			assessment.nis += coordinates(index) * coordinates(index) / eigenvalue;
			assessment.rank++;
		}
		else
		{
			maximumNullInnovation = std::max(
				maximumNullInnovation,
				std::abs(coordinates(index)));
		}
	}
	assessment.nullResidual=maximumNullInnovation;
 assessment.status="DETERMINISTIC_RESIDUAL_CONFLICT";
	if (maximumNullInnovation > 1e-7 || !std::isfinite(assessment.nis))
	{
		return assessment;
	}
	// A posterior that is already fully conditioned on a certified integer has
	// zero covariance in that direction.  It carries no new statistical degree
	// of freedom, but an exactly satisfied affine integer is valid certificate
	// evidence and must not disappear merely because chi-square(0) is undefined.
	// A non-zero innovation in the same null space still fails closed above.
	if (assessment.rank == 0)
	{
		assessment.deterministic = true;
		assessment.threshold = 0;
		assessment.valid = true;
        assessment.status="DETERMINISTIC_ACCEPTED";
		return assessment;
	}

	boost::math::chi_squared distribution(assessment.rank);
	assessment.threshold = quantile(complement(distribution, alpha));
	assessment.valid = std::isfinite(assessment.threshold);
 assessment.status=assessment.valid?(assessment.nis<=assessment.threshold?"ACCEPTED":"WHOLE_NIS_REJECTED"):"INVALID_THRESHOLD";
	return assessment;
}
