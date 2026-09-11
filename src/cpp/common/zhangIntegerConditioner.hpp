#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "common/eigenIncluder.hpp"
#include "common/zhangIarGainAudit.hpp"

struct ZhangIntegerConditionedState
{
	bool     valid = false;
	int      constraintRows = 0;
	int      constraintRank = 0;
	int      covarianceRank = 0;
	double   minimumConstraintEigenvalue =
		std::numeric_limits<double>::quiet_NaN();
	double   minimumSquareRootDiagonal =
		std::numeric_limits<double>::quiet_NaN();
	double   maximumSquareRootDiagonal =
		std::numeric_limits<double>::quiet_NaN();
	double   maximumConstraintResidual =
		std::numeric_limits<double>::quiet_NaN();
	VectorXd mean;
	MatrixXd covariance;
	std::string failureReason;
};

struct ZhangPosteriorEffectiveRowBasis
{
	bool valid = false;
	int inputRows = 0;
	int posteriorRank = 0;
	double varianceTolerance = std::numeric_limits<double>::quiet_NaN();
	std::vector<int> retainedRows;
	std::vector<int> deterministicRows;
	std::vector<double> conditionalVariances;
	std::string failureReason = "NOT_EVALUATED";
};

struct ZhangPosteriorConditionalInnovationAudit
{
	bool valid = false;
	VectorXd residuals;
	double maximumAbsoluteResidual =
		std::numeric_limits<double>::quiet_NaN();
	std::string failureReason = "NOT_EVALUATED";
};

struct ZhangPosteriorEffectiveConditioningResult
{
	bool valid = false;
	bool conditioned = false;
	int inputRows = 0;
	int conditionerRank = 0;
	int deterministicRows = 0;
	double maximumDeterministicResidual = 0;
	VectorXd mean;
	MatrixXd covariance;
	std::string failureReason = "NOT_EVALUATED";
};

/** Select a minimal row basis in the stochastic quotient induced by P.
 *
 * Exact/HNF independence is a certificate property.  Conditioning requires
 * the stronger property rank(A P A')=#rows.  This routine forms the positive
 * eigensquare-root P=L L' and performs deterministic modified Gram-Schmidt on
 * rows of A L.  Rows in the covariance null space remain valid certificates,
 * but are explicitly excluded from the conditioner basis. */
// Call-scoped immutable decomposition. A workspace must never outlive the
// posterior on which it was built; no global cache, pointer key or hash alias.
struct ZhangPosteriorEffectiveWorkspace
{
	bool valid = false;
	MatrixXd squareRoot;
	int covarianceRank = 0;
	double minimumEigenvalue = std::numeric_limits<double>::quiet_NaN();
	double maximumEigenvalue = std::numeric_limits<double>::quiet_NaN();
	double eigenTolerance = std::numeric_limits<double>::quiet_NaN();
	std::string failureReason;
};

/** Audit an exact affine certificate against the same posterior numerical
 * quotient used to build the searchable integer space.
 *
 * Re-diagonalising D P D' with an unrelated tolerance can classify a
 * numerical-null direction as stochastic even though the ambient posterior
 * decomposition has already removed it.  Reusing the immutable workspace
 * makes the two rank decisions identical by construction. */
struct ZhangPosteriorDeterministicCertificateAudit
{
	bool valid = false;
	bool deterministic = false;
	bool affineConsistent = false;
	int inputRows = 0;
	int projectedStochasticRank = 0;
	double posteriorEigenTolerance = std::numeric_limits<double>::quiet_NaN();
	double projectedVarianceTolerance = std::numeric_limits<double>::quiet_NaN();
	double maximumAbsoluteAffineResidual =
		std::numeric_limits<double>::quiet_NaN();
	std::string failureReason = "SNAPSHOT_NOT_EVALUATED";
};

inline ZhangPosteriorEffectiveRowBasis zhangPosteriorEffectiveIntegerRowBasis(
	const ZhangPosteriorEffectiveWorkspace& workspace,
	const MatrixXd& rows,
	double relativeTolerance);

inline ZhangPosteriorEffectiveWorkspace zhangBuildPosteriorEffectiveWorkspace(
	const MatrixXd& covariance)
{
	ZhangPosteriorEffectiveWorkspace result;
	if (covariance.rows() <= 0 || covariance.rows() != covariance.cols() ||
		!covariance.allFinite())
	{
		result.failureReason = "POSTERIOR_EFFECTIVE_BASIS_INPUT_INVALID";
		return result;
	}
	const MatrixXd symmetric = 0.5 * (covariance + covariance.transpose());
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(symmetric);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
	{
		result.failureReason = "POSTERIOR_EFFECTIVE_BASIS_EIGEN_FAILED";
		return result;
	}
	const double maximumEigenvalue = eigen.eigenvalues().maxCoeff();
	const double eigenTolerance = std::max(1e-14,
		std::max(1.0, std::abs(maximumEigenvalue)) *
		std::numeric_limits<double>::epsilon() * covariance.rows() * 64.0);
	if (eigen.eigenvalues().minCoeff() < -eigenTolerance)
	{
		result.failureReason = "POSTERIOR_EFFECTIVE_BASIS_COVARIANCE_NOT_PSD";
		return result;
	}
	result.minimumEigenvalue = eigen.eigenvalues().minCoeff();
	result.maximumEigenvalue = maximumEigenvalue;
	result.eigenTolerance = eigenTolerance;
	std::vector<int> retainedEigenvalues;
	for (int index = 0; index < eigen.eigenvalues().size(); index++)
		if (eigen.eigenvalues()(index) > eigenTolerance)
			retainedEigenvalues.push_back(index);
	MatrixXd squareRoot(covariance.rows(), retainedEigenvalues.size());
	for (int column = 0; column < static_cast<int>(retainedEigenvalues.size()); column++)
	{
		const int index = retainedEigenvalues[column];
		squareRoot.col(column) = eigen.eigenvectors().col(index) *
			std::sqrt(eigen.eigenvalues()(index));
	}
	result.squareRoot = std::move(squareRoot);
	result.covarianceRank = retainedEigenvalues.size();
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangPosteriorDeterministicCertificateAudit
zhangAuditPosteriorDeterministicCertificate(
	const ZhangPosteriorEffectiveWorkspace& workspace,
	const MatrixXd& rows,
	const VectorXd& innovation,
	double affineResidualTolerance = 1e-7,
	double relativeRankTolerance = 1e-10)
{
	ZhangPosteriorDeterministicCertificateAudit result;
	result.inputRows = rows.rows();
	result.posteriorEigenTolerance = workspace.eigenTolerance;
	if (!workspace.valid || rows.rows() <= 0 || rows.rows() != innovation.size() ||
		rows.cols() != workspace.squareRoot.rows() || !rows.allFinite() ||
		!innovation.allFinite() || !(affineResidualTolerance >= 0) ||
		!(relativeRankTolerance > 0))
	{
		result.failureReason = "SNAPSHOT_NUMERICAL_INVALID";
		return result;
	}
	result.maximumAbsoluteAffineResidual = innovation.cwiseAbs().maxCoeff();
	const auto effective = zhangPosteriorEffectiveIntegerRowBasis(
		workspace, rows, relativeRankTolerance);
	if (!effective.valid)
	{
		result.failureReason = "SNAPSHOT_NUMERICAL_INVALID";
		return result;
	}
	result.projectedStochasticRank = effective.posteriorRank;
	result.projectedVarianceTolerance = effective.varianceTolerance;
	result.valid = true;
	if (result.projectedStochasticRank > 0)
	{
		result.failureReason = "SNAPSHOT_HAS_STOCHASTIC_COMPONENT";
		return result;
	}
	result.deterministic = true;
	if (result.maximumAbsoluteAffineResidual > affineResidualTolerance)
	{
		result.failureReason = "SNAPSHOT_AFFINE_RESIDUAL";
		return result;
	}
	result.affineConsistent = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangPosteriorEffectiveRowBasis zhangPosteriorEffectiveIntegerRowBasis(
	const ZhangPosteriorEffectiveWorkspace& workspace,
	const MatrixXd& rows,
	double relativeTolerance = 1e-10)
{
	ZhangPosteriorEffectiveRowBasis result;
	result.inputRows = rows.rows();
	if (!workspace.valid)
	{
		result.failureReason = workspace.failureReason;
		return result;
	}
	if (rows.cols() != workspace.squareRoot.rows() || !rows.allFinite() ||
		!(relativeTolerance > 0))
	{
		result.failureReason = "POSTERIOR_EFFECTIVE_BASIS_INPUT_INVALID";
		return result;
	}
	const MatrixXd whitened = rows * workspace.squareRoot;
	double maximumNorm = 0;
	for (int row = 0; row < whitened.rows(); row++)
		maximumNorm = std::max(maximumNorm, whitened.row(row).norm());
	const double normTolerance = std::max(1e-12,
		relativeTolerance * std::max(1.0, maximumNorm) *
		std::max(rows.rows(), rows.cols()));
	result.varianceTolerance = normTolerance * normTolerance;
	result.conditionalVariances.resize(whitened.rows(),
		std::numeric_limits<double>::quiet_NaN());
	std::vector<VectorXd> orthonormalRows;
	for (int row = 0; row < whitened.rows(); row++)
	{
		VectorXd residual = whitened.row(row).transpose();
		// Re-orthogonalise once to keep the rank decision stable for highly
		// correlated integer rows.
		for (int pass = 0; pass < 2; pass++)
		for (const auto& basis : orthonormalRows)
			residual -= basis.dot(residual) * basis;
		const double variance = residual.squaredNorm();
		result.conditionalVariances[row] = variance;
		if (std::isfinite(variance) && variance > result.varianceTolerance)
		{
			orthonormalRows.push_back(residual / std::sqrt(variance));
			result.retainedRows.push_back(row);
		}
		else result.deterministicRows.push_back(row);
	}
	result.posteriorRank = result.retainedRows.size();
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline ZhangPosteriorEffectiveRowBasis zhangPosteriorEffectiveIntegerRowBasis(
	const MatrixXd& covariance,
	const MatrixXd& rows,
	double relativeTolerance = 1e-10)
{
	return zhangPosteriorEffectiveIntegerRowBasis(
		zhangBuildPosteriorEffectiveWorkspace(covariance), rows, relativeTolerance);
}

/** Audit affine consistency in the covariance quotient.
 *
 * For retained stochastic rows R and every row i this returns
 *
 *   v_i - S_iR S_RR^+ v_R,  S=A P A'.
 *
 * Only residuals of rows classified deterministic by the matching effective
 * basis are consistency constraints.  Testing their raw innovations is wrong
 * whenever a row is deterministic only conditionally on R.
 */
inline ZhangPosteriorConditionalInnovationAudit
zhangPosteriorConditionalInnovationResiduals(
	const MatrixXd& covariance,
	const MatrixXd& rows,
	const VectorXd& innovation,
	const ZhangPosteriorEffectiveRowBasis& effective)
{
	ZhangPosteriorConditionalInnovationAudit result;
	if (!effective.valid || rows.rows() != innovation.size() ||
		rows.cols() != covariance.rows() || covariance.rows() != covariance.cols() ||
		!rows.allFinite() || !covariance.allFinite() || !innovation.allFinite())
	{
		result.failureReason = "CONDITIONAL_INNOVATION_INPUT_INVALID";
		return result;
	}
	result.residuals = innovation;
	if (!effective.retainedRows.empty())
	{
		MatrixXd stochasticRows(effective.retainedRows.size(), rows.cols());
		VectorXd stochasticInnovation(effective.retainedRows.size());
		for (int local = 0;
			local < static_cast<int>(effective.retainedRows.size()); local++)
		{
			stochasticRows.row(local) = rows.row(effective.retainedRows[local]);
			stochasticInnovation(local) = innovation(effective.retainedRows[local]);
		}
		const MatrixXd cross = rows * covariance * stochasticRows.transpose();
		const MatrixXd stochasticCovariance = stochasticRows * covariance *
			stochasticRows.transpose();
		Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(
			0.5 * (stochasticCovariance + stochasticCovariance.transpose()));
		if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
		{
			result.failureReason = "CONDITIONAL_INNOVATION_EIGEN_FAILED";
			return result;
		}
		const double largest = std::max(0.0, eigen.eigenvalues().maxCoeff());
		const double tolerance = std::max(effective.varianceTolerance,
			std::max(1.0, largest) * std::numeric_limits<double>::epsilon() *
				stochasticCovariance.rows() * 64.0);
		VectorXd inverse = VectorXd::Zero(eigen.eigenvalues().size());
		for (int index = 0; index < inverse.size(); index++)
			if (eigen.eigenvalues()(index) > tolerance)
				inverse(index) = 1 / eigen.eigenvalues()(index);
		const MatrixXd pseudoInverse = eigen.eigenvectors() *
			inverse.asDiagonal() * eigen.eigenvectors().transpose();
		result.residuals -= cross * pseudoInverse * stochasticInnovation;
	}
	result.maximumAbsoluteResidual = result.residuals.size() > 0
		? result.residuals.cwiseAbs().maxCoeff() : 0;
	result.valid = result.residuals.allFinite();
	result.failureReason = result.valid ? "NONE" :
		"CONDITIONAL_INNOVATION_NONFINITE";
	return result;
}

inline bool zhangIntegerConditioningInputsValid(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const ZhangIarFunctional& constraints,
	const VectorXd& integers);

/** Independent square-root equality conditioning without P^-1 or normal
 * equations.
 *
 * A rank-revealing eigensquare-root P=L*L' is followed by an orthogonal QR of
 * (A*L)'.  For B=A*L and B'=Q*R*Pi', the minimum-norm constrained increment
 * is L*Q1*R^{-T}*Pi'*innovation, and the conditional covariance is
 * P-(L*Q1)*(L*Q1)'.  This path is deliberately independent of the analytical
 * LDLT solve above and remains defined for positive-semidefinite covariances
 * produced by earlier exact WL constraints. */
inline ZhangIntegerConditionedState
zhangConditionIntegersSquareRootOrthogonal(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const ZhangIarFunctional& constraints,
	const VectorXd& integers)
{
	ZhangIntegerConditionedState result;
	result.constraintRows = constraints.rows();
	if (!zhangIntegerConditioningInputsValid(
			mean, covariance, constraints, integers))
	{
		result.failureReason = "INVALID_INTEGER_SQUARE_ROOT_INPUT";
		return result;
	}

	const MatrixXd symmetric =
		0.5 * (covariance + covariance.transpose());
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(symmetric);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
	{
		result.failureReason = "INTEGER_SQUARE_ROOT_EIGEN_FAILED";
		return result;
	}
	const double maximumEigenvalue = eigen.eigenvalues().maxCoeff();
	const double eigenTolerance = std::max(
		1e-14,
		std::max(1.0, std::abs(maximumEigenvalue))
			* std::numeric_limits<double>::epsilon()
			* covariance.rows() * 64.0);
	if (eigen.eigenvalues().minCoeff() < -eigenTolerance)
	{
		result.failureReason = "INTEGER_SQUARE_ROOT_COVARIANCE_NOT_PSD";
		return result;
	}
	std::vector<int> retained;
	retained.reserve(covariance.rows());
	for (int index = 0; index < eigen.eigenvalues().size(); index++)
	{
		if (eigen.eigenvalues()(index) > eigenTolerance)
		{
			retained.push_back(index);
		}
	}
	result.covarianceRank = retained.size();
	if (result.covarianceRank < constraints.rows())
	{
		result.failureReason = "INTEGER_SQUARE_ROOT_INSUFFICIENT_COVARIANCE_RANK";
		return result;
	}
	MatrixXd squareRoot(covariance.rows(), result.covarianceRank);
	for (int column = 0; column < result.covarianceRank; column++)
	{
		const int eigenIndex = retained[column];
		squareRoot.col(column) = eigen.eigenvectors().col(eigenIndex)
			* std::sqrt(eigen.eigenvalues()(eigenIndex));
	}

	const MatrixXd whitenedConstraints = constraints * squareRoot;
	Eigen::ColPivHouseholderQR<MatrixXd> qr(
		whitenedConstraints.transpose());
	qr.setThreshold(std::max(
		1e-14,
		std::numeric_limits<double>::epsilon()
			* std::max(whitenedConstraints.rows(), whitenedConstraints.cols())
			* 64.0));
	result.constraintRank = qr.rank();
	if (result.constraintRank != constraints.rows())
	{
		result.failureReason = "INTEGER_SQUARE_ROOT_CONSTRAINT_NOT_FULL_ROW_RANK";
		return result;
	}

	const int rows = constraints.rows();
	const MatrixXd upper = qr.matrixR().topLeftCorner(rows, rows)
		.template triangularView<Eigen::Upper>();
	const VectorXd absoluteDiagonal = upper.diagonal().cwiseAbs();
	result.minimumSquareRootDiagonal = absoluteDiagonal.minCoeff();
	result.maximumSquareRootDiagonal = absoluteDiagonal.maxCoeff();
	if (!upper.allFinite() || result.minimumSquareRootDiagonal <= 0)
	{
		result.failureReason = "INTEGER_SQUARE_ROOT_QR_SINGULAR";
		return result;
	}
	MatrixXd selector = MatrixXd::Zero(result.covarianceRank, rows);
	selector.topRows(rows).setIdentity();
	const MatrixXd q1 = qr.householderQ() * selector;
	const VectorXd innovation = integers - constraints * mean;
	const VectorXd permutedInnovation =
		qr.colsPermutation().transpose() * innovation;
	const VectorXd orthogonalIncrement = upper.transpose()
		.template triangularView<Eigen::Lower>()
		.solve(permutedInnovation);
	const MatrixXd constrainedDirections = squareRoot * q1;
	result.mean = mean + constrainedDirections * orthogonalIncrement;
	result.covariance = symmetric
		- constrainedDirections * constrainedDirections.transpose();
	result.covariance = 0.5
		* (result.covariance + result.covariance.transpose());
	result.maximumConstraintResidual =
		(constraints * result.mean - integers).cwiseAbs().maxCoeff();
	if (!result.mean.allFinite() || !result.covariance.allFinite()
	 || !std::isfinite(result.maximumConstraintResidual))
	{
		result.failureReason = "NONFINITE_INTEGER_SQUARE_ROOT_STATE";
		return result;
	}
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

inline bool zhangIntegerConditioningInputsValid(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const ZhangIarFunctional& constraints,
	const VectorXd& integers)
{
	return mean.size() > 0
		&& covariance.rows() == mean.size()
		&& covariance.cols() == mean.size()
		&& constraints.rows() > 0
		&& constraints.cols() == mean.size()
		&& integers.size() == constraints.rows()
		&& mean.allFinite() && covariance.allFinite()
		&& integers.allFinite() && zhangIarSparseAllFinite(constraints);
}

inline bool zhangFactorConstraintCovariance(
	const MatrixXd& covariance,
	const ZhangIarFunctional& constraints,
	MatrixXd& cross,
	MatrixXd& constraintCovariance,
	Eigen::LDLT<MatrixXd>& factor,
	ZhangIntegerConditionedState& result,
	double diagonalVariance = 0)
{
	cross = covariance * constraints.transpose();
	constraintCovariance = constraints * cross;
	constraintCovariance = 0.5
		* (constraintCovariance + constraintCovariance.transpose());
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(constraintCovariance);
	if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
	{
		result.failureReason = "INTEGER_CONSTRAINT_EIGEN_FAILED";
		return false;
	}
	const double maximum = eigen.eigenvalues().maxCoeff();
	const double tolerance = std::max(
		1e-14,
		std::max(1.0, std::abs(maximum))
			* std::numeric_limits<double>::epsilon()
			* constraintCovariance.rows() * 64.0);
	result.constraintRank =
		(eigen.eigenvalues().array() > tolerance).count();
	if (result.constraintRank != constraints.rows())
	{
		result.failureReason = "INTEGER_CONSTRAINT_NOT_FULL_ROW_RANK";
		return false;
	}
	result.minimumConstraintEigenvalue = eigen.eigenvalues().minCoeff();
	MatrixXd innovationCovariance = constraintCovariance;
	innovationCovariance.diagonal().array() += diagonalVariance;
	factor.compute(innovationCovariance);
	if (factor.info() != Eigen::Success || !factor.isPositive())
	{
		result.failureReason = "INTEGER_CONSTRAINT_LDLT_FAILED";
		return false;
	}
	return true;
}

/** Exact one-shot Gaussian conditioning for primitive full-row-rank A*x=z. */
inline ZhangIntegerConditionedState zhangConditionIntegersExact(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const ZhangIarFunctional& constraints,
	const VectorXd& integers)
{
	ZhangIntegerConditionedState result;
	result.constraintRows = constraints.rows();
	if (!zhangIntegerConditioningInputsValid(
			mean, covariance, constraints, integers))
	{
		result.failureReason = "INVALID_INTEGER_CONDITIONING_INPUT";
		return result;
	}
	MatrixXd cross;
	MatrixXd constraintCovariance;
	Eigen::LDLT<MatrixXd> factor;
	if (!zhangFactorConstraintCovariance(
			covariance, constraints, cross, constraintCovariance,
			factor, result))
	{
		return result;
	}
	const VectorXd innovation = integers - constraints * mean;
	result.mean = mean + cross * factor.solve(innovation);
	result.covariance = covariance
		- cross * factor.solve(cross.transpose());
	result.covariance = 0.5
		* (result.covariance + result.covariance.transpose());
	result.maximumConstraintResidual =
		(constraints * result.mean - integers).cwiseAbs().maxCoeff();
	if (!result.mean.allFinite() || !result.covariance.allFinite())
	{
		result.failureReason = "NONFINITE_INTEGER_CONDITIONED_STATE";
		return result;
	}
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

/** Condition only the stochastic quotient of a certified affine row set.
 *
 * A relation that has already been imposed on the posterior has zero
 * conditional variance.  Passing it again to the full-row-rank conditioner
 * incorrectly fails with INTEGER_CONSTRAINT_NOT_FULL_ROW_RANK.  This wrapper
 * keeps such rows as certificates, audits their conditional affine residual,
 * and conditions only the remaining posterior-effective basis. */
inline ZhangPosteriorEffectiveConditioningResult
zhangConditionPosteriorEffectiveIntegers(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const MatrixXd& rows,
	const VectorXd& integers,
	double deterministicRelativeTolerance = 1e-7)
{
	ZhangPosteriorEffectiveConditioningResult result;
	result.inputRows = rows.rows();
	if (mean.size() <= 0 || covariance.rows() != mean.size() ||
		covariance.cols() != mean.size() || rows.cols() != mean.size() ||
		rows.rows() != integers.size() || !mean.allFinite() ||
		!covariance.allFinite() || !rows.allFinite() || !integers.allFinite() ||
		!(deterministicRelativeTolerance > 0))
	{
		result.failureReason = "POSTERIOR_EFFECTIVE_CONDITIONING_INPUT_INVALID";
		return result;
	}
	const auto effective = zhangPosteriorEffectiveIntegerRowBasis(
		covariance, rows);
	if (!effective.valid)
	{
		result.failureReason = effective.failureReason;
		return result;
	}
	result.conditionerRank = effective.posteriorRank;
	result.deterministicRows = effective.deterministicRows.size();
	const auto conditional = zhangPosteriorConditionalInnovationResiduals(
		covariance, rows, integers - rows * mean, effective);
	if (!conditional.valid)
	{
		result.failureReason = conditional.failureReason;
		return result;
	}
	for (const int row : effective.deterministicRows)
	{
		const double residual = std::abs(conditional.residuals(row));
		result.maximumDeterministicResidual = std::max(
			result.maximumDeterministicResidual, residual);
		const double scale = std::max(1.0, std::abs(integers(row)));
		if (residual > deterministicRelativeTolerance * scale)
		{
			result.failureReason = "POSTERIOR_DETERMINISTIC_AFFINE_CONFLICT";
			return result;
		}
	}
	if (effective.retainedRows.empty())
	{
		result.mean = mean;
		result.covariance = covariance;
		result.valid = true;
		result.failureReason = "NONE";
		return result;
	}
	ZhangIarFunctional constraints(
		effective.retainedRows.size(), mean.size());
	VectorXd selected(effective.retainedRows.size());
	for (int local = 0;
		 local < static_cast<int>(effective.retainedRows.size()); local++)
	{
		const int row = effective.retainedRows[local];
		for (int column = 0; column < rows.cols(); column++)
			if (rows(row, column) != 0)
				constraints.insert(local, column) = rows(row, column);
		selected(local) = integers(row);
	}
	constraints.makeCompressed();
	const auto conditioned = zhangConditionIntegersExact(
		mean, covariance, constraints, selected);
	if (!conditioned.valid)
	{
		result.failureReason = conditioned.failureReason;
		return result;
	}
	result.mean = conditioned.mean;
	result.covariance = conditioned.covariance;
	result.conditioned = true;
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}

/** Independent near-zero-noise pseudo-observation re-solve.  Joseph covariance
 * form is deliberately used instead of the exact subtraction above. */
inline ZhangIntegerConditionedState zhangConditionIntegersPseudoObservation(
	const VectorXd& mean,
	const MatrixXd& covariance,
	const ZhangIarFunctional& constraints,
	const VectorXd& integers,
	double sigmaCycles = 1e-8)
{
	ZhangIntegerConditionedState result;
	result.constraintRows = constraints.rows();
	if (!zhangIntegerConditioningInputsValid(
			mean, covariance, constraints, integers)
	 || !std::isfinite(sigmaCycles) || sigmaCycles <= 0)
	{
		result.failureReason = "INVALID_INTEGER_PSEUDO_OBSERVATION_INPUT";
		return result;
	}
	MatrixXd cross;
	MatrixXd constraintCovariance;
	Eigen::LDLT<MatrixXd> factor;
	const double variance = sigmaCycles * sigmaCycles;
	if (!zhangFactorConstraintCovariance(
			covariance, constraints, cross, constraintCovariance,
			factor, result, variance))
	{
		return result;
	}
	const MatrixXd gain = factor.solve(cross.transpose()).transpose();
	const VectorXd innovation = integers - constraints * mean;
	result.mean = mean + gain * innovation;
	const MatrixXd identity = MatrixXd::Identity(mean.size(), mean.size());
	const MatrixXd residualTransform = identity - gain * constraints;
	result.covariance = residualTransform * covariance
		* residualTransform.transpose()
		+ variance * gain * gain.transpose();
	result.covariance = 0.5
		* (result.covariance + result.covariance.transpose());
	result.maximumConstraintResidual =
		(constraints * result.mean - integers).cwiseAbs().maxCoeff();
	if (!result.mean.allFinite() || !result.covariance.allFinite())
	{
		result.failureReason = "NONFINITE_INTEGER_PSEUDO_OBSERVATION_STATE";
		return result;
	}
	result.valid = true;
	result.failureReason = "NONE";
	return result;
}
