#pragma once
#include "common/eigenIncluder.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

struct ZhangUserHeldClosure
{
    bool valid = false;
    double residual = 0, covarianceResidual = 0, tolerance = 0;
    std::string reason = "INVALID_INPUT";
};

// Only for a previously committed integer on a continuous physical/product arc.
// This verifies that the posterior still contains the constraint. It does not
// turn a nearly integral FLOAT estimate into a new integer decision.
inline ZhangUserHeldClosure zhangAssessUserHeldClosure(
    const VectorXd& mean, const MatrixXd& covariance,
    const MatrixXd& rows, const VectorXd& integers, bool provenanceValid)
{
    ZhangUserHeldClosure out;
    if (!provenanceValid) { out.reason="MISSING_OR_CHANGED_PROVENANCE"; return out; }
    if (rows.rows()==0 || rows.cols()!=mean.size() || rows.rows()!=integers.size() ||
        covariance.rows()!=mean.size() || covariance.cols()!=mean.size() ||
        !mean.allFinite() || !covariance.allFinite() || !rows.allFinite() || !integers.allFinite()) return out;
    if ((rows.array()-rows.array().round()).abs().maxCoeff()>1e-8 ||
        (integers.array()-integers.array().round()).abs().maxCoeff()>1e-8)
    { out.reason="NON_INTEGER_STATEMENT"; return out; }
    const MatrixXd ap=rows*covariance;
    const MatrixXd bounds=rows.cwiseAbs()*covariance.cwiseAbs();
    // Account for cancellation in A*P using the unprojected summands. Cap the
    // absolute allowance so large initial variances cannot certify a FLOAT row.
    // A propagated posterior also contains accumulated conditioning/Joseph
    // roundoff, not just this final matrix multiply. Use a 1e-10 cycles^2
    // absolute floor (stricter than the exact-update 1e-8 closure contract).
    out.tolerance=std::min(1e-8,std::max(1e-10,
        64*std::numeric_limits<double>::epsilon()*mean.size()*std::max(1.0,bounds.maxCoeff())));
    out.residual=(rows*mean-integers).cwiseAbs().maxCoeff();
    out.covarianceResidual=ap.cwiseAbs().maxCoeff();
    if (out.residual>1e-7) { out.reason="INTEGER_VALUE_CONFLICT"; return out; }
    if (out.covarianceResidual>out.tolerance) { out.reason="POSTERIOR_NOT_CONDITIONED"; return out; }
    out.valid=true; out.reason="COMMITTED_DETERMINISTIC_CLOSURE";
    return out;
}
