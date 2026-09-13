#pragma once
#include "common/zhangIntegerCandidateNis.hpp"

// An empty condition domain has no statistical test and grants no certificate.
// Check before evaluating Eigen products: BLAS requires LDA >= 1 even for M=0.
inline ZhangIntegerCandidateNis zhangR49ConstraintNis(
 const VectorXd& mean, const MatrixXd& covariance,
 const MatrixXd& rows, const VectorXd& values, double alpha)
{
 ZhangIntegerCandidateNis out;
 if(rows.rows()==0) {out.status="NOT_EVALUATED_EMPTY_DOMAIN";return out;}
 if(rows.cols()!=mean.size() || rows.rows()!=values.size() ||
    covariance.rows()!=mean.size() || covariance.cols()!=mean.size())return out;
 return assessZhangIntegerCandidateNis(values-rows*mean,
     rows*covariance*rows.transpose(),alpha);
}
