#pragma once
#include "common/eigenIncluder.hpp"
#include <vector>
#include <string>
// A one-entry covariance factor cache: exact values, state ordering and epoch
// identify the root. Mean/RHS are not part of a covariance decomposition and
// are ALWAYS recomputed by the caller. No same-dimension reuse is permitted.
struct ZhangR49PosteriorWorkspace {
 MatrixXd covariance,squareRoot;
 std::vector<std::string> stateOrder;
 std::string epoch,reason="NOT_EVALUATED";
 bool valid=false;
 std::uint64_t generation=0,hits=0,decompositions=0;
 bool factor(const MatrixXd& p,const std::vector<std::string>& order,const std::string& time) {
  if(valid && epoch==time && stateOrder==order && p.rows()==covariance.rows() &&
     p.cols()==covariance.cols() && (p.array()==covariance.array()).all()) {++hits;return true;}
  valid=false;covariance=p;stateOrder=order;epoch=time;++generation;++decompositions;
  if(p.rows()==0 || p.rows()!=p.cols() || !p.allFinite()){reason="PRIOR_COVARIANCE_EIGENSOLVER_FAILED";return false;}
  Eigen::SelfAdjointEigenSolver<MatrixXd> e((p+p.transpose())*0.5);
  if(e.info()!=Eigen::Success || !e.eigenvalues().allFinite()){reason="PRIOR_COVARIANCE_EIGENSOLVER_FAILED";return false;}
  const double scale=std::max(1.0,e.eigenvalues().cwiseAbs().maxCoeff());
  if(e.eigenvalues().minCoeff() < -1e-9*scale){reason="PRIOR_COVARIANCE_NOT_PSD";return false;}
  squareRoot=e.eigenvectors()*e.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
  valid=true;reason="NONE";return true;
 }
};
