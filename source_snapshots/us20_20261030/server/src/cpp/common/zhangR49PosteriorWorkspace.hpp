#pragma once
#include "common/eigenIncluder.hpp"
#include <vector>
#include <string>
#include "common/lapackWrapper.hpp"
#include "common/zhangArThreading.hpp"

// Lifetime is exactly one AR call, not the observation history. A separate
// immutable FLOAT slot survives rolling WL/L1 Gram updates. Exact covariance
// equality, complete state order and epoch are all required for reuse.
struct ZhangR49FloatRootEntry {
 const MatrixXd& source;
 std::vector<std::string> order;
 std::string epoch;
 MatrixXd covariance,root;
 bool valid=false;
 bool matches(const MatrixXd& p,const std::vector<std::string>& keys,const std::string& time) const {
  return valid && epoch==time && order==keys && p.rows()==covariance.rows() &&
   p.cols()==covariance.cols() && (p.array()==covariance.array()).all();
 }
 bool isSource(const MatrixXd& p,const std::vector<std::string>& keys,const std::string& time) const {
  return epoch==time && order==keys && p.rows()==source.rows() && p.cols()==source.cols() &&
   (p.array()==source.array()).all();
 }
};
inline thread_local ZhangR49FloatRootEntry* zhangR49PinnedFloatRoot=nullptr;
class ZhangR49FloatRootScope {
 ZhangR49FloatRootEntry entry;
 ZhangR49FloatRootEntry* previous=zhangR49PinnedFloatRoot;
public:
 ZhangR49FloatRootScope(const MatrixXd& p,std::vector<std::string> order,std::string time)
  :entry{p,std::move(order),std::move(time)} {zhangR49PinnedFloatRoot=&entry;}
 ~ZhangR49FloatRootScope(){zhangR49PinnedFloatRoot=previous;}
 ZhangR49FloatRootScope(const ZhangR49FloatRootScope&)=delete;
 ZhangR49FloatRootScope& operator=(const ZhangR49FloatRootScope&)=delete;
};
// A one-entry covariance factor cache: exact values, state ordering and epoch
// identify the root. Mean/RHS are not part of a covariance decomposition and
// are ALWAYS recomputed by the caller. No same-dimension reuse is permitted.
struct ZhangR49PosteriorWorkspace {
 MatrixXd covariance,squareRoot;
 std::vector<std::string> stateOrder;
 std::string epoch,reason="NOT_EVALUATED";
 bool valid=false;
 std::uint64_t generation=0,hits=0,decompositions=0;
 std::uint64_t protectedFloatHits=0;
 std::string backend="NOT_EVALUATED";
 // Retain a root we have just conditioned, keyed by the exact covariance
 // formed here. Do not accept a caller-supplied (possibly mismatched) P/root
 // pair. The next factor() still requires identical P, ordering and epoch.
 bool formGram(MatrixXd root,const std::vector<std::string>& order,const std::string& time) {
  valid=false;
  if(root.rows()==0 || root.rows()!=root.cols() || order.size()!=root.rows() || !root.allFinite()) {
   reason="CONDITIONED_ROOT_INVALID";return false;
  }
  MatrixXd gram=root*root.transpose();
  gram=(0.5*(gram+gram.transpose())).eval();
  if(!gram.allFinite()){reason="CONDITIONED_GRAM_INVALID";return false;}
  covariance=std::move(gram);squareRoot=std::move(root);stateOrder=order;epoch=time;
  ++generation;valid=true;reason="NONE";return true;
 }
 bool factor(const MatrixXd& p,const std::vector<std::string>& order,const std::string& time) {
  if(valid && epoch==time && stateOrder==order && p.rows()==covariance.rows() &&
     p.cols()==covariance.cols() && (p.array()==covariance.array()).all()) {++hits;return true;}
  auto* pinned=zhangR49PinnedFloatRoot;
  if(pinned && pinned->matches(p,order,time)) {
   covariance=p;squareRoot=pinned->root;stateOrder=order;epoch=time;
   ++generation;++hits;++protectedFloatHits;valid=true;reason="NONE";backend="PINNED_FLOAT_ROOT";return true;
  }
  valid=false;covariance=p;stateOrder=order;epoch=time;++generation;++decompositions;
  if(p.rows()==0 || p.rows()!=p.cols() || !p.allFinite()){reason="PRIOR_COVARIANCE_EIGENSOLVER_FAILED";return false;}
  VectorXd values;
  MatrixXd vectors;
  bool solved=false;
#ifdef GINAN_USE_OPENBLAS
  if(p.rows()>=64) {
   BlasThreading::ScopedOpenBlasThreadLimit limit(ZhangArThreading::eigensolverThreads());
   vectors=(p+p.transpose())*0.5;values.resize(p.rows());
   solved=LapackWrapper::dsyevd(p.rows(),vectors.data(),values.data())==0 &&
    values.allFinite() && vectors.allFinite();
   backend="LAPACK_DSYEVD";
  }
#endif
  if(!solved) {
   Eigen::SelfAdjointEigenSolver<MatrixXd> e((p+p.transpose())*0.5);
   if(e.info()!=Eigen::Success || !e.eigenvalues().allFinite()){reason="PRIOR_COVARIANCE_EIGENSOLVER_FAILED";return false;}
   values=e.eigenvalues();vectors=e.eigenvectors();backend="EIGEN_SELFADJOINT";
  }
  const double scale=std::max(1.0,values.cwiseAbs().maxCoeff());
  if(values.minCoeff() < -1e-9*scale){reason="PRIOR_COVARIANCE_NOT_PSD";return false;}
  squareRoot=vectors*values.cwiseMax(0).cwiseSqrt().asDiagonal();
  if(pinned && pinned->isSource(p,order,time)) {
   pinned->covariance=p;pinned->root=squareRoot;pinned->valid=true;
  }
  valid=true;reason="NONE";return true;
 }
};
