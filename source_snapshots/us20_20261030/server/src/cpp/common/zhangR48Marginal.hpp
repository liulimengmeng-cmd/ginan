#pragma once
#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangIntegerCandidateNis.hpp"
#include "common/zhangPhaseTimer.hpp"

// Owner has an immutable actual posterior. Nothing is keyed by size alone.
// H decomposition is reusable for a different RHS; residual/null-space checks
// and target means are recalculated for that RHS. Proofs are never cached here.
struct ZhangR48MarginalWorkspace {
 bool valid=false;
 Eigen::VectorXd mean;
 Eigen::MatrixXd squareRoot;
 struct Conditioner {
  Eigen::MatrixXd rows,inverse,crossRoot,rightBasis,eigenvectors;
  Eigen::VectorXd eigenvalues;
  double tolerance=0;
  bool valid=false;
  std::string reason;
 };
 struct Target {
  Eigen::MatrixXd j,h,covariance;
  Eigen::VectorXd values,mean;
 };
 std::vector<Conditioner> cache;
 std::vector<Target> targets;
 int hits=0,decompositions=0,targetHits=0;
 explicit ZhangR48MarginalWorkspace(const Eigen::VectorXd& m,const Eigen::MatrixXd& p):mean(m) {
  auto w=zhangBuildPosteriorEffectiveWorkspace(p);
  valid=w.valid && m.allFinite() && m.size()==p.rows();
  squareRoot=std::move(w.squareRoot);
 }
 static bool same(const Eigen::MatrixXd& a,const Eigen::MatrixXd& b) {
  return a.rows()==b.rows() && a.cols()==b.cols() && (a.array()==b.array()).all();
 }
 ZhangPosteriorEffectiveConditioningResult project(const Eigen::MatrixXd& j,
  const Eigen::MatrixXd& h,const Eigen::VectorXd& v) {
  ZhangPhaseTimer totalTimer("R48_MARGINAL_TOTAL");
  ZhangPosteriorEffectiveConditioningResult out;
  if(!valid || j.cols()!=mean.size() || h.cols()!=mean.size() || h.rows()!=v.size() ||
     !j.allFinite() || !h.allFinite() || !v.allFinite()) {
   out.failureReason="MARGINAL_COORDINATE_MISMATCH";return out;
  }
  out.inputRows=h.rows();
  for(const auto& t:targets)if(same(t.j,j) && same(t.h,h) && same(t.values,v)) {
   ++targetHits;out.mean=t.mean;out.covariance=t.covariance;
   out.valid=true;out.conditioned=h.rows()>0;out.failureReason="NONE";return out;
  }
  Eigen::MatrixXd projectedRoot=[&](){ZhangPhaseTimer timer("R48_MARGINAL_TARGET_ROOT");return Eigen::MatrixXd(j*squareRoot);}();
  out.mean=j*mean;
  if(h.rows()>0) {
   Conditioner* found=nullptr;
   for(auto& c:cache)if(same(c.rows,h)){found=&c;++hits;break;}
   if(!found) {
    Conditioner c;c.rows=h;
    const auto s=[&](){ZhangPhaseTimer timer("R48_MARGINAL_CONDITIONER_GRAM");
     c.crossRoot=h*squareRoot;return Eigen::MatrixXd(c.crossRoot*c.crossRoot.transpose());}();
    const auto e=[&](){ZhangPhaseTimer timer("R48_MARGINAL_SPECTRUM");
     return ZhangArThreading::eigenSelfAdjoint((s+s.transpose())*0.5);}();++decompositions;
    if(e.info()!=Eigen::Success || !e.eigenvalues().allFinite())c.reason="CONDITIONER_EIGEN_FAILED";
    else {
     c.eigenvalues=e.eigenvalues();c.eigenvectors=e.eigenvectors();
     c.tolerance=std::max(1e-14,1e-12*std::max(0.0,c.eigenvalues.maxCoeff()));
     c.valid=c.eigenvalues.minCoeff()>=-c.tolerance;
     c.reason=c.valid?"NONE":"DETERMINISTIC_OR_NON_PSD_CONDITIONER";
     Eigen::VectorXd inv=Eigen::VectorXd::Zero(h.rows());
     int rank=0;for(int i=0;i<h.rows();++i)if(c.eigenvalues(i)>c.tolerance){inv(i)=1/c.eigenvalues(i);++rank;}
     ZhangPhaseTimer inverseTimer("R48_MARGINAL_INVERSE_AND_RIGHT_BASIS");
     c.inverse=c.eigenvectors*inv.asDiagonal()*c.eigenvectors.transpose();
     // The same retained eigenvectors, order and threshold as before. One
     // level-3 GEMM replaces rank separate level-2 GEMVs/BLAS team launches.
     Eigen::MatrixXd retained(h.rows(),rank);Eigen::VectorXd scale(rank);int col=0;
     for(int i=0;i<h.rows();++i)if(c.eigenvalues(i)>c.tolerance) {
      retained.col(col)=c.eigenvectors.col(i);scale(col++)=1/std::sqrt(c.eigenvalues(i));
     }
     c.rightBasis.noalias()=c.crossRoot.transpose()*retained;
     c.rightBasis=c.rightBasis*scale.asDiagonal();
    }
    if(cache.size()>=8)cache.erase(cache.begin());
    cache.push_back(std::move(c));found=&cache.back();
   }
   if(!found->valid){out.failureReason=found->reason;return out;}
   ZhangPhaseTimer applicationTimer("R48_MARGINAL_APPLY");
   const Eigen::VectorXd residual=v-h*mean;
   const Eigen::VectorXd inBasis=found->eigenvectors.transpose()*residual;
   for(int i=0;i<v.size();++i)if(found->eigenvalues(i)<=found->tolerance && std::abs(inBasis(i))>1e-7) {
    out.failureReason="DETERMINISTIC_OR_NON_PSD_CONDITIONER";return out;
   }
   const Eigen::MatrixXd cross=projectedRoot*found->crossRoot.transpose();
   out.mean+=cross*found->inverse*residual;
   if(found->rightBasis.cols()>0)
    projectedRoot-=(projectedRoot*found->rightBasis)*found->rightBasis.transpose();
  }
  out.covariance=[&](){ZhangPhaseTimer timer("R48_MARGINAL_TARGET_COVARIANCE");return Eigen::MatrixXd(projectedRoot*projectedRoot.transpose());}();
  out.covariance=((out.covariance+out.covariance.transpose())*0.5).eval();
  if(j.rows()==0){out.valid=true;out.failureReason="NONE";return out;}
  const auto check=assessZhangIntegerCandidateNis(Eigen::VectorXd::Zero(j.rows()),out.covariance,1e-6);
  out.valid=check.valid && out.mean.allFinite();out.conditioned=out.valid && h.rows()>0;
  out.failureReason=out.valid?"NONE":check.status;
  if(out.valid) {
   if(targets.size()>=16)targets.erase(targets.begin());
   targets.push_back({j,h,out.covariance,v,out.mean});
  }
  return out;
 }
};
