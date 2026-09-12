#pragma once
#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangIntegerCandidateNis.hpp"

// Immutable posterior ownership is the cache boundary. No reuse across epochs,
// branches, column orders or noise models is possible through this API.
struct ZhangR48MarginalWorkspace {
 bool valid=false;
 Eigen::VectorXd mean;
 Eigen::MatrixXd squareRoot;
 struct Conditioner {
  Eigen::MatrixXd rows, inverse, crossRoot;
  Eigen::VectorXd values, residual;
  bool valid=false;
  std::string reason;
 };
 std::vector<Conditioner> cache;
 int hits=0,decompositions=0;
 explicit ZhangR48MarginalWorkspace(const Eigen::VectorXd& m,const Eigen::MatrixXd& p):mean(m) {
  auto w=zhangBuildPosteriorEffectiveWorkspace(p);
  valid=w.valid && m.allFinite() && m.size()==p.rows();
  squareRoot=std::move(w.squareRoot);
 }
 ZhangPosteriorEffectiveConditioningResult project(const Eigen::MatrixXd& j,
  const Eigen::MatrixXd& h,const Eigen::VectorXd& v) {
  ZhangPosteriorEffectiveConditioningResult out;
  if(!valid || j.cols()!=mean.size() || h.cols()!=mean.size() || h.rows()!=v.size() ||
     !j.allFinite() || !h.allFinite() || !v.allFinite()) {
   out.failureReason="MARGINAL_COORDINATE_MISMATCH";return out;
  }
  const Eigen::MatrixXd root=j*squareRoot;
  out.mean=j*mean;out.covariance=root*root.transpose();out.inputRows=h.rows();
  if(h.rows()==0){out.valid=true;out.failureReason="NONE";return out;}
  Conditioner* found=nullptr;
  for(auto& c:cache)if(c.rows.rows()==h.rows() && c.rows.cols()==h.cols() &&
     (c.rows.array()==h.array()).all() && (c.values.array()==v.array()).all()) {
   found=&c;++hits;break;
  }
  if(!found) {
   Conditioner c;c.rows=h;c.values=v;c.residual=v-h*mean;c.crossRoot=h*squareRoot;
   Eigen::MatrixXd s=c.crossRoot*c.crossRoot.transpose();
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> e((s+s.transpose())*0.5);++decompositions;
   if(e.info()!=Eigen::Success || !e.eigenvalues().allFinite())c.reason="CONDITIONER_EIGEN_FAILED";
   else {
    const double tolerance=std::max(1e-14,1e-12*std::max(0.0,e.eigenvalues().maxCoeff()));
    auto residual=(e.eigenvectors().transpose()*c.residual).eval();
    Eigen::VectorXd inv=Eigen::VectorXd::Zero(v.size());bool ok=e.eigenvalues().minCoeff()>=-tolerance;
    for(int i=0;i<v.size();++i) {
     if(e.eigenvalues()(i)>tolerance)inv(i)=1/e.eigenvalues()(i);
     else if(std::abs(residual(i))>1e-7)ok=false;
    }
    c.valid=ok;c.reason=ok?"NONE":"DETERMINISTIC_OR_NON_PSD_CONDITIONER";
    c.inverse=e.eigenvectors()*inv.asDiagonal()*e.eigenvectors().transpose();
   }
   if(cache.size()>=8)cache.erase(cache.begin());
   cache.push_back(std::move(c));found=&cache.back();
  }
  if(!found->valid){out.failureReason=found->reason;return out;}
  const Eigen::MatrixXd cross=root*found->crossRoot.transpose();
  out.mean+=cross*found->inverse*found->residual;
  out.covariance-=cross*found->inverse*cross.transpose();
  out.covariance=((out.covariance+out.covariance.transpose())*0.5).eval();
  auto check=assessZhangIntegerCandidateNis(Eigen::VectorXd::Zero(j.rows()),out.covariance,1e-6);
  out.valid=check.valid && out.mean.allFinite();
  out.conditioned=out.valid;out.failureReason=out.valid?"NONE":check.status;
  return out;
 }
};

