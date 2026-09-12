#pragma once
#include "common/zhangSequentialQuotientShadow.hpp"
#include "common/zhangR47ProductDomain.hpp"

struct ZhangR48BridgeResult {
 bool accepted=false;
 int rank=0,wlRank=0,l1Rank=0;
 double spent=0;
 std::string status="UNREPRESENTABLE";
 ZhangR47ProductSearchFrame frame;
 ZhangExactMatrix rows;
 ZhangExactVector values;
 Eigen::VectorXd mean;
 Eigen::MatrixXd covariance;
 ZhangIntegerCandidateNis nis;
};
inline Eigen::MatrixXd zhangR48Numeric(const ZhangExactMatrix& rows,int columns) {
 Eigen::MatrixXd m(rows.size(),columns);
 for(int i=0;i<m.rows();++i)m.row(i)=zhangExactRowToDouble(rows[i]).transpose();return m;
}
inline ZhangR48BridgeResult zhangR48SearchBridge(
 const Eigen::VectorXd& mean,const Eigen::MatrixXd& covariance,
 const ZhangExactMatrix& targets,const ZhangExactMatrix& held,const ZhangExactVector& hv,
 double allocation,double alpha,
 const std::function<ZhangSequentialShadowProposal(const Eigen::VectorXd&,const Eigen::MatrixXd&,double,bool)>& search)
{
 ZhangR48BridgeResult out;const int n=mean.size();
 if(targets.size()!=2 || !zhangExactRectangularMatrix(targets,n))return out;
 out.frame=zhangR47CompileProductSearchFrame(targets,held,hv,n);
 if(!out.frame.valid){out.status=out.frame.reason;return out;}
 out.rank=out.frame.searchRank;
 auto residualRank=[&](int row) {
  for(const auto& k:out.frame.affine.kernelBasis) {
   ZhangExactInteger sum=0;for(int c=0;c<out.frame.columns.size();++c)sum+=targets[row][out.frame.columns[c]]*k[c];
   if(sum!=0)return 1;
  }return 0;
 };
 out.wlRank=residualRank(0);out.l1Rank=residualRank(1);
 if(out.rank==0){out.status="ALREADY_PROVEN";return out;}
 if(out.rank>2){out.status="INCONSISTENT_COMPONENT";return out;}
 if(allocation<=0){out.status="RISK_EXHAUSTED";return out;}
 auto mu=mean;auto q=covariance;
 if(!held.empty()) {
  auto c=zhangConditionPosteriorEffectiveIntegers(mean,covariance,zhangR48Numeric(held,n),zhangExactRowToDouble(hv));
  if(!c.valid){out.status=c.failureReason;return out;}mu=c.mean;q=c.covariance;
 }
 auto p=zhangR48Numeric(out.frame.projector,n);
 out.mean=p*mu+zhangExactRowToDouble(out.frame.offsets);out.covariance=p*q*p.transpose();
 out.spent=allocation;
 auto candidate=search(out.mean,out.covariance,allocation,true);
 if(!candidate.valid || candidate.rows.size()!=out.rank || candidate.values.size()!=out.rank ||
    !zhangExactRectangularMatrix(candidate.rows,out.rank) ||
    zhangExactRowHermiteNormalForm(candidate.rows).basis.size()!=out.rank ||
    !std::isfinite(candidate.failureProbability) || candidate.failureProbability<0 ||
    candidate.failureProbability>allocation) {out.status="BRIDGE_SEARCH_REJECTED";return out;}
 out.rows=zhangExactMultiply(candidate.rows,out.frame.projector);out.values=candidate.values;
 auto offsets=zhangExactMatrixTimesColumn(candidate.rows,out.frame.offsets);
 for(int i=0;i<out.values.size();++i)out.values[i]-=offsets[i];
 auto combined=held;auto cv=hv;
 combined.insert(combined.end(),out.rows.begin(),out.rows.end());cv.insert(cv.end(),out.values.begin(),out.values.end());
 if(!zhangR47AffineIntegerFeasible(combined,cv,n)){out.status="AFFINE_INFEASIBLE";return out;}
 auto h=zhangExactRowHermiteNormalForm(combined,cv);
 auto a=zhangR48Numeric(h.basis,n);
 out.nis=assessZhangIntegerCandidateNis(zhangExactRowToDouble(h.values)-a*mean,a*covariance*a.transpose(),alpha);
 out.accepted=out.nis.valid && out.nis.nis<=out.nis.threshold;
 out.status=out.accepted?"BRIDGE_ACCEPTED":out.nis.status;
 return out;
}

