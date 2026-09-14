#pragma once
#include "common/zhangR50Validation.hpp"
#include <memory>
#include "common/zhangSequentialQuotientShadow.hpp"
#include "common/zhangR47ProductDomain.hpp"
#include "common/zhangR48Marginal.hpp"

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
 ZhangIntegerCandidateNis nis,baseNis,incrementNis;
 double schurIdentityError=0;
};
inline Eigen::MatrixXd zhangR48Numeric(const ZhangExactMatrix& rows,int columns) {
 Eigen::MatrixXd m(rows.size(),columns);
 for(int i=0;i<m.rows();++i)m.row(i)=zhangExactRowToDouble(rows[i]).transpose();return m;
}
inline ZhangR48BridgeResult zhangR48SearchBridge(
 const Eigen::VectorXd& mean,const Eigen::MatrixXd& covariance,
 const ZhangExactMatrix& targets,const ZhangExactMatrix& held,const ZhangExactVector& hv,
 double allocation,double alpha,
 const std::function<ZhangSequentialShadowProposal(const Eigen::VectorXd&,const Eigen::MatrixXd&,double,bool)>& search,
 ZhangR48MarginalWorkspace* sharedWorkspace=nullptr,
 const ZhangR47ProductSearchFrame* domainFrame=nullptr, bool partialShadow=false)
{
 ZhangR48BridgeResult out;const int n=mean.size();
 if(!(allocation>0) || !std::isfinite(allocation)) {
  out.status="NOT_EVALUATED_NO_SEARCH_BUDGET";return out;
 }
 if(targets.size()!=2 || !zhangExactRectangularMatrix(targets,n))return out;
 out.frame=domainFrame?zhangR49CompileTargetOnDomain(targets,*domainFrame,n):
  zhangR47CompileProductSearchFrame(targets,held,hv,n);
 if(!out.frame.valid && out.frame.reason=="TARGET_OUTSIDE_CACHED_DOMAIN")
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

 std::unique_ptr<ZhangR48MarginalWorkspace> local;
 if(!sharedWorkspace){local=std::make_unique<ZhangR48MarginalWorkspace>(mean,covariance);sharedWorkspace=local.get();}
 auto p=zhangR48Numeric(out.frame.projector,n);
 auto marginal=sharedWorkspace->project(p,zhangR48Numeric(held,n),zhangExactRowToDouble(hv));
 if(!marginal.valid){out.status=marginal.failureReason;return out;}
 out.mean=marginal.mean+zhangExactRowToDouble(out.frame.offsets);out.covariance=marginal.covariance;
 out.spent=allocation;
 auto candidate=search(out.mean,out.covariance,allocation,true);
 const int fixedRank=candidate.rows.size();
 if(!candidate.valid || fixedRank<=0 || fixedRank>out.rank || (!partialShadow && fixedRank!=out.rank) || candidate.values.size()!=fixedRank ||
    !zhangExactRectangularMatrix(candidate.rows,out.rank) ||
    zhangExactRowHermiteNormalForm(candidate.rows).basis.size()!=fixedRank ||
    !std::isfinite(candidate.failureProbability) || candidate.failureProbability<0 ||
    zhangR50StatisticalReject(candidate.failureProbability>allocation)) {out.status="BRIDGE_SEARCH_REJECTED";return out;}
 out.rows=zhangExactMultiply(candidate.rows,out.frame.projector);out.values=candidate.values;
 auto offsets=zhangExactMatrixTimesColumn(candidate.rows,out.frame.offsets);
 for(int i=0;i<out.values.size();++i)out.values[i]-=offsets[i];
 auto combined=held;auto cv=hv;
 combined.insert(combined.end(),out.rows.begin(),out.rows.end());cv.insert(cv.end(),out.values.begin(),out.values.end());
 if(!zhangR47AffineIntegerFeasible(combined,cv,n)){out.status="AFFINE_INFEASIBLE";return out;}
 auto h=zhangExactRowHermiteNormalForm(combined,cv);
 auto a=zhangR48Numeric(h.basis,n);
 out.nis=assessZhangIntegerCandidateNis(zhangExactRowToDouble(h.values)-a*mean,a*covariance*a.transpose(),alpha);
 auto hn=zhangR48Numeric(held,n),dn=zhangR48Numeric(out.rows,n);
 if(held.empty()){out.baseNis.valid=true;out.baseNis.nis=0;out.baseNis.rank=0;}
 else out.baseNis=assessZhangIntegerCandidateNis(zhangExactRowToDouble(hv)-hn*mean,hn*covariance*hn.transpose(),alpha);
 const auto dm=sharedWorkspace->project(dn,hn,zhangExactRowToDouble(hv));
 if(dm.valid)out.incrementNis=assessZhangIntegerCandidateNis(zhangExactRowToDouble(out.values)-dm.mean,dm.covariance,alpha);
 if(out.nis.valid && out.baseNis.valid && out.incrementNis.valid)
  out.schurIdentityError=std::abs(out.nis.nis-out.baseNis.nis-out.incrementNis.nis);
 out.accepted=out.nis.valid && zhangR50StatisticalAccept(out.nis.nis<=out.nis.threshold);
 out.status=out.accepted?(partialShadow && fixedRank<out.rank?"PARTIAL_INTEGER_PROGRESS":"BRIDGE_ACCEPTED"):out.nis.status;
 return out;
}

