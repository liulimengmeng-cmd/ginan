#pragma once
#include "common/zhangRatioOnly.hpp"
#include "common/zhangSequentialQuotientShadow.hpp"
#include "common/zhangR47Candidate.hpp"
#include <iomanip>

// One deterministic pass per round, one optional common-prior rescue. A failed
// block cannot erase a prefix that already passed the complete joint gate.
inline ZhangSequentialShadowResult zhangR48SafePrefix(
 const Eigen::VectorXd& mean,const Eigen::MatrixXd& covariance,
 const ZhangExactMatrix& baseline,const ZhangExactVector& baselineValues,
 double budget,double nisAlpha,
 const std::function<ZhangSequentialShadowProposal(
  const Eigen::VectorXd&,const Eigen::MatrixXd&,double,bool)>& search,
 const std::function<void(const std::string&)>& emit,
 int maxRounds=4,int blockSize=8,int maxBlocks=8)
{
 ZhangSequentialShadowResult out;out.rows=baseline;out.values=baselineValues;
 const int n=mean.size();
 if(n<=0 || covariance.rows()!=n || covariance.cols()!=n || !mean.allFinite() ||
    !covariance.allFinite() || budget<=0 || maxRounds<=0 || blockSize<2 || maxBlocks<=0) {
  out.status="COORDINATE_MISMATCH";return out;
 }
 auto numeric=[](const ZhangExactMatrix& a,int cols) {
  Eigen::MatrixXd m(a.size(),cols);
  for(int i=0;i<m.rows();++i)m.row(i)=zhangExactRowToDouble(a[i]).transpose();return m;
 };
 const double allocation=budget/(maxRounds*(maxBlocks+1));
 auto gate=[&](const ZhangExactMatrix& rows,const ZhangExactVector& values,int oldRank,
               ZhangExactMatrix& accepted,ZhangExactVector& rhs,const std::string& id) {
  out.provisionalRows=rows.size();out.previousAcceptedRank=oldRank;
  if(!zhangR47AffineIntegerFeasible(rows,values,n)) {
   emit(id+" status=AFFINE_INFEASIBLE");return false;
  }
  auto h=zhangExactRowHermiteNormalForm(rows,values);out.wholeUnionRank=h.basis.size();
  if(!h.consistent || h.basis.size()<=oldRank) {emit(id+" status=NO_NEW_RANK");return false;}
  auto a=numeric(h.basis,n);
  auto nis=assessZhangIntegerCandidateNis(zhangExactRowToDouble(h.values)-a*mean,a*covariance*a.transpose(),nisAlpha);
  out.lastNis=nis;
  std::ostringstream s;s<<std::setprecision(17)<<id<<" status="<<nis.status
   <<" whole_union_rank="<<h.basis.size()<<" previous_accepted_rank="<<oldRank
   <<" nis="<<nis.nis<<" threshold="<<nis.threshold<<" stochastic_rank="<<nis.rank
   <<" null_residual="<<nis.nullResidual<<" min_eigenvalue="<<nis.minEigenvalue
   <<" max_eigenvalue="<<nis.maxEigenvalue<<" tolerance="<<nis.rankTolerance;
  emit(s.str());
  if(!nis.valid || zhangRatioStatisticalReject(nis.nis>nis.threshold))return false;
  accepted=std::move(h.basis);rhs=std::move(h.values);return true;
 };
 if(!baseline.empty()) {
  ZhangExactMatrix checked;ZhangExactVector rhs;
  if(!gate(baseline,baselineValues,-1,checked,rhs,"baseline")) {
   out.rows.clear();out.values.clear();out.status="INVALID_BASELINE";return out;
  }
 }
 for(int round=0;round<maxRounds;++round) {
  ++out.rounds;const auto roundRows=out.rows;const auto roundValues=out.values;
  const double riskBefore=out.reservedRisk;
  auto quotient=zhangExactAffineIntegerQuotient(roundRows,roundValues,n);
  if(!quotient.valid){out.status="AFFINE_INFEASIBLE";break;}
  if(!quotient.quotientRank){out.status="COMPLETE";break;}
  Eigen::VectorXd mu=mean;Eigen::MatrixXd q=covariance;
  if(!roundRows.empty()) {
   auto c=zhangConditionPosteriorEffectiveIntegers(mean,covariance,numeric(roundRows,n),zhangExactRowToDouble(roundValues));
   if(!c.valid){out.status=c.failureReason;break;}mu=c.mean;q=c.covariance;
  }
  auto p=numeric(quotient.quotientProjector,n);
  const Eigen::VectorXd qm=p*(mu-zhangExactRowToDouble(quotient.particularSolution));
  const Eigen::MatrixXd qq=p*q*p.transpose();
  std::vector<int> order(qm.size());std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int a,int b){return qq(a,a)<qq(b,b);});
  ZhangExactMatrix safe=roundRows;ZhangExactVector safeValues=roundValues;
  std::vector<int> acceptedCoordinates;bool rescueUsed=false;
  auto attempt=[&](std::vector<int> coords,bool merge,ZhangExactMatrix& rows,ZhangExactVector& values) {
   if(out.reservedRisk+allocation>budget*(1+1e-12))return false;
   out.reservedRisk+=allocation;++out.attempts;
   Eigen::VectorXd bm(coords.size());Eigen::MatrixXd bq(coords.size(),coords.size());
   for(int i=0;i<coords.size();++i){bm(i)=qm(coords[i]);for(int j=0;j<coords.size();++j)bq(i,j)=qq(coords[i],coords[j]);}
   auto candidate=search(bm,bq,allocation,merge);
   std::ostringstream e;e<<std::setprecision(17)<<"round="<<round<<" attempt="<<out.attempts
    <<" merge="<<merge<<" allocation="<<allocation<<" nominal_perr="<<candidate.failureProbability
    <<" proposal_rows="<<candidate.rows.size()<<" mu="<<bm.transpose()<<" covariance=\n"<<bq;
   for(int i=0;i<candidate.rows.size() && i<candidate.values.size();++i) {
    e<<"\nproposal["<<i<<"]=";for(const auto& v:candidate.rows[i])e<<v<<",";e<<" rhs="<<candidate.values[i];
   }
   emit(e.str());
   if(!candidate.valid || candidate.rows.empty() || candidate.rows.size()!=candidate.values.size() ||
      !std::isfinite(candidate.failureProbability) || candidate.failureProbability<0 ||
      zhangRatioStatisticalReject(candidate.failureProbability>allocation) || (merge && candidate.rows.size()!=coords.size()))return false;
   ZhangExactMatrix selector(coords.size(),ZhangExactVector(qm.size()));
   for(int i=0;i<coords.size();++i)selector[i][coords[i]]=1;
   for(const auto& row:candidate.rows)if(row.size()!=coords.size())return false;
   rows=zhangExactMultiply(zhangExactMultiply(candidate.rows,selector),quotient.quotientProjector);
   values=candidate.values;auto offsets=zhangExactMatrixTimesColumn(rows,quotient.particularSolution);
   for(int i=0;i<values.size();++i)values[i]+=offsets[i];
   return true;
  };
  for(int start=0,block=0;start<order.size() && block<maxBlocks;start+=blockSize/2,++block) {
   const int stop=std::min<int>(order.size(),start+blockSize);
   std::vector<int> coords(order.begin()+start,order.begin()+stop);
   ZhangExactMatrix proposed;ZhangExactVector pv;
   if(!attempt(coords,false,proposed,pv))continue;
   auto combined=safe;auto cv=safeValues;
   combined.insert(combined.end(),proposed.begin(),proposed.end());cv.insert(cv.end(),pv.begin(),pv.end());
   const int common=zhangExactRowHermiteNormalForm(safe).basis.size()+
    zhangExactRowHermiteNormalForm(proposed).basis.size()-zhangExactRowHermiteNormalForm(combined).basis.size();
   if(common>0)++out.overlapChecks;out.overlapCommonRankTotal+=common;
   const std::string id="round="+std::to_string(round)+" block="+std::to_string(block);
   ZhangExactMatrix accepted;ZhangExactVector av;
   if(gate(combined,cv,safe.size(),accepted,av,id)) {
    safe=std::move(accepted);safeValues=std::move(av);
    acceptedCoordinates.insert(acceptedCoordinates.end(),coords.begin(),coords.end());
   } else {
    emit(id+" rejected_block_id="+std::to_string(block)+" safe_prefix_retained=1");
    if(!zhangR47AffineIntegerFeasible(combined,cv,n))++out.overlapConflicts;
    // Conditional innovation locates the disputed directions; never accepts.
    if(!safe.empty()) {
     auto c=zhangConditionPosteriorEffectiveIntegers(mean,covariance,numeric(safe,n),zhangExactRowToDouble(safeValues));
     if(c.valid) {auto a=numeric(proposed,n);auto v=zhangExactRowToDouble(pv)-a*c.mean;
      auto cq=(a*c.covariance*a.transpose()).eval();std::ostringstream e;
      e<<id<<" conditional_innovation="<<v.transpose()<<" conditional_covariance=\n"<<cq;emit(e.str());}
    }
    if(!rescueUsed && !acceptedCoordinates.empty()) {
     rescueUsed=true;++out.mergeAttempts;
     auto merged=acceptedCoordinates;merged.insert(merged.end(),coords.begin(),coords.end());
     std::sort(merged.begin(),merged.end());merged.erase(std::unique(merged.begin(),merged.end()),merged.end());
     // Do not grow an unbounded rescue ILS. No answers from this round are
     // conditioned into qm/qq; replacing the prefix also drops its children.
     if(merged.size()<=2*blockSize && attempt(merged,true,proposed,pv)) {
      auto trial=roundRows;auto tv=roundValues;
      trial.insert(trial.end(),proposed.begin(),proposed.end());tv.insert(tv.end(),pv.begin(),pv.end());
      if(gate(trial,tv,safe.size(),accepted,av,id+" rescue=1")) {
       safe=std::move(accepted);safeValues=std::move(av);acceptedCoordinates=merged;++out.mergeAccepted;
      }
     }
    }
   }
   if(stop==order.size())break;
  }
  if(safe.size()<=roundRows.size()){out.status="NO_NEW_RANK";break;}
  out.rows=std::move(safe);out.values=std::move(safeValues);
  out.acceptedRoundRows.push_back(out.rows);out.acceptedRoundValues.push_back(out.values);
  out.acceptedRoundRisk.push_back(out.reservedRisk-riskBefore);
  out.status="SEARCH_LIMIT";
 }
 return out;
}

