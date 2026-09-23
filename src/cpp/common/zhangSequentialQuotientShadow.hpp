#pragma once
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangProductRelationSolver.hpp"
#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangIntegerCandidateNis.hpp"
#include <functional>
#include <numeric>
#include <sstream>

struct ZhangSequentialShadowProposal {
    ZhangExactMatrix rows;
    ZhangExactVector values;
    double failureProbability = 1;
    bool valid = false;
};
struct ZhangSequentialShadowResult {
    ZhangExactMatrix rows;
    ZhangExactVector values;
    std::vector<ZhangExactMatrix> acceptedRoundRows;
    std::vector<ZhangExactVector> acceptedRoundValues;
    std::vector<double> acceptedRoundRisk;
    double reservedRisk = 0;
    int rounds = 0, attempts = 0, overlapChecks = 0, overlapConflicts = 0;
    int mergeAttempts = 0, mergeAccepted = 0, overlapCommonRankTotal = 0;
    std::string status = "NOT_STARTED";
    ZhangIntegerCandidateNis lastNis;
    int provisionalRows=0, previousAcceptedRank=0, wholeUnionRank=0;
};

// Pure private workspace: no KF state, writer, ledger, or static registry is
// reachable here. The search callback must likewise operate on local copies.
// Blocks within a round share its pre-round conditional prior, so overlapping
// decisions are checked independently rather than being forced to agree by
// conditioning on the answer under test. The next round conditions the union.
inline ZhangSequentialShadowResult zhangSequentialQuotientShadow(
    const Eigen::VectorXd& mean, const Eigen::MatrixXd& covariance,
    const ZhangExactMatrix& baseline, const ZhangExactVector& baselineValues,
    double budget, double nisAlpha,
    const std::function<ZhangSequentialShadowProposal(
        const Eigen::VectorXd&, const Eigen::MatrixXd&, double, bool)>& search,
    const std::function<void(const std::string&)>& emit,
    int maxRounds=4, int blockSize=8, int maxBlocks=4)
{
    ZhangSequentialShadowResult out;
    out.rows=baseline; out.values=baselineValues;
    const int ambient=mean.size();
    if(ambient<=0 || covariance.rows()!=ambient || covariance.cols()!=ambient ||
       !mean.allFinite() || !covariance.allFinite() || !std::isfinite(budget) ||
       budget<=0 || maxRounds<=0 || blockSize<2 || maxBlocks<=0) {
        out.status="INVALID_INPUT";return out;
    }
    const double allocation=budget/(maxRounds*(maxBlocks+1));
    auto numeric=[](const ZhangExactMatrix& rows,int cols) {
        Eigen::MatrixXd m(rows.size(),cols);
        for(int i=0;i<m.rows();++i) m.row(i)=zhangExactRowToDouble(rows[i]).transpose();
        return m;
    };
    for(int round=0;round<maxRounds;++round) {
        ++out.rounds;
        const double riskBeforeRound=out.reservedRisk;
        const auto quotient=zhangExactAffineIntegerQuotient(out.rows,out.values,ambient);
        if(!quotient.valid) {out.status="INVALID_EXACT_QUOTIENT";break;}
        if(!quotient.quotientRank) {out.status="COMPLETE";break;}
        Eigen::VectorXd conditionalMean=mean;
        Eigen::MatrixXd conditionalCov=covariance;
        if(!out.rows.empty()) {
            const auto condition=zhangConditionPosteriorEffectiveIntegers(mean,covariance,
                numeric(out.rows,ambient),zhangExactRowToDouble(out.values));
            if(!condition.valid) {out.status="CONDITIONING_FAILED";break;}
            conditionalMean=condition.mean;conditionalCov=condition.covariance;
        }
        const auto projector=numeric(quotient.quotientProjector,ambient);
        const Eigen::VectorXd qmean=projector*(conditionalMean-
            zhangExactRowToDouble(quotient.particularSolution));
        const Eigen::MatrixXd qcov=projector*conditionalCov*projector.transpose();
        const int dimension=quotient.quotientRank;
        std::vector<int> order(dimension);std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](int a,int b){return qcov(a,a)<qcov(b,b);});
        std::vector<int> attemptedCoordinates;
        ZhangExactMatrix provisional;ZhangExactVector provisionalValues;
        bool conflict=false;
        auto attempt=[&](const std::vector<int>& coordinates,bool merge,
                         ZhangExactMatrix& rows,ZhangExactVector& values) {
            if(out.reservedRisk+allocation>budget*(1+1e-12)) return false;
            out.reservedRisk+=allocation;++out.attempts;
            Eigen::VectorXd blockMean(coordinates.size());
            Eigen::MatrixXd blockCov(coordinates.size(),coordinates.size());
            for(int i=0;i<coordinates.size();++i) {
                blockMean(i)=qmean(coordinates[i]);
                for(int j=0;j<coordinates.size();++j) blockCov(i,j)=qcov(coordinates[i],coordinates[j]);
            }
            const auto proposal=search(blockMean,blockCov,allocation,merge);
            std::ostringstream event;
            event<<"round="<<round<<" quotient_dimension="<<dimension
                 <<" block_dimension="<<coordinates.size()<<" merge="<<merge
                 <<" proposal_rows="<<proposal.rows.size()<<" allocation="<<allocation
                 <<" nominal_perr="<<proposal.failureProbability;
            emit(event.str());
            if(!proposal.valid || proposal.rows.empty() ||
               proposal.rows.size()!=proposal.values.size() ||
               !std::isfinite(proposal.failureProbability) || proposal.failureProbability<0 ||
               proposal.failureProbability>allocation) return false;
            // A conflict is resolved only by a full joint answer on the merged
            // coordinates. A partial answer may omit the disputed overlap.
            if(merge && proposal.rows.size()!=coordinates.size()) return false;
            ZhangExactMatrix selector(coordinates.size(),ZhangExactVector(dimension));
            for(int i=0;i<coordinates.size();++i) selector[i][coordinates[i]]=1;
            for(const auto& row:proposal.rows) if(row.size()!=coordinates.size()) return false;
            rows=zhangExactMultiply(zhangExactMultiply(proposal.rows,selector),quotient.quotientProjector);
            values=proposal.values;
            const auto offset=zhangExactMatrixTimesColumn(rows,quotient.particularSolution);
            for(int i=0;i<values.size();++i) values[i]+=offset[i];
            return true;
        };
        for(int start=0,block=0;start<dimension && block<maxBlocks;
            start+=blockSize/2,++block) {
            const int stop=std::min(dimension,start+blockSize);
            std::vector<int> coordinates(order.begin()+start,order.begin()+stop);
            attemptedCoordinates.insert(attemptedCoordinates.end(),coordinates.begin(),coordinates.end());
            ZhangExactMatrix rows;ZhangExactVector values;
            if(!attempt(coordinates,false,rows,values)) continue;
            if(!provisional.empty()) {

                auto united=provisional;auto rhs=provisionalValues;
                united.insert(united.end(),rows.begin(),rows.end());rhs.insert(rhs.end(),values.begin(),values.end());
                const int commonRank=zhangExactRowHermiteNormalForm(provisional).basis.size()+
                    zhangExactRowHermiteNormalForm(rows).basis.size()-
                    zhangExactRowHermiteNormalForm(united).basis.size();
                if(commonRank>0) ++out.overlapChecks;
                out.overlapCommonRankTotal+=commonRank;
                emit("round="+std::to_string(round)+" overlap_common_rank="+std::to_string(commonRank));
                // Exact HNF checks common affine functionals in one ambient
                // integer lattice, including mixed rows and opposite signs.
                if(!zhangExactAffineIntegerQuotient(united,rhs,ambient,ZhangExactQuotientWork::FEASIBILITY_ONLY).valid) {
                    ++out.overlapConflicts;conflict=true;break;
                }
            }
            provisional.insert(provisional.end(),rows.begin(),rows.end());
            provisionalValues.insert(provisionalValues.end(),values.begin(),values.end());
            if(stop==dimension) break;
        }
        if(conflict) {
            std::sort(attemptedCoordinates.begin(),attemptedCoordinates.end());
            attemptedCoordinates.erase(std::unique(attemptedCoordinates.begin(),attemptedCoordinates.end()),attemptedCoordinates.end());
            provisional.clear();provisionalValues.clear();++out.mergeAttempts;
            if(!attempt(attemptedCoordinates,true,provisional,provisionalValues)) {
                out.status="OVERLAP_CONFLICT_MERGE_REJECTED";break;
            }
        }
        auto rows=out.rows;auto values=out.values;
        rows.insert(rows.end(),provisional.begin(),provisional.end());
        values.insert(values.end(),provisionalValues.begin(),provisionalValues.end());
        out.provisionalRows=provisional.size();out.previousAcceptedRank=out.rows.size();
        const auto united=zhangExactRowHermiteNormalForm(rows,values);
        out.wholeUnionRank=united.basis.size();
        if(!united.consistent || !zhangExactAffineIntegerQuotient(united.basis,united.values,ambient,ZhangExactQuotientWork::FEASIBILITY_ONLY).valid) {
            out.status="WHOLE_LATTICE_AFFINE_CONFLICT";break;
        }
        if(united.basis.size()<=out.rows.size()) {out.status="NO_RANK_GAIN";break;}
        const auto matrix=numeric(united.basis,ambient);
        const auto nis=assessZhangIntegerCandidateNis(
            zhangExactRowToDouble(united.values)-matrix*mean,
            matrix*covariance*matrix.transpose(),nisAlpha);
        out.lastNis=nis;
        if(!zhangIntegerCandidateAdmissible(nis)) {out.status=nis.valid?"WHOLE_LATTICE_NIS_REJECTED":nis.status;break;}
        out.rows=united.basis;out.values=united.values;
        out.acceptedRoundRows.push_back(out.rows);
        out.acceptedRoundValues.push_back(out.values);
        out.acceptedRoundRisk.push_back(out.reservedRisk-riskBeforeRound);
        if(conflict) ++out.mergeAccepted;
        std::ostringstream event;
        event<<"round="<<round<<" accepted_rank="<<out.rows.size()
             <<" removed_fixed_dimensions="<<ambient-quotient.quotientRank
             <<" joint_nis="<<nis.nis<<" threshold="<<nis.threshold;
        emit(event.str());out.status="ROUND_LIMIT";
    }
    return out;
}
