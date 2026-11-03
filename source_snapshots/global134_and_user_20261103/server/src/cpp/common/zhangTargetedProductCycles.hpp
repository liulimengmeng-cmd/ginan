#pragma once

#include "common/zhangR47ProductDomain.hpp"
#include "common/zhangR48Marginal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

/** A proposal, not an integer decision.  rows are the ORIGINAL root-coordinate
 * integer rows: neither primitive saturation nor gcd normalisation is allowed.
 * In particular 2*n is not silently replaced by n.  No candidate RHS, proof,
 * ratio result, physical parent or publication authority is produced here. */
struct ZhangTargetedCycleBlock
{
    ZhangExactMatrix rows;
    std::vector<int> candidateIndices;
    int exactRankGain=0;
    int targetRankBefore=0, targetRankAfter=0, targetRankReduction=0;
    double conditionalTargetTraceBefore=0, conditionalTargetTraceAfter=0;
    double varianceGain=0, score=0;
    std::string status="NOT_EVALUATED";
};

struct ZhangTargetedCyclePlan
{
    bool valid=false;
    int targetRank=0, heldRank=0, candidateRowsExamined=0;
    int independentCandidates=0, shortlistedCandidates=0;
    // Diagnostics for the hard memory bound: before shortlisting every
    // numerical target has only T0,T1,c; the final joint target has <=26 rows.
    int jointCandidateRows=0, maximumNumericalTargetRows=0;
    int skippedZeroOrDependent=0, skippedDuplicate=0, skippedNumericNull=0;
    std::vector<ZhangTargetedCycleBlock> blocks;
    std::string status="NOT_EVALUATED";
};

namespace ZhangTargetedCycleDetail
{
using Rational=boost::multiprecision::cpp_rational;

/** Rational row-space calculations are used ONLY for independence/dimension.
 * They never supply an integer search basis or change an integer image.  The
 * caller must still compile each proposed block on its actual R47 domain. */
struct RowSpace
{
    std::size_t dimension=0;
    std::vector<std::vector<Rational>> rows;
    std::vector<std::size_t> pivots;
    explicit RowSpace(std::size_t n):dimension(n) {}

    bool add(const ZhangExactVector& input)
    {
        std::vector<Rational> row(dimension);
        for(std::size_t c=0;c<dimension;++c) row[c]=Rational(input[c]);
        for(std::size_t r=0;r<rows.size();++r)
        {
            zhangExactWorkCheckpoint("TARGETED_CYCLE_EXACT_ROW_REDUCTION");
            const Rational multiplier=row[pivots[r]];
            if(multiplier==0) continue;
            for(std::size_t c=pivots[r];c<dimension;++c)
                row[c]-=multiplier*rows[r][c];
        }
        std::size_t pivot=0;
        while(pivot<dimension && row[pivot]==0) ++pivot;
        if(pivot==dimension) return false;
        const Rational divisor=row[pivot];
        for(std::size_t c=pivot;c<dimension;++c) row[c]/=divisor;
        // Insert in pivot order: later inserted rows may start before a
        // previous pivot, but they have already been reduced at every pivot.
        const auto position=std::lower_bound(pivots.begin(),pivots.end(),pivot);
        const auto index=static_cast<std::size_t>(position-pivots.begin());
        pivots.insert(position,pivot);
        rows.insert(rows.begin()+index,std::move(row));
        return true;
    }
};

inline int remainingTargetRank(RowSpace space,const ZhangExactMatrix& targets)
{
    const int before=static_cast<int>(space.rows.size());
    for(const auto& row:targets) space.add(row);
    return static_cast<int>(space.rows.size())-before;
}

inline bool opposite(const ZhangExactVector& a,const ZhangExactVector& b)
{
    for(std::size_t c=0;c<a.size();++c) if(a[c]!=-b[c]) return false;
    return true;
}

inline bool numericRows(const ZhangExactMatrix& input,int dimension,
    Eigen::MatrixXd& output)
{
    output.resize(input.size(),dimension);
    for(std::size_t r=0;r<input.size();++r)
    for(int c=0;c<dimension;++c)
    {
        output(static_cast<Eigen::Index>(r),c)=input[r][c].convert_to<double>();
        if(!std::isfinite(output(static_cast<Eigen::Index>(r),c))) return false;
    }
    return true;
}

/** Match a shared verified affine owner to THIS H,b, not merely dimension or
 * epoch.  A column absent from the owner is not constrained to zero: target
 * and proposal rows retain all their ambient root columns. */
inline bool sameHeldDomain(const ZhangR47ProductSearchFrame& domain,
    const ZhangExactRowHnf& held,int dimension)
{
    if(!domain.valid || !domain.affine || !domain.affine->valid) return false;
    const auto& affine=*domain.affine;
    if(affine.ambientDimension!=static_cast<int>(domain.columns.size()) ||
       affine.deterministicRank!=static_cast<int>(affine.deterministicBasis.size()) ||
       affine.quotientRank!=affine.ambientDimension-affine.deterministicRank ||
       affine.deterministicBasis.size()!=affine.deterministicValues.size() ||
       affine.particularSolution.size()!=domain.columns.size() ||
       !zhangExactRectangularMatrix(affine.deterministicBasis,domain.columns.size()))
        return false;
    std::vector<bool> seen(dimension,false);
    for(int column:domain.columns)
    {
        if(column<0 || column>=dimension || seen[column]) return false;
        seen[column]=true;
    }
    ZhangExactMatrix expanded(affine.deterministicBasis.size(),ZhangExactVector(dimension));
    for(std::size_t r=0;r<expanded.size();++r)
    {
        ZhangExactInteger value=0;
        for(std::size_t c=0;c<domain.columns.size();++c)
        {
            expanded[r][domain.columns[c]]=affine.deterministicBasis[r][c];
            value+=affine.deterministicBasis[r][c]*affine.particularSolution[c];
        }
        if(value!=affine.deterministicValues[r]) return false;
    }
    const auto actual=zhangExactRowHermiteNormalForm(expanded,affine.deterministicValues);
    return actual.consistent && actual.basis==held.basis && actual.values==held.values;
}

struct Candidate
{
    int source=0, numericIndex=0, targetReduction=0;
    double gain=0;
};

/** Full conditional covariance, including all cross terms.  The small block
 * eigensolve scores a proposed block; its eigenvectors NEVER become integer
 * rows.  Singular stochastic directions do not manufacture exact integers. */
inline bool blockGain(const Eigen::MatrixXd& covariance,
    const std::vector<int>& selected,double& gain,double& traceAfter)
{
    Eigen::MatrixXd variance(selected.size(),selected.size());
    Eigen::MatrixXd cross(2,selected.size());
    for(std::size_t i=0;i<selected.size();++i)
    {
        cross.col(i)=covariance.block(0,2+selected[i],2,1);
        for(std::size_t j=0;j<selected.size();++j)
            variance(i,j)=covariance(2+selected[i],2+selected[j]);
    }
    const auto eigen=ZhangArThreading::eigenSelfAdjoint(
        Eigen::MatrixXd((variance+variance.transpose())*0.5));
    if(eigen.info()!=Eigen::Success || !eigen.eigenvalues().allFinite()) return false;
    const double maximum=eigen.eigenvalues().maxCoeff();
    const double tolerance=std::max(1e-14,std::max(0.0,maximum)*1e-12);
    if(eigen.eigenvalues().minCoeff() < -tolerance) return false;
    Eigen::VectorXd inverse=Eigen::VectorXd::Zero(selected.size());
    int rank=0;
    for(int i=0;i<inverse.size();++i)
    {
        if(eigen.eigenvalues()(i)>tolerance)
        {inverse(i)=1/eigen.eigenvalues()(i);++rank;}
        else if((cross*eigen.eigenvectors().col(i)).norm()>
            1e-7*std::max(1.0,cross.norm())) return false;
    }
    if(!rank) return false;
    const Eigen::MatrixXd reduction=cross*eigen.eigenvectors()*inverse.asDiagonal()*
        eigen.eigenvectors().transpose()*cross.transpose();
    const Eigen::MatrixXd before=covariance.topLeftCorner(2,2);
    const Eigen::MatrixXd after=((before-reduction)+(before-reduction).transpose())*0.5;
    if(!after.allFinite()) return false;
    const auto check=ZhangArThreading::eigenSelfAdjoint(after);
    const double targetTolerance=1e-9*std::max(1.0,before.cwiseAbs().maxCoeff());
    if(check.info()!=Eigen::Success || check.eigenvalues().minCoeff() < -targetTolerance)
        return false;
    gain=reduction.trace();traceAfter=std::max(0.0,after.trace());
    return std::isfinite(gain) && gain>=-targetTolerance;
}
} // namespace ZhangTargetedCycleDetail

/** Select bounded auxiliary CURRENT root-coordinate cycles for an uncovered
 * dual-frequency product target.  All candidate sources (unit cycles, WL,
 * accepted-current-network row coefficients, existing product projector,
 * product targets) are supplied by the caller with their physical catalogue.
 * Candidate accepted-network RHS are intentionally NOT an input: proposals
 * must be searched/validated on the actual route and receive current parents.
 *
 * targetRank is the exact dimension rank([H;T])-rank(H).  It does not describe
 * the integer image's invariant factors and cannot authorise a product.  The
 * caller MUST use R47 exact image compilation, unchanged ratio/PAR gates,
 * physical mapping, full affine feasibility, full posterior conditioning and
 * atomic joint route commit for a successful proposal.  On commit, increment
 * the domain version and retry the original target in that new domain.
 *
 * If ws is supplied, its owner must be the immutable root posterior.  Relevant
 * joint Gram blocks are checked against rootQ to reject a stale workspace.
 * Candidate screening uses only T0,T1,c marginals; only a shortlist of AT
 * MOST 24 rows can enter a candidate-by-candidate Gram/marginal.  The full
 * held H,b is never truncated.  No state/covariance or caller-owned exact
 * domain is mutated here. */
inline ZhangTargetedCyclePlan zhangPlanTargetedProductCycles(
    const Eigen::VectorXd& rootMean,const Eigen::MatrixXd& rootQ,
    const ZhangExactMatrix& targets2,const ZhangExactMatrix& held,
    const ZhangExactVector& heldValues,const ZhangExactMatrix& candidateRows,
    int maxBlock=8,int maxCandidates=24,ZhangR48MarginalWorkspace* ws=nullptr,
    const ZhangR47ProductSearchFrame* domain=nullptr)
{
    using namespace ZhangTargetedCycleDetail;
    ZhangTargetedCyclePlan out;
    const int dimension=static_cast<int>(rootMean.size());
    auto fail=[&](const char* reason) {out.status=reason;return out;};
    if(dimension<=0 || rootQ.rows()!=dimension || rootQ.cols()!=dimension ||
       !rootMean.allFinite() || !rootQ.allFinite() || targets2.size()!=2 ||
       held.size()!=heldValues.size() ||
       !zhangExactRectangularMatrix(targets2,dimension) ||
       !zhangExactRectangularMatrix(held,dimension) ||
       !zhangExactRectangularMatrix(candidateRows,dimension) ||
       maxBlock<=0 || maxBlock>32 || maxCandidates<=0 || maxCandidates>24)
        return fail("TARGETED_CYCLE_INPUT_INVALID");
    const double covarianceScale=std::max(1.0,rootQ.cwiseAbs().maxCoeff());
    if((rootQ-rootQ.transpose()).cwiseAbs().maxCoeff()>1e-10*covarianceScale)
        return fail("TARGETED_CYCLE_COVARIANCE_ASYMMETRIC");
    const auto canonical=zhangExactRowHermiteNormalForm(held,heldValues);
    if(!canonical.consistent) return fail("TARGETED_CYCLE_HELD_CONFLICT");
    if(domain)
    {
        if(!sameHeldDomain(*domain,canonical,dimension))
            return fail("TARGETED_CYCLE_DOMAIN_MISMATCH");
    }
    else if(!zhangExactAffineIntegerQuotient(canonical.basis,canonical.values,
        dimension,ZhangExactQuotientWork::FEASIBILITY_ONLY).valid)
        return fail("TARGETED_CYCLE_HELD_INTEGER_INFEASIBLE");
    RowSpace heldSpace(dimension);
    for(const auto& row:canonical.basis) heldSpace.add(row);
    out.heldRank=static_cast<int>(heldSpace.rows.size());
    out.targetRank=remainingTargetRank(heldSpace,targets2);
    if(!out.targetRank)
    {out.valid=true;out.status="TARGET_ALREADY_DETERMINED";return out;}

    std::vector<Candidate> candidates;
    for(std::size_t index=0;index<candidateRows.size();++index)
    {
        ++out.candidateRowsExamined;
        auto independent=heldSpace;
        if(!independent.add(candidateRows[index]))
        {++out.skippedZeroOrDependent;continue;}
        bool duplicate=false;
        for(const auto& previous:candidates)
            if(candidateRows[previous.source]==candidateRows[index] ||
               opposite(candidateRows[previous.source],candidateRows[index]))
            {duplicate=true;break;}
        if(duplicate) {++out.skippedDuplicate;continue;}
        Candidate candidate;
        candidate.source=static_cast<int>(index);
        candidate.targetReduction=out.targetRank-remainingTargetRank(independent,targets2);
        candidates.push_back(candidate);
    }
    out.independentCandidates=static_cast<int>(candidates.size());
    if(candidates.empty())
    {out.valid=true;out.status="NO_INDEPENDENT_AUXILIARY_CYCLES";return out;}
    Eigen::MatrixXd numericTargets,numericHeld;
    if(!numericRows(targets2,dimension,numericTargets) ||
       !numericRows(canonical.basis,dimension,numericHeld))
        return fail("TARGETED_CYCLE_NUMERIC_COEFFICIENT_OVERFLOW");
    Eigen::VectorXd values(canonical.values.size());
    for(std::size_t i=0;i<canonical.values.size();++i)
    {
        values(i)=canonical.values[i].convert_to<double>();
        if(!std::isfinite(values(i))) return fail("TARGETED_CYCLE_NUMERIC_RHS_OVERFLOW");
    }
    std::unique_ptr<ZhangR48MarginalWorkspace> owned;
    if(!ws) {owned=std::make_unique<ZhangR48MarginalWorkspace>(rootMean,rootQ);ws=owned.get();}
    if(!ws->valid || ws->mean.size()!=dimension || ws->squareRoot.rows()!=dimension ||
       !(ws->mean.array()==rootMean.array()).all())
        return fail("TARGETED_CYCLE_WORKSPACE_MISMATCH");
    // This base includes the COMPLETE held system, but NO candidate pool.
    // Its H/T Gram is built once.  A candidate's H/T cross terms and scalar
    // variance are checked individually below, without a pool-size square.
    Eigen::MatrixXd base(numericTargets.rows()+numericHeld.rows(),dimension);
    base.topRows(numericTargets.rows())=numericTargets;
    base.bottomRows(numericHeld.rows())=numericHeld;
    const Eigen::MatrixXd baseTimesQ=base*rootQ;
    const Eigen::MatrixXd baseGram=baseTimesQ*base.transpose();
    const Eigen::MatrixXd baseRoot=base*ws->squareRoot;
    const Eigen::MatrixXd cachedBaseGram=baseRoot*baseRoot.transpose();
    const double baseScale=std::max(1.0,baseGram.cwiseAbs().maxCoeff());
    if(!baseGram.allFinite() || !cachedBaseGram.allFinite() ||
       (baseGram-cachedBaseGram).cwiseAbs().maxCoeff()>1e-8*baseScale)
        return fail("TARGETED_CYCLE_WORKSPACE_COVARIANCE_MISMATCH");
    const auto targetMarginal=ws->project(numericTargets,numericHeld,values);
    out.maximumNumericalTargetRows=2;
    if(!targetMarginal.valid || targetMarginal.covariance.rows()!=2 ||
       !targetMarginal.covariance.allFinite())
        return fail("TARGETED_CYCLE_CONDITIONAL_MARGINAL_INVALID");
    const double traceBefore=std::max(0.0,targetMarginal.covariance.trace());
    const double gainTolerance=1e-12*std::max(1.0,traceBefore);
    auto better=[](const Candidate& a,const Candidate& b) {
        if(a.targetReduction!=b.targetReduction) return a.targetReduction>b.targetReduction;
        if(a.gain!=b.gain) return a.gain>b.gain;
        return a.source<b.source;
    };
    std::vector<Candidate> useful;
    for(auto candidate:candidates)
    {
        const auto& exactRow=candidateRows[candidate.source];
        Eigen::MatrixXd numericCandidate;
        if(!numericRows({exactRow},dimension,numericCandidate))
            return fail("TARGETED_CYCLE_NUMERIC_COEFFICIENT_OVERFLOW");
        std::vector<int> nonzero;
        for(int c=0;c<dimension;++c) if(exactRow[c]!=0) nonzero.push_back(c);
        Eigen::VectorXd cross=Eigen::VectorXd::Zero(base.rows());
        Eigen::RowVectorXd candidateRoot=Eigen::RowVectorXd::Zero(ws->squareRoot.cols());
        double variance=0;
        for(int c:nonzero)
        {
            const double coefficient=numericCandidate(0,c);
            cross.noalias()+=coefficient*baseTimesQ.col(c);
            candidateRoot.noalias()+=coefficient*ws->squareRoot.row(c);
            for(int d:nonzero)
                variance+=coefficient*rootQ(c,d)*numericCandidate(0,d);
        }
        const Eigen::VectorXd cachedCross=baseRoot*candidateRoot.transpose();
        const double cachedVariance=candidateRoot.squaredNorm();
        const double scalarScale=std::max(baseScale,std::max(std::abs(variance),
            cross.cwiseAbs().maxCoeff()));
        if(!cross.allFinite() || !cachedCross.allFinite() || !std::isfinite(variance) ||
           !std::isfinite(cachedVariance) ||
           (cross-cachedCross).cwiseAbs().maxCoeff()>1e-8*scalarScale ||
           std::abs(variance-cachedVariance)>1e-8*scalarScale)
            return fail("TARGETED_CYCLE_WORKSPACE_COVARIANCE_MISMATCH");
        // Exactly one candidate enters this marginal.  Reuse the complete H
        // decomposition, including target-H and candidate-H covariance.
        Eigen::MatrixXd singleTarget(3,dimension);
        singleTarget.topRows(2)=numericTargets;
        singleTarget.bottomRows(1)=numericCandidate;
        const auto singleMarginal=ws->project(singleTarget,numericHeld,values);
        out.maximumNumericalTargetRows=3;
        if(!singleMarginal.valid || singleMarginal.covariance.rows()!=3 ||
           !singleMarginal.covariance.allFinite())
            return fail("TARGETED_CYCLE_CONDITIONAL_MARGINAL_INVALID");
        double after=0;
        if(!blockGain(singleMarginal.covariance,{0},candidate.gain,after))
        {++out.skippedNumericNull;continue;}
        if(candidate.targetReduction>0 || candidate.gain>gainTolerance)
        {
            useful.push_back(candidate);
            std::stable_sort(useful.begin(),useful.end(),better);
            if(useful.size()>static_cast<std::size_t>(maxCandidates)) useful.pop_back();
        }
    }
    out.shortlistedCandidates=static_cast<int>(useful.size());
    if(useful.empty())
    {out.valid=true;out.status="NO_TARGETED_AUXILIARY_CYCLE_GAIN";return out;}

    // This is the FIRST candidate-by-candidate numeric projection, AFTER the
    // hard shortlist bound.  It contains at most 24 candidates plus two target
    // rows, irrespective of the size of the caller's proposal pool.
    ZhangExactMatrix joint=targets2;
    for(std::size_t i=0;i<useful.size();++i)
    {
        useful[i].numericIndex=static_cast<int>(i);
        joint.push_back(candidateRows[useful[i].source]);
    }
    Eigen::MatrixXd numericJoint;
    if(!numericRows(joint,dimension,numericJoint))
        return fail("TARGETED_CYCLE_NUMERIC_COEFFICIENT_OVERFLOW");
    out.jointCandidateRows=static_cast<int>(useful.size());
    out.maximumNumericalTargetRows=static_cast<int>(numericJoint.rows());
    // Candidate-candidate cross terms are also verified for the SHORTLIST.
    // They were not needed for scalar screening and must not be guessed from
    // diagonal variances when a joint block is scored.
    Eigen::MatrixXd required(numericJoint.rows()+numericHeld.rows(),dimension);
    required.topRows(numericJoint.rows())=numericJoint;
    required.bottomRows(numericHeld.rows())=numericHeld;
    const Eigen::MatrixXd rootGram=required*rootQ*required.transpose();
    const Eigen::MatrixXd projectedRoot=required*ws->squareRoot;
    const Eigen::MatrixXd cachedGram=projectedRoot*projectedRoot.transpose();
    const double gramScale=std::max(1.0,rootGram.cwiseAbs().maxCoeff());
    if(!rootGram.allFinite() || !cachedGram.allFinite() ||
       (rootGram-cachedGram).cwiseAbs().maxCoeff()>1e-8*gramScale)
        return fail("TARGETED_CYCLE_WORKSPACE_COVARIANCE_MISMATCH");
    const auto marginal=ws->project(numericJoint,numericHeld,values);
    if(!marginal.valid || marginal.covariance.rows()!=numericJoint.rows() ||
       !marginal.covariance.allFinite())
        return fail("TARGETED_CYCLE_CONDITIONAL_MARGINAL_INVALID");
    const auto& covariance=marginal.covariance;
    auto makeBlock=[&](const std::vector<int>& positions) {
        ZhangTargetedCycleBlock block;
        RowSpace space=heldSpace;
        std::vector<int> numericIndices;
        for(int position:positions)
        {
            const auto& candidate=useful[position];
            if(!space.add(candidateRows[candidate.source])) return block;
            block.rows.push_back(candidateRows[candidate.source]);
            block.candidateIndices.push_back(candidate.source);
            numericIndices.push_back(candidate.numericIndex);
        }
        if(!blockGain(covariance,numericIndices,block.varianceGain,
            block.conditionalTargetTraceAfter)) {block.rows.clear();return block;}
        block.exactRankGain=static_cast<int>(positions.size());
        block.targetRankBefore=out.targetRank;
        block.targetRankAfter=remainingTargetRank(space,targets2);
        block.targetRankReduction=out.targetRank-block.targetRankAfter;
        block.conditionalTargetTraceBefore=traceBefore;
        block.score=block.targetRankReduction+
            std::max(0.0,block.varianceGain)/std::max(1e-14,traceBefore);
        block.status=block.targetRankReduction>0 ? "EXACT_TARGET_RANK_REDUCTION" :
            "CONDITIONAL_TARGET_VARIANCE_GAIN";
        return block;
    };
    for(int i=0;i<static_cast<int>(useful.size());++i)
    {
        auto block=makeBlock({i});
        if(!block.rows.empty()) out.blocks.push_back(std::move(block));
    }
    std::vector<int> selected;
    int selectedTargetReduction=0;
    double selectedGain=0;
    while(selected.size()<static_cast<std::size_t>(maxBlock))
    {
        int best=-1;
        ZhangTargetedCycleBlock bestBlock;
        for(int i=0;i<static_cast<int>(useful.size());++i)
        {
            if(std::find(selected.begin(),selected.end(),i)!=selected.end()) continue;
            auto proposal=selected;proposal.push_back(i);
            auto block=makeBlock(proposal);
            if(block.rows.empty()) continue;
            if(best<0 || block.targetRankReduction>bestBlock.targetRankReduction ||
               (block.targetRankReduction==bestBlock.targetRankReduction &&
                block.varianceGain>bestBlock.varianceGain))
            {best=i;bestBlock=std::move(block);}
        }
        if(best<0 || (bestBlock.targetRankReduction<=selectedTargetReduction &&
            bestBlock.varianceGain<=selectedGain+gainTolerance)) break;
        selected.push_back(best);
        selectedTargetReduction=bestBlock.targetRankReduction;
        selectedGain=bestBlock.varianceGain;
        const bool targetDetermined=bestBlock.targetRankAfter==0;
        if(selected.size()>1) out.blocks.push_back(std::move(bestBlock));
        if(targetDetermined) break;
    }
    std::stable_sort(out.blocks.begin(),out.blocks.end(),[](const auto& a,const auto& b) {
        if(a.targetRankReduction!=b.targetRankReduction)
            return a.targetRankReduction>b.targetRankReduction;
        if(a.varianceGain!=b.varianceGain) return a.varianceGain>b.varianceGain;
        if(a.rows.size()!=b.rows.size()) return a.rows.size()<b.rows.size();
        return a.candidateIndices<b.candidateIndices;
    });
    out.valid=true;
    out.status=out.blocks.empty()?"NO_TARGETED_AUXILIARY_CYCLE_GAIN":"TARGETED_CYCLE_PROPOSALS_READY";
    return out;
}
