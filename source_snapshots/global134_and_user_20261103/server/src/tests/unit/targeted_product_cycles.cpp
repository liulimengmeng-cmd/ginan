// Pending execution.  Independent literal/rank and Gaussian Schur-complement
// oracles; do not use the production planner/conditioner to compute the oracle.
#include "common/zhangTargetedProductCycles.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
int checks=0;
void require(bool condition,const char* message)
{
    ++checks;
    if(!condition) throw std::runtime_error(message);
}
bool close(double a,long double b,double tolerance=1e-9)
{
    return std::isfinite(a) && std::abs(static_cast<long double>(a)-b)<=
        tolerance*std::max(1.0L,std::abs(b));
}
ZhangExactVector unit(int dimension,int column,int coefficient=1)
{
    ZhangExactVector row(dimension);row[column]=coefficient;return row;
}
const ZhangTargetedCycleBlock* single(const ZhangTargetedCyclePlan& plan,int source)
{
    for(const auto& block:plan.blocks)
        if(block.candidateIndices==std::vector<int>{source}) return &block;
    return nullptr;
}

// Literal long-double scalar Schur update.  No Eigen inverse, production
// marginal workspace, exact row reducer or product projector is used here.
using Literal=std::vector<std::vector<long double>>;
Literal scalarCondition(Literal covariance,int column)
{
    const auto before=covariance;
    const long double variance=before[column][column];
    if(variance<=0) throw std::runtime_error("invalid literal oracle variance");
    for(std::size_t i=0;i<before.size();++i)
    for(std::size_t j=0;j<before.size();++j)
        covariance[i][j]-=before[i][column]*before[column][j]/variance;
    return covariance;
}
Eigen::MatrixXd numeric(const Literal& literal)
{
    Eigen::MatrixXd out(literal.size(),literal.size());
    for(std::size_t r=0;r<literal.size();++r)
    for(std::size_t c=0;c<literal.size();++c) out(r,c)=static_cast<double>(literal[r][c]);
    return out;
}
} // namespace

int main()
{
    try
    {
        const Literal literal={
            {4,.3L,1,.1L,.4L},
            {.3L,3,.8L,.05L,.2L},
            {1,.8L,3,.1L,.3L},
            {.1L,.05L,.1L,1,.05L},
            {.4L,.2L,.3L,.05L,2}
        };
        const Eigen::MatrixXd covariance=numeric(literal);
        Eigen::VectorXd mean(5);mean<<.2,-.1,.5,.4,1.8;
        const ZhangExactMatrix targets={unit(5,0),unit(5,1)};
        const ZhangExactMatrix held={unit(5,4)};
        const ZhangExactVector values={2};
        const ZhangExactMatrix candidates={unit(5,2),unit(5,3)};
        const auto savedMean=mean;
        const auto savedQ=covariance;
        const auto savedTargets=targets,savedHeld=held,savedCandidates=candidates;
        const auto plan=zhangPlanTargetedProductCycles(mean,covariance,targets,held,values,
            candidates,8,24);
        require(plan.valid && plan.targetRank==2 && plan.heldRank==1,
            "literal independent dual target rank differs");
        require(plan.independentCandidates==2 && plan.shortlistedCandidates==2,
            "eligible auxiliary cycle catalogue was truncated");
        const auto first=single(plan,0),second=single(plan,1);
        require(first && second,"cross-correlated independent singleton was not proposed");
        const auto heldOracle=scalarCondition(literal,4);
        const auto firstOracle=scalarCondition(heldOracle,2);
        const auto secondOracle=scalarCondition(heldOracle,3);
        const long double before=heldOracle[0][0]+heldOracle[1][1];
        const long double afterFirst=firstOracle[0][0]+firstOracle[1][1];
        const long double afterSecond=secondOracle[0][0]+secondOracle[1][1];
        require(close(first->conditionalTargetTraceBefore,before) &&
            close(first->conditionalTargetTraceAfter,afterFirst) &&
            close(first->varianceGain,before-afterFirst),
            "full held/candidate/target covariance differs from independent scalar oracle");
        require(close(second->varianceGain,before-afterSecond) &&
            first->varianceGain>second->varianceGain,
            "cross covariance did not rank the useful auxiliary cycle first");
        require(first->targetRankReduction==0 && first->targetRankAfter==2 &&
            first->exactRankGain==1,
            "variance improvement was misreported as an exact product consequence");
        const ZhangTargetedCycleBlock* pair=nullptr;
        for(const auto& block:plan.blocks) if(block.rows.size()==2) pair=&block;
        require(pair!=nullptr,"bounded independent auxiliary pair was not proposed");
        const auto pairOracle=scalarCondition(firstOracle,3);
        require(close(pair->conditionalTargetTraceAfter,pairOracle[0][0]+pairOracle[1][1]),
            "joint block score omitted covariance between auxiliary cycles");
        require((mean.array()==savedMean.array()).all() &&
            (covariance.array()==savedQ.array()).all() && targets==savedTargets &&
            held==savedHeld && candidates==savedCandidates,
            "planning mutated the protected posterior or root-coordinate catalogue");

        // A diagonal surrogate erases these correlations.  It must not be
        // reported as the same useful proposal set as the full covariance.
        const Eigen::MatrixXd diagonal=covariance.diagonal().asDiagonal();
        const auto noCross=zhangPlanTargetedProductCycles(mean,diagonal,targets,held,values,
            candidates,8,24);
        require(noCross.valid && noCross.blocks.empty(),
            "unrelated cycles were proposed without rank or conditional covariance gain");

        // Exact target directions outrank merely correlated directions; the
        // input integer coefficients and the non-saturated image survive.
        const Eigen::VectorXd zero=Eigen::VectorXd::Zero(4);
        const Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(4,4);
        const ZhangExactMatrix evenTargets={unit(4,0,2),unit(4,1,4)};
        const ZhangExactMatrix nonprimitive={unit(4,0,2),unit(4,1,4),unit(4,2)};
        const auto even=zhangPlanTargetedProductCycles(zero,identity,evenTargets,{}, {},
            nonprimitive,8,24);
        require(even.valid && even.targetRank==2 && !even.blocks.empty() &&
            even.blocks.front().targetRankAfter==0,
            "two independent integer product directions were not prioritised");
        const auto twice=single(even,0),fourTimes=single(even,1);
        require(twice && fourTimes && twice->rows.front()==unit(4,0,2) &&
            fourTimes->rows.front()==unit(4,1,4),
            "non-saturated 2Z/4Z proposal was divided by a gcd");
        for(int n0=-3;n0<=3;++n0) for(int n1=-3;n1<=3;++n1)
        {
            const int image0=2*n0,image1=4*n1;
            require(image0%2==0 && image1%4==0,
                "independent literal integer image lost divisibility");
        }
        const auto capped=zhangPlanTargetedProductCycles(zero,identity,evenTargets,{}, {},
            nonprimitive,1,1);
        require(capped.valid && capped.shortlistedCandidates==1 &&
            capped.blocks.size()==1 && capped.blocks.front().rows.size()==1,
            "proposal shortlist or block limit was ignored");

        // Signed duplicates may share a proposal; a scaled row is NEVER
        // primitive-normalised or confused with the original exact row.
        const ZhangExactMatrix duplicates={unit(5,2),unit(5,2,-1),unit(5,2,2),
            unit(5,3),unit(5,4),ZhangExactVector(5)};
        const auto dedup=zhangPlanTargetedProductCycles(mean,covariance,targets,held,values,
            duplicates,8,24);
        require(dedup.valid && dedup.skippedDuplicate==1 &&
            dedup.skippedZeroOrDependent==2,
            "signed duplicate, held row, or zero row accounting differs");
        const auto scaled=single(dedup,2);
        require(scaled && scaled->rows.front()==unit(5,2,2) &&
            close(scaled->varianceGain,first->varianceGain),
            "scaled exact candidate changed its preserved row or Gaussian information");
        for(const auto& block:dedup.blocks)
        {
            bool hasUnit=false,hasTwice=false;
            for(int source:block.candidateIndices)
            {hasUnit|=source==0;hasTwice|=source==2;}
            require(!(hasUnit && hasTwice),
                "joint proposal contains an exact-dependent integer row");
        }

        // Reuse a matching immutable affine owner; mismatched RHS must fail.
        const auto domain=zhangR47CompileProductSearchFrame(targets,held,values,5);
        const auto oldOwner=domain.affine;
        const auto oldBasis=domain.affine->deterministicBasis;
        ZhangR48MarginalWorkspace workspace(mean,covariance);
        const auto reused=zhangPlanTargetedProductCycles(mean,covariance,targets,held,values,
            candidates,8,24,&workspace,&domain);
        require(reused.valid && single(reused,0) &&
            close(single(reused,0)->varianceGain,first->varianceGain),
            "matching shared domain/workspace changed proposal information");
        require(domain.affine==oldOwner && domain.affine->deterministicBasis==oldBasis,
            "planning changed an immutable affine owner");
        require(!zhangPlanTargetedProductCycles(mean,covariance,targets,held,{3},
            candidates,8,24,nullptr,&domain).valid,
            "shared affine owner was reused for another RHS");
        require(!zhangPlanTargetedProductCycles(mean,diagonal,targets,held,values,
            candidates,8,24,&workspace).valid,
            "stale workspace with unchanged variances but changed cross blocks was accepted");
        Eigen::MatrixXd changedCandidateCross=covariance;
        changedCandidateCross(2,3)=changedCandidateCross(3,2)=.2;
        require(!zhangPlanTargetedProductCycles(mean,changedCandidateCross,targets,held,values,
            candidates,8,24,&workspace).valid,
            "shortlist candidate-candidate cross blocks were not checked against the root");

        // The entire pool is intentionally larger than the shortlist.  Only
        // scalar T0,T1,c screening may precede the <=24-row joint marginal.
        // All ten held equations remain active; they are NOT shortlisted.
        const int largeDimension=40;
        const auto largeMean=Eigen::VectorXd::Zero(largeDimension).eval();
        Eigen::MatrixXd largeQ=Eigen::MatrixXd::Identity(largeDimension,largeDimension);
        ZhangExactMatrix largePool,largeHeld;
        for(int c=2;c<largeDimension;++c)
        {
            largeQ(0,c)=largeQ(c,0)=.01;
            largeQ(1,c)=largeQ(c,1)=.005;
            largePool.push_back(unit(largeDimension,c));
            if(c>=30) largeHeld.push_back(unit(largeDimension,c));
        }
        const auto bounded=zhangPlanTargetedProductCycles(largeMean,largeQ,
            {unit(largeDimension,0),unit(largeDimension,1)},largeHeld,
            ZhangExactVector(largeHeld.size()),largePool,8,24);
        require(bounded.valid && bounded.independentCandidates==28 &&
            bounded.heldRank==10 && bounded.shortlistedCandidates==24 &&
            bounded.jointCandidateRows==24 && bounded.maximumNumericalTargetRows==26,
            "large candidate pool was projected jointly before a bounded shortlist");
        require(!zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{},{},{unit(4,2)},8,25).valid,
            "hard 24-candidate joint-projection memory bound was bypassed");

        // A root column absent from an old compact frame is FREE, not zero.
        Eigen::VectorXd mean3=Eigen::VectorXd::Zero(3);
        Eigen::MatrixXd q3=Eigen::MatrixXd::Identity(3,3);q3(1,2)=q3(2,1)=.6;
        const ZhangExactMatrix target3={unit(3,0),unit(3,1)};
        const ZhangExactMatrix held3={unit(3,0)};
        const auto compact=zhangR47CompileProductSearchFrame(target3,held3,{0},3);
        require(compact.valid && compact.columns==std::vector<int>({0,1}),
            "literal compact domain includes the unmentioned free column");
        const auto freeColumn=zhangPlanTargetedProductCycles(mean3,q3,target3,held3,{0},
            {unit(3,2)},8,24,nullptr,&compact);
        require(freeColumn.valid && freeColumn.targetRank==1 && single(freeColumn,0) &&
            single(freeColumn,0)->targetRankReduction==0,
            "unknown compact-domain column was treated as a deterministic zero");

        // Logical consequence is not integer row-lattice membership: 2n=2
        // determines n=1, whereas 2n=1 is INTEGER infeasible, despite having a
        // rational solution.  Neither case may create a new ILS row for n.
        const auto implied=zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{unit(4,0,2)},{2},
            {unit(4,0),unit(4,1)},8,24);
        require(implied.valid && implied.targetRank==1 && !single(implied,0),
            "rational consequence was mistaken for a new independent target direction");
        require(!zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{unit(4,0,2)},{1},
            {unit(4,1)},8,24).valid,
            "rational but non-integer held system was accepted");
        require(!zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{unit(4,0),unit(4,0)},{1,2},
            {unit(4,1)},8,24).valid,
            "inconsistent held RHS was accepted");
        const auto determined=zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{unit(4,0),unit(4,1)},{2,-1},
            {unit(4,2)},8,24);
        require(determined.valid && determined.targetRank==0 && determined.blocks.empty() &&
            determined.status=="TARGET_ALREADY_DETERMINED",
            "already determined target manufactured an auxiliary AR decision");
        require(!zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{},{},{{1,0}},8,24).valid,
            "malformed/unmapped root row was silently padded with zero");
        Eigen::MatrixXd badQ=identity;badQ(3,3)=-1;
        require(!zhangPlanTargetedProductCycles(zero,badQ,
            {unit(4,0),unit(4,1)},{},{},{unit(4,2)},8,24).valid,
            "non-PSD root covariance was accepted");
        Eigen::MatrixXd asymmetric=identity;asymmetric(0,2)=.1;
        require(!zhangPlanTargetedProductCycles(zero,asymmetric,
            {unit(4,0),unit(4,1)},{},{},{unit(4,2)},8,24).valid,
            "asymmetric full covariance was accepted");
        require(!zhangPlanTargetedProductCycles(zero,identity,
            {unit(4,0),unit(4,1)},{},{},{unit(4,2)},0,24).valid,
            "invalid proposal limit was accepted");
        std::cout<<"targeted_product_cycles: "<<checks<<" checks passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"targeted_product_cycles failed: "<<error.what()<<"\n";
        return 1;
    }
}
