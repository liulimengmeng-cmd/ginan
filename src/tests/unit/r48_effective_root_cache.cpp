#include "common/zhangR48Marginal.hpp"
#include "common/zhangR49PosteriorWorkspace.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
int assertions=0;
void require(bool value,const char* reason)
{
    ++assertions;
    if(!value)throw std::runtime_error(reason);
}
bool exact(const MatrixXd& a,const MatrixXd& b)
{
    return a.rows()==b.rows() && a.cols()==b.cols() &&
        (!a.size() || std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);
}
bool near(const MatrixXd& a,const MatrixXd& b,double tolerance=1e-11)
{
    return a.rows()==b.rows() && a.cols()==b.cols() &&
        (a-b).norm()<=tolerance*std::max(1.0,b.norm());
}
const std::vector<std::string> order={"G01-L1","E11-L1","C36-L1"};

void sameRootAndMeans()
{
    MatrixXd p(3,3);p<<4,.5,0,.5,3,0,0,0,2;
    VectorXd mean(3);mean<<1,2,3;
    MatrixXd j=MatrixXd::Identity(3,3),h(1,3);h<<1,0,0;
    VectorXd rhs(1);rhs<<1.25;
    ZhangPosteriorEffectiveWorkspaceCache cache;
    ZhangR48MarginalWorkspace first(mean,p,&cache,order,"t0");
    ZhangR48MarginalWorkspace second(mean,p,&cache,order,"t0");
    ZhangR48MarginalWorkspace legacy(mean,p);
    require(first.valid && second.valid && legacy.valid,"valid SPD root rejected");
    require(cache.misses==1 && cache.builds==1 && cache.hits==1,"same root not reused");
    require(!first.rootIdentity().cacheHit && second.rootIdentity().cacheHit,"wrong hit flag");
    require(first.rootIdentity().reason=="EMPTY_CACHE" &&
        second.rootIdentity().reason=="EXACT_ROOT_HIT","wrong first/hit reason");
    require(first.rootIdentity().ownerId==cache.id() && cache.id()!=0,"owner identity missing");
    require(first.rootIdentity().generation==second.rootIdentity().generation,"hit changed generation");
    require(first.squareRoot.data()==second.squareRoot.data(),"cached root was copied");
    require(exact(first.squareRoot,legacy.squareRoot),"effective eig path/direction/order changed");
    auto got=second.project(j,h,rhs),baseline=legacy.project(j,h,rhs);
    require(got.valid && baseline.valid && exact(got.mean,baseline.mean) &&
        exact(got.covariance,baseline.covariance),"cached/uncached moments not byte exact");
    // Independent closed-form Gaussian conditioning (no cached eig oracle).
    const VectorXd expectedMean=mean+p.col(0)*((rhs(0)-mean(0))/p(0,0));
    const MatrixXd expectedP=p-p.col(0)*p.row(0)/p(0,0);
    require(near(got.mean,expectedMean) && near(got.covariance,expectedP),"analytic conditional oracle failed");
    VectorXd otherMean=mean;otherMean(1)+=7;otherMean(2)-=2;
    ZhangR48MarginalWorkspace changed(otherMean,p,&cache,order,"t0");
    auto changedResult=changed.project(j,h,rhs);
    const VectorXd changedExpected=otherMean+p.col(0)*((rhs(0)-otherMean(0))/p(0,0));
    require(changed.rootIdentity().cacheHit && cache.builds==1,"mean change rebuilt covariance factor");
    require(changedResult.valid && near(changedResult.mean,changedExpected) &&
        (changedResult.mean-got.mean).norm()>1,"mean/target cache leaked between owners");
    require(changed.targetHits==0 && changed.decompositions==1 && second.targetHits==0,
        "conditioner/target caches were shared");
    VectorXd invalidMean=mean;invalidMean(0)=std::numeric_limits<double>::quiet_NaN();
    ZhangR48MarginalWorkspace invalid(invalidMean,p,&cache,order,"t0");
    require(!invalid.valid && invalid.rootIdentity().cacheHit,"invalid mean accepted by factor reuse");
    require(cache.retainedBytes==2*9*sizeof(double) &&
        cache.lookupSeconds>=0 && cache.buildSeconds>=0,"cache accounting invalid");
}

void changedIdentityAndLifetime()
{
    MatrixXd p=MatrixXd::Identity(3,3);VectorXd mean=VectorXd::Zero(3);
    ZhangPosteriorEffectiveWorkspaceCache cache;
    ZhangR48MarginalWorkspace first(mean,p,&cache,order,"t0");
    const MatrixXd oldRoot=first.squareRoot;
    auto reordered=order;std::swap(reordered[0],reordered[1]);
    ZhangR48MarginalWorkspace otherOrder(mean,p,&cache,reordered,"t0");
    require(!otherOrder.rootIdentity().cacheHit &&
        otherOrder.rootIdentity().reason=="STATE_ORDER_CHANGED","same-size order incorrectly reused");
    require(cache.builds==2 && cache.invalidations==1,"order invalidation missing");
    ZhangR48MarginalWorkspace otherEpoch(mean,p,&cache,reordered,"t1");
    require(otherEpoch.rootIdentity().reason=="EPOCH_CHANGED" && cache.builds==3,"epoch incorrectly reused");
    p(0,0)=2;
    ZhangR48MarginalWorkspace otherP(mean,p,&cache,reordered,"t1");
    require(otherP.rootIdentity().reason=="COVARIANCE_CHANGED" && cache.builds==4,
        "same-size changed covariance incorrectly reused");
    require(exact(first.squareRoot,oldRoot),"eviction mutated live immutable lease");
    p(0,1)=-0.0;
    ZhangR48MarginalWorkspace signedZero(mean,p,&cache,reordered,"t1");
    require(signedZero.rootIdentity().reason=="COVARIANCE_CHANGED" && cache.builds==5,
        "P comparison silently ignored signed-zero bytes");
    const auto previous=cache.generation;
    cache.clear();
    require(cache.retainedBytes==0 && cache.generation>previous,"clear did not invalidate cache");
    require(exact(first.squareRoot,oldRoot),"clear invalidated a live borrower");
    ZhangR48MarginalWorkspace afterClear(mean,p,&cache,reordered,"t1");
    require(afterClear.rootIdentity().reason=="EMPTY_CACHE" && cache.builds==6,"cleared cache reused old root");
    ZhangPosteriorEffectiveWorkspaceCache otherOwner;
    ZhangR48MarginalWorkspace isolated(mean,p,&otherOwner,reordered,"t1");
    require(otherOwner.builds==1 && !isolated.rootIdentity().cacheHit &&
        isolated.rootIdentity().ownerId!=first.rootIdentity().ownerId,"independent owner reused a root");
}

void deterministicAndNonPsd()
{
    MatrixXd p=MatrixXd::Zero(3,3);p(0,0)=4;p(1,1)=2;
    VectorXd mean=VectorXd::Zero(3);MatrixXd h(1,3);h<<0,0,1;
    ZhangPosteriorEffectiveWorkspaceCache cache;
    ZhangR48MarginalWorkspace first(mean,p,&cache,order,"t0");
    ZhangR48MarginalWorkspace cached(mean,p,&cache,order,"t0");
    require(first.valid && cached.valid && cached.squareRoot.cols()==2,"deterministic rank changed");
    auto consistent=cached.project(MatrixXd::Identity(3,3),h,VectorXd::Zero(1));
    auto contradiction=cached.project(MatrixXd::Identity(3,3),h,VectorXd::Ones(1));
    require(consistent.valid && near(consistent.mean,mean) && near(consistent.covariance,p),
        "deterministic consistent condition rejected");
    require(!contradiction.valid && contradiction.failureReason=="DETERMINISTIC_OR_NON_PSD_CONDITIONER",
        "deterministic contradiction lost on cache hit");
    // The original effective PSD threshold remains authoritative, including
    // tiny negative numerical-null modes and strictly non-PSD rejection.
    p(2,2)=-1e-16;
    ZhangR48MarginalWorkspace numericalNull(mean,p,&cache,order,"t0");
    ZhangR48MarginalWorkspace legacyNull(mean,p);
    require(numericalNull.valid && exact(numericalNull.squareRoot,legacyNull.squareRoot) &&
        numericalNull.squareRoot.cols()==2,"original effective tolerance changed");
    p(2,2)=-1e-6;
    ZhangR48MarginalWorkspace rejected(mean,p,&cache,order,"t0");
    ZhangR48MarginalWorkspace legacyRejected(mean,p);
    require(!rejected.valid && !legacyRejected.valid && cache.retainedBytes==0,
        "non-PSD root accepted/cached");
    const auto builds=cache.builds,hits=cache.hits;
    ZhangR48MarginalWorkspace rejectedAgain(mean,p,&cache,order,"t0");
    require(!rejectedAgain.valid && cache.builds==builds+1 && cache.hits==hits,
        "failed PSD audit reused as a valid root");
    require(rejectedAgain.rootIdentity().reason==
        "BUILD_REJECTED:POSTERIOR_EFFECTIVE_BASIS_COVARIANCE_NOT_PSD","PSD failure reason changed");
    p(0,0)=std::numeric_limits<double>::infinity();
    ZhangR48MarginalWorkspace nonFinite(mean,p,&cache,order,"t0");
    require(!nonFinite.valid && cache.retainedBytes==0,"non-finite covariance cached");
}

void capacityAndMissingIdentity()
{
    const MatrixXd p=MatrixXd::Identity(3,3);const VectorXd mean=VectorXd::Zero(3);
    ZhangPosteriorEffectiveWorkspaceCache tooSmall(2*9*sizeof(double)-1);
    ZhangR48MarginalWorkspace first(mean,p,&tooSmall,order,"t0");
    ZhangR48MarginalWorkspace second(mean,p,&tooSmall,order,"t0");
    require(first.valid && second.valid && tooSmall.retainedBytes==0 &&
        tooSmall.capacityBypasses==2 && tooSmall.builds==2 && tooSmall.hits==0,
        "byte cap changed science or retained oversize root");
    require(first.rootIdentity().reason=="CACHE_BYTE_BUDGET_EXCEEDED" &&
        first.squareRoot.rows()==3 && first.squareRoot.cols()==3,"oversize fallback truncated root");
    ZhangPosteriorEffectiveWorkspaceCache exactBudget(2*9*sizeof(double));
    ZhangR48MarginalWorkspace boundary(mean,p,&exactBudget,order,"t0");
    ZhangR48MarginalWorkspace boundaryHit(mean,p,&exactBudget,order,"t0");
    require(boundary.valid && boundaryHit.rootIdentity().cacheHit &&
        exactBudget.retainedBytes==exactBudget.byteBudget(),"inclusive byte budget boundary wrong");
    ZhangR48MarginalWorkspace missingOrder(mean,p,&exactBudget,{},"t0");
    require(missingOrder.valid && !missingOrder.rootIdentity().cacheHit &&
        missingOrder.rootIdentity().reason=="INCOMPLETE_ROOT_IDENTITY" &&
        exactBudget.retainedBytes==0 && exactBudget.identityBypasses==1,
        "missing identity did not fall back safely");
    ZhangR48MarginalWorkspace missingEpoch(mean,p,&exactBudget,order,"");
    require(missingEpoch.valid && exactBudget.identityBypasses==2 && exactBudget.hits==1,
        "empty epoch reused cached root");
}

void nestedScopeRestoration()
{
    require(zhangR49PinnedFloatRoot==nullptr,"test inherited an unexpected AR owner");
    const MatrixXd p=MatrixXd::Identity(3,3);const VectorXd mean=VectorXd::Zero(3);
    std::shared_ptr<const ZhangPosteriorEffectiveWorkspace> borrowed;
    {
        ZhangR49FloatRootScope outer(p,order,"t0");
        auto* outerOwner=zhangR49PinnedFloatRoot;
        ZhangR48MarginalWorkspace first(mean,p,&outerOwner->effectiveRoots,order,"t0");
        borrowed=first.rootIdentity().workspace;
        const auto id=outerOwner->effectiveRoots.id();
        {
            ZhangR49FloatRootScope inner(p,order,"t0");
            auto* innerOwner=zhangR49PinnedFloatRoot;
            ZhangR48MarginalWorkspace innerFirst(mean,p,&innerOwner->effectiveRoots,order,"t0");
            require(innerOwner!=outerOwner && innerOwner->effectiveRoots.id()!=id &&
                !innerFirst.rootIdentity().cacheHit,"nested owner improperly shared a root");
            require(outerOwner->effectiveRoots.builds==1 && outerOwner->effectiveRoots.hits==0,
                "nested scope modified outer cache");
        }
        require(zhangR49PinnedFloatRoot==outerOwner,"nested scope did not restore outer owner");
        ZhangR48MarginalWorkspace after(mean,p,&outerOwner->effectiveRoots,order,"t0");
        require(after.rootIdentity().cacheHit && after.rootIdentity().ownerId==id,
            "restored outer scope lost original root");
    }
    require(zhangR49PinnedFloatRoot==nullptr && borrowed->valid && borrowed->squareRoot.rows()==3,
        "AR scope leaked owner or invalidated explicit borrower");
    ZhangR49FloatRootScope fresh(p,order,"t0");
    ZhangR48MarginalWorkspace first(mean,p,&zhangR49PinnedFloatRoot->effectiveRoots,order,"t0");
    require(!first.rootIdentity().cacheHit && zhangR49PinnedFloatRoot->effectiveRoots.builds==1,
        "later AR call reused a historical root");
}

void largeKernelContracts()
{
    // The real R49 full-root path switches to LAPACK at >=64. R48 must NOT
    // borrow that differently-thresholded root: it keeps the exact original
    // effective Eigen path, retained columns and ordering at this size too.
    constexpr int n=65;
    MatrixXd p=2*MatrixXd::Identity(n,n);
    const VectorXd common=VectorXd::LinSpaced(n,.01,.09);
    p.noalias()+=common*common.transpose();
    VectorXd mean=VectorXd::LinSpaced(n,-2,3);
    std::vector<std::string> keys(n);
    for(int i=0;i<n;++i)keys[i]="large-"+std::to_string(i);
    ZhangR49FloatRootScope scope(p,keys,"large");
    auto& cache=zhangR49PinnedFloatRoot->effectiveRoots;
    ZhangR49PosteriorWorkspace fullRoot;
    require(fullRoot.factor(p,keys,"large"),"large R49 full-root factor failed");
#ifdef GINAN_USE_OPENBLAS
    require(fullRoot.backend=="LAPACK_DSYEVD","large production LAPACK branch not exercised");
#else
    require(fullRoot.backend=="EIGEN_SELFADJOINT","non-BLAS reference backend unexpected");
#endif
    ZhangR48MarginalWorkspace first(mean,p,&cache,keys,"large");
    ZhangR48MarginalWorkspace cached(mean,p,&cache,keys,"large");
    ZhangR48MarginalWorkspace legacy(mean,p);
    require(first.valid && cached.valid && legacy.valid &&
        cached.rootIdentity().workspace->covarianceRank==n,"large SPD rank changed");
    require(cache.builds==1 && cache.hits==1 && cached.rootIdentity().cacheHit &&
        cached.squareRoot.data()==first.squareRoot.data() &&
        cached.squareRoot.data()!=fullRoot.squareRoot.data(),"large effective/full R49 roots aliased");
    require(exact(cached.squareRoot,legacy.squareRoot),"large effective root direction/order changed");
    require(cached.rootIdentity().workspace->eigenTolerance==legacy.rootIdentity().workspace->eigenTolerance &&
        cached.rootIdentity().workspace->minimumEigenvalue==legacy.rootIdentity().workspace->minimumEigenvalue,
        "large effective spectral policy changed");
    MatrixXd h=MatrixXd::Zero(3,n);h(0,0)=1;h(1,2)=1;h(2,n-1)=1;
    const MatrixXd j=MatrixXd::Identity(n,n);
    const VectorXd values=h*mean+VectorXd::Constant(3,.02);
    auto got=cached.project(j,h,values),baseline=legacy.project(j,h,values);
    require(got.valid && baseline.valid && exact(got.mean,baseline.mean) &&
        exact(got.covariance,baseline.covariance),"large cached/uncached full moments not byte exact");
    const MatrixXd s=h*p*h.transpose(),cross=p*h.transpose();
    const VectorXd expectedMean=mean+cross*s.ldlt().solve(values-h*mean);
    const MatrixXd expectedP=p-cross*s.ldlt().solve(cross.transpose());
    require(near(got.mean,expectedMean) && near(got.covariance,expectedP),
        "large independent LDLT conditional oracle failed");
    require(near(cached.squareRoot*cached.squareRoot.transpose(),p),"large full-rank covariance changed");

    p=MatrixXd::Identity(n,n);p(n-2,n-2)=1e-15;p(n-1,n-1)=0;
    // A positive direction below R48's ORIGINAL effective tolerance is kept
    // by R49's full root. This deliberately distinguishes the two consumers.
    require(fullRoot.factor(p,keys,"large"),"large semidefinite full root failed");
    ZhangR48MarginalWorkspace psdFirst(mean,p,&cache,keys,"large");
    ZhangR48MarginalWorkspace psdCached(mean,p,&cache,keys,"large");
    ZhangR48MarginalWorkspace psdLegacy(mean,p);
    require(psdCached.valid && psdLegacy.valid && psdCached.squareRoot.cols()==n-2 &&
        psdCached.rootIdentity().workspace->covarianceRank==n-2,"large effective PSD rank changed");
    require(cache.builds==2 && cache.hits==2 && psdFirst.rootIdentity().reason=="COVARIANCE_CHANGED" &&
        psdCached.rootIdentity().cacheHit,"large changed P did not invalidate/reuse correctly");
    require(exact(psdCached.squareRoot,psdLegacy.squareRoot),"large PSD eig ordering changed");
    require(fullRoot.squareRoot.row(n-2).squaredNorm()>0 && psdCached.squareRoot.row(n-2).squaredNorm()==0,
        "R49 threshold leaked into effective numerical-null space");
    auto psdGot=psdCached.project(j,MatrixXd(0,n),VectorXd(0));
    auto psdBaseline=psdLegacy.project(j,MatrixXd(0,n),VectorXd(0));
    MatrixXd expectedPsd=MatrixXd::Identity(n,n);expectedPsd(n-2,n-2)=0;expectedPsd(n-1,n-1)=0;
    require(psdGot.valid && psdBaseline.valid && exact(psdGot.mean,psdBaseline.mean) &&
        exact(psdGot.covariance,psdBaseline.covariance),"large PSD full moments not byte exact");
    require(near(psdGot.mean,mean) && near(psdGot.covariance,expectedPsd),
        "large analytic effective-null-space covariance oracle failed");
}
}

int main()
{
    try
    {
        ZhangArThreading::ScopedThreads policy(true);
        sameRootAndMeans();
        changedIdentityAndLifetime();
        deterministicAndNonPsd();
        capacityAndMissingIdentity();
        nestedScopeRestoration();
        largeKernelContracts();
        std::cout<<"R48_EFFECTIVE_ROOT_CACHE_PASS assertions="<<assertions
            <<" same_root_shared=1 different_mean_recomputed=1 exact_order_epoch_P=1"
            <<" original_eig_path=1 deterministic_psd_checked=1 bounded_cache=1 nested_scope=1\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"R48_EFFECTIVE_ROOT_CACHE_FAIL assertions="<<assertions<<" "<<error.what()<<"\n";
        return 1;
    }
}
