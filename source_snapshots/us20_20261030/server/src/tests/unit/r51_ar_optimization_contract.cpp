#include "common/zhangR49PosteriorWorkspace.hpp"
#include "common/zhangR51QuotientCompatibility.hpp"
#include "common/zhangPrivateUnitMembership.hpp"
#include <iostream>
#include <random>
#include <stdexcept>

static void check(bool pass,const char* why)
{if(!pass)throw std::runtime_error(why);}

static void privateUnits()
{
    std::mt19937 rng(20261019);
    for(int trial=0;trial<180;++trial)
    {
        const int m=1+rng()%9,n=m+4;
        ZhangExactMatrix rows(m,ZhangExactVector(n));
        ZhangExactVector witness(m),target(n);
        for(int i=0;i<m;++i)
        {
            rows[i][i]=i%2?-1:1;witness[i]=int(rng()%9)-4;
            for(int j=m;j<n;++j)rows[i][j]=int(rng()%9)-4;
            for(int j=0;j<n;++j)target[j]+=witness[i]*rows[i][j];
        }
        auto bad=target;bad[n-1]+=1;
        auto got=zhangPrivateUnitLatticeMembershipBatch(rows,{target,bad});
        // Independent oracle: the private columns fix the sole possible
        // coefficients; changed non-pivot data cannot be represented.
        check(got[0].contained && got[0].combination==witness && !got[1].contained,
            "private-column exact witness oracle failed");
        auto old=zhangIntegerRowLatticeContainsBatch(rows,{target,bad});
        for(int i=0;i<2;++i)check(got[i].contained==old[i].contained && got[i].rank==old[i].rank &&
            got[i].smithInvariants==old[i].smithInvariants && got[i].combination==old[i].combination,
            "private-column result differs from generic exact backend");
    }
    auto parity=zhangPrivateUnitLatticeMembershipBatch({{2,0},{0,3}},{{1,0},{4,6}});
    check(!parity[0].contained && parity[1].contained && parity[1].combination==ZhangExactVector({2,2}),
        "non-unit fallback lost divisibility");
    auto dependent=zhangPrivateUnitLatticeMembershipBatch({{1,1},{2,2}},{{3,3},{3,4}});
    check(dependent[0].contained && !dependent[1].contained,"dependent fallback failed");
    ZhangExactInteger big=1;big<<=150;
    auto large=zhangPrivateUnitLatticeMembershipBatch({{-1,big}},{{2,-2*big}});
    check(large[0].contained && large[0].combination==ZhangExactVector({-2}),"large integer witness failed");
    for(int maximum:{1,3}) {
        ZhangExactWorkBudget budget;budget.maximumCheckpoints=maximum;bool cancelled=false;
        try {ZhangScopedExactWorkBudget scope(budget);zhangPrivateUnitLatticeMembershipBatch({{1,2}},{{3,6}});}
        catch(const ZhangExactWorkBudgetExceeded&){cancelled=true;}
        check(cancelled && !zhangActiveExactWorkBudget,"exact cancellation disappeared in scan/witness path");
    }
    ZhangExactMatrix rows(80,ZhangExactVector(83)),targets(16,ZhangExactVector(83));
    for(int i=0;i<80;++i)rows[i][i]=1;
    for(int i=0;i<16;++i)targets[i][i]=i;
    auto parallel=zhangPrivateUnitLatticeMembershipBatch(rows,targets);
    check(std::all_of(parallel.begin(),parallel.end(),[](const auto& r){return r.contained;}),
        "parallel exact witnesses failed");
    check(zhangPrivateUnitMembershipMetrics.maximumWorkers==omp_get_max_threads(),
        "actual all-CPU exact target team was not created");
    std::cout<<"PASS private unit witnesses, 180 generic comparisons, parity, dependent rows, 150-bit values, cancellation\n";
}

static void components()
{
    ZhangR51QuotientCompatibility domain(5,true),foreign(5,true);
    auto accept=[&](ZhangExactMatrix rows,ZhangExactVector values,std::vector<int> columns){
        auto t=domain.assess(rows,values,columns);check(t.compatible() && domain.commit(t),"valid component rejected");};
    accept({{1,-1}},{0},{0,1});
    accept({{2}},{2},{3});
    check(domain.maximumCheckedColumns<=2 && domain.disjointChecks==2,"disjoint component grew full domain");
    auto stale=domain.assess({{1}},{3},{2});
    accept({{1,-1}},{0},{1,2}); // Must walk 2 -> row -> 1 -> previous row -> 0.
    check(!domain.commit(stale) && !foreign.commit(stale),"stale or foreign component ticket accepted");
    check(!domain.assess({{1,-1}},{1},{0,2}).compatible(),"transitive conflict was missed");
    check(!domain.assess({{0}},{1},{4}).compatible(),"zero-row contradiction accepted");
    check(domain.assess({{0}},{0},{4}).compatible(),"zero identity rejected");
    check(!domain.assess({{2}},{1},{4}).compatible(),"odd parity accepted");
    check(domain.assess({{2}},{4},{4}).compatible(),"even parity rejected");
    // All variables fixed: independently compute every expected consequence.
    ZhangR51QuotientCompatibility fixed(3,true);
    auto all=fixed.assess({{1,0,0},{0,1,0},{0,0,1}},{2,-3,4},{0,1,2});
    check(all.compatible() && fixed.commit(all),"fixed oracle setup failed");
    for(int a=-2;a<=2;++a)for(int b=-2;b<=2;++b)for(int c=-2;c<=2;++c)
        for(int error:{0,1})
            check(fixed.assess({{a,b,c}},{2*a-3*b+4*c+error},{0,1,2}).compatible()==(error==0),
                "independent affine consequence oracle failed");
    std::cout<<"PASS connected-component feasibility, transitive closure, zero rows, parity, 250 independent consequences, tickets\n";
}

static void roots()
{
    const int n=96;
    MatrixXd l(n,n);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)l(i,j)=0.04*std::sin(i*0.17+j*0.31);
    MatrixXd p=l*l.transpose()+MatrixXd::Identity(n,n);
    const int blasBefore=BlasThreading::openblas_get_num_threads();
    const auto eigen=ZhangArThreading::eigenSelfAdjoint(p);
    check(eigen.info()==Eigen::Success && BlasThreading::openblas_get_num_threads()==blasBefore &&
        (eigen.eigenvectors()*eigen.eigenvalues().asDiagonal()*eigen.eigenvectors().transpose()-p).norm()/p.norm()<1e-12,
        "Eigen tridiagonal phase changed covariance or leaked its serial BLAS setting");
    std::vector<std::string> order;
    for(int i=0;i<n;++i)order.push_back(std::to_string(i));
    ZhangR49PosteriorWorkspace workspace;
    {
        ZhangR49FloatRootScope pin(p,order,"epoch");
        check(workspace.factor(p,order,"epoch"),"root failed");
        const auto initial=workspace.decompositions;
        auto root=workspace.squareRoot;root.col(0)*=.5;
        check(workspace.formGram(std::move(root),order,"epoch"),"Gram failed");
        check(workspace.factor(p,order,"epoch") && workspace.decompositions==initial &&
            workspace.protectedFloatHits==1,"FLOAT root was decomposed again");
        MatrixXd a=MatrixXd::Zero(3,n);a(0,0)=1;a(0,1)=-1;a(1,4)=2;a(2,5)=1;
        auto root2=workspace.squareRoot;
        MatrixXd w=a*root2;
        Eigen::JacobiSVD<MatrixXd> svd(w,Eigen::ComputeThinV);
        MatrixXd v=svd.matrixV();
        MatrixXd projected=root2-(root2*v)*v.transpose();
        MatrixXd actual=projected*projected.transpose();
        MatrixXd ap=a*p;
        // Independent subtractive covariance oracle, small well-conditioned case.
        MatrixXd expected=p-ap.transpose()*(ap*a.transpose()).ldlt().solve(ap);
        check((actual-expected).cwiseAbs().maxCoeff()<2e-12,"full conditioned covariance mismatch");
        MatrixXd changed=p;changed(0,0)+=.2;
        check(workspace.factor(changed,order,"epoch") && workspace.decompositions==initial+1,
            "changed covariance reused stale root");
        auto reordered=order;std::swap(reordered[0],reordered[1]);
        check(workspace.factor(p,reordered,"epoch") && workspace.decompositions==initial+2,
            "changed state order reused root");
        check(workspace.factor(p,order,"next") && workspace.decompositions==initial+3,
            "changed epoch reused root");
    }
    check(!zhangR49PinnedFloatRoot,"epoch root retained beyond AR scope");
    MatrixXd bad=MatrixXd::Identity(n,n);bad(0,0)=-1;
    check(!workspace.factor(bad,order,"bad") && workspace.reason=="PRIOR_COVARIANCE_NOT_PSD",
        "negative covariance accepted");
    MatrixXd singular=p;singular.row(0).setZero();singular.col(0).setZero();
    check(workspace.factor(singular,order,"singular") &&
        (workspace.squareRoot*workspace.squareRoot.transpose()-singular).norm()/singular.norm()<1e-12,
        "PSD covariance reconstruction failed");
    std::cout<<"PASS production root: protected reuse, full covariance oracle, P/order/epoch invalidation, PSD, negative rejection, lifetime\n";
}

int main()
{
    omp_set_num_threads(4);BlasThreading::openblas_set_num_threads(4);
    const int oldEigen=Eigen::nbThreads(),oldOmp=omp_get_max_threads();
    const int oldBlas=BlasThreading::openblas_get_num_threads();
    {
        ZhangArThreading::ScopedThreads threads(true);
        check(omp_get_max_threads()==ZhangFloatThreading::availableThreads(),"AR all-CPU profile missing");
        const char* previous=std::getenv("ZHANG_R51_AR_THREADS");
        const bool hadOverride=previous!=nullptr;
        const std::string saved=previous?previous:"";
        setenv("ZHANG_R51_AR_THREADS","legacy",1);
        {
            const int inherited=BlasThreading::openblas_get_num_threads();
            ZhangArThreading::ScopedThreads legacy(true);
            check(!ZhangArThreading::active && inherited==BlasThreading::openblas_get_num_threads(),
                "legacy policy unexpectedly enabled the managed Eigen thread scope");
        }
        if(hadOverride)setenv("ZHANG_R51_AR_THREADS",saved.c_str(),1);
        else unsetenv("ZHANG_R51_AR_THREADS");
        check(ZhangArThreading::active,"nested legacy scope did not restore managed policy");
        privateUnits();components();roots();
    }
    check(omp_get_max_threads()==oldOmp && Eigen::nbThreads()==oldEigen &&
        BlasThreading::openblas_get_num_threads()==oldBlas && !ZhangArThreading::active,
        "AR profile leaked into FLOAT");
    std::cout<<"PASS AR thread restoration\n";
}
