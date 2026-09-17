#include "common/zhangR51QuotientCompatibility.hpp"
#include "common/zhangR51PhysicalImage.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>

static void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static std::size_t comparisons=0,domains=0,rejections=0,discarded=0;
static double referenceSeconds=0,quotientSeconds=0;
static void domain(const ZhangExactMatrix& history,const ZhangExactVector& values,
                   const ZhangExactMatrix& targets,unsigned seed,int attempts=24)
{
    const int n=targets.front().size(),m=targets.size();
    auto frame=zhangR51PhysicalImage(targets,history,values,zhangExactIdentityMatrix(m),n);
    check(frame.valid,"reference domain must be valid");++domains;
    const int rank=frame.projector.size();if(!rank)return;
    ZhangR51QuotientCompatibility fast(rank,true);
    ZhangExactMatrix oldRows=history;ZhangExactVector oldValues=values;
    std::mt19937 random(seed);
    ZhangExactVector truth(rank);for(auto& x:truth)x=int(random()%11)-5;
    for(int attempt=0;attempt<attempts;++attempt){
        const int k=1+random()%std::min(rank,4),nr=1+random()%k;
        std::vector<int> subset;while(subset.size()<k){int c=random()%rank;if(std::find(subset.begin(),subset.end(),c)==subset.end())subset.push_back(c);}
        ZhangExactMatrix b(nr,ZhangExactVector(k));ZhangExactVector z(nr);
        for(int i=0;i<nr;++i){for(int c=0;c<k;++c){b[i][c]=int(random()%7)-3;z[i]+=b[i][c]*truth[subset[c]];}if(attempt%3==1)z[i]+=1+random()%3;}
        std::vector<ZhangR51RationalRow> p;ZhangR51RationalRow offset;
        for(int c:subset){p.push_back(frame.projector[c]);offset.push_back(frame.offsets[c]);}
        ZhangExactMatrix visible;ZhangExactVector rhs;
        check(zhangR51LiftRationalRows(b,z,p,offset,visible,rhs),"exact lift failed");
        const auto physical=zhangExactMultiply(visible,targets);
        auto all=oldRows;all.insert(all.end(),physical.begin(),physical.end());
        auto allValues=oldValues;allValues.insert(allValues.end(),rhs.begin(),rhs.end());
        auto start=std::chrono::steady_clock::now();
        bool expected=zhangR47AffineIntegerFeasible(all,allValues,n);
        referenceSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        start=std::chrono::steady_clock::now();auto ticket=fast.assess(b,z,subset);
        quotientSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        ++comparisons;check(ticket.compatible()==expected,"physical and quotient compatibility disagree");
        if(!expected){++rejections;check(!fast.commit(ticket),"incompatible ticket committed");}
        else if(attempt%4!=0){check(fast.commit(ticket),"accepted ticket failed");oldRows=std::move(all);oldValues=std::move(allValues);}
        else ++discarded; // Numerical/proof rejection: the next proposal must see the old domain.
    }
}
int main(){
    domain({{1,-2}},{1},{{1,0}},1); // Hidden retired arc leaves odd parity.
    domain({}, {},{{2,0},{1,3}},2); // Nonprimitive coupled image and rational lift.
    domain({{1,-2,0},{0,1,-3}},{1,-2},{{1,0,0},{0,1,0}},3);
    domain({{1,1,0}},{5},{{1,0,0},{0,1,0},{1,1,0},{0,0,2}},4); // Dependent target directions.
    ZhangExactInteger big=1;big<<=130;
    domain({{1,-big,0}},{big+1},{{1,0,0},{0,0,3}},5);
    domain({{2}},{6},{{1}},6); // Completely determined domain: no search.
    check(!zhangR51PhysicalImage({{1}},{{2}},{5},{{1}},1).valid,"infeasible initial domain accepted");
    std::mt19937 random(20260917);
    for(int test=0;test<100;++test){
        int n=2+random()%7,m=1+random()%std::min(n,5),h=random()%n;
        ZhangExactMatrix history(h,ZhangExactVector(n)),targets(m,ZhangExactVector(n));
        ZhangExactVector values(h),witness(n);for(auto& x:witness)x=int(random()%7)-3;
        for(int i=0;i<h;++i)for(int j=0;j<n;++j){history[i][j]=int(random()%5)-2;values[i]+=history[i][j]*witness[j];}
        for(auto& row:targets)for(auto& x:row)x=int(random()%5)-2;
        domain(history,values,targets,1000+test);
    }
    // Actual unconstrained-unit-image shortcut used at startup.
    auto unit=zhangR51UnconstrainedUnitImage({{{"a",1}},{{"b",-1}}},{{1,0},{0,1}});
    check(unit.valid && unit.projector.size()==2,"unit shortcut fixture invalid");
    ZhangR51QuotientCompatibility one(2,true),other(2,true),invalid(2,false),zero(0,true);
    auto stale=one.assess({{1}},{3},{0}),first=one.assess({{1}},{3},{0});
    check(first.compatible() && one.commit(first),"first commit failed");
    check(!one.commit(stale) && !one.commit(first) && !other.commit(stale),"stale/foreign ticket accepted");
    check(!one.assess({{1}},{4},{0}).compatible(),"conflicting RHS accepted");
    check(!one.assess({{2}},{3},{1}).compatible(),"fractional integer accepted");
    check(one.assess({{2}},{4},{1}).compatible(),"divisible integer rejected");
    check(!one.assess({{1,2}},{3},{0,0}).compatible(),"duplicate coordinate accepted");
    check(!one.assess({{1}},{3},{2}).compatible(),"unknown coordinate accepted");
    check(!one.assess({{1}},{3},{-1}).compatible(),"negative coordinate accepted");
    check(!one.assess({{1}}, {},{1}).compatible(),"RHS shape mismatch accepted");
    check(!invalid.assess({{1}},{3},{0}).compatible() && !zero.assess({{1}},{3},{0}).compatible(),"invalid source image accepted");
    const double beforeOld=referenceSeconds,beforeNew=quotientSeconds;
    const int n=256,h=208,m=n-h;
    ZhangExactMatrix hist(h,ZhangExactVector(n)),targets(m,ZhangExactVector(n));ZhangExactVector rhs(h);
    for(int i=0;i<h;++i){hist[i][i]=1;hist[i][h+i%m]=-1;rhs[i]=i%7;}
    for(int i=0;i<m;++i)targets[i][h+i]=1;
    domain(hist,rhs,targets,42000,32);
    check(comparisons>1000 && rejections>100 && discarded>0,"insufficient reference coverage");
    std::cout<<"PASS domains="<<domains<<" comparisons="<<comparisons<<" rejected="<<rejections<<" discarded_valid_trials="<<discarded
             <<" reference_seconds="<<referenceSeconds<<" quotient_seconds="<<quotientSeconds
             <<" benchmark_physical_variables="<<n<<" benchmark_quotient_variables="<<m
             <<" benchmark_old_seconds="<<referenceSeconds-beforeOld<<" benchmark_new_seconds="<<quotientSeconds-beforeNew<<"\n";
}
