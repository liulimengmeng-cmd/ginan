#include "common/zhangR51PhysicalEntailment.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>

using Rows=ZhangR51PhysicalEntailmentSnapshot::Rows;
using Row=ZhangR51PhysicalEntailmentSnapshot::Row;
using Rational=boost::multiprecision::cpp_rational;
static void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}

// Independent Gaussian elimination, not a second call to the production HNF,
// affine quotient, kernel, or integer row-lattice membership.
static std::size_t rationalRank(std::vector<std::vector<Rational>> a){
    if(a.empty())return 0;
    const std::size_t columns=a.front().size();
    std::size_t rank=0;
    for(std::size_t c=0;c<columns && rank<a.size();++c){
        std::size_t pivot=rank;
        while(pivot<a.size() && a[pivot][c]==0)++pivot;
        if(pivot==a.size())continue;
        std::swap(a[rank],a[pivot]);
        const Rational divisor=a[rank][c];
        for(std::size_t j=c;j<columns;++j)a[rank][j]/=divisor;
        for(std::size_t i=rank+1;i<a.size();++i){
            const Rational factor=a[i][c];
            for(std::size_t j=c;j<columns;++j)a[i][j]-=factor*a[rank][j];
        }
        ++rank;
    }
    return rank;
}
// For feasible systems, (t,v) is in the rational span of (H,b) iff t*x=v
// for all integer solutions. Test fixtures know feasibility by construction.
static bool independentConsequence(const Rows& rows,const ZhangExactVector& values,
    bool knownFeasible,const Rows& targets,const ZhangExactVector& rhs){
    if(!knownFeasible || rows.size()!=values.size() || targets.size()!=rhs.size())return false;
    std::map<std::string,std::size_t> columns;
    for(const auto* group:{&rows,&targets})for(const auto& row:*group)
        for(const auto& [id,value]:row)if(value!=0 && !columns.contains(id))
            columns[id]=columns.size();
    std::vector<std::vector<Rational>> augmented;
    for(std::size_t i=0;i<rows.size();++i){
        std::vector<Rational> row(columns.size()+1);
        for(const auto& [id,value]:rows[i])if(value!=0)row[columns.at(id)]=Rational(value);
        row.back()=Rational(values[i]);augmented.push_back(std::move(row));
    }
    const auto historyRank=rationalRank(augmented);
    for(std::size_t i=0;i<targets.size();++i){
        auto trial=augmented;
        std::vector<Rational> row(columns.size()+1);
        for(const auto& [id,value]:targets[i])if(value!=0)row[columns.at(id)]=Rational(value);
        row.back()=Rational(rhs[i]);trial.push_back(std::move(row));
        if(rationalRank(std::move(trial))!=historyRank)return false;
    }
    return true;
}

int main(){
    std::size_t compared=0;
    auto compare=[&](const Rows& rows,const ZhangExactVector& values,bool knownFeasible,
        const std::vector<std::pair<Rows,ZhangExactVector>>& queries){
        ZhangR51PhysicalEntailmentSnapshot snapshot(rows,values,true);
        for(const auto& [target,rhs]:queries){
            check(snapshot.entails(target,rhs)==
                independentConsequence(rows,values,knownFeasible,target,rhs),
                "independent entailment mismatch");
            ++compared;
        }
        check(snapshot.builds<=1,"immutable scan rebuilt integer domain");
    };
    auto query=[](Row row,ZhangExactInteger rhs)
        -> std::pair<Rows,ZhangExactVector>{return {{std::move(row)},{rhs}};};
    compare({{{"x",2}}},{6},true,{query({{"x",1}},3),query({{"x",1}},2),
        query({{"retired",1}},0),query({},0),{{},{}},query({{"x",2}},6)});
    compare({{{"x",2}}},{5},false,{query({{"x",1}},2),
        query({{"x",2}},5),query({},0)});
    compare({{{"x",1},{"hidden",-2}}},{0},true,
        {query({{"x",1}},0),query({{"x",2},{"hidden",-4}},0)});
    compare({{{"x",1}},{{"x",1}}},{0,1},false,
        {query({{"x",2}},0),query({{"x",1}},0)});
    compare({}, {},true,{{{},{}},query({},0),query({{"x",1}},0)});
    ZhangExactInteger big=ZhangExactInteger(1)<<130;
    compare({{{"x",big}}},{big*3},true,
        {query({{"x",1}},3),query({{"x",1}},4)});

    check(!zhangR47AffineIntegerFeasible({{2}},{1},1),"2*x=1 accepted");
    check(zhangR47AffineIntegerFeasible({{2}},{6},1),"2*x=6 rejected");
    check(zhangR47AffineIntegerFeasible({{2,1}},{1},2),"2*x+y=1 rejected");
    check(!zhangR47AffineIntegerFeasible({{1},{1}},{0,1},1),
        "contradictory duplicate accepted");
    const ZhangExactMatrix evenGenerators{{2}};
    const auto noWitness=zhangIntegerRowLatticeContains(evenGenerators,{6},false);
    const auto withWitness=zhangIntegerRowLatticeContains(evenGenerators,{6});
    check(noWitness.contained && noWitness.combination.empty() &&
        withWitness.contained && withWitness.combination==ZhangExactVector{3},
        "feasibility-only membership changed witness semantics");
    check(!zhangExactAffineIntegerQuotient({{2}},{1},1,
        ZhangExactQuotientWork::FEASIBILITY_ONLY).valid,
        "feasibility-only quotient accepted a parity conflict");
    const auto quotient=zhangExactAffineIntegerQuotient({{2}},{6},1,
        ZhangExactQuotientWork::FEASIBILITY_ONLY);
    check(quotient.valid && quotient.particularSolution.empty() &&
        quotient.kernelBasis.empty(),"feasibility-only quotient built a solution");

    std::mt19937 random(7319);
    for(int system=0;system<60;++system){
        const int n=1+random()%5,m=random()%7;
        Rows rows;ZhangExactVector values;std::vector<int> planted(n);
        for(auto& value:planted)value=int(random()%7)-3;
        for(int i=0;i<m;++i){
            Row row;ZhangExactInteger rhs=0;
            for(int j=0;j<n;++j){
                const int coefficient=int(random()%5)-2;
                if(coefficient)row["x"+std::to_string(j)]=coefficient;
                rhs+=coefficient*planted[j];
            }
            rows.push_back(row);values.push_back(rhs);
        }
        const bool conflict=system%3==0 && !rows.empty();
        if(conflict){rows.push_back(rows.back());values.push_back(values.back()+1);}
        std::vector<std::pair<Rows,ZhangExactVector>> queries;
        for(int q=0;q<12;++q){
            Row target;ZhangExactInteger rhs=0;
            if(!rows.empty() && q%2==0){
                for(int j=0;j<m;++j){
                    const int coefficient=int(random()%3)-1;
                    for(const auto& [id,value]:rows[j])target[id]+=coefficient*value;
                    rhs+=coefficient*values[j];
                }
            }else for(int j=0;j<n;++j){
                const int coefficient=int(random()%5)-2;
                if(coefficient)target["x"+std::to_string(j)]=coefficient;
                rhs+=coefficient*planted[j];
            }
            if(q%3==0)rhs+=1;
            queries.push_back({{target},{rhs}});
        }
        compare(rows,values,!conflict,queries);
    }

    Rows mutableRows{{{"x",1}}};ZhangExactVector values{1};
    ZhangR51PhysicalEntailmentSnapshot old(mutableRows,values,true);
    values[0]=2;mutableRows[0]={{"new_arc",1}};
    ZhangR51PhysicalEntailmentSnapshot fresh(mutableRows,values,true);
    check(old.entails({{{"x",2}}},{2}),"snapshot mutated through source");
    check(!fresh.entails({{{"x",2}}},{2}),"retired coordinate leaked into new scan");
    check(fresh.entails({{{"new_arc",2}}},{4}),"new epoch rhs lost");

    Rows chain;ZhangExactVector chainValues;
    const int dimension=64,queryCount=96;
    for(int i=1;i<dimension;++i){
        chain.push_back({{"x"+std::to_string(i),1},{"x0",-1}});
        chainValues.push_back(i);
    }
    std::vector<std::pair<Rows,ZhangExactVector>> queries;
    for(int i=0;i<queryCount;++i){
        const int a=1+i%(dimension-2),b=a+1;
        queries.push_back({{{{"x"+std::to_string(a),1},
            {"x"+std::to_string(b),-1}}},{-1}});
    }
    const auto start=std::chrono::steady_clock::now();
    for(const auto& [target,rhs]:queries)
        check(independentConsequence(chain,chainValues,true,target,rhs),
            "benchmark oracle rejected a consequence");
    const auto middle=std::chrono::steady_clock::now();
    ZhangR51PhysicalEntailmentSnapshot cached(chain,chainValues,true);
    for(const auto& [target,rhs]:queries)
        check(cached.entails(target,rhs),"benchmark entailment failed");
    const auto end=std::chrono::steady_clock::now();
    check(cached.builds==1 && cached.queries==queryCount,"benchmark reuse contract");
    std::cout<<"PASS independent_reference_cases="<<compared
        <<" immutable_epoch_isolation=1 parity_and_conflicts=1\n";
    std::cout<<"BENCH dimension="<<dimension<<" queries="<<queryCount
        <<" independent_reference_seconds="<<std::chrono::duration<double>(middle-start).count()
        <<" consequence_seconds="<<std::chrono::duration<double>(end-middle).count()
        <<" hnf_seconds="<<cached.hnfSeconds
        <<" feasibility_seconds="<<cached.feasibilitySeconds<<"\n";
}
