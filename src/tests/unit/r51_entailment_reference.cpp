#include "common/zhangR51PhysicalEntailment.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
using Rows=ZhangR51PhysicalEntailmentSnapshot::Rows;
using Row=ZhangR51PhysicalEntailmentSnapshot::Row;
static void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
// The original production algorithm, copied without algebraic changes.
static bool reference(const Rows& rows,const ZhangExactVector& values,bool valid,
                      const Rows& target,const ZhangExactVector& rhs){
    if(!valid || target.size()!=rhs.size())return false;
    bool direct=true;
    for(int i=0;i<target.size();++i){bool found=false;for(int j=0;j<rows.size();++j)found|=target[i]==rows[j]&&rhs[i]==values[j];direct&=found;}
    if(direct)return true;
    std::map<std::string,int> columns;
    for(const auto* group:{&rows,&target})for(const auto& row:*group)for(const auto& [id,x]:row)
        if(x!=0 && !columns.contains(id))columns[id]=columns.size();
    auto dense=[&](const auto& group){ZhangExactMatrix out(group.size(),ZhangExactVector(columns.size()));
        for(int i=0;i<group.size();++i)for(const auto& [id,x]:group[i])if(x!=0)out[i][columns.at(id)]=x;return out;};
    const auto targets=dense(target);
    const auto frame=zhangR47CompileProductSearchFrame(targets,dense(rows),values,columns.size());
    if(!frame.valid)return false;
    for(int i=0;i<targets.size();++i){ZhangExactInteger value;if(!zhangR47ProductConsequence(frame,targets[i],0,value)||value!=rhs[i])return false;}
    return true;
}
int main(){
    std::size_t compared=0;
    auto compare=[&](const Rows& rows,const ZhangExactVector& rhs,bool valid,const std::vector<std::pair<Rows,ZhangExactVector>>& queries){
        ZhangR51PhysicalEntailmentSnapshot cached(rows,rhs,valid);
        for(const auto& [t,z]:queries){check(cached.entails(t,z)==reference(rows,rhs,valid,t,z),"cached/reference mismatch");++compared;}
        check(cached.builds<=1,"immutable scan rebuilt affine quotient");
    };
    compare({{{"x",2}}},{6},true,{{{{{"x",1}}},{3}},{{{{"x",1}}},{2}},{{{{"retired",1}}},{0}},{{{}},{0}},{{},{}},{{{{"x",2}}},{6}}});
    compare({{{"x",2}}},{5},true,{{{{{"x",1}}},{2}},{{{{"x",2}}},{5}},{{{}},{0}}});
    compare({{{"x",1},{"hidden",-2}}},{0},true,{{{{{"x",1}}},{0}},{{{{"x",2},{"hidden",-4}}},{0}}});
    compare({{{"x",1}},{{"x",1}}},{0,1},true,{{{{{"x",2}}},{0}},{{{{"x",1}}},{0}}});
    compare({}, {},true,{{{},{}},{{{}},{0}},{{{{"x",1}}},{0}}});
    compare({{{"x",1}}},{2},false,{{{{{"x",1}}},{2}}});
    ZhangExactInteger big=ZhangExactInteger(1)<<130;
    compare({{{"x",big}}},{big*3},true,{{{{{"x",1}}},{3}},{{{{"x",1}}},{4}}});
    std::mt19937 random(7319);
    for(int system=0;system<60;++system){
        const int n=1+random()%5,m=random()%7;
        Rows rows;ZhangExactVector rhs;std::vector<int> x(n);
        for(auto& v:x)v=int(random()%7)-3;
        for(int i=0;i<m;++i){Row row;ZhangExactInteger value=0;
            for(int j=0;j<n;++j){int c=int(random()%5)-2;if(c)row["x"+std::to_string(j)]=c;value+=c*x[j];}
            rows.push_back(row);rhs.push_back(value);}
        if(system%3==0 && !rhs.empty())rhs.back()+=1;
        std::vector<std::pair<Rows,ZhangExactVector>> queries;
        for(int q=0;q<12;++q){Row t;ZhangExactInteger z=0;
            if(!rows.empty() && q%2==0){for(int j=0;j<m;++j){int c=int(random()%3)-1;for(const auto& [id,a]:rows[j])t[id]+=c*a;z+=c*rhs[j];}}
            else {for(int j=0;j<n;++j){int c=int(random()%5)-2;if(c)t["x"+std::to_string(j)]=c;z+=c*x[j];}}
            if(q%3==0)z+=1;
            queries.push_back({{t},{z}});
        }
        compare(rows,rhs,true,queries);
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
    for(int i=1;i<dimension;++i){chain.push_back({{"x"+std::to_string(i),1},{"x0",-1}});chainValues.push_back(i);}
    std::vector<std::pair<Rows,ZhangExactVector>> queries;
    for(int i=0;i<queryCount;++i){int a=1+i%(dimension-2),b=a+1;queries.push_back({{{{"x"+std::to_string(a),1},{"x"+std::to_string(b),-1}}},{-1}});}
    const auto start=std::chrono::steady_clock::now();
    for(const auto& [t,z]:queries)check(reference(chain,chainValues,true,t,z),"benchmark reference");
    const auto middle=std::chrono::steady_clock::now();
    ZhangR51PhysicalEntailmentSnapshot cached(chain,chainValues,true);
    for(const auto& [t,z]:queries)check(cached.entails(t,z),"benchmark cached");
    const auto end=std::chrono::steady_clock::now();
    check(cached.builds==1 && cached.queries==queryCount,"benchmark reuse contract");
    std::cout<<"PASS reference_cases="<<compared<<" immutable_epoch_isolation=1 parity_and_conflicts=1\n";
    std::cout<<"BENCH dimension="<<dimension<<" queries="<<queryCount<<" reference_seconds="<<std::chrono::duration<double>(middle-start).count()
        <<" cached_seconds="<<std::chrono::duration<double>(end-middle).count()<<" cached_affine_builds="<<cached.builds<<"\n";
}
