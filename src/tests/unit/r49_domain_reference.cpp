#include "common/zhangR49ConditionDomain.hpp"
#include "common/zhangR49BridgeScheduler.hpp"
#include "common/zhangR48Bridge.hpp"
#include <iostream>
#include <stdexcept>
void require(bool ok){if(!ok)throw std::runtime_error("R49 domain reference failed");}
int main() {
    ZhangExactMatrix targets{{1,0,0},{0,1,0},{0,0,1}},held{{1,-1,0}};
    ZhangExactVector values{3};
    auto domain=zhangR47CompileProductSearchFrame(targets,held,values,3);
    require(domain.valid);
    require(zhangR49CanonicalTarget(domain,{{1,0,-1}},3)==zhangR49CanonicalTarget(domain,{{0,1,-1}},3));
    require(zhangR49CanonicalTarget(domain,{{2,0,-2}},3)!=zhangR49CanonicalTarget(domain,{{1,0,-1}},3));
    auto cached=zhangR49CompileTargetOnDomain({{1,0,-1},{0,1,-1}},domain,3);
    auto reference=zhangR47CompileProductSearchFrame({{1,0,-1},{0,1,-1}},held,values,3);
    require(cached.valid && reference.valid && cached.searchRank==1 && reference.searchRank==1);
    ZhangExactInteger value;
    require(zhangR47ProductConsequence(domain,{1,-1,0},0,value) && value==3);
    auto safe=values;
    require(!zhangR49CommitTrial(safe,ZhangExactVector{4},[](const auto&){return false;}) && safe==values);
    require(zhangR49CommitTrial(safe,ZhangExactVector{5},[](const auto&){return true;}) && safe[0]==5);
    int calls=0;
    auto skipped=zhangR48SearchBridge(Eigen::VectorXd::Zero(3),Eigen::MatrixXd::Identity(3,3),{}, {},{},0,1e-6,
      [&](const auto&,const auto&,double,bool){++calls;return ZhangSequentialShadowProposal{};});
    require(calls==0 && skipped.status=="NOT_EVALUATED_NO_SEARCH_BUDGET");
    std::cout<<"R49 domain reference PASS: integer translations, no 2Z/Z conflation, cached target image, atomic rollback, budget-before-compile\n";
}
