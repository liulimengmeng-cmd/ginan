#include "common/zhangRegionalJointProduct.hpp"
#include "common/zhangR51Warmup.hpp"
#include <cassert>
#include <iostream>

int main() {
    // Independent dense oracle, not another call to the sparse projection.
    Eigen::MatrixXd l(5,5);l<<2,0,0,0,0, 1,3,0,0,0, -.5,2,1,0,0,
        2,-1,.7,4,0, -.2,1,2,.5,2;
    Eigen::MatrixXd p=l*l.transpose();Eigen::VectorXd x(5);x<<1,2,3,4,5;
    const auto originalP=p;const auto originalX=x;
    std::vector<ZhangRegionalRow> rows={{{0,1}},{{1,1},{3,-.5}},{{2,.16237244751199473}},{{4,.16237244751199473}}};
    Eigen::MatrixXd h=Eigen::MatrixXd::Zero(4,5);h(0,0)=1;h(1,1)=1;h(1,3)=-.5;
    h(2,2)=.16237244751199473;h(3,4)=.16237244751199473;
    auto [m,q]=zhangRegionalProject(x,p,rows);
    assert((m-h*x).norm()<1e-12);assert((q-h*p*h.transpose()).norm()<1e-12);
    assert(q(0,2)!=0 && q(1,2)!=0 && q(2,3)!=0);
    // Full affine frontend transform; cross-ionosphere blocks must move too.
    Eigen::MatrixXd t=Eigen::MatrixXd::Identity(4,4);t(1,0)=-.3;t(0,1)=.2;
    const auto transported=zhangRegionalTransformRows(t,rows);
    auto [tm,tq]=zhangRegionalProject(x,p,transported);
    assert((tm-t*h*x).norm()<1e-12);assert((tq-t*h*p*h.transpose()*t.transpose()).norm()<1e-12);
    zhangRegionalValidate(tm,tq);
    // Hou combination must include DI and BI, not an independent-sigma sum.
    Eigen::VectorXd a(4);a<<1,-1,1.6469444444444445,0;
    const double expected=(h.transpose()*a).dot(p*(h.transpose()*a));
    assert(std::abs(a.dot(q*a)-expected)<1e-12);
    assert(std::abs(a.dot(q*a)-a.cwiseProduct(a).dot(q.diagonal()))>1e-3);
    assert((x-originalX).norm()==0 && (p-originalP).norm()==0);
    bool rejected=false;auto bad=q;bad(0,0)=-100;
    try{zhangRegionalValidate(m,bad);}catch(const std::runtime_error&){rejected=true;}assert(rejected);
    rejected=false;try{zhangRegionalProject(x,p,{{{50,1}}});}catch(const std::runtime_error&){rejected=true;}assert(rejected);
    const auto start=zhangR51ParseArStart("1405216800");
    assert(zhangR51FloatWarmup(1405216770,start));
    assert(!zhangR51FloatWarmup(1405216800,start));
    assert(!zhangR51FloatWarmup(1405216830,start));
    for(const char* s:{"nan","inf","1405216800junk","-1","1405216800.5"}) {
        rejected=false;try{zhangR51ParseArStart(s);}catch(...){rejected=true;}assert(rejected);
    }
    std::cout<<"PASS independent full Q, DI/BI/II, frontend transport, Hou combination, immutable source, reject invalid Q/index, exact 2h boundary\n";
}
