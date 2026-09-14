#include "common/zhangR51PhysicalImage.hpp"
#include "common/zhangProductIntegerLedger.hpp"
#include <iostream>
#include <stdexcept>
void check(bool x,const char* msg){if(!x)throw std::runtime_error(msg);}
int main(){
    // Hidden integer r imposes parity on visible a; it must not disappear.
    auto parity=zhangR51PhysicalImage({{1,0}},{{1,-2}},{0},{{1}},2);
    check(parity.valid && parity.projector.size()==1,"parity image");
    check(parity.projector[0][0]==ZhangR51Rational(1)/2,"parity coordinate a/2");
    ZhangExactMatrix rows;ZhangExactVector rhs;
    check(zhangR51LiftRationalRows({{1}},{3},parity.projector,parity.offsets,rows,rhs),"lift");
    check(rows==ZhangExactMatrix{{1}} && rhs==ZhangExactVector{6},"q=3 implies a=6");
    auto fixed=zhangR51PhysicalImage({{1}},{{2}},{6},{{1}},1);
    check(fixed.valid && fixed.projector.empty() && fixed.particularTarget[0]==3,"2a=6 determined");
    auto bad=zhangR51PhysicalImage({{1}},{{2}},{5},{{1}},1);check(!bad.valid,"2a=5 infeasible");
    auto two=zhangR51PhysicalImage({{2,0},{1,3}}, {},{},{{1,0},{0,1}},2);
    std::cerr<<two.reason<<" rank="<<two.projector.size()<<"\n";for(const auto& row:two.generators){for(const auto& x:row)std::cerr<<x<<",";std::cerr<<"\n";} check(two.valid && two.projector.size()==2,"nonsaturated coupled image");
    // Verify every lattice point in a small cube reconstructs integer image coordinates.
    for(int a=-5;a<=5;++a)for(int b=-5;b<=5;++b)for(int i=0;i<2;++i){
        ZhangR51Rational q=two.projector[i][0]*(2*a)+two.projector[i][1]*(a+3*b)+two.offsets[i];
        check(denominator(q)==1,"congruence coordinate integer");
    }
    std::cout<<"R51 physical affine images, hidden-variable parity, determined directions and rational lifts PASS\n";
}
