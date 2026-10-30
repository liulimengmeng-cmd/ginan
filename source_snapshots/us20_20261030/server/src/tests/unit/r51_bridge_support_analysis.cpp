#include "common/zhangBridgeSupportAnalysis.hpp"
#include <iostream>
#include <stdexcept>

static void check(bool value, const char* why)
{ if (!value) throw std::runtime_error(why); }

int main()
{
    // Independent algebra: 2*old+live=1 => old=(1-live)/2. The parity
    // restriction belongs to the ORIGINAL integer image, not this diagnostic.
    auto a = zhangAnalyzeBridgeSupport({{2,1}}, {{1,0}}, {0}, 2);
    check(a.valid && a.representable && a.missingRank == 0, "rational old rewrite");
    auto b = zhangAnalyzeBridgeSupport({{0,1}}, {{1,-1}}, {0}, 2);
    check(b.valid && !b.representable && b.missingRank == 1, "free retired arc");
    check(b.counterexamples == ZhangExactMatrix{{1,0}}, "independent missing-direction oracle");
    auto c = zhangAnalyzeBridgeSupport({{1,2,1}}, {{0,1,0},{1,2,0}}, {0,1}, 3);
    check(c.valid && !c.representable && c.missingRank == 1, "shared missing direction");
    check(c.counterexamples.front() == ZhangExactVector({-2,1,0}), "nonunit exact null witness");
    check(zhangAnalyzeBridgeSupport({}, {{1,2},{2,4}}, {0,1}, 2).missingRank == 1,
        "dependent targets do not double count deficiency");
    check(!zhangAnalyzeBridgeSupport({{1}}, {{1,2}}, {0}, 2).valid, "malformed reject");
    // Independent bounded integer enumeration with live held fixed. Whenever
    // a witness is returned it must preserve H and alter the target exactly.
    std::size_t cases = 0;
    for (int x=-3;x<=3;++x) for(int y=-3;y<=3;++y)
    for (int p=-2;p<=2;++p) for(int q=-2;q<=2;++q)
    {
        const auto result=zhangAnalyzeBridgeSupport({{x,y,1}}, {{p,q,0}}, {0,1}, 3);
        check(result.valid,"small exact result");
        const bool expected = x==0 && y==0 ? p==0 && q==0 : x*q==y*p;
        check(result.representable==expected,"independent determinant oracle");
        for(const auto& v:result.counterexamples)
            check(x*v[0]+y*v[1]==0 && v[2]==0 && p*v[0]+q*v[1]!=0,"independent witness audit");
        ++cases;
    }
    std::cout<<"bridge support exact null-direction tests PASS cases="<<cases<<"\n";
}
