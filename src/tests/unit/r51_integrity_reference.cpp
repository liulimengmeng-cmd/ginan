#include "common/zhangR51Integrity.hpp"
#include "common/zhangProductIntegerLedger.hpp"
#include <cassert>
#include <iostream>
static void audit(const ZhangExactMatrix& f,const ZhangExactVector& n,bool feasible) {
    const auto w=zhangR51ExactConflictWitness(f,n);
    assert(w.valid && w.feasible==feasible);
    if(feasible) return;
    ZhangExactVector row(f[0].size());ZhangExactInteger value=0;
    for(int i=0;i<f.size();++i) {
        value+=w.combination[i]*n[i];
        for(int j=0;j<row.size();++j) row[j]+=w.combination[i]*f[i][j];
    }
    assert(value==w.rhs);
    for(const auto& x:row) assert(w.modulus==0?x==0:x%w.modulus==0);
    assert(w.modulus==0?value!=0:value%w.modulus!=0);
}
int main() {
    audit({{1,0},{0,1},{1,1}},{1,2,4},false);
    audit({{1,-2},{1,0}},{0,1},false);
    audit({{2}},{6},true);audit({{2}},{5},false);
    audit({{2,4},{4,8}},{6,12},true);
    audit({{2,4},{4,8}},{6,13},false);
    for(int a=-3;a<=3;++a)for(int b=-3;b<=3;++b)
    for(int c=-2;c<=2;++c)for(int d=-2;d<=2;++d) {
        ZhangExactMatrix f={{a,b},{c,d}};
        for(int x=-2;x<=2;++x) {
            ZhangExactVector n={x,1};
            const auto exact=zhangIntegerRowLatticeContains({{a,c},{b,d}},n);
            audit(f,n,exact.contained);
        }
    }
    ProductIntegerLedger l;ProductIntegerLedgerRow row;
    row.system=E_Sys::GPS;row.physicalExpansion={{"a",1}};
    row.physicalExpansionExact=true;row.canonicalProductExpansion={{"p",1}};
    row.phaseSegmentFingerprint="segment1";row.integerValue=3;
    auto receipt=l.preflight(10,{row},1,"root1","physical1");
    assert(receipt.update.valid && l.rows().empty());
    assert(!l.commit(receipt,"root2","physical1").valid && l.rows().empty());
    assert(l.commit(receipt,"root1","physical1").valid && l.rows().size()==1);
    assert(!l.commit(receipt,"root1","physical1").valid);
    row.integerValue=4;auto failure=l.preflight(11,{row},1,"root2","physical1");
    assert(!failure.update.valid && l.rows()[0].integerValue==3);
    std::cout<<"R51 exact witnesses: 6125 exhaustive small systems and 6 examples; immutable preflight PASS\n";
}
