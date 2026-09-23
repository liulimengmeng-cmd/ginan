#pragma once
#include "common/zhangR47ProductDomain.hpp"
using ZhangR51Rational=boost::multiprecision::cpp_rational;
using ZhangR51RationalRow=std::vector<ZhangR51Rational>;
struct ZhangR51PhysicalImage {
    bool valid=false;std::string reason="NOT_EVALUATED";
    std::vector<ZhangR51RationalRow> projector;
    ZhangR51RationalRow offsets;
    ZhangExactMatrix generators;ZhangExactVector particularTarget;
};
// A sparse exact right-inverse witness proves an unconstrained target is
// already primitive: choose a private unit column for each target row.
// No floating rank test and no relaxation of a nonprimitive image is used.
inline ZhangR51PhysicalImage zhangR51UnconstrainedUnitImage(
    const std::vector<std::map<std::string,ZhangExactInteger>>& physical,
    const ZhangExactMatrix& numerical)
{
    ZhangR51PhysicalImage out;
    if(physical.empty() || physical.size()!=numerical.size())return out;
    std::map<std::string,int> counts;
    for(const auto& row:physical)for(const auto& [id,x]:row)if(x!=0)++counts[id];
    for(const auto& row:physical) {
        bool unit=false;for(const auto& [id,x]:row)unit |= (x==1 || x==-1) && counts[id]==1;
        if(!unit){out.reason="NO_PRIVATE_UNIT_COLUMN_WITNESS";return out;}
    }
    const int m=numerical.size(),n=numerical.front().size();
    if(!zhangExactRectangularMatrix(numerical,n))return out;
    out.projector.assign(m,ZhangR51RationalRow(n));out.offsets.resize(m);out.particularTarget.resize(m);
    out.generators=zhangExactIdentityMatrix(m);
    for(int i=0;i<m;++i)for(int c=0;c<n;++c)if(numerical[i][c]!=0)out.projector[i][c]=numerical[i][c];
    out.valid=true;out.reason="EXACT_PRIVATE_UNIT_COLUMN_RIGHT_INVERSE";return out;
}
// The old physical variables are existential integers, including variables no
// longer in the KF. The image of their affine solution is an integer lattice,
// which can be nonprimitive in the visible target coordinates.
inline ZhangR51PhysicalImage zhangR51PhysicalImage(
    const ZhangExactMatrix& physicalTargets,const ZhangExactMatrix& history,
    const ZhangExactVector& values,const ZhangExactMatrix& numericalTargets,int physicalDimension,
    ZhangProductFrameWork work=ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
{
    ZhangR51PhysicalImage out;
    if(physicalTargets.empty() || physicalTargets.size()!=numericalTargets.size())return out;
    const auto frame=zhangR47CompileProductSearchFrame(physicalTargets,history,values,physicalDimension,work);
    if(!frame.valid || !frame.affine || !frame.affine->valid){out.reason=frame.reason;return out;}
    const int m=physicalTargets.size(),n=numericalTargets.front().size(),rank=frame.searchRank;
    out.particularTarget.resize(m);
    for(int i=0;i<m;++i)for(int c=0;c<frame.columns.size();++c)
        out.particularTarget[i]+=physicalTargets[i][frame.columns[c]]*frame.affine->particularSolution[c];
    out.generators=frame.imageGenerators;
    out.projector.assign(rank,ZhangR51RationalRow(n));out.offsets.resize(rank);
    // imageGenerators is the transpose of row HNF and has ordered pivot rows.
    // Its exact triangular left inverse retains divisibility (e.g. a=2q).
    std::vector<ZhangR51RationalRow> inverse(rank,ZhangR51RationalRow(m));
    for(int k=0;k<rank;++k) {
        int pivot=-1;for(int i=0;i<m;++i)if(out.generators[i][k]!=0){pivot=i;break;}
        if(pivot<0){out.reason="IMAGE_PIVOT_MISSING";return out;}
        inverse[k][pivot]=1;
        for(int j=0;j<k;++j)if(out.generators[pivot][j]!=0)for(int i=0;i<m;++i)if(inverse[j][i]!=0)
            inverse[k][i]-=ZhangR51Rational(out.generators[pivot][j])*inverse[j][i];
        for(auto& x:inverse[k])x/=ZhangR51Rational(out.generators[pivot][k]);
    }
    for(int i=0;i<rank;++i) {
        ZhangR51RationalRow check(rank);
        for(int k=0;k<m;++k)if(inverse[i][k]!=0)for(int j=0;j<rank;++j)if(out.generators[k][j]!=0)
            check[j]+=inverse[i][k]*ZhangR51Rational(out.generators[k][j]);
        for(int j=0;j<rank;++j)if(check[j]!=(i==j?1:0)){out.reason="IMAGE_LEFT_INVERSE_AUDIT_FAILED";return out;}
        for(int k=0;k<m;++k)if(inverse[i][k]!=0) {
            out.offsets[i]-=inverse[i][k]*ZhangR51Rational(out.particularTarget[k]);
            for(int c=0;c<n;++c)if(numericalTargets[k][c]!=0)out.projector[i][c]+=inverse[i][k]*ZhangR51Rational(numericalTargets[k][c]);
        }
    }
    out.valid=true;out.reason=rank?"EXACT_PHYSICAL_AFFINE_IMAGE":"PHYSICAL_TARGET_ALREADY_DETERMINED";return out;
}
inline bool zhangR51LiftRationalRows(const ZhangExactMatrix& transform,const ZhangExactVector& integers,
    const std::vector<ZhangR51RationalRow>& projection,const ZhangR51RationalRow& offsets,
    ZhangExactMatrix& rows,ZhangExactVector& values)
{
    if(transform.empty() || projection.empty() || transform.size()!=integers.size())return false;
    const int n=projection.front().size();rows.clear();values.clear();
    for(int i=0;i<transform.size();++i) {
        if(transform[i].size()!=projection.size())return false;
        ZhangR51RationalRow row(n);ZhangR51Rational value=integers[i];
        for(int k=0;k<projection.size();++k)if(transform[i][k]!=0) {
            value-=ZhangR51Rational(transform[i][k])*offsets[k];
            for(int c=0;c<n;++c)if(projection[k][c]!=0)row[c]+=ZhangR51Rational(transform[i][k])*projection[k][c];
        }
        ZhangExactInteger scale=denominator(value);
        for(const auto& x:row) {const ZhangExactInteger d=denominator(x);scale=scale/boost::multiprecision::gcd(scale,d)*d;}
        ZhangExactVector exact(n);for(int c=0;c<n;++c)exact[c]=numerator(ZhangR51Rational(row[c]*scale));
        rows.push_back(exact);values.push_back(numerator(ZhangR51Rational(value*scale)));
    }
    return true;
}
