#pragma once
#include "common/zhangR47ProductDomain.hpp"
using ZhangR51Rational=boost::multiprecision::cpp_rational;
using ZhangR51RationalRow=std::vector<ZhangR51Rational>;
// Kept by the exact image factory. Numerical current-domain rewrites may
// replace W/offset, but cannot silently replace the original COMPLETE integer
// image p+G*Z^r by an over-lattice when final residual admission is checked.
struct ZhangR51IntegerImageSource
{
    const ZhangExactMatrix targets,generators;
    const ZhangExactVector particularTarget;
};
struct ZhangR51PhysicalImage {
    bool valid=false;std::string reason="NOT_EVALUATED";
    std::vector<ZhangR51RationalRow> projector;
    ZhangR51RationalRow offsets;
    ZhangExactMatrix generators;ZhangExactVector particularTarget;
    std::shared_ptr<const ZhangR51IntegerImageSource> integerImageSource;
    ZhangExactColumnPlanDiagnostics columnPlan;
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
        zhangExactWorkCheckpoint("zhangR51PhysicalImage_LOOP");
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
inline ZhangR51PhysicalImage zhangR51PhysicalImageFromFrame(
    const ZhangExactMatrix& physicalTargets,const ZhangExactMatrix& numericalTargets,
    const ZhangR47ProductSearchFrame& frame)
{
    ZhangR51PhysicalImage out;
    out.columnPlan=frame.columnPlan;
    if(physicalTargets.empty() || physicalTargets.size()!=numericalTargets.size())return out;
    if(!frame.valid || !frame.affine || !frame.affine->valid){out.reason=frame.reason;return out;}
    if(!zhangExactRectangularMatrix(numericalTargets,numericalTargets.front().size()))
    {out.reason="NUMERICAL_TARGET_SHAPE_INVALID";return out;}
    const int m=physicalTargets.size(),n=numericalTargets.front().size(),rank=frame.searchRank;
    out.particularTarget.resize(m);
    for(int i=0;i<m;++i)for(int c=0;c<frame.columns.size();++c) {
        zhangExactWorkCheckpoint("PHYSICAL_IMAGE_PARTICULAR_TARGET");
        out.particularTarget[i]+=physicalTargets[i][frame.columns[c]]*frame.affine->particularSolution[c];
    }
    out.generators=frame.imageGenerators;
    out.projector.assign(rank,ZhangR51RationalRow(n));out.offsets.resize(rank);
    // imageGenerators is the transpose of row HNF and has ordered pivot rows.
    // Its exact triangular left inverse retains divisibility (e.g. a=2q).
    std::vector<ZhangR51RationalRow> inverse(rank,ZhangR51RationalRow(m));
    for(int k=0;k<rank;++k) {
        zhangExactWorkCheckpoint("zhangR51PhysicalImage_LOOP");
        int pivot=-1;for(int i=0;i<m;++i)if(out.generators[i][k]!=0){pivot=i;break;}
        if(pivot<0){out.reason="IMAGE_PIVOT_MISSING";return out;}
        inverse[k][pivot]=1;
        for(int j=0;j<k;++j)if(out.generators[pivot][j]!=0)for(int i=0;i<m;++i) {
            zhangExactWorkCheckpoint("PHYSICAL_IMAGE_RATIONAL_INVERSE");
            if(inverse[j][i]!=0)inverse[k][i]-=ZhangR51Rational(out.generators[pivot][j])*inverse[j][i];
        }
        for(auto& x:inverse[k]) {
            zhangExactWorkCheckpoint("PHYSICAL_IMAGE_RATIONAL_DIVISION");
            x/=ZhangR51Rational(out.generators[pivot][k]);
        }
    }
    for(int i=0;i<rank;++i) {
        zhangExactWorkCheckpoint("zhangR51PhysicalImage_LOOP");
        ZhangR51RationalRow check(rank);
        for(int k=0;k<m;++k)if(inverse[i][k]!=0)for(int j=0;j<rank;++j)if(out.generators[k][j]!=0)
            check[j]+=inverse[i][k]*ZhangR51Rational(out.generators[k][j]);
        for(int j=0;j<rank;++j)if(check[j]!=(i==j?1:0)){out.reason="IMAGE_LEFT_INVERSE_AUDIT_FAILED";return out;}
        for(int k=0;k<m;++k)if(inverse[i][k]!=0) {
            out.offsets[i]-=inverse[i][k]*ZhangR51Rational(out.particularTarget[k]);
            for(int c=0;c<n;++c) {
                zhangExactWorkCheckpoint("PHYSICAL_IMAGE_TARGET_PROJECTION");
                if(numericalTargets[k][c]!=0)out.projector[i][c]+=inverse[i][k]*ZhangR51Rational(numericalTargets[k][c]);
            }
        }
    }
    out.integerImageSource=std::make_shared<const ZhangR51IntegerImageSource>(
        ZhangR51IntegerImageSource{physicalTargets,out.generators,out.particularTarget});
    out.valid=true;out.reason=rank?"EXACT_PHYSICAL_AFFINE_IMAGE":"PHYSICAL_TARGET_ALREADY_DETERMINED";return out;
}
inline ZhangR51PhysicalImage zhangR51PhysicalImage(
    const ZhangExactMatrix& physicalTargets,const ZhangExactMatrix& history,
    const ZhangExactVector& values,const ZhangExactMatrix& numericalTargets,int physicalDimension,
    ZhangProductFrameWork work=ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
{
    const auto frame=zhangR47CompileProductSearchFrame(
        physicalTargets,history,values,physicalDimension,work);
    return zhangR51PhysicalImageFromFrame(physicalTargets,numericalTargets,frame);
}
inline ZhangR51PhysicalImage zhangR51PhysicalImageOnVerifiedDomain(
    const ZhangExactMatrix& physicalTargets,
    std::shared_ptr<const ZhangExactAffineIntegerQuotient> domain,
    const ZhangExactMatrix& numericalTargets,int physicalDimension,
    ZhangProductFrameWork work=ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
{
    const auto frame=zhangR47CompileTargetOnVerifiedDomain(
        physicalTargets,std::move(domain),physicalDimension,work);
    return zhangR51PhysicalImageFromFrame(physicalTargets,numericalTargets,frame);
}
inline bool zhangR51LiftRationalRows(const ZhangExactMatrix& transform,const ZhangExactVector& integers,
    const std::vector<ZhangR51RationalRow>& projection,const ZhangR51RationalRow& offsets,
    ZhangExactMatrix& rows,ZhangExactVector& values)
{
    if(transform.empty() || projection.empty() || transform.size()!=integers.size())return false;
    const int n=projection.front().size();rows.clear();values.clear();
    for(int i=0;i<transform.size();++i) {
        zhangExactWorkCheckpoint("zhangR51PhysicalImage_LOOP");
        if(transform[i].size()!=projection.size())return false;
        ZhangR51RationalRow row(n);ZhangR51Rational value=integers[i];
        for(int k=0;k<projection.size();++k)if(transform[i][k]!=0) {
            value-=ZhangR51Rational(transform[i][k])*offsets[k];
            for(int c=0;c<n;++c) {
                zhangExactWorkCheckpoint("PHYSICAL_IMAGE_RATIONAL_LIFT");
                if(projection[k][c]!=0)row[c]+=ZhangR51Rational(transform[i][k])*projection[k][c];
            }
        }
        ZhangExactInteger scale=denominator(value);
        for(const auto& x:row) {
            zhangExactWorkCheckpoint("zhangR51PhysicalImage_LOOP");const ZhangExactInteger d=denominator(x);scale=scale/boost::multiprecision::gcd(scale,d)*d;}
        ZhangExactVector exact(n);for(int c=0;c<n;++c) {
            zhangExactWorkCheckpoint("PHYSICAL_IMAGE_DENOMINATOR_CLEARING");
            exact[c]=numerator(ZhangR51Rational(row[c]*scale));
        }
        rows.push_back(exact);values.push_back(numerator(ZhangR51Rational(value*scale)));
    }
    return true;
}

struct ZhangR51ResidualImageAdmission
{
    bool valid=false;
    std::string reason="NOT_EVALUATED";
    ZhangExactVector residualSolution,targetValues;
};

/** Final admission in the ORIGINAL complete integer target image
 * T*N=p+G*y, y in Z^r. Only the small searched r*r system Z*y=z is solved.
 * The input must come from the trusted complete integer-image factory; the
 * immutable source receipt below guards G/p/T against accidental replacement.
 * The identity T=G*W+constant alone would NOT prove a complete integer image.
 * Before using it, verify the current HNF proves T=G*W+p+G*offset and audit
 * the supplied denominator-cleared lift. This remains valid after a bound
 * current-domain numerical rewrite of W/offset. Non-unimodular Z is allowed
 * only if its UNIQUE solution is integer. No full physical HNF, feasibility
 * factorisation, kernel or particular origin is reconstructed here. */
inline ZhangR51ResidualImageAdmission zhangR51CertifyResidualImageAdmission(
    const ZhangR51PhysicalImage& image,const ZhangExactMatrix& originalTargets,
    const ZhangExactMatrix& currentBasis,const ZhangExactVector& currentValues,
    const ZhangExactMatrix& searchRows,const ZhangExactVector& searchValues,
    const ZhangExactMatrix& liftedRows,const ZhangExactVector& liftedValues)
{
    ZhangR51ResidualImageAdmission out;
    if (!image.integerImageSource ||
        image.integerImageSource->targets!=originalTargets ||
        image.integerImageSource->generators!=image.generators ||
        image.integerImageSource->particularTarget!=image.particularTarget)
    {out.reason="RESIDUAL_ORIGINAL_INTEGER_IMAGE_SOURCE_MISMATCH";return out;}
    if (!image.valid || originalTargets.empty() ||
        originalTargets.size()!=image.particularTarget.size() ||
        (!image.projector.empty() && originalTargets.size()!=image.generators.size()) ||
        image.projector.size()!=image.offsets.size() ||
        currentBasis.size()!=currentValues.size())
    {out.reason="RESIDUAL_IMAGE_SHAPE_INVALID";return out;}
    const auto dimension=originalTargets.front().size();
    const auto rank=image.projector.size();
    if (!zhangExactRectangularMatrix(originalTargets,dimension) ||
        !zhangExactRectangularMatrix(currentBasis,dimension) ||
        !zhangExactRectangularMatrix(image.generators,rank) ||
        searchRows.size()!=rank || searchValues.size()!=rank ||
        !zhangExactRectangularMatrix(searchRows,rank) ||
        liftedRows.size()!=rank || liftedValues.size()!=rank ||
        !zhangExactRectangularMatrix(liftedRows,dimension))
    {out.reason="RESIDUAL_IMAGE_SHAPE_INVALID";return out;}
    for (const auto& row:image.projector) if (row.size()!=dimension)
    {out.reason="RESIDUAL_IMAGE_PROJECTOR_SHAPE_INVALID";return out;}
    std::vector<std::size_t> pivots;
    for (const auto& row:currentBasis)
    {
        const auto entry=std::find_if(row.begin(),row.end(),[](const auto& x){return x!=0;});
        if (entry==row.end() || *entry<0)
        {out.reason="RESIDUAL_SOURCE_BASIS_NOT_HNF";return out;}
        const auto pivot=std::size_t(entry-row.begin());
        if (!pivots.empty() && pivot<=pivots.back())
        {out.reason="RESIDUAL_SOURCE_BASIS_NOT_HNF";return out;}
        pivots.push_back(pivot);
    }
    for (std::size_t t=0;t<originalTargets.size();++t)
    {
        ZhangR51RationalRow remainder(dimension);
        for (std::size_t c=0;c<dimension;++c) remainder[c]=originalTargets[t][c];
        ZhangR51Rational rhs=image.particularTarget[t];
        for (std::size_t k=0;k<rank;++k) if (image.generators[t][k]!=0)
        {
            const ZhangR51Rational g=image.generators[t][k];
            rhs+=g*image.offsets[k];
            for (std::size_t c=0;c<dimension;++c)
            {
                zhangExactWorkCheckpoint("RESIDUAL_ORIGINAL_TARGET_REPRESENTATION");
                remainder[c]-=g*image.projector[k][c];
            }
        }
        for (std::size_t h=0;h<currentBasis.size();++h)
        {
            const auto p=pivots[h];
            const ZhangR51Rational factor=remainder[p]/ZhangR51Rational(currentBasis[h][p]);
            if (factor==0) continue;
            rhs-=factor*ZhangR51Rational(currentValues[h]);
            for (std::size_t c=p;c<dimension;++c)
            {
                zhangExactWorkCheckpoint("RESIDUAL_CURRENT_TARGET_REDUCTION");
                if (currentBasis[h][c]!=0) remainder[c]-=factor*ZhangR51Rational(currentBasis[h][c]);
            }
        }
        if (rhs!=0 || std::any_of(remainder.begin(),remainder.end(),[](const auto& x){return x!=0;}))
        {out.reason="RESIDUAL_ORIGINAL_TARGET_REPRESENTATION_MISMATCH";return out;}
    }
    if (rank==0)
    {out.targetValues=image.particularTarget;out.valid=true;out.reason="EXACT_ALREADY_DETERMINED_IMAGE";return out;}
    ZhangExactMatrix expectedRows;ZhangExactVector expectedValues;
    if (!zhangR51LiftRationalRows(searchRows,searchValues,image.projector,image.offsets,
        expectedRows,expectedValues) || expectedRows!=liftedRows || expectedValues!=liftedValues)
    {out.reason="RESIDUAL_EXACT_LIFT_MISMATCH";return out;}
    // Exact Gaussian elimination in residual coordinates only. It is not a
    // replacement for the source-domain integer feasibility certification.
    std::vector<ZhangR51RationalRow> system(rank,ZhangR51RationalRow(rank+1));
    for (std::size_t r=0;r<rank;++r)
    {
        for (std::size_t c=0;c<rank;++c) system[r][c]=searchRows[r][c];
        system[r][rank]=searchValues[r];
    }
    for (std::size_t c=0;c<rank;++c)
    {
        std::size_t pivot=c;
        while (pivot<rank && system[pivot][c]==0) ++pivot;
        if (pivot==rank) {out.reason="RESIDUAL_SEARCH_NOT_FULL_RANK";return out;}
        std::swap(system[c],system[pivot]);
        const ZhangR51Rational divisor=system[c][c];
        for (std::size_t k=c;k<=rank;++k) system[c][k]/=divisor;
        for (std::size_t r=0;r<rank;++r) if (r!=c && system[r][c]!=0)
        {
            const ZhangR51Rational factor=system[r][c];
            for (std::size_t k=c;k<=rank;++k)
            {
                zhangExactWorkCheckpoint("RESIDUAL_SMALL_EXACT_SYSTEM");
                system[r][k]-=factor*system[c][k];
            }
        }
    }
    for (std::size_t r=0;r<rank;++r)
    {
        if (denominator(system[r][rank])!=1)
        {out.reason="RESIDUAL_SEARCH_NO_INTEGER_SOLUTION";return out;}
        out.residualSolution.push_back(numerator(system[r][rank]));
    }
    if (zhangExactMatrixTimesColumn(searchRows,out.residualSolution)!=searchValues)
    {out.reason="RESIDUAL_SEARCH_SOLUTION_AUDIT_FAILED";return out;}
    out.targetValues=image.particularTarget;
    const auto delta=zhangExactMatrixTimesColumn(image.generators,out.residualSolution);
    for (std::size_t t=0;t<out.targetValues.size();++t) out.targetValues[t]+=delta[t];
    out.valid=true;out.reason="EXACT_RESIDUAL_INTEGER_IMAGE_ADMITTED";
    return out;
}
