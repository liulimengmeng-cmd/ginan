#pragma once
#include "common/zhangR47Candidate.hpp"

struct ZhangR47ProductSearchFrame
{
    bool valid=false;
    int targetRank=0, conditionerRank=0, searchRank=0;
    std::vector<int> columns;
    ZhangExactAffineIntegerQuotient affine;
    ZhangExactMatrix projector, imageGenerators;
    ZhangExactVector offsets;
    std::string reason="NOT_EVALUATED";
};

/** Parameterise the product image on the feasible integer affine set.
 * Only columns touched by the products/conditioners enter exact elimination;
 * ILS dimension is rank(T K), never the number of network cycles. HNF column
 * image generators retain all Smith divisibility without treating 2Z as Z. */
inline ZhangR47ProductSearchFrame zhangR47CompileProductSearchFrame(
    const ZhangExactMatrix& targets,const ZhangExactMatrix& conditioners,
    const ZhangExactVector& values,int dimension)
{
    ZhangR47ProductSearchFrame out;
    if(targets.empty() || !zhangExactRectangularMatrix(targets,dimension) ||
       !zhangExactRectangularMatrix(conditioners,dimension)) return out;
    for(int c=0;c<dimension;++c)
    {
        bool used=false;
        for(const auto& row:targets) used |= row[c]!=0;
        for(const auto& row:conditioners) used |= row[c]!=0;
        if(used) out.columns.push_back(c);
    }
    if(out.columns.empty()) {out.valid=true;out.reason="ZERO_PRODUCT_IMAGE";return out;}
    auto compact=[&](const ZhangExactMatrix& input)
    {
        ZhangExactMatrix result(input.size(),ZhangExactVector(out.columns.size()));
        for(std::size_t r=0;r<input.size();++r) for(std::size_t c=0;c<out.columns.size();++c)
            result[r][c]=input[r][out.columns[c]];
        return result;
    };
    const auto t=compact(targets),h=compact(conditioners);
    out.targetRank=zhangExactRowHermiteNormalForm(t).basis.size();
    out.affine=zhangExactAffineIntegerQuotient(h,values,out.columns.size());
    if(!out.affine.valid) {out.reason=out.affine.failureReason;return out;}
    out.conditionerRank=out.affine.deterministicRank;
    if(!out.affine.quotientRank) {out.valid=true;out.reason="PRODUCT_FULLY_DETERMINED";return out;}
    ZhangExactMatrix k(out.columns.size(),ZhangExactVector(out.affine.quotientRank));
    for(int r=0;r<out.affine.quotientRank;++r) for(int c=0;c<out.columns.size();++c)
        k[c][r]=out.affine.kernelBasis[r][c];
    const auto m=zhangExactMultiply(t,k);
    const auto image=zhangPrimitiveImageCoordinates(m);
    if(!image.valid) {out.reason="PRODUCT_IMAGE_COORDINATES_FAILED";return out;}
    out.imageGenerators=image.imageGenerators;
    out.searchRank=image.primitiveRows.size();
    const auto projected=zhangExactMultiply(image.primitiveRows,out.affine.quotientProjector);
    for(const auto& row:projected)
    {
        ZhangExactVector full(dimension); ZhangExactInteger offset=0;
        for(std::size_t c=0;c<out.columns.size();++c)
        {full[out.columns[c]]=row[c];offset-=row[c]*out.affine.particularSolution[c];}
        out.projector.push_back(std::move(full));out.offsets.push_back(offset);
    }
    out.valid=true;out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
    return out;
}

inline bool zhangR47ProductConsequence(const ZhangR47ProductSearchFrame& frame,
    const ZhangExactVector& target,const ZhangExactInteger& offset,ZhangExactInteger& value)
{
    if(!frame.valid || !frame.affine.valid) return false;
    std::vector<bool> represented(target.size());
    ZhangExactVector compact(frame.columns.size());
    for(std::size_t c=0;c<frame.columns.size();++c)
    {
        if(frame.columns[c]>=target.size()) return false;
        compact[c]=target[frame.columns[c]];represented[frame.columns[c]]=true;
    }
    for(std::size_t c=0;c<target.size();++c) if(target[c]!=0 && !represented[c]) return false;
    for(const auto& kernel:frame.affine.kernelBasis)
    {
        ZhangExactInteger residual=0;
        for(std::size_t c=0;c<compact.size();++c) residual+=compact[c]*kernel[c];
        if(residual!=0) return false;
    }
    value=offset;
    for(std::size_t c=0;c<compact.size();++c) value+=compact[c]*frame.affine.particularSolution[c];
    return true;
}
