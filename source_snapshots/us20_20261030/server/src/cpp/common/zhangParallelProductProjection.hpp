#pragma once

#include "common/zhangQuotientIntegerLattice.hpp"
#include "common/zhangParallelProjection.hpp"
#include <stdexcept>

// Compute the image of each integer-kernel generator on every product target.
// This is an exact projection, not a feasibility or certificate test.  The
// kernel generators own disjoint output rows, so scheduling cannot change the
// integer sums or the subsequent consequence lattice.
inline ZhangExactMatrix zhangProjectTargetsOnKernel(
    const ZhangExactMatrix& targets,const std::vector<int>& columns,
    const ZhangExactMatrix& kernel)
{
    const int targetCount=targets.size();
    const int dimension=targets.empty()?0:targets.front().size();
    if(!zhangExactRectangularMatrix(targets,dimension))
        throw std::invalid_argument("PRODUCT_TARGET_KERNEL_TARGET_DIMENSION");
    for(int c:columns)if(c<0 || c>=dimension)
        throw std::invalid_argument("PRODUCT_TARGET_KERNEL_COLUMN_OUTSIDE_TARGET");
    for(const auto& row:kernel)if(row.size()!=columns.size())
        throw std::invalid_argument("PRODUCT_TARGET_KERNEL_BASIS_DIMENSION");

    std::vector<std::vector<std::pair<int,ZhangExactInteger>>> sparse(targetCount);
    std::size_t nonzeros=0;
    for(int r=0;r<targetCount;++r)
    for(int c=0;c<columns.size();++c)if(targets[r][columns[c]]!=0) {
        sparse[r].emplace_back(c,targets[r][columns[c]]);++nonzeros;
    }
    ZhangExactMatrix result(kernel.size(),ZhangExactVector(targetCount));
    const int threads=zhangR51ParallelThreads();
    // Multiprecision allocation and OpenMP startup dominate small images.
    const bool parallel=threads>1 && kernel.size()>=64 &&
        nonzeros*kernel.size()>=2000000;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
#endif
    for(int k=0;k<static_cast<int>(kernel.size());++k)
    for(int r=0;r<targetCount;++r)
    for(const auto& [c,value]:sparse[r])
        result[k][r]+=value*kernel[k][c];
    return result;
}
