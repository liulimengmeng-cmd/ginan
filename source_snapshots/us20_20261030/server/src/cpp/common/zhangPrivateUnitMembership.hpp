#pragma once
#include "common/zhangIntegerAudit.hpp"
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif

struct ZhangPrivateUnitMembershipMetrics
{
    std::size_t fastBatches=0,fallbackBatches=0,targets=0,maximumBasisRows=0;
    int maximumWorkers=1;
};
inline thread_local ZhangPrivateUnitMembershipMetrics zhangPrivateUnitMembershipMetrics;

// Explicit integer witnesses for generators with distinct private +/-1
// columns. This is a proved structural shortcut, not a real-rank test or a
// saturated approximation. Non-unit/overlapping inputs use the old backend.
inline std::vector<ZhangIntegerLatticeMembership> zhangPrivateUnitLatticeMembershipBatch(
    const ZhangExactMatrix& rows,const ZhangExactMatrix& targets)
{
    if(targets.empty())return {};
    const auto dimension=targets.front().size(),rowCount=rows.size();
    for(const auto& target:targets)if(target.size()!=dimension)
        return std::vector<ZhangIntegerLatticeMembership>(targets.size());
    for(const auto& row:rows)if(row.size()!=dimension)
        return std::vector<ZhangIntegerLatticeMembership>(targets.size());
    if(dimension==0)return zhangIntegerRowLatticeContainsBatch(rows,targets);
    std::vector<int> owner(dimension,-1),pivots(rowCount,-1);
    std::vector<std::vector<std::pair<std::size_t,ZhangExactInteger>>> sparse(rowCount);
    for(std::size_t r=0;r<rowCount;++r)for(std::size_t c=0;c<dimension;++c) {
        zhangExactWorkCheckpoint("PRIVATE_UNIT_MEMBERSHIP_SCAN");
        if(rows[r][c]==0)continue;
        sparse[r].emplace_back(c,rows[r][c]);
        owner[c]=owner[c]==-1?static_cast<int>(r):-2;
    }
    for(std::size_t c=0;c<dimension;++c)if(owner[c]>=0) {
        const auto r=owner[c];
        if(pivots[r]<0 && (rows[r][c]==1 || rows[r][c]==-1))pivots[r]=c;
    }
    if(!std::all_of(pivots.begin(),pivots.end(),[](int c){return c>=0;})) {
        ++zhangPrivateUnitMembershipMetrics.fallbackBatches;
        return zhangIntegerRowLatticeContainsBatch(rows,targets);
    }
    ++zhangPrivateUnitMembershipMetrics.fastBatches;
    zhangPrivateUnitMembershipMetrics.targets+=targets.size();
    zhangPrivateUnitMembershipMetrics.maximumBasisRows=std::max(
        zhangPrivateUnitMembershipMetrics.maximumBasisRows,rowCount);
    std::vector<ZhangIntegerLatticeMembership> results(targets.size());
    auto solveTarget=[&](std::size_t t) {
        ZhangExactVector combination(rowCount),rebuilt(dimension);
        for(std::size_t r=0;r<rowCount;++r) {
            zhangExactWorkCheckpoint("PRIVATE_UNIT_MEMBERSHIP_WITNESS");
            combination[r]=targets[t][pivots[r]]*rows[r][pivots[r]];
            if(combination[r]!=0)for(const auto& [c,value]:sparse[r])
                rebuilt[c]+=combination[r]*value;
        }
        auto& result=results[t];result.rank=rowCount;
        result.smithInvariants=ZhangExactVector(rowCount,1);
        result.contained=rebuilt==targets[t];
        if(result.contained)result.combination=std::move(combination);
    };
#ifdef _OPENMP
    // An OpenMP `if(false)` still creates a serialized region across which a
    // C++ cancellation cannot escape safely. Budgeted calls use no region at
    // all, and unexpected worker exceptions are rethrown on the owning thread.
    const bool parallel=targets.size()>=8 && rowCount>=64 &&
        !zhangActiveExactWorkBudget && !omp_in_parallel();
    if(parallel) {
        int workers=1;
        std::vector<std::exception_ptr> errors(targets.size());
#pragma omp parallel for schedule(static) reduction(max:workers)
        for(std::ptrdiff_t t=0;t<static_cast<std::ptrdiff_t>(targets.size());++t) {
            workers=std::max(workers,omp_get_num_threads());
            try {solveTarget(t);}catch(...) {errors[t]=std::current_exception();}
        }
        zhangPrivateUnitMembershipMetrics.maximumWorkers=std::max(
            zhangPrivateUnitMembershipMetrics.maximumWorkers,workers);
        for(const auto& error:errors)if(error)std::rethrow_exception(error);
    } else
#endif
        for(std::size_t t=0;t<targets.size();++t)solveTarget(t);
    return results;
}
