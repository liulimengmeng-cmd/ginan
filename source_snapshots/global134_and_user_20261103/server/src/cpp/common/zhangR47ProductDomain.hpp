#pragma once
#include "common/zhangR47Candidate.hpp"
#include <memory>
#include <mutex>
#include <cstdlib>
#include "common/zhangExactColumnPlan.hpp"
#include "common/zhangProductFrameReplay.hpp"

struct ZhangExactColumnPlanDiagnostics
{
    bool sharedPlanHit=false,targetReplayHit=false,legacyResourceFallback=false,targetComplete=false,
        legacyOverride=false;
    std::size_t buildCalls=0,columnOperations=0,journalBytes=0;
};

/** Cache by the actual immutable affine owner, never an epoch, pointer value
 * alone, or a hash of a subset of rows. The current-domain caller validates
 * its full runtime/generation identity before handing out this aliasing view.
 * Expired owners are evicted; at most two bounded plans and eight small
 * target replays per plan survive. No historical integer RHS is retained
 * after its domain owner dies. */
inline ZhangExactMatrix zhangProjectedKernelOnOwnedDomain(
    const std::shared_ptr<const ZhangExactAffineIntegerQuotient>& domain,
    const ZhangExactMatrix& targets,ZhangExactColumnPlanDiagnostics* diagnostics=nullptr)
{
    if (!domain || !domain->valid) return {};
    if (const char* mode=std::getenv("ZHANG_R51_BRIDGE_COLUMN_PLAN");
        mode && std::string(mode)=="legacy")
    {
        ZhangPhaseTimer timer("BRIDGE_SOURCE_COLUMN_PLAN_LEGACY");
        zhangExactLibraryBoundary("COLUMN_PLAN_LEGACY_ENTRY");
        auto extended=domain->deterministicBasis;
        const std::size_t width=targets.empty()?domain->ambientDimension:targets.front().size();
        for (auto& row:extended) row.resize(width);
        auto result=zhangExactProjectedKernelColumns(std::move(extended),targets,width);
        zhangExactLibraryBoundary("COLUMN_PLAN_LEGACY_RESULT");
        if(diagnostics) { *diagnostics={};diagnostics->legacyOverride=true;
            diagnostics->targetComplete=true; }
        return result;
    }
    struct ProjectionEntry
    {
        ZhangExactMatrix targets;
        ZhangExactColumnProjection replay;
        std::size_t targetBytes=0;
        explicit ProjectionEntry(const ZhangExactMatrix& t):targets(t),replay(t)
        {
            for(const auto& row:t) for(const auto& x:row)
            {
                targetBytes+=sizeof(x);
                if(x!=0) {ZhangExactInteger magnitude=x<0?-x:x;
                    targetBytes+=(boost::multiprecision::msb(magnitude)+8)/8;}
            }
        }
        std::size_t residentBytes() const {return targetBytes+replay.residentBytes();}
    };
    struct Entry
    {
        std::weak_ptr<const ZhangExactAffineIntegerQuotient> owner;
        std::timed_mutex mutex;
        std::unique_ptr<ZhangExactColumnPlan> plan;
        std::vector<std::shared_ptr<ProjectionEntry>> projections;
        bool resourceLimited=false;
    };
    static std::mutex cacheMutex;
    static std::vector<std::shared_ptr<Entry>> cache;
    std::shared_ptr<Entry> entry;
    bool hit=false;
    {
        std::lock_guard<std::mutex> lock(cacheMutex);
        cache.erase(std::remove_if(cache.begin(),cache.end(),
            [](const auto& item){return item->owner.expired();}),cache.end());
        for (const auto& item:cache)
            if (!item->owner.owner_before(domain) && !domain.owner_before(item->owner))
            {
                const auto existing=item->owner.lock();
                if (existing.get()==domain.get()) {entry=item;hit=true;break;}
            }
        if (!entry)
        {
            if (cache.size()>=2) cache.erase(cache.begin());
            entry=std::make_shared<Entry>();entry->owner=domain;cache.push_back(entry);
        }
    }
    zhangExactWorkCheckpoint("COLUMN_PLAN_OWNER_LOCK_ENTRY",true);
    std::unique_lock<std::timed_mutex> entryLock(entry->mutex,std::defer_lock);
    auto* workBudget=zhangActiveExactWorkBudget;
    if(workBudget && workBudget->deadline!=ZhangExactWorkBudget::Clock::time_point::max())
    {
        if(!entryLock.try_lock_until(workBudget->deadline))
        {
            workBudget->phase="COLUMN_PLAN_OWNER_LOCK";
            workBudget->exhausted=true;workBudget->reason="WALL_TIMEOUT";
            throw ZhangExactWorkBudgetExceeded{};
        }
    }
    else entryLock.lock();
    zhangExactWorkCheckpoint("COLUMN_PLAN_OWNER_LOCK_ACQUIRED",true);
    ZhangExactColumnPlanDiagnostics localDiagnostics;
    auto& diagnostic=diagnostics?*diagnostics:localDiagnostics;
    diagnostic={};diagnostic.sharedPlanHit=hit;
    struct Sample
    {
        Entry& entry;
        const ZhangExactColumnPlanDiagnostics& diagnostic;
        std::ostream* trace;
        std::size_t rows,columns,targetRows,targetColumns;
        ~Sample()
        {
            trimTargetCache();
            if(!trace) return;
            const auto* plan=entry.plan.get();
            *trace<<"\nZHANG_BRIDGE_COLUMN_PLAN_SAMPLE source_rows="<<rows
                <<" source_columns="<<columns
                <<" source_metrics_available="<<static_cast<bool>(plan)
                <<" source_nonzeros="<<(plan?plan->sourceNonzeros:0)
                <<" source_max_bits="<<(plan?plan->sourceMaxBits:0)
                <<" target_rows="<<targetRows<<" target_columns="<<targetColumns
                <<" plan_hit="<<diagnostic.sharedPlanHit
                <<" replay_hit="<<diagnostic.targetReplayHit
                <<" plan_complete="<<(plan && plan->complete())
                <<" target_complete="<<diagnostic.targetComplete
                <<" advance_calls="<<(plan?plan->advanceCalls:0)
                <<" resumes="<<(plan && plan->advanceCalls?plan->advanceCalls-1:0)
                <<" rank="<<(plan?plan->rank:0)
                <<" column_operations="<<(plan?plan->operations.size():0)
                <<" journal_bytes="<<(plan?plan->journalBytes:0)
                <<" resident_matrix_bytes="<<(plan?plan->residentMatrixBytes():0)
                <<" total_resident_matrix_bytes="<<ZhangExactColumnPlan::totalResidentMatrixBytes()
                <<" resource_fallback="<<diagnostic.legacyResourceFallback;
        }
        void trimTargetCache()
        {
            auto bytes=[&](){std::size_t sum=0;
                for(const auto& projection:entry.projections) sum+=projection->residentBytes();
                return sum;};
            while(!entry.projections.empty() && bytes()>64*1024*1024)
                entry.projections.erase(entry.projections.begin());
        }
    } sample{*entry,diagnostic,ZhangPhaseTimer::currentTrace(),
        domain->deterministicBasis.size(),static_cast<std::size_t>(domain->ambientDimension),
        targets.size(),targets.empty()?0:targets.front().size()};
    try
    {
        if (!entry->resourceLimited)
        {
            {
                ZhangPhaseTimer timer("BRIDGE_SOURCE_COLUMN_PLAN");
                if (!entry->plan)
                {
                    zhangExactLibraryBoundary("COLUMN_PLAN_SOURCE_COPY_ENTRY");
                    entry->plan=std::make_unique<ZhangExactColumnPlan>(
                        domain->deterministicBasis,domain->ambientDimension);
                    zhangExactLibraryBoundary("COLUMN_PLAN_SOURCE_COPY_EXIT");
                }
                // The object is installed before the first budgeted advance:
                // cancellation leaves the completed operations reusable.
                if (!entry->plan->complete()) entry->plan->advance();
            }
            std::shared_ptr<ProjectionEntry> projection;
            for (const auto& previous:entry->projections)
                if (previous->targets==targets) {projection=previous;break;}
            diagnostic.targetReplayHit=static_cast<bool>(projection);
            if (!projection)
            {
                projection=std::make_shared<ProjectionEntry>(targets);
                std::size_t cells=0;
                for (const auto& row:targets) cells+=row.size();
                // Large targets remain valid but are replayed transiently.
                if (cells<=131072 && projection->residentBytes()<=64*1024*1024)
                {
                    if (entry->projections.size()>=8) entry->projections.erase(entry->projections.begin());
                    entry->projections.push_back(projection);
                }
            }
            {
                ZhangPhaseTimer timer("BRIDGE_TARGET_COLUMN_REPLAY");
                projection->replay.advance(*entry->plan);
            }
            diagnostic.buildCalls=entry->plan->advanceCalls;
            diagnostic.columnOperations=entry->plan->operations.size();
            diagnostic.journalBytes=entry->plan->journalBytes;
            diagnostic.targetComplete=projection->replay.complete();
            zhangExactLibraryBoundary("COLUMN_PLAN_RESULT_COPY_ENTRY");
            auto result=projection->replay.image;
            zhangExactLibraryBoundary("COLUMN_PLAN_RESULT_COPY_EXIT");
            return result;
        }
    }
    catch (const ZhangExactColumnPlanLimitExceeded&)
    {
        // Resource caps change only the implementation, never authorization.
        // Discard the limited plan and use the old bounded-by-caller path.
        entry->plan.reset();entry->projections.clear();entry->resourceLimited=true;
    }
    diagnostic.legacyResourceFallback=true;
    auto extended=domain->deterministicBasis;
    const std::size_t width=targets.empty()?domain->ambientDimension:targets.front().size();
    for (auto& row:extended) row.resize(width);
    auto result=zhangExactProjectedKernelColumns(std::move(extended),targets,width);
    zhangExactLibraryBoundary("COLUMN_PLAN_RESOURCE_FALLBACK_RESULT");
    diagnostic.targetComplete=true;
    return result;
}

struct ZhangR47ProductSearchFrame
{
    bool valid=false;
    int targetRank=0, conditionerRank=0, searchRank=0;
    std::vector<int> columns;
    std::shared_ptr<const ZhangExactAffineIntegerQuotient> affine;
    ZhangExactMatrix projector, imageGenerators;
    ZhangExactVector offsets;
    ZhangExactColumnPlanDiagnostics columnPlan;
    std::string reason="NOT_EVALUATED";
};

/** Parameterise the product image on the feasible integer affine set.
 * Only columns touched by the products/conditioners enter exact elimination;
 * ILS dimension is rank(T K), never the number of network cycles. HNF column
 * image generators retain all Smith divisibility without treating 2Z as Z. */
enum class ZhangProductFrameWork { INTEGER_PROJECTOR, IMAGE_GENERATORS_ONLY };

inline ZhangR47ProductSearchFrame zhangR47CompileProductSearchFrame(
    const ZhangExactMatrix& targets,const ZhangExactMatrix& conditioners,
    const ZhangExactVector& values,int dimension,
    ZhangProductFrameWork work=ZhangProductFrameWork::INTEGER_PROJECTOR)
{
    ZhangPhaseTimer frameTimer("PRODUCT_FRAME_TOTAL");
    zhangCaptureProductFrame(targets,conditioners,values,dimension,static_cast<int>(work));
    ZhangR48ExactTimer productTimer("R48_PRODUCT_IMAGE",targets,dimension);
    ZhangR47ProductSearchFrame out;
    if(targets.empty() || conditioners.size()!=values.size() ||
       !zhangExactRectangularMatrix(targets,dimension) ||
       !zhangExactRectangularMatrix(conditioners,dimension)) return out;
    for(int c=0;c<dimension;++c)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        bool used=false;
        for(const auto& row:targets) used |= row[c]!=0;
        for(const auto& row:conditioners) used |= row[c]!=0;
        if(used) out.columns.push_back(c);
    }
    if(out.columns.empty()) {
        // No physical column remains.  Zero rows are feasible only when every
        // right-hand side is zero; do not let this fast path certify 0 = b.
        for(const auto& value:values) if(value!=0) {
            out.reason="AFFINE_QUOTIENT_INTEGER_CONFLICT";return out;
        }
        out.affine=std::make_shared<const ZhangExactAffineIntegerQuotient>();
        out.valid=true;out.reason="ZERO_PRODUCT_IMAGE";return out;
    }
    auto compact=[&](const ZhangExactMatrix& input)
    {
        ZhangExactMatrix result(input.size(),ZhangExactVector(out.columns.size()));
        for(std::size_t r=0;r<input.size();++r) for(std::size_t c=0;c<out.columns.size();++c) {
            zhangExactWorkCheckpoint("PRODUCT_DOMAIN_COMPACT_COPY");
            result[r][c]=input[r][out.columns[c]];
        }
        return result;
    };
    // The identity column catalogue is already in the exact legacy order.
    // Borrow both inputs rather than duplicating two cpp_int matrices solely
    // to reproduce that same catalogue. A genuinely compressed input remains
    // owned locally, with the original ascending column order unchanged.
    const bool identityColumns=out.columns.size()==static_cast<std::size_t>(dimension);
    ZhangExactMatrix compactTargets,compactConditioners;
    if(!identityColumns)
    {
        {ZhangPhaseTimer timer("PRODUCT_FRAME_TARGET_COMPACT");compactTargets=compact(targets);}
        {ZhangPhaseTimer timer("PRODUCT_FRAME_CONDITIONER_COMPACT");compactConditioners=compact(conditioners);}
    }
    const auto& t=identityColumns?targets:compactTargets;
    const auto& h=identityColumns?conditioners:compactConditioners;
    if(auto* trace=ZhangPhaseTimer::currentTrace())
        *trace << "\nZHANG_PRODUCT_FRAME_COLUMN_VIEW source_columns=" << dimension
            << " compact_columns=" << out.columns.size()
            << " target_rows=" << targets.size()
            << " conditioner_rows=" << conditioners.size()
            << " identity_borrowed=" << identityColumns
            << " compact_copied_cells=" << (identityColumns?0:
                (targets.size()+conditioners.size())*out.columns.size());
    if(work==ZhangProductFrameWork::INTEGER_PROJECTOR)
    {
        ZhangPhaseTimer timer("PRODUCT_FRAME_TARGET_HNF");
        out.targetRank=zhangExactRowHermiteNormalForm(t).basis.size();
    }
    const auto quotientWork=work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY
        ? ZhangExactQuotientWork::PARTICULAR_ONLY
        : ZhangExactQuotientWork::PARTICULAR_AND_KERNEL;
    out.affine=std::make_shared<const ZhangExactAffineIntegerQuotient>(
        zhangExactAffineIntegerQuotient(h,values,out.columns.size(),quotientWork));
    if(!out.affine->valid) {out.reason=out.affine->failureReason;return out;}
    out.conditionerRank=out.affine->deterministicRank;
    if(!out.affine->quotientRank) {out.valid=true;out.reason="PRODUCT_FULLY_DETERMINED";return out;}
    ZhangExactMatrix k,m;
    if(work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
    {
        ZhangR48ExactTimer timer("R48_PRODUCT_IMAGE_PROJECTED_KERNEL",
            out.affine->deterministicBasis,out.columns.size());
        m=zhangExactProjectedKernelColumns(
            out.affine->deterministicBasis,t,out.columns.size());
        if(m.empty() || m.front().size()!=
            static_cast<std::size_t>(out.affine->quotientRank))
        {out.reason="PRODUCT_PROJECTED_KERNEL_FAILED";return out;}
    }
    else
    {
        ZhangPhaseTimer timer("PRODUCT_FRAME_KERNEL_PROJECT");
        k.assign(out.columns.size(),ZhangExactVector(out.affine->quotientRank));
        for(int r=0;r<out.affine->quotientRank;++r)
        for(int c=0;c<out.columns.size();++c)
            k[c][r]=out.affine->kernelBasis[r][c];
        m=zhangExactMultiply(t,k);
    }
    const auto image=[&](){ZhangPhaseTimer timer("PRODUCT_FRAME_IMAGE_HNF");return zhangPrimitiveImageCoordinates(m);}();
    if(!image.valid) {out.reason="PRODUCT_IMAGE_COORDINATES_FAILED";return out;}
    out.imageGenerators=image.imageGenerators;
    out.searchRank=image.primitiveRows.size();
    // R51 constructs and audits its own rational left inverse of this exact
    // image. This path keeps T*V during column operations and never stores
    // the full physical kernel or the unused network integer lift.
    if(work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY) {
        out.valid=true;
        out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
        return out;
    }
    // Solve only the product-image coordinate RHS, never construct the
    // entire network quotient left inverse. K is saturated, so these rows
    // have exact integer lifts. Audit each lift before exposing it to ILS.
    const auto lifts=[&]() {
        ZhangPhaseTimer liftTimer("PRODUCT_FRAME_INTEGER_LIFT");
        ZhangR48ExactTimer timer("R48_PRODUCT_IMAGE_LIFT",k);
        return zhangIntegerRowLatticeContainsBatch(k,image.primitiveRows);
    }();
    ZhangExactMatrix projected;
    for(std::size_t i=0;i<lifts.size();++i) {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        if(!lifts[i].contained || lifts[i].combination.size()!=out.columns.size()) {
            out.reason="PRODUCT_IMAGE_INTEGER_LIFT_FAILED";return out;
        }
        projected.push_back(lifts[i].combination);
    }
    if(!projected.empty() && zhangExactMultiply(projected,k)!=image.primitiveRows) {
        out.reason="PRODUCT_IMAGE_INTEGER_LIFT_AUDIT_FAILED";return out;
    }
    for(const auto& row:projected)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        ZhangExactVector full(dimension); ZhangExactInteger offset=0;
        for(std::size_t c=0;c<out.columns.size();++c)
        {
            zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");full[out.columns[c]]=row[c];offset-=row[c]*out.affine->particularSolution[c];}
        out.projector.push_back(std::move(full));out.offsets.push_back(offset);
    }
    out.valid=true;out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
    return out;
}

inline ZhangR47ProductSearchFrame zhangR49CompileTargetOnDomain(
 const ZhangExactMatrix& targets,const ZhangR47ProductSearchFrame& domain,int dimension)
{
 ZhangR47ProductSearchFrame out;
 if(!domain.valid || !domain.affine || !domain.affine->valid || targets.empty() || !zhangExactRectangularMatrix(targets,dimension))return out;
 std::vector<bool> represented(dimension);for(int c:domain.columns)represented[c]=true;
 for(const auto& row:targets)for(int c=0;c<dimension;++c)if(row[c]!=0 && !represented[c]) {
  out.reason="TARGET_OUTSIDE_CACHED_DOMAIN";return out;
 }
 out.columns=domain.columns;out.affine=domain.affine;out.conditionerRank=domain.conditionerRank;
 ZhangExactMatrix t(targets.size(),ZhangExactVector(out.columns.size()));
 for(int i=0;i<targets.size();++i)for(int c=0;c<out.columns.size();++c)t[i][c]=targets[i][out.columns[c]];
 out.targetRank=zhangExactRowHermiteNormalForm(t).basis.size();
 if(!out.affine->quotientRank){out.valid=true;out.reason="PRODUCT_FULLY_DETERMINED";return out;}
    ZhangExactMatrix k(out.columns.size(),ZhangExactVector(out.affine->quotientRank));
    for(int r=0;r<out.affine->quotientRank;++r) for(int c=0;c<out.columns.size();++c)
        k[c][r]=out.affine->kernelBasis[r][c];
    const auto m=zhangExactMultiply(t,k);
    const auto image=zhangPrimitiveImageCoordinates(m);
    if(!image.valid) {out.reason="PRODUCT_IMAGE_COORDINATES_FAILED";return out;}
    out.imageGenerators=image.imageGenerators;
    out.searchRank=image.primitiveRows.size();
    // Solve only the product-image coordinate RHS, never construct the
    // entire network quotient left inverse. K is saturated, so these rows
    // have exact integer lifts. Audit each lift before exposing it to ILS.
    const auto lifts=[&]() {
        ZhangR48ExactTimer timer("R48_PRODUCT_IMAGE_LIFT",k);
        return zhangIntegerRowLatticeContainsBatch(k,image.primitiveRows);
    }();
    ZhangExactMatrix projected;
    for(std::size_t i=0;i<lifts.size();++i) {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        if(!lifts[i].contained || lifts[i].combination.size()!=out.columns.size()) {
            out.reason="PRODUCT_IMAGE_INTEGER_LIFT_FAILED";return out;
        }
        projected.push_back(lifts[i].combination);
    }
    if(!projected.empty() && zhangExactMultiply(projected,k)!=image.primitiveRows) {
        out.reason="PRODUCT_IMAGE_INTEGER_LIFT_AUDIT_FAILED";return out;
    }
    for(const auto& row:projected)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        ZhangExactVector full(dimension); ZhangExactInteger offset=0;
        for(std::size_t c=0;c<out.columns.size();++c)
        {
            zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");full[out.columns[c]]=row[c];offset-=row[c]*out.affine->particularSolution[c];}
        out.projector.push_back(std::move(full));out.offsets.push_back(offset);
    }
    out.valid=true;out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
    return out;
}

/** Target-only view of an already verified immutable affine domain.
 * This deliberately does not call HNF, a feasibility/rank-membership query,
 * or a particular-solution solver again. The caller must keep the owner and
 * its current epoch/generation identity alive. Appended target-only columns
 * remain FREE integers; their zero particular value is not a zero constraint.
 * The existing projected-kernel operation sequence and canonical image HNF
 * are retained, so product search coordinates do not change. */
inline ZhangR47ProductSearchFrame zhangR47CompileTargetOnVerifiedDomain(
    const ZhangExactMatrix& targets,
    std::shared_ptr<const ZhangExactAffineIntegerQuotient> domain,
    int dimension,
    ZhangProductFrameWork work=ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
{
    ZhangR48ExactTimer timer("R51_PRODUCT_IMAGE_VERIFIED_DOMAIN",targets,dimension);
    ZhangR47ProductSearchFrame out;
    if (!domain || !domain->valid || targets.empty() || dimension <= 0 ||
        domain->ambientDimension < 0 || domain->ambientDimension > dimension ||
        !zhangExactRectangularMatrix(targets,dimension) ||
        !zhangExactRectangularMatrix(domain->deterministicBasis,domain->ambientDimension) ||
        domain->deterministicBasis.size()!=domain->deterministicValues.size() ||
        domain->particularSolution.size()!=static_cast<std::size_t>(domain->ambientDimension) ||
        domain->deterministicRank!=static_cast<int>(domain->deterministicBasis.size()) ||
        domain->quotientRank!=domain->ambientDimension-domain->deterministicRank)
    {out.reason="VERIFIED_DOMAIN_SHAPE_INVALID";return out;}
    const auto sourceDomain=domain;
    // A current compiled domain owns precisely the nonzero source columns.
    // Require that contract, rather than silently recompact a changed basis.
    for (int c=0;c<domain->ambientDimension;++c)
    {
        bool used=false;
        for (const auto& row:domain->deterministicBasis) used |= row[c]!=0;
        if (!used) {out.reason="VERIFIED_DOMAIN_SOURCE_COLUMN_UNUSED";return out;}
        out.columns.push_back(c);
    }
    for (int c=domain->ambientDimension;c<dimension;++c)
    {
        bool used=false;
        for (const auto& row:targets) used |= row[c]!=0;
        if (used) out.columns.push_back(c);
    }
    if (out.columns.empty())
    {out.affine=std::move(domain);out.valid=true;out.reason="ZERO_PRODUCT_IMAGE";return out;}
    const bool needKernel=work==ZhangProductFrameWork::INTEGER_PROJECTOR &&
        domain->kernelBasis.size()!=static_cast<std::size_t>(domain->quotientRank);
    if (out.columns.size()==static_cast<std::size_t>(domain->ambientDimension) && !needKernel)
        out.affine=std::move(domain);
    else
    {
        auto extended=std::make_shared<ZhangExactAffineIntegerQuotient>(*domain);
        extended->ambientDimension=out.columns.size();
        extended->quotientRank=extended->ambientDimension-extended->deterministicRank;
        for (auto& row:extended->deterministicBasis) row.resize(out.columns.size());
        extended->particularSolution.resize(out.columns.size());
        extended->kernelBasis.clear();extended->quotientProjector.clear();
        if (work==ZhangProductFrameWork::INTEGER_PROJECTOR)
            extended->kernelBasis=zhangExactIntegerKernel(
                extended->deterministicBasis,out.columns.size());
        out.affine=std::move(extended);
    }
    const bool identityColumns=out.columns.size()==static_cast<std::size_t>(dimension);
    ZhangExactMatrix compactTargets;
    if(!identityColumns)
    {
        compactTargets.assign(targets.size(),ZhangExactVector(out.columns.size()));
        for (std::size_t r=0;r<compactTargets.size();++r)
        for (std::size_t c=0;c<out.columns.size();++c)
        {
            zhangExactWorkCheckpoint("VERIFIED_DOMAIN_TARGET_COPY");
            compactTargets[r][c]=targets[r][out.columns[c]];
        }
    }
    const auto& t=identityColumns?targets:compactTargets;
    if(auto* trace=ZhangPhaseTimer::currentTrace())
        *trace << "\nZHANG_VERIFIED_PRODUCT_TARGET_COLUMN_VIEW source_columns=" << dimension
            << " compact_columns=" << out.columns.size()
            << " target_rows=" << targets.size()
            << " identity_borrowed=" << identityColumns
            << " compact_copied_cells=" << (identityColumns?0:targets.size()*out.columns.size());
    out.conditionerRank=out.affine->deterministicRank;
    if (work==ZhangProductFrameWork::INTEGER_PROJECTOR)
        out.targetRank=zhangExactRowHermiteNormalForm(t).basis.size();
    if (!out.affine->quotientRank)
    {out.valid=true;out.reason="PRODUCT_FULLY_DETERMINED";return out;}
    ZhangExactMatrix k,m;
    if (work==ZhangProductFrameWork::IMAGE_GENERATORS_ONLY)
    {
        ZhangR48ExactTimer projectionTimer("R48_PRODUCT_IMAGE_PROJECTED_KERNEL",
            out.affine->deterministicBasis,out.columns.size());
        m=zhangProjectedKernelOnOwnedDomain(sourceDomain,t,&out.columnPlan);
        if (m.empty() || m.front().size()!=static_cast<std::size_t>(out.affine->quotientRank))
        {out.reason="PRODUCT_PROJECTED_KERNEL_FAILED";return out;}
    }
    else
    {
        k.assign(out.columns.size(),ZhangExactVector(out.affine->quotientRank));
        for (int r=0;r<out.affine->quotientRank;++r)
        for (std::size_t c=0;c<out.columns.size();++c) k[c][r]=out.affine->kernelBasis[r][c];
        m=zhangExactMultiply(t,k);
    }
    const auto image=zhangPrimitiveImageCoordinates(m);
    if (!image.valid) {out.reason="PRODUCT_IMAGE_COORDINATES_FAILED";return out;}
    out.imageGenerators=image.imageGenerators;
    out.searchRank=image.primitiveRows.size();
    if (work==ZhangProductFrameWork::INTEGER_PROJECTOR)
    {
        const auto lifts=zhangIntegerRowLatticeContainsBatch(k,image.primitiveRows);
        ZhangExactMatrix projected;
        for (const auto& lift:lifts)
        {
            if (!lift.contained || lift.combination.size()!=out.columns.size())
            {out.reason="PRODUCT_IMAGE_INTEGER_LIFT_FAILED";return out;}
            projected.push_back(lift.combination);
        }
        if (!projected.empty() && zhangExactMultiply(projected,k)!=image.primitiveRows)
        {out.reason="PRODUCT_IMAGE_INTEGER_LIFT_AUDIT_FAILED";return out;}
        for (const auto& row:projected)
        {
            ZhangExactVector full(dimension);ZhangExactInteger offset=0;
            for (std::size_t c=0;c<out.columns.size();++c)
            {full[out.columns[c]]=row[c];offset-=row[c]*out.affine->particularSolution[c];}
            out.projector.push_back(std::move(full));out.offsets.push_back(std::move(offset));
        }
    }
    out.valid=true;out.reason=out.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
    return out;
}

inline bool zhangR47ProductConsequence(const ZhangR47ProductSearchFrame& frame,
    const ZhangExactVector& target,const ZhangExactInteger& offset,ZhangExactInteger& value)
{
    if(!frame.valid || !frame.affine || !frame.affine->valid) return false;
    // A target-only image view deliberately omits the full physical kernel.
    // An empty, UNREQUESTED kernel is not evidence that all directions are
    // fixed. Callers needing consequences must use the compiled HNF query or
    // a genuinely complete-kernel frame.
    if(frame.affine->quotientRank<0 ||
       frame.affine->kernelBasis.size()!=static_cast<std::size_t>(frame.affine->quotientRank) ||
       frame.affine->particularSolution.size()!=frame.columns.size() ||
       !zhangExactRectangularMatrix(frame.affine->kernelBasis,frame.columns.size())) return false;
    std::vector<bool> represented(target.size());
    ZhangExactVector compact(frame.columns.size());
    for(std::size_t c=0;c<frame.columns.size();++c)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        if(frame.columns[c]<0 || frame.columns[c]>=target.size()) return false;
        compact[c]=target[frame.columns[c]];represented[frame.columns[c]]=true;
    }
    for(std::size_t c=0;c<target.size();++c) if(target[c]!=0 && !represented[c]) return false;
    for(const auto& kernel:frame.affine->kernelBasis)
    {
        zhangExactWorkCheckpoint("zhangR47ProductDomain_LOOP");
        ZhangExactInteger residual=0;
        for(std::size_t c=0;c<compact.size();++c) residual+=compact[c]*kernel[c];
        if(residual!=0) return false;
    }
    value=offset;
    for(std::size_t c=0;c<compact.size();++c) value+=compact[c]*frame.affine->particularSolution[c];
    return true;
}
