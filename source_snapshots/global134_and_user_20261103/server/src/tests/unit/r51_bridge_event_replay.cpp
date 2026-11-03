#include "common/zhangBridgeSupportAnalysis.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace fs=std::filesystem;

static void require(bool condition,const std::string& message)
{if(!condition) throw std::runtime_error(message);}

static void setLegacy(bool legacy)
{
#ifdef _WIN32
    _putenv_s("ZHANG_R51_BRIDGE_COLUMN_PLAN",legacy?"legacy":"");
#else
    if(legacy) ::setenv("ZHANG_R51_BRIDGE_COLUMN_PLAN","legacy",1);
    else ::unsetenv("ZHANG_R51_BRIDGE_COLUMN_PLAN");
#endif
}

struct Capture
{
    std::size_t dimension=0;
    std::uint64_t generation=0;
    std::string epoch,satellite,root;
    std::vector<std::string> ids;
    std::set<std::size_t> unavailable;
    ZhangExactMatrix targets,rows;
    ZhangExactVector rhs;
};

static Capture readCaptureStream(std::istream& input)
{
    require(bool(input),"cannot read capture stream");
    Capture out;std::string magic;
    input>>magic>>out.dimension>>out.generation;
    require(magic=="R51_BRIDGE_SUPPORT_V1" && out.dimension>0 && out.dimension<=20000,
        "unsupported bridge capture header");
    input>>std::quoted(out.epoch)>>std::quoted(out.satellite)>>std::quoted(out.root);
    out.ids.resize(out.dimension);
    for(auto& id:out.ids) input>>std::quoted(id);
    std::size_t count=0;input>>count;
    require(count<=out.dimension,"too many unavailable columns");
    for(std::size_t n=0;n<count;++n)
    {
        std::size_t column=0;input>>column;
        require(bool(input) && column<out.dimension && out.unavailable.insert(column).second,
            "invalid/duplicate unavailable column");
    }
    auto matrix=[&]() {
        std::size_t rowCount=0;input>>rowCount;
        require(bool(input) && rowCount<=32000000/out.dimension,
            "capture matrix exceeds driver allocation guard");
        ZhangExactMatrix rows(rowCount,ZhangExactVector(out.dimension));
        for(auto& row:rows)
        {
            std::size_t nonzeros=0;input>>nonzeros;
            require(bool(input) && nonzeros<=out.dimension,"invalid sparse row width");
            std::set<std::size_t> seen;
            for(std::size_t n=0;n<nonzeros;++n)
            {
                std::size_t column=0;ZhangExactInteger value;input>>column>>value;
                require(bool(input) && column<out.dimension && seen.insert(column).second && value!=0,
                    "invalid sparse matrix entry");
                row[column]=std::move(value);
            }
        }
        return rows;
    };
    out.targets=matrix();out.rows=matrix();
    input>>count;require(bool(input) && count==out.rows.size(),"invalid RHS length");
    out.rhs.resize(count);for(auto& value:out.rhs) input>>value;
    require(bool(input) && out.rhs.size()==out.rows.size() && !out.targets.empty(),
        "truncated or inconsistent bridge capture");
    input>>std::ws;require(input.eof(),"unexpected trailing bridge capture data");
    return out;
}

static Capture readCapture(const fs::path& path)
{
    std::ifstream input(path);
    require(bool(input),"cannot read capture: "+path.string());
    return readCaptureStream(input);
}

struct CachedDomain
{
    std::uint64_t generation=0;
    std::string epoch,root;
    std::vector<std::string> sourceIds;
    std::shared_ptr<const ZhangExactAffineIntegerQuotient> domain;
};

static std::shared_ptr<const ZhangExactAffineIntegerQuotient> verifiedDomain(
    const Capture& capture,std::vector<CachedDomain>& cache,bool& reused)
{
    std::vector<bool> used(capture.dimension);
    for(const auto& row:capture.rows)
        for(std::size_t c=0;c<capture.dimension;++c) used[c]=used[c] || row[c]!=0;
    std::size_t sourceDimension=0;
    while(sourceDimension<capture.dimension && used[sourceDimension]) ++sourceDimension;
    for(std::size_t c=sourceDimension;c<capture.dimension;++c)
        require(!used[c],"capture violates runtime source-prefix/free-target-column contract");
    ZhangExactMatrix rows(capture.rows.size(),ZhangExactVector(sourceDimension));
    for(std::size_t r=0;r<rows.size();++r)
        std::copy_n(capture.rows[r].begin(),sourceDimension,rows[r].begin());
    std::vector<std::string> ids(capture.ids.begin(),capture.ids.begin()+sourceDimension);
    for(const auto& previous:cache)
        if(previous.generation==capture.generation && previous.epoch==capture.epoch &&
            previous.root==capture.root && previous.sourceIds==ids &&
            previous.domain->deterministicBasis==rows && previous.domain->deterministicValues==capture.rhs)
        {reused=true;return previous.domain;}
    reused=false;
    std::shared_ptr<const ZhangExactAffineIntegerQuotient> domain;
    if(sourceDimension)
        domain=std::make_shared<const ZhangExactAffineIntegerQuotient>(
            zhangExactAffineIntegerQuotient(rows,capture.rhs,sourceDimension,
                ZhangExactQuotientWork::PARTICULAR_ONLY));
    else
    {
        require(rows.empty() && capture.rhs.empty(),"empty source has noncanonical constraints");
        auto freeDomain=std::make_shared<ZhangExactAffineIntegerQuotient>();
        freeDomain->valid=true;freeDomain->failureReason="NONE";domain=freeDomain;
    }
    require(domain->valid,"capture has no certified feasible integer source domain: "+domain->failureReason);
    // Runtime captures the compiled HNF, not arbitrary insertion-order input.
    require(domain->deterministicBasis==rows && domain->deterministicValues==capture.rhs,
        "captured source is not the canonical runtime HNF/RHS");
    if(cache.size()>=2) cache.erase(cache.begin());
    cache.push_back({capture.generation,capture.epoch,capture.root,std::move(ids),domain});
    return domain;
}

static void replay(const fs::path& path,const Capture& capture,std::vector<CachedDomain>& cache)
{
    using Clock=std::chrono::steady_clock;
    const auto begin=Clock::now();bool domainReused=false;
    const auto domain=verifiedDomain(capture,cache,domainReused);
    const auto ready=Clock::now();
    const auto support=zhangAnalyzeBridgeSupport(capture.rows,capture.targets,
        capture.unavailable,capture.dimension);
    require(support.valid,"support analysis failed its exact witness audit");
    // Independently audit each negative witness against FULL H and targets.
    // This does not certify covariance, ratio acceptance or real-time latency.
    for(const auto& witness:support.counterexamples)
    {
        require(witness.size()==capture.dimension,"negative witness shape mismatch");
        for(std::size_t c=0;c<capture.dimension;++c)
            require(capture.unavailable.contains(c) || witness[c]==0,"negative witness changes a live arc");
        for(const auto& row:capture.rows)
        {
            ZhangExactInteger value=0;
            for(std::size_t c=0;c<capture.dimension;++c) value+=row[c]*witness[c];
            require(value==0,"negative witness does not satisfy H*v=0");
        }
        bool changes=false;
        for(const auto& target:capture.targets)
        {
            ZhangExactInteger value=0;
            for(std::size_t c=0;c<capture.dimension;++c) value+=target[c]*witness[c];
            changes=changes || value!=0;
        }
        require(changes,"negative witness does not change any target");
    }
    setLegacy(true);const auto oldBegin=Clock::now();
    const auto oldFrame=zhangR47CompileTargetOnVerifiedDomain(capture.targets,domain,capture.dimension);
    const auto oldImage=zhangR51PhysicalImageFromFrame(capture.targets,capture.targets,oldFrame);
    const auto oldEnd=Clock::now();setLegacy(false);
    const auto newBegin=Clock::now();
    const auto newFrame=zhangR47CompileTargetOnVerifiedDomain(capture.targets,domain,capture.dimension);
    const auto newImage=zhangR51PhysicalImageFromFrame(capture.targets,capture.targets,newFrame);
    const auto newEnd=Clock::now();
    require(oldFrame.valid && newFrame.valid && oldImage.valid && newImage.valid,
        "old/new image compilation failed");
    require(oldFrame.columns==newFrame.columns && oldFrame.searchRank==newFrame.searchRank &&
        oldFrame.projector==newFrame.projector && oldFrame.offsets==newFrame.offsets &&
        oldFrame.imageGenerators==newFrame.imageGenerators &&
        oldImage.particularTarget==newImage.particularTarget && oldImage.generators==newImage.generators &&
        oldImage.projector==newImage.projector && oldImage.offsets==newImage.offsets,
        "old/new exact bridge coordinates/origin/projector/offset differ");
    const auto millis=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
    std::cout<<"{\"file\":"<<std::quoted(path.string())<<",\"epoch\":"<<std::quoted(capture.epoch)
        <<",\"satellite\":"<<std::quoted(capture.satellite)<<",\"generation\":"<<capture.generation
        <<",\"rows\":"<<capture.rows.size()<<",\"dimension\":"<<capture.dimension
        <<",\"unavailable\":"<<capture.unavailable.size()<<",\"missing_rank\":"<<support.missingRank
        <<",\"negative_witnesses\":"<<support.counterexamples.size()
        <<",\"representable\":"<<(support.representable?"true":"false")
        <<",\"domain_reused\":"<<(domainReused?"true":"false")
        <<",\"domain_prepare_ms\":"<<millis(begin,ready)<<",\"legacy_image_ms\":"<<millis(oldBegin,oldEnd)
        <<",\"new_image_ms\":"<<millis(newBegin,newEnd)<<",\"plan_hit\":"<<(newFrame.columnPlan.sharedPlanHit?"true":"false")
        <<",\"plan_advance_calls\":"<<newFrame.columnPlan.buildCalls
        <<",\"resource_fallback\":"<<(newFrame.columnPlan.legacyResourceFallback?"true":"false")
        <<",\"exact_equal\":true,\"scope\":\"exact_affine_image_only_no_covariance_or_AR_search\"}\n";
}

// Compatibility experiment for EXISTING 02:00 ProductFrame inputs, not a
// captured bridge event. Preserve original column order and compact exactly
// as the generic factory did; free source columns remain free integers.
static void replayFrozenAlgebra(const fs::path& path,const Capture& capture)
{
    using Clock=std::chrono::steady_clock;
    std::vector<int> columns;
    for(std::size_t c=0;c<capture.dimension;++c)
    {
        bool used=false;
        for(const auto& row:capture.rows) used=used || row[c]!=0;
        for(const auto& row:capture.targets) used=used || row[c]!=0;
        if(used) columns.push_back(c);
    }
    require(!columns.empty(),"frozen zero-dimensional product frame");
    const auto compact=[&](const ZhangExactMatrix& original) {
        ZhangExactMatrix result(original.size(),ZhangExactVector(columns.size()));
        for(std::size_t r=0;r<original.size();++r)
        for(std::size_t c=0;c<columns.size();++c) result[r][c]=original[r][columns[c]];
        return result;
    };
    const auto target=compact(capture.targets),source=compact(capture.rows);
    const auto prepareBegin=Clock::now();
    auto domain=std::make_shared<const ZhangExactAffineIntegerQuotient>(
        zhangExactAffineIntegerQuotient(source,capture.rhs,columns.size(),
            ZhangExactQuotientWork::PARTICULAR_ONLY));
    require(domain->valid,"frozen frame has infeasible affine source");
    const auto prepared=Clock::now();
    auto makeFrame=[&](const ZhangExactMatrix& projected) {
        ZhangR47ProductSearchFrame frame;
        const auto image=zhangPrimitiveImageCoordinates(projected);
        require(image.valid,"frozen projected integer image HNF failed");
        frame.columns=columns;frame.affine=domain;
        frame.conditionerRank=domain->deterministicRank;
        frame.imageGenerators=image.imageGenerators;frame.searchRank=image.primitiveRows.size();
        frame.valid=true;frame.reason=frame.searchRank?"EXACT_PRODUCT_IMAGE_QUOTIENT":"PRODUCT_FULLY_DETERMINED";
        return frame;
    };
    const auto oldBegin=Clock::now();
    auto oldFrame=makeFrame(zhangExactProjectedKernelColumns(domain->deterministicBasis,target,columns.size()));
    const auto oldEnd=Clock::now();
    ZhangExactMatrix rawNew;
    ZhangExactColumnPlanDiagnostics statistics;
    const auto newBegin=Clock::now();
    // Large generic search frames contain over 1000 targets, unlike a bridge's
    // two targets. Batch target replay without changing row/order or the final
    // image HNF. This is test-harness adaptation, not a production capability.
    for(std::size_t first=0;first<target.size();first+=32)
    {
        const auto last=std::min(target.size(),first+32);
        ZhangExactMatrix batch(target.begin()+first,target.begin()+last);
        auto image=zhangProjectedKernelOnOwnedDomain(domain,batch,&statistics);
        rawNew.insert(rawNew.end(),std::make_move_iterator(image.begin()),std::make_move_iterator(image.end()));
    }
    auto newFrame=makeFrame(rawNew);
    const auto newEnd=Clock::now();
    require(oldFrame.imageGenerators==newFrame.imageGenerators && oldFrame.searchRank==newFrame.searchRank,
        "frozen integer image or search coordinate ordering changed");
    const auto oldImage=zhangR51PhysicalImageFromFrame(capture.targets,capture.targets,oldFrame);
    const auto newImage=zhangR51PhysicalImageFromFrame(capture.targets,capture.targets,newFrame);
    require(oldImage.valid && newImage.valid && oldImage.generators==newImage.generators &&
        oldImage.particularTarget==newImage.particularTarget && oldImage.projector==newImage.projector &&
        oldImage.offsets==newImage.offsets,"frozen full physical image differs");
    const auto millis=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
    std::cout<<"{\"file\":"<<std::quoted(path.string())<<",\"epoch\":"<<std::quoted(capture.epoch)
        <<",\"input_dimension\":"<<capture.dimension<<",\"compact_dimension\":"<<columns.size()
        <<",\"source_rows\":"<<source.size()<<",\"target_rows\":"<<target.size()
        <<",\"source_prepare_ms\":"<<millis(prepareBegin,prepared)
        <<",\"legacy_projection_image_ms\":"<<millis(oldBegin,oldEnd)
        <<",\"new_batched_projection_image_ms\":"<<millis(newBegin,newEnd)
        <<",\"plan_advance_calls\":"<<statistics.buildCalls<<",\"plan_operations\":"<<statistics.columnOperations
        <<",\"resource_fallback\":"<<(statistics.legacyResourceFallback?"true":"false")
        <<",\"exact_equal\":true,\"scope\":\"generic_real0200_inputs_exact_backend_consistency_not_bridge_event_no_covariance_or_AR_search\"}\n";
}

int main(int argc,char** argv)
{
    try
    {
        std::vector<CachedDomain> cache;
        if(argc==1 || std::string(argv[1])=="--self-test")
        {
            std::istringstream input(
                "R51_BRIDGE_SUPPORT_V1 3 1\n"
                "\"2024-07-17 02:37:30\" \"G10\" \"synthetic_float_root\"\n"
                "\"old\"\n\"live\"\n\"unknown\"\n"
                "1 0\n2\n2 1 1 2 3\n2 0 2 2 -1\n1\n2 0 2 1 1\n1 1\n");
            auto fixture=readCaptureStream(input);
            replay("SYNTHETIC_REPRESENTABLE",fixture,cache);
            fixture.unavailable={0,2};
            replay("SYNTHETIC_MISSING_DIRECTION",fixture,cache);
            return 0;
        }
        std::vector<fs::path> paths;
        const bool frozen=std::string(argv[1])=="--frozen-algebra";
        for(int i=frozen?2:1;i<argc;++i)
        {
            const fs::path path(argv[i]);
            if(fs::is_directory(path))
            {
                for(const auto& item:fs::directory_iterator(path))
                    if(item.is_regular_file() && item.path().extension()==".txt") paths.push_back(item.path());
            }
            else paths.push_back(path);
        }
        std::sort(paths.begin(),paths.end());
        require(!paths.empty(),"no bridge capture files selected");
        setLegacy(false);
        for(const auto& path:paths)
        {
            const auto capture=readCapture(path);
            if(frozen) replayFrozenAlgebra(path,capture);
            else replay(path,capture,cache);
        }
        setLegacy(false);
        return 0;
    }
    catch(const std::exception& error)
    {
        setLegacy(false);std::cerr<<"bridge replay failed: "<<error.what()<<"\n";return 1;
    }
}
