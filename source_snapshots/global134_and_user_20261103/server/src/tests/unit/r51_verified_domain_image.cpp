// Old/new coordinate-equivalence tests plus independent integer oracles.
// These are mathematical regressions, not a real-epoch GNSS replay.
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangCurrentBridgeRepresentation.hpp"

#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
void require(bool value,const char* reason)
{
    if (!value) throw std::runtime_error(reason);
}

constexpr long epoch=1405216800;
ZhangDecisionProofPtr proof()
{
    auto out=std::make_shared<ZhangIntegerDecisionProof>();
    out->id="verified-domain-test";
    out->originalStatement="current exact fixture";
    out->observationProvenance="runtime|PPP_COMMITTED|E1405216800";
    out->conditionalFailureBound=1e-6;
    return out;
}

ZhangCurrentEpochEvidence evidenceFor(const ZhangExactMatrix& rows,const ZhangExactVector& rhs)
{
    const std::size_t dimension=rows.empty()?3:rows.front().size();
    ZhangProductPhysicalCycleChart chart;
    chart.columns=dimension;
    for (std::size_t c=0;c<dimension;++c)
    {
        const auto id="physical-"+std::to_string(c);
        chart.chordIdentities[c]=id;
        chart.expansions[c][id]=1;
        chart.activePhysicalIdentities.insert(id);
    }
    ZhangCurrentEpochEvidence out;
    require(out.bind("runtime",E_Sys::GPS,epoch,"posterior","chart",chart,"root"),"fixture bind");
    const auto decision=proof();
    for (std::size_t r=0;r<rows.size();++r)
    {
        ZhangCurrentEpochPhysicalRow row;
        for (std::size_t c=0;c<dimension;++c) if (rows[r][c]!=0)
            row["physical-"+std::to_string(c)]=rows[r][c];
        if (!row.empty()) require(out.add(row,rhs[r],{decision}),"fixture add");
    }
    return out;
}

void compareImages(const ZhangR51PhysicalImage& old,const ZhangR51PhysicalImage& fresh)
{
    require(old.valid && fresh.valid,"old/new image invalid");
    require(old.reason==fresh.reason,"image reason changed");
    require(old.particularTarget==fresh.particularTarget,"particular target changed");
    require(old.generators==fresh.generators,"integer image generators changed");
    require(old.projector==fresh.projector,"rational projector/search coordinates changed");
    require(old.offsets==fresh.offsets,"rational offsets/search mean changed");
}

void compareFrames(const ZhangR47ProductSearchFrame& old,const ZhangR47ProductSearchFrame& fresh)
{
    require(old.valid && fresh.valid,"old/new frame invalid");
    require(old.columns==fresh.columns && old.targetRank==fresh.targetRank &&
        old.conditionerRank==fresh.conditionerRank && old.searchRank==fresh.searchRank,
        "frame ranks or column order changed");
    require(old.projector==fresh.projector && old.offsets==fresh.offsets &&
        old.imageGenerators==fresh.imageGenerators,"frame search coordinates changed");
    require(old.affine->particularSolution==fresh.affine->particularSolution &&
        old.affine->deterministicBasis==fresh.affine->deterministicBasis &&
        old.affine->deterministicValues==fresh.affine->deterministicValues,
        "frame original affine coordinates changed");
}

ZhangR51ResidualImageAdmission admit(const ZhangR51PhysicalImage& image,
    const ZhangExactMatrix& targets,const ZhangExactMatrix& basis,const ZhangExactVector& rhs,
    const ZhangExactMatrix& z,const ZhangExactVector& integers)
{
    ZhangExactMatrix lifted;ZhangExactVector values;
    require(zhangR51LiftRationalRows(z,integers,image.projector,image.offsets,lifted,values),"lift fixture");
    return zhangR51CertifyResidualImageAdmission(image,targets,basis,rhs,z,integers,lifted,values);
}
}

int main()
{
    std::mt19937 random(199);
    std::size_t witnesses=0,images=0,admissions=0;
    // The compact tape must reproduce the legacy particular witness, not
    // merely another solution of the same integer system. Known integer
    // generators independently establish feasibility; every witness is audited.
    for (int trial=0;trial<500;++trial)
    {
        const int count=1+random()%5,dimension=1+random()%4;
        ZhangExactMatrix generators(count,ZhangExactVector(dimension));
        ZhangExactVector coefficients(count),target(dimension);
        for (int r=0;r<count;++r)
        {
            coefficients[r]=int(random()%7)-3;
            for (int c=0;c<dimension;++c) generators[r][c]=int(random()%9)-4;
        }
        for (int r=0;r<count;++r) for (int c=0;c<dimension;++c)
            target[c]+=coefficients[r]*generators[r][c];
        for (int inconsistent=0;inconsistent<2;++inconsistent)
        {
            const auto old=zhangIntegerRowLatticeContains(generators,target,true);
            const auto fresh=zhangIntegerRowLatticeContains(generators,target,true,true);
            require(old.contained==fresh.contained && old.rank==fresh.rank &&
                old.smithInvariants==fresh.smithInvariants && old.combination==fresh.combination,
                "compact Smith witness/feasibility changed");
            if (inconsistent==0) require(fresh.contained,"known feasible integer generator rejected");
            if (fresh.contained)
            {
                ZhangExactVector actual(dimension);
                for (int r=0;r<count;++r) for (int c=0;c<dimension;++c)
                    actual[c]+=fresh.combination[r]*generators[r][c];
                require(actual==target,"independent compact witness multiply failed");
            }
            ++witnesses;
            ++target[0];
        }
    }
    require(!zhangIntegerRowLatticeContains({{2}},{3},true,true).contained,"2*N=3 accepted");
    require(!zhangIntegerRowLatticeContains({{0,0}},{0,1},true,true).contained,"0=1 accepted");
    auto infeasible=evidenceFor({{2}},{1});
    const auto infeasibleDomain=infeasible.compiledCurrentDomain();
    require(!infeasibleDomain->feasible &&
        !zhangCurrentEpochAffineDomainView(infeasibleDomain,infeasible.currentDomainIdentity()),
        "integer-infeasible 2*N=1 created verified affine authority");
    auto empty=evidenceFor({},{});
    const auto emptyDomain=empty.compiledCurrentDomain();
    require(emptyDomain->feasible && emptyDomain->columns.empty(),"empty current domain infeasible");
    const auto emptyView=zhangCurrentEpochAffineDomainView(emptyDomain,empty.currentDomainIdentity());
    const ZhangExactMatrix freeTargets{{1,0},{0,2}};
    compareImages(zhangR51PhysicalImage(freeTargets,{},{},freeTargets,2),
        zhangR51PhysicalImageOnVerifiedDomain(freeTargets,emptyView,freeTargets,2));
    compareImages(zhangR51PhysicalImage({{2,0},{0,1}}, {}, {}, {{2,0},{0,1}},2),
        zhangR51PhysicalImageOnVerifiedDomain({{2,0},{0,1}},emptyView,{{2,0},{0,1}},2));
    // Include non-unit, rank-deficient and redundant source equations. Each
    // current origin must equal the OLD affine quotient's exact origin.
    for (int trial=0;trial<120;++trial)
    {
        const int count=1+random()%4;
        ZhangExactMatrix source(count,ZhangExactVector(3));
        ZhangExactVector origin(3),rhs(count);
        for (auto& x:origin) x=int(random()%7)-3;
        for (int r=0;r<count;++r)
        {
            for (int c=0;c<3;++c) source[r][c]=int(random()%7)-3;
            if (source[r]==ZhangExactVector(3)) source[r][0]=2;
            for (int c=0;c<3;++c) rhs[r]+=source[r][c]*origin[c];
        }
        source.push_back(source.front());rhs.push_back(rhs.front());
        if (trial%2) {std::reverse(source.begin(),source.end());std::reverse(rhs.begin(),rhs.end());}
        auto evidence=evidenceFor(source,rhs);
        const auto domain=evidence.compiledCurrentDomain();
        require(domain->feasible,"known feasible current domain rejected");
        const auto view=zhangCurrentEpochAffineDomainView(domain,evidence.currentDomainIdentity());
        require(view && view.get()==&domain->affine,"immutable domain owner not shared");
        const int base=domain->columns.size(),dimension=base+2;
        auto conditioners=domain->basis;
        for (auto& row:conditioners) row.resize(dimension);
        ZhangExactMatrix targets(2,ZhangExactVector(dimension));
        for (auto& row:targets) for (int c=0;c<dimension;++c) row[c]=int(random()%7)-3;
        // An unused trailing ambient column tests full projector padding.
        if (trial%3==0) for (auto& row:targets) row.back()=0;
        for (auto work:{ZhangProductFrameWork::IMAGE_GENERATORS_ONLY,ZhangProductFrameWork::INTEGER_PROJECTOR})
        {
            const auto old=zhangR47CompileProductSearchFrame(targets,conditioners,domain->basisValues,dimension,work);
            const auto fresh=zhangR47CompileTargetOnVerifiedDomain(targets,view,dimension,work);
            compareFrames(old,fresh);
            compareImages(zhangR51PhysicalImage(targets,conditioners,domain->basisValues,targets,dimension,work),
                zhangR51PhysicalImageOnVerifiedDomain(targets,view,targets,dimension,work));
            ++images;
        }
        // Replace the old final whole-physical-domain rebuilding by the
        // residual-only admission and compare the ORIGINAL target values.
        const auto projected=zhangR51PhysicalImageOnVerifiedDomain(targets,view,targets,dimension);
        const std::size_t rank=projected.projector.size();
        if (rank)
        {
            ZhangExactMatrix z(rank,ZhangExactVector(rank));
            ZhangExactVector y(rank);
            for (std::size_t r=0;r<rank;++r)
            {
                z[r][r]=1+random()%3;
                y[r]=int(random()%7)-3;
                for (std::size_t c=r+1;c<rank;++c) z[r][c]=int(random()%5)-2;
            }
            const auto integers=zhangExactMatrixTimesColumn(z,y);
            ZhangExactMatrix lifted;ZhangExactVector liftedValues;
            require(zhangR51LiftRationalRows(z,integers,projected.projector,projected.offsets,
                lifted,liftedValues),"random lift fixture");
            const auto admission=zhangR51CertifyResidualImageAdmission(projected,targets,
                conditioners,domain->basisValues,z,integers,lifted,liftedValues);
            require(admission.valid && admission.residualSolution==y,"random residual solution changed");
            auto full=conditioners;auto fullRhs=domain->basisValues;
            full.insert(full.end(),lifted.begin(),lifted.end());
            fullRhs.insert(fullRhs.end(),liftedValues.begin(),liftedValues.end());
            const auto oldFinal=zhangR51PhysicalImage(targets,full,fullRhs,targets,dimension);
            require(oldFinal.valid && oldFinal.projector.empty() &&
                oldFinal.particularTarget==admission.targetValues,
                "residual-only final admission differs from old whole-physical reconstruction");
            ++admissions;
        }
        // Repeated image queries must not allocate another source origin or
        // run the old full affine quotient; target projection remains allowed.
        zhangR48ExactMetrics.clear();
        const ZhangExactMatrix baseTargets(1,ZhangExactVector(base,1));
        const auto shared=zhangR47CompileTargetOnVerifiedDomain(baseTargets,view,base);
        require(shared.valid && shared.affine==view,"same-domain target copied large affine owner");
        require(zhangR48ExactMetrics["R48_EXACT_FEASIBILITY"].calls==0,
            "verified-domain target rebuilt source affine feasibility");
        require(evidence.compiledCurrentDomain()==domain && evidence.domainCache->builds==1,
            "immutable evidence domain rebuilt");
    }
    auto evidence=evidenceFor({{2,1,0}},{1});
    const auto domain=evidence.compiledCurrentDomain();
    auto identity=evidence.currentDomainIdentity();
    const auto view=zhangCurrentEpochAffineDomainView(domain,identity);
    for (int mismatch=0;mismatch<7;++mismatch)
    {
        auto wrong=identity;
        if (mismatch==0) ++wrong.epochGpstSeconds;
        if (mismatch==1) ++wrong.constraintGeneration;
        if (mismatch==2) wrong.posteriorId="other";
        if (mismatch==3) wrong.rootId="other";
        if (mismatch==4) wrong.chartId="other";
        if (mismatch==5) wrong.system=E_Sys::GAL;
        if (mismatch==6) wrong.physicalColumnsHash="other";
        require(!zhangCurrentEpochAffineDomainView(domain,wrong),"cross-scope affine cache reused");
    }
    require(evidence.add({{"physical-0",1}},2,{proof()}),"generation fixture add");
    const auto next=evidence.compiledCurrentDomain();
    require(next!=domain && next->identity.constraintGeneration>identity.constraintGeneration &&
        !zhangCurrentEpochAffineDomainView(domain,next->identity),"stale generation accepted");
    ZhangExactInteger value;
    require(!domain->consequence({{"unknown",1}},value),"unknown current variable treated as zero");
    // Parity oracle independent of HNF: 2*old+live=1 permits EVERY integer
    // old, and 2*old takes precisely 2Z, never the saturated Z.
    const auto parity=zhangR51PhysicalImageOnVerifiedDomain({{2,0}},view,{{2,0}},2);
    require(parity.valid && parity.generators==ZhangExactMatrix{{2}},"2Z image saturated");
    const auto onlyImage=zhangR47CompileTargetOnVerifiedDomain({{2,0}},view,2);
    require(!zhangR47ProductConsequence(onlyImage,{2,0},0,value),
        "unrequested empty physical kernel certified a free target");
    const auto fullKernel=zhangR47CompileTargetOnVerifiedDomain({{2,0}},view,2,
        ZhangProductFrameWork::INTEGER_PROJECTOR);
    require(!zhangR47ProductConsequence(fullKernel,{2,0},0,value) &&
        zhangR47ProductConsequence(fullKernel,{2,1},0,value) && value==1,
        "full-kernel consequence behavior changed");
    for (int y=-6;y<=6;++y)
    {
        const auto result=admit(parity,{{2,0}},domain->basis,domain->basisValues,{{2}},{2*y});
        require(result.valid && result.targetValues==ZhangExactVector{2*y},"parity residual target changed");
        bool witness=false;
        for (int old=-10;old<=10;++old) for (int live=-25;live<=25;++live)
            witness |= 2*old+live==1 && 2*old==2*y;
        require(witness,"independent physical parity witness missing");
        ++admissions;
    }
    require(!admit(parity,{{2,0}},domain->basis,domain->basisValues,{{2}},{3}).valid,
        "noninteger residual solution allowed through nonprimitive image");
    auto saturated=parity;
    saturated.generators={{1}};
    saturated.projector={{2,0}};
    saturated.offsets={0};
    require(!admit(saturated,{{2,0}},domain->basis,domain->basisValues,{{1}},{1}).valid,
        "fake saturated image with correct T=G*W identity accepted");
    // Target-only variables are explicit free coordinates, not missing zeros.
    const auto unknown=zhangR51PhysicalImageOnVerifiedDomain({{0,0,1}},view,{{0,0,1}},3);
    require(unknown.valid && unknown.generators==ZhangExactMatrix{{1}} &&
        unknown.projector==std::vector<ZhangR51RationalRow>{{0,0,1}},"unknown target column lost free integer");
    // Independent small-system oracle: Z=[2 1;0 3], y=(-2,3), z=(-1,9).
    const ZhangExactMatrix targets{{1,0},{0,2}};
    const auto free=zhangR51PhysicalImage(targets,{}, {},targets,2);
    const auto result=admit(free,targets,{}, {},{{2,1},{0,3}},{-1,9});
    require(result.valid && result.residualSolution==ZhangExactVector{-2,3} &&
        result.targetValues==ZhangExactVector{-2,6},"nonunimodular exact residual admission failed");
    require(!admit(free,targets,{}, {},{{2,1},{0,3}},{-1,10}).valid,"fractional residual accepted");
    require(!admit(free,targets,{}, {},{{1,0},{2,0}},{1,2}).valid,"rank deficient residual accepted");
    auto corrupt=free;corrupt.offsets[0]=1;
    require(!admit(corrupt,targets,{}, {},{{1,0},{0,1}},{1,2}).valid,"mismatched target offset accepted");
    corrupt=free;corrupt.projector[0][1]=1;
    require(!admit(corrupt,targets,{}, {},{{1,0},{0,1}},{1,2}).valid,"mismatched original target accepted");
    ZhangExactMatrix lift;ZhangExactVector liftValues;
    require(zhangR51LiftRationalRows({{1,0},{0,1}},{1,2},free.projector,free.offsets,lift,liftValues),"bad lift fixture");
    ++liftValues[0];
    require(!zhangR51CertifyResidualImageAdmission(free,targets,{}, {},{{1,0},{0,1}},{1,2},lift,liftValues).valid,
        "wrong lifted RHS accepted");
    // A current-domain rewrite changes W/offset, NOT the original image.
    // old=live+7, target old. The legacy origin is (7,0), hence q=old-7;
    // rewriting gives q=live, not a new target image or a changed origin.
    auto translated=zhangR51PhysicalImage({{1,0}},{{1,-1}},{7},{{1,0}},2);
    require(translated.particularTarget==ZhangExactVector{7},"legacy translation origin changed");
    const auto rewrite=zhangRewriteCurrentBridgeRepresentation(
        translated.projector,translated.offsets,{{1,-1}},{7},{0});
    require(rewrite.valid,"rewrite fixture");
    translated.projector=rewrite.projector;translated.offsets=rewrite.offsets;
    const auto translatedAdmission=admit(translated,{{1,0}},{{1,-1}},{7},{{1}},{2});
    require(translatedAdmission.valid && translatedAdmission.targetValues==ZhangExactVector{9},
        "current rewrite invalidated original target admission");
    const auto fullyFixed=zhangR51PhysicalImage({{1}},{{2}},{6},{{1}},1);
    const auto already=zhangR51CertifyResidualImageAdmission(fullyFixed,{{1}},{{2}},{6},{},{},{},{});
    require(already.valid && already.targetValues==ZhangExactVector{3},"already-fixed image rejected");
    std::cout<<"PASS verified_domain_image witnesses="<<witnesses<<" images="<<images
        <<" independent_admissions="<<admissions<<'\n';
}
