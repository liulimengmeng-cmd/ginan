#include "common/zhangExactFloatTargetSpan.hpp"
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace zhang_float_target_q_span;
static void check(bool value,const char* why){if(!value)throw std::runtime_error(why);}

static void smallLiteralCases()
{
    {
        const Rows rows={{{0,1},{1,1}},{{1,1}}};ImmutableComponentSpan span(rows,2);
        const auto r=span.resolve({{0,1},{1,2}});
        check(r.losslesslyEncodable && r.coefficients==std::map<std::size_t,double>{{0,1},{1,1}},"literal shear witness differs");
        const auto before=span.metrics();span.resolve({{0,2},{1,3}});
        check(span.metrics().componentsBuilt==before.componentsBuilt && span.metrics().rowsFactored==before.rowsFactored,
              "same immutable component was refactored");
    }
    for(const int multiplier:{2,3})
    {
        const Rows rows={{{0,multiplier}}};ImmutableComponentSpan span(rows,1);const auto r=span.resolve({{0,1}});
        check(r.inRationalSpan,"nonunit FLOAT scale was treated as integer nonmembership");
        if(multiplier==2)check(r.losslesslyEncodable && r.coefficients.at(0)==.5,"dyadic half witness rejected");
        else check(!r.losslesslyEncodable && r.reason=="CHOSEN_Q_WITNESS_NOT_LOSSLESSLY_ENCODABLE" &&
            !r.alternativeDyadicWitnessSearchComplete,"one-third witness was rounded or called mathematically impossible");
    }
    {
        const Rows rows={{{0,3}},{{0,1}}};ImmutableComponentSpan span(rows,1);const auto r=span.resolve({{0,1}});
        check(r.losslesslyEncodable && r.coefficients==std::map<std::size_t,double>{{1,1}},"unit row did not take priority over3");
    }
    {
        const Rows rows={{{0,3}},{{0,5}}};ImmutableComponentSpan span(rows,1);const auto r=span.resolve({{0,1}});
        check(r.inRationalSpan && !r.losslesslyEncodable && !r.alternativeDyadicWitnessSearchComplete,
            "free-variable non-dyadic witness refusal falsely claimed completeness");
        check(Q(3)*2+Q(5)*(-1)==1,"independent alternative dyadic Bezout witness lost");
    }
    {
        const Rows redundant={{{0,1}},{{0,1}}};ImmutableComponentSpan span(redundant,2);
        check(span.resolve({{0,1}}).losslesslyEncodable,"redundant Gaussian rows changed representability");
        check(!span.resolve({{1,1}}).inRationalSpan,"unknown column was silently zero");
        const Rows sum={{{0,1},{1,1}}};ImmutableComponentSpan onlySum(sum,2);
        check(!onlySum.resolve({{0,1}}).inRationalSpan,"sum recovered an independent coordinate");
    }
}

static void exactNoiseAndWholeNumericTerms()
{
    const Q eps=*zhangFloatCarrierExactDouble(1e-12);
    {
        // Actual A=[T G], with the same random source on two outputs.
        const Rows rows={{{0,1},{2,eps}},{{1,1},{2,eps}}};ImmutableComponentSpan span(rows,3);
        const auto cancel=span.resolve({{0,1},{1,-1}});
        check(cancel.losslesslyEncodable && cancel.coefficients==std::map<std::size_t,double>{{0,1},{1,-1}},
              "actual shared tiny-noise cancellation was rejected");
        check(!span.resolve({{0,1}}).inRationalSpan,"1e-12 new noise was discarded");
    }
    {
        const Rows rows={{{0,1},{1,eps}},{{0,1},{1,eps}}};ImmutableComponentSpan span(rows,2);
        check(!span.resolve({{0,1}}).inRationalSpan,"duplicated noisy rows invented a noise-free old state");
    }
    {
        // Source columns: ONE, phase, ambiguity, clock; no symbolic shortcut.
        const Rows rows={{{0,1}},{{1,1},{2,1}},{{2,1}},{{3,1}}};ImmutableComponentSpan span(rows,4);
        const auto r=span.resolve({{0,1},{1,1},{2,2},{3,Q(1)/4}});
        check(r.losslesslyEncodable && r.coefficients==std::map<std::size_t,double>{{0,1},{1,1},{2,1},{3,.25}},
              "full ONE/phase/ambiguity/clock row was not preserved");
        const Rows lostClock={rows[0],rows[1],rows[2]};ImmutableComponentSpan lost(lostClock,4);
        check(!lost.resolve({{0,1},{1,1},{2,2},{3,Q(1)/4}}).inRationalSpan,"lost clock nuisance was omitted");
    }
}

static void lazyDisconnectedComponents()
{
    Rows rows;for(std::size_t i=0;i<7000;++i)rows.push_back({{i,1}});
    rows.push_back({{7000,1},{7001,1}});rows.push_back({{7001,1}});
    ImmutableComponentSpan span(rows,7002);
    check(span.metrics().componentsBuilt==0,"unqueried components were eagerly factored");
    check(span.resolve({{7000,1},{7001,2}}).losslesslyEncodable,"large disconnected shear failed");
    check(span.metrics().totalRows==7002 && span.metrics().componentsBuilt==1 && span.metrics().rowsFactored==2 &&
        span.metrics().largestComponent==2,"7000 untouched identity rows entered a dense/global decomposition");
    check(span.resolve({{12,1}}).losslesslyEncodable && span.metrics().rowsFactored==3,"second tiny component was not independently lazy");
    // One appended random source connects rows that otherwise have disjoint
    // source support. Closure must include BOTH before deciding cancellation.
    const Q eps=*zhangFloatCarrierExactDouble(1e-12);
    const Rows noisy={{{0,1},{2,eps}},{{1,1},{2,eps}}};ImmutableComponentSpan noise(noisy,3);
    check(noise.resolve({{0,1},{1,-1}}).losslesslyEncodable && noise.metrics().largestComponent==2,
          "noise column failed to join its actual connected component");
}

static void independentTwoDimensionalOracle()
{
    std::mt19937 random(20261102);
    for(int fixture=0;fixture<256;++fixture)
    {
        std::vector<std::array<int,2>> original;Rows rows;
        for(int i=0,n=random()%5;i<n;++i)
        {
            const std::array<int,2> a={int(random()%7)-3,int(random()%7)-3};original.push_back(a);
            std::vector<std::pair<std::size_t,Q>> row;
            for(int c=0;c<2;++c)if(a[c])row.emplace_back(c,a[c]);rows.push_back(std::move(row));
        }
        const std::array<int,2> target={int(random()%7)-3,int(random()%7)-3};
        std::array<int,2> first{};bool have=false,rank2=false;
        for(const auto& row:original)if(row[0] || row[1])
        {
            if(!have){first=row;have=true;}
            else rank2|=first[0]*row[1]-first[1]*row[0]!=0;
        }
        const bool expected=rank2 || (have?first[0]*target[1]-first[1]*target[0]==0:target[0]==0 && target[1]==0);
        ImmutableComponentSpan span(rows,2);Sparse t;
        for(int c=0;c<2;++c)if(target[c])t.emplace(c,target[c]);
        const auto r=span.resolve(t);check(r.inRationalSpan==expected,"exact sparse solver disagrees with independent determinant-rank oracle");
        if(r.inRationalSpan)for(int c=0;c<2;++c)
        {
            Q literal=0;for(const auto&[row,value]:r.exactCoefficients)literal+=value*original.at(row)[c];
            check(literal==target[c],"independent literal witness dot product differs");
        }
    }
}

int main()
{
    smallLiteralCases();exactNoiseAndWholeNumericTerms();lazyDisconnectedComponents();independentTwoDimensionalOracle();
    std::cout<<"PASS exact_FLOAT_Q_rowspan shear=1 dyadic_half=1 nondyadic_conservative=1 unit_priority=1 "
        "tiny_noise_retained=1 real_noise_cancellation=1 full_ONE_clock_terms=1 "
        "7000_unqueried_components_not_factored=1 independent_2D_determinant_oracle=256 integer_authority=0\n";
}
