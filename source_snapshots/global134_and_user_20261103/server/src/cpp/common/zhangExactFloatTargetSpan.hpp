#pragma once
#include "common/zhangProductFloatCarriers.hpp"
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>

// Gaussian FLOAT row-space algebra over Q. This is NOT an integer lattice,
// integer feasibility test, historical constraint, or publication authority.
// The caller owns immutable actual [T G] rows for this object's entire life.
// No cache survives that actual transform call or uses pointer/epoch identity.
namespace zhang_float_target_q_span
{
using Q=ZhangFloatCarrierRational;
using Sparse=ZhangFloatCarrierSparseRow;
using Rows=std::vector<std::vector<std::pair<std::size_t,Q>>>;

struct Result
{
    bool inRationalSpan=false,losslesslyEncodable=false;
    std::string reason;
    Sparse exactCoefficients;
    std::map<std::size_t,double> coefficients;
    // One Q witness need not find every possible dyadic solution when free
    // coefficients exist. An unencodable witness is a conservative refusal.
    bool alternativeDyadicWitnessSearchComplete=false;
};
struct Metrics
{
    std::size_t totalRows=0,componentsBuilt=0,rowsFactored=0,largestComponent=0;
};

class ImmutableComponentSpan
{
    struct BasisRow {Sparse row,witness;};
    using Basis=std::map<std::size_t,BasisRow>;
    const Rows& rows;
    std::size_t columns;
    std::vector<std::size_t> parent;
    std::map<std::size_t,std::size_t> columnOwner;
    std::map<std::size_t,std::vector<std::size_t>> members;
    std::map<std::size_t,Basis> cached;
    Metrics observed;

    std::size_t root(std::size_t row)
    {
        while(parent[row]!=row){parent[row]=parent[parent[row]];row=parent[row];}
        return row;
    }
    void unite(std::size_t a,std::size_t b)
    {a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);}
    static void add(Sparse& to,const Sparse& from,const Q& weight)
    {
        if(weight==0)return;
        for(const auto&[column,value]:from)
        {
            zhangExactWorkCheckpoint("FLOAT_Q_SPAN_SPARSE_ADD");
            const Q next=to[column]+weight*value;
            if(next==0)to.erase(column);else to[column]=next;
        }
    }
    static bool powerOfTwo(const ZhangExactInteger& value)
    {return value>0 && (value & (value-1))==0;}
    static int pivotPreference(const std::vector<std::pair<std::size_t,Q>>& row)
    {
        int best=3;
        for(const auto&[column,q]:row)
        {
            if(q==1 || q==-1)return 0;
            ZhangExactInteger numeratorAbs=numerator(q);
            if(numeratorAbs<0)numeratorAbs=-numeratorAbs;
            if(q!=0 && powerOfTwo(numeratorAbs) && powerOfTwo(denominator(q)))best=1;
            else if(q!=0)best=std::min(best,2);
        }
        return best;
    }
    const Basis& basis(std::size_t component)
    {
        if(const auto old=cached.find(component);old!=cached.end())return old->second;
        ZhangPhaseTimer timer("FLOAT_TARGET_EXACT_COMPONENT_FACTOR");
        auto ordered=members.at(component);
        std::sort(ordered.begin(),ordered.end(),[&](auto a,auto b)
        {
            return std::tuple{pivotPreference(rows[a]),rows[a].size(),a}<
                   std::tuple{pivotPreference(rows[b]),rows[b].size(),b};
        });
        Basis built;
        for(const auto index:ordered)
        {
            zhangExactWorkCheckpoint("FLOAT_Q_SPAN_FACTOR_ROW");
            Sparse row,witness{{index,1}};
            for(const auto&[column,value]:rows[index])if(value!=0)row[column]+=value;
            for(auto i=row.begin();i!=row.end();)if(i->second==0)i=row.erase(i);else ++i;
            for(const auto&[pivot,source]:built)
                if(const auto found=row.find(pivot);found!=row.end())
                {
                    const Q weight=found->second;
                    add(row,source.row,-weight);add(witness,source.witness,-weight);
                }
            if(row.empty())continue;
            const auto pivot=row.begin()->first;const Q divisor=row.begin()->second;
            for(auto&[column,value]:row){zhangExactWorkCheckpoint("FLOAT_Q_SPAN_NORMALIZE");value/=divisor;}
            for(auto&[column,value]:witness){zhangExactWorkCheckpoint("FLOAT_Q_SPAN_WITNESS_NORMALIZE");value/=divisor;}
            built.emplace(pivot,BasisRow{std::move(row),std::move(witness)});
        }
        ++observed.componentsBuilt;observed.rowsFactored+=ordered.size();
        observed.largestComponent=std::max(observed.largestComponent,ordered.size());
        return cached.emplace(component,std::move(built)).first->second;
    }

public:
    ImmutableComponentSpan(const Rows& actualRows,std::size_t actualColumns)
        :rows(actualRows),columns(actualColumns),parent(rows.size())
    {
        ZhangPhaseTimer timer("FLOAT_TARGET_EXACT_COMPONENT_INDEX");
        observed.totalRows=rows.size();std::iota(parent.begin(),parent.end(),0);
        for(std::size_t row=0;row<rows.size();++row)
            for(const auto&[column,value]:rows[row])if(value!=0)
            {
                zhangExactWorkCheckpoint("FLOAT_Q_SPAN_INDEX_NONZERO");
                if(column>=columns)throw std::invalid_argument("FLOAT_Q_SPAN_COLUMN_OUT_OF_RANGE");
                const auto [found,inserted]=columnOwner.emplace(column,row);
                if(!inserted)unite(row,found->second);
            }
        for(std::size_t row=0;row<rows.size();++row)members[root(row)].push_back(row);
        for(auto&[column,row]:columnOwner)row=root(row);
    }
    ImmutableComponentSpan(const ImmutableComponentSpan&)=delete;
    ImmutableComponentSpan& operator=(const ImmutableComponentSpan&)=delete;
    const Metrics& metrics() const noexcept{return observed;}

    Result resolve(const Sparse& target)
    {
        ZhangPhaseTimer timer("FLOAT_TARGET_EXACT_COMPONENT_QUERY");
        Result result;
        std::map<std::size_t,Sparse> partitions;
        for(const auto&[column,value]:target)if(value!=0)
        {
            if(column>=columns || !columnOwner.contains(column))
            {result.reason="TARGET_OUTSIDE_RATIONAL_ROW_SPAN";return result;}
            partitions[columnOwner.at(column)].emplace(column,value);
        }
        Sparse witness;
        for(auto&[component,remaining]:partitions)
        {
            const auto& factor=basis(component);
            for(const auto&[pivot,source]:factor)
                if(const auto found=remaining.find(pivot);found!=remaining.end())
                {
                    const Q weight=found->second;
                    add(remaining,source.row,-weight);add(witness,source.witness,weight);
                }
            if(!remaining.empty())
            {result.reason="TARGET_OUTSIDE_RATIONAL_ROW_SPAN";return result;}
        }
        // Audit against original actual rows, including ALL appended noise,
        // constants, phase, clock and other Gaussian/nuisance coordinates.
        Sparse pulled,canonicalTarget;
        for(const auto&[column,value]:target)if(value!=0)canonicalTarget[column]=value;
        for(const auto&[row,weight]:witness)
            for(const auto&[column,value]:rows.at(row))
            {
                zhangExactWorkCheckpoint("FLOAT_Q_SPAN_ORIGINAL_WITNESS_AUDIT");
                const Q next=pulled[column]+weight*value;
                if(next==0)pulled.erase(column);else pulled[column]=next;
            }
        if(pulled!=canonicalTarget)
        {result.reason="INTERNAL_EXACT_Q_WITNESS_AUDIT_FAILED";return result;}
        result.inRationalSpan=true;result.exactCoefficients=std::move(witness);
        for(const auto&[row,weight]:result.exactCoefficients)
        {
            const double value=weight.convert_to<double>();
            const auto roundtrip=zhangFloatCarrierExactDouble(value);
            if(!roundtrip || *roundtrip!=weight)
            {result.reason="CHOSEN_Q_WITNESS_NOT_LOSSLESSLY_ENCODABLE";return result;}
            if(value!=0)result.coefficients.emplace(row,value);
        }
        result.losslesslyEncodable=true;result.reason="EXACT_Q_WITNESS_LOSSLESS_IEEE_COEFFICIENTS";
        return result;
    }
};
} // namespace zhang_float_target_q_span
