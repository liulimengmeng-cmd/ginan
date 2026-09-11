#pragma once
#include "common/zhangProductPhysicalCycleChart.hpp"
#include "common/zhangIntegerDecisionProof.hpp"
#include "common/zhangR47Candidate.hpp"

struct ZhangR47TransportedHistory
{
    bool valid=false;
    ZhangExactMatrix rows;
    ZhangExactVector values;
    std::vector<ZhangDecisionProofs> parents;
    int inputRank=0, survivingRank=0;
    std::string reason="NOT_EVALUATED";
};
inline ZhangR47TransportedHistory zhangR47TransportHistory(
    const std::vector<std::map<std::string,ZhangExactInteger>>& physical,
    const ZhangExactVector& values,const std::vector<ZhangDecisionProofs>& parents,
    const ZhangProductPhysicalCycleChart& chart)
{
    ZhangR47TransportedHistory out;
    if(physical.size()!=values.size() || physical.size()!=parents.size()) return out;
    std::set<std::string> supported;
    for(const auto& [column,expansion]:chart.expansions)
        for(const auto& [id,value]:expansion) if(value!=0) supported.insert(id);
    std::map<std::string,int> columns;
    for(const auto& row:physical) for(const auto& [id,value]:row)
        if(value!=0 && !columns.contains(id)) columns[id]=columns.size();
    ZhangExactMatrix rows(physical.size(),ZhangExactVector(columns.size()));
    std::vector<bool> available(columns.size());
    std::vector<std::string> names(columns.size());
    for(const auto& [id,column]:columns) {available[column]=supported.contains(id);names[column]=id;}
    for(std::size_t r=0;r<physical.size();++r)
        for(const auto& [id,value]:physical[r]) if(value!=0) rows[r][columns.at(id)]=value;
    const auto surviving=zhangExactSurvivingLattice(rows,values,available,true);
    if(!surviving.consistent) {out.reason="WHOLE_HISTORY_AFFINE_CONFLICT";return out;}
    out.inputRank=zhangExactRowHermiteNormalForm(rows).basis.size();
    for(std::size_t r=0;r<surviving.basis.size();++r)
    {
        std::map<std::string,ZhangExactInteger> expansion;
        int compact=0;
        for(std::size_t c=0;c<available.size();++c) if(available[c])
        { const auto value=surviving.basis[r][compact++]; if(value!=0) expansion[names[c]]=value; }
        ZhangExactVector projected;
        if(!chart.project(expansion,projected,&out.reason)) return out;
        ZhangDecisionProofs dependencies;
        for(std::size_t p=0;p<parents.size();++p)
            if(surviving.rowTransform[r][p]!=0) dependencies=zhangMergeDecisionProofs(dependencies,parents[p]);
        if(dependencies.empty() || !zhangDecisionRiskClosure(dependencies).valid)
        {out.reason="TRANSPORTED_HISTORY_PARENT_INCOMPLETE";return out;}
        out.rows.push_back(std::move(projected));out.values.push_back(surviving.values[r]);
        out.parents.push_back(std::move(dependencies));
    }
    out.survivingRank=out.rows.size();out.valid=true;
    out.reason="WHOLE_HISTORY_EXACT_CANCELLATION";return out;
}

struct ZhangR47HistorySubset
{
    ZhangExactMatrix rows;
    ZhangExactVector values;
    std::vector<int> selected;
    std::vector<std::string> reasons;
    ZhangDecisionProofs parents;
    bool valid=false;
    double risk=1;
};
/** Greedy budget-feasible selection runs before conditioning. The gain is a
 * deterministic score supplied from one posterior, never a risk waiver. */
inline ZhangR47HistorySubset zhangR47SelectHistorySubset(
    const ZhangExactMatrix& rows,const ZhangExactVector& values,
    const std::vector<ZhangDecisionProofs>& parents,const ZhangDecisionProofs& baseline,
    const std::vector<double>& gain,int dimension,double ceiling,double reserve)
{
    ZhangR47HistorySubset out;out.parents=baseline;
    out.reasons.assign(rows.size(),"NOT_SELECTED");
    if(rows.size()!=values.size() || rows.size()!=parents.size() || rows.size()!=gain.size() ||
       !std::isfinite(ceiling) || !std::isfinite(reserve) || reserve<0 || reserve>ceiling) return out;
    const auto base=zhangDecisionRiskClosure(baseline);
    if(!base.valid || base.bound>ceiling-reserve) return out;
    std::vector<int> order(rows.size());std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return gain[a]>gain[b];});
    for(const auto index:order)
    {
        if(parents[index].empty() || !zhangDecisionRiskClosure(parents[index]).valid)
        {out.reasons[index]="RISK_PARENT_INCOMPLETE";continue;}
        const auto merged=zhangMergeDecisionProofs(out.parents,parents[index]);
        const auto closure=zhangDecisionRiskClosure(merged);
        if(!closure.valid || closure.bound>ceiling-reserve)
        {out.reasons[index]="BUDGET_RESERVED_FOR_NEW_SEARCH";continue;}
        auto candidateRows=out.rows; auto candidateValues=out.values;
        candidateRows.push_back(rows[index]);candidateValues.push_back(values[index]);
        if(!zhangR47AffineIntegerFeasible(candidateRows,candidateValues,dimension))
        {out.reasons[index]="EXACT_AFFINE_CONFLICT";continue;}
        if(zhangExactRowHermiteNormalForm(candidateRows).basis.size()==out.rows.size())
        {out.reasons[index]="EXACT_REDUNDANCY_NOT_CONDITIONED";continue;}
        out.rows=std::move(candidateRows);out.values=std::move(candidateValues);
        out.parents=merged;out.selected.push_back(index);out.reasons[index]="ADMITTED_BEFORE_CONDITIONING";
    }
    out.risk=zhangDecisionRiskClosure(out.parents).bound;out.valid=true;return out;
}
