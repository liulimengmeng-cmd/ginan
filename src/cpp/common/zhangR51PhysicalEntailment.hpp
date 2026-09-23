#pragma once
#include <chrono>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include "common/zhangR47ProductDomain.hpp"

/** One immutable physical H*x=b snapshot, scoped to a single certificate scan.
 * Integer feasibility is checked once. Logical consequences are then tested
 * against the rational row span of H with exact arithmetic; no particular
 * integer solution or saturated physical kernel is needed for this query.
 * This is not an integer row-lattice membership or certificate-witness API.
 */
class ZhangR51PhysicalEntailmentSnapshot
{
public:
    using Row = std::map<std::string, ZhangExactInteger>;
    using Rows = std::vector<Row>;
    ZhangR51PhysicalEntailmentSnapshot(const Rows& rows,
        const ZhangExactVector& values, bool valid)
        : rows_(rows), values_(values), valid_(valid && rows.size()==values.size())
    {
        for (const auto& row : rows_)
            for (const auto& [id,value] : row)
                if (value!=0)
                {
                    ++inputNonzeros;
                    if (!originalColumns_.contains(id)) originalColumns_[id]=originalColumns_.size();
                }
    }
    bool entails(const Rows& targets, const ZhangExactVector& rhs)
    {
        ++queries;
        if (!valid_ || targets.size()!=rhs.size()) return false;
        for (const auto& row : targets)
            for (const auto& [id,value] : row)
                if (value!=0 && !originalColumns_.contains(id)) return false;
        if (!prepared_)
        {
            prepared_=true;
            const auto start=std::chrono::steady_clock::now();
            Rows sparse=rows_;
            ZhangExactVector sparseValues=values_;
            const auto sparseStart=std::chrono::steady_clock::now();
            const bool sparseConsistent=eliminateUnits(sparse,sparseValues);
            sparseSeconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-sparseStart).count();
            if (sparseConsistent)
            {
                for (const auto& row:sparse)
                    for (const auto& [id,value]:row)
                        if (value!=0 && !columns_.contains(id)) columns_[id]=columns_.size();
                residualRows=sparse.size();
                residualColumns=columns_.size();
                for (const auto& row:sparse) residualNonzeros+=row.size();
                ZhangExactMatrix h(sparse.size(),ZhangExactVector(columns_.size()));
                for (std::size_t i=0;i<sparse.size();++i)
                    for (const auto& [id,value]:sparse[i])
                        h[i][columns_.at(id)]=value;
                const auto hnfStart=std::chrono::steady_clock::now();
                auto hnf=zhangExactRowHermiteNormalForm(std::move(h),sparseValues);
                hnfSeconds=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-hnfStart).count();
                if (hnf.consistent)
                {
                    basis_=std::move(hnf.basis);
                    basisValues_=std::move(hnf.values);
                    basisRows=basis_.size();
                    for (const auto& row : basis_)
                    {
                        const auto pivot=std::find_if(row.begin(),row.end(),
                            [](const auto& value){return value!=0;});
                        pivots_.push_back(std::distance(row.begin(),pivot));
                        for (const auto& value : row) basisNonzeros+=value!=0;
                    }
                    const auto feasibilityStart=std::chrono::steady_clock::now();
                    feasible_=zhangR47AffineIntegerFeasible(
                        basis_,basisValues_,columns_.size());
                    feasibilitySeconds=std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-feasibilityStart).count();
                }
            }
            ++builds;
            buildSeconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count();
        }
        if (!feasible_) return false;
        bool direct=true;
        for (std::size_t i=0;i<targets.size();++i)
        {
            bool found=false;
            for (std::size_t j=0;j<rows_.size();++j)
                found |= targets[i]==rows_[j] && rhs[i]==values_[j];
            direct &= found;
        }
        if (direct) { ++directQueries; return true; }
        using Rational=boost::multiprecision::cpp_rational;
        for (std::size_t i=0;i<targets.size();++i)
        {
            Row reduced=targets[i];
            ZhangExactInteger reducedRhs=rhs[i];
            for (const auto& substitution:substitutions_)
            {
                const auto pivot=reduced.find(substitution.variable);
                if (pivot==reduced.end() || pivot->second==0) continue;
                const ZhangExactInteger factor=pivot->second;
                reduced.erase(pivot);
                reducedRhs-=factor*substitution.constant;
                for (const auto& [id,coefficient]:substitution.coefficients)
                {
                    reduced[id]+=factor*coefficient;
                    if (reduced[id]==0) reduced.erase(id);
                }
            }
            std::vector<Rational> remaining(columns_.size());
            for (const auto& [id,value] : reduced)
            {
                if (value!=0 && !columns_.contains(id)) return false;
                if (value!=0) remaining[columns_.at(id)]=Rational(value);
            }
            Rational remainingRhs(reducedRhs);
            for (std::size_t row=0;row<basis_.size();++row)
            {
                const auto pivot=pivots_[row];
                if (remaining[pivot]==0) continue;
                const Rational factor=remaining[pivot]/Rational(basis_[row][pivot]);
                for (std::size_t column=pivot;column<remaining.size();++column)
                    if (basis_[row][column]!=0)
                        remaining[column]-=factor*Rational(basis_[row][column]);
                remainingRhs-=factor*Rational(basisValues_[row]);
            }
            if (remainingRhs!=0 ||
                std::any_of(remaining.begin(),remaining.end(),
                    [](const auto& value){return value!=0;})) return false;
        }
        return true;
    }
    std::size_t rowCount() const { return rows_.size(); }
    std::size_t columnCount() const { return originalColumns_.size(); }
    std::size_t queries=0, directQueries=0, builds=0;
    std::size_t inputNonzeros=0, unitPivots=0, fillAdded=0;
    std::size_t residualRows=0, residualColumns=0, residualNonzeros=0;
    std::size_t basisRows=0, basisNonzeros=0;
    double buildSeconds=0, sparseSeconds=0, hnfSeconds=0, feasibilitySeconds=0;
private:
    struct UnitSubstitution
    {
        std::string variable;
        ZhangExactInteger constant;
        Row coefficients;
    };
    bool eliminateUnits(Rows& rows,ZhangExactVector& values)
    {
        std::map<std::string,std::set<std::size_t>> adjacency;
        std::vector<bool> active(rows.size(),true);
        std::deque<std::size_t> pending;
        for (std::size_t row=0;row<rows.size();++row)
        {
            for (auto it=rows[row].begin();it!=rows[row].end();)
            {
                if (it->second==0) it=rows[row].erase(it);
                else {adjacency[it->first].insert(row);++it;}
            }
            if (rows[row].empty())
            {
                if (values[row]!=0) return false;
                active[row]=false;
            }
            else pending.push_back(row);
        }
        while (!pending.empty())
        {
            const auto row=pending.front();pending.pop_front();
            if (!active[row]) continue;
            std::optional<std::string> variable;
            std::size_t bestCost=std::numeric_limits<std::size_t>::max();
            for (const auto& [id,coefficient]:rows[row])
                if (coefficient==1 || coefficient==-1)
                {
                    const auto cost=adjacency[id].size()*rows[row].size();
                    if (cost<bestCost) {variable=id;bestCost=cost;}
                }
            if (!variable) continue;
            const ZhangExactInteger sign=rows[row].at(*variable);
            UnitSubstitution substitution;
            substitution.variable=*variable;
            substitution.constant=sign*values[row];
            for (const auto& [id,coefficient]:rows[row])
                if (id!=*variable) substitution.coefficients[id]=-sign*coefficient;
            for (const auto& [id,coefficient]:rows[row])
                if (id!=*variable) adjacency[id].erase(row);
            active[row]=false;
            auto neighbors=std::move(adjacency[*variable]);
            adjacency.erase(*variable);
            for (const auto other:neighbors)
            {
                if (other==row || !active[other]) continue;
                auto pivot=rows[other].find(*variable);
                if (pivot==rows[other].end()) continue;
                const ZhangExactInteger factor=pivot->second;
                rows[other].erase(pivot);
                values[other]-=factor*substitution.constant;
                for (const auto& [id,coefficient]:substitution.coefficients)
                {
                    auto existing=rows[other].find(id);
                    const ZhangExactInteger previous=existing==rows[other].end()
                        ? ZhangExactInteger(0) : existing->second;
                    const ZhangExactInteger updated=previous+factor*coefficient;
                    if (updated==0)
                    {
                        if (existing!=rows[other].end())
                        {
                            rows[other].erase(existing);
                            adjacency[id].erase(other);
                        }
                    }
                    else
                    {
                        if (existing==rows[other].end())
                        {
                            ++fillAdded;
                            adjacency[id].insert(other);
                        }
                        rows[other][id]=updated;
                    }
                }
                if (rows[other].empty())
                {
                    if (values[other]!=0) return false;
                    active[other]=false;
                }
                else pending.push_back(other);
            }
            substitutions_.push_back(std::move(substitution));
            ++unitPivots;
        }
        Rows residual;
        ZhangExactVector residualValues;
        for (std::size_t row=0;row<rows.size();++row)
            if (active[row])
            {
                residual.push_back(std::move(rows[row]));
                residualValues.push_back(std::move(values[row]));
            }
        rows=std::move(residual);
        values=std::move(residualValues);
        return true;
    }
    const Rows rows_;
    const ZhangExactVector values_;
    const bool valid_;
    std::map<std::string,int> originalColumns_, columns_;
    bool prepared_=false;
    bool feasible_=false;
    std::vector<UnitSubstitution> substitutions_;
    ZhangExactMatrix basis_;
    ZhangExactVector basisValues_;
    std::vector<std::size_t> pivots_;
};
