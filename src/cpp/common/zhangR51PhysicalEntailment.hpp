#pragma once
#include <chrono>
#include <map>
#include <string>
#include <vector>
#include "common/zhangR47ProductDomain.hpp"

/** One immutable physical H*x=b snapshot, scoped to a single certificate scan.
 * No epoch/global cache and no pointer to mutable runtime state. Exact integer
 * feasibility, a particular integer solution and the saturated kernel are
 * computed once. Entailment is t*K=0 and t*x0=rhs, using exact arithmetic.
 * Unknown physical variables are independent free integers, never zeroed.
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
                if (value!=0 && !columns_.contains(id)) columns_[id]=columns_.size();
    }
    bool entails(const Rows& targets, const ZhangExactVector& rhs)
    {
        ++queries;
        if (!valid_ || targets.size()!=rhs.size()) return false;
        // Preserve the existing direct-receipt fast path exactly.
        bool direct=true;
        for (std::size_t i=0;i<targets.size();++i)
        {
            bool found=false;
            for (std::size_t j=0;j<rows_.size();++j)
                found |= targets[i]==rows_[j] && rhs[i]==values_[j];
            direct &= found;
        }
        if (direct) { ++directQueries; return true; }
        for (const auto& row : targets)
            for (const auto& [id,value] : row)
                if (value!=0 && !columns_.contains(id)) return false;
        if (columns_.empty()) return false;
        if (!prepared_)
        {
            prepared_=true;
            const auto start=std::chrono::steady_clock::now();
            ZhangExactMatrix h(rows_.size(),ZhangExactVector(columns_.size()));
            for (std::size_t i=0;i<rows_.size();++i)
                for (const auto& [id,value] : rows_[i])
                    if (value!=0) h[i][columns_.at(id)]=value;
            affine_=zhangExactAffineIntegerQuotient(h,values_,columns_.size(),
                ZhangExactQuotientWork::PARTICULAR_AND_KERNEL);
            ++builds;
            buildSeconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count();
        }
        if (!affine_.valid) return false;
        for (std::size_t i=0;i<targets.size();++i)
        {
            for (const auto& kernel : affine_.kernelBasis)
            {
                ZhangExactInteger residual=0;
                for (const auto& [id,value] : targets[i])
                    if (value!=0) residual+=value*kernel[columns_.at(id)];
                if (residual!=0) return false;
            }
            ZhangExactInteger value=0;
            for (const auto& [id,coefficient] : targets[i])
                if (coefficient!=0) value+=coefficient*affine_.particularSolution[columns_.at(id)];
            if (value!=rhs[i]) return false;
        }
        return true;
    }
    std::size_t rowCount() const { return rows_.size(); }
    std::size_t columnCount() const { return columns_.size(); }
    std::size_t queries=0, directQueries=0, builds=0;
    double buildSeconds=0;
private:
    const Rows rows_;
    const ZhangExactVector values_;
    const bool valid_;
    std::map<std::string,int> columns_;
    bool prepared_=false;
    ZhangExactAffineIntegerQuotient affine_;
};
