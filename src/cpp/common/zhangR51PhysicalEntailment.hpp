#pragma once
#include <chrono>
#include <map>
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
                    if (!columns_.contains(id)) columns_[id]=columns_.size();
                }
    }
    bool entails(const Rows& targets, const ZhangExactVector& rhs)
    {
        ++queries;
        if (!valid_ || targets.size()!=rhs.size()) return false;
        for (const auto& row : targets)
            for (const auto& [id,value] : row)
                if (value!=0 && !columns_.contains(id)) return false;
        if (!prepared_)
        {
            prepared_=true;
            const auto start=std::chrono::steady_clock::now();
            ZhangExactMatrix h(rows_.size(),ZhangExactVector(columns_.size()));
            for (std::size_t i=0;i<rows_.size();++i)
                for (const auto& [id,value] : rows_[i])
                    if (value!=0) h[i][columns_.at(id)]=value;
            const auto hnfStart=std::chrono::steady_clock::now();
            auto hnf=zhangExactRowHermiteNormalForm(std::move(h),values_);
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
            std::vector<Rational> remaining(columns_.size());
            for (const auto& [id,value] : targets[i])
                if (value!=0) remaining[columns_.at(id)]=Rational(value);
            Rational remainingRhs(rhs[i]);
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
    std::size_t columnCount() const { return columns_.size(); }
    std::size_t queries=0, directQueries=0, builds=0;
    std::size_t inputNonzeros=0, basisRows=0, basisNonzeros=0;
    double buildSeconds=0, hnfSeconds=0, feasibilitySeconds=0;
private:
    const Rows rows_;
    const ZhangExactVector values_;
    const bool valid_;
    std::map<std::string,int> columns_;
    bool prepared_=false;
    bool feasible_=false;
    ZhangExactMatrix basis_;
    ZhangExactVector basisValues_;
    std::vector<std::size_t> pivots_;
};
