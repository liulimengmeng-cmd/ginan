#pragma once
#include <set>
#include <vector>
#include <cstddef>
#include "common/zhangR47Candidate.hpp"

/** Exact compatibility in an already certified affine integer image.
 * The caller has proved t=t0+G*q, q in Z^rank, with G an integer image
 * basis. It lifts B*q=z into physical rows with zhangR51LiftRationalRows.
 * Thus physical feasibility is exactly feasibility of the accepted B rows
 * plus the proposed B rows. Retired physical arcs and nonprimitive images
 * are represented by the initial image, never discarded by this class.
 * One instance belongs to one immutable image in one block-search call.
 */
class ZhangR51QuotientCompatibility
{
public:
    class Trial {
        friend class ZhangR51QuotientCompatibility;
        const ZhangR51QuotientCompatibility* owner_=nullptr;
        std::size_t generation_=0;
        ZhangExactMatrix rows_;
        ZhangExactVector values_;
        bool compatible_=false;
    public:
        bool compatible() const { return compatible_; }
    };
    explicit ZhangR51QuotientCompatibility(int rank,bool sourceValid)
        : rank_(rank), valid_(sourceValid && rank>0), columnRows_(std::max(0,rank)) {}
    ZhangR51QuotientCompatibility(const ZhangR51QuotientCompatibility&)=delete;
    ZhangR51QuotientCompatibility& operator=(const ZhangR51QuotientCompatibility&)=delete;

    Trial assess(const ZhangExactMatrix& reduced,const ZhangExactVector& integers,
                 const std::vector<int>& subset)
    {
        ++queries;
        Trial trial;trial.owner_=this;trial.generation_=generation_;
        if(!valid_ || reduced.empty() || reduced.size()!=integers.size() ||
           !zhangExactRectangularMatrix(reduced,subset.size()))return trial;
        std::set<int> selected;
        for(int c:subset)if(c<0 || c>=rank_ || !selected.insert(c).second)return trial;
        trial.rows_=ZhangExactMatrix(reduced.size(),ZhangExactVector(rank_));
        trial.values_=integers;
        for(std::size_t i=0;i<reduced.size();++i)
            for(std::size_t j=0;j<subset.size();++j)trial.rows_[i][subset[j]]=reduced[i][j];
        // Accepted components are already feasible. Only the connected closure
        // touched by a new row can change feasibility. Disjoint new blocks must
        // not repeatedly diagonalise every previously fixed quotient variable.
        std::set<int> columns;
        std::vector<int> pending;
        for(const auto& row:trial.rows_)for(int c:subset)
            if(row[c]!=0 && columns.insert(c).second)pending.push_back(c);
        std::set<std::size_t> touched;
        for(std::size_t next=0;next<pending.size();++next)
            for(const auto row:columnRows_[pending[next]])if(touched.insert(row).second)
                for(const int c:acceptedSupport_[row])
                    if(columns.insert(c).second)pending.push_back(c);
        ZhangExactMatrix rows;
        ZhangExactVector values;
        auto append=[&](const ZhangExactVector& source,const ZhangExactInteger& value) {
            ZhangExactVector compact;
            compact.reserve(columns.size());
            for(const int c:columns)compact.push_back(source[c]);
            rows.push_back(std::move(compact));values.push_back(value);
        };
        for(const auto row:touched)append(acceptedRows_[row],acceptedValues_[row]);
        for(std::size_t row=0;row<trial.rows_.size();++row)append(trial.rows_[row],integers[row]);
        ++componentChecks;
        if(touched.empty())++disjointChecks;
        checkedRows+=rows.size();checkedColumns+=columns.size();
        maximumCheckedColumns=std::max(maximumCheckedColumns,columns.size());
        if(columns.empty())
            trial.compatible_=std::all_of(values.begin(),values.end(),[](const auto& v){return v==0;});
        else trial.compatible_=zhangR47AffineIntegerFeasible(rows,values,columns.size());
        return trial;
    }
    // Commit only after the caller's numerical conditioning and proof gates.
    // Rejected proposals never constrain later candidates. Tickets cannot be
    // reused after another commit or transferred to another image instance.
    bool commit(const Trial& trial)
    {
        if(!valid_ || !trial.compatible_ || trial.owner_!=this || trial.generation_!=generation_)return false;
        acceptedRows_.insert(acceptedRows_.end(),trial.rows_.begin(),trial.rows_.end());
        acceptedValues_.insert(acceptedValues_.end(),trial.values_.begin(),trial.values_.end());
        for(const auto& row:trial.rows_) {
            const auto rowId=acceptedSupport_.size();
            std::vector<int> support;
            for(int c=0;c<rank_;++c)if(row[c]!=0) {
                support.push_back(c);columnRows_[c].push_back(rowId);
            }
            acceptedSupport_.push_back(std::move(support));
        }
        ++generation_;return true;
    }
    std::size_t acceptedRows() const { return acceptedRows_.size(); }
    std::size_t queries=0;
    std::size_t componentChecks=0,disjointChecks=0,checkedRows=0,checkedColumns=0,maximumCheckedColumns=0;
private:
    const int rank_;
    const bool valid_;
    std::size_t generation_=0;
    ZhangExactMatrix acceptedRows_;
    ZhangExactVector acceptedValues_;
    std::vector<std::vector<std::size_t>> columnRows_;
    std::vector<std::vector<int>> acceptedSupport_;
};
