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
        : rank_(rank), valid_(sourceValid && rank>0) {}
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
        auto rows=acceptedRows_;rows.insert(rows.end(),trial.rows_.begin(),trial.rows_.end());
        auto values=acceptedValues_;values.insert(values.end(),integers.begin(),integers.end());
        trial.compatible_=zhangR47AffineIntegerFeasible(rows,values,rank_);
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
        ++generation_;return true;
    }
    std::size_t acceptedRows() const { return acceptedRows_.size(); }
    std::size_t queries=0;
private:
    const int rank_;
    const bool valid_;
    std::size_t generation_=0;
    ZhangExactMatrix acceptedRows_;
    ZhangExactVector acceptedValues_;
};
