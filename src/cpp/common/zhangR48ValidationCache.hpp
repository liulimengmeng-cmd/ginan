#pragma once
#include "common/zhangIntegerCandidateNis.hpp"
#include "common/zhangPhaseTimer.hpp"

// Owned by ONE immutable FLOAT posterior call, never by epoch number or
// matrix dimension alone. Only numerical validation is cached; proof/risk
// closure and exact integer feasibility remain the caller's responsibility.
class ZhangR48ValidationCache
{
    const Eigen::VectorXd& mean;
    const Eigen::MatrixXd& covariance;
    struct Entry {
        Eigen::MatrixXd rows;
        Eigen::VectorXd values;
        double alpha;
        ZhangIntegerCandidateNis assessment;
    };
    std::vector<Entry> entries;
public:
    int hits=0,evaluations=0;
    ZhangR48ValidationCache(const Eigen::VectorXd& m,const Eigen::MatrixXd& p)
        :mean(m),covariance(p){}
    ZhangIntegerCandidateNis assess(const Eigen::MatrixXd& rows,
        const Eigen::VectorXd& values,double alpha)
    {
        ZhangPhaseTimer timer("R48_ROUTE_NUMERICAL_VALIDATION");
        for(const auto& entry:entries)
            if(entry.alpha==alpha && entry.rows.rows()==rows.rows() &&
               entry.rows.cols()==rows.cols() && entry.values.size()==values.size() &&
               (entry.rows.array()==rows.array()).all() &&
               (entry.values.array()==values.array()).all())
            {++hits;return entry.assessment;}
        ++evaluations;
        const auto assessment=assessZhangIntegerCandidateNis(values-rows*mean,
            rows*covariance*rows.transpose(),alpha);
        // Two competing routes only. Rejected trial systems cannot grow a
        // historical cache; RHS changes must re-evaluate null-space checks.
        if(entries.size()>=2)entries.erase(entries.begin());
        entries.push_back({rows,values,alpha,assessment});
        return assessment;
    }
    ZhangR48ValidationCache(const ZhangR48ValidationCache&)=delete;
};
