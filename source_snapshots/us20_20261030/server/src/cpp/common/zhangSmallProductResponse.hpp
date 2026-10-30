#pragma once

#include "common/zhangProductConditioningReceipt.hpp"
#include "common/zhangExactWorkBudget.hpp"

/** Exact-set recipe difference. Every row actually applied to the cached
 * posterior must remain in the current recipe with the SAME RHS. This is not
 * rank equivalence, proportionality, or an approximate coefficient match. */
inline bool zhangProductResponseExtraConditions(
    const std::map<ZhangCoordinateRow,ZhangExactInteger>& applied,
    const std::vector<ZhangCurrentEpochNumericConstraint>& current,
    std::vector<ZhangCurrentEpochNumericConstraint>& extras,std::string& reason)
{
    extras.clear();
    std::vector<ZhangCurrentEpochNumericConstraint> normalized;
    if(!zhangNormalizeProductNumericConstraints(current,normalized,reason))return false;
    std::map<ZhangCoordinateRow,ZhangExactInteger> rows;
    for(const auto& row:normalized) {
        zhangExactWorkCheckpoint("PRODUCT_RESPONSE_RECIPE_CURRENT");
        rows.emplace(row.rowCycles,row.value);
    }
    for(const auto& [row,rhs]:applied) {
        zhangExactWorkCheckpoint("PRODUCT_RESPONSE_RECIPE_BASE");
        const auto found=rows.find(row);
        if(found==rows.end() || found->second!=rhs) {
            reason="PRODUCT_RESPONSE_APPLIED_ROW_OR_RHS_NOT_CURRENT";return false;
        }
        rows.erase(found);
    }
    if(rows.size()>64) {reason="PRODUCT_RESPONSE_EXTRA_CONDITION_CAP";return false;}
    for(const auto& [row,rhs]:rows)extras.push_back({row,rhs});
    reason="PRODUCT_RESPONSE_COMPLETE_BASE_RECIPE_PRESENT";return true;
}

struct ZhangSmallProductResponse
{
    bool valid=false;
    std::string reason="NOT_EVALUATED";
    Eigen::VectorXd mean,extraValues;
    Eigen::MatrixXd covariance,targetProjector,extraSelector;
    std::uint64_t covarianceTerms=0;
};

/** Form U*x0 and U*P0*U^T directly from a VERIFIED immutable, already
 * conditioned full-state posterior. U contains the <=2 target rows followed
 * by <=64 newly admitted conditions. No cross covariance is omitted and no
 * thousands-dimensional submatrix is copied or factorised. */
inline ZhangSmallProductResponse zhangBuildSmallProductResponse(
    const KFState& base,const std::vector<ZhangCoordinateRow>& targets,
    const std::vector<ZhangCurrentEpochNumericConstraint>& extras)
{
    ZhangSmallProductResponse out;
    const int r=targets.size(),k=extras.size(),n=r+k;
    if(r<1 || r>2 || k>64 || base.P.rows()!=base.x.size() ||
        base.P.cols()!=base.x.size()) {out.reason="PRODUCT_RESPONSE_DIMENSION_INVALID";return out;}
    std::vector<std::vector<std::pair<int,double>>> rows(n);
    auto mapRow=[&](int row,const ZhangCoordinateRow& coordinates) {
        for(const auto& [key,value]:coordinates) {
            zhangExactWorkCheckpoint("PRODUCT_RESPONSE_SPARSE_MAP");
            if(!std::isfinite(value))return false;
            if(value==0)continue;
            const auto found=base.kfIndexMap.find(key);
            if(found==base.kfIndexMap.end() || found->second<0 || found->second>=base.x.size())return false;
            rows[row].emplace_back(found->second,value);
        }
        return true;
    };
    for(int i=0;i<r;++i)if(!mapRow(i,targets[i])) {out.reason="PRODUCT_RESPONSE_TARGET_KEY_INVALID";return out;}
    out.extraValues.resize(k);
    for(int i=0;i<k;++i) {
        if(!mapRow(r+i,extras[i].rowCycles) || extras[i].value < -9007199254740991LL ||
            extras[i].value > 9007199254740991LL) {out.reason="PRODUCT_RESPONSE_EXTRA_ROW_OR_RHS_INVALID";return out;}
        out.extraValues[i]=extras[i].value.convert_to<double>();
    }
    out.mean=Eigen::VectorXd::Zero(n);out.covariance=Eigen::MatrixXd::Zero(n,n);
    for(int i=0;i<n;++i) {
        long double value=0;
        for(const auto& [index,coefficient]:rows[i]) {
            zhangExactWorkCheckpoint("PRODUCT_RESPONSE_MEAN");
            if(!std::isfinite(base.x[index])) {out.reason="PRODUCT_RESPONSE_NONFINITE_MEAN";return out;}
            value+=static_cast<long double>(coefficient)*base.x[index];
        }
        out.mean[i]=static_cast<double>(value);
        for(int j=i;j<n;++j) {
            long double covariance=0;
            for(const auto& [a,ca]:rows[i])for(const auto& [b,cb]:rows[j]) {
                zhangExactWorkCheckpoint("PRODUCT_RESPONSE_COMPLETE_CROSS_COVARIANCE");
                ++out.covarianceTerms;
                const double pab=base.P(a,b),pba=base.P(b,a);
                if(!std::isfinite(pab) || !std::isfinite(pba) ||
                    std::abs(pab-pba)>1e-10*std::max({1.0,std::abs(pab),std::abs(pba)})) {
                    out.reason="PRODUCT_RESPONSE_NONFINITE_OR_ASYMMETRIC_COVARIANCE";return out;
                }
                covariance+=static_cast<long double>(ca)*cb*(static_cast<long double>(pab)+pba)*0.5L;
            }
            out.covariance(i,j)=out.covariance(j,i)=static_cast<double>(covariance);
        }
    }
    if(!out.mean.allFinite() || !out.covariance.allFinite()) {out.reason="PRODUCT_RESPONSE_NONFINITE_PROJECTION";return out;}
    out.targetProjector=Eigen::MatrixXd::Zero(r,n);
    out.extraSelector=Eigen::MatrixXd::Zero(k,n);
    for(int i=0;i<r;++i)out.targetProjector(i,i)=1;
    for(int i=0;i<k;++i)out.extraSelector(i,r+i)=1;
    out.valid=true;out.reason="PRODUCT_RESPONSE_FULL_CROSS_COVARIANCE";return out;
}
