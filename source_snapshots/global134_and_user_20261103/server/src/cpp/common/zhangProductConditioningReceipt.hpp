#pragma once

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "common/zhangCurrentEpochEvidence.hpp"

// A damaged authoritative FLOAT is not an eligible fallback. Latch the
// violation so that a caught exception cannot make a later commit eligible.
inline void zhangRequireUnchangedFloatAuthority(bool coreEqual,bool& violation)
{
    violation |= !coreEqual;
    if(violation)throw std::runtime_error("ZHANG_FLOAT_AUTHORITY_CHANGED_NO_PUBLICATION_OR_CHECKPOINT");
}

/** Canonical equality of the complete numerical conditioning recipe. This
 * does NOT replace exact integer feasibility, proof closure, or a successful
 * conditioner. No rank-only, gcd/saturation, or approximate row comparison is
 * allowed. Floating coefficients must match exactly after removing zeros. */
inline bool zhangNormalizeProductNumericConstraints(
    const std::vector<ZhangCurrentEpochNumericConstraint>& input,
    std::vector<ZhangCurrentEpochNumericConstraint>& output,
    std::string& reason)
{
    output.clear();
    std::map<ZhangCoordinateRow,ZhangExactInteger> rows;
    for(const auto& constraint:input) {
        ZhangCoordinateRow normalized;
        for(const auto& [key,coefficient]:constraint.rowCycles) {
            if(!std::isfinite(coefficient)) {reason="CONDITIONING_RECEIPT_NONFINITE_ROW";return false;}
            if(coefficient!=0)normalized.emplace(key,coefficient);
        }
        if(normalized.empty()) {
            if(constraint.value!=0) {reason="CONDITIONING_RECEIPT_ZERO_ROW_CONFLICT";return false;}
            continue; // Structural identity: not a numerical fixed row.
        }
        const auto [prior,inserted]=rows.emplace(std::move(normalized),constraint.value);
        if(!inserted && prior->second!=constraint.value) {
            reason="CONDITIONING_RECEIPT_DUPLICATE_RHS_CONFLICT";return false;
        }
    }
    for(auto& [row,value]:rows)output.push_back({row,value});
    reason="NONE";return true;
}

struct ZhangProductConditioningIdentity
{
    std::string runtimeId;
    long int epochGpstSeconds=0;
    std::string posteriorId,chartId;
    std::string sourceRoot,sourceStateOrder;
    std::string resultRoot,resultStateOrder;
    bool operator==(const ZhangProductConditioningIdentity&) const = default;
};

/** Issued ONLY after a successful full-state numerical conditioner. The
 * caller must supply roots computed from the actual protected FLOAT and the
 * actual successful private result. The writer checks those actual objects
 * again after every current-epoch event admission. This receipt is ephemeral
 * and cannot carry integer authority into another epoch. */
struct ZhangProductConditioningReceipt
{
    ZhangProductConditioningIdentity identity;
    std::vector<ZhangCurrentEpochNumericConstraint> constraints;
    bool valid=false;

    bool bind(const ZhangProductConditioningIdentity& binding,
        const std::vector<ZhangCurrentEpochNumericConstraint>& applied,
        std::string& reason)
    {
        valid=false;constraints.clear();identity={};
        if(binding.runtimeId.empty() || binding.epochGpstSeconds<=0 ||
            binding.posteriorId.empty() || binding.chartId.empty() ||
            binding.sourceRoot.empty() || binding.resultRoot.empty() ||
            binding.sourceStateOrder.empty() ||
            binding.sourceStateOrder!=binding.resultStateOrder) {
            reason="CONDITIONING_RECEIPT_INCOMPLETE_IDENTITY";return false;
        }
        if(!zhangNormalizeProductNumericConstraints(applied,constraints,reason) || constraints.empty()) {
            if(reason=="NONE")reason="CONDITIONING_RECEIPT_EMPTY_RECIPE";
            return false;
        }
        identity=binding;valid=true;reason="CONDITIONING_RECEIPT_ISSUED";return true;
    }

    bool matches(const ZhangProductConditioningIdentity& actual,
        const std::vector<ZhangCurrentEpochNumericConstraint>& required,
        std::string& reason) const
    {
        if(!valid || !(identity==actual)) {
            reason="CONDITIONING_RECEIPT_ROOT_ORDER_SCOPE_MISMATCH";return false;
        }
        std::vector<ZhangCurrentEpochNumericConstraint> normalized;
        if(!zhangNormalizeProductNumericConstraints(required,normalized,reason))return false;
        if(normalized.size()!=constraints.size()) {
            reason="CONDITIONING_RECEIPT_RECIPE_SIZE_MISMATCH";return false;
        }
        for(std::size_t row=0;row<normalized.size();++row)
            if(normalized[row].rowCycles!=constraints[row].rowCycles ||
                normalized[row].value!=constraints[row].value) {
                reason="CONDITIONING_RECEIPT_ROW_OR_RHS_MISMATCH";return false;
            }
        reason="CONDITIONING_RECEIPT_COMPLETE_RECIPE_MATCH";return true;
    }
};

inline std::map<std::string,ZhangProductConditioningReceipt>&
zhangProductConditioningReceiptRegistry()
{
    static std::map<std::string,ZhangProductConditioningReceipt> receipts;
    return receipts;
}
