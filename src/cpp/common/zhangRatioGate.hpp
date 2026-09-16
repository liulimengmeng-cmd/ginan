#pragma once
#include <cmath>
#include <limits>
#include <string>
struct ZhangRatioGateResult {
    bool executed=false, accepted=false;
    double ratio=std::numeric_limits<double>::quiet_NaN();
    std::string reason;
};
inline ZhangRatioGateResult zhangEvaluateRatioGate(bool complete, bool twoDistinct,
    double first, double second, double threshold)
{
    ZhangRatioGateResult r;
    if (!complete) { r.reason="ILS_SEARCH_INCOMPLETE"; return r; }
    if (!twoDistinct) { r.reason="TWO_DISTINCT_CANDIDATES_REQUIRED"; return r; }
    if (!std::isfinite(first) || !std::isfinite(second) || first<0 || second<first ||
        !std::isfinite(threshold) || threshold<=1)
    { r.reason="INVALID_RATIO_INPUT"; return r; }
    r.executed=true;
    r.ratio=first>0 ? second/first : (second>0 ? std::numeric_limits<double>::infinity() : 1);
    r.accepted=r.ratio>=threshold;
    r.reason=r.accepted ? "ACCEPTED" : "RATIO_REJECTED";
    return r;
}
