#pragma once
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

// R50 acceptance ablation only. Raw probabilities/NIS remain diagnostics;
// ratio-only receipts MUST NOT be interpreted as a calibrated risk bound.
inline bool zhangR50RatioOnly()
{
    const char* value=std::getenv("ZHANG_R50_RATIO_ONLY");
    return value && std::strcmp(value,"1")==0;
}
inline bool zhangR50StatisticalAccept(bool originalDecision)
{ return zhangR50RatioOnly() || originalDecision; }
inline bool zhangR50StatisticalReject(bool originalDecision)
{ return !zhangR50RatioOnly() && originalDecision; }

struct ZhangR50RatioResult
{
    bool valid=false, accepted=false;
    double ratio=std::numeric_limits<double>::quiet_NaN();
};
inline ZhangR50RatioResult zhangR50AssessRatio(double best,double second,double threshold)
{
    ZhangR50RatioResult out;
    if(!std::isfinite(best) || !std::isfinite(second) || best<0 || second<best ||
       !std::isfinite(threshold) || threshold<=1) return out;
    out.valid=true;
    out.ratio=best==0 ? (second>0 ? std::numeric_limits<double>::infinity() : 1) : second/best;
    out.accepted=out.ratio>=threshold;
    return out;
}
