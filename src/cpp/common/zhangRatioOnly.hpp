#pragma once
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

// R51 acceptance ablation only. Integer/product NIS is not evaluated in this
// mode; ratio-only receipts MUST NOT be interpreted as calibrated risk bounds.
inline bool zhangRatioOnly()
{
    const char* value=std::getenv("ZHANG_R51_RATIO_ONLY");
    return value && std::strcmp(value,"1")==0;
}
// Product runs avoid diagnostic-only decompositions by default. This switch
// must never bypass a solver, integer-validity check, or writer authorization.
inline bool zhangExpensiveDiagnosticsEnabled()
{
    if (zhangRatioOnly()) return false;
    const char* value=std::getenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS");
    return value ? std::strcmp(value,"1")==0 : true;
}
inline bool zhangRatioStatisticalAccept(bool originalDecision)
{ return zhangRatioOnly() || originalDecision; }
inline bool zhangRatioStatisticalReject(bool originalDecision)
{ return !zhangRatioOnly() && originalDecision; }

// In ratio-only mode lambda_candidate_nis_valid records that LAMBDA returned
// a finite candidate distance; it does not certify an evaluated NIS.
template<class Candidate>
bool zhangLambdaCandidateAdmissible(const Candidate& candidate)
{
    return candidate.lambda_candidate_nis_valid &&
        (zhangRatioOnly() || candidate.lambda_candidate_nis <=
            candidate.lambda_candidate_nis_threshold);
}

// Do not construct a diagnostic covariance or run its eigensolver when the
// statistical NIS decision is disabled.  The caller must keep its separate
// exact-integer and posterior-conditioning validity checks.  A skipped result
// is deliberately invalid/NaN, never a fabricated NIS pass.
template<class Compute>
auto zhangRatioOnlySkipNis(Compute&& compute)
{
    if (!zhangRatioOnly()) return compute();
    auto skipped = decltype(compute()){};
    skipped.status = "SKIPPED_RATIO_ONLY";
    return skipped;
}

struct ZhangR50RatioResult
{
    bool valid=false, accepted=false;
    double ratio=std::numeric_limits<double>::quiet_NaN();
};
inline ZhangR50RatioResult zhangRatioAssessRatio(double best,double second,double threshold)
{
    ZhangR50RatioResult out;
    if(!std::isfinite(best) || !std::isfinite(second) || best<0 || second<best ||
       !std::isfinite(threshold) || threshold<=1) return out;
    out.valid=true;
    out.ratio=best==0 ? (second>0 ? std::numeric_limits<double>::infinity() : 1) : second/best;
    out.accepted=out.ratio>=threshold;
    return out;
}
