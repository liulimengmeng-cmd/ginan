#pragma once
#include <array>
#include <cstdlib>
#include <string>
#include <vector>
#include <boost/serialization/version.hpp>

inline bool zhangEvidenceSlipQcEnabled()
{
    const char* value=std::getenv("ZHANG_R51_EVIDENCE_SLIP_QC");
    return value && std::string(value)=="1";
}

// Local receiver-link QC, NOT network AR/product joint NIS. Combination units
// are [metres, metres, cycles, metres] for IF, GF, MW and the triple combination.
struct ZhangSlipSample
{
    double time=0, rawIf=0, clock=0;
    int signature=0, thirdSignature=0;
    bool geometry=false, third=false, code=false;
    std::array<double,4> value{};
    std::array<double,16> covariance{};
    std::array<double,3> wavelength{};
    double gfRate=0; // -lambda1*D1 + lambda2*D2, metres/second
    bool gfRateValid=false;
    template<class A> void serialize(A& a, const unsigned int& version)
    { a & time & rawIf & clock & signature & thirdSignature & geometry & third & code
        & value & covariance & wavelength;
      if(version>=1)a & gfRate & gfRateValid; }
};
BOOST_CLASS_VERSION(ZhangSlipSample, 1)

struct ZhangSlipEvidenceState
{
    std::vector<ZhangSlipSample> accepted; // bounded at 32, no suspect learning
    // Screening is not trusted prediction history or an integer certificate.
    std::vector<ZhangSlipSample> screening; // bounded at 96, includes suspects
    bool phaseSuspect=false;
    template<class A> void serialize(A& a, const unsigned int&) { a & accepted & screening & phaseSuspect; }
};

enum class ZhangSlipAction { ACCEPT, IONO_TREND, WARMUP, CODE_OUTLIER,
    PHASE_SUSPECT, INTEGER_STEP, UNRESOLVED_STEP, NO_EVIDENCE };

inline const char* zhangSlipActionName(ZhangSlipAction action)
{
    switch(action)
    {
        case ZhangSlipAction::ACCEPT:return "ACCEPT";
        case ZhangSlipAction::IONO_TREND:return "IONO_TREND";
        case ZhangSlipAction::WARMUP:return "WARMUP";
        case ZhangSlipAction::CODE_OUTLIER:return "CODE_OUTLIER";
        case ZhangSlipAction::PHASE_SUSPECT:return "PHASE_SUSPECT";
        case ZhangSlipAction::INTEGER_STEP:return "INTEGER_STEP";
        case ZhangSlipAction::UNRESOLVED_STEP:return "CONTINUITY_UNKNOWN_PHASE_STEP";
        default:return "NO_EVIDENCE";
    }
}

struct ZhangSlipEvidenceDecision
{
    ZhangSlipAction action=ZhangSlipAction::NO_EVIDENCE;
    bool usePhase=false, useCode=true, breakArc=false;
    std::array<int,3> integer{};
    std::array<bool,3> retire{};
    std::array<double,4> innovation{}, sigma{};
    double phaseScore=0, bestScore=0, ratio=0, changePointScore=0, changePointTime=0;
    bool dopplerAvailable=false,dopplerAssisted=false;
    double dopplerInnovation=0,dopplerSigma=0,dopplerScore=0,dopplerIntegerFit=0;
    int clockPartners=0;
};
