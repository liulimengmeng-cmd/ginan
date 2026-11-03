#include "common/zhangUserPhaseQc.hpp"
#include "common/zhangPhaseArcDecision.hpp"

#include <iostream>
#include <limits>

static unsigned assertions = 0;
#define CHECK(condition) do { ++assertions; if (!(condition)) { \
    std::cerr << "FAIL line " << __LINE__ << ": " #condition "\n"; return EXIT_FAILURE; \
} } while (false)

static ZhangUserTrendDecision inspect(ZhangUserQcTrend& state, double time, double gf,
                                     bool reset=false, double sigma=.004)
{
    return zhangInspectUserTrend(state, time, gf, sigma, .05, 4, 1e-6, 90, reset);
}

int main()
{
    // Controlled smooth ionosphere, including a >0.05m raw inter-epoch rate.
    // Predictor uses past accepted samples only, not future observations.
    ZhangUserQcTrend linear;
    for (int i=0; i<40; ++i)
    {
        const double t=i*30;
        const auto d=inspect(linear,t,2+.0019*t);
        CHECK(d.accepted);
        if (i>=3) CHECK(std::abs(d.innovation)<1e-10);
    }
    CHECK(linear.accepted.size()==8);
    CHECK(linear.consecutiveSuspects==0);
    ZhangUserQcTrend quadratic;
    for (int i=0; i<80; ++i)
    {
        const double t=i*30;
        const auto d=inspect(quadratic,t,2+.0007*t+.00000022*t*t);
        CHECK(d.accepted);
        if (i>=5) CHECK(std::abs(d.innovation)<1e-8);
    }

    // Warmup is explicit, and gross acquisition spikes cannot enter history.
    ZhangUserQcTrend warmup;
    CHECK(inspect(warmup,0,1).warmup);
    CHECK(!inspect(warmup,30,6).accepted);
    CHECK(warmup.accepted.size()==1);
    CHECK(inspect(warmup,60,1.02).accepted);

    // Confirmed physical evidence has a per-band serial. Quarantine alone
    // never increments it or writes a guessed new physical integer identity.
    ZhangSignalTracking tracking;
    auto tr=zhangDecidePhaseTracking(tracking,0,true,1,0,45);
    CHECK(!tr.breakArc && tr.useCurrentPhase && tracking.eventSequence==0);
    tr=zhangDecidePhaseTracking(tracking,1,true,1,30,45);
    CHECK(tr.breakArc && tr.useCurrentPhase && tracking.eventSequence==1);
    tr=zhangDecidePhaseTracking(tracking,0,true,1,60,45);
    CHECK(!tr.breakArc && tracking.eventSequence==1);
    tr=zhangDecidePhaseTracking(tracking,2,true,1,90,45);
    CHECK(tr.breakArc && !tr.useCurrentPhase && tracking.eventSequence==2);
    tr=zhangDecidePhaseTracking(tracking,2,true,1,120,45);
    CHECK(!tr.breakArc && !tr.useCurrentPhase && tracking.eventSequence==2);
    tr=zhangDecidePhaseTracking(tracking,0,true,1,150,45);
    CHECK(tr.breakArc && tr.useCurrentPhase && tracking.eventSequence==3);
    tr=zhangDecidePhaseTracking(tracking,0,false,1,180,45);
    CHECK(!tr.breakArc && !tr.useCurrentPhase && tracking.eventSequence==3);
    tr=zhangDecidePhaseTracking(tracking,0,true,1,210,45);
    CHECK(tr.breakArc && tr.gap && tracking.eventSequence==4);
    tr=zhangDecidePhaseTracking(tracking,0,true,2,240,45);
    CHECK(tr.breakArc && tr.gap && tracking.eventSequence==5);
    tr=zhangDecidePhaseTracking(tracking,0,true,2,230,45);
    CHECK(tr.breakArc && tr.gap && tracking.eventSequence==6);

    // L1-only and equal-frequency-count jumps produce GF evidence after
    // warmup; neither is labelled a proven physical slip without tracking.
    constexpr double l1=.190293672798365, l2=.244210213424568;
    for (const auto jump : {l1, -l2, l1-l2})
    {
        ZhangUserQcTrend trend;
        for (int i=0; i<12; ++i) CHECK(inspect(trend,i*30,1).accepted);
        const auto count=trend.accepted.size();
        auto d=inspect(trend,360,1+jump);
        CHECK(!d.accepted);
        CHECK(trend.accepted.size()==count);
        auto e=zhangUserCombinationEligibility(false,false,true,d.accepted,true,true);
        CHECK(!e.usePhase && e.useCode && e.state==ZhangUserQcState::SUSPECT);
        d=inspect(trend,390,1+jump);
        CHECK(!d.accepted && trend.consecutiveSuspects==2);
        CHECK(!inspect(trend,420,1+jump).accepted);
        CHECK(!inspect(trend,450,1+jump).accepted);
        CHECK(!inspect(trend,480,1+jump).accepted);
        CHECK(trend.accepted.size()==count);
        CHECK(inspect(trend,510,1,true).accepted);
    }
    // (9,7) is almost GF-null. At realistic code noise its 2-cycle MW change
    // can also be below a four-sigma gate: the test documents this blind case,
    // rather than claiming that GF+MW guarantees complete slip detection.
    CHECK(std::abs(9*l1-7*l2)<.005);
    ZhangUserQcTrend blindGf,blindMw;
    for(int i=0;i<12;++i)
    {
        CHECK(inspect(blindGf,i*30,1).accepted);
        CHECK(zhangInspectUserTrend(blindMw,i*30,0,.6,1,4,1e-4,90).accepted);
    }
    const auto bGF=inspect(blindGf,360,1+9*l1-7*l2);
    const auto bMW=zhangInspectUserTrend(blindMw,360,2,.6,1,4,1e-4,90);
    CHECK(bGF.accepted && bMW.accepted);
    CHECK(zhangUserCombinationEligibility(true,false,true,true,true,true).state==
          ZhangUserQcState::CONFIRMED_BREAK); // independent LLI still retires

    // MW cannot identify whether code or phase caused a spike: without
    // independent code evidence it preserves good pseudorange.
    auto e=zhangUserCombinationEligibility(false,false,true,true,true,false);
    CHECK(e.useCode && !e.usePhase && e.state==ZhangUserQcState::SUSPECT);
    e=zhangUserCombinationEligibility(false,false,true,false,true,true);
    CHECK(e.useCode && !e.usePhase && e.state==ZhangUserQcState::SUSPECT);
    e=zhangUserCombinationEligibility(false,false,true,true,false,false);
    CHECK(e.useCode && e.usePhase); // missing code must not disable GF
    e=zhangUserCombinationEligibility(false,true,false,false,false,false);
    CHECK(e.useCode && !e.usePhase);
    e=zhangUserCombinationEligibility(true,false,true,true,true,true);
    CHECK(e.useCode && e.usePhase && e.state==ZhangUserQcState::CONFIRMED_BREAK);

    // Controlled code-only spike: the code GF cancels range/clock, independently
    // identifies a single-band 20m code anomaly, and does not retire an integer.
    ZhangUserQcTrend codeGf;
    for (int i=0;i<12;++i)
        CHECK(zhangInspectUserTrend(codeGf,i*30,2+.002*i*30,.425,3,4,4e-6,90).accepted);
    auto code=zhangInspectUserTrend(codeGf,360,2+.002*360+20,.425,3,4,4e-6,90);
    CHECK(!code.accepted);
    e=zhangUserCombinationEligibility(false,false,true,true,true,false,!code.accepted);
    CHECK(!e.useCode && !e.usePhase && e.state==ZhangUserQcState::SUSPECT);
    // A carrier (9,7) slip changes MW but not code GF. Good code remains usable.
    ZhangUserQcTrend unchangedCode;
    for (int i=0;i<12;++i)
        CHECK(zhangInspectUserTrend(unchangedCode,i*30,2,.425,3,4,4e-6,90).accepted);
    code=zhangInspectUserTrend(unchangedCode,360,2,.425,3,4,4e-6,90);
    CHECK(code.accepted);
    e=zhangUserCombinationEligibility(false,false,true,true,true,false,!code.accepted);
    CHECK(e.useCode && !e.usePhase);

    // Small sampling irregularities are handled by actual dt, not an assumed
    // 30-second grid. The first sample after a confirmed gap is a new baseline.
    ZhangUserQcTrend irregular;
    for (double t : {0.,29.,61.,89.,122.,150.,181.,211.,242.})
        CHECK(inspect(irregular,t,1+.0017*t).accepted);
    auto d=inspect(irregular,500,4);
    CHECK(d.accepted && d.gap && d.warmup && irregular.accepted.size()==1);
    const auto count=irregular.accepted.size();
    CHECK(!inspect(irregular,530,std::numeric_limits<double>::quiet_NaN()).accepted);
    CHECK(!inspect(irregular,530,4,false,0).accepted);
    CHECK(irregular.accepted.size()==count);

    // A large absolute GPST origin must not destroy local trend numerics.
    ZhangUserQcTrend gpsTime;
    for (int i=0;i<50;++i)
        CHECK(inspect(gpsTime,1405216800.+i*30,5+.0016*i*30).accepted);
    std::cout << "USER_PHASE_QC_CONTRACT PASS assertions=" << assertions
              << " nis_calculated=0 true_slip_probability=UNASSESSED\n";
    return EXIT_SUCCESS;
}
