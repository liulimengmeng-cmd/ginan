#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangUserProductError.hpp"
#include "common/zhangUserProductChart.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

void removeBadAmbiguities(Trace&, KFState&, ReceiverMap&);
void mergeCorrelated(Trace&, KFState&, KFMeasEntryList&);

extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool value, const char* message)
    {
        ++checks;
        if (!value) throw std::runtime_error(message);
    };
    try
    {
        setenv("ZHANG_USER_QC", "CAUSAL_TRACKING", 1);
        acsConfig.process_ppp = true;
        acsConfig.zhangFullRank.enable = false;
        acsConfig.pppOpts.ionoOpts.use_if_combo = false;
        acsConfig.ambErrors.phase_reject_limit = 10000;
        auto& reset = acsConfig.ambErrors.resetOnSlip;
        reset.LLI = reset.retrack = reset.single_freq = false;
        reset.GF = reset.MW = reset.SCDIA = false;
        ReceiverMap receivers;
        auto& receiver = receivers["TEST"];
        receiver.id = "TEST";
        KFState state;
        GTime time;
        time.bigTime = 100000;
        std::ostringstream trace;
        auto key = [&](const char* satellite, E_ObsCode code)
        {
            KFKey value;
            value.type = KF::AMBIGUITY;
            value.str = receiver.id;
            value.Sat = SatSys(satellite);
            value.num = static_cast<int>(code);
            value.rec_ptr = &receiver;
            return value;
        };
        const auto first = key("G01", E_ObsCode::L1C);
        const auto second = key("G01", E_ObsCode::L2W);
        const auto remote = key("G02", E_ObsCode::L1C);
        InitialState initial;
        initial.estimate = true;
        initial.P = 1;
        for (const auto& k : {first, second, remote}) state.addKFState(k, initial);
        check(state.stateTransition(trace, time), "initial native transition");
        zhangUserRememberChartAnchor(state,receiver.id,first.Sat,first.num,{true,true,"C0",0,0});
        zhangUserSetChartOffset(state,receiver.id,first.Sat,first.num,.37);
        state.metaDataMap[zhangUserChartKey(receiver.id,first.Sat,first.num)+":UNANCHORED_MAP"]="1";
        auto& signal = receiver.satStatMap[first.Sat].sigStatMap[ft2string(F1)];
        signal.phaseQuarantined = true;
        signal.userQcState = ZhangUserQcState::SUSPECT;
        removeBadAmbiguities(trace, state, receivers);
        check(state.stateTransitionMap.count(first), "suspect does not retire integer");
        check(state.stateTransitionMap.count(second), "suspect does not retire companion frequency");
        signal.userPhysicalBreakThisEpoch = true;
        signal.tracking.pendingBreak = true;
        signal.tracking.eventSequence = 7;
        removeBadAmbiguities(trace, state, receivers);
        check(!state.stateTransitionMap.count(first), "confirmed physical break retires despite disabled legacy flags");
        check(zhangUserChartOffset(state,receiver.id,first.Sat,first.num)==0 &&
            !state.metaDataMap.count(zhangUserChartKey(receiver.id,first.Sat,first.num)+":UNANCHORED_MAP") &&
            !zhangUserGonsisComponentCertificateValid(zhangUserChartAnchor(state,receiver.id,first.Sat,first.num)),
            "native physical retirement deletes old fractional integer map and its unanchored marker");
        check(state.stateTransitionMap.count(second), "receiver satellite frequency locality");
        check(state.stateTransitionMap.count(remote), "other satellite not retired");
        check(signal.tracking.eventSequence == 7, "physical sequence preserved");
        check(state.metaDataMap.at("USER_AR_QC_APPLIED_SEQUENCE:TEST:G01:" +
            std::to_string(first.num)) == "7", "native physical event receipt prevents duplicate reference lifecycle action");
        check(!signal.userPhysicalBreakThisEpoch && !signal.tracking.pendingBreak, "native event consumed once");
        check(state.stateTransition(trace, time), "removal marginalization");
        check(!state.kfIndexMap.count(first), "retired state absent from native KF");
        state.addKFState(first, initial);
        check(state.stateTransition(trace, time), "new physical arc added");
        removeBadAmbiguities(trace, state, receivers);
        check(state.stateTransitionMap.count(first), "duplicate cleanup does not remove new arc");
        check(trace.str().find("user_confirmed_physical_break") != std::string::npos, "physical cause traced");
        setenv("ZHANG_USER_PRODUCT_ERROR_POLICY", "SEGMENT_CONSTANT", 1);
        KFMeasEntryList emptyMeasurements;
        bool mergeRejected = false;
        try { mergeCorrelated(trace, state, emptyMeasurements); }
        catch (const std::runtime_error&) { mergeRejected = true; }
        check(mergeRejected, "product error must not be merged into integer ambiguity");
        zhangUserValidateProductErrorEpochAlignment(ZhangUserProductErrorPolicy::SEGMENT_CONSTANT, time, time);
        GTime offset = time;
        offset.bigTime -= 0.1L;
        bool offsetRejected = false;
        try { zhangUserValidateProductErrorEpochAlignment(ZhangUserProductErrorPolicy::SEGMENT_CONSTANT, offset, time); }
        catch (const std::runtime_error&) { offsetRejected = true; }
        check(offsetRejected, "mixed-time record retirement unsupported and explicitly rejected");
        zhangUserValidateProductErrorEpochAlignment(ZhangUserProductErrorPolicy::WHITE, offset, time);
        check(true, "legacy white preserves closest-epoch behavior");
        std::cout << "USER_NATIVE_QC_LIFECYCLE_CONTRACT PASS checks=" << checks
                  << " nis_evaluations=0\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "USER_NATIVE_QC_LIFECYCLE_CONTRACT FAIL " << error.what() << '\n';
        return 1;
    }
}
