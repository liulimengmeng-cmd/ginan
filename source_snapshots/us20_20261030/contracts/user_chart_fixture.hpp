#pragma once
#include "common/algebra.hpp"
#include "common/receiver.hpp"
#include "pea/zhangPppAr.hpp"

// Fixture construction must follow the actual order: establish an owner on
// an empty filter, then create observation-backed phase states. Existing
// moments below are the synthetic observation result, not unowned history.
inline void initializeFixtureChart(Trace& trace,ReceiverMap& receivers,KFState& state)
{
    // Keep the actual bound KF object address: a copied object is a different
    // checkpoint runtime owner and is correctly rejected by the real hook.
    const auto rts=state.rts_basename;
    const auto ownerOverride=state.factorCallbackOwnerOverride;
    KFState measured=state;
    state.kfIndexMap.clear();state.stateTransitionMap.clear();
    state.kfIndexMap[state.oneKey]=0;
    state.stateTransitionMap[state.oneKey][state.oneKey][0]=1;
    state.x=VectorXd::Ones(1);state.dx=VectorXd::Zero(1);state.P=MatrixXd::Zero(1,1);
    updateZhangPppArUserReferences(trace,receivers,state);
    auto metadata=std::move(state.metaDataMap);
    state=std::move(measured);state.metaDataMap=std::move(metadata);
    state.rts_basename=rts;state.factorCallbackOwnerOverride=ownerOverride;
    persistZhangUserAppliedChartReceipts(state);
}
