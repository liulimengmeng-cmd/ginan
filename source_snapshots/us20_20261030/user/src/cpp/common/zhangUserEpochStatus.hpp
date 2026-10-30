#ifndef ZHANG_USER_EPOCH_STATUS_HPP
#define ZHANG_USER_EPOCH_STATUS_HPP

#include "common/algebra.hpp"

// Sideband availability, not a substitute for the authoritative KF transaction.
inline constexpr const char* ZHANG_USER_LAST_UPDATE = "ZHANG_USER_LAST_MEASUREMENT_UPDATE";
inline constexpr const char* ZHANG_USER_EPOCH_STATUS = "ZHANG_USER_EPOCH_STATUS";

inline bool zhangUserHasMeasurementHistory(const KFState& state)
{
    const auto it = state.metaDataMap.find(ZHANG_USER_LAST_UPDATE);
    return (it != state.metaDataMap.end() && !it->second.empty()) ||
        (state.lastPppTransactionCommitted && state.lastPppQcConverged);
}

inline bool zhangUserMayContinueEmptyEpoch(
    bool standaloneUser, Eigen::Index rows, const KFState& predicted)
{
    return standaloneUser && rows == 0 && predicted.x.size() > 0 &&
        predicted.P.rows() == predicted.x.size() &&
        predicted.P.cols() == predicted.x.size() &&
        predicted.x.allFinite() && predicted.P.allFinite();
}

inline bool zhangUserMayRunEpochAr(const KFState& state)
{
    return state.lastPppTransactionCommitted && state.lastPppQcConverged;
}

inline void zhangUserMarkEmptyEpoch(KFState& state)
{
    if (state.lastPppTransactionCommitted && state.lastPppQcConverged &&
        !state.metaDataMap.contains(ZHANG_USER_LAST_UPDATE))
        state.metaDataMap[ZHANG_USER_LAST_UPDATE] = state.time.to_string(0);
    state.lastPppTransactionCommitted = false;
    state.lastPppQcConverged = false;
    state.metaDataMap[ZHANG_USER_EPOCH_STATUS] = zhangUserHasMeasurementHistory(state)
        ? "PREDICTION_ONLY" : "UNINITIALIZED";
}

#endif
