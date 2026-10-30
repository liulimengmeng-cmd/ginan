#pragma once
#include "common/algebra.hpp"
#include "common/zhangRegionalJointProduct.hpp"

inline const std::string zhangRegionalTropReceiptPrefix="REGIONAL_TROP_FINAL_PPP_DESIGN:";
// Bookkeeping only. Unlike outageLimitMap, this does NOT introduce an outage
// reset policy or change x/P, Q, model weights, integer authority or the solve.
// It is copied with the authoritative state into private conditional branches.
inline void zhangCaptureRegionalTropDesign(KFState& state,const KFMeas& measurement)
{
    if (!zhangRegionalTropEnabled()) return;
    for(auto it=state.metaDataMap.begin();it!=state.metaDataMap.end();)
        it=it->first.starts_with(zhangRegionalTropReceiptPrefix)?state.metaDataMap.erase(it):std::next(it);
    if (!state.lastPppTransactionCommitted || !state.lastPppQcConverged ||
        measurement.time!=state.time || measurement.H.cols()!=state.x.size() ||
        measurement.H.rows()!=measurement.obsKeys.size() || measurement.R.rows()!=measurement.H.rows() ||
        measurement.R.cols()!=measurement.H.rows()) return;
    for(const auto& [key,column]:state.kfIndexMap) {
        if(key.type!=KF::TROP || key.num!=0 || key.str.empty()) continue;
        bool supported=false;
        for(int row=0;row<measurement.H.rows();++row) {
            const auto& obs=measurement.obsKeys[row];
            if(obs.str!=key.str || (obs.type!=KF::CODE_MEAS && obs.type!=KF::PHAS_MEAS)) continue;
            if(std::isfinite(measurement.H(row,column)) && measurement.H(row,column)!=0 &&
                std::isfinite(measurement.R(row,row)) && measurement.R(row,row)>0) {supported=true;break;}
        }
        if(supported) state.metaDataMap[zhangRegionalTropReceiptPrefix+key.str]=state.time.to_string(9);
    }
}
