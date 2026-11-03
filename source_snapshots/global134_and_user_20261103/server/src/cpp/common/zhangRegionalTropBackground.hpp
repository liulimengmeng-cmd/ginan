#pragma once
#include "common/zhangRegionalTropModel.hpp"
#include "common/acsConfig.hpp"
#include "common/zhangCheckpoint.hpp"
#include "trop/tropModels.hpp"
#include "orbprop/coordinates.hpp"
#include <sstream>
#include <stdexcept>

// Both ends use the existing GINAN GPT2 pressure/water-vapour -> ZHD/ZWD
// formulas. Input is ECEF -> WGS84 ellipsoidal height; GPT2 performs its
// internal geoid/grid height reduction exactly once. Regional H means
// ellipsoidal height, NOT local ENU U and NOT meteorological orthometric H.
inline double zhangRegionalTropBackground(GTime time,const Eigen::Vector3d& xyz)
{
    VectorPos pos=ecef2pos(xyz);
    std::ostringstream trace;
    double dry=0,wet=0,dm=0,wm=0,var=0;
    const double delay=tropGPT2(trace,time,pos,PI/2,dry,dm,wet,wm,var);
    const double ztd=dry+wet;
    if (!(delay>0) || !std::isfinite(ztd) || ztd<0.5 || ztd>4)
        throw std::runtime_error("REGIONAL_TROP_GPT2_BACKGROUND_UNAVAILABLE");
    return ztd;
}
inline std::string zhangRegionalTropGridHash()
{
    if (acsConfig.gpt2grid_files.size()!=1)
        throw std::runtime_error("REGIONAL_TROP_EXACTLY_ONE_GPT2_GRID_REQUIRED");
    std::string reason;
    const auto hash=zhangCheckpointFileSha256(acsConfig.gpt2grid_files.front(),&reason);
    if (hash.size()!=64) throw std::runtime_error("REGIONAL_TROP_GRID_HASH:"+reason);
    return hash;
}
