#pragma once

#include "common/algebra.hpp"
#include "common/zhangUserArcGeneration.hpp"
#include "common/zhangUserGonsisComponents.hpp"
#include "common/zhangUserDatumPolicy.hpp"
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// A legal integer coordinate is n = a - offset, not necessarily the current
// effective ambiguity a. The certificate belongs to that coordinate's chart.
// No fixed integer value, search result, or HOLD constraint is stored here.
inline std::string zhangUserChartKey(const std::string& receiver,
    const SatSys& satellite, int code)
{
    return "USER_INTEGER_CHART:" + std::to_string(receiver.size()) + ":" + receiver +
        ":" + satellite.id() + ":" + std::to_string(code);
}

inline double zhangUserChartOffset(const KFState& state, const std::string& receiver,
    const SatSys& satellite, int code)
{
    const auto key=zhangUserChartKey(receiver,satellite,code);
    const auto generation=state.metaDataMap.find(key+":GENERATION");
    if (generation!=state.metaDataMap.end())
    {
        const auto current=std::to_string(zhangUserArcGeneration(state.metaDataMap,
            zhangUserArcGenerationKey(true,receiver,satellite.id(),code)))+":"+
            std::to_string(zhangUserArcGeneration(state.metaDataMap,
            zhangUserArcGenerationKey(false,receiver,satellite.id(),code)));
        if (generation->second!=current) return 0; // re-created current-chart state, not the old a
    }
    const auto found = state.metaDataMap.find(zhangUserChartKey(receiver,satellite,code)+":OFFSET");
    if (found == state.metaDataMap.end()) return 0;
    std::size_t used = 0;
    const double value = std::stod(found->second, &used);
    if (used != found->second.size() || !std::isfinite(value))
        throw std::runtime_error("INVALID_USER_INTEGER_CHART_OFFSET");
    return value;
}

inline void zhangUserSetChartOffset(KFState& state, const std::string& receiver,
    const SatSys& satellite, int code, double value)
{
    if (!std::isfinite(value)) throw std::runtime_error("INVALID_USER_INTEGER_CHART_OFFSET");
    std::ostringstream out; out << std::setprecision(17) << value;
    state.metaDataMap[zhangUserChartKey(receiver,satellite,code)+":OFFSET"] = out.str();
}

inline std::string zhangUserChartGeneration(const KFState& state,
    const std::string& receiver, const SatSys& satellite, int code)
{
    return std::to_string(zhangUserArcGeneration(state.metaDataMap,
        zhangUserArcGenerationKey(true,receiver,satellite.id(),code))) + ":" +
        std::to_string(zhangUserArcGeneration(state.metaDataMap,
        zhangUserArcGenerationKey(false,receiver,satellite.id(),code)));
}

inline ZhangUserGonsisComponentCertificate zhangUserChartAnchor(const KFState& state,
    const std::string& receiver, const SatSys& satellite, int code)
{
    const auto key = zhangUserChartKey(receiver,satellite,code);
    auto get = [&](const std::string& suffix) -> std::string
    { const auto it=state.metaDataMap.find(key+suffix);
      return it==state.metaDataMap.end()?std::string{}:it->second; };
    if (get(":GENERATION") != zhangUserChartGeneration(state,receiver,satellite,code) ||
        get(":COMPONENT").empty() || get(":VERSION").empty() || get(":ALIGNMENT").empty()) return {};
    return {true,true,get(":COMPONENT"),std::stoi(get(":VERSION")),std::stoi(get(":ALIGNMENT"))};
}

inline void zhangUserRememberChartAnchor(KFState& state, const std::string& receiver,
    const SatSys& satellite, int code, const ZhangUserGonsisComponentCertificate& certificate)
{
    if (zhangUserGonsisComponentCertificateValid(zhangUserChartAnchor(state,receiver,satellite,code)) ||
        !zhangUserGonsisComponentCertificateValid(certificate)) return;
    const auto key=zhangUserChartKey(receiver,satellite,code);
    if (state.metaDataMap.count(key+":UNANCHORED_MAP")) return;
    state.metaDataMap[key+":GENERATION"]=zhangUserChartGeneration(state,receiver,satellite,code);
    state.metaDataMap[key+":COMPONENT"]=certificate.component;
    state.metaDataMap[key+":VERSION"]=std::to_string(certificate.version);
    state.metaDataMap[key+":ALIGNMENT"]=std::to_string(certificate.alignmentGeneration);
    if(zhangUserDatumPolicyRequested())state.metaDataMap[key+":AUTHORIZATION"]=std::to_string(zhangUserArcGeneration(
        state.metaDataMap,zhangUserIntegerAuthorizationGenerationKey(receiver,satellite.id(),code)));
    zhangUserSetChartOffset(state,receiver,satellite,code,0);
}

inline void zhangUserRetireChartAnchor(KFState& state, const std::string& receiver,
    const SatSys& satellite, int code)
{
    const auto key=zhangUserChartKey(receiver,satellite,code);
    for (const auto suffix : {":GENERATION",":COMPONENT",":VERSION",":ALIGNMENT",
                             ":OFFSET",":UNANCHORED_MAP",":AUTHORIZATION"}) state.metaDataMap.erase(key+suffix);
}

inline std::string zhangUserChartIsolationKey(const std::string& receiver,E_Sys system,int code)
{
    return "USER_PRODUCT_CHART_ISOLATION:"+std::to_string(receiver.size())+":"+receiver+":"+
        std::to_string(static_cast<int>(system))+":"+std::to_string(code);
}

inline bool zhangUserChartIsolated(const KFState& state,const std::string& receiver,
    E_Sys system,int code)
{
    return state.metaDataMap.count(zhangUserChartIsolationKey(receiver,system,code));
}

enum class ZhangUserUnknownChartPolicy { RECOVER, ISOLATE, RELEASE };

inline ZhangUserUnknownChartPolicy zhangUserUnknownChartPolicy()
{
    const char* value=std::getenv("ZHANG_USER_UNKNOWN_CHART");
    // Recovery is a coordinate-lifecycle operation, not permission to invent
    // an old->new chart map. ISOLATE remains an explicit regression control.
    if (!value || !*value || std::string(value)=="RECOVER") return ZhangUserUnknownChartPolicy::RECOVER;
    if (std::string(value)=="ISOLATE") return ZhangUserUnknownChartPolicy::ISOLATE;
    if (std::string(value)=="RELEASE") return ZhangUserUnknownChartPolicy::RELEASE;
    throw std::runtime_error("UNKNOWN_USER_PRODUCT_CHART_POLICY");
}

inline bool zhangUserUnknownChartIsolationEnabled()
{
    return zhangUserUnknownChartPolicy()==ZhangUserUnknownChartPolicy::ISOLATE;
}

inline bool zhangUserCanonicalAmbiguityRow(const KFState& state,
    const std::string& receiver,const SatSys& satellite,int code,VectorXd& row)
{
    const double offset=zhangUserChartOffset(state,receiver,satellite,code);
    if (offset==0) return true;
    const auto one=state.kfIndexMap.find(state.oneKey);
    if (one==state.kfIndexMap.end() || state.x(one->second)!=1 ||
        state.dx(one->second)!=0 || state.P.row(one->second).cwiseAbs().maxCoeff()!=0 ||
        state.P.col(one->second).cwiseAbs().maxCoeff()!=0) return false;
    row(one->second)-=offset;
    return true;
}

// History is attached to the anchored integer variable n, not its current
// representation a=n+offset. Revoked integer authorization still changes its
// history identity; physical/state generations are independently in the caller.
inline bool zhangUserChartHistoryAuthorized(const KFState& state,const std::string& receiver,const SatSys& satellite,int code)
{
    const auto key=zhangUserChartKey(receiver,satellite,code);
    if(!zhangUserGonsisComponentCertificateValid(zhangUserChartAnchor(state,receiver,satellite,code)) ||
        zhangUserChartIsolated(state,receiver,satellite.sys,code) || state.metaDataMap.count(key+":UNANCHORED_MAP") ||
        !state.metaDataMap.count(key+":OFFSET"))return false;
    const auto authorization=state.metaDataMap.find(key+":AUTHORIZATION");
    return authorization!=state.metaDataMap.end() && authorization->second==std::to_string(zhangUserArcGeneration(
        state.metaDataMap,zhangUserIntegerAuthorizationGenerationKey(receiver,satellite.id(),code)));
}

inline std::string zhangUserProductHistoryIdentity(const KFState& state,const std::string& receiver,
    const SatSys& satellite,int code,const ZhangUserGonsisComponentCertificate& certificate)
{
    if(!zhangUserGonsisComponentCertificateValid(certificate) ||
        (zhangUserDatumPolicyRequested() && !zhangUserChartHistoryAuthorized(state,receiver,satellite,code)))return {};
    std::ostringstream out;out<<std::setprecision(17)<<"INTERNAL_IDENTIFIED_INTEGER_CHART:"
        <<satellite.id()<<":"<<code<<":"<<zhangUserGonsisComponentToken(certificate);
    if(zhangUserDatumPolicyRequested())
        out<<":CANONICAL_HISTORY_V1:IA:"<<zhangUserArcGeneration(state.metaDataMap,
            zhangUserIntegerAuthorizationGenerationKey(receiver,satellite.id(),code));
    else out<<":OFFSET:"<<zhangUserChartOffset(state,receiver,satellite,code);
    return out.str();
}
