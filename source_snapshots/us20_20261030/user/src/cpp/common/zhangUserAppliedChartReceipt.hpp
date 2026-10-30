#pragma once

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// A stochastic-coordinate ownership receipt, not an integer decision or an
// estimate. Stored inside KF metadata, so copies/transactions/checkpoints own
// it together with x/P. Runtime catalogues may be rebuilt only from this exact
// applied identity; current product values are never used to infer an owner.
struct ZhangAppliedChartEndpoint
{
    std::string satellite, segment, event;
    int counter=-1, datum=-1;
    long double appliedEpoch=0;
    std::uint64_t physical=0, state=0, integerAuthorization=0;
};

struct ZhangAppliedChartReceipt
{
    std::string receiver, reference;
    int system=0, signal=0;
    long double appliedEpoch=0;
    std::set<std::string> ownedKeys;
    std::vector<ZhangAppliedChartEndpoint> endpoints;
};

inline std::string zhangAppliedChartReceiptKey(const std::string& receiver,int system,int signal)
{
    return "USER_APPLIED_CHART_V1:"+std::to_string(receiver.size())+":"+receiver+":"+
        std::to_string(system)+":"+std::to_string(signal);
}

inline std::string zhangSerializeAppliedChartReceipt(const ZhangAppliedChartReceipt& r)
{
    std::ostringstream s;s<<std::setprecision(20)<<1<<' '<<std::quoted(r.receiver)<<' '<<r.system<<' '
        <<r.signal<<' '<<std::quoted(r.reference)<<' '<<r.appliedEpoch<<' '<<r.ownedKeys.size();
    for (const auto& key:r.ownedKeys)s<<' '<<std::quoted(key);
    s<<' '<<r.endpoints.size();
    for (const auto& e:r.endpoints)s<<' '<<std::quoted(e.satellite)<<' '<<e.counter<<' '<<e.datum<<' '
        <<std::quoted(e.segment)<<' '<<std::quoted(e.event)<<' '<<e.appliedEpoch<<' '
        <<e.physical<<' '<<e.state<<' '<<e.integerAuthorization;
    return s.str();
}

inline bool zhangDeserializeAppliedChartReceipt(const std::string& text,ZhangAppliedChartReceipt& out)
{
    ZhangAppliedChartReceipt r;std::istringstream s(text);int version=0;std::size_t keys=0,ends=0;
    if (!(s>>version>>std::quoted(r.receiver)>>r.system>>r.signal>>std::quoted(r.reference)>>r.appliedEpoch>>keys) ||
        version!=1 || r.receiver.empty() || r.reference.empty() || !std::isfinite(r.appliedEpoch) || keys>10000)return false;
    for(std::size_t i=0;i<keys;++i){std::string key;if(!(s>>std::quoted(key)) || key.empty() || !r.ownedKeys.insert(key).second)return false;}
    if(!(s>>ends) || ends>10000)return false;
    std::set<std::string> satellites;
    for(std::size_t i=0;i<ends;++i)
    {
        ZhangAppliedChartEndpoint e;
        if(!(s>>std::quoted(e.satellite)>>e.counter>>e.datum>>std::quoted(e.segment)>>std::quoted(e.event)>>
            e.appliedEpoch>>e.physical>>e.state>>e.integerAuthorization) || e.satellite.empty() ||
            e.counter<0 || e.datum<0 || !std::isfinite(e.appliedEpoch) || e.event.empty() ||
            !satellites.insert(e.satellite).second)return false;
        r.endpoints.push_back(std::move(e));
    }
    s>>std::ws;if(!s.eof() || !satellites.contains(r.reference))return false;
    out=std::move(r);return true;
}

inline bool zhangAppliedChartOwns(const ZhangAppliedChartReceipt& r,
    const std::string& receiver,int system,int signal,const std::set<std::string>& actual,long double now)
{
    return r.receiver==receiver && r.system==system && r.signal==signal &&
        r.ownedKeys==actual && !actual.empty() && std::isfinite(now) && r.appliedEpoch<=now+1e-6L;
}
