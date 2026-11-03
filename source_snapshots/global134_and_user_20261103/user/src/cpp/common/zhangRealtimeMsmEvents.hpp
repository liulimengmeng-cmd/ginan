#pragma once
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

// Loss-of-lock possibilities are evidence, not confirmed physical slips.
// Normal reports cannot erase an unconsumed event. Future observations never
// annotate an earlier processing epoch. Each new generation has empty state.
namespace zhangRealtimeMsmEvents
{
using Key=std::tuple<std::string,std::string,int>;
struct Signal
{
    double time=0,lower=0,upper=0,consumed=-std::numeric_limits<double>::infinity();
    bool seen=false,half=false;
    std::deque<double> events;
};
struct Tracker
{
    std::mutex mutex;
    std::map<Key,Signal> signals;
    bool observe(const std::string& receiver,double time,const std::string& satellite,int code,
                 double lower,double upper,bool half)
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (!std::isfinite(time)||lower<0||!(upper>lower))
            throw std::runtime_error("REALTIME_MSM_INVALID_LOCK_INTERVAL");
        auto& signal=signals[{receiver,satellite,code}];
        bool possible=false;
        if (signal.seen)
        {
            if(time<=signal.time) return false; // Duplicate/past input cannot rewrite history.
            const double dt=time-signal.time;
            // Compare true report spacing with the quantised lock interval.
            // A reset may fit inside dt even if the intervals also permit an
            // uninterrupted history. This uncertainty must remain visible.
            possible=dt>5 || lower+1e-8<dt ||
                lower>signal.upper+dt+1e-8 || upper+1e-8<signal.lower+dt || half!=signal.half;
        }
        possible|=half;
        if (possible && time>signal.consumed)
        {
            if(signal.events.size()>=128) throw std::runtime_error("REALTIME_MSM_EVENT_QUEUE_LIMIT");
            signal.events.push_back(time);
        }
        signal.time=time;signal.lower=lower;signal.upper=upper;signal.half=half;signal.seen=true;
        if(signals.size()>8192)
        {
            for(auto it=signals.begin();it!=signals.end();)
                if(it->second.time<time-120) it=signals.erase(it); else ++it;
            if(signals.size()>8192) throw std::runtime_error("REALTIME_MSM_SIGNAL_LIMIT");
        }
        return possible;
    }
    unsigned consume(const std::string& receiver,double time,const std::string& satellite,int code)
    {
        std::lock_guard<std::mutex> lock(mutex);
        const auto it=signals.find({receiver,satellite,code});
        if(it==signals.end()) return 0;
        auto& signal=it->second;
        if(time<signal.consumed) throw std::runtime_error("REALTIME_MSM_PROCESSING_TIME_REVERSED");
        const double after=std::isfinite(signal.consumed)?signal.consumed:time-30;
        unsigned count=0;
        while(!signal.events.empty() && signal.events.front()<=time+1e-8)
        {
            if(signal.events.front()>after+1e-8) ++count;
            signal.events.pop_front();
        }
        signal.consumed=time;
        return count;
    }
};
inline Tracker tracker;
}
