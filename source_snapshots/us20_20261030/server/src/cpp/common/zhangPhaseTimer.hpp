#pragma once
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <ostream>
#include <string>

// One per-thread hierarchy. Header-only kernels inherit the caller's trace;
// standalone mathematical tests need no logger and pay no output cost.
class ZhangPhaseTimer
{
    std::ostream* trace;
    std::string phase,epoch;
    std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    double children=0;
    ZhangPhaseTimer* parent;
    inline static thread_local ZhangPhaseTimer* current=nullptr;
public:
    static std::ostream* currentTrace() noexcept {return current?current->trace:nullptr;}
    ZhangPhaseTimer(std::ostream& t,std::string p,std::string e={})
        :trace(&t),phase(std::move(p)),epoch(std::move(e)),parent(current){current=this;}
    explicit ZhangPhaseTimer(std::string p)
        :trace(current?current->trace:nullptr),phase(std::move(p)),parent(current){current=this;}
    ~ZhangPhaseTimer()
    {
        const double ms=std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        current=parent;if(parent)parent->children+=ms;
        if(!trace)return;
        *trace<<"\nZHANG_AR_PHASE_TIMER phase="<<phase;
        if(!epoch.empty())*trace<<" epoch="<<std::quoted(epoch);
        *trace<<" inclusive_ms="<<ms<<" exclusive_ms="<<std::max(0.0,ms-children)
            <<" root="<<(parent==nullptr)<<" parent="<<(parent?parent->phase:"NONE");
    }
    ZhangPhaseTimer(const ZhangPhaseTimer&)=delete;
    ZhangPhaseTimer& operator=(const ZhangPhaseTimer&)=delete;
};
