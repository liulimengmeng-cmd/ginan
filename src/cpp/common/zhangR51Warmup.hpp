#pragma once
#include <cmath>
#include <stdexcept>
#include <string>

// GPST seconds on the same Unix calendar origin as GTime::bigTime.
// A missing setting preserves the original R51 execution path.
inline long double zhangR51ParseArStart(const char* text) {
    if (!text) return 0;
    std::string value(text); std::size_t used=0;
    long double start=std::stold(value,&used);
    if (used!=value.size() || !std::isfinite(start) || start<=0 || std::floor(start)!=start)
        throw std::invalid_argument("Invalid ZHANG_R51_AR_START_GPST_SECONDS");
    return start;
}
inline bool zhangR51FloatWarmup(long double now,long double start) {
    return start>0 && now<start;
}
