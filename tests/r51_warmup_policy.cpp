#include "common/zhangR51Warmup.hpp"
#include <cassert>
#include <iostream>
int main() {
 const auto start=zhangR51ParseArStart("1405296000");
 assert(start==1405296000.L);
 assert(zhangR51FloatWarmup(start-30,start));
 assert(zhangR51FloatWarmup(start-0.01L,start));
 assert(!zhangR51FloatWarmup(start,start));
 assert(!zhangR51FloatWarmup(start+30,start));
 assert(!zhangR51FloatWarmup(start,zhangR51ParseArStart(nullptr)));
 for(auto invalid:{"","0","-1","NaN","inf","1405296000x","1405296000.5"}) {
  bool threw=false;try{zhangR51ParseArStart(invalid);}catch(...){threw=true;}assert(threw);
 }
 std::cout << "PASS disabled, before/at/after GPST boundary, invalid settings rejected\n";
}
