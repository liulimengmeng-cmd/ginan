#pragma once
#include "common/zhangIntegerAudit.hpp"
#include <chrono>
#include <map>
#include <string>
#include "common/zhangPhaseTimer.hpp"
struct ZhangR48ExactMetric {
 double seconds=0;std::size_t calls=0,rows=0,columns=0,nonzeros=0,bits=0;
};
inline thread_local std::map<std::string,ZhangR48ExactMetric> zhangR48ExactMetrics;
struct ZhangR48ExactSample {
 std::string phase;double seconds=0;std::size_t rows=0,columns=0,nonzeros=0,bits=0;
};
// Cleared per official-product call. Actual per-input samples, NOT a density
// assembled from unrelated maxima. Cap diagnostics independently of history.
inline thread_local std::vector<ZhangR48ExactSample> zhangR48ExactSamples;
struct ZhangR48ExactTimer {
 std::string phase;std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
 bool done=false;
 ZhangR48ExactSample sample;
 ZhangR48ExactTimer(std::string name,const ZhangExactMatrix& matrix,int columns=0):phase(std::move(name)) {
  auto& m=zhangR48ExactMetrics[phase];++m.calls;m.rows=std::max(m.rows,matrix.size());
  m.columns=std::max(m.columns,matrix.empty()?std::size_t(columns):matrix.front().size());
  std::size_t nnz=0,bits=0;
  for(const auto& row:matrix)for(const auto& value:row)if(value!=0) {
   ++nnz;ZhangExactInteger a=value<0?-value:value;
   m.bits=std::max(m.bits,std::size_t(boost::multiprecision::msb(a)+1));
   bits=std::max(bits,std::size_t(boost::multiprecision::msb(a)+1));
  }
  m.nonzeros=std::max(m.nonzeros,nnz);
  sample={phase,0,matrix.size(),matrix.empty()?std::size_t(columns):matrix.front().size(),nnz,bits};
 }
 void finish(){if(!done){sample.seconds=
  std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
  zhangR48ExactMetrics[phase].seconds+=sample.seconds;
  if(zhangR48ExactSamples.size()<512)zhangR48ExactSamples.push_back(sample);done=true;}}
 ~ZhangR48ExactTimer(){finish();}
};
