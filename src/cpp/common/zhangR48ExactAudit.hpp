#pragma once
#include "common/zhangIntegerAudit.hpp"
#include <chrono>
#include <map>
#include <string>
struct ZhangR48ExactMetric {
 double seconds=0;std::size_t calls=0,rows=0,columns=0,nonzeros=0,bits=0;
};
inline thread_local std::map<std::string,ZhangR48ExactMetric> zhangR48ExactMetrics;
struct ZhangR48ExactTimer {
 std::string phase;std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
 bool done=false;
 ZhangR48ExactTimer(std::string name,const ZhangExactMatrix& matrix,int columns=0):phase(std::move(name)) {
  auto& m=zhangR48ExactMetrics[phase];++m.calls;m.rows=std::max(m.rows,matrix.size());
  m.columns=std::max(m.columns,matrix.empty()?std::size_t(columns):matrix.front().size());
  std::size_t nnz=0;
  for(const auto& row:matrix)for(const auto& value:row)if(value!=0) {
   ++nnz;ZhangExactInteger a=value<0?-value:value;
   m.bits=std::max(m.bits,std::size_t(boost::multiprecision::msb(a)+1));
  }
  m.nonzeros=std::max(m.nonzeros,nnz);
 }
 void finish(){if(!done){zhangR48ExactMetrics[phase].seconds+=
  std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();done=true;}}
 ~ZhangR48ExactTimer(){finish();}
};
