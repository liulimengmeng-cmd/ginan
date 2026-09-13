#pragma once
#include <cstdlib>
#include <cmath>
#include <stdexcept>
struct ZhangR49RouteFusionPolicy {
 double weight=0;
 int maximumSearches=4;
 ZhangR49RouteFusionPolicy() {
  if(const char* value=std::getenv("ZHANG_R49_FUSION_WEIGHT")) {
   char* end=nullptr;weight=std::strtod(value,&end);
   if(end==value || *end!='\0' || !std::isfinite(weight) || weight<0 || weight>0.5)
    throw std::runtime_error("INVALID_R49_FUSION_WEIGHT");
  }
 }
};
