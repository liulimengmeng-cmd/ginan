#include "common/zhangR51FloatReuse.hpp"
#include <cassert>
#include <iostream>
int main(){
std::map<std::string,std::string> env{{"ZHANG_R51_REUSE_RATIO_20240717","1"},{"ZHANG_R51_ENABLE","1"},{"ZHANG_R51_RATIO_ONLY","1"},{"ZHANG_R51_AR_START_GPST_SECONDS","1405216800"},{"ZHANG_R49_FUSION_WEIGHT","0.10"},{"OPENBLAS_NUM_THREADS","4"},{"OMP_NUM_THREADS","4"},{"MKL_NUM_THREADS","1"}};
std::string h[]={"9e5ea901db26708ecb2f772b7ad5efde886d919491811a536302a915ec776c3d","4463f3a4c97f0e213c15b995d9e215a873178428a0cf7bcab7e436c98f96b80a","6ee5be509848128c71ff23e267b9e10df839c746ec06cd363fe3462ccbbbe841","4f6351a7db3a2af79ae728ffc1d8d6e9b19986cafe49d6c47c6a703cfeac4874"};
auto pass=[&](){return zhangR51AuditedRatioReuseAllowed(h[0],h[1],h[2],h[3],env);};assert(pass());
for(int i=0;i<4;++i){auto saved=h[i];h[i]="WRONG";assert(!pass());h[i]=saved;}
auto original=env;for(const auto& [key,value]:original){env[key]="WRONG";assert(!pass());env=original;env.erase(key);assert(!pass());env=original;}
assert(pass());std::cout<<"PASS ratio-only checkpoint and 20 mismatch/missing rejection cases\n";}
