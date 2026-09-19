#include "common/zhangR51FloatReuse.hpp"
#include <cassert>
#include <iostream>
int main(){
std::map<std::string,std::string> env{{"ZHANG_R51_REUSE_FLOAT_20240717","1"},{"ZHANG_R51_ENABLE","1"},{"ZHANG_R51_RATIO_ONLY","1"},{"ZHANG_R51_AR_START_GPST_SECONDS","1405216800"},{"ZHANG_R49_FUSION_WEIGHT","0.10"},{"OPENBLAS_NUM_THREADS","4"},{"OMP_NUM_THREADS","4"},{"MKL_NUM_THREADS","1"}};
std::string h[]={"ab397e655774da5c0a9632becd08fe171508b93ffeaf2dd7b04de899bae2182a","b7b0d35686204154b44448eb097fb6990ca07453cecf9c1fc20f2e80f192cc38","02bbaabbe2a836ec3dfc08d2ad206326c9dbddf41237692e0d3806b0043b764f","4f6351a7db3a2af79ae728ffc1d8d6e9b19986cafe49d6c47c6a703cfeac4874"};
auto pass=[&](){return zhangR51AuditedFloatReuseAllowed(h[0],h[1],h[2],h[3],env);};assert(pass());
for(int i=0;i<4;++i){auto saved=h[i];h[i]="WRONG";assert(!pass());h[i]=saved;}
auto original=env;for(const auto& [key,value]:original){env[key]="WRONG";assert(!pass());env=original;env.erase(key);assert(!pass());env=original;}
assert(pass());std::cout<<"PASS pinned checkpoint and 20 mismatch/missing rejection cases\n";}
