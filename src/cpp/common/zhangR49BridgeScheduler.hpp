#pragma once
#include <set>
#include <string>
#include <utility>
struct ZhangR49BridgeScheduler {
    int attempts=0;
    std::set<std::string> testedRequests;
};
template<class Domain, class Validator>
bool zhangR49CommitTrial(Domain& safe, Domain trial, Validator validate) {
    if(!validate(trial)) return false;
    safe=std::move(trial);
    return true;
}
