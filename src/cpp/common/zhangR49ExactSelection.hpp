#pragma once
#include <map>
#include <cstdint>
#include "common/zhangIntegerDecisionProof.hpp"
// Fraction-free sparse elimination only determines rational independence.
// It never substitutes for integer feasibility or saturated kernel transport.
struct ZhangR49IncrementalRank {
    using Row=std::map<int,ZhangExactInteger>;
    std::map<int,Row> basis;
    static ZhangExactInteger gcd(ZhangExactInteger a,ZhangExactInteger b) {
        if(a<0)a=-a;if(b<0)b=-b;
        while(b!=0){ZhangExactInteger c=a%b;a=b;b=c;}return a;
    }
    static void normalize(Row& row) {
        ZhangExactInteger g=0;
        for(auto i=row.begin();i!=row.end();) {
            if(i->second==0)i=row.erase(i);else {g=gcd(g,i->second);++i;}
        }
        if(row.empty())return;if(row.begin()->second<0)g=-g;
        if(g!=1)for(auto& [c,v]:row)v/=g;
    }
    bool append(const ZhangExactVector& input) {
        Row row;for(int c=0;c<input.size();++c)if(input[c]!=0)row[c]=input[c];
        normalize(row);
        while(!row.empty()) {
            const int pivot=row.begin()->first;auto old=basis.find(pivot);
            if(old==basis.end()){basis.emplace(pivot,std::move(row));return true;}
            const ZhangExactInteger a=old->second.begin()->second,b=row.begin()->second,g=gcd(a,b);
            for(auto& [c,v]:row)v*=a/g;
            for(const auto& [c,v]:old->second)row[c]-=(b/g)*v;
            normalize(row);
        }
        return false;
    }
};
struct ZhangR49ProofCache {
    std::map<std::vector<std::uintptr_t>,ZhangDecisionRiskClosure> cache;
    const ZhangDecisionRiskClosure& closure(const ZhangDecisionProofs& roots) {
        std::vector<std::uintptr_t> key;
        for(const auto& p:roots)key.push_back(reinterpret_cast<std::uintptr_t>(p.get()));
        auto i=cache.find(key);if(i!=cache.end())return i->second;
        return cache.emplace(std::move(key),zhangDecisionRiskClosure(roots)).first->second;
    }
};
