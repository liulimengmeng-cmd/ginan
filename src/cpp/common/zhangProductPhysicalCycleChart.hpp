#pragma once
#include "common/zhangFullRank.hpp"
#include "common/zhangIntegerAudit.hpp"
#include <map>
#include <set>
#include <string>

// Current cycle coordinates a=C*N. Column names identify chords, not physical
// scalar ambiguities. Every retained row must expand through the complete C.
struct ZhangProductPhysicalCycleChart
{
    int columns = 0;
    std::map<int, std::string> chordIdentities;
    std::map<int, std::map<std::string, ZhangExactInteger>> expansions;
    std::set<std::string> activePhysicalIdentities;
    std::set<std::string> cataloguedSignals;
    std::string failureReason = "NONE";

    bool add(int column, const std::string& signal, const ZhangGraphEdge& chord,
        const ZhangGraphBasis& basis, const std::map<ZhangGraphEdge,int>& versions)
    {
        auto identity = [&](const ZhangGraphEdge& edge) {
            const auto version=versions.find(edge);
            return version==versions.end() ? std::string{} : signal+"|"+
                edge.receiver+"|"+edge.satellite.id()+"|V"+std::to_string(version->second);
        };
        if(column<0 || column>=columns || basis.treeEdges.contains(chord)) return false;
        const auto cycle=zhangFundamentalCycle(basis,chord);
        if(cycle.empty()) return false;
        std::map<std::string,ZhangExactInteger> row;
        std::map<std::string,ZhangExactInteger> nuisance;
        for(const auto& [edge,value]:cycle)
        {
            const auto id=identity(edge);
            if(id.empty()) return false;
            row[id]+=value;
            nuisance["R|"+edge.receiver]+=value;
            nuisance["S|"+edge.satellite.id()]+=value;
        }
        for(const auto& [key,value]:nuisance) if(value!=0) return false;
        const auto id=identity(chord);
        if(id.empty() || row[id]!=1) return false;
        chordIdentities[column]=id;
        expansions[column]=std::move(row);
        if(cataloguedSignals.insert(signal).second) for(const auto& edge:basis.edges) {
            const auto name=identity(edge);
            if(!name.empty()) activePhysicalIdentities.insert(name);
        }
        return true;
    }
    bool expand(const ZhangExactVector& cycleRow,
        std::map<std::string,ZhangExactInteger>& physical) const
    {
        physical.clear();
        if(cycleRow.size()!=static_cast<std::size_t>(columns)) return false;
        for(int c=0;c<columns;++c) {
            if(cycleRow[c]==0) continue;
            const auto it=expansions.find(c);
            if(it==expansions.end()) return false;
            for(const auto& [id,value]:it->second) physical[id]+=cycleRow[c]*value;
        }
        for(auto it=physical.begin();it!=physical.end();)
            if(it->second==0) it=physical.erase(it); else ++it;
        return !physical.empty();
    }
    bool project(const std::map<std::string,ZhangExactInteger>& physical,
        ZhangExactVector& current, std::string* reason=nullptr) const
    {
        auto fail=[&](const char* why) {if(reason) *reason=why;return false;};
        current=ZhangExactVector(columns);
        if(physical.empty()) return fail("EMPTY_PHYSICAL_ROW");
        for(const auto& [id,value]:physical)
            if(value!=0 && !activePhysicalIdentities.contains(id)) return fail("ARC_RETIRED");
        for(const auto& [c,id]:chordIdentities) {
            auto it=physical.find(id);if(it!=physical.end()) current[c]=it->second;
        }
        std::map<std::string,ZhangExactInteger> rebuilt;
        if(!expand(current,rebuilt) || rebuilt!=physical)
            return fail("MISSING_POSTERIOR_COLUMN_OR_NO_CURRENT_REPRESENTATION");
        if(reason) *reason="EXACT_PHYSICAL_ROUNDTRIP";
        return true;
    }
};
