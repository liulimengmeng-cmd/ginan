#pragma once
#include "common/zhangFullRank.hpp"
#include "common/zhangIntegerAudit.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>

struct ZhangPhysicalProjectionDiagnostic
{
	std::map<std::string, ZhangExactInteger> requested;
	std::map<std::string, ZhangExactInteger> rebuilt;
	std::map<std::string, ZhangExactInteger> delta;
};

// Current cycle coordinates a=C*N. Column names identify chords, not physical
// scalar ambiguities. Every retained row must expand through the complete C.
struct ZhangProductPhysicalCycleChart
{
    int columns = 0;
    std::map<int, std::string> chordIdentities;
    std::map<int, std::map<std::string, ZhangExactInteger>> expansions;
    // The numeric coordinate/order remains intact. A missing signal version
    // is not a zero physical coefficient and is never inserted in expansions.
    std::map<int, std::string> unexpandedColumns;
    std::set<std::string> activePhysicalIdentities;
    std::set<std::string> cataloguedSignals;
    std::string failureReason = "NONE";

    bool add(int column, const std::string& signal, const ZhangGraphEdge& chord,
        const ZhangGraphBasis& basis, const std::map<ZhangGraphEdge,int>& versions,
        bool retainUnknownColumn=false)
    {
        auto identity = [&](const ZhangGraphEdge& edge) {
            const auto version=versions.find(edge);
            return version==versions.end() || version->second<0 ? std::string{} : signal+"|"+
                edge.receiver+"|"+edge.satellite.id()+"|V"+std::to_string(version->second);
        };
        if(column<0 || column>=columns || basis.treeEdges.contains(chord)) return false;
        if(retainUnknownColumn && !basis.edges.contains(chord))return false;
        const auto cycle=zhangFundamentalCycle(basis,chord);
        if(cycle.empty()) return false;
        std::map<std::string,ZhangExactInteger> row;
        std::map<std::string,ZhangExactInteger> nuisance;
        std::string unknown;
        for(const auto& [edge,value]:cycle)
        {
            if(value==0)continue;
            const auto id=identity(edge);
            if(id.empty()) {
                if(unknown.empty())unknown="UNKNOWN_SIGNAL_ARC_VERSION|SIGNAL="+signal+
                    "|EDGE="+edge.receiver+":"+edge.satellite.id();
            } else row[id]+=value;
            nuisance["R|"+edge.receiver]+=value;
            nuisance["S|"+edge.satellite.id()]+=value;
        }
        for(const auto& [key,value]:nuisance) if(value!=0) return false;
        if(!unknown.empty()) {
            if(!retainUnknownColumn)return false;
            expansions.erase(column);chordIdentities.erase(column);
            unexpandedColumns[column]=std::move(unknown);
            return true;
        }
        const auto id=identity(chord);
        if(id.empty() || row[id]!=1) return false;
        chordIdentities[column]=id;
        expansions[column]=std::move(row);
        unexpandedColumns.erase(column);
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
        std::map<std::string,ZhangExactInteger> candidate;
        for(int c=0;c<columns;++c) {
            if(cycleRow[c]==0) continue;
            const auto it=expansions.find(c);
            if(it==expansions.end()) return false;
            for(const auto& [id,value]:it->second) candidate[id]+=cycleRow[c]*value;
        }
        for(auto it=candidate.begin();it!=candidate.end();)
            if(it->second==0) it=candidate.erase(it); else ++it;
        if(candidate.empty())return false;
        physical=std::move(candidate);return true;
    }
    bool project(const std::map<std::string,ZhangExactInteger>& physical,
		ZhangExactVector& current,
		std::string* reason=nullptr,
		ZhangPhysicalProjectionDiagnostic* diagnostic=nullptr) const
    {
		if(diagnostic)
		{
			*diagnostic={};
			diagnostic->requested=physical;
			diagnostic->delta=physical;
		}
        auto fail=[&](const char* why) {if(reason) *reason=why;return false;};
        current=ZhangExactVector(columns);
        if(physical.empty()) return fail("EMPTY_PHYSICAL_ROW");
        for(const auto& [id,value]:physical)
            if(value!=0 && !activePhysicalIdentities.contains(id)) return fail("ARC_RETIRED");
        for(const auto& [c,id]:chordIdentities) {
            auto it=physical.find(id);if(it!=physical.end()) current[c]=it->second;
        }
        std::map<std::string,ZhangExactInteger> rebuilt;
		const bool rebuiltValid=expand(current,rebuilt);
		if(diagnostic)
		{
			diagnostic->rebuilt=rebuilt;
			diagnostic->delta=physical;
			for(const auto& [identity,coefficient]:rebuilt)
				diagnostic->delta[identity]-=coefficient;
			for(auto it=diagnostic->delta.begin();it!=diagnostic->delta.end();)
				if(it->second==0) it=diagnostic->delta.erase(it); else ++it;
		}
        if(!rebuiltValid || rebuilt!=physical)
            return fail("MISSING_POSTERIOR_COLUMN_OR_NO_CURRENT_REPRESENTATION");
        if(reason) *reason="EXACT_PHYSICAL_ROUNDTRIP";
        return true;
    }
};

/** Filter only current search proposals, never conditioning equations or a
 * mandatory product catalogue. Both target arrays retain the same ordered
 * original row indices; the full numeric coordinate dimension is unchanged.
 * A malformed chart/row is a global failure, not a proposal to silently drop. */
struct ZhangPhysicalTargetProposalSubset
{
    bool valid=false;
    std::string reason="NOT_EVALUATED";
    ZhangExactMatrix numericTargets;
    std::vector<std::map<std::string,ZhangExactInteger>> physicalTargets;
    std::vector<int> retainedIndices,rejectedIndices;
};

// A partial identity catalogue is usable only for explicitly supplied rows
// whose complete physical expansion and inverse pullback are exact. This is
// not a certificate, and never removes a conditioning equation. Structural
// zeros can be allowed for coordinate maps, not as new integer evidence.
inline bool zhangProductPhysicalRowsRepresentable(
    const ZhangProductPhysicalCycleChart& chart, const ZhangExactMatrix& rows,
    bool allowStructuralZero=false)
{
    if(chart.failureReason!="NONE" || chart.columns<0)return false;
    for(const auto& row:rows)if(row.size()!=static_cast<std::size_t>(chart.columns))return false;
    for(const auto& row:rows) {
        const bool zero=std::all_of(row.begin(),row.end(),[](const auto& value){return value==0;});
        if(zero) {if(!allowStructuralZero)return false;continue;}
        std::map<std::string,ZhangExactInteger> physical;
        ZhangExactVector rebuilt;
        if(!chart.expand(row,physical) || !chart.project(physical,rebuilt) || rebuilt!=row)return false;
    }
    return true;
}

inline ZhangPhysicalTargetProposalSubset zhangPhysicalTargetProposalSubset(
    const ZhangProductPhysicalCycleChart& chart,const ZhangExactMatrix& targets,
    bool filterUnexpandedProposals)
{
    ZhangPhysicalTargetProposalSubset out;
    if(chart.failureReason!="NONE") {out.reason=chart.failureReason;return out;}
    if(chart.columns<0) {out.reason="TARGET_COORDINATE_DIMENSION_MISMATCH";return out;}
    for(const auto& row:targets)if(row.size()!=static_cast<std::size_t>(chart.columns)) {
        out.reason="TARGET_COORDINATE_DIMENSION_MISMATCH";return out;
    }
    for(int i=0;i<targets.size();++i) {
        std::map<std::string,ZhangExactInteger> physical;
        if(!chart.expand(targets[i],physical)) {
            bool unknown=false;
            for(int c=0;c<chart.columns;++c)
                unknown |= targets[i][c]!=0 && chart.unexpandedColumns.contains(c);
            out.rejectedIndices.push_back(i);
            if(!filterUnexpandedProposals || !unknown) {
                out.numericTargets.clear();out.physicalTargets.clear();out.retainedIndices.clear();
                out.reason="TARGET_PHYSICAL_EXPANSION_REQUIRED";return out;
            }
            continue;
        }
        out.numericTargets.push_back(targets[i]);
        out.physicalTargets.push_back(std::move(physical));
        out.retainedIndices.push_back(i);
    }
    out.valid=true;out.reason="NONE";return out;
}
