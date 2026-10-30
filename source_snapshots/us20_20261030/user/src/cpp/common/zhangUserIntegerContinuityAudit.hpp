#pragma once

#include "common/zhangUserResidualQuotient.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

struct ZhangUserIntegerAuditSnapshot
{
    std::vector<std::string> columns;
    std::vector<std::string> owners;
    ZhangExactMatrix rows;
    ZhangExactVector values;
    double time = 0;
};

struct ZhangUserIntegerContinuityDecision
{
    bool valid = false;
    bool contradiction = false;
    bool confirmationEligible = false;
    int confirmationCount = 0;
    int isolatedColumns = 0;
    int isolatedOwners = 0;
    int retainedRank = 0;
    int auditEvidenceRank = 0;
    int auditEvidenceColumns = 0;
    ZhangUserIntegerAuditSnapshot output;
    std::string reason = "NOT_EVALUATED";
};

// Bounded audit metadata only: one HNF-compressed domain per stable-owner
// catalogue and one quarantined identity per owner. It NEVER conditions FLOAT,
// never declares an old candidate correct, and never produces HOLD factors.
// Targets are exact physical/product SD equations, so chart/reference changes
// are reindexed algebraically instead of being called integer contradictions.
class ZhangUserIntegerContinuityAudit
{
public:
    ZhangUserIntegerContinuityDecision observe(const std::string& group,
        double time, const std::vector<std::string>& columns,
        const ZhangExactMatrix& rows, const ZhangExactVector& values,
        int confirmations, bool isolateContradictions,
        const std::vector<std::string>& columnOwners = {})
    {
        ZhangUserIntegerContinuityDecision out;
        if (group.empty() || !std::isfinite(time) || confirmations < 1)
        { out.reason = "INVALID_CONTINUITY_AUDIT_INPUT"; return out; }
        ZhangUserIntegerAuditSnapshot current;
        current.columns = columns; current.time = time;
        std::set<std::string> unique(columns.begin(), columns.end());
        if (unique.size() != columns.size() || unique.count(""))
        { out.reason = "INVALID_PHYSICAL_COLUMN_IDENTITIES"; return out; }
        const auto& owners = columnOwners.empty() && !isolateContradictions ? columns : columnOwners;
        std::set<std::string> uniqueOwners(owners.begin(), owners.end());
        if (owners.size() != columns.size() || uniqueOwners.size() != owners.size() || uniqueOwners.count(""))
        { out.reason = "EXPLICIT_STABLE_COLUMN_OWNERS_REQUIRED"; return out; }
        current.owners = owners;
        const auto canonical = zhangAuditUserIntegerUnion(rows, values, columns.size());
        if (!canonical.valid) { out.reason = canonical.reason; return out; }
        current.rows = canonical.rows; current.values = canonical.values;
        auto& state = states_[group];
        if (state.observed && time < state.lastObservedTime - 1e-6)
        { state = {}; out.reason = "TIME_REVERSED_AUDIT_RESET"; return out; }
        // Absence is not retirement evidence: a temporarily QC-excluded arc
        // can reappear with exactly the same physical/product identity. Only
        // an observed replacement identity for the same stable owner releases
        // its quarantine. One identity per stable owner bounds metadata size.
        for (int c = 0; c < columns.size(); ++c)
        {
            state.latestIdentities[owners[c]] = columns[c];
            const auto q = state.quarantined.find(owners[c]);
            if (q != state.quarantined.end() && q->second != columns[c])
                state.quarantined.erase(q);
        }

        // Keep absent same-identity variables in AUDIT evidence. A missing
        // target or an empty/rejected search is not physical retirement proof.
        // Owner identity replacement is the only catalogue retirement below.
        auto auditColumns = columns;
        auto auditOwners = owners;
        for (int c = 0; c < state.previous.columns.size(); ++c)
            if (state.latestIdentities.at(state.previous.owners[c]) == state.previous.columns[c] &&
                !unique.count(state.previous.columns[c]))
            {
                auditColumns.push_back(state.previous.columns[c]);
                auditOwners.push_back(state.previous.owners[c]);
            }
        auto previous = reindex(state.previous, auditColumns);
        const auto proposed = reindex(current, auditColumns);
        if (!previous.rows.empty() && !proposed.rows.empty())
        {
            auto combinedRows = previous.rows;
            auto combinedValues = previous.values;
            combinedRows.insert(combinedRows.end(), proposed.rows.begin(), proposed.rows.end());
            combinedValues.insert(combinedValues.end(), proposed.values.begin(), proposed.values.end());
            if (!zhangAuditUserIntegerUnion(combinedRows, combinedValues, auditColumns.size()).valid)
            {
                out.contradiction = true;
                if (isolateContradictions)
                {
                    bool localized = false;
                    // First localize rows individually inconsistent with the
                    // old feasible domain. No assumption that the old is true.
                    for (std::size_t r = 0; r < rows.size(); ++r)
                    {
                        auto testRows = previous.rows; auto testValues = previous.values;
                        auto expandedRow = rows[r]; expandedRow.resize(auditColumns.size());
                        testRows.push_back(std::move(expandedRow)); testValues.push_back(values[r]);
                        if (zhangAuditUserIntegerUnion(testRows, testValues, auditColumns.size()).valid) continue;
                        localized = true;
                        for (std::size_t c = 0; c < columns.size(); ++c)
                            if (rows[r][c] != 0) state.quarantined[owners[c]] = columns[c];
                    }
                    if (!localized)
                        for (const auto& row : current.rows)
                            for (std::size_t c = 0; c < columns.size(); ++c)
                                if (row[c] != 0) state.quarantined[owners[c]] = columns[c];
                }
            }
        }
        if (!state.quarantined.empty())
        {
            std::vector<int> safe;
            for (int c = 0; c < columns.size(); ++c)
                if (!state.quarantined.count(owners[c])) safe.push_back(c);
            const auto restricted = zhangRestrictUserIntegerLedger(
                current.rows, current.values, safe, columns.size());
            if (!restricted.valid) { out.reason = restricted.reason; return out; }
            current.rows.clear(); current.values = restricted.values;
            for (const auto& row : restricted.rows)
            {
                ZhangExactVector full(columns.size());
                for (int c = 0; c < safe.size(); ++c) full[safe[c]] = row[c];
                current.rows.push_back(std::move(full));
            }
        }
        for (int c = 0; c < columns.size(); ++c)
            if (state.quarantined.count(owners[c])) ++out.isolatedColumns;
        out.isolatedOwners = state.quarantined.size();
        // Compare equations in a sorted physical catalogue; HNF removes row
        // insertion order and unimodular reference-chart reparameterization.
        const std::string signature = canonicalSignature(current);
        const double delta = state.observed ? time - state.lastObservedTime : 0;
        const bool sameCandidate = !signature.empty() && signature == state.signature &&
            delta >= -1e-6 && delta <= 45 + 1e-6;
        if (!sameCandidate) state.confirmationCount = signature.empty() ? 0 : 1;
        else if (delta > 1e-6) ++state.confirmationCount;
        // Duplicate invocation in one epoch does not manufacture evidence.
        out.confirmationCount = state.confirmationCount;
        out.confirmationEligible = !signature.empty() && state.confirmationCount >= confirmations;
        state.signature = signature;
        state.lastObservedTime = time;
        state.observed = true;
        // Compress metadata evidence only. Current output/confirmation above
        // never receives these previous equations. In a conflict neither side
        // is treated as truth: safe mode removes implicated owners; legacy
        // audit-only retains the new candidate plus unaffected old consequences.
        std::vector<int> historySafe;
        for (int c = 0; c < auditColumns.size(); ++c)
        {
            bool remove = state.quarantined.count(auditOwners[c]);
            if (out.contradiction && !isolateContradictions && c < columns.size())
                remove = std::any_of(rows.begin(), rows.end(), [&](const auto& row) { return row[c] != 0; });
            if (!remove) historySafe.push_back(c);
        }
        const auto historyRestriction = zhangRestrictUserIntegerLedger(previous.rows,
            previous.values, historySafe, auditColumns.size());
        if (!historyRestriction.valid)
        { out.reason = historyRestriction.reason; return out; }
        ZhangExactMatrix historyRows;
        for (const auto& row : historyRestriction.rows)
        {
            ZhangExactVector full(auditColumns.size());
            for (int c = 0; c < historySafe.size(); ++c) full[historySafe[c]] = row[c];
            historyRows.push_back(std::move(full));
        }
        auto historyValues = historyRestriction.values;
        const auto cleanCurrent = reindex(current, auditColumns);
        historyRows.insert(historyRows.end(), cleanCurrent.rows.begin(), cleanCurrent.rows.end());
        historyValues.insert(historyValues.end(), cleanCurrent.values.begin(), cleanCurrent.values.end());
        const auto history = zhangAuditUserIntegerUnion(historyRows, historyValues, auditColumns.size());
        if (!history.valid) { out.reason = "AUDIT_METADATA_DOMAIN_CONFLICT:" + history.reason; return out; }
        state.previous.columns = auditColumns; state.previous.owners = auditOwners;
        state.previous.rows = history.rows; state.previous.values = history.values;
        state.previous.time = time;
        out.auditEvidenceRank = history.rank;
        out.auditEvidenceColumns = auditColumns.size();
        out.output = std::move(current);
        out.retainedRank = out.output.rows.size();
        out.valid = true;
        out.reason = out.contradiction ? "SAME_PHYSICAL_DOMAIN_INTEGER_CONTRADICTION" :
            (out.isolatedColumns ? "LIVE_CONTRADICTORY_COLUMNS_ISOLATED" :
             (out.confirmationEligible ? "CORRELATED_CANDIDATE_CONFIRMATION_PASSED" :
              "CORRELATED_CANDIDATE_CONFIRMATION_PENDING"));
        return out;
    }

    void clear() { states_.clear(); }

private:
    struct State
    {
        bool observed = false;
        double lastObservedTime = 0;
        ZhangUserIntegerAuditSnapshot previous;
        std::map<std::string,std::string> latestIdentities;
        std::map<std::string,std::string> quarantined;
        std::string signature;
        int confirmationCount = 0;
    };
    std::map<std::string, State> states_;

    static ZhangUserIntegerAuditSnapshot reindex(const ZhangUserIntegerAuditSnapshot& previous,
        const std::vector<std::string>& columns)
    {
        ZhangUserIntegerAuditSnapshot out; out.columns = columns;
        if (previous.rows.empty()) return out;
        std::map<std::string,int> directory;
        for (int c = 0; c < columns.size(); ++c) directory[columns[c]] = c;
        std::vector<int> retained;
        for (int c = 0; c < previous.columns.size(); ++c)
            if (directory.count(previous.columns[c])) retained.push_back(c);
        const auto projected = zhangRestrictUserIntegerLedger(previous.rows,
            previous.values, retained, previous.columns.size());
        if (!projected.valid) throw std::runtime_error("CONTINUITY_AUDIT_REINDEX_FAILED:" + projected.reason);
        out.values = projected.values;
        for (const auto& row : projected.rows)
        {
            ZhangExactVector full(columns.size());
            for (int c = 0; c < retained.size(); ++c)
                full[directory.at(previous.columns[retained[c]])] = row[c];
            out.rows.push_back(std::move(full));
        }
        return out;
    }

    static std::string canonicalSignature(const ZhangUserIntegerAuditSnapshot& candidate)
    {
        if (candidate.rows.empty()) return {};
        std::vector<std::string> columns;
        std::vector<int> active;
        for (int c = 0; c < candidate.columns.size(); ++c)
            if (std::any_of(candidate.rows.begin(), candidate.rows.end(),
                [&](const auto& row) {return row[c] != 0;}))
            { columns.push_back(candidate.columns[c]); active.push_back(c); }
        std::sort(columns.begin(), columns.end());
        std::map<std::string,int> directory;
        for (int c = 0; c < columns.size(); ++c) directory[columns[c]] = c;
        ZhangExactMatrix rows(candidate.rows.size(), ZhangExactVector(columns.size()));
        for (int r = 0; r < rows.size(); ++r)
            for (int c : active)
                rows[r][directory.at(candidate.columns[c])] = candidate.rows[r][c];
        const auto hnf = zhangAuditUserIntegerUnion(rows, candidate.values, columns.size());
        if (!hnf.valid) throw std::runtime_error("CONTINUITY_SIGNATURE_INFEASIBLE");
        std::string signature;
        for (const auto& column : columns)
            signature += std::to_string(column.size()) + ":" + column;
        for (int r = 0; r < hnf.rows.size(); ++r)
        {
            for (const auto& value : hnf.rows[r]) signature += value.convert_to<std::string>() + ",";
            signature += "=" + hnf.values[r].convert_to<std::string>() + ";";
        }
        return signature;
    }
};
