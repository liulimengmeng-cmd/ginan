#pragma once

#include "common/zhangFullRank.hpp"
#include "common/zhangIntegerAudit.hpp"
#include <iomanip>
#include <sstream>

// Implemented by the checkpoint module; the validator is read-only.
std::string zhangCheckpointSha256(const std::string& bytes);

// This namespace deliberately has no conversion to a CURRENT AR recipe,
// conditioning equation, decision proof, or current product certificate.
// Archived RHSs authorize chart transport only, never present-epoch fixing.
namespace zhang_transport_certificate
{
using Row = std::map<std::string, ZhangExactInteger>;
inline constexpr const char* metadataKey = "ZHANG_TRANSPORT_ONLY_COMMITTED_V1";
inline constexpr std::size_t maximumBytes = 16ULL * 1024 * 1024;
inline constexpr std::size_t maximumFacts = 16384;

inline bool physicalIdentityMatchesSystem(const std::string& id, int system)
{
    const auto p1 = id.find('|'), p2 = id.find('|', p1 == std::string::npos ? 0 : p1 + 1);
    const auto p3 = id.find("|V", p2 == std::string::npos ? 0 : p2 + 1);
    if (p1 == std::string::npos || p1 == 0 || p2 == std::string::npos || p2 == p1 + 1 ||
        p3 == std::string::npos || p3 <= p2 + 1 || p3 + 2 == id.size() || id.find('|', p3 + 1) != std::string::npos)
        return false;
    for (std::size_t i = p3 + 2; i < id.size(); ++i) if (id[i] < '0' || id[i] > '9') return false;
    const char constellation = id[p2 + 1];
    switch (static_cast<E_Sys>(system))
    {
        case E_Sys::GPS: return constellation == 'G';
        case E_Sys::GAL: return constellation == 'E';
        case E_Sys::BDS: return constellation == 'C';
        case E_Sys::GLO: return constellation == 'R';
        case E_Sys::QZS: return constellation == 'J';
        case E_Sys::SBS: return constellation == 'S';
        case E_Sys::LEO: return constellation == 'L';
        default: return false;
    }
}

struct SourceScope
{
    long int committedEpoch = 0;
    std::string posterior, chart, floatRoot, stateOrder;
    bool operator==(const SourceScope&) const = default;
};

struct Fact
{
    Row physical;
    ZhangExactInteger value = 0;
    std::string sourceScopeId;
    std::vector<std::string> proofIds;
};

struct Directory
{
    std::string authoritativeOwner;
    int system = 0;
    long int lastCommittedEpoch = 0;
    std::vector<Fact> facts;
    // Intern full provenance/order and proof DAG payloads once. Full KF order
    // can be large: copying it into every sparse equation is forbidden.
    std::map<std::string, SourceScope> sourceScopes;
    std::map<std::string, std::string> proofAtoms;

    std::string encode() const
    {
        std::ostringstream out;
        out << "TRANSPORT_ONLY_COMMITTED_V1 " << std::quoted(authoritativeOwner)
            << ' ' << system << ' ' << lastCommittedEpoch << ' ' << facts.size()
            << ' ' << sourceScopes.size() << ' ' << proofAtoms.size() << '\n';
        for (const auto& [id, scope] : sourceScopes)
            out << std::quoted(id) << ' ' << scope.committedEpoch << ' ' << std::quoted(scope.posterior) << ' '
                << std::quoted(scope.chart) << ' ' << std::quoted(scope.floatRoot) << ' '
                << std::quoted(scope.stateOrder) << '\n';
        for (const auto& [id, atom] : proofAtoms)
            out << std::quoted(id) << ' ' << std::quoted(atom) << '\n';
        for (const auto& fact : facts)
        {
            out << std::quoted(fact.sourceScopeId) << ' ' << fact.value << ' '
                << fact.physical.size() << ' ' << fact.proofIds.size() << '\n';
            for (const auto& [id, coefficient] : fact.physical)
                out << std::quoted(id) << ' ' << coefficient << '\n';
            for (const auto& id : fact.proofIds) out << std::quoted(id) << '\n';
        }
        return out.str();
    }

    static bool decode(const std::string& bytes, const std::string& expectedOwner,
        int expectedSystem, Directory& result, std::string& reason)
    {
        if (bytes.size() > maximumBytes)
        { reason = "TRANSPORT_DIRECTORY_BYTE_LIMIT"; return false; }
        std::istringstream in(bytes);
        Directory candidate;
        std::string schema;
        std::map<Row, ZhangExactInteger> seenRows;
        std::size_t count = 0, scopes = 0, proofCount = 0;
        if (!(in >> schema >> std::quoted(candidate.authoritativeOwner) >> candidate.system
              >> candidate.lastCommittedEpoch >> count >> scopes >> proofCount) || schema != "TRANSPORT_ONLY_COMMITTED_V1" ||
            candidate.authoritativeOwner != expectedOwner || candidate.system != expectedSystem ||
            count > maximumFacts || scopes > maximumFacts || proofCount > 65536)
        { reason = "TRANSPORT_DIRECTORY_OWNER_OR_SCHEMA_MISMATCH"; return false; }
        for (std::size_t i = 0; i < scopes; ++i)
        {
            SourceScope scope;
            std::string id;
            if (!(in >> std::quoted(id) >> scope.committedEpoch >> std::quoted(scope.posterior) >>
                std::quoted(scope.chart) >> std::quoted(scope.floatRoot) >> std::quoted(scope.stateOrder)) ||
                id.empty() || scope.committedEpoch > candidate.lastCommittedEpoch || scope.posterior.empty() ||
                scope.chart.empty() || scope.floatRoot.empty() || scope.stateOrder.empty() ||
                !candidate.sourceScopes.emplace(id, std::move(scope)).second)
            { reason = "TRANSPORT_SOURCE_SCOPE_INVALID"; return false; }
        }
        for (std::size_t i = 0; i < proofCount; ++i)
        {
            std::string id, payload;
            if (!(in >> std::quoted(id) >> std::quoted(payload)) || id.empty() || payload.empty() ||
                !candidate.proofAtoms.emplace(id, payload).second)
            { reason = "TRANSPORT_PROOF_PAYLOAD_INVALID_OR_DUPLICATE"; return false; }
        }
        for (std::size_t index = 0; index < count; ++index)
        {
            Fact fact;
            std::size_t terms = 0, atoms = 0;
            if (!(in >> std::quoted(fact.sourceScopeId) >> fact.value >> terms >> atoms) ||
                terms == 0 || terms > 65536 || atoms == 0 || atoms > 65536 ||
                !candidate.sourceScopes.contains(fact.sourceScopeId))
            { reason = "TRANSPORT_FACT_TRUNCATED_OR_UNBOUND"; return false; }
            for (std::size_t term = 0; term < terms; ++term)
            {
                std::string id;
                ZhangExactInteger coefficient;
                if (!(in >> std::quoted(id) >> coefficient) || !physicalIdentityMatchesSystem(id, expectedSystem) || coefficient == 0 ||
                    !fact.physical.emplace(id, coefficient).second)
                { reason = "TRANSPORT_PHYSICAL_IDENTITY_INVALID"; return false; }
            }
            for (std::size_t atom = 0; atom < atoms; ++atom)
            {
                std::string id;
                if (!(in >> std::quoted(id)) || !candidate.proofAtoms.contains(id) ||
                    std::find(fact.proofIds.begin(), fact.proofIds.end(), id) != fact.proofIds.end())
                { reason = "TRANSPORT_PROOF_PAYLOAD_INVALID"; return false; }
                fact.proofIds.push_back(std::move(id));
            }
            const auto [row, inserted] = seenRows.emplace(fact.physical, fact.value);
            if (!inserted)
            { reason = row->second == fact.value ? "TRANSPORT_DUPLICATE_ROW" : "TRANSPORT_ROW_VALUE_CONFLICT"; return false; }
            candidate.facts.push_back(std::move(fact));
        }
        in >> std::ws;
        if (!in.eof()) { reason = "TRANSPORT_DIRECTORY_TRAILING_DATA"; return false; }
        result = std::move(candidate);
        reason = "TRANSPORT_DIRECTORY_OWNER_BOUND";
        return true;
    }
};

// Every row refers to its OWN signal/version identities. WL-only information
// can be a known joint direction, while both scalar L1 and L2 remain unknown.
// Empty/unknown coordinates are never supplied as zeros to another consumer.
struct DirectionContract
{
    std::string owner, sourceEpoch, targetEpoch, reason = "NOT_EVALUATED";
    std::vector<Row> requestedDirections;
    std::vector<Row> knownDirections;
    ZhangExactVector knownValues;
    std::vector<Row> unknownDirections;
    // Rows act on the ORDERED requested target vector, not on physical KF
    // coordinates. They can certify WL while neither scalar phase is known.
    ZhangExactMatrix knownTargetCombinations;
    ZhangExactVector knownCombinationValues;
    std::size_t knownRank = 0, unknownRank = 0;
    bool ranksEvaluated = false;
    bool currentArAuthority = false;
};

// One immutable build for all candidates of a graph event. Integer feasibility
// is checked even though the inputs carry accepted-decision provenance.
struct Query
{
    struct UnitEquation { std::string pivot; Row row; ZhangExactInteger value; };
    std::vector<UnitEquation> units;
    std::map<std::string, std::size_t> columns;
    ZhangExactMatrix basis;
    ZhangExactVector values;
    std::vector<std::size_t> pivots;
    bool feasible = false;
    std::string reason = "NO_COMMITTED_TRANSPORT_CERTIFICATES";

    void reduceUnits(Row& row, ZhangExactInteger& constant) const
    {
        for (const auto& unit : units)
        {
            const auto found = row.find(unit.pivot);
            if (found == row.end()) continue;
            const ZhangExactInteger multiplier = found->second;
            for (const auto& [id, coefficient] : unit.row)
            {
                auto& remaining = row[id];
                remaining -= multiplier * coefficient;
                if (remaining == 0) row.erase(id);
            }
            constant += multiplier * unit.value;
        }
    }

    explicit Query(const Directory& directory)
    {
        for (const auto& fact : directory.facts)
            for (const auto& [id, coefficient] : fact.physical)
                if (coefficient != 0 && !columns.contains(id)) columns.emplace(id, columns.size());
        if (columns.empty()) return;
        // Physical cycle equations usually have sparse unit pivots. Exact
        // affine substitution eliminates them before any dense integer
        // decomposition; a coefficient 2 is NEVER divided as if it were a
        // unit, and its congruence survives in the residual integer backend.
        std::vector<std::pair<Row, ZhangExactInteger>> pending;
        for (const auto& fact : directory.facts) pending.emplace_back(fact.physical, fact.value);
        while (!pending.empty())
        {
            const auto before = units.size();
            std::vector<std::pair<Row, ZhangExactInteger>> residual;
            for (auto& [row, rhs] : pending)
            {
                ZhangExactInteger removed = 0;
                reduceUnits(row, removed); rhs -= removed;
                if (row.empty())
                {
                    if (rhs != 0) { reason = "TRANSPORT_HISTORY_INCONSISTENT"; return; }
                    continue;
                }
                const auto unit = std::find_if(row.begin(), row.end(),
                    [](const auto& term) { return term.second == 1 || term.second == -1; });
                if (unit == row.end()) { residual.emplace_back(std::move(row), std::move(rhs)); continue; }
                const auto pivot = unit->first;
                if (unit->second == -1) { for (auto& [id, coefficient] : row) coefficient = -coefficient; rhs = -rhs; }
                units.push_back({pivot, std::move(row), std::move(rhs)});
            }
            pending = std::move(residual);
            if (units.size() == before) break;
        }
        if (pending.empty()) { feasible = true; reason = "TRANSPORT_HISTORY_SPARSE_UNIT_FEASIBLE"; return; }
        std::set<std::size_t> residualPositions;
        for (const auto& [row, rhs] : pending)
            for (const auto& [id, coefficient] : row) residualPositions.insert(columns.at(id));
        std::vector<std::size_t> positions(residualPositions.begin(), residualPositions.end());
        ZhangExactMatrix rows(pending.size(), ZhangExactVector(positions.size()));
        ZhangExactVector rhs;
        for (std::size_t r = 0; r < pending.size(); ++r)
        {
            for (const auto& [id, coefficient] : pending[r].first)
                rows[r][std::lower_bound(positions.begin(), positions.end(), columns.at(id)) - positions.begin()] = coefficient;
            rhs.push_back(pending[r].second);
        }
        auto hnf = zhangExactRowHermiteNormalForm(std::move(rows), rhs);
        if (!hnf.consistent) { reason = "TRANSPORT_HISTORY_INCONSISTENT"; return; }
        ZhangExactMatrix transpose(positions.size(), ZhangExactVector(hnf.basis.size()));
        for (std::size_t r = 0; r < hnf.basis.size(); ++r)
            for (std::size_t c = 0; c < positions.size(); ++c) transpose[c][r] = hnf.basis[r][c];
        const auto origin = zhangIntegerRowLatticeContains(transpose, hnf.values, true, true);
        if (!origin.contained || origin.combination.size() != positions.size() ||
            zhangExactMatrixTimesColumn(hnf.basis, origin.combination) != hnf.values)
        { reason = "TRANSPORT_HISTORY_INTEGER_INFEASIBLE"; return; }
        basis.resize(hnf.basis.size(), ZhangExactVector(columns.size()));
        for (std::size_t r = 0; r < hnf.basis.size(); ++r)
            for (std::size_t c = 0; c < positions.size(); ++c) basis[r][positions[c]] = std::move(hnf.basis[r][c]);
        values = std::move(hnf.values);
        for (const auto& row : basis)
            pivots.push_back(std::distance(row.begin(), std::find_if(row.begin(), row.end(),
                [](const auto& coefficient) { return coefficient != 0; })));
        feasible = true;
        reason = "TRANSPORT_HISTORY_INTEGER_FEASIBLE";
    }

    bool consequence(const Row& target, ZhangExactInteger& value) const
    {
        if (!feasible || target.empty()) return false;
        using Rational = boost::multiprecision::cpp_rational;
        Row reduced = target;
        ZhangExactInteger unitConstant = 0;
        reduceUnits(reduced, unitConstant);
        std::vector<Rational> remaining(columns.size());
        for (const auto& [id, coefficient] : reduced)
        {
            if (coefficient == 0) continue;
            const auto found = columns.find(id);
            if (found == columns.end()) return false;
            remaining[found->second] = Rational(coefficient);
        }
        Rational rhs = Rational(unitConstant);
        for (std::size_t r = 0; r < basis.size(); ++r)
        {
            const auto p = pivots[r];
            const Rational multiplier = remaining[p] / Rational(basis[r][p]);
            if (multiplier == 0) continue;
            for (std::size_t c = p; c < remaining.size(); ++c)
                remaining[c] -= multiplier * Rational(basis[r][c]);
            rhs += multiplier * Rational(values[r]);
        }
        if (std::any_of(remaining.begin(), remaining.end(), [](const auto& c) { return c != 0; }) ||
            denominator(rhs) != 1) return false;
        value = numerator(rhs);
        return true;
    }

    DirectionContract classify(const std::vector<Row>& targets) const
    {
        DirectionContract result;
        result.requestedDirections = targets;
        std::vector<std::size_t> unknownIndices;
        for (std::size_t index = 0; index < targets.size(); ++index)
        {
            const auto& target = targets[index];
            ZhangExactInteger value;
            if (consequence(target, value))
            {
                result.knownDirections.push_back(target); result.knownValues.push_back(value);
                ZhangExactVector unit(targets.size()); unit[index] = 1;
                result.knownTargetCombinations.push_back(std::move(unit));
                result.knownCombinationValues.push_back(value);
            }
            else { result.unknownDirections.push_back(target); unknownIndices.push_back(index); }
        }
        if (feasible && !targets.empty())
        {
            result.ranksEvaluated = true;
            using Rational = boost::multiprecision::cpp_rational;
            auto allColumns = columns;
            for (const auto& target : result.unknownDirections)
                for (const auto& [id, coefficient] : target)
                    if (coefficient != 0 && !allColumns.contains(id)) allColumns.emplace(id, allColumns.size());
            std::vector<std::vector<Rational>> residual(unknownIndices.size(), std::vector<Rational>(allColumns.size()));
            for (std::size_t t = 0; t < unknownIndices.size(); ++t)
            {
                Row reduced = targets[unknownIndices[t]];
                ZhangExactInteger unusedConstant = 0;
                reduceUnits(reduced, unusedConstant);
                for (const auto& [id, coefficient] : reduced) residual[t][allColumns.at(id)] = Rational(coefficient);
                for (std::size_t r = 0; r < basis.size(); ++r)
                {
                    const Rational multiplier = residual[t][pivots[r]] / Rational(basis[r][pivots[r]]);
                    for (std::size_t c = pivots[r]; c < columns.size(); ++c)
                        residual[t][c] -= multiplier * Rational(basis[r][c]);
                }
            }
            ZhangExactMatrix transpose;
            for (std::size_t c = 0; c < allColumns.size(); ++c)
            {
                ZhangExactInteger denominatorProduct = 1;
                for (const auto& row : residual) denominatorProduct *= denominator(row[c]);
                ZhangExactVector row(unknownIndices.size());
                bool nonzero = false;
                for (std::size_t t = 0; t < unknownIndices.size(); ++t)
                {
                    row[t] = numerator(residual[t][c]) * (denominatorProduct / denominator(residual[t][c]));
                    nonzero |= row[t] != 0;
                }
                if (nonzero) transpose.push_back(std::move(row));
            }
            const auto residualKernel = zhangExactIntegerKernel(std::move(transpose), unknownIndices.size());
            for (const auto& residualCombination : residualKernel)
            {
                ZhangExactVector combination(targets.size());
                for (std::size_t u = 0; u < unknownIndices.size(); ++u)
                    combination[unknownIndices[u]] = residualCombination[u];
                Row combined;
                for (std::size_t t = 0; t < targets.size(); ++t)
                    for (const auto& [id, coefficient] : targets[t]) combined[id] += combination[t] * coefficient;
                for (auto it = combined.begin(); it != combined.end();)
                    if (it->second == 0) it = combined.erase(it); else ++it;
                ZhangExactInteger constant = 0;
                if (!combined.empty() && !consequence(combined, constant))
                {
                    result.ranksEvaluated = false;
                    result.knownTargetCombinations.clear();
                    result.knownCombinationValues.clear();
                    result.reason = "TRANSPORT_DIRECTION_COMBINATION_RECHECK_FAILED";
                    return result;
                }
                result.knownTargetCombinations.push_back(std::move(combination));
                result.knownCombinationValues.push_back(constant);
            }
            result.knownRank = result.knownTargetCombinations.size();
            result.unknownRank = targets.size() - result.knownRank;
        }
        result.reason = result.unknownRank == 0 && feasible
            ? "ALL_REQUESTED_TRANSPORT_DIRECTIONS_CERTIFIED" : "PARTIAL_TRANSPORT_DIRECTIONS_UNKNOWN";
        return result;
    }
};

struct TreePlan
{
    ZhangGraphBasis basis;
    std::set<ZhangGraphEdge> certifiedReplacements, fallbackReplacements, newEndpointExtensions;
    std::map<int, DirectionContract> perSignal;
    DirectionContract jointSignalDirections;
    bool allOldEndpointsCertified = false;
};

// Admission is lexicographic: old continuous forest, then edges certified in
// ALL required signals, then the legal FLOAT fallback. Quality only orders
// edges within one class. No certificate fabricates an old endpoint for a leaf.
template<class Certify>
TreePlan select(const ZhangGraphBasis& oldBasis, const std::set<ZhangGraphEdge>& available,
    const std::set<ZhangGraphEdge>& broken, const std::string& root,
    const std::map<ZhangGraphEdge, double>& quality, Certify certify)
{
    std::set<ZhangGraphEdge> retained, certified, legal;
    for (const auto& edge : available) if (!broken.contains(edge)) legal.insert(edge);
    for (const auto& edge : oldBasis.treeEdges)
        if (available.contains(edge) && !broken.contains(edge)) retained.insert(edge);
    zhang_graph_detail::DisjointSet certifiedComponents;
    for (const auto& edge : retained)
        certifiedComponents.unite(zhang_graph_detail::receiverNode(edge.receiver),
            zhang_graph_detail::satelliteNode(edge.satellite));
    std::vector<ZhangGraphEdge> ordered(legal.begin(), legal.end());
    std::sort(ordered.begin(), ordered.end(), [&](const auto& a, const auto& b)
    {
        const double qa = quality.contains(a) ? quality.at(a) : 0;
        const double qb = quality.contains(b) ? quality.at(b) : 0;
        return qa != qb ? qa > qb : a < b;
    });
    for (const auto& edge : ordered)
        if (!retained.contains(edge) && !broken.contains(edge) && oldBasis.edges.contains(edge) &&
            oldBasis.receivers.contains(edge.receiver) && oldBasis.satellites.contains(edge.satellite) &&
            certifiedComponents.find(zhang_graph_detail::receiverNode(edge.receiver)) !=
                certifiedComponents.find(zhang_graph_detail::satelliteNode(edge.satellite)) && certify(edge))
        {
            certified.insert(edge);
            certifiedComponents.unite(zhang_graph_detail::receiverNode(edge.receiver),
                zhang_graph_detail::satelliteNode(edge.satellite));
        }
    auto admitted = retained;
    admitted.insert(certified.begin(), certified.end());
    const auto strict = zhangBuildSpanningTree(admitted, root, retained, quality);
    TreePlan result;
    result.allOldEndpointsCertified = !oldBasis.receivers.empty() &&
        strict.connected && strict.receivers == oldBasis.receivers && strict.satellites == oldBasis.satellites;
    // Passing strict.treeEdges as preferred enforces an admitted forest;
    // unproved edges cannot displace it by a better elevation score.
    result.basis = zhangBuildSpanningTree(legal, root, strict.treeEdges, quality);
    for (const auto& edge : result.basis.treeEdges)
        if (!retained.contains(edge))
        {
            if (!oldBasis.receivers.contains(edge.receiver) || !oldBasis.satellites.contains(edge.satellite))
                result.newEndpointExtensions.insert(edge);
            else (certified.contains(edge) ? result.certifiedReplacements : result.fallbackReplacements).insert(edge);
        }
    return result;
}
}

/** Checkpoint validation only: no graph lookup, no live-column requirement,
 * no historical RHS application, and no mutation. Absence is valid. */
inline bool zhangValidateCertifiedTreeTransportMetadata(
    const std::map<std::string, std::string>& metadata, const std::string& expectedOwner,
    std::string& reason)
{
    using namespace zhang_transport_certificate;
    const std::string prefix = std::string(metadataKey) + "/";
    for (const auto& [key, envelope] : metadata)
    {
        if (key.rfind(metadataKey, 0) != 0) continue;
        if (key.rfind(prefix, 0) != 0 || envelope.size() < 65 || envelope[64] != '\n' ||
            zhangCheckpointSha256(envelope.substr(65)) != envelope.substr(0, 64))
        { reason = "TRANSPORT_CHECKPOINT_KEY_OR_SEAL_INVALID"; return false; }
        std::string schema, owner;
        int system = 0;
        std::istringstream header(envelope.substr(65));
        if (!(header >> schema >> std::quoted(owner) >> system) || system <= static_cast<int>(E_Sys::NONE) ||
            system >= static_cast<int>(E_Sys::SUPPORTED) ||
            key != prefix + enum_to_string(static_cast<E_Sys>(system)))
        { reason = "TRANSPORT_CHECKPOINT_SYSTEM_MISMATCH"; return false; }
        Directory directory;
        if (!Directory::decode(envelope.substr(65), expectedOwner, system, directory, reason)) return false;
    }
    reason = "TRANSPORT_CHECKPOINT_METADATA_VALID_OR_ABSENT";
    return true;
}
