#pragma once
#include <Eigen/Core>

#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include "common/zhangIntegerAudit.hpp"

/** One named satellite-minus-reference product relation, independently
 * expanded from current fundamental cycles to physical ambiguity arcs.
 *
 * The representation deliberately contains no ZhangProductIntegerFunctional:
 * it is structural evidence derived from the graph/S-system itself, not a
 * self-attestation by the product coordinate that will later consume it.
 */
struct ZhangProductRelationRow
{
    SatSys                                  satellite;
    SatSys                                  referenceSatellite;
    ZhangExactVector                        currentCycleCoefficients;
    std::map<ZhangGraphEdge, ZhangExactInteger> physicalArcCoefficients;
    ZhangExactVector                        nuisanceCoefficients;
};

enum class ZhangCanonicalProductDirectionStatus
{
    SPATIALLY_MAPPABLE,
    REQUIRES_COMPONENT_GAUGE,
    REQUIRES_BESD,
    PHYSICALLY_UNSUPPORTED
};

inline const char* zhangCanonicalProductDirectionStatusName(
    ZhangCanonicalProductDirectionStatus status)
{
    switch (status)
    {
        case ZhangCanonicalProductDirectionStatus::SPATIALLY_MAPPABLE:
            return "SPATIALLY_MAPPABLE";
        case ZhangCanonicalProductDirectionStatus::REQUIRES_COMPONENT_GAUGE:
            return "REQUIRES_COMPONENT_GAUGE";
        case ZhangCanonicalProductDirectionStatus::REQUIRES_BESD:
            return "REQUIRES_BESD";
        case ZhangCanonicalProductDirectionStatus::PHYSICALLY_UNSUPPORTED:
            return "PHYSICALLY_UNSUPPORTED";
    }
    return "PHYSICALLY_UNSUPPORTED";
}

struct ZhangUnmappableProductRelationAudit
{
    SatSys          satellite;
    SatSys          referenceSatellite;
    E_ObsCode       observable = E_ObsCode::NONE;
    ZhangGraphEdge  missingChord;
    int             missingArcVersion = -1;
    std::string     missingReason = "NOT_EVALUATED";
    std::string     treeOrChord = "UNKNOWN";
    bool            temporalRecoverable = false;
    int             alternativeReceiverCount = 0;
    ZhangCanonicalProductDirectionStatus status =
        ZhangCanonicalProductDirectionStatus::PHYSICALLY_UNSUPPORTED;
};

struct ZhangProductRelationBasis
{
    E_Sys                                   system = E_Sys::NONE;
    E_ObsCode                               observable = E_ObsCode::NONE;
    SatSys                                  referenceSatellite;
    std::vector<SatSys>                     satellites;
    std::vector<ZhangGraphEdge>             currentChords;
    std::vector<ZhangProductRelationRow>     namedRelations;
    // Complete canonical satellite-product target in physical ambiguity-arc
    // coordinates.  Unlike namedRelations, these rows are not current-cycle
    // corrections and therefore are never truncated by a product core.
    std::vector<ZhangProductRelationRow>     canonicalPhysicalRelations;
    std::vector<int>                        canonicalPhysicalNamedIndices;
    std::vector<int>                        independentNamedIndices;
    std::vector<int>                        mappableNamedIndices;
    std::vector<ZhangUnmappableProductRelationAudit> unmappableNamedRelations;
    std::vector<ZhangGraphEdge>             physicalArcColumns;
    ZhangExactMatrix                        networkIntegerBasis;
    ZhangExactMatrix                        canonicalTargetBasis;
    std::vector<ZhangExactInteger>          canonicalSmithInvariants;
    std::string                             canonicalTargetHnf;
    int                                     canonicalTargetRank = 0;
    bool                                    canonicalTargetPrimitive = false;
    ZhangExactMatrix                        exactRowBasis;
    ZhangExactMatrix                        networkContainmentTransform;
    Eigen::MatrixXd                         transform;
    ZhangExactVector                        affineOffsets;
    int                                     fullTargetRank = 0;
    int                                     mappableTargetRank = 0;
    int                                     primitiveRank = 0;
    int                                     unmappableTargetRank = 0;
    int                                     namedRelationCount = 0;
    int                                     exactRank = 0;
    ZhangExactInteger                       saturationIndex = 0;
    std::string                             exactHnf;
    bool                                    primitive = false;
    bool                                    admissibleCompletionProven = false;
    bool                                    networkLatticeContained = false;
    bool                                    networkClosureExactZero = false;
    bool                                    temporalRecoveryRequired = false;
    bool                                    nuisanceOrthogonal = false;
    bool                                    physicalExpansionValid = false;
    // The authoritative graph/product trees remain untouched.  A private
    // product-search basis may be rebuilt from ambiguity arcs that actually
    // have current posterior columns when that exact reparameterisation
    // exposes a larger canonical product rank.
    bool                                    availabilityRebased = false;
    bool                                    productTreeTransported = false;
    std::vector<ZhangGraphEdge>             staleProductTreeEdges;
    int                                     legacyFilteredMappableTargetRank = 0;
    int                                     authoritativeMappableTargetRank = 0;
    int                                     posteriorColumnEdgeCount = 0;
    int                                     availableStateEdgeCount = 0;
    int                                     ignoredStalePosteriorEdgeCount = 0;
    std::string                             compilationRootReceiver;
    bool                                    valid = false;
    std::string                             failureReason;
};

/** Compare only the physical satellite-minus-reference semantics of two
 * ordered relation catalogues.  Numeric row positions, current-cycle
 * coefficients and physical preimages are representation-specific and must
 * not be used to claim that L1 and L2 describe the same canonical target. */
inline bool zhangProductRelationSemanticOrderingMatches(
    const std::vector<ZhangProductRelationRow>& first,
    const std::vector<ZhangProductRelationRow>& second)
{
    if (first.size() != second.size()) return false;
    for (std::size_t index = 0; index < first.size(); index++)
    {
        if (first[index].satellite != second[index].satellite ||
            first[index].referenceSatellite !=
                second[index].referenceSatellite)
        {
            return false;
        }
    }
    return true;
}

/** Select an exact independent named product subset after posterior-column
 * availability is known.
 *
 * The builder's independentNamedIndices are independent in the complete
 * structural graph.  Filtering that fixed subset afterwards is incorrect: a
 * row that was redundant in the full graph can replace an unavailable pivot
 * row and increase the rank visible in the current posterior. */
inline std::vector<int> zhangIndependentMappableProductRelationIndices(
    const ZhangProductRelationBasis& basis,
    const std::set<ZhangGraphEdge>&  availableEdges)
{
    std::vector<int> indices;
    ZhangExactMatrix independentRows;
    int rank = 0;
    for (int named = 0;
         named < static_cast<int>(basis.namedRelations.size()); named++)
    {
        const auto& row = basis.namedRelations[named].currentCycleCoefficients;
        if (row.size() != basis.currentChords.size()) continue;
        bool mappable = true;
        for (int chord = 0;
             chord < static_cast<int>(basis.currentChords.size()); chord++)
        {
            if (row[chord] != 0 &&
                !availableEdges.contains(basis.currentChords[chord]))
            {
                mappable = false;
                break;
            }
        }
        if (!mappable) continue;
        ZhangExactMatrix candidate = independentRows;
        candidate.push_back(row);
        const int candidateRank = static_cast<int>(
            zhangExactRowHermiteNormalForm(candidate).basis.size());
        if (candidateRank <= rank) continue;
        indices.push_back(named);
        independentRows.push_back(row);
        rank = candidateRank;
    }
    return indices;
}

/** Rank exposed by the legacy "select globally, then filter" policy.
 *
 * independentNamedIndices is an exact independent set in the complete graph,
 * so every surviving subset remains independent.  This diagnostic therefore
 * counts the globally selected rows whose non-zero current chords all have a
 * posterior column.  Comparing it with the post-availability reselection rank
 * proves whether a redundant global row supplied a missing local pivot. */
inline int zhangLegacyFilteredMappableProductRelationRank(
    const ZhangProductRelationBasis& basis,
    const std::set<ZhangGraphEdge>&  availableEdges)
{
    int rank = 0;
    for (const int named : basis.independentNamedIndices)
    {
        if (named < 0 || named >= static_cast<int>(basis.namedRelations.size()))
            continue;
        const auto& row = basis.namedRelations[named].currentCycleCoefficients;
        if (row.size() != basis.currentChords.size()) continue;
        bool mappable = true;
        for (int chord = 0;
             chord < static_cast<int>(basis.currentChords.size()); chord++)
        {
            if (row[chord] != 0 &&
                !availableEdges.contains(basis.currentChords[chord]))
            {
                mappable = false;
                break;
            }
        }
        rank += mappable;
    }
    return rank;
}

/** Construct one deterministic physical arc preimage for every canonical
 * satellite-minus-reference direction.  A path alternates +1/-1 from the
 * target satellite to the reference satellite, so every receiver and every
 * intermediate satellite cancels while the endpoint satellite divergence is
 * exactly e_s-e_ref.  Different legal paths are alternative representations
 * of the same canonical q and are compared modulo the certified closure
 * lattice downstream. */
inline bool zhangBuildCanonicalProductPhysicalRelations(
    const ZhangGraphBasis& representedBasis,
    const std::set<SatSys>& canonicalSatellites,
    const SatSys& requestedReference,
    std::vector<ZhangProductRelationRow>& rows,
    std::string& failureReason)
{
    rows.clear();
    if (canonicalSatellites.empty())
    {
        failureReason = "EMPTY_CANONICAL_SATELLITE_TARGET";
        return false;
    }
    SatSys reference = requestedReference;
    if (reference.prn == 0) reference = *canonicalSatellites.begin();
    if (!canonicalSatellites.contains(reference))
    {
        failureReason = "CANONICAL_REFERENCE_OUTSIDE_TARGET";
        return false;
    }
    using Step = std::pair<std::string, ZhangGraphEdge>;
    std::map<std::string, std::vector<Step>> adjacency;
    auto receiverNode = [](const std::string& receiver)
        { return std::string("R:") + receiver; };
    auto satelliteNode = [](const SatSys& satellite)
        {
            return std::string("S:") +
                std::to_string(static_cast<int>(satellite.sys)) + ":" +
                std::to_string(satellite.prn);
        };
    for (const auto& edge : representedBasis.edges)
    {
        const auto receiver = receiverNode(edge.receiver);
        const auto satellite = satelliteNode(edge.satellite);
        adjacency[receiver].push_back({satellite, edge});
        adjacency[satellite].push_back({receiver, edge});
    }
    for (auto& [node, neighbours] : adjacency)
        std::sort(neighbours.begin(), neighbours.end(),
            [](const auto& left, const auto& right)
            {
                if (left.first != right.first) return left.first < right.first;
                return left.second < right.second;
            });
    const std::string destination = satelliteNode(reference);
    bool complete = true;
    for (const auto& satellite : canonicalSatellites)
    {
        if (satellite == reference) continue;
        ZhangProductRelationRow relation;
        relation.satellite = satellite;
        relation.referenceSatellite = reference;
        const std::string source = satelliteNode(satellite);
        std::queue<std::string> pending;
        std::set<std::string> visited{source};
        std::map<std::string, Step> parent;
        pending.push(source);
        while (!pending.empty() && !visited.contains(destination))
        {
            const auto node = pending.front();
            pending.pop();
            for (const auto& next : adjacency[node])
            {
                if (!visited.insert(next.first).second) continue;
                parent[next.first] = {node, next.second};
                pending.push(next.first);
            }
        }
        if (!visited.contains(destination))
        {
            // Retain the named canonical coordinate even when no current
            // spatial preimage exists.  Downstream exact-intersection code
            // excludes this empty physical row and embeds the supported
            // result back into the complete canonical coordinate space.
            rows.push_back(std::move(relation));
            complete = false;
            continue;
        }
        std::vector<ZhangGraphEdge> reversePath;
        std::string node = destination;
        while (node != source)
        {
            auto found = parent.find(node);
            if (found == parent.end())
            {
                failureReason = "CANONICAL_PHYSICAL_PATH_RECONSTRUCTION_FAILED";
                rows.clear();
                return false;
            }
            reversePath.push_back(found->second.second);
            node = found->second.first;
        }
        std::reverse(reversePath.begin(), reversePath.end());
        for (size_t edge = 0; edge < reversePath.size(); edge++)
            relation.physicalArcCoefficients[reversePath[edge]] +=
                edge % 2 == 0 ? 1 : -1;
        rows.push_back(std::move(relation));
    }
    failureReason = complete
        ? "NONE" : "PARTIAL_CANONICAL_PHYSICAL_MAPPING";
    return true;
}

/** Build the product-relevant satellite integer lattice without estimator or
 * ProductRelationManager state.
 *
 * Each named row is first obtained from the exact current-graph/product-tree
 * incidence map.  It is then independently expanded through the current
 * fundamental cycles to original receiver-satellite ambiguity arcs.  The
 * nuisance block is explicit and identically zero.  Exact row HNF and Smith
 * invariants prove rank and primitivity; any non-primitive or malformed
 * result fails closed instead of being silently treated as an integer basis.
 */
struct ProductRelationBasisBuilder
{
    static ZhangProductRelationBasis build(
        const ZhangGraphBasis& currentBasis,
        const ZhangGraphBasis& productBasis,
        const SatSys& requestedReference = SatSys(),
        E_Sys system = E_Sys::NONE,
        E_ObsCode observable = E_ObsCode::NONE)
    {
        ZhangProductRelationBasis result;
        result.system = system;
        result.observable = observable;
        const ZhangSatelliteProductTarget target =
            zhangBuildSatelliteProductTarget(
                currentBasis, productBasis, requestedReference);
        if (!target.valid)
        {
            result.failureReason = target.failureReason;
            return result;
        }

        result.referenceSatellite = target.referenceSatellite;
        result.satellites = target.targetSatellites;
        result.satellites.push_back(target.referenceSatellite);
        std::sort(result.satellites.begin(), result.satellites.end());
        result.satellites.erase(
            std::unique(result.satellites.begin(), result.satellites.end()),
            result.satellites.end());
        result.currentChords = target.currentChords;
        result.canonicalTargetBasis = target.canonicalMatrix;
        result.canonicalSmithInvariants = target.canonicalSmithInvariants;
        result.canonicalTargetHnf = target.canonicalHnf;
        result.canonicalTargetRank = target.canonicalRank;
        result.canonicalTargetPrimitive = target.canonicalPrimitive;
        result.productTreeTransported = target.productTreeTransported;
        result.staleProductTreeEdges = target.staleProductTreeEdges;
        result.namedRelationCount = target.matrix.size();
        // fullTargetRank denotes the physical product target, not the
        // transient rank that happens to pull back through today's chord
        // columns.  Mappability is reported separately below.
        result.fullTargetRank = result.canonicalTargetRank;
        if (target.matrix.size() != target.targetSatellites.size() ||
            target.matrix.empty())
        {
            result.failureReason = "PRODUCT_RELATION_TARGET_DIMENSION_MISMATCH";
            return result;
        }

        std::map<std::string, std::size_t> receiverIndex;
        std::map<SatSys, std::size_t> satelliteIndex;
        std::size_t nuisanceIndex = 0;
        for (const auto& receiver : currentBasis.receivers)
        {
            receiverIndex[receiver] = nuisanceIndex++;
        }
        for (const auto& satellite : currentBasis.satellites)
        {
            satelliteIndex[satellite] = nuisanceIndex++;
        }
        result.physicalExpansionValid = true;
        result.nuisanceOrthogonal = true;
        for (std::size_t row = 0; row < target.matrix.size(); row++)
        {
            if (target.matrix[row].size() != target.currentChords.size())
            {
                result.failureReason =
                    "PRODUCT_RELATION_CYCLE_DIMENSION_MISMATCH";
                result.physicalExpansionValid = false;
                return result;
            }
            ZhangProductRelationRow relation;
            relation.satellite = target.targetSatellites[row];
            relation.referenceSatellite = target.referenceSatellite;
            relation.currentCycleCoefficients = target.matrix[row];
            relation.nuisanceCoefficients = ZhangExactVector(nuisanceIndex);

            for (std::size_t chord = 0;
                 chord < target.currentChords.size(); chord++)
            {
                const ZhangExactInteger multiplier = target.matrix[row][chord];
                if (multiplier == 0)
                {
                    continue;
                }
                const auto cycle = zhangFundamentalCycle(
                    currentBasis, target.currentChords[chord]);
                if (cycle.empty())
                {
                    result.failureReason =
                        "PRODUCT_RELATION_PHYSICAL_CYCLE_EXPANSION_FAILED";
                    result.physicalExpansionValid = false;
                    return result;
                }
                for (const auto& [edge, coefficient] : cycle)
                {
                    relation.physicalArcCoefficients[edge] +=
                        multiplier * coefficient;
                }
            }
            for (auto iterator = relation.physicalArcCoefficients.begin();
                 iterator != relation.physicalArcCoefficients.end();)
            {
                if (iterator->second == 0)
                {
                    iterator = relation.physicalArcCoefficients.erase(iterator);
                }
                else
                {
                    const auto& edge = iterator->first;
                    const auto& coefficient = iterator->second;
                    // A physical ambiguity arc carries one additive receiver
                    // and one additive satellite nuisance datum.  Exact cycle
                    // relations must annihilate both incidences.  Compute the
                    // coefficients from the expanded physical row instead of
                    // declaring them zero by construction.
                    relation.nuisanceCoefficients[
                        receiverIndex.at(edge.receiver)] += coefficient;
                    relation.nuisanceCoefficients[
                        satelliteIndex.at(edge.satellite)] += coefficient;
                    ++iterator;
                }
            }
            result.nuisanceOrthogonal &= std::all_of(
                relation.nuisanceCoefficients.begin(),
                relation.nuisanceCoefficients.end(),
                [](const auto& coefficient) { return coefficient == 0; });
            result.namedRelations.push_back(std::move(relation));
        }

        // Build the physical ambient from every exact row that will be
        // compared.  After an exact pivot, a legal network chord may be
        // represented by the current tree even when it is not present in the
        // instantaneous observation-edge set.  Using currentBasis.edges alone
        // therefore creates a representation-only
        // NETWORK_INTEGER_BASIS_ARC_OUTSIDE_AMBIENT_SPACE failure.
        std::set<ZhangGraphEdge> physicalAmbient = currentBasis.edges;
        for (const auto& relation : result.namedRelations)
            for (const auto& [edge, coefficient] :
                 relation.physicalArcCoefficients)
                if (coefficient != 0) physicalAmbient.insert(edge);
        std::vector<std::map<ZhangGraphEdge, int>> networkCycles;
        for (const auto& chord : target.currentChords)
        {
            const auto cycle = zhangFundamentalCycle(currentBasis, chord);
            if (cycle.empty())
            {
                result.failureReason =
                    "NETWORK_INTEGER_BASIS_CYCLE_EXPANSION_FAILED";
                return result;
            }
            for (const auto& [edge, coefficient] : cycle)
                if (coefficient != 0) physicalAmbient.insert(edge);
            networkCycles.push_back(cycle);
        }
        result.physicalArcColumns.assign(
            physicalAmbient.begin(), physicalAmbient.end());
        std::map<ZhangGraphEdge, std::size_t> physicalColumnIndex;
        for (std::size_t column = 0;
             column < result.physicalArcColumns.size(); column++)
        {
            physicalColumnIndex[result.physicalArcColumns[column]] = column;
        }
        ZhangExactMatrix physicalRows;
        for (const auto& relation : result.namedRelations)
        {
            ZhangExactVector row(result.physicalArcColumns.size());
            for (const auto& [edge, coefficient] :
                 relation.physicalArcCoefficients)
            {
                row[physicalColumnIndex.at(edge)] = coefficient;
            }
            physicalRows.push_back(std::move(row));
        }

        // D_N^T: the already legal network integer-estimable basis, expressed
        // as exact fundamental-cycle rows in the same physical-arc ambient
        // coordinate as H_P.
        for (const auto& cycle : networkCycles)
        {
            ZhangExactVector row(result.physicalArcColumns.size());
            for (const auto& [edge, coefficient] : cycle)
            {
                auto column = physicalColumnIndex.find(edge);
                if (column == physicalColumnIndex.end())
                {
                    // This is now an internal construction invariant rather
                    // than an expected runtime representation failure.
                    result.failureReason =
                        "NETWORK_INTEGER_BASIS_ARC_OUTSIDE_AMBIENT_SPACE";
                    return result;
                }
                row[column->second] = coefficient;
            }
            result.networkIntegerBasis.push_back(std::move(row));
        }
        if (!result.nuisanceOrthogonal)
        {
            result.failureReason =
                "PRODUCT_RELATION_REAL_NUISANCE_NOT_ANNIHILATED";
            return result;
        }

        ZhangExactRowHnf hnf = zhangExactRowHermiteNormalForm(physicalRows);
        if (!hnf.consistent)
        {
            result.failureReason = "PRODUCT_RELATION_HNF_FAILED";
            return result;
        }
        result.exactRowBasis = hnf.basis;
        result.exactRank = hnf.basis.size();
        result.exactHnf = zhangExactMatrixFingerprint(result.exactRowBasis);

        // Machine closure of H_P = U D_N^T.  Membership returns each exact
        // integer coefficient row of U; equality below is cpp_int exact zero,
        // never a floating tolerance check.
        result.networkLatticeContained = true;
        const auto networkMemberships = zhangIntegerRowLatticeContainsBatch(
            result.networkIntegerBasis, result.exactRowBasis);
        for (const auto& membership : networkMemberships)
        {
            if (!membership.contained ||
                membership.combination.size() !=
                    result.networkIntegerBasis.size())
            {
                result.networkLatticeContained = false;
                break;
            }
            result.networkContainmentTransform.push_back(
                membership.combination);
        }
        result.networkClosureExactZero = result.networkLatticeContained &&
            zhangExactMultiply(
                result.networkContainmentTransform,
                result.networkIntegerBasis) == result.exactRowBasis;
        if (!result.networkClosureExactZero)
        {
            result.failureReason =
                "PRODUCT_RELATION_NOT_IN_NETWORK_INTEGER_LATTICE";
            return result;
        }
        ZhangExactMatrix independentRows;
        int independentRank = 0;
        for (int row = 0; row < static_cast<int>(physicalRows.size()); row++)
        {
            ZhangExactMatrix candidate = independentRows;
            candidate.push_back(physicalRows[row]);
            const int candidateRank = zhangExactRowHermiteNormalForm(
                candidate).basis.size();
            if (candidateRank > independentRank)
            {
                result.independentNamedIndices.push_back(row);
                independentRows.push_back(physicalRows[row]);
                independentRank = candidateRank;
            }
        }
        if (independentRank != result.exactRank)
        {
            result.failureReason =
                "PRODUCT_RELATION_INDEPENDENT_NAMED_RANK_MISMATCH";
            return result;
        }
        const ZhangIntegerLatticeMembership smith =
            zhangIntegerRowLatticeContains(
                physicalRows,
                ZhangExactVector(result.physicalArcColumns.size()));
        if (smith.rank != result.exactRank)
        {
            result.failureReason = "PRODUCT_RELATION_HNF_SNF_RANK_MISMATCH";
            return result;
        }
        result.saturationIndex = 1;
        result.primitive = true;
        for (const auto& invariant : smith.smithInvariants)
        {
            const ZhangExactInteger magnitude = zhangExactAbs(invariant);
            result.saturationIndex *= magnitude;
            result.primitive &= magnitude == 1;
        }
        if (!result.primitive)
        {
            result.failureReason =
                "PRODUCT_RELATION_EXACT_LATTICE_NOT_PRIMITIVE";
            return result;
        }
        result.primitiveRank = result.primitive ? result.exactRank : 0;
        // A saturated (index-one) primitive sublattice of Z^E is a direct
        // summand, hence it admits an integer unimodular completion.
        result.admissibleCompletionProven = result.primitive &&
            result.saturationIndex == 1 &&
            result.primitiveRank == result.exactRank;
        result.mappableNamedIndices = result.independentNamedIndices;
        result.mappableTargetRank = result.exactRank;
        result.unmappableTargetRank =
            result.fullTargetRank - result.mappableTargetRank;
        result.temporalRecoveryRequired = result.unmappableTargetRank > 0;
        result.affineOffsets = ZhangExactVector(result.exactRank);
        result.valid = result.exactRank > 0 &&
            result.physicalExpansionValid && result.nuisanceOrthogonal &&
            result.networkClosureExactZero &&
            result.admissibleCompletionProven;
        if (!result.valid && result.failureReason.empty())
        {
            result.failureReason = "EMPTY_PRODUCT_RELATION_EXACT_LATTICE";
        }
        return result;
    }
};
