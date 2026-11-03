#pragma once

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include "common/zhangFullRank.hpp"

/** A graph event is not a physical slip.  Several flags can be true together
 * (e.g. a signal slip can also split a component); the primary name is only
 * a diagnostic classification and never authorizes integer continuity. */
struct ZhangCommittedSGraphEvent
{
    std::string kind = "TOPOLOGY_UNCHANGED";
    bool initialized = false;
    bool physicalArcBreak = false;
    bool chordCatalogueChanged = false;
    bool nodeAugmentation = false;
    bool treeReplacement = false;
    bool componentSplit = false;
    bool componentReconnect = false;
    std::set<ZhangGraphEdge> addedEdges;
    std::set<ZhangGraphEdge> removedEdges;
};

inline ZhangCommittedSGraphEvent zhangClassifyCommittedSGraphEvent(
    const ZhangGraphBasis& before, const ZhangGraphBasis& after,
    bool wasInitialized, bool physicalArcBreak)
{
    ZhangCommittedSGraphEvent out;
    out.initialized = !wasInitialized;
    out.physicalArcBreak = physicalArcBreak;
    std::set_difference(after.edges.begin(), after.edges.end(),
        before.edges.begin(), before.edges.end(),
        std::inserter(out.addedEdges, out.addedEdges.end()));
    std::set_difference(before.edges.begin(), before.edges.end(),
        after.edges.begin(), after.edges.end(),
        std::inserter(out.removedEdges, out.removedEdges.end()));
    out.nodeAugmentation = wasInitialized &&
        (std::includes(after.receivers.begin(), after.receivers.end(),
             before.receivers.begin(), before.receivers.end()) &&
         std::includes(after.satellites.begin(), after.satellites.end(),
             before.satellites.begin(), before.satellites.end())) &&
        (after.receivers != before.receivers || after.satellites != before.satellites);
    out.treeReplacement = wasInitialized &&
        (before.rootReceiver != after.rootReceiver ||
         before.treeEdges != after.treeEdges);
    out.componentSplit = wasInitialized && after.componentCount > before.componentCount;
    out.componentReconnect = wasInitialized && after.componentCount < before.componentCount;
    out.chordCatalogueChanged = wasInitialized &&
        before.rootReceiver == after.rootReceiver && before.treeEdges == after.treeEdges &&
        before.edges != after.edges;
    if (out.initialized) out.kind = "INITIALISE";
    else if (out.physicalArcBreak) out.kind = "PHYSICAL_ARC_BREAK";
    else if (out.componentSplit) out.kind = "COMPONENT_SPLIT";
    else if (out.componentReconnect) out.kind = "COMPONENT_RECONNECT";
    else if (out.nodeAugmentation) out.kind = "NODE_AUGMENTATION";
    else if (out.treeReplacement) out.kind = "TREE_REPLACEMENT";
    else if (out.chordCatalogueChanged)
        out.kind = out.removedEdges.empty() ? "CHORD_ONLY_INSERTION" : "CHORD_CATALOGUE_CHANGE";
    return out;
}

/** Commit extra chord coordinates without choosing a tree or touching x/P.
 * Existing KF chords may be retained while not currently observed.  They
 * remain in the coordinate catalogue, NOT in current signal admission.
 * Unknown endpoints are rejected: node augmentation is a separate event. */
inline bool zhangAugmentCommittedSChords(
    const ZhangGraphBasis& committed,
    const ZhangGraphBasis& currentCandidate,
    const std::set<ZhangGraphEdge>& representedEdges,
    ZhangGraphBasis& output)
{
    if (committed.rootReceiver != currentCandidate.rootReceiver ||
        committed.treeEdges != currentCandidate.treeEdges ||
        committed.receivers != currentCandidate.receivers ||
        committed.satellites != currentCandidate.satellites)
        return false;
    output = currentCandidate;
    for (const auto& edge : representedEdges)
    {
        if (!output.receivers.contains(edge.receiver) ||
            !output.satellites.contains(edge.satellite)) return false;
        output.edges.insert(edge);
    }
    return true;
}
