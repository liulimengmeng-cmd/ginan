#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>
#include "common/zhangFullRank.hpp"

struct ZhangPreselectedProductTargetAudit
{
    bool valid = false;
    bool hasEdge = false;
    ZhangGraphEdge edge;
    std::string reason = "NOT_EVALUATED";
};

/** A pre-factor product selector is not evidence that its new nodes/arcs
 * entered the committed FLOAT.  This read-only audit checks only that
 * selector and its own metadata; unrelated graph failures must not be
 * "repaired" by discarding product state.  Version/history tests agree with
 * the graph checkpoint exporter, and never synthesize an observation. */
inline ZhangPreselectedProductTargetAudit zhangAuditPreselectedProductTarget(
    const ZhangGraphBasis& committed,
    const ZhangGraphBasis& product,
    const std::set<ZhangGraphEdge>& representedTreeEdges,
    const std::map<ZhangGraphEdge, int>& selectedArcs,
    const std::map<ZhangGraphEdge, int>& currentArcs,
    const std::set<std::pair<ZhangGraphEdge, int>>& knownArcHistories,
    const std::map<ZhangGraphEdge, std::map<int, int>>& selectedSignals,
    const std::map<ZhangGraphEdge, std::map<int, int>>& currentSignals,
    const std::vector<int>& requiredSignals)
{
    ZhangPreselectedProductTargetAudit result;
    auto reject = [&](const char* reason)
    {
        result.reason = reason;
        return result;
    };
    if (requiredSignals.empty()) return reject("PRODUCT_BASELINE_SIGNALS_NOT_CONFIGURED");
    if (!product.connected || !product.receivers.contains(product.rootReceiver) ||
        product.rootReceiver != committed.rootReceiver)
        return reject("PRODUCT_ROOT_NOT_IN_COMMITTED_S_CHART");
    for (const auto& edge : product.treeEdges)
    {
        result.edge = edge;
        result.hasEdge = true;
        if (!committed.edges.contains(edge) || !representedTreeEdges.contains(edge))
            return reject("PRODUCT_TREE_EDGE_NOT_REPRESENTED_BY_COMMITTED_S");
        const auto selected = selectedArcs.find(edge);
        const auto current = currentArcs.find(edge);
        if (selected == selectedArcs.end() || selected->second < 0 ||
            current == currentArcs.end() || selected->second > current->second)
            return reject("PRODUCT_TOPOLOGY_ARC_VERSION_INVALID");
        if (!knownArcHistories.contains({edge, selected->second}))
            return reject("PRODUCT_TOPOLOGY_ARC_OBSERVATION_HISTORY_MISSING");
        const auto signals = selectedSignals.find(edge);
        const auto physical = currentSignals.find(edge);
        if (signals == selectedSignals.end() || signals->second.empty() ||
            physical == currentSignals.end())
            return reject("PRODUCT_SIGNAL_ARC_VERSIONS_MISSING");
        for (const auto signal : requiredSignals)
            if (!signals->second.contains(signal) || !physical->second.contains(signal))
                return reject("PRODUCT_BASELINE_SIGNAL_ARC_INCOMPLETE");
        for (const auto& [signal, version] : signals->second)
            if (signal <= 0 || version < 0 || !physical->second.contains(signal) ||
                version > physical->second.at(signal))
                return reject("PRODUCT_SIGNAL_ARC_VERSION_INVALID");
    }
    result.hasEdge = false;
    if (selectedArcs.size() != product.treeEdges.size() ||
        selectedSignals.size() != product.treeEdges.size())
        return reject("PRODUCT_TREE_ARC_METADATA_NOT_BIJECTIVE");
    if (!std::includes(committed.receivers.begin(), committed.receivers.end(),
            product.receivers.begin(), product.receivers.end()) ||
        !std::includes(committed.satellites.begin(), committed.satellites.end(),
            product.satellites.begin(), product.satellites.end()))
        return reject("PRODUCT_NODES_NOT_REPRESENTED_BY_COMMITTED_S");
    result.valid = true;
    result.reason = "CURRENT_PRESELECTED_PRODUCT_TARGET_VALID";
    return result;
}
