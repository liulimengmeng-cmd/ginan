#pragma once

#include <map>
#include <set>
#include <string>

#include "common/zhangGraphCoordinateTransport.hpp"

struct ZhangProductChartIdentityPrecondition
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
};

/** An unchanged KF S chart is an identity map only for the same physical
 * arc versions.  In particular, a missing old P-tree version must not be
 * replaced with the current edgeHistory version.  Physical retirement needs
 * a separate correlated auxiliary transform, never this identity shortcut. */
inline ZhangProductChartIdentityPrecondition
zhangCheckProductChartIdentityReceiptPreconditions(
    const ZhangGraphBasis& stateBasis,
    const ZhangGraphBasis& candidate,
    const ZhangGraphBasis& oldProduct,
    const ZhangGraphBasis& newProduct,
    const std::map<ZhangGraphEdge, int>& oldProductTreeVersions,
    const std::map<ZhangGraphEdge, int>& currentVersions,
    const std::set<ZhangGraphEdge>& representedEdges,
    const std::set<ZhangGraphEdge>& discontinuousEdges)
{
    ZhangProductChartIdentityPrecondition result;
    if (!discontinuousEdges.empty())
        result.reason = "PHYSICAL_ARC_BREAK_REQUIRES_CORRELATED_TRANSFORM";
    else if (!stateBasis.connected || !candidate.connected ||
        !oldProduct.connected || !newProduct.connected ||
        stateBasis.rootReceiver != candidate.rootReceiver ||
        stateBasis.rootReceiver != oldProduct.rootReceiver ||
        stateBasis.rootReceiver != newProduct.rootReceiver)
        result.reason = "PRODUCT_REFERENCE_ROOT_CHANGED";
    else if (stateBasis.treeEdges != candidate.treeEdges ||
        stateBasis.receivers != candidate.receivers ||
        stateBasis.satellites != candidate.satellites)
        result.reason = "KF_S_CHART_CHANGED";
    else if (oldProduct.satellites.empty() ||
        newProduct.satellites.empty() ||
        *oldProduct.satellites.begin() !=
            *newProduct.satellites.begin())
        result.reason = "PRODUCT_REFERENCE_SATELLITE_CHANGED";
    else
    {
        std::set<std::string> receivers;
        std::set<SatSys> satellites;
        for (const auto& edge : representedEdges)
        {
            receivers.insert(edge.receiver);
            satellites.insert(edge.satellite);
        }
        if (receivers != stateBasis.receivers ||
            satellites != stateBasis.satellites)
            result.reason = "S_NODE_UNIVERSE_CHANGED";
        else
        {
            auto currentPhysical = [&](const ZhangGraphEdge& edge)
            {
                return representedEdges.contains(edge) &&
                    currentVersions.contains(edge);
            };
            for (const auto& edge : stateBasis.treeEdges)
                if (!currentPhysical(edge))
                {
                    result.reason = "OLD_S_TREE_ARC_NOT_CURRENT";
                    return result;
                }
            for (const auto& edge : oldProduct.treeEdges)
            {
                const auto old = oldProductTreeVersions.find(edge);
                const auto current = currentVersions.find(edge);
                if (old == oldProductTreeVersions.end())
                {
                    result.reason = "OLD_PRODUCT_TREE_ARC_VERSION_MISSING";
                    return result;
                }
                if (!currentPhysical(edge) || current == currentVersions.end()
                    || old->second != current->second)
                {
                    result.reason = "OLD_PRODUCT_TREE_ARC_VERSION_CHANGED";
                    return result;
                }
            }
            for (const auto& edge : newProduct.treeEdges)
                if (!currentPhysical(edge))
                {
                    result.reason = "NEW_PRODUCT_TREE_ARC_NOT_CURRENT";
                    return result;
                }
            result.valid = true;
            result.reason = "S_CHART_AND_PHYSICAL_ARC_VERSIONS_IDENTICAL";
        }
    }
    return result;
}
