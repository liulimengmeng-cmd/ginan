#pragma once

#include "common/zhangProductBasisMode.hpp"
#include "common/zhangZeroDatumNativeProduct.hpp"

/** The request is frozen by the authoritative FLOAT owner, not by a cloned
 * disposable AR branch. Equality of a chord map, product rank zero, or a
 * conditioning receipt alone does not prove which committed S graph defines
 * the native products. The caller supplies the ACTUAL authoritative snapshot.
 */
template<class GraphContext>
bool zhangSingleSBasisNativeViewValid(const GraphContext& graph, std::string& reason)
{
    auto fail = [&](const char* why) { reason = why; return false; };
    if (graph.productBasisMode != "SINGLE_S_BASIS_NATIVE" ||
        graph.committedSBasisIdentity.empty())
        return fail("SINGLE_S_BASIS_NATIVE_GRAPH_MODE_OR_SNAPSHOT_MISSING");
    if (graph.productDatumVersion != graph.committedSGaugeVersion)
        return fail("SINGLE_S_BASIS_NATIVE_PRODUCT_DATUM_NOT_COMMITTED_S_GAUGE");
    const auto& current = graph.basis;
    const auto& view = graph.productBasis;
    if (!graph.initialized || current.rootReceiver != view.rootReceiver ||
        current.edges != view.edges || current.treeEdges != view.treeEdges ||
        current.receivers != view.receivers || current.satellites != view.satellites ||
        current.componentCount != view.componentCount || current.connected != view.connected)
        return fail("SINGLE_S_BASIS_NATIVE_PRODUCT_VIEW_NOT_COMMITTED_S");
    reason = "SINGLE_S_BASIS_NATIVE_PRODUCT_VIEW_IS_COMMITTED_S";
    return true;
}

template<class GraphContext>
bool zhangSingleSBasisNativeSnapshotMatches(
    const ZhangZeroDatumNativeRequest& request,
    const GraphContext& graph,
    std::string& reason)
{
    auto fail = [&](const char* why) { reason = why; return false; };
    if (!request.requested || !request.singleSBasisNative ||
        request.committedSBasisIdentity.empty() ||
        graph.committedSBasisIdentity != request.committedSBasisIdentity)
        return fail("SINGLE_S_BASIS_NATIVE_AUTHORITATIVE_SNAPSHOT_MISMATCH");
    if (!zhangSingleSBasisNativeViewValid(graph, reason)) return false;
    const std::set<SatSys> requested(request.satellites.begin(), request.satellites.end());
    if (requested.size() != request.satellites.size() || requested.size() < 2 ||
        requested != graph.basis.satellites)
        return fail("SINGLE_S_BASIS_NATIVE_AUTHORITATIVE_CATALOGUE_MISMATCH");
    reason = "SINGLE_S_BASIS_NATIVE_UNIQUE_COMMITTED_SNAPSHOT_VERIFIED";
    return true;
}

// A disposable conditioned branch has its own runtime lease. Do not compare
// runtime strings across owners; compare every coordinate/physical fact that
// defines its S chart to the authoritative snapshot instead.
template<class GraphContext>
bool zhangSingleSBasisNativeGraphCoordinatesMatch(
    const GraphContext& source, const GraphContext& branch, std::string& reason)
{
    if (!zhangSingleSBasisNativeViewValid(source, reason) ||
        !zhangSingleSBasisNativeViewValid(branch, reason)) return false;
    const auto& left = source.basis;
    const auto& right = branch.basis;
    if (left.rootReceiver != right.rootReceiver || left.edges != right.edges ||
        left.treeEdges != right.treeEdges || left.receivers != right.receivers ||
        left.satellites != right.satellites || left.componentCount != right.componentCount ||
        left.connected != right.connected ||
        source.signalArcVersions != branch.signalArcVersions ||
        source.signalObservationEdges != branch.signalObservationEdges ||
        source.signalObservationFinalized != branch.signalObservationFinalized ||
        source.representationVersion != branch.representationVersion ||
        source.floatGaugeVersion != branch.floatGaugeVersion ||
        source.integerComponentVersion != branch.integerComponentVersion ||
        source.committedSGaugeVersion != branch.committedSGaugeVersion ||
        source.productDatumVersion != branch.productDatumVersion)
    {
        reason = "SINGLE_S_BASIS_NATIVE_CONDITIONED_BRANCH_S_CHART_MISMATCH";
        return false;
    }
    reason = "SINGLE_S_BASIS_NATIVE_CONDITIONED_BRANCH_S_CHART_VERIFIED";
    return true;
}

/** Current support is per satellite AND signal. Retained arc-version metadata
 * is not an observation. In the native mode the physical datum functional
 * must be the actual committed S-tree path and every nonzero path term must
 * have been admitted by the current final-factor selector for this signal.
 * This does not invent a local integer-rank or covariance-gain threshold.
 * The actual complete posterior receipt and numerical publication gates are
 * separate requirements and are not issued or relaxed here.
 */
template<class GraphContext, class Functional>
bool zhangSingleSBasisNativeSignalSupport(
    const GraphContext& graph, SatSys satellite, E_ObsCode observable,
    const Functional& functional, std::string& reason)
{
    auto fail = [&](const char* why) { reason = why; return false; };
    if (!graph.signalObservationFinalized)
        return fail("SINGLE_S_BASIS_NATIVE_CURRENT_SIGNAL_SUPPORT_NOT_FINALIZED");
    const int signal = static_cast<int>(observable);
    const auto accepted = graph.signalObservationEdges.find(signal);
    if (accepted == graph.signalObservationEdges.end())
        return fail("SINGLE_S_BASIS_NATIVE_CURRENT_SIGNAL_SUPPORT_MISSING");
    bool satelliteObserved = false;
    for (const auto& edge : accepted->second)
    {
        if (edge.satellite != satellite || !graph.basis.edges.contains(edge)) continue;
        const auto versions = graph.signalArcVersions.find(edge);
        if (versions != graph.signalArcVersions.end())
        {
            const auto version = versions->second.find(signal);
            if (version != versions->second.end() && version->second >= 0)
                satelliteObserved = true;
        }
    }
    if (!satelliteObserved)
        return fail("SINGLE_S_BASIS_NATIVE_SATELLITE_SIGNAL_NOT_CURRENTLY_ADMITTED");
    if (!functional.valid || functional.satellite != satellite ||
        functional.observable != observable ||
        functional.affineOffsetCycles != 0 ||
        functional.physicalEdges.size() != functional.networkCoefficients.size() ||
        functional.physicalEdges.size() != functional.physicalArcVersions.size())
        return fail("SINGLE_S_BASIS_NATIVE_PHYSICAL_FUNCTIONAL_INVALID");
    for (std::size_t term = 0; term < functional.networkCoefficients.size(); ++term)
    {
        if (functional.networkCoefficients[term] == 0) continue;
        const auto& edge = functional.physicalEdges[term];
        if (!graph.basis.treeEdges.contains(edge))
            return fail("SINGLE_S_BASIS_NATIVE_PATH_OUTSIDE_COMMITTED_S_TREE");
        if (!accepted->second.contains(edge))
            return fail("SINGLE_S_BASIS_NATIVE_DATUM_PATH_SIGNAL_NOT_CURRENTLY_ADMITTED");
        const auto versions = graph.signalArcVersions.find(edge);
        if (versions == graph.signalArcVersions.end())
            return fail("SINGLE_S_BASIS_NATIVE_DATUM_PATH_PHYSICAL_ARC_MISSING");
        const auto version = versions->second.find(signal);
        if (version == versions->second.end() || version->second < 0 ||
            version->second != functional.physicalArcVersions[term])
            return fail("SINGLE_S_BASIS_NATIVE_DATUM_PATH_PHYSICAL_ARC_VERSION_MISMATCH");
    }
    reason = "SINGLE_S_BASIS_NATIVE_CURRENT_SIGNAL_AND_COMPLETE_S_PATH_VERIFIED";
    return true;
}
