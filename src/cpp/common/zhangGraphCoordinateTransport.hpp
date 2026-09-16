#pragma once

#include "common/algebra.hpp"
#include "common/zhangFullRank.hpp"

using ZhangCoordinateRow = std::map<KFKey, double>;

inline KFKey zhangTransportReceiverKey(E_Sys system, E_ObsCode code, const std::string& receiver)
{
    KFKey key;
    key.type = KF::PHASE_BIAS; key.Sat = SatSys(system, 0);
    key.str = receiver; key.num = static_cast<int>(code);
    return key;
}
inline KFKey zhangTransportSatelliteKey(E_ObsCode code, const SatSys& satellite)
{
    KFKey key; key.type = KF::PHASE_BIAS; key.Sat = satellite;
    key.num = static_cast<int>(code); return key;
}
inline KFKey zhangTransportAmbiguityKey(E_ObsCode code, const ZhangGraphEdge& edge)
{
    KFKey key; key.type = KF::AMBIGUITY; key.Sat = edge.satellite;
    key.str = edge.receiver; key.num = static_cast<int>(code); return key;
}
inline void zhangAddCoordinateRow(
    ZhangCoordinateRow& destination, const ZhangCoordinateRow& source, double scale = 1)
{
    for (const auto& [key, value] : source)
    {
        auto& coefficient = destination[key];
        coefficient += scale * value;
        if (coefficient == 0) destination.erase(key);
    }
}

struct ZhangGraphCoordinatePlan
{
    bool valid = false;
    std::string failureReason;
    std::map<KFKey, ZhangCoordinateRow> transform;
    std::map<KFKey, double> independentSourceVariances;
    std::map<std::pair<E_ObsCode, ZhangGraphEdge>, ZhangCoordinateRow> survivingPhysicalRows;
    std::set<ZhangGraphEdge> freshEdges;
    std::set<std::string> componentRoots;
};

// A rectangular change of coordinates on physical phase edge functionals.
// Surviving rows are copied exactly. New arcs get independent, finite prior
// sources, never values estimated from this epoch's measurements. Disconnected
// components keep an uncertain root (and its old cross-covariances), not a
// second hard zero datum. No integer certification is performed here.
inline ZhangGraphCoordinatePlan zhangPlanGraphCoordinateTransport(
    const KFState& state, E_Sys system,
    const ZhangGraphBasis& oldBasis, const ZhangGraphBasis& destination,
    const std::set<ZhangGraphEdge>& retiredArcs,
    const std::map<E_ObsCode, double>& wavelengths,
    const std::function<double(const ZhangGraphEdge&, E_ObsCode)>& freshArcVariance,
    const std::function<double(const std::string&, E_ObsCode)>& freshGaugeVariance,
    const std::map<ZhangGraphEdge, std::set<E_ObsCode>>* retiredSignals = nullptr)
{
    ZhangGraphCoordinatePlan plan;
    auto target = [&](const KFKey& key)
    {
        return key.Sat.sys == system && wavelengths.contains(static_cast<E_ObsCode>(key.num)) &&
            (key.type == KF::PHASE_BIAS || key.type == KF::AMBIGUITY);
    };
    for (const auto& [key, index] : state.kfIndexMap)
    {
        if (!target(key)) plan.transform[key][key] = 1;
        else
        {
            const auto q = state.procNoiseMap.find(key);
            const auto tau = state.gaussMarkovTauMap.find(key);
            if ((q != state.procNoiseMap.end() && q->second != 0) ||
                (tau != state.gaussMarkovTauMap.end() && tau->second > 0))
            {
                plan.failureReason = "NONCONSTANT_PHASE_PROCESS_REQUIRES_FULL_MODEL_TRANSPORT";
                return plan;
            }
        }
    }
    if (destination.edges.empty() || wavelengths.empty())
    {
        plan.failureReason = "EMPTY_DESTINATION_PHASE_GRAPH";
        return plan;
    }
    auto sourcePresent = [&](const KFKey& key)
    {
        return state.kfIndexMap.contains(key) && state.stateTransitionMap.contains(key);
    };
    for (const auto& [code, wavelength] : wavelengths)
    {
        if (!(wavelength > 0) || !std::isfinite(wavelength))
        {
            plan.failureReason = "INVALID_WAVELENGTH"; return plan;
        }
        std::map<ZhangGraphEdge, ZhangCoordinateRow> physical;
        for (const auto& edge : destination.edges)
        {
            const auto receiver = zhangTransportReceiverKey(system, code, edge.receiver);
            const auto satellite = zhangTransportSatelliteKey(code, edge.satellite);
            const auto ambiguity = zhangTransportAmbiguityKey(code, edge);
            const bool receiverKnown = edge.receiver == oldBasis.rootReceiver || sourcePresent(receiver);
            const bool satelliteKnown = sourcePresent(satellite);
            const bool represented = oldBasis.edges.contains(edge) &&
                (!retiredArcs.contains(edge) || (retiredSignals && retiredSignals->contains(edge) &&
                 !retiredSignals->at(edge).contains(code))) && receiverKnown && satelliteKnown &&
                (oldBasis.treeEdges.contains(edge) || sourcePresent(ambiguity));
            auto& row = physical[edge];
            if (edge.receiver != oldBasis.rootReceiver && receiverKnown) row[receiver] = 1;
            if (satelliteKnown) row[satellite] = 1;
            if (represented)
            {
                if (!oldBasis.treeEdges.contains(edge)) row[ambiguity] = wavelength;
                plan.survivingPhysicalRows[{code, edge}] = row;
            }
            else
            {
                auto independent = ambiguity;
                independent.str = "__ZHANG_NEW_ARC__:" + edge.receiver;
                const double variance = freshArcVariance(edge, code);
                if (!(variance > 0) || !std::isfinite(variance) || state.kfIndexMap.contains(independent))
                {
                    plan.failureReason = "INVALID_FRESH_ARC_PRIOR"; return plan;
                }
                row[independent] = 1; // metres, not cycles
                plan.independentSourceVariances[independent] = variance;
                plan.freshEdges.insert(edge);
            }
        }

        std::map<std::string, ZhangCoordinateRow> receivers;
        std::map<SatSys, ZhangCoordinateRow> satellites;
        std::set<std::string> knownReceivers;
        std::set<SatSys> knownSatellites;
        auto seedRoot = [&](const std::string& root) -> bool
        {
            knownReceivers.insert(root);
            plan.componentRoots.insert(root);
            if (root == destination.rootReceiver) { receivers[root] = {}; return true; }
            const auto key = zhangTransportReceiverKey(system, code, root);
            if (sourcePresent(key)) receivers[root][key] = 1;
            else
            {
                auto independent = key;
                independent.str = "__ZHANG_NEW_GAUGE__:" + root;
                const double variance = freshGaugeVariance(root, code);
                if (!(variance > 0) || !std::isfinite(variance)) return false;
                receivers[root][independent] = 1;
                plan.independentSourceVariances[independent] = variance;
            }
            return true;
        };
        if (destination.receivers.contains(destination.rootReceiver)) seedRoot(destination.rootReceiver);
        while (knownReceivers.size() < destination.receivers.size() ||
               knownSatellites.size() < destination.satellites.size())
        {
            bool progressed = false;
            for (const auto& edge : destination.treeEdges)
            {
                const bool r = knownReceivers.contains(edge.receiver);
                const bool s = knownSatellites.contains(edge.satellite);
                if (r && !s)
                {
                    satellites[edge.satellite] = physical.at(edge);
                    zhangAddCoordinateRow(satellites[edge.satellite], receivers.at(edge.receiver), -1);
                    knownSatellites.insert(edge.satellite); progressed = true;
                }
                else if (s && !r)
                {
                    receivers[edge.receiver] = physical.at(edge);
                    zhangAddCoordinateRow(receivers[edge.receiver], satellites.at(edge.satellite), -1);
                    knownReceivers.insert(edge.receiver); progressed = true;
                }
            }
            if (!progressed)
            {
                auto root = std::find_if(destination.receivers.begin(), destination.receivers.end(),
                    [&](const auto& receiver) { return !knownReceivers.contains(receiver); });
                if (root == destination.receivers.end() || !seedRoot(*root))
                {
                    plan.failureReason = "INVALID_COMPONENT_FOREST_OR_GAUGE_PRIOR"; return plan;
                }
            }
        }
        for (const auto& [receiver, row] : receivers)
            if (receiver != destination.rootReceiver)
                plan.transform[zhangTransportReceiverKey(system, code, receiver)] = row;
        for (const auto& [satellite, row] : satellites)
            plan.transform[zhangTransportSatelliteKey(code, satellite)] = row;
        for (const auto& edge : destination.edges)
        {
            if (destination.treeEdges.contains(edge)) continue;
            auto row = physical.at(edge);
            zhangAddCoordinateRow(row, receivers.at(edge.receiver), -1);
            zhangAddCoordinateRow(row, satellites.at(edge.satellite), -1);
            for (auto& [key, coefficient] : row) coefficient /= wavelength;
            plan.transform[zhangTransportAmbiguityKey(code, edge)] = std::move(row);
        }
        // Prove H_new [T G] == [H_old 0] symbolically for every surviving
        // physical row. This also proves all their covariance/cross-covariance
        // identities, without forming a large H P H' during every graph event.
        for (const auto& [identity, before] : plan.survivingPhysicalRows)
        {
            if (identity.first != code) continue;
            const auto& edge = identity.second;
            ZhangCoordinateRow difference;
            zhangAddCoordinateRow(difference, receivers.at(edge.receiver));
            zhangAddCoordinateRow(difference, satellites.at(edge.satellite));
            if (!destination.treeEdges.contains(edge))
                zhangAddCoordinateRow(difference,
                    plan.transform.at(zhangTransportAmbiguityKey(code, edge)), wavelength);
            zhangAddCoordinateRow(difference, before, -1);
            for (const auto& [key, value] : difference)
                if (std::abs(value) > 1e-10)
                {
                    plan.failureReason = "SURVIVING_PHYSICAL_ROW_IDENTITY_FAILED"; return plan;
                }
        }
    }
    plan.valid = true;
    return plan;
}
