#pragma once

#include "common/algebra.hpp"
#include "common/zhangFullRank.hpp"
#include "common/zhangProductHandoffFloatRegistry.hpp"

#include <iomanip>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>

using ZhangCoordinateRow = std::map<KFKey, double>;

// A broken old arc may still have a valid observation of its new physical
// arc.  Current product support tests observation/QC availability here; arc
// continuity and integer inheritance are separate versioned proof gates.
inline bool zhangCurrentObservedProductEdge(
    const ZhangGraphEdge& edge,
    const std::set<ZhangGraphEdge>& observed,
    const std::set<ZhangGraphEdge>& qcExcluded,
    const std::set<ZhangGraphEdge>& elevationExcluded)
{
    return observed.contains(edge) && !qcExcluded.contains(edge) &&
        !elevationExcluded.contains(edge);
}

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
    std::map<std::pair<E_ObsCode, ZhangGraphEdge>, KFKey>
        retainedOldArcAuxiliaryKeys;
};

inline KFKey zhangGraphOldArcAuxiliaryKey(
    const KFState& state,
    E_ObsCode code,
    const ZhangGraphEdge& edge,
    int oldArcVersion)
{
    // Hexadecimal floating-point text preserves every bit of the finite
    // long-double epoch value without an out-of-line GTime formatter.
    std::ostringstream epochId;
    epochId.imbue(std::locale::classic());
    epochId << std::hexfloat
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << state.time.bigTime;
    KFKey key;
    key.type = KF::PRODUCT_HANDOFF;
    key.Sat = edge.satellite;
    key.num = static_cast<int>(code);
    key.str = "ZHANG_OLD_ARC:" + edge.receiver + "@"
        + std::to_string(oldArcVersion) + "#" + epochId.str();
    return key;
}

/** Pull a public linear functional on the post-event coordinate back through
 * x_new = T x_old + G u.  u contains independent zero-mean fresh-arc/gauge
 * priors; a nonzero coefficient on u forbids claiming same-time continuity.
 * Affine offsets are measured in metres, not ambiguity cycles. */
struct ZhangGraphFunctionalPullback
{
    bool valid = false;
    std::string failureReason;
    ZhangCoordinateRow oldStateTerms;
    ZhangCoordinateRow independentPriorTerms;
    double affineOffsetMetres = 0;
};

inline ZhangGraphFunctionalPullback zhangPullbackGraphFunctional(
    const ZhangGraphCoordinatePlan& plan,
    const ZhangCoordinateRow& newRow,
    double newOffsetMetres = 0)
{
    ZhangGraphFunctionalPullback result;
    if (!plan.valid || !std::isfinite(newOffsetMetres))
    {
        result.failureReason = !plan.valid
            ? "INVALID_GRAPH_COORDINATE_PLAN" : "NONFINITE_FUNCTIONAL_OFFSET";
        return result;
    }
    result.affineOffsetMetres = newOffsetMetres;
    for (const auto& [key, coefficient] : newRow)
    {
        if (!std::isfinite(coefficient))
        {
            result.failureReason = "NONFINITE_FUNCTIONAL_COEFFICIENT";
            return result;
        }
        const auto transformed = plan.transform.find(key);
        if (transformed == plan.transform.end())
        {
            result.failureReason = "FUNCTIONAL_DESTINATION_MISSING";
            return result;
        }
        for (const auto& [source, loading] : transformed->second)
        {
            if (!std::isfinite(loading))
            {
                result.failureReason = "NONFINITE_COORDINATE_LOADING";
                return result;
            }
            auto& terms = plan.independentSourceVariances.contains(source)
                ? result.independentPriorTerms : result.oldStateTerms;
            terms[source] += coefficient * loading;
            if (terms[source] == 0) terms.erase(source);
        }
    }
    result.valid = true;
    result.failureReason = "NONE";
    return result;
}

struct ZhangGraphFunctionalProof
{
    bool sameTimeMapProven = false;
    std::string failureReason;
    ZhangGraphFunctionalPullback pullback;
};

/** Prove M_new T=M_old, M_new G=0 and equal affine offsets.  This is a
 * state-functional proof only.  It does not certify a physical integer bridge
 * or authorize a user AR product by itself. */
inline ZhangGraphFunctionalProof zhangVerifyGraphFunctionalTransport(
    const ZhangGraphCoordinatePlan& plan,
    const ZhangCoordinateRow& oldRow,
    const ZhangCoordinateRow& newRow,
    double oldOffsetMetres,
    double newOffsetMetres,
    double tolerance = 1e-10)
{
    ZhangGraphFunctionalProof proof;
    proof.pullback = zhangPullbackGraphFunctional(
        plan, newRow, newOffsetMetres);
    if (!proof.pullback.valid || !std::isfinite(oldOffsetMetres)
        || !(tolerance >= 0) || !std::isfinite(tolerance))
    {
        proof.failureReason = proof.pullback.valid
            ? "INVALID_OLD_OFFSET_OR_TOLERANCE" : proof.pullback.failureReason;
        return proof;
    }
    if (std::abs(newOffsetMetres - oldOffsetMetres) > tolerance)
    {
        proof.failureReason = "AFFINE_OFFSET_MISMATCH";
        return proof;
    }
    for (const auto& [source, coefficient] :
         proof.pullback.independentPriorTerms)
    {
        if (std::abs(coefficient) > tolerance)
        {
            proof.failureReason = "NEW_PRIOR_ENTERED_PUBLIC_FUNCTIONAL";
            return proof;
        }
    }
    ZhangCoordinateRow residual = proof.pullback.oldStateTerms;
    zhangAddCoordinateRow(residual, oldRow, -1);
    for (const auto& [source, coefficient] : residual)
    {
        if (!std::isfinite(coefficient) || std::abs(coefficient) > tolerance)
        {
            proof.failureReason = "PRODUCT_FUNCTIONAL_PULLBACK_MISMATCH";
            return proof;
        }
    }
    proof.sameTimeMapProven = true;
    proof.failureReason = "NONE";
    return proof;
}

struct ZhangGraphFunctionalMoments
{
    bool valid = false;
    std::string failureReason;
    VectorXd mean;
    MatrixXd covariance;
};

/** Project an arbitrary set of post-event product rows directly from the
 * frozen pre-event FLOAT posterior, retaining every cross-covariance. */
inline ZhangGraphFunctionalMoments zhangGraphFunctionalMoments(
    const KFState& oldState,
    const ZhangGraphCoordinatePlan& plan,
    const std::vector<std::pair<ZhangCoordinateRow, double>>& newTargets)
{
    ZhangGraphFunctionalMoments result;
    const int count = static_cast<int>(newTargets.size());
    result.mean = VectorXd::Zero(count);
    result.covariance = MatrixXd::Zero(count, count);
    std::vector<ZhangGraphFunctionalPullback> rows;
    rows.reserve(count);
    for (int i = 0; i < count; i++)
    {
        auto row = zhangPullbackGraphFunctional(
            plan, newTargets[i].first, newTargets[i].second);
        if (!row.valid)
        {
            result.failureReason = row.failureReason;
            return result;
        }
        result.mean(i) = row.affineOffsetMetres;
        for (const auto& [key, coefficient] : row.oldStateTerms)
        {
            auto stateIndex = oldState.kfIndexMap.find(key);
            if (stateIndex == oldState.kfIndexMap.end())
            {
                result.failureReason = "FUNCTIONAL_OLD_STATE_MISSING";
                return result;
            }
            result.mean(i) += coefficient * oldState.x(stateIndex->second);
        }
        rows.push_back(std::move(row));
    }
    for (int i = 0; i < count; i++)
    for (int j = i; j < count; j++)
    {
        double covariance = 0;
        for (const auto& [firstKey, firstCoefficient] : rows[i].oldStateTerms)
        for (const auto& [secondKey, secondCoefficient] : rows[j].oldStateTerms)
        {
            const int first = oldState.kfIndexMap.at(firstKey);
            const int second = oldState.kfIndexMap.at(secondKey);
            covariance += firstCoefficient * secondCoefficient
                * oldState.P(first, second);
        }
        for (const auto& [source, firstCoefficient] :
             rows[i].independentPriorTerms)
        {
            const auto second = rows[j].independentPriorTerms.find(source);
            if (second != rows[j].independentPriorTerms.end())
                covariance += firstCoefficient * second->second
                    * plan.independentSourceVariances.at(source);
        }
        result.covariance(i, j) = result.covariance(j, i) = covariance;
    }
    if (!result.mean.allFinite() || !result.covariance.allFinite())
    {
        result.failureReason = "FUNCTIONAL_MOMENT_NONFINITE";
        return result;
    }
    result.valid = true;
    result.failureReason = "NONE";
    return result;
}

struct ZhangGraphJointMomentAgreement
{
    bool proven = false;
    std::string failureReason = "NOT_EVALUATED";
    double maximumMeanError = 0;
    double maximumCovarianceError = 0;
};

/** Compare every entry of one product block's mean and full covariance,
 * including off-diagonal satellite, frequency and clock cross terms. */
inline ZhangGraphJointMomentAgreement zhangVerifyGraphJointMoments(
    const VectorXd& expectedMean,
    const MatrixXd& expectedCovariance,
    const VectorXd& actualMean,
    const MatrixXd& actualCovariance,
    double relativeTolerance = 1e-8)
{
    ZhangGraphJointMomentAgreement result;
    const int size = expectedMean.size();
    if (size == 0 || actualMean.size() != size
        || expectedCovariance.rows() != size
        || expectedCovariance.cols() != size
        || actualCovariance.rows() != size
        || actualCovariance.cols() != size
        || !(std::isfinite(relativeTolerance)
            && relativeTolerance > 0))
    {
        result.failureReason = "JOINT_PRODUCT_MOMENT_DIMENSION_INVALID";
        return result;
    }
    if (!expectedMean.allFinite() || !actualMean.allFinite()
        || !expectedCovariance.allFinite()
        || !actualCovariance.allFinite())
    {
        result.failureReason = "JOINT_PRODUCT_MOMENT_NONFINITE";
        return result;
    }
    for (int i = 0; i < size; i++)
    {
        const double before = expectedMean(i);
        const double after = actualMean(i);
        const double error = std::abs(after - before);
        result.maximumMeanError = std::max(result.maximumMeanError, error);
        if (error > relativeTolerance *
                std::max({1.0, std::abs(before), std::abs(after)}))
        {
            result.failureReason = "JOINT_PRODUCT_MEAN_TRANSPORT_MISMATCH";
            return result;
        }
    }
    for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
    {
        const double before = expectedCovariance(i, j);
        const double after = actualCovariance(i, j);
        const double error = std::abs(after - before);
        result.maximumCovarianceError = std::max(
            result.maximumCovarianceError, error);
        if (error > relativeTolerance *
                std::max({1.0, std::abs(before), std::abs(after)}))
        {
            result.failureReason =
                "JOINT_PRODUCT_COVARIANCE_TRANSPORT_MISMATCH";
            return result;
        }
    }
    result.proven = true;
    result.failureReason = "COMPLETE_PHASE_CLOCK_MEAN_COVARIANCE_MATCH";
    return result;
}

// A rectangular change of coordinates on physical phase edge functionals.
// Surviving rows are copied exactly. New arcs get independent, finite prior
// sources, never values estimated from this epoch's measurements. Disconnected
// components keep an uncertain root (and its old cross-covariances), not a
// second hard zero datum. No integer certification is performed here.
// A signal admitted on only one frequency may add a cycle coordinate, but it
// cannot supply a shared receiver/satellite datum or join two common-tree
// components. In particular, the other frequency must not acquire a fictive
// physical arc merely because the edge appears in this union basis.
inline ZhangGraphBasis zhangAugmentCommonTreeWithSignalChords(
    ZhangGraphBasis common,
    const std::set<ZhangGraphEdge>& signalCandidates)
{
    for (const auto& edge : signalCandidates)
    {
        if (common.edges.contains(edge) ||
            !common.receivers.contains(edge.receiver) ||
            !common.satellites.contains(edge.satellite) ||
            zhangFundamentalCycle(common, edge).empty())
            continue;
        common.edges.insert(edge);
    }
    return common;
}

// The raw observation inventory is not the committed physical state graph.
// A phase row can remain in the KF measurement after its edge has left the
// active graph component; it must not enter checkpointed graph admissions or
// current integer support. An empty old-history catalogue is allowed here:
// signal arc versions are assigned separately from this support selection.
inline std::map<int, std::set<ZhangGraphEdge>>
zhangGraphSupportedSignalEdges(
    const std::map<int, std::set<ZhangGraphEdge>>& signalEdges,
    const std::set<ZhangGraphEdge>& stateEdges,
    const std::set<ZhangGraphEdge>& basisEdges)
{
    std::map<int, std::set<ZhangGraphEdge>> supported;
    for (const auto& [code, edges] : signalEdges)
        for (const auto& edge : edges)
            if (stateEdges.contains(edge) && basisEdges.contains(edge))
                supported[code].insert(edge);
    return supported;
}

inline std::map<int, std::set<ZhangGraphEdge>>
zhangFinalAdmittedSignalFactorRows(
    E_Sys system,
    const KFMeas& finalMeasurement,
    const std::map<int, std::set<ZhangGraphEdge>>& preQcSignalEdges,
    const std::map<ZhangGraphEdge, std::map<int, int>>& preQcAdmissions)
{
    std::map<int, std::set<ZhangGraphEdge>> retained;
    if (finalMeasurement.obsKeys.size() !=
            size_t(finalMeasurement.H.rows()) ||
        finalMeasurement.metaDataMaps.size() !=
            finalMeasurement.obsKeys.size() ||
        finalMeasurement.R.rows() != finalMeasurement.H.rows() ||
        finalMeasurement.R.cols() != finalMeasurement.H.rows())
        return retained;
    for (size_t index = 0;
         index < finalMeasurement.obsKeys.size(); index++)
    {
        const auto& key = finalMeasurement.obsKeys[index];
        if (key.type != KF::PHAS_MEAS ||
            key.Sat.sys != system || key.Sat.prn <= 0 ||
            key.str.empty()) continue;
        const auto& meta = finalMeasurement.metaDataMaps[index];
        const auto deweighted = meta.find("zhangPppDeweighted");
        const auto phaseRejectCount = meta.find("phaseRejectCount");
        if ((deweighted != meta.end() &&
             deweighted->second != nullptr) ||
            phaseRejectCount == meta.end() ||
            phaseRejectCount->second == nullptr ||
            !(finalMeasurement.R(index, index) > 0) ||
            !std::isfinite(finalMeasurement.R(index, index)))
            continue;
        const ZhangGraphEdge edge{key.str, key.Sat};
        const auto signal = preQcSignalEdges.find(key.num);
        const auto admission = preQcAdmissions.find(edge);
        if (signal == preQcSignalEdges.end() ||
            !signal->second.contains(edge) ||
            admission == preQcAdmissions.end() ||
            !admission->second.contains(key.num)) continue;
        retained[key.num].insert(edge);
    }
    return retained;
}

inline ZhangGraphCoordinatePlan zhangPlanGraphCoordinateTransport(
    const KFState& state, E_Sys system,
    const ZhangGraphBasis& oldBasis, const ZhangGraphBasis& destination,
    const std::set<ZhangGraphEdge>& retiredArcs,
    const std::map<E_ObsCode, double>& wavelengths,
    const std::function<double(const ZhangGraphEdge&, E_ObsCode)>& freshArcVariance,
    const std::function<double(const std::string&, E_ObsCode)>& freshGaugeVariance,
    const std::map<ZhangGraphEdge, std::set<E_ObsCode>>* retiredSignals = nullptr,
    bool retainRetiredOldTreeArcs = false,
    const std::map<ZhangGraphEdge, int>* oldArcVersions = nullptr,
    const std::set<ZhangGraphEdge>* oldProductTargetTreeEdges = nullptr,
    const std::map<ZhangGraphEdge, std::map<int, int>>*
        oldSignalArcVersions = nullptr,
    const std::map<int, std::set<ZhangGraphEdge>>*
        destinationSignalEdges = nullptr)
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
        if (destinationSignalEdges)
            for (const auto& edge : destination.treeEdges)
                if (!destinationSignalEdges->contains(int(code)) ||
                    !destinationSignalEdges->at(int(code)).contains(edge))
                {
                    plan.failureReason =
                        "COMMON_TREE_SIGNAL_SUPPORT_MISSING";
                    return plan;
                }
        if (retainRetiredOldTreeArcs)
        {
            if (!oldArcVersions && !oldSignalArcVersions)
            {
                plan.failureReason = "OLD_ARC_VERSION_MAP_MISSING";
                return plan;
            }
            std::set<ZhangGraphEdge> auxiliaryCandidates =
                oldBasis.treeEdges;
            if (oldProductTargetTreeEdges)
                auxiliaryCandidates.insert(
                    oldProductTargetTreeEdges->begin(),
                    oldProductTargetTreeEdges->end());
            for (const auto& edge : auxiliaryCandidates)
            {
                const bool retiredForSignal = retiredArcs.contains(edge)
                    && (!retiredSignals || !retiredSignals->contains(edge)
                        || retiredSignals->at(edge).contains(code));
                if (!retiredForSignal && destination.edges.contains(edge))
                    continue;
                std::optional<int> oldVersion;
                if (oldSignalArcVersions)
                {
                    const auto identity = oldSignalArcVersions->find(edge);
                    if (identity != oldSignalArcVersions->end())
                    {
                        const auto signal = identity->second.find(int(code));
                        if (signal != identity->second.end())
                            oldVersion = signal->second;
                    }
                }
                else if (oldArcVersions)
                {
                    const auto topology = oldArcVersions->find(edge);
                    if (topology != oldArcVersions->end())
                        oldVersion = topology->second;
                }
                const auto receiver = zhangTransportReceiverKey(
                    system, code, edge.receiver);
                const auto satellite = zhangTransportSatelliteKey(
                    code, edge.satellite);
                const auto ambiguity = zhangTransportAmbiguityKey(
                    code, edge);
                if (!oldBasis.edges.contains(edge)
                    || !oldVersion || *oldVersion < 0
                    || (edge.receiver != oldBasis.rootReceiver
                        && !sourcePresent(receiver))
                    || !sourcePresent(satellite)
                    || (!oldBasis.treeEdges.contains(edge)
                        && !sourcePresent(ambiguity)))
                {
                    // Complete product support can contain signal-specific
                    // chords absent on this frequency. They are not source
                    // FLOAT coordinates: preserve the filter transaction but
                    // leave that optional bridge support explicitly missing.
                    // A source estimation-tree edge must still be complete.
                    if (!oldBasis.treeEdges.contains(edge)) continue;
                    plan.failureReason =
                        "OLD_PRODUCT_OR_TREE_ARC_AUXILIARY_SOURCE_MISSING";
                    return plan;
                }
                const KFKey auxiliary = zhangGraphOldArcAuxiliaryKey(
                    state, code, edge, *oldVersion);
                KFKey retained;
                if (zhangHandoffFloatFind(state,
                        {edge, int(code), *oldVersion}, retained))
                {
                    plan.retainedOldArcAuxiliaryKeys[{code, edge}] = retained;
                    continue;
                }
                if (state.kfIndexMap.contains(auxiliary)
                    || plan.transform.contains(auxiliary))
                {
                    plan.failureReason = "OLD_TREE_ARC_AUXILIARY_DUPLICATE";
                    return plan;
                }
                ZhangCoordinateRow oldPhysical;
                if (edge.receiver != oldBasis.rootReceiver)
                    oldPhysical[receiver] = 1;
                oldPhysical[satellite] = 1;
                if (!oldBasis.treeEdges.contains(edge))
                    oldPhysical[ambiguity] = wavelength;
                plan.transform[auxiliary] = std::move(oldPhysical);
                plan.retainedOldArcAuxiliaryKeys[{code, edge}] = auxiliary;
            }
        }
        std::map<ZhangGraphEdge, ZhangCoordinateRow> physical;
        for (const auto& edge : destination.edges)
        {
            if (destinationSignalEdges &&
                (!destinationSignalEdges->contains(int(code)) ||
                 !destinationSignalEdges->at(int(code)).contains(edge)))
                continue;
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
            if (!physical.contains(edge)) continue;
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
