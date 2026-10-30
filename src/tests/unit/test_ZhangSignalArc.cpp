#include <boost/test/unit_test.hpp>

#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductPhysicalCycleChart.hpp"
#include "common/zhangProductConditionalMap.hpp"
#include "common/zhangUserIntegerFunctional.hpp"

#include <algorithm>
#include <array>

namespace
{
KFKey signalTestKey(KF type, const char* name)
{
    KFKey key;
    key.type = type;
    key.str = name;
    return key;
}
}

BOOST_AUTO_TEST_CASE(signal_arc_l2_break_does_not_change_l1_product_identity)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g13(E_Sys::GPS, 13);
    const ZhangGraphEdge a01{"A", g01};
    const ZhangGraphEdge b01{"B", g01};
    const ZhangGraphEdge b13{"B", g13};
    const auto basis = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, b13}, "A");
    BOOST_REQUIRE(basis.connected);

    std::map<ZhangGraphEdge, std::map<int, int>> before;
    for (const auto& edge : basis.edges)
    {
        before[edge][int(E_ObsCode::L1C)] = 0;
        before[edge][int(E_ObsCode::L2W)] = 0;
    }
    auto after = before;
    after[b13][int(E_ObsCode::L2W)] = 1;

    const auto l1Before = zhangBuildProductIntegerFunctionals(
        basis, before, E_ObsCode::L1C, g01);
    const auto l1After = zhangBuildProductIntegerFunctionals(
        basis, after, E_ObsCode::L1C, g01);
    const auto l2Before = zhangBuildProductIntegerFunctionals(
        basis, before, E_ObsCode::L2W, g01);
    const auto l2After = zhangBuildProductIntegerFunctionals(
        basis, after, E_ObsCode::L2W, g01);
    BOOST_REQUIRE(l1Before.contains(g13));
    BOOST_REQUIRE(l1After.contains(g13));
    BOOST_REQUIRE(l2Before.contains(g13));
    BOOST_REQUIRE(l2After.contains(g13));
    BOOST_CHECK_EQUAL(
        zhangProductPhysicalFunctionalFingerprint(l1Before.at(g13)),
        zhangProductPhysicalFunctionalFingerprint(l1After.at(g13)));
    BOOST_CHECK_NE(
        zhangProductPhysicalFunctionalFingerprint(l2Before.at(g13)),
        zhangProductPhysicalFunctionalFingerprint(l2After.at(g13)));
    const auto l1Difference = zhangProductIntegerFunctionalDifference(
        l1Before.at(g13), l1After.at(g13));
    const auto l2Difference = zhangProductIntegerFunctionalDifference(
        l2Before.at(g13), l2After.at(g13));
    BOOST_CHECK(l1Difference.valid);
    BOOST_CHECK(l1Difference.coefficients.empty());
    BOOST_CHECK(l2Difference.valid);
    BOOST_CHECK(!l2Difference.coefficients.empty());

    auto unknown = after;
    unknown[b13].erase(int(E_ObsCode::L2W));
    const auto l2Unknown = zhangBuildProductIntegerFunctionals(
        basis, unknown, E_ObsCode::L2W, g01);
    BOOST_CHECK(!l2Unknown.contains(g13));
    BOOST_CHECK(!zhangBuildProductIntegerFunctionals(
        basis, unknown, E_ObsCode::L1C, g01).empty());
	// The flattened compatibility overload must not make a missing signal
	// identity look like version zero either.
	std::map<ZhangGraphEdge, int> partialL2;
	partialL2[a01] = 0;
	partialL2[b01] = 0;
	BOOST_CHECK(!zhangBuildProductIntegerFunctionals(
		basis, partialL2, g01, 0, E_ObsCode::L2W).contains(g13));
}

BOOST_AUTO_TEST_CASE(signal_arc_unknown_unused_chord_keeps_other_products)
{
	const SatSys g01(E_Sys::GPS, 1);
	const SatSys g13(E_Sys::GPS, 13);
	const SatSys g30(E_Sys::GPS, 30);
	const ZhangGraphEdge a01{"A", g01};
	const ZhangGraphEdge b01{"B", g01};
	const ZhangGraphEdge a13{"A", g13};
	const ZhangGraphEdge b13{"B", g13}; // unused product chord
	const ZhangGraphEdge b30{"B", g30};
	const auto basis = zhangBuildSpanningTree(
		std::set<ZhangGraphEdge>{a01, b01, a13, b13, b30}, "A",
		{a01, b01, a13, b30});
	BOOST_REQUIRE(basis.connected);
	BOOST_REQUIRE(!basis.treeEdges.contains(b13));
	std::map<ZhangGraphEdge, std::map<int, int>> versions;
	for (const auto& edge : basis.edges)
		versions[edge][int(E_ObsCode::L2W)] = 4;
	const auto fullyKnown = zhangBuildProductIntegerFunctionals(
		basis, versions, E_ObsCode::L2W, g01);
	versions[b13].erase(int(E_ObsCode::L2W));
	const auto withoutUnusedChord = zhangBuildProductIntegerFunctionals(
		basis, versions, E_ObsCode::L2W, g01);
	BOOST_CHECK(withoutUnusedChord.contains(g13));
	BOOST_CHECK(withoutUnusedChord.contains(g30));
	BOOST_CHECK_EQUAL(
		zhangProductPhysicalFunctionalFingerprint(fullyKnown.at(g13)),
		zhangProductPhysicalFunctionalFingerprint(withoutUnusedChord.at(g13)));
	for (const auto& [satellite, functional] : withoutUnusedChord)
		for (std::size_t i = 0; i < functional.networkCoefficients.size(); ++i)
			if (functional.networkCoefficients[i] != 0)
				BOOST_CHECK_GE(functional.physicalArcVersions[i], 0);
	versions[a13].erase(int(E_ObsCode::L2W));
	const auto missingOnlyG13Path = zhangBuildProductIntegerFunctionals(
		basis, versions, E_ObsCode::L2W, g01);
	BOOST_CHECK(!missingOnlyG13Path.contains(g13));
	BOOST_CHECK(missingOnlyG13Path.contains(g30));
}

BOOST_AUTO_TEST_CASE(signal_arc_absolute_reference_bridge_retains_full_joint_covariance)
{
    const KFKey reference = signalTestKey(KF::PHASE_BIAS, "reference");
    const KFKey commonBridge = signalTestKey(KF::PRODUCT_HANDOFF, "old-gauge");
    const KFKey subject = signalTestKey(KF::PHASE_BIAS, "subject");
    KFState state;
    state.kfIndexMap.clear();
    state.kfIndexMap[reference] = 0;
    state.kfIndexMap[commonBridge] = 1;
    state.kfIndexMap[subject] = 2;
    state.x = VectorXd(3);
    state.x << 1.0, 2.0, 3.0;
    state.P = MatrixXd(3, 3);
    state.P << 4.0, 0.6, -0.3,
               0.6, 9.0, 1.2,
              -0.3, 1.2, 16.0;
    ZhangGraphCoordinatePlan identity;
    identity.valid = true;
    identity.transform[reference][reference] = 1;
    identity.transform[commonBridge][commonBridge] = 1;
    identity.transform[subject][subject] = 1;

    const ZhangCoordinateRow oldReference{{reference, 1}};
    const ZhangCoordinateRow oldSubject{{reference, 1}, {subject, 1}};
    const auto oldMoments = zhangGraphFunctionalMoments(state, identity,
        {{oldReference, 0}, {oldSubject, 0}});
    BOOST_REQUIRE(oldMoments.valid);
    BOOST_CHECK_CLOSE(oldMoments.covariance(0, 1), 3.7, 1e-9);
    const auto valid = zhangVerifyGraphJointMoments(
        oldMoments.mean, oldMoments.covariance,
        oldMoments.mean, oldMoments.covariance);
    BOOST_CHECK(valid.proven);
    auto brokenCovariance = oldMoments.covariance;
    brokenCovariance(0, 1) += 0.5;
    brokenCovariance(1, 0) += 0.5;
    BOOST_CHECK(!zhangVerifyGraphJointMoments(
        oldMoments.mean, oldMoments.covariance,
        oldMoments.mean, brokenCovariance).proven);

    const ZhangGraphEdge commonEdge{"A", SatSys(E_Sys::GPS, 1)};
    ZhangProductIntegerTransition absoluteBridge;
    absoluteBridge.observable = E_ObsCode::L1C;
    absoluteBridge.physicalEdges = {commonEdge};
    absoluteBridge.physicalArcVersions = {1};
    absoluteBridge.coefficients = {1};
    absoluteBridge.valid = true;
    const ZhangCoordinateRow newReference{
        {reference, 1}, {commonBridge, 1}};
    const ZhangCoordinateRow newSubject{
        {reference, 1}, {commonBridge, 1}, {subject, 1}};
    const auto physicalRow = [&](const ZhangGraphEdge& edge, int version,
        ZhangCoordinateRow& row, std::string& reason)
    {
        if (edge != commonEdge || version != 1) return false;
        row = {{commonBridge, 1}};
        reason = "CORRELATED_OLD_REFERENCE";
        return true;
    };
    const auto referenceProof = zhangProveProductConditionalMap(
        oldReference, newReference, absoluteBridge, 1.0, true, physicalRow);
    const auto subjectProof = zhangProveProductConditionalMap(
        oldSubject, newSubject, absoluteBridge, 1.0, true, physicalRow);
    BOOST_CHECK(referenceProof.valid);
    BOOST_CHECK(subjectProof.valid);
    BOOST_CHECK(!referenceProof.bridgeStateRowCycles.empty());
    BOOST_CHECK(!subjectProof.bridgeStateRowCycles.empty());

    const auto offset = zhangVerifyGraphFunctionalTransport(
        identity, oldReference, oldReference, 0.25, 0.0);
    BOOST_CHECK(!offset.sameTimeMapProven);
    BOOST_CHECK_EQUAL(offset.failureReason, "AFFINE_OFFSET_MISMATCH");
    absoluteBridge.affineOffsetCycles = 1;
    BOOST_CHECK(!zhangProveProductConditionalMap(
        oldReference, newReference, absoluteBridge, 1.0, true,
        physicalRow).valid);
}

BOOST_AUTO_TEST_CASE(signal_arc_auxiliary_key_includes_signal_and_version)
{
    KFState state;
    const ZhangGraphEdge edge{"A", SatSys(E_Sys::GPS, 13)};
    const auto l1 = zhangGraphOldArcAuxiliaryKey(
        state, E_ObsCode::L1C, edge, 0);
    const auto l2 = zhangGraphOldArcAuxiliaryKey(
        state, E_ObsCode::L2W, edge, 1);
    const auto nextL1 = zhangGraphOldArcAuxiliaryKey(
        state, E_ObsCode::L1C, edge, 1);
    BOOST_CHECK(l1 != l2);
    BOOST_CHECK(l1 != nextL1);
    BOOST_CHECK(l2 != nextL1);
}

BOOST_AUTO_TEST_CASE(signal_only_arc_is_a_chord_not_a_shared_tree_edge)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g13(E_Sys::GPS, 13);
    const SatSys g30(E_Sys::GPS, 30);
    const ZhangGraphEdge a01{"A", g01};
    const ZhangGraphEdge b01{"B", g01};
    const ZhangGraphEdge a13{"A", g13};
    const ZhangGraphEdge b13{"B", g13};
    const auto old = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, a13, b13}, "A",
        {a01, b01, b13});
    const auto common = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, a13}, "A");
    BOOST_REQUIRE(old.connected && common.connected);
    BOOST_REQUIRE(old.treeEdges.contains(b13));
    const auto next = zhangAugmentCommonTreeWithSignalChords(
        common, {b13, {"A", g30}});
    BOOST_CHECK(next.edges.contains(b13));
    BOOST_CHECK(!next.treeEdges.contains(b13));
    BOOST_CHECK(next.treeEdges == common.treeEdges);
    BOOST_CHECK(!next.satellites.contains(g30));

    KFState state;
    state.kfIndexMap.clear();
    auto add = [&](const KFKey& key)
    {
        state.kfIndexMap[key] = int(state.kfIndexMap.size());
        state.stateTransitionMap[key][key][0] = 1;
    };
    for (const auto code : {E_ObsCode::L1C, E_ObsCode::L2W})
    {
        add(zhangTransportReceiverKey(E_Sys::GPS, code, "B"));
        add(zhangTransportSatelliteKey(code, g01));
        add(zhangTransportSatelliteKey(code, g13));
        add(zhangTransportAmbiguityKey(code, a13));
    }
    const int n = int(state.kfIndexMap.size());
    state.x = Eigen::VectorXd::Zero(n);
    state.P = Eigen::MatrixXd::Identity(n, n);
    std::map<ZhangGraphEdge, std::set<E_ObsCode>> retired{
        {b13, {E_ObsCode::L2W}}};
    const std::map<int, std::set<ZhangGraphEdge>> destinationSignalEdges{
        {int(E_ObsCode::L1C), next.edges},
        {int(E_ObsCode::L2W), common.edges}};
    const auto plan = zhangPlanGraphCoordinateTransport(
        state, E_Sys::GPS, old, next, {b13},
        {{E_ObsCode::L1C, 0.1902936728},
         {E_ObsCode::L2W, 0.2442102134}},
        [](const auto&, auto) { return 100.; },
        [](const auto&, auto) { return 25.; },
        &retired, false, nullptr, nullptr, nullptr,
        &destinationSignalEdges);
    BOOST_REQUIRE_MESSAGE(plan.valid, plan.failureReason);
    BOOST_CHECK(plan.survivingPhysicalRows.contains(
        {E_ObsCode::L1C, b13}));
    BOOST_CHECK(!plan.survivingPhysicalRows.contains(
        {E_ObsCode::L2W, b13}));
    BOOST_CHECK(plan.transform.contains(
        zhangTransportAmbiguityKey(E_ObsCode::L1C, b13)));
    BOOST_CHECK(!plan.transform.contains(
        zhangTransportAmbiguityKey(E_ObsCode::L2W, b13)));

    const ZhangGraphEdge c30{"C", g30};
    const ZhangGraphEdge b30{"B", g30};
    const auto forest = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, c30}, "A");
    const auto isolated = zhangAugmentCommonTreeWithSignalChords(
        forest, {b30});
    BOOST_CHECK(!isolated.edges.contains(b30));
    BOOST_CHECK(isolated.treeEdges == forest.treeEdges);
}

BOOST_AUTO_TEST_CASE(signal_arc_current_support_uses_final_accepted_phase_rows)
{
    const ZhangGraphEdge a01{"A", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge b01{"B", SatSys(E_Sys::GPS, 1)};
    const std::map<int, std::set<ZhangGraphEdge>> preQc{
        {int(E_ObsCode::L1C), {a01, b01}},
        {int(E_ObsCode::L2W), {a01}}};
    const std::map<ZhangGraphEdge, std::map<int, int>> admissions{
        {a01, {{int(E_ObsCode::L1C), 4},
                {int(E_ObsCode::L2W), 2}}},
        {b01, {{int(E_ObsCode::L1C), 1}}}};
    KFMeas finalMeasurement;
    finalMeasurement.H = Eigen::MatrixXd::Ones(4, 1);
    finalMeasurement.R = Eigen::MatrixXd::Identity(4, 4);
    finalMeasurement.metaDataMaps.resize(4);
    unsigned int count = 0;
    auto phaseKey = [](const ZhangGraphEdge& edge,
                       E_ObsCode observable, KF type)
    {
        KFKey key;
        key.type = type;
        key.str = edge.receiver;
        key.Sat = edge.satellite;
        key.num = int(observable);
        return key;
    };
    finalMeasurement.obsKeys = {
        phaseKey(a01, E_ObsCode::L1C, KF::PHAS_MEAS),
        phaseKey(a01, E_ObsCode::L2W, KF::PHAS_MEAS),
        phaseKey(b01, E_ObsCode::L1C, KF::PHAS_MEAS),
        phaseKey(a01, E_ObsCode::L2W, KF::CODE_MEAS)};
    finalMeasurement.metaDataMaps[0]["phaseRejectCount"] = &count;
    finalMeasurement.metaDataMaps[1]["phaseRejectCount"] = &count;
    finalMeasurement.metaDataMaps[1]["zhangPppDeweighted"] = &count;
    finalMeasurement.metaDataMaps[2]["phaseRejectCount"] = nullptr;
    finalMeasurement.metaDataMaps[3]["phaseRejectCount"] = &count;
    const auto accepted = zhangFinalAdmittedSignalFactorRows(
        E_Sys::GPS, finalMeasurement, preQc, admissions);
    BOOST_REQUIRE(accepted.contains(int(E_ObsCode::L1C)));
    BOOST_CHECK(accepted.at(int(E_ObsCode::L1C)) ==
        std::set<ZhangGraphEdge>{a01});
    BOOST_CHECK(!accepted.contains(int(E_ObsCode::L2W)));

    finalMeasurement.metaDataMaps[0]["zhangPppDeweighted"] = &count;
    BOOST_CHECK(zhangFinalAdmittedSignalFactorRows(
        E_Sys::GPS, finalMeasurement, preQc, admissions).empty());
}

BOOST_AUTO_TEST_CASE(schema3_float_restart_first_graph_admission_is_checkpointable)
{
    const ZhangGraphEdge a01{"A", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge b01{"B", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge b13{"B", SatSys(E_Sys::GPS, 13)};
    const ZhangGraphEdge c30{"C", SatSys(E_Sys::GPS, 30)};
    const std::set<ZhangGraphEdge> oldBasisEdges{a01, b01, b13, c30};
    const std::set<ZhangGraphEdge> activeStateEdges{a01, b01, b13};
    const std::map<int, std::set<ZhangGraphEdge>> rawPhaseCatalogue{
        {int(E_ObsCode::L1C), oldBasisEdges},
        {int(E_ObsCode::L2W), oldBasisEdges}};
    // Schema 3 supplies valid FLOAT coordinates but no authoritative
    // per-signal integer versions. Only active current rows start new V0 arcs.
    std::map<ZhangGraphEdge, std::map<int, int>> signalVersions;
    BOOST_CHECK(signalVersions.empty());
    const auto preQc = zhangGraphSupportedSignalEdges(
        rawPhaseCatalogue, activeStateEdges, oldBasisEdges);
    for (const auto& [code, edges] : preQc)
        for (const auto& edge : edges)
            signalVersions[edge][code] = 0;
    BOOST_CHECK(!signalVersions.contains(c30));
    BOOST_REQUIRE(preQc.at(int(E_ObsCode::L1C)).size() == 3);

    KFMeas finalMeasurement;
    finalMeasurement.H = Eigen::MatrixXd::Ones(4, 1);
    finalMeasurement.R = Eigen::MatrixXd::Identity(4, 4);
    finalMeasurement.metaDataMaps.resize(4);
    unsigned int rejectCount = 0;
    const std::array<ZhangGraphEdge, 4> rows{a01, b01, b13, c30};
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        KFKey key;
        key.type = KF::PHAS_MEAS;
        key.str = rows[i].receiver;
        key.Sat = rows[i].satellite;
        key.num = int(E_ObsCode::L1C);
        finalMeasurement.obsKeys.push_back(key);
        finalMeasurement.metaDataMaps[i]["phaseRejectCount"] = &rejectCount;
    }
    std::map<ZhangGraphEdge, std::map<int, int>> admissions;
    for (const auto& edge : preQc.at(int(E_ObsCode::L1C)))
        admissions[edge][int(E_ObsCode::L1C)] = 0;
    const auto finalQc = zhangFinalAdmittedSignalFactorRows(
        E_Sys::GPS, finalMeasurement, preQc, admissions);
    const auto persisted = zhangGraphSupportedSignalEdges(
        finalQc, activeStateEdges, oldBasisEdges);
    std::set<ZhangGraphEdge> observationEdges;
    for (const auto& [code, edges] : persisted)
        observationEdges.insert(edges.begin(), edges.end());
    BOOST_CHECK(std::includes(activeStateEdges.begin(),
        activeStateEdges.end(), observationEdges.begin(),
        observationEdges.end()));
    BOOST_CHECK(!observationEdges.contains(c30));
    const auto productBasis = zhangBuildSpanningTree(activeStateEdges, "A");
    const auto products = zhangBuildProductIntegerFunctionals(
        productBasis, signalVersions, E_ObsCode::L1C,
        SatSys(E_Sys::GPS, 1));
    BOOST_CHECK(products.contains(SatSys(E_Sys::GPS, 13)));
}

BOOST_AUTO_TEST_CASE(final_qc_graph_admission_rechecks_retired_state_edges)
{
    const ZhangGraphEdge a01{"A", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge b01{"B", SatSys(E_Sys::GPS, 1)};
    const std::map<int, std::set<ZhangGraphEdge>> preQc{
        {int(E_ObsCode::L1C), {a01, b01}}};
    std::map<ZhangGraphEdge, std::map<int, int>> admissions{
        {a01, {{int(E_ObsCode::L1C), 0}}},
        {b01, {{int(E_ObsCode::L1C), 0}}}};
    KFMeas finalMeasurement;
    finalMeasurement.H = Eigen::MatrixXd::Ones(2, 1);
    finalMeasurement.R = Eigen::MatrixXd::Identity(2, 2);
    finalMeasurement.metaDataMaps.resize(2);
    unsigned int rejectCount = 0;
    for (const auto& edge : {a01, b01})
    {
        KFKey key;
        key.type = KF::PHAS_MEAS;
        key.str = edge.receiver;
        key.Sat = edge.satellite;
        key.num = int(E_ObsCode::L1C);
        finalMeasurement.obsKeys.push_back(key);
    }
    for (auto& meta : finalMeasurement.metaDataMaps)
        meta["phaseRejectCount"] = &rejectCount;
    const auto rawFinal = zhangFinalAdmittedSignalFactorRows(
        E_Sys::GPS, finalMeasurement, preQc, admissions);
    BOOST_REQUIRE(rawFinal.at(int(E_ObsCode::L1C)).contains(b01));
    const auto persisted = zhangGraphSupportedSignalEdges(
        rawFinal, {a01}, {a01, b01});
    BOOST_CHECK(persisted.at(int(E_ObsCode::L1C)) ==
        std::set<ZhangGraphEdge>{a01});
}

BOOST_AUTO_TEST_CASE(schema3_float_tree_rebase_needs_no_old_integer_arc_certificate)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g13(E_Sys::GPS, 13);
    const ZhangGraphEdge a01{"A", g01};
    const ZhangGraphEdge b01{"B", g01};
    const ZhangGraphEdge a13{"A", g13};
    const ZhangGraphEdge b13{"B", g13};
    const auto oldBasis = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, a13, b13}, "A",
        {a01, b01, b13});
    const auto destination = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, a13}, "A",
        {a01, b01, a13});
    BOOST_REQUIRE(oldBasis.connected && destination.connected);
    BOOST_REQUIRE(oldBasis.treeEdges.contains(b13));
    BOOST_REQUIRE(destination.treeEdges.contains(a13));

    KFState state;
    state.kfIndexMap.clear();
    for (const auto code : {E_ObsCode::L1C, E_ObsCode::L2W})
        for (const auto& key : {
             zhangTransportReceiverKey(E_Sys::GPS, code, "B"),
             zhangTransportSatelliteKey(code, g01),
             zhangTransportSatelliteKey(code, g13),
             zhangTransportAmbiguityKey(code, a13)})
        {
            state.kfIndexMap[key] = int(state.kfIndexMap.size());
            state.stateTransitionMap[key][key][0] = 1;
        }
    const int n = int(state.kfIndexMap.size());
    state.x = Eigen::VectorXd::LinSpaced(n, 0.2, 1.6);
    state.P = Eigen::MatrixXd::Identity(n, n) * 4;
    state.P(2, 6) = state.P(6, 2) = 0.35;
    state.P(3, 7) = state.P(7, 3) = -0.2;
    const std::map<E_ObsCode, double> wavelengths{
        {E_ObsCode::L1C, 0.1902936728},
        {E_ObsCode::L2W, 0.2442102134}};
    const std::map<ZhangGraphEdge, std::set<E_ObsCode>> retiredSignals{
        {b13, {E_ObsCode::L1C, E_ObsCode::L2W}}};
    const std::map<int, std::set<ZhangGraphEdge>> destinationSignals{
        {int(E_ObsCode::L1C), destination.edges},
        {int(E_ObsCode::L2W), destination.edges}};
    const auto plan = zhangPlanGraphCoordinateTransport(
        state, E_Sys::GPS, oldBasis, destination, {b13}, wavelengths,
        [](const auto&, auto) { return 100.; },
        [](const auto&, auto) { return 25.; },
        &retiredSignals, false, nullptr, nullptr, nullptr,
        &destinationSignals);
    BOOST_REQUIRE_MESSAGE(plan.valid, plan.failureReason);
    BOOST_CHECK(plan.retainedOldArcAuxiliaryKeys.empty());
    for (const auto& [destinationKey, row] : plan.transform)
        BOOST_CHECK(destinationKey.type != KF::PRODUCT_HANDOFF);
    const auto impossibleLegacyCertificate = zhangPlanGraphCoordinateTransport(
        state, E_Sys::GPS, oldBasis, destination, {b13}, wavelengths,
        [](const auto&, auto) { return 100.; },
        [](const auto&, auto) { return 25.; },
        &retiredSignals, true, nullptr, nullptr, nullptr,
        &destinationSignals);
    BOOST_CHECK(!impossibleLegacyCertificate.valid);
    BOOST_CHECK_EQUAL(impossibleLegacyCertificate.failureReason,
        "OLD_ARC_VERSION_MAP_MISSING");

    std::vector<std::pair<ZhangCoordinateRow, double>> newRows;
    Eigen::MatrixXd oldRows = Eigen::MatrixXd::Zero(2, n);
    int signal = 0;
    for (const auto& [code, wavelength] : wavelengths)
    {
        const KFKey satellite = zhangTransportSatelliteKey(code, g13);
        const KFKey ambiguity = zhangTransportAmbiguityKey(code, a13);
        const ZhangCoordinateRow oldRow{{satellite, 1},
                                        {ambiguity, wavelength}};
        const ZhangCoordinateRow newRow{{satellite, 1}};
        BOOST_CHECK(zhangVerifyGraphFunctionalTransport(
            plan, oldRow, newRow, 0, 0).sameTimeMapProven);
        oldRows(signal, state.kfIndexMap.at(satellite)) = 1;
        oldRows(signal, state.kfIndexMap.at(ambiguity)) = wavelength;
        newRows.emplace_back(newRow, 0);
        ++signal;
    }
    const auto moments = zhangGraphFunctionalMoments(state, plan, newRows);
    BOOST_REQUIRE(moments.valid);
    BOOST_CHECK_SMALL((moments.mean - oldRows * state.x).norm(), 1e-10);
    BOOST_CHECK_SMALL((moments.covariance -
        oldRows * state.P * oldRows.transpose()).norm(), 1e-10);
}

BOOST_AUTO_TEST_CASE(full_float_cycle_needs_known_arc_outside_product_core)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g13(E_Sys::GPS, 13);
    const ZhangGraphEdge a01{"A", g01};
    const ZhangGraphEdge b01{"B", g01};
    const ZhangGraphEdge a13{"A", g13};
    const ZhangGraphEdge b13{"B", g13};
    const auto floatBasis = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, a13, b13}, "A",
        {a01, b01, a13});
    BOOST_REQUIRE(floatBasis.connected);
    BOOST_REQUIRE(!floatBasis.treeEdges.contains(b13));
    const auto productCore = zhangBuildSpanningTree(
        std::set<ZhangGraphEdge>{a01, b01, a13}, "A");
    BOOST_REQUIRE(!productCore.edges.contains(b13));
    std::map<ZhangGraphEdge, int> productOnlyVersions;
    for (const auto& edge : productCore.edges)
        productOnlyVersions[edge] = 2;
    auto fullKnownVersions = productOnlyVersions;
    fullKnownVersions[b13] = 3;
    ZhangProductPhysicalCycleChart incomplete;
    incomplete.columns = 1;
    BOOST_CHECK(!incomplete.add(0, "L1C", b13, floatBasis,
        productOnlyVersions));
    ZhangProductPhysicalCycleChart complete;
    complete.columns = 1;
    BOOST_REQUIRE(complete.add(0, "L1C", b13, floatBasis,
        fullKnownVersions));
    BOOST_CHECK_EQUAL(complete.chordIdentities.at(0),
        "L1C|B|G13|V3");
    BOOST_CHECK(complete.activePhysicalIdentities.contains(
        "L1C|B|G13|V3"));
}
