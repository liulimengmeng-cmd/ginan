#include <boost/test/unit_test.hpp>
#include "common/zhangProductPhysicalPullback.hpp"

namespace {
ZhangGraphBasis r47Graph(bool pivot)
{
    ZhangGraphBasis b;
    b.rootReceiver = "A";
    b.receivers = {"A", "B"};
    b.satellites = {SatSys("G01"), SatSys("G02"), SatSys("G03")};
    for (const auto& r : b.receivers)
        for (const auto& s : b.satellites) b.edges.insert({r, s});
    b.treeEdges = {{"A",SatSys("G01")}, {"A",SatSys("G02")},
        {"A",SatSys("G03")}, {"B",SatSys(pivot ? "G02" : "G01")}};
    b.connected = true;
    b.componentCount = 1;
    return b;
}
ZhangProductRelationBasis r47Target(const ZhangGraphBasis& privateGraph)
{
    ZhangProductRelationBasis target;
    target.observable = E_ObsCode::L1C;
    ZhangProductRelationRow row;
    row.satellite = SatSys("G03");
    row.referenceSatellite = SatSys("G01");
    for (auto [edge, value] : zhangFundamentalCycle(privateGraph, {"B",SatSys("G03")}))
        row.physicalArcCoefficients[edge] = value;
    target.namedRelations.push_back(row);
    return target;
}
}

BOOST_AUTO_TEST_CASE(r47_a_private_target_pulls_back_before_search_and_constraint_binding)
{
    auto target = r47Target(r47Graph(false));
    const auto intended = target.namedRelations[0].physicalArcCoefficients;
    const auto authoritative = r47Graph(true);
    std::map<ZhangGraphEdge,int> versions;
    for (auto edge : authoritative.edges) versions[edge] = 7;
    std::vector<ZhangGraphEdge> chords = {{"B",SatSys("G01")}, {"B",SatSys("G03")}};
    const auto result = zhangPullbackProductPhysicalRelations(target,
        authoritative, versions, chords, {{chords[0],0},{chords[1],1}}, 2);
    BOOST_REQUIRE(result.valid);
    BOOST_REQUIRE(result.available[0]);
    BOOST_CHECK(result.posteriorRows[0] == ZhangExactVector({-1,1}));
    BOOST_CHECK(target.namedRelations[0].currentCycleCoefficients == result.posteriorRows[0]);
    BOOST_CHECK(target.currentChords == chords);
    BOOST_CHECK(target.namedRelations[0].physicalArcCoefficients == intended);
    Eigen::VectorXd row;
    BOOST_REQUIRE(zhangExactPosteriorRowToDouble(result.posteriorRows[0], row));
    Eigen::Vector2d mu(-4,5);
    Eigen::Matrix2d covariance;
    covariance << 4,1,1,9;
    BOOST_CHECK_SMALL(row.dot(mu) - 9, 1e-12);
    BOOST_CHECK_SMALL((row.transpose()*covariance*row)(0,0) - 11, 1e-12);
    // The former chord-name lookup would use (0,1), yielding 5 and 9.
    BOOST_CHECK_NE(row.dot(mu), mu(1));
    // Rz=k with d=3 gives RTa=k-Rd, including non-unit R.
    const auto rhs = ZhangExactInteger(24) - ZhangExactInteger(2)*3;
    BOOST_CHECK_EQUAL(rhs, 18);
    BOOST_CHECK_SMALL((2*row).dot(mu) - rhs.convert_to<double>(), 1e-12);
}

BOOST_AUTO_TEST_CASE(r47_a_rectangular_private_subgraph_does_not_require_unimodular_inverse)
{
    auto subgraph = r47Graph(false);
    subgraph.edges.erase({"A",SatSys("G02")});
    subgraph.edges.erase({"B",SatSys("G02")});
    subgraph.treeEdges.erase({"A",SatSys("G02")});
    subgraph.satellites.erase(SatSys("G02"));
    auto target = r47Target(subgraph);
    auto graph = r47Graph(true);
    std::map<ZhangGraphEdge,int> versions;
    for (auto edge : graph.edges) versions[edge] = 0;
    std::vector<ZhangGraphEdge> chords = {{"B",SatSys("G01")}, {"B",SatSys("G03")}};
    auto result = zhangPullbackProductPhysicalRelations(target, graph, versions,
        chords, {{chords[0],0},{chords[1],1}}, 2);
    BOOST_REQUIRE(result.valid);
    BOOST_REQUIRE(result.available[0]);
    BOOST_CHECK_EQUAL(result.posteriorRows.size(), 1);
    BOOST_CHECK(result.posteriorRows[0] == ZhangExactVector({-1,1}));
    // Missing G01 cannot be hidden by relabelling private G03 as current G03.
    result = zhangPullbackProductPhysicalRelations(target, graph, versions,
        chords, {{chords[1],1}}, 2);
    BOOST_REQUIRE(result.valid);
    BOOST_CHECK(!result.available[0]);
    target.namedRelations[0].physicalArcCoefficients.begin()->second += 1;
    result = zhangPullbackProductPhysicalRelations(target, graph, versions,
        chords, {{chords[0],0},{chords[1],1}}, 2);
    BOOST_CHECK(!result.valid);
}

BOOST_AUTO_TEST_CASE(r47_a_numeric_integer_coefficients_must_be_exactly_representable)
{
    Eigen::VectorXd numeric;
    const ZhangExactInteger huge = (ZhangExactInteger(1) << 53) + 1;
    BOOST_CHECK(!zhangExactPosteriorRowToDouble({huge}, numeric));
    BOOST_REQUIRE(zhangExactPosteriorRowToDouble({-7,0,9}, numeric));
    BOOST_CHECK_EQUAL(numeric(0), -7);
}
