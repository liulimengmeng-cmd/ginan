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

#include <filesystem>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>
namespace {
struct R47Fixture {
    std::vector<std::string> names, columns;
    ZhangExactMatrix rows;
};
R47Fixture r47ReadFixture(const std::string& file, const std::string& key)
{
    boost::property_tree::ptree input;
    boost::property_tree::read_json((std::filesystem::path(__FILE__).parent_path().parent_path()
        / "fixtures" / "r47" / file).string(), input);
    R47Fixture f;
    std::map<std::string,int> columns;
    std::vector<std::map<std::string,int>> sparse;
    for (const auto& [name, node] : input.get_child(key))
    {
        f.names.push_back(name);
        std::map<std::string,int> row;
        for (const auto& entry : node)
        {
            auto it = entry.second.begin();
            const auto receiver = (it++)->second.get_value<std::string>();
            const auto sat = (it++)->second.get_value<std::string>();
            const int value = it->second.get_value<int>();
            const auto arc = receiver+"/"+sat;
            row[arc] = value;
            if (!columns.count(arc)) { columns[arc] = f.columns.size(); f.columns.push_back(arc); }
        }
        sparse.push_back(row);
    }
    for (const auto& sparseRow : sparse)
    {
        ZhangExactVector row(columns.size());
        for (const auto& [arc,value] : sparseRow) row[columns.at(arc)] = value;
        f.rows.push_back(row);
    }
    return f;
}
}
BOOST_AUTO_TEST_CASE(r47_b_actual_2100_whole_lattice_recovers_six_cancelled_directions)
{
    const auto f = r47ReadFixture("availability_2100_exact_cancellation.json", "full_targets");
    std::vector<bool> available(f.columns.size(), true);
    const auto missing = std::find(f.columns.begin(), f.columns.end(), "KIRI/G27") - f.columns.begin();
    available.at(missing) = false;
    const auto domain = zhangCompileWholeProductLattice(f.rows, ZhangExactVector(f.rows.size()), available);
    BOOST_REQUIRE(domain.valid);
    BOOST_CHECK_EQUAL(domain.directlyAvailableRank, 16);
    BOOST_CHECK_EQUAL(domain.wholeAvailableRank, 22);
    BOOST_CHECK_EQUAL(domain.namedExpressions.size(), 28);
    BOOST_CHECK(domain.primitiveSearchImage);
    const auto base = std::find(f.names.begin(), f.names.end(), "G10") - f.names.begin();
    for (const std::string sat : {"G16","G18","G23","G26","G27","G29"})
    {
        const auto n = std::find(f.names.begin(), f.names.end(), sat) - f.names.begin();
        BOOST_CHECK(!domain.namedRecoverable.at(n));
        ZhangExactVector difference;
        for (std::size_t c = 0; c < available.size(); ++c)
            if (available[c]) difference.push_back(f.rows[n][c]-f.rows[base][c]);
        BOOST_CHECK_EQUAL(f.rows[n][missing]-f.rows[base][missing], 0);
        BOOST_CHECK(zhangIntegerRowLatticeContains(domain.searchPosteriorRows, difference).contained);
    }
}
BOOST_AUTO_TEST_CASE(r47_b_actual_2730_catalogue_retains_all_dependent_targets)
{
    const auto f = r47ReadFixture("r46_2730_target_catalogue.json", "current_target_rows");
    const auto domain = zhangCompileWholeProductLattice(f.rows, ZhangExactVector(f.rows.size()),
        std::vector<bool>(f.columns.size(), true));
    BOOST_REQUIRE(domain.valid);
    BOOST_CHECK_EQUAL(domain.targetIndependentRank, 22);
    BOOST_CHECK_EQUAL(domain.namedExpressions.size(), 28);
    BOOST_CHECK_EQUAL(domain.structuralIdentityRows.size(), 6);
    BOOST_CHECK(zhangExactMultiply(domain.namedRecovery, domain.searchPosteriorRows) == f.rows);
    for (const std::string sat : {"G12","G19","G24","G25","G26","G32"})
    {
        const auto n = std::find(f.names.begin(), f.names.end(), sat) - f.names.begin();
        BOOST_CHECK(domain.namedRecoverable.at(n));
        BOOST_CHECK_EQUAL(domain.namedStatus.at(n), "REPRESENTABLE_DEPENDENT");
    }
}
BOOST_AUTO_TEST_CASE(r47_b_affine_zero_identities_and_divisibility_are_not_ar_proofs)
{
    const auto d = zhangCompileWholeProductLattice({{2,1},{2,1},{0,0}}, {3,7,5}, {true,false});
    BOOST_REQUIRE(d.valid);
    BOOST_CHECK_EQUAL(d.wholeAvailableRank, 0);
    BOOST_CHECK_EQUAL(d.structuralIdentityRows.size(), 2);
    BOOST_CHECK(d.namedRecoverable[2]);
    BOOST_CHECK_EQUAL(d.namedRecoveryOffsets[2], 5);
    auto q = zhangCompileWholeProductLattice({{2,0}}, {3}, {true,true});
    BOOST_REQUIRE(q.valid);
    BOOST_CHECK(!q.primitiveSearchImage);
    BOOST_CHECK_EQUAL(q.searchSmithInvariants.front(), 2);
}
