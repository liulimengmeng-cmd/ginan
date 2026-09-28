#include <boost/test/unit_test.hpp>
#include <numeric>
#include "common/zhangProductPhysicalPullback.hpp"
#include "common/zhangProductAffineRecovery.hpp"
#include "common/zhangParallelProjection.hpp"
#include "common/zhangParallelProductProjection.hpp"

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
    BOOST_REQUIRE_EQUAL(result.diagnostics.size(), 1);
    BOOST_CHECK(!result.diagnostics[0].requested.empty());
    BOOST_CHECK(result.diagnostics[0].rebuilt != result.diagnostics[0].requested);
    BOOST_CHECK(!result.diagnostics[0].delta.empty());
    BOOST_REQUIRE_EQUAL(result.missingRequiredChords[0].size(), 1);
    BOOST_CHECK(result.missingRequiredChords[0][0] == chords[0]);
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

#include "common/zhangR47Candidate.hpp"
BOOST_AUTO_TEST_CASE(r47_c_mixed_conditioner_changes_product_without_certifying_it)
{
    Eigen::Vector2d mean(3.20,2.98);
    Eigen::Matrix2d covariance=Eigen::Vector2d(.04,.0025).asDiagonal();
    Eigen::Matrix<double,1,2> h; h<<1,1;
    Eigen::VectorXd value(1); value<<6;
    const auto conditioned=zhangConditionPosteriorEffectiveIntegers(mean,covariance,h,value);
    BOOST_REQUIRE(conditioned.valid);
    BOOST_CHECK_SMALL(conditioned.mean(0)-1288.0/425,1e-12);
    BOOST_CHECK_SMALL(conditioned.covariance(0,0)-1.0/425,1e-12);
    const auto image=zhangAppliedLatticeProductImage({{1,1}},{6},{{1,0}},{0});
    BOOST_REQUIRE(image.consistent);
    BOOST_CHECK(image.basis.empty());
    Eigen::Matrix2d uncoupled=Eigen::Vector2d(.04,.0025).asDiagonal();
    h<<0,1; value<<3;
    const auto noGain=zhangConditionPosteriorEffectiveIntegers(mean,uncoupled,h,value);
    BOOST_REQUIRE(noGain.valid);
    BOOST_CHECK_EQUAL(noGain.mean(0),mean(0));
    BOOST_CHECK_EQUAL(noGain.covariance(0,0),uncoupled(0,0));
}
BOOST_AUTO_TEST_CASE(r47_c_nonprimitive_equalities_require_integer_feasibility)
{
    BOOST_CHECK(zhangR47AffineIntegerFeasible({{2}},{6},1));
    BOOST_CHECK(!zhangR47AffineIntegerFeasible({{2}},{5},1));
    const auto q=zhangExactAffineIntegerQuotient({{2}},{6},1);
    BOOST_REQUIRE(q.valid);
    BOOST_CHECK_EQUAL(q.particularSolution[0],3);
    BOOST_CHECK_EQUAL(q.quotientRank,0);
    BOOST_CHECK(!zhangExactAffineIntegerQuotient({{2}},{5},1).valid);
    BOOST_CHECK(!zhangR47AffineIntegerFeasible({{1,1},{1,-1}},{1,0},2));
}
BOOST_AUTO_TEST_CASE(r47_c_full_state_conditioning_retains_cross_covariance_and_rolls_back)
{
    Eigen::Vector3d mean(3.20,2.98,10);
    Eigen::Matrix3d covariance;
    covariance<<.04,0,.01, 0,.0025,.001, .01,.001,.1;
    Eigen::Matrix<double,1,3> h; h<<1,1,0;
    Eigen::VectorXd value(1); value<<6;
    const auto result=zhangConditionPosteriorEffectiveIntegers(mean,covariance,h,value);
    BOOST_REQUIRE(result.valid);
    BOOST_CHECK_SMALL(result.mean(2)-(10-.011*.18/.0425),1e-12);
    BOOST_CHECK_LT(result.covariance(2,2),covariance(2,2));
    Eigen::Matrix<double,2,3> conflict; conflict<<1,1,0,1,1,0;
    Eigen::Vector2d bad(6,7);
    BOOST_CHECK(!zhangConditionPosteriorEffectiveIntegers(mean,covariance,conflict,bad).valid);
    BOOST_CHECK_EQUAL(mean(2),10); // root never mutated, no inverse update.
}
BOOST_AUTO_TEST_CASE(r47_c_immutable_candidate_rejects_dropped_conditioner_or_parent)
{
    ZhangR47Candidate c;
    c.sourcePosteriorId="root"; c.authoritativeCycleChartId="generation|columns";
    c.frontendSemanticId="beta+kappa"; c.sourceIntegerParentCount=0; c.stateIndices={0,1};
    c.admittedStateConditioners={{1,1}}; c.admittedStateValues={6};
    c.newIntegerConstraints={{1,0}}; c.newIntegerValues={3};
    c.jointRows={{1,0},{0,1}}; c.jointValues={3,3};
    c.physicalFunctionals={{{"arc1",1}},{{"arc2",1}}};
    c.allDecisionParents={std::make_shared<const ZhangIntegerDecisionProof>(
        ZhangIntegerDecisionProof{"decision","integer","epoch",1e-4,{}})};
    c.riskBound=1e-4;
    BOOST_REQUIRE(zhangR47CandidateContractValid(c));
    auto dropped=c; dropped.jointRows={{1,0}}; dropped.jointValues={3};
    dropped.physicalFunctionals.resize(1);
    BOOST_CHECK(!zhangR47CandidateContractValid(dropped));
    dropped=c; dropped.allDecisionParents.clear();
    BOOST_CHECK(!zhangR47CandidateContractValid(dropped));
    dropped=c; dropped.sourceIntegerParentCount=1;
    BOOST_CHECK(!zhangR47CandidateContractValid(dropped));

    c.productConsequences.fullJointProductMappingExact=true;
    c.productConsequences.fullJointProductNetworkRows={{1,0},{0,1}};
    c.productConsequences.fullJointProductAffineOffsets={0,7};
    c.productConsequences.networkRows=c.jointRows;
    c.productConsequences.networkIntegers=c.jointValues;
    c.productConsequences.jointProductRows={{1,0},{0,1}};
    c.productConsequences.decisionProofs=c.allDecisionParents;
    BOOST_REQUIRE(zhangR47CandidateProductChartMatches(
        c,c.productConsequences));
    auto changedChart=c.productConsequences;
    changedChart.fullJointProductAffineOffsets[1]=8;
    BOOST_CHECK(!zhangR47CandidateProductChartMatches(c,changedChart));
    changedChart=c.productConsequences;
    changedChart.decisionProofs={std::make_shared<const ZhangIntegerDecisionProof>(
        ZhangIntegerDecisionProof{"other","integer","epoch",1e-4,{}})};
    BOOST_CHECK(!zhangR47CandidateProductChartMatches(c,changedChart));
}

#include "common/zhangR47History.hpp"
BOOST_AUTO_TEST_CASE(r47_d_retired_physical_arcs_cancel_before_projection_with_both_parents)
{
    const auto graph=r47Graph(true);
    std::map<ZhangGraphEdge,int> versions; for(auto edge:graph.edges) versions[edge]=2;
    ZhangProductPhysicalCycleChart chart; chart.columns=2;
    BOOST_REQUIRE(chart.add(0,"L1C",{"B",SatSys("G01")},graph,versions));
    BOOST_REQUIRE(chart.add(1,"L1C",{"B",SatSys("G03")},graph,versions));
    auto first=chart.expansions.at(0), second=chart.expansions.at(1);
    first["L1C|OLD|G04|V1"]=1;second["L1C|OLD|G04|V1"]=1;
    ZhangExactVector ignored;
    BOOST_CHECK(!chart.project(first,ignored));BOOST_CHECK(!chart.project(second,ignored));
    auto a=std::make_shared<const ZhangIntegerDecisionProof>(ZhangIntegerDecisionProof{"a","first","old",1e-4,{}});
    auto b=std::make_shared<const ZhangIntegerDecisionProof>(ZhangIntegerDecisionProof{"b","second","old",2e-4,{}});
    const auto result=zhangR47TransportHistory({first,second},{4,1},{{a},{b}},chart);
    BOOST_REQUIRE(result.valid);BOOST_REQUIRE_EQUAL(result.rows.size(),1);
    const auto membership=zhangIntegerRowLatticeContains(result.rows,{1,-1});
    BOOST_REQUIRE(membership.contained);
    BOOST_CHECK_EQUAL(membership.combination[0]*result.values[0],3);
    const auto risk=zhangDecisionRiskClosure(result.parents[0]);
    BOOST_REQUIRE(risk.valid);BOOST_CHECK_EQUAL(risk.atoms.size(),2);
    BOOST_CHECK_SMALL(risk.bound-3e-4,1e-15);
}
BOOST_AUTO_TEST_CASE(r47_d_budget_subset_deduplicates_ancestors_and_reserves_search_risk)
{
    auto a=std::make_shared<const ZhangIntegerDecisionProof>(ZhangIntegerDecisionProof{"a","first","old",6e-4,{}});
    auto b=std::make_shared<const ZhangIntegerDecisionProof>(ZhangIntegerDecisionProof{"b","second","old",1e-4,{a}});
    auto c=std::make_shared<const ZhangIntegerDecisionProof>(ZhangIntegerDecisionProof{"c","third","old",2e-4,{}});
    const auto choice=zhangR47SelectHistorySubset({{1,0,0},{0,1,0},{0,0,1}},
        {2,3,4},{{a},{b},{c}},{},{3,2,1},3,1e-3,2.5e-4);
    BOOST_REQUIRE(choice.valid);BOOST_REQUIRE_EQUAL(choice.selected.size(),2);
    BOOST_CHECK_SMALL(choice.risk-7e-4,1e-15);
    BOOST_CHECK_EQUAL(choice.reasons[2],"BUDGET_RESERVED_FOR_NEW_SEARCH");
    BOOST_CHECK_EQUAL(zhangDecisionRiskClosure(choice.parents).atoms.size(),2);
    const auto conflict=zhangR47SelectHistorySubset({{1,0},{1,0},{0,1}},
        {2,5,3},{{a},{a},{a}},{},{3,2,1},2,1e-3,2.5e-4);
    BOOST_REQUIRE(conflict.valid);
    BOOST_CHECK_EQUAL(conflict.reasons[1],"EXACT_AFFINE_CONFLICT");
    BOOST_CHECK_EQUAL(conflict.selected.size(),2); // shared risk is zero increment, not automatic acceptance.
}

#include "common/zhangR47ProductDomain.hpp"
#include "common/zhangSequentialQuotientShadow.hpp"
BOOST_AUTO_TEST_CASE(r47_e_nonprimitive_product_image_preserves_even_and_odd_cosets)
{
    const auto even=zhangR47CompileProductSearchFrame({{2,0}}, {},{},2);
    BOOST_REQUIRE(even.valid);BOOST_CHECK_EQUAL(even.searchRank,1);
    BOOST_CHECK_EQUAL(zhangExactAbs(even.imageGenerators[0][0]),2);
    const auto odd=zhangR47CompileProductSearchFrame({{1,-1}},{{1,1}},{1},2);
    BOOST_REQUIRE(odd.valid);BOOST_REQUIRE_EQUAL(odd.searchRank,1);
    BOOST_CHECK_EQUAL(zhangExactAbs(odd.imageGenerators[0][0]),2);
    auto rows=ZhangExactMatrix{{1,1}};auto values=ZhangExactVector{1};
    rows.push_back(odd.projector[0]);values.push_back(2-odd.offsets[0]);
    const auto final=zhangR47CompileProductSearchFrame({{1,-1}},rows,values,2);
    BOOST_REQUIRE(final.valid);BOOST_CHECK_EQUAL(final.searchRank,0);
    ZhangExactInteger consequence;
    BOOST_REQUIRE(zhangR47ProductConsequence(final,{1,-1},0,consequence));
    BOOST_CHECK_EQUAL(zhangExactAbs(consequence)%2,1);
    BOOST_CHECK(!zhangR47CompileProductSearchFrame({{1,0}},{{2,0}},{5},2).valid);
}
BOOST_AUTO_TEST_CASE(r47_e_mixed_history_reduces_product_uncertainty_without_searching_network_nullspace)
{
    const auto f=zhangR47CompileProductSearchFrame({{1,0}},{{1,1}},{6},2);
    BOOST_REQUIRE(f.valid);BOOST_REQUIRE_EQUAL(f.searchRank,1);
    ZhangExactInteger value;
    BOOST_CHECK(!zhangR47ProductConsequence(f,{1,0},0,value));
    ZhangExactVector target(1000);target[517]=3;
    const auto sparse=zhangR47CompileProductSearchFrame({target},{},{},1000);
    BOOST_REQUIRE(sparse.valid);BOOST_CHECK_EQUAL(sparse.columns.size(),1);
    BOOST_CHECK_EQUAL(sparse.searchRank,1);BOOST_CHECK_EQUAL(sparse.columns[0],517);
    BOOST_CHECK_EQUAL(zhangExactAbs(sparse.imageGenerators[0][0]),3);
}
BOOST_AUTO_TEST_CASE(r47_e_direct_projection_reuses_only_the_same_feasible_domain)
{
    const ZhangExactMatrix held{{1,-1,0},{0,0,2}};
    const ZhangExactVector values{1,4};
    const auto domain=zhangR47CompileProductSearchFrame({{1,1,1}},held,values,3);
    BOOST_REQUIRE(domain.valid);
    const ZhangExactMatrix targets{{1,0,0},{0,1,0}};
    const auto projected=zhangR49CompileTargetOnDomain(targets,domain,3);
    const auto rebuilt=zhangR47CompileProductSearchFrame(targets,held,values,3);
    BOOST_REQUIRE(projected.valid && rebuilt.valid);
    BOOST_CHECK(projected.affine.get()==domain.affine.get());
    BOOST_CHECK_EQUAL(projected.searchRank,rebuilt.searchRank);
    BOOST_CHECK(projected.imageGenerators==rebuilt.imageGenerators);
    BOOST_CHECK(projected.projector==rebuilt.projector);
    BOOST_CHECK(projected.offsets==rebuilt.offsets);

    // Appending an integer decision changes H and b. The old projected
    // frame must not stand in for a newly solved affine domain.
    auto extended=held;extended.push_back({1,0,0});
    auto rhs=values;rhs.push_back(5);
    const auto afterFix=zhangR47CompileProductSearchFrame(targets,extended,rhs,3);
    BOOST_REQUIRE(afterFix.valid);
    BOOST_CHECK_EQUAL(afterFix.searchRank,0);
    BOOST_CHECK(afterFix.affine.get()!=domain.affine.get());
    extended.push_back({1,0,0});rhs.push_back(6);
    BOOST_CHECK(!zhangR47CompileProductSearchFrame(targets,extended,rhs,3).valid);

    // Even the zero-column shortcut must reject an inconsistent RHS.
    BOOST_CHECK(!zhangR47CompileProductSearchFrame({{0}},{{0}},{1},1).valid);
    BOOST_CHECK(zhangR47CompileProductSearchFrame({{0}},{{0}},{0},1).valid);
}
BOOST_AUTO_TEST_CASE(r51_parallel_sparse_moments_equal_dense_reference)
{
    Eigen::MatrixXd a=Eigen::MatrixXd::Zero(96,160);
    for(int r=0;r<a.rows();++r) {
        a(r,(7*r)%a.cols())=1;
        a(r,(7*r+19)%a.cols())=-2;
    }
    Eigen::MatrixXd root(160,160);
    for(int r=0;r<root.rows();++r)
    for(int c=0;c<root.cols();++c)
        root(r,c)=r==c?2.0:(r+c)%13==0?0.01:0.0;
    const Eigen::MatrixXd p=(root*root.transpose()).eval();
    const ZhangParallelRowProjection projected(a);
    BOOST_REQUIRE(projected.sparse);
    const Eigen::MatrixXd cross=projected.multiply(p);
    const Eigen::MatrixXd covariance=projected.covariance(cross);
    const Eigen::VectorXd mean=Eigen::VectorXd::LinSpaced(160,-1,1);
    BOOST_CHECK_SMALL((cross-a*p).cwiseAbs().maxCoeff(),1e-11);
    BOOST_CHECK_SMALL((covariance-a*p*a.transpose()).cwiseAbs().maxCoeff(),1e-10);
    BOOST_CHECK_SMALL((projected.multiply(mean)-a*mean).cwiseAbs().maxCoeff(),1e-12);

    const Eigen::MatrixXd dense=Eigen::MatrixXd::Constant(4,4,1);
    const ZhangParallelRowProjection fallback(dense);
    const Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(4,4);
    BOOST_CHECK(!fallback.sparse);
    BOOST_CHECK_SMALL((fallback.covariance(fallback.multiply(identity))
        -dense*dense.transpose()).cwiseAbs().maxCoeff(),1e-12);

    const Eigen::MatrixXd whitened=(a*root).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(
        whitened*whitened.transpose());
    BOOST_REQUIRE(eigen.info()==Eigen::Success);
    const double tolerance=std::max(1e-14,1e-12*eigen.eigenvalues().maxCoeff());
    const auto batched=zhangR51ConstraintRightBasis(whitened,
        eigen.eigenvectors(),eigen.eigenvalues(),tolerance);
    Eigen::MatrixXd serialBasis(whitened.cols(),batched.cols());
    int column=0;
    for(int i=0;i<eigen.eigenvalues().size();++i)
        if(eigen.eigenvalues()(i)>tolerance)
            serialBasis.col(column++)=whitened.transpose()*eigen.eigenvectors().col(i)
                /std::sqrt(eigen.eigenvalues()(i));
    BOOST_REQUIRE_EQUAL(column,batched.cols());
    BOOST_CHECK_SMALL((batched-serialBasis).cwiseAbs().maxCoeff(),1e-11);
}
BOOST_AUTO_TEST_CASE(r51_parallel_exact_product_projection_equals_serial_order)
{
    const int products=48,dimension=64,generators=128;
    ZhangExactMatrix targets(products,ZhangExactVector(dimension));
    ZhangExactMatrix kernel(generators,ZhangExactVector(dimension));
    std::vector<int> columns(dimension);
    std::iota(columns.begin(),columns.end(),0);
    for(int r=0;r<products;++r)
    for(int c=0;c<dimension;++c)
        if((r+3*c)%4!=0)targets[r][c]=(r+c)%7-3;
    for(int k=0;k<generators;++k)
    for(int c=0;c<dimension;++c)kernel[k][c]=(k+2*c)%5-2;
    const auto projected=zhangProjectTargetsOnKernel(targets,columns,kernel);
    ZhangExactMatrix serial(generators,ZhangExactVector(products));
    for(int r=0;r<products;++r)
    for(int k=0;k<generators;++k)
    for(int c=0;c<dimension;++c)
        serial[k][r]+=targets[r][columns[c]]*kernel[k][c];
    BOOST_CHECK(projected==serial);
    BOOST_CHECK_THROW(zhangProjectTargetsOnKernel(targets,{dimension},kernel),
        std::invalid_argument);
}
BOOST_AUTO_TEST_CASE(r47_e_actual_2100_search_domain_has_22_dimensions_after_missing_arc_elimination)
{
    const auto f=r47ReadFixture("availability_2100_exact_cancellation.json","full_targets");
    std::vector<bool> available(f.columns.size(),true);
    available[std::find(f.columns.begin(),f.columns.end(),"KIRI/G27")-f.columns.begin()]=false;
    const auto domain=zhangCompileWholeProductLattice(f.rows,ZhangExactVector(f.rows.size()),available);
    const auto frame=zhangR47CompileProductSearchFrame(domain.searchPosteriorRows,{}, {},
        std::count(available.begin(),available.end(),true));
    BOOST_REQUIRE(frame.valid);BOOST_CHECK_EQUAL(frame.searchRank,22);
}
BOOST_AUTO_TEST_CASE(r47_e_official_round_receipts_account_for_all_attempts)
{
    Eigen::Vector4d mean(1,2,3,4);Eigen::Matrix4d q=Eigen::Matrix4d::Identity()*1e-4;
    const auto result=zhangSequentialQuotientShadow(mean,q,{}, {},1e-4,1e-6,
        [](const Eigen::VectorXd& mean,const Eigen::MatrixXd&,double budget,bool)
        {
            ZhangSequentialShadowProposal p;p.rows=zhangExactIdentityMatrix(mean.size());
            for(int i=0;i<mean.size();++i) p.values.push_back(std::llround(mean[i]));
            p.valid=true;p.failureProbability=budget/2;return p;
        },[](const std::string&){},2,3,4);
    BOOST_CHECK_EQUAL(result.rows.size(),4);
    BOOST_REQUIRE_EQUAL(result.acceptedRoundRows.size(),1);
    BOOST_CHECK_EQUAL(result.acceptedRoundValues.size(),1);
    BOOST_CHECK_SMALL(result.acceptedRoundRisk[0]-result.reservedRisk,1e-15);
    BOOST_CHECK_LE(result.reservedRisk,1e-4);
    BOOST_CHECK_EQUAL(result.status,"COMPLETE");
}

BOOST_AUTO_TEST_CASE(r51_product_affine_recovery_preserves_original_target)
{
    const ZhangExactMatrix current{{1,-1,0}};
    const ZhangExactMatrix history{{0,1,-1}};
    const ZhangExactVector values{7};
    const ZhangExactVector target{1,0,-1};
    const auto recovery=zhangRecoverProductAffineTarget(
        current,history,values,target);
    BOOST_REQUIRE_MESSAGE(recovery.valid,recovery.reason);
    BOOST_CHECK(recovery.posteriorRow==ZhangExactVector({1}));
    BOOST_CHECK_EQUAL(recovery.offset,7);
    BOOST_CHECK(zhangVerifyProductAffineRecovery(
        current,history,values,target,recovery));
    auto corrupted=recovery;
    corrupted.offset=8;
    BOOST_CHECK(!zhangVerifyProductAffineRecovery(
        current,history,values,target,corrupted));
}

BOOST_AUTO_TEST_CASE(r51_product_affine_recovery_checks_integer_domain_and_rational_history)
{
    const auto rational=zhangRecoverProductAffineTarget(
        {},{{2}},{2},{1});
    BOOST_REQUIRE_MESSAGE(rational.valid,rational.reason);
    BOOST_CHECK_EQUAL(rational.offset,1);
    BOOST_CHECK_EQUAL(rational.historyDenominator,2);
    BOOST_CHECK_EQUAL(rational.historyNumerators[0],1);
    BOOST_CHECK(zhangVerifyProductAffineRecovery(
        {},{{2}},{2},{1},rational));
    const auto twoPivots=zhangRecoverProductAffineTarget(
        {},{{2,0},{0,3}},{4,6},{1,1});
    BOOST_REQUIRE_MESSAGE(twoPivots.valid,twoPivots.reason);
    BOOST_CHECK_EQUAL(twoPivots.offset,4);
    BOOST_CHECK(zhangVerifyProductAffineRecovery(
        {},{{2,0},{0,3}},{4,6},{1,1},twoPivots));
    const auto redundantCurrent=zhangRecoverProductAffineTarget(
        {{1,0}},{{2,0},{0,3}},{4,6},{1,1});
    BOOST_REQUIRE_MESSAGE(redundantCurrent.valid,redundantCurrent.reason);
    BOOST_CHECK(redundantCurrent.posteriorRow==ZhangExactVector({0}));
    BOOST_CHECK_EQUAL(redundantCurrent.offset,4);
    BOOST_CHECK(!zhangRecoverProductAffineTarget({},{{2}},{1},{1}).valid);
    BOOST_CHECK(!zhangRecoverProductAffineTarget({{2}},{},{},{1}).valid);
    BOOST_CHECK(!zhangRecoverProductAffineTarget(
        {{1,0,0}},{{0,1,0}},{7},{0,0,1}).valid);
}

BOOST_AUTO_TEST_CASE(r51_alternative_support_requires_a_constant_physical_bridge)
{
    // The alternative support observes the same satellite target only after
    // the old integer relation supplies a constant gauge offset.  A live
    // but unfixed direction is not allowed to masquerade as that offset.
    const ZhangExactVector desired{1,0,-1,0};
    const ZhangExactVector alternative{1,-1,0,0};
    ZhangExactVector difference(desired.size());
    for (int i=0;i<difference.size();++i)
        difference[i]=desired[i]-alternative[i];
    const ZhangExactMatrix history{{0,1,-1,0}};
    const auto bridge=zhangRecoverProductAffineTarget(
        {},history,{7},difference);
    BOOST_REQUIRE_MESSAGE(bridge.valid,bridge.reason);
    BOOST_CHECK(bridge.posteriorRow.empty());
    BOOST_CHECK_EQUAL(bridge.offset,7);
    BOOST_CHECK(zhangVerifyProductAffineRecovery(
        {},history,{7},difference,bridge));
    BOOST_CHECK(!zhangRecoverProductAffineTarget(
        {},{}, {},difference).valid);
    auto crossComponent=difference;
    crossComponent[3]=1;
    BOOST_CHECK(!zhangRecoverProductAffineTarget(
        {},history,{7},crossComponent).valid);
}
