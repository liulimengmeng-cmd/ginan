#include <boost/test/unit_test.hpp>
#include "common/zhangProductPhysicalCycleChart.hpp"
#include "common/zhangProductIntegerLedger.hpp"
namespace {
ZhangGraphBasis r46Basis(bool pivot) {
    ZhangGraphBasis b; b.rootReceiver="A";b.receivers={"A","B"};
    b.satellites={SatSys("G01"),SatSys("G02"),SatSys("G03")};
    for(const auto& r:b.receivers) for(const auto& s:b.satellites) b.edges.insert({r,s});
    b.treeEdges={{"A",SatSys("G01")},{"A",SatSys("G02")},{"A",SatSys("G03")},
        {"B",SatSys(pivot?"G02":"G01")}};b.connected=true;b.componentCount=1;return b;
}
ZhangProductPhysicalCycleChart r46Chart(bool pivot,int version=0,bool omit=false) {
    auto b=r46Basis(pivot);std::map<ZhangGraphEdge,int> v;
    for(const auto& e:b.edges) v[e]=(e.receiver=="A" && e.satellite==SatSys("G01"))?version:0;
    ZhangProductPhysicalCycleChart chart;chart.columns=2;
    if(!omit) BOOST_REQUIRE(chart.add(0,"L1C",{"B",SatSys(pivot?"G01":"G02")},b,v));
    BOOST_REQUIRE(chart.add(1,"L1C",{"B",SatSys("G03")},b,v));return chart;
}
ProductIntegerLedgerRow r46Row(std::map<std::string,ZhangExactInteger> f,int n) {
    ProductIntegerLedgerRow r;r.system=E_Sys::GPS;r.firstObservable=E_ObsCode::L1C;
    r.secondObservable=E_ObsCode::L2W;r.productRow={1,-1};r.physicalExpansion=std::move(f);
    r.physicalExpansionExact=true;r.integerValue=n;r.canonicalProductExpansion={{"1|G01",1},{"1|G02",-1}};
    r.phaseSegmentFingerprint="G01|L1C|SEG0;G02|L1C|SEG0;";
    r.admissionFailureProbabilityBound=1e-6;return r;
}
}
BOOST_AUTO_TEST_CASE(r46_tree_pivot_expands_and_roundtrips_true_physical_functional) {
    const auto old=r46Chart(false), current=r46Chart(true);
    std::map<std::string,ZhangExactInteger> physical;
    BOOST_REQUIRE(old.expand({0,1},physical));
    BOOST_CHECK_EQUAL(physical.size(),4);
    ZhangExactVector row;BOOST_REQUIRE(current.project(physical,row));
    BOOST_CHECK_EQUAL(row[0],-1);BOOST_CHECK_EQUAL(row[1],1);
    BOOST_CHECK_EQUAL(row[0]*(-4)+row[1]*5,9);
    std::map<std::string,ZhangExactInteger> fake={{"L1C|B|G03|V0",1}};
    BOOST_CHECK(!current.project(fake,row));
}
BOOST_AUTO_TEST_CASE(r46_retired_arc_and_missing_posterior_column_fail_closed) {
    std::map<std::string,ZhangExactInteger> f;BOOST_REQUIRE(r46Chart(false).expand({0,1},f));
    ZhangExactVector row;std::string reason;
    BOOST_CHECK(!r46Chart(true,1).project(f,row,&reason));BOOST_CHECK_EQUAL(reason,"ARC_RETIRED");
    BOOST_CHECK(!r46Chart(true,0,true).project(f,row,&reason));
    BOOST_CHECK_EQUAL(reason,"MISSING_POSTERIOR_COLUMN_OR_NO_CURRENT_REPRESENTATION");
}
BOOST_AUTO_TEST_CASE(r46_distinct_physical_support_requires_joint_affine_consistency) {
    ProductIntegerLedger ledger;
    auto a=r46Row({{"N1",1}},27), b=r46Row({{"N2",1}},11);
    BOOST_REQUIRE(ledger.observe(100,{a},1).valid);
    BOOST_REQUIRE(ledger.observe(101,{b},1).valid);
    BOOST_CHECK_EQUAL(ledger.rows().size(),2);
    auto consistent=r46Row({{"N1",1},{"N2",-1}},16);
    BOOST_REQUIRE(ledger.observe(102,{consistent},1).valid);
    const auto before=ledger.rows().size();
    auto contradiction=r46Row({{"N1",1},{"N2",1}},39);
    const auto rejected=ledger.observe(103,{contradiction},1);
    BOOST_CHECK(!rejected.valid);
    BOOST_CHECK_EQUAL(rejected.failureReason,"PRODUCT_LEDGER_TRUE_PHYSICAL_AFFINE_CONFLICT");
    BOOST_CHECK_EQUAL(ledger.rows().size(),before);
}
BOOST_AUTO_TEST_CASE(r46_same_true_physical_row_confirms_across_chart_labels) {
    ProductIntegerLedger ledger;auto a=r46Row({{"N1",1},{"N2",-1}},16);
    BOOST_REQUIRE(ledger.observe(100,{a},2).valid);
    a.backendBasisGeneration=44;a.canonicalProductExpansion={{"1|G04",1},{"1|G05",-1}};
    BOOST_REQUIRE(ledger.observe(101,{a},2).valid);
    BOOST_CHECK_EQUAL(ledger.rows().size(),1);BOOST_CHECK(ledger.rows()[0].certified);
    a.integerValue=17;auto rejected=ledger.observe(102,{a},2);
    BOOST_CHECK(!rejected.valid);BOOST_CHECK(rejected.conflictPhysicalRowsEqual);
}

#include "common/zhangSequentialQuotientShadow.hpp"
namespace {
ZhangSequentialShadowProposal r46ShadowIdentity(const Eigen::VectorXd& mean,double allocation) {
    ZhangSequentialShadowProposal p;p.valid=true;p.failureProbability=allocation/2;
    p.rows=ZhangExactMatrix(mean.size(),ZhangExactVector(mean.size()));
    for(int i=0;i<mean.size();++i) {p.rows[i][i]=1;p.values.push_back(std::llround(mean(i)));}
    return p;
}
}
BOOST_AUTO_TEST_CASE(r46_sequential_shadow_removes_fixed_mixed_dimensions_and_keeps_input) {
    const ZhangExactMatrix baseline={{1,-1,0,0}};
    const ZhangExactVector rhs={-1};
    Eigen::VectorXd mean(4);mean<<1,2,3,4;
    std::vector<int> dimensions;
    const auto out=zhangSequentialQuotientShadow(mean,0.01*Eigen::MatrixXd::Identity(4,4),
        baseline,rhs,1e-4,1e-6,
        [&](const Eigen::VectorXd& m,const Eigen::MatrixXd& c,double allocation,bool) {
            dimensions.push_back(m.size());BOOST_CHECK(c.diagonal().minCoeff()>0);
            return r46ShadowIdentity(m,allocation);
        },[](const std::string&){},4,2,1);
    BOOST_CHECK_EQUAL(out.rows.size(),4);BOOST_CHECK_GE(out.rounds,2);
    BOOST_REQUIRE(dimensions.size()>=2);BOOST_CHECK_EQUAL(dimensions[0],2);
    BOOST_CHECK_EQUAL(dimensions[1],1);BOOST_CHECK_LE(out.reservedRisk,1e-4);
    BOOST_CHECK_EQUAL(baseline.size(),1);BOOST_CHECK_EQUAL(rhs[0],-1);
}
BOOST_AUTO_TEST_CASE(r46_overlap_conflict_requires_joint_merged_resolve) {
    int calls=0;
    const auto out=zhangSequentialQuotientShadow(Eigen::VectorXd::Zero(4),
        Eigen::MatrixXd::Identity(4,4),{},{},1e-4,1e-6,
        [&](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double allocation,bool merge) {
            auto p=r46ShadowIdentity(m,allocation);++calls;
            if(calls==2 && !merge) p.values[0]=1;
            return p;
        },[](const std::string&){},4,2,4);
    BOOST_CHECK_EQUAL(out.overlapConflicts,1);BOOST_CHECK_EQUAL(out.mergeAttempts,1);
    BOOST_CHECK_EQUAL(out.mergeAccepted,1);BOOST_CHECK_EQUAL(out.rows.size(),4);
    BOOST_CHECK_GT(out.overlapCommonRankTotal,0);
    for(const auto& v:out.values) BOOST_CHECK_EQUAL(v,0);
    BOOST_CHECK_LE(out.reservedRisk,1e-4);
}
BOOST_AUTO_TEST_CASE(r46_failed_overlap_merge_rolls_back_round_and_charges_discarded_attempts) {
    int calls=0;
    const auto out=zhangSequentialQuotientShadow(Eigen::VectorXd::Zero(4),
        Eigen::MatrixXd::Identity(4,4),{},{},1e-4,1e-6,
        [&](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double allocation,bool merge) {
            auto p=r46ShadowIdentity(m,allocation);++calls;
            if(calls==2) p.values[0]=1;
            if(merge) p.valid=false;
            return p;
        },[](const std::string&){},4,2,4);
    BOOST_CHECK(out.rows.empty());BOOST_CHECK_EQUAL(out.attempts,3);
    BOOST_CHECK_EQUAL(out.mergeAccepted,0);BOOST_CHECK_GT(out.reservedRisk,0);
    BOOST_CHECK_EQUAL(out.status,"OVERLAP_CONFLICT_MERGE_REJECTED");
}
BOOST_AUTO_TEST_CASE(r46_shadow_rejects_overbudget_and_joint_nis_outlier) {
    auto run=[&](bool overbudget) {
        return zhangSequentialQuotientShadow(Eigen::VectorXd::Zero(2),
            0.0001*Eigen::MatrixXd::Identity(2,2),{},{},1e-4,1e-6,
            [&](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double allocation,bool) {
                auto p=r46ShadowIdentity(m,allocation);
                if(overbudget) p.failureProbability=allocation*2;
                else p.values[0]=100;
                return p;
            },[](const std::string&){},2,2,1);
    };
    const auto budget=run(true),nis=run(false);
    BOOST_CHECK(budget.rows.empty());BOOST_CHECK(nis.rows.empty());
    BOOST_CHECK_EQUAL(budget.status,"NO_RANK_GAIN");
    BOOST_CHECK_EQUAL(nis.status,"WHOLE_LATTICE_NIS_REJECTED");
}

BOOST_AUTO_TEST_CASE(r46_true_physical_sign_is_independent_of_canonical_chart_orientation) {
    ProductIntegerLedger ledger;
    auto a=r46Row({{"N1",1},{"N2",-1}},16);
    auto b=r46Row({{"N1",-1},{"N2",1}},-16);
    // Same canonical sign but opposite physical orientation: the true
    // physical row, rather than a transient named view, owns normalization.
    BOOST_REQUIRE(ledger.observe(100,{a},2).valid);
    BOOST_REQUIRE(ledger.observe(130,{b},2).valid);
    BOOST_CHECK_EQUAL(ledger.rows().size(),1);
    BOOST_CHECK_EQUAL(ledger.rows()[0].confirmationEpochs,2);
}
BOOST_AUTO_TEST_CASE(r46_joint_physical_integer_feasibility_rejects_parity_conflict) {
    ProductIntegerLedger ledger;
    auto a=r46Row({{"N1",1},{"N2",1}},0);
    auto b=r46Row({{"N1",1},{"N2",-1}},1);
    BOOST_REQUIRE(ledger.observe(100,{a},1).valid);
    const auto result=ledger.observe(130,{b},1);
    BOOST_CHECK(!result.valid);BOOST_CHECK_EQUAL(ledger.rows().size(),1);
    BOOST_CHECK_EQUAL(result.failureReason,"PRODUCT_LEDGER_TRUE_PHYSICAL_AFFINE_CONFLICT");
}

BOOST_AUTO_TEST_CASE(r46_partial_merged_answer_cannot_hide_disputed_overlap) {
    int calls=0;
    const auto out=zhangSequentialQuotientShadow(Eigen::VectorXd::Zero(4),
        Eigen::MatrixXd::Identity(4,4),{},{},1e-4,1e-6,
        [&](const Eigen::VectorXd& m,const Eigen::MatrixXd&,double allocation,bool merge) {
            auto p=r46ShadowIdentity(m,allocation);++calls;
            if(calls==2) p.values[0]=1;
            if(merge) {p.rows.resize(1);p.values.resize(1);}
            return p;
        },[](const std::string&){},4,2,4);
    BOOST_CHECK(out.rows.empty());BOOST_CHECK_EQUAL(out.mergeAttempts,1);
    BOOST_CHECK_EQUAL(out.mergeAccepted,0);
    BOOST_CHECK_EQUAL(out.status,"OVERLAP_CONFLICT_MERGE_REJECTED");
}
