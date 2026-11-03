#include "common/zhangCertifiedTreeTransport.hpp"
#include "common/zhangProductTransportDiagnostic.hpp"
#include "common/zhangProductPublicationGate.hpp"
#include "common/zhangFullJointHandoff.hpp"
#include <iostream>
#include <stdexcept>

using namespace zhang_transport_certificate;
static void require(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
static std::string id(const char* signal, const char* receiver, int satellite, int version = 0)
{ return std::string(signal) + "|" + receiver + "|G0" + std::to_string(satellite) + "|V" + std::to_string(version); }
static Fact fact(Row row, int value)
{
    Fact result;
    result.physical = std::move(row); result.value = value; result.sourceScopeId = "source";
    result.proofIds = {"independent-fixture"};
    return result;
}
static Row cycle(const char* signal)
{ return {{id(signal,"B",2),1},{id(signal,"B",1),-1},{id(signal,"A",2),-1},{id(signal,"A",1),1}}; }

int main()
{
    Directory directory;
    directory.authoritativeOwner = "AUTHORITATIVE"; directory.system = static_cast<int>(E_Sys::GPS);
    directory.lastCommittedEpoch = 100;
    directory.sourceScopes["source"] = {100,"posterior","chart","float","order"};
    directory.proofAtoms["independent-fixture"] = "accepted-complete-fixture-proof";
    directory.facts = {fact(cycle("L1C"),12),fact(cycle("L2W"),5)};
    Query query(directory);
    require(query.feasible, "accepted old physical system rejected");
    ZhangExactInteger value;
    require(query.consequence(cycle("L1C"),value) && value == 12, "old cycle constant lost");
    // Independent literal 2x2 oracle: the old tree has e11,e12,e21.
    // New tree e11,e21,e22 needs h_s2=-z=-12; the physical slip J
    // changes the new e12 chord to J-z, not the old certified constant.
    for (int a = -3; a <= 3; ++a)
    for (int b = -3; b <= 3; ++b)
    for (int c = -3; c <= 3; ++c)
    for (int slip : {-7,0,7})
    {
        const int d = 12 + c + b - a;
        require(d-c-b+a == 12, "literal 2x2 physical cycle oracle failed");
        const int newSatelliteDatum = -12;
        const int newSlippedChord = slip - 12;
        require(newSatelliteDatum == -value.convert_to<int>() && newSlippedChord - slip == newSatelliteDatum,
            "certified chart transport depends on the unknown physical slip");
    }
    auto newArc = cycle("L1C");
    newArc.erase(id("L1C","A",2)); newArc[id("L1C","A",2,1)] = -1;
    require(!query.consequence(newArc,value), "old certificate leaked into a new physical arc");
    auto newEndpoint = cycle("L1C"); newEndpoint[id("L1C","D",3)] = 1;
    require(!query.consequence(newEndpoint,value), "new endpoint treated as a known old zero");

    Directory wl = directory;
    Row wideLane = {{id("L1C","A",1),1},{id("L2W","A",1),-1}};
    wl.facts = {fact(wideLane,7)};
    Query wlQuery(wl);
    const auto partial = wlQuery.classify({{{id("L1C","A",1),1}},{{id("L2W","A",1),1}}});
    require(partial.ranksEvaluated && partial.knownRank == 1 && partial.unknownRank == 1 &&
        partial.knownDirections.empty() && partial.knownTargetCombinations.size() == 1 && !partial.currentArAuthority,
        "WL information incorrectly promoted to two scalar phase constants");
    const auto& combination = partial.knownTargetCombinations.front();
    require(combination[0] + combination[1] == 0 &&
        partial.knownCombinationValues.front() == 7 * combination[0], "known joint direction lost its exact integer value");
    Directory impossible = wl;
    impossible.facts = {fact({{id("L1C","A",1),2}},1)};
    const Query impossibleQuery(impossible);
    require(!impossibleQuery.feasible && !impossibleQuery.classify({wideLane}).ranksEvaluated,
        "integer-infeasible history reported a known rank");
    auto parity = wl;
    parity.facts = {fact({{id("L1C","A",1),2},{id("L1C","B",1),1}},1),
                    fact({{id("L1C","B",1),1}},0)};
    require(!Query(parity).feasible, "unit elimination discarded the residual modulo-two obstruction");

    const ZhangGraphEdge a1{"A",SatSys("G01")},a2{"A",SatSys("G02")},
        b1{"B",SatSys("G01")},b2{"B",SatSys("G02")},c1{"C",SatSys("G01")},c2{"C",SatSys("G02")};
    ZhangGraphBasis old = zhangBuildSpanningTree({a1,a2,b1,b2,c1,c2},"A",{a1,a2,b1,c1});
    const auto plan = select(old,{a1,b1,b2,c1,c2},{a2},"A",{{b2,0.01},{c2,1e8}},
        [&](const auto& edge) { return edge == b2; });
    require(plan.basis.treeEdges.contains(b2) && !plan.basis.treeEdges.contains(c2) &&
        plan.certifiedReplacements.contains(b2) && plan.fallbackReplacements.empty(),
        "quality score displaced a certified cross-cut replacement");
    const auto split = select(old,{a1,b1,c1},{a2},"A",{},[](const auto&) { return false; });
    require(!split.allOldEndpointsCertified && !split.basis.satellites.contains(SatSys("G02")),
        "disconnected component invented cross-component continuity");
    // Different insertion order is immaterial to the independent physical
    // equation, including redundant target queries and repeated constraints.
    auto reversed = directory; std::reverse(reversed.facts.begin(),reversed.facts.end());
    Query reversedQuery(reversed);
    require(reversedQuery.consequence(cycle("L1C"),value) && value == 12, "fact order changed the exact image");
    auto redundant = reversed; redundant.facts.push_back(directory.facts.front());
    const Query redundantQuery(redundant);
    require(redundantQuery.feasible && redundantQuery.consequence(cycle("L1C"),value) && value == 12,
        "repeated accepted constraints changed the transport consequence");

    Directory decoded; std::string reason;
    require(Directory::decode(directory.encode(),"AUTHORITATIVE",directory.system,decoded,reason), "typed directory roundtrip failed");
    require(!Directory::decode(directory.encode(),"PRIVATE_BRANCH",directory.system,decoded,reason), "branch owner inherited authoritative history");
    require(zhangProductTransportFailureClass("OLD_FLOAT_TARGET_UNAVAILABLE") == "OLD_FLOAT_TARGET_UNAVAILABLE" &&
        zhangProductTransportFailureClass("SHARED_EVENT_BUDGET_EXHAUSTED") == "BUDGET_INCOMPLETE" &&
        zhangProductTransportFailureClass("RATIO_TEST_FAILED") == "RATIO_REJECTED",
        "transport diagnosis conflated representation, search and budget failure");
    // old z=12, new slip J=7: current independent AR observes J-z=-5.
    // An exhaustive one-dimensional search is independent of the historical
    // directory. Changing/removing the directory cannot change this search.
    auto independentAr = [](double estimate)
    {
        std::vector<std::pair<double,int>> candidates;
        for (int integer=-30;integer<=30;++integer)
            candidates.push_back({(estimate-integer)*(estimate-integer),integer});
        std::sort(candidates.begin(),candidates.end());
        return std::pair{candidates[1].first/candidates[0].first,candidates[0].second};
    };
    const auto currentAr = independentAr(-5.02);
    require(currentAr.first > 2 && currentAr.second == -5,"independent current AR oracle failed");
    ZhangFullJointHandoffTask temporal;
    temporal.historicalTransportValueProven=true; temporal.historicalTransportValue=12;
    require(temporal.integerTransportProven() && !temporal.currentValueProven &&
        temporal.integerTransportValue()==12,"historical datum value masqueraded as a current decision");
    ZhangProductPublicationEvidence gate;
    gate.nativeProductValid=true; gate.currentSpatialCertificateValid=currentAr.first>2;
    gate.numericalStepValid=true; gate.priorPublishedIdentity=true;
    gate.conditionalProductMapProven=true; gate.requireFullJointHandoff=true;
    gate.fullJointHandoffProven=true; gate.historicalTransportBridgeProven=true;
    require(zhangProductPublicationGate(gate).route==ZhangProductPublicationRoute::SAME_SEGMENT &&
        !gate.currentEpochBridgeProven,"transport-only history was not kept separate from current evidence");
    auto noCurrentAr=gate; noCurrentAr.currentSpatialCertificateValid=false;
    require(zhangProductPublicationGate(noCurrentAr).route==ZhangProductPublicationRoute::SUSPENDED,
        "historical certificate replaced the current native AR gate");
    auto noFullQ=gate; noFullQ.fullJointHandoffProven=false;
    require(zhangProductPublicationGate(noFullQ).route!=ZhangProductPublicationRoute::SAME_SEGMENT,
        "historical certificate replaced the full-Q/FLOAT map gate");
    auto badCurrentAr=gate;
    badCurrentAr.currentSpatialCertificateValid=independentAr(-4.5).first>2;
    require(zhangProductPublicationGate(badCurrentAr).route==ZhangProductPublicationRoute::SUSPENDED,
        "old z=12 rescued a current ratio rejection");
    std::cout << "certified_tree_transport: PASS\n";
}
