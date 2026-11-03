#include "common/zhangCurrentBridgeWorkIdentity.hpp"
#include "common/zhangHandoffIncompleteTicket.hpp"

#include <iostream>
#include <set>
#include <stdexcept>

static void require(bool ok,const char* message)
{ if(!ok)throw std::runtime_error(message); }

int main()
{
    std::map<ZhangCurrentEpochEvidenceKey,ZhangCurrentEpochEvidence> registry;
    auto makeEvidence=[](E_Sys sys) {
        ZhangCurrentEpochEvidence e;
        e.valid=true;e.runtimeId="runtime";e.system=sys;e.epochGpstSeconds=1405216800;
        e.posteriorId="posterior";e.chartId="chart";e.floatStateRoot="root";
        e.floatPosteriorPsdStateOrderId="order";e.constraintGeneration=7;
        KFKey key;key.type=KF::AMBIGUITY;key.num=1;key.str="receiver";
        e.numericConstraints.push_back({{{key,1}},3});
        return e;
    };
    registry[{"runtime",E_Sys::GPS}]=makeEvidence(E_Sys::GPS);
    registry[{"runtime",E_Sys::GAL}]=makeEvidence(E_Sys::GAL);
    const auto key=[&] {
        return zhangCurrentBridgeConditioningKey(registry,"runtime",1405216800,
            "posterior","chart","root","order");
    };
    const auto original=key();
    require(original.size()==2 && original==key(),"work identity unstable");
    auto& other=registry.at({"runtime",E_Sys::GAL});
    ++other.constraintGeneration;
    require(key()!=original,"other constellation mutation reused stale pair");
    using PairKey=std::tuple<E_Sys,int,int,std::uint64_t,ZhangCurrentBridgeConditioningKey>;
    const PairKey initialPair{E_Sys::GPS,1,2,7,original};
    const auto finalPair=zhangCurrentBridgeRekey(initialPair,7,key());
    std::map<PairKey,ZhangHandoffPairOutcome> outcomes;
    outcomes.emplace(initialPair,ZhangHandoffPairOutcome{
        "CURRENT_EVENT_OVERLAP_WALL_TIMEOUT","initial-representative",false});
    outcomes.emplace(finalPair,ZhangHandoffPairOutcome{
        "CURRENT_EVENT_OVERLAP_RATIO_REJECTED","retry-representative",false,2});
    require(finalPair!=initialPair && outcomes.at(finalPair).reason=="CURRENT_EVENT_OVERLAP_RATIO_REJECTED",
        "alias retained old timeout instead of final generation outcome");
    require(!outcomes.contains(zhangCurrentBridgeRekey(initialPair,8,key())),
        "unsearched generation falsely inherited old terminal result");
    --other.constraintGeneration;
    other.numericConstraints[0].value=4;
    require(key()!=original,"same-size numerical RHS change reused stale pair");
    other.numericConstraints[0].value=3;
    other.numericConstraints[0].rowCycles.begin()->second=2;
    require(key()!=original,"same-size numerical row change reused stale pair");
    other.numericConstraints[0].rowCycles.begin()->second=1;
    registry[{"unrelated",E_Sys::GPS}]=makeEvidence(E_Sys::GPS);
    require(key()==original,"unrelated runtime contaminated key");
    other.floatStateRoot="different-root";
    require(key().size()==1,"unmatched FLOAT root entered conditioner key");

    std::set<int> tickets;
    std::string reason;
    ZhangExactWorkBudget budget;
    auto insertAndExit=[&](const char* finalReason,bool exhausted,const char* budgetReason) {
        tickets.insert(1);reason=finalReason;budget.exhausted=exhausted;budget.reason=budgetReason;
        ZhangHandoffIncompleteTicket guard([&] { tickets.erase(1); },reason,budget);
    };
    insertAndExit("CURRENT_EVENT_OVERLAP_WALL_TIMEOUT",false,"NONE");
    require(tickets.empty(),"reported wall timeout retained candidate ticket");
    insertAndExit("CURRENT_EVENT_OVERLAP_RATIO_REJECTED",false,"NONE");
    require(tickets.contains(1),"terminal ratio decision was retried");
    tickets.clear();
    insertAndExit("CURRENT_EVENT_OVERLAP_EXACT_OPERATION_LIMIT",true,"EXACT_OPERATION_LIMIT");
    require(tickets.contains(1),"operation limit was incorrectly treated as renewable time");
    tickets.clear();
    try {
        tickets.insert(1);reason="NOT_EVALUATED";budget.exhausted=true;budget.reason="WALL_TIMEOUT";
        ZhangHandoffIncompleteTicket guard([&] { tickets.erase(1); },reason,budget);
        throw ZhangExactWorkBudgetExceeded{};
    } catch(const ZhangExactWorkBudgetExceeded&) {}
    require(tickets.empty(),"cooperative unwinding retained incomplete ticket");
    std::cout<<"complete conditioner work identity and timeout ticket contracts passed\n";
}
