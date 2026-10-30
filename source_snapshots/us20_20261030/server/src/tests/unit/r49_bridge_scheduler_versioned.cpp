#include "common/zhangR49BridgeScheduler.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
int checks = 0;
void require(bool condition, const char* message)
{
    ++checks;
    if (!condition) throw std::runtime_error(message);
}

// This scheduler contract intentionally has no dependency on the production
// exact backend. The caller supplies already computed exact union ranks and
// the validator verdict; it cannot hide a shared algebra error in this test.
struct Domain
{
    std::string posteriorId = "FLOAT_ROOT_A";
    std::string chartId = "CURRENT_PHYSICAL_ORDER_A";
    std::string policyId = "RATIO_2_NO_NIS";
    std::uint64_t version = 0;
};
struct Route
{
    Domain domain;
    int acceptedPayload = 0;
};
}

int main()
{
    Domain initial;
    const auto g18 = zhangR49VersionedBridgeRequest(initial, "G18_WL_L1_EXACT_IMAGE");
    ZhangR49BridgeScheduler scheduler;
    require(g18.valid(), "Complete initial version-zero identity rejected");
    require(scheduler.claimRequest(g18), "First target request rejected");
    require(!scheduler.claimRequest(g18), "Same domain and target retried");
    require(scheduler.attempts == 0, "Target preparation consumed a nonexistent ILS call");
    require(scheduler.consumeSearchCall(3), "First real ILS call rejected");
    require(scheduler.chargeRisk(0.0001), "First spent allocation rejected");

    // Simulate a failed ratio and a failed append. No domain commit occurs,
    // but the actual call and spent risk stay consumed outside the safe route.
    Route safe, invalid = safe;
    invalid.acceptedPayload = 7;
    require(!zhangR49CommitGrowingTrial(safe, invalid, 4, 5,
        [](const Route&) { return false; }), "Invalid union committed");
    require(safe.domain.version == 0 && safe.acceptedPayload == 0,
        "Rejected union mutated safe route");
    require(scheduler.attempts == 1 && scheduler.spentRisk == 0.0001,
        "Rejected trial refunded search accounting");
    require(!scheduler.claimRequest(g18), "Failed same-domain request repeated");

    Route dependent = safe;
    dependent.acceptedPayload = 8;
    require(!zhangR49CommitGrowingTrial(safe, dependent, 4, 4,
        [](const Route&) { return true; }), "Entailed row advanced domain version");
    require(safe.domain.version == 0 && safe.acceptedPayload == 0,
        "Same-rank union mutated safe route");

    Route growing = safe;
    growing.acceptedPayload = 9;
    require(zhangR49CommitGrowingTrial(safe, growing, 4, 5,
        [](const Route&) { return true; }), "Valid rank-growing union rejected");
    require(safe.domain.version == 1 && safe.acceptedPayload == 9,
        "Valid extension did not commit exactly one version");
    const auto retryG18 = zhangR49VersionedBridgeRequest(safe.domain,
        "G18_WL_L1_EXACT_IMAGE");
    require(scheduler.claimRequest(retryG18), "Changed-domain G18 retry suppressed");
    require(!scheduler.claimRequest(retryG18), "Changed-domain request repeated twice");
    require(scheduler.consumeSearchCall(3), "Changed-domain real call rejected");
    require(scheduler.chargeRisk(0.0002), "Changed-domain spent allocation rejected");
    require(scheduler.attempts == 2 && scheduler.spentRisk > 0.00029,
        "Changed-domain retry discarded previous spent budget");

    // A different target in the SAME domain is still a different request.
    require(scheduler.claimRequest(zhangR49VersionedBridgeRequest(safe.domain,
        "G18_TARGETED_INDEPENDENT_CYCLE")), "Different target incorrectly deduplicated");
    require(scheduler.consumeSearchCall(3), "Last budgeted call rejected");
    require(!scheduler.consumeSearchCall(3) && scheduler.attempts == 3,
        "Changed-domain retry exceeded the fixed global call cap");
    require(!scheduler.consumeSearchCall(0), "Nonpositive budget accepted");

    // Complete scope and optional exact H/b identity are key material.
    Domain different = initial;
    different.posteriorId = "FLOAT_ROOT_B";
    require(scheduler.claimRequest(zhangR49VersionedBridgeRequest(different,
        "G18_WL_L1_EXACT_IMAGE")), "Different posterior deduplicated");
    different = initial; different.chartId += "_B";
    require(scheduler.claimRequest(zhangR49VersionedBridgeRequest(different,
        "G18_WL_L1_EXACT_IMAGE")), "Different physical chart deduplicated");
    different = initial; different.policyId += "_B";
    require(scheduler.claimRequest(zhangR49VersionedBridgeRequest(different,
        "G18_WL_L1_EXACT_IMAGE")), "Different search policy deduplicated");
    const auto h0 = zhangR49VersionedBridgeRequest(initial, "TARGET", "FULL_H0_AND_B0");
    const auto h1 = zhangR49VersionedBridgeRequest(initial, "TARGET", "FULL_H1_AND_B1");
    require(scheduler.claimRequest(h0) && scheduler.claimRequest(h1),
        "Different complete conditioner fingerprints deduplicated");
    require(!scheduler.claimRequest(h0), "Identical complete conditioner repeated");

    Domain delimitersA = initial, delimitersB = initial;
    delimitersA.posteriorId = "a|b"; delimitersA.chartId = "c";
    delimitersB.posteriorId = "a"; delimitersB.chartId = "b|c";
    require(zhangR49VersionedBridgeRequest(delimitersA, "TARGET").key() !=
        zhangR49VersionedBridgeRequest(delimitersB, "TARGET").key(),
        "Delimiter-bearing posterior/chart names collided");
    const auto requestCount = scheduler.testedRequests.size();
    Domain incomplete = initial; incomplete.posteriorId.clear();
    require(!scheduler.claimRequest(zhangR49VersionedBridgeRequest(incomplete, "TARGET")) &&
        !scheduler.claimRequest(zhangR49VersionedBridgeRequest(initial, "")) &&
        scheduler.testedRequests.size() == requestCount,
        "Invalid request admitted or changed dedup state");

    const double consumed = scheduler.spentRisk;
    require(!scheduler.chargeRisk(-0.0001) &&
        !scheduler.chargeRisk(std::numeric_limits<double>::quiet_NaN()) &&
        !scheduler.chargeRisk(std::numeric_limits<double>::infinity()) &&
        scheduler.spentRisk == consumed,
        "Invalid allocation refunded or corrupted charged risk");

    Route stale = safe;
    stale.domain.version = 0;
    require(!zhangR49CommitGrowingTrial(safe, stale, 5, 6,
        [](const Route&) { return true; }) && safe.domain.version == 1,
        "Stale affine domain committed");
    Route changedRoot = safe;
    changedRoot.domain.posteriorId += "_OTHER";
    require(!zhangR49CommitGrowingTrial(safe, changedRoot, 5, 6,
        [](const Route&) { return true; }), "Trial changed immutable posterior scope");
    Route manuallyAdvanced = safe;
    ++manuallyAdvanced.domain.version;
    require(!zhangR49CommitGrowingTrial(safe, manuallyAdvanced, 5, 6,
        [](const Route&) { return true; }), "Caller manufactured a new domain version");
    Route exhausted = safe;
    exhausted.domain.version = std::numeric_limits<std::uint64_t>::max();
    require(!zhangR49CommitGrowingTrial(exhausted, exhausted, 5, 6,
        [](const Route&) { return true; }), "Domain version wrapped around");
    require(scheduler.attempts == 3 && scheduler.spentRisk == consumed,
        "Rejected domain operation reset call/risk accounting");

    // Keep the original generic atomic-commit API compatible for existing
    // non-versioned contracts. It must still leave a rejected value unchanged.
    int atomic = 1;
    require(!zhangR49CommitTrial(atomic, 2, [](int) { return false; }) && atomic == 1,
        "Legacy rejected trial lost atomicity");
    require(zhangR49CommitTrial(atomic, 3, [](int) { return true; }) && atomic == 3,
        "Legacy valid trial rejected");
    std::cout << "R49 versioned bridge scheduler: " << checks
        << " source contract checks; changed-domain retry, exact-rank progression,"
           " nonrefundable call/risk budgets, immutable scope\n";
}
