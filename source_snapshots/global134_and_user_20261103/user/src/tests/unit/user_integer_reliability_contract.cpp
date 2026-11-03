#include "common/zhangUserArPolicy.hpp"
#include "common/zhangUserArExperiment.hpp"
#include "common/zhangUserIntegerContinuityAudit.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main()
{
    int checks = 0;
    auto check = [&](bool value, const char* reason)
    { ++checks; if (!value) throw std::runtime_error(reason); };
    try
    {
        unsetenv("ZHANG_USER_AR_RELIABILITY_PROFILE");
        check(!zhangUserArReliabilityPolicy().safe, "legacy profile default");
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "SAFE_RESEARCH", 1);
        const auto safe = zhangUserArReliabilityPolicy();
        check(safe.enforceBootstrapSuccess && safe.minimumArcSeconds == 120 &&
            safe.minimumArcSamples == 5 && safe.minimumParRank == 2 &&
            safe.maximumParTrials == 4 && safe.confirmationEpochs == 3,
            "safe policy defaults explicit");
        setenv("ZHANG_USER_AR_MIN_ARC_SECONDS", "0", 1);
        bool refused = false;
        try { zhangUserArReliabilityPolicy(); } catch (const std::exception&) { refused = true; }
        check(refused, "safe profile cannot secretly weaken maturity");
        unsetenv("ZHANG_USER_AR_MIN_ARC_SECONDS");
        unsetenv("ZHANG_USER_QC");
        refused = false;
        try { zhangUserSharedSearchRequested(); } catch (const std::exception&) { refused = true; }
        check(refused, "safe application cannot silently retain legacy observation QC");
        setenv("ZHANG_USER_QC", "CAUSAL_TRACKING", 1);
        unsetenv("ZHANG_USER_AR_SEARCH");
        check(zhangUserSharedSearchRequested(), "safe application shared search selected");
        setenv("ZHANG_USER_AR_SEARCH", "LEGACY_EPOCH", 1);
        refused = false;
        try { zhangUserSharedSearchRequested(); } catch (const std::exception&) { refused = true; }
        check(refused, "safe application cannot bypass new acceptance layer");
        unsetenv("ZHANG_USER_AR_SEARCH");
        setenv("ZHANG_USER_AR_FEEDBACK", "QUOTIENT_HOLD", 1);
        setenv("ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH", "1", 1);
        refused = false;
        try { zhangUserQuotientHoldRequested(); } catch (const std::exception&) { refused = true; }
        check(refused, "safe research never deploys HOLD");
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "LEGACY_RATIO", 1);
        check(zhangUserQuotientHoldRequested(), "legacy HOLD requires explicit research opt-in");
        unsetenv("ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH");
        refused = false;
        try { zhangUserQuotientHoldRequested(); } catch (const std::exception&) { refused = true; }
        check(refused, "default HOLD deployment refused");
        unsetenv("ZHANG_USER_AR_FEEDBACK");

        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "GONSIS_UDUC", 1);
        const auto gonsis = zhangUserArReliabilityPolicy();
        check(gonsis.gonsis && !gonsis.safe && !gonsis.enforceBootstrapSuccess &&
            gonsis.minimumArcSeconds == 0 && gonsis.minimumArcSamples == 1 &&
            gonsis.confirmationEpochs == 1 && gonsis.minimumParRank == 4 &&
            gonsis.maximumConditionalSigmaCycles == 0,
            "GONSIS UC strategy is separate not relabeled SAFE");
        setenv("ZHANG_USER_AR_SEARCH", "GONSIS_UDUC", 1);
        check(zhangUserGonsisSearchRequested() && !zhangUserSharedSearchRequested(),
            "GONSIS explicit strategy dispatch");
        setenv("ZHANG_USER_AR_SEARCH", "SHARED_WL_L1_PAR", 1);
        refused = false;
        try { zhangUserGonsisSearchRequested(); } catch (const std::exception&) { refused = true; }
        check(refused, "mixed GONSIS and shared strategy refused");
        setenv("ZHANG_USER_AR_SEARCH", "GONSIS_UDUC", 1);
        setenv("ZHANG_USER_AR_FEEDBACK", "QUOTIENT_HOLD", 1);
        setenv("ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH", "1", 1);
        refused = false;
        try { zhangUserQuotientHoldRequested(); } catch (const std::exception&) { refused = true; }
        check(refused, "GONSIS does not inherit historical integer HOLD");
        unsetenv("ZHANG_USER_AR_FEEDBACK");
        unsetenv("ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH");
        setenv("ZHANG_USER_AR_CONFIRM_EPOCHS", "3", 1);
        refused = false;
        try { zhangUserArReliabilityPolicy(); } catch (const std::exception&) { refused = true; }
        check(refused, "stale SAFE confirmation environment refused for GONSIS");
        unsetenv("ZHANG_USER_AR_CONFIRM_EPOCHS");
        unsetenv("ZHANG_USER_AR_SEARCH");
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "LEGACY_RATIO", 1);

        ZhangUserArExperimentPolicy maturity;
        maturity.minimumArcSeconds = 120; maturity.minimumArcSamples = 5;
        ZhangUserArSeparatedArcTracker ages(maturity);
        ZhangUserArExperimentKey key{"run", "site", "G13"};
        for (int t = 0; t <= 120; t += 30)
        {
            const auto a = ages.update(key, t, "physical-1", true, "product-1", true, true);
            check(a.effective.eligible == (t == 120), "young arc remains candidate-only rejected");
        }
        const auto suspect = ages.update(key, 150, "physical-1", true, "product-1", true, false);
        check(suspect.physical.age == 150 && suspect.product.age == 150 && !suspect.effective.valid,
            "suspect quarantine does not invent physical or product breaks");
        const auto resumed = ages.update(key, 180, "physical-1", true, "product-1", true, true);
        check(resumed.physical.age == 180 && resumed.effective.age == 0,
            "effective evidence restarts while physical arc retained");
        const auto product = ages.update(key, 210, "physical-1", true, "product-2", true, true);
        check(product.physical.age == 210 && product.product.age == 0 && product.effective.age == 0,
            "product identity change never relabels a physical slip");
        const auto physical = ages.update(key, 240, "physical-2", true, "product-2", true, true);
        check(physical.physical.age == 0 && physical.product.age == 30 && physical.effective.age == 0,
            "true arc break resets physical maturity only relevant domain");
        ZhangUserArSeparatedArcTracker numericResetAges(maturity);
        numericResetAges.update(key, 0, "tracking-T1", true, "product-P1", true, true, "state-G1-T1", true);
        const auto beforeReset = numericResetAges.update(key, 30, "tracking-T1", true,
            "product-P1", true, true, "state-G1-T1", true);
        check(beforeReset.physical.age == 30 && beforeReset.state.age == 30,
            "physical and state ages initially match independent identities");
        const auto numericReset = numericResetAges.update(key, 60, "tracking-T1", true,
            "product-P1", true, true, "state-G2-T1", true);
        check(numericReset.physical.age == 60 && numericReset.state.age == 0 &&
            numericReset.product.age == 60 && numericReset.effective.age == 0,
            "numerical KF replacement does not become a physical tracking break");
        const auto missingState = numericResetAges.update(key, 90, "tracking-T1", true,
            "product-P1", true, false, "state-G2-T1", false);
        check(missingState.physical.age == 90 && !missingState.state.valid && !missingState.effective.valid,
            "missing numeric state preserves observed physical tracking age");
        const auto trackingBreak = numericResetAges.update(key, 120, "tracking-T2", true,
            "product-P1", true, true, "state-G3-T2", true);
        check(trackingBreak.physical.age == 0 && trackingBreak.state.age == 0,
            "actual tracking sequence change resets physical and state ages");

        ZhangUserIntegerContinuityAudit audit;
        const std::vector<std::string> columns{"A/1", "B/1", "C/1", "D/1"};
        const std::vector<std::string> owners{"A", "B", "C", "D"};
        // Independent known truth: B-A=2, C-A=5, D-C=4. New reference B
        // yields A-B=-2, C-B=3, D-C=4; no estimator or old candidate is applied.
        const ZhangExactMatrix refA{{-1,1,0,0},{-1,0,1,0},{0,0,-1,1}};
        const ZhangExactVector valA{2,5,4};
        auto first = audit.observe("group", 0, columns, refA, valA, 3, true, owners);
        check(first.valid && !first.confirmationEligible && first.confirmationCount == 1,
            "first candidate cannot satisfy three confirmations");
        auto duplicate = audit.observe("group", 0, columns, refA, valA, 3, true, owners);
        check(duplicate.confirmationCount == 1, "duplicate same epoch is not new evidence");
        const ZhangExactMatrix refB{{1,-1,0,0},{0,-1,1,0},{0,0,-1,1}};
        auto second = audit.observe("group", 30, columns, refB, {-2,3,4}, 3, true, owners);
        check(!second.contradiction && second.confirmationCount == 2,
            "exact reference change does not reset confirmation");
        auto third = audit.observe("group", 60, {"D/1","C/1","B/1","A/1"},
            {{0,0,-1,1},{0,1,-1,0},{1,-1,0,0}}, {-2,3,4}, 3, true, {"D","C","B","A"});
        check(third.confirmationEligible && !third.contradiction,
            "physical reindex insertion order irrelevant");
        auto changed = audit.observe("group", 90, columns,
            {{-1,1,0,0},{-1,0,1,0},{0,0,-1,1}}, {3,5,4}, 3, true, owners);
        check(changed.contradiction && changed.isolatedColumns > 0,
            "same physical SD integer change isolated rather than trusting history");
        check(changed.retainedRank >= 1, "unrelated C-D relation preserved through exact restriction");
        for (const auto& row : changed.output.rows)
            check(row[0] == 0 && row[1] == 0, "quarantined physical variables cannot publish equations");
        auto missing = audit.observe("group", 120, {"C/1","D/1"}, {{-1,1}}, {4}, 3, true, {"C","D"});
        check(missing.isolatedColumns == 0 && missing.isolatedOwners == 2,
            "temporarily absent isolated arcs remain remembered by stable owner");
        auto still = audit.observe("group", 150, columns, refA, valA, 3, true, owners);
        check(still.isolatedColumns > 0, "next apparent old candidate does not decide truth");
        check(still.isolatedOwners == 2, "same identity reappearance does not release quarantine");
        auto newArcs = audit.observe("group", 180, {"A/2","B/2","C/1","D/1"},
            refA, {3,5,4}, 3, true, owners);
        check(newArcs.isolatedColumns == 0 && !newArcs.contradiction,
            "actual arc identity replacement retires quarantine no false inheritance");
        auto productChange = audit.observe("group", 210,
            {"A/2/p2","B/2/p2","C/1/p2","D/1/p2"}, refA, {9,4,3}, 3, true, owners);
        check(!productChange.contradiction && productChange.confirmationCount == 1,
            "unproved product datum change resets evidence no raw integer comparison");
        audit.observe("support", 0, {"A/1","B/1"}, {{-1,1}}, {2}, 3, true, {"A","B"});
        auto youngUnused = audit.observe("support", 30, {"A/1","B/1","X/1"},
            {{-1,1,0}}, {2}, 3, true, {"A","B","X"});
        check(youngUnused.confirmationCount == 2 && !youngUnused.contradiction,
            "unused young physical column does not reset active relation confirmation");
        auto youngLeaves = audit.observe("support", 60, {"A/1","B/1"},
            {{-1,1}}, {2}, 3, true, {"A","B"});
        check(youngLeaves.confirmationEligible && youngLeaves.confirmationCount == 3,
            "unused column disappearance does not reset active evidence");
        audit.observe("empty-gap", 0, {"A/1","B/1"}, {{-1,1}}, {2}, 2, true, {"A","B"});
        const auto emptyGap = audit.observe("empty-gap", 30, {}, {}, {}, 2, true, {});
        check(emptyGap.confirmationCount == 0 && emptyGap.auditEvidenceRank == 1,
            "empty candidate resets confirmation but preserves independent audit evidence");
        const auto afterEmpty = audit.observe("empty-gap", 60, {"A/1","B/1"},
            {{-1,1}}, {3}, 2, true, {"A","B"});
        check(afterEmpty.contradiction && afterEmpty.isolatedColumns == 2,
            "empty/rejected epoch cannot hide same-identity integer contradiction");
        audit.observe("partial-gap", 0, columns, refA, valA, 2, true, owners);
        const auto partialGap = audit.observe("partial-gap", 30, {"C/1","D/1"},
            {{-1,1}}, {4}, 2, true, {"C","D"});
        check(partialGap.auditEvidenceRank == 3 && partialGap.auditEvidenceColumns == 4,
            "nonempty subset cannot retire temporarily absent A-B audit relation");
        const auto afterPartial = audit.observe("partial-gap", 60, columns, refA,
            {3,5,4}, 2, true, owners);
        check(afterPartial.contradiction && afterPartial.isolatedColumns >= 2,
            "partial C-D-only epoch cannot hide later A-B contradiction");
        for (int epoch = 0; epoch < 100; ++epoch)
        {
            const auto bounded = audit.observe("bounded", epoch*30,
                {"A/"+std::to_string(epoch), "B/"+std::to_string(epoch)},
                {{-1,1}}, {epoch}, 1, true, {"A","B"});
            check(bounded.auditEvidenceColumns == 2 && bounded.auditEvidenceRank == 1,
                "audit metadata rank/columns bounded by owners not historical length");
        }
        const std::vector<std::string> frequencyColumns{"A/L1/G1","A/L2/G1","B/L1/G1","B/L2/G1"};
        const std::vector<std::string> frequencyOwners{"A/L1","A/L2","B/L1","B/L2"};
        audit.observe("frequency", 0, frequencyColumns, {{0,-1,0,1}}, {2}, 1, true, frequencyOwners);
        const auto l2Conflict = audit.observe("frequency", 30, frequencyColumns,
            {{0,-1,0,1}}, {3}, 1, true, frequencyOwners);
        check(l2Conflict.contradiction && l2Conflict.isolatedColumns == 2,
            "L2-only contradiction isolates its actual frequency owners");
        const auto l1Replacement = audit.observe("frequency", 60,
            {"A/L1/G2","A/L2/G1","B/L1/G1","B/L2/G1"},
            {{0,-1,0,1}}, {3}, 1, true, frequencyOwners);
        check(l1Replacement.isolatedColumns == 2 && l1Replacement.isolatedOwners == 2,
            "L1 replacement never releases unchanged L2 contradiction quarantine");
        std::cout << "USER_INTEGER_RELIABILITY_CONTRACT checks=" << checks << " status=PASS\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    { std::cerr << "FAIL " << error.what() << '\n'; return EXIT_FAILURE; }
}
