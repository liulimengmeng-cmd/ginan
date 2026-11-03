#include "common/zhangCurrentEpochEvidence.hpp"

#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

using Row = ZhangCurrentEpochPhysicalRow;
using Rational = boost::multiprecision::cpp_rational;

static void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}

// Independent Gaussian elimination with synchronized RHS. The reference does
// not call production HNF, integer quotient/kernel or membership routines.
// Random systems below are known integer-feasible by their construction.
static bool referenceConsequence(std::vector<std::vector<Rational>> h,
    std::vector<Rational> rhs, std::vector<Rational> target, Rational& value)
{
    std::size_t rank = 0;
    value = 0;
    for (std::size_t column = 0; column < target.size() && rank < h.size(); ++column)
    {
        std::size_t pivot = rank;
        while (pivot < h.size() && h[pivot][column] == 0) ++pivot;
        if (pivot == h.size()) continue;
        std::swap(h[rank], h[pivot]);
        std::swap(rhs[rank], rhs[pivot]);
        const Rational divisor = h[rank][column];
        for (auto& coefficient : h[rank]) coefficient /= divisor;
        rhs[rank] /= divisor;
        for (std::size_t row = rank + 1; row < h.size(); ++row)
        {
            const Rational factor = h[row][column];
            for (std::size_t c = 0; c < target.size(); ++c)
                h[row][c] -= factor * h[rank][c];
            rhs[row] -= factor * rhs[rank];
        }
        const Rational factor = target[column];
        for (std::size_t c = 0; c < target.size(); ++c)
            target[c] -= factor * h[rank][c];
        value += factor * rhs[rank];
        ++rank;
    }
    return std::all_of(target.begin(), target.end(),
        [](const Rational& coefficient) { return coefficient == 0; });
}

static ZhangProductPhysicalCycleChart chartFor(const std::vector<std::string>& ids)
{
    ZhangProductPhysicalCycleChart chart;
    chart.columns = static_cast<int>(ids.size());
    for (std::size_t i = 0; i < ids.size(); ++i)
    {
        chart.chordIdentities[i] = ids[i];
        chart.expansions[i][ids[i]] = 1;
        chart.activePhysicalIdentities.insert(ids[i]);
    }
    return chart;
}

static ZhangDecisionProofPtr proofFor(const std::string& id)
{
    auto proof = std::make_shared<ZhangIntegerDecisionProof>();
    proof->id = id;
    proof->originalStatement = "current fixture integer constraint";
    proof->observationProvenance = "runtime|PPP_COMMITTED|E1405216800";
    proof->conditionalFailureBound = 1e-6;
    return proof;
}

int main()
{
    constexpr long int epoch = 1405216800;
    const std::string arc = "L1C|R|G02|V1";
    const auto chart = chartFor({arc});
    const auto proof = proofFor("one");
    ZhangCurrentEpochEvidence evidence;
    require(evidence.bind("runtime", E_Sys::GPS, epoch, "posterior", "chart", chart, "float-root"),
        "current fixture bind failed");
    require(evidence.add({{arc, 2}}, 6, {proof}), "current fixture admission failed");
    ZhangProductIntegerTransition target;
    target.physicalEdges = {{"R", SatSys(E_Sys::GPS, 2)}};
    target.physicalArcVersions = {1};
    target.coefficients = {1};
    target.valid = true;
    auto result = evidence.query(E_ObsCode::L1C, epoch, "posterior", "chart", target, false);
    require(result.proven && result.value == 3, "2a=6 did not entail a=3");
    require(result.proofs.empty() && result.proofsView().size() == 1,
        "shared proof query copied roots or lost proof ancestry");
    const auto domain = evidence.compiledCurrentDomain();
    require(!evidence.floatPosteriorPsdValidated && evidence.floatPosteriorPsdRootId.empty() &&
        evidence.floatPosteriorPsdStateOrderId.empty(),
        "whole-FLOAT PSD certificate fabricated by current evidence binding");
    require(domain->identity.epochGpstSeconds == epoch &&
        domain->identity.rootId == "float-root" &&
        !domain->identity.physicalColumnsHash.empty(), "cache identity incomplete");
    for (int i = 0; i < 25; ++i)
    {
        const auto repeated = evidence.query(E_ObsCode::L1C, epoch, "posterior", "chart", target, false);
        require(repeated.proven && repeated.sharedProofs == result.sharedProofs,
            "current domain/proof roots not reused");
    }
    require(evidence.domainCache->builds == 1 && evidence.domainCache->hits >= 25,
        "one immutable current domain repeatedly compiled");
    const auto generation = evidence.constraintGeneration;
    require(evidence.add({{arc, 2}}, 6, {proof}) && evidence.constraintGeneration == generation,
        "duplicate current row unnecessarily invalidated domain");
    require(!evidence.query(E_ObsCode::L1C, epoch + 30, "posterior", "chart", target).proven,
        "cross-epoch cache authorized bridge");
    require(!evidence.query(E_ObsCode::L1C, epoch, "other-root", "chart", target).proven,
        "different posterior cache authorized bridge");
    require(!evidence.query(E_ObsCode::L1C, epoch, "posterior", "other-chart", target).proven,
        "different chart cache authorized bridge");
    const auto membership = zhangIntegerRowLatticeContains(ZhangExactMatrix{{2}}, ZhangExactVector{1});
    require(!membership.contained, "logical consequence confused with integer row membership");

    // A transactional copy can add a constraint without mutating its parent's
    // cached domain. The newly contradictory system must fail even for a
    // direct exact-row query, not just for rationally derived targets.
    auto proposed = evidence;
    require(proposed.add({{arc, 3}}, 10, {proof}), "contradictory admission fixture failed");
    require(proposed.constraintGeneration > generation &&
        proposed.compiledCurrentDomain() != domain, "new integer row kept stale domain");
    require(!proposed.query(E_ObsCode::L1C, epoch, "posterior", "chart", target).proven,
        "infeasible integer system authorized consequence");
    require(evidence.compiledCurrentDomain() == domain && domain->feasible,
        "proposed mutation changed parent domain");
    ZhangCurrentEpochEvidence parity;
    require(parity.bind("runtime", E_Sys::GPS, epoch, "posterior", "chart", chart, "float-root"),
        "parity fixture bind failed");
    require(parity.add({{arc, 2}}, 3, {proof}), "parity fixture admission failed");
    target.coefficients = {2};
    require(!parity.query(E_ObsCode::L1C, epoch, "posterior", "chart", target).proven,
        "direct 2a=3 row bypassed integer feasibility");
    target.coefficients.clear();
    target.physicalEdges.clear();
    target.physicalArcVersions.clear();
    require(!evidence.query(E_ObsCode::L1C, epoch, "posterior", "chart", target).proven,
        "structural zero manufactured integer evidence");
    auto changedChart = chart;
    changedChart.activePhysicalIdentities.insert("L1C|other|G03|V1");
    require(!evidence.bind("runtime", E_Sys::GPS, epoch, "posterior", "chart", changedChart, "float-root"),
        "changed physical catalog rebound old cache");
    auto samePointerChangedRoot = evidence;
    samePointerChangedRoot.floatStateRoot = "different-float-root";
    require(samePointerChangedRoot.compiledCurrentDomain() != domain &&
        samePointerChangedRoot.compiledCurrentDomain()->identity.rootId == "different-float-root",
        "cache relied on shared-pointer identity instead of FLOAT root");
    auto samePointerChangedColumns = evidence;
    samePointerChangedColumns.chart.cataloguedSignals.insert("L2W");
    require(samePointerChangedColumns.compiledCurrentDomain()->identity.physicalColumnsHash !=
        domain->identity.physicalColumnsHash, "changed physical catalog kept cache content hash");

    // Unknown physical columns remain unknown even if overlap evidence admits
    // querying outside the frozen current chart.
    auto overlap = evidence;
    require(overlap.addEventOverlap({{"old-arc", 1}}, 7, {proof}, "current-event"),
        "overlap fixture admission failed");
    ZhangExactInteger value;
    require(overlap.compiledCurrentDomain()->consequence({{"old-arc", 1}}, value) && value == 7,
        "current overlap consequence missing");
    require(!overlap.compiledCurrentDomain()->consequence({{"unknown", 1}}, value),
        "unknown physical variable treated as zero");

    KFKey key;
    key.type = KF::AMBIGUITY;
    key.Sat = SatSys(E_Sys::GPS, 2);
    key.str = "R";
    require(evidence.addNumericConstraint({{key, 2}}, 6), "numeric current row rejected");
    require(evidence.addNumericConstraint({{key, 2}}, 6) && evidence.numericConstraints.size() == 1,
        "duplicate numeric constraint applied twice");
    require(!evidence.addNumericConstraint({{key, 2}}, 7), "numeric RHS conflict accepted");
    require(!evidence.addNumericConstraint({{key, 0}}, 0), "structural zero became numeric AR row");
    require(!evidence.addNumericConstraint({{key, std::numeric_limits<double>::infinity()}}, 0),
        "nonfinite numeric constraint accepted");

    // Logical evidence queries remain backward compatible, but the stronger
    // missing-arc representation operation needs an explicit admission-time
    // physical/numerical link for EVERY raw H row, including event rows.
    std::string bindingReason;
    require(!evidence.validatePhysicalNumericBindings(bindingReason),
        "unbound legacy mathematical evidence authorized numerical rewrite");
    auto bound = evidence;
    bound.floatPosteriorPsdStateOrderId = "state-order";
    require(bound.addNumericConstraint({{key, 2}}, 6, {{arc, 2}}),
        "same-root physical/numeric admission binding rejected");
    require(bound.validatePhysicalNumericBindings(bindingReason) &&
        bindingReason == "CURRENT_PHYSICAL_NUMERIC_BINDINGS_VERIFIED" &&
        bound.numericConstraints.size() == 1 && bound.physicalNumericBindings.size() == 1,
        "bound current raw row absent from numerical recipe");
    require(bound.addNumericConstraint({{key, 2}}, 6, {{arc, 2}}) &&
        bound.physicalNumericBindings.size() == 1, "duplicate binding grew current recipe");
    require(!bound.addNumericConstraint({{key, 3}}, 6, {{arc, 2}}) &&
        bound.numericConstraints.size() == 1 && bound.physicalNumericBindings.size() == 1,
        "conflicting binding mutated recipe before rejection");
    require(!bound.addNumericConstraint({{key, 2}}, 6, {{"not-admitted", 2}}),
        "binding admitted a physical equation not in current raw H");
    auto missingRecipe = bound;
    missingRecipe.numericConstraints.clear();
    require(!missingRecipe.validatePhysicalNumericBindings(bindingReason) &&
        bindingReason == "CURRENT_PHYSICAL_NUMERIC_BINDING_RECIPE_MISSING_OR_RHS_MISMATCH",
        "isolated physical proof substituted for numerical conditioning");
    auto wrongNumericRhs = bound;
    wrongNumericRhs.numericConstraints[0].value = 7;
    require(!wrongNumericRhs.validatePhysicalNumericBindings(bindingReason),
        "binding ignored numerical RHS mismatch");
    auto wrongPhysicalRhs = bound;
    wrongPhysicalRhs.values[0] = 7;
    require(!wrongPhysicalRhs.validatePhysicalNumericBindings(bindingReason) &&
        bindingReason == "CURRENT_PHYSICAL_NUMERIC_RAW_ROW_UNBOUND_OR_RHS_MISMATCH",
        "binding ignored raw physical RHS mismatch");
    for (int mismatch = 0; mismatch < 4; ++mismatch)
    {
        auto wrongScope = bound;
        if (mismatch == 0) wrongScope.epochGpstSeconds += 30;
        if (mismatch == 1) wrongScope.floatStateRoot = "different-float-root";
        if (mismatch == 2) wrongScope.floatPosteriorPsdStateOrderId = "different-state-order";
        if (mismatch == 3) wrongScope.posteriorId = "different-posterior";
        require(!wrongScope.validatePhysicalNumericBindings(bindingReason) &&
            bindingReason == "CURRENT_PHYSICAL_NUMERIC_BINDING_SCOPE_MISMATCH",
            "binding crossed epoch/FLOAT root/state order/posterior scope");
    }
    auto unboundEvent = bound;
    require(unboundEvent.addEventOverlap({{"old-arc", 1}}, 7, {proof}, "current-event"),
        "event binding fixture rejected");
    require(!unboundEvent.validatePhysicalNumericBindings(bindingReason),
        "unbound event H row authorized numerical rewrite");
    KFKey eventKey = key;
    eventKey.str = "event-auxiliary";
    require(unboundEvent.addNumericConstraint({{eventKey, 1}}, 7, {{"old-arc", 1}}) &&
        unboundEvent.validatePhysicalNumericBindings(bindingReason),
        "event raw H binding absent from complete numerical recipe");
    require(bound.validatePhysicalNumericBindings(bindingReason) &&
        evidence.physicalNumericBindings.empty(), "provisional binding changed parent evidence");
    ZhangCurrentEpochEvidence noRoot;
    require(noRoot.bind("runtime", E_Sys::GPS, epoch, "posterior", "chart", chart),
        "legacy root-free bind fixture failed");
    require(!noRoot.addNumericConstraint({{key, 1}}, 0), "numeric constraint lacked FLOAT root");

    std::mt19937 random(719);
    std::size_t comparisons = 0;
    for (int trial = 0; trial < 80; ++trial)
    {
        const std::vector<std::string> ids{"a", "b", "c"};
        ZhangCurrentEpochEvidence current;
        require(current.bind("runtime", E_Sys::GPS, epoch, "posterior", "chart",
            chartFor(ids), "float-root"), "random bind failed");
        std::vector<std::vector<Rational>> h;
        std::vector<Rational> rhs;
        const int witness[3]{static_cast<int>(random() % 7) - 3,
            static_cast<int>(random() % 7) - 3, static_cast<int>(random() % 7) - 3};
        const int count = 1 + random() % 3;
        for (int i = 0; i < count; ++i)
        {
            Row row;
            ZhangExactInteger right = 0;
            std::vector<Rational> dense(3);
            for (int c = 0; c < 3; ++c)
            {
                const int coefficient = static_cast<int>(random() % 7) - 3;
                if (coefficient) row[ids[c]] = coefficient;
                dense[c] = coefficient;
                right += coefficient * witness[c];
            }
            if (row.empty()) { --i; continue; }
            require(current.add(row, right, {proof}), "random feasible row rejected");
            h.push_back(dense);
            rhs.emplace_back(right);
        }
        const auto compiled = current.compiledCurrentDomain();
        require(compiled->feasible, "known feasible random domain rejected");
        for (int q = 0; q < 12; ++q)
        {
            Row sparse;
            std::vector<Rational> dense(3);
            for (int c = 0; c < 3; ++c)
            {
                const int coefficient = static_cast<int>(random() % 7) - 3;
                if (coefficient) sparse[ids[c]] = coefficient;
                dense[c] = coefficient;
            }
            Rational expected;
            const bool expectedKnown = referenceConsequence(h, rhs, dense, expected);
            const bool actualKnown = compiled->consequence(sparse, value);
            require(actualKnown == expectedKnown, "independent Gaussian consequence mismatch");
            if (actualKnown) require(Rational(value) == expected, "independent consequence value mismatch");
            ++comparisons;
        }
    }
    std::cout << "PASS current_epoch_cache comparisons=" << comparisons << '\n';
}
