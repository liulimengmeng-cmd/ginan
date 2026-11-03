#include "common/zhangCurrentEpochEvidence.hpp"
#include <iostream>
#include <stdexcept>

namespace {
int checks = 0;
void require(bool ok, const char* why)
{ ++checks; if (!ok) throw std::runtime_error(why); }
ZhangDecisionProofs proof(const std::string& id, const std::string& source = "runtime|PPP_COMMITTED|E123")
{
    auto p = std::make_shared<ZhangIntegerDecisionProof>();
    p->id = id; p->originalStatement = "controlled current overlap";
    p->observationProvenance = source; p->conditionalFailureBound = 1e-7;
    return {p};
}
KFKey key(int satellite)
{
    KFKey k; k.type = KF::AMBIGUITY; k.str = "R";
    k.Sat = SatSys(E_Sys::GPS, satellite); return k;
}
ZhangCurrentEpochEvidence parent()
{
    ZhangCurrentEpochEvidence e;
    e.runtimeId = "runtime"; e.system = E_Sys::GPS; e.epochGpstSeconds = 123;
    e.posteriorId = "float-posterior"; e.floatStateRoot = "immutable-root";
    e.chartId = "physical-chart"; e.floatPosteriorPsdStateOrderId = "full-order";
    e.floatPosteriorPsdRootId = "immutable-root"; e.floatPosteriorPsdValidated = true;
    e.valid = true;
    e.rows = {{{"base", 1}}}; e.values = {4}; e.rowProofs = {proof("base")};
    require(e.addNumericConstraint({{key(2), 1}}, 4, {{"base", 1}}), "base binding");
    return e;
}
bool same(const ZhangCurrentEpochEvidence& a, const ZhangCurrentEpochEvidence& b)
{
    if (a.currentDomainIdentity() != b.currentDomainIdentity() || a.rows != b.rows ||
        a.values != b.values || a.rowProofs != b.rowProofs ||
        a.eventOverlapRows != b.eventOverlapRows || a.eventOverlapValues != b.eventOverlapValues ||
        a.eventOverlapProofs != b.eventOverlapProofs || a.eventOverlapWitnesses != b.eventOverlapWitnesses ||
        a.numericConstraints.size() != b.numericConstraints.size() ||
        a.physicalNumericBindings.size() != b.physicalNumericBindings.size()) return false;
    for (std::size_t i = 0; i < a.numericConstraints.size(); ++i)
        if (a.numericConstraints[i].rowCycles != b.numericConstraints[i].rowCycles ||
            a.numericConstraints[i].value != b.numericConstraints[i].value) return false;
    for (std::size_t i = 0; i < a.physicalNumericBindings.size(); ++i)
    {
        const auto& x = a.physicalNumericBindings[i]; const auto& y = b.physicalNumericBindings[i];
        if (x.physical != y.physical || x.rowCycles != y.rowCycles || x.value != y.value ||
            x.runtimeId != y.runtimeId || x.system != y.system || x.epochGpstSeconds != y.epochGpstSeconds ||
            x.posteriorId != y.posteriorId || x.rootId != y.rootId || x.chartId != y.chartId ||
            x.stateOrderId != y.stateOrderId) return false;
    }
    return true;
}
}
int main()
{
    const ZhangCurrentEpochPhysicalRow a{{"old-A", 1}}, b{{"old-B", 1}};
    const ZhangCoordinateRow na{{key(3), 1}}, nb{{key(4), 1}};
    auto base = parent(); auto old = base; auto untouched = base;
    const auto pa = proof("A"), pb = proof("B");
    require(old.addEventOverlap(a, 7, pa, "event") && old.addNumericConstraint(na, 7, a), "old A");
    require(old.addEventOverlap(b, 9, pb, "event") && old.addNumericConstraint(nb, 9, b), "old B");
    std::string reason;
    ZhangCurrentEpochEventTransaction tx(base);
    require(tx.stage(a, 7, pa, "event", na, reason), "stage A");
    require(tx.stage(b, 9, pb, "event", nb, reason), "stage B");
    require(same(base, untouched), "staging mutated parent");
    require(tx.prepare(reason) && tx.prepare(reason), "repeat prepare");
    require(same(base, untouched) && !tx.committed(), "prepare published evidence");
    require(tx.commit(base, reason) && tx.committed(), "atomic commit");
    require(same(base, old), "delta differs from independent old copy+add result");
    require(base.validatePhysicalNumericBindings(reason), "full binding audit");
    const auto size = base.eventOverlapRows.size(); const auto gen = base.constraintGeneration;
    require(tx.commit(base, reason) && base.eventOverlapRows.size() == size &&
        base.constraintGeneration == gen, "commit not idempotent");
    auto anotherOwner = base;
    require(!tx.commit(anotherOwner, reason), "different owner reused committed transaction");
    require(!tx.stage(a, 7, pa, "event", na, reason), "closed transaction staged row");

    ZhangCurrentEpochEventTransaction duplicate(base);
    require(duplicate.stage(a, 7, pa, "event", na, reason) && duplicate.commit(base, reason), "duplicate stage");
    require(same(base, old), "duplicate changed generation or bindings");
    for (int failure = 0; failure < 7; ++failure)
    {
        auto e = parent(); auto before = e; ZhangCurrentEpochEventTransaction bad(e);
        bool accepted = false;
        if (failure == 0) accepted = bad.stage(a, 7, proof("wrong", "another|PPP_COMMITTED|E123"), "event", na, reason);
        if (failure == 1) accepted = bad.stage(a, 7, proof("stale", "runtime|PPP_COMMITTED|E122"), "event", na, reason);
        if (failure == 2) accepted = bad.stage(a, 7, {}, "event", na, reason);
        if (failure == 3) accepted = bad.stage(a, 7, pa, "", na, reason);
        if (failure == 4) accepted = bad.stage(a, 7, pa, "event", {}, reason);
        if (failure == 5) accepted = bad.stage({{"base", 1}}, 5, pa, "event", na, reason);
        if (failure == 6) accepted = bad.stage(a, 7, pa, "event", {{key(2), 1}}, reason);
        require(!accepted && same(e, before), "invalid stage admitted or mutated parent");
    }
    auto e = parent(); ZhangCurrentEpochEventTransaction conflict(e);
    require(conflict.stage(a, 7, pa, "event", na, reason), "first valid stage");
    require(!conflict.stage(a, 8, pa, "event", na, reason), "intra-delta RHS conflict");
    require(!conflict.stage(a, 7, pa, "other-event", na, reason), "intra-delta witness conflict");
    require(!conflict.stage(a, 7, pa, "event", nb, reason), "intra-delta binding conflict");
    require(conflict.commit(e, reason) && e.eventOverlapRows.size() == 1 &&
        e.eventOverlapValues.front() == 7, "failed stage polluted prior valid delta");
    auto collisionBase = parent(); const auto collisionBefore = collisionBase;
    ZhangCurrentEpochEventTransaction collision(collisionBase);
    require(collision.stage(a, 7, proof("base"), "event", na, reason), "collision staged fixture");
    require(!collision.commit(collisionBase, reason) &&
        reason == "EVENT_DELTA_BASE_UNION_DELTA_PROOF_CONFLICT" && same(collisionBase, collisionBefore),
        "base UNION delta proof ID collision or its reason was lost");

    for (int field = 0; field < 9; ++field)
    {
        auto scoped = parent(); ZhangCurrentEpochEventTransaction stale(scoped);
        require(stale.stage(a, 7, pa, "event", na, reason), "scope stage");
        if (field == 0) scoped.epochGpstSeconds++;
        if (field == 1) scoped.runtimeId += "other";
        if (field == 2) scoped.posteriorId += "other";
        if (field == 3) scoped.floatStateRoot += "other";
        if (field == 4) scoped.chartId += "other";
        if (field == 5) scoped.floatPosteriorPsdStateOrderId += "other";
        if (field == 6) scoped.floatPosteriorPsdValidated = false;
        if (field == 7) scoped.invalidateCurrentDomain();
        if (field == 8) scoped.chart.cataloguedSignals.insert("new-signal");
        const auto before = scoped;
        require(!stale.commit(scoped, reason) && same(scoped, before), "changed parent committed");
    }
    // Interrupt every fallible preparation/reserve boundary deterministically.
    // Capacity changes are harmless; logical registry and generation must match.
    int cancellations = 0;
    for (int limit = 1; limit < 50; ++limit)
    {
        auto cancelBase = parent(); const auto before = cancelBase;
        ZhangCurrentEpochEventTransaction cancel(cancelBase);
        require(cancel.stage(a, 7, pa, "event", na, reason), "cancel stage");
        ZhangExactWorkBudget budget; budget.maximumCheckpoints = limit;
        bool done = false;
        try { ZhangScopedExactWorkBudget scope(budget); done = cancel.commit(cancelBase, reason); }
        catch (const ZhangExactWorkBudgetExceeded&) { ++cancellations; }
        if (done) break;
        require(same(cancelBase, before) && !cancel.committed(), "cancel partially published delta");
        // Same frozen-parent retry reuses private validated work safely.
        require(cancel.commit(cancelBase, reason), "canceled same-parent transaction not resumable");
    }
    require(cancellations >= 7, "cancellation boundaries not exercised");
    auto large = parent();
    large.floatPosteriorPsdStateOrderId = std::string(16384, 's');
    large.physicalNumericBindings.front().stateOrderId = large.floatPosteriorPsdStateOrderId;
    const auto oldStateOrderAllocation = large.physicalNumericBindings.front().stateOrderId.data();
    ZhangCurrentEpochEventTransaction shared(large);
    require(shared.stage(a, 7, pa, "event", na, reason) && shared.commit(large, reason), "large-scope commit");
    require(large.physicalNumericBindings.front().stateOrderId.data() == oldStateOrderAllocation,
        "commit deep-copied the existing state-order binding");
    auto overflow = parent(); overflow.constraintGeneration = std::numeric_limits<std::uint64_t>::max();
    const auto overflowBefore = overflow; ZhangCurrentEpochEventTransaction overflowTx(overflow);
    require(overflowTx.stage(a, 7, pa, "event", na, reason) && !overflowTx.commit(overflow, reason) &&
        reason == "EVENT_DELTA_CONSTRAINT_GENERATION_OVERFLOW" && same(overflow, overflowBefore),
        "constraint generation wrapped or partially published");
    std::cout << "current_event_evidence_delta PASS checks=" << checks
        << " cancellation_boundaries=" << cancellations << '\n';
}
