#include "common/zhangPrivateGraphState.hpp"
#include <iostream>
#include <map>
#include <stdexcept>

static std::map<std::string, std::string> registry;
static std::string rootId = "fixture-owner";
static void require(bool condition, const char* reason)
{ if (!condition) throw std::runtime_error(reason); }
static void eraseFixture(const KFState& state)
{ registry.erase(zhangPrivateGraphIdentity(state)); }

int main()
{
    KFState source;
    source.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA] = rootId;
    source.x.resize(3); source.x << 1, 76.2739575, 63.3698112;
    source.P.resize(3, 3); source.P << 1, .1, .2, .1, 2, .3, .2, .3, 3;
    const auto meanBefore = source.x;
    const auto covarianceBefore = source.P;
    const std::string expectedGraph = "old-S-basis:nonzero-cycle:signal-arcs=17,19";
    registry[rootId] = expectedGraph;

    // The registry fixture is an independent ownership oracle. It deliberately
    // uses only IDs and literal graph payloads, not the production graph clone,
    // graph serializer, or the conditioner under test in other executables.
    auto clone = [](const KFState& from, const KFState& to, std::string& why)
    {
        const auto sourceId = zhangPrivateGraphIdentity(from);
        const auto destinationId = zhangPrivateGraphIdentity(to);
        require(sourceId != destinationId, "clone observed inherited authoritative ID");
        registry[destinationId] = registry.at(sourceId);
        why = "CLONED";
        return true;
    };
    std::string reason;
    {
        ZhangPrivateGraphState first, second;
        require(zhangMakePrivateGraphState(source, rootId, "final-product", clone,
            eraseFixture, first, reason), "private factory failed");
        const auto firstId = zhangPrivateGraphIdentity(*first);
        require(firstId != rootId && registry.at(firstId) == expectedGraph,
            "private graph was not isolated");
        require(zhangMakePrivateGraphState(source, rootId, "final-product", clone,
            eraseFixture, second, reason), "second live private factory failed");
        require(firstId != zhangPrivateGraphIdentity(*second), "live leases reused same ID");
        first->x(1) = 0; first->P.setZero();
        first.reset();
        require(registry.count(rootId) == 1 && registry.at(rootId) == expectedGraph &&
            registry.count(zhangPrivateGraphIdentity(*second)) == 1,
            "cleanup erased owner or another live branch");
    }
    require(registry.size() == 1 && registry.at(rootId) == expectedGraph,
        "normal destruction lost authoritative graph");

    auto failedClone = [](const KFState&, const KFState& to, std::string& why)
    {
        // Even if a failed callback staged some private data, cleanup cannot
        // touch its source and must reclaim only this newly minted lease.
        registry[zhangPrivateGraphIdentity(to)] = "partial-private";
        why = "INJECTED_CLONE_FAILURE";
        return false;
    };
    ZhangPrivateGraphState rejected;
    require(!zhangMakePrivateGraphState(source, rootId, "final-product", failedClone,
        eraseFixture, rejected, reason) && !rejected && reason == "INJECTED_CLONE_FAILURE",
        "clone failure was not propagated");
    require(registry.size() == 1 && registry.at(rootId) == expectedGraph,
        "clone failure cleanup erased owner or leaked private state");
    {
        ZhangPrivateGraphState existing;
        require(zhangMakePrivateGraphState(source, rootId, "final-product", clone,
            eraseFixture, existing, reason), "existing-destination fixture failed");
        const auto existingId = zhangPrivateGraphIdentity(*existing);
        const auto existingPointer = existing.get();
        require(!zhangMakePrivateGraphState(source, rootId, "final-product", failedClone,
            eraseFixture, existing, reason) && existing.get() == existingPointer &&
            registry.at(existingId) == expectedGraph && registry.at(rootId) == expectedGraph,
            "failed replacement destroyed previous owned destination");
    }

    bool conditioningEntered = false;
    try
    {
        ZhangPrivateGraphState branch;
        require(zhangMakePrivateGraphState(source, rootId, "final-product", clone,
            eraseFixture, branch, reason), "exception fixture creation failed");
        conditioningEntered = true;
        throw std::runtime_error("injected conditioning failure");
    }
    catch (const std::runtime_error&) {}
    require(conditioningEntered && registry.size() == 1 && registry.at(rootId) == expectedGraph,
        "conditioning exception erased owner or leaked branch");

    require(!zhangMakePrivateGraphState(source, "wrong-owner", "final-product", clone,
        eraseFixture, rejected, reason) && registry.size() == 1,
        "wrong owner acquired cleanup rights");
    require(!zhangMakePrivateGraphState(source, rootId, "bad/role", clone,
        eraseFixture, rejected, reason), "invalid role acquired cleanup rights");

    // A copied state without an explicit lease has no graph-cleanup authority.
    { ZhangPrivateGraphState plain(new KFState(source), ZhangPrivateGraphStateDeleter{}); }
    { ZhangPrivateGraphState alias(new KFState(source),
        ZhangPrivateGraphStateDeleter{rootId, rootId, eraseFixture}); }
    require(registry.count(rootId) && registry.at(rootId) == expectedGraph,
        "plain copied metadata acquired cleanup ownership");
    std::string abandonedId;
    {
        ZhangPrivateGraphState tampered;
        require(zhangMakePrivateGraphState(source, rootId, "final-product", clone,
            eraseFixture, tampered, reason), "metadata-tamper fixture failed");
        abandonedId = zhangPrivateGraphIdentity(*tampered);
        tampered->metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA] = rootId;
    }
    require(registry.at(rootId) == expectedGraph && registry.count(abandonedId),
        "changed branch metadata redirected cleanup to owner");
    registry.erase(abandonedId); // Fixture cleanup of the refused, orphaned lease.
    require((source.x.array() == meanBefore.array()).all() &&
        (source.P.array() == covarianceBefore.array()).all(),
        "authoritative numeric posterior changed");
    std::cout << "PASS private_graph_state independent_registry_lifecycle_oracle\n";
}
