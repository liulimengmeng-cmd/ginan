#pragma once

#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Identity belongs to the accepted decision, NEVER to a transported row hash.
// Bounds are conditional on all parents being correct; union accounting visits
// every original decision once. This is a model-bound contract, not empirical
// calibration of the GNSS receiver or of an adaptive PAR procedure.
struct ZhangIntegerDecisionProof;
using ZhangDecisionProofPtr = std::shared_ptr<const ZhangIntegerDecisionProof>;
using ZhangDecisionProofs = std::vector<ZhangDecisionProofPtr>;
struct ZhangIntegerDecisionProof
{
    std::string id;
    std::string originalStatement;
    std::string observationProvenance;
    double conditionalFailureBound = 1;
    ZhangDecisionProofs parents;
};

struct ZhangDecisionRiskClosure
{
    bool valid = false;
    double bound = 0;
    std::string reason = "MISSING_DECISION_PROOF";
    std::map<std::string, ZhangDecisionProofPtr> atoms;
};

struct ZhangIncrementalDecisionRisk
{
    bool valid = false;
    int baselineAtoms = 0;
    int unionAtoms = 0;
    int addedAtoms = 0;
    double bound = 0;
    std::string reason = "NOT_EVALUATED";
};

inline ZhangDecisionRiskClosure zhangDecisionRiskClosure(const ZhangDecisionProofs& roots)
{
    ZhangDecisionRiskClosure result;
    std::set<std::string> visiting;
    auto visit = [&](auto&& self, const ZhangDecisionProofPtr& proof) -> bool
    {
        if (!proof || proof->id.empty() || proof->originalStatement.empty() ||
            proof->observationProvenance.empty() ||
            !std::isfinite(proof->conditionalFailureBound) ||
            proof->conditionalFailureBound < 0 || proof->conditionalFailureBound > 1)
            return false;
        if (visiting.count(proof->id))
        { result.reason = "DECISION_DEPENDENCY_CYCLE"; return false; }
        auto prior = result.atoms.find(proof->id);
        if (prior != result.atoms.end())
        {
            // Two different immutable objects claiming an ID are not silently
            // merged, even if their current numerical integer values agree.
            if (prior->second.get() != proof.get())
            { result.reason = "DECISION_ID_COLLISION"; return false; }
            return true;
        }
        visiting.insert(proof->id);
        for (const auto& parent : proof->parents)
            if (!self(self, parent)) return false;
        visiting.erase(proof->id);
        result.atoms.emplace(proof->id, proof);
        result.bound += proof->conditionalFailureBound;
        return true;
    };
    for (const auto& root : roots) if (!visit(visit, root)) return result;
    result.valid = true;
    result.reason = "DECISION_DEPENDENCIES_CLOSED";
    return result;
}

inline ZhangDecisionProofs zhangMergeDecisionProofs(
    const ZhangDecisionProofs& first, const ZhangDecisionProofs& second)
{
    auto result = first;
    for (const auto& proof : second)
    {
        bool present = false;
        for (const auto& old : result) present |= old.get() == proof.get();
        if (!present) result.push_back(proof);
    }
    return result;
}

/** Risk contributed by candidate decisions that are not already ancestors of
 * the baseline certificate.  Comparing two accumulated floating-point sums is
 * unsafe here: a shared proof can be rounded differently after a view or HNF
 * change.  Immutable proof IDs define the statistical decision exactly once;
 * an ID collision remains fail-closed through zhangDecisionRiskClosure. */
inline ZhangIncrementalDecisionRisk zhangIncrementalDecisionRisk(
    const ZhangDecisionProofs& baseline,
    const ZhangDecisionProofs& candidate)
{
    ZhangIncrementalDecisionRisk result;
    const auto baselineClosure = zhangDecisionRiskClosure(baseline);
    if (!baselineClosure.valid)
    {
        result.reason = baselineClosure.reason;
        return result;
    }
    const auto unionClosure = zhangDecisionRiskClosure(
        zhangMergeDecisionProofs(baseline, candidate));
    if (!unionClosure.valid)
    {
        result.reason = unionClosure.reason;
        return result;
    }
    result.baselineAtoms = static_cast<int>(baselineClosure.atoms.size());
    result.unionAtoms = static_cast<int>(unionClosure.atoms.size());
    for (const auto& [id, proof] : unionClosure.atoms)
    {
        if (baselineClosure.atoms.contains(id)) continue;
        result.bound = std::min(1.0,
            result.bound + proof->conditionalFailureBound);
        result.addedAtoms++;
    }
    result.valid = true;
    result.reason = "NONE";
    return result;
}
