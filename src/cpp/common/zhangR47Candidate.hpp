#pragma once
#include "common/zhangRatioOnly.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangProductRelationSolver.hpp"
#include "common/zhangIntegerConditioner.hpp"

/** A proposed integer transaction owns its complete conditioning domain.
 * Product consequences are a separate view and cannot replace this domain. */
struct ZhangR47Candidate
{
    std::string sourcePosteriorId, authoritativeCycleChartId, frontendSemanticId;
    std::vector<int> stateIndices;
    ZhangExactMatrix admittedStateConditioners, newIntegerConstraints, jointRows;
    ZhangExactVector admittedStateValues, newIntegerValues, jointValues;
    std::vector<std::map<std::string,ZhangExactInteger>> physicalFunctionals;
    ZhangDecisionProofs allDecisionParents;
    double riskBound = 1;
    std::string riskScope = "CURRENT_CERTIFICATE_ANCESTOR_UNION";
    int sourceIntegerParentCount = -1;
    ZhangProductIntegerConstraintSet productConsequences;
};

inline bool zhangR47AffineIntegerFeasible(const ZhangExactMatrix& rows,
    const ZhangExactVector& values, int dimension)
{
    if (rows.size()!=values.size() || !zhangExactRectangularMatrix(rows,dimension)) return false;
    if (rows.empty()) return true;
    ZhangExactMatrix columns(dimension,ZhangExactVector(rows.size()));
    for (int c=0;c<dimension;++c) for (std::size_t r=0;r<rows.size();++r) columns[c][r]=rows[r][c];
    return zhangIntegerRowLatticeContains(columns,values,false).contained;
}

inline bool zhangR47CandidateContractValid(const ZhangR47Candidate& c)
{
    const int n=c.stateIndices.size();
    if (c.sourcePosteriorId.empty() || c.authoritativeCycleChartId.empty() ||
        c.frontendSemanticId.empty() || c.sourceIntegerParentCount!=0 || n<=0 ||
        c.allDecisionParents.empty() || c.physicalFunctionals.size()!=c.jointRows.size() ||
        !zhangR47AffineIntegerFeasible(c.jointRows,c.jointValues,n)) return false;
    const auto risk=zhangDecisionRiskClosure(c.allDecisionParents);
    if (!risk.valid || zhangRatioStatisticalReject(risk.bound>1e-3) || c.riskBound!=risk.bound) return false;
    auto rows=c.admittedStateConditioners, valuesRows=c.newIntegerConstraints;
    rows.insert(rows.end(),valuesRows.begin(),valuesRows.end());
    auto values=c.admittedStateValues;
    values.insert(values.end(),c.newIntegerValues.begin(),c.newIntegerValues.end());
    // The whole affine union is immutable; removing a conditioner is a new
    // candidate from the clean root, never an update of a conditional state.
    const auto expected=zhangExactRowHermiteNormalForm(rows,values);
    const auto actual=zhangExactRowHermiteNormalForm(c.jointRows,c.jointValues);
    return expected.consistent and actual.consistent and expected.basis==actual.basis and expected.values==actual.values;
}

/** A writer may attach physical identities to a frozen candidate, but it may
 * not substitute another frontend product chart or silently change an affine
 * offset.  The exact q=Mx+c chart (including its completeness flag) is part
 * of the delivery contract. */
inline bool zhangR47CandidateProductChartMatches(
    const ZhangR47Candidate& candidate,
    const ZhangProductIntegerConstraintSet& delivered)
{
    const auto& frozen = candidate.productConsequences;
    const auto ancestors = zhangDecisionRiskClosure(candidate.allDecisionParents);
    const auto required = zhangDecisionRiskClosure(delivered.decisionProofs);
    if (!ancestors.valid || !required.valid ||
        frozen.fullJointProductMappingExact !=
            delivered.fullJointProductMappingExact ||
        frozen.system != delivered.system ||
        frozen.firstObservable != delivered.firstObservable ||
        frozen.secondObservable != delivered.secondObservable ||
        frozen.referenceSatellite != delivered.referenceSatellite ||
        frozen.coordinateSatellites != delivered.coordinateSatellites ||
        frozen.fullJointProductNetworkRows !=
            delivered.fullJointProductNetworkRows ||
        frozen.fullJointProductAffineOffsets !=
            delivered.fullJointProductAffineOffsets ||
        frozen.networkRows != delivered.networkRows ||
        frozen.networkIntegers != delivered.networkIntegers ||
        frozen.jointProductRows != delivered.jointProductRows)
        return false;
    for (const auto& [id, proof] : required.atoms)
    {
        const auto found = ancestors.atoms.find(id);
        if (found == ancestors.atoms.end() ||
            found->second.get() != proof.get()) return false;
    }
    return true;
}
