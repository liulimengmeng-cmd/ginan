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
