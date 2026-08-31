#include <iostream>
#include <sstream>

#include "ambres/GNSSambres.hpp"

int traceLevel = 0;

namespace
{
bool check(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
    }
    return condition;
}

bool sameAmbiguityMap(const map<int, KFKey>& left, const map<int, KFKey>& right)
{
    if (left.size() != right.size())
    {
        return false;
    }

    for (const auto& [index, leftKey] : left)
    {
        const auto rightEntry = right.find(index);
        if (rightEntry == right.end())
        {
            return false;
        }

        const auto& rightKey = rightEntry->second;
        if (leftKey.type != rightKey.type || leftKey.str != rightKey.str ||
            leftKey.Sat.sys != rightKey.Sat.sys || leftKey.Sat.prn != rightKey.Sat.prn ||
            leftKey.num != rightKey.num)
        {
            return false;
        }
    }

    return true;
}
}  // namespace

int main()
{
    bool passed = true;
    std::ostringstream trace;

    GinAR_mtx offMatrix;
    GinAR_opt offOptions;
    offOptions.mode = E_ARmode::OFF;
    passed &= check(GNSS_AR(trace, offMatrix, offOptions) == 0, "off mode resolves nothing");
    passed &= check(offMatrix.diagnosticStatus == "MODE_OFF", "off mode diagnostic");

    GinAR_mtx weakMatrix;
    weakMatrix.aflt = VectorXd::Zero(4);
    weakMatrix.Paflt = MatrixXd::Identity(4, 4) * 100;
    GinAR_opt lambdaOptions;
    lambdaOptions.mode = E_ARmode::LAMBDA_ALT;
    lambdaOptions.sucthr = 0.9999;
    lambdaOptions.ratthr = 3;
    passed &= check(
        GNSS_AR(trace, weakMatrix, lambdaOptions) == 0,
        "weak lambda covariance resolves nothing"
    );
    passed &= check(
        weakMatrix.diagnosticStatus == "SUCCESS_RATE_BELOW_THRESHOLD",
        "weak lambda success-rate diagnostic"
    );
    passed &= check(
        weakMatrix.bootstrappedSuccessRate >= 0 &&
            weakMatrix.bootstrappedSuccessRate < lambdaOptions.sucthr,
        "weak lambda actual success rate"
    );
    passed &= check(
        weakMatrix.selectedDecorrelatedAmbiguityCount == 1,
        "weak lambda selected ambiguity count"
    );

    GinAR_mtx oneDimensionalMatrix;
    oneDimensionalMatrix.aflt = VectorXd::Constant(1, 0.1);
    oneDimensionalMatrix.Paflt = MatrixXd::Identity(1, 1) * 1e-4;
    GinAR_opt oneDimensionalOptions;
    oneDimensionalOptions.mode = E_ARmode::LAMBDA_ALT;
    oneDimensionalOptions.sucthr = 0.9999;
    oneDimensionalOptions.ratthr = 3;
    passed &= check(
        GNSS_AR(trace, oneDimensionalMatrix, oneDimensionalOptions) == 0,
        "production lambda retains the three-dimension minimum"
    );
    passed &= check(
        oneDimensionalMatrix.diagnosticStatus ==
            "INSUFFICIENT_DECORRELATED_AMBIGUITIES",
        "production minimum reports insufficient dimension"
    );
    oneDimensionalOptions.minimumDecorrelatedAmbiguityCount = 1;
    passed &= check(
        GNSS_AR(trace, oneDimensionalMatrix, oneDimensionalOptions) == 1,
        "diagnostic lambda can search a one-dimensional remainder"
    );
    passed &= check(
        oneDimensionalMatrix.diagnosticStatus == "RESOLVED_RATIO_ACCEPTED",
        "one-dimensional diagnostic remainder passes the ratio test"
    );

    GinAR_mtx roundMatrix;
    roundMatrix.aflt = VectorXd::Constant(2, 0.01);
    roundMatrix.Paflt = MatrixXd::Identity(2, 2) * 1e-4;
    GinAR_opt roundOptions;
    roundOptions.mode = E_ARmode::ROUND;
    roundOptions.sucthr = 0.99;
    roundOptions.ratthr = 3;
    passed &= check(GNSS_AR(trace, roundMatrix, roundOptions) == 2, "round resolves both values");
    passed &= check(
        roundMatrix.diagnosticStatus == "RESOLVED_BY_ROUNDING",
        "round success diagnostic"
    );
    passed &= check(
        roundMatrix.selectedDecorrelatedAmbiguityCount == 2,
        "round selected ambiguity count"
    );

    // Regression for a false ratio-test acceptance in the LAMBDA candidate
    // search.  The depth-first search encounters squared norms 1.371, 11.171
    // and 11.171 before the actual second-best norm 3.971.  Stopping when the
    // third candidate is encountered therefore inflates the ratio from
    // 3.971 / 1.371 = 2.896 to 11.171 / 1.371 = 8.148.
    GinAR_mtx ratioMatrix;
    ratioMatrix.aflt.resize(3);
    ratioMatrix.aflt << 0.01, 0.01, 0.37;
    ratioMatrix.Paflt = MatrixXd::Identity(3, 3) * 0.1;
    GinAR_opt ratioOptions;
    ratioOptions.mode = E_ARmode::LAMBDA_ALT;
    ratioOptions.sucthr = 0;
    ratioOptions.ratthr = 3;
    ratioOptions.nset = 2;
    passed &= check(
        GNSS_AR(trace, ratioMatrix, ratioOptions) == 0,
        "lambda rejects when the true second-best candidate fails the ratio test"
    );
    passed &= check(
        ratioMatrix.diagnosticStatus == "RATIO_BELOW_THRESHOLD",
        "lambda true second-best ratio diagnostic"
    );
    passed &= check(
        ratioMatrix.integerCandidateCount == 2,
        "lambda candidate pool respects the configured size"
    );
    passed &= check(
        std::abs(ratioMatrix.bestSquaredNorm - 1.371) < 1e-12,
        "lambda finds the true best squared norm"
    );
    passed &= check(
        std::abs(ratioMatrix.secondSquaredNorm - 3.971) < 1e-12,
        "lambda finds the true second-best squared norm"
    );
    passed &= check(
        ratioMatrix.solutionRatio < ratioOptions.ratthr,
        "lambda uses the non-inflated solution ratio"
    );

    GinAR_mtx tiedRatioMatrix;
    tiedRatioMatrix.aflt.resize(3);
    tiedRatioMatrix.aflt << 0.5, 0.01, 0.01;
    tiedRatioMatrix.Paflt = MatrixXd::Identity(3, 3) * 0.1;
    GinAR_opt tiedRatioOptions = ratioOptions;
    tiedRatioOptions.nset = 8;
    passed &= check(
        GNSS_AR(trace, tiedRatioMatrix, tiedRatioOptions) == 0,
        "lambda rejects equal-distance best candidates"
    );
    passed &= check(
        tiedRatioMatrix.diagnosticStatus == "RATIO_BELOW_THRESHOLD",
        "lambda equal-distance ratio diagnostic"
    );
    passed &= check(
        tiedRatioMatrix.integerCandidateCount == 2,
        "lambda retains both equal-distance candidates"
    );
    passed &= check(
        std::abs(tiedRatioMatrix.bestSquaredNorm - tiedRatioMatrix.secondSquaredNorm) < 1e-12,
        "lambda equal-distance candidates remain distinct"
    );
    passed &= check(
        std::abs(tiedRatioMatrix.solutionRatio - 1) < 1e-12,
        "lambda tie produces unit ratio"
    );

    GinAR_mtx commonSetTieMatrix;
    commonSetTieMatrix.aflt.resize(3);
    commonSetTieMatrix.aflt << 0.5, 0.5, 0.01;
    commonSetTieMatrix.Paflt = MatrixXd::Identity(3, 3) * 0.1;
    GinAR_opt commonSetTieOptions = ratioOptions;
    commonSetTieOptions.mode = E_ARmode::LAMBDA_AL2;
    passed &= check(
        GNSS_AR(trace, commonSetTieMatrix, commonSetTieOptions) == 1,
        "lambda common-set retains only the coordinate shared by cutoff ties"
    );
    passed &= check(
        commonSetTieMatrix.integerCandidateCount == 4,
        "lambda common-set retains every candidate tied at the cutoff"
    );
    MatrixXd expectedCommonSetTransform(1, 3);
    expectedCommonSetTransform << 0, 0, 1;
    passed &= check(
        commonSetTieMatrix.Ztrs.isApprox(expectedCommonSetTransform, 1e-12) &&
            commonSetTieMatrix.zfix.size() == 1 &&
            std::abs(commonSetTieMatrix.zfix(0)) < 1e-12,
        "lambda common-set identifies only the third coordinate as common"
    );

    GinAR_mtx bieTieMatrix;
    bieTieMatrix.aflt.resize(3);
    bieTieMatrix.aflt << 0.5, 0.5, 0.01;
    bieTieMatrix.Paflt = MatrixXd::Identity(3, 3) * 0.1;
    GinAR_opt bieTieOptions = ratioOptions;
    bieTieOptions.mode = E_ARmode::LAMBDA_BIE;
    passed &= check(
        GNSS_AR(trace, bieTieMatrix, bieTieOptions) == 3,
        "lambda BIE returns all ambiguity coordinates for cutoff ties"
    );
    passed &= check(
        bieTieMatrix.integerCandidateCount == 4,
        "lambda BIE retains every candidate tied at the cutoff"
    );
    VectorXd expectedBieTie(3);
    expectedBieTie << 0.5, 0.5, 0;
    passed &= check(
        bieTieMatrix.zfix.isApprox(expectedBieTie, 1e-12),
        "lambda BIE preserves the symmetric cutoff-tie mean"
    );

    GinAR_mtx ambiguityMatrix;
    ambiguityMatrix.aflt.resize(7);
    ambiguityMatrix.aflt << 10.25, 15.25, -3.75, 20.4, 22.4, 7.1, 9.2;
    ambiguityMatrix.Paflt = MatrixXd::Identity(7, 7);
    ambiguityMatrix.Paflt.diagonal() << 4, 1, 2, 3, 0.5, 1, 1;

    auto addAmbiguity = [&](int index, const char* receiver, E_Sys system, int prn, E_ObsCode code)
    {
        KFKey key;
        key.type = KF::AMBIGUITY;
        key.str = receiver;
        key.Sat = SatSys(system, prn);
        key.num = static_cast<int>(code);
        ambiguityMatrix.ambmap[index] = key;
    };
    addAmbiguity(0, "A", E_Sys::GPS, 1, E_ObsCode::L1C);
    addAmbiguity(1, "A", E_Sys::GPS, 3, E_ObsCode::L1C);
    addAmbiguity(2, "A", E_Sys::GPS, 5, E_ObsCode::L1C);
    addAmbiguity(3, "A", E_Sys::GPS, 1, E_ObsCode::L2W);
    addAmbiguity(4, "A", E_Sys::GPS, 8, E_ObsCode::L2W);
    addAmbiguity(5, "B", E_Sys::GAL, 2, E_ObsCode::L1C);
    addAmbiguity(6, "B", E_Sys::GAL, 7, E_ObsCode::L1C);

    map<E_Sys, bool> receiverPivot = {
        {E_Sys::GPS, true},
        {E_Sys::GAL, false},
    };
    const auto receiverTransform =
        buildReceiverAmbiguityIntegerTransform(ambiguityMatrix, receiverPivot);
    passed &= check(receiverTransform.matrix.rows() == 5, "receiver transform row count");
    passed &= check(receiverTransform.matrix.cols() == 7, "receiver transform column count");
    passed &= check(
        receiverTransform.singleDifferencedGroupCount == 2,
        "receiver transform differenced group count"
    );
    passed &= check(receiverTransform.identityGroupCount == 1, "receiver transform identity group");
    passed &= check(
        receiverTransform.droppedSingletonGroupCount == 0,
        "receiver transform no dropped singleton"
    );

    VectorXd expectedCoordinates(5);
    expectedCoordinates << -5, -19, -2, 7.1, 9.2;
    const VectorXd integerCoordinates = receiverTransform.matrix * ambiguityMatrix.aflt;
    passed &= check(
        integerCoordinates.isApprox(expectedCoordinates, 1e-12),
        "receiver transform produces satellite single differences"
    );

    VectorXd datumShift = VectorXd::Zero(7);
    datumShift << 0.37, 0.37, 0.37, -0.22, -0.22, 0, 0;
    passed &= check(
        (receiverTransform.matrix * (ambiguityMatrix.aflt + datumShift))
            .isApprox(integerCoordinates, 1e-12),
        "receiver transform removes receiver signal datum"
    );

    const MatrixXd transformedCovariance =
        receiverTransform.matrix * ambiguityMatrix.Paflt * receiverTransform.matrix.transpose();
    VectorXd expectedVariance(5);
    expectedVariance << 5, 3, 3.5, 1, 1;
    passed &= check(
        transformedCovariance.diagonal().isApprox(expectedVariance, 1e-12),
        "receiver transform propagates covariance"
    );
    passed &= check(
        std::abs(transformedCovariance(0, 1) - 1) < 1e-12,
        "receiver transform retains covariance induced by the shared pivot"
    );
    passed &= check(
        std::abs(transformedCovariance(0, 2)) < 1e-12,
        "receiver transform does not couple independent signal groups"
    );
    passed &= check(
        transformedCovariance.isApprox(transformedCovariance.transpose(), 1e-12),
        "receiver transformed covariance is symmetric"
    );

    GinAR_mtx integerFeedbackMatrix;
    integerFeedbackMatrix.Ztrs.resize(2, 5);
    integerFeedbackMatrix.Ztrs << 2, -1, 3, 0, 0, 0, 1, -2, 3, -1;
    integerFeedbackMatrix.zfix.resize(2);
    integerFeedbackMatrix.zfix << 7, -4;
    const VectorXd originalFixedIntegers = integerFeedbackMatrix.zfix;
    passed &= check(
        mapIntegerAmbiguityConstraintsToOriginalState(
            integerFeedbackMatrix,
            receiverTransform,
            ambiguityMatrix.ambmap
        ),
        "integer feedback maps resolved constraints to original ambiguity states"
    );
    MatrixXd expectedOriginalFeedback(2, 7);
    expectedOriginalFeedback << 2, -1, -1, 3, -3, 0, 0, 0, -1, 1, -2, 2, 3, -1;
    passed &= check(
        integerFeedbackMatrix.Ztrs.isApprox(expectedOriginalFeedback, 1e-12),
        "integer feedback applies non-identity Z times receiver transform"
    );
    passed &= check(
        sameAmbiguityMap(integerFeedbackMatrix.ambmap, ambiguityMatrix.ambmap),
        "integer feedback restores the original ambiguity map"
    );
    passed &= check(
        integerFeedbackMatrix.zfix.isApprox(originalFixedIntegers, 1e-12),
        "integer feedback preserves fixed integer values"
    );

    GinAR_mtx invalidFeedbackMatrix;
    invalidFeedbackMatrix.Ztrs = MatrixXd::Identity(2, 2);
    invalidFeedbackMatrix.zfix = VectorXd::Zero(2);
    const MatrixXd invalidFeedbackBefore = invalidFeedbackMatrix.Ztrs;
    passed &= check(
        !mapIntegerAmbiguityConstraintsToOriginalState(
            invalidFeedbackMatrix,
            receiverTransform,
            ambiguityMatrix.ambmap
        ),
        "integer feedback rejects incompatible Z and receiver-transform dimensions"
    );
    passed &= check(
        invalidFeedbackMatrix.Ztrs.isApprox(invalidFeedbackBefore, 1e-12) &&
            invalidFeedbackMatrix.ambmap.empty(),
        "failed integer feedback mapping leaves its input unchanged"
    );

    GinAR_mtx singletonMatrix;
    singletonMatrix.aflt = VectorXd::Constant(1, 12.25);
    singletonMatrix.Paflt = MatrixXd::Identity(1, 1);
    KFKey singletonKey;
    singletonKey.type = KF::AMBIGUITY;
    singletonKey.str = "C";
    singletonKey.Sat = SatSys(E_Sys::GPS, 9);
    singletonKey.num = static_cast<int>(E_ObsCode::L1C);
    singletonMatrix.ambmap[0] = singletonKey;
    const auto singletonTransform =
        buildReceiverAmbiguityIntegerTransform(singletonMatrix, receiverPivot);
    passed &= check(singletonTransform.matrix.rows() == 0, "singleton datum has no integer row");
    passed &= check(
        singletonTransform.droppedSingletonGroupCount == 1,
        "singleton datum is reported"
    );

    GinAR_mtx complementInput;
    complementInput.aflt.resize(3);
    complementInput.aflt << 10.2, 20.3, 30.4;
    complementInput.Paflt = MatrixXd::Identity(3, 3);
    GinAR_mtx acceptedPartial;
    acceptedPartial.fullDecorrelatedTransform.resize(3, 3);
    acceptedPartial.fullDecorrelatedTransform << 1, 0, 0, -1, 1, 0, 0, -1, 1;
    acceptedPartial.Ztrs = acceptedPartial.fullDecorrelatedTransform.bottomRows(1);
    acceptedPartial.zfix = VectorXd::Constant(1, 10);
    const auto complement =
        buildConditionalIntegerComplement(complementInput, acceptedPartial);
    passed &= check(
        complement.diagnosticStatus == "CONDITIONAL_INTEGER_COMPLEMENT_READY",
        "partial integer solution exposes a conditional integer complement"
    );
    VectorXd expectedConditionalFloat(2);
    expectedConditionalFloat << 10.2, 10.15;
    passed &= check(
        complement.ambiguityResolution.aflt.isApprox(expectedConditionalFloat, 1e-12),
        "conditional complement mean uses the accepted integer innovation"
    );
    MatrixXd expectedConditionalCovariance(2, 2);
    expectedConditionalCovariance << 1, -1, -1, 1.5;
    passed &= check(
        complement.ambiguityResolution.Paflt.isApprox(
            expectedConditionalCovariance,
            1e-12
        ),
        "conditional complement covariance uses the Schur complement"
    );
    passed &= check(
        complement.transformToInputCoordinates.isApprox(
            acceptedPartial.fullDecorrelatedTransform.topRows(2),
            1e-12
        ),
        "conditional complement preserves the unused unimodular rows"
    );
    MatrixXd stackedIntegerBasis(3, 3);
    stackedIntegerBasis << acceptedPartial.Ztrs,
                           complement.transformToInputCoordinates;
    passed &= check(
        std::abs(stackedIntegerBasis.determinant()) == 1,
        "accepted rows and complement rows retain a full unimodular integer basis"
    );

    GinAR_mtx invalidPartial = acceptedPartial;
    invalidPartial.Ztrs = acceptedPartial.fullDecorrelatedTransform.middleRows(1, 1);
    const auto invalidComplement =
        buildConditionalIntegerComplement(complementInput, invalidPartial);
    passed &= check(
        invalidComplement.diagnosticStatus == "FIXED_ROWS_NOT_DECORRELATED_TAIL",
        "conditional complement rejects an unsupported non-tail fixed-row subset"
    );

    if (passed)
    {
        std::cout << "PASS: gnss_ambiguity_diagnostics_tests\n";
        return 0;
    }
    return 1;
}
