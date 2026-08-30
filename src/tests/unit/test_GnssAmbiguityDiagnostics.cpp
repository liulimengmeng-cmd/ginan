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

    if (passed)
    {
        std::cout << "PASS: gnss_ambiguity_diagnostics_tests\n";
        return 0;
    }
    return 1;
}
