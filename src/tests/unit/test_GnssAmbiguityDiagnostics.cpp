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

    if (passed)
    {
        std::cout << "PASS: gnss_ambiguity_diagnostics_tests\n";
        return 0;
    }
    return 1;
}
