#include "common/zhangPublishedContinuity.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

static int checks = 0;
static void require(bool result, const char* reason)
{
    ++checks;
    if (!result) throw std::runtime_error(reason);
}

int main()
{
    ZhangPublishedContinuityContract contract;
    contract.bundleId = "1405216830:PRODUCT_FIXED";
    require(zhangPublishedContinuityValid(contract, 1405216830, "new"), "honest UNKNOWN invalid");
    std::ostringstream unknown;
    zhangWritePublishedContinuityCsv(unknown, contract);
    require(unknown.str().find(",0,nan,0,nan,") != std::string::npos,
        "unknown transport became a zero shift");
    contract.previousEpochGpstSeconds = 1405216800;
    contract.previousCounter = 4;
    contract.previousDatumVersion = 7;
    contract.previousPublicSegment = "public-A";
    contract.jointCovarianceComplete = true;
    require(!zhangAuthorizePublishedIdentity(contract, 1405216830, "public-A", true, false, true, "proof-A"),
        "same segment label replaced full moment-map evidence");
    require(!zhangAuthorizePublishedIdentity(contract, 1405216830, "public-B", true, true, true, "proof-A"),
        "new public segment fabricated identity");
    require(!zhangAuthorizePublishedIdentity(contract, 1405216800, "public-A", true, true, true, "proof-A"),
        "same-epoch/stale endpoint certificate accepted");
    require(zhangAuthorizePublishedIdentity(contract, 1405216830, "public-A", true, true, true, "proof-A"),
        "certified aligned public chart rejected");
    require(contract.floatCorrectionShiftMetres == 0 && contract.integerTransportShiftCycles == 0,
        "native gamma was applied twice at the consumer");
    require(contract.floatTransportValid && contract.integerTransportValid, "aligned chart lost authorization");
    auto invalid = contract;
    invalid.jointCovarianceComplete = false;
    require(!zhangPublishedContinuityValid(invalid, 1405216830, "public-A"), "missing Qp accepted");
    invalid = contract;
    invalid.proofId = "UNAVAILABLE";
    require(!zhangPublishedContinuityValid(invalid, 1405216830, "public-A"), "unbound proof accepted");
    invalid = contract;
    invalid.floatCorrectionShiftMetres = 1;
    require(!zhangPublishedContinuityValid(invalid, 1405216830, "public-A"), "nonzero identity accepted");
    invalid = contract;
    invalid.kind = "UNKNOWN";
    require(!zhangPublishedContinuityValid(invalid, 1405216830, "public-A"), "contradictory UNKNOWN flags accepted");
    auto realOnly = contract;
    realOnly.kind = "CERTIFIED_AFFINE";
    realOnly.floatCorrectionShiftMetres = .125;
    realOnly.integerTransportValid = false;
    require(zhangPublishedContinuityValid(realOnly, 1405216830, "public-B"), "real chart map requires fictitious integer authorization");
    std::ostringstream realWire;
    zhangWritePublishedContinuityCsv(realWire, realOnly);
    require(realWire.str().find(",1,0.125,0,nan,") != std::string::npos,
        "real-only transport fabricated integer shift");
    require(!zhangPublishedProductHasNoiseParameters(false, true), "G10 phase without clock Qp accepted");
    require(!zhangPublishedProductHasNoiseParameters(true, false), "G10 clock without phase Qp accepted");
    require(zhangPublishedProductHasNoiseParameters(true, true), "PPP-only pair lost Qp coverage");
    require(zhangPublishedProductRequiresJointCovariance("PRODUCT_FIXED", true),
        "PPP-only G10 row was excluded from the published Qp parameter set");
    require(!zhangPublishedProductRequiresJointCovariance("PRODUCT_FIXED", false),
        "unusable row unnecessarily expanded the published Qp");
    require(!zhangPublishedProductRequiresJointCovariance("NETWORK_FIXED_DIAGNOSTIC", true),
        "diagnostic branch was presented as a formal fixed product");
    // Independent affine oracle, no call to solver/canonicalizer. A known
    // offset changes the mean, never covariance; linear mixing transports
    // every cross block, not just marginal sigmas.
    const double Q[3][3] = {{4, 1, 2}, {1, 9, 3}, {2, 3, 16}};
    const double T[3][3] = {{1, 0, 0}, {-1, 1, 0}, {0, 0, 1}};
    double transported[3][3] = {};
    for (int a=0; a<3; ++a) for (int b=0; b<3; ++b)
        for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
            transported[a][b] += T[a][i]*Q[i][j]*T[b][j];
    const double oracle[3][3] = {{4,-3,2},{-3,11,1},{2,1,16}};
    for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
        require(transported[i][j] == oracle[i][j], "complete affine covariance oracle mismatch");
    std::cout << "PASS user_published_continuity_contract checks=" << checks << '\n';
}

