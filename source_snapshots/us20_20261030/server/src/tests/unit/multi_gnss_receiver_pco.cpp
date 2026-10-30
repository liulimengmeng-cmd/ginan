// Exercise the production receiver PCO lookup, including the legacy fallback.
#include "common/acsConfig.hpp"
#include "common/antenna.hpp"
#include "common/constants.hpp"
#include "common/navigation.hpp"
#include <iostream>
#include <stdexcept>

namespace
{
void check(bool pass, const char* reason)
{
    if (!pass) throw std::runtime_error(reason);
}
}

int main()
{
    const std::string id = "PCO_TEST_ONLY";
    const double ep[] = {2026, 10, 1, 0, 0, 0};
    const auto time = epoch2time(ep);
    const Vector3d first(0.01, -0.02, 0.10);
    const Vector3d second(-0.02, 0.01, 0.14);
    auto& low = nav.pcoMap[id][E_Sys::GPS][F1][time];
    auto& high = nav.pcoMap[id][E_Sys::GPS][F2][time];
    low.recPco = first;
    high.recPco = second;
    acsConfig.auto_fill_pco = true;
    // Independent wavelength-domain linear model evaluated at E5a/B2a.
    const double weight = (CLIGHT / FREQ5 - CLIGHT / FREQ1) /
                          (CLIGHT / FREQ2 - CLIGHT / FREQ1);
    const Vector3d expected = first + weight * (second - first);
    for (const auto sys : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        double variance = -1;
        const auto actual = antPco(id, sys, F5, time, variance, E_Radio::RECEIVER, true);
        check(actual.allFinite() && (actual - expected).norm() < 1e-12,
              "L2 wavelength is incorrect in receiver PCO interpolation");
        check(variance == 0, "Unexpected PCO variance contract");
    }
    double variance = -1;
    check((antPco(id, E_Sys::GPS, F1, time, variance, E_Radio::RECEIVER, true) - first).norm() == 0,
          "Exact GPS L1 calibration changed");
    check((antPco(id, E_Sys::GPS, F2, time, variance, E_Radio::RECEIVER, true) - second).norm() == 0,
          "Exact GPS L2 calibration changed");
    const Vector3d measured(0.03, -0.01, 0.12);
    nav.pcoMap[id][E_Sys::GAL][F5][time].recPco = measured;
    check((antPco(id, E_Sys::GAL, F5, time, variance, E_Radio::RECEIVER, true) - measured).norm() == 0,
          "Measured Galileo calibration was replaced by interpolation");
    nav.pcoMap.erase(id);
    std::cout << "MULTI_GNSS_RECEIVER_PCO_PASS exact_calibration_preserved=1 interpolation_verified=1\n";
}
