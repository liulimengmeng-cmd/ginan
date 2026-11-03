#include "common/acsConfig.hpp"
#include "common/acsQC.hpp"
#include "common/constants.hpp"
#include "common/linearCombo.hpp"
#include "common/navigation.hpp"
#include "common/observations.hpp"
#include "common/satStat.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

void obs2lc(Trace&, GObs&, lc_t&);

static bool finite(const S_LC& lc)
{
    const double values[] = {lc.GF_Phas_m, lc.GF_Code_m, lc.WL_Phas_m, lc.WL_Phas_c,
        lc.WL_Code_m, lc.WL_Code_c, lc.NL_Phas_m, lc.NL_Phas_c, lc.NL_Code_m,
        lc.NL_Code_c, lc.IF_Phas_m, lc.IF_Code_m, lc.MW_m, lc.MW_c,
        lc.lam_A, lc.lam_B, lc.lam_WL, lc.lam_NL};
    for (const auto value : values) if (!std::isfinite(value)) return false;
    return true;
}

int main(int argc, char** argv)
{
    if (argc != 2 || (std::string(argv[1]) != "--expect-old-bug" &&
                      std::string(argv[1]) != "--expect-fixed"))
        throw std::runtime_error("Usage: --expect-old-bug|--expect-fixed");
    const bool old = std::string(argv[1]) == "--expect-old-bug";
    for (const auto system : {E_Sys::GAL, E_Sys::BDS})
    {
        acsConfig.code_priorities[system] = system == E_Sys::GAL
            ? std::vector<E_ObsCode>{E_ObsCode::L1C, E_ObsCode::L5Q}
            : std::vector<E_ObsCode>{E_ObsCode::L1P, E_ObsCode::L5P};
        SatNav navigation;
        navigation.lamMap[F1] = CLIGHT / 1575.42e6;
        navigation.lamMap[F5] = CLIGHT / 1176.45e6;
        SatStat status;
        GObs observation;
        observation.Sat = SatSys(system, system == E_Sys::GAL ? 2 : 19);
        observation.time = GTime(boost::posix_time::time_from_string("2024-07-17 00:00:30"));
        observation.satNav_ptr = &navigation;
        observation.satStat_ptr = &status;
        observation.sigs[F1].L = 21000000. / navigation.lamMap[F1];
        observation.sigs[F1].P = 21000010.;
        observation.sigs[F5].L = 21000001. / navigation.lamMap[F5];
        observation.sigs[F5].P = 21000012.;
        E_FType first, second, third;
        const bool selected = satFreqs(system, first, second, third);
        lc_t base;
        std::ostringstream trace;
        obs2lc(trace, observation, base);
        bool poisoned = false;
        for (const auto& [pair, lc] : base.lcMap)
        {
            const bool bad = !finite(lc);
            poisoned |= bad;
            std::cout << "ACTUAL_PRODUCTION_LC system=" << enum_to_string(system)
                << " pair=" << enum_to_string(pair.first) << ":" << enum_to_string(pair.second)
                << " valid=" << lc.valid << " finite_all_18=" << !bad
                << " wl_wavelength=" << lc.lam_WL << " if_phase=" << lc.IF_Phas_m
                << " mw_cycles=" << lc.MW_c << '\n';
        }
        std::cout << "ACTUAL_PRODUCTION_FREQS system=" << enum_to_string(system)
            << " selected=" << selected << " first=" << enum_to_string(first)
            << " second=" << enum_to_string(second) << " third=" << enum_to_string(third)
            << " poisoned_lc_present=" << poisoned << '\n';
        if (old != poisoned)
            throw std::runtime_error(old ? "OLD_DUPLICATE_FREQUENCY_BUG_NOT_REPRODUCED" : "FIXED_PRODUCTION_STILL_HAS_NONFINITE_LC");
    }
    std::cout << (old ? "CONFIRMED_OLD_PRODUCTION_DUPLICATE_F5_LC_BUG" : "PASS_FIXED_PRODUCTION_TWO_FREQUENCY_GAL_BDS_FINITE") << '\n';
}
