#include "common/acsConfig.hpp"
#include "common/acsQC.hpp"
#include "common/observations.hpp"
#include "common/satStat.hpp"
#include <array>
#include <iostream>
#include <stdexcept>

void require(bool value, const char* reason)
{
    if (!value) throw std::runtime_error(reason);
}

int main()
{
    for (auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
        acsConfig.process_sys[system] = true;
    acsConfig.process_sys[E_Sys::GLO] = false;
    ObsList observations;
    auto excluded = std::make_shared<GObs>();
    excluded->Sat = SatSys("C06");
    excluded->excludeConfig = true;
    require(excluded->satStat_ptr == nullptr, "Expected unprepared observation");
    observations.push_back(excluded);
    auto unattached = std::make_shared<GObs>();
    unattached->Sat = SatSys("C19");
    observations.push_back(unattached);
    std::array<SatStat, 5> status;
    const std::array<const char*, 5> satellites{"G01", "E11", "C20", "C07", "R01"};
    for (size_t i = 0; i < status.size(); ++i)
    {
        auto obs = std::make_shared<GObs>();
        obs->Sat = SatSys(satellites[i]);
        obs->satStat_ptr = &status[i];
        obs->excludeConfig = i == 3;
        auto& sig = status[i].sigStatMap["test"];
        sig.slip.any = 3;
        sig.phaseQuarantined = true;
        sig.codeQuarantined = true;
        sig.phaseContinuityUnproven = true;
        sig.tracking.pendingBreak = i == 2;
        observations.push_back(obs);
    }
    clearSlips(observations);
    for (int i = 0; i < 3; ++i)
    {
        const auto& sig = status[i].sigStatMap.at("test");
        require(sig.slip.any == 0 && sig.previousPhaseQuarantined &&
                !sig.phaseQuarantined && !sig.codeQuarantined,
                "Active G/E/C per-epoch slip reset changed");
        require(sig.phaseContinuityUnproven == (i == 2),
                "Pending phase-break evidence was erased");
    }
    for (int i : {3, 4})
    {
        const auto& sig = status[i].sigStatMap.at("test");
        require(sig.slip.any == 3 && sig.phaseQuarantined && sig.codeQuarantined,
                "Excluded observation changed persistent state");
    }
    std::cout << "MULTI_GNSS_EXCLUDED_SLIPS_PASS excluded_untouched=1 active_reset_preserved=1\n";
}
