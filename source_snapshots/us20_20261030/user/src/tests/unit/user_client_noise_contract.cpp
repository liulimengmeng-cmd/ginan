#include "common/zhangUserNoiseAudit.hpp"
#include "common/zhangUserArExperiment.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>

// Linked against the real PEA objects with --wrap=main. This deliberately
// exercises KFMeas's actual keyed-noise constructor, not a reimplementation.
extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool ok, const char* label)
    {
        ++checks;
        if (!ok) throw std::runtime_error(std::string("FAIL ") + label);
    };
    try
    {
        KFState state;
        KFMeasEntryList collided, independent;
        constexpr double epsilon = 1e-6;
        for (int row = 0; row < 3; ++row)
        {
            KFMeasEntry oldEntry(&state), entry(&state);
            oldEntry.obsKey.type = entry.obsKey.type = KF::Z_AMB;
            entry.obsKey.num = row;
            oldEntry.addNoiseEntry(oldEntry.obsKey, 1, epsilon);
            entry.addNoiseEntry(entry.obsKey, 1, epsilon);
            collided.push_back(oldEntry);
            independent.push_back(entry);
        }
        KFMeas oldMeas(state, collided, GTime::noTime());
        KFMeas meas(state, independent, GTime::noTime());
        check(oldMeas.noiseIndexMap.size() == 1, "old keys collide");
        check((oldMeas.R - epsilon * MatrixXd::Ones(3, 3)).norm() == 0, "old epsilon 11 transpose");
        check(meas.noiseIndexMap.size() == 3, "new keys independent");
        check((meas.R - epsilon * MatrixXd::Identity(3, 3)).norm() == 0, "new epsilon I");

        // Independent expected covariance: diag(observation) plus one common
        // clock plus one common signed phase-bias source. SD cancels the clock
        // but doubles opposite-sign bias errors. No diagonal-only surrogate.
        KFMeasEntryList observations;
        KFKey clock, bias;
        clock.type = KF::SAT_CLOCK;
        bias.type = KF::PHASE_BIAS;
        for (int row = 0; row < 2; ++row)
        {
            KFMeasEntry entry(&state);
            entry.obsKey.type = KF::PHAS_MEAS;
            entry.obsKey.str = "TEST";
            entry.obsKey.num = 1;
            entry.obsKey.Sat = SatSys(row ? "G02" : "G01");
            entry.addNoiseEntry(entry.obsKey, 1, row ? 9e-6 : 4e-6);
            entry.addNoiseEntry(clock, 1, 0.04);
            entry.addNoiseEntry(bias, row ? -1 : 1, 0.0025);
            observations.push_back(entry);
        }
        KFMeas full(state, observations, GTime::noTime());
        MatrixXd expected(2, 2);
        expected << 0.042504, 0.0375, 0.0375, 0.042509;
        check((full.R - expected).cwiseAbs().maxCoeff() < 1e-15, "full correlated covariance");
        const auto beforeR = full.R;
        const auto beforeH = full.H_star;
        const auto parts = zhangUserNoiseParts(full, full.noiseIndexMap);
        check(parts.untracked.norm() < 1e-15, "source closure");
        VectorXd d(2);
        d << 1, -1;
        check(std::abs((d.transpose() * full.R * d)(0) - 0.010013) < 1e-15, "DRDt with correlation");
        check(std::abs((d.transpose() * parts.covariance[1] * d)(0)) < 1e-15, "common clock cancels");
        check(std::abs((d.transpose() * parts.covariance[2] * d)(0) - 0.01) < 1e-15, "opposite bias signs");
        check((full.R - beforeR).norm() == 0 && (full.H_star - beforeH).norm() == 0, "audit read only");
        std::cout << "PASS user_client_noise_contract checks=" << checks << '\n';
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
