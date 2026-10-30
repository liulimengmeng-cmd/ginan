#include "common/zhangSystemClockDatum.hpp"

#include <Eigen/Dense>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

using Identity = std::tuple<int, int, int, std::string, int>;
Identity identity(const KFKey& key)
{
    return {static_cast<int>(key.type), static_cast<int>(key.Sat.sys),
            key.Sat.prn, key.str, key.num};
}

void enable(ACSConfig& config, E_Sys system)
{
    config.process_sys[system] = true;
    config.zhangFullRank.sysOpts[system].baseline_observables =
        {E_ObsCode::L1C, E_ObsCode::L2W};
    config.zhangPppAr.baseline_observables[system] =
        {E_ObsCode::L1C, E_ObsCode::L2W};
}

void keyContracts()
{
    ACSConfig config;
    enable(config, E_Sys::GPS);
    auto pointer = reinterpret_cast<Receiver*>(0x1234);
    require(!zhangUsesSystemReceiverClocks(config), "traditional mode changed");
    config.zhangFullRank.enable = true;
    auto gps = zhangReceiverClockKey(config, "A", E_Sys::GPS, 2,
                                     KF::REC_CLOCK_RATE, pointer);
    require(gps.Sat.sys == E_Sys::NONE && gps.Sat.prn == 0,
            "GPS-only checkpoint key changed");
    require(gps.rec_ptr == pointer && gps.num == 2 &&
            gps.type == KF::REC_CLOCK_RATE && gps.str == "A",
            "clock component or receiver metadata lost");
    require(zhangReceiverClockMatchesSystem(config, gps, E_Sys::GPS),
            "legacy GPS clock not matched");
    require(!zhangReceiverClockMatchesSystem(config, gps, E_Sys::GAL),
            "legacy GPS clock matched absent constellation");

    enable(config, E_Sys::GAL);
    enable(config, E_Sys::BDS);
    require(zhangUsesSystemReceiverClocks(config), "joint clock mode not selected");
    std::set<Identity> keys;
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        auto key = zhangReceiverClockKey(config, "A", system, 2,
                                         KF::REC_CLOCK_RATE, pointer);
        keys.insert(identity(key));
        require(key.Sat.sys == system && key.rec_ptr == pointer,
                "joint effective clock has wrong system or receiver");
        require(zhangReceiverClockMatchesSystem(config, key, system),
                "joint clock cannot find its own system");
        require(!zhangReceiverClockMatchesSystem(config, gps, system),
                "legacy NONE clock silently imported into joint model");
    }
    require(keys.size() == 3, "G/E/C effective clocks collapsed");
    KFKey satellite;
    satellite.type = KF::SAT_CLOCK;
    satellite.Sat = SatSys(E_Sys::GAL, 11);
    auto gPivot = zhangReceiverClockKey(config, "A", E_Sys::GPS);
    auto ePivot = zhangReceiverClockKey(config, "A", E_Sys::GAL);
    require(!zhangClockDatumPivotApplies(config, satellite, gPivot),
            "GPS reference subtracted from Galileo clock");
    require(zhangClockDatumPivotApplies(config, satellite, ePivot),
            "same-system reference excluded");

    config.zhangFullRank.enable = false;
    config.zhangPppAr.user_adapter = true;
    require(zhangUsesSystemReceiverClocks(config), "held-out joint user not scoped");
    config.zhangPppAr.user_adapter = false;
    require(!zhangUsesSystemReceiverClocks(config), "non-Zhang Ginan keys changed");
}

void independentRankModel()
{
    // Independent observation construction: three constellations, three
    // receivers, eight satellites/system. Receiver R0 fixes three independent
    // clock gauges. XYZ/ZTD of R1/R2 are shared by all their constellations.
    constexpr int receivers = 3, perSystem = 8, satellites = 24;
    constexpr int common = 8, clocks = 6;
    Eigen::MatrixXd gauged = Eigen::MatrixXd::Zero(72, common + clocks + satellites);
    Eigen::MatrixXd free = Eigen::MatrixXd::Zero(72, common + 9 + satellites);
    Eigen::MatrixXd collapsed = Eigen::MatrixXd::Zero(72, common + 2 + satellites);
    for (int receiver = 0; receiver < receivers; ++receiver)
        for (int sat = 0; sat < satellites; ++sat)
        {
            const int row = receiver * satellites + sat;
            const int system = sat / perSystem;
            if (receiver)
                for (int component = 0; component < 4; ++component)
                {
                    const double angle = 0.17 * (sat + 1) * (component + 1) +
                                         0.13 * receiver * (sat + 3);
                    const double value = component == 3
                        ? 1.4 + 0.2 * std::cos(angle)
                        : std::sin(angle);
                    const int column = (receiver - 1) * 4 + component;
                    gauged(row, column) = free(row, column) =
                        collapsed(row, column) = value;
                }
            free(row, common + receiver * 3 + system) = 1;
            free(row, common + 9 + sat) = -1;
            gauged(row, common + clocks + sat) = -1;
            collapsed(row, common + 2 + sat) = -1;
            if (receiver)
            {
                gauged(row, common + (receiver - 1) * 3 + system) = 1;
                collapsed(row, common + receiver - 1) = 1;
            }
        }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> freeQr(free), gaugeQr(gauged);
    freeQr.setThreshold(1e-10);
    gaugeQr.setThreshold(1e-10);
    require(free.cols() - freeQr.rank() == 3,
            "independent joint model does not have exactly three clock gauges");
    require(gaugeQr.rank() == gauged.cols(),
            "per-system reference clocks do not restore full column rank");
    Eigen::VectorXd truth = Eigen::VectorXd::Zero(gauged.cols());
    truth.segment(common, clocks) << 2, 12, -7, -1, 4, 9;
    const Eigen::VectorXd observation = gauged * truth;
    const double independentResidual = (observation - gauged * gaugeQr.solve(observation)).norm();
    const double collapsedResidual =
        (observation - collapsed * collapsed.colPivHouseholderQr().solve(observation)).norm();
    require(independentResidual < 1e-10, "independent clock observations not representable");
    require(collapsedResidual > 1,
            "one receiver clock incorrectly represents nonzero inter-system offsets");
    std::cout << "rank independent=38 nullity=3 collapsed_residual="
              << collapsedResidual << '\n';
}

void independentGaugeScope()
{
    ACSConfig config;
    config.zhangFullRank.enable = true;
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS}) enable(config, system);
    // [common XYZ/ZTD scalar, G receiver, G satellite, E receiver,
    //  E satellite, C receiver, C satellite]. Re-reference E only.
    std::vector<KFKey> keys(7);
    keys[0].type = KF::TROP;
    for (int system = 0; system < 3; ++system)
    {
        const E_Sys sys = std::array<E_Sys, 3>{E_Sys::GPS, E_Sys::GAL, E_Sys::BDS}[system];
        keys[1 + 2 * system] = zhangReceiverClockKey(config, "A", sys);
        keys[2 + 2 * system].type = KF::SAT_CLOCK;
        keys[2 + 2 * system].Sat = SatSys(sys, 10);
    }
    Eigen::MatrixXd transport = Eigen::MatrixXd::Identity(7, 7);
    const KFKey pivot = keys[3];
    for (int row = 1; row < 7; ++row)
        if (zhangClockDatumPivotApplies(config, keys[row], pivot))
            transport(row, 3) -= 1;
    Eigen::VectorXd mean(7);
    mean << 1, 2, 3, 5, 8, 13, 21;
    Eigen::MatrixXd factor(7, 7);
    for (int row = 0; row < 7; ++row)
        for (int column = 0; column < 7; ++column)
            factor(row, column) = (row == column ? 3.0 : 0.0) +
                                  0.03 * (row + 1) * (column + 2);
    const Eigen::MatrixXd covariance = factor * factor.transpose();
    const Eigen::VectorXd actualMean = transport * mean;
    const Eigen::MatrixXd actualCovariance = transport * covariance * transport.transpose();
    Eigen::VectorXd expectedMean = mean;
    expectedMean[3] = 0;
    expectedMean[4] = mean[4] - mean[3];
    require((actualMean - expectedMean).norm() < 1e-13,
            "gauge scope changes another constellation or common state");
    // Independently expand Cov(X_i-a_i X_E, X_j-a_j X_E) rather than
    // comparing a second call to the transformation under test.
    for (int row = 0; row < 7; ++row)
        for (int column = 0; column < 7; ++column)
        {
            const int a = row == 3 || row == 4;
            const int b = column == 3 || column == 4;
            const double expected = covariance(row, column) -
                a * covariance(3, column) - b * covariance(row, 3) +
                a * b * covariance(3, 3);
            require(std::abs(actualCovariance(row, column) - expected) < 1e-12,
                    "full cross-covariance not transported with scoped gauge");
        }
}

void independentJointNoise()
{
    // Explicitly nonzero cross-G/E/C Q. Reference entries are scalar dot
    // products, independent of the key aggregation implementation below.
    const double rows[3][2] = {{1, 0.5}, {-0.4, 2}, {0.6, 0.3}};
    std::map<Identity, std::array<double, 3>> randomSources;
    for (int measurement = 0; measurement < 3; ++measurement)
        for (int factor = 0; factor < 2; ++factor)
        {
            const auto key = zhangProductJointNoiseKey(12345, "PRODUCT_FIXED", factor);
            require(key.Sat.sys == E_Sys::NONE, "joint Q split by constellation");
            randomSources[identity(key)][measurement] = rows[measurement][factor];
        }
    Eigen::Matrix3d actual = Eigen::Matrix3d::Zero();
    for (const auto& [key, coefficients] : randomSources)
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 3; ++column)
                actual(row, column) += coefficients[row] * coefficients[column];
    const double expected[3][3] = {{1.25, 0.6, 0.75},
                                  {0.6, 4.16, 0.36},
                                  {0.75, 0.36, 0.45}};
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            require(std::abs(actual(row, column) - expected[row][column]) < 1e-13,
                    "cross-constellation product covariance lost");
    require(identity(zhangProductJointNoiseKey(12345, "PRODUCT_FIXED", 0)) !=
            identity(zhangProductJointNoiseKey(12346, "PRODUCT_FIXED", 0)),
            "different epochs share product noise source");
    require(identity(zhangProductJointNoiseKey(12345, "PRODUCT_FIXED", 0)) !=
            identity(zhangProductJointNoiseKey(12345, "FLOAT", 0)),
            "different product solutions share noise source");
}
}

int main()
{
    try
    {
        keyContracts();
        independentRankModel();
        independentGaugeScope();
        independentJointNoise();
        std::cout << "multi_gnss_clock_contract: PASS\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "multi_gnss_clock_contract: FAIL " << error.what() << '\n';
        return 1;
    }
}
