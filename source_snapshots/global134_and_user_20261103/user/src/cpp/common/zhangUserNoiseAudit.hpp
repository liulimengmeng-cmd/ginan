#pragma once

#include <array>
#include <cstdlib>
#include <iomanip>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include "common/algebra.hpp"

// Read-only, opt-in diagnostics. Bound by distinct epochs, not by function
// calls or receivers. The default performs no covariance products or logging.
inline bool zhangUserNoiseAuditEpoch(GTime time)
{
    static const unsigned limit = []()
    {
        const char* value = std::getenv("ZHANG_USER_NOISE_AUDIT_EPOCHS");
        if (!value) return 0u;
        const std::string text(value);
        if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("INVALID_USER_NOISE_AUDIT_EPOCH_BUDGET");
        const unsigned long parsed = std::stoul(text);
        if (parsed > 100) throw std::runtime_error("USER_NOISE_AUDIT_LIMIT_EXCEEDS_100_EPOCHS");
        return static_cast<unsigned>(parsed);
    }();
    if (!limit || time == GTime::noTime()) return false;
    static std::mutex mutex;
    static std::set<long double> epochs;
    const std::lock_guard<std::mutex> lock(mutex);
    if (epochs.count(time.bigTime)) return true;
    if (epochs.size() >= limit) return false;
    epochs.insert(time.bigTime);
    return true;
}

inline unsigned zhangUserNoiseCategory(KF type)
{
    if (type == KF::PHAS_MEAS || type == KF::CODE_MEAS) return 0;
    if (type == KF::SAT_CLOCK || type == KF::REC_CLOCK) return 1;
    if (type == KF::PHASE_BIAS || type == KF::CODE_BIAS) return 2;
    return 3; // Includes inseparable joint internal-product covariance factors.
}

struct ZhangUserNoiseParts
{
    std::array<MatrixXd, 4> covariance;
    MatrixXd untracked;
};

inline ZhangUserNoiseParts zhangUserNoiseParts(
    const KFMeas& meas, const std::map<KFKey, int>& originalNoiseKeys)
{
    ZhangUserNoiseParts out;
    for (auto& part : out.covariance)
        part = MatrixXd::Zero(meas.R.rows(), meas.R.cols());
    if (meas.H_star.rows() != meas.R.rows() ||
        meas.H_star.cols() != meas.uncorrelatedNoise.size())
        throw std::runtime_error("NOISE_AUDIT_DIMENSION_MISMATCH");
    for (const auto& [key, column] : originalNoiseKeys)
    {
        if (column < 0 || column >= meas.H_star.cols())
            throw std::runtime_error("NOISE_AUDIT_COLUMN_MISSING");
        const auto a = meas.H_star.col(column);
        out.covariance[zhangUserNoiseCategory(key.type)].noalias() +=
            meas.uncorrelatedNoise(column) * a * a.transpose();
    }
    out.untracked = meas.R;
    for (const auto& part : out.covariance) out.untracked -= part;
    return out;
}

inline void traceZhangUserMeasurementNoise(
    Trace& trace, const KFMeas& meas, const std::map<KFKey, int>& originalNoiseKeys,
    const char* stage)
{
    if (!zhangUserNoiseAuditEpoch(meas.time)) return;
    const auto savedPrecision = trace.precision();
    const auto savedFlags = trace.flags();
    trace << std::setprecision(17);
    const auto parts = zhangUserNoiseParts(meas, originalNoiseKeys);
    trace << "\nZHANG_USER_NOISE_AUDIT time=" << meas.time.to_string(0)
          << " stage=" << stage << " rows=" << meas.R.rows()
          << " sources=" << originalNoiseKeys.size()
          << " untracked_covariance_max_m2="
          << (parts.untracked.size() ? parts.untracked.cwiseAbs().maxCoeff() : 0)
          << " categories=OBSERVATION,CLOCK,BIAS,OTHER"
          << " joint_product_factors_category=OTHER statistic=NIS_NOT_COMPUTED";
    for (int row = 0; row < meas.R.rows(); ++row)
    {
        const auto& obs = meas.obsKeys.at(row);
        trace << "\nZHANG_USER_NOISE_ROW time=" << meas.time.to_string(0)
              << " stage=" << stage << " row=" << row << " obs=" << std::quoted(std::string(obs))
              << " R_total_m2=" << meas.R(row, row)
              << " R_observation_m2=" << parts.covariance[0](row, row)
              << " R_clock_m2=" << parts.covariance[1](row, row)
              << " R_bias_m2=" << parts.covariance[2](row, row)
              << " R_other_m2=" << parts.covariance[3](row, row)
              << " R_untracked_or_QC_delta_m2=" << parts.untracked(row, row);
        for (const auto& [key, column] : originalNoiseKeys)
        {
            const double coefficient = meas.H_star(row, column);
            if (coefficient == 0) continue;
            trace << "\nZHANG_USER_NOISE_SOURCE time=" << meas.time.to_string(0)
                  << " stage=" << stage << " row=" << row
                  << " key=" << std::quoted(std::string(key))
                  << " noise_column=" << column << " category=" << zhangUserNoiseCategory(key.type)
                  << " coefficient=" << coefficient
                  << " variance_native_unit2=" << meas.uncorrelatedNoise(column)
                  << " contribution_m2=" << coefficient * coefficient * meas.uncorrelatedNoise(column);
        }
    }
    // Measurement-level diagnostic D=[1,-1]. This is NOT the AR search's
    // posterior ambiguity covariance. Keep the common-clock/bias cross terms.
    std::map<std::tuple<std::string, int, int>, int> references;
    for (int row = 0; row < meas.R.rows(); ++row)
    {
        const auto& obs = meas.obsKeys.at(row);
        if (obs.type != KF::PHAS_MEAS && obs.type != KF::CODE_MEAS) continue;
        const auto group = std::make_tuple(obs.str, static_cast<int>(obs.type), obs.num);
        const auto [it, inserted] = references.emplace(group, row);
        if (inserted) continue;
        const int ref = it->second;
        trace << "\nZHANG_USER_NOISE_SD time=" << meas.time.to_string(0)
              << " stage=" << stage << " row=" << row << " reference_row=" << ref
              << " reference_policy=FIRST_DIAGNOSTIC_ROW_NOT_AR_PIVOT"
              << " R_SD_m2=" << meas.R(row, row) + meas.R(ref, ref) - meas.R(row, ref) - meas.R(ref, row)
              << " cross_covariance_m2=" << meas.R(row, ref);
        for (unsigned category = 0; category < 4; ++category)
        {
            const auto& part = parts.covariance[category];
            trace << " category_" << category << "_SD_m2="
                  << part(row, row) + part(ref, ref) - part(row, ref) - part(ref, row);
        }
    }
    trace.precision(savedPrecision);
    trace.flags(savedFlags);
}
