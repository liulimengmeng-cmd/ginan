// Link against the real PEA objects (with main renamed), never a KF stub.
// --snapshot PATH records every complete production x/dx/P and QC workspace
// for an external byte-exact baseline/new comparison. No real GNSS data runs.
#include "common/algebra.hpp"
#include "common/zhangFloatThreading.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Float = ZhangFloatThreading;
namespace
{
int checks = 0, helperCases = 0, productionCases = 0;
int actualOpenmpTeam = 0, memoryWorkers = 0, parallelHelperCases = 0, productionMaximumDimension = 0;
std::ostringstream snapshots(std::ios::out | std::ios::binary);

void require(bool ok, const std::string& reason)
{
    ++checks;
    if (!ok) throw std::runtime_error(reason);
}

struct Environment
{
    std::string name;
    std::optional<std::string> old;
    Environment(const char* key, const char* value) : name(key)
    {
        if (const auto* setting = std::getenv(key)) old = setting;
        setenv(key, value, 1);
    }
    ~Environment()
    {
        if (old) setenv(name.c_str(), old->c_str(), 1);
        else unsetenv(name.c_str());
    }
};

struct Context
{
    const bool before = Float::inFloat;
    explicit Context(bool enabled) { Float::inFloat = enabled; }
    ~Context() { Float::inFloat = before; }
};

template<class A, class B> void exact(const A& a, const B& b, const char* reason)
{
    require(a.rows() == b.rows() && a.cols() == b.cols(), std::string(reason) + " dimensions");
    require(std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0, reason);
}

void appendString(const std::string& value)
{
    const std::uint64_t size = value.size();
    snapshots.write(reinterpret_cast<const char*>(&size), sizeof(size));
    snapshots.write(value.data(), value.size());
}

template<class A> void appendDense(const A& value)
{
    const std::int64_t rows = value.rows(), columns = value.cols();
    snapshots.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    snapshots.write(reinterpret_cast<const char*>(&columns), sizeof(columns));
    snapshots.write(reinterpret_cast<const char*>(value.data()), value.size() * sizeof(double));
}

void parallelRuntimeProbe()
{
    Context context(true);Environment policy("ZHANG_FLOAT_MEMORY_MODE", "parallel");
    require(Float::parallelMemoryWork(512, 512), "512x512 helper threshold is not enabled in FLOAT");
    memoryWorkers = Float::memoryThreads();
    require(memoryWorkers >= 2, "native execution must allow multiple physical memory workers");
#ifdef _OPENMP
    require(omp_get_max_threads() >= 2 && omp_get_thread_limit() >= 2,
        "native OpenMP execution is still pinned to one worker");
#pragma omp parallel num_threads(memoryWorkers)
    {
#pragma omp single
        actualOpenmpTeam = omp_get_num_threads();
    }
#endif
    require(actualOpenmpTeam >= 2, "actual OpenMP team did not execute with multiple workers");
}

#ifndef FLOAT_POSTERIOR_WRITEBACK_BASELINE
void helperContracts()
{
    const std::uint64_t special[] = {
        0x0000000000000000ULL, 0x8000000000000000ULL,
        0x0000000000000001ULL, 0x8000000000000001ULL,
        0x000fffffffffffffULL, 0x800fffffffffffffULL,
        0x0010000000000000ULL, 0x8010000000000000ULL,
        0x7fefffffffffffffULL, 0xffefffffffffffffULL,
        0x7ff0000000000000ULL, 0xfff0000000000000ULL,
        0x7ff8000000000123ULL, 0xfff8000000000456ULL,
        0x3ff0000000000000ULL, 0xc000000000000000ULL};
    for (bool floating : {false, true})
    for (const char* mode : {"serial", "parallel"})
    for (int n : {0, 1, 2, 7, 33, 511, 512, 513})
    for (bool padded : {false, true})
    for (bool exceptional : {false, true})
    {
        Environment policy("ZHANG_FLOAT_MEMORY_MODE", mode);
        Context context(floating);
        const int sourceLeading = n + (padded ? 5 : 0);
        const int targetLeading = n + (padded ? 9 : 0);
        const int offset = padded ? 11 : 0;
        std::vector<double> source(offset + sourceLeading * n + 13, 5.25);
        std::vector<double> target(offset + targetLeading * n + 17, -91.25);
        for (int column = 0; column < n; ++column)
        for (int row = 0; row < n; ++row)
        {
            double value = std::sin(row * 0.17 + column * 0.31) * (row + 1);
            if (exceptional)
            {
                const auto bits = special[(row * 5 + column * 7) % 16];
                std::memcpy(&value, &bits, sizeof(value));
            }
            source[offset + row + column * sourceLeading] = value;
        }
        const auto unchanged = source;
        auto expected = target;
        if (n)
        {
            using Stride = Eigen::OuterStride<Eigen::Dynamic>;
            Eigen::Map<const MatrixXd, 0, Stride> oldSource(
                source.data() + offset, n, n, Stride(sourceLeading));
            // Independent literal old expression: do not call the new helper
            // to construct its expected result. Include NaN payloads, signed
            // zero and diagonal overflow in the full byte comparison.
            const MatrixXd old = (oldSource + oldSource.transpose()).eval() / 2;
            for (int column = 0; column < n; ++column)
            for (int row = 0; row < n; ++row)
                expected[offset + row + column * targetLeading] = old(row, column);
        }
        const auto eigenThreads = Eigen::nbThreads();
#ifdef _OPENMP
        const auto ompThreads = omp_get_max_threads();
        const auto ompDynamic = omp_get_dynamic();
#endif
#ifdef GINAN_USE_OPENBLAS
        const auto blasThreads = BlasThreading::openblas_get_num_threads();
#endif
        Float::symmetricAverageCopy(source.data() + offset, sourceLeading,
            target.data() + offset, targetLeading, n);
        if (Float::parallelMemoryWork(n, n))
        {
            require(Float::memoryThreads() >= 2, "eligible helper call has no parallel worker budget");
            ++parallelHelperCases;
        }
        const std::string tag = std::string("helper n=") + std::to_string(n) +
            " mode=" + mode + " exceptional=" + std::to_string(exceptional);
        require(std::memcmp(expected.data(), target.data(), target.size() * sizeof(double)) == 0,
            tag + " differs from the old expression or overwrites padding");
        require(std::memcmp(source.data(), unchanged.data(), source.size() * sizeof(double)) == 0,
            tag + " mutates the source");
        require(Float::inFloat == floating && Eigen::nbThreads() == eigenThreads,
            tag + " changes FLOAT/Eigen context");
#ifdef _OPENMP
        require(omp_get_max_threads() == ompThreads && omp_get_dynamic() == ompDynamic,
            tag + " changes OpenMP policy");
#endif
#ifdef GINAN_USE_OPENBLAS
        require(BlasThreading::openblas_get_num_threads() == blasThreads,
            tag + " changes BLAS policy");
#endif
        ++helperCases;
    }
    for (const auto bits : special)
    {
        Context context(true);
        MatrixXd source(1, 1);std::memcpy(source.data(), &bits, sizeof(bits));
        const MatrixXd old = (source + source.transpose()).eval() / 2;
        MatrixXd changed(1, 1);
        Float::symmetricAverageCopy(source.data(), 1, changed.data(), 1, 1);
        exact(old, changed, "scalar signed-zero/subnormal/overflow/NaN diagonal differs");
        ++helperCases;
    }
}
#endif

struct Fixture
{
    KFState state;
    KFMeas meas;
    int begX = 0, numX = 7, begH = 0, numH = 5;
};

Fixture fixture(bool correlated, bool joseph, bool illConditioned, bool advanced, bool chunk, int dimension = 7)
{
    Fixture f;
    const int n = dimension, m = n == 7 ? 5 : n - 4;
    f.numX = n;f.numH = m;
    f.state.kfIndexMap.clear();f.state.stateTransitionMap.clear();
    for (int i = 0; i < n; ++i)
    {
        KFKey key;key.type = KF::IONO_STEC;key.str = "WRITEBACK";key.num = i;
        f.state.kfIndexMap[key] = i;f.state.stateTransitionMap[key][key][0] = 1;
    }
    f.state.initNoiseMap.clear();
    f.state.x = VectorXd::LinSpaced(n, -0.7, 1.1);
    f.state.dx = VectorXd::LinSpaced(n, 0.25, 0.5);
    MatrixXd a(n, n);
    for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)
        a(i, j) = 0.05 * std::sin(i * 3 + j * 7 + 0.5);
    f.state.P = a * a.transpose();f.state.P.diagonal().array() += 1.0;
    f.meas.H.resize(m, n);
    for (int i = 0; i < m; ++i)for (int j = 0; j < n; ++j)
        f.meas.H(i, j) = 0.13 * std::cos(i * 7 + j * 3 + 0.5);
    f.meas.R = MatrixXd::Identity(m, m) * 0.7;
    if (correlated)
    {
        MatrixXd b(m, m);
        for (int i = 0; i < m; ++i)for (int j = 0; j < m; ++j)
            b(i, j) = 0.07 * std::cos(i * 2 + j * 5 + 0.3);
        f.meas.R += b * b.transpose();
    }
    if (illConditioned)
    {
        require(n == 7, "small analytic ill-conditioned fixture dimension changed");
        f.state.P.setZero();f.meas.H.setZero();f.meas.R.setZero();
        const double variances[] = {1e-8, 1e-4, 1, 3, 1e4, 2, 7};
        const double noise[] = {1e-12, 1e-6, 0.5, 4, 1e3};
        for (int i = 0; i < n; ++i)f.state.P(i, i) = variances[i];
        for (int i = 0; i < m; ++i){f.meas.H(i, i) = 1;f.meas.R(i, i) = noise[i];}
    }
    f.meas.V = VectorXd::LinSpaced(m, -0.02, 0.03);
    f.meas.Y = f.meas.V;f.meas.VV = f.meas.V;
    if (advanced)
    {
        f.meas.H_star = MatrixXd::Identity(m, m);
        f.meas.uncorrelatedNoise = f.meas.R.diagonal();
    }
    else
    {
        f.meas.H_star = MatrixXd::Zero(m, 0);
        f.meas.uncorrelatedNoise = VectorXd::Zero(0);
    }
    f.meas.prefitRatios = f.meas.postfitRatios = VectorXd::Zero(m);
    f.meas.metaDataMaps.resize(m);f.meas.componentsMaps.resize(m);
    for (int i = 0; i < m; ++i)f.meas.obsKeys.push_back(f.state.kfIndexMap.begin()->first);
    f.state.joseph_stabilisation = joseph;f.state.advanced_postfits = advanced;
    f.state.prefitOpts.sigma_check = f.state.prefitOpts.omega_test = false;
    f.state.postfitOpts.sigma_check = f.state.postfitOpts.omega_test = false;
    f.state.postfitOpts.max_iterations = 2;
    f.state.chiSquareTest.enable = false;
    f.state.factorCommitSequence = 42;
    f.state.metaDataMap["WRITEBACK_SENTINEL"] = "before";
    if (chunk){f.begX = 1;f.numX = 5;f.begH = 1;f.numH = 3;}
    return f;
}

struct Output
{
    bool accepted = false;
    VectorXd x, dx, vv;
    MatrixXd p, qinv, qinvh;
};

Output solve(Fixture f, bool floating, bool omega = false)
{
    Context context(floating);Environment policy("ZHANG_R51_FLOAT_MATRIX_MODE", "legacy");
    f.state.postfitOpts.omega_test = omega;
    const auto beforeX = f.state.x;const auto beforeP = f.state.P;
    Output result;result.x = f.state.x;result.p = f.state.P;
    result.dx = VectorXd::Constant(f.state.x.size(), 0.125);
    result.qinv = MatrixXd::Constant(f.numH, f.numH, -17.0);
    result.qinvh = MatrixXd::Constant(f.numH, f.numX, -19.0);
    std::ostringstream trace;
    result.accepted = f.state.kFilter(trace, f.meas, result.x, result.p, result.dx,
        result.qinv, result.qinvh, f.begX, f.numX, f.begH, f.numH, true);
    result.vv = f.meas.VV;
    exact(beforeX, f.state.x, "kernel modifies prior x");
    exact(beforeP, f.state.P, "kernel modifies prior full P");
    return result;
}

void equalOutput(const Output& left, const Output& right)
{
    require(left.accepted == right.accepted, "writeback changes kernel acceptance");
    exact(left.x, right.x, "legacy/new complete x differs");
    exact(left.dx, right.dx, "legacy/new complete dx differs");
    exact(left.p, right.p, "legacy/new complete covariance differs");
    exact(left.qinv, right.qinv, "legacy/new omega Q inverse differs");
    exact(left.qinvh, right.qinvh, "legacy/new omega QinvH differs");
    exact(left.vv, right.vv, "legacy/new advanced postfits differ");
}

void record(const std::string& name, const Output& out)
{
    appendString(name);const std::uint8_t accepted = out.accepted;
    snapshots.write(reinterpret_cast<const char*>(&accepted), sizeof(accepted));
    appendDense(out.x);appendDense(out.dx);appendDense(out.p);
    appendDense(out.qinv);appendDense(out.qinvh);appendDense(out.vv);++productionCases;
    productionMaximumDimension = std::max(productionMaximumDimension, static_cast<int>(out.p.rows()));
}

void independentDenseOracle(const Fixture& f, const Output& actual)
{
    // Independent long-double dense equations, not the production LAPACK
    // chain, symmetric helper, or a re-call of kFilter as the only oracle.
    using LMatrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;
    using LVector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
    const LMatrix p = f.state.P.block(f.begX, f.begX, f.numX, f.numX).cast<long double>();
    const LMatrix h = f.meas.H.block(f.begH, f.begX, f.numH, f.numX).cast<long double>();
    const LMatrix r = f.meas.R.block(f.begH, f.begH, f.numH, f.numH).cast<long double>();
    const LMatrix innovation = h * p * h.transpose() + r;
    const Eigen::FullPivLU<LMatrix> decomposition(innovation);
    require(decomposition.rank() == f.numH, "independent fixture innovation is singular");
    const LMatrix k = decomposition.solve(LMatrix(h * p)).transpose();
    const LVector delta = k * f.meas.V.segment(f.begH, f.numH).cast<long double>();
    LMatrix posterior;
    if (f.state.joseph_stabilisation)
    {
        const LMatrix ikh = LMatrix::Identity(f.numX, f.numX) - k * h;
        posterior = ikh * p * ikh.transpose() + k * r * k.transpose();
    }
    else posterior = p - k * h * p;
    posterior = (posterior + posterior.transpose()).eval() / 2;
    VectorXd wantedX = f.state.x, wantedDx = VectorXd::Constant(f.state.x.size(), 0.125);
    MatrixXd wantedP = f.state.P;
    wantedDx.segment(f.begX, f.numX) = delta.cast<double>();
    wantedX.segment(f.begX, f.numX) += delta.cast<double>();
    wantedP.block(f.begX, f.begX, f.numX, f.numX) = posterior.cast<double>();
    require(actual.accepted, "well-defined independent fixture rejected");
    require((actual.x - wantedX).cwiseAbs().maxCoeff() <= 1e-12 * std::max(1.0, wantedX.cwiseAbs().maxCoeff()),
        "complete x differs from independent dense equations");
    require((actual.dx - wantedDx).cwiseAbs().maxCoeff() <= 1e-12 * std::max(1.0, wantedDx.cwiseAbs().maxCoeff()),
        "complete dx differs from independent dense equations");
    require((actual.p - wantedP).cwiseAbs().maxCoeff() <= 1e-12 * std::max(1.0, wantedP.cwiseAbs().maxCoeff()),
        "full P differs from independent dense equations");
}

void productionContracts()
{
    for (bool correlated : {false, true})for (bool joseph : {false, true})
    for (bool ill : {false, true})for (bool advanced : {false, true})for (bool chunk : {false, true})
    {
        const auto f = fixture(correlated, joseph, ill, advanced, chunk);
        for (bool omega : {false, true})
        {
            const auto old = solve(f, false, omega), changed = solve(f, true, omega);
            equalOutput(old, changed);independentDenseOracle(f, changed);
            record(std::string("dense_") + std::to_string(correlated) + std::to_string(joseph) +
                std::to_string(ill) + std::to_string(advanced) + std::to_string(chunk) + std::to_string(omega), changed);
        }
    }
    for (int defect = 0; defect < 4; ++defect)
    {
        auto f = fixture(true, false, false, false, false);
        if (defect == 0){f.state.P.setZero();f.meas.R.setZero();}
        if (defect == 1)f.meas.R(0, 0) = std::numeric_limits<double>::quiet_NaN();
        if (defect == 2)f.meas.R(0, 0) = std::numeric_limits<double>::infinity();
        if (defect == 3)
        {
            f.state.P = MatrixXd::Identity(7, 7) * (std::numeric_limits<double>::max() * 0.75);
            f.meas.H.setZero();f.meas.R.setIdentity();
        }
        const auto old = solve(f, false), changed = solve(f, true);
        equalOutput(old, changed);require(!changed.accepted, "singular/nonfinite/overflow fixture accepted");
        record("rejected_" + std::to_string(defect), changed);
    }
    // Cross the common 64-element BLAS/vectorisation scale with actual
    // production solves, correlated R, Joseph, omega and a 129-state advanced
    // postfit case. Complete output covariance/unchanged blocks remain in the
    // byte-exact old-path and external frozen-binary snapshots.
    for (int dimension : {65, 129})for (bool joseph : {false, true})for (bool omega : {false, true})
    {
        const auto f = fixture(true, joseph, false, dimension == 129, false, dimension);
        const auto old = solve(f, false, omega), changed = solve(f, true, omega);
        equalOutput(old, changed);independentDenseOracle(f, changed);
        record("large_correlated_" + std::to_string(dimension) + std::to_string(joseph) + std::to_string(omega), changed);
    }
    // Preserve the existing successful symmetric-indefinite fallback too.
    // This is a numerical-kernel fixture, not an acceptable physical noise
    // model or permission to weaken the caller's measurement-noise QC.
    auto indefinite = fixture(false, false, false, false, false);
    indefinite.state.P.setIdentity();indefinite.meas.H.setZero();
    indefinite.meas.R.setIdentity();indefinite.meas.R(0, 0) = -2;
    for (int i = 0; i < 5; ++i)indefinite.meas.H(i, i) = 1;
    const auto old = solve(indefinite, false), changed = solve(indefinite, true);
    equalOutput(old, changed);independentDenseOracle(indefinite, changed);
    record("indefinite_innovation_legacy_fallback", changed);
}

void rejectedTransactionContract()
{
    auto f = fixture(true, false, false, false, false);
    const auto entry = f.state;const auto measurement = f.meas;
    int packets = 0, committed = 0;
    f.state.authoritativeFactorTransactionCallback =
        [&](const KFState& state, const KFMeasurementTransaction& packet, std::string& reason)
        {
            ++packets;require(packet.qcConverged, "factor packet lacks ordinary passing QC");
            require(state.activeMeasurementTransaction == &packet, "factor callback escaped its provisional transaction");
            state.activeMeasurementTransaction->commitSideEffects.push_back([&]{++committed;});
            reason = "WRITEBACK_TEST_REJECTED";return false;
        };
    Context context(true);Environment policy("ZHANG_R51_FLOAT_MATRIX_MODE", "legacy");
    std::ostringstream trace;
    const auto result = f.state.filterKalman(trace, f.meas, "/PPP", true);
    require(result == KFFilterResult::FACTOR_TRANSACTION_REJECTED && packets == 1 && committed == 0,
        "rejected transaction commits a posterior or deferred side effect");
    exact(f.state.x, entry.x, "rejected transaction changes complete x");
    exact(f.state.dx, entry.dx, "rejected transaction changes complete dx");
    exact(f.state.P, entry.P, "rejected transaction changes full P");
    exact(f.meas.H, measurement.H, "rejected transaction changes H");
    exact(f.meas.R, measurement.R, "rejected transaction changes correlated R");
    exact(f.meas.VV, measurement.VV, "rejected transaction changes postfits");
    require(f.state.kfIndexMap == entry.kfIndexMap && f.state.stateTransitionMap == entry.stateTransitionMap &&
        f.state.factorCommitSequence == entry.factorCommitSequence && f.state.metaDataMap == entry.metaDataMap &&
        !f.state.lastPppTransactionCommitted && !f.state.lastPppQcConverged && !f.state.activeMeasurementTransaction,
        "rejected transaction changes ordering/metadata/sequence or retains authority");
    Output out;out.x = f.state.x;out.dx = f.state.dx;out.p = f.state.P;out.vv = f.meas.VV;
    record("factor_transaction_rejected", out);
}
}

int main(int argc, char** argv)
{
    try
    {
        std::optional<std::string> snapshot;
        if (argc == 3 && std::string(argv[1]) == "--snapshot") snapshot = argv[2];
        else require(argc == 1, "usage: float_posterior_writeback_contract [--snapshot NEW_PATH]");
        BlasThreading::ScopedOpenBlasThreadLimit deterministicBlas(1);
        parallelRuntimeProbe();
#ifndef FLOAT_POSTERIOR_WRITEBACK_BASELINE
        helperContracts();
#endif
        productionContracts();rejectedTransactionContract();
        const auto payload = snapshots.str();
        if (snapshot)
        {
            std::ifstream previous(*snapshot, std::ios::binary);
            require(!previous.good(), "preserve existing snapshot; do not overwrite");
            std::ofstream stream(*snapshot, std::ios::binary | std::ios::out);
            require(stream.good(), "cannot open new snapshot");
            stream.write(payload.data(), payload.size());stream.close();
            require(stream.good(), "complete snapshot write failed");
        }
        std::cout << "FLOAT_POSTERIOR_WRITEBACK_CONTRACT_PASS\n"
            << "{\"schema\":\"FLOAT_POSTERIOR_WRITEBACK_CONTRACT_V1\",\"passed\":true,\"checks\":" << checks
            << ",\"helper_cases\":" << helperCases << ",\"production_cases\":" << productionCases
            << ",\"actual_openmp_team\":" << actualOpenmpTeam << ",\"memory_workers\":" << memoryWorkers
            << ",\"parallel_helper_cases\":" << parallelHelperCases
            << ",\"production_max_dimension\":" << productionMaximumDimension
            << ",\"parallel_team_evidence\":\"RUNTIME_PROBE_WITH_SAME_HELPER_WORKER_POLICY\""
            << ",\"full_state_snapshot_bytes\":" << payload.size()
            << ",\"real_epoch_equivalence_claimed\":false}\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FLOAT_POSTERIOR_WRITEBACK_CONTRACT_FAILED: " << error.what() << '\n';
        return 1;
    }
}
