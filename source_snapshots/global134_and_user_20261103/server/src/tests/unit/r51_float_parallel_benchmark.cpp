// Link against the real PEA objects, not a second implementation of kFilter.
#include "common/algebra.hpp"
#include "common/blasThreading.hpp"
#ifdef FLOAT_THREAD_PROFILE
#include "common/zhangFloatThreading.hpp"
#endif
#include <omp.h>
#include <sys/resource.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <ctime>

static void require(bool pass, const char* reason)
{
    if (!pass) throw std::runtime_error(reason);
}

static void input(KFState& state, KFMeas& meas, int n, int m)
{
    state.x = VectorXd::Zero(n);
    state.dx = state.x;
    VectorXd u(n);
    for (int j = 0; j < n; ++j) u(j) = 0.02 * std::sin(j + 0.5);
    state.P = u * u.transpose();
    state.P.diagonal().array() += 1.0;
    meas.H = MatrixXd::Zero(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < 6; ++j)
            meas.H(i, (i * 17 + j * 109) % n) += std::cos(i * 0.01 + j) / 3;
    meas.R = MatrixXd::Identity(m, m) * 0.7;
    meas.V.resize(m);
    for (int i = 0; i < m; ++i) meas.V(i) = 0.01 * std::sin(i * 0.2);
    meas.Y = meas.V;
    meas.VV = meas.V;
    meas.H_star = MatrixXd::Zero(m, 0);
    meas.uncorrelatedNoise.resize(0);
    state.advanced_postfits = false;
    state.postfitOpts.omega_test = false;
}

#ifdef FLOAT_THREAD_PROFILE
class ScopedUnsetTestEnvironment
{
    const std::string name;
    const bool existed;
    const std::string previous;
public:
    explicit ScopedUnsetTestEnvironment(const char* variable)
        : name(variable), existed(std::getenv(variable) != nullptr),
          previous(existed ? std::getenv(variable) : "")
    {
        unsetenv(name.c_str());
    }
    ~ScopedUnsetTestEnvironment()
    {
        if (existed) setenv(name.c_str(), previous.c_str(), 1);
        else unsetenv(name.c_str());
    }
    ScopedUnsetTestEnvironment(const ScopedUnsetTestEnvironment&) = delete;
    ScopedUnsetTestEnvironment& operator=(const ScopedUnsetTestEnvironment&) = delete;
};

static void threadContract()
{
    using namespace ZhangFloatThreading;
    // Exercise the original inherited all-CPU contract independently of an
    // external benchmark-only BLAS budget, then restore that budget for BENCH.
    // The separate helper test covers explicit matrix=1 and exceptional exits.
    ScopedUnsetTestEnvironment matrixSetting("ZHANG_FLOAT_MATRIX_THREADS");
    const char* original = std::getenv("ZHANG_R51_FLOAT_THREADS");
    const bool hadSetting = original != nullptr;
    const std::string saved = original ? original : "";
    require(requestedThreads(nullptr, 12) == 12 && requestedThreads("all", 12) == 12 &&
            requestedThreads("legacy", 12) == 0 && requestedThreads("6", 12) == 6,
            "thread profile parsing failed");
    for (const char* invalid : {"", "0", "13", "-1", "2x", " 2", "ALL"})
    {
        bool rejected = false;
        try { requestedThreads(invalid, 12); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid thread profile accepted");
    }
    const int oldOmp = omp_get_max_threads();
    const int oldDynamic = omp_get_dynamic();
    const int oldLevels = omp_get_max_active_levels();
    const int oldEigen = Eigen::nbThreads();
    const int oldBlas = BlasThreading::openblas_get_num_threads();
    require(physicalCores() >= 1 && physicalCores() <= availableThreads(),
        "invalid physical core inventory");
    setenv("ZHANG_R51_FLOAT_THREADS", "all", 1);
    try
    {
        ScopedFloatThreads profile;
        require(omp_get_max_threads() == availableThreads() && Eigen::nbThreads() == 1 &&
            BlasThreading::openblas_get_num_threads() == availableThreads(),
            "all-CPU settings not applied");
        const int expectedFactor = std::min(availableThreads(), physicalCores());
        require(factorThreads() == expectedFactor, "factorisation thread cap mismatch");
        {
            BlasThreading::ScopedOpenBlasThreadLimit factor(factorThreads());
            require(BlasThreading::openblas_get_num_threads() == expectedFactor,
                "factorisation cap not applied");
        }
        require(BlasThreading::openblas_get_num_threads() == availableThreads(),
            "factorisation leaked its thread cap into GEMM");
        {
            ScopedOuterWork outer;
            int unsafeWorkers = 0;
#pragma omp parallel reduction(+:unsafeWorkers)
            {
                unsafeWorkers += BlasThreading::openblas_get_num_threads() != 1 ||
                    Eigen::nbThreads() != 1;
            }
            require(unsafeWorkers == 0, "nested BLAS/Eigen pools enabled in station workers");
        }
        require(BlasThreading::openblas_get_num_threads() == availableThreads(),
            "outer-work guard did not restore dense BLAS threads");
        throw std::runtime_error("injected FLOAT exit");
    }
    catch (const std::runtime_error& error)
    {
        require(std::string(error.what()) == "injected FLOAT exit", error.what());
    }
    require(omp_get_max_threads() == oldOmp && omp_get_dynamic() == oldDynamic &&
        omp_get_max_active_levels() == oldLevels && Eigen::nbThreads() == oldEigen &&
        BlasThreading::openblas_get_num_threads() == oldBlas && !inFloat,
        "FLOAT exit leaked thread settings into AR");
    if (hadSetting) setenv("ZHANG_R51_FLOAT_THREADS", saved.c_str(), 1);
    else unsetenv("ZHANG_R51_FLOAT_THREADS");
    std::cout << "PASS thread contract: parsing, actual all-CPU setting, no nested pools, exception restoration\n";
}

#ifndef FLOAT_FROZEN_PROFILE
static void snapshotContract()
{
    KFState state;KFMeas source;input(state,source,600,700);
    source.W=VectorXd::Ones(700)*3;
    source.prefitRatios=VectorXd::Ones(700)*4;
    source.postfitRatios=VectorXd::Ones(700)*5;
    source.H_star=MatrixXd::Constant(700,600,7);
    source.uncorrelatedNoise=VectorXd::Constant(600,8);
    KFKey key;key.str="SNAPSHOT";source.noiseIndexMap[key]=11;
    source.obsKeys={key};int marker=19;
    source.metaDataMaps={{{"MARKER",&marker}}};
    source.componentsMaps.resize(1);
    source.componentsMaps[0][E_Component::RANGE].value=9;
    const KFMeas copy=source.snapshotForFilter();
    require(copy.time==source.time && copy.noiseIndexMap==source.noiseIndexMap &&
        copy.obsKeys==source.obsKeys && copy.metaDataMaps==source.metaDataMaps &&
        copy.componentsMaps.size()==source.componentsMaps.size() &&
        copy.componentsMaps[0].at(E_Component::RANGE).value==9,"snapshot metadata lost");
    for(auto pair:{std::pair<const MatrixXd*,const MatrixXd*>{&source.H,&copy.H},
                  {&source.R,&copy.R},{&source.H_star,&copy.H_star}})
        require(pair.first->rows()==pair.second->rows()&&pair.first->cols()==pair.second->cols()&&
            std::memcmp(pair.first->data(),pair.second->data(),pair.first->size()*sizeof(double))==0,
            "snapshot matrix bytes changed");
    for(auto pair:{std::pair<const VectorXd*,const VectorXd*>{&source.Y,&copy.Y},
                  {&source.V,&copy.V},{&source.VV,&copy.VV},{&source.W,&copy.W},
                  {&source.uncorrelatedNoise,&copy.uncorrelatedNoise},
                  {&source.prefitRatios,&copy.prefitRatios},{&source.postfitRatios,&copy.postfitRatios}})
        require(pair.first->size()==pair.second->size()&&
            std::memcmp(pair.first->data(),pair.second->data(),pair.first->size()*sizeof(double))==0,
            "snapshot vector bytes changed");
    source.H.setZero();source.R.setZero();source.H_star.setZero();source.metaDataMaps.clear();
    require(copy.H.norm()>0&&copy.R.norm()>0&&copy.H_star(0,0)==7&&
        copy.metaDataMaps[0].at("MARKER")==&marker,"snapshot aliases mutable source storage");
    std::cout<<"PASS complete measurement snapshot: dense bytes, metadata, independent storage\n";
}
#endif
#endif

static void assemblyContract()
{
    KFState state;
    state.x = VectorXd::Zero(3);
    std::vector<KFKey> keys(3);
    for (int j = 0; j < 3; ++j)
    {
        keys[j].type = KF::IONO_STEC;
        keys[j].str = "ASSEMBLY";
        keys[j].num = j;
        state.kfIndexMap[keys[j]] = j;
    }
    KFKey noiseKey;
    noiseKey.type = KF::CODE_MEAS;
    noiseKey.str = "COMMON_NOISE";
    for (bool correlated : {false, true})
    {
        KFMeasEntryList entries;
        const int count = 129;
        std::vector<int> identity(count);
        for (int i = 0; i < count; ++i)
        {
            KFMeasEntry row;
            row.value = i * 0.5;
            row.innov = i * 0.25;
            row.noise = i + 1;
            row.designEntryMap[keys[i % 3]] = i + 2;
            row.obsKey = keys[i % 3];
            identity[i] = i;
            row.metaDataMap["ROW"] = &identity[i];
            if (correlated)
            {
                row.noiseElementMap[noiseKey] = 0.3;
                row.noiseEntryMap[noiseKey] = i % 2 ? -1 : 1;
            }
            entries.push_back(std::move(row));
        }
        const int previousEigen = Eigen::nbThreads();
        const int previousBlas = BlasThreading::openblas_get_num_threads();
        KFMeas meas(state, entries, GTime::noTime());
#ifdef FLOAT_THREAD_PROFILE
        require(Eigen::nbThreads() == previousEigen &&
            BlasThreading::openblas_get_num_threads() == previousBlas,
            "matrix assembly leaked thread settings");
#endif
        for (int i = 0; i < count; ++i)
        {
            require(meas.Y(i) == i * 0.5 && meas.V(i) == i * 0.25 &&
                meas.obsKeys[i] == keys[i % 3] &&
                meas.metaDataMaps[i].at("ROW") == &identity[i], "measurement order/value mismatch");
            for (int j = 0; j < 3; ++j)
                require(meas.H(i, j) == (j == i % 3 ? i + 2 : 0), "design row mismatch");
            for (int j = 0; j < count; ++j)
            {
                const double expected = correlated ? (i % 2 == j % 2 ? 0.3 : -0.3) :
                    (i == j ? i + 1 : 0);
                require(std::abs(meas.R(i, j) - expected) < 1e-14, "noise covariance mismatch");
            }
        }
    }
    std::cout << "PASS production matrix assembly: row order, design, diagonal/correlated noise\n";
}

static void oracle()
{
    for (bool joseph : {false, true})
    for (bool advanced : {false, true})
    for (bool omega : {false, true})
    {
        KFState state;
        KFMeas meas;
        input(state, meas, 7, 11);
        state.joseph_stabilisation = joseph;
        state.advanced_postfits = advanced;
        state.postfitOpts.omega_test = omega;
        if (advanced)
        {
            meas.H_star = MatrixXd::Identity(11, 11);
            meas.uncorrelatedNoise = meas.R.diagonal();
        }
        // Independent information-form oracle: no production Kalman core call.
        MatrixXd information = state.P.inverse() + meas.H.transpose() *
            meas.R.ldlt().solve(meas.H);
        MatrixXd expectedP = information.ldlt().solve(MatrixXd::Identity(7, 7));
        VectorXd expectedX = information.ldlt().solve(
            meas.H.transpose() * meas.R.ldlt().solve(meas.V));
        VectorXd xp = state.x, dx = state.dx;
        MatrixXd pp = state.P, qi, qh;
        std::ostringstream trace;
        require(state.kFilter(trace, meas, xp, pp, dx, qi, qh, 0, 7, 0, 11),
            "production kFilter failed");
        require((xp - expectedX).cwiseAbs().maxCoeff() < 2e-12 &&
                (pp - expectedP).cwiseAbs().maxCoeff() < 2e-12,
            "independent information-form posterior mismatch");
        if (omega)
        {
            MatrixXd q = meas.H * state.P * meas.H.transpose() + meas.R;
            require((q * qi - MatrixXd::Identity(11, 11)).norm() < 1e-11 &&
                    (qh - qi * meas.H).norm() < 1e-11, "omega matrix mismatch");
        }
        if (advanced)
        {
            MatrixXd q = meas.H * state.P * meas.H.transpose() + meas.R;
            require((meas.VV - q.ldlt().solve(meas.R) * meas.V).norm() < 1e-11,
                "advanced residual mismatch");
        }
    }
    std::cout << "PASS independent information-form oracle: 8 combinations\n";
}

static void edgeOracle()
{
    // Full cross-covariance outside a partial block must remain untouched.
    KFState state;KFMeas meas;input(state,meas,7,9);
    VectorXd xp=state.x,dx=state.dx;MatrixXd pp=state.P,qi,qh;
    const MatrixXd prior=state.P;
    const MatrixXd h=meas.H.block(2,1,4,3),p=prior.block(1,1,3,3);
    const MatrixXd s=h*p*h.transpose()+meas.R.block(2,2,4,4);
    const MatrixXd k=p*h.transpose()*s.ldlt().solve(MatrixXd::Identity(4,4));
    MatrixXd expected=prior;expected.block(1,1,3,3)=p-k*h*p;
    std::ostringstream trace;
    require(state.kFilter(trace,meas,xp,pp,dx,qi,qh,1,3,2,4),"partial block failed");
    require((pp-expected).cwiseAbs().maxCoeff()<2e-12,"partial block cross covariance changed");
    // Exercise the unchanged symmetric-indefinite fallback, independently.
    input(state,meas,3,3);meas.H.setIdentity();meas.R.setIdentity();meas.R(0,0)=-2;
    state.P.setIdentity();xp=state.x;pp=state.P;dx=state.dx;
    const MatrixXd indefinite=meas.H*state.P*meas.H.transpose()+meas.R;
    const MatrixXd gain=state.P*meas.H.transpose()*indefinite.inverse();
    const MatrixXd expectedFallback=state.P-gain*meas.H*state.P;
    require(state.kFilter(trace,meas,xp,pp,dx,qi,qh,0,3,0,3),"indefinite fallback failed");
    require((pp-expectedFallback).cwiseAbs().maxCoeff()<2e-12,"fallback posterior changed");
    meas.H.setZero();meas.R.setZero();xp=state.x;pp=state.P;
    require(!state.kFilter(trace,meas,xp,pp,dx,qi,qh,0,3,0,3),"singular innovation accepted");
    require((pp-state.P).isZero(),"failed solve did not restore covariance");
    std::cout<<"PASS partial-state cross covariance, indefinite fallback, singular rejection\n";
}

#if defined(FLOAT_THREAD_PROFILE) && !defined(FLOAT_FROZEN_PROFILE)
static void numericalFallbackContract()
{
    const char* previous=std::getenv("ZHANG_R51_FLOAT_MATRIX_MODE");
    const bool hadSetting=previous!=nullptr;const std::string saved=previous?previous:"";
    KFState state;KFMeas meas;input(state,meas,3,3);
    state.P=MatrixXd::Identity(3,3)*1e9;meas.H.setIdentity();
    meas.R=MatrixXd::Identity(3,3)*1e-5;
    VectorXd a=state.x,da=state.dx,b=state.x,db=state.dx;
    MatrixXd pa=state.P,pb=state.P,qi,qh;std::ostringstream legacy,guarded;
    setenv("ZHANG_R51_FLOAT_MATRIX_MODE","legacy",1);
    require(state.kFilter(legacy,meas,a,pa,da,qi,qh,0,3,0,3),"legacy cancellation fixture failed");
    setenv("ZHANG_R51_FLOAT_MATRIX_MODE","symmetric",1);
    require(state.kFilter(guarded,meas,b,pb,db,qi,qh,0,3,0,3),"guarded cancellation fixture failed");
    require((a-b).isZero()&&(pa-pb).isZero(),"cancellation fallback changed legacy moments");
    require(guarded.str().find("reason=LEGACY_CANCELLATION_RISK")!=std::string::npos,
        "unsafe symmetric update was not bypassed");
    if(hadSetting)setenv("ZHANG_R51_FLOAT_MATRIX_MODE",saved.c_str(),1);
    else unsetenv("ZHANG_R51_FLOAT_MATRIX_MODE");
    std::cout<<"PASS cancellation-risk dispatch: exact legacy moments, no statistical rejection\n";
}
#endif

int main(int argc, char** argv)
{
    const int threads = argc > 1 ? (std::string(argv[1]) == "all" ?
        omp_get_num_procs() : std::stoi(argv[1])) : 4;
    const int n = argc > 2 ? std::stoi(argv[2]) : 1400;
    const int m = argc > 3 ? std::stoi(argv[3]) : 2000;
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    Eigen::setNbThreads(0);
    BlasThreading::openblas_set_num_threads(4);
#ifdef FLOAT_THREAD_PROFILE
    threadContract();
    const std::string setting = threads == 0 ? "legacy" : std::to_string(threads);
    setenv("ZHANG_R51_FLOAT_THREADS", setting.c_str(), 1);
    ZhangFloatThreading::ScopedFloatThreads profile;
#ifndef FLOAT_FROZEN_PROFILE
    snapshotContract();
#endif
#else
    omp_set_num_threads(threads);
    BlasThreading::openblas_set_num_threads(threads);
#endif
    assemblyContract();
    oracle();
    edgeOracle();
#if defined(FLOAT_THREAD_PROFILE) && !defined(FLOAT_FROZEN_PROFILE)
    numericalFallbackContract();
#endif
    int team = 1;
#pragma omp parallel
    {
#pragma omp single
        team = omp_get_num_threads();
    }
    KFState state;
    KFMeas meas;
    input(state, meas, n, m);
    VectorXd xp = state.x, dx = state.dx;
    MatrixXd pp = state.P, qi, qh;
    std::ostringstream trace;
    const auto started = std::chrono::steady_clock::now();
    const auto cpuStarted = std::clock();
    require(state.kFilter(trace, meas, xp, pp, dx, qi, qh, 0, n, 0, m), "benchmark failed");
    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    const double cpu = double(std::clock() - cpuStarted) / CLOCKS_PER_SEC;
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout << std::setprecision(15) << "BENCH n=" << n << " m=" << m
        << " omp_team=" << team << " blas=" << BlasThreading::openblas_get_num_threads()
#ifdef FLOAT_THREAD_PROFILE
        << " factor_blas=" << ZhangFloatThreading::factorThreads()
#endif
        << " wall_s=" << wall << " cpu_s=" << cpu << " average_cores=" << cpu / wall
        << " peak_rss_kib=" << usage.ru_maxrss << " x_norm=" << xp.norm()
        << " p_norm=" << pp.norm() << '\n' << trace.str() << '\n';
    if (argc > 4)
    {
        std::ofstream stream(argv[4], std::ios::binary);
        stream.write(reinterpret_cast<const char*>(xp.data()), xp.size() * sizeof(double));
        stream.write(reinterpret_cast<const char*>(pp.data()), pp.size() * sizeof(double));
        require(bool(stream), "posterior output failed");
    }
}
