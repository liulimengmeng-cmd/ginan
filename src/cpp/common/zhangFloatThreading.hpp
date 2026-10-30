#pragma once
#include <cstring>
#include "common/blasThreading.hpp"
#include "common/eigenIncluder.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <set>
#include <utility>
#ifdef __linux__
#include <sched.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ZhangFloatThreading
{
inline thread_local bool inFloat = false;

inline int availableThreads()
{
#ifdef _OPENMP
    return std::max(1, std::min(omp_get_num_procs(), omp_get_thread_limit()));
#else
    return 1;
#endif
}

// Count physical cores within the current Linux CPU affinity mask. SMT pairs
// help independent station work but repeatedly synchronised POTRF panels can
// be slower on both siblings. Do not assume every machine has SMT or divide
// its logical CPU count by two blindly.
inline int physicalCores()
{
#ifdef __linux__
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
    {
        std::set<std::pair<int, int>> cores;
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        {
            if (!CPU_ISSET(cpu, &mask)) continue;
            const std::string base = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/";
            int package, core;
            if (!(std::ifstream(base + "physical_package_id") >> package) ||
                !(std::ifstream(base + "core_id") >> core)) return availableThreads();
            cores.emplace(package, core);
        }
        if (!cores.empty()) return static_cast<int>(cores.size());
    }
#endif
    return availableThreads();
}

inline int factorThreads()
{
#ifdef GINAN_USE_OPENBLAS
    if (!inFloat) return 0;
#ifdef _OPENMP
    if (omp_in_parallel()) return 0; // Never change global BLAS control in a worker.
#endif
    return std::min(BlasThreading::openblas_get_num_threads(), physicalCores());
#else
    return 0;
#endif
}

// Unset/all selects every CPU visible to this process, not the inherited
// OMP_NUM_THREADS (which deliberately remains the AR-stage setting).
inline int requestedThreads(const char* setting, int available)
{
    if (!setting || std::string_view(setting) == "all") return available;
    if (std::string_view(setting) == "legacy") return 0;
    int threads = 0;
    const std::string_view text(setting);
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), threads);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
        threads < 1 || threads > available)
        throw std::invalid_argument("ZHANG_R51_FLOAT_THREADS must be all, legacy, or 1..available CPUs");
    return threads;
}

// A/B switch for bit-preserving memory kernels only. It does not change the
// BLAS arithmetic, filter model, or any acceptance criterion. Set before run.
inline bool parallelMemoryRequested()
{
    const char* value = std::getenv("ZHANG_FLOAT_MEMORY_MODE");
    if (!value || std::string_view(value) == "parallel") return true;
    if (std::string_view(value) == "serial") return false;
    throw std::invalid_argument("ZHANG_FLOAT_MEMORY_MODE must be parallel or serial");
}

// Matrix workers are independent of the outer receiver/copy team. More
// virtual CPUs need not give faster BLAS; leave the old setting unchanged
// unless a measured run profile explicitly selects a different count.
inline int requestedMatrixThreads(int outerThreads)
{
    const char* value = std::getenv("ZHANG_FLOAT_MATRIX_THREADS");
    if (!value) return outerThreads;
    const int threads = requestedThreads(value, availableThreads());
    if (threads == 0)
        throw std::invalid_argument("ZHANG_FLOAT_MATRIX_THREADS must be all or 1..available CPUs");
    return threads;
}

inline int memoryThreads()
{
    int threads = std::min(availableThreads(), physicalCores());
#ifdef _OPENMP
    threads = std::min(threads, omp_get_max_threads());
#endif
    return threads;
}

inline bool parallelMemoryWork(Eigen::Index rows, Eigen::Index columns)
{
    // Division avoids overflowing rows*columns for untrusted dimensions.
    return inFloat && rows > 0 && columns > 0 &&
        rows >= 262144 / columns + (262144 % columns != 0) && parallelMemoryRequested()
#ifdef _OPENMP
        && !omp_in_parallel()
#endif
        ;
}

// Non-overlapping, column-major blocks, including padded leading dimensions.
// Callers own valid storage; no arithmetic or floating-point reduction occurs.
inline void copyBlock(const double* source, Eigen::Index sourceLeading,
    double* destination, Eigen::Index destinationLeading,
    Eigen::Index rows, Eigen::Index columns)
{
    if (rows <= 0 || columns <= 0 ||
        (source == destination && sourceLeading == destinationLeading)) return;
    const bool parallel = parallelMemoryWork(rows, columns);
    const int threads = parallel ? memoryThreads() : 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
#endif
    for (Eigen::Index column = 0; column < columns; ++column)
        std::memcpy(destination + column * destinationLeading,
            source + column * sourceLeading, rows * sizeof(double));
}

inline void transposeBlock(const double* source, Eigen::Index sourceLeading,
    double* destination, Eigen::Index destinationLeading,
    Eigen::Index rows, Eigen::Index columns)
{
    if (rows <= 0 || columns <= 0) return;
    const bool parallel = parallelMemoryWork(rows, columns);
    const int threads = parallel ? memoryThreads() : 1;
    constexpr Eigen::Index tile = 32;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) num_threads(threads) if(parallel)
#endif
    for (Eigen::Index column = 0; column < columns; column += tile)
        for (Eigen::Index row = 0; row < rows; row += tile)
            for (Eigen::Index j = column; j < std::min(column + tile, columns); ++j)
                for (Eigen::Index i = row; i < std::min(row + tile, rows); ++i)
                    destination[j + i * destinationLeading] = source[i + j * sourceLeading];
}

// Preserve the computed upper triangle, including signed zeros and NaN bits.
// Each lower entry has exactly one owner and never aliases an upper read.
inline void mirrorUpper(double* data, Eigen::Index leading, Eigen::Index size)
{
    if (size <= 1) return;
    const bool parallel = parallelMemoryWork(size, size);
    const int threads = parallel ? memoryThreads() : 1;
    constexpr Eigen::Index tile = 32;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) num_threads(threads) if(parallel)
#endif
    for (Eigen::Index column = 0; column < size; column += tile)
        for (Eigen::Index row = 0; row < size; row += tile)
            for (Eigen::Index j = column; j < std::min(column + tile, size); ++j)
                for (Eigen::Index i = std::max(row, j + 1); i < std::min(row + tile, size); ++i)
                    data[i + j * leading] = data[j + i * leading];
}

// Fuse the legacy (source + source.transpose()).eval() / 2 writeback into
// non-overlapping, possibly padded column-major destination storage. Each
// output preserves its own original addition order, including the diagonal
// (a+a)/2: replacing that diagonal by a would change overflow/rejection cases.
// There is no matrix-sized intermediate, reduction, or BLAS thread change.
inline void symmetricAverageCopy(const double* source, Eigen::Index sourceLeading,
    double* destination, Eigen::Index destinationLeading, Eigen::Index size)
{
    if (size <= 0) return;
    const bool parallel = parallelMemoryWork(size, size);
    const int threads = parallel ? memoryThreads() : 1;
    constexpr Eigen::Index tile = 32;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) num_threads(threads) if(parallel)
#endif
    for (Eigen::Index column = 0; column < size; column += tile)
        for (Eigen::Index row = 0; row < size; row += tile)
            for (Eigen::Index j = column; j < std::min(column + tile, size); ++j)
                for (Eigen::Index i = row; i < std::min(row + tile, size); ++i)
                    destination[i + j * destinationLeading] =
                        (source[i + j * sourceLeading] + source[j + i * sourceLeading]) / 2.0;
}

// Set process-global OpenBLAS control only on the caller thread, outside any
// worker region. All exit paths restore the effective pre-FLOAT settings.
class ScopedFloatThreads
{
    const bool previousFloat = inFloat;
    const int previousEigen = Eigen::nbThreads();
#ifdef _OPENMP
    const int previousOmp = omp_get_max_threads();
    const int previousDynamic = omp_get_dynamic();
    const int previousLevels = omp_get_max_active_levels();
#endif
    const int requested = requestedThreads(std::getenv("ZHANG_R51_FLOAT_THREADS"), availableThreads());
    const int matrixRequested = requestedMatrixThreads(requested);
    BlasThreading::ScopedOpenBlasThreadLimit blas{matrixRequested};
public:
    ScopedFloatThreads()
    {
        (void)parallelMemoryRequested(); // Reject bad policy before changing state.
        inFloat = true;
        if (!requested) return;
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(requested);
#endif
        // Dense GEMM uses the one BLAS pool. Never run an Eigen OpenMP team
        // around a multi-threaded BLAS call.
#ifdef GINAN_USE_OPENBLAS
        Eigen::setNbThreads(1);
#else
        Eigen::setNbThreads(requested);
#endif
    }
    ~ScopedFloatThreads()
    {
        inFloat = previousFloat;
        if (!requested) return;
#ifdef _OPENMP
        omp_set_num_threads(previousOmp);
        omp_set_dynamic(previousDynamic);
        omp_set_max_active_levels(previousLevels);
#endif
        Eigen::setNbThreads(previousEigen);
    }
    int threads() const { return requested; }
    ScopedFloatThreads(const ScopedFloatThreads&) = delete;
    ScopedFloatThreads& operator=(const ScopedFloatThreads&) = delete;
};

// Outer station/measurement loops use OpenMP; each worker must not start a
// second BLAS or Eigen pool. Construct before entering, destroy after joining.
class ScopedOuterWork
{
    const int previousEigen = Eigen::nbThreads();
    BlasThreading::ScopedOpenBlasThreadLimit blas{1};
public:
    ScopedOuterWork() { Eigen::setNbThreads(1); }
    ~ScopedOuterWork() { Eigen::setNbThreads(previousEigen); }
    ScopedOuterWork(const ScopedOuterWork&) = delete;
    ScopedOuterWork& operator=(const ScopedOuterWork&) = delete;
};

// Independent contiguous chunks; bit-preserving copies, no BLAS workers.
template<class Plain> inline void copyDense(const Plain& source, Plain& destination)
{
    if(&source==&destination) return;
    destination.resize(source.rows(), source.cols());
    const Eigen::Index size = source.size(), chunk = 65536;
    int threads = std::min(availableThreads(), physicalCores());
#ifdef _OPENMP
    threads = std::min(threads,omp_get_max_threads());
#endif
    const bool parallel = inFloat && size >= 262144
#ifdef _OPENMP
        && !omp_in_parallel()
#endif
        ;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
#endif
    for (Eigen::Index begin = 0; begin < size; begin += chunk)
        std::memcpy(destination.data()+begin, source.data()+begin,
            std::min(chunk,size-begin)*sizeof(double));
}

class ScopedKernelContext
{
    bool previous = inFloat;
public:
    explicit ScopedKernelContext(bool enabled) { inFloat = enabled || previous; }
    ~ScopedKernelContext() { inFloat = previous; }
};

inline bool symmetricUpdateEnabled()
{
    const char* value = std::getenv("ZHANG_R51_FLOAT_MATRIX_MODE");
    if (!value || std::string(value) == "symmetric") return inFloat;
    if (std::string(value) == "legacy") return false;
    throw std::invalid_argument("ZHANG_R51_FLOAT_MATRIX_MODE must be symmetric or legacy");
}

class KernelTimer
{
    std::ostream& trace;
    const bool enabled = inFloat;
    const int states, measurements;
    std::chrono::steady_clock::time_point last = std::chrono::steady_clock::now();
    std::clock_t lastCpu = std::clock();
    std::ostringstream records;
public:
    KernelTimer(std::ostream& output, int n, int m) : trace(output), states(n), measurements(m) {}
    void mark(const char* phase)
    {
        if (!enabled) return;
        const auto now = std::chrono::steady_clock::now();
        const auto cpu = std::clock();
        const auto wallMs = std::chrono::duration<double, std::milli>(now - last).count();
        const double cpuMs = cpu != std::clock_t(-1) && lastCpu != std::clock_t(-1)
            ? 1000.0 * (cpu - lastCpu) / CLOCKS_PER_SEC : -1.0;
        records << "\nZHANG_FLOAT_KERNEL phase=" << phase << " states=" << states
            << " measurements=" << measurements << " elapsed_ms="
            << wallMs << " process_cpu_ms=" << cpuMs
            << " process_cpu_core_equivalent=" << (cpuMs >= 0 && wallMs > 0 ? cpuMs / wallMs : -1.0)
            << " cpu_scope=PROCESS_INCLUDING_OTHER_THREADS";
        last = now;
        lastCpu = cpu;
    }
    ~KernelTimer()
    {
        if (!enabled) return;
        static std::mutex outputMutex;
        const std::lock_guard<std::mutex> lock(outputMutex);
        trace << records.str();
    }
};
}
