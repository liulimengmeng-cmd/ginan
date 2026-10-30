// Independent byte-loop oracles for production matrix movement helpers.
// No KF reconstruction or numerical solver is invoked by this standalone test.
#include "common/zhangFloatThreading.hpp"

#include <atomic>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Float = ZhangFloatThreading;
namespace
{
std::atomic<int> checks{0};

void require(bool pass, const char* reason)
{
    ++checks;
    if (!pass) throw std::runtime_error(reason);
}

std::vector<double> pattern(std::size_t size, std::uint64_t salt)
{
    std::vector<double> values(size);
    // -0, distinct quiet NaN payloads and infinities must survive movement
    // byte-for-byte. Never evaluate arithmetic on these copied values.
    const std::uint64_t special[] = {0, 0x8000000000000000ULL,
        0x7ff8000000000123ULL, 0xfff8000000000456ULL,
        0x7ff0000000000000ULL, 0xfff0000000000000ULL,
        0x3ff0000000000000ULL, 0xbff0000000000000ULL};
    for (std::size_t index = 0; index < size; ++index)
    {
        const std::uint64_t bits = index % 11 < 8 ? special[index % 11]
            : (0x3fe0000000000000ULL | ((index * 0x9e3779b97f4a7c15ULL + salt) & 0x000fffffffffffffULL));
        std::memcpy(values.data() + index, &bits, sizeof(bits));
    }
    return values;
}

void exact(const std::vector<double>& first, const std::vector<double>& second, const char* reason)
{
    require(first.size() == second.size(), "byte-comparison size mismatch");
    require(std::memcmp(first.data(), second.data(), first.size() * sizeof(double)) == 0, reason);
}

struct ThreadState
{
    bool inFloat;
    int eigen;
    int blas = -1;
#ifdef _OPENMP
    int ompMaximum, ompDynamic, ompLevels, activeLevel, insideParallel;
#endif

    static ThreadState capture()
    {
        ThreadState state{};
        state.inFloat = Float::inFloat;
        state.eigen = Eigen::nbThreads();
#ifdef GINAN_USE_OPENBLAS
        state.blas = BlasThreading::openblas_get_num_threads();
#endif
#ifdef _OPENMP
        state.ompMaximum = omp_get_max_threads();
        state.ompDynamic = omp_get_dynamic();
        state.ompLevels = omp_get_max_active_levels();
        state.activeLevel = omp_get_active_level();
        state.insideParallel = omp_in_parallel();
#endif
        return state;
    }

    void unchanged() const
    {
        const auto now = capture();
        require(now.inFloat == inFloat, "matrix helper changed FLOAT context");
        require(now.eigen == eigen, "matrix helper changed global Eigen threads");
        require(now.blas == blas, "matrix helper changed global OpenBLAS threads");
#ifdef _OPENMP
        require(now.ompMaximum == ompMaximum && now.ompDynamic == ompDynamic
            && now.ompLevels == ompLevels && now.activeLevel == activeLevel
            && now.insideParallel == insideParallel, "matrix helper leaked an OpenMP setting or active region");
#endif
    }
};

struct RestoreRuntime
{
    ThreadState before = ThreadState::capture();
    ~RestoreRuntime()
    {
        Float::inFloat = before.inFloat;
        Eigen::setNbThreads(before.eigen);
#ifdef GINAN_USE_OPENBLAS
        BlasThreading::openblas_set_num_threads(before.blas);
#endif
#ifdef _OPENMP
        omp_set_num_threads(before.ompMaximum);
        omp_set_dynamic(before.ompDynamic);
        omp_set_max_active_levels(before.ompLevels);
#endif
    }
};

struct RestoreMemoryMode
{
    std::optional<std::string> previous;
    RestoreMemoryMode()
    {
        if (const char* value = std::getenv("ZHANG_FLOAT_MEMORY_MODE")) previous = value;
    }
    ~RestoreMemoryMode()
    {
        if (previous) setenv("ZHANG_FLOAT_MEMORY_MODE", previous->c_str(), 1);
        else unsetenv("ZHANG_FLOAT_MEMORY_MODE");
    }
};

struct RestoreSetting
{
    std::string name;
    std::optional<std::string> previous;
    explicit RestoreSetting(const char* setting) : name(setting)
    {
        if (const char* value = std::getenv(setting)) previous = value;
    }
    ~RestoreSetting()
    {
        if (previous) setenv(name.c_str(), previous->c_str(), 1);
        else unsetenv(name.c_str());
    }
};

void independentMatrixThreadBudgetContracts()
{
    RestoreSetting restoreOuter("ZHANG_R51_FLOAT_THREADS");
    RestoreSetting restoreMatrix("ZHANG_FLOAT_MATRIX_THREADS");
    RestoreMemoryMode restoreMode;
    setenv("ZHANG_FLOAT_MEMORY_MODE", "parallel", 1);
    unsetenv("ZHANG_FLOAT_MATRIX_THREADS");
    const int available = Float::availableThreads();
    require(Float::requestedMatrixThreads(available) == available,
        "unset matrix setting changed inherited outer thread budget");
    require(Float::requestedMatrixThreads(0) == 0, "unset matrix setting changed legacy disabled override");
    setenv("ZHANG_FLOAT_MATRIX_THREADS", "all", 1);
    require(Float::requestedMatrixThreads(1) == available, "matrix all does not select available CPUs");
    setenv("ZHANG_FLOAT_MATRIX_THREADS", "1", 1);
    require(Float::requestedMatrixThreads(available) == 1, "explicit one-thread matrix budget ignored");
    setenv("ZHANG_R51_FLOAT_THREADS", "all", 1);
    const auto before = ThreadState::capture();
    {
        Float::ScopedFloatThreads scope;
        require(Float::inFloat && scope.threads() == available, "outer FLOAT all policy changed");
#ifdef _OPENMP
        require(omp_get_max_threads() == available,
            "matrix=1 incorrectly reduced the independent outer OpenMP team");
#endif
#ifdef GINAN_USE_OPENBLAS
        require(BlasThreading::openblas_get_num_threads() == 1 && Float::factorThreads() == 1,
            "matrix=1 did not constrain BLAS/factor workers");
        require(Eigen::nbThreads() == 1, "FLOAT Eigen pool no longer suppressed behind BLAS");
#endif
    }
    before.unchanged();
    bool injectedCaught = false;
    try
    {
        Float::ScopedFloatThreads scope;
        throw std::runtime_error("INJECTED_MATRIX_BUDGET_EXIT");
    }
    catch (const std::runtime_error& error)
    {
        injectedCaught = std::string(error.what()) == "INJECTED_MATRIX_BUDGET_EXIT";
    }
    require(injectedCaught, "exception fixture did not execute scoped matrix override");
    before.unchanged();
    for (const char* value : {"legacy", "", "0", "-1", "1x", "ALL"})
    {
        setenv("ZHANG_FLOAT_MATRIX_THREADS", value, 1);
        bool rejected = false;
        try { Float::ScopedFloatThreads scope; }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid independent matrix thread policy accepted");
        before.unchanged();
    }
    unsetenv("ZHANG_FLOAT_MATRIX_THREADS");
    {
        Float::ScopedFloatThreads scope;
#ifdef GINAN_USE_OPENBLAS
        require(BlasThreading::openblas_get_num_threads() == available,
            "unset matrix setting no longer preserves the original all-CPU BLAS policy");
#endif
    }
    before.unchanged();
}

void memoryModeContracts()
{
    RestoreMemoryMode restore;
    unsetenv("ZHANG_FLOAT_MEMORY_MODE");
    require(Float::parallelMemoryRequested(), "unset memory policy must select parallel eligibility");
    Float::inFloat = false;
    require(!Float::parallelMemoryWork(512, 512), "non-FLOAT work incorrectly eligible for new OpenMP helper");
    Float::inFloat = true;
    require(!Float::parallelMemoryWork(512, 511), "below-threshold work was parallelised");
    require(Float::parallelMemoryWork(512, 512), "exact-threshold work not eligible");
    require(Float::parallelMemoryWork(513, 513), "above-threshold work not eligible");
    require(!Float::parallelMemoryWork(0, 1000) && !Float::parallelMemoryWork(1000, 0),
        "zero-dimensional parallel eligibility");
    require(!Float::parallelMemoryWork(-1, 1000) && !Float::parallelMemoryWork(1000, -1),
        "negative dimensions admitted by helper predicate");
    require(Float::parallelMemoryWork(std::numeric_limits<Eigen::Index>::max(), 2),
        "eligibility arithmetic overflowed rows*columns");
    setenv("ZHANG_FLOAT_MEMORY_MODE", "serial", 1);
    require(!Float::parallelMemoryRequested() && !Float::parallelMemoryWork(513, 513),
        "serial control arm still admits parallel memory helpers");
    setenv("ZHANG_FLOAT_MEMORY_MODE", "parallel", 1);
    require(Float::parallelMemoryRequested(), "explicit parallel policy not selected");
    for (const char* value : {"", "0", "1", "on", "Serial", "parallel ", " serial"})
    {
        setenv("ZHANG_FLOAT_MEMORY_MODE", value, 1);
        bool rejected = false;
        try { (void)Float::parallelMemoryRequested(); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid memory policy silently accepted");
        const auto before = ThreadState::capture();
        rejected = false;
        try { Float::ScopedFloatThreads scope; }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "FLOAT scope did not validate an invalid memory mode");
        before.unchanged();
    }
}

void copyCase(Eigen::Index rows, Eigen::Index columns)
{
    const Eigen::Index sourceLeading = rows + 19, targetLeading = rows + 31;
    const Eigen::Index sourceOffset = 3 * sourceLeading + 5, targetOffset = 2 * targetLeading + 7;
    auto source = pattern(static_cast<std::size_t>((columns + 4) * sourceLeading + 43), 17);
    const auto originalSource = source;
    auto actual = pattern(static_cast<std::size_t>((columns + 4) * targetLeading + 47), 31);
    auto expected = actual;
    // Independent column-major scalar loop. Expected padding/prefix/suffix
    // remains exactly as initialised; compare entire allocations afterwards.
    for (Eigen::Index column = 0; column < columns; ++column)
        for (Eigen::Index row = 0; row < rows; ++row)
            std::memcpy(expected.data() + targetOffset + column * targetLeading + row,
                source.data() + sourceOffset + column * sourceLeading + row, sizeof(double));
    const auto threads = ThreadState::capture();
    Float::copyBlock(source.data() + sourceOffset, sourceLeading,
        actual.data() + targetOffset, targetLeading, rows, columns);
    threads.unchanged();
    exact(source, originalSource, "copyBlock changed source storage");
    exact(actual, expected, "copyBlock bytes/order/padding mismatch");
    const auto beforeNoOp = actual;
    Float::copyBlock(actual.data() + targetOffset, targetLeading,
        actual.data() + targetOffset, targetLeading, rows, columns);
    threads.unchanged();
    exact(actual, beforeNoOp, "copyBlock same-pointer/same-leading alias must be a no-op");
}

void transposeCase(Eigen::Index rows, Eigen::Index columns)
{
    const Eigen::Index sourceLeading = rows + 23, targetLeading = columns + 29;
    const Eigen::Index sourceOffset = 2 * sourceLeading + 7, targetOffset = 3 * targetLeading + 11;
    auto source = pattern(static_cast<std::size_t>((columns + 4) * sourceLeading + 47), 53);
    const auto originalSource = source;
    auto actual = pattern(static_cast<std::size_t>((rows + 5) * targetLeading + 59), 71);
    auto expected = actual;
    for (Eigen::Index column = 0; column < columns; ++column)
        for (Eigen::Index row = 0; row < rows; ++row)
            std::memcpy(expected.data() + targetOffset + row * targetLeading + column,
                source.data() + sourceOffset + column * sourceLeading + row, sizeof(double));
    const auto threads = ThreadState::capture();
    Float::transposeBlock(source.data() + sourceOffset, sourceLeading,
        actual.data() + targetOffset, targetLeading, rows, columns);
    threads.unchanged();
    exact(source, originalSource, "transposeBlock changed source storage");
    exact(actual, expected, "transposeBlock byte transpose/edge-tile/padding mismatch");
}

void mirrorCase(Eigen::Index size)
{
    const Eigen::Index leading = size + 17, offset = 2 * leading + 5;
    auto actual = pattern(static_cast<std::size_t>((size + 5) * leading + 29), 89);
    const auto original = actual;
    auto expected = actual;
    for (Eigen::Index column = 0; column < size; ++column)
        for (Eigen::Index row = column + 1; row < size; ++row)
            std::memcpy(expected.data() + offset + column * leading + row,
                original.data() + offset + row * leading + column, sizeof(double));
    const auto threads = ThreadState::capture();
    Float::mirrorUpper(actual.data() + offset, leading, size);
    threads.unchanged();
    exact(actual, expected, "mirrorUpper changed upper/diagonal/padding or failed lower byte copy");
}

void zeroContracts()
{
    const auto threads = ThreadState::capture();
    Float::copyBlock(nullptr, 0, nullptr, 0, 0, 0);
    Float::copyBlock(nullptr, 0, nullptr, 0, 0, 1000);
    Float::copyBlock(nullptr, 0, nullptr, 0, 1000, 0);
    Float::transposeBlock(nullptr, 0, nullptr, 0, 0, 0);
    Float::transposeBlock(nullptr, 0, nullptr, 0, 0, 1000);
    Float::transposeBlock(nullptr, 0, nullptr, 0, 1000, 0);
    Float::mirrorUpper(nullptr, 0, 0);
    threads.unchanged();
}

void scalarAndParallelContracts()
{
    // Below, exactly on and above 262144; both sides of each 32-wide tile.
    const std::pair<Eigen::Index, Eigen::Index> dimensions[] = {
        {0, 0}, {0, 7}, {7, 0}, {1, 1}, {1, 33}, {33, 1},
        {31, 33}, {32, 32}, {33, 31}, {512, 511}, {512, 512},
        {513, 513}, {257, 1031}, {1031, 257}};
    for (const bool inFloat : {false, true})
    {
        Float::inFloat = inFloat;
        zeroContracts();
        for (const auto& [rows, columns] : dimensions)
        {
            copyCase(rows, columns);
            transposeCase(rows, columns);
        }
        for (const Eigen::Index size : {0, 1, 31, 32, 33, 511, 512, 513}) mirrorCase(size);
    }
}

void existingOuterRegionContracts()
{
#ifdef _OPENMP
    const int workers = std::min(2, Float::availableThreads());
    std::vector<std::string> failures(workers);
    const auto caller = ThreadState::capture();
#pragma omp parallel num_threads(workers)
    {
        const int worker = omp_get_thread_num();
        const bool previous = Float::inFloat;
        // thread_local context must be explicit in a worker; this exercises
        // omp_in_parallel suppression even with inFloat=true on that worker.
        Float::inFloat = true;
        try
        {
            const auto before = ThreadState::capture();
            require(omp_in_parallel() || workers == 1, "outer OpenMP fixture did not enter a region");
            if (workers > 1)
                require(!Float::parallelMemoryWork(513, 513), "memory helper would create a nested OpenMP team");
            copyCase(513, 513);
            transposeCase(513, 517);
            mirrorCase(513);
            before.unchanged();
        }
        catch (const std::exception& error) { failures[worker] = error.what(); }
        Float::inFloat = previous;
    }
    caller.unchanged();
    for (const auto& failure : failures) require(failure.empty(), failure.c_str());
#endif
}
}

int main()
{
    try
    {
        RestoreRuntime restore;
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(std::min(4, Float::availableThreads()));
        // Allow nesting in the fixture: helpers must avoid starting it,
        // rather than relying on a process-level nesting prohibition.
        omp_set_max_active_levels(2);
#endif
        Eigen::setNbThreads(3);
#ifdef GINAN_USE_OPENBLAS
        BlasThreading::openblas_set_num_threads(2);
#endif
        RestoreMemoryMode restoreMode;
        independentMatrixThreadBudgetContracts();
        memoryModeContracts();
        for (const char* mode : {"serial", "parallel"})
        {
            setenv("ZHANG_FLOAT_MEMORY_MODE", mode, 1);
            scalarAndParallelContracts();
            existingOuterRegionContracts();
        }
        std::cout << "FLOAT_MATRIX_PARALLEL_HELPERS checks=" << checks.load()
            << " exact_bytes=true padded_blocks=true zero_and_edge_tiles=true outer_omp=true passed=true\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FLOAT_MATRIX_PARALLEL_HELPERS failed=" << error.what() << '\n';
        return 1;
    }
}
