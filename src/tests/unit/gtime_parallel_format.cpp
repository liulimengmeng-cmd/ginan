// Production-linked shared-const time formatting contract. No solver is run.
#include "common/gTime.hpp"
#include <atomic>
#include <barrier>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace
{
std::uint64_t assertions = 0;
void require(bool value, const char* name)
{
    ++assertions;
    if (!value) throw std::runtime_error(name);
}

struct Known
{
    long double seconds;
    int precision;
    const char* expected;
};

// Constants and expected strings are independent of production formatters.
const Known known[] = {
    {0.L, 0, "1980-01-06 00:00:00"},
    {86400.L, 0, "1980-01-07 00:00:00"},
    {0.125L, 3, "1980-01-06 00:00:00.125"},
    {0.125L, 9, "1980-01-06 00:00:00.125000000"},
    {0.125L, 12, "1980-01-06 00:00:00.125000000000"},
    {1405209600.L, 0, "2024-07-17 00:00:00"},
    {1405209600.L, 9, "2024-07-17 00:00:00.000000000"},
    {1405209600.L, 12, "2024-07-17 00:00:00.000000000000"},
    {1405209600.125L, 9, "2024-07-17 00:00:00.125000000"},
    {1405209600.125L, 12, "2024-07-17 00:00:00.125000000000"},
    {1405216800.L, 0, "2024-07-17 02:00:00"},
    {1405209600.L, -4, "2024-07-17 00:00:00"},
    {1405209600.L, 40, "2024-07-17 00:00:00.000000000000"}
};

std::string iso(std::string value)
{
    value[10] = 'T';
    return value;
}

std::uint64_t serialSweep(bool dump = false)
{
    std::uint64_t hash = 14695981039346656037ULL;
    auto append = [&](const std::string& text)
    {
        if (dump) std::cout << text << '\n';
        for (unsigned char value : text)
        {
            hash ^= value;
            hash *= 1099511628211ULL;
        }
        hash ^= 255;
        hash *= 1099511628211ULL;
    };
    const long double bases[] = {0.L, 86399.L, 86400.L, 1261872000.L,
                                  1405209600.L, 1405216799.L};
    const int precisions[] = {-7, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 19};
    for (long double base : bases)
    for (int tick = 0; tick < 257; ++tick)
    {
        GTime time;
        time.bigTime = base + static_cast<long double>(tick) / 256;
        for (int precision : precisions)
        {
            append(time.to_string(precision));
            append(time.to_ISOstring(precision));
        }
    }
    return hash;
}
}

int main(int argc, char** argv)
{
    try
    {
        for (const auto& item : known)
        {
            GTime time;
            time.bigTime = item.seconds;
            require(time.to_string(item.precision) == item.expected, "known formatted value");
            require(time.to_ISOstring(item.precision) == iso(item.expected), "known ISO value");
        }
        const bool dump = argc == 2 && std::string(argv[1]) == "--dump-serial";
        const auto fingerprint = serialSweep(dump);
        if (dump || (argc == 2 && std::string(argv[1]) == "--serial-only"))
        {
            std::cout << "{\"passed\":true,\"known_assertions\":" << assertions
                      << ",\"serial_calls\":46260,\"serial_fingerprint\":\""
                      << std::hex << fingerprint << "\"}\n";
            return 0;
        }

        // Poison the ABI-preserved but now unused object cache. It must not be
        // consulted or changed, even if it pretends to be a matching entry.
        GTime shared;
        shared.bigTime = 1405209600.L;
        shared.cacheTime = shared.bigTime;
        shared.cacheN = 0;
        shared.cacheString = "DO_NOT_CONSUME_OBJECT_CACHE";
        const GTime& immutable = shared;
        require(immutable.to_string(0) == "2024-07-17 00:00:00", "poisoned cache ignored");
        shared.cacheN = 9;
        require(immutable.to_string(9) == "2024-07-17 00:00:00.000000000", "poisoned precision ignored");
        shared.cacheN = 0;

        GTime other;
        other.bigTime = 0.125L;
        const GTime& otherImmutable = other;
        constexpr int workers = 12;
        constexpr int iterations = 100000;
        std::barrier start(workers);
        std::atomic<std::uint64_t> failures{0};
        std::vector<std::thread> threads;
        for (int worker = 0; worker < workers; ++worker)
            threads.emplace_back([&, worker]()
            {
                start.arrive_and_wait();
                for (int iteration = 0; iteration < iterations; ++iteration)
                {
                    const int select = (iteration + worker) % 6;
                    const int precision = select % 3 == 0 ? 0 : select % 3 == 1 ? 9 : 12;
                    std::string expected = precision == 0 ? "2024-07-17 00:00:00" :
                        precision == 9 ? "2024-07-17 00:00:00.000000000" :
                                         "2024-07-17 00:00:00.000000000000";
                    const std::string result = select < 3 ? immutable.to_string(precision)
                                                         : immutable.to_ISOstring(precision);
                    if (select >= 3) expected[10] = 'T';
                    if (result != expected) ++failures;
                    // Switch the per-thread time key too, not only precision.
                    if (otherImmutable.to_string(9) != "1980-01-06 00:00:00.125000000")
                        ++failures;
                }
            });
        for (auto& thread : threads) thread.join();
        require(failures == 0, "shared const mixed-precision parallel format");
        require(shared.cacheTime == shared.bigTime && shared.cacheN == 0 &&
                shared.cacheString == "DO_NOT_CONSUME_OBJECT_CACHE", "legacy cache untouched");
        std::cout << "{\"passed\":true,\"assertions\":" << assertions
                  << ",\"parallel_calls\":" << 2ULL * workers * iterations
                  << ",\"workers\":" << workers << ",\"failures\":" << failures
                  << ",\"serial_calls\":46260,\"serial_fingerprint\":\""
                  << std::hex << fingerprint << "\"}\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GTIME_CONTRACT_FAILED: " << error.what() << '\n';
        return 1;
    }
}
