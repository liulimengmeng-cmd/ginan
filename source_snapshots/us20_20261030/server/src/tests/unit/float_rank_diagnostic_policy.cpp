#include "common/zhangFloatRankDiagnosticPolicy.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Policy = ZhangFloatRankDiagnostic;

namespace
{
int checks = 0;

void require(bool condition, const char* reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}

void environment(const char* value)
{
#ifdef _WIN32
    require(_putenv_s(Policy::environmentVariable, value ? value : "") == 0,
        "cannot set diagnostic test environment");
#else
    require((value ? setenv(Policy::environmentVariable, value, 1)
                   : unsetenv(Policy::environmentVariable)) == 0,
        "cannot set diagnostic test environment");
#endif
}

struct RestoreEnvironment
{
    std::optional<std::string> value;
    RestoreEnvironment()
    {
        if (const char* current = std::getenv(Policy::environmentVariable)) value = current;
    }
    ~RestoreEnvironment()
    {
#ifdef _WIN32
        _putenv_s(Policy::environmentVariable, value ? value->c_str() : "");
#else
        if (value) setenv(Policy::environmentVariable, value->c_str(), 1);
        else unsetenv(Policy::environmentVariable);
#endif
    }
};

void parseContracts()
{
    require(std::string_view(Policy::environmentVariable) == "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC",
        "gate is not narrowly scoped to the pure-observation diagnostic");
    require(Policy::parseMode({}) == Policy::Mode::Off, "unset/empty default must be off");
    require(Policy::parseMode("off") == Policy::Mode::Off, "off parse");
    require(Policy::parseMode("bounded") == Policy::Mode::Bounded, "bounded parse");
    require(Policy::parseMode("full") == Policy::Mode::Full, "full parse");
    for (const auto value : {"on", "1", "true", "FALSE", "OFF", "Full", "full ", " full", "\n", "bounded\n"})
    {
        bool rejected = false;
        try { (void)Policy::parseMode(value); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "ambiguous/invalid diagnostic mode was silently accepted");
    }
    bool embeddedNullRejected = false;
    try { (void)Policy::parseMode(std::string_view("full\0off", 8)); }
    catch (const std::invalid_argument&) { embeddedNullRejected = true; }
    require(embeddedNullRejected, "parser ignores string_view length");
    require(std::string_view(Policy::modeName(Policy::Mode::Off)) == "off", "off name");
    require(std::string_view(Policy::modeName(Policy::Mode::Bounded)) == "bounded", "bounded name");
    require(std::string_view(Policy::modeName(Policy::Mode::Full)) == "full", "full name");
}

void gateContracts()
{
    require(Policy::boundedMaximumRows == 1024 && Policy::boundedMaximumColumns == 1024,
        "bounded rows/columns contract changed without review");
    for (const auto mode : {Policy::Mode::Off, Policy::Mode::Bounded, Policy::Mode::Full})
    {
        require(!Policy::requested(mode, false), "runtime gate bypasses original parent diagnostics gate");
        require(Policy::requested(mode, true) == (mode != Policy::Mode::Off), "parent/runtime intersection");
    }
    for (const std::size_t size : {std::size_t{0}, std::size_t{1}, std::size_t{1023},
            std::size_t{1024}, std::size_t{1025}, std::numeric_limits<std::size_t>::max()})
    {
        require(!Policy::allowsRows(Policy::Mode::Off, size), "off admits diagnostic row allocation");
        require(!Policy::allowsColumns(Policy::Mode::Off, size), "off admits diagnostic column allocation");
        require(Policy::allowsRows(Policy::Mode::Full, size), "full no longer retains legacy unbounded input");
        require(Policy::allowsColumns(Policy::Mode::Full, size), "full no longer retains legacy unbounded state width");
        require(Policy::allowsRows(Policy::Mode::Bounded, size) == (size <= 1024), "bounded row boundary");
        require(Policy::allowsColumns(Policy::Mode::Bounded, size) == (size <= 1024), "bounded column boundary");
    }
    // Checking rows independently cannot accidentally approve a wide state,
    // or vice versa. No product of dimensions/overflow is involved.
    require(Policy::allowsRows(Policy::Mode::Bounded, 1)
        && !Policy::allowsColumns(Policy::Mode::Bounded, 1025), "wide diagnostic admitted");
    require(!Policy::allowsRows(Policy::Mode::Bounded, 1025)
        && Policy::allowsColumns(Policy::Mode::Bounded, 1), "tall diagnostic admitted");
}

void environmentContracts()
{
    RestoreEnvironment restore;
    environment(nullptr);
    require(Policy::modeFromEnvironment() == Policy::Mode::Off, "missing environment must not run SVD");
    environment("");
    require(Policy::modeFromEnvironment() == Policy::Mode::Off, "empty environment must not run SVD");
    for (const auto name : {"off", "bounded", "full"})
    {
        environment(name);
        require(Policy::modeFromEnvironment() == Policy::parseMode(name), "runtime mode not parsed exactly");
    }
    environment("on");
    bool rejected = false;
    try { (void)Policy::modeFromEnvironment(); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid environment does not fail explicitly");
}
}

int main()
{
    try
    {
        parseContracts();
        gateContracts();
        environmentContracts();
        std::cout << "PURE_OBSERVATION_RANK_POLICY checks=" << checks << " passed=true\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "PURE_OBSERVATION_RANK_POLICY failed=" << error.what() << '\n';
        return 1;
    }
}
