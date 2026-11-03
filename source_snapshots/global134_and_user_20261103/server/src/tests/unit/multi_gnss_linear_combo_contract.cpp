// Production-linked LC guards; independent frequency-domain long-double oracle.
// Does not include producer implementation or use getLC() as a reference.
#include "common/acsConfig.hpp"
#include "common/acsQC.hpp"
#include "common/constants.hpp"
#include "common/linearCombo.hpp"
#include "common/navigation.hpp"
#include "common/observations.hpp"
#include "common/satStat.hpp"
#include <boost/json.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

void obs2lc(Trace&, GObs&, lc_t&);
namespace {
unsigned checks = 0, negativeCases = 0, validCases = 0, duplicateDefaultCases = 0;
double maximumOracleDifference = 0;
boost::json::array oracleRows;
void require(bool pass, const std::string& why)
{ ++checks; if (!pass) throw std::runtime_error(why); }
struct Field { const char* name; double S_LC::* member; };
const std::array<Field, 18> fields{{
    {"GF_Phas_m", &S_LC::GF_Phas_m}, {"GF_Code_m", &S_LC::GF_Code_m},
    {"WL_Phas_m", &S_LC::WL_Phas_m}, {"WL_Phas_c", &S_LC::WL_Phas_c},
    {"WL_Code_m", &S_LC::WL_Code_m}, {"WL_Code_c", &S_LC::WL_Code_c},
    {"NL_Phas_m", &S_LC::NL_Phas_m}, {"NL_Phas_c", &S_LC::NL_Phas_c},
    {"NL_Code_m", &S_LC::NL_Code_m}, {"NL_Code_c", &S_LC::NL_Code_c},
    {"IF_Phas_m", &S_LC::IF_Phas_m}, {"IF_Code_m", &S_LC::IF_Code_m},
    {"MW_m", &S_LC::MW_m}, {"MW_c", &S_LC::MW_c},
    {"lam_A", &S_LC::lam_A}, {"lam_B", &S_LC::lam_B},
    {"lam_WL", &S_LC::lam_WL}, {"lam_NL", &S_LC::lam_NL}}};
using Inputs = std::array<double, 6>;
const Inputs normal{21000000.375, 20999999.875, 21000001.25, 21000000.75,
                    299792458.0 / 1575420000.0, 299792458.0 / 1176450000.0};

std::array<long double, 18> reference(const Inputs& in, bool reverseMw)
{
    // Observation laws written in frequency coordinates, not the producer's
    // wavelength-coefficient arithmetic. Preserve legacy NL_Phas_m's WL
    // denominator intentionally: changing that physical definition is out of scope.
    const long double c = 299792458.L;
    const long double f = c / in[4], g = c / in[5];
    const long double a = in[0], b = in[1], p = in[2], q = in[3];
    const long double wl = c / (f-g), nl = c / (f+g);
    const long double phaseW = (f*a-g*b)/(f-g), codeN = (f*p+g*q)/(f+g);
    const long double mw = phaseW-codeN;
    return {a-b, p-q, phaseW, (f*a-g*b)/c, (f*p-g*q)/(f+g), (f*p-g*q)/c,
        (f*a+g*b)/(f-g), (f*a+g*b)/c, codeN, (f*p+g*q)/c,
        (f*f*a-g*g*b)/(f*f-g*g), (f*f*p-g*g*q)/(f*f-g*g), mw,
        (reverseMw ? -1 : 1)*mw/wl, in[4], in[5], wl, nl};
}
void finiteFields(const S_LC& value, const std::string& where)
{
    for (const auto& field : fields)
        require(std::isfinite(value.*field.member), where + ":NONFINITE:" + field.name);
}
void match(const S_LC& value, const Inputs& in, bool reverseMw, const std::string& name)
{
    ++validCases; require(value.valid, name + ":EXPECTED_VALID"); finiteFields(value, name);
    const auto expected = reference(in, reverseMw);
    boost::json::array measured, independentlyExpected;
    for (std::size_t i = 0; i < fields.size(); ++i) {
        const double difference = std::abs((long double)(value.*fields[i].member) - expected[i]);
        maximumOracleDifference = std::max(maximumOracleDifference, difference);
        // Cancellation at ~21 Mm is tested in absolute native units, not a
        // metre-scale relative tolerance hiding GF/MW errors. 0.2 microunit.
        require(difference <= 2e-7, name + ":ORACLE_MISMATCH:" + fields[i].name +
            ":difference=" + std::to_string(difference));
        measured.push_back(value.*fields[i].member);
        independentlyExpected.push_back(static_cast<double>(expected[i]));
    }
    oracleRows.push_back(boost::json::object{{"case", name}, {"MW_sign_reversed", reverseMw},
        {"actual_18_fields", measured}, {"independent_expected_18_fields", independentlyExpected}});
}
S_LC scalar(const Inputs& in)
{
    double c1 = 0, c2 = 0;
    auto value = getLC(in[0], in[1], in[2], in[3], in[4], in[5], &c1, &c2);
    require(std::isfinite(c1) && std::isfinite(c2), "NONFINITE_IF_COEFFICIENT_OUTPUT");
    return value;
}
void invalidScalar(const Inputs& in, const std::string& name)
{
    ++negativeCases; const auto value = scalar(in);
    require(!value.valid, name + ":INVALID_SCALAR_ACCEPTED"); finiteFields(value, name);
}

struct Fixture {
    SatNav satellite;
    SatStat status;
    GObs obs;
    Fixture(E_Sys system, E_FType second, Inputs in = normal) {
        obs.Sat = SatSys(system, system == E_Sys::BDS ? 19 : 13);
        obs.time.bigTime = 1405216800;
        obs.satNav_ptr = &satellite; obs.satStat_ptr = &status;
        satellite.lamMap[F1] = in[4]; satellite.lamMap[second] = in[5];
        Sig first{}, other{};
        first.code = system == E_Sys::BDS ? E_ObsCode::L1P : E_ObsCode::L1C;
        other.code = second == F2 ? E_ObsCode::L2W :
            system == E_Sys::BDS ? E_ObsCode::L5P : E_ObsCode::L5Q;
        first.L = in[0]/in[4]; first.P = in[2];
        other.L = in[1]/in[5]; other.P = in[3];
        obs.sigs[F1] = first; obs.sigs[second] = other;
        acsConfig.code_priorities[system] = {first.code, other.code};
    }
    Inputs observed(E_FType second) const {
        return {obs.sigs.at(F1).L * satellite.lamMap.at(F1),
                obs.sigs.at(second).L * satellite.lamMap.at(second),
                obs.sigs.at(F1).P, obs.sigs.at(second).P,
                satellite.lamMap.at(F1), satellite.lamMap.at(second)};
    }
};
void finiteBase(const lc_t& base, const std::string& name)
{
    for (const auto* values : {&base.L_m, &base.P, &base.mp})
        for (const auto& [frequency, value] : *values)
            require(std::isfinite(value), name + ":NONFINITE_DERIVED_BASE_MAP");
    for (const auto& [pair, value] : base.lcMap) {
        finiteFields(value, name);
        if (pair.first == pair.second) require(!value.valid, name + ":VALID_SAME_FREQUENCY");
    }
}
}

int main()
try {
    for (const E_FType second : {F2, F5}) {
        Inputs in = normal; if (second == F2) in[5] = 299792458.0 / 1227600000.0;
        match(scalar(in), in, false, "scalar_normal");
        in[0] = -in[0]; in[1] = -in[1];
        match(scalar(in), in, false, "scalar_negative_phase_legal");
    }
    for (int component = 0; component < 6; ++component)
        for (double value : {0., std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()}) {
            auto in = normal; in[component] = value;
            invalidScalar(in, "scalar_bad_component_" + std::to_string(component));
        }
    for (int component : {2, 3, 4, 5}) {
        auto in = normal; in[component] = -1;
        invalidScalar(in, "negative_code_or_wavelength");
    }
    auto duplicate = normal; duplicate[5] = duplicate[4]; invalidScalar(duplicate, "duplicate_wavelength");
    auto tiny = normal; tiny[4] = std::numeric_limits<double>::denorm_min(); invalidScalar(tiny, "tiny_wave_overflow");
    auto overflow = normal; overflow[0] = 1e308; overflow[1] = -1e308; invalidScalar(overflow, "finite_observation_overflow");
    auto squareOverflow = normal; squareOverflow[4] = 1e200; squareOverflow[5] = 2e200;
    invalidScalar(squareOverflow, "finite_wavelength_squared_overflow");

    std::ostringstream trace;
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS}) {
        for (const E_FType second : {F2, F5}) {
            if (second == F2 && system != E_Sys::GPS) continue;
            auto in = normal; if (second == F2) in[5] = 299792458.0 / 1227600000.0;
            Fixture fixture(system, second, in);
            lc_t base{}; obs2lc(trace, fixture.obs, base);
            finiteBase(base, "obs2lc_normal");
            const bool flip = system != E_Sys::GPS && second == F5;
            match(getLC(base, F1, second), fixture.observed(second), flip, "obs2lc_pair");
            unsigned validCount = 0; for (const auto& [pair, value] : base.lcMap) validCount += value.valid;
            require(validCount == 1, "TWO_FREQ_FIXTURE_MORE_THAN_ONE_VALID_PAIR");
            E_FType f1, f2, f3; require(satFreqs(system, f1, f2, f3), "SAT_FREQS_FAILED");
            if (second == F5) { ++duplicateDefaultCases; require(f2 == F5 && f3 == F5, "DUPLICATE_DEFAULT_F5_NOT_EXERCISED"); }
            const auto& same = getLC(fixture.obs, base, second, second);
            require(!same.valid, "SAME_FREQUENCY_OBS_ACCEPTED"); finiteFields(same, "same_freq");

            lc_t missing{}; fixture.obs.sigs[second].P = 0;
            const auto& noP = getLC(fixture.obs, missing, F1, second);
            require(!noP.valid, "MISSING_CODE_OUTER_VALID_TRUE"); finiteBase(missing, "missing_code");
        }
        Fixture fixture(system, F5);
        for (int bad = 0; bad < 5; ++bad) {
            lc_t base{};
            auto& cached = base.lcMap[{F1, F5}]; cached = scalar(normal);
            if (bad == 0) cached.MW_c = std::numeric_limits<double>::quiet_NaN();
            if (bad == 1) cached.lam_B = cached.lam_A;
            if (bad == 2) cached.lam_A = 0;
            if (bad == 3) cached.lam_B = -1;
            if (bad == 4) { cached.valid = false; cached.MW_c = std::numeric_limits<double>::quiet_NaN(); }
            const auto& lookup = getLC(base, F1, F5);
            require(!lookup.valid, "INVALID_CACHED_TRUE_PRESERVED_BY_LOOKUP");
            finiteFields(lookup, "invalid_cache_lookup");
            const auto& rebuilt = getLC(fixture.obs, base, F1, F5);
            finiteBase(base, "invalid_cache_obs");
            if (rebuilt.valid) match(rebuilt, fixture.observed(F5), system != E_Sys::GPS, "cache_recomputed");
        }
        for (int bad = 0; bad < 6; ++bad) {
            Fixture invalid(system, F5);
            if (bad == 0) invalid.satellite.lamMap[F5] = 0;
            if (bad == 1) invalid.satellite.lamMap[F5] = -1;
            if (bad == 2) invalid.obs.sigs[F5].L = std::numeric_limits<double>::quiet_NaN();
            if (bad == 3) invalid.obs.sigs[F5].P = std::numeric_limits<double>::infinity();
            if (bad == 4) invalid.obs.sigs[F5].P = -1;
            if (bad == 5) invalid.satellite.lamMap[F5] = invalid.satellite.lamMap[F1];
            lc_t base{}; obs2lc(trace, invalid.obs, base);
            finiteBase(base, "bad_obs2lc_input");
            for (const auto& [pair, value] : base.lcMap) require(!value.valid, "BAD_OBS2LC_PAIR_VALID");
        }
    }
    boost::json::object result{{"status", "PASS"}, {"assertions", checks},
        {"scalar_negative_cases", negativeCases}, {"independent_18_field_oracle_cases", validCases},
        {"duplicate_default_F5_system_cases", duplicateDefaultCases},
        {"maximum_absolute_oracle_difference_native", maximumOracleDifference},
        {"oracle", "Independent long-double frequency-domain equations; no production core call in reference"},
        {"negative_phase_is_legal", true}, {"legacy_NL_phase_definition_preserved", true},
        {"producer_recompiled_in_driver", false},
        {"oracle_rows", oracleRows},
        {"scope", "PRODUCTION_LINKED_CONTROLLED_LC_CONTRACT_NOT_REAL_NETWORK_QUALITY"}};
    boost::json::array names; for (const auto& field : fields) names.push_back(field.name);
    result["field_order"] = std::move(names);
    std::cout << "LINEAR_COMBO_CONTRACT_RESULT " << boost::json::serialize(result) << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "LINEAR_COMBO_CONTRACT_FAILED assertions=" << checks << " reason=" << error.what() << '\n'; return 1;
}
