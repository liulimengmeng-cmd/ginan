// Mathematical planning/reference tests. These do not replace the required
// production registration/transform/checkpoint tests after P2 integration.
#include "common/zhangProductFloatCarriers.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangR47ProductDomain.hpp"
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>

static void check(bool condition, const char* why)
{ if (!condition) throw std::runtime_error(why); }
using Rational = ZhangFloatCarrierRational;
using Row = ZhangFloatCarrierSparseRow;

static std::vector<std::vector<long double>> denseTransform(const ZhangFloatCarrierRowPlan& plan)
{
    std::vector<std::vector<long double>> t(plan.transformRows.size(),
        std::vector<long double>(plan.sourceDimension));
    for (std::size_t i = 0; i < plan.transformRows.size(); ++i)
        for (const auto& [j, value] : plan.transformRows[i]) t[i][j] = value.convert_to<long double>();
    return t;
}

// Independent literal loop oracle; no Eigen or production transform call.
static std::vector<std::vector<long double>> projectedCovariance(
    const std::vector<std::vector<long double>>& t,
    const std::vector<std::vector<long double>>& p)
{
    std::vector<std::vector<long double>> result(t.size(), std::vector<long double>(t.size()));
    for (std::size_t i = 0; i < t.size(); ++i)
    for (std::size_t j = 0; j < t.size(); ++j)
    for (std::size_t a = 0; a < p.size(); ++a)
    for (std::size_t b = 0; b < p.size(); ++b)
        result[i][j] += t[i][a] * p[a][b] * t[j][b];
    return result;
}

static void carrierRowsAndFullQ()
{
    std::vector<ZhangFloatCarrierIndexedTarget> targets{
        {"phase", "m", {{0, 1}, {1, 1}, {2, 1}}, 0, true},
        {"scaled", "m", {{1, 2}, {2, 2}}, 3, true},
        {"live", "m", {{0, -3}}, 0, true}};
    const auto plan = planZhangProductFloatCarrierRows(3, {0}, targets);
    check(plan.valid && plan.reducesDimension && plan.carrierRank == 1,
        "whole-target rank compression failed");
    check(plan.transformRows == std::vector<Row>{{{0, 1}}, {{1, 1}, {2, 1}}},
        "carrier is not the retained whole raw target");
    std::reverse(targets.begin(), targets.end());
    const auto reversed = planZhangProductFloatCarrierRows(3, {0}, targets);
    check(reversed.valid && reversed.transformRows == plan.transformRows &&
        reversed.carrierSourceTargetIds == plan.carrierSourceTargetIds,
        "input registration order changed the carrier basis");
    const std::vector<std::vector<long double>> p{{4, 1, 2}, {1, 3, .5L}, {2, .5L, 5}};
    const auto q = projectedCovariance(denseTransform(plan), p);
    check(q == std::vector<std::vector<long double>>{{4, 3}, {3, 9}},
        "clock/carrier covariance was lost");
    // Independent scalar-observation Gaussian oracle: conditioning using
    // retained coordinates commutes with prior compression/marginalization.
    auto condition = [](std::vector<long double> mean,
        std::vector<std::vector<long double>> covariance,
        const std::vector<long double>& h)
    {
        std::vector<long double> ph(mean.size());
        for (std::size_t i = 0; i < mean.size(); ++i)
        for (std::size_t j = 0; j < mean.size(); ++j) ph[i] += covariance[i][j] * h[j];
        long double variance = 2;
        for (std::size_t i = 0; i < mean.size(); ++i) variance += h[i] * ph[i];
        for (std::size_t i = 0; i < mean.size(); ++i)
        {
            mean[i] += ph[i] * .8L / variance;
            for (std::size_t j = 0; j < mean.size(); ++j)
                covariance[i][j] -= ph[i] * ph[j] / variance;
        }
        return std::pair{mean, covariance};
    };
    const auto [fullMean, fullCovariance] = condition({10, 2, 5}, p, {1, .25L, .25L});
    const auto [smallMean, smallCovariance] = condition({10, 7}, q, {1, .25L});
    const auto transportedPosterior = projectedCovariance(denseTransform(plan), fullCovariance);
    check(std::abs(fullMean[0] - smallMean[0]) < 1e-16L &&
        std::abs(fullMean[1] + fullMean[2] - smallMean[1]) < 1e-16L,
        "future live-observation conditional mean changed after compression");
    for (std::size_t i = 0; i < smallCovariance.size(); ++i)
    for (std::size_t j = 0; j < smallCovariance.size(); ++j)
        check(std::abs(transportedPosterior[i][j] - smallCovariance[i][j]) < 1e-16L,
            "future live-observation conditional full covariance changed");
    const auto whole = analyzeZhangBridgeSupportOnActualSpan({}, {{1, 1}}, {{1, 1}, {1, 0}}, 2);
    check(whole.valid && !whole.representable && whole.unresolvedTargets == std::vector<std::size_t>{1},
        "individual arc failure vetoed a representable whole sum");
    // A newly requested independent direction cannot be invented after only
    // a sum was retained. The planner can add a direction only while its
    // original source coordinates still exist.
    targets.push_back({"difference", "m", {{1, 1}, {2, -1}}, 0, true});
    const auto extended = planZhangProductFloatCarrierRows(3, {0}, targets);
    check(extended.valid && extended.carrierRank == 2 && !extended.reducesDimension,
        "independent active target was silently discarded");
    targets.front().complete = false;
    check(!planZhangProductFloatCarrierRows(3, {0}, targets).valid,
        "unbound active consumer was omitted");
    check(!planZhangProductFloatCarrierRows(3, {0, 0}, {}).valid,
        "duplicate protected index accepted");

    // Nonzero inter-system and inter-signal blocks must survive global TPT'.
    const std::vector<long double> u{1, .1L, .2L, .3L, .4L, .5L, .6L};
    std::vector<std::vector<long double>> globalP(7, std::vector<long double>(7));
    for (std::size_t i = 0; i < 7; ++i)
    for (std::size_t j = 0; j < 7; ++j) globalP[i][j] = (i == j ? 1 : 0) + u[i] * u[j];
    const auto global = planZhangProductFloatCarrierRows(7, {0}, {
        {"G", "m", {{1, 1}, {2, 1}}, 0, true},
        {"E", "m", {{3, 1}, {4, 1}}, 0, true},
        {"C", "m", {{5, 1}, {6, -1}}, 0, true}});
    check(global.valid && global.carrierRank == 3, "global carrier plan invalid");
    const auto globalQ = projectedCovariance(denseTransform(global), globalP);
    // Carrier ordering is canonical by target ID: C, E, G.
    const std::vector<std::vector<long double>> expected{
        {2, -.1L, .7L, .3L}, {-.1L, 2.01L, -.07L, -.03L},
        {.7L, -.07L, 2.49L, .21L}, {.3L, -.03L, .21L, 2.09L}};
    for (std::size_t i = 0; i < expected.size(); ++i)
    for (std::size_t j = 0; j < expected.size(); ++j)
        check(std::abs(globalQ[i][j] - expected[i][j]) < 1e-17L,
            "global complete covariance differs from independent oracle");

    const auto before = zhangR47CompileProductSearchFrame({{2}}, {}, {}, 1,
        ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);
    const auto floating = planZhangProductFloatCarrierRows(2, {0}, {
        {"half", "m", {{1, Rational(1) / 2}}, 0, true}});
    const auto after = zhangR47CompileProductSearchFrame({{2}}, {}, {}, 1,
        ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);
    check(floating.valid && before.valid && after.valid &&
        before.imageGenerators == after.imageGenerators &&
        !zhangIntegerRowLatticeContains(before.imageGenerators, {1}).contained,
        "FLOAT carrier changed the independent exact integer image 2Z");
}

// For three variables with coefficients <=2, a null basis has integer
// directions bounded by the 2x2 minors (<=8); lower ranks need only <=2.
// Enumerating [-8,8]^3 is therefore complete for these homogeneous fixtures.
static bool enumeratedObservable(const std::vector<std::array<int, 3>>& rows,
    const std::array<int, 3>& target)
{
    for (int x = -8; x <= 8; ++x)
    for (int y = -8; y <= 8; ++y)
    for (int z = -8; z <= 8; ++z)
    {
        bool null = true;
        for (const auto& row : rows)
            if (row[0] * x + row[1] * y + row[2] * z != 0) { null = false; break; }
        if (null && target[0] * x + target[1] * y + target[2] * z != 0) return false;
    }
    return true;
}

static void generalSpanAndGauge()
{
    const auto exact = analyzeZhangBridgeSupportOnActualSpan({{1, -1}}, {{1, 1}}, {{1, 0}}, 2);
    check(exact.valid && exact.representable &&
        exact.currentCombinations[0] == ZhangBridgeSpanRow{Rational(1) / 2} &&
        exact.numericCombinations[0] == ZhangBridgeSpanRow{Rational(1) / 2},
        "current H and actual S rational witnesses differ");
    const auto gauge = analyzeZhangBridgeSupportOnActualSpan({}, {{1, 1, 1}}, {{1, 1, 0}}, 3);
    check(gauge.valid && !gauge.representable && !gauge.counterexamples.empty(),
        "a nuisance-bearing phase carrier was interpreted as an integer sum");
    const auto cancelled = analyzeZhangBridgeSupportOnActualSpan({},
        {{1, 0, 1}, {0, 1, 1}}, {{1, -1, 0}}, 3);
    check(cancelled.valid && cancelled.representable,
        "explicit nuisance cancellation was not resolved");
    check(!analyzeZhangBridgeSupportOnActualSpan({{1}}, {{1, 2}}, {}, 1).valid,
        "malformed general span accepted");
    std::mt19937 random(202611021);
    std::size_t checked = 0;
    for (int fixture = 0; fixture < 192; ++fixture)
    {
        ZhangExactMatrix h;
        ZhangBridgeSpanMatrix s;
        std::vector<std::array<int, 3>> independentRows;
        const int hCount = random() % 4, sCount = random() % 4;
        for (int row = 0; row < hCount + sCount; ++row)
        {
            const std::array<int, 3> values{int(random() % 5) - 2,
                int(random() % 5) - 2, int(random() % 5) - 2};
            independentRows.push_back(values);
            if (row < hCount) h.push_back({values[0], values[1], values[2]});
            else s.push_back({values[0], values[1], values[2]});
        }
        for (int query = 0; query < 3; ++query)
        {
            const std::array<int, 3> t{int(random() % 5) - 2,
                int(random() % 5) - 2, int(random() % 5) - 2};
            const auto result = analyzeZhangBridgeSupportOnActualSpan(h, s,
                {{t[0], t[1], t[2]}}, 3);
            check(result.valid && result.representable == enumeratedObservable(independentRows, t),
                "general S disagrees with independent complete integer-null enumeration");
            ++checked;
        }
    }
    for (const double value : {0., -0., 1., .1,
            std::numeric_limits<double>::min(), std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::max()})
    {
        const auto rational = zhangFloatCarrierExactDouble(value);
        check(rational && rational->convert_to<double>() == value,
            "dyadic coefficient conversion changed its binary value");
    }
    check(!zhangFloatCarrierExactDouble(std::numeric_limits<double>::infinity()),
        "nonfinite carrier coefficient accepted");
    std::cout << "PASS general_actual_S independent_null_enumeration=" << checked
        << " explicit_nuisance=1 exact_binary_coefficients=1\n";
}

int main()
{
    carrierRowsAndFullQ(); generalSpanAndGauge();
    std::cout << "PASS float_target_carrier_math rank_modulo_live=1 whole_targets=1 "
        "clock_crossQ=1 cross_system_crossQ=1 integer_image_unchanged=1 "
        "future_gaussian_update_commutes=1 "
        "production_integration_not_yet_claimed=1\n";
}
