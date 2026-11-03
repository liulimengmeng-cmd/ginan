// Link against the real PEA objects, including ppp_ambres.cpp. The test calls
// the non-inline entry point used by ZhangR51PhysicalStore, not a replacement
// adapter or a second copy of its new implementation.
#include "common/zhangR51PhysicalConsequence.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangProductRelationSolver.hpp"
#include "common/zhangR48ExactAudit.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>

using Row = ZhangR51PhysicalRow;
using Rows = ZhangR51PhysicalRows;
static void require(bool condition, const char* why)
{ if (!condition) throw std::runtime_error(why); }

// Old full affine backend is the second reference. Deliberately request the
// complete particular solution and saturated physical kernel. No new
// consequence snapshot/helper participates in this reference.
static bool oldFullConsequence(const Rows& rows, const ZhangExactVector& values,
    bool valid, const Rows& targets, const ZhangExactVector& rhs)
{
    if (!valid || rows.size() != values.size() || targets.size() != rhs.size())
        return false;
    std::map<std::string, std::size_t> columns;
    for (const auto* group : {&rows, &targets})
        for (const auto& row : *group)
            for (const auto& [id, value] : row)
                if (value != 0 && !columns.contains(id)) columns[id] = columns.size();
    // A harmless free dummy coordinate lets the old API describe the empty
    // domain; it does not turn an unknown target coordinate into zero.
    const std::size_t dimension = std::max<std::size_t>(1, columns.size());
    ZhangExactMatrix h(rows.size(), ZhangExactVector(dimension));
    for (std::size_t r = 0; r < rows.size(); ++r)
        for (const auto& [id, value] : rows[r])
            if (value != 0) h[r][columns.at(id)] = value;
    const auto affine = zhangExactAffineIntegerQuotient(h, values, dimension,
        ZhangExactQuotientWork::PARTICULAR_AND_KERNEL);
    if (!affine.valid) return false;
    for (std::size_t r = 0; r < targets.size(); ++r)
    {
        ZhangExactInteger constant = 0;
        for (const auto& [id, value] : targets[r])
            if (value != 0) constant += value * affine.particularSolution[columns.at(id)];
        if (constant != rhs[r]) return false;
        for (const auto& direction : affine.kernelBasis)
        {
            ZhangExactInteger change = 0;
            for (const auto& [id, value] : targets[r])
                if (value != 0) change += value * direction[columns.at(id)];
            if (change != 0) return false;
        }
    }
    return true;
}

// Independent exhaustive integer oracle for two-variable fixtures. All input
// coefficients are <=3 and RHS <=9. Rank-two integral solutions have |x|,|y|
// <=54 by Cramer's rule. A feasible rank-one equation has a Bezout particular
// solution <=27 and a nonzero integer null step <=3, so [-64,64]^2 contains
// both a particular point and an independent free-direction step. Thus this
// finite enumeration decides feasibility AND constant linear consequences for
// this bounded *fixture family*, not just the sampled residuals.
static bool enumeratedConsequence(const std::vector<std::array<int, 2>>& h,
    const std::vector<int>& b, const std::array<int, 2>& target, int rhs)
{
    bool feasible = false;
    for (int x = -64; x <= 64; ++x)
    for (int y = -64; y <= 64; ++y)
    {
        bool accepted = true;
        for (std::size_t row = 0; row < h.size(); ++row)
            if (h[row][0] * x + h[row][1] * y != b[row])
            { accepted = false; break; }
        if (!accepted) continue;
        feasible = true;
        if (target[0] * x + target[1] * y != rhs) return false;
    }
    return feasible;
}

static std::size_t kernelCalls()
{
    const auto found = zhangR48ExactMetrics.find("R48_EXACT_KERNEL");
    return found == zhangR48ExactMetrics.end() ? 0 : found->second.calls;
}

int main()
{
    std::size_t compared = 0;
    double productionSeconds = 0, oldSeconds = 0;
    auto compare = [&](const Rows& rows, const ZhangExactVector& values,
        bool domainValid, const Rows& targets, const ZhangExactVector& rhs,
        bool independentExpected)
    {
        ZhangR51PhysicalConsequenceMetrics metrics;
        const auto kernelsBefore = kernelCalls();
        const auto start = std::chrono::steady_clock::now();
        const bool actual = zhangR51PhysicalConsequence(
            rows, values, domainValid, targets, rhs, &metrics);
        const auto middle = std::chrono::steady_clock::now();
        require(kernelCalls() == kernelsBefore,
            "production consequence requested a full physical kernel");
        require(metrics.builds <= 1, "one query rebuilt the domain");
        const bool legacy = oldFullConsequence(rows, values, domainValid, targets, rhs);
        const auto end = std::chrono::steady_clock::now();
        require(actual == independentExpected, "production disagrees with independent oracle");
        require(legacy == independentExpected, "old full backend disagrees with independent oracle");
        productionSeconds += std::chrono::duration<double>(middle - start).count();
        oldSeconds += std::chrono::duration<double>(end - middle).count();
        ++compared;
    };

    // Feasibility is mandatory even for a literal source row or empty query.
    compare({{{"x", 2}}}, {1}, true, {{{"x", 2}}}, {1}, false);
    compare({{{"x", 2}}}, {1}, true, {}, {}, false);
    compare({{{"x", 2}}}, {2}, true, {{{"x", 1}}}, {1}, true);
    require(!zhangIntegerRowLatticeContains({{2}}, {1}).contained,
        "test lost the distinction between consequence and integer row membership");
    compare({{{"x", 1}, {"old", -2}}}, {0}, true, {{{"x", 1}}}, {0}, false);
    compare({{{"old", 2}, {"live", 1}}}, {1}, true,
        {{{"old", 2}, {"live", 1}}}, {1}, true);
    compare({{{"old", 2}, {"live", 1}}, {{"live", 1}}}, {1, 0}, true,
        {{{"old", 2}}}, {1}, false);
    compare({{{"x", 1}}}, {3}, true, {{{"unknown_arc", 1}}}, {0}, false);
    compare({{{"x", 1}}}, {3}, true, {{{"x", 1}, {"unknown_arc", 0}}}, {3}, true);
    compare({{{"x", 1}}, {{"x", 1}}}, {0, 1}, true, {{{"x", 1}}}, {0}, false);
    compare({{}}, {1}, true, {{}}, {0}, false);
    compare({}, {}, true, {}, {}, true);
    compare({}, {}, true, {{}}, {0}, true);
    compare({}, {}, true, {{{"x", 1}}}, {0}, false);
    compare({{{"x", 1}}}, {}, true, {{{"x", 1}}}, {0}, false);
    compare({{{"x", 1}}}, {1}, false, {{{"x", 1}}}, {1}, false);
    compare({{{"x", 1}}}, {1}, true, {{{"x", 1}}}, {}, false);
    const ZhangExactInteger large = ZhangExactInteger(1) << 130;
    compare({{{"x", large}}}, {3 * large}, true, {{{"x", 1}}}, {3}, true);
    compare({{{"x", large}}}, {3 * large + 1}, true,
        {{{"x", large}}}, {3 * large + 1}, false);

    std::mt19937 random(20261102);
    for (int fixture = 0; fixture < 192; ++fixture)
    {
        const int count = random() % 6;
        std::vector<std::array<int, 2>> h;
        std::vector<int> b;
        Rows rows;
        ZhangExactVector values;
        const int plantedX = int(random() % 3) - 1;
        const int plantedY = int(random() % 3) - 1;
        for (int i = 0; i < count; ++i)
        {
            const std::array<int, 2> row{int(random() % 7) - 3, int(random() % 7) - 3};
            const int value = fixture % 2
                ? row[0] * plantedX + row[1] * plantedY : int(random() % 19) - 9;
            h.push_back(row); b.push_back(value);
            Row sparse;
            if (row[0]) sparse["x"] = row[0];
            if (row[1]) sparse["y"] = row[1];
            rows.push_back(std::move(sparse)); values.push_back(value);
        }
        const std::array<int, 2> target{int(random() % 7) - 3, int(random() % 7) - 3};
        const int rhs = fixture % 3
            ? target[0] * plantedX + target[1] * plantedY : int(random() % 19) - 9;
        Row sparseTarget;
        if (target[0]) sparseTarget["x"] = target[0];
        if (target[1]) sparseTarget["y"] = target[1];
        const bool expected = enumeratedConsequence(h, b, target, rhs);
        compare(rows, values, true, {sparseTarget}, {rhs}, expected);
        // Redundant source rows and insertion order cannot change the result.
        if (!rows.empty()) { rows.push_back(rows.front()); values.push_back(values.front()); }
        std::reverse(rows.begin(), rows.end());
        std::reverse(values.begin(), values.end());
        compare(rows, values, true, {sparseTarget}, {rhs}, expected);
    }
    std::cout << "PASS production_physical_consequence cases=" << compared
        << " independent_integer_enumerator=1 old_full_affine_reference=1"
        << " infeasible_direct_source_rejected=1 new_full_kernel_calls=0\n"
        << "BENCH production_seconds=" << productionSeconds
        << " old_full_backend_seconds=" << oldSeconds << '\n';
}
