// Pure exact-representation / physical-image regression. This does not claim
// to exercise the production graph admission entry or a real GNSS replay.
#include "common/zhangCurrentBridgeRepresentation.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
void check(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}

void near(double actual, double expected, const char* reason)
{
    check(std::isfinite(actual) && std::abs(actual - expected) < 1e-12, reason);
}

ZhangR51Rational evaluate(const ZhangR51RationalRow& row,
    const ZhangR51Rational& offset, const std::vector<int>& physical)
{
    check(row.size() == physical.size(), "oracle argument shape");
    ZhangR51Rational result = offset;
    for (std::size_t c = 0; c < row.size(); ++c) result += row[c] * physical[c];
    return result;
}
}

int main()
{
    // Independent affine oracle: old=live+7, q=old+2=live+9.
    const auto simple = zhangRewriteCurrentBridgeRepresentation(
        {{1,0}}, {2}, {{1,-1}}, {7}, {0});
    check(simple.valid && simple.projector == std::vector<ZhangR51RationalRow>{{0,1}}
        && simple.offsets == ZhangR51RationalRow{9}, "translation and offset rewrite");
    for (int live = -20; live <= 20; ++live)
        check(evaluate(simple.projector[0], simple.offsets[0], {live+7,live}) == live+9,
            "translation independent integer oracle");

    // Complete dual target oracle, not a second invocation of the rewriter:
    // old1=2*a+b+7; old2=a-b-3; WL=old1-old2=a+2*b+10;
    // L1=old1=2*a+b+7. Their full cross covariance is not diagonalized.
    const auto dual = zhangRewriteCurrentBridgeRepresentation(
        {{1,-1,0,0},{1,0,0,0}}, {0,0},
        {{1,0,-2,-1},{0,1,-1,1}}, {7,-3}, {0,1});
    check(dual.valid && dual.projector == std::vector<ZhangR51RationalRow>{{0,0,1,2},{0,0,2,1}}
        && dual.offsets == ZhangR51RationalRow{10,7}, "dual full target rewrite");
    for (int a = -7; a <= 7; ++a) for (int b = -7; b <= 7; ++b)
    {
        const std::vector<int> physical{2*a+b+7,a-b-3,a,b};
        check(evaluate(dual.projector[0], dual.offsets[0], physical) == a+2*b+10,
            "WL independent oracle");
        check(evaluate(dual.projector[1], dual.offsets[1], physical) == 2*a+b+7,
            "L1 independent oracle");
    }
    const std::array<double,2> liveMean{0.1,0.2};
    const std::array<std::array<double,2>,2> liveCovariance{{{4,1.5},{1.5,9}}};
    std::array<std::array<double,2>,2> numeric{};
    for (int r = 0; r < 2; ++r)
    {
        double mean = dual.offsets[r].convert_to<double>();
        for (int c = 0; c < 2; ++c)
        {
            numeric[r][c] = dual.projector[r][c+2].convert_to<double>();
            mean += numeric[r][c] * liveMean[c];
        }
        near(mean, r == 0 ? 10.5 : 7.4, "dual closed-form mean");
    }
    const std::array<std::array<double,2>,2> covarianceOracle{{{46,33.5},{33.5,31}}};
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b)
    {
        double covariance = 0;
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j)
            covariance += numeric[a][i] * liveCovariance[i][j] * numeric[b][j];
        near(covariance, covarianceOracle[a][b], "dual closed-form full covariance");
    }

    // Nonprimitive visible lattice: 2*old+live=1. The original target 2*old
    // has image 2Z. q=old is integer, but live must be odd. Never replace this
    // by an unconstrained/saturated real live direction or divide the lift.
    auto parityFrame = zhangR51PhysicalImage({{2,0}}, {{2,1}}, {1}, {{2,0}}, 2);
    check(parityFrame.valid && parityFrame.generators == ZhangExactMatrix{{2}}
        && parityFrame.particularTarget == ZhangExactVector{0}
        && parityFrame.projector == std::vector<ZhangR51RationalRow>{{1,0}}
        && parityFrame.offsets == ZhangR51RationalRow{0}, "original nonprimitive image oracle");
    const auto originalGenerators = parityFrame.generators;
    const auto originalParticular = parityFrame.particularTarget;
    const std::array<std::string,2> physicalOrder{"old","live"};
    const auto originalOrder = physicalOrder;
    const auto parity = zhangRewriteCurrentBridgeRepresentation(parityFrame.projector,
        parityFrame.offsets, {{2,1}}, {1}, {0});
    check(parity.valid && parity.projector == std::vector<ZhangR51RationalRow>{{0,ZhangR51Rational(-1)/2}}
        && parity.offsets == ZhangR51RationalRow{ZhangR51Rational(1)/2},
        "nonprimitive rational representation, constant retained");
    parityFrame.projector = parity.projector;
    parityFrame.offsets = parity.offsets;
    check(parityFrame.generators == originalGenerators && parityFrame.particularTarget == originalParticular
        && physicalOrder == originalOrder, "rewrite cannot alter image generators/particular/order");
    for (int old = -12; old <= 12; ++old)
    {
        const int live = 1-2*old;
        check(evaluate(parity.projector[0], parity.offsets[0], {old,live}) == old,
            "parity exact q image independent oracle");
        check(live % 2 != 0, "visible support remains odd");
    }
    for (int even = -20; even <= 20; even += 2)
        check(denominator(evaluate(parity.projector[0], parity.offsets[0], {0,even})) == 2,
            "even live is not in integer q image");
    near(parity.offsets[0].convert_to<double>()
        + parity.projector[0][1].convert_to<double>() * -5.2, 3.1,
        "nonprimitive closed-form q mean");
    near(std::pow(parity.projector[0][1].convert_to<double>(),2) * 4.0, 1.0,
        "nonprimitive closed-form q variance");
    ZhangExactMatrix lift;
    ZhangExactVector liftedValues;
    check(zhangR51LiftRationalRows({{1}}, {3}, parity.projector, parity.offsets, lift, liftedValues),
        "nonprimitive exact lift");
    check(lift == ZhangExactMatrix{{0,-1}} && liftedValues == ZhangExactVector{5},
        "q=3 implies live=-5, not an unconstrained rounded old");
    // Check the ORIGINAL physical target, not the rewritten numerical target.
    auto allConditions = ZhangExactMatrix{{2,1}};
    auto allValues = ZhangExactVector{1};
    allConditions.insert(allConditions.end(), lift.begin(), lift.end());
    allValues.insert(allValues.end(), liftedValues.begin(), liftedValues.end());
    const auto finalTarget = zhangR51PhysicalImage({{2,0}}, allConditions, allValues, {{2,0}}, 2);
    check(finalTarget.valid && finalTarget.projector.empty()
        && finalTarget.particularTarget == ZhangExactVector{6}, "original target lift oracle 2*old=6");

    // Discover another unavailable variable introduced by a first rewrite.
    const auto first = zhangRewriteCurrentBridgeRepresentation({{1,0,0}}, {0},
        {{1,-1,0},{0,1,-1}}, {2,3}, {0});
    check(first.valid && first.projector == std::vector<ZhangR51RationalRow>{{0,1,0}}
        && first.offsets == ZhangR51RationalRow{2}, "first missing support discovery");
    const auto second = zhangRewriteCurrentBridgeRepresentation({{1,0,0}}, {0},
        {{1,-1,0},{0,1,-1}}, {2,3}, {0,1});
    check(second.valid && second.projector == std::vector<ZhangR51RationalRow>{{0,0,1}}
        && second.offsets == ZhangR51RationalRow{5}, "second original-coordinate rewrite");

    const auto unavailable = zhangRewriteCurrentBridgeRepresentation({{1,0}}, {0}, {}, {}, {0});
    check(!unavailable.valid && unavailable.unresolvedColumns == std::set<std::size_t>{0},
        "missing current constraint cannot assign old zero");
    const auto hiddenFreedom = zhangRewriteCurrentBridgeRepresentation({{1,0}}, {0}, {{1,-1}}, {7}, {0,1});
    check(!hiddenFreedom.valid && hiddenFreedom.unresolvedColumns == std::set<std::size_t>{1},
        "unobservable hidden freedom rejects instead of history restoration");
    check(!zhangRewriteCurrentBridgeRepresentation({{1,0}}, {0}, {{1}}, {1}, {0}).valid,
        "malformed condition cannot authorize substitution");
    check(!zhangRewriteCurrentBridgeRepresentation({{1,0}}, {0}, {{1,-1}}, {7}, {2}).valid,
        "out-of-range missing column rejects");
    const auto unchanged = zhangRewriteCurrentBridgeRepresentation({{1,0}}, {2}, {}, {}, {});
    check(unchanged.valid && unchanged.projector == std::vector<ZhangR51RationalRow>{{1,0}}
        && unchanged.offsets == ZhangR51RationalRow{2}, "healthy original support unchanged");

    ZhangExactWorkBudget budget;
    budget.maximumCheckpoints = 3;
    ZhangCurrentBridgeRepresentation cancelled;
    check(!zhangRunWithinExactWorkBudget(budget, [&] {
        return zhangRewriteCurrentBridgeRepresentation({{1,0}}, {0}, {{2,1}}, {1}, {0});
    }, cancelled) && budget.exhausted && !cancelled.valid && cancelled.projector.empty(),
        "cancelled rewrite cannot expose partial success");

    std::cout << "R51 current bridge exact representation, nonprimitive image, original lift, "
        "closed-form full covariance and fail-closed support PASS\n";
}
