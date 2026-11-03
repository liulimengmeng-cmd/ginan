// Structural-zero regression with literal graph, integer, and covariance
// oracles. The oracle does not call the production target/pullback builder.
#include "common/zhangProductRelationBasis.hpp"
#include "common/zhangR47Candidate.hpp"
#include "common/zhangR47ProductDomain.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

static int checks = 0;
static void require(bool condition, const char* why)
{
    ++checks;
    if (!condition) throw std::runtime_error(why);
}

int main()
{
    const SatSys g01(E_Sys::GPS, 1), g02(E_Sys::GPS, 2);
    const std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    const auto current = zhangBuildSpanningTree(edges, "R0",
        {{"R0", g01}, {"R0", g02}, {"R1", g01}});
    require(current.connected, "literal K2,2 graph is disconnected");
    const auto zero = ProductRelationBasisBuilder::build(
        current, current, g01, E_Sys::GPS, E_ObsCode::L1C);
    require(zero.valid && zero.structuralMapValid,
        "identical S-bases rejected their exact structural map");
    require(zero.zeroDatumCorrection && zero.exactRank == 0 &&
        zero.datumCorrectionRank == 0, "zero G was confused with nonzero rank");
    require(zero.canonicalTargetRank == 1 && zero.namedRelationCount == 1 &&
        zero.namedRelations.size() == 1 && zero.satellites.size() == 2,
        "zero G truncated the physical product catalogue");
    require(zero.structuralZeroNamedIndices == std::vector<int>{0} &&
        zero.namedRelations[0].currentCycleCoefficients == ZhangExactVector{0},
        "zero correction was not recorded as a structural named identity");
    require(zero.exactRowBasis.empty() && zero.independentNamedIndices.empty() &&
        zero.affineRecoveryParents.empty(),
        "zero structural correction manufactured an integer decision");

    // Literal zero expressions preserve all named targets and affine shifts;
    // they do not become positive-rank searches or disappear from the domain.
    const auto whole = zhangCompileWholeProductLattice({{0}, {0}}, {0, 7}, {true});
    require(whole.valid && !whole.hasNonzeroSearchDirections &&
        whole.wholeAvailableRank == 0 && whole.searchPosteriorRows.empty(),
        "structural zeros were advertised as search directions");
    require(whole.namedExpressions.size() == 2 &&
        whole.namedRecoveryOffsets == ZhangExactVector({0, 7}) &&
        whole.namedRecoverable == std::vector<bool>({true, true}),
        "zero search domain lost named constants or affine shifts");
    require(whole.structuralIdentityRows.size() == 2,
        "zero domain lost exact structural identities");
    require(!zhangR47CandidateContractValid(ZhangR47Candidate{}),
        "a structural zero without a current decision passed the AR contract");

    // An unknown nonzero coordinate must not silently turn into zero.
    const auto missing = zhangCompileWholeProductLattice({{1}}, {0}, {false});
    require(missing.valid && !missing.namedRecoverable[0] &&
        missing.namedStatus[0] == "NOT_CURRENTLY_REPRESENTABLE",
        "missing nonzero physical support was accepted as a zero identity");

    const auto product = zhangBuildSpanningTree(edges, "R0",
        {{"R0", g01}, {"R1", g01}, {"R1", g02}});
    const auto nonzero = ProductRelationBasisBuilder::build(
        current, product, g01, E_Sys::GPS, E_ObsCode::L1C);
    require(nonzero.valid && nonzero.structuralMapValid &&
        !nonzero.zeroDatumCorrection && nonzero.datumCorrectionRank == 1,
        "nonzero datum correction changed its previous rank contract");
    require(nonzero.namedRelations[0].currentCycleCoefficients == ZhangExactVector{1},
        "literal K2,2 product correction differs from +k");
    // Independent incidence arithmetic: N00=2,N01=5,N10=7,N11=14.
    const long long nativePotential = 5 - 2;
    const long long cycle = 14 - 7 - 5 + 2;
    const long long productPotential = 14 - 7;
    require(nativePotential + cycle == productPotential,
        "literal complete product potential is not invariant");
    require(nativePotential != cycle && cycle != productPotential,
        "test failed to distinguish G*k from the complete product potential");

    // A nonzero target can have no residual image ONLY AFTER a certified
    // conditioner fixes it. This is different from an identically zero G.
    const auto fixedImage = zhangR47CompileProductSearchFrame({{1}}, {{1}}, {5}, 1);
    const auto freeImage = zhangR47CompileProductSearchFrame({{1}}, {}, {}, 1);
    require(fixedImage.valid && fixedImage.searchRank == 0 &&
        freeImage.valid && freeImage.searchRank == 1,
        "zero datum correction was confused with zero conditional image rank");

    // G=0 still leaves a native phase state with nonzero uncertainty. A
    // relevant current ambiguity constraint updates it through cross blocks.
    // Literal covariance P=[[4,1],[1,9]], mean=(10,3), constraint a=5.
    Eigen::Vector2d mean; mean << 10, 3;
    Eigen::Matrix2d covariance; covariance << 4, 1, 1, 9;
    Eigen::RowVector2d conditioner; conditioner << 0, 1;
    const double innovationVariance =
        (conditioner * covariance * conditioner.transpose())(0, 0);
    const Eigen::Vector2d gain = covariance * conditioner.transpose() /
        innovationVariance;
    const Eigen::Vector2d conditionalMean = mean + gain *
        (5 - (conditioner * mean)(0, 0));
    const Eigen::Matrix2d conditionalCovariance = covariance -
        gain * conditioner * covariance;
    require(std::abs(conditionalMean[0] - (10.0 + 2.0 / 9.0)) < 1e-12 &&
        std::abs(conditionalCovariance(0, 0) - 35.0 / 9.0) < 1e-12,
        "zero G native state lost its nonzero conditional variance");
    require(conditionalCovariance(0, 0) > 0 &&
        conditionalCovariance(0, 0) < covariance(0, 0),
        "native projection did not retain a genuine cross-covariance gain");
    Eigen::VectorXd fixedValue(1); fixedValue << 5;
    const auto production = zhangConditionPosteriorEffectiveIntegers(
        mean, covariance, conditioner, fixedValue);
    require(production.valid && production.conditioned &&
        (production.mean - conditionalMean).norm() < 1e-11 &&
        (production.covariance - conditionalCovariance).norm() < 1e-11,
        "production native-state conditioning differs from independent full-P oracle");
    const auto repeated = zhangConditionPosteriorEffectiveIntegers(
        production.mean, production.covariance, conditioner, fixedValue);
    require(repeated.valid && !repeated.conditioned &&
        (repeated.mean - production.mean).norm() < 1e-11 &&
        (repeated.covariance - production.covariance).norm() < 1e-11,
        "repeated current constraint shrank the native product covariance twice");

    const std::set<ZhangGraphEdge> otherSatellites = {
        {"R0", g01}, {"R0", SatSys(E_Sys::GPS, 3)}
    };
    const auto mismatch = ProductRelationBasisBuilder::build(current,
        zhangBuildSpanningTree(otherSatellites, "R0"), g01);
    require(!mismatch.valid && !mismatch.structuralMapValid &&
        mismatch.failureReason == "product_tree_target_node_or_root_mismatch",
        "malformed product map was rescued as a structural zero");

    std::cout << "epoch_product_zero_datum: " << checks << " checks passed\n";
}
