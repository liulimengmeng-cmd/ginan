// This test source is an unexecuted validation artifact until the user
// authorizes compilation/testing. Its reference is scalar Gaussian algebra,
// not another invocation of the production conditioner or target compiler.
#include "common/zhangProductConditioningReceipt.hpp"
#include "common/zhangZeroDatumNativeProduct.hpp"
#include "common/zhangR47Candidate.hpp"

#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

static int checks = 0;
static void require(bool condition, const char* reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}

// Literal serialization binds the actual fixture numbers without reusing a
// production root hasher. A mutation of ANY mean/covariance entry changes it.
static std::string rootOf(const Eigen::VectorXd& x, const Eigen::MatrixXd& p)
{
    std::ostringstream s;
    s << std::setprecision(17) << x.size() << ':' << p.rows() << ':' << p.cols();
    for (int i = 0; i < x.size(); ++i) s << '|' << x[i];
    for (int i = 0; i < p.rows(); ++i)
        for (int j = 0; j < p.cols(); ++j) s << '|' << p(i,j);
    return s.str();
}

struct GaussianOracle
{
    std::array<long double, 6> x;
    std::array<std::array<long double, 6>, 6> p;

    // Exact scalar conditional formula; no matrix inverse, eigen/SVD,
    // production conditioner, or production covariance projector is called.
    void condition(int column, long double integer)
    {
        const auto beforeX = x;
        const auto beforeP = p;
        const long double variance = beforeP[column][column];
        require(variance > 0, "independent scalar oracle has zero variance");
        for (int i = 0; i < 6; ++i)
        {
            x[i] = beforeX[i] + beforeP[i][column] / variance *
                (integer - beforeX[column]);
            for (int j = 0; j < 6; ++j)
                p[i][j] = beforeP[i][j] - beforeP[i][column] *
                    beforeP[column][j] / variance;
        }
    }

    long double mean(const std::array<long double, 6>& row) const
    {
        long double out = 0;
        for (int i = 0; i < 6; ++i) out += row[i] * x[i];
        return out;
    }
    long double covariance(const std::array<long double, 6>& left,
        const std::array<long double, 6>& right) const
    {
        long double out = 0;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) out += left[i] * p[i][j] * right[j];
        return out;
    }
};

int main()
{
    try
    {
#ifdef _WIN32
        _putenv_s("ZHANG_R51_ENABLE", "1");
        _putenv_s("ZHANG_R51_RATIO_ONLY", "1");
#else
        setenv("ZHANG_R51_ENABLE", "1", 1);
        setenv("ZHANG_R51_RATIO_ONLY", "1", 1);
#endif
        const SatSys g01(E_Sys::GPS, 1), g02(E_Sys::GPS, 2);
        const std::set<ZhangGraphEdge> edges = {
            {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
        };
        const auto basis = zhangBuildSpanningTree(edges, "R0",
            {{"R0", g01}, {"R0", g02}, {"R1", g01}});
        require(basis.connected && basis.treeEdges.size() == 3,
            "literal K2,2 S-basis is not connected");
        const ZhangGraphEdge chord{"R1", g02};
        const std::map<ZhangGraphEdge, int> versions = {
            {{"R0", g01}, 1}, {{"R0", g02}, 2},
            {{"R1", g01}, 3}, {{"R1", g02}, 4}
        };
        ZhangProductPhysicalCycleChart physicalChart;
        physicalChart.columns = 2;
        require(physicalChart.add(0, "L1C", chord, basis, versions) &&
            physicalChart.add(1, "L2W", chord, basis, versions),
            "literal current dual-frequency physical chart was rejected");
        std::map<std::string, ZhangExactInteger> physicalL1, physicalL2;
        require(physicalChart.expand({1, 0}, physicalL1) &&
            physicalChart.expand({0, 1}, physicalL2),
            "literal current physical rows failed expansion");

        auto first = ProductRelationBasisBuilder::build(
            basis, basis, g01, E_Sys::GPS, E_ObsCode::L1C);
        auto second = ProductRelationBasisBuilder::build(
            basis, basis, g01, E_Sys::GPS, E_ObsCode::L2W);
        for (auto* view : {&first, &second})
        {
            // Production compilation binds G to the full ambiguity order.
            // It has zero independent rows but preserves its named catalogue.
            view->transform = Eigen::MatrixXd::Zero(0, 2);
            view->wholeLattice = zhangCompileWholeProductLattice({{0}}, {0}, {true});
        }
        std::vector<SatSys> nativeSatellites;
        std::string reason;
        require(zhangZeroDatumNativeStructureValid(first, second, 2, nativeSatellites, reason),
            "complete dual-frequency zero correction structure was rejected");
        require(nativeSatellites == std::vector<SatSys>({g01, g02}) &&
            first.exactRank == 0 && second.exactRank == 0,
            "zero correction invented a positive integer rank");

        const long int epoch = 1405216800;
        const std::string runtime = "native-zero-test";
        auto atom = std::make_shared<ZhangIntegerDecisionProof>();
        atom->id = "CURRENT_NETWORK_NONZERO_TWO_ROWS";
        atom->originalStatement = "literal current cycles k1=2 and k2=-1";
        atom->observationProvenance = runtime + "|PPP_COMMITTED|E" + std::to_string(epoch);
        atom->conditionalFailureBound = 0;
        const ZhangDecisionProofs parents{atom};
        KFKey k1;
        k1.type = KF::AMBIGUITY; k1.str = "R1"; k1.Sat = g02;
        k1.num = static_cast<int>(E_ObsCode::L1C);
        KFKey k2 = k1; k2.num = static_cast<int>(E_ObsCode::L2W);
        const std::vector<ZhangCurrentEpochNumericConstraint> recipe = {
            {{{k1, 1}}, 2}, {{{k2, 1}}, -1}
        };

        // Order is clock,phase-L1,phase-L2,iono,k1,k2. Strict diagonal
        // dominance makes the literal symmetric prior positive definite.
        GaussianOracle oracle{
            {0.2L, -0.1L, 0.3L, 0.8L, 1.7L, -0.6L},
            {{{4,.4L,.6L,.2L,.5L,.7L},
              {.4L,3,.3L,.1L,.8L,.2L},
              {.6L,.3L,5,.5L,.4L,1},
              {.2L,.1L,.5L,2,.3L,.6L},
              {.5L,.8L,.4L,.3L,9,2},
              {.7L,.2L,1,.6L,2,16}}}
        };
        Eigen::VectorXd mean(6);
        Eigen::MatrixXd covariance(6,6);
        for (int i = 0; i < 6; ++i)
        {
            mean[i] = static_cast<double>(oracle.x[i]);
            for (int j = 0; j < 6; ++j) covariance(i,j) = static_cast<double>(oracle.p[i][j]);
        }
        const auto protectedMean = mean;
        const auto protectedCovariance = covariance;
        const std::string sourceRoot = rootOf(mean, covariance);
        const std::string ambiguityMoment = rootOf(mean.tail(2), covariance.bottomRightCorner(2,2));
        const std::string chartId = "K22_CURRENT_ORDER_L1C_L2W_V1_2_3_4";
        const std::string stateOrder = "CLOCK|L1|L2|IONO|K1|K2";
        ZhangCurrentEpochEvidence current;
        require(current.bind(runtime, E_Sys::GPS, epoch, ambiguityMoment,
            chartId, physicalChart, sourceRoot), "current physical evidence binding failed");
        current.floatPosteriorPsdValidated = true;
        current.floatPosteriorPsdRootId = sourceRoot;
        current.floatPosteriorPsdStateOrderId = stateOrder;
        require(current.add(physicalL1, 2, parents) && current.add(physicalL2, -1, parents) &&
            current.addNumericConstraint(recipe[0].rowCycles, 2, physicalL1) &&
            current.addNumericConstraint(recipe[1].rowCycles, -1, physicalL2) &&
            current.validatePhysicalNumericBindings(reason),
            "nonzero current H,b did not bind to the numerical recipe");
        const auto exactDomain = current.compiledCurrentDomain();
        require(exactDomain && exactDomain->feasible && exactDomain->basis.size() == 2 &&
            !exactDomain->proofClosure.atoms.empty(),
            "nonempty current physical integer domain or parents were lost");

        ZhangZeroDatumNativeRequest request;
        request.requested = true; request.satellites = nativeSatellites; request.conditioningRank = 2;
        ZhangR47Candidate candidate;
        candidate.sourcePosteriorId = ambiguityMoment;
        candidate.authoritativeCycleChartId = chartId;
        candidate.frontendSemanticId = "CURRENT_NATIVE_G_ZERO_COMPLETE_TWO_SIGNALS";
        candidate.stateIndices = {4,5};
        candidate.admittedStateConditioners = candidate.jointRows = {{1,0},{0,1}};
        candidate.admittedStateValues = candidate.jointValues = {2,-1};
        candidate.physicalFunctionals = {physicalL1, physicalL2};
        candidate.allDecisionParents = parents;
        candidate.riskBound = 0;
        candidate.sourceIntegerParentCount = 0;
        auto& nativeProduct = candidate.productConsequences;
        nativeProduct.zeroDatumNative = request;
        nativeProduct.reliable = nativeProduct.exactNetworkMapping =
            nativeProduct.fullJointProductMappingExact = true;
        nativeProduct.system = E_Sys::GPS;
        nativeProduct.firstObservable = E_ObsCode::L1C;
        nativeProduct.secondObservable = E_ObsCode::L2W;
        nativeProduct.referenceSatellite = g01;
        nativeProduct.coordinateSatellites = {g02};
        nativeProduct.productCoordinateDimension = 1;
        nativeProduct.networkAmbiguityDimension = 2;
        nativeProduct.networkRows = candidate.jointRows;
        nativeProduct.networkIntegers = candidate.jointValues;
        nativeProduct.fullJointProductNetworkRows = {{0,0},{0,0}};
        nativeProduct.fullJointProductAffineOffsets = {0,0};
        nativeProduct.decisionProofs = parents;
        nativeProduct.conditioningRank = 2;
        require(zhangR47CandidateContractValid(candidate),
            "native-zero current network candidate failed its native immutable contract");
        require(candidate.productConsequences.zeroDatumNative.requested &&
            zhangR47CandidateProductChartMatches(candidate, nativeProduct),
            "native producer/delivery chart did not preserve its explicit authority route");
        Eigen::MatrixXd h = Eigen::MatrixXd::Zero(2,6);
        h(0,4) = 1; h(1,5) = 1;
        Eigen::VectorXd integers(2); integers << 2, -1;
        const auto result = zhangConditionPosteriorEffectiveIntegers(mean, covariance, h, integers);
        oracle.condition(4, 2);
        oracle.condition(5, -1);
        require(result.valid && result.conditioned && result.conditionerRank == 2,
            "full current nonzero integer conditioning failed");
        for (int i = 0; i < 6; ++i)
        {
            require(std::abs(result.mean[i] - static_cast<double>(oracle.x[i])) < 1e-11,
                "full conditioned mean differs from independent scalar oracle");
            for (int j = 0; j < 6; ++j)
                require(std::abs(result.covariance(i,j) - static_cast<double>(oracle.p[i][j])) < 1e-11,
                    "full conditioned covariance/cross block differs from scalar oracle");
        }
        require((mean - protectedMean).norm() == 0 &&
            (covariance - protectedCovariance).norm() == 0,
            "private native conditioning mutated the protected FLOAT");

        // Check the complete native correction plus ionosphere, including
        // cross-frequency, clock-phase and phase-ionosphere covariance.
        const std::array<long double,6> l1{1,-1,0,-1,0,0};
        const std::array<long double,6> l2{1,0,-1,-1.6469444444444444L,0,0};
        Eigen::MatrixXd projection(2,6);
        for (int c = 0; c < 6; ++c) { projection(0,c) = l1[c]; projection(1,c) = l2[c]; }
        const Eigen::VectorXd productMean = projection * result.mean;
        const Eigen::MatrixXd productCovariance = projection * result.covariance * projection.transpose();
        const std::array<std::array<long double,6>,2> productRows{l1,l2};
        for (int r = 0; r < 2; ++r)
        {
            require(std::abs(productMean[r] - static_cast<double>(oracle.mean(productRows[r]))) < 1e-11,
                "native product mean lost a non-ambiguity state term");
            for (int s = 0; s < 2; ++s)
                require(std::abs(productCovariance(r,s) - static_cast<double>(
                    oracle.covariance(productRows[r], productRows[s]))) < 1e-11,
                    "native full joint product covariance differs from scalar oracle");
        }
        require(std::abs(result.covariance(0,1)) > 1e-4 &&
            std::abs(result.covariance(1,3)) > 1e-4 &&
            std::abs(productCovariance(0,1)) > 1e-4 && productCovariance(0,0) > 0,
            "native zero datum incorrectly zeroed real product/cross covariance");
        require(result.covariance(1,1) < covariance(1,1) &&
            result.covariance(2,2) < covariance(2,2),
            "fixture has no genuine native phase conditioning gain");

        ZhangProductConditioningIdentity actual;
        actual.runtimeId = runtime; actual.epochGpstSeconds = epoch;
        actual.posteriorId = ambiguityMoment; actual.chartId = chartId;
        actual.sourceRoot = sourceRoot; actual.resultRoot = rootOf(result.mean, result.covariance);
        actual.sourceStateOrder = actual.resultStateOrder = stateOrder;
        ZhangProductConditioningReceipt receipt;
        require(receipt.bind(actual, recipe, reason),
            "successful full-state conditioner could not issue a nonempty receipt");
        require(zhangZeroDatumNativePosteriorVerified(request, receipt, actual, recipe, reason),
            "verified native posterior was rejected");
        require(first.exactRank == 0 && second.exactRank == 0 &&
            candidate.productConsequences.certifiedPairRank == 0 &&
            candidate.productConsequences.dualFrequencyCertifiedPairs.empty() &&
            candidate.jointRows == ZhangExactMatrix({{1,0},{0,1}}),
            "native admission manufactured 0=0 pair rows or changed network H");

        for (int field = 0; field < 8; ++field)
        {
            auto changed = actual;
            switch (field)
            {
                case 0: changed.runtimeId += "-other"; break;
                case 1: changed.epochGpstSeconds += 30; break;
                case 2: changed.posteriorId += "-other"; break;
                case 3: changed.chartId += "-other-arc"; break;
                case 4: changed.sourceRoot += "-other"; break;
                case 5: changed.resultRoot += "-other"; break;
                case 6: changed.sourceStateOrder += "-other"; break;
                case 7: changed.resultStateOrder += "-other"; break;
            }
            require(!zhangZeroDatumNativePosteriorVerified(request, receipt, changed, recipe, reason),
                "wrong current epoch/root/chart/order reused a native receipt");
        }
        auto mutatedCovariance = result.covariance;
        mutatedCovariance(0,3) += .125; mutatedCovariance(3,0) += .125;
        auto changedActual = actual;
        changedActual.resultRoot = rootOf(result.mean, mutatedCovariance);
        require(!zhangZeroDatumNativePosteriorVerified(request, receipt, changedActual, recipe, reason),
            "mutated full clock-ionosphere cross block reused the old receipt");
        for (int kind = 0; kind < 7; ++kind)
        {
            auto changedRequest = request;
            auto changedRecipe = recipe;
            switch (kind)
            {
                case 0: changedRequest.requested = false; break;
                case 1: changedRequest.conditioningRank = 0; break;
                case 2: changedRequest.conditioningRank = 3; break;
                case 3: changedRecipe.clear(); break;
                case 4: changedRecipe = {{{{k1,0}},0}}; break;
                case 5: changedRecipe.push_back({recipe[0].rowCycles,3}); break;
                case 6: changedRecipe[0].rowCycles[k1] = std::numeric_limits<double>::quiet_NaN(); break;
            }
            require(!zhangZeroDatumNativePosteriorVerified(changedRequest, receipt, actual, changedRecipe, reason),
                "missing/zero/conflicting/nonfinite current network recipe authorized native AR");
        }
        auto inconsistent = candidate;
        inconsistent.jointRows.push_back({1,0}); inconsistent.jointValues.push_back(3);
        require(!zhangR47CandidateContractValid(inconsistent),
            "inconsistent integer network was accepted for native products");
        auto noParents = candidate; noParents.allDecisionParents.clear();
        require(!zhangR47CandidateContractValid(noParents),
            "nonzero rows without current decision parents authorized native products");
        for (int kind = 0; kind < 14; ++kind)
        {
            auto badCandidate = candidate;
            auto& badProduct = badCandidate.productConsequences;
            switch (kind)
            {
                case 0: badProduct.networkIntegers[0] = 3; break;
                case 1: badProduct.networkRows[0][0] = 2; break;
                case 2: badProduct.zeroDatumNative.conditioningRank = 0; break;
                case 3: badProduct.fullJointProductNetworkRows[0][0] = 1; break;
                case 4: badProduct.fullJointProductAffineOffsets[0] = 1; break;
                case 5: badProduct.fullJointProductMappingExact = false; break;
                case 6: badProduct.coordinateSatellites.clear(); break;
                case 7: badProduct.zeroDatumNative.satellites.push_back(g02); break;
                case 8: badProduct.secondObservable = badProduct.firstObservable; break;
                case 9: badProduct.jointProductRows = {{0,0}}; break;
                case 10: badProduct.conditioningOnlyRows = {{0}}; break;
                case 11: badProduct.certifiedPairRank = 1; break;
                case 12: badCandidate.newIntegerConstraints = {{1,0}};
                    badCandidate.newIntegerValues = {2}; break;
                case 13: badProduct.wideLaneProductRows = {{0}}; break;
            }
            require(!zhangR47CandidateContractValid(badCandidate),
                "malformed native metadata or fictitious pair equality passed native candidate contract");
        }
        for (int kind = 0; kind < 3; ++kind)
        {
            auto badCandidate = candidate;
            if (kind == 0)
            {
                badCandidate.jointRows.clear(); badCandidate.jointValues.clear();
            }
            else if (kind == 1)
            {
                badCandidate.jointRows = {{0,0}}; badCandidate.jointValues = {0};
            }
            else
            {
                badCandidate.jointRows = {{1,0},{1,0}}; badCandidate.jointValues = {2,2};
            }
            badCandidate.admittedStateConditioners = badCandidate.jointRows;
            badCandidate.admittedStateValues = badCandidate.jointValues;
            badCandidate.physicalFunctionals.resize(badCandidate.jointRows.size(), physicalL1);
            badCandidate.productConsequences.networkRows = badCandidate.jointRows;
            badCandidate.productConsequences.networkIntegers = badCandidate.jointValues;
            badCandidate.productConsequences.zeroDatumNative.conditioningRank = badCandidate.jointRows.size();
            require(!zhangR47CandidateContractValid(badCandidate),
                "empty, all-zero, or dependent rows were counted as native AR conditioning rank");
        }
        for (int field = 0; field < 3; ++field)
        {
            auto delivered = nativeProduct;
            if (field == 0) delivered.zeroDatumNative.requested = false;
            if (field == 1) delivered.zeroDatumNative.satellites.clear();
            if (field == 2) --delivered.zeroDatumNative.conditioningRank;
            require(!zhangR47CandidateProductChartMatches(candidate, delivered),
                "writer substituted native authority metadata after candidate freeze");
        }
        auto oldParent = std::make_shared<ZhangIntegerDecisionProof>(*atom);
        oldParent->id = "OLD_PARENT";
        oldParent->observationProvenance = runtime + "|PPP_COMMITTED|E" + std::to_string(epoch-30);
        auto wrongParentEvidence = current;
        require(!wrongParentEvidence.add(physicalL1,2,{oldParent}),
            "a previous-epoch parent was accepted as current native authority");
        auto wrongPhysical = physicalL1;
        const auto firstArc = wrongPhysical.begin()->first;
        const auto coefficient = wrongPhysical.begin()->second;
        wrongPhysical.erase(firstArc); wrongPhysical[firstArc + "_OTHER_ARC"] = coefficient;
        require(!current.add(wrongPhysical,2,parents),
            "unknown or retired physical arc was accepted as current native authority");
        for (int kind = 0; kind < 11; ++kind)
        {
            auto badFirst = first, badSecond = second;
            switch (kind)
            {
                case 0: badFirst.zeroDatumCorrection = false; break;
                case 1: badFirst.namedRelations[0].currentCycleCoefficients[0] = 1; break;
                case 2: badFirst.wholePosteriorRows = {{1,0}}; break;
                case 3: badFirst.wholeLattice.wholeAvailableRank = 1; break;
                case 4: badFirst.wholeLattice.namedRecoverable[0] = false; break;
                case 5: badSecond.observable = badFirst.observable; break;
                case 6: badFirst.structuralZeroNamedIndices.clear(); break;
                case 7: badFirst.transform = Eigen::MatrixXd::Zero(1,2); break;
                case 8: badFirst.exactRowBasis = {{1}}; break;
                case 9: badFirst.wholeLattice.namedOffsets[0] = 1; break;
                case 10: badFirst.wholeLattice.structuralIdentityValues[0] = 1; break;
            }
            std::vector<SatSys> rejectedCatalogue;
            require(!zhangZeroDatumNativeStructureValid(badFirst,badSecond,2,rejectedCatalogue,reason) &&
                rejectedCatalogue.empty(),
                "nonzero/residual/incomplete product target entered native-zero authority");
        }
        const auto repeated = zhangConditionPosteriorEffectiveIntegers(
            result.mean, result.covariance, h, integers);
        require(repeated.valid && !repeated.conditioned &&
            (repeated.mean - result.mean).norm() < 1e-11 &&
            (repeated.covariance - result.covariance).norm() < 1e-11,
            "duplicate native delivery conditioned the same network twice");
        std::cout << "native_zero_datum_posterior: " << checks << " checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
