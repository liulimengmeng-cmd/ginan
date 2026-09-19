#define BOOST_TEST_MODULE ZhangFullRankTests
#include <boost/test/unit_test.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <tuple>
#include "ambres/GNSSambres.hpp"
#include "common/algebra.hpp"
#include "common/eigenIncluder.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangIntegerAudit.hpp"
#include "common/zhangIarGainAudit.hpp"
#include "common/zhangProductRelationBasis.hpp"
#include "common/zhangProductRelationAdmission.hpp"
#include "common/zhangConflictAwareProductForest.hpp"
#include "common/zhangProductIntegerLedger.hpp"
#include "common/zhangProductGaugeCertificateLedger.hpp"
#include "common/zhangProductIntegerCandidateGenerator.hpp"
#include "common/zhangFullProductLatticeOracle.hpp"
#include "common/zhangProductRelationSolver.hpp"
#include "common/zhangIntegerSupportQuality.hpp"
#include "common/zhangIntegerSupportResidualAudit.hpp"
#include "common/zhangTargetedBesdTracker.hpp"
#include "common/zhangTheoryRegression.hpp"
#include "common/zhangPhaseContinuity.hpp"
#include "common/zhangSatelliteDatum.hpp"
#include "common/zhangPersistentProductDatum.hpp"
#include "common/zhangFullRank.hpp"
#include "common/zhangFixedLagSquareRoot.hpp"
#include "common/zhangIncrementalFixedLag.hpp"
#include "common/zhangIncrementalRawSquareRoot.hpp"
#include "common/zhangPersistentRawTargetWindow.hpp"
#include "common/zhangIntegerTargets.hpp"
#include "common/zhangLambdaBeam.hpp"
#include "common/zhangResidualStatistics.hpp"
#include "common/zhangFactorCapture.hpp"
#include "common/zhangRawFactorWindow.hpp"
#include "common/zhangUserTarget.hpp"
#include "common/zhangUserIntegerFunctional.hpp"
#include "common/zhangIfUser.hpp"
#include "common/zhangIfWideLane.hpp"
#include "common/zhangProductGaugeCompiler.hpp"
#include "common/zhangQuotientIntegerLattice.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangIncrementalSubsetNis.hpp"
#include "common/zhangHybridUserModel.hpp"
#include "common/zhangHybridService.hpp"
#include "pea/zhangPppAr.hpp"

namespace
{
double zhangTestRoundPerr(double fractional, double variance)
{
	if (variance < 1e-20) return 0;
	double probability = 0;
	const double factor = -0.25 / variance;
	for (int integer = 1; integer < 10; integer++)
	{
		probability += std::exp(
			(integer + 2 * fractional) * integer * factor);
		probability += std::exp(
			(integer - 2 * fractional) * integer * factor);
	}
	return probability / (probability + 1);
}

void makeProductGaugeEvidenceComplete(ProductGaugeCertificate& certificate,
	bool noResidualDof = false)
{
	certificate.wideLaneReliable = true;
	certificate.firstSignalReliable = true;
	certificate.exactProductLatticeMembership = true;
	certificate.jointNisPassed = true;
	certificate.cycleClosurePassed = true;
	certificate.temporalAlignmentCertified = true;
	certificate.noResidualDof = noResidualDof;
	certificate.independentSupportPaths = noResidualDof ? 0 : 1;
}

void makeFixedLagTemporalEvidence(ProductGaugeCertificate& certificate,
	const std::vector<std::size_t>& factorSequences,
	const Matrix2d& covariance)
{
	certificate.requiresIndependentTemporalEvidence = true;
	certificate.factorWindowIdentity = "TEST_RAW_FACTOR_WINDOW";
	certificate.factorSequences = factorSequences;
	certificate.factorMeasurementRows = 100 * factorSequences.size();
	certificate.temporalEvidenceCovariance = covariance;
}

BOOST_AUTO_TEST_CASE(
	routine_temporal_snapshot_capture_covers_persistent_product_coordinates)
{
	BOOST_CHECK(zhangShouldCaptureRoutineTemporalProductSnapshot(
		true, false, true, true, true, true));
	BOOST_CHECK(zhangShouldCaptureRoutineTemporalProductSnapshot(
		true, true, true, true, true, true));
	BOOST_CHECK(!zhangShouldCaptureRoutineTemporalProductSnapshot(
		false, true, true, true, true, true));
	BOOST_CHECK(!zhangShouldCaptureRoutineTemporalProductSnapshot(
		true, true, false, true, true, true));
	BOOST_CHECK(!zhangShouldCaptureRoutineTemporalProductSnapshot(
		true, true, true, false, true, true));
	BOOST_CHECK(!zhangShouldCaptureRoutineTemporalProductSnapshot(
		true, true, true, true, false, true));
	BOOST_CHECK(!zhangShouldCaptureRoutineTemporalProductSnapshot(
		true, true, true, true, true, false));
}

BOOST_AUTO_TEST_CASE(
	product_integer_pending_row_rechecks_on_later_posterior)
{
	const auto accepted = zhangRecheckProductIntegerOnPosterior(
		12.0001, 1e-4, ZhangExactInteger(12), 1e-3, 1e-6);
	BOOST_CHECK(accepted.valid);
	BOOST_CHECK(accepted.sameInteger);
	BOOST_CHECK(accepted.reliable);
	BOOST_CHECK_EQUAL(accepted.failureReason, "NONE");

	const auto changed = zhangRecheckProductIntegerOnPosterior(
		13.0001, 1e-4, ZhangExactInteger(12), 1e-3, 1e-6);
	BOOST_CHECK(changed.valid);
	BOOST_CHECK(!changed.sameInteger);
	BOOST_CHECK(!changed.reliable);
	BOOST_CHECK_EQUAL(changed.failureReason, "INTEGER_CHANGED");

	const auto uncertain = zhangRecheckProductIntegerOnPosterior(
		12.01, 0.25, ZhangExactInteger(12), 1e-3, 1e-6);
	BOOST_CHECK(uncertain.valid);
	BOOST_CHECK(uncertain.sameInteger);
	BOOST_CHECK(!uncertain.reliable);
	BOOST_CHECK_EQUAL(
		uncertain.failureReason, "FAILURE_PROBABILITY_EXCEEDED");

	const auto deterministic = zhangRecheckProductIntegerOnPosterior(
		12.0, 0.0, ZhangExactInteger(12), 1e-3, 1e-6);
	BOOST_CHECK(deterministic.valid);
	BOOST_CHECK(deterministic.deterministic);
	BOOST_CHECK(deterministic.reliable);
}

BOOST_AUTO_TEST_CASE(
	product_integer_pending_recheck_preserves_original_joint_admission_budget)
{
	ProductIntegerLedgerRow row;
	BOOST_CHECK_CLOSE(
		zhangProductLedgerPosteriorFailureBudget(row, 1e-3), 1e-3, 1e-9);

	row.admissionFailureProbabilityBound = 2e-5;
	BOOST_CHECK_CLOSE(
		zhangProductLedgerPosteriorFailureBudget(row, 1e-3), 2e-5, 1e-9);

	row.admissionFailureProbabilityBound = 0;
	BOOST_CHECK_CLOSE(
		zhangProductLedgerPosteriorFailureBudget(row, 1e-3), 1e-12, 1e-9);

	row.admissionFailureProbabilityBound =
		std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK(!std::isfinite(
		zhangProductLedgerPosteriorFailureBudget(row, 1e-3)));
}

BOOST_AUTO_TEST_CASE(
	product_integer_pending_row_matures_without_lambda_reselection)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.firstObservable = E_ObsCode::L1C;
	row.secondObservable = E_ObsCode::L2W;
	row.productRow = {1, -1};
	row.integerValue = 12;
	row.physicalExpansion = {{"GPS|L1C|R0-G01|ARC4", 1}};
	row.canonicalProductExpansion = {{"GPS|L1C|G01", 1}};
	row.phaseSegmentFingerprint = "G01|L1C|SEG7";
	row.backendBasisGeneration = 3;
	row.admissionFailureProbabilityBound = 2e-5;
	row.conditioningOnly = true;
	row.pairCertificate = false;

	const auto first = ledger.observe(100, {row}, 2);
	BOOST_REQUIRE(first.valid);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(!ledger.rows().front().certified);
	BOOST_CHECK_EQUAL(ledger.rows().front().confirmationEpochs, 1);
	BOOST_CHECK_CLOSE(ledger.rows().front().admissionFailureProbabilityBound,
		2e-5, 1e-9);

	// The next LAMBDA basis may not contain `row`; the stored physical row is
	// instead evaluated directly on the later posterior and then re-observed.
	const auto recheck = zhangRecheckProductIntegerOnPosterior(
		12.0001, 1e-4, row.integerValue, 1e-3, 1e-6);
	BOOST_REQUIRE(recheck.reliable);
	const auto second = ledger.observe(130, {ledger.rows().front()}, 2);
	BOOST_REQUIRE(second.valid);
	BOOST_CHECK(ledger.rows().front().certified);
	BOOST_CHECK_EQUAL(ledger.rows().front().confirmationEpochs, 2);
	BOOST_CHECK_CLOSE(ledger.rows().front().admissionFailureProbabilityBound,
		2e-5, 1e-9);
}

Receiver& checkpointTestReceiver()
{
	static Receiver receiver;
	receiver.id = "R0";
	return receiver;
}

struct TemporaryCheckpointFile
{
	std::filesystem::path path;

	explicit TemporaryCheckpointFile(const std::string& suffix)
	{
		const auto nonce =
			std::chrono::steady_clock::now().time_since_epoch().count();
		path = std::filesystem::temp_directory_path() /
			("ginan_e29_checkpoint_" + std::to_string(nonce) + suffix);
	}

	~TemporaryCheckpointFile()
	{
		std::error_code error;
		std::filesystem::remove(path, error);
	}
};

KFState makeCheckpointTestState()
{
	KFState state;
	KFKey receiverClock;
	receiverClock.type = KF::REC_CLOCK;
	receiverClock.str = "R0";
	receiverClock.comment = "receiver datum";
	receiverClock.estimatedTime.bigTime = 123456700.5L;
	receiverClock.rec_ptr = &checkpointTestReceiver();
	KFKey satelliteClock;
	satelliteClock.type = KF::SAT_CLOCK;
	satelliteClock.Sat = SatSys(E_Sys::GPS, 7);
	satelliteClock.comment = "satellite datum";
	satelliteClock.estimatedTime.bigTime = 123456710.25L;

	state.time.bigTime = 123456789.25L;
	state.x = VectorXd(3);
	state.x << 1, -2.5, 3.75;
	state.P = MatrixXd(3, 3);
	state.P <<
		0, 0, 0,
		0, 4, -0.25,
		0, -0.25, 9;
	state.dx = VectorXd(3);
	state.dx << 0, 0.125, -0.5;
	state.prefitRatios = VectorXd(2);
	state.prefitRatios << 1.5, 2.5;
	state.postfitRatios = VectorXd(2);
	state.postfitRatios << 0.75, 1.25;
	state.kfIndexMap = {
		{KFState::oneKey, 0},
		{receiverClock, 1},
		{satelliteClock, 2}};
	state.stateTransitionMap[satelliteClock][receiverClock][0] = -1.25;
	state.gaussMarkovTauMap[receiverClock] = 3600;
	state.gaussMarkovMuMap[receiverClock] = 0.25;
	state.procNoiseMap[satelliteClock] = 0.01;
	state.initNoiseMap[satelliteClock] = 4;
	state.sigmaMaxMap[satelliteClock] = 20;
	state.outageLimitMap[satelliteClock] = 120;
	state.exponentialNoiseMap[receiverClock] = {0.75, 45};
	state.pseudoStateMap[satelliteClock][receiverClock] = -1;
	state.pseudoParentMap[receiverClock] = satelliteClock;
	state.errorCountMap[satelliteClock] = 2;
	FilterChunk chunk;
	chunk.id = "connected-product-core";
	chunk.begX = 1;
	chunk.numX = 2;
	chunk.begH = 7;
	chunk.numH = 11;
	state.filterChunkMap["zhang"] = chunk;
	state.metaDataMap["zhang_checkpoint_runtime_id"] = "runtime-00";
	state.lsqRequired = true;
	state.sigmaPass = true;
	state.chiQCPass = true;
	state.chi2 = 12.5;
	state.dof = 7;
	state.chi2PerDof = 12.5 / 7;
	state.qc = 0.875;
	state.id = "E29-test-state";
	state.rts_basename = "e29-test";
	state.output_residuals = true;
	state.outputMongoMeasurements = true;
	state.statisticsMap["accepted"] = 17;
	state.statisticsMapSum["accepted"] = 41;
	return state;
}

ZhangCheckpointBundle makeCheckpointTestBundle()
{
	ZhangCheckpointBundle bundle;
	bundle.manifest.runtimeId = "runtime-00";
	bundle.manifest.checkpointId = "seed-00";
	bundle.manifest.parentCheckpointId = "cold-start";
	bundle.manifest.epoch = "2019-07-18T00:00:00Z";
	bundle.manifest.binarySha256 = std::string(64, 'a');
	bundle.manifest.configText = "frozen-e29-config";
	bundle.manifest.inputManifestText = "frozen-e29-input-manifest";
	bundle.manifest.configSha256 =
		zhangCheckpointSha256(bundle.manifest.configText);
	bundle.manifest.inputManifestSha256 =
		zhangCheckpointSha256(bundle.manifest.inputManifestText);
	bundle.manifest.platformFingerprint = "x86_64-linux";
	bundle.manifest.compilerFingerprint = "gcc-11.4-cxx20";
	bundle.manifest.linearAlgebraFingerprint = "eigen-3.4.1-openblas";
	bundle.manifest.endianness = "LITTLE";
	bundle.manifest.createdUtc = "2026-08-10T00:00:00Z";
	bundle.kfCore = captureZhangCheckpointKfCore(makeCheckpointTestState());
	bundle.sections["zhang.graph"] = {
		1, "pointer-free-graph-runtime", ""};
	bundle.sections["zhang.graph"].sha256 =
		zhangCheckpointSha256(bundle.sections["zhang.graph"].payload);
	return bundle;
}

ZhangCheckpointExpectations checkpointTestExpectations()
{
	ZhangCheckpointExpectations expectations;
	expectations.experimentMode = "E29_GPS_L1C_L2W_ZHANG_FULL_RANK";
	expectations.binarySha256 = std::string(64, 'a');
	expectations.configSha256 =
		zhangCheckpointSha256("frozen-e29-config");
	expectations.inputManifestSha256 =
		zhangCheckpointSha256("frozen-e29-input-manifest");
	expectations.platformFingerprint = "x86_64-linux";
	expectations.compilerFingerprint = "gcc-11.4-cxx20";
	expectations.linearAlgebraFingerprint = "eigen-3.4.1-openblas";
	expectations.endianness = "LITTLE";
	return expectations;
}

struct ZhangFormalDesign
{
    MatrixXd raw;
    MatrixXd full;
    MatrixXd nullSpace;
    MatrixXd rawToFull;
};

ZhangFormalDesign buildZhangFormalDesign()
{
    constexpr int receiverCount = 3;
    constexpr int satelliteCount = 4;
    constexpr int frequencyCount = 2;

    constexpr int rowCount = 2 * frequencyCount * receiverCount * satelliteCount;

    constexpr int rawTroposphereOffset = 0;
    constexpr int rawReceiverClockOffset = rawTroposphereOffset + receiverCount;
    constexpr int rawSatelliteClockOffset = rawReceiverClockOffset + receiverCount;
    constexpr int rawIonosphereOffset = rawSatelliteClockOffset + satelliteCount;
    constexpr int rawReceiverIfOffset =
        rawIonosphereOffset + receiverCount * satelliteCount;
    constexpr int rawSatelliteIfOffset = rawReceiverIfOffset + receiverCount;
    constexpr int rawReceiverGfOffset = rawSatelliteIfOffset + satelliteCount;
    constexpr int rawSatelliteGfOffset = rawReceiverGfOffset + receiverCount;
    constexpr int rawReceiverPhaseOffset = rawSatelliteGfOffset + satelliteCount;
    constexpr int rawSatellitePhaseOffset =
        rawReceiverPhaseOffset + frequencyCount * receiverCount;
    constexpr int rawAmbiguityOffset =
        rawSatellitePhaseOffset + frequencyCount * satelliteCount;
    constexpr int rawColumnCount =
        rawAmbiguityOffset + frequencyCount * receiverCount * satelliteCount;

    constexpr int fullTroposphereOffset = 0;
    constexpr int fullReceiverClockOffset = fullTroposphereOffset + receiverCount;
    constexpr int fullSatelliteClockOffset =
        fullReceiverClockOffset + receiverCount - 1;
    constexpr int fullIonosphereOffset = fullSatelliteClockOffset + satelliteCount;
    constexpr int fullReceiverPhaseOffset =
        fullIonosphereOffset + receiverCount * satelliteCount;
    constexpr int fullSatellitePhaseOffset =
        fullReceiverPhaseOffset + frequencyCount * (receiverCount - 1);
    constexpr int fullAmbiguityOffset =
        fullSatellitePhaseOffset + frequencyCount * satelliteCount;
    constexpr int fullColumnCount =
        fullAmbiguityOffset +
        frequencyCount * (receiverCount - 1) * (satelliteCount - 1);

    constexpr int nullity =
        2 * receiverCount + 2 * satelliteCount + 1 +
        frequencyCount * (receiverCount + satelliteCount);

    const double mu[frequencyCount] = {1.0, 1.6469444444444445};
    const double lambda[frequencyCount] = {0.190293672798365, 0.244210213424568};
    const double tropMapping[receiverCount][satelliteCount] = {
        {1.10, 1.43, 2.04, 3.12},
        {1.22, 1.67, 2.31, 2.86},
        {1.35, 1.82, 2.18, 3.57}
    };

    auto rawIonosphere = [&](int receiver, int satellite)
    {
        return rawIonosphereOffset + receiver * satelliteCount + satellite;
    };
    auto rawReceiverPhase = [&](int frequency, int receiver)
    {
        return rawReceiverPhaseOffset + frequency * receiverCount + receiver;
    };
    auto rawSatellitePhase = [&](int frequency, int satellite)
    {
        return rawSatellitePhaseOffset + frequency * satelliteCount + satellite;
    };
    auto rawAmbiguity = [&](int frequency, int receiver, int satellite)
    {
        return rawAmbiguityOffset +
               frequency * receiverCount * satelliteCount +
               receiver * satelliteCount + satellite;
    };
    auto fullReceiverClock = [&](int receiver)
    {
        return receiver == 0 ? -1 : fullReceiverClockOffset + receiver - 1;
    };
    auto fullIonosphere = [&](int receiver, int satellite)
    {
        return fullIonosphereOffset + receiver * satelliteCount + satellite;
    };
    auto fullReceiverPhase = [&](int frequency, int receiver)
    {
        return receiver == 0
                   ? -1
                   : fullReceiverPhaseOffset +
                         frequency * (receiverCount - 1) + receiver - 1;
    };
    auto fullSatellitePhase = [&](int frequency, int satellite)
    {
        return fullSatellitePhaseOffset + frequency * satelliteCount + satellite;
    };
    auto fullAmbiguity = [&](int frequency, int receiver, int satellite)
    {
        if (receiver == 0 || satellite == 0)
        {
            return -1;
        }

        return fullAmbiguityOffset +
               frequency * (receiverCount - 1) * (satelliteCount - 1) +
               (receiver - 1) * (satelliteCount - 1) + satellite - 1;
    };

    ZhangFormalDesign result;
    result.raw = MatrixXd::Zero(rowCount, rawColumnCount);
    result.full = MatrixXd::Zero(rowCount, fullColumnCount);
    result.nullSpace = MatrixXd::Zero(rawColumnCount, nullity);
    result.rawToFull = MatrixXd::Zero(fullColumnCount, rawColumnCount);

    int row = 0;
    for (int frequency = 0; frequency < frequencyCount; frequency++)
    {
        for (int receiver = 0; receiver < receiverCount; receiver++)
        {
            for (int satellite = 0; satellite < satelliteCount; satellite++)
            {
                // Raw code equation.
                result.raw(row, rawTroposphereOffset + receiver) =
                    tropMapping[receiver][satellite];
                result.raw(row, rawReceiverClockOffset + receiver) = +1;
                result.raw(row, rawSatelliteClockOffset + satellite) = -1;
                result.raw(row, rawIonosphere(receiver, satellite)) = +mu[frequency];
                result.raw(row, rawReceiverIfOffset + receiver) = +1;
                result.raw(row, rawSatelliteIfOffset + satellite) = -1;
                result.raw(row, rawReceiverGfOffset + receiver) = +mu[frequency];
                result.raw(row, rawSatelliteGfOffset + satellite) = -mu[frequency];

                // Full-rank code equation.
                result.full(row, fullTroposphereOffset + receiver) =
                    tropMapping[receiver][satellite];
                if (int column = fullReceiverClock(receiver); column >= 0)
                {
                    result.full(row, column) = +1;
                }
                result.full(row, fullSatelliteClockOffset + satellite) = -1;
                result.full(row, fullIonosphere(receiver, satellite)) = +mu[frequency];
                row++;

                // Raw phase equation in Ginan's satellite-phase sign convention.
                result.raw(row, rawTroposphereOffset + receiver) =
                    tropMapping[receiver][satellite];
                result.raw(row, rawReceiverClockOffset + receiver) = +1;
                result.raw(row, rawSatelliteClockOffset + satellite) = -1;
                result.raw(row, rawIonosphere(receiver, satellite)) = -mu[frequency];
                result.raw(row, rawReceiverPhase(frequency, receiver)) = +1;
                result.raw(row, rawSatellitePhase(frequency, satellite)) = +1;
                result.raw(row, rawAmbiguity(frequency, receiver, satellite)) =
                    lambda[frequency];

                // Full-rank phase equation.
                result.full(row, fullTroposphereOffset + receiver) =
                    tropMapping[receiver][satellite];
                if (int column = fullReceiverClock(receiver); column >= 0)
                {
                    result.full(row, column) = +1;
                }
                result.full(row, fullSatelliteClockOffset + satellite) = -1;
                result.full(row, fullIonosphere(receiver, satellite)) = -mu[frequency];
                if (int column = fullReceiverPhase(frequency, receiver); column >= 0)
                {
                    result.full(row, column) = +1;
                }
                result.full(row, fullSatellitePhase(frequency, satellite)) = +1;
                if (int column = fullAmbiguity(frequency, receiver, satellite);
                    column >= 0)
                {
                    result.full(row, column) = lambda[frequency];
                }
                row++;
            }
        }
    }

    int nullColumn = 0;

    // Receiver IF code-bias directions.
    for (int receiver = 0; receiver < receiverCount; receiver++)
    {
        result.nullSpace(rawReceiverClockOffset + receiver, nullColumn) = -1;
        result.nullSpace(rawReceiverIfOffset + receiver, nullColumn) = +1;
        for (int frequency = 0; frequency < frequencyCount; frequency++)
        {
            result.nullSpace(rawReceiverPhase(frequency, receiver), nullColumn) = +1;
        }
        nullColumn++;
    }

    // Satellite IF code-bias directions.
    for (int satellite = 0; satellite < satelliteCount; satellite++)
    {
        result.nullSpace(rawSatelliteClockOffset + satellite, nullColumn) = -1;
        result.nullSpace(rawSatelliteIfOffset + satellite, nullColumn) = +1;
        for (int frequency = 0; frequency < frequencyCount; frequency++)
        {
            result.nullSpace(rawSatellitePhase(frequency, satellite), nullColumn) = -1;
        }
        nullColumn++;
    }

    // Receiver GF code-bias directions.
    for (int receiver = 0; receiver < receiverCount; receiver++)
    {
        result.nullSpace(rawReceiverGfOffset + receiver, nullColumn) = +1;
        for (int satellite = 0; satellite < satelliteCount; satellite++)
        {
            result.nullSpace(rawIonosphere(receiver, satellite), nullColumn) = -1;
        }
        for (int frequency = 0; frequency < frequencyCount; frequency++)
        {
            result.nullSpace(rawReceiverPhase(frequency, receiver), nullColumn) =
                -mu[frequency];
        }
        nullColumn++;
    }

    // Satellite GF code-bias directions.
    for (int satellite = 0; satellite < satelliteCount; satellite++)
    {
        result.nullSpace(rawSatelliteGfOffset + satellite, nullColumn) = +1;
        for (int receiver = 0; receiver < receiverCount; receiver++)
        {
            result.nullSpace(rawIonosphere(receiver, satellite), nullColumn) = +1;
        }
        for (int frequency = 0; frequency < frequencyCount; frequency++)
        {
            result.nullSpace(rawSatellitePhase(frequency, satellite), nullColumn) =
                +mu[frequency];
        }
        nullColumn++;
    }

    // Common receiver/satellite clock datum.
    for (int receiver = 0; receiver < receiverCount; receiver++)
    {
        result.nullSpace(rawReceiverClockOffset + receiver, nullColumn) = +1;
    }
    for (int satellite = 0; satellite < satelliteCount; satellite++)
    {
        result.nullSpace(rawSatelliteClockOffset + satellite, nullColumn) = +1;
    }
    nullColumn++;

    // Receiver phase-bias/ambiguity directions.
    for (int frequency = 0; frequency < frequencyCount; frequency++)
    {
        for (int receiver = 0; receiver < receiverCount; receiver++)
        {
            result.nullSpace(rawReceiverPhase(frequency, receiver), nullColumn) =
                lambda[frequency];
            for (int satellite = 0; satellite < satelliteCount; satellite++)
            {
                result.nullSpace(
                    rawAmbiguity(frequency, receiver, satellite),
                    nullColumn
                ) = -1;
            }
            nullColumn++;
        }
    }

    // Satellite phase-bias/ambiguity directions.
    for (int frequency = 0; frequency < frequencyCount; frequency++)
    {
        for (int satellite = 0; satellite < satelliteCount; satellite++)
        {
            result.nullSpace(rawSatellitePhase(frequency, satellite), nullColumn) =
                lambda[frequency];
            for (int receiver = 0; receiver < receiverCount; receiver++)
            {
                result.nullSpace(
                    rawAmbiguity(frequency, receiver, satellite),
                    nullColumn
                ) = -1;
            }
            nullColumn++;
        }
    }

    BOOST_REQUIRE_EQUAL(nullColumn, nullity);

    auto addTransform = [&](int fullRow, int rawColumn, double coefficient)
    {
        if (fullRow >= 0)
        {
            result.rawToFull(fullRow, rawColumn) += coefficient;
        }
    };

    for (int receiver = 0; receiver < receiverCount; receiver++)
    {
        addTransform(
            fullTroposphereOffset + receiver,
            rawTroposphereOffset + receiver,
            +1
        );
    }

    for (int receiver = 1; receiver < receiverCount; receiver++)
    {
        int fullRow = fullReceiverClock(receiver);
        addTransform(fullRow, rawReceiverClockOffset + receiver, +1);
        addTransform(fullRow, rawReceiverIfOffset + receiver, +1);
        addTransform(fullRow, rawReceiverClockOffset, -1);
        addTransform(fullRow, rawReceiverIfOffset, -1);
    }

    for (int satellite = 0; satellite < satelliteCount; satellite++)
    {
        int fullRow = fullSatelliteClockOffset + satellite;
        addTransform(fullRow, rawSatelliteClockOffset + satellite, +1);
        addTransform(fullRow, rawSatelliteIfOffset + satellite, +1);
        addTransform(fullRow, rawReceiverClockOffset, -1);
        addTransform(fullRow, rawReceiverIfOffset, -1);
    }

    for (int receiver = 0; receiver < receiverCount; receiver++)
    {
        for (int satellite = 0; satellite < satelliteCount; satellite++)
        {
            int fullRow = fullIonosphere(receiver, satellite);
            addTransform(fullRow, rawIonosphere(receiver, satellite), +1);
            addTransform(fullRow, rawReceiverGfOffset + receiver, +1);
            addTransform(fullRow, rawSatelliteGfOffset + satellite, -1);
        }
    }

    for (int frequency = 0; frequency < frequencyCount; frequency++)
    {
        for (int receiver = 1; receiver < receiverCount; receiver++)
        {
            int fullRow = fullReceiverPhase(frequency, receiver);
            addTransform(fullRow, rawReceiverPhase(frequency, receiver), +1);
            addTransform(fullRow, rawReceiverIfOffset + receiver, -1);
            addTransform(fullRow, rawReceiverGfOffset + receiver, +mu[frequency]);
            addTransform(fullRow, rawReceiverPhase(frequency, 0), -1);
            addTransform(fullRow, rawReceiverIfOffset, +1);
            addTransform(fullRow, rawReceiverGfOffset, -mu[frequency]);
            addTransform(
                fullRow,
                rawAmbiguity(frequency, receiver, 0),
                +lambda[frequency]
            );
            addTransform(
                fullRow,
                rawAmbiguity(frequency, 0, 0),
                -lambda[frequency]
            );
        }

        for (int satellite = 0; satellite < satelliteCount; satellite++)
        {
            int fullRow = fullSatellitePhase(frequency, satellite);
            addTransform(fullRow, rawSatellitePhase(frequency, satellite), +1);
            addTransform(fullRow, rawSatelliteIfOffset + satellite, +1);
            addTransform(
                fullRow,
                rawSatelliteGfOffset + satellite,
                -mu[frequency]
            );
            addTransform(fullRow, rawReceiverPhase(frequency, 0), +1);
            addTransform(fullRow, rawReceiverIfOffset, -1);
            addTransform(fullRow, rawReceiverGfOffset, +mu[frequency]);
            addTransform(
                fullRow,
                rawAmbiguity(frequency, 0, satellite),
                +lambda[frequency]
            );
        }

        for (int receiver = 1; receiver < receiverCount; receiver++)
        {
            for (int satellite = 1; satellite < satelliteCount; satellite++)
            {
                int fullRow = fullAmbiguity(frequency, receiver, satellite);
                addTransform(
                    fullRow,
                    rawAmbiguity(frequency, receiver, satellite),
                    +1
                );
                addTransform(
                    fullRow,
                    rawAmbiguity(frequency, 0, satellite),
                    -1
                );
                addTransform(
                    fullRow,
                    rawAmbiguity(frequency, receiver, 0),
                    -1
                );
                addTransform(fullRow, rawAmbiguity(frequency, 0, 0), +1);
            }
        }
    }

    return result;
}
}  // namespace

BOOST_AUTO_TEST_CASE(reference_row_and_column_are_ambiguity_s_bases)
{
    std::vector<E_ObsCode> baseline = {E_ObsCode::L1C, E_ObsCode::L2W};

    BOOST_CHECK(!zhangFullRankRetainsAmbiguity(
        "ZIM2",
        SatSys(E_Sys::GPS, 8),
        E_ObsCode::L1C,
        baseline,
        "ZIM2",
        "G08"
    ));
    BOOST_CHECK(!zhangFullRankRetainsAmbiguity(
        "FFMJ",
        SatSys(E_Sys::GPS, 8),
        E_ObsCode::L1C,
        baseline,
        "ZIM2",
        "G08"
    ));
    BOOST_CHECK(zhangFullRankRetainsAmbiguity(
        "FFMJ",
        SatSys(E_Sys::GPS, 12),
        E_ObsCode::L2W,
        baseline,
        "ZIM2",
        "G08"
    ));
    BOOST_CHECK(!zhangFullRankRetainsAmbiguity(
        "FFMJ",
        SatSys(E_Sys::GPS, 12),
        E_ObsCode::L5Q,
        baseline,
        "ZIM2",
        "G08"
    ));
}

BOOST_AUTO_TEST_CASE(raw_model_null_space_matches_the_29_zhang_s_basis_directions)
{
    ZhangFormalDesign design = buildZhangFormalDesign();

    Eigen::FullPivLU<MatrixXd> rawDecomposition(design.raw);
    rawDecomposition.setThreshold(1e-11);
    Eigen::FullPivLU<MatrixXd> nullDecomposition(design.nullSpace);
    nullDecomposition.setThreshold(1e-11);

    BOOST_CHECK_EQUAL(design.raw.rows(), 48);
    BOOST_CHECK_EQUAL(design.raw.cols(), 74);
    BOOST_CHECK_EQUAL(rawDecomposition.rank(), 45);
    BOOST_CHECK_EQUAL(rawDecomposition.dimensionOfKernel(), 29);
    BOOST_CHECK_EQUAL(nullDecomposition.rank(), 29);
    BOOST_CHECK_SMALL((design.raw * design.nullSpace).norm(), 1e-12);
    BOOST_CHECK_SMALL((design.rawToFull * design.nullSpace).norm(), 1e-12);
    BOOST_CHECK_SMALL(
        (design.full * design.rawToFull - design.raw).norm(),
        1e-12
    );
}

BOOST_AUTO_TEST_CASE(raw_generalized_inverse_and_full_rank_solution_are_equivalent)
{
    ZhangFormalDesign design = buildZhangFormalDesign();

    VectorXd rawState = VectorXd::LinSpaced(design.raw.cols(), -2.5, 3.5);
    VectorXd observations = design.raw * rawState;
    for (int row = 0; row < observations.rows(); row++)
    {
        observations(row) += 1e-3 * std::sin(0.37 * row);
    }

    Eigen::CompleteOrthogonalDecomposition<MatrixXd> rawSolver(design.raw);
    rawSolver.setThreshold(1e-11);
    Eigen::CompleteOrthogonalDecomposition<MatrixXd> fullSolver(design.full);
    fullSolver.setThreshold(1e-11);

    VectorXd rawSolution = rawSolver.solve(observations);
    VectorXd fullSolution = fullSolver.solve(observations);
    VectorXd transformedRawSolution = design.rawToFull * rawSolution;

    VectorXd rawPrediction = design.raw * rawSolution;
    VectorXd fullPrediction = design.full * fullSolution;

    BOOST_CHECK_SMALL((transformedRawSolution - fullSolution).norm(), 1e-9);
    BOOST_CHECK_SMALL((rawPrediction - fullPrediction).norm(), 1e-10);
    BOOST_CHECK_SMALL(
        ((observations - rawPrediction) - (observations - fullPrediction)).norm(),
        1e-10
    );
}

BOOST_AUTO_TEST_CASE(code_phase_ionosphere_float_design_is_full_column_rank)
{
    constexpr int receiverCount = 3;
    constexpr int satelliteCount = 4;
    constexpr int frequencyCount = 2;

    constexpr int troposphereCount = receiverCount;
    constexpr int receiverClockCount = receiverCount - 1;
    constexpr int satelliteClockCount = satelliteCount;
    constexpr int ionosphereCount = receiverCount * satelliteCount;
    constexpr int receiverPhaseCount = frequencyCount * (receiverCount - 1);
    constexpr int satellitePhaseCount = frequencyCount * satelliteCount;
    constexpr int ambiguityCount =
        frequencyCount * (receiverCount - 1) * (satelliteCount - 1);

    constexpr int columnCount =
        troposphereCount + receiverClockCount + satelliteClockCount + ionosphereCount +
        receiverPhaseCount + satellitePhaseCount + ambiguityCount;
    constexpr int rowCount = 2 * frequencyCount * receiverCount * satelliteCount;

    MatrixXd design = MatrixXd::Zero(rowCount, columnCount);

    int receiverClockOffset = troposphereCount;
    int satelliteClockOffset = receiverClockOffset + receiverClockCount;
    int ionosphereOffset = satelliteClockOffset + satelliteClockCount;
    int receiverPhaseOffset = ionosphereOffset + ionosphereCount;
    int satellitePhaseOffset = receiverPhaseOffset + receiverPhaseCount;
    int ambiguityOffset = satellitePhaseOffset + satellitePhaseCount;

    const double mu[frequencyCount] = {1.0, 1.6469444444444445};
    const double lambda[frequencyCount] = {0.190293672798365, 0.244210213424568};
    const double tropMapping[receiverCount][satelliteCount] = {
        {1.10, 1.43, 2.04, 3.12},
        {1.22, 1.67, 2.31, 2.86},
        {1.35, 1.82, 2.18, 3.57}
    };

    int row = 0;
    for (int frequency = 0; frequency < frequencyCount; frequency++)
    {
        for (int receiver = 0; receiver < receiverCount; receiver++)
        {
            for (int satellite = 0; satellite < satelliteCount; satellite++)
            {
                int ionosphereColumn =
                    ionosphereOffset + receiver * satelliteCount + satellite;

                design(row, receiver) = tropMapping[receiver][satellite];
                if (receiver > 0)
                {
                    design(row, receiverClockOffset + receiver - 1) = 1;
                }
                design(row, satelliteClockOffset + satellite) = -1;
                design(row, ionosphereColumn) = mu[frequency];
                row++;

                design(row, receiver) = tropMapping[receiver][satellite];
                if (receiver > 0)
                {
                    design(row, receiverClockOffset + receiver - 1) = 1;
                }
                design(row, satelliteClockOffset + satellite) = -1;
                design(row, ionosphereColumn) = -mu[frequency];

                if (receiver > 0)
                {
                    int receiverPhaseColumn =
                        receiverPhaseOffset + frequency * (receiverCount - 1) + receiver - 1;
                    design(row, receiverPhaseColumn) = 1;
                }

                int satellitePhaseColumn =
                    satellitePhaseOffset + frequency * satelliteCount + satellite;
                design(row, satellitePhaseColumn) = 1;

                if (receiver > 0 && satellite > 0)
                {
                    int ambiguityColumn =
                        ambiguityOffset +
                        frequency * (receiverCount - 1) * (satelliteCount - 1) +
                        (receiver - 1) * (satelliteCount - 1) + satellite - 1;
                    design(row, ambiguityColumn) = lambda[frequency];
                }
                row++;
            }
        }
    }

    BOOST_REQUIRE_EQUAL(row, rowCount);
    BOOST_REQUIRE_EQUAL(design.cols(), columnCount);

    Eigen::FullPivLU<MatrixXd> decomposition(design);
    decomposition.setThreshold(1e-11);

    BOOST_CHECK_EQUAL(decomposition.rank(), columnCount);
    BOOST_CHECK_EQUAL(decomposition.dimensionOfKernel(), 0);
}

BOOST_AUTO_TEST_CASE(reference_change_preserves_code_and_phase_observables)
{
    constexpr int receiverCount = 3;
    constexpr int satelliteCount = 4;
    constexpr int frequencyCount = 2;

    const double receiverClock[receiverCount] = {13.2, -7.1, 4.3};
    const double satelliteClock[satelliteCount] = {-2.4, 8.7, 1.1, -5.2};
    const double ionosphere[receiverCount][satelliteCount] = {
        {2.1, 3.2, 1.7, 4.0},
        {2.5, 3.6, 2.0, 4.4},
        {1.9, 3.0, 1.5, 3.8}
    };
    const double receiverPhase[frequencyCount][receiverCount] = {
        {0.12, -0.31, 0.25},
        {-0.21, 0.17, 0.38}
    };
    const double satellitePhase[frequencyCount][satelliteCount] = {
        {0.45, -0.22, 0.31, -0.16},
        {-0.37, 0.28, 0.19, -0.42}
    };
    const int ambiguity[frequencyCount][receiverCount][satelliteCount] = {
        {
            {10, 14, 21, 8},
            {17, 5, 12, 23},
            {6, 19, 15, 11}
        },
        {
            {32, 27, 18, 41},
            {25, 39, 31, 16},
            {44, 22, 36, 29}
        }
    };
    const double mu[frequencyCount] = {1.0, 1.6469444444444445};
    const double lambda[frequencyCount] = {0.190293672798365, 0.244210213424568};

    auto reconstructed = [&](int referenceReceiver, int referenceSatellite)
    {
        VectorXd values(2 * frequencyCount * receiverCount * satelliteCount);
        int row = 0;

        for (int frequency = 0; frequency < frequencyCount; frequency++)
        {
            for (int receiver = 0; receiver < receiverCount; receiver++)
            {
                for (int satellite = 0; satellite < satelliteCount; satellite++)
                {
                    double receiverClockEstimate =
                        receiverClock[receiver] - receiverClock[referenceReceiver];
                    double satelliteClockEstimate =
                        satelliteClock[satellite] - receiverClock[referenceReceiver];

                    values(row++) =
                        receiverClockEstimate -
                        satelliteClockEstimate +
                        mu[frequency] * ionosphere[receiver][satellite];

                    double receiverPhaseEstimate = 0;
                    if (receiver != referenceReceiver)
                    {
                        receiverPhaseEstimate =
                            receiverPhase[frequency][receiver] -
                            receiverPhase[frequency][referenceReceiver] +
                            lambda[frequency] *
                                (ambiguity[frequency][receiver][referenceSatellite] -
                                 ambiguity[frequency][referenceReceiver][referenceSatellite]);
                    }

                    double satellitePhaseEstimate =
                        satellitePhase[frequency][satellite] +
                        receiverPhase[frequency][referenceReceiver] +
                        lambda[frequency] *
                            ambiguity[frequency][referenceReceiver][satellite];

                    double ambiguityEstimate = 0;
                    if (receiver != referenceReceiver && satellite != referenceSatellite)
                    {
                        ambiguityEstimate =
                            ambiguity[frequency][receiver][satellite] -
                            ambiguity[frequency][referenceReceiver][satellite] -
                            ambiguity[frequency][receiver][referenceSatellite] +
                            ambiguity[frequency][referenceReceiver][referenceSatellite];
                    }

                    values(row++) =
                        receiverClockEstimate -
                        satelliteClockEstimate -
                        mu[frequency] * ionosphere[receiver][satellite] +
                        receiverPhaseEstimate +
                        satellitePhaseEstimate +
                        lambda[frequency] * ambiguityEstimate;
                }
            }
        }

        return values;
    };

    VectorXd datumA = reconstructed(0, 0);
    VectorXd datumB = reconstructed(1, 2);

    BOOST_CHECK_SMALL((datumA - datumB).norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(reference_change_transforms_full_covariance_without_information_loss)
{
    constexpr int receiverCount = 3;
    constexpr int satelliteCount = 3;
    constexpr double wavelength = 0.190293672798365;

    // Per datum: two receiver clocks, three satellite clocks, two receiver phase states,
    // three satellite phase states, and four DD ambiguities.
    constexpr int stateCount =
        (receiverCount - 1) + satelliteCount +
        (receiverCount - 1) + satelliteCount +
        (receiverCount - 1) * (satelliteCount - 1);

    auto transform = [&](int oldReceiver, int oldSatellite, int newReceiver, int newSatellite)
    {
        MatrixXd T = MatrixXd::Zero(stateCount, stateCount);

        auto recClockIndex = [&](int receiver, int reference)
        {
            int index = 0;
            for (int r = 0; r < receiverCount; r++)
            {
                if (r == reference)
                    continue;
                if (r == receiver)
                    return index;
                index++;
            }
            return -1;
        };
        auto satClockIndex = [&](int satellite)
        {
            return receiverCount - 1 + satellite;
        };
        auto recPhaseIndex = [&](int receiver, int reference)
        {
            int offset = receiverCount - 1 + satelliteCount;
            int index = 0;
            for (int r = 0; r < receiverCount; r++)
            {
                if (r == reference)
                    continue;
                if (r == receiver)
                    return offset + index;
                index++;
            }
            return -1;
        };
        auto satPhaseIndex = [&](int satellite)
        {
            return 2 * (receiverCount - 1) + satelliteCount + satellite;
        };
        auto ambiguityIndex = [&](int receiver, int satellite, int referenceReceiver, int referenceSatellite)
        {
            if (receiver == referenceReceiver || satellite == referenceSatellite)
                return -1;

            int offset = 2 * (receiverCount - 1) + 2 * satelliteCount;
            int index = 0;
            for (int r = 0; r < receiverCount; r++)
                for (int s = 0; s < satelliteCount; s++)
                {
                    if (r == referenceReceiver || s == referenceSatellite)
                        continue;
                    if (r == receiver && s == satellite)
                        return offset + index;
                    index++;
                }
            return -1;
        };
        auto add = [&](int row, int column, double value)
        {
            if (column >= 0)
                T(row, column) += value;
        };

        for (int receiver = 0; receiver < receiverCount; receiver++)
        {
            if (receiver == newReceiver)
                continue;

            int row = recClockIndex(receiver, newReceiver);
            add(row, recClockIndex(receiver, oldReceiver), +1);
            add(row, recClockIndex(newReceiver, oldReceiver), -1);
        }

        for (int satellite = 0; satellite < satelliteCount; satellite++)
        {
            int row = satClockIndex(satellite);
            add(row, satClockIndex(satellite), +1);
            add(row, recClockIndex(newReceiver, oldReceiver), -1);
        }

        for (int receiver = 0; receiver < receiverCount; receiver++)
        {
            if (receiver == newReceiver)
                continue;

            int row = recPhaseIndex(receiver, newReceiver);
            add(row, recPhaseIndex(receiver, oldReceiver), +1);
            add(row, recPhaseIndex(newReceiver, oldReceiver), -1);
            add(
                row,
                ambiguityIndex(receiver, newSatellite, oldReceiver, oldSatellite),
                +wavelength
            );
            add(
                row,
                ambiguityIndex(newReceiver, newSatellite, oldReceiver, oldSatellite),
                -wavelength
            );
        }

        for (int satellite = 0; satellite < satelliteCount; satellite++)
        {
            int row = satPhaseIndex(satellite);
            add(row, satPhaseIndex(satellite), +1);
            add(row, recPhaseIndex(newReceiver, oldReceiver), +1);
            add(
                row,
                ambiguityIndex(newReceiver, satellite, oldReceiver, oldSatellite),
                +wavelength
            );
        }

        for (int receiver = 0; receiver < receiverCount; receiver++)
            for (int satellite = 0; satellite < satelliteCount; satellite++)
            {
                if (receiver == newReceiver || satellite == newSatellite)
                    continue;

                int row = ambiguityIndex(receiver, satellite, newReceiver, newSatellite);
                add(
                    row,
                    ambiguityIndex(receiver, satellite, oldReceiver, oldSatellite),
                    +1
                );
                add(
                    row,
                    ambiguityIndex(newReceiver, satellite, oldReceiver, oldSatellite),
                    -1
                );
                add(
                    row,
                    ambiguityIndex(receiver, newSatellite, oldReceiver, oldSatellite),
                    -1
                );
                add(
                    row,
                    ambiguityIndex(newReceiver, newSatellite, oldReceiver, oldSatellite),
                    +1
                );
            }

        return T;
    };

    MatrixXd forward = transform(0, 0, 1, 2);
    MatrixXd reverse = transform(1, 2, 0, 0);

    BOOST_CHECK_SMALL((reverse * forward - MatrixXd::Identity(stateCount, stateCount)).norm(), 1e-12);

    auto observationDesign = [&](int referenceReceiver, int referenceSatellite)
    {
        MatrixXd design =
            MatrixXd::Zero(2 * receiverCount * satelliteCount, stateCount);

        auto recClockIndex = [&](int receiver)
        {
            int index = 0;
            for (int r = 0; r < receiverCount; r++)
            {
                if (r == referenceReceiver)
                    continue;
                if (r == receiver)
                    return index;
                index++;
            }
            return -1;
        };
        auto satClockIndex = [&](int satellite)
        {
            return receiverCount - 1 + satellite;
        };
        auto recPhaseIndex = [&](int receiver)
        {
            int offset = receiverCount - 1 + satelliteCount;
            int index = 0;
            for (int r = 0; r < receiverCount; r++)
            {
                if (r == referenceReceiver)
                    continue;
                if (r == receiver)
                    return offset + index;
                index++;
            }
            return -1;
        };
        auto satPhaseIndex = [&](int satellite)
        {
            return 2 * (receiverCount - 1) + satelliteCount + satellite;
        };
        auto ambiguityIndex = [&](int receiver, int satellite)
        {
            if (receiver == referenceReceiver || satellite == referenceSatellite)
                return -1;

            int offset = 2 * (receiverCount - 1) + 2 * satelliteCount;
            int index = 0;
            for (int r = 0; r < receiverCount; r++)
                for (int s = 0; s < satelliteCount; s++)
                {
                    if (r == referenceReceiver || s == referenceSatellite)
                        continue;
                    if (r == receiver && s == satellite)
                        return offset + index;
                    index++;
                }
            return -1;
        };

        int row = 0;
        for (int receiver = 0; receiver < receiverCount; receiver++)
            for (int satellite = 0; satellite < satelliteCount; satellite++)
            {
                if (int column = recClockIndex(receiver); column >= 0)
                {
                    design(row, column) = +1;
                }
                design(row, satClockIndex(satellite)) = -1;
                row++;

                if (int column = recClockIndex(receiver); column >= 0)
                {
                    design(row, column) = +1;
                }
                design(row, satClockIndex(satellite)) = -1;
                if (int column = recPhaseIndex(receiver); column >= 0)
                {
                    design(row, column) = +1;
                }
                design(row, satPhaseIndex(satellite)) = +1;
                if (int column = ambiguityIndex(receiver, satellite); column >= 0)
                {
                    design(row, column) = wavelength;
                }
                row++;
            }

        return design;
    };

    MatrixXd oldDesign = observationDesign(0, 0);
    MatrixXd newDesign = observationDesign(1, 2);
    BOOST_CHECK_SMALL((newDesign * forward - oldDesign).norm(), 1e-12);

    MatrixXd generator = MatrixXd::Random(stateCount, stateCount);
    MatrixXd oldCovariance =
        generator * generator.transpose() + 0.1 * MatrixXd::Identity(stateCount, stateCount);
    MatrixXd newCovariance = forward * oldCovariance * forward.transpose();
    MatrixXd recoveredCovariance = reverse * newCovariance * reverse.transpose();

    BOOST_CHECK_SMALL((newCovariance - newCovariance.transpose()).norm(), 1e-12);
    BOOST_CHECK_SMALL((recoveredCovariance - oldCovariance).norm(), 1e-10);

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(newCovariance);
    BOOST_REQUIRE_EQUAL(eigenSolver.info(), Eigen::Success);
    BOOST_CHECK_GT(eigenSolver.eigenvalues().minCoeff(), 0);

    VectorXd oldState = VectorXd::LinSpaced(stateCount, -0.8, 1.3);
    VectorXd newState = forward * oldState;
    VectorXd measurementNoise =
        VectorXd::LinSpaced(oldDesign.rows(), -2e-3, 3e-3);
    VectorXd observations = oldDesign * oldState + measurementNoise;

    VectorXd oldPrediction = oldDesign * oldState;
    VectorXd newPrediction = newDesign * newState;
    VectorXd oldInnovation = observations - oldPrediction;
    VectorXd newInnovation = observations - newPrediction;

    MatrixXd measurementCovariance =
        0.04 * MatrixXd::Identity(oldDesign.rows(), oldDesign.rows());
    MatrixXd oldInnovationCovariance =
        oldDesign * oldCovariance * oldDesign.transpose() + measurementCovariance;
    MatrixXd newInnovationCovariance =
        newDesign * newCovariance * newDesign.transpose() + measurementCovariance;

    BOOST_CHECK_SMALL((newPrediction - oldPrediction).norm(), 1e-12);
    BOOST_CHECK_SMALL((newInnovation - oldInnovation).norm(), 1e-12);
    BOOST_CHECK_SMALL(
        (newInnovationCovariance - oldInnovationCovariance).norm() /
            oldInnovationCovariance.norm(),
        1e-12
    );
}

namespace
{
MatrixXd zhangGraphPhaseDesign(const ZhangGraphBasis& basis, double wavelength = 0.190293672798365)
{
    std::vector<std::string> receivers;
    for (const auto& receiver : basis.receivers)
    {
        if (receiver != basis.rootReceiver)
        {
            receivers.push_back(receiver);
        }
    }
    std::vector<SatSys> satellites(basis.satellites.begin(), basis.satellites.end());

    std::vector<ZhangGraphEdge> cycleEdges;
    for (const auto& edge : basis.edges)
    {
        if (!basis.isTreeEdge(edge.receiver, edge.satellite))
        {
            cycleEdges.push_back(edge);
        }
    }

    std::map<std::string, int> receiverColumn;
    std::map<SatSys, int> satelliteColumn;
    std::map<ZhangGraphEdge, int> cycleColumn;

    int column = 0;
    for (const auto& receiver : receivers)
    {
        receiverColumn[receiver] = column++;
    }
    for (const auto& satellite : satellites)
    {
        satelliteColumn[satellite] = column++;
    }
    for (const auto& edge : cycleEdges)
    {
        cycleColumn[edge] = column++;
    }

    MatrixXd design = MatrixXd::Zero(basis.edges.size(), column);
    int row = 0;
    for (const auto& edge : basis.edges)
    {
        if (edge.receiver != basis.rootReceiver)
        {
            design(row, receiverColumn.at(edge.receiver)) = +1;
        }
        design(row, satelliteColumn.at(edge.satellite)) = +1;

        auto cycleIt = cycleColumn.find(edge);
        if (cycleIt != cycleColumn.end())
        {
            design(row, cycleIt->second) = wavelength;
        }
        row++;
    }
    return design;
}

MatrixXd zhangCycleMatrix(const ZhangGraphBasis& basis)
{
    std::vector<ZhangGraphEdge> edges(basis.edges.begin(), basis.edges.end());
    std::map<ZhangGraphEdge, int> edgeColumn;
    for (int index = 0; index < edges.size(); index++)
    {
        edgeColumn[edges[index]] = index;
    }

    std::vector<ZhangGraphEdge> nonTreeEdges;
    for (const auto& edge : edges)
    {
        if (!basis.isTreeEdge(edge.receiver, edge.satellite))
        {
            nonTreeEdges.push_back(edge);
        }
    }

    MatrixXd cycles = MatrixXd::Zero(nonTreeEdges.size(), edges.size());
    for (int row = 0; row < nonTreeEdges.size(); row++)
    {
        auto cycle = zhangFundamentalCycle(basis, nonTreeEdges[row]);
        for (const auto& [edge, coefficient] : cycle)
        {
            cycles(row, edgeColumn.at(edge)) = coefficient;
        }
    }
    return cycles;
}

std::set<ZhangGraphEdge> sparseConnectedGraph(int receiverCount, int satelliteCount, int seed)
{
    std::set<ZhangGraphEdge> edges;

    edges.insert({"R0", SatSys(E_Sys::GPS, 1)});
    for (int receiver = 1; receiver < receiverCount; receiver++)
    {
        int previousSatellite = 1 + (receiver - 1) % satelliteCount;
        int nextSatellite     = 1 + receiver % satelliteCount;
        edges.insert(
            {"R" + std::to_string(receiver), SatSys(E_Sys::GPS, previousSatellite)}
        );
        edges.insert(
            {"R" + std::to_string(receiver), SatSys(E_Sys::GPS, nextSatellite)}
        );
    }
    for (int satellite = receiverCount + 1; satellite <= satelliteCount; satellite++)
    {
        int receiver = (satellite + seed) % receiverCount;
        edges.insert(
            {"R" + std::to_string(receiver), SatSys(E_Sys::GPS, satellite)}
        );
    }

    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> receiverDistribution(0, receiverCount - 1);
    std::uniform_int_distribution<int> satelliteDistribution(1, satelliteCount);
    for (int extra = 0; extra < receiverCount + satelliteCount; extra++)
    {
        edges.insert(
            {"R" + std::to_string(receiverDistribution(generator)),
             SatSys(E_Sys::GPS, satelliteDistribution(generator))}
        );
    }
    return edges;
}
}  // namespace

BOOST_AUTO_TEST_CASE(random_sparse_connected_graphs_have_full_rank_tree_coordinates)
{
    for (int seed = 1; seed <= 30; seed++)
    {
        int receiverCount  = 3 + seed % 6;
        int satelliteCount = 4 + seed % 9;
        auto edges = sparseConnectedGraph(receiverCount, satelliteCount, seed);

        ZhangGraphBasis basis = zhangBuildSpanningTree(edges, "R0");
        BOOST_REQUIRE(basis.connected);
        BOOST_CHECK_EQUAL(
            basis.treeEdges.size(),
            basis.receivers.size() + basis.satellites.size() - 1
        );

        MatrixXd design = zhangGraphPhaseDesign(basis);
        Eigen::FullPivLU<MatrixXd> decomposition(design);
        decomposition.setThreshold(1e-11);
        BOOST_CHECK_EQUAL(decomposition.rank(), design.cols());
        BOOST_CHECK_EQUAL(design.rows(), design.cols());
    }
}

BOOST_AUTO_TEST_CASE(product_receiver_core_spans_all_satellites_with_redundancy)
{
    const SatSys g1(E_Sys::GPS, 1);
    const SatSys g2(E_Sys::GPS, 2);
    const SatSys g3(E_Sys::GPS, 3);
    const SatSys g4(E_Sys::GPS, 4);
    const SatSys g5(E_Sys::GPS, 5);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g1}, {"R0", g2},
        {"R1", g1}, {"R1", g3}, {"R1", g4},
        {"R2", g2}, {"R2", g4}, {"R2", g5},
        {"R3", g3}, {"R3", g5},
        {"R4", g1}, {"R4", g2}, {"R4", g3}, {"R4", g4}, {"R4", g5},
        {"R5", g1}, {"R5", g2}
    };
    const auto core = zhangBuildProductReceiverCore(
        edges, "R0", {"R3", "R5"}, 2);
    BOOST_REQUIRE(core.connected);
    BOOST_CHECK_EQUAL(core.satellites.size(), 5);
    BOOST_CHECK_EQUAL(core.minimumSatelliteSupport, 2);
    BOOST_CHECK(core.receivers.find("R0") != core.receivers.end());
    BOOST_CHECK(core.receivers.find("R3") != core.receivers.end());
    // A prior receiver that contributes no remaining support deficit is not
    // forced into the new core; this is the controlled-retirement invariant.
    BOOST_CHECK(core.receivers.find("R5") == core.receivers.end());
    BOOST_CHECK_LT(core.receivers.size(), 6);
    std::map<SatSys, int> support;
    for (const auto& edge : core.edges)
    {
        support[edge.satellite]++;
    }
    for (const auto& satellite : core.satellites)
    {
        BOOST_CHECK_GE(support[satellite], 2);
    }
    BOOST_CHECK(zhangBuildSpanningTree(core.edges, "R0").connected);
}

BOOST_AUTO_TEST_CASE(integer_support_core_excludes_unqualified_arcs_but_keeps_float_graph)
{
	const SatSys g1(E_Sys::GPS, 1);
	const SatSys g2(E_Sys::GPS, 2);
	const std::set<ZhangGraphEdge> fullFloatEdges = {
		{"R0", g1}, {"R0", g2}, {"R1", g1}, {"R1", g2}};
	ZhangIntegerArcQuality good;
	good.ageEpochs = 20;
	good.observations = 20;
	good.phaseResidualRms = 0.01;
	good.codeResidualRms = 1;
	good.phaseResidualMad = 0.01;
	good.codeResidualMad = 1;
	good.elevationScore = 0.5;
	good.whitenedResidualScore = 1;
	std::map<ZhangGraphEdge, ZhangIntegerArcQuality> quality;
	for (const auto& edge : fullFloatEdges)
		quality[edge] = good;
	quality[{"R1", g2}].phaseResidualMad = 1;

	const auto core = zhangBuildIntegerSupportCore(
		fullFloatEdges, "R0", {}, 1, quality);
	BOOST_REQUIRE(core.valid);
	BOOST_CHECK_EQUAL(core.qualifiedEdges.size(), 3);
	BOOST_CHECK_EQUAL(core.rejectedEdges.at({"R1", g2}),
		"PHASE_MAD_GATE_FAILED");
	BOOST_CHECK(core.receiverCore.edges.find({"R1", g2}) ==
		core.receiverCore.edges.end());
	// The filtering result is product-only; it must never mutate or redefine
	// the authoritative FLOAT observation graph.
	BOOST_CHECK_EQUAL(fullFloatEdges.size(), 4);
	BOOST_CHECK(zhangBuildSpanningTree(fullFloatEdges, "R0").connected);
}

BOOST_AUTO_TEST_CASE(risk_aware_s_basis_rejects_bad_integer_support_arc)
{
	const SatSys g1(E_Sys::GPS, 1);
	const SatSys g2(E_Sys::GPS, 2);
	const std::set<ZhangGraphEdge> edges = {
		{"R0", g1}, {"R1", g1}, {"R0", g2}, {"R1", g2}};
	ZhangIntegerArcQuality good;
	good.ageEpochs = 20; good.observations = 20;
	good.phaseResidualRms = 0.01; good.codeResidualRms = 0.5;
	good.phaseResidualMad = 0.01; good.codeResidualMad = 0.5;
	good.elevationScore = 0.5; good.whitenedResidualScore = 1;
	auto poor = good;
	poor.phaseResidualRms = 1; // fails the strict support gate.
	std::map<ZhangGraphEdge, ZhangIntegerArcQuality> quality;
	quality[{"R0", g1}] = good;
	quality[{"R1", g1}] = poor;
	quality[{"R0", g2}] = good;
	quality[{"R1", g2}] = good;
	const auto basis = zhangBuildRiskAwareSpanningTree(
		edges, "R0", {{"R1", g1}}, quality);
	BOOST_REQUIRE(basis.connected);
	BOOST_CHECK(basis.treeEdges.find({"R1", g1}) == basis.treeEdges.end());
}

BOOST_AUTO_TEST_CASE(rooted_product_tree_limits_nonroot_satellite_path_load)
{
    const SatSys g1(E_Sys::GPS, 1);
    const SatSys g2(E_Sys::GPS, 2);
    const SatSys g3(E_Sys::GPS, 3);
    const SatSys g4(E_Sys::GPS, 4);
    const std::set<ZhangGraphEdge> edges = {
        {"R0", g1}, {"R0", g2}, {"R0", g3},
        {"R1", g1}, {"R1", g2},
        {"R2", g2}, {"R2", g3}, {"R2", g4}
    };
    const std::set<ZhangGraphEdge> chainPreferred = {
        {"R0", g1}, {"R1", g1}, {"R1", g2},
        {"R2", g2}, {"R2", g3}, {"R2", g4}
    };
    const auto kruskal = zhangBuildSpanningTree(
        edges, "R0", chainPreferred);
    const auto rooted = zhangBuildRootedProductTree(
        edges, "R0", chainPreferred);
    BOOST_REQUIRE(kruskal.connected);
    BOOST_REQUIRE(rooted.connected);

    auto vulnerableMaximum = [](const ZhangGraphBasis& basis)
    {
        int maximum = 0;
        for (const auto& [edge, load] :
             zhangProductTreeSatellitePathLoads(basis))
        {
            if (edge.receiver != basis.rootReceiver)
            {
                maximum = std::max(maximum, load);
            }
        }
        return maximum;
    };
    BOOST_CHECK_EQUAL(vulnerableMaximum(kruskal), 3);
    BOOST_CHECK_EQUAL(vulnerableMaximum(rooted), 1);
    BOOST_CHECK_LT(vulnerableMaximum(rooted), vulnerableMaximum(kruskal));
}

BOOST_AUTO_TEST_CASE(two_spanning_trees_preserve_observations_and_integer_cycle_lattice)
{
    std::set<ZhangGraphEdge> edges = {
        {"R0", SatSys(E_Sys::GPS, 1)},
        {"R0", SatSys(E_Sys::GPS, 2)},
        {"R1", SatSys(E_Sys::GPS, 1)},
        {"R1", SatSys(E_Sys::GPS, 3)},
        {"R1", SatSys(E_Sys::GPS, 4)},
        {"R2", SatSys(E_Sys::GPS, 2)},
        {"R2", SatSys(E_Sys::GPS, 4)},
        {"R2", SatSys(E_Sys::GPS, 5)},
        {"R3", SatSys(E_Sys::GPS, 3)},
        {"R3", SatSys(E_Sys::GPS, 5)}
    };

    std::set<ZhangGraphEdge> preferredA = {
        {"R0", SatSys(E_Sys::GPS, 1)},
        {"R0", SatSys(E_Sys::GPS, 2)},
        {"R1", SatSys(E_Sys::GPS, 1)},
        {"R1", SatSys(E_Sys::GPS, 3)},
        {"R1", SatSys(E_Sys::GPS, 4)},
        {"R2", SatSys(E_Sys::GPS, 2)},
        {"R2", SatSys(E_Sys::GPS, 5)},
        {"R3", SatSys(E_Sys::GPS, 3)}
    };
    std::set<ZhangGraphEdge> preferredB = {
        {"R0", SatSys(E_Sys::GPS, 2)},
        {"R1", SatSys(E_Sys::GPS, 1)},
        {"R1", SatSys(E_Sys::GPS, 3)},
        {"R2", SatSys(E_Sys::GPS, 2)},
        {"R2", SatSys(E_Sys::GPS, 4)},
        {"R2", SatSys(E_Sys::GPS, 5)},
        {"R3", SatSys(E_Sys::GPS, 3)},
        {"R3", SatSys(E_Sys::GPS, 5)}
    };

    ZhangGraphBasis basisA = zhangBuildSpanningTree(edges, "R0", preferredA);
    ZhangGraphBasis basisB = zhangBuildSpanningTree(edges, "R0", preferredB);
    BOOST_REQUIRE(basisA.connected);
    BOOST_REQUIRE(basisB.connected);
    BOOST_CHECK(basisA.treeEdges != basisB.treeEdges);

    MatrixXd designA = zhangGraphPhaseDesign(basisA);
    MatrixXd designB = zhangGraphPhaseDesign(basisB);
    VectorXd observations = VectorXd::LinSpaced(edges.size(), -1.7, 2.3);
    VectorXd stateA = designA.fullPivLu().solve(observations);
    VectorXd stateB = designB.fullPivLu().solve(observations);
    BOOST_CHECK_SMALL((designA * stateA - observations).norm(), 1e-12);
    BOOST_CHECK_SMALL((designB * stateB - observations).norm(), 1e-12);

    MatrixXd treeTransform = designB.inverse() * designA;
    BOOST_CHECK_SMALL((designB * treeTransform - designA).norm(), 1e-12);

    MatrixXd generator = MatrixXd::Random(designA.cols(), designA.cols());
    MatrixXd covarianceA =
        generator * generator.transpose() +
        0.1 * MatrixXd::Identity(designA.cols(), designA.cols());
    MatrixXd covarianceB = treeTransform * covarianceA * treeTransform.transpose();
    VectorXd transformedState = treeTransform * stateA;
    VectorXd noise = VectorXd::LinSpaced(observations.size(), -2e-3, 3e-3);
    VectorXd measured = observations + noise;
    MatrixXd measurementCovariance =
        0.02 * MatrixXd::Identity(observations.size(), observations.size());

    BOOST_CHECK_SMALL((designB * transformedState - designA * stateA).norm(), 1e-12);
    BOOST_CHECK_SMALL(
        ((measured - designB * transformedState) - (measured - designA * stateA)).norm(),
        1e-12
    );
    MatrixXd innovationCovarianceA =
        designA * covarianceA * designA.transpose() + measurementCovariance;
    MatrixXd innovationCovarianceB =
        designB * covarianceB * designB.transpose() + measurementCovariance;
    BOOST_CHECK_SMALL(
        (innovationCovarianceB - innovationCovarianceA).norm() /
            innovationCovarianceA.norm(),
        1e-12
    );

    MatrixXd cyclesA = zhangCycleMatrix(basisA);
    MatrixXd cyclesB = zhangCycleMatrix(basisB);
    BOOST_REQUIRE_EQUAL(cyclesA.rows(), cyclesB.rows());

    MatrixXd basisTransform =
        cyclesB * cyclesA.transpose() * (cyclesA * cyclesA.transpose()).inverse();
    MatrixXd integerTransform = basisTransform.array().round().matrix();
    BOOST_CHECK_SMALL((basisTransform - integerTransform).norm(), 1e-12);
    BOOST_CHECK_SMALL((cyclesB - integerTransform * cyclesA).norm(), 1e-12);
    BOOST_CHECK_CLOSE(std::abs(integerTransform.determinant()), 1.0, 1e-10);

    VectorXd rawIntegers(edges.size());
    for (int index = 0; index < rawIntegers.size(); index++)
    {
        rawIntegers(index) = (7 * index + 3) % 19 - 9;
    }
    VectorXd integersA = cyclesA * rawIntegers;
    VectorXd integersB = cyclesB * rawIntegers;
    BOOST_CHECK_SMALL((integersA.array() - integersA.array().round()).matrix().norm(), 1e-12);
    BOOST_CHECK_SMALL((integersB - integerTransform * integersA).norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(tree_edge_failure_exchanges_basis_without_losing_rank)
{
    auto edges = sparseConnectedGraph(5, 7, 42);
    ZhangGraphBasis oldBasis = zhangBuildSpanningTree(edges, "R0");
    BOOST_REQUIRE(oldBasis.connected);

    ZhangGraphEdge failed = *oldBasis.treeEdges.rbegin();
    std::set<ZhangGraphEdge> remaining = edges;
    remaining.erase(failed);
    ZhangGraphBasis newBasis =
        zhangBuildSpanningTree(remaining, "R0", oldBasis.treeEdges);

    BOOST_REQUIRE(newBasis.connected);
    BOOST_CHECK(newBasis.treeEdges.find(failed) == newBasis.treeEdges.end());
    BOOST_CHECK(newBasis.treeEdges != oldBasis.treeEdges);

    MatrixXd design = zhangGraphPhaseDesign(newBasis);
    Eigen::FullPivLU<MatrixXd> decomposition(design);
    decomposition.setThreshold(1e-11);
    BOOST_CHECK_EQUAL(decomposition.rank(), design.cols());
}

BOOST_AUTO_TEST_CASE(
    tree_slip_pivots_on_surviving_represented_graph_before_retire)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g02(E_Sys::GPS, 2);
    const std::set<ZhangGraphEdge> represented = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}};
    const std::set<ZhangGraphEdge> preferred = {
        {"R0", g01}, {"R1", g01}, {"R1", g02}};
    const ZhangGraphBasis oldBasis = zhangBuildSpanningTree(
        represented, "R0", preferred);
    BOOST_REQUIRE(oldBasis.connected);
    const ZhangGraphEdge slipped{"R1", g01};
    BOOST_REQUIRE(oldBasis.treeEdges.count(slipped));

    const auto plan = zhangPlanPivotBeforeRetire(
        oldBasis, represented, {slipped});
    BOOST_REQUIRE_MESSAGE(plan.connected, plan.failureReason);
    BOOST_CHECK(plan.required);
    BOOST_CHECK(!plan.survivingRepresentedEdges.count(slipped));
    BOOST_CHECK(!plan.replacementBasis.treeEdges.count(slipped));
    BOOST_CHECK_EQUAL(plan.replacementEdges.size(), 1);
    BOOST_CHECK(plan.replacementEdges.count({"R0", g02}));
    BOOST_CHECK(plan.replacementBasis.receivers == oldBasis.receivers);
    BOOST_CHECK(plan.replacementBasis.satellites == oldBasis.satellites);

    // The rectangular post-event transform is determined entirely from the
    // surviving pre-event observables.  It preserves their means and full
    // covariance while intentionally dropping the slipped arc coordinate.
    const MatrixXd oldDesign = zhangGraphPhaseDesign(oldBasis);
    const MatrixXd newDesign = zhangGraphPhaseDesign(plan.replacementBasis);
    BOOST_REQUIRE_EQUAL(oldDesign.rows(), 4);
    BOOST_REQUIRE_EQUAL(newDesign.rows(), 3);
    MatrixXd survivingOldDesign(3, oldDesign.cols());
    int sourceRow = 0;
    int destinationRow = 0;
    for (const auto& edge : oldBasis.edges)
    {
        if (!(edge == slipped))
            survivingOldDesign.row(destinationRow++) = oldDesign.row(sourceRow);
        sourceRow++;
    }
    BOOST_REQUIRE_EQUAL(destinationRow, 3);
    const MatrixXd transform = newDesign.inverse() * survivingOldDesign;
    const VectorXd oldState = VectorXd::LinSpaced(oldDesign.cols(), -2.0, 3.0);
    MatrixXd generator = MatrixXd::Random(oldDesign.cols(), oldDesign.cols());
    const MatrixXd oldCovariance = generator * generator.transpose() +
        0.2 * MatrixXd::Identity(oldDesign.cols(), oldDesign.cols());
    const VectorXd newState = transform * oldState;
    const MatrixXd newCovariance =
        transform * oldCovariance * transform.transpose();
    BOOST_CHECK_SMALL(
        (newDesign * newState - survivingOldDesign * oldState).norm(), 1e-12);
    const MatrixXd oldObservableCovariance =
        survivingOldDesign * oldCovariance * survivingOldDesign.transpose();
    const MatrixXd newObservableCovariance =
        newDesign * newCovariance * newDesign.transpose();
    BOOST_CHECK_SMALL(
        (newObservableCovariance - oldObservableCovariance).norm() /
            oldObservableCovariance.norm(),
        1e-12);
}

BOOST_AUTO_TEST_CASE(
    tree_slip_allows_component_split_only_when_surviving_graph_disconnects)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g02(E_Sys::GPS, 2);
    const std::set<ZhangGraphEdge> represented = {
        {"R0", g01}, {"R1", g01}, {"R1", g02}};
    const ZhangGraphBasis oldBasis =
        zhangBuildSpanningTree(represented, "R0");
    BOOST_REQUIRE(oldBasis.connected);
    const ZhangGraphEdge bridge{"R1", g01};
    BOOST_REQUIRE(oldBasis.treeEdges.count(bridge));

    const auto plan = zhangPlanPivotBeforeRetire(
        oldBasis, represented, {bridge});
    BOOST_CHECK(plan.required);
    BOOST_CHECK(!plan.connected);
    BOOST_CHECK_EQUAL(
        plan.failureReason,
        "SURVIVING_REPRESENTED_GRAPH_DISCONNECTED");
    BOOST_CHECK(plan.replacementBasis.receivers != oldBasis.receivers ||
                plan.replacementBasis.satellites != oldBasis.satellites);
}

BOOST_AUTO_TEST_CASE(
    tree_slip_joint_pivot_rejects_single_frequency_replacement)
{
    const ZhangGraphEdge rootG01{"R0", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge remoteG01{"R1", SatSys(E_Sys::GPS, 1)};
    const ZhangGraphEdge remoteG02{"R1", SatSys(E_Sys::GPS, 2)};
    const ZhangGraphEdge rootG02{"R0", SatSys(E_Sys::GPS, 2)};
    const std::set<ZhangGraphEdge> l1Edges = {
        rootG01, remoteG01, remoteG02, rootG02};
    const std::set<ZhangGraphEdge> l2Edges = {
        rootG01, remoteG01, remoteG02};

    const auto joint = zhangJointObservableEdgeIntersection({l1Edges, l2Edges});
    BOOST_CHECK_EQUAL(joint.size(), 3);
    BOOST_CHECK(!joint.count(rootG02));

    const ZhangGraphBasis oldBasis = zhangBuildSpanningTree(
        {rootG01, remoteG01, remoteG02}, "R0");
    BOOST_REQUIRE(oldBasis.connected);

    const auto abstractPlan = zhangPlanPivotBeforeRetire(
        oldBasis, l1Edges, {remoteG01});
    BOOST_REQUIRE_MESSAGE(abstractPlan.connected, abstractPlan.failureReason);
    BOOST_CHECK(abstractPlan.replacementEdges.count(rootG02));

    const auto statePlan = zhangPlanPivotBeforeRetire(
        oldBasis, joint, {remoteG01});
    BOOST_CHECK(!statePlan.connected);
    BOOST_CHECK_EQUAL(
        statePlan.failureReason,
        "SURVIVING_REPRESENTED_GRAPH_DISCONNECTED");
}

BOOST_AUTO_TEST_CASE(
    p3_experiment3_fault_injection_preserves_connected_float_posterior)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g02(E_Sys::GPS, 2);
    const std::set<ZhangGraphEdge> represented = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}};
    const std::set<ZhangGraphEdge> preferred = {
        {"R0", g01}, {"R1", g01}, {"R1", g02}};
    const ZhangGraphBasis oldBasis = zhangBuildSpanningTree(
        represented, "R0", preferred);
    BOOST_REQUIRE(oldBasis.connected);
    const ZhangGraphEdge failed{"R1", g01};
    BOOST_REQUIRE(oldBasis.treeEdges.count(failed));

    // These labels exercise the unified transaction contract.  Event type
    // changes arc-retirement semantics, not the exact S-transform on the
    // common surviving represented graph.
    const std::vector<std::string> injectedEvents = {
        "L1_ONLY_SLIP", "L2_ONLY_SLIP", "DUAL_FREQUENCY_SLIP",
        "STATION_QC_REMOVAL", "OBSERVATION_LOSS_30S",
        "OBSERVATION_LOSS_60S", "OBSERVATION_LOSS_120S",
        "RECEIVER_EXIT", "SATELLITE_SET", "SATELLITE_RISE",
        "SIGNAL_UNAVAILABLE"};
    const MatrixXd oldDesign = zhangGraphPhaseDesign(oldBasis);
    MatrixXd survivingOldDesign(3, oldDesign.cols());
    int oldRow = 0;
    int survivingRow = 0;
    for (const auto& edge : oldBasis.edges)
    {
        if (!(edge == failed))
            survivingOldDesign.row(survivingRow++) = oldDesign.row(oldRow);
        oldRow++;
    }
    BOOST_REQUIRE_EQUAL(survivingRow, 3);
    const VectorXd oldState = VectorXd::LinSpaced(
        oldDesign.cols(), -1.75, 2.25);
    MatrixXd generator = MatrixXd::Random(oldDesign.cols(), oldDesign.cols());
    const MatrixXd oldCovariance = generator * generator.transpose() +
        0.1 * MatrixXd::Identity(oldDesign.cols(), oldDesign.cols());

    for (const auto& event : injectedEvents)
    {
        BOOST_TEST_CONTEXT("fault=" << event)
        {
            const auto plan = zhangPlanPivotBeforeRetire(
                oldBasis, represented, {failed});
            BOOST_REQUIRE_MESSAGE(plan.connected, plan.failureReason);
            const MatrixXd newDesign =
                zhangGraphPhaseDesign(plan.replacementBasis);
            const MatrixXd transform =
                newDesign.inverse() * survivingOldDesign;
            const VectorXd newState = transform * oldState;
            const MatrixXd newCovariance =
                transform * oldCovariance * transform.transpose();
            BOOST_CHECK_SMALL(
                (newDesign * newState - survivingOldDesign * oldState).norm(),
                1e-12);
            const MatrixXd before = survivingOldDesign * oldCovariance *
                survivingOldDesign.transpose();
            const MatrixXd after = newDesign * newCovariance *
                newDesign.transpose();
            BOOST_CHECK_SMALL((after - before).norm() / before.norm(), 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE(
    p3_experiment3_multi_edge_split_loses_only_cross_component_rank)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g02(E_Sys::GPS, 2);
    const SatSys g03(E_Sys::GPS, 3);
    const SatSys g04(E_Sys::GPS, 4);
    const ZhangGraphEdge bridge{"R1", g03};
    const std::set<ZhangGraphEdge> represented = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02},
        bridge,
        {"R2", g03}, {"R2", g04}, {"R3", g03}, {"R3", g04}};
    const ZhangGraphBasis oldBasis = zhangBuildSpanningTree(
        represented, "R0");
    BOOST_REQUIRE(oldBasis.connected);
    BOOST_REQUIRE(oldBasis.treeEdges.count(bridge));

    const auto plan = zhangPlanPivotBeforeRetire(
        oldBasis, represented, {bridge});
    BOOST_CHECK(!plan.connected);
    BOOST_CHECK_EQUAL(
        plan.failureReason,
        "SURVIVING_REPRESENTED_GRAPH_DISCONNECTED");

    // Both 2x2 components retain their own fundamental cycle.  Only a
    // cross-component gauge direction is unavailable after the bridge loss.
    const int oldCycleRank = static_cast<int>(represented.size()) - 8 + 1;
    const int survivingCycleRank =
        static_cast<int>(plan.survivingRepresentedEdges.size()) - 8 + 2;
    BOOST_CHECK_EQUAL(oldCycleRank, 2);
    BOOST_CHECK_EQUAL(survivingCycleRank, oldCycleRank);
}

BOOST_AUTO_TEST_CASE(represented_edges_are_preferred_for_replacement_tree)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    std::set<ZhangGraphEdge> represented = {
        {"R1", g01}, {"R1", g02}
    };

    ZhangGraphBasis basis = zhangBuildSpanningTree(
        edges,
        "R0",
        {},
        {},
        represented
    );

    BOOST_REQUIRE(basis.connected);
    BOOST_CHECK(basis.treeEdges.find({"R1", g01}) != basis.treeEdges.end());
    BOOST_CHECK(basis.treeEdges.find({"R1", g02}) != basis.treeEdges.end());
}

BOOST_AUTO_TEST_CASE(longer_continuous_arc_breaks_historical_edge_ties)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    std::map<ZhangGraphEdge, int> persistence = {
        {{"R0", g01}, 5},
        {{"R0", g02}, 4},
        {{"R1", g01}, 100},
        {{"R1", g02}, 90}
    };

    ZhangGraphBasis basis = zhangBuildSpanningTree(
        edges,
        "R0",
        {},
        {},
        edges,
        persistence
    );

    BOOST_REQUIRE(basis.connected);
    BOOST_CHECK(basis.treeEdges.find({"R1", g01}) != basis.treeEdges.end());
    BOOST_CHECK(basis.treeEdges.find({"R1", g02}) != basis.treeEdges.end());
}

BOOST_AUTO_TEST_CASE(disconnected_graph_is_detected_and_root_component_isolated)
{
    std::set<ZhangGraphEdge> edges = {
        {"R0", SatSys(E_Sys::GPS, 1)},
        {"R1", SatSys(E_Sys::GPS, 1)},
        {"R2", SatSys(E_Sys::GPS, 8)},
        {"R2", SatSys(E_Sys::GPS, 9)}
    };

    ZhangGraphBasis forest = zhangBuildSpanningTree(edges, "R0");
    BOOST_CHECK(!forest.connected);
    BOOST_CHECK_EQUAL(forest.componentCount, 2);

    auto rootEdges = zhangRootComponentEdges(edges, "R0");
    BOOST_CHECK_EQUAL(rootEdges.size(), 2);
    BOOST_CHECK(rootEdges.find({"R2", SatSys(E_Sys::GPS, 8)}) == rootEdges.end());

    ZhangGraphBasis rootBasis = zhangBuildSpanningTree(rootEdges, "R0");
    BOOST_CHECK(rootBasis.connected);
}

BOOST_AUTO_TEST_CASE(canonical_integer_coordinates_close_exactly_on_two_by_two_graph)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    ZhangGraphBasis basis = zhangBuildSpanningTree(edges, "R0");
    BOOST_REQUIRE(basis.connected);

    ZhangCanonicalIntegerAudit audit = zhangCanonicalIntegerAudit(basis);
    BOOST_REQUIRE(audit.valid);
    BOOST_REQUIRE_EQUAL(audit.treeEdges.size(), 3);
    BOOST_REQUIRE_EQUAL(audit.chordEdges.size(), 1);
    BOOST_CHECK(zhangExactAbs(zhangExactDeterminant(audit.canonicalToArc)) == 1);
    BOOST_REQUIRE_EQUAL(audit.satelliteDatumSingleDifferences.size(), 1);
    BOOST_REQUIRE_EQUAL(audit.satelliteFixQuotient.size(), 1);
    BOOST_CHECK(
        std::all_of(
            audit.satelliteFixQuotient.front().begin(),
            audit.satelliteFixQuotient.front().end(),
            [](const auto& value) { return value == 0; }
        )
    );

    ZhangExactVector canonical = {2, -1, 3, 5};
    ZhangExactVector raw =
        zhangExactMatrixTimesColumn(audit.canonicalToArc, canonical);
    std::vector<ZhangGraphEdge> arcs = audit.treeEdges;
    arcs.insert(arcs.end(), audit.chordEdges.begin(), audit.chordEdges.end());
    std::map<ZhangGraphEdge, std::size_t> arcIndex;
    for (std::size_t index = 0; index < arcs.size(); index++)
    {
        arcIndex[arcs[index]] = index;
    }
    ZhangExactInteger recoveredCycle = 0;
    for (const auto& [edge, coefficient] :
         zhangFundamentalCycle(basis, audit.chordEdges.front()))
    {
        recoveredCycle += coefficient * raw[arcIndex.at(edge)];
    }
    BOOST_CHECK(recoveredCycle == canonical.back());
}

BOOST_AUTO_TEST_CASE(canonical_audit_rejects_stale_arcs_outside_active_tree)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    std::set<ZhangGraphEdge> activeEdges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    ZhangGraphBasis activeBasis = zhangBuildSpanningTree(activeEdges, "R0");
    BOOST_REQUIRE(activeBasis.connected);

    ZhangGraphBasis stateTransformBasis = activeBasis;
    stateTransformBasis.edges.insert({"STALE", SatSys(E_Sys::GPS, 9)});
    ZhangCanonicalIntegerAudit audit =
        zhangCanonicalIntegerAudit(stateTransformBasis);
    BOOST_CHECK(!audit.valid);
    BOOST_CHECK(audit.failureReason.find("missing_chord_endpoint") == 0);

    ZhangCanonicalIntegerAudit activeAudit =
        zhangCanonicalIntegerAudit(activeBasis);
    BOOST_CHECK(activeAudit.valid);
}

BOOST_AUTO_TEST_CASE(global_scale_canonical_audit_uses_sparse_structure)
{
    std::set<ZhangGraphEdge> edges;
    for (int receiver = 0; receiver < 17; receiver++)
    {
        for (int satellite = 1; satellite <= 17; satellite++)
        {
            edges.insert(
                {"R" + std::to_string(receiver), SatSys(E_Sys::GPS, satellite)}
            );
        }
    }
    ZhangGraphBasis basis = zhangBuildSpanningTree(edges, "R0");
    BOOST_REQUIRE(basis.connected);

    ZhangCanonicalIntegerAudit audit = zhangCanonicalIntegerAudit(basis);
    BOOST_REQUIRE(audit.valid);
    BOOST_CHECK(!audit.denseCanonicalMaterialised);
    BOOST_CHECK(audit.canonicalToArc.empty());
    BOOST_CHECK(!audit.canonicalToArcFingerprint.empty());
    BOOST_CHECK_EQUAL(audit.treeEdges.size(), 33);
    BOOST_CHECK_EQUAL(audit.chordEdges.size(), 256);
}

BOOST_AUTO_TEST_CASE(sparse_tree_exchange_is_an_exact_unimodular_integer_transition)
{
    std::set<ZhangGraphEdge> edges = {
        {"R0", SatSys(E_Sys::GPS, 1)}, {"R0", SatSys(E_Sys::GPS, 2)},
        {"R1", SatSys(E_Sys::GPS, 1)}, {"R1", SatSys(E_Sys::GPS, 3)},
        {"R2", SatSys(E_Sys::GPS, 2)}, {"R2", SatSys(E_Sys::GPS, 3)},
        {"R0", SatSys(E_Sys::GPS, 3)}, {"R2", SatSys(E_Sys::GPS, 1)}
    };
    std::set<ZhangGraphEdge> preferredA = {
        {"R0", SatSys(E_Sys::GPS, 1)}, {"R0", SatSys(E_Sys::GPS, 2)},
        {"R0", SatSys(E_Sys::GPS, 3)}, {"R1", SatSys(E_Sys::GPS, 1)},
        {"R2", SatSys(E_Sys::GPS, 2)}
    };
    std::set<ZhangGraphEdge> preferredB = {
        {"R0", SatSys(E_Sys::GPS, 1)}, {"R1", SatSys(E_Sys::GPS, 1)},
        {"R1", SatSys(E_Sys::GPS, 3)}, {"R2", SatSys(E_Sys::GPS, 2)},
        {"R2", SatSys(E_Sys::GPS, 3)}
    };
    ZhangGraphBasis basisA = zhangBuildSpanningTree(edges, "R0", preferredA);
    ZhangGraphBasis basisB = zhangBuildSpanningTree(edges, "R0", preferredB);
    BOOST_REQUIRE(basisA.connected);
    BOOST_REQUIRE(basisB.connected);
    BOOST_REQUIRE(basisA.treeEdges != basisB.treeEdges);

    ZhangExactMatrix forward = zhangCanonicalTransition(basisA, basisB);
    ZhangExactMatrix reverse = zhangCanonicalTransition(basisB, basisA);
    BOOST_REQUIRE(!forward.empty());
    BOOST_CHECK(zhangExactAbs(zhangExactDeterminant(forward)) == 1);
    BOOST_CHECK(
        zhangExactMultiply(reverse, forward) ==
        zhangExactIdentityMatrix(forward.size())
    );
}

BOOST_AUTO_TEST_CASE(e27_if_user_integer_and_covariance_algebra_closes)
{
    constexpr double lambda1 = 0.190293672798365;
    constexpr double lambda2 = 0.244210213424568;
    const auto coefficients = zhangIfUserCoefficients(lambda1, lambda2);
    BOOST_REQUIRE(coefficients.valid);
    BOOST_CHECK_SMALL(coefficients.alpha + coefficients.beta - 1, 1e-14);
    BOOST_CHECK_SMALL(
        coefficients.narrowLaneWavelength -
            lambda1 * lambda2 / (lambda1 + lambda2),
        1e-14);

    constexpr double firstInteger = 123456;
    constexpr double secondInteger = 123411;
    constexpr double wideLaneInteger = firstInteger - secondInteger;
    const double ifAmbiguity = zhangIfAmbiguityMetres(
        coefficients, lambda1, lambda2, firstInteger, secondInteger);
    const double recoveredFirst = zhangIfConditionedFirstInteger(
        coefficients, lambda2, ifAmbiguity, wideLaneInteger);
    BOOST_CHECK_SMALL(recoveredFirst - firstInteger, 1e-10);

    // Per-frequency correction precedes IF construction.  The combined
    // correction must be exactly alpha*c1+beta*c2.
    constexpr double rawFirst = 12.125;
    constexpr double rawSecond = 15.750;
    constexpr double correctionFirst = 0.237;
    constexpr double correctionSecond = -0.119;
    const double correctedIf =
        coefficients.alpha * (rawFirst + correctionFirst) +
        coefficients.beta * (rawSecond + correctionSecond);
    const double rawIf = coefficients.alpha * rawFirst +
        coefficients.beta * rawSecond;
    BOOST_CHECK_SMALL(
        (correctedIf - rawIf) -
            (coefficients.alpha * correctionFirst +
             coefficients.beta * correctionSecond),
        1e-12);

    // Three satellites with [clock, phase1, phase2] parameters.  The exact
    // functional includes cross-frequency and cross-satellite covariance.
    MatrixXd factor(9, 4);
    factor <<
        0.30,  0.02,  0.00,  0.00,
        0.10,  0.15,  0.01,  0.00,
        0.12, -0.03,  0.18,  0.00,
        0.28,  0.01,  0.00,  0.02,
        0.09,  0.14, -0.02,  0.01,
        0.11, -0.04,  0.17, -0.01,
        0.31,  0.03,  0.01, -0.02,
        0.08,  0.16,  0.00,  0.03,
        0.13, -0.02,  0.19,  0.01;
    const MatrixXd covariance = factor * factor.transpose();
    const MatrixXd transform = zhangIfProductSdFunctional(
        3, 0, coefficients, true);
    const MatrixXd propagated = zhangPropagateIfProductSdCovariance(
        covariance, 3, 0, coefficients, true);
    BOOST_REQUIRE_EQUAL(propagated.rows(), 2);
    BOOST_CHECK_SMALL(
        (propagated - transform * covariance * transform.transpose()).norm(),
        1e-14);
    BOOST_CHECK((propagated.diagonal().array() >= 0).all());

    // Removing cross-satellite covariance is not an equivalent stochastic
    // model and must be observable in this regression.
    MatrixXd blockDiagonal = covariance;
    for (int left = 0; left < 3; left++)
    for (int right = 0; right < 3; right++)
    {
        if (left != right)
        {
            blockDiagonal.block<3, 3>(3 * left, 3 * right).setZero();
        }
    }
    const MatrixXd independentApproximation =
        transform * blockDiagonal * transform.transpose();
    BOOST_CHECK_GT((propagated - independentApproximation).norm(), 1e-6);
}

BOOST_AUTO_TEST_CASE(e29_hybrid_user_model_follows_frozen_document_equations)
{
	constexpr double lambda1 = 0.190293672798365;
	constexpr double mu1 = 1;
	constexpr double mu2 = 1.646944444444444;
	// The existing server stores B^phi.  The supplied user equations use the
	// correction-side phase bias delta^G=-B^phi, so the adapter must map signs
	// explicitly instead of silently renaming the internal coordinate.
	const auto product =
		zhangDualFrequencyHybridProductsFromInternalPhaseStates(
			12.5, -0.08, 0.11);
	BOOST_REQUIRE(product.valid);
	// Equations (24)--(27): both baseline code signals consume exactly the
	// same satellite clock.  Only phase has signal-specific products.
	BOOST_CHECK_SMALL(
		product.codeCorrectionToAddMetres(0)
		- product.codeCorrectionToAddMetres(1), 1e-14);
	BOOST_CHECK_SMALL(
		product.phaseCorrectionToAddMetres(0) - 12.58, 1e-14);
	BOOST_CHECK_SMALL(
		product.phaseCorrectionToAddMetres(1) - 12.39, 1e-14);
	constexpr double rawCode = 23456789.25;
	constexpr double rawPhase = 23456780.75;
	BOOST_CHECK_SMALL(
		zhangHybridApplyLeftCorrection(
			rawCode, product.codeCorrectionToAddMetres(0))
		- (rawCode + 12.5), 1e-12);
	BOOST_CHECK_SMALL(
		zhangHybridApplyLeftCorrection(
			rawPhase, product.phaseCorrectionToAddMetres(0))
		- (rawPhase + 12.58), 1e-12);

	// Equations (3)--(13): physical receiver/satellite code biases are not
	// discarded; their two IF/GF directions are absorbed by the estimable
	// clocks, ionosphere and phase biases.
	constexpr double receiverCodeIf = 2.75;
	constexpr double receiverCodeGf = -0.625;
	constexpr double satelliteCodeIf = -1.2;
	constexpr double satelliteCodeGf = 0.35;
	const double receiverCode1 = receiverCodeIf + mu1 * receiverCodeGf;
	const double receiverCode2 = receiverCodeIf + mu2 * receiverCodeGf;
	const double satelliteCode1 = satelliteCodeIf + mu1 * satelliteCodeGf;
	const double satelliteCode2 = satelliteCodeIf + mu2 * satelliteCodeGf;
	const auto receiverDatum = zhangHybridCodeIfGfDatum(
		receiverCode1, receiverCode2, mu1, mu2);
	const auto satelliteDatum = zhangHybridCodeIfGfDatum(
		satelliteCode1, satelliteCode2, mu1, mu2);
	BOOST_REQUIRE(receiverDatum.valid);
	BOOST_REQUIRE(satelliteDatum.valid);
	BOOST_CHECK_SMALL(
		receiverDatum.ifBiasMetres - receiverCodeIf, 1e-13);
	BOOST_CHECK_SMALL(
		receiverDatum.gfBiasMetres - receiverCodeGf, 1e-13);
	BOOST_CHECK_SMALL(
		satelliteDatum.ifBiasMetres - satelliteCodeIf, 1e-13);
	BOOST_CHECK_SMALL(
		satelliteDatum.gfBiasMetres - satelliteCodeGf, 1e-13);

	constexpr double receiverClock = 8.5;
	constexpr double satelliteClock = -3.25;
	constexpr double ionosphere = 4.2;
	const double estimableReceiverClock = receiverClock + receiverCodeIf;
	const double estimableSatelliteClock = satelliteClock + satelliteCodeIf;
	const double estimableIonosphere = ionosphere
		+ receiverCodeGf - satelliteCodeGf;
	for (const auto& [mu, receiverCode, satelliteCode] :
		{std::tuple{mu1, receiverCode1, satelliteCode1},
		 std::tuple{mu2, receiverCode2, satelliteCode2}})
	{
		BOOST_CHECK_SMALL(
			zhangHybridOriginalCodePrediction(
				receiverClock, satelliteClock, ionosphere, mu,
				receiverCode, satelliteCode)
			- zhangHybridFullRankCodePrediction(
				estimableReceiverClock, estimableSatelliteClock,
				estimableIonosphere, mu),
			1e-13);
	}

	constexpr double receiverPhaseBias = 0.14;
	constexpr double satellitePhaseBias = -0.09;
	constexpr long long ambiguity = 123456;
	const double estimableReceiverPhaseBias = receiverPhaseBias
		- receiverCodeIf + mu1 * receiverCodeGf;
	const double estimableSatellitePhaseBias = satellitePhaseBias
		- satelliteCodeIf + mu1 * satelliteCodeGf;
	BOOST_CHECK_SMALL(
		zhangHybridOriginalPhasePrediction(
			receiverClock, satelliteClock, ionosphere, mu1,
			receiverPhaseBias, satellitePhaseBias,
			lambda1, ambiguity)
		- zhangHybridFullRankPhasePrediction(
			estimableReceiverClock, estimableSatelliteClock,
			estimableIonosphere, mu1,
			estimableReceiverPhaseBias, estimableSatellitePhaseBias,
			lambda1, ambiguity),
		1e-12);

	const MatrixXd sd = zhangHybridSatelliteSingleDifferenceTransform(4, 1);
	BOOST_REQUIRE_EQUAL(sd.rows(), 3);
	BOOST_REQUIRE_EQUAL(sd.cols(), 4);
	BOOST_CHECK_SMALL((sd * Vector4d::Ones()).norm(), 1e-14);
	const Vector4d ambiguities(17, -4, 8, 21);
	const Vector3d expectedSd(21, 12, 25);
	BOOST_CHECK_SMALL((sd * ambiguities - expectedSd).norm(), 1e-14);

	const Matrix2d integerTransform =
		zhangHybridWideLaneFirstIntegerTransform();
	BOOST_CHECK_EQUAL(std::llround(integerTransform.determinant()), 1);
	const Vector2d integerPair(31, 24);
	const Vector2d wlFirst = integerTransform * integerPair;
	BOOST_CHECK_SMALL(wlFirst(0) - 7, 1e-14);
	BOOST_CHECK_SMALL(wlFirst(1) - 31, 1e-14);
	BOOST_CHECK_SMALL(
		(integerTransform.inverse() * wlFirst - integerPair).norm(), 1e-14);

	BOOST_CHECK(
		zhangHybridIntegerUsability(true, false)
		== ZhangHybridIntegerUsability::FLOAT_ONLY);
	BOOST_CHECK(
		zhangHybridIntegerUsability(true, true)
		== ZhangHybridIntegerUsability::PPP_AR_USABLE);
	BOOST_CHECK(
		zhangHybridIntegerUsability(false, true)
		== ZhangHybridIntegerUsability::UNUSABLE);
	BOOST_CHECK(zhangHybridRelativeIntegerPairCertified(
		true, "GPS-L1C-COMP-A", true, "GPS-L1C-COMP-A"));
	BOOST_CHECK(!zhangHybridRelativeIntegerPairCertified(
		true, "GPS-L1C-COMP-A", true, "GPS-L1C-COMP-B"));
	BOOST_CHECK(!zhangHybridRelativeIntegerPairCertified(
		true, "NONE", true, "NONE"));
	BOOST_CHECK(!zhangHybridRelativeIntegerPairCertified(
		true, "GPS-L1C-COMP-A", false, "GPS-L1C-COMP-A"));

	Matrix4d userNoise = Matrix4d::Identity() * 0.04;
	Matrix4d factor;
	factor <<
		0.30,  0.02, 0.00, 0.00,
		0.25, -0.01, 0.04, 0.00,
		0.28,  0.03, 0.00, 0.02,
		0.22, -0.02, 0.05, 0.01;
	const Matrix4d productCovariance = factor * factor.transpose();
	const MatrixXd corrected = zhangHybridCorrectedObservationCovariance(
		userNoise, productCovariance);
	BOOST_CHECK_SMALL(
		(corrected - userNoise - productCovariance).norm(), 1e-14);
	const MatrixXd propagated = zhangHybridSingleDifferenceCovariance(
		userNoise, productCovariance, sd);
	BOOST_CHECK_SMALL(
		(propagated - sd * corrected * sd.transpose()).norm(), 1e-14);

	Matrix4d diagonalProduct = productCovariance.diagonal().asDiagonal();
	const MatrixXd diagonalApproximation = zhangHybridSingleDifferenceCovariance(
		userNoise, diagonalProduct, sd);
	BOOST_CHECK_GT((propagated - diagonalApproximation).norm(), 1e-6);
}

BOOST_AUTO_TEST_CASE(e29_hybrid_real_gauge_is_overlap_gls_and_not_epoch_zero_mean)
{
	const std::vector<SatSys> firstSatellites = {
		SatSys(E_Sys::GPS, 1), SatSys(E_Sys::GPS, 2),
		SatSys(E_Sys::GPS, 3)};
	const std::vector<std::string> firstSegments = {"A", "B", "C"};
	Vector3d firstRaw(4.0, 7.0, 13.0);
	Matrix3d firstCovariance = Matrix3d::Zero();
	firstCovariance.diagonal() << 1.0, 4.0, 9.0;
	ZhangHybridRealGaugeTransport gauge;
	const auto first = gauge.transport(
		firstSatellites, firstSegments, firstRaw, firstCovariance);
	BOOST_REQUIRE_MESSAGE(first.valid, first.failureReason);
	BOOST_CHECK(first.newGeneration);
	BOOST_CHECK_EQUAL(first.generation, 0);
	const Vector3d inverseVarianceWeights(1.0, 0.25, 1.0 / 9.0);
	BOOST_CHECK_SMALL(
		inverseVarianceWeights.dot(first.values), 1e-12);
	BOOST_CHECK_SMALL(
		(first.covariance
		 - first.transform * firstCovariance * first.transform.transpose()).norm(),
		1e-13);

	// G04 joins, but the old overlap remains.  A single common offset aligns
	// G01--G03 to their previous gauge; the four-satellite result is therefore
	// not re-zeroed over the changed membership.
	const std::vector<SatSys> secondSatellites = {
		SatSys(E_Sys::GPS, 1), SatSys(E_Sys::GPS, 2),
		SatSys(E_Sys::GPS, 3), SatSys(E_Sys::GPS, 4)};
	const std::vector<std::string> secondSegments = {"A", "B", "C", "D"};
	Vector4d secondRaw;
	secondRaw.head<3>() = firstRaw.array() + 2.5;
	secondRaw(3) = 100.0;
	Matrix4d secondCovariance = Matrix4d::Identity();
	const auto second = gauge.transport(
		secondSatellites, secondSegments, secondRaw, secondCovariance);
	BOOST_REQUIRE_MESSAGE(second.valid, second.failureReason);
	BOOST_CHECK(!second.newGeneration);
	BOOST_CHECK_EQUAL(second.overlapCount, 3);
	BOOST_CHECK_SMALL((second.values.head<3>() - first.values).norm(), 1e-12);
	BOOST_CHECK_GT(std::abs(second.values.sum()), 1.0);

	// No unchanged physical segment remains: continuity cannot be fabricated.
	const std::vector<std::string> thirdSegments = {"A2", "B2", "C2", "D2"};
	const auto third = gauge.transport(
		secondSatellites, thirdSegments, secondRaw, secondCovariance);
	BOOST_REQUIRE_MESSAGE(third.valid, third.failureReason);
	BOOST_CHECK(third.newGeneration);
	BOOST_CHECK_EQUAL(third.generation, 1);
	BOOST_CHECK_EQUAL(third.overlapCount, 0);
	BOOST_CHECK_SMALL(third.values.sum(), 1e-12);
}

BOOST_AUTO_TEST_CASE(e29_persistent_dynamic_gate_uses_manager_not_product_tree_proof)
{
	const auto persistent = zhangHybridInitialIntegerGate(
		true,
		true,   // backend graph valid
		true,   // independently audited product functional valid
		false,  // PRODUCT_TREE runtime alignment deliberately unavailable
		false,  // named row deliberately absent from held lattice
		true,   // persistent kappa datum continuous
		true);  // persistent component precision valid
	BOOST_CHECK(persistent.structureValid);
	BOOST_CHECK(persistent.datumContinuous);
	BOOST_CHECK(persistent.precisionValid);

	const auto productTree = zhangHybridInitialIntegerGate(
		false, true, true, false, false, true, true);
	BOOST_CHECK(productTree.structureValid);
	BOOST_CHECK(!productTree.datumContinuous);
	BOOST_CHECK(!productTree.precisionValid);

	const auto invalidStructure = zhangHybridInitialIntegerGate(
		true, true, false, true, true, true, true);
	BOOST_CHECK(!invalidStructure.structureValid);
}

BOOST_AUTO_TEST_CASE(e29_hybrid_real_gauge_checkpoint_preserves_transport)
{
	const std::vector<SatSys> satellites = {
		SatSys(E_Sys::GPS, 1), SatSys(E_Sys::GPS, 3)};
	const std::vector<std::string> segments = {"G01-S0", "G03-S0"};
	Vector2d raw(2.0, 5.0);
	Matrix2d covariance;
	covariance << 1.0, 0.2, 0.2, 2.0;
	ZhangHybridRealGaugeTransport original;
	BOOST_REQUIRE(original.transport(
		satellites, segments, raw, covariance).valid);
	const auto checkpoint = original.checkpointState();
	ZhangHybridRealGaugeTransport restored;
	std::string failureReason;
	BOOST_REQUIRE_MESSAGE(
		restored.restoreCheckpointState(checkpoint, &failureReason),
		failureReason);
	Vector2d shifted = raw.array() + 7.0;
	const auto expected = original.transport(
		satellites, segments, shifted, covariance);
	const auto actual = restored.transport(
		satellites, segments, shifted, covariance);
	BOOST_REQUIRE(expected.valid);
	BOOST_REQUIRE(actual.valid);
	BOOST_CHECK_SMALL((actual.values - expected.values).norm(), 1e-13);
	BOOST_CHECK_SMALL((actual.covariance - expected.covariance).norm(), 1e-13);
	BOOST_CHECK_EQUAL(actual.generation, expected.generation);
}

BOOST_AUTO_TEST_CASE(e29_hybrid_real_gauge_candidate_copy_is_transactional)
{
	const std::vector<SatSys> satellites = {
		SatSys(E_Sys::GPS, 1), SatSys(E_Sys::GPS, 3)};
	const std::vector<std::string> segments = {"G01-S0", "G03-S0"};
	Vector2d raw(2.0, 5.0);
	const Matrix2d covariance = Matrix2d::Identity();
	ZhangHybridRealGaugeTransport committed;
	const auto initial = committed.transport(
		satellites, segments, raw, covariance);
	BOOST_REQUIRE(initial.valid);
	const auto before = committed.checkpointState();

	// appendProductCovariance evaluates every block on a copy and assigns the
	// copies only after every block succeeds.  Aborting the epoch must therefore
	// leave the committed history byte-for-byte equivalent to its checkpoint.
	ZhangHybridRealGaugeTransport candidate = committed;
	// Include a relative change as well as a common change.  A pure common
	// shift is deliberately removed by the gauge and may produce a committed
	// checkpoint numerically identical to the old one; that is not evidence
	// of a failed transaction.
	Vector2d shifted;
	shifted << raw(0) + 11.0, raw(1) + 12.0;
	const auto evaluated = candidate.transport(
		satellites, segments, shifted, covariance);
	BOOST_REQUIRE(evaluated.valid);
	const auto afterAbort = committed.checkpointState();
	BOOST_CHECK_EQUAL(afterAbort.initialized, before.initialized);
	BOOST_CHECK_EQUAL(afterAbort.generation, before.generation);
	BOOST_CHECK(afterAbort.previousValues == before.previousValues);
	BOOST_CHECK(afterAbort.previousSegments == before.previousSegments);

	const auto candidateCheckpoint = candidate.checkpointState();
	committed = std::move(candidate);
	const auto afterCommit = committed.checkpointState();
	BOOST_CHECK(afterCommit.previousValues != before.previousValues);
	BOOST_CHECK(afterCommit.previousValues == candidateCheckpoint.previousValues);
	BOOST_CHECK(afterCommit.previousSegments == candidateCheckpoint.previousSegments);
	BOOST_CHECK_SMALL(
		afterCommit.previousValues.at(satellites.front())
			- evaluated.values(0),
		1e-13);
}

BOOST_AUTO_TEST_CASE(hybrid_stable_frontend_controller_is_prepare_commit_atomic)
{
	using State = std::map<std::string, ZhangHybridRealGaugeTransport>;
	const std::vector<SatSys> satellites = {
		SatSys(E_Sys::GPS, 1), SatSys(E_Sys::GPS, 3)};
	const std::vector<std::string> segments = {"G01-S0", "G03-S0"};
	const Matrix2d covariance = Matrix2d::Identity();
	Vector2d raw(2.0, 5.0);
	State persistent;
	BOOST_REQUIRE(persistent["L1"].transport(
		satellites, segments, raw, covariance).valid);
	const auto before = persistent.at("L1").checkpointState();

	ZhangHybridStableFrontend controller;
	auto rejected = controller.prepare(persistent);
	Vector2d changed(14.0, 19.0);
	BOOST_REQUIRE(rejected.preparedState["L1"].transport(
		satellites, segments, changed, covariance).valid);
	controller.validateIntegerAlignment(rejected, true);
	controller.validateRealGauge(rejected, true);
	controller.validateComponentConsistency(rejected, false);
	controller.validateMetadata(rejected, true);
	BOOST_CHECK(!controller.commit(persistent, rejected));
	controller.rollback(rejected);
	const auto afterRollback = persistent.at("L1").checkpointState();
	BOOST_CHECK(afterRollback.previousValues == before.previousValues);
	BOOST_CHECK(afterRollback.previousSegments == before.previousSegments);

	auto accepted = controller.prepare(persistent);
	BOOST_REQUIRE(accepted.preparedState["L1"].transport(
		satellites, segments, changed, covariance).valid);
	controller.validateIntegerAlignment(accepted, true);
	controller.validateRealGauge(accepted, true);
	controller.validateComponentConsistency(accepted, true);
	controller.validateMetadata(accepted, true);
	const auto expected = accepted.preparedState.at("L1").checkpointState();
	BOOST_REQUIRE(controller.commit(persistent, accepted));
	const auto committed = persistent.at("L1").checkpointState();
	BOOST_CHECK(committed.previousValues == expected.previousValues);
	BOOST_CHECK(committed.previousSegments == expected.previousSegments);
}

BOOST_AUTO_TEST_CASE(e29_hybrid_real_gauge_maps_transform_cross_block_covariance)
{
	const std::vector<SatSys> satellites = {
		SatSys(E_Sys::GPS, 1), SatSys(E_Sys::GPS, 3)};
	const std::vector<std::string> clockSegments = {"CLOCK-G01", "CLOCK-G03"};
	const std::vector<std::string> phaseSegments = {"G01-L1-S0", "G03-L1-S0"};
	Matrix4d factor;
	factor <<
		0.8, 0.0, 0.0, 0.0,
		0.2, 0.7, 0.0, 0.0,
		0.3, 0.1, 0.6, 0.0,
		0.1, 0.2, 0.2, 0.5;
	const Matrix4d rawCovariance = factor * factor.transpose();
	Vector2d clockMean(4.0, 7.0);
	Vector2d phaseMean(-2.0, 3.0);
	ZhangHybridRealGaugeTransport clockGauge;
	ZhangHybridRealGaugeTransport phaseGauge;
	const auto clock = clockGauge.transport(
		satellites, clockSegments, clockMean,
		rawCovariance.block<2, 2>(0, 0));
	const auto phase = phaseGauge.transport(
		satellites, phaseSegments, phaseMean,
		rawCovariance.block<2, 2>(2, 2));
	BOOST_REQUIRE(clock.valid);
	BOOST_REQUIRE(phase.valid);

	Matrix4d fullTransform = Matrix4d::Zero();
	fullTransform.block<2, 2>(0, 0) = clock.transform;
	fullTransform.block<2, 2>(2, 2) = phase.transform;
	const Matrix4d propagated =
		fullTransform * rawCovariance * fullTransform.transpose();
	const Matrix2d expectedCross = clock.transform
		* rawCovariance.block<2, 2>(0, 2)
		* phase.transform.transpose();
	BOOST_CHECK_SMALL(
		(propagated.block<2, 2>(0, 2) - expectedCross).norm(), 1e-14);
	BOOST_CHECK_GT(
		(propagated.block<2, 2>(0, 2)
			- rawCovariance.block<2, 2>(0, 2)).norm(),
		1e-6);
}

BOOST_AUTO_TEST_CASE(e29_dual_frequency_partition_is_component_intersection)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	SatSys g03(E_Sys::GPS, 3);
	SatSys g04(E_Sys::GPS, 4);
	const std::map<SatSys, std::string> l1 = {
		{g01, "L1-A"}, {g02, "L1-A"}, {g03, "L1-A"}, {g04, "L1-B"}};
	const std::map<SatSys, std::string> l2 = {
		{g01, "L2-X"}, {g02, "L2-X"}, {g03, "L2-Y"}, {g04, "L2-X"}};
	const auto components = zhangHybridDualFrequencyComponents(l1, l2);
	BOOST_REQUIRE_EQUAL(components.size(), 1);
	const auto& members = components.begin()->second;
	BOOST_CHECK_EQUAL(members.size(), 2);
	BOOST_CHECK(members.count(g01));
	BOOST_CHECK(members.count(g02));
	BOOST_CHECK(!members.count(g03));
	BOOST_CHECK(!members.count(g04));
}

BOOST_AUTO_TEST_CASE(e27_if_wl_conditioning_uses_the_full_cross_covariance)
{
    const Vector2d ifMean(17.35, -4.20);
    const Vector2d wideLaneMean(3.08, -1.12);
    const Vector2d fixedWideLane(3, -1);
    Matrix2d ifCovariance;
    ifCovariance << 0.40, 0.08,
                    0.08, 0.30;
    Matrix2d wideLaneCovariance;
    wideLaneCovariance << 0.20, 0.03,
                          0.03, 0.16;
    Matrix2d crossCovariance;
    crossCovariance << 0.050, -0.010,
                       0.015,  0.040;
    constexpr double coefficient = -0.73;

    const auto conditioned = zhangConditionFirstIntegerGivenWideLane(
        ifMean, ifCovariance, wideLaneMean, wideLaneCovariance,
        crossCovariance, fixedWideLane, coefficient);
    BOOST_REQUIRE_MESSAGE(conditioned.valid, conditioned.failureReason);
    const Matrix2d inverse = wideLaneCovariance.inverse();
    const Vector2d expectedMean = ifMean + coefficient * fixedWideLane +
        crossCovariance * inverse * (fixedWideLane - wideLaneMean);
    const Matrix2d expectedCovariance = ifCovariance -
        crossCovariance * inverse * crossCovariance.transpose();
    BOOST_CHECK_SMALL((conditioned.mean - expectedMean).norm(), 1e-13);
    BOOST_CHECK_SMALL(
        (conditioned.covariance - expectedCovariance).norm(), 1e-13);
    BOOST_CHECK_LT(conditioned.covariance.trace(), ifCovariance.trace());

    const auto independentApproximation =
        zhangConditionFirstIntegerGivenWideLane(
            ifMean, ifCovariance, wideLaneMean, wideLaneCovariance,
            Matrix2d::Zero(), fixedWideLane, coefficient);
    BOOST_REQUIRE(independentApproximation.valid);
    BOOST_CHECK_SMALL(
        (independentApproximation.covariance - ifCovariance).norm(), 1e-14);
    BOOST_CHECK_GT(
        (independentApproximation.mean - conditioned.mean).norm(), 1e-4);
}

BOOST_AUTO_TEST_CASE(e27_wide_lane_raw_noise_sensitivity_reconstructs_covariance)
{
    ZhangIfWideLaneAccumulator accumulator(3600, 60, 32);
    const std::vector<int> satellites = {3, 7, 11};
    for (int satellite : satellites)
    {
        accumulator.setArcVersion(satellite, 1);
    }
    const Vector3d physical(10.2, 15.2, 7.2);
    const Vector3d variances(0.04, 0.09, 0.16);
    const Matrix3d rawDesign = Matrix3d::Identity();
    const Matrix3d covariance = variances.asDiagonal();
    for (int epoch = 0; epoch < 4; epoch++)
    {
        const std::vector<std::string> keys = {
            "E" + std::to_string(epoch) + "-G03",
            "E" + std::to_string(epoch) + "-G07",
            "E" + std::to_string(epoch) + "-G11"};
        accumulator.addEpoch(
            epoch * 60, satellites, physical, covariance,
            keys, variances, rawDesign);
    }
    const auto estimate = accumulator.estimate(satellites, 3, 240);
    BOOST_REQUIRE_MESSAGE(estimate.valid, estimate.failureReason);
    BOOST_REQUIRE_EQUAL(estimate.noiseSensitivity.size(), 12);
    Matrix2d reconstructed = Matrix2d::Zero();
    for (const auto& [key, sensitivity] : estimate.noiseSensitivity)
    {
        reconstructed += estimate.noiseVariance.at(key) *
            sensitivity * sensitivity.transpose();
    }
    BOOST_CHECK_SMALL(
        (reconstructed - estimate.covariance).norm(), 1e-12);

    const auto exchanged = accumulator.estimate(satellites, 11, 240);
    BOOST_REQUIRE(exchanged.valid);
    Matrix2d transform;
    // ref=3 targets [7,11]; ref=11 targets [3,7].
    transform << 0, -1,
                 1, -1;
    for (const auto& [key, sensitivity] : estimate.noiseSensitivity)
    {
        BOOST_REQUIRE(exchanged.noiseSensitivity.count(key));
        BOOST_CHECK_SMALL(
            (exchanged.noiseSensitivity.at(key) -
             transform * sensitivity).norm(), 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(e27_wide_lane_window_is_reference_invariant_and_arc_safe)
{
    ZhangIfWideLaneAccumulator accumulator(3600, 60, 32);
    for (int satellite : {3, 7, 11, 19})
    {
        accumulator.setArcVersion(satellite, 1);
    }
    const std::vector<int> satellites = {3, 7, 11, 19};
    const VectorXd physical = (VectorXd(4) << 10.2, 15.2, 7.2, 21.2).finished();
    MatrixXd covariance = MatrixXd::Constant(4, 4, 0.01);
    covariance.diagonal().array() += 0.04;
    for (int epoch = 0; epoch < 10; epoch++)
    {
        accumulator.addEpoch(epoch * 60, satellites, physical, covariance);
    }
    const auto reference3 = accumulator.estimate(satellites, 3, 600);
    const auto reference11 = accumulator.estimate(satellites, 11, 600);
    BOOST_REQUIRE(reference3.valid && reference11.valid);

    // Transform SDs relative to G03 into SDs relative to G11 exactly.
    MatrixXd exchange = MatrixXd::Zero(3, 3);
    // ref=3 targets [7,11,19]; ref=11 targets [3,7,19].
    exchange.row(0) << 0, -1, 0;
    exchange.row(1) << 1, -1, 0;
    exchange.row(2) << 0, -1, 1;
    BOOST_CHECK_SMALL(
        (reference11.mean - exchange * reference3.mean).norm(), 1e-12);
    BOOST_CHECK_SMALL(
        (reference11.covariance -
         exchange * reference3.covariance * exchange.transpose()).norm(),
        1e-12);

    // A real physical arc change invalidates old factors involving G11.
    accumulator.setArcVersion(11, 2);
    const auto changedArc = accumulator.estimate(satellites, 3, 600);
    BOOST_CHECK(!changedArc.valid);
}

BOOST_AUTO_TEST_CASE(satellite_product_target_is_exact_across_tree_exchange)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);

    // K2,2 has one current fundamental cycle.  The persistent product tree
    // uses that current chord, so the G02-G01 product correction is +k.
    std::set<ZhangGraphEdge> k22Edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    ZhangGraphBasis current = zhangBuildSpanningTree(k22Edges, "R0");
    ZhangGraphBasis product = zhangBuildSpanningTree(
        k22Edges,
        "R0",
        {{"R0", g01}, {"R1", g01}, {"R1", g02}}
    );
    ZhangSatelliteProductTarget k22 =
        zhangBuildSatelliteProductTarget(current, product, g01);
    BOOST_REQUIRE(k22.valid);
    BOOST_REQUIRE_EQUAL(k22.matrix.size(), 1);
    BOOST_CHECK(k22.matrix.front() == ZhangExactVector({1}));
    BOOST_CHECK(k22.targetSatellites == std::vector<SatSys>({g02}));

    // Three stations/three satellites: the correction G*k changes when the
    // dynamic tree changes, but z_T + G*k must equal the same persistent
    // product datum exactly.  Comparing G*k alone would be mathematically
    // wrong because the dynamic-tree node integer potential changes too.
    SatSys g03(E_Sys::GPS, 3);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R0", g03},
        {"R1", g01}, {"R1", g02}, {"R1", g03},
        {"R2", g01}, {"R2", g02}, {"R2", g03},
    };
    ZhangGraphBasis basisA = zhangBuildSpanningTree(
        edges, "R0",
        {{"R0", g01}, {"R0", g02}, {"R0", g03},
         {"R1", g01}, {"R2", g02}}
    );
    ZhangGraphBasis basisB = zhangBuildSpanningTree(
        edges, "R0",
        {{"R0", g01}, {"R1", g01}, {"R1", g03},
         {"R2", g02}, {"R2", g03}}
    );
    ZhangGraphBasis productBasis = zhangBuildSpanningTree(
        edges, "R0",
        {{"R0", g02}, {"R1", g02}, {"R1", g03},
         {"R2", g01}, {"R2", g03}}
    );
    BOOST_REQUIRE(basisA.connected && basisB.connected && productBasis.connected);

    map<ZhangGraphEdge, ZhangExactInteger> physical;
    int value = 1;
    for (const auto& edge : edges)
    {
        physical[edge] = value++;
    }
    auto cycleValues = [&](const ZhangGraphBasis& basis)
    {
        ZhangCanonicalIntegerAudit audit = zhangCanonicalIntegerAudit(basis);
        ZhangExactVector cycles;
        for (const auto& chord : audit.chordEdges)
        {
            ZhangExactInteger cycle = 0;
            for (const auto& [edge, coefficient] :
                 zhangFundamentalCycle(basis, chord))
            {
                cycle += coefficient * physical.at(edge);
            }
            cycles.push_back(cycle);
        }
        return cycles;
    };
    auto satelliteDatum = [&](const ZhangGraphBasis& basis)
    {
        ZhangCanonicalIntegerAudit audit = zhangCanonicalIntegerAudit(basis);
        ZhangExactVector treeValues;
        for (const auto& edge : audit.treeEdges)
        {
            treeValues.push_back(physical.at(edge));
        }
        ZhangExactVector nodes =
            zhangExactMatrixTimesColumn(audit.treeInverse, treeValues);
        const std::size_t satelliteOffset = basis.receivers.size() - 1;
        map<SatSys, ZhangExactInteger> satelliteValues;
        std::size_t row = satelliteOffset;
        for (const auto& satellite : basis.satellites)
        {
            satelliteValues[satellite] = nodes[row++];
        }
        ZhangExactVector differences;
        for (const auto& satellite : basis.satellites)
        {
            if (satellite != g01)
            {
                differences.push_back(
                    satelliteValues[satellite] - satelliteValues[g01]
                );
            }
        }
        return differences;
    };
    ZhangExactVector productDatum = satelliteDatum(productBasis);
    for (const auto& basis : {basisA, basisB})
    {
        ZhangSatelliteProductTarget target =
            zhangBuildSatelliteProductTarget(basis, productBasis, g01);
        BOOST_REQUIRE(target.valid);
		ZhangIntegerLatticeMembership smith =
			zhangIntegerRowLatticeContains(
				target.matrix,
				ZhangExactVector(target.currentChords.size()));
		BOOST_CHECK_EQUAL(smith.rank, target.matrix.size());
		for (const auto& invariant : smith.smithInvariants)
		{
			BOOST_CHECK(zhangExactAbs(invariant) == 1);
		}
        ZhangExactVector corrected = satelliteDatum(basis);
        ZhangExactVector correction =
            zhangExactMatrixTimesColumn(target.matrix, cycleValues(basis));
        for (std::size_t row = 0; row < corrected.size(); row++)
        {
            corrected[row] += correction[row];
        }
        BOOST_CHECK(corrected == productDatum);
    }

    // A product tree needs all target satellites but not every estimation
    // receiver.  The same exact projection must close when R2 remains in the
    // current state graph but is deliberately absent from the product core.
    std::set<ZhangGraphEdge> coreEdges;
    for (const auto& edge : edges)
    {
        if (edge.receiver != "R2")
        {
            coreEdges.insert(edge);
        }
    }
    ZhangGraphBasis coreProduct = zhangBuildSpanningTree(
        coreEdges, "R0",
        {{"R0", g02}, {"R1", g01}, {"R1", g03}});
    BOOST_REQUIRE(coreProduct.connected);
    BOOST_CHECK_LT(coreProduct.receivers.size(), basisA.receivers.size());
    const ZhangExactVector coreProductDatum = satelliteDatum(coreProduct);
    for (const auto& basis : {basisA, basisB})
    {
        ZhangSatelliteProductTarget target =
            zhangBuildSatelliteProductTarget(basis, coreProduct, g01);
        BOOST_REQUIRE_MESSAGE(target.valid, target.failureReason);
        ZhangExactVector corrected = satelliteDatum(basis);
        ZhangExactVector correction = zhangExactMatrixTimesColumn(
            target.matrix, cycleValues(basis));
        for (std::size_t row = 0; row < corrected.size(); row++)
        {
            corrected[row] += correction[row];
        }
        BOOST_CHECK(corrected == coreProductDatum);
    }
}

BOOST_AUTO_TEST_CASE(product_relation_basis_expands_independently_to_physical_arcs)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    const std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    const ZhangGraphBasis current = zhangBuildSpanningTree(edges, "R0");
    const ZhangGraphBasis product = zhangBuildSpanningTree(
        edges,
        "R0",
        {{"R0", g01}, {"R1", g01}, {"R1", g02}});
    const ZhangProductRelationBasis relationBasis =
        ProductRelationBasisBuilder::build(current, product, g01);

    BOOST_REQUIRE_MESSAGE(relationBasis.valid, relationBasis.failureReason);
    BOOST_CHECK_EQUAL(relationBasis.namedRelationCount, 1);
    BOOST_CHECK_EQUAL(relationBasis.exactRank, 1);
    BOOST_CHECK(relationBasis.independentNamedIndices == std::vector<int>({0}));
    BOOST_CHECK(relationBasis.primitive);
    BOOST_CHECK(relationBasis.saturationIndex == 1);
    BOOST_CHECK(relationBasis.admissibleCompletionProven);
    BOOST_CHECK(relationBasis.networkLatticeContained);
    BOOST_CHECK(relationBasis.networkClosureExactZero);
    BOOST_CHECK(
        zhangExactMultiply(
            relationBasis.networkContainmentTransform,
            relationBasis.networkIntegerBasis) ==
        relationBasis.exactRowBasis);
    BOOST_CHECK(relationBasis.nuisanceOrthogonal);
    BOOST_CHECK(relationBasis.physicalExpansionValid);
    BOOST_REQUIRE_EQUAL(relationBasis.namedRelations.size(), 1);
    const auto& relation = relationBasis.namedRelations.front();
    BOOST_CHECK(relation.satellite == g02);
    BOOST_CHECK(relation.referenceSatellite == g01);
    BOOST_CHECK(relation.currentCycleCoefficients == ZhangExactVector({1}));
    BOOST_CHECK(std::all_of(
        relation.nuisanceCoefficients.begin(),
        relation.nuisanceCoefficients.end(),
        [](const auto& coefficient) { return coefficient == 0; }));

    BOOST_REQUIRE_EQUAL(relationBasis.currentChords.size(), 1);
    const auto expectedCycle = zhangFundamentalCycle(
        current, relationBasis.currentChords.front());
    std::map<ZhangGraphEdge, ZhangExactInteger> expectedPhysical;
    for (const auto& [edge, coefficient] : expectedCycle)
    {
        expectedPhysical[edge] += coefficient;
    }
    BOOST_CHECK(relation.physicalArcCoefficients == expectedPhysical);
}

BOOST_AUTO_TEST_CASE(
	product_relation_independence_is_reselected_after_column_availability)
{
	const SatSys g01(E_Sys::GPS, 1);
	const SatSys g02(E_Sys::GPS, 2);
	const SatSys g03(E_Sys::GPS, 3);
	ZhangProductRelationBasis basis;
	basis.currentChords = {
		{"R0", g01}, {"R0", g02}, {"R0", g03}};
	for (const ZhangExactVector& row : ZhangExactMatrix{
			{1, 0, 1}, {0, 1, 1}, {1, -1, 0}})
	{
		ZhangProductRelationRow relation;
		relation.currentCycleCoefficients = row;
		basis.namedRelations.push_back(std::move(relation));
	}
	// Rows 0 and 1 are a valid global pivot subset, but both use the missing
	// third chord.  Their dependent difference (row 2) remains exactly
	// representable and must replace them after availability filtering.
	basis.independentNamedIndices = {0, 1};
	const std::set<ZhangGraphEdge> partial = {
		basis.currentChords[0], basis.currentChords[1]};
	BOOST_CHECK_EQUAL(
		zhangLegacyFilteredMappableProductRelationRank(basis, partial), 0);
	BOOST_CHECK(
		zhangIndependentMappableProductRelationIndices(basis, partial) ==
		std::vector<int>({2}));
	BOOST_CHECK_EQUAL(
		zhangLegacyFilteredMappableProductRelationRank(
			basis, std::set<ZhangGraphEdge>(
				basis.currentChords.begin(), basis.currentChords.end())), 2);
	BOOST_CHECK(
		zhangIndependentMappableProductRelationIndices(
			basis, std::set<ZhangGraphEdge>(
				basis.currentChords.begin(), basis.currentChords.end())) ==
		std::vector<int>({0, 1}));
	BOOST_CHECK(
		zhangIndependentMappableProductRelationIndices(basis, {}).empty());
}

BOOST_AUTO_TEST_CASE(product_relation_basis_is_reference_invariant)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g02(E_Sys::GPS, 2);
    const SatSys g03(E_Sys::GPS, 3);
    const std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R0", g03},
        {"R1", g01}, {"R1", g02}, {"R1", g03}
    };
    const ZhangGraphBasis current = zhangBuildSpanningTree(edges, "R0");
    const ZhangGraphBasis product = zhangBuildSpanningTree(
        edges,
        "R0",
        {{"R0", g01}, {"R1", g01}, {"R1", g02}, {"R1", g03}});

    const auto relativeToG01 = ProductRelationBasisBuilder::build(
        current, product, g01, E_Sys::GPS, E_ObsCode::L1C);
    const auto relativeToG02 = ProductRelationBasisBuilder::build(
        current, product, g02, E_Sys::GPS, E_ObsCode::L1C);
    BOOST_REQUIRE_MESSAGE(relativeToG01.valid, relativeToG01.failureReason);
    BOOST_REQUIRE_MESSAGE(relativeToG02.valid, relativeToG02.failureReason);
    BOOST_CHECK_EQUAL(relativeToG01.fullTargetRank, 2);
    BOOST_CHECK_EQUAL(relativeToG02.fullTargetRank, 2);
    BOOST_CHECK(relativeToG01.referenceSatellite == g01);
    BOOST_CHECK(relativeToG02.referenceSatellite == g02);

    // The named coordinate matrices differ with the reference, but their
    // canonical physical row lattice must be exactly identical.
    BOOST_CHECK(relativeToG01.exactRowBasis == relativeToG02.exactRowBasis);
    BOOST_CHECK_EQUAL(relativeToG01.exactHnf, relativeToG02.exactHnf);
    BOOST_CHECK(relativeToG01.networkClosureExactZero);
    BOOST_CHECK(relativeToG02.networkClosureExactZero);
    BOOST_CHECK(
        zhangExactMultiply(
            relativeToG01.networkContainmentTransform,
            relativeToG01.networkIntegerBasis) ==
        relativeToG01.exactRowBasis);
    BOOST_CHECK(
        zhangExactMultiply(
            relativeToG02.networkContainmentTransform,
            relativeToG02.networkIntegerBasis) ==
        relativeToG02.exactRowBasis);
}

BOOST_AUTO_TEST_CASE(product_relation_basis_fails_closed_on_graph_mismatch)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    const std::set<ZhangGraphEdge> currentEdges = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    const std::set<ZhangGraphEdge> productEdges = {
        {"R0", g01}, {"R0", g03}, {"R1", g01}, {"R1", g03}
    };
    const auto relationBasis = ProductRelationBasisBuilder::build(
        zhangBuildSpanningTree(currentEdges, "R0"),
        zhangBuildSpanningTree(productEdges, "R0"),
        g01);
    BOOST_CHECK(!relationBasis.valid);
    BOOST_CHECK_EQUAL(
        relationBasis.failureReason,
        "product_tree_target_node_or_root_mismatch");
}

BOOST_AUTO_TEST_CASE(product_transition_transport_includes_exact_affine_offset)
{
	ZhangProductIntegerFunctional previous;
	ZhangProductIntegerFunctional current;
	previous.satellite = SatSys(E_Sys::GPS, 7);
	current.satellite = previous.satellite;
	previous.valid = true;
	current.valid = true;
	previous.affineOffsetCycles = -4;
	current.affineOffsetCycles = 9;
	const auto transition = zhangProductIntegerFunctionalDifference(
		previous, current);
	BOOST_REQUIRE(transition.valid);
	BOOST_CHECK(transition.coefficients.empty());
	BOOST_CHECK(transition.affineOffsetCycles == 13);
	BOOST_CHECK(zhangCompleteProductTransitionInteger(transition, 21) == 34);
}

BOOST_AUTO_TEST_CASE(product_support_metrics_distinguish_paths_bridges_and_capacity)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);

    std::set<ZhangGraphEdge> k22 = {
        {"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}
    };
    BOOST_CHECK_EQUAL(
        zhangAlternativePhysicalPathCount(k22, {"R0", g01}), 1
    );
    ZhangSatelliteSupportMetrics redundant =
        zhangSatelliteSupportMetrics(k22);
    BOOST_REQUIRE_EQUAL(redundant.supportCounts.size(), 1);
    BOOST_CHECK_EQUAL(redundant.supportCounts.at({g01, g02}), 2);
    BOOST_CHECK_EQUAL(redundant.bridgeEdges.size(), 1);
    BOOST_CHECK_EQUAL(redundant.edgeConnectivity, 2);

    std::set<ZhangGraphEdge> chain = k22;
    chain.insert({"R2", g02});
    chain.insert({"R2", g03});
    ZhangSatelliteSupportMetrics metrics =
        zhangSatelliteSupportMetrics(chain);
    BOOST_CHECK_EQUAL(metrics.satellites.size(), 3);
    BOOST_CHECK_EQUAL(metrics.supportCounts.at({g01, g02}), 2);
    BOOST_CHECK_EQUAL(metrics.supportCounts.at({g02, g03}), 1);
    BOOST_CHECK_EQUAL(metrics.bridgeEdges.size(), 2);
    BOOST_CHECK_EQUAL(metrics.minimumSupport, 1);
    BOOST_CHECK_EQUAL(metrics.maximumSupport, 2);
    BOOST_CHECK_EQUAL(metrics.edgeConnectivity, 1);

    std::set<ZhangGraphEdge> tree = {
        {"R0", g01}, {"R0", g02}, {"R1", g02}
    };
    BOOST_CHECK_EQUAL(
        zhangAlternativePhysicalPathCount(tree, {"R0", g02}), 0
    );
}

BOOST_AUTO_TEST_CASE(promoted_satellite_relation_survives_source_arc_retirement)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(
        g01, g03, 2, "four_physical_arcs", true
    ));
    manager.retireUnprovedBridges({g03});

    long long difference = 0;
    BOOST_CHECK(manager.relation(g01, g03, difference));
    BOOST_CHECK_EQUAL(difference, 2);
    BOOST_CHECK_EQUAL(manager.relationCount(), 1);
}

BOOST_AUTO_TEST_CASE(frontend_integer_gauge_birth_defines_zero_kappa_component)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1C);
	// Estimator warm-up may replace a physical segment before the broadcast
	// frontend has a fixed integer from which to define its t0.
	manager.recordSatelliteDiscontinuity(g03);

    const auto initial = manager.initialiseFrontendGaugeComponent(
        {g03, g01, g02});
    BOOST_REQUIRE_MESSAGE(initial.accepted, initial.reason);
    BOOST_CHECK_EQUAL(initial.satelliteCount, 3);
    BOOST_CHECK_EQUAL(initial.relationCount, 2);

    long long difference = 99;
    BOOST_REQUIRE(manager.relation(g01, g02, difference));
    BOOST_CHECK_EQUAL(difference, 0);
    BOOST_REQUIRE(manager.relation(g01, g03, difference));
    BOOST_CHECK_EQUAL(difference, 0);
    const auto status = manager.status(g03, true);
    BOOST_CHECK_EQUAL(status.componentSize, 3);
    BOOST_CHECK_EQUAL(status.componentRank, 2);
    BOOST_CHECK(status.integerDatumContinuous);
    BOOST_CHECK(status.integerPrecisionValid);
    BOOST_CHECK(status.integerValid);

    const auto repeated = manager.initialiseFrontendGaugeComponent(
        {g01, g02, g03});
    BOOST_CHECK(!repeated.accepted);
    BOOST_CHECK_EQUAL(
        repeated.reason, "FRONTEND_GAUGE_ALREADY_INITIALISED");
    BOOST_CHECK_EQUAL(manager.relationCount(), 2);
}

BOOST_AUTO_TEST_CASE(product_support_path_switch_preserves_component_and_value)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g02, 5, "path_p1", true));
    auto before = manager.status(g02, true);
    double rawBefore = 0.37;
    double productBefore = rawBefore + 0.19 * before.alignmentCycles;

    // Retiring p1 is a provenance event only; p2 proves the same relation.
    manager.retireUnprovedBridges({g02});
    BOOST_REQUIRE(manager.promoteRelation(g01, g02, 5, "path_p2", true));
    auto after = manager.status(g02, true);
    double productAfter = rawBefore + 0.19 * after.alignmentCycles;

    BOOST_CHECK_EQUAL(after.datumVersion, before.datumVersion);
    BOOST_CHECK_EQUAL(after.componentId, before.componentId);
    BOOST_CHECK_SMALL(productAfter - productBefore, 1e-15);
}

BOOST_AUTO_TEST_CASE(detached_subtree_keeps_internal_promoted_relations)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(
        g01, g02, 7, "unproved_dynamic_bridge", false
    ));
    BOOST_REQUIRE(manager.promoteRelation(
        g02, g03, -3, "promoted_subtree_relation", true
    ));

    manager.retireUnprovedBridges({g02, g03});
    long long difference = 0;
    BOOST_CHECK(!manager.relation(g01, g02, difference));
    BOOST_CHECK(manager.relation(g02, g03, difference));
    BOOST_CHECK_EQUAL(difference, -3);
    BOOST_CHECK_EQUAL(manager.status(g01, true).componentSize, 1);
    BOOST_CHECK_EQUAL(manager.status(g02, true).componentSize, 2);
}

BOOST_AUTO_TEST_CASE(inconsistent_satellite_integer_bridge_is_rejected)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g02, 2, "edge_12"));
    BOOST_REQUIRE(manager.promoteRelation(g02, g03, 4, "edge_23"));
    BOOST_CHECK(!manager.promoteRelation(g01, g03, 7, "bad_cycle"));
    BOOST_CHECK_EQUAL(manager.conflicts(), 1);
    long long difference = 0;
    BOOST_REQUIRE(manager.relation(g01, g03, difference));
    BOOST_CHECK_EQUAL(difference, 6);
}

BOOST_AUTO_TEST_CASE(satellite_product_events_distinguish_topology_progress)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    SatSys g12(E_Sys::GPS, 12);
    SatSys g25(E_Sys::GPS, 25);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);

    auto edgeA = manager.promoteRelationDetailed(
        g01, g03, 9, "component_a"
    );
    BOOST_REQUIRE(edgeA.accepted);
    BOOST_CHECK(edgeA.type ==
        ZhangProductRelationEventType::NEW_COMPONENT_EDGE);

    auto edgeB = manager.promoteRelationDetailed(
        g12, g25, 150, "component_b"
    );
    BOOST_REQUIRE(edgeB.accepted);
    BOOST_CHECK(edgeB.type ==
        ZhangProductRelationEventType::NEW_COMPONENT_EDGE);

    auto merge = manager.promoteRelationDetailed(
        g03, g12, -21, "bridge"
    );
    BOOST_REQUIRE(merge.accepted);
    BOOST_CHECK(merge.type ==
        ZhangProductRelationEventType::COMPONENT_MERGE);
    BOOST_CHECK_EQUAL(merge.oldComponentSizeA, 2);
    BOOST_CHECK_EQUAL(merge.oldComponentSizeB, 2);
    BOOST_CHECK_EQUAL(merge.newComponentSize, 4);

    auto confirmation = manager.promoteRelationDetailed(
        g01, g25, 138, "redundant_path"
    );
    BOOST_REQUIRE(confirmation.accepted);
    BOOST_CHECK(confirmation.type ==
        ZhangProductRelationEventType::REDUNDANT_CONFIRMATION);

    auto conflict = manager.promoteRelationDetailed(
        g01, g25, 139, "bad_cycle"
    );
    BOOST_CHECK(!conflict.accepted);
    BOOST_CHECK(conflict.type ==
        ZhangProductRelationEventType::CONFLICT_REJECTED);
    BOOST_CHECK_EQUAL(manager.eventCount(
        ZhangProductRelationEventType::NEW_COMPONENT_EDGE), 2);
    BOOST_CHECK_EQUAL(manager.eventCount(
        ZhangProductRelationEventType::COMPONENT_MERGE), 1);
    BOOST_CHECK_EQUAL(manager.eventCount(
        ZhangProductRelationEventType::REDUNDANT_CONFIRMATION), 1);
}

BOOST_AUTO_TEST_CASE(hybrid_broadcast_component_metadata_is_persistent_and_exact)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g03(E_Sys::GPS, 3);
	SatSys g12(E_Sys::GPS, 12);
	ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1C);
	BOOST_CHECK_EQUAL(
		zhangHybridPhaseProductSegmentId(g03, E_ObsCode::L1C, 0),
		"G03-L1C-SEG0");
	BOOST_REQUIRE(manager.promoteRelation(g01, g03, 4, "edge_13"));
	BOOST_REQUIRE(manager.promoteRelation(g03, g12, -7, "edge_312"));
	BOOST_REQUIRE(manager.promoteRelation(g01, g12, -3, "cycle_112"));
	const auto before = manager.status(g03, true);
	BOOST_CHECK(before.integerValid);
	BOOST_CHECK_EQUAL(before.componentSize, 3);
	BOOST_CHECK_EQUAL(before.componentRank, 2);
	BOOST_CHECK_EQUAL(before.certifiedRelationCount, 3);
	BOOST_CHECK_EQUAL(before.redundantRelationCount, 1);
	BOOST_CHECK(before.cycleClosureValid);
	BOOST_CHECK_GT(before.componentVersion, 0);

	// An exact backend coordinate shift transports integer potentials but is
	// not a frontend physical-segment event.
	manager.applyDynamicTreeShift(g03, 2);
	const auto transported = manager.status(g03, true);
	BOOST_CHECK_EQUAL(
		transported.componentVersion, before.componentVersion);
	BOOST_CHECK_EQUAL(
		transported.alignmentGeneration, before.alignmentGeneration);
	BOOST_CHECK_EQUAL(transported.phaseSegment, before.phaseSegment);

	manager.recordSatelliteDiscontinuity(g03);
	const auto after = manager.status(g03, true);
	BOOST_CHECK_GT(after.componentVersion, before.componentVersion);
	BOOST_CHECK_GT(after.alignmentGeneration, before.alignmentGeneration);
	BOOST_CHECK_EQUAL(after.phaseSegment, before.phaseSegment + 1);
	BOOST_CHECK_EQUAL(after.componentSize, 1);
	BOOST_CHECK(!after.integerValid);
}

BOOST_AUTO_TEST_CASE(local_fractional_alignment_loss_can_relink_to_component_anchor)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    SatSys g22(E_Sys::GPS, 22);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g03, 9, "fixed_13"));
    BOOST_REQUIRE(manager.promoteRelation(g01, g22, -3, "fixed_122"));

    auto before = manager.status(g03, true);
    auto preserved = manager.applyDynamicTreeTransform({
        {g01, 10.0}, {g03, 10.25}, {g22, 7.0}
    });
    BOOST_CHECK(preserved.at(g01));
    BOOST_CHECK(!preserved.at(g03));
    BOOST_CHECK(preserved.at(g22));
    BOOST_CHECK(manager.alignmentState(g03) ==
        ZhangCurrentAlignmentState::CURRENT_ALIGNMENT_PENDING);
    BOOST_CHECK(manager.status(g01, true).integerDatumContinuous);
    BOOST_CHECK(!manager.status(g03, true).integerDatumContinuous);

    auto relink = manager.realignRelation(
        g01, g03, 11, "same_component_relink"
    );
    BOOST_REQUIRE(relink.accepted);
    BOOST_CHECK(relink.type ==
        ZhangProductRelationEventType::CURRENT_REALIGNMENT);
    BOOST_CHECK(manager.alignmentState(g03) ==
        ZhangCurrentAlignmentState::CURRENT_ALIGNMENT_VALID);
    auto after = manager.status(g03, true);
    BOOST_CHECK(after.integerDatumContinuous);
    BOOST_CHECK_EQUAL(after.datumVersion, before.datumVersion);
    BOOST_CHECK_EQUAL(after.discontinuityCounter, before.discontinuityCounter);
    long long currentDifference = 0;
    BOOST_REQUIRE(manager.relation(g01, g03, currentDifference));
    BOOST_CHECK_EQUAL(currentDifference, 11);
}

BOOST_AUTO_TEST_CASE(attaching_left_singleton_preserves_established_component_alignment)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    SatSys g11(E_Sys::GPS, 11);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g03, 9, "established"));
    auto g01Before = manager.status(g01, true);
    auto g03Before = manager.status(g03, true);

    BOOST_REQUIRE(manager.promoteRelation(
        g11, g01, -12, "left_singleton_attachment"
    ));
    auto g01After = manager.status(g01, true);
    auto g03After = manager.status(g03, true);
    auto g11After = manager.status(g11, true);
    BOOST_CHECK_EQUAL(
        g01After.alignmentCycles, g01Before.alignmentCycles
    );
    BOOST_CHECK_EQUAL(
        g03After.alignmentCycles, g03Before.alignmentCycles
    );
    BOOST_CHECK_EQUAL(g11After.alignmentCycles, 12);
}

BOOST_AUTO_TEST_CASE(conflicting_current_relation_is_quarantined_then_relinked)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g07(E_Sys::GPS, 7);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g07, 12, "initial"));
    auto before = manager.status(g07, true);

    auto quarantined = manager.quarantineCurrentAlignment(
        g01, g07, g01
    );
    BOOST_CHECK(quarantined.type ==
        ZhangProductRelationEventType::CURRENT_ALIGNMENT_QUARANTINED);
    BOOST_CHECK(quarantined.quarantinedSatellite == g07);
    BOOST_CHECK(!manager.status(g07, true).integerValid);
    BOOST_CHECK(manager.status(g01, true).integerValid);

    auto relink = manager.realignRelation(
        g01, g07, 7559, "confirmed_current_coordinate"
    );
    BOOST_REQUIRE(relink.accepted);
    BOOST_CHECK(relink.type ==
        ZhangProductRelationEventType::CURRENT_REALIGNMENT);
    auto after = manager.status(g07, true);
    BOOST_CHECK(after.integerValid);
    BOOST_CHECK_EQUAL(after.datumVersion, before.datumVersion);
    BOOST_CHECK_EQUAL(after.discontinuityCounter, before.discontinuityCounter);
    long long currentDifference = 0;
    BOOST_REQUIRE(manager.relation(g01, g07, currentDifference));
    BOOST_CHECK_EQUAL(currentDifference, 7559);
}

BOOST_AUTO_TEST_CASE(held_support_quarantine_preserves_trusted_anchor)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g07(E_Sys::GPS, 7);
    SatSys g23(E_Sys::GPS, 23);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g07, 12, "initial_07"));
    BOOST_REQUIRE(manager.promoteRelation(g01, g23, -94, "initial_23"));

    std::set<SatSys> support{g01, g07, g23};
    BOOST_CHECK_EQUAL(
        manager.quarantineCurrentAlignments(support, g01), 2
    );
    BOOST_CHECK(manager.status(g01, true).integerValid);
    BOOST_CHECK(!manager.status(g07, true).integerValid);
    BOOST_CHECK(!manager.status(g23, true).integerValid);
    long long persistentDifference = 0;
    BOOST_REQUIRE(manager.relation(g01, g23, persistentDifference));
    BOOST_CHECK_EQUAL(persistentDifference, -94);
}

BOOST_AUTO_TEST_CASE(certified_temporal_batch_restores_quarantined_frontend)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g07(E_Sys::GPS, 7);
	SatSys g23(E_Sys::GPS, 23);
	ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
	BOOST_REQUIRE(manager.promoteRelation(g01, g07, 12, "initial_07"));
	BOOST_REQUIRE(manager.promoteRelation(g01, g23, -94, "initial_23"));
	BOOST_REQUIRE_EQUAL(
		manager.quarantineCurrentAlignments({g07, g23}, g01), 2);

	// The common +2 cycle raw shift is an unobservable component gauge.
	// Relative changes are +5 for G07 and -3 for G23, so kappa must change
	// by -5 and +3 respectively to keep raw+lambda*kappa invariant.
	const auto restored = manager.applyCertifiedTemporalTransform({
		{g01, 2}, {g07, 7}, {g23, -1}});
	BOOST_REQUIRE_MESSAGE(restored.accepted, restored.reason);
	BOOST_CHECK_EQUAL(restored.restoredSatellites, 2);
	BOOST_CHECK(manager.status(g07, true).integerValid);
	BOOST_CHECK(manager.status(g23, true).integerValid);
	long long difference = 0;
	BOOST_REQUIRE(manager.relation(g01, g07, difference));
	BOOST_CHECK_EQUAL(difference, 7);
	BOOST_REQUIRE(manager.relation(g01, g23, difference));
	BOOST_CHECK_EQUAL(difference, -91);
}

BOOST_AUTO_TEST_CASE(certified_temporal_batch_fails_without_aligned_anchor)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g07(E_Sys::GPS, 7);
	ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
	BOOST_REQUIRE(manager.promoteRelation(g01, g07, 12, "initial"));
	BOOST_REQUIRE_EQUAL(
		manager.quarantineCurrentAlignments({g01, g07}), 2);
	const auto before = manager.checkpointState();
	const auto rejected = manager.applyCertifiedTemporalTransform({{g07, 5}});
	BOOST_CHECK(!rejected.accepted);
	BOOST_CHECK_EQUAL(rejected.reason, "NO_ALIGNED_COMPONENT_ANCHOR");
	BOOST_CHECK(manager.checkpointState().alignmentCycles == before.alignmentCycles);
	BOOST_CHECK(manager.checkpointState().alignmentKnown == before.alignmentKnown);
}

BOOST_AUTO_TEST_CASE(dynamic_tree_integer_changes_leave_product_invariant)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g02, 3, "initial_target"));
    constexpr double lambda = 0.190293672798365;
    std::map<SatSys, double> raw = {{g01, 0.2}, {g02, -0.4}};
    std::map<SatSys, double> product;
    for (const auto& [satellite, value] : raw)
    {
        product[satellite] = value +
            lambda * manager.status(satellite, true).alignmentCycles;
    }
    auto component = manager.status(g02, true).componentId;

    for (const auto& [satellite, stateJump] :
         std::map<SatSys, long long>{{g01, 4}, {g02, -5}})
    {
        raw[satellite] += lambda * stateJump;
        manager.applyDynamicTreeShift(satellite, -stateJump);
    }
    for (const auto& [satellite, value] : raw)
    {
        double transformedProduct = value +
            lambda * manager.status(satellite, true).alignmentCycles;
        BOOST_CHECK_SMALL(transformedProduct - product.at(satellite), 1e-14);
        BOOST_CHECK_EQUAL(manager.status(satellite, true).datumVersion, 0);
        BOOST_CHECK_EQUAL(manager.status(satellite, true).componentId, component);
    }
}

BOOST_AUTO_TEST_CASE(component_common_fractional_gauge_preserves_integer_datum)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    SatSys g22(E_Sys::GPS, 22);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g03, 9, "fixed_13"));
    BOOST_REQUIRE(manager.promoteRelation(g01, g22, -3, "fixed_122"));

    auto before = manager.status(g03, true);
    auto preserved = manager.applyDynamicTreeTransform({
        {g01, 111.76080459977078},
        {g03, 114.76080459977078},
        {g22, 107.76080459977078},
    });
    BOOST_CHECK(preserved.at(g01));
    BOOST_CHECK(preserved.at(g03));
    BOOST_CHECK(preserved.at(g22));
    auto after = manager.status(g03, true);
    BOOST_CHECK(after.integerDatumContinuous);
    BOOST_CHECK_EQUAL(after.datumVersion, before.datumVersion);
    BOOST_CHECK_EQUAL(after.componentId, before.componentId);
    BOOST_CHECK_EQUAL(after.alignmentCycles - before.alignmentCycles, 3);
}

BOOST_AUTO_TEST_CASE(hybrid_tree_invariance_is_relative_then_real_gauge)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    SatSys g22(E_Sys::GPS, 22);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g03, 9, "fixed_13"));
    BOOST_REQUIRE(manager.promoteRelation(g01, g22, -3, "fixed_122"));
    constexpr double lambda = 0.190293672798365;
    const std::map<SatSys, double> cycleChanges = {
        {g01, 111.76080459977078},
        {g03, 114.76080459977078},
        {g22, 107.76080459977078}
    };
    std::map<SatSys, ZhangSatelliteDatumStatus> before;
    for (const auto& [satellite, ignored] : cycleChanges)
    {
        before[satellite] = manager.status(satellite, true);
    }
    const auto preserved = manager.applyDynamicTreeTransform(cycleChanges);
    std::vector<ZhangHybridTreeTransformSample> samples;
    for (const auto& [satellite, change] : cycleChanges)
    {
        const auto after = manager.status(satellite, true);
        const auto& old = before.at(satellite);
        samples.push_back({
            satellite, old.componentId, after.componentId, -lambda * change,
            old.alignmentCycles, after.alignmentCycles,
            old.phaseSegment, after.phaseSegment,
            old.datumVersion, after.datumVersion,
            old.componentVersion, after.componentVersion,
            old.alignmentGeneration, after.alignmentGeneration,
            preserved.at(satellite)
        });
    }
    const auto audit = zhangHybridTreeTransformInvariance(samples, lambda);
    BOOST_REQUIRE_EQUAL(audit.size(), 3);
    const double expectedCommon = -lambda * cycleChanges.at(g01);
    for (const auto& row : audit)
    {
        BOOST_CHECK(row.valid);
        BOOST_CHECK(row.invariant);
        BOOST_CHECK_EQUAL(row.reason, "INVARIANT");
        BOOST_CHECK_EQUAL(row.componentSupportCount, 3);
        BOOST_CHECK_SMALL(row.relativeFrontendDeltaMetres, 1e-12);
        BOOST_CHECK_SMALL(
            row.componentCommonDeltaMetres - expectedCommon, 1e-12);
        BOOST_CHECK_SMALL(
            row.expectedRealGaugeShiftMetres + expectedCommon, 1e-12);
		BOOST_CHECK(row.hybridClosureMachineZero);
		BOOST_CHECK_SMALL(row.hybridClosureResidualMetres, 1e-12);
    }

    // The second half of the frontend absorbs the one remaining component
    // common mode and leaves the broadcast products exactly unchanged.
    ZhangHybridRealGaugeTransport gauge;
    const std::vector<SatSys> satellites = {g01, g03, g22};
    const std::vector<std::string> segments = {
        "G01-L1W-SEG0", "G03-L1W-SEG0", "G22-L1W-SEG0"};
    VectorXd oldFrontend(3);
    oldFrontend << 0.2, -0.4, 0.2;
    const MatrixXd covariance = MatrixXd::Identity(3, 3);
    const auto first = gauge.transport(
        satellites, segments, oldFrontend, covariance);
    BOOST_REQUIRE(first.valid);
    const VectorXd shiftedFrontend =
        oldFrontend + VectorXd::Constant(3, expectedCommon);
    const auto second = gauge.transport(
        satellites, segments, shiftedFrontend, covariance);
    BOOST_REQUIRE(second.valid);
    BOOST_CHECK_EQUAL(second.overlapCount, 3);
    BOOST_CHECK_SMALL(
        second.commonShiftMetres + expectedCommon, 1e-12);
    BOOST_CHECK_SMALL((second.values - first.values).norm(), 1e-12);

    // One member losing integer alignment must fail closed even when the
    // remaining members still share a valid common real gauge.
    samples[1].alignmentPreserved = false;
    const auto suspended = zhangHybridTreeTransformInvariance(samples, lambda);
    BOOST_CHECK(!suspended[1].invariant);
    BOOST_CHECK_EQUAL(suspended[1].reason, "ALIGNMENT_SUSPENDED");
}

BOOST_AUTO_TEST_CASE(hybrid_pure_s_basis_event_closes_integer_and_real_gauges)
{
	constexpr double lambda = 0.190293672798365;
	constexpr long long backendIntegerGauge = 3;
	constexpr double gamma = -0.047;
	const double backendDelta = lambda * backendIntegerGauge + gamma;
	const auto closure = zhangHybridPureCoordinateClosure(
		lambda, backendIntegerGauge, gamma, backendDelta, 1e-12);
	BOOST_REQUIRE(closure.valid);
	BOOST_CHECK(closure.machineZero);
	BOOST_CHECK_EQUAL(closure.reason, "PURE_COORDINATE_INVARIANT");
	BOOST_CHECK_SMALL(closure.hybridResidualMetres, 1e-12);
	BOOST_CHECK_SMALL(
		closure.integerCompensationMetres + lambda * backendIntegerGauge,
		1e-12);
	BOOST_CHECK_SMALL(closure.realGaugeCompensationMetres + gamma, 1e-12);

	const auto rejected = zhangHybridPureCoordinateClosure(
		lambda, backendIntegerGauge, gamma, backendDelta + 1e-3, 1e-12);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(
		rejected.reason, "BACKEND_INTEGER_REAL_DECOMPOSITION_MISMATCH");
}

BOOST_AUTO_TEST_CASE(real_gauge_transport_audit_separates_coordinate_and_time)
{
	Vector2d raw;
	raw << 0.31, 0.29;
	Vector2d integerRemoved = Vector2d::Zero();
	Matrix2d posterior;
	posterior << 4, 1, 1, 2;
	auto pure = zhangAuditRealGaugeTransport(
		ZhangRealGaugeTransportEventKind::
			SAME_POSTERIOR_COORDINATE_TRANSFORM,
		raw, integerRemoved, posterior, posterior, posterior);
	BOOST_REQUIRE_MESSAGE(pure.valid, pure.failureReason);
	BOOST_CHECK(pure.samePosteriorEvent);
	BOOST_CHECK(pure.differenceCovarianceMachineZero);
	BOOST_CHECK_SMALL(pure.glsShiftVariance, 1e-15);
	BOOST_CHECK_SMALL(pure.realShiftMetres - 0.30, 1e-15);

	Matrix2d oldCovariance;
	oldCovariance << 4, 1, 1, 3;
	Matrix2d newCovariance;
	newCovariance << 5, 1.5, 1.5, 4;
	Matrix2d crossCovariance;
	crossCovariance << 3, 0.5, 0.75, 2;
	auto temporal = zhangAuditRealGaugeTransport(
		ZhangRealGaugeTransportEventKind::CROSS_EPOCH_TRANSPORT,
		raw, integerRemoved, oldCovariance, newCovariance,
		crossCovariance);
	BOOST_REQUIRE_MESSAGE(temporal.valid, temporal.failureReason);
	const Matrix2d expectedDifference = newCovariance + oldCovariance
		- crossCovariance - crossCovariance.transpose();
	BOOST_CHECK_SMALL(
		(temporal.differenceCovariance - expectedDifference).norm(), 1e-14);
	BOOST_CHECK_GT(temporal.glsShiftVariance, 0);

	// Marginals without Q-+ are not an admissible temporal audit input.
	auto missingCross = zhangAuditRealGaugeTransport(
		ZhangRealGaugeTransportEventKind::CROSS_EPOCH_TRANSPORT,
		raw, integerRemoved, oldCovariance, newCovariance, MatrixXd());
	BOOST_CHECK(!missingCross.valid);
	BOOST_CHECK_EQUAL(missingCross.failureReason,
		"REAL_GAUGE_AUDIT_REQUIRES_FULL_JOINT_MARGINAL");
}

BOOST_AUTO_TEST_CASE(hybrid_server_to_user_dual_frequency_integer_theorem_closes)
{
	auto closure = zhangHybridUserIntegerClosure(
		ZhangExactInteger(105), ZhangExactInteger(77),
		ZhangExactInteger(88), ZhangExactInteger(61),
		ZhangExactInteger(12), ZhangExactInteger(9),
		ZhangExactInteger(-4), ZhangExactInteger(-8),
		true, true);
	BOOST_REQUIRE_MESSAGE(closure.valid, closure.failureReason);
	BOOST_CHECK_EQUAL(closure.firstSignalSatelliteSd, 31);
	BOOST_CHECK_EQUAL(closure.secondSignalSatelliteSd, 31);
	BOOST_CHECK_EQUAL(closure.wideLaneSatelliteSd, 0);
	BOOST_CHECK(closure.admissibleDualFrequencyTransform);
	BOOST_CHECK(closure.exactInverseClosure);

	auto disconnected = zhangHybridUserIntegerClosure(
		1, 0, 1, 0, 0, 0, 0, 0, false, true);
	BOOST_CHECK(!disconnected.valid);
	BOOST_CHECK_EQUAL(disconnected.failureReason,
		"USER_SATELLITES_NOT_IN_SAME_INTEGER_COMPONENT");

	auto uncertified = zhangHybridUserIntegerClosure(
		1, 0, 1, 0, 0, 0, 0, 0, true, false);
	BOOST_CHECK(!uncertified.valid);
	BOOST_CHECK_EQUAL(uncertified.failureReason,
		"SERVER_INTEGER_RELATION_NOT_CERTIFIED");
}

BOOST_AUTO_TEST_CASE(product_reference_exchange_preserves_integer_relations)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    SatSys g22(E_Sys::GPS, 22);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g03, 9, "fixed_13"));
    BOOST_REQUIRE(manager.promoteRelation(g01, g22, -3, "fixed_122"));

    long long fromG01ToG22 = 0;
    long long fromG03ToG22 = 0;
    BOOST_REQUIRE(manager.relation(g01, g22, fromG01ToG22));
    BOOST_REQUIRE(manager.relation(g03, g22, fromG03ToG22));
    BOOST_CHECK_EQUAL(fromG01ToG22, -3);
    BOOST_CHECK_EQUAL(fromG03ToG22, -12);
    BOOST_CHECK_EQUAL(fromG01ToG22 - 9, fromG03ToG22);
}

BOOST_AUTO_TEST_CASE(only_satellite_discontinuity_changes_product_version)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    ZhangSatelliteDatumManager manager(E_Sys::GPS, E_ObsCode::L1W);
    BOOST_REQUIRE(manager.promoteRelation(g01, g03, 2, "fixed_relation"));
    auto initial = manager.status(g03, true);

    manager.applyDynamicTreeTransform({{g01, 4.25}, {g03, 7.25}});
    manager.markDynamicAlignmentUnknown({g03});
    auto dynamic = manager.status(g03, true);
    BOOST_CHECK_EQUAL(dynamic.datumVersion, initial.datumVersion);
    BOOST_CHECK_EQUAL(dynamic.phaseSegment, initial.phaseSegment);
    BOOST_CHECK_EQUAL(dynamic.discontinuityCounter, initial.discontinuityCounter);
    BOOST_CHECK(dynamic.integerDatumContinuous);

    manager.recordSatelliteDiscontinuity(g03);
    auto discontinuous = manager.status(g03, true);
    BOOST_CHECK_EQUAL(discontinuous.datumVersion, initial.datumVersion + 1);
    BOOST_CHECK_EQUAL(discontinuous.phaseSegment, initial.phaseSegment + 1);
    BOOST_CHECK_EQUAL(
        discontinuous.discontinuityCounter,
        initial.discontinuityCounter + 1
    );
    BOOST_CHECK(!discontinuous.integerDatumContinuous);
}

BOOST_AUTO_TEST_CASE(canonical_product_relations_reject_silent_satellite_substitution)
{
	const SatSys g01(E_Sys::GPS, 1);
	const SatSys g02(E_Sys::GPS, 2);
	const SatSys g03(E_Sys::GPS, 3);
	const SatSys g05(E_Sys::GPS, 5);
	const SatSys g07(E_Sys::GPS, 7);
	ZhangPersistentProductDatumRegistry registry;
	const auto initial = registry.selectRelations(
		E_Sys::GPS,
		{{g01, g02}, {g01, g03}, {g01, g05}},
		{g01, g02, g03, g05}, 3);
	BOOST_REQUIRE(initial.established);
	BOOST_REQUIRE_EQUAL(initial.selected.size(), 3);
	BOOST_CHECK_EQUAL(initial.selected[0].id(), "G01->G02");

	// The current graph proposes G07 after G02 becomes unavailable.  The
	// persistent registry keeps G01-G02 missing and must not replace it.
	const auto changedGraph = registry.selectRelations(
		E_Sys::GPS,
		{{g01, g03}, {g01, g05}, {g01, g07}},
		{g01, g03, g05, g07}, 3);
	BOOST_CHECK(changedGraph.silentSubstitutionRejected);
	BOOST_REQUIRE_EQUAL(changedGraph.selected.size(), 2);
	BOOST_REQUIRE_EQUAL(changedGraph.missing.size(), 1);
	BOOST_CHECK_EQUAL(changedGraph.missing[0].id(), "G01->G02");
	BOOST_REQUIRE_EQUAL(changedGraph.ignoredSubstitutes.size(), 1);
	BOOST_CHECK_EQUAL(changedGraph.ignoredSubstitutes[0].id(), "G01->G07");

	// Reversed reference orientation is the same canonical coordinate.
	const auto restored = registry.selectRelations(
		E_Sys::GPS,
		{{g02, g01}, {g03, g01}, {g05, g01}},
		{g01, g02, g03, g05}, 3);
	BOOST_CHECK(!restored.silentSubstitutionRejected);
	BOOST_CHECK(restored.selected == initial.selected);
}

BOOST_AUTO_TEST_CASE(l1c_and_l2w_product_datum_versions_are_independent)
{
	const SatSys g01(E_Sys::GPS, 1);
	const SatSys g03(E_Sys::GPS, 3);
	const auto relation = ZhangCanonicalSatelliteRelation::ordered(g01, g03);
	ZhangPersistentProductDatumRegistry registry;
	const auto l1Initial = registry.observe(
		E_Sys::GPS, E_ObsCode::L1C, relation, 0, 0, 0, 0, true);
	const auto l2Initial = registry.observe(
		E_Sys::GPS, E_ObsCode::L2W, relation, 0, 0, 0, 0, true);
	BOOST_REQUIRE(l1Initial.valid && l2Initial.valid);
	BOOST_CHECK_EQUAL(l1Initial.version, 0);
	BOOST_CHECK_EQUAL(l2Initial.version, 0);
	BOOST_CHECK_NE(l1Initial.productDatumId, l2Initial.productDatumId);

	// Temporary loss of absolute observability is quotient-only, not a new
	// datum version.
	const auto l1Quotient = registry.observe(
		E_Sys::GPS, E_ObsCode::L1C, relation, 0, 0, 0, 0, false);
	BOOST_CHECK(l1Quotient.quotientOnly);
	BOOST_CHECK(!l1Quotient.absoluteValid);
	BOOST_CHECK(!l1Quotient.versionChanged);
	BOOST_CHECK_EQUAL(l1Quotient.version, 0);

	// A real L1C endpoint discontinuity advances L1C only.
	const auto l1Changed = registry.observe(
		E_Sys::GPS, E_ObsCode::L1C, relation, 0, 1, 0, 1, false);
	const auto l2Unchanged = registry.observe(
		E_Sys::GPS, E_ObsCode::L2W, relation, 0, 0, 0, 0, true);
	BOOST_CHECK(l1Changed.versionChanged);
	BOOST_CHECK_EQUAL(l1Changed.version, 1);
	BOOST_CHECK_EQUAL(l2Unchanged.version, 0);
}

BOOST_AUTO_TEST_CASE(wide_lane_alignment_transport_is_invariant_across_s_basis_changes)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g03(E_Sys::GPS, 3);
    ZhangSatelliteDatumManager l1(E_Sys::GPS, E_ObsCode::L1W);
    ZhangSatelliteDatumManager l2(E_Sys::GPS, E_ObsCode::L2W);
    BOOST_REQUIRE(l1.promoteRelation(g01, g03, 5, "l1_relation"));
    BOOST_REQUIRE(l2.promoteRelation(g01, g03, 2, "l2_relation"));

    long long l1Before = 0;
    long long l2Before = 0;
    BOOST_REQUIRE(l1.relation(g01, g03, l1Before));
    BOOST_REQUIRE(l2.relation(g01, g03, l2Before));
    const auto l1AnchorBefore = l1.status(g01, true);
    const auto l1SatelliteBefore = l1.status(g03, true);
    const auto l2AnchorBefore = l2.status(g01, true);
    const auto l2SatelliteBefore = l2.status(g03, true);
    const long long persistentWideLaneBefore =
        (l1Before - l2Before) -
        (l1SatelliteBefore.alignmentCycles -
            l1AnchorBefore.alignmentCycles) +
        (l2SatelliteBefore.alignmentCycles -
            l2AnchorBefore.alignmentCycles);

    l1.applyDynamicTreeTransform({{g01, 4.25}, {g03, 7.25}});
    l2.applyDynamicTreeTransform({{g01, -1.4}, {g03, 0.6}});
    long long l1After = 0;
    long long l2After = 0;
    BOOST_REQUIRE(l1.relation(g01, g03, l1After));
    BOOST_REQUIRE(l2.relation(g01, g03, l2After));
    const auto l1Anchor = l1.status(g01, true);
    const auto l1Satellite = l1.status(g03, true);
    const auto l2Anchor = l2.status(g01, true);
    const auto l2Satellite = l2.status(g03, true);
    const long long persistentWideLaneAfter =
        (l1After - l2After) -
        (l1Satellite.alignmentCycles - l1Anchor.alignmentCycles) +
        (l2Satellite.alignmentCycles - l2Anchor.alignmentCycles);

    BOOST_CHECK_NE(l1After - l2After, persistentWideLaneBefore);
    BOOST_CHECK_EQUAL(persistentWideLaneAfter, persistentWideLaneBefore);
    BOOST_CHECK_EQUAL(l1Satellite.phaseSegment, 0);
    BOOST_CHECK_EQUAL(l2Satellite.phaseSegment, 0);
}

BOOST_AUTO_TEST_CASE(product_constraint_promotion_requires_exact_named_membership)
{
    auto recovered = ProductConstraintPromotion::recoverNamedTargets(
        {{1, 1}, {0, 1}}, {7, 4}, 2
    );
    BOOST_REQUIRE_EQUAL(recovered.size(), 2);
    BOOST_CHECK(recovered.at(0) == 3);
    BOOST_CHECK(recovered.at(1) == 4);

    auto unsaturated = ProductConstraintPromotion::recoverNamedTargets(
        {{2}}, {6}, 1
    );
    BOOST_CHECK(unsaturated.empty());
}

BOOST_AUTO_TEST_CASE(local_subtree_break_preserves_unaffected_product_continuity)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R1", g02}, {"R1", g03}
    };
    ZhangGraphBasis oldBasis = zhangBuildSpanningTree(edges, "R0");
    BOOST_REQUIRE(oldBasis.connected);

    edges.erase({"R0", g02});
    std::set<ZhangGraphEdge> retained = zhangRootComponentEdges(edges, "R0");
    ZhangGraphBasis rootBasis = zhangBuildSpanningTree(retained, "R0");
    BOOST_REQUIRE(rootBasis.connected);
    BOOST_CHECK_EQUAL(rootBasis.satellites.size(), 1);
    BOOST_CHECK(rootBasis.satellites.find(g01) != rootBasis.satellites.end());

    ZhangPhaseContinuityState unaffected;
    ZhangPhaseContinuityState detached;
    unaffected.markFixed();
    detached.markFixed();
    GTime time;
    detached.reinitialise(time, "local_subtree_break", 2);
    BOOST_CHECK(unaffected.integerValid());
    BOOST_CHECK_EQUAL(unaffected.counter, 0);
    BOOST_CHECK(!detached.integerValid());
    BOOST_CHECK_EQUAL(detached.counter, 1);
}

BOOST_AUTO_TEST_CASE(wide_lane_only_lattice_does_not_validate_individual_signals)
{
    ZhangExactMatrix held = {{1, -1}};
    ZhangDualSignalLatticeValidity validity =
        zhangClassifyDualSignalLattice(held, 1);
    BOOST_CHECK(!validity.l1);
    BOOST_CHECK(!validity.l2);
    BOOST_CHECK(validity.wideLane);

    ZhangIntegerLatticeMembership unsaturated =
        zhangIntegerRowLatticeContains({{2}}, {1});
    BOOST_CHECK(!unsaturated.contained);
    BOOST_REQUIRE_EQUAL(unsaturated.smithInvariants.size(), 1);
    BOOST_CHECK(unsaturated.smithInvariants.front() == 2);
}

BOOST_AUTO_TEST_CASE(wide_lane_plus_l1_lattice_recovers_both_signals)
{
    ZhangExactMatrix held = {{1, -1}, {1, 0}};
    ZhangDualSignalLatticeValidity validity =
        zhangClassifyDualSignalLattice(held, 1);
    BOOST_CHECK(validity.l1);
    BOOST_CHECK(validity.l2);
    BOOST_CHECK(validity.wideLane);
}

BOOST_AUTO_TEST_CASE(exact_row_hnf_removes_redundancy_and_tracks_integer_values)
{
    ZhangExactMatrix rows = {
        {2, 0},
        {0, 3},
        {2, 3},
        {4, 0}
    };
    ZhangExactVector values = {4, 6, 10, 8};
    ZhangExactRowHnf hnf = zhangExactRowHermiteNormalForm(rows, values);
    BOOST_REQUIRE(hnf.consistent);
    BOOST_REQUIRE_EQUAL(hnf.basis.size(), 2);
    BOOST_CHECK(zhangIntegerRowLatticeContains(hnf.basis, {2, 0}).contained);
    BOOST_CHECK(zhangIntegerRowLatticeContains(hnf.basis, {0, 3}).contained);
    BOOST_CHECK(!zhangIntegerRowLatticeContains(hnf.basis, {1, 0}).contained);

    ZhangExactRowHnf inconsistent = zhangExactRowHermiteNormalForm(
        {{1, 0}, {1, 0}},
        {2, 3}
    );
    BOOST_CHECK(!inconsistent.consistent);

    // Equation (17): membership must return the actual integer row
    // combination so the persistent-product shift can be evaluated exactly.
    ZhangIntegerLatticeMembership represented =
        zhangIntegerRowLatticeContains({{2, 1}, {1, 1}}, {3, 2});
    BOOST_REQUIRE(represented.contained);
    BOOST_REQUIRE_EQUAL(represented.combination.size(), 2);
    BOOST_CHECK(represented.combination == ZhangExactVector({1, 1}));
    ZhangExactVector heldValues = {5, 7};
    ZhangExactInteger shift = 0;
    for (std::size_t row = 0; row < represented.combination.size(); row++)
    {
        shift += represented.combination[row] * heldValues[row];
    }
    BOOST_CHECK(shift == 12);
}

BOOST_AUTO_TEST_CASE(exact_physical_row_hnf_tracks_a_basis_invariant_search_frame)
{
    const ZhangExactMatrix physicalRows = {
        { 1, -1,  0,  1,  0},
        { 0,  1, -1,  0,  1},
        { 0,  0,  0,  1, -1},
    };
    const ZhangExactMatrix unimodular = {
        {1,  1, 0},
        {0,  1, 1},
        {0, -1, 0},
    };
    const ZhangExactMatrix rebasedRows = zhangExactMultiply(
        unimodular, physicalRows);

    const ZhangExactRowHnf base = zhangExactRowHermiteNormalForm(
        physicalRows, {}, true);
    const ZhangExactRowHnf rebased = zhangExactRowHermiteNormalForm(
        rebasedRows, {}, true);
    BOOST_REQUIRE(base.consistent);
    BOOST_REQUIRE(rebased.consistent);
    BOOST_REQUIRE_EQUAL(base.basis.size(), physicalRows.size());
    BOOST_REQUIRE_EQUAL(rebased.basis.size(), physicalRows.size());
    BOOST_CHECK(base.basis == rebased.basis);
    BOOST_CHECK(
        zhangExactMultiply(base.rowTransform, physicalRows) == base.basis);
    BOOST_CHECK(
        zhangExactMultiply(rebased.rowTransform, rebasedRows) ==
        rebased.basis);
    BOOST_CHECK(
        zhangExactMultiply(rebased.rowTransform, unimodular) ==
        base.rowTransform);

    const ZhangExactVector currentState = {7, -2, 5};
    const ZhangExactVector rebasedState =
        zhangExactMatrixTimesColumn(unimodular, currentState);
    BOOST_CHECK(
        zhangExactMatrixTimesColumn(base.rowTransform, currentState) ==
        zhangExactMatrixTimesColumn(
            rebased.rowTransform, rebasedState));
}

BOOST_AUTO_TEST_CASE(iar_gain_low_rank_covariance_matches_dense_conditioning)
{
    Matrix4d covariance;
    covariance <<
        4.0, 0.8, 0.3, 0.1,
        0.8, 3.0, 0.4, 0.2,
        0.3, 0.4, 2.0, 0.5,
        0.1, 0.2, 0.5, 1.5;
    Matrix<double, 2, 4> denseConstraints;
    denseConstraints <<
        1, -1, 0, 0,
        0,  1, 1, -1;
    ZhangIarFunctional constraints = denseConstraints.sparseView();
    ZhangIarCovarianceCondition condition =
        zhangIarCovarianceCondition(covariance, constraints);
    BOOST_REQUIRE(condition.valid);
    BOOST_CHECK_EQUAL(condition.rank, 2);

    const Matrix2d constraintCovariance =
        denseConstraints * covariance * denseConstraints.transpose();
    const Matrix4d densePosterior = covariance -
        covariance * denseConstraints.transpose() *
        constraintCovariance.inverse() * denseConstraints * covariance;
    const Matrix4d factorPosterior = covariance -
        condition.reductionFactor * condition.reductionFactor.transpose();
    BOOST_CHECK_SMALL(
        (densePosterior - factorPosterior).norm(), 1e-11);

    Matrix<double, 2, 4> denseTarget;
    denseTarget <<
        1, 0, -1, 0,
        0, 1,  0, -1;
    ZhangIarFunctional target = denseTarget.sparseView();
    const double auditedTrace = zhangIarProjectedCovarianceTrace(
        covariance, condition, target);
    const double denseTrace =
        (denseTarget * densePosterior * denseTarget.transpose()).trace();
    BOOST_CHECK_SMALL(auditedTrace - denseTrace, 1e-11);

    Matrix2d unimodular;
    unimodular << 1, 1, 0, 1;
    Matrix<double, 2, 4> denseRebased =
        unimodular * denseConstraints;
    ZhangIarFunctional rebased = denseRebased.sparseView();
    ZhangIarCovarianceCondition rebasedCondition =
        zhangIarCovarianceCondition(covariance, rebased);
    BOOST_REQUIRE(rebasedCondition.valid);
    const double rebasedTrace = zhangIarProjectedCovarianceTrace(
        covariance, rebasedCondition, target);
    BOOST_CHECK_SMALL(rebasedTrace - auditedTrace, 1e-11);
}

BOOST_AUTO_TEST_CASE(theory_regression_physical_dd_row_closes_in_cycle_basis)
{
    const SatSys g01(E_Sys::GPS, 1);
    const SatSys g02(E_Sys::GPS, 2);
    const SatSys g03(E_Sys::GPS, 3);
    const std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R0", g03},
        {"R1", g01}, {"R1", g02}, {"R1", g03},
        {"R2", g01}, {"R2", g02}, {"R2", g03}
    };
    const ZhangGraphBasis basis = zhangBuildSpanningTree(edges, "R0");
    BOOST_REQUIRE(basis.connected);
    std::map<ZhangGraphEdge, int> chordColumns;
    for (const ZhangGraphEdge& edge : basis.edges)
    {
        if (!basis.isTreeEdge(edge.receiver, edge.satellite))
        {
            chordColumns[edge] = chordColumns.size();
        }
    }
    BOOST_REQUIRE_EQUAL(chordColumns.size(), 4);

    VectorXd row;
    BOOST_REQUIRE(zhangDdCycleCoordinateRow(
        basis, chordColumns, "R0", "R2", g01, g03, row));
    BOOST_CHECK_EQUAL(row.size(), chordColumns.size());
    BOOST_CHECK(row.allFinite());
    for (double value : row)
    {
        BOOST_CHECK_SMALL(value - std::round(value), 1e-14);
    }
}

BOOST_AUTO_TEST_CASE(theory_regression_detects_stronger_satellite_sd_correlation)
{
    Matrix3d covariance;
    covariance <<
         1.00, -0.30, 0.30,
        -0.30,  1.13, 0.91,
         0.30,  0.91, 1.13;
    Matrix<double, 1, 3> ambiguityDense;
    ambiguityDense << 1, 0, 0;
    Matrix<double, 1, 3> undifferencedDense;
    undifferencedDense << 0, 0, 1;
    Matrix<double, 1, 3> satelliteDifferenceDense;
    satelliteDifferenceDense << 0, -1, 1;
    const ZhangIarFunctional ambiguity = ambiguityDense.sparseView();
    const ZhangIarFunctional undifferenced =
        undifferencedDense.sparseView();
    const ZhangIarFunctional satelliteDifference =
        satelliteDifferenceDense.sparseView();

    const ZhangPairedCorrelationSummary ud = zhangPairedCorrelations(
        covariance, ambiguity, undifferenced);
    const ZhangPairedCorrelationSummary sd = zhangPairedCorrelations(
        covariance, ambiguity, satelliteDifference);
    BOOST_REQUIRE(ud.valid);
    BOOST_REQUIRE(sd.valid);
    BOOST_REQUIRE_EQUAL(ud.pairs, 1);
    BOOST_REQUIRE_EQUAL(sd.pairs, 1);
    BOOST_CHECK_SMALL(ud.coefficients.front() -
        0.30 / std::sqrt(1.13), 1e-12);
    BOOST_CHECK_SMALL(sd.coefficients.front() -
        0.60 / std::sqrt(0.44), 1e-12);
    BOOST_CHECK_SMALL(
        ud.pooledCorrelation - ud.coefficients.front(), 1e-12);
    BOOST_CHECK_SMALL(
        sd.pooledCorrelation - sd.coefficients.front(), 1e-12);
    BOOST_CHECK_LT(ud.rmsAbsolute, sd.rmsAbsolute);
}

BOOST_AUTO_TEST_CASE(exact_surviving_lattice_eliminates_removed_arcs_without_rounding)
{
    ZhangExactMatrix kernel = zhangExactIntegerKernel({
        {1, 1, 0},
        {0, 1, 1},
    });
    BOOST_REQUIRE_EQUAL(kernel.size(), 1);
    BOOST_REQUIRE_EQUAL(kernel.front().size(), 3);
    BOOST_CHECK(
        zhangExactMatrixTimesColumn(
            {{1, 1, 0}, {0, 1, 1}},
            kernel.front()
        ) == ZhangExactVector({0, 0})
    );

    // Every input row touches the removed third arc.  Their exact integer
    // combination r1-r2 survives as n1-n2=-2; row 3 is redundant.
    ZhangExactSurvivingLattice surviving = zhangExactSurvivingLattice(
        {
            {1, 0, 1},
            {0, 1, 1},
            {1, 1, 2},
        },
        {5, 7, 12},
        {true, true, false}
    );
    BOOST_CHECK(surviving.consistent);
    BOOST_CHECK_EQUAL(surviving.touchedRows, 3);
    BOOST_REQUIRE_EQUAL(surviving.basis.size(), 1);
    BOOST_CHECK(surviving.basis.front() == ZhangExactVector({1, -1}));
    BOOST_CHECK(surviving.values.front() == -2);

    ZhangExactSurvivingLattice none = zhangExactSurvivingLattice(
        {{1, 1}},
        {3},
        {true, false}
    );
    BOOST_CHECK(none.consistent);
    BOOST_CHECK(none.basis.empty());
}

BOOST_AUTO_TEST_CASE(integer_cycle_fixing_feedback_updates_full_state_and_covariance)
{
    constexpr int stateCount     = 9;
    constexpr int ambiguityCount = 3;

    MatrixXd generator = MatrixXd::Random(stateCount, stateCount);
    MatrixXd covariance =
        generator * generator.transpose() +
        0.5 * MatrixXd::Identity(stateCount, stateCount);
    VectorXd state = VectorXd::LinSpaced(stateCount, -1.2, 2.4);

    MatrixXd selector = MatrixXd::Zero(ambiguityCount, stateCount);
    selector(0, 4) = 1;
    selector(1, 6) = 1;
    selector(2, 8) = 1;

    VectorXd floatAmbiguities = selector * state;
    VectorXd fixedAmbiguities = floatAmbiguities.array().round().matrix();
    MatrixXd ambiguityCovariance = selector * covariance * selector.transpose();
    MatrixXd crossCovariance = covariance * selector.transpose();

    VectorXd fixedState =
        state -
        crossCovariance * ambiguityCovariance.inverse() *
            (floatAmbiguities - fixedAmbiguities);
    MatrixXd fixedCovariance =
        covariance -
        crossCovariance * ambiguityCovariance.inverse() * crossCovariance.transpose();
    fixedCovariance = 0.5 * (fixedCovariance + fixedCovariance.transpose());

    BOOST_CHECK_SMALL((selector * fixedState - fixedAmbiguities).norm(), 1e-12);
    BOOST_CHECK_SMALL((selector * fixedCovariance * selector.transpose()).norm(), 1e-10);
    BOOST_CHECK_LE(fixedCovariance.diagonal().sum(), covariance.diagonal().sum());

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(fixedCovariance);
    BOOST_REQUIRE_EQUAL(eigenSolver.info(), Eigen::Success);
    BOOST_CHECK_GE(eigenSolver.eigenvalues().minCoeff(), -1e-10);

    // Satellite clock/phase products occupy non-ambiguity states and must receive the correlated
    // conditional update rather than leaving the float products unchanged.
    BOOST_CHECK_GT((fixedState.head(4) - state.head(4)).norm(), 1e-6);
}

BOOST_AUTO_TEST_CASE(leave_one_out_internal_phase_products_restore_user_integer_differences)
{
    constexpr double wavelength = 0.190293672798365;

    std::set<ZhangGraphEdge> networkEdges = {
        {"R0", SatSys(E_Sys::GPS, 1)},
        {"R0", SatSys(E_Sys::GPS, 2)},
        {"R1", SatSys(E_Sys::GPS, 1)},
        {"R1", SatSys(E_Sys::GPS, 3)},
        {"R1", SatSys(E_Sys::GPS, 4)},
        {"R2", SatSys(E_Sys::GPS, 2)},
        {"R2", SatSys(E_Sys::GPS, 3)},
        {"R2", SatSys(E_Sys::GPS, 4)}
    };
    ZhangGraphBasis basis = zhangBuildSpanningTree(networkEdges, "R0");
    BOOST_REQUIRE(basis.connected);

    std::map<std::string, double> receiverBias = {
        {"R0", 0.13}, {"R1", -0.27}, {"R2", 0.41}
    };
    std::map<SatSys, double> satelliteBias;
    for (const auto& satellite : basis.satellites)
    {
        satelliteBias[satellite] = 0.07 * satellite.prn - 0.19;
    }

    VectorXd networkObservations(networkEdges.size());
    int edgeIndex = 0;
    for (const auto& edge : networkEdges)
    {
        int integerAmbiguity = 11 * (edgeIndex + 1) - 17;
        networkObservations(edgeIndex) =
            receiverBias.at(edge.receiver) +
            satelliteBias.at(edge.satellite) +
            wavelength * integerAmbiguity;
        edgeIndex++;
    }

    MatrixXd networkDesign = zhangGraphPhaseDesign(basis, wavelength);
    VectorXd networkState = networkDesign.fullPivLu().solve(networkObservations);
    BOOST_CHECK_SMALL(
        (networkDesign * networkState - networkObservations).norm(),
        1e-12
    );

    const int receiverPhaseCount = basis.receivers.size() - 1;
    std::vector<SatSys> satellites(basis.satellites.begin(), basis.satellites.end());
    std::map<SatSys, double> internalSatelliteProduct;
    for (int satellite = 0; satellite < satellites.size(); satellite++)
    {
        internalSatelliteProduct[satellites[satellite]] =
            networkState(receiverPhaseCount + satellite);
    }

    const double userBias = -0.33;
    std::map<SatSys, double> correctedUserPhase;
    for (int satellite = 0; satellite < satellites.size(); satellite++)
    {
        int userAmbiguity = 23 + 7 * satellite;
        double rawUserPhase =
            userBias +
            satelliteBias.at(satellites[satellite]) +
            wavelength * userAmbiguity;
        correctedUserPhase[satellites[satellite]] =
            rawUserPhase - internalSatelliteProduct.at(satellites[satellite]);
    }

    SatSys referenceSatellite = satellites.front();
    for (const auto& satellite : satellites)
    {
        double userIntegerDifference =
            (correctedUserPhase.at(satellite) -
             correctedUserPhase.at(referenceSatellite)) /
            wavelength;
        BOOST_CHECK_SMALL(
            userIntegerDifference - std::round(userIntegerDifference),
            1e-11
        );
    }
}

BOOST_AUTO_TEST_CASE(phase_continuity_integer_branch_change_preserves_validity)
{
    ZhangPhaseContinuityState state;
    state.markFixed();

    GTime time;
    time.bigTime = 1000;
    auto change = state.applyExactTransform(time, 4.0, 2);

    BOOST_CHECK(change == ZhangPhaseContinuityChange::EXACT_INTEGER);
    BOOST_CHECK_EQUAL(state.integerShiftCycles, 4);
    BOOST_CHECK_EQUAL(state.counter, 0);
    BOOST_CHECK_EQUAL(state.datumVersion, 0);
    BOOST_CHECK(state.integerValid());
}

BOOST_AUTO_TEST_CASE(phase_continuity_fractional_change_forces_user_reinitialisation)
{
    ZhangPhaseContinuityState state;
    state.markFixed();

    GTime time;
    time.bigTime = 2000;
    auto change = state.applyExactTransform(time, -2.25, 2);

    BOOST_CHECK(change == ZhangPhaseContinuityChange::EXACT_FRACTIONAL);
    BOOST_CHECK_EQUAL(state.counter, 1);
    BOOST_CHECK_EQUAL(state.datumVersion, 1);
    BOOST_CHECK_EQUAL(state.iod, 1);
    BOOST_CHECK_CLOSE(state.fractionalShiftCycles, -2.25, 1e-12);
    BOOST_CHECK(!state.integerValid());
    BOOST_CHECK_EQUAL(state.stabilizationRemaining, 2);
}

BOOST_AUTO_TEST_CASE(auxiliary_product_tree_generation_is_not_a_global_product_reset)
{
    ZhangProductDatumVersionTracker tracker;
    BOOST_CHECK(!tracker.observe(17));
    BOOST_CHECK(!tracker.observe(17));
    BOOST_CHECK(tracker.observe(18));
    BOOST_CHECK(!tracker.observe(18));

    ZhangPhaseContinuityState continuity;
    continuity.markFixed();
    BOOST_CHECK(tracker.observe(19));
    BOOST_CHECK_EQUAL(continuity.counter, 0);
    BOOST_CHECK_EQUAL(continuity.datumVersion, 0);
    BOOST_CHECK_EQUAL(continuity.iod, 0);
    BOOST_CHECK(continuity.integerValid());
}

BOOST_AUTO_TEST_CASE(hou_product_coordinate_absorbs_fractional_tree_transform)
{
    constexpr double wavelength = 0.190293672798365;
    ZhangPhaseContinuityState continuity;
    continuity.markFixed();

    Vector2d before;
    before << 12.4, -3.7; // satellite clock and internal phase node [m]
    auto beforeTarget = zhangHouOsbLikePhaseCorrectionTarget(
        2, 0, 1, wavelength, 0.0);
    const double referenceCorrection = beforeTarget.value(before);

    constexpr double treeOffsetCycles = -2.25;
    Vector2d after = before;
    after(1) -= treeOffsetCycles * wavelength;
    continuity.applyHouProductTransform(treeOffsetCycles);
    const double alignmentCycles =
        continuity.integerShiftCycles + continuity.fractionalShiftCycles;
    auto afterTarget = zhangHouOsbLikePhaseCorrectionTarget(
        2, 0, 1, wavelength, alignmentCycles);

    BOOST_CHECK_SMALL(
        afterTarget.value(after) - referenceCorrection,
        1e-12
    );
    BOOST_CHECK_EQUAL(continuity.counter, 0);
    BOOST_CHECK_EQUAL(continuity.datumVersion, 0);
    BOOST_CHECK_EQUAL(continuity.iod, 0);
    BOOST_CHECK(continuity.integerValid());
    BOOST_CHECK_EQUAL(
        continuity.resetReason,
        "hou_exact_affine_s_transform"
    );
}

BOOST_AUTO_TEST_CASE(phase_continuity_reinitialisation_resets_branch_and_stabilises)
{
    ZhangPhaseContinuityState state;
    state.integerShiftCycles = 7;
    state.fractionalShiftCycles = 0.3;
    state.markFixed();

    GTime resetTime;
    resetTime.bigTime = 3000;
    auto change = state.reinitialise(resetTime, "tree_state_missing", 2);

    BOOST_CHECK(change == ZhangPhaseContinuityChange::REINITIALISED);
    BOOST_CHECK_EQUAL(state.counter, 1);
    BOOST_CHECK_EQUAL(state.datumVersion, 1);
    BOOST_CHECK_EQUAL(state.integerShiftCycles, 0);
    BOOST_CHECK_SMALL(state.fractionalShiftCycles, 1e-15);
    BOOST_CHECK_EQUAL(state.resetReason, "tree_state_missing");
    BOOST_CHECK(!state.integerValid());

    state.markFixed();
    state.advanceEpoch(resetTime);
    BOOST_CHECK(!state.integerValid());
    GTime next = resetTime + 300;
    state.advanceEpoch(next);
    BOOST_CHECK(state.integerValid());
}

BOOST_AUTO_TEST_CASE(phase_continuity_lost_integer_rank_invalidates_once)
{
    ZhangPhaseContinuityState state;
    state.markFixed();
    BOOST_CHECK(state.integerValid());

    GTime changeTime;
    changeTime.bigTime = 1247461200;
    BOOST_CHECK(state.invalidateIntegerDatum(
        changeTime,
        "integer_datum_incomplete",
        2
    ));
    BOOST_CHECK(!state.integerValid());
    BOOST_CHECK_EQUAL(state.counter, 1);
    BOOST_CHECK_EQUAL(state.datumVersion, 1);
    BOOST_CHECK_EQUAL(state.iod, 1);
    BOOST_CHECK_EQUAL(state.resetReason, "integer_datum_incomplete");

    BOOST_CHECK(!state.invalidateIntegerDatum(
        changeTime,
        "integer_datum_incomplete",
        2
    ));
    BOOST_CHECK_EQUAL(state.counter, 1);
}

BOOST_AUTO_TEST_CASE(held_out_user_reference_exchange_preserves_phase_and_covariance)
{
    constexpr double lambda = 0.190293672798365;

    // Old coordinates: receiver phase relative to satellite A, followed by
    // single-difference ambiguities B-A and C-A.
    Vector3d oldState(0.27, 13.4, -7.2);
    Matrix3d generator;
    generator << 0.7, -0.2, 0.1,
                 0.3,  1.1, 0.4,
                -0.1,  0.5, 0.9;
    Matrix3d oldCovariance =
        generator * generator.transpose() +
        0.1 * Matrix3d::Identity();

    // New reference B: b_B=b_A+lambda*D_BA,
    // D_AB=-D_BA and D_CB=D_CA-D_BA.
    Matrix3d transform;
    transform << 1, lambda, 0,
                 0,     -1, 0,
                 0,     -1, 1;
    Vector3d newState = transform * oldState;
    Matrix3d newCovariance =
        transform * oldCovariance * transform.transpose();

    Vector3d oldPredictions(
        oldState(0),
        oldState(0) + lambda * oldState(1),
        oldState(0) + lambda * oldState(2)
    );
    Vector3d newPredictions(
        newState(0) + lambda * newState(1),
        newState(0),
        newState(0) + lambda * newState(2)
    );
    BOOST_CHECK_SMALL((oldPredictions - newPredictions).norm(), 1e-12);

    Matrix3d inverse = transform.inverse();
    BOOST_CHECK_SMALL(
        (inverse * newState - oldState).norm(),
        1e-12
    );
    BOOST_CHECK_SMALL(
        (
            inverse * newCovariance * inverse.transpose() -
            oldCovariance
        ).norm(),
        1e-12
    );
}

namespace
{
struct ZhangProjectedPhysicalTarget
{
    double mean = 0;
    double variance = 0;
    ZhangAffineUserTarget coordinateTarget;
};

VectorXd satelliteMeanContrast(
    const std::set<ZhangGraphEdge>& edges,
    const SatSys& positive,
    const SatSys& negative)
{
    int positiveCount = 0;
    int negativeCount = 0;
    for (const auto& edge : edges)
    {
        positiveCount += edge.satellite == positive;
        negativeCount += edge.satellite == negative;
    }
    VectorXd contrast = VectorXd::Zero(edges.size());
    if (positiveCount == 0 || negativeCount == 0)
    {
        return VectorXd();
    }
    int row = 0;
    for (const auto& edge : edges)
    {
        if (edge.satellite == positive)
        {
            contrast(row) = +1.0 / positiveCount;
        }
        if (edge.satellite == negative)
        {
            contrast(row) = -1.0 / negativeCount;
        }
        row++;
    }
    return contrast;
}

ZhangProjectedPhysicalTarget projectPhysicalTarget(
    const ZhangGraphBasis& basis,
    const VectorXd& physicalEdgeState,
    const MatrixXd& physicalEdgeCovariance,
    const VectorXd& physicalTargetRow,
    double wavelength = 0.190293672798365)
{
    MatrixXd design = zhangGraphPhaseDesign(basis, wavelength);
    MatrixXd inverse = design.inverse();
    VectorXd coordinateState = inverse * physicalEdgeState;
    MatrixXd coordinateCovariance =
        inverse * physicalEdgeCovariance * inverse.transpose();

    ZhangProjectedPhysicalTarget result;
    result.coordinateTarget.row = design.transpose() * physicalTargetRow;
    result.coordinateTarget.offset = 0;
    result.coordinateTarget.units = "cycle";
    result.mean = result.coordinateTarget.value(coordinateState);
    result.variance =
        result.coordinateTarget.variance(coordinateCovariance);
    return result;
}

double testScalarRoundErrorProbability(double fractional, double variance)
{
	if (variance < 1e-20)
	{
		return 0;
	}
	double alternateMass = 0;
	const double exponentScale = -0.25 / variance;
	for (int offset = 1; offset < 10; offset++)
	{
		alternateMass += std::exp(
			(offset + 2 * fractional) * offset * exponentScale);
		alternateMass += std::exp(
			(offset - 2 * fractional) * offset * exponentScale);
	}
	return alternateMass / (1 + alternateMass);
}

void checkProjectedTargetInvariant(
    const ZhangProjectedPhysicalTarget& first,
    const ZhangProjectedPhysicalTarget& second)
{
    BOOST_CHECK_SMALL(first.mean - second.mean, 1e-10);
    BOOST_CHECK_LT(
        zhangProtectedRelativeVarianceDifference(
            first.variance, second.variance
        ),
        1e-10
    );
	BOOST_REQUIRE_GT(first.variance, 0);
	BOOST_REQUIRE_GT(second.variance, 0);
	const long long firstCandidate = std::llround(first.mean);
	const long long secondCandidate = std::llround(second.mean);
	BOOST_CHECK_EQUAL(firstCandidate, secondCandidate);
	const double firstFractional = first.mean - firstCandidate;
	const double secondFractional = second.mean - secondCandidate;
	BOOST_CHECK_SMALL(
		testScalarRoundErrorProbability(firstFractional, first.variance)
		- testScalarRoundErrorProbability(secondFractional, second.variance),
		1e-10);
	BOOST_CHECK_SMALL(
		firstFractional * firstFractional / first.variance
		- secondFractional * secondFractional / second.variance,
		1e-10);
}
}

BOOST_AUTO_TEST_CASE(user_phase_and_wl_targets_survive_three_affine_s_bases)
{
    constexpr double lambda1 = 0.190293672798365;
    constexpr double lambda2 = 0.244210213424568;
    VectorXd physicalState(6);
    // C_s, C_r, B_s1, B_r1, B_s2, B_r2, all in metres.
    physicalState << 0.18, -0.07, 1.22, -0.41, 0.87, -1.31;
    MatrixXd generator(6, 6);
    generator <<
        1.0,  0.1,  0.0,  0.2, -0.1,  0.0,
        0.2,  0.9, -0.1,  0.0,  0.1,  0.2,
        0.1, -0.2,  1.1,  0.3,  0.0, -0.1,
        0.0,  0.1,  0.2,  0.8, -0.2,  0.0,
       -0.1,  0.0,  0.1, -0.1,  1.2,  0.2,
        0.2, -0.1,  0.0,  0.1,  0.3,  0.9;
    MatrixXd physicalCovariance =
        generator * generator.transpose()
        + 0.2 * MatrixXd::Identity(6, 6);

    auto s1 = zhangUserPhaseCorrectionTarget(6, 0, 2, lambda1, 5);
	auto houS1 = zhangHouOsbLikePhaseCorrectionTarget(
		6, 0, 2, lambda1, 5);
    auto r1 = zhangUserPhaseCorrectionTarget(6, 1, 3, lambda1, -2);
    auto s2 = zhangUserPhaseCorrectionTarget(6, 0, 4, lambda2, 3);
    auto r2 = zhangUserPhaseCorrectionTarget(6, 1, 5, lambda2, 1);
    auto l1Difference = zhangLinearCombination(s1, +1, r1, -1, "metre");
    auto l2Difference = zhangLinearCombination(s2, +1, r2, -1, "metre");
    auto wideLane = zhangLinearCombination(
        l1Difference, +1 / lambda1,
        l2Difference, -1 / lambda2,
        "cycle"
    );
    BOOST_CHECK_EQUAL(
        zhangUserPhaseCorrectionValue(
            physicalState(0), physicalState(2), lambda1, 5
        ),
        physicalState(0) - (physicalState(2) + 5 * lambda1)
    );
	BOOST_CHECK_SMALL((houS1.row - s1.row).norm(), 1e-15);
	BOOST_CHECK_SMALL(houS1.offset - s1.offset, 1e-15);
	VectorXd commonPhaseDatumState = physicalState;
	commonPhaseDatumState(0) += 37.25;
	commonPhaseDatumState(2) += 37.25;
	BOOST_CHECK_SMALL(
		houS1.value(commonPhaseDatumState) - houS1.value(physicalState),
		1e-12);
    const double referenceMean = wideLane.value(physicalState);
    const double referenceVariance = wideLane.variance(physicalCovariance);

    std::vector<MatrixXd> transforms;
    MatrixXd first = MatrixXd::Identity(6, 6);
    first(2, 0) = 1;
    first(3, 1) = -1;
    first(4, 2) = 1;
    transforms.push_back(first);
    MatrixXd second = MatrixXd::Identity(6, 6);
    second.row(0).swap(second.row(1));
    second(4, 0) = -2;
    second(5, 1) = 1;
    transforms.push_back(second);
    MatrixXd third = MatrixXd::Identity(6, 6);
    third(0, 2) = 1;
    third(1, 3) = 1;
    third(4, 5) = -1;
    transforms.push_back(third);

    int transformNumber = 1;
    for (const auto& transform : transforms)
    {
        VectorXd translation = VectorXd::LinSpaced(
            6, -0.03 * transformNumber, 0.02 * transformNumber
        );
        VectorXd transformedState = transform * physicalState + translation;
        MatrixXd transformedCovariance =
            transform * physicalCovariance * transform.transpose();
        ZhangAffineUserTarget transformedTarget;
        BOOST_REQUIRE(zhangTransportAffineUserTarget(
            wideLane, transform, translation, transformedTarget
        ));
        BOOST_CHECK_SMALL(
            transformedTarget.value(transformedState) - referenceMean,
            1e-10
        );
        BOOST_CHECK_LT(
            zhangProtectedRelativeVarianceDifference(
                transformedTarget.variance(transformedCovariance),
                referenceVariance
            ),
            1e-10
        );
        transformNumber++;
    }
}

BOOST_AUTO_TEST_CASE(
    user_target_is_invariant_for_tree_receiver_and_satellite_reference_changes)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R0", g03},
        {"R1", g01}, {"R1", g02}, {"R1", g03},
        {"R2", g01}, {"R2", g02}, {"R2", g03}
    };
    ZhangGraphBasis t1 = zhangBuildSpanningTree(
        edges, "R0",
        {{"R0", g01}, {"R0", g02}, {"R0", g03},
         {"R1", g01}, {"R2", g02}}
    );
    ZhangGraphBasis t2 = zhangBuildSpanningTree(
        edges, "R0",
        {{"R0", g01}, {"R1", g01}, {"R1", g03},
         {"R2", g02}, {"R2", g03}}
    );
    ZhangGraphBasis t3 = zhangBuildSpanningTree(
        edges, "R1",
        {{"R1", g01}, {"R1", g02}, {"R2", g02},
         {"R2", g03}, {"R0", g03}}
    );
    BOOST_REQUIRE(t1.connected && t2.connected && t3.connected);
    BOOST_REQUIRE(t1.treeEdges != t2.treeEdges);
    BOOST_CHECK_NE(t1.rootReceiver, t3.rootReceiver);

    VectorXd physicalState = VectorXd::LinSpaced(edges.size(), -2.1, 1.7);
    MatrixXd generator = MatrixXd::Zero(edges.size(), edges.size());
    for (int row = 0; row < generator.rows(); row++)
    for (int column = 0; column <= row; column++)
    {
        generator(row, column) =
            row == column ? 1.0 + 0.1 * row : 0.01 * (row + column + 1);
    }
    MatrixXd physicalCovariance =
        generator * generator.transpose()
        + 0.1 * MatrixXd::Identity(edges.size(), edges.size());
    VectorXd g03MinusG01 = satelliteMeanContrast(edges, g03, g01);
    VectorXd g03MinusG02 = satelliteMeanContrast(edges, g03, g02);
    VectorXd g02MinusG01 = satelliteMeanContrast(edges, g02, g01);
    BOOST_CHECK_SMALL(
        (g03MinusG01 - g03MinusG02 - g02MinusG01).norm(),
        1e-15
    );

    auto p1 = projectPhysicalTarget(
        t1, physicalState, physicalCovariance, g03MinusG01
    );
    auto p2 = projectPhysicalTarget(
        t2, physicalState, physicalCovariance, g03MinusG01
    );
    auto p3 = projectPhysicalTarget(
        t3, physicalState, physicalCovariance, g03MinusG01
    );
    checkProjectedTargetInvariant(p1, p2);
    checkProjectedTargetInvariant(p1, p3);

    // Re-express the same G03-G01 product through reference satellite G02.
    auto p32 = projectPhysicalTarget(
        t3, physicalState, physicalCovariance, g03MinusG02
    );
    auto p21 = projectPhysicalTarget(
        t3, physicalState, physicalCovariance, g02MinusG01
    );
    BOOST_CHECK_SMALL(p3.mean - p32.mean - p21.mean, 1e-10);
    BOOST_CHECK_SMALL(
        (p3.coordinateTarget.row
            - p32.coordinateTarget.row
            - p21.coordinateTarget.row).norm(),
        1e-10
    );
}

BOOST_AUTO_TEST_CASE(
    satellite_join_and_leaf_exit_preserve_common_physical_target_subspace)
{
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    std::set<ZhangGraphEdge> oldEdges = {
        {"R0", g01}, {"R0", g02},
        {"R1", g01}, {"R1", g02}
    };
    std::set<ZhangGraphEdge> joinedEdges = oldEdges;
    joinedEdges.insert({"R0", g03});
    joinedEdges.insert({"R1", g03});
    ZhangGraphBasis oldBasis = zhangBuildSpanningTree(oldEdges, "R0");
    ZhangGraphBasis joinedBasis = zhangBuildSpanningTree(
        joinedEdges, "R1",
        {{"R1", g01}, {"R1", g02}, {"R1", g03},
         {"R0", g02}}
    );
    BOOST_REQUIRE(oldBasis.connected && joinedBasis.connected);

    std::map<ZhangGraphEdge, double> value;
    std::map<ZhangGraphEdge, double> variance;
    int index = 1;
    for (const auto& edge : joinedEdges)
    {
        value[edge] = -0.8 + 0.17 * index;
        variance[edge] = 0.05 + 0.01 * index;
        index++;
    }
    auto assembleState = [&](const std::set<ZhangGraphEdge>& selected)
    {
        VectorXd state(selected.size());
        int row = 0;
        for (const auto& edge : selected)
        {
            state(row++) = value.at(edge);
        }
        return state;
    };
    auto assembleCovariance = [&](const std::set<ZhangGraphEdge>& selected)
    {
        MatrixXd covariance = MatrixXd::Zero(selected.size(), selected.size());
        int row = 0;
        for (const auto& edge : selected)
        {
            covariance(row, row) = variance.at(edge);
            row++;
        }
        return covariance;
    };

    auto beforeJoin = projectPhysicalTarget(
        oldBasis,
        assembleState(oldEdges),
        assembleCovariance(oldEdges),
        satelliteMeanContrast(oldEdges, g02, g01)
    );
    auto afterJoin = projectPhysicalTarget(
        joinedBasis,
        assembleState(joinedEdges),
        assembleCovariance(joinedEdges),
        satelliteMeanContrast(joinedEdges, g02, g01)
    );
    checkProjectedTargetInvariant(beforeJoin, afterJoin);

    // The reverse comparison is the leaf-exit case.  G03 targets are retired,
    // while the common G02-G01 target remains unchanged.
    checkProjectedTargetInvariant(afterJoin, beforeJoin);
}

BOOST_AUTO_TEST_CASE(dual_frequency_common_target_survives_independent_backbone_changes)
{
    constexpr double lambda1 = 0.190293672798365;
    constexpr double lambda2 = 0.244210213424568;
    SatSys g01(E_Sys::GPS, 1);
    SatSys g02(E_Sys::GPS, 2);
    SatSys g03(E_Sys::GPS, 3);
    std::set<ZhangGraphEdge> edges = {
        {"R0", g01}, {"R0", g02}, {"R0", g03},
        {"R1", g01}, {"R1", g02}, {"R1", g03},
        {"R2", g01}, {"R2", g02}, {"R2", g03}
    };
    ZhangGraphBasis a = zhangBuildSpanningTree(edges, "R0");
    ZhangGraphBasis b = zhangBuildSpanningTree(
        edges, "R1",
        {{"R1", g01}, {"R1", g02}, {"R1", g03},
         {"R0", g01}, {"R2", g03}}
    );
    ZhangGraphBasis c = zhangBuildSpanningTree(
        edges, "R2",
        {{"R2", g01}, {"R2", g02}, {"R2", g03},
         {"R0", g02}, {"R1", g01}}
    );
    BOOST_REQUIRE(a.connected && b.connected && c.connected);

    const int edgeCount = edges.size();
    VectorXd physicalState = VectorXd::LinSpaced(2 * edgeCount, -1.9, 2.4);
    MatrixXd generator = MatrixXd::Identity(2 * edgeCount, 2 * edgeCount);
    for (int row = 1; row < generator.rows(); row++)
    {
        generator(row, row - 1) = 0.13;
    }
    MatrixXd physicalCovariance =
        generator * generator.transpose()
        + 0.1 * MatrixXd::Identity(2 * edgeCount, 2 * edgeCount);
    VectorXd satelliteContrast = satelliteMeanContrast(edges, g03, g01);
    VectorXd wideLanePhysicalRow(2 * edgeCount);
    wideLanePhysicalRow.head(edgeCount) = satelliteContrast / lambda1;
    wideLanePhysicalRow.tail(edgeCount) = -satelliteContrast / lambda2;

    auto projectDual = [&](const ZhangGraphBasis& first,
                           const ZhangGraphBasis& second)
    {
        MatrixXd design = MatrixXd::Zero(2 * edgeCount, 2 * edgeCount);
        design.topLeftCorner(edgeCount, edgeCount) =
            zhangGraphPhaseDesign(first, lambda1);
        design.bottomRightCorner(edgeCount, edgeCount) =
            zhangGraphPhaseDesign(second, lambda2);
        MatrixXd inverse = design.inverse();
        VectorXd state = inverse * physicalState;
        MatrixXd covariance =
            inverse * physicalCovariance * inverse.transpose();
        ZhangProjectedPhysicalTarget result;
        result.coordinateTarget.row = design.transpose() * wideLanePhysicalRow;
        result.coordinateTarget.units = "cycle";
        result.mean = result.coordinateTarget.value(state);
        result.variance = result.coordinateTarget.variance(covariance);
        return result;
    };

    auto first = projectDual(a, b);
    auto second = projectDual(c, a);
    auto third = projectDual(b, c);
    checkProjectedTargetInvariant(first, second);
    checkProjectedTargetInvariant(first, third);
}

BOOST_AUTO_TEST_CASE(fixed_lag_identity_resets_only_for_real_physical_change)
{
    using Transition = ZhangFixedLagIdentityTransition;
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(true, true, false, false, false, false)
        == Transition::CONTINUE
    );
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(true, true, false, false, true, true)
        == Transition::CONTINUE_EXACT_TRANSFORM
    );
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(true, true, false, false, true, false)
        == Transition::RESET_EXACT_TRANSFORM_UNAVAILABLE
    );
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(true, true, true, false, true, true)
        == Transition::RESET_PHYSICAL_IDENTITY
    );
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(true, true, false, true, true, true)
        == Transition::RESET_PHYSICAL_IDENTITY
    );
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(false, true, false, false, true, true)
        == Transition::START_NEW_TARGET
    );
    BOOST_CHECK(
        zhangClassifyFixedLagTransition(true, false, false, false, true, true)
        == Transition::RETIRE_TARGET
    );
}

BOOST_AUTO_TEST_CASE(square_root_window_matches_dense_schur_marginal)
{
    // Columns are [two epoch-local nuisance, two retained physical targets].
    MatrixXd denseFactor(8, 4);
    denseFactor <<
        1.0,  0.0,  0.4, -0.1,
        0.0,  1.0, -0.2,  0.3,
        0.7, -0.1,  1.0,  0.0,
       -0.3,  0.8,  0.0,  1.0,
        0.5,  0.2,  0.3,  0.7,
       -0.2,  0.4,  0.8, -0.5,
        0.1, -0.6,  0.2,  0.9,
        0.9,  0.3, -0.4,  0.2;
    VectorXd rhs(8);
    rhs << 0.7, -0.2, 1.1, -0.4, 0.5, 0.9, -0.7, 0.3;
    SparseMatrix<double> sparseFactor = denseFactor.sparseView();
    ZhangSquareRootMarginal marginal =
        zhangMarginaliseSquareRootFactors(sparseFactor, rhs, 2);
    BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
    BOOST_CHECK_EQUAL(marginal.nuisanceRank, 2);
    BOOST_CHECK_EQUAL(marginal.targetRank, 2);

    MatrixXd normal = denseFactor.transpose() * denseFactor;
    VectorXd natural = denseFactor.transpose() * rhs;
    Matrix2d Nnn = normal.topLeftCorner<2, 2>();
    Matrix2d Nnt = normal.topRightCorner<2, 2>();
    Matrix2d Ntn = normal.bottomLeftCorner<2, 2>();
    Matrix2d Ntt = normal.bottomRightCorner<2, 2>();
    Vector2d hn = natural.head<2>();
    Vector2d ht = natural.tail<2>();
    Matrix2d denseSchur = Ntt - Ntn * Nnn.inverse() * Nnt;
    Vector2d denseNatural = ht - Ntn * Nnn.inverse() * hn;
    Vector2d denseMean = denseSchur.inverse() * denseNatural;
    Matrix2d denseCovariance = denseSchur.inverse();
    BOOST_CHECK_SMALL((marginal.mean - denseMean).norm(), 1e-11);
    BOOST_CHECK_SMALL(
        (marginal.covariance - denseCovariance).norm()
            / denseCovariance.norm(),
        1e-11
    );
}

BOOST_AUTO_TEST_CASE(square_root_window_is_invariant_to_retained_coordinate_change)
{
    MatrixXd originalFactor(7, 4);
    originalFactor <<
        1.0,  0.0,  0.2, -0.3,
        0.0,  1.0,  0.4,  0.1,
        0.6, -0.2,  1.0,  0.0,
        0.1,  0.7,  0.0,  1.0,
       -0.4,  0.3,  0.5,  0.8,
        0.8,  0.1, -0.2,  0.6,
        0.2, -0.5,  0.7, -0.4;
    VectorXd originalRhs(7);
    originalRhs << 0.3, -0.6, 1.2, 0.4, -0.1, 0.8, -0.5;
    auto original = zhangMarginaliseSquareRootFactors(
        originalFactor.sparseView(), originalRhs, 2
    );
    BOOST_REQUIRE_MESSAGE(original.valid, original.failureReason);

    // z_new = T z_old + b.  Substitute z_old=T^-1(z_new-b) in every factor.
    Matrix2d transform;
    transform << 1, 1, 0, 1;
    Vector2d translation(3, -2);
    Matrix2d inverse = transform.inverse();
    MatrixXd changedFactor = originalFactor;
    changedFactor.rightCols(2) = originalFactor.rightCols(2) * inverse;
    VectorXd changedRhs = originalRhs
        + originalFactor.rightCols(2) * inverse * translation;
    auto changed = zhangMarginaliseSquareRootFactors(
        changedFactor.sparseView(), changedRhs, 2
    );
    BOOST_REQUIRE_MESSAGE(changed.valid, changed.failureReason);
    Vector2d expectedMean = transform * original.mean + translation;
    Matrix2d expectedCovariance =
        transform * original.covariance * transform.transpose();
    BOOST_CHECK_SMALL((changed.mean - expectedMean).norm(), 1e-11);
    BOOST_CHECK_SMALL(
        (changed.covariance - expectedCovariance).norm()
            / expectedCovariance.norm(),
        1e-11
    );
}

BOOST_AUTO_TEST_CASE(whitening_reveals_rank_and_rejects_negative_variance)
{
    Vector3d residual(0.3, -0.5, 0.2);
    Matrix3d semidefinite;
    semidefinite <<
        2.0, 0.4, 0.0,
        0.4, 1.0, 0.0,
        0.0, 0.0, 1e-16;
    ZhangWhitenedBlock whitened =
        zhangWhitenRetainedResidual(residual, semidefinite, 1e-12);
    BOOST_REQUIRE_MESSAGE(whitened.valid, whitened.failureReason);
    BOOST_CHECK_EQUAL(whitened.rank, 2);
    Matrix2d leading = semidefinite.topLeftCorner<2, 2>();
    Vector2d leadingResidual = residual.head<2>();
    BOOST_CHECK_SMALL(
        whitened.squaredNorm
            - leadingResidual.dot(leading.ldlt().solve(leadingResidual)),
        1e-12
    );

    Matrix3d invalid = semidefinite;
    invalid(2, 2) = -1e-3;
    auto rejected = zhangWhitenRetainedResidual(residual, invalid, 1e-12);
    BOOST_CHECK(!rejected.valid);
    BOOST_CHECK_EQUAL(rejected.failureReason, "NEGATIVE_COVARIANCE_DIRECTION");
}

BOOST_AUTO_TEST_CASE(final_accepted_factor_capture_preserves_order_and_linearisation)
{
    KFKey clock;
    clock.type = KF::SAT_CLOCK;
    clock.Sat = SatSys(E_Sys::GPS, 1);
    KFKey phase;
    phase.type = KF::PHASE_BIAS;
    phase.Sat = SatSys(E_Sys::GPS, 1);
    phase.num = static_cast<int>(E_ObsCode::L1W);
    std::vector<ZhangCapturedStateKey> keys = {
        zhangCapturedStateKey(clock), zhangCapturedStateKey(phase)
    };

    Vector2d priorMean(0.4, -0.7);
    Matrix2d priorCovariance;
    priorCovariance << 0.8, 0.1, 0.1, 1.2;
    KFMeas firstMeasurement;
    firstMeasurement.time.bigTime = 1000;
    firstMeasurement.H = MatrixXd::Zero(2, 2);
    firstMeasurement.H << 1, -1, 0.5, 0.2;
    firstMeasurement.V = Vector2d(0.03, -0.04);
    firstMeasurement.R = Matrix2d::Zero();
    firstMeasurement.R << 0.01, 0.002, 0.002, 0.02;
    firstMeasurement.obsKeys = {clock, phase};
	firstMeasurement.prefitRatios = Vector2d(0.5, -2.0);

    ZhangFactorCaptureBuffer capture;
    capture.setMaximumEvents(10);
	const Matrix2d firstInnovationCovariance = firstMeasurement.H
		* priorCovariance * firstMeasurement.H.transpose()
		+ firstMeasurement.R;
	const Matrix2d firstGain = priorCovariance * firstMeasurement.H.transpose()
		* firstInnovationCovariance.inverse();
	const Vector2d firstPosteriorMean = priorMean
		+ firstGain * firstMeasurement.V;
	Matrix2d firstPosteriorCovariance = priorCovariance
		- firstGain * firstMeasurement.H * priorCovariance;
	firstPosteriorCovariance = 0.5
		* (firstPosteriorCovariance + firstPosteriorCovariance.transpose());
    BOOST_REQUIRE(capture.recordMeasurement(
        firstMeasurement.time,
        keys,
        priorMean,
        priorCovariance,
        firstMeasurement,
        "/PPP",
        firstPosteriorMean,
        firstPosteriorCovariance
    ));
    BOOST_REQUIRE_EQUAL(capture.capturedEvents().size(), 1);
    const auto& first = capture.capturedEvents().front();
    BOOST_CHECK_SMALL(
        (first.rightHandSide
            - (firstMeasurement.V + firstMeasurement.H * priorMean)).norm(),
        1e-15
    );
	Vector2d firstTargetRow(1, 0);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		firstMeasurement.time,
		"GPS:WL:G01:G02",
		"L1W:REC:G01:A0=1;L2W:REC:G01:A0=-1;",
		"G01:0:0->G02:0:0",
		{{"L1W:REC:G01", 0}, {"L2W:REC:G01", 0}},
		keys,
		firstTargetRow,
		2,
		firstPosteriorMean,
		firstPosteriorCovariance
	));

    SparseMatrix<double> transition(2, 2);
    transition.insert(0, 0) = 1;
    transition.insert(1, 1) = 1;
    Matrix2d processCovariance = Matrix2d::Zero();
    processCovariance(1, 1) = 0.04;
    GTime transitionTime;
    transitionTime.bigTime = 1030;
    BOOST_REQUIRE(capture.recordTransition(
        transitionTime,
        keys,
        keys,
        transition,
        processCovariance,
        "KF_STATE_TRANSITION"
    ));

    SparseMatrix<double> exactTransform(2, 2);
    exactTransform.insert(0, 0) = 1;
    exactTransform.insert(0, 1) = 1;
    exactTransform.insert(1, 1) = 1;
    BOOST_REQUIRE(capture.recordCoordinateTransform(
        transitionTime,
        keys,
        keys,
        exactTransform,
        "synthetic S-basis exchange"
    ));

    KFMeas secondMeasurement = firstMeasurement;
    secondMeasurement.time.bigTime = 1030;
    Matrix2d denseExactTransform = MatrixXd(exactTransform);
    Vector2d transitionedPrior = firstPosteriorMean;
    Matrix2d transitionedCovariance =
        firstPosteriorCovariance + processCovariance;
    Vector2d transformedPrior = denseExactTransform * transitionedPrior;
    Matrix2d transformedPriorCovariance =
        denseExactTransform * transitionedCovariance
            * denseExactTransform.transpose();
	const Matrix2d secondInnovationCovariance = secondMeasurement.H
		* transformedPriorCovariance * secondMeasurement.H.transpose()
		+ secondMeasurement.R;
	const Matrix2d secondGain = transformedPriorCovariance
		* secondMeasurement.H.transpose()
		* secondInnovationCovariance.inverse();
	const Vector2d secondPosteriorMean = transformedPrior
		+ secondGain * secondMeasurement.V;
	Matrix2d secondPosteriorCovariance = transformedPriorCovariance
		- secondGain * secondMeasurement.H * transformedPriorCovariance;
	secondPosteriorCovariance = 0.5
		* (secondPosteriorCovariance + secondPosteriorCovariance.transpose());
    BOOST_REQUIRE(capture.recordMeasurement(
        secondMeasurement.time,
        keys,
        transformedPrior,
        transformedPriorCovariance,
        secondMeasurement,
        "/PPP",
		secondPosteriorMean,
		secondPosteriorCovariance
    ));
	Vector2d transformedTargetRow(1, -1);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		secondMeasurement.time,
		"GPS:WL:G01:G02",
		"L1W:REC:G02:A0=1;L2W:REC:G02:A0=-1;",
		"G01:0:0->G02:0:0",
		{{"L1W:REC:G01", 0}, {"L2W:REC:G01", 0}},
		keys,
		transformedTargetRow,
		2,
		secondPosteriorMean,
		secondPosteriorCovariance
	));
	const auto rawIntegerDatum = capture.currentRawIntegerDatumMarginal();
	BOOST_REQUIRE_MESSAGE(rawIntegerDatum.valid, rawIntegerDatum.failureReason);
	BOOST_CHECK_EQUAL(rawIntegerDatum.targetRank, 1);
	BOOST_CHECK(!capture.recordPhysicalTarget(
		secondMeasurement.time,
		"GPS:WL:G01:G02",
		"L1W:REC:G02:A1=1;L2W:REC:G02:A0=-1;",
		"G01:0:0->G02:0:0",
		{{"L1W:REC:G01", 1}, {"L2W:REC:G01", 0}},
		keys,
		transformedTargetRow,
		2,
		secondPosteriorMean,
		secondPosteriorCovariance
	));
	BOOST_CHECK_EQUAL(
		capture.lastTargetReason(),
		"PERSISTENT_RAW_TARGET_PHYSICAL_VERSION_CHANGED");
    ZhangFactorCaptureSummary summary = capture.summary();
    BOOST_REQUIRE_MESSAGE(summary.valid, summary.failureReason);
    BOOST_CHECK_EQUAL(summary.events, 4);
    BOOST_CHECK_EQUAL(summary.measurements, 2);
    BOOST_CHECK_EQUAL(summary.transitions, 1);
    BOOST_CHECK_EQUAL(summary.coordinateTransforms, 1);
	BOOST_CHECK_EQUAL(summary.physicalTargets, 2);
	BOOST_CHECK_EQUAL(summary.physicalTargetCoordinateContinuations, 1);
	BOOST_CHECK_EQUAL(summary.physicalTargetIdentityResets, 0);
    BOOST_CHECK_EQUAL(summary.measurementRows, 4);
    BOOST_CHECK_SMALL(summary.maximumReplayPriorMeanRelativeError, 1e-15);
    BOOST_CHECK_SMALL(
        summary.maximumReplayPriorCovarianceRelativeError, 1e-15
    );
    BOOST_CHECK(
        capture.capturedEvents()[2].nonsingularCoordinateTransform
    );
	const auto scaleDiagnostics = capture.innovationScaleDiagnostics();
	BOOST_REQUIRE_EQUAL(scaleDiagnostics.size(), 2);
	std::vector<double> diagnosticScales;
	for (const auto& scale : scaleDiagnostics)
	{
		BOOST_CHECK_EQUAL(scale.blocks, 2);
		BOOST_CHECK_EQUAL(scale.samples, 2);
		diagnosticScales.push_back(scale.predictiveCovarianceScaleMle());
	}
	std::sort(diagnosticScales.begin(), diagnosticScales.end());
	BOOST_CHECK_CLOSE(diagnosticScales[0], 0.25, 1e-12);
	BOOST_CHECK_CLOSE(diagnosticScales[1], 4.0, 1e-12);
	BOOST_REQUIRE_EQUAL(capture.capturedPhysicalTargets().size(), 2);
	BOOST_CHECK_SMALL(
		capture.capturedPhysicalTargets()[1].mean
			- (2 + transformedTargetRow.dot(secondPosteriorMean)),
		1e-15
	);
	BOOST_CHECK_SMALL(
		capture.capturedPhysicalTargets()[1].variance
			- (transformedTargetRow.transpose()
				* secondPosteriorCovariance * transformedTargetRow)(0, 0),
		1e-15
	);
	BOOST_CHECK_SMALL(summary.maximumTargetMeanRelativeError, 1e-15);
	BOOST_CHECK_SMALL(summary.maximumTargetVarianceRelativeError, 1e-15);
}

BOOST_AUTO_TEST_CASE(factor_capture_fails_closed_on_state_key_chain_mismatch)
{
    KFKey firstKey;
    firstKey.type = KF::SAT_CLOCK;
    firstKey.Sat = SatSys(E_Sys::GPS, 1);
    KFKey secondKey = firstKey;
    secondKey.Sat = SatSys(E_Sys::GPS, 2);
    std::vector<ZhangCapturedStateKey> firstMap = {
        zhangCapturedStateKey(firstKey)
    };
    std::vector<ZhangCapturedStateKey> secondMap = {
        zhangCapturedStateKey(secondKey)
    };
    KFMeas measurement;
    measurement.H = MatrixXd::Ones(1, 1);
    measurement.V = VectorXd::Zero(1);
    measurement.R = MatrixXd::Identity(1, 1);
    measurement.obsKeys = {firstKey};
    VectorXd mean = VectorXd::Zero(1);
    MatrixXd covariance = MatrixXd::Identity(1, 1);

    ZhangFactorCaptureBuffer capture;
	MatrixXd posteriorCovariance = MatrixXd::Constant(1, 1, 0.5);
    BOOST_REQUIRE(capture.recordMeasurement(
        GTime(),
        firstMap,
        mean,
        covariance,
        measurement,
        "/PPP",
		mean,
		posteriorCovariance
    ));
    SparseMatrix<double> identity(1, 1);
    identity.insert(0, 0) = 1;
    BOOST_CHECK(!capture.recordTransition(
        GTime(), secondMap, secondMap, identity, covariance, "bad chain"
    ));
    BOOST_CHECK_EQUAL(
        capture.summary().failureReason,
        "INVALID_TRANSITION_CAPTURE_OR_KEY_CHAIN"
    );
}

BOOST_AUTO_TEST_CASE(authoritative_factor_commit_packet_preserves_affine_chain)
{
	KFKey state;
	state.type = KF::SAT_CLOCK;
	state.Sat = SatSys(E_Sys::GPS, 3);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(state)};
	const VectorXd priorMean = VectorXd::Constant(1, 2.0);
	const MatrixXd priorCovariance = MatrixXd::Constant(1, 1, 4.0);
	KFMeas measurement;
	measurement.H = MatrixXd::Zero(1, 1);
	measurement.V = VectorXd::Zero(1);
	measurement.R = MatrixXd::Identity(1, 1);
	measurement.obsKeys = {state};

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		priorMean, priorCovariance, 40, 41));
	BOOST_REQUIRE(capture.bindPersistentSnapshot(
		"X", "X", VectorXd::Ones(1), 0));

	SparseMatrix<double> transition(1, 1);
	transition.insert(0, 0) = 1.5;
	const MatrixXd processCovariance = MatrixXd::Constant(1, 1, 0.25);
	const VectorXd transitionOffset = VectorXd::Constant(1, -0.75);
	const VectorXd transitionedMean = transition * priorMean + transitionOffset;
	const MatrixXd transitionedCovariance =
		transition * priorCovariance * transition.transpose()
		+ processCovariance;
	BOOST_REQUIRE(capture.recordTransition(
		GTime(), keys, keys, transition, processCovariance,
		transitionOffset, "AFFINE_STATE_TRANSITION", transitionedMean,
		transitionedCovariance, 41, 42));

	SparseMatrix<double> transform(1, 1);
	transform.insert(0, 0) = -2.0;
	const VectorXd transformOffset = VectorXd::Constant(1, 1.25);
	const VectorXd transformedMean = transform * transitionedMean
		+ transformOffset;
	const MatrixXd transformedCovariance =
		transform * transitionedCovariance * transform.transpose();
	BOOST_REQUIRE(capture.recordCoordinateTransform(
		GTime(), keys, keys, transform, "AFFINE_EXACT_TRANSFORM", false,
		transformOffset, transformedMean, transformedCovariance, 42, 43));

	const auto summary = capture.summary();
	BOOST_REQUIRE_MESSAGE(summary.valid, summary.failureReason);
	BOOST_CHECK_EQUAL(summary.lastCommitSequence, 43);
	BOOST_CHECK_EQUAL(summary.uniqueMeasurementFactorIds, 1);
	BOOST_CHECK_EQUAL(summary.duplicateMeasurementFactorIds, 0);
	BOOST_CHECK_EQUAL(summary.factorChronologyHash.size(), 16);
	BOOST_CHECK_SMALL(summary.maximumCommitMeanRelativeError, 1e-15);
	BOOST_CHECK_SMALL(summary.maximumCommitCovarianceRelativeError, 1e-15);
	BOOST_REQUIRE_EQUAL(capture.capturedEvents().size(), 3);
	BOOST_CHECK_SMALL(
		(capture.capturedEvents()[1].rightHandSide - transitionOffset).norm(),
		1e-15);
	BOOST_CHECK_SMALL(
		(capture.capturedEvents()[2].rightHandSide - transformOffset).norm(),
		1e-15);
	BOOST_CHECK(!capture.capturedEvents()[0].factorId.empty());
	BOOST_CHECK(!capture.capturedEvents()[1].factorId.empty());
	BOOST_CHECK(!capture.capturedEvents()[2].factorId.empty());
	const auto replayA = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; });
	const auto replayB = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; });
	BOOST_REQUIRE_MESSAGE(replayA.valid, replayA.failureReason);
	BOOST_REQUIRE_MESSAGE(replayB.valid, replayB.failureReason);
	BOOST_REQUIRE_EQUAL(replayA.mean.size(), 1);
	BOOST_CHECK_SMALL(replayA.mean(0) - 2.0, 1e-12);
	BOOST_CHECK_SMALL((replayA.mean - replayB.mean).norm(), 1e-15);
	BOOST_CHECK_SMALL(
		(replayA.covariance - replayB.covariance).norm(), 1e-15);
}

BOOST_AUTO_TEST_CASE(authoritative_factor_commit_packet_fails_closed)
{
	KFKey state;
	state.type = KF::SAT_CLOCK;
	state.Sat = SatSys(E_Sys::GPS, 4);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(state)};
	const VectorXd mean = VectorXd::Zero(1);
	const MatrixXd covariance = MatrixXd::Identity(1, 1);
	KFMeas measurement;
	measurement.H = MatrixXd::Zero(1, 1);
	measurement.V = VectorXd::Zero(1);
	measurement.R = MatrixXd::Identity(1, 1);
	measurement.obsKeys = {state};
	SparseMatrix<double> identity(1, 1);
	identity.insert(0, 0) = 1;

	ZhangFactorCaptureBuffer sequenceCapture;
	BOOST_REQUIRE(sequenceCapture.recordMeasurement(
		GTime(), keys, mean, covariance, measurement, "/PPP",
		mean, covariance, 7, 8));
	BOOST_CHECK(!sequenceCapture.recordTransition(
		GTime(), keys, keys, identity, MatrixXd::Zero(1, 1),
		VectorXd::Zero(1), "SEQUENCE_BREAK", mean, covariance, 9, 10));
	BOOST_CHECK_EQUAL(
		sequenceCapture.summary().failureReason,
		"AUTHORITATIVE_COMMIT_SEQUENCE_BREAK");

	ZhangFactorCaptureBuffer momentCapture;
	BOOST_REQUIRE(momentCapture.recordMeasurement(
		GTime(), keys, mean, covariance, measurement, "/PPP",
		mean, covariance, 11, 12));
	BOOST_CHECK(!momentCapture.recordTransition(
		GTime(), keys, keys, identity, MatrixXd::Zero(1, 1),
		VectorXd::Zero(1), "MOMENT_BREAK", mean,
		MatrixXd::Constant(1, 1, 2.0), 12, 13));
	BOOST_CHECK_EQUAL(
		momentCapture.summary().failureReason,
		"AUTHORITATIVE_TRANSITION_COMMIT_MISMATCH");
}

BOOST_AUTO_TEST_CASE(
	authoritative_postfit_relaxation_moments_preserve_mean_and_cross_covariance)
{
	Vector3d sourceMean;
	sourceMean << 0.4, -0.7, 1.2;
	Matrix3d sourceCovariance;
	sourceCovariance <<
		4.0, 0.5, -0.25,
		0.5, 3.0, 0.4,
		-0.25, 0.4, 2.0;
	const Matrix3d transition = Matrix3d::Identity();
	Matrix3d processCovariance = Matrix3d::Zero();
	processCovariance(1, 1) = 2.5;

	VectorXd destinationMean;
	MatrixXd destinationCovariance;
	BOOST_REQUIRE(kfApplyStochasticTransitionMoments(
		sourceMean, sourceCovariance, transition, processCovariance,
		destinationMean, destinationCovariance));

	BOOST_CHECK_SMALL((destinationMean - sourceMean).norm(), 1e-15);
	const string sourceMomentHash =
		kfMomentHash(sourceMean, sourceCovariance);
	const string destinationMomentHash =
		kfMomentHash(destinationMean, destinationCovariance);
	const Matrix3d covarianceIncrement =
		destinationCovariance - sourceCovariance;
	BOOST_CHECK_SMALL(covarianceIncrement(1, 1) - 2.5, 1e-15);
	for (int row = 0; row < covarianceIncrement.rows(); row++)
	for (int column = 0; column < covarianceIncrement.cols(); column++)
	{
		if (row == 1 && column == 1)
		{
			continue;
		}
		BOOST_CHECK_SMALL(covarianceIncrement(row, column), 1e-15);
	}
	BOOST_CHECK_NE(sourceMomentHash, destinationMomentHash);
	BOOST_CHECK_EQUAL(
		destinationMomentHash,
		kfMomentHash(destinationMean, destinationCovariance));

	Matrix3d invalidProcess = processCovariance;
	invalidProcess(0, 0) = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK(!kfApplyStochasticTransitionMoments(
		sourceMean, sourceCovariance, transition, invalidProcess,
		destinationMean, destinationCovariance));
}

BOOST_AUTO_TEST_CASE(
	authoritative_postfit_noncommuting_staged_transitions_preserve_order)
{
	Vector2d sourceMean(1.0, -0.5);
	Matrix2d sourceCovariance;
	sourceCovariance << 2.0, 0.3, 0.3, 1.0;
	Matrix2d firstTransition;
	firstTransition << 1.0, 1.0, 0.0, 1.0;
	Matrix2d secondTransition;
	secondTransition << 1.0, 0.0, 1.0, 1.0;
	Matrix2d firstProcess = Matrix2d::Zero();
	firstProcess.diagonal() << 0.1, 0.2;
	Matrix2d secondProcess = Matrix2d::Zero();
	secondProcess.diagonal() << 0.05, 0.07;

	KFMeasurementTransaction transaction;
	transaction.entryCommitSequence = 40;
	transaction.provisionalCommitSequence = 40;
	const std::vector<string> expectedRelaxationLabels = {
		"KF_POSTFIT_STATE_RELAXATION/Postfit/REC_CLOCK/MDVJ/NONE",
		"KF_POSTFIT_STATE_RELAXATION/Postfit/REC_CLOCK/KOKV/NONE",
		"KF_POSTFIT_STATE_RELAXATION/Postfit/REC_CLOCK/MKEA/NONE",
		"KF_POSTFIT_STATE_RELAXATION/Postfit/REC_CLOCK/CPVG/NONE",
		"KF_POSTFIT_STATE_RELAXATION/Postfit/REC_CLOCK/UTQI/NONE",
		"KF_POSTFIT_STATE_RELAXATION/Postfit/REC_CLOCK/ALBH/NONE"};
	auto stage = [&](const Matrix2d& transition,
					 const Matrix2d& processCovariance,
					 const string& label,
					 const VectorXd& mean,
					 const MatrixXd& covariance,
					 VectorXd& destinationMean,
					 MatrixXd& destinationCovariance)
	{
		BOOST_REQUIRE(kfApplyStochasticTransitionMoments(
			mean, covariance, transition, processCovariance,
			destinationMean, destinationCovariance));
		KFStagedFactor factor;
		factor.kind = KFStagedFactorKind::STATE_TRANSITION;
		factor.label = label;
		factor.adaptiveQcTransition = true;
		factor.sourceMean = mean;
		factor.sourceCovariance = covariance;
		factor.transition = transition.sparseView();
		factor.processCovariance = processCovariance.sparseView();
		factor.destinationMean = destinationMean;
		factor.destinationCovariance = destinationCovariance;
		factor.beforeCommitSequence =
			transaction.provisionalCommitSequence;
		factor.afterCommitSequence =
			++transaction.provisionalCommitSequence;
		transaction.stagedFactors.push_back(std::move(factor));
		transaction.modelGeneration++;
	};

	VectorXd afterFirstMean;
	MatrixXd afterFirstCovariance;
	stage(firstTransition, firstProcess,
		expectedRelaxationLabels[0],
		sourceMean, sourceCovariance,
		afterFirstMean, afterFirstCovariance);
	VectorXd afterSecondMean;
	MatrixXd afterSecondCovariance;
	stage(secondTransition, secondProcess,
		expectedRelaxationLabels[1],
		afterFirstMean, afterFirstCovariance,
		afterSecondMean, afterSecondCovariance);

	VectorXd orderedMean = afterSecondMean;
	MatrixXd orderedCovariance = afterSecondCovariance;
	for (std::size_t index = 2; index < expectedRelaxationLabels.size(); index++)
	{
		Matrix2d process = Matrix2d::Zero();
		process(index % 2, index % 2) = 0.01 * (index + 1);
		VectorXd nextMean;
		MatrixXd nextCovariance;
		stage(Matrix2d::Identity(), process,
			expectedRelaxationLabels[index],
			orderedMean, orderedCovariance, nextMean, nextCovariance);
		orderedMean = std::move(nextMean);
		orderedCovariance = std::move(nextCovariance);
	}

	VectorXd reverseFirstMean;
	MatrixXd reverseFirstCovariance;
	BOOST_REQUIRE(kfApplyStochasticTransitionMoments(
		sourceMean, sourceCovariance, secondTransition, secondProcess,
		reverseFirstMean, reverseFirstCovariance));
	VectorXd reversedMean;
	MatrixXd reversedCovariance;
	BOOST_REQUIRE(kfApplyStochasticTransitionMoments(
		reverseFirstMean, reverseFirstCovariance,
		firstTransition, firstProcess,
		reversedMean, reversedCovariance));
	for (std::size_t index = 2; index < expectedRelaxationLabels.size(); index++)
	{
		Matrix2d process = Matrix2d::Zero();
		process(index % 2, index % 2) = 0.01 * (index + 1);
		VectorXd nextMean;
		MatrixXd nextCovariance;
		BOOST_REQUIRE(kfApplyStochasticTransitionMoments(
			reversedMean, reversedCovariance, Matrix2d::Identity(), process,
			nextMean, nextCovariance));
		reversedMean = std::move(nextMean);
		reversedCovariance = std::move(nextCovariance);
	}

	KFStagedFactor measurementFactor;
	measurementFactor.kind = KFStagedFactorKind::MEASUREMENT;
	measurementFactor.label = "/PPP";
	measurementFactor.sourceMean = orderedMean;
	measurementFactor.sourceCovariance = orderedCovariance;
	measurementFactor.destinationMean = orderedMean;
	measurementFactor.destinationCovariance = orderedCovariance;
	measurementFactor.beforeCommitSequence =
		transaction.provisionalCommitSequence;
	measurementFactor.afterCommitSequence =
		++transaction.provisionalCommitSequence;
	transaction.stagedFactors.push_back(std::move(measurementFactor));

	BOOST_REQUIRE_EQUAL(transaction.stagedFactors.size(), 7U);
	for (std::size_t index = 0; index < expectedRelaxationLabels.size(); index++)
	{
		BOOST_CHECK_EQUAL(
			transaction.stagedFactors[index].label,
			expectedRelaxationLabels[index]);
		BOOST_CHECK(transaction.stagedFactors[index].adaptiveQcTransition);
		BOOST_CHECK_EQUAL(
			transaction.stagedFactors[index].beforeCommitSequence,
			40 + index);
		BOOST_CHECK_EQUAL(
			transaction.stagedFactors[index].afterCommitSequence,
			41 + index);
	}
	BOOST_CHECK_EQUAL(transaction.stagedFactors.back().label, "/PPP");
	BOOST_CHECK(!transaction.stagedFactors.back().adaptiveQcTransition);
	BOOST_CHECK_EQUAL(
		transaction.stagedFactors.back().beforeCommitSequence, 46);
	BOOST_CHECK_EQUAL(
		transaction.stagedFactors.back().afterCommitSequence, 47);
	BOOST_CHECK_EQUAL(transaction.provisionalCommitSequence, 47);
	BOOST_CHECK_EQUAL(transaction.modelGeneration, 6U);
	BOOST_CHECK_SMALL(
		(transaction.stagedFactors[5].destinationMean - orderedMean).norm(),
		1e-15);
	BOOST_CHECK_SMALL(
		(transaction.stagedFactors[5].destinationCovariance
			- orderedCovariance).norm(), 1e-15);
	BOOST_CHECK_GT((orderedMean - reversedMean).norm(), 1e-6);
	BOOST_CHECK_GT((orderedCovariance - reversedCovariance).norm(), 1e-6);
}

BOOST_AUTO_TEST_CASE(
	authoritative_postfit_transaction_candidate_rejection_rolls_back_live_state_and_capture)
{
	KFKey stateKey;
	stateKey.type = KF::SAT_CLOCK;
	stateKey.Sat = SatSys(E_Sys::GPS, 3);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(stateKey)};
	VectorXd liveMean = VectorXd::Constant(1, 2.0);
	MatrixXd liveCovariance = MatrixXd::Constant(1, 1, 4.0);
	std::uint64_t liveSequence = 71;
	int liveSidecarGeneration = 17;

	KFMeas anchorMeasurement;
	anchorMeasurement.H = MatrixXd::Zero(1, 1);
	anchorMeasurement.V = VectorXd::Zero(1);
	anchorMeasurement.R = MatrixXd::Identity(1, 1);
	anchorMeasurement.obsKeys = {stateKey};
	ZhangFactorCaptureBuffer liveCapture;
	BOOST_REQUIRE(liveCapture.recordMeasurement(
		GTime(), keys, liveMean, liveCovariance, anchorMeasurement,
		"/PPP_ANCHOR", liveMean, liveCovariance, 70, liveSequence));
	const auto liveCheckpoint = liveCapture.checkpointReplay();
	const string liveChronologyHash =
		liveCapture.summary().factorChronologyHash;
	const string liveMomentIdentity = kfMomentHash(liveMean, liveCovariance);

	ZhangFactorCaptureBuffer candidateCapture = liveCapture;
	std::uint64_t candidateSequence = liveSequence;
	int candidateSidecarGeneration = liveSidecarGeneration;
	MatrixXd transition = MatrixXd::Identity(1, 1);
	MatrixXd processCovariance = MatrixXd::Constant(1, 1, 0.25);
	VectorXd candidateMean;
	MatrixXd candidateCovariance;
	BOOST_REQUIRE(kfApplyStochasticTransitionMoments(
		liveMean, liveCovariance, transition, processCovariance,
		candidateMean, candidateCovariance));
	const std::uint64_t transitionBeforeSequence = candidateSequence;
	const std::uint64_t transitionAfterSequence = ++candidateSequence;
	BOOST_REQUIRE(candidateCapture.recordTransition(
		GTime(), keys, keys, transition.sparseView(), processCovariance,
		VectorXd::Zero(1), "KF_POSTFIT_STATE_RELAXATION/REJECTED",
		candidateMean, candidateCovariance,
		transitionBeforeSequence, transitionAfterSequence));
	candidateSidecarGeneration++;

	KFMeas rejectedMeasurement = anchorMeasurement;
	// Deliberately skip sequence 72.  The candidate is allowed to fail and be
	// discarded, but the live capture/estimator packet must remain unchanged.
	BOOST_CHECK(!candidateCapture.recordMeasurement(
		GTime(), keys, candidateMean, candidateCovariance,
		rejectedMeasurement, "/PPP_REJECTED",
		candidateMean, candidateCovariance, candidateSequence + 1,
		candidateSequence + 2));
	BOOST_CHECK_EQUAL(candidateCapture.summary().failureReason,
		"AUTHORITATIVE_COMMIT_SEQUENCE_BREAK");

	const auto liveAfterReject = liveCapture.checkpointReplay();
	BOOST_CHECK_EQUAL(kfMomentHash(liveMean, liveCovariance),
		liveMomentIdentity);
	BOOST_CHECK_EQUAL(liveSequence, 71);
	BOOST_CHECK_EQUAL(liveSidecarGeneration, 17);
	BOOST_CHECK_EQUAL(candidateSequence, 72);
	BOOST_CHECK_EQUAL(candidateSidecarGeneration, 18);
	BOOST_CHECK_EQUAL(liveAfterReject.lastCommitSequence,
		liveCheckpoint.lastCommitSequence);
	BOOST_CHECK(liveAfterReject.currentKeys == liveCheckpoint.currentKeys);
	BOOST_CHECK_SMALL(
		(liveAfterReject.replayMean - liveCheckpoint.replayMean).norm(),
		1e-15);
	BOOST_CHECK_SMALL(
		(liveAfterReject.replayCovariance
			- liveCheckpoint.replayCovariance).norm(), 1e-15);
	BOOST_REQUIRE_EQUAL(liveAfterReject.events.size(),
		liveCheckpoint.events.size());
	for (std::size_t index = 0; index < liveCheckpoint.events.size(); index++)
	{
		BOOST_CHECK_EQUAL(liveAfterReject.events[index].factorId,
			liveCheckpoint.events[index].factorId);
	}
	BOOST_CHECK_EQUAL(liveCapture.summary().factorChronologyHash,
		liveChronologyHash);
	BOOST_CHECK_EQUAL(liveCapture.summary().events, 1U);
	BOOST_CHECK_EQUAL(candidateCapture.summary().events, 2U);
}

BOOST_AUTO_TEST_CASE(
	terminal_reconciliation_closes_final_model_generation_without_new_factor)
{
	Vector3d priorMean;
	priorMean << 0.2, -0.4, 0.7;
	Matrix3d priorCovariance;
	priorCovariance <<
		2.0, 0.3, -0.1,
		0.3, 1.5, 0.2,
		-0.1, 0.2, 1.0;
	MatrixXd design(1, 3);
	design << 1.0, -0.5, 0.25;
	VectorXd absoluteObservation = VectorXd::Constant(1, 3.5);
	MatrixXd initialNoise = MatrixXd::Constant(1, 1, 0.25);
	MatrixXd finalNoise = MatrixXd::Constant(1, 1, 4.0);

	auto solve = [&](const VectorXd& mean,
					 const MatrixXd& covariance,
					 const MatrixXd& noise)
	{
		const MatrixXd innovation =
			design * covariance * design.transpose() + noise;
		const MatrixXd gain = covariance * design.transpose()
			* innovation.inverse();
		const VectorXd posteriorMean = mean + gain
			* (absoluteObservation - design * mean);
		const MatrixXd identity = MatrixXd::Identity(mean.size(), mean.size());
		const MatrixXd attenuation = identity - gain * design;
		MatrixXd posteriorCovariance = attenuation * covariance
			* attenuation.transpose() + gain * noise * gain.transpose();
		posteriorCovariance = 0.5
			* (posteriorCovariance + posteriorCovariance.transpose());
		return std::pair<VectorXd, MatrixXd>{
			posteriorMean, posteriorCovariance};
	};

	KFMeasurementTransaction transaction;
	transaction.entryCommitSequence = 10;
	transaction.provisionalCommitSequence = 10;
	transaction.modelGeneration = 0;
	transaction.solveGeneration = 0;
	transaction.qcConverged = false;
	transaction.maxIterationsReached = true;
	const auto initialSolve = solve(priorMean, priorCovariance, initialNoise);
	BOOST_CHECK(initialSolve.first.allFinite());

	// The final deweight happens after the last robust solve.  Until a terminal
	// no-callback solve is performed, the measurement packet is not stageable.
	transaction.modelGeneration++;
	BOOST_CHECK_NE(transaction.modelGeneration, transaction.solveGeneration);
	const std::size_t factorsBeforeReconciliation =
		transaction.stagedFactors.size();
	const std::uint64_t sequenceBeforeReconciliation =
		transaction.provisionalCommitSequence;

	const auto reconciled = solve(priorMean, priorCovariance, finalNoise);
	transaction.solveGeneration = transaction.modelGeneration;
	transaction.terminalReconciliationPerformed = true;
	BOOST_CHECK_EQUAL(transaction.modelGeneration, transaction.solveGeneration);
	BOOST_CHECK_EQUAL(transaction.stagedFactors.size(),
		factorsBeforeReconciliation);
	BOOST_CHECK_EQUAL(transaction.provisionalCommitSequence,
		sequenceBeforeReconciliation);

	KFStagedFactor measurementFactor;
	measurementFactor.kind = KFStagedFactorKind::MEASUREMENT;
	measurementFactor.label = "/PPP";
	measurementFactor.sourceMean = priorMean;
	measurementFactor.sourceCovariance = priorCovariance;
	measurementFactor.destinationMean = reconciled.first;
	measurementFactor.destinationCovariance = reconciled.second;
	measurementFactor.measurement.H = design;
	measurementFactor.measurement.R = finalNoise;
	measurementFactor.measurement.V =
		absoluteObservation - design * priorMean;
	measurementFactor.beforeCommitSequence =
		transaction.provisionalCommitSequence;
	measurementFactor.afterCommitSequence =
		++transaction.provisionalCommitSequence;
	transaction.stagedFactors.push_back(std::move(measurementFactor));

	std::size_t measurementFactors = 0;
	for (const auto& factor : transaction.stagedFactors)
	{
		measurementFactors += factor.kind == KFStagedFactorKind::MEASUREMENT;
	}
	BOOST_CHECK_EQUAL(measurementFactors, 1U);
	BOOST_CHECK(transaction.terminalReconciliationPerformed);
	BOOST_CHECK(!transaction.qcConverged);
	BOOST_CHECK(transaction.maxIterationsReached);
	BOOST_CHECK_EQUAL(transaction.provisionalCommitSequence, 11);
	BOOST_CHECK_SMALL(
		(transaction.stagedFactors.back().destinationMean
			- reconciled.first).norm(), 1e-15);
	BOOST_CHECK_SMALL(
		(transaction.stagedFactors.back().destinationCovariance
			- reconciled.second).norm(), 1e-15);

	const auto doubleUpdated = solve(
		reconciled.first, reconciled.second, finalNoise);
	BOOST_CHECK_GT((doubleUpdated.first - reconciled.first).norm(), 1e-6);
	BOOST_CHECK_GT((doubleUpdated.second - reconciled.second).norm(), 1e-6);
}

BOOST_AUTO_TEST_CASE(least_square_state_initialisation_is_stochastic_reset_commit)
{
	KFKey retainedKey;
	retainedKey.type = KF::SAT_CLOCK;
	retainedKey.Sat = SatSys(E_Sys::GPS, 5);
	KFKey freshArcKey;
	freshArcKey.type = KF::AMBIGUITY;
	freshArcKey.str = "R0";
	freshArcKey.Sat = SatSys(E_Sys::GPS, 5);
	freshArcKey.num = static_cast<int>(E_ObsCode::L1C);

	std::vector<ZhangCapturedStateKey> keys(2);
	keys[0] = zhangCapturedStateKey(retainedKey);
	keys[1] = zhangCapturedStateKey(freshArcKey);
	Vector2d previousMean(5.0, 3.0);
	Matrix2d previousCovariance = Matrix2d::Zero();
	previousCovariance.diagonal() << 4.0, 9.0;
	KFMeas zeroMeasurement;
	zeroMeasurement.H = MatrixXd::Zero(1, 2);
	zeroMeasurement.V = VectorXd::Zero(1);
	zeroMeasurement.R = MatrixXd::Identity(1, 1);
	zeroMeasurement.obsKeys = {retainedKey};

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, previousMean, previousCovariance, zeroMeasurement,
		"/PPP_FIRST", previousMean, previousCovariance, 10, 11));

	const Vector2d sentinelMean(5.0, 0.0);
	Matrix2d sentinelCovariance = Matrix2d::Zero();
	sentinelCovariance.diagonal() << 4.0, -1.0;
	const Vector2d initialisedMean(5.0, 2.0);
	Matrix2d initialisedCovariance = Matrix2d::Zero();
	initialisedCovariance.diagonal() << 4.0, 1.0;
	SparseMatrix<double> transition;
	MatrixXd processCovariance;
	BOOST_REQUIRE(makeKFStateInitialisationTransition(
		sentinelCovariance, initialisedCovariance, {1},
		transition, processCovariance));
	BOOST_CHECK_SMALL(transition.coeff(0, 0) - 1.0, 1e-15);
	BOOST_CHECK_SMALL(transition.coeff(1, 1), 1e-15);
	BOOST_CHECK_SMALL(processCovariance(0, 0), 1e-15);
	BOOST_CHECK_SMALL(processCovariance(1, 1) - 1.0, 1e-12);
	BOOST_REQUIRE(capture.recordTransition(
		GTime(), keys, keys, transition, processCovariance,
		initialisedMean - transition * sentinelMean,
		"KF_LEAST_SQUARE_STATE_INITIALISATION/LSQ",
		initialisedMean, initialisedCovariance, 11, 12));

	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, initialisedMean, initialisedCovariance,
		zeroMeasurement, "/PPP_SECOND",
		initialisedMean, initialisedCovariance, 12, 13));
	const auto summary = capture.summary();
	BOOST_REQUIRE_MESSAGE(summary.valid, summary.failureReason);
	BOOST_CHECK_EQUAL(summary.measurements, 2);
	BOOST_CHECK_EQUAL(summary.transitions, 1);
	BOOST_CHECK_EQUAL(summary.lastCommitSequence, 13);
	BOOST_CHECK_EQUAL(summary.uniqueMeasurementFactorIds, 2);
	BOOST_REQUIRE_EQUAL(capture.capturedEvents().size(), 3);
	BOOST_CHECK_SMALL(capture.capturedEvents()[1].design.coeff(1, 1), 1e-15);
	BOOST_CHECK_SMALL(
		capture.capturedEvents()[1].covariance.coeff(1, 1) - 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(temporal_product_snapshots_survive_rectangular_arc_reinitialisation)
{
	KFKey oldArc;
	oldArc.type = KF::AMBIGUITY;
	oldArc.str = "R0";
	oldArc.Sat = SatSys(E_Sys::GPS, 1);
	oldArc.num = static_cast<int>(E_ObsCode::L1C);
	KFKey surviving = oldArc;
	surviving.str = "R1";
	std::vector<ZhangCapturedStateKey> source = {
		zhangCapturedStateKey(surviving), zhangCapturedStateKey(oldArc)};
	std::vector<ZhangCapturedStateKey> destination = {
		zhangCapturedStateKey(surviving)};

	Vector2d priorMean(3.2, -1.4);
	Matrix2d priorCovariance;
	priorCovariance << 0.4, 0.1, 0.1, 0.7;
	KFMeas measurement;
	measurement.H = MatrixXd::Identity(2, 2);
	measurement.V = Vector2d(0.05, -0.02);
	measurement.R = Matrix2d::Identity() * 0.2;
	measurement.obsKeys = {surviving, oldArc};
	const Matrix2d innovation = priorCovariance + measurement.R;
	const Matrix2d gain = priorCovariance * innovation.inverse();
	const Vector2d posteriorMean = priorMean + gain * measurement.V;
	Matrix2d posteriorCovariance = priorCovariance
		- gain * priorCovariance;
	posteriorCovariance = 0.5 *
		(posteriorCovariance + posteriorCovariance.transpose());

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), source, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	BOOST_REQUIRE(capture.bindPersistentSnapshot(
		"OLD", "OLD@A0", Vector2d(0, 1), 0));

	SparseMatrix<double> projection(1, 2);
	projection.insert(0, 0) = 1;
	BOOST_REQUIRE(capture.recordCoordinateTransform(
		GTime(), source, destination, projection,
		"local phase-coordinate reinitialisation", true));
	BOOST_REQUIRE(capture.bindPersistentSnapshot(
		"NEW", "NEW@A1", VectorXd::Ones(1), 0));
	const auto marginal = capture.persistentSnapshotMarginal();
	BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
	const auto replayed = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; });
	BOOST_REQUIRE_MESSAGE(replayed.valid, replayed.failureReason);
	BOOST_CHECK(marginal.identities == replayed.identities);
	BOOST_CHECK_SMALL((marginal.mean - replayed.mean).norm(), 1e-12);
	BOOST_CHECK_SMALL(
		(marginal.covariance - replayed.covariance).norm(), 1e-12);
	BOOST_REQUIRE_EQUAL(marginal.identities.size(), 2);
	BOOST_CHECK_EQUAL(marginal.identities[0], "OLD");
	BOOST_CHECK_EQUAL(marginal.identities[1], "NEW");
	const double expectedDifference = posteriorMean(0) - posteriorMean(1);
	const double expectedVariance = posteriorCovariance(0, 0)
		+ posteriorCovariance(1, 1) - 2 * posteriorCovariance(0, 1);
	BOOST_CHECK_SMALL(
		(marginal.mean(1) - marginal.mean(0)) - expectedDifference, 1e-12);
	BOOST_CHECK_SMALL(
		(marginal.covariance(1, 1) + marginal.covariance(0, 0)
			 - 2 * marginal.covariance(0, 1)) - expectedVariance, 1e-12);
}

BOOST_AUTO_TEST_CASE(candidate_product_snapshot_is_bound_before_its_coordinate_is_dropped)
{
	KFKey retained;
	retained.type = KF::AMBIGUITY;
	retained.str = "R0";
	retained.Sat = SatSys(E_Sys::GPS, 1);
	retained.num = static_cast<int>(E_ObsCode::L1C);
	KFKey candidateChord = retained;
	candidateChord.str = "R1";
	std::vector<ZhangCapturedStateKey> source = {
		zhangCapturedStateKey(retained),
		zhangCapturedStateKey(candidateChord)};
	std::vector<ZhangCapturedStateKey> destination = {
		zhangCapturedStateKey(retained)};

	Vector2d priorMean(2.4, -0.8);
	Matrix2d priorCovariance;
	priorCovariance << 0.5, 0.12, 0.12, 0.9;
	KFMeas measurement;
	measurement.H = MatrixXd::Identity(2, 2);
	measurement.V = Vector2d(0.03, -0.04);
	measurement.R = Matrix2d::Identity() * 0.25;
	measurement.obsKeys = {retained, candidateChord};
	const Matrix2d innovation = priorCovariance + measurement.R;
	const Matrix2d gain = priorCovariance * innovation.inverse();
	const Vector2d posteriorMean = priorMean + gain * measurement.V;
	Matrix2d posteriorCovariance = priorCovariance
		- gain * priorCovariance;
	posteriorCovariance = 0.5 *
		(posteriorCovariance + posteriorCovariance.transpose());

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), source, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	BOOST_REQUIRE(capture.bindPersistentSnapshots({
		{"OLD_PRODUCT", "OLD_PRODUCT@A0", Vector2d(1, 0), 0},
		{"CANDIDATE_PRODUCT", "CANDIDATE_PRODUCT@A0", Vector2d(0, 1), 0}}));

	SparseMatrix<double> projection(1, 2);
	projection.insert(0, 0) = 1;
	BOOST_REQUIRE(capture.recordCoordinateTransform(
		GTime(), source, destination, projection,
		"candidate chord local reinitialisation", true));
	const auto marginal = capture.persistentSnapshotMarginal();
	BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
	const auto replayed = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; });
	BOOST_REQUIRE_MESSAGE(replayed.valid, replayed.failureReason);
	BOOST_CHECK(marginal.identities == replayed.identities);
	BOOST_CHECK_SMALL((marginal.mean - replayed.mean).norm(), 1e-12);
	BOOST_CHECK_SMALL(
		(marginal.covariance - replayed.covariance).norm(), 1e-12);
	BOOST_REQUIRE_EQUAL(marginal.identities.size(), 2);
	BOOST_CHECK_EQUAL(marginal.identities[0], "OLD_PRODUCT");
	BOOST_CHECK_EQUAL(marginal.identities[1], "CANDIDATE_PRODUCT");
	const double expectedDifference = posteriorMean(1) - posteriorMean(0);
	const double expectedVariance = posteriorCovariance(0, 0)
		+ posteriorCovariance(1, 1) - 2 * posteriorCovariance(0, 1);
	const double projectedDifference =
		marginal.mean(1) - marginal.mean(0);
	const double projectedVariance =
		marginal.covariance(1, 1) + marginal.covariance(0, 0)
		- 2 * marginal.covariance(0, 1);
	BOOST_CHECK_SMALL(
		projectedDifference - expectedDifference, 1e-12);
	BOOST_CHECK_SMALL(
		projectedVariance - expectedVariance, 1e-12);

	const long long expectedCandidate = std::llround(expectedDifference);
	const long long projectedCandidate = std::llround(projectedDifference);
	const double expectedFractional = expectedDifference - expectedCandidate;
	const double projectedFractional =
		projectedDifference - projectedCandidate;
	const double expectedPerr = testScalarRoundErrorProbability(
		expectedFractional, expectedVariance);
	const double projectedPerr = testScalarRoundErrorProbability(
		projectedFractional, projectedVariance);
	const double expectedNis =
		expectedFractional * expectedFractional / expectedVariance;
	const double projectedNis =
		projectedFractional * projectedFractional / projectedVariance;
	BOOST_CHECK_EQUAL(projectedCandidate, expectedCandidate);
	BOOST_CHECK_SMALL(projectedPerr - expectedPerr, 1e-10);
	BOOST_CHECK_SMALL(projectedNis - expectedNis, 1e-10);
	BOOST_CHECK_EQUAL(
		projectedPerr <= 1e-3 && projectedNis <= 23.9281,
		expectedPerr <= 1e-3 && expectedNis <= 23.9281);
}

BOOST_AUTO_TEST_CASE(temporal_zero_row_snapshot_is_a_valid_besd_endpoint)
{
	KFKey ambiguity;
	ambiguity.type = KF::AMBIGUITY;
	ambiguity.str = "R0";
	ambiguity.Sat = SatSys(E_Sys::GPS, 2);
	ambiguity.num = static_cast<int>(E_ObsCode::L1C);
	std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(ambiguity)};

	VectorXd priorMean = VectorXd::Constant(1, 3.25);
	MatrixXd priorCovariance = MatrixXd::Constant(1, 1, 0.4);
	KFMeas measurement;
	measurement.H = MatrixXd::Identity(1, 1);
	measurement.V = VectorXd::Constant(1, 0.1);
	measurement.R = MatrixXd::Constant(1, 1, 0.2);
	measurement.obsKeys = {ambiguity};
	const double gain = 0.4 / 0.6;
	VectorXd posteriorMean = VectorXd::Constant(1, 3.25 + gain * 0.1);
	MatrixXd posteriorCovariance = MatrixXd::Constant(
		1, 1, 0.4 - gain * 0.4);

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	BOOST_REQUIRE(capture.bindPersistentSnapshot(
		"ZERO", "ZERO@A0", VectorXd::Zero(1), 0));
	BOOST_REQUIRE(capture.bindPersistentSnapshot(
		"CURRENT", "CURRENT@A0", VectorXd::Ones(1), 0));

	const auto marginal = capture.persistentSnapshotMarginal();
	BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
	BOOST_REQUIRE_EQUAL(marginal.identities.size(), 2);
	BOOST_CHECK_SMALL(marginal.mean(0), 1e-12);
	BOOST_CHECK_SMALL(marginal.covariance(0, 0), 1e-12);
	BOOST_CHECK_SMALL(
		(marginal.mean(1) - marginal.mean(0)) - posteriorMean(0), 1e-12);
	BOOST_CHECK_SMALL(
		(marginal.covariance(1, 1) + marginal.covariance(0, 0)
			- 2 * marginal.covariance(0, 1)) - posteriorCovariance(0, 0),
		1e-12);
}

BOOST_AUTO_TEST_CASE(
	persistent_snapshot_chronology_supports_correlated_measurement_row_ablation)
{
	KFKey ambiguity;
	ambiguity.type = KF::AMBIGUITY;
	ambiguity.str = "R0";
	ambiguity.Sat = SatSys(E_Sys::GPS, 4);
	ambiguity.num = static_cast<int>(E_ObsCode::L1C);
	std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(ambiguity)};

	VectorXd priorMean = VectorXd::Zero(1);
	MatrixXd priorCovariance = MatrixXd::Ones(1, 1);
	KFMeas measurement;
	measurement.H = MatrixXd::Ones(2, 1);
	measurement.V = Vector2d(1.0, 0.0);
	measurement.R = Matrix2d::Zero();
	measurement.R(0, 0) = 1.0;
	measurement.R(1, 1) = 0.01;
	measurement.R(0, 1) = 0.05;
	measurement.R(1, 0) = 0.05;
	KFKey codeObservation = ambiguity;
	codeObservation.type = KF::CODE_MEAS;
	KFKey phaseObservation = ambiguity;
	phaseObservation.type = KF::PHAS_MEAS;
	measurement.obsKeys = {codeObservation, phaseObservation};
	const Matrix2d innovation = measurement.H * priorCovariance
		* measurement.H.transpose() + measurement.R;
	const MatrixXd gain = priorCovariance * measurement.H.transpose()
		* innovation.inverse();
	const VectorXd posteriorMean = priorMean + gain * measurement.V;
	MatrixXd posteriorCovariance = priorCovariance
		- gain * measurement.H * priorCovariance;
	posteriorCovariance = 0.5
		* (posteriorCovariance + posteriorCovariance.transpose());

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	BOOST_REQUIRE(capture.bindPersistentSnapshots({
		{"ZERO", "ZERO@A0", VectorXd::Zero(1), 0},
		{"X", "X@A0", VectorXd::Ones(1), 0}}));

	const auto full = capture.persistentSnapshotMarginal();
	const auto replayedFull = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; });
	BOOST_REQUIRE_MESSAGE(full.valid, full.failureReason);
	BOOST_REQUIRE_MESSAGE(replayedFull.valid, replayedFull.failureReason);
	BOOST_CHECK(full.identities == replayedFull.identities);
	BOOST_CHECK_SMALL((full.mean - replayedFull.mean).norm(), 1e-12);
	BOOST_CHECK_SMALL(
		(full.covariance - replayedFull.covariance).norm(), 1e-12);

	BOOST_REQUIRE(capture.retainPersistentSnapshots({"X"}));
	BOOST_REQUIRE_EQUAL(capture.capturedSnapshotOperations().size(), 2);
	const auto codeOnly = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent& event, int row)
		{
			return event.observationKeys[row].type ==
				static_cast<int>(KF::CODE_MEAS);
		});
	BOOST_REQUIRE_MESSAGE(codeOnly.valid, codeOnly.failureReason);
	BOOST_REQUIRE_EQUAL(codeOnly.identities.size(), 1);
	BOOST_CHECK_EQUAL(codeOnly.identities.front(), "X");
	// N(0,1) updated by y=x+e, y=1, Var(e)=1.
	BOOST_CHECK_SMALL(codeOnly.mean(0) - 0.5, 1e-12);
	BOOST_CHECK_SMALL(codeOnly.covariance(0, 0) - 0.5, 1e-12);
	const auto phaseOnly = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent& event, int row)
		{
			return event.observationKeys[row].type ==
				static_cast<int>(KF::PHAS_MEAS);
		});
	BOOST_REQUIRE_MESSAGE(phaseOnly.valid, phaseOnly.failureReason);
	BOOST_REQUIRE_EQUAL(phaseOnly.identities.size(), 1);
	BOOST_CHECK_SMALL(phaseOnly.mean(0), 1e-12);
	BOOST_CHECK_SMALL(
		phaseOnly.covariance(0, 0) - 1.0 / 101.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(
	persistent_snapshot_replay_has_separate_zero_process_noise_control)
{
	KFKey state;
	state.type = KF::SAT_CLOCK;
	state.Sat = SatSys(E_Sys::GPS, 5);
	std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(state)};
	VectorXd priorMean = VectorXd::Zero(1);
	MatrixXd priorCovariance = MatrixXd::Ones(1, 1);
	KFMeas anchorMeasurement;
	anchorMeasurement.H = MatrixXd::Zero(1, 1);
	anchorMeasurement.V = VectorXd::Zero(1);
	anchorMeasurement.R = MatrixXd::Ones(1, 1);
	anchorMeasurement.obsKeys = {state};

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, anchorMeasurement,
		"/PPP", priorMean, priorCovariance));
	SparseMatrix<double> transition(1, 1);
	transition.insert(0, 0) = 1;
	MatrixXd processCovariance = MatrixXd::Constant(1, 1, 4.0);
	BOOST_REQUIRE(capture.recordTransition(
		GTime(), keys, keys, transition, processCovariance, "random walk"));
	BOOST_REQUIRE(capture.bindPersistentSnapshot(
		"X", "X@A0", VectorXd::Ones(1), 0));

	const auto full = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; });
	const auto zeroProcess = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; }, 0);
	BOOST_REQUIRE_MESSAGE(full.valid, full.failureReason);
	BOOST_REQUIRE_MESSAGE(zeroProcess.valid, zeroProcess.failureReason);
	BOOST_CHECK_SMALL(full.covariance(0, 0) - 5.0, 1e-12);
	BOOST_CHECK_SMALL(zeroProcess.covariance(0, 0) - 1.0, 1e-12);
	auto invalid = capture.replayPersistentSnapshotsKeepingRows(
		[](const ZhangCapturedFactorEvent&, int) { return true; }, -1);
	BOOST_CHECK(!invalid.valid);
	BOOST_CHECK_EQUAL(invalid.failureReason, "INVALID_PROCESS_NOISE_SCALE");
}

BOOST_AUTO_TEST_CASE(
	e29_product_gauge_compiler_closes_ten_exact_three_by_three_by_two_s_bases)
{
	constexpr int receivers = 3;
	constexpr int satellites = 3;
	constexpr int frequencies = 2;
	constexpr int observations =
		2 * receivers * satellites * frequencies;
	constexpr int tauOffset = 0;
	constexpr int receiverClockOffset = tauOffset + receivers;
	constexpr int satelliteClockOffset = receiverClockOffset + receivers;
	constexpr int ionosphereOffset = satelliteClockOffset + satellites;
	constexpr int receiverCodeOffset =
		ionosphereOffset + receivers * satellites;
	constexpr int satelliteCodeOffset =
		receiverCodeOffset + receivers * frequencies;
	constexpr int receiverPhaseOffset =
		satelliteCodeOffset + satellites * frequencies;
	constexpr int satellitePhaseOffset =
		receiverPhaseOffset + receivers * frequencies;
	constexpr int ambiguityOffset =
		satellitePhaseOffset + satellites * frequencies;
	constexpr int parameters =
		ambiguityOffset + receivers * satellites * frequencies;
	ZhangExactMatrix raw = zhangExactZeroMatrix(observations, parameters);
	const int mu[frequencies] = {1, 2};
	const int wavelength[frequencies] = {1, 2};
	auto ionosphere = [](int receiver, int satellite)
	{
		return receiver * satellites + satellite;
	};
	auto receiverSignal = [](int receiver, int frequency)
	{
		return receiver * frequencies + frequency;
	};
	auto satelliteSignal = [](int satellite, int frequency)
	{
		return satellite * frequencies + frequency;
	};
	auto ambiguity = [](int receiver, int satellite, int frequency)
	{
		return (receiver * satellites + satellite) * frequencies + frequency;
	};
	int row = 0;
	for (int receiver = 0; receiver < receivers; receiver++)
	for (int satellite = 0; satellite < satellites; satellite++)
	for (int frequency = 0; frequency < frequencies; frequency++)
	{
		auto common = [&](ZhangExactVector& design)
		{
			design[tauOffset + receiver] = 1;
			design[receiverClockOffset + receiver] = 1;
			design[satelliteClockOffset + satellite] = -1;
		};
		auto& code = raw[row++];
		common(code);
		code[ionosphereOffset + ionosphere(receiver, satellite)] =
			mu[frequency];
		code[receiverCodeOffset + receiverSignal(receiver, frequency)] = 1;
		code[satelliteCodeOffset + satelliteSignal(satellite, frequency)] = -1;
		auto& phase = raw[row++];
		common(phase);
		phase[ionosphereOffset + ionosphere(receiver, satellite)] =
			-mu[frequency];
		phase[receiverPhaseOffset + receiverSignal(receiver, frequency)] =
			wavelength[frequency];
		phase[satellitePhaseOffset + satelliteSignal(satellite, frequency)] =
			-wavelength[frequency];
		phase[ambiguityOffset + ambiguity(receiver, satellite, frequency)] =
			wavelength[frequency];
	}
	BOOST_REQUIRE_EQUAL(row, observations);

	MatrixXd rawDouble(observations, parameters);
	for (int r = 0; r < observations; r++)
	for (int c = 0; c < parameters; c++)
	{
		rawDouble(r, c) = raw[r][c].convert_to<double>();
	}
	Eigen::FullPivLU<MatrixXd> rawLu(rawDouble);
	rawLu.setThreshold(1e-12);
	const int estimableRank = rawLu.rank();
	BOOST_REQUIRE_GT(estimableRank, 0);

	auto chooseBasis = [&](std::vector<int> order)
	{
		std::vector<int> selected;
		int rank = 0;
		for (int column : order)
		{
			MatrixXd candidate(observations, selected.size() + 1);
			for (int existing = 0;
				 existing < static_cast<int>(selected.size()); existing++)
			{
				candidate.col(existing) = rawDouble.col(selected[existing]);
			}
			candidate.col(selected.size()) = rawDouble.col(column);
			Eigen::FullPivLU<MatrixXd> lu(candidate);
			lu.setThreshold(1e-12);
			if (lu.rank() > rank)
			{
				selected.push_back(column);
				rank++;
				if (rank == estimableRank)
				{
					break;
				}
			}
		}
		return selected;
	};
	auto selectExactColumns = [&](const std::vector<int>& columns)
	{
		ZhangExactMatrix design = zhangExactZeroMatrix(
			observations, columns.size());
		for (int r = 0; r < observations; r++)
		for (int c = 0; c < static_cast<int>(columns.size()); c++)
		{
			design[r][c] = raw[r][columns[c]];
		}
		return design;
	};
	auto exactToDouble = [](const ZhangExactMatrix& design)
	{
		MatrixXd result(design.size(), design.front().size());
		for (int r = 0; r < result.rows(); r++)
		for (int c = 0; c < result.cols(); c++)
		{
			result(r, c) = design[r][c].convert_to<double>();
		}
		return result;
	};

	std::vector<int> natural(parameters);
	std::iota(natural.begin(), natural.end(), 0);
	const auto frontendColumns = chooseBasis(natural);
	BOOST_REQUIRE_EQUAL(frontendColumns.size(), estimableRank);
	const auto frontendExact = selectExactColumns(frontendColumns);
	const MatrixXd frontendDense = exactToDouble(frontendExact);
	std::set<std::vector<int>> distinctBases;
	std::mt19937 generator(29001);
	for (int trial = 0; trial < 200 && distinctBases.size() < 10; trial++)
	{
		auto order = natural;
		std::shuffle(order.begin(), order.end(), generator);
		auto columns = chooseBasis(order);
		if (columns.size() == static_cast<std::size_t>(estimableRank))
		{
			distinctBases.insert(std::move(columns));
		}
	}
	BOOST_REQUIRE_GE(distinctBases.size(), 10);

	std::normal_distribution<double> normal(0, 1);
	int audited = 0;
	for (const auto& backendColumns : distinctBases)
	{
		if (audited++ == 10)
		{
			break;
		}
		const auto backendExact = selectExactColumns(backendColumns);
		const auto exact = zhangCompileExactProductGaugeTransform(
			frontendExact, backendExact);
		BOOST_REQUIRE_MESSAGE(exact.valid, exact.failureReason);
		BOOST_REQUIRE_EQUAL(exact.rank, estimableRank);
		for (int r = 0; r < observations; r++)
		for (int c = 0; c < estimableRank; c++)
		{
			ZhangExactRational predicted = 0;
			for (int k = 0; k < estimableRank; k++)
			{
				predicted += ZhangExactRational(frontendExact[r][k])
					* exact.transform[k][c];
			}
			BOOST_CHECK(predicted == ZhangExactRational(backendExact[r][c]));
		}

		const MatrixXd backendDense = exactToDouble(backendExact);
		const auto compiled = zhangCompileProductGaugeTransform(
			frontendDense.sparseView(0, 0),
			backendDense.sparseView(0, 0), 1e-12);
		BOOST_REQUIRE_MESSAGE(compiled.valid, compiled.failureReason);
		BOOST_CHECK_SMALL(compiled.maximumClosureError, 1e-12);
		VectorXd backendState(estimableRank);
		MatrixXd squareRoot(estimableRank, estimableRank);
		for (int r = 0; r < estimableRank; r++)
		{
			backendState(r) = normal(generator);
			for (int c = 0; c < estimableRank; c++)
			{
				squareRoot(r, c) = normal(generator);
			}
		}
		const MatrixXd backendCovariance = squareRoot * squareRoot.transpose()
			+ 0.1 * MatrixXd::Identity(estimableRank, estimableRank);
		const VectorXd frontendState = compiled.transform * backendState;
		const MatrixXd frontendCovariance =
			zhangProjectProductGaugeCovariance(
				backendCovariance, compiled.transform);
		BOOST_REQUIRE_EQUAL(frontendCovariance.rows(), estimableRank);
		BOOST_CHECK_SMALL(
			(backendDense * backendState
				- frontendDense * frontendState).cwiseAbs().maxCoeff(),
			1e-10);
		const MatrixXd backendPredictionCovariance =
			backendDense * backendCovariance * backendDense.transpose();
		const MatrixXd frontendPredictionCovariance =
			frontendDense * frontendCovariance * frontendDense.transpose();
		BOOST_CHECK_SMALL(
			(backendPredictionCovariance - frontendPredictionCovariance)
				.cwiseAbs().maxCoeff(),
			1e-10);
	}
	BOOST_CHECK_EQUAL(audited, 10);
}

BOOST_AUTO_TEST_CASE(
	e29_integer_conditioner_matches_near_zero_noise_fixed_resolve)
{
	constexpr int dimension = 12;
	constexpr int constraintsCount = 4;
	std::mt19937 generator(29002);
	std::normal_distribution<double> normal(0, 1);
	VectorXd mean(dimension);
	MatrixXd squareRoot(dimension, dimension);
	for (int row = 0; row < dimension; row++)
	{
		mean(row) = normal(generator);
		for (int column = 0; column < dimension; column++)
		{
			squareRoot(row, column) = normal(generator);
		}
	}
	const MatrixXd covariance = squareRoot * squareRoot.transpose()
		+ 0.5 * MatrixXd::Identity(dimension, dimension);
	std::vector<Eigen::Triplet<double>> triplets = {
		{0, 0, 1}, {0, 1, -1},
		{1, 2, 1}, {1, 3, 1}, {1, 4, -1},
		{2, 5, 1},
		{3, 6, 1}, {3, 7, -1}, {3, 8, 1}};
	ZhangIarFunctional constraints(constraintsCount, dimension);
	constraints.setFromTriplets(triplets.begin(), triplets.end());
	constraints.makeCompressed();
	VectorXd integers = constraints * mean;
	for (int row = 0; row < integers.size(); row++)
	{
		integers(row) = std::round(integers(row));
	}

	const auto exact = zhangConditionIntegersExact(
		mean, covariance, constraints, integers);
	const auto squareRootConditioned =
		zhangConditionIntegersSquareRootOrthogonal(
		mean, covariance, constraints, integers);
	const auto pseudo = zhangConditionIntegersPseudoObservation(
		mean, covariance, constraints, integers, 1e-8);
	BOOST_REQUIRE_MESSAGE(exact.valid, exact.failureReason);
	BOOST_REQUIRE_MESSAGE(
		squareRootConditioned.valid,
		squareRootConditioned.failureReason);
	BOOST_REQUIRE_MESSAGE(pseudo.valid, pseudo.failureReason);
	BOOST_CHECK_EQUAL(exact.constraintRank, constraintsCount);
	BOOST_CHECK_EQUAL(
		squareRootConditioned.constraintRank,
		constraintsCount);
	BOOST_CHECK_EQUAL(pseudo.constraintRank, constraintsCount);
	BOOST_CHECK_SMALL(exact.maximumConstraintResidual, 1e-10);
	BOOST_CHECK_SMALL(
		squareRootConditioned.maximumConstraintResidual,
		1e-10);
	BOOST_CHECK_SMALL(
		(exact.mean - squareRootConditioned.mean)
			.cwiseAbs().maxCoeff(),
		1e-9);
	BOOST_CHECK_SMALL(
		(exact.covariance - squareRootConditioned.covariance)
			.cwiseAbs().maxCoeff(),
		1e-9);
	BOOST_CHECK_SMALL(
		(exact.mean - pseudo.mean).cwiseAbs().maxCoeff(), 1e-9);
	BOOST_CHECK_SMALL(
		(exact.covariance - pseudo.covariance).cwiseAbs().maxCoeff(), 1e-9);

	std::vector<Eigen::Triplet<double>> redundantTriplets = triplets;
	redundantTriplets.emplace_back(4, 0, 2);
	redundantTriplets.emplace_back(4, 1, -2);
	ZhangIarFunctional redundant(constraintsCount + 1, dimension);
	redundant.setFromTriplets(
		redundantTriplets.begin(), redundantTriplets.end());
	redundant.makeCompressed();
	VectorXd redundantIntegers(constraintsCount + 1);
	redundantIntegers.head(constraintsCount) = integers;
	redundantIntegers(4) = 2 * integers(0);
	const auto rejected = zhangConditionIntegersExact(
		mean, covariance, redundant, redundantIntegers);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(
		rejected.failureReason,
		"INTEGER_CONSTRAINT_NOT_FULL_ROW_RANK");
}

BOOST_AUTO_TEST_CASE(
	accepted_measurement_families_use_obs_key_and_actual_state_support)
{
	KFKey clock;
	clock.type = KF::SAT_CLOCK;
	clock.Sat = SatSys(E_Sys::GPS, 4);
	KFKey ionosphere = clock;
	ionosphere.type = KF::IONO_STEC;
	KFKey phaseDatum = clock;
	phaseDatum.type = KF::PHASE_BIAS;
	KFKey position = clock;
	position.type = KF::REC_POS;

	ZhangCapturedFactorEvent event;
	event.kind = ZhangCapturedFactorKind::MEASUREMENT;
	event.destinationKeys = {
		zhangCapturedStateKey(clock),
		zhangCapturedStateKey(ionosphere),
		zhangCapturedStateKey(phaseDatum),
		zhangCapturedStateKey(position)};
	event.design.resize(6, 4);
	event.design.insert(0, 3) = 1;
	event.design.insert(1, 3) = 1;
	event.design.insert(2, 0) = 1;
	event.design.insert(3, 1) = 1;
	event.design.insert(4, 2) = 1;
	event.design.insert(5, 0) = 1;
	event.design.insert(5, 1) = -1;
	KFKey phaseObservation = phaseDatum;
	phaseObservation.type = KF::PHAS_MEAS;
	KFKey codeObservation = phaseDatum;
	codeObservation.type = KF::CODE_MEAS;
	KFKey pseudo = phaseDatum;
	pseudo.type = KF::PSEUDO_MEAS;
	event.observationKeys = {
		zhangCapturedStateKey(phaseObservation),
		zhangCapturedStateKey(codeObservation),
		zhangCapturedStateKey(pseudo),
		zhangCapturedStateKey(pseudo),
		zhangCapturedStateKey(pseudo),
		zhangCapturedStateKey(pseudo)};

	BOOST_CHECK(
		zhangCapturedMeasurementFamily(event, 0) ==
		ZhangCapturedMeasurementFamily::PHASE_OBSERVATION);
	BOOST_CHECK(
		zhangCapturedMeasurementFamily(event, 1) ==
		ZhangCapturedMeasurementFamily::CODE_OBSERVATION);
	BOOST_CHECK(
		zhangCapturedMeasurementFamily(event, 2) ==
		ZhangCapturedMeasurementFamily::CLOCK_FACTOR);
	BOOST_CHECK(
		zhangCapturedMeasurementFamily(event, 3) ==
		ZhangCapturedMeasurementFamily::IONOSPHERE_FACTOR);
	BOOST_CHECK(
		zhangCapturedMeasurementFamily(event, 4) ==
		ZhangCapturedMeasurementFamily::PHASE_DATUM_FACTOR);
	BOOST_CHECK(
		zhangCapturedMeasurementFamily(event, 5) ==
		ZhangCapturedMeasurementFamily::MIXED_PSEUDO_FACTOR);
}

BOOST_AUTO_TEST_CASE(
	persistent_snapshot_is_not_rebound_after_unrepresentable_rectangular_reset)
{
	Vector2d priorMean(0.4, -1.2);
	Matrix2d priorCovariance;
	priorCovariance << 0.8, 0.25, 0.25, 1.1;
	ZhangPersistentRawTargetWindow window;
	BOOST_REQUIRE(window.initialise(priorMean, priorCovariance));
	Vector2d oldPhysicalRow(0, 1);
	BOOST_REQUIRE(window.bindTarget(
		"OLD", "OLD:R0/G01@0", oldPhysicalRow, 0, 1));
	const auto beforeReset = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(beforeReset.valid, beforeReset.failureReason);
	const long long beforeCandidate = std::llround(beforeReset.mean(0));
	const double beforeFractional = beforeReset.mean(0) - beforeCandidate;
	const double beforePerr = testScalarRoundErrorProbability(
		beforeFractional, beforeReset.covariance(0, 0));
	const double beforeNis = beforeFractional * beforeFractional /
		beforeReset.covariance(0, 0);

	// The current-state projection removes x1, but the explicit immutable
	// target variable is carried with identity.  Its marginal must survive.
	MatrixXd projection(1, 2);
	projection << 1, 0;
	BOOST_REQUIRE(window.applyExactCoordinateTransform(projection));
	const auto afterReset = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(afterReset.valid, afterReset.failureReason);
	BOOST_CHECK_SMALL(afterReset.mean(0) - beforeReset.mean(0), 1e-12);
	BOOST_CHECK_SMALL(
		afterReset.covariance(0, 0) - beforeReset.covariance(0, 0), 1e-12);
	const long long afterCandidate = std::llround(afterReset.mean(0));
	const double afterFractional = afterReset.mean(0) - afterCandidate;
	const double afterPerr = testScalarRoundErrorProbability(
		afterFractional, afterReset.covariance(0, 0));
	const double afterNis = afterFractional * afterFractional /
		afterReset.covariance(0, 0);
	BOOST_CHECK_EQUAL(afterCandidate, beforeCandidate);
	BOOST_CHECK_SMALL(afterPerr - beforePerr, 1e-12);
	BOOST_CHECK_SMALL(afterNis - beforeNis, 1e-12);
	BOOST_CHECK_EQUAL(
		beforePerr <= 1e-3 && beforeNis <= 23.9281,
		afterPerr <= 1e-3 && afterNis <= 23.9281);

	// A newly initialised current row is not an exact transport proof for the
	// removed covector.  Reusing the old snapshot identity must therefore not
	// inject a zero-noise constraint.
	VectorXd newlyInitialisedRow = VectorXd::Ones(1);
	BOOST_REQUIRE(window.bindTarget(
		"OLD", "OLD:R0/G01@0", newlyInitialisedRow, 0, 2));
	const auto afterFirstRebind = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(
		afterFirstRebind.valid, afterFirstRebind.failureReason);
	BOOST_CHECK_SMALL(
		afterFirstRebind.mean(0) - afterReset.mean(0), 1e-12);
	BOOST_CHECK_SMALL(
		afterFirstRebind.covariance(0, 0)
			- afterReset.covariance(0, 0), 1e-12);

	MatrixXd measurementDesign = MatrixXd::Ones(1, 1);
	MatrixXd measurementCovariance = MatrixXd::Constant(1, 1, 0.2);
	VectorXd observation = VectorXd::Constant(1, 0.1);
	BOOST_REQUIRE(window.addAcceptedMeasurement(
		measurementDesign, measurementCovariance, observation));
	const auto beforeSecondRebind = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(
		beforeSecondRebind.valid, beforeSecondRebind.failureReason);
	BOOST_REQUIRE(window.bindTarget(
		"OLD", "OLD:R0/G01@0", newlyInitialisedRow, 0, 3));
	const auto afterSecondRebind = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(
		afterSecondRebind.valid, afterSecondRebind.failureReason);
	BOOST_CHECK_SMALL(
		afterSecondRebind.mean(0) - beforeSecondRebind.mean(0), 1e-12);
	BOOST_CHECK_SMALL(
		afterSecondRebind.covariance(0, 0)
			- beforeSecondRebind.covariance(0, 0), 1e-12);
	BOOST_CHECK_EQUAL(window.summary().exactConstraintsApplied, 0);
}

BOOST_AUTO_TEST_CASE(
	persistent_snapshot_batch_augmentation_matches_sequential_bindings)
{
	Vector3d priorMean(0.4, -0.7, 1.2);
	Matrix3d priorCovariance;
	priorCovariance <<
		0.8, 0.1, -0.05,
		0.1, 1.2, 0.2,
		-0.05, 0.2, 0.6;
	MatrixXd design(2, 3);
	design << 1, 0.2, -0.1, -0.3, 1, 0.4;
	Matrix2d measurementCovariance = 0.1 * Matrix2d::Identity();
	Vector2d observation(0.1, -0.2);

	ZhangPersistentRawTargetWindow sequential;
	ZhangPersistentRawTargetWindow batch;
	BOOST_REQUIRE(sequential.initialise(priorMean, priorCovariance));
	BOOST_REQUIRE(batch.initialise(priorMean, priorCovariance));
	BOOST_REQUIRE(sequential.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	BOOST_REQUIRE(batch.addAcceptedMeasurement(
		design, measurementCovariance, observation));

	Vector3d firstRow(1, -1, 0.5);
	Vector3d secondRow(-0.2, 0.4, 1);
	BOOST_REQUIRE(sequential.bindTarget(
		"FIRST", "FIRST@0", firstRow, 2.0, 1));
	BOOST_REQUIRE(sequential.bindTarget(
		"SECOND", "SECOND@0", secondRow, -0.3, 1));
	MatrixXd rows(2, 3);
	rows.row(0) = firstRow.transpose();
	rows.row(1) = secondRow.transpose();
	Vector2d offsets(2.0, -0.3);
	BOOST_REQUIRE(batch.bindNewTargets(
		{"FIRST", "SECOND"}, {"FIRST@0", "SECOND@0"},
		rows, offsets, 1));

	const auto sequentialMarginal = sequential.targetMarginal();
	const auto batchMarginal = batch.targetMarginal();
	BOOST_REQUIRE_MESSAGE(
		sequentialMarginal.valid, sequentialMarginal.failureReason);
	BOOST_REQUIRE_MESSAGE(batchMarginal.valid, batchMarginal.failureReason);
	BOOST_CHECK_EQUAL_COLLECTIONS(
		sequentialMarginal.identities.begin(),
		sequentialMarginal.identities.end(),
		batchMarginal.identities.begin(), batchMarginal.identities.end());
	BOOST_CHECK_SMALL(
		(sequentialMarginal.mean - batchMarginal.mean).norm(), 1e-11);
	BOOST_CHECK_SMALL(
		(sequentialMarginal.covariance - batchMarginal.covariance).norm(),
		1e-11);

	Matrix3d transition = Matrix3d::Identity();
	transition(0, 1) = 0.1;
	Matrix3d processCovariance = 0.02 * Matrix3d::Identity();
	BOOST_REQUIRE(sequential.advance(transition, processCovariance));
	BOOST_REQUIRE(batch.advance(transition, processCovariance));
	BOOST_REQUIRE(sequential.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	BOOST_REQUIRE(batch.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	const auto sequentialFuture = sequential.targetMarginal();
	const auto batchFuture = batch.targetMarginal();
	BOOST_REQUIRE_MESSAGE(sequentialFuture.valid, sequentialFuture.failureReason);
	BOOST_REQUIRE_MESSAGE(batchFuture.valid, batchFuture.failureReason);
	BOOST_CHECK_SMALL((sequentialFuture.mean - batchFuture.mean).norm(), 1e-10);
	BOOST_CHECK_SMALL(
		(sequentialFuture.covariance - batchFuture.covariance).norm(), 1e-10);
}

BOOST_AUTO_TEST_CASE(
	persistent_snapshot_lifecycle_marginalises_only_released_targets)
{
	Vector2d priorMean(0.3, -0.4);
	Matrix2d priorCovariance;
	priorCovariance << 0.8, 0.2, 0.2, 1.1;
	ZhangPersistentRawTargetWindow window;
	BOOST_REQUIRE(window.initialise(priorMean, priorCovariance));
	MatrixXd rows(3, 2);
	rows << 1, 0, 0, 1, 1, -1;
	BOOST_REQUIRE(window.bindNewTargets(
		{"A", "B", "C"}, {"A@0", "B@0", "C@0"}, rows,
		Vector3d(0.1, -0.2, 0.3), 1));
	const auto before = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(before.valid, before.failureReason);
	BOOST_REQUIRE(window.retainTargets({"A", "C"}));
	const auto after = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(after.valid, after.failureReason);
	BOOST_REQUIRE_EQUAL(after.identities.size(), 2);
	BOOST_CHECK_EQUAL(after.identities[0], "A");
	BOOST_CHECK_EQUAL(after.identities[1], "C");
	Vector2d expectedMean(before.mean(0), before.mean(2));
	Matrix2d expectedCovariance;
	expectedCovariance <<
		before.covariance(0, 0), before.covariance(0, 2),
		before.covariance(2, 0), before.covariance(2, 2);
	BOOST_CHECK_SMALL((after.mean - expectedMean).norm(), 1e-11);
	BOOST_CHECK_SMALL(
		(after.covariance - expectedCovariance).norm(), 1e-11);
}

BOOST_AUTO_TEST_CASE(retained_target_information_block_matches_scalar_schur_update)
{
	KFKey first;
	first.type = KF::SAT_CLOCK;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)
	};
	Vector2d priorMean(0.2, -0.1);
	Matrix2d priorCovariance = Matrix2d::Zero();
	priorCovariance.diagonal() << 2.0, 3.0;
	KFMeas measurement;
	measurement.H = MatrixXd::Zero(1, 2);
	measurement.H(0, 0) = 1;
	measurement.V = VectorXd::Constant(1, 0.3);
	measurement.R = MatrixXd::Constant(1, 1, 0.5);
	measurement.obsKeys = {first};
	Vector2d posteriorMean(0.44, -0.1);
	Matrix2d posteriorCovariance = Matrix2d::Zero();
	posteriorCovariance.diagonal() << 0.4, 3.0;

	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance
	));
	Vector2d targetRow(1, 0);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "GPS:WL:L1W:L2W:G01:G02", "basis-a",
		"G01:0:0->G02:0:0", {{"L1W:REC:G01", 0}}, keys,
		targetRow, 0, posteriorMean, posteriorCovariance
	));
	const auto& block = capture.currentRetainedBlock();
	BOOST_REQUIRE_MESSAGE(block.valid, block.failureReason);
	BOOST_CHECK_EQUAL(block.targetCount, 1);
	BOOST_CHECK_EQUAL(block.informationRank, 1);
	BOOST_CHECK_SMALL(block.whitenedSquaredNorm - 0.036, 1e-12);
	const auto incremental = capture.currentIncrementalTargetMarginal();
	BOOST_REQUIRE_MESSAGE(incremental.valid, incremental.failureReason);
	BOOST_CHECK_EQUAL(incremental.informationRank, 1);
	BOOST_CHECK_EQUAL(incremental.quotientValidRank, 1);
	BOOST_CHECK_EQUAL(incremental.absoluteValidRank, 1);
	BOOST_CHECK_EQUAL(incremental.storedRows, 1);
	BOOST_CHECK_EQUAL(incremental.storedColumns, 1);
}

BOOST_AUTO_TEST_CASE(raw_multi_epoch_window_matches_gaussian_kalman_posterior)
{
	// x = [C_s, C_r, B_r1, B_r2, I, N1, N2].  These signs match
	// ppp_obs.cpp: +C_r-C_s, phase -alpha*I, code +alpha*I,
	// receiver phase bias +B_rj and ambiguity +lambda_j*N_j.
	constexpr double lambda1 = 0.190293672798365;
	constexpr double lambda2 = 0.244210213424568;
	constexpr double alpha1 = 1.0;
	constexpr double alpha2 = 1.646944444444444;
	MatrixXd H = MatrixXd::Zero(4, 7);
	H.row(0) << -1, +1, +1,  0, -alpha1, lambda1, 0;
	H.row(1) << -1, +1,  0, +1, -alpha2, 0, lambda2;
	H.row(2) << -1, +1,  0,  0, +alpha1, 0, 0;
	H.row(3) << -1, +1,  0,  0, +alpha2, 0, 0;
	MatrixXd R = MatrixXd::Zero(4, 4);
	R.diagonal() << 4e-4, 6e-4, 0.16, 0.25;
	R(0, 1) = R(1, 0) = 8e-5;

	VectorXd priorMean(7);
	priorMean << 0.3, -0.1, 0.04, -0.02, 1.4, 12.2, -3.7;
	MatrixXd priorCovariance = MatrixXd::Zero(7, 7);
	priorCovariance.diagonal() << 4, 4, 1, 1, 9, 16, 16;
	VectorXd truth(7);
	truth << 0.25, -0.08, 0.03, -0.01, 1.2, 12, -4;
	VectorXd y1 = H * truth;
	y1 << y1(0) + 0.006, y1(1) - 0.004,
		y1(2) + 0.05, y1(3) - 0.08;

	MatrixXd F = MatrixXd::Identity(7, 7);
	MatrixXd Q = MatrixXd::Zero(7, 7);
	Q.diagonal() << 0.04, 0.09, 0.0025, 0.0025, 0.16, 0, 0;
	VectorXd y2 = H * truth;
	y2 << y2(0) - 0.003, y2(1) + 0.005,
		y2(2) - 0.03, y2(3) + 0.04;

	ZhangRawFactorWindow window;
	BOOST_REQUIRE_MESSAGE(
		window.initialise(priorMean, priorCovariance),
		window.lastFailureReason());
	BOOST_REQUIRE_MESSAGE(
		window.addAcceptedMeasurement(H, R, y1),
		window.lastFailureReason());
	BOOST_REQUIRE_MESSAGE(
		window.addStateTransition(F, Q), window.lastFailureReason());
	BOOST_REQUIRE_MESSAGE(
		window.addAcceptedMeasurement(H, R, y2),
		window.lastFailureReason());
	MatrixXd integerDatum = MatrixXd::Zero(1, 7);
	integerDatum(0, 5) = 1;
	integerDatum(0, 6) = -1;
	auto marginal = window.marginaliseToIntegerDatum(
		integerDatum, VectorXd::Zero(1));
	BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);

	// Independent covariance-form control.  Agreement proves numerical
	// equivalence to the finite-prior Gaussian posterior; it is deliberately
	// not used as a proof of classical observation estimability.
	VectorXd controlMean = priorMean;
	MatrixXd controlCovariance = priorCovariance;
	auto update = [&](const VectorXd& observation)
	{
		MatrixXd innovation = H * controlCovariance * H.transpose() + R;
		MatrixXd gain = controlCovariance * H.transpose()
			* innovation.ldlt().solve(MatrixXd::Identity(4, 4));
		controlMean += gain * (observation - H * controlMean);
		MatrixXd I = MatrixXd::Identity(7, 7);
		controlCovariance = (I - gain * H) * controlCovariance
			* (I - gain * H).transpose() + gain * R * gain.transpose();
		controlCovariance = 0.5
			* (controlCovariance + controlCovariance.transpose());
	};
	update(y1);
	controlMean = F * controlMean;
	controlCovariance = F * controlCovariance * F.transpose() + Q;
	update(y2);
	const double expectedMean = (integerDatum * controlMean)(0);
	const double expectedVariance =
		(integerDatum * controlCovariance * integerDatum.transpose())(0, 0);
	BOOST_TEST_MESSAGE(
		"raw-window complete-equation mean_error="
		<< std::abs(marginal.mean(0) - expectedMean)
		<< " variance_error="
		<< std::abs(marginal.covariance(0, 0) - expectedVariance));
	BOOST_CHECK_SMALL(marginal.mean(0) - expectedMean, 2e-10);
	BOOST_CHECK_SMALL(marginal.covariance(0, 0) - expectedVariance, 2e-10);
	BOOST_CHECK_EQUAL(marginal.targetRank, 1);
}

BOOST_AUTO_TEST_CASE(single_receiver_wl_is_not_observation_estimable_with_free_phase_biases)
{
	constexpr double lambda1 = 0.190293672798365;
	constexpr double lambda2 = 0.244210213424568;
	constexpr double alpha1 = 1.0;
	constexpr double alpha2 = 1.646944444444444;
	MatrixXd H = MatrixXd::Zero(4, 7);
	H.row(0) << -1, +1, +1,  0, -alpha1, lambda1, 0;
	H.row(1) << -1, +1,  0, +1, -alpha2, 0, lambda2;
	H.row(2) << -1, +1,  0,  0, +alpha1, 0, 0;
	H.row(3) << -1, +1,  0,  0, +alpha2, 0, 0;
	VectorXd target = VectorXd::Zero(7);
	target(5) = 1;
	target(6) = -1;
	MatrixXd augmented(H.rows() + 1, H.cols());
	augmented.topRows(H.rows()) = H;
	augmented.bottomRows(1) = target.transpose();
	Eigen::FullPivLU<MatrixXd> observationRows(H);
	Eigen::FullPivLU<MatrixXd> augmentedRows(augmented);
	BOOST_CHECK_EQUAL(augmentedRows.rank(), observationRows.rank() + 1);
}

BOOST_AUTO_TEST_CASE(network_fundamental_cycle_is_primitive_and_nuisance_orthogonal)
{
	// Edge order: R0-S0, R0-S1, R1-S0, R1-S1.  Removing one node datum
	// gives the real-valued receiver/satellite phase-bias incidence B.
	MatrixXd B = MatrixXd::Zero(4, 3);
	B.row(0) << 1, -1,  0;
	B.row(1) << 1,  0, -1;
	B.row(2) << 0, -1,  0;
	B.row(3) << 0,  0, -1;
	Vector4d cycle;
	cycle << 1, -1, -1, 1;
	BOOST_CHECK_SMALL((B.transpose() * cycle).norm(), 1e-14);
	int coefficientGcd = 0;
	for (int index = 0; index < cycle.size(); index++)
	{
		coefficientGcd = std::gcd(
			coefficientGcd, std::abs(static_cast<int>(cycle(index))));
	}
	BOOST_CHECK_EQUAL(coefficientGcd, 1);
	BOOST_CHECK_EQUAL(Eigen::FullPivLU<MatrixXd>(B).rank(), 3);
	MatrixXd complete(4, 4);
	complete.leftCols(3) = B;
	complete.col(3) = cycle;
	BOOST_CHECK_EQUAL(Eigen::FullPivLU<MatrixXd>(complete).rank(), 4);
}

BOOST_AUTO_TEST_CASE(six_legal_trees_and_receiver_roots_replay_identical_raw_factors)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	SatSys g03(E_Sys::GPS, 3);
	std::set<ZhangGraphEdge> edges = {
		{"R0", g01}, {"R0", g02}, {"R0", g03},
		{"R1", g01}, {"R1", g02}, {"R1", g03},
		{"R2", g01}, {"R2", g02}, {"R2", g03}
	};
	std::vector<ZhangGraphBasis> trees;
	trees.push_back(zhangBuildSpanningTree(edges, "R0"));
	trees.push_back(zhangBuildSpanningTree(
		edges, "R0", {{"R0", g01}, {"R1", g01}, {"R1", g02},
			{"R2", g02}, {"R2", g03}}));
	trees.push_back(zhangBuildSpanningTree(
		edges, "R0", {{"R0", g03}, {"R1", g02}, {"R1", g03},
			{"R2", g01}, {"R2", g03}}));
	// Explicit receiver-root changes.  Preferred edges make the physical tree
	// sets distinct, so this is not merely the same tree traversed from another
	// root.
	trees.push_back(zhangBuildSpanningTree(
		edges, "R1", {{"R1", g01}, {"R1", g02}, {"R1", g03},
			{"R0", g01}, {"R2", g02}}));
	trees.push_back(zhangBuildSpanningTree(
		edges, "R2", {{"R2", g01}, {"R2", g02}, {"R2", g03},
			{"R0", g03}, {"R1", g02}}));
	std::mt19937 generator(20260806);
	std::uniform_real_distribution<double> quality(0.0, 1.0);
	for (int attempt = 0; attempt < 100 && trees.size() < 6; attempt++)
	{
		std::map<ZhangGraphEdge, double> seededQuality;
		for (const auto& edge : edges)
		{
			seededQuality[edge] = quality(generator);
		}
		ZhangGraphBasis candidate = zhangBuildSpanningTree(
			edges, "R0", {}, seededQuality);
		bool duplicate = std::any_of(
			trees.begin(), trees.end(), [&](const auto& existing)
			{
				return existing.treeEdges == candidate.treeEdges;
			});
		if (!duplicate)
		{
			trees.push_back(std::move(candidate));
		}
	}
	BOOST_REQUIRE_EQUAL(trees.size(), 6);
	for (const auto& tree : trees)
	{
		BOOST_REQUIRE(tree.connected);
	}
	std::set<std::set<ZhangGraphEdge>> distinctTrees;
	for (const auto& tree : trees)
	{
		distinctTrees.insert(tree.treeEdges);
	}
	BOOST_REQUIRE_EQUAL(distinctTrees.size(), 6);

	// The four coordinate systems all map to one physical state containing
	// four continuous nuisance terms followed by the nine physical arc
	// ambiguities.  Thus every replay receives byte-identical observations.
	const int nuisanceCount = 4; // satellite clock, receiver clock, phase, iono
	const int arcCount = edges.size();
	const int stateCount = nuisanceCount + arcCount;
	VectorXd physicalMean = VectorXd::Zero(stateCount);
	physicalMean.head(nuisanceCount) << 0.3, -0.2, 0.05, 1.1;
	for (int index = 0; index < arcCount; index++)
	{
		physicalMean(nuisanceCount + index) = index - 3.25;
	}
	MatrixXd physicalCovariance = MatrixXd::Identity(stateCount, stateCount);
	physicalCovariance.topLeftCorner(nuisanceCount, nuisanceCount) *= 4;
	MatrixXd physicalDesign = MatrixXd::Zero(arcCount, stateCount);
	int edgeRow = 0;
	for (const auto& ignored : edges)
	{
		physicalDesign(edgeRow, 0) = -1;
		physicalDesign(edgeRow, 1) = +1;
		physicalDesign(edgeRow, 2) = (edgeRow % 2 == 0) ? +1 : -1;
		physicalDesign(edgeRow, 3) = 0.2 + 0.03 * edgeRow;
		physicalDesign(edgeRow, nuisanceCount + edgeRow) = 1;
		edgeRow++;
	}
	MatrixXd observationCovariance = 0.04
		* MatrixXd::Identity(arcCount, arcCount);
	VectorXd observation = physicalDesign * physicalMean;
	for (int row = 0; row < observation.size(); row++)
	{
		observation(row) += 0.002 * (row - 4);
	}
	VectorXd secondObservation = physicalDesign * physicalMean;
	for (int row = 0; row < secondObservation.size(); row++)
	{
		secondObservation(row) -= 0.0015 * (row - 3);
	}
	MatrixXd physicalProcess = MatrixXd::Zero(stateCount, stateCount);
	physicalProcess.topLeftCorner(nuisanceCount, nuisanceCount).diagonal()
		<< 0.04, 0.09, 0.0025, 0.16;

	// Physical product datum: G02-G01 in a fixed product tree.
	ZhangGraphBasis productTree = trees[1];
	ZhangCanonicalIntegerAudit productAudit =
		zhangCanonicalIntegerAudit(productTree);
	BOOST_REQUIRE(productAudit.valid);
	BOOST_REQUIRE(productAudit.denseCanonicalMaterialised);
	VectorXd productDatumCoordinate = VectorXd::Zero(arcCount);
	const auto& productDifference =
		productAudit.satelliteDatumSingleDifferences.front();
	for (int treeColumn = 0;
		 treeColumn < static_cast<int>(productAudit.treeEdges.size());
		 treeColumn++)
	{
		productDatumCoordinate(treeColumn) =
			productDifference[treeColumn].convert_to<double>();
	}
	std::map<ZhangGraphEdge, int> physicalEdgeIndex;
	edgeRow = 0;
	for (const auto& edge : edges)
	{
		physicalEdgeIndex[edge] = edgeRow++;
	}
	auto canonicalToPhysicalArcs = [&](const ZhangCanonicalIntegerAudit& audit)
	{
		std::vector<ZhangGraphEdge> canonicalArcs = audit.treeEdges;
		canonicalArcs.insert(
			canonicalArcs.end(), audit.chordEdges.begin(), audit.chordEdges.end());
		MatrixXd result = MatrixXd::Zero(arcCount, arcCount);
		for (int auditRow = 0; auditRow < arcCount; auditRow++)
		for (int column = 0; column < arcCount; column++)
		{
			result(physicalEdgeIndex.at(canonicalArcs[auditRow]), column) =
				audit.canonicalToArc[auditRow][column].convert_to<double>();
		}
		return result;
	};
	const MatrixXd productCanonicalToPhysical =
		canonicalToPhysicalArcs(productAudit);
	const VectorXd physicalDatum = productCanonicalToPhysical.transpose()
		.fullPivLu().solve(productDatumCoordinate);

	std::vector<double> means;
	std::vector<double> variances;
	std::vector<double> roundErrorProbabilities;
	std::vector<double> integerNis;
	std::vector<long long> integerCandidates;
	std::vector<bool> reliableStates;
	for (const auto& tree : trees)
	{
		ZhangCanonicalIntegerAudit audit = zhangCanonicalIntegerAudit(tree);
		BOOST_REQUIRE(audit.valid);
		BOOST_REQUIRE(audit.denseCanonicalMaterialised);
		MatrixXd canonicalToPhysical = MatrixXd::Identity(
			stateCount, stateCount);
		canonicalToPhysical.bottomRightCorner(arcCount, arcCount) =
			canonicalToPhysicalArcs(audit);
		MatrixXd physicalToCanonical = canonicalToPhysical.inverse();
		VectorXd coordinateMean = physicalToCanonical * physicalMean;
		MatrixXd coordinateCovariance = physicalToCanonical
			* physicalCovariance * physicalToCanonical.transpose();
		MatrixXd coordinateDesign = physicalDesign * canonicalToPhysical;
		MatrixXd target = MatrixXd::Zero(1, stateCount);
		for (int column = 0; column < arcCount; column++)
		{
			target(0, nuisanceCount + column) = physicalDatum(column);
		}
		target *= canonicalToPhysical;

		ZhangRawFactorWindow replay;
		BOOST_REQUIRE_MESSAGE(
			replay.initialise(coordinateMean, coordinateCovariance),
			replay.lastFailureReason());
		BOOST_REQUIRE_MESSAGE(
			replay.addAcceptedMeasurement(
				coordinateDesign, observationCovariance, observation),
			replay.lastFailureReason());
		const MatrixXd coordinateProcess = physicalToCanonical
			* physicalProcess * physicalToCanonical.transpose();
		BOOST_REQUIRE_MESSAGE(
			replay.addStateTransition(
				MatrixXd::Identity(stateCount, stateCount), coordinateProcess),
			replay.lastFailureReason());
		BOOST_REQUIRE_MESSAGE(
			replay.addAcceptedMeasurement(
				coordinateDesign, observationCovariance, secondObservation),
			replay.lastFailureReason());
		auto marginal = replay.marginaliseToIntegerDatum(
			target, VectorXd::Zero(1));
		BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
		means.push_back(marginal.mean(0));
		variances.push_back(marginal.covariance(0, 0));
		const double fractional = marginal.mean(0)
			- std::round(marginal.mean(0));
		integerCandidates.push_back(std::llround(marginal.mean(0)));
		roundErrorProbabilities.push_back(
			testScalarRoundErrorProbability(
				fractional, marginal.covariance(0, 0)));
		integerNis.push_back(
			fractional * fractional / marginal.covariance(0, 0));
		reliableStates.push_back(
			roundErrorProbabilities.back() <= 1e-3
			&& integerNis.back() <= 23.9281);
	}
	double maximumMeanDifference = 0;
	double maximumVarianceDifference = 0;
	double maximumPerrDifference = 0;
	double maximumNisDifference = 0;
	for (int strategy = 1; strategy < static_cast<int>(trees.size()); strategy++)
	{
		maximumMeanDifference = std::max(
			maximumMeanDifference, std::abs(means[strategy] - means[0]));
		maximumVarianceDifference = std::max(
			maximumVarianceDifference,
			std::abs(variances[strategy] - variances[0]));
		maximumPerrDifference = std::max(
			maximumPerrDifference,
			std::abs(
				roundErrorProbabilities[strategy]
				- roundErrorProbabilities[0]));
		maximumNisDifference = std::max(
			maximumNisDifference,
			std::abs(integerNis[strategy] - integerNis[0]));
		BOOST_CHECK_SMALL(means[strategy] - means[0], 1e-10);
		BOOST_CHECK_SMALL(variances[strategy] - variances[0], 1e-10);
		BOOST_CHECK_SMALL(
			roundErrorProbabilities[strategy]
			- roundErrorProbabilities[0], 1e-10);
		BOOST_CHECK_SMALL(
			integerNis[strategy] - integerNis[0], 1e-10);
		BOOST_CHECK_EQUAL(integerCandidates[strategy], integerCandidates[0]);
		BOOST_CHECK_EQUAL(reliableStates[strategy], reliableStates[0]);
	}
	BOOST_TEST_MESSAGE(
		"six-tree/root same-factor maximum_mean_difference="
		<< maximumMeanDifference
		<< " maximum_variance_difference="
		<< maximumVarianceDifference
		<< " maximum_perr_difference="
		<< maximumPerrDifference
		<< " maximum_nis_difference="
		<< maximumNisDifference);
}

BOOST_AUTO_TEST_CASE(
	satellite_reference_change_preserves_raw_integer_mean_variance_perr_and_nis)
{
	// Physical satellite potentials are only observed through differences.
	// Ref-G01 coordinates are [s2-s1,s3-s1], while Ref-G02 coordinates are
	// [s1-s2,s3-s2].  Both describe the same target s3-s1.
	Vector3d physicalMean(0.35, -1.2, 2.45);
	Matrix3d physicalCovariance;
	physicalCovariance <<
		0.9, 0.12, -0.04,
		0.12, 0.7, 0.08,
		-0.04, 0.08, 1.1;
	MatrixXd refG01(2, 3);
	refG01 << -1, 1, 0, -1, 0, 1;
	MatrixXd refG02(2, 3);
	refG02 << 1, -1, 0, 0, -1, 1;
	MatrixXd designG01(3, 2);
	designG01 << 1, 0, 0, 1, -1, 1;
	MatrixXd designG02(3, 2);
	designG02 << -1, 0, -1, 1, 0, 1;
	MatrixXd physicalDifferenceDesign(3, 3);
	physicalDifferenceDesign << -1, 1, 0, -1, 0, 1, 0, -1, 1;
	Vector3d observation = physicalDifferenceDesign * physicalMean;
	observation += Vector3d(0.006, -0.004, 0.002);
	Matrix3d observationCovariance = 0.04 * Matrix3d::Identity();

	std::vector<double> means;
	std::vector<double> variances;
	std::vector<double> perrs;
	std::vector<double> nises;
	std::vector<long long> candidates;
	for (int reference = 0; reference < 2; reference++)
	{
		const MatrixXd& coordinate = reference == 0 ? refG01 : refG02;
		const MatrixXd& design = reference == 0 ? designG01 : designG02;
		VectorXd target(2);
		target = reference == 0 ? Vector2d(0, 1) : Vector2d(-1, 1);
		ZhangRawFactorWindow replay;
		BOOST_REQUIRE(replay.initialise(
			coordinate * physicalMean,
			coordinate * physicalCovariance * coordinate.transpose()));
		BOOST_REQUIRE(replay.addAcceptedMeasurement(
			design, observationCovariance, observation));
		const auto marginal = replay.marginaliseToIntegerDatum(
			target.transpose(), VectorXd::Zero(1));
		BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
		const long long candidate = std::llround(marginal.mean(0));
		const double fractional = marginal.mean(0) - candidate;
		means.push_back(marginal.mean(0));
		variances.push_back(marginal.covariance(0, 0));
		candidates.push_back(candidate);
		perrs.push_back(testScalarRoundErrorProbability(
			fractional, marginal.covariance(0, 0)));
		nises.push_back(fractional * fractional / marginal.covariance(0, 0));
	}
	BOOST_CHECK_SMALL(means[1] - means[0], 1e-10);
	BOOST_CHECK_SMALL(variances[1] - variances[0], 1e-10);
	BOOST_CHECK_EQUAL(candidates[1], candidates[0]);
	BOOST_CHECK_SMALL(perrs[1] - perrs[0], 1e-10);
	BOOST_CHECK_SMALL(nises[1] - nises[0], 1e-10);
	BOOST_CHECK_EQUAL(
		perrs[0] <= 1e-3 && nises[0] <= 23.9281,
		perrs[1] <= 1e-3 && nises[1] <= 23.9281);
}

BOOST_AUTO_TEST_CASE(incremental_fixed_lag_matches_dense_batch_and_kalman_at_2_5_10_epochs)
{
	const Vector2d boundaryMean(0.2, -0.1);
	Matrix2d boundaryCovariance;
	boundaryCovariance << 1.0, 0.2, 0.2, 0.8;
	Matrix2d transition;
	transition << 1.0, 0.1, 0.0, 1.0;
	Matrix2d processCovariance;
	processCovariance << 0.04, 0.006, 0.006, 0.02;
	Matrix2d design;
	design << 1.0, 0.3, -0.2, 1.0;
	Matrix2d measurementCovariance;
	measurementCovariance << 0.09, 0.01, 0.01, 0.16;
	const MatrixXd noLocal = MatrixXd::Zero(2, 0);

	for (int epochCount : {2, 5, 10})
	{
		std::vector<Vector2d> observations;
		Vector2d truth(0.35, -0.25);
		for (int epoch = 0; epoch < epochCount; epoch++)
		{
			if (epoch > 0)
			{
				truth = transition * truth;
			}
			Vector2d observation = design * truth;
			observation(0) += 0.004 * (epoch - 2);
			observation(1) -= 0.003 * (epoch + 1);
			observations.push_back(observation);
		}

		Vector2d kalmanMean = boundaryMean;
		Matrix2d kalmanCovariance = boundaryCovariance;
		auto kalmanUpdate = [&](const Vector2d& observation)
		{
			const Matrix2d innovation = design * kalmanCovariance
				* design.transpose() + measurementCovariance;
			const Matrix2d gain = kalmanCovariance * design.transpose()
				* innovation.inverse();
			kalmanMean += gain * (observation - design * kalmanMean);
			const Matrix2d identity = Matrix2d::Identity();
			kalmanCovariance = (identity - gain * design) * kalmanCovariance
				* (identity - gain * design).transpose()
				+ gain * measurementCovariance * gain.transpose();
			kalmanCovariance = 0.5
				* (kalmanCovariance + kalmanCovariance.transpose());
		};

		ZhangIncrementalFixedLagSquareRoot incremental(3);
		BOOST_REQUIRE(incremental.initialise(
			boundaryMean, boundaryCovariance));
		BOOST_REQUIRE(incremental.addLatestMeasurement(
			design, noLocal, measurementCovariance, observations[0]));
		kalmanUpdate(observations[0]);
		for (int epoch = 1; epoch < epochCount; epoch++)
		{
			kalmanMean = transition * kalmanMean;
			kalmanCovariance = transition * kalmanCovariance
				* transition.transpose() + processCovariance;
			BOOST_REQUIRE(incremental.advance(
				transition, processCovariance));
			BOOST_REQUIRE(incremental.addLatestMeasurement(
				design, noLocal, measurementCovariance,
				observations[epoch]));
			kalmanUpdate(observations[epoch]);
		}

		const int stateSize = 2 * epochCount;
		const int rows = 2 + 2 * epochCount + 2 * (epochCount - 1);
		MatrixXd batchFactor = MatrixXd::Zero(rows, stateSize);
		VectorXd batchRhs = VectorXd::Zero(rows);
		const Matrix2d boundaryWeight = boundaryCovariance.llt()
			.matrixL().solve(Matrix2d::Identity());
		const Matrix2d processWeight = processCovariance.llt()
			.matrixL().solve(Matrix2d::Identity());
		const Matrix2d measurementWeight = measurementCovariance.llt()
			.matrixL().solve(Matrix2d::Identity());
		int row = 0;
		batchFactor.block<2, 2>(row, 0) = boundaryWeight;
		batchRhs.segment<2>(row) = boundaryWeight * boundaryMean;
		row += 2;
		for (int epoch = 0; epoch < epochCount; epoch++)
		{
			batchFactor.block<2, 2>(row, 2 * epoch) =
				measurementWeight * design;
			batchRhs.segment<2>(row) =
				measurementWeight * observations[epoch];
			row += 2;
			if (epoch + 1 < epochCount)
			{
				batchFactor.block<2, 2>(row, 2 * epoch) =
					-processWeight * transition;
				batchFactor.block<2, 2>(row, 2 * (epoch + 1)) =
					processWeight;
				row += 2;
			}
		}
		const MatrixXd batchInformation = batchFactor.transpose() * batchFactor;
		const MatrixXd batchCovariance = batchInformation.inverse();
		const VectorXd batchMean = batchCovariance
			* batchFactor.transpose() * batchRhs;
		const Vector2d batchLatestMean = batchMean.tail<2>();
		const Matrix2d batchLatestCovariance =
			batchCovariance.bottomRightCorner<2, 2>();

		const auto qrLatest = incremental.latestMarginal();
		BOOST_REQUIRE_MESSAGE(qrLatest.valid, qrLatest.failureReason);
		BOOST_CHECK_SMALL((qrLatest.mean - batchLatestMean).norm(), 1e-10);
		BOOST_CHECK_SMALL(
			(qrLatest.covariance - batchLatestCovariance).norm()
				/ batchLatestCovariance.norm(),
			1e-10);
		BOOST_CHECK_SMALL((qrLatest.mean - kalmanMean).norm(), 1e-10);
		BOOST_CHECK_SMALL(
			(qrLatest.covariance - kalmanCovariance).norm()
				/ kalmanCovariance.norm(),
			1e-10);
		const auto summary = incremental.summary();
		BOOST_REQUIRE(summary.valid);
		BOOST_CHECK_LE(summary.activeEpochs, 3);
		BOOST_CHECK_LE(summary.storedColumns, 6);
		BOOST_CHECK_LE(summary.storedRows, 6);
		BOOST_TEST_MESSAGE(
			"incremental-equivalence epochs=" << epochCount
			<< " mean_error=" << (qrLatest.mean - batchLatestMean).norm()
			<< " covariance_relative_error="
			<< (qrLatest.covariance - batchLatestCovariance).norm()
				/ batchLatestCovariance.norm()
			<< " stored=" << summary.storedRows << "x"
			<< summary.storedColumns);
	}
}

BOOST_AUTO_TEST_CASE(raw_square_root_boundary_matches_kalman_with_semidefinite_process)
{
	for (int epochCount : {2, 5, 10})
	{
		Vector2d mean(0.2, -0.4);
		Matrix2d covariance;
		covariance << 0.8, 0.15, 0.15, 0.6;
		ZhangIncrementalRawSquareRoot incremental;
		BOOST_REQUIRE(incremental.initialise(mean, covariance));
		double innovationSquaredNorm = 0;
		for (int epoch = 0; epoch < epochCount; epoch++)
		{
			Matrix2d design;
			design << 1.0, 0.35, -0.2, 1.0;
			Matrix2d measurementCovariance;
			measurementCovariance << 0.09, 0.01, 0.01, 0.16;
			Vector2d observation(
				0.3 + 0.04 * epoch, -0.1 + 0.02 * epoch);
			const Vector2d innovation = observation - design * mean;
			const Matrix2d innovationCovariance = design * covariance
				* design.transpose() + measurementCovariance;
			innovationSquaredNorm += innovation.dot(
				innovationCovariance.ldlt().solve(innovation));
			const Matrix2d gain = covariance * design.transpose()
				* innovationCovariance.inverse();
			mean += gain * innovation;
			covariance -= gain * design * covariance;
			covariance = 0.5 * (covariance + covariance.transpose());
			BOOST_REQUIRE(incremental.addAcceptedMeasurement(
				design, measurementCovariance, observation));
			VectorXd actualMean;
			MatrixXd actualCovariance;
			BOOST_REQUIRE(incremental.currentMarginal(
				actualMean, actualCovariance));
			BOOST_CHECK_SMALL((actualMean - mean).norm(), 1e-11);
			BOOST_CHECK_SMALL(
				(actualCovariance - covariance).norm() / covariance.norm(),
				1e-11);

			MatrixXd targetRow(1, 2);
			targetRow << 1, -1;
			VectorXd targetOffset(1);
			targetOffset << 3;
			const auto target = incremental.marginaliseTargets(
				targetRow, targetOffset);
			BOOST_REQUIRE_MESSAGE(target.valid, target.failureReason);
			BOOST_CHECK_SMALL(
				target.mean(0) - (mean(0) - mean(1) + 3), 1e-11);
			BOOST_CHECK_SMALL(
				target.covariance(0, 0)
					- (targetRow * covariance * targetRow.transpose())(0, 0),
				1e-11);

			if (epoch + 1 == epochCount)
			{
				continue;
			}
			Matrix2d transition;
			transition << 1, 0.1, 0, 1;
			Matrix2d processCovariance = Matrix2d::Zero();
			processCovariance(0, 0) = 0.01;
			mean = transition * mean;
			covariance = transition * covariance * transition.transpose()
				+ processCovariance;
			BOOST_REQUIRE(incremental.advance(
				transition, processCovariance));
			if (epoch == 2)
			{
				Matrix2d exact;
				exact << 1, 1, 0, 1;
				Vector2d shift(2, -1);
				mean = exact * mean + shift;
				covariance = exact * covariance * exact.transpose();
				BOOST_REQUIRE(incremental.applyExactCoordinateTransform(
					exact, shift));
			}
		}
		const auto summary = incremental.summary();
		BOOST_REQUIRE_MESSAGE(summary.valid, summary.failureReason);
		BOOST_CHECK_EQUAL(summary.batchOrthogonalDof, 2 * epochCount);
		BOOST_CHECK_SMALL(
			summary.batchOrthogonalSquaredNorm - innovationSquaredNorm, 1e-10);
		BOOST_CHECK_EQUAL(summary.storedRows, 2);
		BOOST_CHECK_EQUAL(summary.storedColumns, 2);
		BOOST_TEST_MESSAGE(
			"raw-square-root epochs=" << epochCount
			<< " orthogonal=" << summary.batchOrthogonalSquaredNorm
			<< "/" << summary.batchOrthogonalDof
			<< " stored=" << summary.storedRows << "x"
			<< summary.storedColumns);
	}
}

BOOST_AUTO_TEST_CASE(raw_square_root_boundary_preserves_deterministic_prior_subspace)
{
	Vector2d mean(0.5, -2.0);
	Matrix2d covariance = Matrix2d::Zero();
	covariance(0, 0) = 0.4;
	ZhangIncrementalRawSquareRoot incremental;
	BOOST_REQUIRE(incremental.initialise(mean, covariance));

	Matrix2d design;
	design << 1, 0.25, -0.4, 1;
	Matrix2d measurementCovariance = Matrix2d::Identity() * 0.1;
	Vector2d observation(0.2, -1.8);
	const Vector2d innovation = observation - design * mean;
	const Matrix2d innovationCovariance = design * covariance
		* design.transpose() + measurementCovariance;
	const Matrix2d gain = covariance * design.transpose()
		* innovationCovariance.inverse();
	mean += gain * innovation;
	covariance -= gain * design * covariance;
	covariance = 0.5 * (covariance + covariance.transpose());
	BOOST_REQUIRE(incremental.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	VectorXd actualMean;
	MatrixXd actualCovariance;
	BOOST_REQUIRE(incremental.currentMarginal(actualMean, actualCovariance));
	BOOST_CHECK_SMALL((actualMean - mean).norm(), 1e-12);
	BOOST_CHECK_SMALL((actualCovariance - covariance).norm(), 1e-12);
	BOOST_CHECK_SMALL(actualMean(1) + 2.0, 1e-15);
	BOOST_CHECK_SMALL(actualCovariance.row(1).norm(), 1e-15);

	Matrix2d transition = Matrix2d::Identity();
	Matrix2d processCovariance = Matrix2d::Zero();
	processCovariance(1, 1) = 0.03;
	mean = transition * mean;
	covariance = transition * covariance * transition.transpose()
		+ processCovariance;
	BOOST_REQUIRE(incremental.advance(transition, processCovariance));
	BOOST_REQUIRE(incremental.currentMarginal(actualMean, actualCovariance));
	BOOST_CHECK_SMALL((actualMean - mean).norm(), 1e-12);
	BOOST_CHECK_SMALL((actualCovariance - covariance).norm(), 1e-12);
	BOOST_CHECK_EQUAL(incremental.summary().storedColumns, 2);
}

BOOST_AUTO_TEST_CASE(rejected_measurement_key_chain_is_transactional)
{
	KFKey firstKey;
	firstKey.type = KF::PHASE_BIAS;
	firstKey.Sat = SatSys(E_Sys::GPS, 1);
	firstKey.num = static_cast<int>(E_ObsCode::L1W);
	KFKey wrongKey = firstKey;
	wrongKey.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(firstKey)};
	const std::vector<ZhangCapturedStateKey> wrongKeys = {
		zhangCapturedStateKey(wrongKey)};

	auto makeMeasurement = [](long int epoch, double innovation)
	{
		KFMeas measurement;
		measurement.time.bigTime = epoch;
		measurement.H = MatrixXd::Ones(1, 1);
		measurement.V = VectorXd::Constant(1, innovation);
		measurement.R = MatrixXd::Constant(1, 1, 0.5);
		return measurement;
	};
	auto posterior = [](const VectorXd& mean, const MatrixXd& covariance,
		const KFMeas& measurement)
	{
		const MatrixXd innovationCovariance = measurement.H * covariance *
			measurement.H.transpose() + measurement.R;
		const MatrixXd gain = covariance * measurement.H.transpose() *
			innovationCovariance.inverse();
		VectorXd nextMean = mean + gain * measurement.V;
		MatrixXd nextCovariance = covariance - gain * measurement.H * covariance;
		nextCovariance = 0.5 * (nextCovariance + nextCovariance.transpose());
		return std::pair{nextMean, nextCovariance};
	};

	ZhangFactorCaptureBuffer capture;
	const VectorXd priorMean = VectorXd::Zero(1);
	const MatrixXd priorCovariance = MatrixXd::Identity(1, 1);
	const auto firstMeasurement = makeMeasurement(1000, 0.2);
	const auto [firstMean, firstCovariance] = posterior(
		priorMean, priorCovariance, firstMeasurement);
	BOOST_REQUIRE(capture.recordMeasurement(
		firstMeasurement.time, keys, priorMean, priorCovariance,
		firstMeasurement, "/PPP", firstMean, firstCovariance));

	const auto secondMeasurement = makeMeasurement(1030, -0.1);
	const auto [secondMean, secondCovariance] = posterior(
		firstMean, firstCovariance, secondMeasurement);
	BOOST_CHECK(!capture.recordMeasurement(
		secondMeasurement.time, wrongKeys, firstMean, firstCovariance,
		secondMeasurement, "/PPP", secondMean, secondCovariance));
	BOOST_CHECK_EQUAL(
		capture.summary().failureReason, "MEASUREMENT_KEY_CHAIN_MISMATCH");
	BOOST_REQUIRE(capture.recordMeasurement(
		secondMeasurement.time, keys, firstMean, firstCovariance,
		secondMeasurement, "/PPP", secondMean, secondCovariance));
	BOOST_CHECK_EQUAL(capture.capturedEvents().size(), 2);
	BOOST_CHECK_EQUAL(capture.summary().measurements, 2);
}

BOOST_AUTO_TEST_CASE(canonical_affine_ledger_unions_before_current_projection)
{
	// Neither certificate is individually representable because both touch the
	// unavailable canonical coordinate 2.  Their exact difference eliminates it
	// and must survive in the permuted current coordinate order.
	const auto projected = zhangExactProjectCanonicalAffineLattice(
		{
			{1, 0, 1},
			{0, 1, 1},
		},
		{5, 7},
		{1, 0, -1},
		2);
	BOOST_REQUIRE(projected.valid);
	BOOST_CHECK_EQUAL(projected.canonicalInputRank, 2);
	BOOST_CHECK_EQUAL(projected.directlyMappableRank, 0);
	BOOST_CHECK_EQUAL(projected.recoveredCombinationRank, 1);
	BOOST_CHECK_EQUAL(projected.touchedRows, 2);
	BOOST_REQUIRE_EQUAL(projected.basis.size(), 1);
	BOOST_REQUIRE_EQUAL(projected.values.size(), 1);
	BOOST_CHECK(projected.basis.front() == ZhangExactVector({1, -1}));
	BOOST_CHECK_EQUAL(projected.values.front(), 2);
}

BOOST_AUTO_TEST_CASE(
	raw_square_root_selected_linear_marginal_matches_full_rank_deficient_projection)
{
	Vector3d mean(0.4, -0.7, 1.1);
	Matrix3d covariance;
	covariance <<
		0.9, 0.2, -0.1,
		0.2, 1.3, 0.15,
		-0.1, 0.15, 0.7;
	ZhangIncrementalRawSquareRoot incremental;
	BOOST_REQUIRE(incremental.initialise(mean, covariance));
	MatrixXd design(2, 3);
	design << 1, -0.2, 0.1, 0.3, 1, -0.4;
	Matrix2d measurementCovariance;
	measurementCovariance << 0.12, 0.02, 0.02, 0.18;
	Vector2d observation(0.2, -0.5);
	BOOST_REQUIRE(incremental.addAcceptedMeasurement(
		design, measurementCovariance, observation));

	VectorXd fullMean;
	MatrixXd fullCovariance;
	BOOST_REQUIRE(incremental.currentMarginal(fullMean, fullCovariance));
	MatrixXd targetRows(4, 3);
	targetRows << 1, -1, 0,
		2, -2, 0,
		0, 1, 1,
		0, 0, 0;
	Vector4d offsets(0.3, 0.6, -1.2, 2.25);
	const auto selected = incremental.linearMarginal(targetRows, offsets);
	BOOST_REQUIRE_MESSAGE(selected.valid, selected.failureReason);
	BOOST_CHECK_EQUAL(selected.targetRank, 2);
	BOOST_CHECK_SMALL(
		(selected.mean - (offsets + targetRows * fullMean)).norm(), 1e-11);
	const MatrixXd expectedCovariance =
		targetRows * fullCovariance * targetRows.transpose();
	BOOST_CHECK_SMALL(
		(selected.covariance - expectedCovariance).norm(), 1e-11);
	BOOST_CHECK_SMALL((selected.covariance.row(1)
		- 2 * selected.covariance.row(0)).norm(), 1e-11);
	BOOST_CHECK_SMALL(selected.mean(3) - offsets(3), 1e-12);
	BOOST_CHECK_SMALL(selected.covariance.row(3).norm(), 1e-12);
	BOOST_CHECK_SMALL(selected.covariance.col(3).norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(
	raw_square_root_qr_permutation_preserves_exact_condition_and_future_update)
{
	Vector3d mean(0.35, -0.8, 1.25);
	Matrix3d covariance;
	covariance <<
		1.2, 0.25, -0.15,
		0.25, 0.7, 0.18,
		-0.15, 0.18, 1.5;
	ZhangIncrementalRawSquareRoot incremental;
	BOOST_REQUIRE(incremental.initialise(mean, covariance));

	MatrixXd exactRow(1, 3);
	exactRow << 1, 3, -2;
	VectorXd exactValue = VectorXd::Constant(1, 0.45);
	const MatrixXd exactInnovationCovariance =
		exactRow * covariance * exactRow.transpose();
	const MatrixXd exactGain = covariance * exactRow.transpose()
		* exactInnovationCovariance.inverse();
	mean += exactGain * (exactValue - exactRow * mean);
	covariance -= exactGain * exactRow * covariance;
	covariance = 0.5 * (covariance + covariance.transpose());
	BOOST_REQUIRE(incremental.applyExactConstraint(exactRow, exactValue));

	VectorXd actualMean;
	MatrixXd actualCovariance;
	BOOST_REQUIRE(incremental.currentMarginal(actualMean, actualCovariance));
	BOOST_CHECK_SMALL((actualMean - mean).norm(), 1e-10);
	BOOST_CHECK_SMALL((actualCovariance - covariance).norm(), 1e-10);
	BOOST_CHECK_SMALL((exactRow * actualMean - exactValue).norm(), 1e-11);
	BOOST_CHECK_SMALL((exactRow * actualCovariance).norm(), 1e-11);

	MatrixXd design(2, 3);
	design << 0.1, 4, -0.3, -2, 0.25, 1;
	Matrix2d measurementCovariance;
	measurementCovariance << 0.08, 0.01, 0.01, 0.14;
	Vector2d observation(-0.7, 1.4);
	const Matrix2d innovationCovariance = design * covariance
		* design.transpose() + measurementCovariance;
	const MatrixXd gain = covariance * design.transpose()
		* innovationCovariance.inverse();
	mean += gain * (observation - design * mean);
	covariance -= gain * design * covariance;
	covariance = 0.5 * (covariance + covariance.transpose());
	BOOST_REQUIRE(incremental.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	BOOST_REQUIRE(incremental.currentMarginal(actualMean, actualCovariance));
	BOOST_CHECK_SMALL((actualMean - mean).norm(), 1e-10);
	BOOST_CHECK_SMALL((actualCovariance - covariance).norm(), 1e-10);
	BOOST_CHECK_SMALL((exactRow * actualMean - exactValue).norm(), 1e-10);
	BOOST_CHECK_SMALL((exactRow * actualCovariance).norm(), 1e-10);

	MatrixXd targetRows(2, 3);
	targetRows << 1, -1, 0.5, 0, 2, -1;
	Vector2d offsets(0.2, -0.4);
	const auto target = incremental.linearMarginal(targetRows, offsets);
	BOOST_REQUIRE_MESSAGE(target.valid, target.failureReason);
	BOOST_CHECK_SMALL(
		(target.mean - (offsets + targetRows * mean)).norm(), 1e-10);
	BOOST_CHECK_SMALL(
		(target.covariance - targetRows * covariance
			* targetRows.transpose()).norm(),
		1e-10);
}

BOOST_AUTO_TEST_CASE(
	persistent_target_joint_difference_cancels_common_gauge_and_rejects_stale_identity)
{
	Vector3d latentMean(10.0, 0.25, -0.4);
	Matrix3d latentCovariance;
	latentCovariance <<
		4.0, 0.3, -0.2,
		0.3, 0.5, 0.1,
		-0.2, 0.1, 0.8;
	Matrix3d endpoints;
	endpoints << 1, 1, 0, 1, 0, 1, 1, 0, 0;
	ZhangPersistentRawTargetMarginal marginal;
	marginal.valid = true;
	marginal.mean = endpoints * latentMean;
	marginal.covariance = endpoints * latentCovariance * endpoints.transpose();
	marginal.identities = {"A@V0", "B@V0", "REF@V0"};
	marginal.physicalVersions = marginal.identities;

	const auto projected = zhangPersistentTargetDifferenceMarginal(
		marginal,
		{{"A@V0", "B@V0"}, {"REF@V0", "A@V0"}}, true);
	BOOST_REQUIRE_MESSAGE(projected.valid, projected.failureReason);
	BOOST_CHECK_EQUAL(projected.availableCount, 2);
	BOOST_CHECK_EQUAL(projected.stochasticRank, 2);
	MatrixXd expectedTransform(2, 3);
	expectedTransform << 0, -1, 1, 0, 1, 0;
	const Vector2d expectedMean = expectedTransform * latentMean;
	const Matrix2d expectedCovariance =
		expectedTransform * latentCovariance * expectedTransform.transpose();
	BOOST_CHECK_SMALL((projected.mean - expectedMean).norm(), 1e-12);
	BOOST_CHECK_SMALL(
		(projected.covariance - expectedCovariance).norm(), 1e-12);
	BOOST_CHECK_LT(projected.covariance(0, 0),
		marginal.covariance(0, 0) + marginal.covariance(1, 1));

	const auto stale = zhangPersistentTargetDifferenceMarginal(
		marginal, {{"A@V0", "B@V1"}}, true);
	BOOST_CHECK(!stale.valid);
	BOOST_CHECK_EQUAL(
		stale.failureReason,
		"PERSISTENT_TARGET_DIFFERENCE_IDENTITIES_MISSING");
}

BOOST_AUTO_TEST_CASE(
	persistent_target_direct_difference_matches_full_target_projection)
{
	Vector3d mean(0.6, -1.1, 0.25);
	Matrix3d covariance;
	covariance <<
		0.9, 0.15, -0.08,
		0.15, 1.4, 0.22,
		-0.08, 0.22, 0.65;
	ZhangPersistentRawTargetWindow window;
	BOOST_REQUIRE(window.initialise(mean, covariance));
	MatrixXd rows(3, 3);
	rows << 1, 1, 0, 1, 0, 1, 0, 1, -1;
	Vector3d offsets(0.2, -0.4, 1.3);
	BOOST_REQUIRE(window.bindNewTargets(
		{"A", "B", "C"}, {"A", "B", "C"}, rows, offsets, 1));

	MatrixXd design(2, 3);
	design << 1, -0.2, 0.3, -0.4, 1, 0.1;
	Matrix2d measurementCovariance;
	measurementCovariance << 0.12, 0.015, 0.015, 0.2;
	Vector2d observation(0.1, -0.35);
	BOOST_REQUIRE(window.addAcceptedMeasurement(
		design, measurementCovariance, observation));

	const std::vector<std::pair<std::string, std::string>> pairs = {
		{"A", "B"}, {"C", "A"}};
	const auto full = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(full.valid, full.failureReason);
	const auto projected = zhangPersistentTargetDifferenceMarginal(
		full, pairs, true);
	const auto direct = window.targetDifferenceMarginal(pairs, true);
	BOOST_REQUIRE_MESSAGE(projected.valid, projected.failureReason);
	BOOST_REQUIRE_MESSAGE(direct.valid, direct.failureReason);
	BOOST_CHECK_EQUAL(direct.availableCount, 2);
	BOOST_CHECK_EQUAL(direct.stochasticRank, projected.stochasticRank);
	BOOST_CHECK_SMALL((direct.mean - projected.mean).norm(), 1e-11);
	BOOST_CHECK_SMALL(
		(direct.covariance - projected.covariance).norm(), 1e-11);

	const auto partial = window.targetDifferenceMarginal(
		{{"A", "B"}, {"A", "MISSING"}}, false);
	BOOST_REQUIRE_MESSAGE(partial.valid, partial.failureReason);
	BOOST_CHECK_EQUAL(partial.availableCount, 1);
	BOOST_CHECK(partial.availableRows[0]);
	BOOST_CHECK(!partial.availableRows[1]);
	BOOST_CHECK_SMALL(partial.mean(1), 1e-15);
	BOOST_CHECK_SMALL(partial.covariance.row(1).norm(), 1e-15);
	BOOST_CHECK_SMALL(partial.covariance.col(1).norm(), 1e-15);
}

BOOST_AUTO_TEST_CASE(
	persistent_target_full_canonical_catalogue_preserves_all_28_directions)
{
	constexpr int satelliteCount = 29;
	constexpr int canonicalRank = satelliteCount - 1;
	VectorXd mean = VectorXd::LinSpaced(satelliteCount, -0.7, 1.1);
	MatrixXd seed(satelliteCount, 5);
	for (int row = 0; row < satelliteCount; row++)
	for (int column = 0; column < seed.cols(); column++)
		seed(row, column) =
			0.03 * (row + 1) * (column + 2) +
			0.01 * ((row + 3 * column) % 7);
	MatrixXd covariance = seed * seed.transpose() +
		0.2 * MatrixXd::Identity(satelliteCount, satelliteCount);
	ZhangPersistentRawTargetWindow window;
	BOOST_REQUIRE(window.initialise(mean, covariance));
	std::vector<std::string> identities;
	for (int satellite = 0; satellite < satelliteCount; satellite++)
		identities.push_back("GPS_CANONICAL_" + std::to_string(satellite));
	BOOST_REQUIRE(window.bindNewTargets(
		identities, identities,
		MatrixXd::Identity(satelliteCount, satelliteCount),
		VectorXd::Zero(satelliteCount), 1));

	std::vector<std::pair<std::string, std::string>> pairs;
	MatrixXd difference = MatrixXd::Zero(canonicalRank, satelliteCount);
	for (int satellite = 1; satellite < satelliteCount; satellite++)
	{
		pairs.emplace_back(identities.front(), identities[satellite]);
		difference(satellite - 1, 0) = -1;
		difference(satellite - 1, satellite) = 1;
	}
	const auto marginal = window.targetDifferenceMarginal(pairs, true);
	BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
	BOOST_CHECK_EQUAL(marginal.availableCount, canonicalRank);
	BOOST_CHECK_EQUAL(marginal.mean.size(), canonicalRank);
	BOOST_CHECK_EQUAL(marginal.covariance.rows(), canonicalRank);
	BOOST_CHECK_EQUAL(marginal.covariance.cols(), canonicalRank);
	BOOST_CHECK_EQUAL(marginal.stochasticRank, canonicalRank);
	BOOST_CHECK_SMALL((marginal.mean - difference * mean).norm(), 1e-11);
	BOOST_CHECK_SMALL((marginal.covariance - difference * covariance *
		difference.transpose()).norm(), 1e-10);
	BOOST_CHECK_GT(std::abs(marginal.covariance(0, 1)), 1e-6);
}

BOOST_AUTO_TEST_CASE(
	canonical_product_relation_semantic_order_ignores_representation)
{
	const SatSys reference(E_Sys::GPS, 2);
	std::vector<ZhangProductRelationRow> first(2);
	std::vector<ZhangProductRelationRow> second(2);
	first[0].satellite = second[0].satellite = SatSys(E_Sys::GPS, 3);
	first[1].satellite = second[1].satellite = SatSys(E_Sys::GPS, 5);
	for (int row = 0; row < 2; row++)
		first[row].referenceSatellite = second[row].referenceSatellite = reference;
	first[0].currentCycleCoefficients = {1, -1};
	second[0].currentCycleCoefficients = {3, 0, -2};
	BOOST_CHECK(zhangProductRelationSemanticOrderingMatches(first, second));

	std::swap(second[0], second[1]);
	BOOST_CHECK(!zhangProductRelationSemanticOrderingMatches(first, second));
	std::swap(second[0], second[1]);
	second[1].referenceSatellite = SatSys(E_Sys::GPS, 7);
	BOOST_CHECK(!zhangProductRelationSemanticOrderingMatches(first, second));
}

BOOST_AUTO_TEST_CASE(
	persistent_canonical_snapshot_never_crosses_phase_segment)
{
	BOOST_CHECK(zhangPersistentProductSnapshotPhaseCompatible(
		"SNAPSHOT", "G03:L1C:SEG7", "G03:L1C:SEG7"));
	BOOST_CHECK(!zhangPersistentProductSnapshotPhaseCompatible(
		"SNAPSHOT", "G03:L1C:SEG7", "G03:L1C:SEG8"));
	BOOST_CHECK(!zhangPersistentProductSnapshotPhaseCompatible(
		"SNAPSHOT", "G03:L1C:SEG7", "UNRESOLVED"));
	BOOST_CHECK(!zhangPersistentProductSnapshotPhaseCompatible(
		"", "G03:L1C:SEG7", "G03:L1C:SEG7"));
}

BOOST_AUTO_TEST_CASE(
	persistent_target_exact_constraint_survives_process_and_s_basis_change)
{
	Vector2d stateMean(0.4, -0.7);
	Matrix2d stateCovariance;
	stateCovariance << 0.8, 0.1, 0.1, 1.2;
	ZhangIncrementalRawSquareRoot window;
	BOOST_REQUIRE(window.initialise(stateMean, stateCovariance));

	// Augment a persistent target a=x0-x1+2 without adding stochastic rank.
	MatrixXd augment = MatrixXd::Zero(3, 2);
	augment.topRows(2) = Matrix2d::Identity();
	augment.row(2) << 1, -1;
	Vector3d translation(0, 0, 2);
	BOOST_REQUIRE(window.applyExactCoordinateTransform(augment, translation));
	VectorXd mean;
	MatrixXd covariance;
	BOOST_REQUIRE(window.currentMarginal(mean, covariance));
	Vector3d relation(-1, 1, 1);
	BOOST_CHECK_SMALL(relation.dot(mean) - 2, 1e-12);
	BOOST_CHECK_SMALL(
		(relation.transpose() * covariance * relation)(0, 0), 1e-12);

	// State process noise would let the current state functional drift, while
	// the persistent target coordinate itself is carried with zero process
	// noise.  The next exact physical-functional factor reconnects them.
	Matrix3d transition = Matrix3d::Identity();
	Matrix3d processCovariance = Matrix3d::Zero();
	processCovariance(0, 0) = 0.1;
	BOOST_REQUIRE(window.advance(transition, processCovariance));
	BOOST_REQUIRE(window.currentMarginal(mean, covariance));
	BOOST_CHECK_GT(
		(relation.transpose() * covariance * relation)(0, 0), 0.09);
	BOOST_REQUIRE(window.applyExactConstraint(
		relation.transpose(), VectorXd::Constant(1, 2)));
	BOOST_REQUIRE(window.currentMarginal(mean, covariance));
	BOOST_CHECK_SMALL(relation.dot(mean) - 2, 1e-11);
	BOOST_CHECK_SMALL(
		std::abs((relation.transpose() * covariance * relation)(0, 0)),
		1e-11);

	// Pure S-basis coordinate change: x0'=x0+x1, x1'=x1, a'=a.
	Matrix3d basisTransform = Matrix3d::Identity();
	basisTransform(0, 1) = 1;
	BOOST_REQUIRE(window.applyExactCoordinateTransform(basisTransform));
	Vector3d transformedRelation(-1, 2, 1);
	BOOST_REQUIRE(window.applyExactConstraint(
		transformedRelation.transpose(), VectorXd::Constant(1, 2)));
	BOOST_REQUIRE(window.currentMarginal(mean, covariance));
	BOOST_CHECK_SMALL(transformedRelation.dot(mean) - 2, 1e-11);
	BOOST_CHECK_SMALL(
		std::abs((transformedRelation.transpose()
			* covariance * transformedRelation)(0, 0)), 1e-11);
	BOOST_CHECK_EQUAL(window.summary().exactConstraintsApplied, 2);
}

BOOST_AUTO_TEST_CASE(
	persistent_raw_target_variable_is_constant_across_s_basis_and_versions_reset)
{
	Vector2d mean(0.4, -0.7);
	Matrix2d covariance;
	covariance << 0.8, 0.1, 0.1, 1.2;
	ZhangPersistentRawTargetWindow window;
	BOOST_REQUIRE(window.initialise(mean, covariance));
	Matrix2d design;
	design << 1, 0.2, -0.3, 1;
	Matrix2d measurementCovariance = 0.1 * Matrix2d::Identity();
	Vector2d observation(0.1, -0.2);
	BOOST_REQUIRE(window.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	Vector2d targetRow(1, -1);
	BOOST_REQUIRE(window.bindTarget(
		"GPS:WL:G01->G03", "G01:0->G03:0", targetRow, 2, 1));
	const auto initial = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(initial.valid, initial.failureReason);
	BOOST_REQUIRE_EQUAL(initial.targetCount, 1);

	Matrix2d transition = Matrix2d::Identity();
	Matrix2d processCovariance = Matrix2d::Zero();
	processCovariance(0, 0) = 0.05;
	BOOST_REQUIRE(window.advance(transition, processCovariance));
	BOOST_REQUIRE(window.addAcceptedMeasurement(
		design, measurementCovariance, observation));
	BOOST_REQUIRE(window.bindTarget(
		"GPS:WL:G01->G03", "G01:0->G03:0", targetRow, 2, 2));
	const auto constrained = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(constrained.valid, constrained.failureReason);

	Matrix2d basisTransform;
	basisTransform << 1, 1, 0, 1;
	BOOST_REQUIRE(window.applyExactCoordinateTransform(basisTransform));
	Vector2d transformedTargetRow(1, -2);
	BOOST_REQUIRE(window.bindTarget(
		"GPS:WL:G01->G03", "G01:0->G03:0",
		transformedTargetRow, 2, 3));
	const auto transformed = window.targetMarginal();
	BOOST_REQUIRE_MESSAGE(transformed.valid, transformed.failureReason);
	BOOST_CHECK_SMALL(
		transformed.mean(0) - constrained.mean(0), 1e-11);
	BOOST_CHECK_SMALL(
		transformed.covariance(0, 0) - constrained.covariance(0, 0),
		1e-11);
	BOOST_CHECK_EQUAL(transformed.identities[0], "GPS:WL:G01->G03");
	BOOST_CHECK_EQUAL(transformed.physicalVersions[0], "G01:0->G03:0");

	BOOST_CHECK(!window.bindTarget(
		"GPS:WL:G01->G03", "G01:1->G03:0",
		transformedTargetRow, 2, 4));
	BOOST_CHECK_EQUAL(
		window.lastFailureReason(),
		"PERSISTENT_RAW_TARGET_PHYSICAL_VERSION_CHANGED");
}

BOOST_AUTO_TEST_CASE(incremental_fixed_lag_eliminates_epoch_local_nuisance_immediately)
{
	VectorXd boundaryMean = VectorXd::Constant(1, 0.2);
	MatrixXd boundaryCovariance = MatrixXd::Constant(1, 1, 0.8);
	MatrixXd separatorDesign(3, 1);
	separatorDesign << 1.0, 0.7, -0.4;
	MatrixXd localDesign(3, 1);
	localDesign << 1.0, -1.0, 0.5;
	MatrixXd measurementCovariance = MatrixXd::Zero(3, 3);
	measurementCovariance.diagonal() << 0.09, 0.16, 0.25;
	VectorXd observation(3);
	observation << 0.31, -0.08, 0.17;

	ZhangIncrementalFixedLagSquareRoot incremental(2);
	BOOST_REQUIRE(incremental.initialise(boundaryMean, boundaryCovariance));
	BOOST_REQUIRE(incremental.addLatestMeasurement(
		separatorDesign, localDesign, measurementCovariance, observation));
	const auto summary = incremental.summary();
	BOOST_REQUIRE(summary.valid);
	BOOST_CHECK_EQUAL(summary.storedColumns, 1);
	BOOST_CHECK_EQUAL(summary.storedRows, 1);

	const double priorWeight = 1 / std::sqrt(boundaryCovariance(0, 0));
	MatrixXd denseFactor = MatrixXd::Zero(4, 2); // [local, separator]
	denseFactor(0, 1) = priorWeight;
	VectorXd denseRhs = VectorXd::Zero(4);
	denseRhs(0) = priorWeight * boundaryMean(0);
	const MatrixXd measurementWeight = measurementCovariance.llt()
		.matrixL().solve(MatrixXd::Identity(3, 3));
	denseFactor.bottomLeftCorner(3, 1) = measurementWeight * localDesign;
	denseFactor.bottomRightCorner(3, 1) = measurementWeight * separatorDesign;
	denseRhs.tail(3) = measurementWeight * observation;
	const auto dense = zhangMarginaliseSquareRootFactors(
		denseFactor.sparseView(), denseRhs, 1);
	const auto current = incremental.latestMarginal();
	BOOST_REQUIRE_MESSAGE(dense.valid, dense.failureReason);
	BOOST_REQUIRE_MESSAGE(current.valid, current.failureReason);
	BOOST_CHECK_SMALL((current.mean - dense.mean).norm(), 1e-11);
	BOOST_CHECK_SMALL(
		(current.covariance - dense.covariance).norm()
			/ dense.covariance.norm(),
		1e-11);
}

BOOST_AUTO_TEST_CASE(incremental_exact_s_basis_change_preserves_physical_posterior)
{
	Vector2d mean(0.4, -0.3);
	Matrix2d covariance;
	covariance << 0.7, 0.1, 0.1, 0.5;
	ZhangIncrementalFixedLagSquareRoot incremental(2);
	BOOST_REQUIRE(incremental.initialise(mean, covariance));
	const auto before = incremental.latestMarginal();
	BOOST_REQUIRE_MESSAGE(before.valid, before.failureReason);
	Matrix2d transform;
	transform << 1, 1, 0, 1;
	Vector2d translation(3, -2);
	BOOST_REQUIRE(incremental.applyExactLatestCoordinateTransform(
		transform, translation));
	const auto after = incremental.latestMarginal();
	BOOST_REQUIRE_MESSAGE(after.valid, after.failureReason);
	BOOST_CHECK_SMALL(
		(after.mean - (transform * before.mean + translation)).norm(), 1e-11);
	const Matrix2d expectedCovariance = transform * before.covariance
		* transform.transpose();
	BOOST_CHECK_SMALL(
		(after.covariance - expectedCovariance).norm()
			/ expectedCovariance.norm(),
		1e-11);
}

BOOST_AUTO_TEST_CASE(residual_domains_keep_independent_dof_and_project_integer_gauge)
{
	Vector2d priorMean(0.1, -0.2);
	Matrix2d priorCovariance;
	priorCovariance << 0.8, 0.1, 0.1, 0.6;
	Matrix2d design;
	design << 1.0, 0.2, -0.3, 1.0;
	Matrix2d measurementCovariance;
	measurementCovariance << 0.09, 0.01, 0.01, 0.16;
	Vector2d observation(0.24, -0.17);
	const auto prefit = zhangPrefitInnovationStatistic(
		observation, design, measurementCovariance,
		priorMean, priorCovariance);
	BOOST_REQUIRE_MESSAGE(prefit.valid, prefit.failureReason);
	BOOST_CHECK_EQUAL(prefit.dof, 2);

	MatrixXd batchDesign(5, 2);
	batchDesign << 1, 0, 0, 1, 1, 1, 0.5, -0.2, -0.3, 0.7;
	VectorXd batchRhs(5);
	batchRhs << 0.2, -0.1, 0.13, 0.09, -0.04;
	const auto orthogonal = zhangBatchOrthogonalResidualStatistic(
		batchDesign, batchRhs);
	BOOST_REQUIRE_MESSAGE(orthogonal.valid, orthogonal.failureReason);
	BOOST_CHECK_EQUAL(orthogonal.dof, 3);

	const Matrix2d innovation = design * priorCovariance
		* design.transpose() + measurementCovariance;
	const Matrix2d gain = priorCovariance * design.transpose()
		* innovation.inverse();
	const Vector2d posteriorMean = priorMean
		+ gain * (observation - design * priorMean);
	const Matrix2d posteriorCovariance = priorCovariance
		- gain * design * priorCovariance;
	Vector2d heldOutObservation(0.18, -0.11);
	const auto heldOut = zhangHeldOutPredictionStatistic(
		heldOutObservation, design, measurementCovariance,
		posteriorMean, posteriorCovariance);
	BOOST_REQUIRE_MESSAGE(heldOut.valid, heldOut.failureReason);
	BOOST_CHECK_EQUAL(heldOut.dof, 2);

	Vector2d targetMean(12.18, -3.81);
	Matrix2d targetCovariance;
	targetCovariance << 0.04, 0.01, 0.01, 0.09;
	Vector2d candidate(12, -4);
	MatrixXd quotientDirection(2, 1);
	quotientDirection << 1, 1;
	const auto integerDistance = zhangTargetToIntegerStatistic(
		targetMean, targetCovariance, candidate, quotientDirection);
	BOOST_REQUIRE_MESSAGE(
		integerDistance.valid, integerDistance.failureReason);
	BOOST_CHECK_EQUAL(integerDistance.removedGaugeRank, 1);
	BOOST_CHECK_EQUAL(integerDistance.dof, 1);
	const auto gaugeShifted = zhangTargetToIntegerStatistic(
		targetMean + 7 * quotientDirection.col(0),
		targetCovariance, candidate, quotientDirection);
	BOOST_REQUIRE(gaugeShifted.valid);
	BOOST_CHECK_SMALL(
		integerDistance.squaredNorm - gaugeShifted.squaredNorm, 1e-12);
	BOOST_TEST_MESSAGE(
		"residual-dof prefit=" << prefit.dof
		<< " batch_orthogonal=" << orthogonal.dof
		<< " held_out=" << heldOut.dof
		<< " integer_quotient=" << integerDistance.dof
		<< " target_distance=" << integerDistance.squaredNorm);
}

BOOST_AUTO_TEST_CASE(generic_primitive_integer_targets_do_not_require_wide_lane)
{
	MatrixXd baseRows = MatrixXd::Identity(2, 2);
	Vector2d baseMean(12.08, -3.96);
	Matrix2d baseCovariance;
	baseCovariance << 0.018, 0.006, 0.006, 0.025;

	const ZhangIntegerMatrix direct = zhangDirectJointIntegerTransform(2);
	const ZhangIntegerMatrix wideLaneL1 = zhangWideLaneL1IntegerTransform();
	const auto directTargets = zhangTransformIntegerTargets(
		baseRows, baseMean, baseCovariance, direct);
	const auto wlTargets = zhangTransformIntegerTargets(
		baseRows, baseMean, baseCovariance, wideLaneL1);
	BOOST_REQUIRE_MESSAGE(directTargets.valid, directTargets.failureReason);
	BOOST_REQUIRE_MESSAGE(wlTargets.valid, wlTargets.failureReason);
	BOOST_CHECK(directTargets.audit.unimodular);
	BOOST_CHECK(wlTargets.audit.unimodular);
	BOOST_CHECK_SMALL(wlTargets.mean(0) - (baseMean(0) - baseMean(1)), 1e-15);
	BOOST_CHECK_SMALL(wlTargets.mean(1) - baseMean(0), 1e-15);

	// Representative unimodular decorrelation returned by an integer solver.
	ZhangIntegerMatrix lambdaTransform(2, 2);
	lambdaTransform << 1, 0, -2, 1;
	const auto lambdaTargets = zhangTransformIntegerTargets(
		baseRows, baseMean, baseCovariance, lambdaTransform);
	BOOST_REQUIRE_MESSAGE(lambdaTargets.valid, lambdaTargets.failureReason);
	BOOST_CHECK(lambdaTargets.audit.unimodular);

	ZhangIntegerMatrix nonPrimitive = ZhangIntegerMatrix::Zero(2, 2);
	nonPrimitive.diagonal() << 2, 1;
	const auto rejected = zhangAuditPrimitiveIntegerTransform(nonPrimitive);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK(!rejected.primitive);

	// All complete primitive coordinates recover the same base pair exactly.
	const Matrix2d wlInverse = wideLaneL1.cast<double>().transpose().inverse();
	BOOST_CHECK_SMALL(
		(wlInverse * wlTargets.mean - directTargets.mean).norm(), 1e-14);
	BOOST_CHECK_SMALL(
		(wlInverse * wlTargets.covariance * wlInverse.transpose()
			- directTargets.covariance).norm(),
		1e-14);
}

BOOST_AUTO_TEST_CASE(integer_quotient_and_wide_lane_coordinates_are_primitive)
{
	const std::vector<std::string> identities = {
		"GPS:K1_L1C:G01:G03|arc-a",
		"GPS:K2_L2W:G01:G03|arc-b",
		"GPS:K1_L1C:G01:G02|arc-c",
		"GPS:K2_L2W:G01:G02|arc-d",
		"GPS:K1_L1C:G01:G05|arc-e",
		"GPS:K2_L2W:G01:G05|arc-f"};
	const std::vector<std::string> gauges = {
		"GPS:K1", "GPS:K2", "GPS:K1", "GPS:K2", "GPS:K1", "GPS:K2"};
	const std::vector<bool> absolute(6, false);
	VectorXd mean(6);
	mean << 9.4, 0.6, -177.6, -24.1, -40.7, 40.0;
	MatrixXd covariance = MatrixXd::Identity(6, 6);
	const auto quotient = zhangBuildIntegerQuotientCoordinates(
		identities, gauges, absolute, mean, covariance);
	BOOST_REQUIRE_MESSAGE(quotient.valid, quotient.failureReason);
	BOOST_CHECK_EQUAL(quotient.transform.rows(), 6);
	BOOST_CHECK_EQUAL(quotient.transform.cols(), 4);
	BOOST_CHECK(zhangAuditPrimitiveIntegerTransform(quotient.transform).valid);
	BOOST_CHECK_SMALL(quotient.mean(0) - (mean(2) - mean(0)), 1e-14);
	BOOST_CHECK_SMALL(quotient.mean(1) - (mean(4) - mean(0)), 1e-14);
	BOOST_CHECK_SMALL(quotient.mean(2) - (mean(3) - mean(1)), 1e-14);
	BOOST_CHECK_SMALL(quotient.mean(3) - (mean(5) - mean(1)), 1e-14);

	const auto wideLane = zhangBuildWideLaneL1BlockCoordinates(quotient);
	BOOST_REQUIRE_MESSAGE(wideLane.valid, wideLane.failureReason);
	const auto audit = zhangAuditPrimitiveIntegerTransform(wideLane.transform);
	BOOST_CHECK(audit.valid);
	BOOST_CHECK(audit.unimodular);
	const MatrixXd transform = wideLane.transform.cast<double>();
	const VectorXd transformed = transform.transpose() * quotient.mean;
	BOOST_CHECK_SMALL(transformed(0) - (quotient.mean(0) - quotient.mean(2)), 1e-14);
	BOOST_CHECK_SMALL(transformed(1) - quotient.mean(0), 1e-14);
}

BOOST_AUTO_TEST_CASE(production_canonical_separator_identities_pair_l1c_l2w)
{
	const std::vector<std::string> identities = {
		"GPS:L1C:G01->G03|datum=GPS:L1C:G01->G03:V0|phase=G01:0->G03:0",
		"GPS:L2W:G01->G03|datum=GPS:L2W:G01->G03:V0|phase=G01:0->G03:0",
		"GPS:L1C:G01->G02|datum=GPS:L1C:G01->G02:V0|phase=G01:0->G02:0",
		"GPS:L2W:G01->G02|datum=GPS:L2W:G01->G02:V0|phase=G01:0->G02:0",
		"GPS:L1C:G01->G05|datum=GPS:L1C:G01->G05:V0|phase=G01:0->G05:0",
		"GPS:L2W:G01->G05|datum=GPS:L2W:G01->G05:V0|phase=G01:0->G05:0"};
	const std::vector<std::string> gauges = {
		"GPS:K1", "GPS:K2", "GPS:K1", "GPS:K2", "GPS:K1", "GPS:K2"};
	const std::vector<bool> absolute(6, false);
	VectorXd mean(6);
	mean << 9.4, 0.6, -177.6, -24.1, -40.7, 40.0;
	const MatrixXd covariance = MatrixXd::Identity(6, 6);
	const auto quotient = zhangBuildIntegerQuotientCoordinates(
		identities, gauges, absolute, mean, covariance);
	BOOST_REQUIRE_MESSAGE(quotient.valid, quotient.failureReason);
	BOOST_REQUIRE_EQUAL(quotient.relations.size(), 4);
	BOOST_CHECK_EQUAL(quotient.relations[0], "G03->G02");
	BOOST_CHECK_EQUAL(quotient.relations[1], "G03->G05");
	BOOST_CHECK_EQUAL(quotient.relations[2], "G03->G02");
	BOOST_CHECK_EQUAL(quotient.relations[3], "G03->G05");
	const auto wideLane = zhangBuildWideLaneL1BlockCoordinates(quotient);
	BOOST_REQUIRE_MESSAGE(wideLane.valid, wideLane.failureReason);
	BOOST_CHECK(zhangAuditPrimitiveIntegerTransform(wideLane.transform).unimodular);
}

BOOST_AUTO_TEST_CASE(par_subset_selection_uses_joint_covariance)
{
	Matrix3d covariance;
	covariance <<
		0.004, 0.001, 0,
		0.001, 0.006, 0,
		0, 0, 4.0;
	double success = 0;
	const auto subset = zhangSelectParSubset(covariance, 0.99, &success);
	BOOST_REQUIRE_EQUAL(subset.size(), 2);
	BOOST_CHECK_EQUAL(subset[0], 0);
	BOOST_CHECK_EQUAL(subset[1], 1);
	BOOST_CHECK_GE(success, 0.99);
}

BOOST_AUTO_TEST_CASE(lambda_reduction_diagnostics_require_conditional_variances)
{
	Matrix2d covariance;
	covariance << 0.84, 0.24, 0.24, 1.38;
	Matrix2d transform;
	transform << 1, 1, 0, 1;
	const Matrix2d reduced =
		transform.transpose() * covariance * transform;
	Eigen::LDLT<Matrix2d> conditionalFactor(reduced);
	BOOST_REQUIRE_EQUAL(conditionalFactor.info(), Eigen::Success);
	Vector2d best(9, 73);
	Vector2d second(9, 72);
	const auto diagnostics = zhangAuditLambdaReduction(
		covariance, transform, reduced, conditionalFactor.vectorD(),
		best, second);
	BOOST_REQUIRE_MESSAGE(diagnostics.valid, diagnostics.failureReason);
	BOOST_CHECK(diagnostics.transformUnimodular);
	BOOST_CHECK(diagnostics.candidateBackTransformConsistent);
	BOOST_CHECK_SMALL(diagnostics.covarianceTransformMaximumError, 1e-12);
	BOOST_CHECK_SMALL(diagnostics.conditionalDeterminantLogError, 1e-12);
	BOOST_CHECK_SMALL(diagnostics.bestCandidateBackTransformMaximumError, 1e-12);
	BOOST_CHECK_SMALL(diagnostics.reducedCandidateIntegerMaximumError, 1e-12);
	BOOST_CHECK_CLOSE(
		diagnostics.ambiguityDilutionOfPrecision,
		std::pow(covariance.determinant(), 0.25), 1e-10);
	BOOST_CHECK_CLOSE(
		diagnostics.jointBootstrappedSuccessRate,
		diagnostics.conditionalSuccessRates.prod(), 1e-10);

	// The marginal diagonal is not the LAMBDA conditional D for correlated
	// coordinates and must fail the determinant-consistency audit.
	const auto marginalMisuse = zhangAuditLambdaReduction(
		covariance, transform, reduced, reduced.diagonal(), best, second);
	BOOST_CHECK(!marginalMisuse.valid);
	BOOST_CHECK_EQUAL(
		marginalMisuse.failureReason, "INCONSISTENT_LAMBDA_REDUCTION");
}

BOOST_AUTO_TEST_CASE(lambda_par_diagnostics_report_joint_ranks_candidates_and_closure)
{
	Vector3d floatMean(5.04, -1.97, -3.02);
	Matrix3d covariance;
	covariance <<
		0.006, 0.001, -0.0005,
		0.001, 0.009, 0.0015,
		-0.0005, 0.0015, 0.012;
	ZhangIntegerVector best(3);
	best << 5, -2, -3;
	ZhangIntegerVector second(3);
	second << 5, -2, -2;
	MatrixXd relationDesign(2, 3);
	relationDesign << 1, -1, 0, 0, 1, -1;
	MatrixXd closureDesign(1, 3);
	closureDesign << 1, 1, 1;
	const auto diagnostics = zhangEvaluateLambdaParCandidates(
		floatMean, covariance, best, second,
		3, 2, relationDesign, closureDesign, 0.99);
	BOOST_REQUIRE_MESSAGE(diagnostics.valid, diagnostics.failureReason);
	BOOST_CHECK_EQUAL(diagnostics.quotientValidRank, 3);
	BOOST_CHECK_EQUAL(diagnostics.absoluteValidRank, 2);
	BOOST_CHECK_EQUAL(diagnostics.productRelationGraphRank, 2);
	BOOST_CHECK_EQUAL(diagnostics.recoverableSatelliteCount, 3);
	BOOST_CHECK_GT(diagnostics.secondCandidateDistance,
		diagnostics.bestCandidateDistance);
	BOOST_CHECK_GT(diagnostics.distanceRatio, 1);
	BOOST_CHECK_EQUAL(diagnostics.maximumCycleClosureError, 0);
	BOOST_CHECK_GT(diagnostics.parTargetCount, 0);
	BOOST_CHECK_LE(diagnostics.parTargetCount, 3);
	BOOST_TEST_MESSAGE(
		"integer-diagnostics quotient_rank="
		<< diagnostics.quotientValidRank
		<< " absolute_rank=" << diagnostics.absoluteValidRank
		<< " graph_rank=" << diagnostics.productRelationGraphRank
		<< " best=" << diagnostics.bestCandidateDistance
		<< " second=" << diagnostics.secondCandidateDistance
		<< " ratio=" << diagnostics.distanceRatio
		<< " joint_success=" << diagnostics.jointBootstrappedSuccessRate
		<< " par_targets=" << diagnostics.parTargetCount
		<< " par_success=" << diagnostics.parBootstrappedSuccessRate
		<< " recoverable_satellites="
		<< diagnostics.recoverableSatelliteCount
		<< " closure=" << diagnostics.maximumCycleClosureError);
}

BOOST_AUTO_TEST_CASE(retained_target_whitening_projects_shared_quotient_direction)
{
	KFKey first;
	first.type = KF::SAT_CLOCK;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)};
	Vector2d priorMean = Vector2d::Zero();
	Matrix2d priorCovariance = Matrix2d::Identity();
	KFMeas measurement;
	measurement.H = Matrix2d::Identity();
	measurement.V = Vector2d(0.2, -0.1);
	measurement.R = Matrix2d::Identity();
	measurement.obsKeys = {first, second};
	const Vector2d posteriorMean = 0.5 * measurement.V;
	const Matrix2d posteriorCovariance = 0.5 * Matrix2d::Identity();
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	Vector2d firstRow(1, 0);
	Vector2d secondRow(0, 1);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K1", "arc-a", "segment-a", {{"arc-a", 0}},
		keys, firstRow, 0, posteriorMean, posteriorCovariance,
		1, "GPS:L1L2:COMPONENT-0"));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "arc-b", "segment-b", {{"arc-b", 0}},
		keys, secondRow, 0, posteriorMean, posteriorCovariance,
		1, "GPS:L1L2:COMPONENT-0"));
	const auto& block = capture.currentRetainedBlock();
	BOOST_REQUIRE_MESSAGE(block.valid, block.failureReason);
	BOOST_CHECK_EQUAL(block.targetCount, 2);
	BOOST_CHECK_EQUAL(block.informationRank, 2);
	BOOST_CHECK_EQUAL(block.projectedGaugeRank, 1);
	BOOST_CHECK_EQUAL(block.residualDof, 1);
	BOOST_CHECK_EQUAL(block.whitenedResidual.size(), 1);
	const auto originalSeparatorIdentities = block.separatorIdentities;

	// A disjoint exact raw-arc representation with the same physical phase
	// segments is a coordinate continuation, not a new separator identity.
	KFMeas secondMeasurement = measurement;
	secondMeasurement.V = Vector2d::Zero();
	const Matrix2d secondPosteriorCovariance =
		(1.0 / 3.0) * Matrix2d::Identity();
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, posteriorMean, posteriorCovariance, secondMeasurement,
		"/PPP", posteriorMean, secondPosteriorCovariance));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K1", "coordinate-c", "segment-a", {{"arc-c", 0}},
		keys, firstRow, 0, posteriorMean, secondPosteriorCovariance,
		1, "GPS:L1L2:COMPONENT-0"));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-d", "segment-b", {{"arc-d", 0}},
		keys, secondRow, 0, posteriorMean, secondPosteriorCovariance,
		1, "GPS:L1L2:COMPONENT-0"));
	const auto continuedIdentities =
		capture.currentRetainedBlock().separatorIdentities;
	BOOST_CHECK(continuedIdentities == originalSeparatorIdentities);

	// A version change on the same physical arc is a hard window boundary.  It
	// must not be appended as a new separator to the old chronology.
	KFMeas thirdMeasurement = measurement;
	thirdMeasurement.V = Vector2d::Zero();
	const Matrix2d thirdPosteriorCovariance =
		0.25 * Matrix2d::Identity();
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, posteriorMean, secondPosteriorCovariance,
		thirdMeasurement, "/PPP", posteriorMean, thirdPosteriorCovariance));
	BOOST_CHECK(!capture.recordPhysicalTarget(
		GTime(), "K1", "coordinate-c", "segment-a", {{"arc-c", 1}},
		keys, firstRow, 0, posteriorMean, thirdPosteriorCovariance,
		1, "GPS:L1L2:COMPONENT-0"));
	BOOST_CHECK_EQUAL(
		capture.lastTargetReason(),
		"PERSISTENT_RAW_TARGET_PHYSICAL_VERSION_CHANGED");
	BOOST_CHECK_EQUAL(capture.capturedPhysicalTargets().size(), 4);
	capture.resetForPhysicalArcChange();
	BOOST_CHECK(!capture.summary().valid);
}

BOOST_AUTO_TEST_CASE(factor_evidence_window_tracks_every_accepted_measurement_block)
{
	KFKey first;
	first.type = KF::SAT_CLOCK;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)};
	const Vector2d priorMean = Vector2d::Zero();
	const Matrix2d priorCovariance = Matrix2d::Identity();
	KFMeas firstMeasurement;
	firstMeasurement.H = Matrix2d::Identity();
	firstMeasurement.V = Vector2d(0.2, -0.1);
	firstMeasurement.R = Matrix2d::Identity();
	firstMeasurement.obsKeys = {first, second};
	const Vector2d firstPosteriorMean = 0.5 * firstMeasurement.V;
	const Matrix2d firstPosteriorCovariance =
		0.5 * Matrix2d::Identity();
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, firstMeasurement, "/PPP",
		firstPosteriorMean, firstPosteriorCovariance));

	const auto firstWindow = capture.factorEvidenceWindow();
	BOOST_REQUIRE_MESSAGE(firstWindow.valid, firstWindow.failureReason);
	BOOST_REQUIRE_EQUAL(firstWindow.measurementSequences.size(), 1);
	BOOST_CHECK_EQUAL(firstWindow.measurementSequences[0], 0);
	BOOST_CHECK_EQUAL(firstWindow.firstSequence, 0);
	BOOST_CHECK_EQUAL(firstWindow.lastSequence, 0);
	BOOST_CHECK_EQUAL(firstWindow.measurementRows, 2);
	const auto windowIdentity = firstWindow.identity;

	KFMeas secondMeasurement = firstMeasurement;
	secondMeasurement.V = Vector2d::Zero();
	const Matrix2d secondPosteriorCovariance =
		(1.0 / 3.0) * Matrix2d::Identity();
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, firstPosteriorMean, firstPosteriorCovariance,
		secondMeasurement, "/PPP", firstPosteriorMean,
		secondPosteriorCovariance));

	const auto secondWindow = capture.factorEvidenceWindow();
	BOOST_REQUIRE_MESSAGE(secondWindow.valid, secondWindow.failureReason);
	BOOST_REQUIRE_EQUAL(secondWindow.measurementSequences.size(), 2);
	BOOST_CHECK_EQUAL(secondWindow.measurementSequences[0], 0);
	BOOST_CHECK_EQUAL(secondWindow.measurementSequences[1], 1);
	BOOST_CHECK_EQUAL(secondWindow.firstSequence, 0);
	BOOST_CHECK_EQUAL(secondWindow.lastSequence, 1);
	BOOST_CHECK_EQUAL(secondWindow.measurementRows, 4);
	BOOST_CHECK_EQUAL(secondWindow.identity, windowIdentity);
}

BOOST_AUTO_TEST_CASE(raw_window_never_promotes_unresolved_peer_to_absolute_datum)
{
	KFKey first;
	first.type = KF::SAT_CLOCK;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)};
	const Vector2d priorMean = Vector2d::Zero();
	const Matrix2d priorCovariance = Matrix2d::Identity();
	KFMeas measurement;
	measurement.H = Matrix2d::Identity();
	measurement.V = Vector2d(0.1, -0.2);
	measurement.R = Matrix2d::Identity();
	measurement.obsKeys = {first, second};
	const Vector2d posteriorMean = 0.5 * measurement.V;
	const Matrix2d posteriorCovariance = 0.5 * Matrix2d::Identity();
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	const Vector2d firstRow(1, 0);
	const Vector2d secondRow(0, 1);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "GPS:K1_L1C:G01:G02", "arc-a", "G01:0->G02:0",
		{{"arc-a", 0}}, keys, firstRow, 0,
		posteriorMean, posteriorCovariance,
		0, "", "GPS:L1C:G01->G02", "GPS:L1C:G01->G02:V0", 0));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "GPS:K1_L1C:G01:G03", "arc-b", "G01:0->G03:0",
		{{"arc-b", 0}}, keys, secondRow, 0,
		posteriorMean, posteriorCovariance,
		1, "GPS:K1_L1C:CANONICAL", "GPS:L1C:G01->G03",
		"GPS:L1C:G01->G03:V0", 0));
	const auto marginal = capture.currentRawSquareRootTargetMarginal();
	BOOST_REQUIRE_MESSAGE(marginal.valid, marginal.failureReason);
	BOOST_REQUIRE_EQUAL(marginal.absoluteValidity.size(), 2);
	BOOST_CHECK(marginal.absoluteValidity[0]);
	BOOST_CHECK(!marginal.absoluteValidity[1]);
	BOOST_CHECK_EQUAL(marginal.absoluteValidRank, 1);
	BOOST_CHECK_EQUAL(marginal.unresolvedGaugeRank, 1);
	BOOST_CHECK_EQUAL(marginal.quotientValidRank, 1);
}

BOOST_AUTO_TEST_CASE(single_unresolved_target_is_valid_zero_dof_block)
{
	KFKey key;
	key.type = KF::PHASE_BIAS;
	key.Sat = SatSys(E_Sys::GPS, 1);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(key)};
	KFMeas measurement;
	measurement.H = MatrixXd::Identity(1, 1);
	measurement.V = VectorXd::Constant(1, 0.1);
	measurement.R = MatrixXd::Identity(1, 1);
	measurement.obsKeys = {key};
	const VectorXd priorMean = VectorXd::Zero(1);
	const MatrixXd priorCovariance = MatrixXd::Identity(1, 1);
	const VectorXd posteriorMean = VectorXd::Constant(1, 0.05);
	const MatrixXd posteriorCovariance = MatrixXd::Constant(1, 1, 0.5);
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K1", "arc-a", "G01:0->G02:0", {{"arc-a", 0}},
		keys, VectorXd::Ones(1), 0, posteriorMean, posteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G02",
		"GPS:L1C:G01->G02:V0", 0));
	const auto& block = capture.currentRetainedBlock();
	BOOST_REQUIRE(block.valid);
	BOOST_CHECK_EQUAL(block.targetCount, 1);
	BOOST_CHECK_EQUAL(block.informationRank, 1);
	BOOST_CHECK_EQUAL(block.residualDof, 0);
	BOOST_CHECK_EQUAL(block.projectedGaugeRank, 1);
	BOOST_CHECK(!block.likelihoodValid);
	BOOST_CHECK_SMALL(block.whitenedSquaredNorm, 1e-15);
}

BOOST_AUTO_TEST_CASE(persistent_canonical_functional_survives_temporary_target_loss)
{
	KFKey first;
	first.type = KF::PHASE_BIAS;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> keys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)};
	KFMeas measurement;
	measurement.H = Matrix2d::Identity();
	measurement.V = Vector2d(0.2, -0.1);
	measurement.R = Matrix2d::Identity();
	measurement.obsKeys = {first, second};
	const Vector2d priorMean = Vector2d::Zero();
	const Matrix2d priorCovariance = Matrix2d::Identity();
	const Vector2d posteriorMean = 0.5 * measurement.V;
	const Matrix2d posteriorCovariance = 0.5 * Matrix2d::Identity();
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	const Vector2d firstRow(1, 0);
	const Vector2d secondRow(0, 1);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K1", "coordinate-a", "G01:0->G02:0",
		{{"arc-a", 0}}, keys, firstRow, 0,
		posteriorMean, posteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G02",
		"GPS:L1C:G01->G02:V0", 0));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-b", "G01:0->G03:0",
		{{"arc-b", 0}}, keys, secondRow, 0,
		posteriorMean, posteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G03",
		"GPS:L1C:G01->G03:V0", 0));

	SparseMatrix<double> transform(2, 2);
	transform.insert(0, 0) = 1;
	transform.insert(0, 1) = 1;
	transform.insert(1, 1) = 1;
	BOOST_REQUIRE(capture.recordCoordinateTransform(
		GTime(), keys, keys, transform, "synthetic tree exchange"));
	const Matrix2d denseTransform = MatrixXd(transform);
	const Vector2d transformedPrior = denseTransform * posteriorMean;
	const Matrix2d transformedPriorCovariance = denseTransform
		* posteriorCovariance * denseTransform.transpose();
	KFMeas secondMeasurement = measurement;
	secondMeasurement.V = Vector2d(0.01, -0.02);
	const Matrix2d innovationCovariance = transformedPriorCovariance
		+ secondMeasurement.R;
	const Matrix2d gain = transformedPriorCovariance
		* innovationCovariance.inverse();
	const Vector2d secondPosteriorMean = transformedPrior
		+ gain * secondMeasurement.V;
	Matrix2d secondPosteriorCovariance = transformedPriorCovariance
		- gain * transformedPriorCovariance;
	secondPosteriorCovariance = 0.5
		* (secondPosteriorCovariance + secondPosteriorCovariance.transpose());
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), keys, transformedPrior, transformedPriorCovariance,
		secondMeasurement, "/PPP", secondPosteriorMean,
		secondPosteriorCovariance));

	// Only K1 is rebuilt in the new S-basis.  K2 must remain active through its
	// transported canonical functional; it must not be retired or replaced.
	const Vector2d transformedFirstRow(1, -1);
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K1", "coordinate-c", "G01:0->G02:0",
		{{"arc-c", 0}}, keys, transformedFirstRow, 0,
		secondPosteriorMean, secondPosteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G02",
		"GPS:L1C:G01->G02:V0", 0));
	const auto raw = capture.currentRawSquareRootTargetMarginal();
	BOOST_REQUIRE_MESSAGE(raw.valid, raw.failureReason);
	BOOST_CHECK_EQUAL(raw.requestedTargetCount, 2);
	BOOST_CHECK_EQUAL(raw.unresolvedGaugeRank, 1);
	BOOST_CHECK_EQUAL(raw.quotientValidRank, 1);
	const auto& retained = capture.currentRetainedBlock();
	BOOST_REQUIRE_EQUAL(retained.targetCount, 2);
	const auto incremental = capture.currentIncrementalTargetMarginal();
	BOOST_REQUIRE_MESSAGE(incremental.valid, incremental.failureReason);
	BOOST_CHECK_EQUAL(incremental.requestedTargetCount, 2);
	BOOST_CHECK_EQUAL(incremental.unresolvedGaugeRank, 1);
}

BOOST_AUTO_TEST_CASE(persistent_canonical_functional_rejects_unrepresentable_s_transform)
{
	KFKey first;
	first.type = KF::PHASE_BIAS;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> sourceKeys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)};
	const std::vector<ZhangCapturedStateKey> destinationKeys = {
		zhangCapturedStateKey(first)};
	KFMeas measurement;
	measurement.H = Matrix2d::Identity();
	measurement.V = Vector2d::Zero();
	measurement.R = Matrix2d::Identity();
	measurement.obsKeys = {first, second};
	const Vector2d priorMean = Vector2d::Zero();
	const Matrix2d priorCovariance = Matrix2d::Identity();
	const Vector2d posteriorMean = Vector2d::Zero();
	const Matrix2d posteriorCovariance = 0.5 * Matrix2d::Identity();
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), sourceKeys, priorMean, priorCovariance, measurement, "/PPP",
		posteriorMean, posteriorCovariance));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-b", "G01:0->G03:0",
		{{"arc-b", 0}}, sourceKeys, Vector2d(0, 1), 0,
		posteriorMean, posteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G03",
		"GPS:L1C:G01->G03:V0", 0));
	SparseMatrix<double> projection(1, 2);
	projection.insert(0, 0) = 1;
	BOOST_CHECK(!capture.recordCoordinateTransform(
		GTime(), sourceKeys, destinationKeys, projection,
		"unrepresentable tree exchange"));
	const auto summary = capture.summary();
	BOOST_CHECK_EQUAL(
		summary.failureReason,
		"PERSISTENT_FUNCTIONAL_NOT_TRANSPORTABLE_EXACT_COORDINATE_TRANSFORM");

	// The same loss is a legal boundary only when the caller has classified it
	// as a real physical-arc reinitialisation.  The old chronology is then
	// closed and the next accepted measurement establishes a fresh anchor.
	capture.resetForPhysicalArcChange();
	BOOST_CHECK(capture.summary().failureReason.empty());
	KFMeas restartedMeasurement;
	restartedMeasurement.H = MatrixXd::Identity(1, 1);
	restartedMeasurement.V = VectorXd::Zero(1);
	restartedMeasurement.R = MatrixXd::Identity(1, 1);
	restartedMeasurement.obsKeys = {first};
	const VectorXd restartedPrior = VectorXd::Zero(1);
	const MatrixXd restartedPriorCovariance = MatrixXd::Identity(1, 1);
	const VectorXd restartedPosterior = VectorXd::Zero(1);
	const MatrixXd restartedPosteriorCovariance =
		0.5 * MatrixXd::Identity(1, 1);
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), destinationKeys, restartedPrior, restartedPriorCovariance,
		restartedMeasurement, "/PPP", restartedPosterior,
		restartedPosteriorCovariance));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-c", "G01:1->G03:0",
		{{"arc-b", 1}}, destinationKeys, VectorXd::Ones(1), 0,
		restartedPosterior, restartedPosteriorCovariance,
		1, "GPS:L1C:CANONICAL:V1", "GPS:L1C:G01->G03",
		"GPS:L1C:G01->G03:V0", 0));
	BOOST_CHECK(capture.summary().valid);
}

BOOST_AUTO_TEST_CASE(persistent_quotient_is_held_when_state_transition_drops_its_row)
{
	KFKey first;
	first.type = KF::PHASE_BIAS;
	first.Sat = SatSys(E_Sys::GPS, 1);
	KFKey second = first;
	second.Sat = SatSys(E_Sys::GPS, 2);
	const std::vector<ZhangCapturedStateKey> sourceKeys = {
		zhangCapturedStateKey(first), zhangCapturedStateKey(second)};
	const std::vector<ZhangCapturedStateKey> destinationKeys = sourceKeys;
	KFMeas firstMeasurement;
	firstMeasurement.H = Matrix2d::Identity();
	firstMeasurement.V = Vector2d::Zero();
	firstMeasurement.R = Matrix2d::Identity();
	firstMeasurement.obsKeys = {first, second};
	const Vector2d sourceMean = Vector2d::Zero();
	const Matrix2d sourceCovariance = Matrix2d::Identity();
	const Matrix2d sourcePosteriorCovariance = 0.5 * Matrix2d::Identity();
	ZhangFactorCaptureBuffer capture;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), sourceKeys, sourceMean, sourceCovariance, firstMeasurement,
		"/PPP", sourceMean, sourcePosteriorCovariance));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K1", "coordinate-a", "G01:0->G02:0",
		{{"arc-a", 0}}, sourceKeys, Vector2d(1, 0), 0,
		sourceMean, sourcePosteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G02",
		"GPS:L1C:G01->G02:V0", 0));
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-b", "G01:0->G03:0",
		{{"arc-b", 0}}, sourceKeys, Vector2d(0, 1), 0,
		sourceMean, sourcePosteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G03",
		"GPS:L1C:G01->G03:V0", 0));
	SparseMatrix<double> projection(2, 2);
	projection.insert(0, 0) = 1;
	Matrix2d processCovariance = Matrix2d::Zero();
	processCovariance(1, 1) = 1;
	BOOST_REQUIRE(capture.recordTransition(
		GTime(), sourceKeys, destinationKeys, projection,
		processCovariance, "state retirement"));

	KFMeas secondMeasurement = firstMeasurement;
	const Vector2d destinationMean = Vector2d::Zero();
	Matrix2d destinationPriorCovariance = Matrix2d::Zero();
	destinationPriorCovariance(0, 0) = 0.5;
	destinationPriorCovariance(1, 1) = 1;
	Matrix2d destinationPosteriorCovariance = Matrix2d::Zero();
	destinationPosteriorCovariance(0, 0) = 1.0 / 3.0;
	destinationPosteriorCovariance(1, 1) = 0.5;
	BOOST_REQUIRE(capture.recordMeasurement(
		GTime(), destinationKeys, destinationMean, destinationPriorCovariance,
		secondMeasurement, "/PPP", destinationMean,
		destinationPosteriorCovariance));
	BOOST_CHECK(!capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-c", "G01:0->G03:0",
		{{"arc-c", 0}}, destinationKeys, Vector2d(0, 1), 0,
		destinationMean, destinationPosteriorCovariance,
		1, "GPS:L1C:CANONICAL", "GPS:L1C:G01->G03",
		"GPS:L1C:G01->G03:V0", 0));
	BOOST_CHECK_EQUAL(
		capture.lastTargetReason(),
		"PERSISTENT_QUOTIENT_FUNCTIONAL_NOT_TRANSPORTABLE");
	BOOST_CHECK(capture.summary().valid);
	const auto held = capture.currentIncrementalTargetMarginal();
	BOOST_REQUIRE_MESSAGE(held.valid, held.failureReason);
	BOOST_CHECK_EQUAL(held.requestedTargetCount, 2);
	BOOST_CHECK_EQUAL(held.quotientValidRank, 1);

	// An exact absolute datum for the same version may rebind the unavailable
	// row without changing the canonical product identity.
	BOOST_REQUIRE(capture.recordPhysicalTarget(
		GTime(), "K2", "coordinate-c", "G01:0->G03:0",
		{{"arc-c", 0}}, destinationKeys, Vector2d(0, 1), 4,
		destinationMean, destinationPosteriorCovariance,
		0, "", "GPS:L1C:G01->G03", "GPS:L1C:G01->G03:V0", 0));
}

BOOST_AUTO_TEST_CASE(incremental_target_separator_accumulates_and_retires_without_history_rows)
{
	ZhangIncrementalTargetSeparator separator;
	const std::vector<std::string> firstKeys = {"A@0", "B@0"};
	const Matrix2d design = Matrix2d::Identity();
	Matrix2d firstCovariance;
	firstCovariance << 0.04, 0.01, 0.01, 0.09;
	Vector2d firstObservation(2.1, -0.9);
	BOOST_REQUIRE(separator.addLikelihood(
		firstKeys, design, firstCovariance, firstObservation,
		{"G0", "G0"}, {false, false}));
	Matrix2d secondCovariance;
	secondCovariance << 0.03, -0.004, -0.004, 0.05;
	Vector2d secondObservation(1.95, -1.05);
	BOOST_REQUIRE(separator.addLikelihood(
		firstKeys, design, secondCovariance, secondObservation,
		{"G0", "G0"}, {false, false}));
	const auto accumulated = separator.marginal();
	BOOST_REQUIRE_MESSAGE(accumulated.valid, accumulated.failureReason);
	const Matrix2d expectedInformation = firstCovariance.inverse()
		+ secondCovariance.inverse();
	const Matrix2d expectedCovariance = expectedInformation.inverse();
	const Vector2d expectedMean = expectedCovariance
		* (firstCovariance.inverse() * firstObservation
			+ secondCovariance.inverse() * secondObservation);
	BOOST_CHECK_SMALL((accumulated.mean - expectedMean).norm(), 1e-11);
	BOOST_CHECK_SMALL(
		(accumulated.covariance - expectedCovariance).norm()
			/ expectedCovariance.norm(),
		1e-11);
	BOOST_CHECK_EQUAL(accumulated.informationRank, 2);
	BOOST_CHECK_EQUAL(accumulated.quotientValidRank, 1);
	BOOST_CHECK_EQUAL(accumulated.absoluteValidRank, 0);

	// An exact datum promotes only that canonical relation.  Other relations in
	// the quotient group retain their unresolved integer gauge.
	ZhangIncrementalTargetSeparator promoted = separator;
	BOOST_REQUIRE(promoted.addLikelihood(
		firstKeys, design, 0.02 * Matrix2d::Identity(),
		Vector2d(2.0, 2.0), {"G0", ""}, {false, true}, {0, 3}));
	const auto promotedMarginal = promoted.marginal();
	BOOST_REQUIRE_MESSAGE(promotedMarginal.valid, promotedMarginal.failureReason);
	BOOST_CHECK_EQUAL(promotedMarginal.unresolvedGaugeRank, 1);
	BOOST_CHECK_EQUAL(promotedMarginal.quotientValidRank, 1);
	BOOST_CHECK_EQUAL(promotedMarginal.absoluteValidRank, 1);
	BOOST_REQUIRE_EQUAL(promotedMarginal.coordinateOffsets.size(), 2);
	BOOST_CHECK_SMALL(promotedMarginal.coordinateOffsets[1] - 3, 1e-15);

	const std::vector<std::string> nextKeys = {"B@0", "C@1"};
	BOOST_REQUIRE(separator.addLikelihood(
		nextKeys, design, 0.02 * Matrix2d::Identity(),
		Vector2d(-1.0, 3.02), {"G0", ""}, {true, true}));
	BOOST_REQUIRE(separator.retainOnly({"B@0", "C@1"}));
	const auto retired = separator.marginal();
	BOOST_REQUIRE_MESSAGE(retired.valid, retired.failureReason);
	BOOST_CHECK_EQUAL(retired.requestedTargetCount, 2);
	BOOST_CHECK_EQUAL(retired.informationRank, 2);
	BOOST_CHECK_EQUAL(retired.quotientValidRank, 2);
	BOOST_CHECK_EQUAL(retired.absoluteValidRank, 2);
	BOOST_CHECK_LE(retired.storedRows, 2);
	BOOST_CHECK_LE(retired.storedColumns, 2);
	BOOST_CHECK_LE(retired.maximumStoredColumns, 3);
}

BOOST_AUTO_TEST_CASE(lambda_beam_single_row_leverage_matches_direct_schur_deletion)
{
	Matrix3d covariance;
	covariance <<
		0.08, 0.02, -0.01,
		0.02, 0.11,  0.03,
		-0.01, 0.03, 0.09;
	Vector3d innovation(0.22, -0.31, 0.17);
	const auto leverage = zhangConstraintNisLeverage(
		innovation, covariance);
	BOOST_REQUIRE(leverage.valid);

	for (int removed = 0; removed < 3; removed++)
	{
		vector<int> retained;
		for (int row = 0; row < 3; row++)
		{
			if (row != removed)
			{
				retained.push_back(row);
			}
		}
		VectorXd reducedInnovation = innovation(retained);
		MatrixXd reducedCovariance = covariance(retained, retained);
		const double reducedNis = reducedInnovation.dot(
			reducedCovariance.ldlt().solve(reducedInnovation));
		BOOST_CHECK_SMALL(
			leverage.nis - reducedNis
				- leverage.deletionReduction(removed),
			1e-12);
	}
}

BOOST_AUTO_TEST_CASE(lambda_beam_product_gain_and_hnf_are_basis_invariant)
{
	Matrix3d ambiguityCovariance;
	ambiguityCovariance <<
		0.12, 0.03, 0.01,
		0.03, 0.15, 0.02,
		0.01, 0.02, 0.09;
	Matrix<double, 2, 3> productCross;
	productCross <<
		0.04, -0.02, 0.01,
		0.01,  0.03, 0.02;
	Matrix<double, 2, 3> rows;
	rows << 1, 0, -1,
		0, 1, 1;
	Matrix2d unimodular;
	unimodular << 1, 2, 0, 1;
	Matrix<double, 2, 3> changedBasis = unimodular * rows;

	const double gain = zhangConstraintProductInformationGain(
		productCross, 0.8, ambiguityCovariance, rows);
	const double changedGain = zhangConstraintProductInformationGain(
		productCross, 0.8, ambiguityCovariance, changedBasis);
	BOOST_REQUIRE(std::isfinite(gain));
	BOOST_CHECK_SMALL(gain - changedGain, 1e-12);
	BOOST_CHECK_EQUAL(
		zhangIntegerRowHnfFingerprint(rows),
		zhangIntegerRowHnfFingerprint(changedBasis));
	BOOST_CHECK_EQUAL(
		zhangIntegerRowHnfCanonicalKey(rows),
		zhangIntegerRowHnfCanonicalKey(changedBasis));

	Matrix<double, 1, 3> reduced = rows.topRows(1);
	const double reducedGain = zhangConstraintProductInformationGain(
		productCross, 0.8, ambiguityCovariance, reduced);
	BOOST_CHECK_GE(gain + 1e-12, reducedGain);
	BOOST_CHECK_NE(
		zhangIntegerRowHnfFingerprint(rows),
		zhangIntegerRowHnfFingerprint(reduced));
	BOOST_CHECK_NE(
		zhangIntegerRowHnfCanonicalKey(rows),
		zhangIntegerRowHnfCanonicalKey(reduced));

	Vector2d rhs(7, -3);
	Vector2d changedRhs = unimodular * rhs;
	BOOST_CHECK_EQUAL(
		zhangIntegerAffineHnfCanonicalKey(rows, rhs),
		zhangIntegerAffineHnfCanonicalKey(changedBasis, changedRhs));
	BOOST_CHECK_EQUAL(
		zhangIntegerAffineHnfFingerprint(rows, rhs),
		zhangIntegerAffineHnfFingerprint(changedBasis, changedRhs));
	Vector2d inconsistentRhs = changedRhs;
	inconsistentRhs(0) += 1;
	BOOST_CHECK_NE(
		zhangIntegerAffineHnfCanonicalKey(rows, rhs),
		zhangIntegerAffineHnfCanonicalKey(changedBasis, inconsistentRhs));
}

BOOST_AUTO_TEST_CASE(iar_product_gain_spectrum_matches_known_real_mode_ceiling)
{
	const Matrix3d ambiguityCovariance = Matrix3d::Identity();
	Matrix3d ambiguityProductCross = Matrix3d::Zero();
	ambiguityProductCross.diagonal() << 3, 2, 1;
	const auto spectrum = zhangIarProductGainSpectrum(
		ambiguityCovariance,
		ambiguityProductCross,
		Matrix3d::Identity());
	BOOST_REQUIRE_MESSAGE(spectrum.valid, spectrum.failureReason);
	BOOST_REQUIRE_EQUAL(spectrum.ambiguityRank, 3);
	BOOST_REQUIRE_EQUAL(spectrum.eigenvaluesDescending.size(), 3);
	BOOST_CHECK_SMALL(spectrum.eigenvaluesDescending(0) - 9, 1e-12);
	BOOST_CHECK_SMALL(spectrum.eigenvaluesDescending(1) - 4, 1e-12);
	BOOST_CHECK_SMALL(spectrum.eigenvaluesDescending(2) - 1, 1e-12);
	BOOST_CHECK_SMALL(spectrum.totalWeightedGain - 14, 1e-12);
	BOOST_CHECK_SMALL(spectrum.rho(1) - 9.0 / 14, 1e-12);
	BOOST_CHECK_SMALL(spectrum.rho(2) - 13.0 / 14, 1e-12);
	BOOST_CHECK_SMALL(spectrum.rho(20) - 1, 1e-12);
	BOOST_CHECK_EQUAL(spectrum.minimumRankForRho(0.80), 2);
	BOOST_CHECK_EQUAL(spectrum.minimumRankForRho(0.95), 3);
}

BOOST_AUTO_TEST_CASE(iar_product_gain_spectrum_is_ambiguity_basis_invariant)
{
	Matrix3d ambiguityCovariance;
	ambiguityCovariance <<
		0.12, 0.03, 0.01,
		0.03, 0.15, 0.02,
		0.01, 0.02, 0.09;
	Matrix<double, 3, 2> ambiguityProductCross;
	ambiguityProductCross <<
		0.04, -0.02,
		0.01,  0.03,
		0.02,  0.01;
	Matrix2d productWeight;
	productWeight << 2.0, 0.2, 0.2, 0.7;
	Matrix3d transform;
	transform <<
		1, 2, 0,
		0, 1, 1,
		1, 0, 1;
	const auto baseline = zhangIarProductGainSpectrum(
		ambiguityCovariance, ambiguityProductCross, productWeight);
	const auto changed = zhangIarProductGainSpectrum(
		transform * ambiguityCovariance * transform.transpose(),
		transform * ambiguityProductCross,
		productWeight);
	BOOST_REQUIRE_MESSAGE(baseline.valid, baseline.failureReason);
	BOOST_REQUIRE_MESSAGE(changed.valid, changed.failureReason);
	BOOST_REQUIRE_EQUAL(
		baseline.eigenvaluesDescending.size(),
		changed.eigenvaluesDescending.size());
	BOOST_CHECK_SMALL(
		(baseline.eigenvaluesDescending -
		 changed.eigenvaluesDescending).norm(),
		1e-11);
	BOOST_CHECK_SMALL(
		baseline.totalWeightedGain - changed.totalWeightedGain,
		1e-11);
	BOOST_CHECK_SMALL(baseline.rho(2) - changed.rho(2), 1e-11);
}

BOOST_AUTO_TEST_CASE(iar_product_gain_spectrum_rejects_nullspace_cross_covariance)
{
	Matrix2d singularCovariance = Matrix2d::Zero();
	singularCovariance(0, 0) = 4;
	Vector2d validCross(2, 0);
	const auto valid = zhangIarProductGainSpectrum(
		singularCovariance,
		validCross,
		Matrix<double, 1, 1>::Identity());
	BOOST_REQUIRE_MESSAGE(valid.valid, valid.failureReason);
	BOOST_CHECK_EQUAL(valid.ambiguityRank, 1);
	BOOST_CHECK_SMALL(valid.totalWeightedGain - 1, 1e-12);

	Vector2d invalidCross(2, 0.1);
	const auto invalid = zhangIarProductGainSpectrum(
		singularCovariance,
		invalidCross,
		Matrix<double, 1, 1>::Identity());
	BOOST_CHECK(!invalid.valid);
	BOOST_CHECK_EQUAL(
		invalid.failureReason,
		"CROSS_COVARIANCE_OUTSIDE_AMBIGUITY_RANGE");
}

BOOST_AUTO_TEST_CASE(lambda_beam_bootstrap_log_failure_does_not_saturate)
{
	Vector3d conditionalVariances(1e-4, 2e-4, 5e-4);
	const double logFailure =
		zhangBootstrapLogFailure(conditionalVariances);
	BOOST_CHECK(std::isfinite(logFailure));
	BOOST_CHECK_LT(logFailure, -100);

	Vector2d moderate(0.08, 0.12);
	const double moderateLogFailure = zhangBootstrapLogFailure(moderate);
	const double directSuccess =
		std::erf(std::sqrt(1 / (8 * moderate(0)))) *
		std::erf(std::sqrt(1 / (8 * moderate(1))));
	BOOST_CHECK_SMALL(
		std::exp(moderateLogFailure) - (1 - directSuccess), 1e-14);
}

BOOST_AUTO_TEST_CASE(e25b_joint_user_integer_functional_s0_s1_closes_exactly)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	SatSys g03(E_Sys::GPS, 3);
	std::set<ZhangGraphEdge> edges = {
		{"R0", g01}, {"R0", g02}, {"R0", g03},
		{"R1", g01}, {"R1", g02}, {"R1", g03},
		{"R2", g01}, {"R2", g02}, {"R2", g03}
	};
	ZhangGraphBasis productBasis = zhangBuildSpanningTree(
		edges, "R0",
		{{"R0", g01}, {"R1", g01}, {"R1", g02},
		 {"R2", g02}, {"R2", g03}});
	BOOST_REQUIRE(productBasis.connected);
	std::map<ZhangGraphEdge, int> versions;
	for (const auto& edge : edges)
	{
		versions[edge] = 4;
	}
	auto products = zhangBuildProductIntegerFunctionals(
		productBasis, versions, g01, 7);
	BOOST_REQUIRE_EQUAL(products.size(), 3);
	for (const auto& [satellite, product] : products)
	{
		BOOST_REQUIRE(product.valid);
		BOOST_CHECK_NE(
			zhangProductIntegerFunctionalFingerprint(product), "INVALID");
	}

	auto functional = zhangBuildJointUserIntegerFunctional(
		products, g01, 5);
	BOOST_REQUIRE(functional.valid);
	auto audit = zhangAuditUserIntegerLattice(functional);
	BOOST_REQUIRE_MESSAGE(audit.valid, audit.failureReason);
	BOOST_CHECK(audit.nuisanceOrthogonal);
	BOOST_CHECK(audit.affineInteger);
	BOOST_CHECK(audit.primitiveAdmissible);
	BOOST_CHECK_SMALL(audit.maximumNuisanceCoefficient, 1e-15);
	BOOST_CHECK_SMALL(audit.maximumAffineIntegerError, 1e-15);

	std::mt19937 generator(20260808);
	std::uniform_int_distribution<int> integerDistribution(-100, 100);
	std::uniform_real_distribution<double> nuisanceDistribution(-100, 100);
	ZhangExactVector integers(functional.integerRows.front().size());
	for (auto& integer : integers)
	{
		integer = integerDistribution(generator);
	}
	VectorXd nuisanceA(5);
	VectorXd nuisanceB(5);
	for (int index = 0; index < 5; index++)
	{
		nuisanceA(index) = nuisanceDistribution(generator);
		nuisanceB(index) = nuisanceDistribution(generator);
	}
	VectorXd valueA = functional.value(integers, nuisanceA);
	VectorXd valueB = functional.value(integers, nuisanceB);
	BOOST_REQUIRE_EQUAL(valueA.size(), 2);
	BOOST_CHECK_SMALL((valueA - valueB).norm(), 1e-12);
	for (int row = 0; row < valueA.size(); row++)
	{
		BOOST_CHECK_SMALL(valueA(row) - std::round(valueA(row)), 1e-12);
	}
}

BOOST_AUTO_TEST_CASE(e25b_random_tree_exchange_preserves_fixed_product_lattice)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	SatSys g03(E_Sys::GPS, 3);
	std::set<ZhangGraphEdge> edges = {
		{"R0", g01}, {"R0", g02}, {"R0", g03},
		{"R1", g01}, {"R1", g02}, {"R1", g03},
		{"R2", g01}, {"R2", g02}, {"R2", g03}
	};
	ZhangGraphBasis currentA = zhangBuildSpanningTree(
		edges, "R0",
		{{"R0", g01}, {"R0", g02}, {"R0", g03},
		 {"R1", g01}, {"R2", g02}});
	ZhangGraphBasis currentB = zhangBuildSpanningTree(
		edges, "R0",
		{{"R0", g01}, {"R1", g01}, {"R1", g03},
		 {"R2", g02}, {"R2", g03}});
	ZhangGraphBasis fixedProduct = zhangBuildSpanningTree(
		edges, "R0",
		{{"R0", g02}, {"R1", g02}, {"R1", g03},
		 {"R2", g01}, {"R2", g03}});
	BOOST_REQUIRE(currentA.connected && currentB.connected && fixedProduct.connected);
	BOOST_REQUIRE(currentA.treeEdges != currentB.treeEdges);

	ZhangSatelliteProductTarget representationA =
		zhangBuildSatelliteProductTarget(currentA, fixedProduct, g01);
	ZhangSatelliteProductTarget representationB =
		zhangBuildSatelliteProductTarget(currentB, fixedProduct, g01);
	BOOST_REQUIRE(representationA.valid && representationB.valid);
	BOOST_CHECK(
		zhangExactAbs(zhangExactDeterminant(
			zhangCanonicalTransition(currentA, currentB))) == 1);

	std::map<ZhangGraphEdge, int> versions;
	for (const auto& edge : edges)
	{
		versions[edge] = 9;
	}
	auto fixedProductsA = zhangBuildProductIntegerFunctionals(
		fixedProduct, versions, g01, 12);
	auto fixedProductsB = zhangBuildProductIntegerFunctionals(
		fixedProduct, versions, g01, 12);
	BOOST_REQUIRE_EQUAL(fixedProductsA.size(), fixedProductsB.size());
	for (const auto& [satellite, product] : fixedProductsA)
	{
		BOOST_CHECK_EQUAL(
			zhangProductIntegerFunctionalFingerprint(product),
			zhangProductIntegerFunctionalFingerprint(fixedProductsB.at(satellite)));
	}
	auto userA = zhangBuildJointUserIntegerFunctional(fixedProductsA, g01, 3);
	auto userB = zhangBuildJointUserIntegerFunctional(fixedProductsB, g01, 3);
	BOOST_REQUIRE(zhangAuditUserIntegerLattice(userA).valid);
	BOOST_REQUIRE(zhangAuditUserIntegerLattice(userB).valid);
	BOOST_CHECK(userA.integerRows == userB.integerRows);
	BOOST_CHECK(
		zhangClassifyTemporalIntegerDatumAction(false, true, false) ==
		ZhangTemporalIntegerDatumAction::EXACT_TRANSPORT_NO_BESD);
}

BOOST_AUTO_TEST_CASE(e25b_cycle_slip_requires_besd_or_datum_reset)
{
	BOOST_CHECK(
		zhangClassifyTemporalIntegerDatumAction(true, true, true) ==
		ZhangTemporalIntegerDatumAction::ESTIMATE_BESD);
	BOOST_CHECK(
		zhangClassifyTemporalIntegerDatumAction(true, true, false) ==
		ZhangTemporalIntegerDatumAction::RESET_PRODUCT_DATUM);
	BOOST_CHECK(
		zhangClassifyTemporalIntegerDatumAction(false, false, false) ==
		ZhangTemporalIntegerDatumAction::RESET_PRODUCT_DATUM);
}

BOOST_AUTO_TEST_CASE(targeted_besd_selector_routes_only_nonheld_retired_arcs)
{
	const SatSys g08(E_Sys::GPS, 8);
	const ZhangGraphEdge active{"R0", g08};
	const ZhangGraphEdge retired{"R1", g08};
	ZhangProductIntegerTransition transition;
	transition.physicalEdges = {active, retired};
	transition.physicalArcVersions = {3, 7};
	transition.coefficients = {1, -1};
	transition.valid = true;
	const std::set<ZhangGraphEdge> postEventEdges = {active};
	const std::map<ZhangGraphEdge, int> postEventVersions = {{active, 3}};

	const auto selected = zhangSelectTargetedBesdTransition(
		transition, postEventEdges, postEventVersions, false, false);
	BOOST_CHECK(selected.selected);
	BOOST_CHECK_EQUAL(selected.reason, "REQUIRES_BESD_RETIRED_ARC");
	BOOST_CHECK_EQUAL(selected.physicalTerms, 2);

	const auto held = zhangSelectTargetedBesdTransition(
		transition, postEventEdges, postEventVersions, true, false);
	BOOST_CHECK(!held.selected);
	BOOST_CHECK_EQUAL(held.reason, "EXACT_HELD_TRANSPORT");

	const auto reset = zhangSelectTargetedBesdTransition(
		transition, postEventEdges, postEventVersions, false, true);
	BOOST_CHECK(!reset.selected);
	BOOST_CHECK_EQUAL(reset.reason, "PHASE_SEGMENT_RESET");

	const std::set<ZhangGraphEdge> allEdges = {active, retired};
	const std::map<ZhangGraphEdge, int> allVersions = {
		{active, 3}, {retired, 7}};
	const auto current = zhangSelectTargetedBesdTransition(
		transition, allEdges, allVersions, false, false);
	BOOST_CHECK(!current.selected);
	BOOST_CHECK_EQUAL(current.reason, "CURRENT_PHYSICAL_GRAPH_RELATION");

	auto versionChanged = allVersions;
	versionChanged[retired] = 8;
	const auto changed = zhangSelectTargetedBesdTransition(
		transition, allEdges, versionChanged, false, false);
	BOOST_CHECK(changed.selected);
	BOOST_CHECK_EQUAL(changed.reason, "REQUIRES_BESD_RETIRED_ARC");
}

BOOST_AUTO_TEST_CASE(e25b_product_phase_transport_is_current_node_plus_cycle_target)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	ZhangGraphEdge first{"R0", g01};
	ZhangGraphEdge second{"R0", g02};
	ZhangProductIntegerFunctional reference;
	reference.satellite = g01;
	reference.referenceSatellite = g01;
	reference.physicalEdges = {first, second};
	reference.networkCoefficients = {0, 0};
	reference.physicalArcVersions = {3, 3};
	reference.affineOffsetCycles = 2;
	reference.valid = true;
	auto satellite = reference;
	satellite.satellite = g02;
	satellite.networkCoefficients = {2, -1};
	satellite.affineOffsetCycles = 5;
	std::map<SatSys, ZhangProductIntegerFunctional> products = {
		{g01, reference}, {g02, satellite}
	};
	auto joint = zhangBuildJointUserIntegerFunctional(products, g01, 0);
	BOOST_REQUIRE(zhangAuditUserIntegerLattice(joint).valid);

	// Primitive order: two network arcs, then user G01/G02 ambiguities.
	ZhangExactVector primitive = {4, -3, 11, 24};
	const long long networkPath = 2 * 4 - (-3);
	const long long userSd = 24 - 11;
	const long long productAffineDifference = 5 - 2;
	VectorXd value = joint.value(primitive, VectorXd());
	BOOST_REQUIRE_EQUAL(value.size(), 1);
	BOOST_CHECK_EQUAL(
		value(0), userSd - networkPath - productAffineDifference);
	BOOST_CHECK_EQUAL(
		joint.affineOffsetsCycles(0), -productAffineDifference);

	// The service state is the satellite phase bias B_P.  The current-tree
	// state already contains z_T, so exact product transport must add G*k.
	// Since the user model applies -(C-B_P) = -C+B_P, the resulting ambiguity
	// is the user SD minus the complete product-tree node potential.
	const long long currentTreeNode = -5;
	const long long productTreeNode =
		currentTreeNode + networkPath + productAffineDifference;
	BOOST_CHECK_EQUAL(
		productTreeNode,
		currentTreeNode + networkPath + productAffineDifference);
	BOOST_CHECK_EQUAL(userSd - productTreeNode,
		userSd - currentTreeNode - networkPath - productAffineDifference);
}

BOOST_AUTO_TEST_CASE(product_physical_identity_ignores_unrelated_tree_generation)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	ZhangProductIntegerFunctional functional;
	functional.satellite = g02;
	functional.referenceSatellite = g01;
	functional.physicalEdges = {{"R0", g01}, {"R0", g02}, {"R1", g02}};
	functional.networkCoefficients = {1, -1, 0};
	functional.physicalArcVersions = {4, 7, 99};
	functional.temporalBasisVersion = 12;
	functional.valid = true;

	auto nextGeneration = functional;
	nextGeneration.temporalBasisVersion = 13;
	nextGeneration.physicalArcVersions[2] = 100;
	BOOST_CHECK_EQUAL(
		zhangProductPhysicalFunctionalFingerprint(functional),
		zhangProductPhysicalFunctionalFingerprint(nextGeneration));
	BOOST_CHECK_NE(
		zhangProductIntegerFunctionalFingerprint(functional),
		zhangProductIntegerFunctionalFingerprint(nextGeneration));

	auto changedSupport = nextGeneration;
	changedSupport.physicalArcVersions[1] = 8;
	BOOST_CHECK_NE(
		zhangProductPhysicalFunctionalFingerprint(functional),
		zhangProductPhysicalFunctionalFingerprint(changedSupport));
}

BOOST_AUTO_TEST_CASE(product_functional_difference_preserves_arc_versions_exactly)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	ZhangGraphEdge a{"R0", g01};
	ZhangGraphEdge b{"R0", g02};
	ZhangProductIntegerFunctional previous;
	previous.satellite = g02;
	previous.referenceSatellite = g01;
	previous.physicalEdges = {a, b};
	previous.networkCoefficients = {1, -1};
	previous.physicalArcVersions = {4, 7};
	previous.affineOffsetCycles = 2;
	previous.valid = true;

	auto generationOnly = previous;
	generationOnly.temporalBasisVersion = 99;
	auto zero = zhangProductIntegerFunctionalDifference(
		previous, generationOnly);
	BOOST_REQUIRE(zero.valid);
	BOOST_CHECK(zero.coefficients.empty());
	BOOST_CHECK_EQUAL(zero.affineOffsetCycles, 0);

	auto current = previous;
	current.physicalArcVersions[1] = 8;
	current.affineOffsetCycles = 5;
	auto changed = zhangProductIntegerFunctionalDifference(previous, current);
	BOOST_REQUIRE(changed.valid);
	BOOST_REQUIRE_EQUAL(changed.coefficients.size(), 2);
	BOOST_CHECK(changed.physicalEdges[0] == b);
	BOOST_CHECK(changed.physicalEdges[1] == b);
	BOOST_CHECK_EQUAL(changed.physicalArcVersions[0], 7);
	BOOST_CHECK_EQUAL(changed.physicalArcVersions[1], 8);
	BOOST_CHECK_EQUAL(changed.coefficients[0], 1);
	BOOST_CHECK_EQUAL(changed.coefficients[1], -1);
	BOOST_CHECK_EQUAL(changed.affineOffsetCycles, 3);
}

BOOST_AUTO_TEST_CASE(product_pair_difference_accepts_distinct_satellites_exactly)
{
	SatSys g01(E_Sys::GPS, 1);
	SatSys g02(E_Sys::GPS, 2);
	SatSys g03(E_Sys::GPS, 3);
	ZhangGraphEdge a{"R0", g01};
	ZhangGraphEdge b{"R0", g02};
	ZhangGraphEdge c{"R0", g03};
	ZhangProductIntegerFunctional first;
	first.satellite = g02;
	first.referenceSatellite = g01;
	first.physicalEdges = {a, b};
	first.networkCoefficients = {1, -1};
	first.physicalArcVersions = {4, 7};
	first.affineOffsetCycles = 2;
	first.valid = true;
	ZhangProductIntegerFunctional second;
	second.satellite = g03;
	second.referenceSatellite = g01;
	second.physicalEdges = {b, c};
	second.networkCoefficients = {1, -1};
	second.physicalArcVersions = {7, 9};
	second.affineOffsetCycles = 7;
	second.valid = true;

	BOOST_CHECK(!zhangProductIntegerFunctionalDifference(first, second).valid);
	const auto pair = zhangProductIntegerFunctionalPairDifference(first, second);
	BOOST_REQUIRE_MESSAGE(pair.valid, pair.failureReason);
	BOOST_REQUIRE_EQUAL(pair.coefficients.size(), 3);
	BOOST_CHECK_EQUAL(pair.coefficients[0], -1);
	BOOST_CHECK_EQUAL(pair.coefficients[1], 2);
	BOOST_CHECK_EQUAL(pair.coefficients[2], -1);
	BOOST_CHECK_EQUAL(pair.affineOffsetCycles, 5);
}

BOOST_AUTO_TEST_CASE(e25b_rejects_nuisance_fractional_offset_and_unsaturated_rows)
{
	ZhangJointUserIntegerFunctional valid;
	valid.integerRows = {{1, -1, 0}, {0, 1, -1}};
	valid.nuisanceRows = MatrixXd::Zero(2, 2);
	valid.affineOffsetsCycles = VectorXd::Zero(2);
	valid.valid = true;
	BOOST_REQUIRE(zhangAuditUserIntegerLattice(valid).valid);

	auto nuisanceLeak = valid;
	nuisanceLeak.nuisanceRows(0, 1) = 0.01;
	auto nuisanceAudit = zhangAuditUserIntegerLattice(nuisanceLeak);
	BOOST_CHECK(!nuisanceAudit.valid);
	BOOST_CHECK_EQUAL(
		nuisanceAudit.failureReason,
		"REAL_NUISANCE_LEAKS_INTO_INTEGER_FUNCTIONAL");

	auto fractionalOffset = valid;
	fractionalOffset.affineOffsetsCycles(1) = 0.25;
	auto offsetAudit = zhangAuditUserIntegerLattice(fractionalOffset);
	BOOST_CHECK(!offsetAudit.valid);
	BOOST_CHECK_EQUAL(offsetAudit.failureReason, "NON_INTEGER_AFFINE_OFFSET");

	auto unsaturated = valid;
	for (auto& row : unsaturated.integerRows)
	for (auto& coefficient : row)
	{
		coefficient *= 2;
	}
	auto saturationAudit = zhangAuditUserIntegerLattice(unsaturated);
	BOOST_CHECK(!saturationAudit.valid);
	BOOST_CHECK(!saturationAudit.primitiveAdmissible);
}

BOOST_AUTO_TEST_CASE(e29_checkpoint_core_roundtrip_is_bitwise_and_preserves_callbacks)
{
	TemporaryCheckpointFile file("_roundtrip.bin");
	auto bundle = makeCheckpointTestBundle();
	auto writeResult = writeZhangCheckpointBundle(
		file.path.string(), bundle);
	BOOST_REQUIRE_MESSAGE(writeResult.valid, writeResult.failureReason);
	BOOST_CHECK_EQUAL(writeResult.payloadSha256.size(), 64);
	BOOST_CHECK_GT(writeResult.payloadBytes, 0);
	std::string fileHashFailure;
	BOOST_CHECK_EQUAL(
		zhangCheckpointFileSha256(file.path.string(), &fileHashFailure).size(),
		64);
	BOOST_CHECK_EQUAL(fileHashFailure, "NONE");

	ZhangCheckpointBundle restoredBundle;
	auto readResult = readZhangCheckpointBundle(
		file.path.string(), checkpointTestExpectations(), restoredBundle);
	BOOST_REQUIRE_MESSAGE(readResult.valid, readResult.failureReason);
	BOOST_CHECK_EQUAL(readResult.payloadSha256, writeResult.payloadSha256);
	BOOST_REQUIRE_EQUAL(restoredBundle.sections.count("zhang.graph"), 1);
	const auto& graphSection = restoredBundle.sections.at("zhang.graph");
	BOOST_CHECK_EQUAL(graphSection.schemaVersion, 1);
	BOOST_CHECK_EQUAL(graphSection.payload, "pointer-free-graph-runtime");
	BOOST_CHECK_EQUAL(
		graphSection.sha256,
		zhangCheckpointSha256(graphSection.payload));

	KFState unresolvedDestination;
	const VectorXd unresolvedBefore = unresolvedDestination.x;
	std::string unresolvedFailure;
	BOOST_CHECK(!restoreZhangCheckpointKfCoreWithReceiverResolver(
		restoredBundle.kfCore,
		unresolvedDestination,
		[](const std::string&) -> Receiver* { return nullptr; },
		&unresolvedFailure));
	BOOST_CHECK_EQUAL(
		unresolvedFailure,
		"CHECKPOINT_CORE_RECEIVER_POINTER_REBIND_FAILED:R0");
	BOOST_CHECK(
		(unresolvedDestination.x.array() == unresolvedBefore.array()).all());

	KFState destination;
	destination.acceptedMeasurementFactorCallback = [](
		const KFState&,
		const KFMeas&,
		const std::string&,
		const VectorXd&,
		const MatrixXd&,
		std::uint64_t,
		std::uint64_t) { return true; };
	std::string restoreFailure;
	BOOST_REQUIRE(restoreZhangCheckpointKfCoreWithReceiverResolver(
		restoredBundle.kfCore,
		destination,
		[](const std::string& id) -> Receiver*
		{
			return id == "R0" ? &checkpointTestReceiver() : nullptr;
		},
		&restoreFailure));
	BOOST_CHECK_EQUAL(restoreFailure, "NONE");
	BOOST_CHECK(static_cast<bool>(
		destination.acceptedMeasurementFactorCallback));
	BOOST_CHECK_EQUAL(zhangCheckpointRuntimeId(destination), "runtime-00");
	std::string bindFailure;
	BOOST_CHECK(bindZhangCheckpointRuntimeId(
		destination, "runtime-00", &bindFailure));
	BOOST_CHECK_EQUAL(bindFailure, "NONE");
	BOOST_CHECK(!bindZhangCheckpointRuntimeId(
		destination, "different-runtime", &bindFailure));
	BOOST_CHECK_EQUAL(bindFailure, "CHECKPOINT_RUNTIME_ID_ALREADY_BOUND");
	BOOST_CHECK(zhangCheckpointKfCoreBitwiseEqual(
		restoredBundle.kfCore,
		captureZhangCheckpointKfCore(destination)));
	BOOST_REQUIRE_EQUAL(destination.filterChunkMap.count("zhang"), 1);
	const auto& chunk = destination.filterChunkMap.at("zhang");
	BOOST_CHECK_EQUAL(chunk.begH, 7);
	BOOST_CHECK_EQUAL(chunk.numH, 11);
	bool estimatedTimeRestored = false;
	bool receiverPointerRestored = false;
	for (const auto& [key, index] : destination.kfIndexMap)
	{
		if (key.type == KF::REC_CLOCK && key.str == "R0")
		{
			receiverPointerRestored =
				key.rec_ptr == &checkpointTestReceiver();
		}
		if (key.type == KF::SAT_CLOCK && key.Sat == SatSys(E_Sys::GPS, 7))
		{
			estimatedTimeRestored =
				key.estimatedTime.bigTime == 123456710.25L;
		}
	}
	BOOST_CHECK(receiverPointerRestored);
	BOOST_CHECK(estimatedTimeRestored);
}

BOOST_AUTO_TEST_CASE(e29_checkpoint_rejects_corruption_and_provenance_drift)
{
	TemporaryCheckpointFile file("_corrupt.bin");
	auto bundle = makeCheckpointTestBundle();
	auto writeResult = writeZhangCheckpointBundle(
		file.path.string(), bundle);
	BOOST_REQUIRE_MESSAGE(writeResult.valid, writeResult.failureReason);

	auto wrongExpectations = checkpointTestExpectations();
	wrongExpectations.configSha256 = std::string(64, 'd');
	ZhangCheckpointBundle ignored;
	auto provenanceResult = readZhangCheckpointBundle(
		file.path.string(), wrongExpectations, ignored);
	BOOST_CHECK(!provenanceResult.valid);
	BOOST_CHECK_EQUAL(
		provenanceResult.failureReason,
		"CHECKPOINT_PROVENANCE_MISMATCH");

	std::fstream stream(
		file.path, std::ios::binary | std::ios::in | std::ios::out);
	BOOST_REQUIRE(stream);
	stream.seekg(-1, std::ios::end);
	char byte = 0;
	stream.read(&byte, 1);
	BOOST_REQUIRE(stream);
	byte ^= 0x5a;
	stream.seekp(-1, std::ios::end);
	stream.write(&byte, 1);
	stream.flush();
	BOOST_REQUIRE(stream);
	stream.close();

	auto corruptResult = readZhangCheckpointBundle(
		file.path.string(), checkpointTestExpectations(), ignored);
	BOOST_CHECK(!corruptResult.valid);
	BOOST_CHECK_EQUAL(
		corruptResult.failureReason,
		"CHECKPOINT_PAYLOAD_SHA256_MISMATCH");
}

BOOST_AUTO_TEST_CASE(e29_checkpoint_writer_fails_closed_on_invalid_identity_and_index)
{
	TemporaryCheckpointFile missingIdentityFile("_identity.bin");
	auto missingIdentity = makeCheckpointTestBundle();
	missingIdentity.manifest.runtimeId.clear();
	auto identityResult = writeZhangCheckpointBundle(
		missingIdentityFile.path.string(), missingIdentity);
	BOOST_CHECK(!identityResult.valid);
	BOOST_CHECK_EQUAL(
		identityResult.failureReason,
		"CHECKPOINT_MANIFEST_IDENTITY_MISSING");
	BOOST_CHECK(!std::filesystem::exists(missingIdentityFile.path));

	TemporaryCheckpointFile badIndexFile("_index.bin");
	auto badIndex = makeCheckpointTestBundle();
	auto firstIndex = badIndex.kfCore.kfIndexMap.begin();
	auto secondIndex = std::next(firstIndex);
	secondIndex->second = firstIndex->second;
	auto indexResult = writeZhangCheckpointBundle(
		badIndexFile.path.string(), badIndex);
	BOOST_CHECK(!indexResult.valid);
	BOOST_CHECK_EQUAL(
		indexResult.failureReason,
		"CHECKPOINT_CORE_INDEX_NOT_BIJECTIVE");
	BOOST_CHECK(!std::filesystem::exists(badIndexFile.path));

	TemporaryCheckpointFile badContentFile("_content.bin");
	auto badContent = makeCheckpointTestBundle();
	badContent.manifest.configText += "-tampered";
	auto contentResult = writeZhangCheckpointBundle(
		badContentFile.path.string(), badContent);
	BOOST_CHECK(!contentResult.valid);
	BOOST_CHECK_EQUAL(
		contentResult.failureReason,
		"CHECKPOINT_MANIFEST_CONTENT_HASH_MISMATCH");
}

BOOST_AUTO_TEST_CASE(e29_checkpoint_required_sections_are_strictly_validated)
{
	auto bundle = makeCheckpointTestBundle();
	const std::vector<ZhangCheckpointSectionRequirement> requirements = {
		{"zhang.graph", 1}};
	std::string failure;
	BOOST_CHECK(validateZhangCheckpointRequiredSections(
		bundle, requirements, &failure));
	BOOST_CHECK_EQUAL(failure, "NONE");

	auto missing = bundle;
	missing.sections.clear();
	BOOST_CHECK(!validateZhangCheckpointRequiredSections(
		missing, requirements, &failure));
	BOOST_CHECK_EQUAL(
		failure, "CHECKPOINT_REQUIRED_SECTION_MISSING:zhang.graph");

	auto wrongVersion = bundle;
	wrongVersion.sections.at("zhang.graph").schemaVersion = 2;
	BOOST_CHECK(!validateZhangCheckpointRequiredSections(
		wrongVersion, requirements, &failure));
	BOOST_CHECK_EQUAL(
		failure,
		"CHECKPOINT_REQUIRED_SECTION_VERSION_MISMATCH:zhang.graph");

	auto corrupt = bundle;
	corrupt.sections.at("zhang.graph").payload += "-corrupt";
	BOOST_CHECK(!validateZhangCheckpointRequiredSections(
		corrupt, requirements, &failure));
	BOOST_CHECK_EQUAL(
		failure, "CHECKPOINT_REQUIRED_SECTION_HASH_MISMATCH:zhang.graph");
}

BOOST_AUTO_TEST_CASE(e29_checkpoint_manifest_json_is_atomic_and_auditable)
{
	TemporaryCheckpointFile file("_manifest.json");
	auto bundle = makeCheckpointTestBundle();
	auto writeResult = writeZhangCheckpointManifestJson(
		file.path.string(), bundle);
	BOOST_REQUIRE_MESSAGE(writeResult.valid, writeResult.failureReason);
	BOOST_CHECK_EQUAL(writeResult.payloadSha256.size(), 64);
	BOOST_CHECK_GT(writeResult.payloadBytes, 0);

	std::ifstream input(file.path, std::ios::binary);
	BOOST_REQUIRE(input);
	const std::string json(
		(std::istreambuf_iterator<char>(input)),
		std::istreambuf_iterator<char>());
	BOOST_CHECK_NE(json.find("\"runtime_id\": \"runtime-00\""),
		std::string::npos);
	BOOST_CHECK_NE(json.find("\"state_dimension\": 3"),
		std::string::npos);
	BOOST_CHECK_NE(json.find("\"name\": \"zhang.graph\""),
		std::string::npos);
	BOOST_CHECK_NE(json.find(
		bundle.sections.at("zhang.graph").sha256), std::string::npos);

	auto secondWrite = writeZhangCheckpointManifestJson(
		file.path.string(), bundle);
	BOOST_CHECK(!secondWrite.valid);
	BOOST_CHECK_EQUAL(
		secondWrite.failureReason, "CHECKPOINT_TARGET_ALREADY_EXISTS");
}

BOOST_AUTO_TEST_CASE(product_relation_admission_commits_exact_redundant_batch)
{
	ZhangProductRelationAdmissionState state;
	auto candidate = [](
		const std::string& id,
		const std::string& satellite,
		const std::string& observable,
		const std::map<std::string, ZhangExactInteger>& row,
		long long value)
	{
		ZhangProductRelationAdmissionCandidate result;
		result.relationId = id;
		result.satellite = satellite;
		result.observable = observable;
		result.physicalCoefficients = row;
		result.integerValue = value;
		result.exactIntegerEstimable = true;
		result.phaseSegmentCompatible = true;
		result.scalarReliabilityPassed = true;
		result.jointNisPassed = true;
		return result;
	};
	const std::vector<ZhangProductRelationAdmissionCandidate> rows = {
		candidate("G02-L1-a", "G02", "L1C", {{"a", 1}}, 3),
		candidate("G02-L1-a-repeat", "G02", "L1C", {{"a", 1}}, 3),
		candidate("G02-L2-b", "G02", "L2W", {{"b", 1}}, -2),
		candidate("G02-L2-b-repeat", "G02", "L2W", {{"b", 1}}, -2),
	};
	const auto result = ProductRelationAdmission::admit(state, rows);
	BOOST_CHECK(result.committed);
	BOOST_CHECK_EQUAL(result.status, "CERTIFIED_NEW_RELATION");
	BOOST_CHECK_EQUAL(result.candidateRows, 4);
	BOOST_CHECK_EQUAL(result.candidateExactRank, 2);
	BOOST_CHECK_EQUAL(result.candidateRedundantRows, 2);
	BOOST_CHECK(result.candidateCycleClosureConsistent);
	BOOST_CHECK(result.persistentCycleClosureConsistent);
	BOOST_CHECK_EQUAL(result.persistentRankAfter, 2);
	BOOST_CHECK_EQUAL(result.restoredSatellites, 1);
	BOOST_CHECK_EQUAL(state.certifiedSatellites.count("G02"), 1);
}

BOOST_AUTO_TEST_CASE(product_relation_admission_aborts_conflicting_cycle_atomically)
{
	ZhangProductRelationAdmissionState state;
	ZhangProductRelationAdmissionCandidate first;
	first.relationId = "first";
	first.satellite = "G02";
	first.observable = "L1C";
	first.physicalCoefficients = {{"a", 1}};
	first.integerValue = 3;
	first.exactIntegerEstimable = true;
	first.phaseSegmentCompatible = true;
	first.scalarReliabilityPassed = true;
	first.jointNisPassed = true;
	auto duplicate = first;
	duplicate.relationId = "conflict";
	duplicate.integerValue = 4;
	const auto result = ProductRelationAdmission::admit(
		state, {first, duplicate});
	BOOST_CHECK(!result.committed);
	BOOST_CHECK_EQUAL(
		result.status, "ABORT_INCONSISTENT_CANDIDATE_CYCLE_CLOSURE");
	BOOST_CHECK(state.certifiedRows.empty());
	BOOST_CHECK(state.certifiedSatellites.empty());
}

BOOST_AUTO_TEST_CASE(product_relation_admission_waits_for_redundancy)
{
	ZhangProductRelationAdmissionState state;
	ZhangProductRelationAdmissionCandidate candidate;
	candidate.relationId = "single-bridge";
	candidate.satellite = "G02";
	candidate.observable = "L1C";
	candidate.physicalCoefficients = {{"a", 1}};
	candidate.integerValue = 3;
	candidate.exactIntegerEstimable = true;
	candidate.phaseSegmentCompatible = true;
	candidate.scalarReliabilityPassed = true;
	candidate.jointNisPassed = true;
	const auto result = ProductRelationAdmission::admit(state, {candidate});
	BOOST_CHECK(!result.committed);
	BOOST_CHECK_EQUAL(
		result.status, "PREPARE_MERGE_AWAITING_REDUNDANCY");
	BOOST_CHECK(state.certifiedRows.empty());
	BOOST_CHECK_EQUAL(state.pendingCandidates.size(), 1);
	auto redundant = candidate;
	redundant.relationId = "redundant-bridge";
	const auto committed = ProductRelationAdmission::admit(
		state, {redundant});
	BOOST_CHECK(committed.committed);
	BOOST_CHECK_EQUAL(committed.candidateRows, 2);
	BOOST_CHECK_EQUAL(committed.candidateExactRank, 1);
	BOOST_CHECK_EQUAL(committed.candidateRedundantRows, 1);
	BOOST_CHECK(state.pendingCandidates.empty());
	BOOST_CHECK_EQUAL(state.certifiedRows.size(), 1);
}

BOOST_AUTO_TEST_CASE(product_relation_admission_fails_closed_at_every_gate)
{
	auto validCandidate = []()
	{
		ZhangProductRelationAdmissionCandidate candidate;
		candidate.relationId = "gate-a";
		candidate.satellite = "G02";
		candidate.observable = "L1C";
		candidate.physicalCoefficients = {{"a", 1}};
		candidate.integerValue = 3;
		candidate.exactIntegerEstimable = true;
		candidate.phaseSegmentCompatible = true;
		candidate.scalarReliabilityPassed = true;
		candidate.jointNisPassed = true;
		return candidate;
	};
	struct GateCase
	{
		std::string expected;
		std::function<void(ZhangProductRelationAdmissionCandidate&)> fail;
	};
	const std::vector<GateCase> cases = {
		{"REJECTED_NOT_EXACT_INTEGER_ESTIMABLE", [](auto& row)
			{ row.exactIntegerEstimable = false; }},
		{"REJECTED_PHASE_SEGMENT_INCOMPATIBLE", [](auto& row)
			{ row.phaseSegmentCompatible = false; }},
		{"REJECTED_SCALAR_RELIABILITY", [](auto& row)
			{ row.scalarReliabilityPassed = false; }},
		{"REJECTED_JOINT_NIS", [](auto& row)
			{ row.jointNisPassed = false; }},
	};
	for (const auto& test : cases)
	{
		ZhangProductRelationAdmissionState state;
		auto first = validCandidate();
		auto second = first;
		second.relationId = "gate-b";
		test.fail(first);
		const auto result = ProductRelationAdmission::admit(
			state, {first, second});
		BOOST_CHECK(!result.committed);
		BOOST_CHECK_EQUAL(result.status, test.expected);
		BOOST_CHECK(state.certifiedRows.empty());
		BOOST_CHECK(state.pendingCandidates.empty());
	}
}

BOOST_AUTO_TEST_CASE(product_relation_admission_requires_redundancy_per_signal)
{
	auto make = [](const std::string& id, const std::string& observable,
		const std::string& column, long long value)
	{
		ZhangProductRelationAdmissionCandidate candidate;
		candidate.relationId = id;
		candidate.satellite = "G02";
		candidate.observable = observable;
		candidate.physicalCoefficients = {{column, 1}};
		candidate.integerValue = value;
		candidate.exactIntegerEstimable = true;
		candidate.phaseSegmentCompatible = true;
		candidate.scalarReliabilityPassed = true;
		candidate.jointNisPassed = true;
		return candidate;
	};
	ZhangProductRelationAdmissionState state;
	const auto result = ProductRelationAdmission::admit(state, {
		make("l1-a", "L1C", "a", 3),
		make("l1-b", "L1C", "a", 3),
		make("l2-a", "L2W", "b", -2),
	});
	BOOST_CHECK(!result.committed);
	BOOST_CHECK_EQUAL(
		result.status, "PREPARE_MERGE_AWAITING_REDUNDANCY");
	BOOST_CHECK_EQUAL(result.observableGroups, 2);
	BOOST_CHECK_EQUAL(result.redundancyCheckedGroups, 1);
	BOOST_CHECK_EQUAL(state.pendingCandidates.size(), 3);
}

BOOST_AUTO_TEST_CASE(temporal_certificate_confirmation_is_value_gap_and_kind_safe)
{
	TemporalCertificateConfirmationState state;
	auto first = zhangConfirmTemporalCertificate(
		state, "L1C:G03-G02", ZhangExactInteger(7), 100,
		"path-a", false, 3, 30, false);
	BOOST_CHECK(!first.accepted);
	BOOST_CHECK_EQUAL(first.consistentEpochs, 1);

	// Re-evaluating the same epoch must not manufacture confirmation count.
	auto duplicate = zhangConfirmTemporalCertificate(
		state, "L1C:G03-G02", ZhangExactInteger(7), 100,
		"path-a", false, 3, 30, false);
	BOOST_CHECK_EQUAL(duplicate.consistentEpochs, 1);

	zhangConfirmTemporalCertificate(
		state, "L1C:G03-G02", ZhangExactInteger(7), 110,
		"path-a", false, 3, 30, false);
	auto third = zhangConfirmTemporalCertificate(
		state, "L1C:G03-G02", ZhangExactInteger(7), 120,
		"path-a", false, 3, 30, false);
	BOOST_CHECK(third.accepted);

	// An integer change is a new hypothesis, not a continuation.
	auto changed = zhangConfirmTemporalCertificate(
		state, "L1C:G03-G02", ZhangExactInteger(8), 130,
		"path-a", false, 3, 30, false);
	BOOST_CHECK(changed.reset);
	BOOST_CHECK_EQUAL(changed.consistentEpochs, 1);

	// A bridge additionally requires independent redundant support.
	TemporalCertificateConfirmationState bridge;
	auto noRedundancy = zhangConfirmTemporalCertificate(
		bridge, "L1C:G05-G02", ZhangExactInteger(2), 200,
		"path-a", false, 1, 30, true);
	BOOST_CHECK(!noRedundancy.accepted);
	BOOST_CHECK_EQUAL(noRedundancy.reason, "AWAITING_REDUNDANCY");
	auto redundant = zhangConfirmTemporalCertificate(
		bridge, "L1C:G05-G02", ZhangExactInteger(2), 210,
		"path-b", true, 1, 30, true);
	BOOST_CHECK(redundant.accepted);

	BOOST_CHECK_EQUAL(
		zhangTemporalCertificateKindName(
			TemporalCertificateKind::SELF_GAUGE_SHIFT),
		"SELF_GAUGE_SHIFT");
	BOOST_CHECK_EQUAL(
		zhangTemporalCertificateKindName(
			TemporalCertificateKind::INTER_SATELLITE_BRIDGE),
		"INTER_SATELLITE_BRIDGE");
}

BOOST_AUTO_TEST_CASE(targeted_besd_tracker_matches_augmented_kalman_update)
{
	VectorXd stateMean(2);
	stateMean << 1.2, -0.7;
	MatrixXd stateCovariance(2, 2);
	stateCovariance << 4.0, 0.6,
		0.6, 2.0;
	MatrixXd targets(2, 2);
	targets << 1, -1,
		2, 1;
	VectorXd offsets(2);
	offsets << 3, -2;

	ZhangTargetedBesdTracker tracker;
	BOOST_REQUIRE(tracker.initialise(
		{"old", "new"}, targets, offsets, stateMean, stateCovariance));

	MatrixXd design(1, 2);
	design << 0.5, 1.5;
	MatrixXd noise(1, 1);
	noise << 0.25;
	VectorXd residual(1);
	residual << -0.8;

	// Independent augmented-state reference update for [x, f].
	VectorXd jointMean(4);
	jointMean << stateMean, targets * stateMean + offsets;
	MatrixXd jointCovariance(4, 4);
	jointCovariance.topLeftCorner(2, 2) = stateCovariance;
	jointCovariance.bottomLeftCorner(2, 2) = targets * stateCovariance;
	jointCovariance.topRightCorner(2, 2) =
		jointCovariance.bottomLeftCorner(2, 2).transpose();
	jointCovariance.bottomRightCorner(2, 2) =
		targets * stateCovariance * targets.transpose();
	MatrixXd jointDesign = MatrixXd::Zero(1, 4);
	jointDesign.leftCols(2) = design;
	MatrixXd innovation = jointDesign * jointCovariance
		* jointDesign.transpose() + noise;
	MatrixXd gain = jointCovariance * jointDesign.transpose()
		* innovation.inverse();
	jointMean += gain * residual;
	jointCovariance -= gain * innovation * gain.transpose();
	jointCovariance = 0.5
		* (jointCovariance + jointCovariance.transpose());

	BOOST_REQUIRE(tracker.updateAcceptedMeasurement(
		stateCovariance, design, noise, residual));
	const auto marginal = tracker.marginal();
	BOOST_REQUIRE(marginal.valid);
	BOOST_CHECK_SMALL(
		(marginal.mean - jointMean.tail(2)).norm(), 1e-12);
	BOOST_CHECK_SMALL(
		(marginal.covariance
			- jointCovariance.bottomRightCorner(2, 2)).norm(), 1e-12);
	BOOST_CHECK_SMALL(
		(tracker.crossCovariance()
			- jointCovariance.bottomLeftCorner(2, 2)).norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(targeted_besd_tracker_carries_only_target_schur_boundary)
{
	VectorXd stateMean(3);
	stateMean << 1, 2, 3;
	MatrixXd stateCovariance = MatrixXd::Identity(3, 3);
	MatrixXd targets(2, 3);
	targets << 1, 0, -1,
		0, 2, 1;
	ZhangTargetedBesdTracker tracker;
	BOOST_REQUIRE(tracker.initialise(
		{"old", "new"}, targets, VectorXd::Zero(2),
		stateMean, stateCovariance));

	MatrixXd transition(2, 3);
	transition << 1, 0, 0,
		0, 1, 1;
	const MatrixXd expectedCross =
		tracker.crossCovariance() * transition.transpose();
	const auto before = tracker.marginal();
	BOOST_REQUIRE(tracker.advanceState(transition));
	const auto after = tracker.marginal();
	BOOST_REQUIRE(after.valid);
	BOOST_CHECK_EQUAL(tracker.targetCount(), 2);
	BOOST_CHECK_EQUAL(tracker.currentStateDimension(), 2);
	BOOST_CHECK_SMALL((after.mean - before.mean).norm(), 1e-15);
	BOOST_CHECK_SMALL(
		(after.covariance - before.covariance).norm(), 1e-15);
	BOOST_CHECK_SMALL(
		(tracker.crossCovariance() - expectedCross).norm(), 1e-15);

	// A malformed update fails closed and leaves no usable marginal.
	BOOST_CHECK(!tracker.updateAcceptedMeasurement(
		MatrixXd::Identity(3, 3), MatrixXd::Zero(1, 3),
		MatrixXd::Identity(1, 1), VectorXd::Zero(1)));
	BOOST_CHECK(!tracker.isActive());
	BOOST_CHECK(!tracker.marginal().valid);
}

BOOST_AUTO_TEST_CASE(product_relation_score_is_reliability_lexicographic)
{
	ZhangProductRelationLexicographicScore unreliable;
	unreliable.componentCoverageGain = 100;
	unreliable.productInformationGain = 1;
	ZhangProductRelationLexicographicScore reliable;
	reliable.reliabilityPassed = true;
	reliable.componentCoverageGain = 1;
	reliable.productInformationGain = 1e-6;
	BOOST_CHECK(unreliable < reliable);

	ZhangProductRelationLexicographicScore moreCoverage = reliable;
	moreCoverage.componentCoverageGain = 2;
	BOOST_CHECK(reliable < moreCoverage);

	ProductParBranch lowerGain;
	lowerGain.reliabilityPassed = true;
	lowerGain.componentCoverageGain = 2;
	lowerGain.productInformationGain = 0.2;
	ProductParBranch higherGain = lowerGain;
	higherGain.productInformationGain = 0.3;
	BOOST_CHECK(zhangProductParScore(lowerGain) <
		zhangProductParScore(higherGain));

	ProductParBranch weakUnreliable;
	weakUnreliable.integerRank = 20;
	weakUnreliable.rawPartialFixedRank = 3;
	weakUnreliable.partialFixFraction = 0.15;
	weakUnreliable.componentCoverageGain = 20;
	weakUnreliable.productInformationGain = 1;
	ProductParBranch fixableUnreliable = weakUnreliable;
	fixableUnreliable.integerRank = 6;
	fixableUnreliable.rawPartialFixedRank = 5;
	fixableUnreliable.partialFixFraction = 5.0 / 6.0;
	fixableUnreliable.componentCoverageGain = 6;
	fixableUnreliable.productInformationGain = 0.2;
	BOOST_CHECK(zhangProductParScore(weakUnreliable) <
		zhangProductParScore(fixableUnreliable));
}

BOOST_AUTO_TEST_CASE(product_relation_named_ordering_is_semantic)
{
	const SatSys g01("G01");
	const SatSys g02("G02");
	const SatSys g03("G03");
	ZhangProductRelationBasis first;
	first.mappableNamedIndices = {0, 1};
	first.namedRelations.resize(2);
	first.namedRelations[0].satellite = g02;
	first.namedRelations[0].referenceSatellite = g01;
	first.namedRelations[1].satellite = g03;
	first.namedRelations[1].referenceSatellite = g01;
	ZhangProductRelationBasis second = first;
	BOOST_CHECK(zhangProductNamedOrderingMatches(first, second));

	// Equal numeric indices do not rescue a semantic L1/L2 row swap.
	std::swap(second.namedRelations[0].satellite,
		second.namedRelations[1].satellite);
	BOOST_CHECK(!zhangProductNamedOrderingMatches(first, second));
}

BOOST_AUTO_TEST_CASE(product_relation_wl_gain_uses_joint_cross_covariance)
{
	MatrixXd jointCovariance(2, 2);
	jointCovariance << 4, 3,
		3, 4;
	ZhangIarFunctional wideLane(1, 2);
	wideLane.insert(0, 0) = 1;
	wideLane.insert(0, 1) = -1;
	wideLane.makeCompressed();
	const double gain = zhangNamedProductInformationGain(
		jointCovariance, wideLane);
	BOOST_CHECK_CLOSE(gain, 0.125, 1e-10);

	// This is Q11+Q22-Q12-Q21 = 2, not the cross-covariance-free value 8.
	const MatrixXd wideLaneCovariance = wideLane * jointCovariance *
		wideLane.transpose();
	BOOST_CHECK_CLOSE(wideLaneCovariance(0, 0), 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(product_gain_spectrum_separates_rank_and_search_failures)
{
	BOOST_CHECK_EQUAL(
		zhangProductGainSpectrumDiagnosis(0.82, 0.03),
		"INTEGER_CANDIDATE_SUBSPACE_MISALIGNED");
	BOOST_CHECK_EQUAL(
		zhangProductGainSpectrumDiagnosis(0.08, 0.03),
		"REAL_RANK_CEILING_LOW_INCREASE_RANK");
	BOOST_CHECK_EQUAL(
		zhangProductGainSpectrumDiagnosis(0.82, 0.50),
		"INTEGER_SUBSET_USES_REAL_CEILING_EFFICIENTLY");
}

BOOST_AUTO_TEST_CASE(product_relation_wl_l1_transform_is_unimodular)
{
	ZhangExactMatrix transform = {
		{1, 0, -1, 0},
		{0, 1, 0, -1},
		{1, 0, 0, 0},
		{0, 1, 0, 0}
	};
	const auto smith = zhangIntegerRowLatticeContains(
		transform, ZhangExactVector(4));
	BOOST_REQUIRE_EQUAL(smith.smithInvariants.size(), 4);
	for (const auto& invariant : smith.smithInvariants)
	{
		BOOST_CHECK_EQUAL(invariant, 1);
	}
}

BOOST_AUTO_TEST_CASE(product_relation_partial_decorrelated_rows_are_not_certificate)
{
	const auto none = zhangRecoverCertifiedNamedProductSubset(
		{{1, 1}}, {7}, 2);
	BOOST_CHECK(none.empty());

	const auto certifiedPartial = zhangRecoverCertifiedNamedProductSubset(
		{{1, 1, 0}, {0, 1, 0}}, {7, 3}, 3);
	BOOST_REQUIRE_EQUAL(certifiedPartial.size(), 2);
	BOOST_CHECK_EQUAL(certifiedPartial.at(0), 4);
	BOOST_CHECK_EQUAL(certifiedPartial.at(1), 3);
	BOOST_CHECK(!certifiedPartial.contains(2));

	const auto partial = zhangRecoverCompleteNamedProductSubset(
		{{1, 1}}, {7}, 2);
	BOOST_CHECK(partial.empty());

	const auto complete = zhangRecoverCompleteNamedProductSubset(
		{{1, 1}, {0, 1}}, {7, 3}, 2);
	BOOST_REQUIRE_EQUAL(complete.size(), 2);
	BOOST_CHECK_EQUAL(complete.at(0), 4);
	BOOST_CHECK_EQUAL(complete.at(1), 3);
}

BOOST_AUTO_TEST_CASE(canonical_physical_product_rows_ignore_restricted_product_core)
{
    std::set<SatSys> satellites;
    std::set<ZhangGraphEdge> edges;
    for (int prn = 1; prn <= 29; prn++)
    {
        const SatSys satellite(E_Sys::GPS, prn);
        satellites.insert(satellite);
        edges.insert({"R0", satellite});
        edges.insert({"R1", satellite});
    }
    const ZhangGraphBasis current = zhangBuildSpanningTree(edges, "R0");
    std::set<ZhangGraphEdge> restrictedEdges;
    for (const auto& edge : edges)
        if (edge.satellite.prn <= 23) restrictedEdges.insert(edge);
    const ZhangGraphBasis restrictedProduct =
        zhangBuildRootedProductTree(restrictedEdges, "R0");
    BOOST_REQUIRE(current.connected);
    BOOST_REQUIRE(restrictedProduct.connected);
    BOOST_CHECK_EQUAL(restrictedProduct.satellites.size(), 23);

    std::vector<ZhangProductRelationRow> physical;
    std::string failure;
    BOOST_REQUIRE(zhangBuildCanonicalProductPhysicalRelations(
        current, satellites, SatSys(E_Sys::GPS, 1), physical, failure));
    BOOST_CHECK_EQUAL(physical.size(), 28);
    for (const auto& row : physical)
    {
        BOOST_CHECK(row.referenceSatellite == SatSys(E_Sys::GPS, 1));
        BOOST_CHECK_EQUAL(row.physicalArcCoefficients.size(), 2);
        ZhangExactInteger receiverSum = 0;
        for (const auto& [edge, coefficient] : row.physicalArcCoefficients)
            receiverSum += coefficient;
        BOOST_CHECK_EQUAL(receiverSum, 0);
    }
}

BOOST_AUTO_TEST_CASE(canonical_physical_catalogue_retains_unrepresented_directions)
{
	std::set<SatSys> canonical;
	std::set<ZhangGraphEdge> representedEdges;
	for (int prn = 1; prn <= 29; prn++)
	{
		const SatSys satellite(E_Sys::GPS, prn);
		canonical.insert(satellite);
		if (prn <= 23)
		{
			representedEdges.insert({"R0", satellite});
			representedEdges.insert({"R1", satellite});
		}
	}
	const auto represented = zhangBuildSpanningTree(representedEdges, "R0");
	BOOST_REQUIRE(represented.connected);
	std::vector<ZhangProductRelationRow> physical;
	std::string failure;
	BOOST_REQUIRE(zhangBuildCanonicalProductPhysicalRelations(
		represented, canonical, SatSys(E_Sys::GPS, 1), physical, failure));
	BOOST_CHECK_EQUAL(failure, "PARTIAL_CANONICAL_PHYSICAL_MAPPING");
	BOOST_REQUIRE_EQUAL(physical.size(), 28);
	int mapped = 0;
	int unresolved = 0;
	for (const auto& row : physical)
	{
		BOOST_CHECK(row.referenceSatellite == SatSys(E_Sys::GPS, 1));
		if (row.physicalArcCoefficients.empty()) unresolved++;
		else mapped++;
	}
	BOOST_CHECK_EQUAL(mapped, 22);
	BOOST_CHECK_EQUAL(unresolved, 6);
}

BOOST_AUTO_TEST_CASE(canonical_29_satellite_product_target_is_reference_invariant)
{
	std::set<SatSys> satellites;
	for (int prn = 1; prn <= 29; prn++)
		satellites.insert(SatSys(E_Sys::GPS, prn));
	const auto reference1 = zhangCanonicalSatelliteProductTarget(
		satellites, SatSys(E_Sys::GPS, 1));
	const auto reference14 = zhangCanonicalSatelliteProductTarget(
		satellites, SatSys(E_Sys::GPS, 14));
	const auto reference29 = zhangCanonicalSatelliteProductTarget(
		satellites, SatSys(E_Sys::GPS, 29));
	BOOST_REQUIRE(reference1.valid && reference14.valid && reference29.valid);
	BOOST_CHECK_EQUAL(reference1.canonicalRank, 28);
	BOOST_CHECK_EQUAL(reference14.canonicalRank, 28);
	BOOST_CHECK_EQUAL(reference29.canonicalRank, 28);
	BOOST_CHECK(reference1.canonicalPrimitive);
	BOOST_CHECK(reference1.canonicalHnf == reference14.canonicalHnf);
	BOOST_CHECK(reference1.canonicalHnf == reference29.canonicalHnf);
	BOOST_CHECK(reference1.canonicalSmithInvariants ==
		reference14.canonicalSmithInvariants);
	BOOST_CHECK(reference1.canonicalSmithInvariants ==
		reference29.canonicalSmithInvariants);
	BOOST_REQUIRE_EQUAL(reference1.canonicalSmithInvariants.size(), 28);
	for (const auto& invariant : reference1.canonicalSmithInvariants)
		BOOST_CHECK_EQUAL(zhangExactAbs(invariant), 1);
	// Every primitive satellite difference is a member, including a direction
	// that uses neither of the three construction references.
	ZhangExactVector g08MinusG09(29);
	g08MinusG09[7] = 1;
	g08MinusG09[8] = -1;
	BOOST_CHECK(zhangIntegerRowLatticeContains(
		reference1.canonicalMatrix, g08MinusG09).contained);
}

BOOST_AUTO_TEST_CASE(only_product_fixed_is_a_formal_pppar_product)
{
	BOOST_CHECK(zhangFormalPppArProductSolution("PRODUCT_FIXED"));
	BOOST_CHECK(!zhangFormalPppArProductSolution("FIXED"));
	BOOST_CHECK(!zhangFormalPppArProductSolution(
		"NETWORK_FIXED_DIAGNOSTIC"));
	BOOST_CHECK(!zhangFormalPppArProductSolution("NETWORK_WL"));
	BOOST_CHECK(!zhangFormalPppArProductSolution("FLOAT"));
}

BOOST_AUTO_TEST_CASE(user_rejects_legacy_network_fixed_ar_claim)
{
	ZhangInternalProduct legacy;
	legacy.solution = "FIXED";
	legacy.ppp_usable = true;
	legacy.pppar_usable = true;
	legacy.ar_valid = true;
	legacy.dual_frequency_ar_valid = true;
	BOOST_CHECK(zhangRejectNonFormalPppArClaim(legacy));
	BOOST_CHECK(legacy.ppp_usable);
	BOOST_CHECK(!legacy.pppar_usable);
	BOOST_CHECK(!legacy.ar_valid);
	BOOST_CHECK(!legacy.dual_frequency_ar_valid);
	BOOST_CHECK_EQUAL(legacy.invalid_reason,
		"NON_PRODUCT_FIXED_AR_CLAIM_REJECTED");

	ZhangInternalProduct formal;
	formal.solution = "PRODUCT_FIXED";
	formal.pppar_usable = true;
	formal.ar_valid = true;
	formal.dual_frequency_ar_valid = true;
	BOOST_CHECK(!zhangRejectNonFormalPppArClaim(formal));
	BOOST_CHECK(formal.pppar_usable);
	BOOST_CHECK(formal.ar_valid);
	BOOST_CHECK(formal.dual_frequency_ar_valid);
}

BOOST_AUTO_TEST_CASE(
	product_relation_named_rows_inherit_accepted_parent_lattice_certificate)
{
	const ZhangExactMatrix parentRows = {
		{1, 1, 0},
		{0, 1, 0}
	};
	const ZhangExactVector parentValues = {7, 3};
	const auto accepted = zhangPromoteNamedCertificateFromAcceptedParent(
		parentRows, parentValues, 3, true);
	BOOST_REQUIRE(accepted.exact);
	BOOST_CHECK_EQUAL(accepted.parentFixedRank, 2);
	BOOST_REQUIRE_EQUAL(accepted.values.size(), 2);
	BOOST_CHECK(accepted.values.at(0) == 4);
	BOOST_CHECK(accepted.values.at(1) == 3);
	BOOST_CHECK(!accepted.values.contains(2));

	// Exact algebra alone cannot bypass the statistical parent gate.
	const auto rejected = zhangPromoteNamedCertificateFromAcceptedParent(
		parentRows, parentValues, 3, false);
	BOOST_CHECK(!rejected.exact);
	BOOST_CHECK(rejected.values.empty());

	// A non-primitive row does not determine a named integer coordinate.
	const auto unsaturated = zhangPromoteNamedCertificateFromAcceptedParent(
		{{2}}, {6}, 1, true);
	BOOST_CHECK(!unsaturated.exact);
	BOOST_CHECK(unsaturated.values.empty());
}

BOOST_AUTO_TEST_CASE(
	product_relation_mixed_lattice_recovers_pair_without_star_coordinate)
{
	// u=z0-z1 is a directly named satellite-pair edge, although neither z0
	// nor z1 relative to the canonical reference is determined.
	int batchTargets = 0;
	int scalarDecompositionsAvoided = 0;
	const auto pairs = zhangRecoverCertifiedPairRelations(
		{{1, -1, 0}}, {4}, 3, true,
		&batchTargets, &scalarDecompositionsAvoided);
	BOOST_REQUIRE_EQUAL(pairs.size(), 1);
	BOOST_CHECK_EQUAL(batchTargets, 6);
	BOOST_CHECK_EQUAL(scalarDecompositionsAvoided, 5);
	BOOST_CHECK_EQUAL(pairs.front().firstNode, 0);
	BOOST_CHECK_EQUAL(pairs.front().secondNode, 1);
	BOOST_CHECK(pairs.front().value == 4);
	const auto stars = zhangPromoteNamedCertificateFromAcceptedParent(
		{{1, -1, 0}}, {4}, 3, true);
	BOOST_CHECK(stars.values.empty());

	// Statistical rejection dominates exact membership.
	BOOST_CHECK(zhangRecoverCertifiedPairRelations(
		{{1, -1, 0}}, {4}, 3, false).empty());
	// A higher-order combination is conditioning evidence, not a pair edge.
	BOOST_CHECK(zhangRecoverCertifiedPairRelations(
		{{1, 1, 1}}, {9}, 3, true).empty());
}

BOOST_AUTO_TEST_CASE(
	product_relation_reliability_forest_uses_only_passed_independent_edges)
{
	std::vector<ZhangPairReliabilityEdge> edges = {
		{0, 1, 1e-5, 0.01},
		{1, 2, 2e-5, 0.02},
		{0, 2, 3e-5, 0.03}, // reliable but closes a cycle
		{2, 3, 2e-3, 0.001}, // precise-looking but fails Perr
		{0, 3, 5e-4, 0.04}
	};
	const auto forest = zhangPairReliabilityForest(4, edges, 1e-3);
	BOOST_REQUIRE_EQUAL(forest.size(), 3);
	BOOST_CHECK_EQUAL(forest[0].firstNode, 0);
	BOOST_CHECK_EQUAL(forest[0].secondNode, 1);
	BOOST_CHECK_EQUAL(forest[1].firstNode, 1);
	BOOST_CHECK_EQUAL(forest[1].secondNode, 2);
	BOOST_CHECK_EQUAL(forest[2].firstNode, 0);
	BOOST_CHECK_EQUAL(forest[2].secondNode, 3);
}

BOOST_AUTO_TEST_CASE(
	product_relation_all_pair_gain_is_reference_invariant)
{
	MatrixXd q(3, 3);
	q << 0.4, 0.1, 0.05,
		 0.1, 0.3, 0.02,
		 0.05, 0.02, 0.2;
	const MatrixXd d = zhangAllPairIncidence(3);
	BOOST_REQUIRE_EQUAL(d.rows(), 6);
	BOOST_REQUIRE_EQUAL(d.cols(), 3);
	const double trace = zhangReferenceInvariantPairTrace(q);
	BOOST_CHECK_CLOSE(trace, (d * q * d.transpose()).trace(), 1e-10);

	// Change star reference from implicit node 3 to node 0.  The new named
	// coordinates are [K1-K0,K2-K0,K3-K0].
	MatrixXd transform(3, 3);
	transform << -1, 1, 0,
		-1, 0, 1,
		-1, 0, 0;
	const MatrixXd changed = transform * q * transform.transpose();
	BOOST_CHECK_CLOSE(trace, zhangReferenceInvariantPairTrace(changed), 1e-9);
}

BOOST_AUTO_TEST_CASE(
	product_relation_exact_conditioning_reports_reference_free_gain)
{
	VectorXd mean(3); mean << 1.1, 2.2, 3.3;
	MatrixXd q = MatrixXd::Identity(3, 3);
	MatrixXd rows(1, 3); rows << 1, -1, 0;
	VectorXd integer(1); integer << -1;
	const auto conditioned = zhangConditionExactProductRows(
		mean, q, rows, integer);
	BOOST_REQUIRE(conditioned.valid);
	BOOST_CHECK_EQUAL(conditioned.effectiveRank, 1);
	BOOST_CHECK_SMALL((rows * conditioned.covariance).norm(), 1e-12);
	BOOST_CHECK_SMALL((rows * conditioned.mean - integer).norm(), 1e-12);
	BOOST_CHECK(zhangReferenceInvariantPairTrace(conditioned.covariance) <
		zhangReferenceInvariantPairTrace(q));
}

BOOST_AUTO_TEST_CASE(component_bridge_gls_aggregates_correlated_edges)
{
	VectorXd edges(3); edges << 4.15, 3.90, 4.05;
	MatrixXd q = MatrixXd::Identity(3, 3) * 0.09;
	const auto bridge = zhangComponentBridgeGls(edges, q);
	BOOST_REQUIRE(bridge.valid);
	BOOST_CHECK_CLOSE(bridge.mean, edges.mean(), 1e-10);
	BOOST_CHECK_CLOSE(bridge.variance, 0.03, 1e-10);
	BOOST_CHECK_EQUAL(bridge.effectiveRank, 3);
	BOOST_CHECK(bridge.residualNis > 0);
}

BOOST_AUTO_TEST_CASE(component_gauge_gls_jointly_recovers_all_component_gauges)
{
	// Three certified components, component zero is the integer datum.  The
	// observations contain redundant correlated edges for c1, c2 and c1-c2.
	MatrixXd design(5, 2);
	design << 1, 0,
		1, 0,
		0, 1,
		0, 1,
		1, -1;
	VectorXd measurements(5);
	measurements << 3.02, 2.98, -1.01, -0.99, 4.01;
	MatrixXd covariance = MatrixXd::Identity(5, 5) * 0.01;
	covariance(0, 1) = covariance(1, 0) = 0.002;
	covariance(2, 3) = covariance(3, 2) = 0.002;
	const auto result = zhangComponentGaugeGls(
		measurements, covariance, design);
	BOOST_REQUIRE(result.valid);
	BOOST_CHECK_EQUAL(result.gaugeRank, 2);
	BOOST_CHECK_EQUAL(result.measurementRank, 5);
	BOOST_REQUIRE_EQUAL(result.mean.size(), 2);
	BOOST_CHECK_SMALL(result.mean(0) - 3, 0.03);
	BOOST_CHECK_SMALL(result.mean(1) + 1, 0.03);
	BOOST_CHECK(result.covariance(0, 0) > 0);
	BOOST_CHECK(result.covariance(1, 1) > 0);
	BOOST_CHECK(result.residualNis >= 0);
}

BOOST_AUTO_TEST_CASE(component_gauge_gls_rejects_unestimable_gauge)
{
	VectorXd measurements(2); measurements << 2, 2;
	MatrixXd covariance = MatrixXd::Identity(2, 2);
	MatrixXd design = MatrixXd::Zero(2, 2);
	design.col(0).setOnes();
	const auto result = zhangComponentGaugeGls(
		measurements, covariance, design);
	BOOST_CHECK(!result.valid);
	BOOST_CHECK_EQUAL(result.gaugeRank, 1);
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_requires_segment_stable_confirmation)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.wideLaneInteger = -4;
	certificate.firstSignalInteger = 17;
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	const auto first = ledger.observe(100, {certificate}, 2);
	BOOST_REQUIRE(first.valid);
	BOOST_CHECK_EQUAL(first.confirmedCertificates, 0);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	BOOST_CHECK(!ledger.certificates().front().active);
	const auto second = ledger.observe(130, {certificate}, 2);
	BOOST_REQUIRE(second.valid);
	BOOST_CHECK_EQUAL(second.confirmedCertificates, 1);
	BOOST_CHECK(ledger.certificates().front().active);
	BOOST_CHECK(ledger.certificates().front().currentAlignmentValid);
}

BOOST_AUTO_TEST_CASE(
	fixed_lag_gauge_confirmation_requires_new_conditional_information)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.wideLaneInteger = -4;
	certificate.firstSignalInteger = 17;
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate, true);
	Matrix2d firstCovariance;
	firstCovariance << 4.0, 0.2,
		0.2, 9.0;
	makeFixedLagTemporalEvidence(certificate, {0, 2, 4}, firstCovariance);

	ProductGaugeCertificateLedger ledger;
	const auto first = ledger.observe(100, {certificate}, 2);
	BOOST_REQUIRE(first.valid);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	BOOST_CHECK_EQUAL(ledger.certificates().front().confirmationEpochs, 1);
	BOOST_CHECK(!ledger.certificates().front().active);

	// Merely evaluating the same cumulative window at a later epoch must not
	// manufacture an independent confirmation.
	const auto duplicate = ledger.observe(130, {certificate}, 2);
	BOOST_REQUIRE(duplicate.valid);
	BOOST_CHECK_EQUAL(duplicate.duplicateEvidenceCandidates, 1);
	BOOST_CHECK_EQUAL(duplicate.independentEvidenceConfirmations, 0);
	BOOST_CHECK_EQUAL(ledger.certificates().front().confirmationEpochs, 1);

	// A new raw factor without a positive relation-space information increment
	// is still not confirmation evidence.
	auto noGain = certificate;
	noGain.factorSequences.push_back(6);
	noGain.factorMeasurementRows += 100;
	const auto zeroGain = ledger.observe(160, {noGain}, 2);
	BOOST_REQUIRE(zeroGain.valid);
	BOOST_CHECK_EQUAL(zeroGain.insufficientInformationGainCandidates, 1);
	BOOST_CHECK_EQUAL(ledger.certificates().front().confirmationEpochs, 1);

	// The same integer may mature only after a genuinely new measurement block
	// produces positive conditional precision in the joint [WL,L1] space.
	auto informative = noGain;
	informative.temporalEvidenceCovariance = 0.5 * firstCovariance;
	const auto confirmed = ledger.observe(190, {informative}, 2);
	BOOST_REQUIRE(confirmed.valid);
	BOOST_CHECK_EQUAL(confirmed.independentEvidenceConfirmations, 1);
	BOOST_CHECK_EQUAL(confirmed.maximumNewFactorCount, 1);
	BOOST_CHECK_EQUAL(confirmed.maximumConditionalInformationRank, 2);
	BOOST_CHECK_CLOSE(confirmed.maximumOverlapFraction, 0.75, 1e-9);
	BOOST_CHECK_GT(confirmed.maximumConditionalInformationGain, 0);
	BOOST_CHECK(ledger.certificates().front().active);
	BOOST_CHECK_EQUAL(ledger.certificates().front().confirmationEpochs, 2);
	BOOST_CHECK_EQUAL(ledger.certificates().front().factorSequences.size(), 4);
}

BOOST_AUTO_TEST_CASE(
	fixed_lag_gauge_confirmation_rejects_non_nested_factor_windows)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 8);
	certificate.reference = SatSys(E_Sys::GPS, 9);
	certificate.satellitePhaseSegments = {"G08-L1-SEG0", "G08-L2-SEG0"};
	certificate.referencePhaseSegments = {"G09-L1-SEG0", "G09-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate, true);
	Matrix2d covariance = Matrix2d::Identity();
	makeFixedLagTemporalEvidence(certificate, {0, 2, 4}, covariance);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 2).valid);

	auto shifted = certificate;
	shifted.factorSequences = {2, 4, 6};
	shifted.temporalEvidenceCovariance = 0.5 * covariance;
	const auto rejected = ledger.observe(130, {shifted}, 2);
	BOOST_REQUIRE(rejected.valid);
	BOOST_CHECK_EQUAL(rejected.nonNestedEvidenceCandidates, 1);
	BOOST_CHECK_EQUAL(rejected.independentEvidenceConfirmations, 0);
	BOOST_CHECK(!ledger.certificates().front().active);
	BOOST_CHECK_EQUAL(ledger.certificates().front().confirmationEpochs, 1);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_pending_recheck_accepts_same_correlated_dual_integer)
{
	Vector2d mean; mean << 7.0001, -3.9998;
	Matrix2d covariance;
	covariance << 1e-4, 5e-5,
		5e-5, 2e-4;
	const auto recheck = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, ZhangExactInteger(7), ZhangExactInteger(-4),
		1e-3, 1e-6);
	BOOST_REQUIRE(recheck.valid);
	BOOST_CHECK(recheck.sameInteger);
	BOOST_CHECK(recheck.reliable);
	BOOST_CHECK(!recheck.alternativeReliable);
	BOOST_CHECK_EQUAL(recheck.nearestWideLaneInteger, 7);
	BOOST_CHECK_EQUAL(recheck.nearestFirstSignalInteger, -4);
	BOOST_CHECK_LE(recheck.jointFailureProbability, 1e-3);
	BOOST_CHECK_LE(recheck.nis, recheck.nisThreshold);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_projection_maps_reference_invariant_q_to_wl_l1_posterior)
{
	const SatSys g02(E_Sys::GPS, 2);
	const SatSys g03(E_Sys::GPS, 3);
	const SatSys g04(E_Sys::GPS, 4);
	ProductGaugeCertificate certificate;
	certificate.productRelationQ = {{g03, 1}, {g04, -1}};
	const std::map<SatSys, int> coordinate = {{g03, 0}, {g04, 1}};
	VectorXd mean(4); mean << 10, 4, 7, 1;
	const MatrixXd covariance = MatrixXd::Identity(4, 4);
	const auto projection = zhangProjectProductGaugeToPosterior(
		certificate, g02, coordinate, 2, mean, covariance);
	BOOST_REQUIRE(projection.valid);
	BOOST_REQUIRE_EQUAL(projection.productRow.size(), 2);
	BOOST_CHECK_EQUAL(projection.productRow[0], 1);
	BOOST_CHECK_EQUAL(projection.productRow[1], -1);
	BOOST_CHECK_SMALL(projection.floatingValue(0), 1e-12);
	BOOST_CHECK_CLOSE(projection.floatingValue(1), 6.0, 1e-12);
	BOOST_CHECK_CLOSE(projection.covariance(0, 0), 4.0, 1e-12);
	BOOST_CHECK_CLOSE(projection.covariance(0, 1), 2.0, 1e-12);
	BOOST_CHECK_CLOSE(projection.covariance(1, 1), 2.0, 1e-12);

	certificate.productRelationQ = {{g03, 1}, {g02, -1}};
	const auto implicitReference = zhangProjectProductGaugeToPosterior(
		certificate, g02, coordinate, 2, mean, covariance);
	BOOST_REQUIRE(implicitReference.valid);
	BOOST_CHECK_EQUAL(implicitReference.productRow[0], 1);
	BOOST_CHECK_EQUAL(implicitReference.productRow[1], 0);
	BOOST_CHECK_CLOSE(implicitReference.floatingValue(0), 3.0, 1e-12);
	BOOST_CHECK_CLOSE(implicitReference.floatingValue(1), 10.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_projection_rejects_noninvariant_or_unavailable_q)
{
	const SatSys g02(E_Sys::GPS, 2);
	const SatSys g03(E_Sys::GPS, 3);
	const SatSys g05(E_Sys::GPS, 5);
	const std::map<SatSys, int> coordinate = {{g03, 0}};
	const VectorXd mean = VectorXd::Zero(2);
	const MatrixXd covariance = MatrixXd::Identity(2, 2);
	ProductGaugeCertificate certificate;
	certificate.productRelationQ = {{g03, 1}};
	const auto nonInvariant = zhangProjectProductGaugeToPosterior(
		certificate, g02, coordinate, 1, mean, covariance);
	BOOST_CHECK(!nonInvariant.valid);
	BOOST_CHECK_EQUAL(nonInvariant.failureReason,
		"PRODUCT_GAUGE_RELATION_NOT_REFERENCE_INVARIANT");

	certificate.productRelationQ = {{g05, 1}, {g02, -1}};
	const auto unavailable = zhangProjectProductGaugeToPosterior(
		certificate, g02, coordinate, 1, mean, covariance);
	BOOST_CHECK(!unavailable.valid);
	BOOST_CHECK_EQUAL(unavailable.failureReason,
		"PRODUCT_GAUGE_SATELLITE_NOT_IN_CURRENT_BASIS");
}

BOOST_AUTO_TEST_CASE(
	product_gauge_pending_recheck_honours_semidefinite_null_space)
{
	Vector2d mean; mean << 7.00001, -3.99999;
	Matrix2d covariance;
	covariance << 1e-4, 1e-4,
		1e-4, 1e-4;
	const auto consistent = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, ZhangExactInteger(7), ZhangExactInteger(-4),
		1e-3, 1e-6);
	BOOST_REQUIRE(consistent.valid);
	BOOST_CHECK(consistent.reliable);

	mean << 7.00001, -4.00001;
	const auto inconsistent = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, ZhangExactInteger(7), ZhangExactInteger(-4),
		1e-3, 1e-6);
	BOOST_CHECK(!inconsistent.valid);
	BOOST_CHECK(!inconsistent.reliable);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_pending_recheck_rejects_reliable_changed_integer)
{
	Vector2d mean; mean << 8.0001, -3.9998;
	Matrix2d covariance;
	covariance << 1e-4, 5e-5,
		5e-5, 2e-4;
	const auto recheck = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, ZhangExactInteger(7), ZhangExactInteger(-4),
		1e-3, 1e-6);
	BOOST_REQUIRE(recheck.valid);
	BOOST_CHECK(!recheck.sameInteger);
	BOOST_CHECK(!recheck.reliable);
	BOOST_CHECK(recheck.alternativeReliable);
	BOOST_CHECK_EQUAL(recheck.nearestWideLaneInteger, 8);
	BOOST_CHECK_EQUAL(recheck.nearestFirstSignalInteger, -4);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_pending_recheck_handles_deterministic_posterior)
{
	Vector2d mean; mean << 7, -4;
	const Matrix2d covariance = Matrix2d::Zero();
	const auto accepted = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, ZhangExactInteger(7), ZhangExactInteger(-4),
		1e-3, 1e-6);
	BOOST_REQUIRE(accepted.valid);
	BOOST_CHECK(accepted.deterministic);
	BOOST_CHECK(accepted.reliable);

	mean(0) = 8;
	const auto changed = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, ZhangExactInteger(7), ZhangExactInteger(-4),
		1e-3, 1e-6);
	BOOST_REQUIRE(changed.valid);
	BOOST_CHECK(changed.deterministic);
	BOOST_CHECK(!changed.reliable);
	BOOST_CHECK(changed.alternativeReliable);
}

BOOST_AUTO_TEST_CASE(
	integer_candidate_nis_accepts_only_consistent_deterministic_certificate)
{
	const Matrix2d covariance = Matrix2d::Zero();
	const auto accepted = assessZhangIntegerCandidateNis(
		Vector2d::Zero(), covariance, 1e-6);
	BOOST_REQUIRE(accepted.valid);
	BOOST_CHECK(accepted.deterministic);
	BOOST_CHECK_EQUAL(accepted.rank, 0);
	BOOST_CHECK_SMALL(accepted.nis, 1e-15);
	BOOST_CHECK_SMALL(accepted.threshold, 1e-15);

	Vector2d conflicting = Vector2d::Zero();
	conflicting(1) = 1e-4;
	const auto rejected = assessZhangIntegerCandidateNis(
		conflicting, covariance, 1e-6);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK(!rejected.deterministic);
}

BOOST_AUTO_TEST_CASE(
	fixed_lag_raw_integer_proposal_requires_current_posterior_reauthorization)
{
	// The smoothed raw-factor marginal can be extremely confident in an old
	// integer.  That makes it useful for proposal generation, but never gives it
	// authority over the current physical posterior.
	const Vector2d proposedInteger(12, -7);
	const Matrix2d rawCovariance = 1e-4 * Matrix2d::Identity();
	const Vector2d rawMean(12.001, -6.999);
	const auto rawEvidence = assessZhangIntegerCandidateNis(
		proposedInteger - rawMean, rawCovariance, 1e-6);
	BOOST_REQUIRE(rawEvidence.valid);
	BOOST_CHECK_LE(rawEvidence.nis, rawEvidence.threshold);

	const Matrix2d currentCovariance = 0.01 * Matrix2d::Identity();
	const Vector2d currentCompatibleMean(12.02, -7.01);
	const auto reauthorised = assessZhangIntegerCandidateNis(
		proposedInteger - currentCompatibleMean, currentCovariance, 1e-6);
	BOOST_REQUIRE(reauthorised.valid);
	BOOST_CHECK_LE(reauthorised.nis, reauthorised.threshold);

	// A one-cycle identity/version error remains a hard rejection even though
	// the raw window itself still strongly supports the stale integer.
	const Vector2d currentChangedMean(13.0, -7.0);
	const auto rejected = assessZhangIntegerCandidateNis(
		proposedInteger - currentChangedMean, currentCovariance, 1e-6);
	BOOST_REQUIRE(rejected.valid);
	BOOST_CHECK_GT(rejected.nis, rejected.threshold);
}

BOOST_AUTO_TEST_CASE(
	product_pair_provenance_survives_exact_union_catalogue_rebuild)
{
	ZhangProductIntegerConstraintSet temporal;
	ZhangCertifiedPairRelation temporalPair;
	temporalPair.firstNode = 2;
	temporalPair.secondNode = 0;
	temporalPair.fromTemporalLedger = true;
	temporal.dualFrequencyCertifiedPairs.push_back(temporalPair);

	ZhangProductIntegerConstraintSet gauge;
	ZhangCertifiedPairRelation gaugePair;
	gaugePair.firstNode = 3;
	gaugePair.secondNode = 4;
	gaugePair.fromProductGaugeLedger = true;
	gauge.dualFrequencyCertifiedPairs.push_back(gaugePair);

	ZhangProductIntegerConstraintSet rebuilt;
	ZhangCertifiedPairRelation reversedTemporal;
	reversedTemporal.firstNode = 0;
	reversedTemporal.secondNode = 2;
	rebuilt.dualFrequencyCertifiedPairs.push_back(reversedTemporal);
	ZhangCertifiedPairRelation rebuiltGauge;
	rebuiltGauge.firstNode = 4;
	rebuiltGauge.secondNode = 3;
	rebuilt.dualFrequencyCertifiedPairs.push_back(rebuiltGauge);
	ZhangCertifiedPairRelation unrelated;
	unrelated.firstNode = 1;
	unrelated.secondNode = 4;
	rebuilt.dualFrequencyCertifiedPairs.push_back(unrelated);

	zhangPropagateProductPairProvenance(rebuilt, {&temporal, &gauge});
	BOOST_CHECK(rebuilt.dualFrequencyCertifiedPairs[0].fromTemporalLedger);
	BOOST_CHECK(!rebuilt.dualFrequencyCertifiedPairs[0].fromProductGaugeLedger);
	BOOST_CHECK(!rebuilt.dualFrequencyCertifiedPairs[1].fromTemporalLedger);
	BOOST_CHECK(rebuilt.dualFrequencyCertifiedPairs[1].fromProductGaugeLedger);
	BOOST_CHECK(!rebuilt.dualFrequencyCertifiedPairs[2].fromTemporalLedger);
	BOOST_CHECK(!rebuilt.dualFrequencyCertifiedPairs[2].fromProductGaugeLedger);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_pending_recheck_preserves_original_admission_budget)
{
	ProductGaugeCertificate certificate;
	BOOST_CHECK_CLOSE(zhangProductGaugePosteriorFailureBudget(
		certificate, 1e-3), 1e-3, 1e-9);
	certificate.failureProbabilityBudget = 2e-5;
	BOOST_CHECK_CLOSE(zhangProductGaugePosteriorFailureBudget(
		certificate, 1e-3), 2e-5, 1e-9);
	certificate.failureProbabilityBudget = 0;
	BOOST_CHECK_CLOSE(zhangProductGaugePosteriorFailureBudget(
		certificate, 1e-3), 1e-12, 1e-9);
	certificate.failureProbabilityBudget =
		std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_CLOSE(zhangProductGaugePosteriorFailureBudget(
		certificate, 1e-3), 1e-3, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_pending_certificate_matures_without_candidate_reselection)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.wideLaneInteger = 7;
	certificate.firstSignalInteger = -4;
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	certificate.failureProbabilityBudget = 2e-5;
	makeProductGaugeEvidenceComplete(certificate, true);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 2).valid);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	BOOST_CHECK(ledger.certificates().front().state ==
		ProductGaugeCertificateState::CONFIRMING);

	const auto pending = ledger.certificates().front();
	Vector2d mean;
	mean << pending.wideLaneInteger.convert_to<double>() + 1e-5,
		pending.firstSignalInteger.convert_to<double>() - 1e-5;
	Matrix2d covariance;
	covariance << 1e-6, 2e-7,
		2e-7, 1e-6;
	const auto recheck = zhangRecheckProductGaugeOnPosterior(
		mean, covariance, pending.wideLaneInteger,
		pending.firstSignalInteger,
		zhangProductGaugePosteriorFailureBudget(pending, 1e-3), 1e-6);
	BOOST_REQUIRE(recheck.reliable);
	BOOST_REQUIRE(ledger.observe(130, {pending}, 2).valid);
	BOOST_REQUIRE(ledger.certificates().front().active);
	BOOST_CHECK_EQUAL(ledger.certificates().front().confirmationEpochs, 2);
	BOOST_CHECK_CLOSE(ledger.certificates().front().failureProbabilityBudget,
		2e-5, 1e-9);

	BOOST_CHECK(ledger.invalidate(ledger.certificates().front(), 160));
	BOOST_CHECK(!ledger.certificates().front().active);
	BOOST_CHECK(ledger.certificates().front().state ==
		ProductGaugeCertificateState::SUSPENDED);
	BOOST_CHECK(ledger.certificates().front().inactiveReason ==
		ProductGaugeCertificateInactiveReason::CONFLICT);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_admission_family_budget_counts_joint_family_once)
{
	ProductGaugeCertificate first;
	first.system = E_Sys::GPS;
	first.firstObservable = E_ObsCode::L1C;
	first.secondObservable = E_ObsCode::L2W;
	first.satellite = SatSys(E_Sys::GPS, 3);
	first.reference = SatSys(E_Sys::GPS, 2);
	first.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	first.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	first.admissionFamilyId = "FAMILY-A";
	first.admissionFamilyFailureProbabilityBudget = 2e-4;
	first.failureProbabilityBudget = 5e-5;
	zhangCanonicaliseProductGaugeCertificate(first);

	ProductGaugeCertificate second = first;
	second.satellite = SatSys(E_Sys::GPS, 4);
	second.satellitePhaseSegments = {"G04-L1-SEG0", "G04-L2-SEG0"};
	second.productRelationQ.clear();
	second.phaseSegmentsBySatellite.clear();
	zhangCanonicaliseProductGaugeCertificate(second);

	ProductGaugeCertificate third = first;
	third.satellite = SatSys(E_Sys::GPS, 5);
	third.satellitePhaseSegments = {"G05-L1-SEG0", "G05-L2-SEG0"};
	third.productRelationQ.clear();
	third.phaseSegmentsBySatellite.clear();
	third.admissionFamilyId = "FAMILY-B";
	third.admissionFamilyFailureProbabilityBudget = 3e-4;
	zhangCanonicaliseProductGaugeCertificate(third);

	const auto grouped = zhangProductGaugeAdmissionFamilyBudgetAudit(
		{&first, &second, &third});
	BOOST_REQUIRE(grouped.valid);
	BOOST_CHECK_EQUAL(grouped.certificates, 3);
	BOOST_CHECK_EQUAL(grouped.families, 2);
	BOOST_CHECK_EQUAL(grouped.legacyFallbackFamilies, 0);
	BOOST_CHECK_CLOSE(grouped.failureProbability, 5e-4, 1e-9);

	second.admissionFamilyFailureProbabilityBudget = 4e-4;
	const auto conservativeMaximum =
		zhangProductGaugeAdmissionFamilyBudgetAudit({&first, &second});
	BOOST_REQUIRE(conservativeMaximum.valid);
	BOOST_CHECK_EQUAL(conservativeMaximum.families, 1);
	BOOST_CHECK_CLOSE(
		conservativeMaximum.failureProbability, 4e-4, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_subset_maximises_coverage_before_redundant_rank)
{
	const std::vector<double> probabilities{6e-4, 4e-4, 3e-4};
	const std::vector<int> rankGain{4, 3, 2};
	const std::vector<int> satelliteGain{5, 3, 4};
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		probabilities, 1e-3, [&](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid = true;
			for (const int family : families)
			{
				score.dualGraphRank += rankGain[family];
				score.certifiedSatelliteCount += satelliteGain[family];
				score.conditioningRank += 2 * rankGain[family];
			}
			return score;
		});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(selected.exactEnumeration);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{0, 2}));
	BOOST_CHECK_EQUAL(selected.score.dualGraphRank, 6);
	BOOST_CHECK_EQUAL(selected.score.certifiedSatelliteCount, 9);
	BOOST_CHECK_CLOSE(selected.failureProbability, 9e-4, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_subset_uses_safe_temporal_overlap_as_tie_breaker)
{
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>{4e-4, 4e-4}, 4e-4,
		[](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid = families.size() == 1;
			if (!score.valid) return score;
			score.largestCertifiedComponentSize = 4;
			score.certifiedSatelliteCount = 4;
			score.dualGraphRank = 3;
			score.certifiedComponentCount = 1;
			score.conditioningRank = 6;
			score.continuityPriority = families.front() == 1 ? 3 : 1;
			return score;
		});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{1}));
	BOOST_CHECK_EQUAL(selected.score.continuityPriority, 3);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_subset_can_replace_expensive_current_family)
{
	// Family zero models the atomic current-epoch certificate set.  The two
	// historical families jointly certify more product directions at the same
	// global Perr and must therefore be allowed to replace, not merely append to,
	// the current set.
	const std::vector<double> probabilities{9e-4, 5e-4, 5e-4};
	const std::vector<int> rankGain{2, 3, 3};
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		probabilities, 1e-3, [&](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid = true;
			for (const int family : families)
				score.dualGraphRank += rankGain[family];
			score.certifiedSatelliteCount = score.dualGraphRank + 1;
			score.conditioningRank = 2 * score.dualGraphRank;
			return score;
		});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{1, 2}));
	BOOST_CHECK_EQUAL(selected.score.dualGraphRank, 6);
	BOOST_CHECK_CLOSE(selected.failureProbability, 1e-3, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_subset_prefers_one_larger_common_component)
{
	// Family 0 represents a connected four-satellite product component.  Family
	// 1 has more total rank and touches more satellites, but only as two
	// disconnected three-satellite islands.  PPP-AR users cannot share an
	// integer datum across those islands, so the connected component must win.
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>{5e-4, 5e-4}, 5e-4,
		[](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			if (families == std::vector<int>{0})
			{
				score.valid = true;
				score.largestCertifiedComponentSize = 4;
				score.certifiedComponentCount = 1;
				score.dualGraphRank = 3;
				score.certifiedSatelliteCount = 4;
			}
			else if (families == std::vector<int>{1})
			{
				score.valid = true;
				score.largestCertifiedComponentSize = 3;
				score.certifiedComponentCount = 2;
				score.dualGraphRank = 4;
				score.certifiedSatelliteCount = 6;
			}
			return score;
		});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{0}));
	BOOST_CHECK_EQUAL(selected.score.largestCertifiedComponentSize, 4);
	BOOST_CHECK_EQUAL(selected.score.certifiedComponentCount, 1);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_subset_evaluates_complementary_joint_family_sets)
{
	const std::vector<double> probabilities{4e-4, 4e-4};
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		probabilities, 1e-3, [&](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			// Neither signal-only family is broadcastable alone, but their joint
			// exact union closes a dual-frequency product direction.
			score.valid = families.size() == 2;
			score.dualGraphRank = score.valid ? 1 : 0;
			score.certifiedSatelliteCount = score.valid ? 2 : 0;
			score.conditioningRank = score.valid ? 2 : 0;
			return score;
		});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{0, 1}));
	BOOST_CHECK_EQUAL(selected.evaluatedSubsets, 3);
	BOOST_CHECK_EQUAL(selected.budgetFeasibleSubsets, 3);
	BOOST_CHECK_EQUAL(selected.evaluatedUniqueSubsets, 3);
	BOOST_CHECK_EQUAL(selected.searchableNonEmptySubsets, 3);
	BOOST_CHECK_EQUAL(selected.prunedNotTestedSubsets, 0);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_bounded_selector_distinguishes_pruned_from_tested_failure)
{
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>(10, 1e-4), 5e-4,
		[](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid = !families.empty() && families.front() == 0;
			score.dualGraphRank = score.valid ? families.size() : 0;
			score.certifiedSatelliteCount = score.valid
				? score.dualGraphRank + 1 : 0;
			return score;
		}, 0, 2, {0});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(!selected.exactEnumeration);
	BOOST_CHECK_GT(selected.evaluatedUniqueSubsets, 0);
	BOOST_CHECK_EQUAL(selected.evaluatedSubsets,
		selected.evaluatedUniqueSubsets);
	BOOST_CHECK_EQUAL(selected.searchableNonEmptySubsets, 1023);
	BOOST_CHECK_GT(selected.prunedNotTestedSubsets, 0);
	BOOST_CHECK_LT(selected.prunedNotTestedSubsets,
		selected.searchableNonEmptySubsets);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_beam_preserves_late_complementary_frontier)
{
	const std::vector<double> probabilities(10, 1e-4);
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		probabilities, 5e-4, [&](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			const bool wide = std::find(families.begin(), families.end(), 8)
				!= families.end();
			const bool first = std::find(families.begin(), families.end(), 9)
				!= families.end();
			const int distractors = std::count_if(
				families.begin(), families.end(), [](const int family)
				{
					return family < 8;
				});
			score.wideLaneRank = wide;
			score.firstSignalRank = first;
			score.valid = distractors > 0 || (wide && first);
			score.dualGraphRank = wide && first
				? 20 + distractors : distractors;
			score.certifiedSatelliteCount = score.valid
				? score.dualGraphRank + 1 : 0;
			score.conditioningRank = score.dualGraphRank * 2;
			return score;
		}, 0, 4);
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(!selected.exactEnumeration);
	BOOST_CHECK(std::find(selected.selectedFamilyIndices.begin(),
		selected.selectedFamilyIndices.end(), 8) !=
		selected.selectedFamilyIndices.end());
	BOOST_CHECK(std::find(selected.selectedFamilyIndices.begin(),
		selected.selectedFamilyIndices.end(), 9) !=
		selected.selectedFamilyIndices.end());
	BOOST_CHECK_GE(selected.score.dualGraphRank, 20);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_beam_width_one_keeps_one_nonempty_candidate)
{
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>{1e-4}, 1e-3,
		[](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid = families == std::vector<int>{0};
			score.dualGraphRank = score.valid;
			score.certifiedSatelliteCount = score.valid ? 2 : 0;
			score.conditioningRank = score.valid ? 2 : 0;
			return score;
		}, 0, 1);
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(!selected.exactEnumeration);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{0}));
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_beam_width_one_prefers_verified_output)
{
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>{1e-4, 1e-4}, 1e-3,
		[](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid = std::find(families.begin(), families.end(), 0) !=
				families.end();
			score.wideLaneRank = std::find(
				families.begin(), families.end(), 1) != families.end();
			score.dualGraphRank = score.valid;
			score.certifiedSatelliteCount = score.valid ? 2 : 0;
			score.conditioningRank = score.valid ? 2 : 0;
			return score;
		}, 0, 1);
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(!selected.exactEnumeration);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{0}));
}

BOOST_AUTO_TEST_CASE(
	product_gauge_ledger_keeps_first_joint_admission_family_immutable)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	certificate.failureProbabilityBudget = 3e-4;
	certificate.admissionFamilyId = "FIRST-JOINT-FAMILY";
	certificate.admissionFamilyFailureProbabilityBudget = 3e-4;
	makeProductGaugeEvidenceComplete(certificate, true);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 2).valid);

	auto later = certificate;
	later.failureProbabilityBudget = 1e-4;
	later.admissionFamilyId = "LATER-JOINT-FAMILY";
	later.admissionFamilyFailureProbabilityBudget = 1e-4;
	BOOST_REQUIRE(ledger.observe(130, {later}, 2).valid);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	const auto& stored = ledger.certificates().front();
	BOOST_CHECK_EQUAL(stored.admissionFamilyId, "FIRST-JOINT-FAMILY");
	BOOST_CHECK_CLOSE(
		stored.admissionFamilyFailureProbabilityBudget, 3e-4, 1e-9);
	BOOST_CHECK_CLOSE(stored.failureProbabilityBudget, 1e-4, 1e-9);
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_rejects_unresolved_segment_sentinels)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 8);
	certificate.reference = SatSys(E_Sys::GPS, 9);
	certificate.satellitePhaseSegments = {"UNRESOLVED", "G08-L2W-SEG0"};
	certificate.referencePhaseSegments = {"G09-L1C-SEG0", "G09-L2W-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	const auto unresolved = ledger.observe(100, {certificate}, 1);
	BOOST_CHECK(!unresolved.valid);
	BOOST_CHECK_EQUAL(unresolved.rejectedCandidates, 1);
	BOOST_CHECK_EQUAL(unresolved.failureReason,
		"PRODUCT_GAUGE_CERTIFICATE_SEGMENT_UNRESOLVED");
	BOOST_CHECK(ledger.certificates().empty());

	certificate.satellitePhaseSegments[0] = "G08-L1C-A";
	const auto malformed = ledger.observe(130, {certificate}, 1);
	BOOST_CHECK(!malformed.valid);
	BOOST_CHECK_EQUAL(malformed.failureReason,
		"PRODUCT_GAUGE_CERTIFICATE_SEGMENT_UNRESOLVED");
	BOOST_CHECK(ledger.certificates().empty());
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_identity_includes_canonical_segments)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 9);
	certificate.reference = SatSys(E_Sys::GPS, 8);
	certificate.productRelationQ[certificate.satellite] = 1;
	certificate.productRelationQ[certificate.reference] = -1;
	certificate.satellitePhaseSegments = {"G09-L1C-SEG0", "G09-L2W-SEG0"};
	certificate.referencePhaseSegments = {"G08-L1C-SEG0", "G08-L2W-SEG0"};
	certificate.phaseSegmentsBySatellite[certificate.satellite] =
		certificate.satellitePhaseSegments;
	certificate.phaseSegmentsBySatellite[certificate.reference] =
		certificate.referencePhaseSegments;
	zhangCanonicaliseProductGaugeCertificate(certificate);
	const auto identity = zhangProductGaugeCertificateIdentity(certificate);
	BOOST_CHECK(identity.find("G08-L1C-SEG0") != std::string::npos);
	BOOST_CHECK(identity.find("G09-L2W-SEG0") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_retires_on_phase_segment_change)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.wideLaneInteger = 2;
	certificate.firstSignalInteger = 9;
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	BOOST_REQUIRE(ledger.certificates().front().active);
	certificate.satellitePhaseSegments[0] = "G03-L1-SEG1";
	const auto changed = ledger.observe(130, {certificate}, 1);
	BOOST_REQUIRE(changed.valid);
	BOOST_CHECK_EQUAL(changed.retiredSegmentCertificates, 1);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 2);
	BOOST_CHECK(!ledger.certificates().front().active);
	BOOST_CHECK(ledger.certificates().back().active);
}

BOOST_AUTO_TEST_CASE(product_gauge_ledger_persists_arbitrary_canonical_q)
{
	const SatSys g02(E_Sys::GPS, 2);
	const SatSys g03(E_Sys::GPS, 3);
	const SatSys g04(E_Sys::GPS, 4);
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = g02;
	certificate.reference = g04;
	certificate.productRelationQ = {{g02, 1}, {g03, 1}, {g04, -2}};
	certificate.wideLaneInteger = 5;
	certificate.firstSignalInteger = 19;
	certificate.phaseSegmentsBySatellite = {
		{g02, {"G02-L1-SEG0", "G02-L2-SEG0"}},
		{g03, {"G03-L1-SEG0", "G03-L2-SEG0"}},
		{g04, {"G04-L1-SEG0", "G04-L2-SEG0"}}};
	certificate.satellitePhaseSegments =
		certificate.phaseSegmentsBySatellite[g02];
	certificate.referencePhaseSegments =
		certificate.phaseSegmentsBySatellite[g04];
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	BOOST_CHECK_EQUAL(ledger.activeRank(), 1);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	const auto& active = ledger.certificates().front();
	BOOST_CHECK(active.productRelationQ == certificate.productRelationQ);
	BOOST_CHECK(zhangProductGaugeCertificateMatchesCurrentSegments(
		active, E_ObsCode::L1C, E_ObsCode::L2W,
		certificate.phaseSegmentsBySatellite));
	auto changed = certificate.phaseSegmentsBySatellite;
	changed[g03][1] = "G03-L2-SEG1";
	BOOST_CHECK(!zhangProductGaugeCertificateMatchesCurrentSegments(
		active, E_ObsCode::L1C, E_ObsCode::L2W, changed));
}

BOOST_AUTO_TEST_CASE(component_gauge_product_row_includes_certified_offsets)
{
	// Component anchors have internal potentials [2,-3,5] relative to their
	// own component gauges.  Fixed row 2*(c1-c0)-(c2-c0)=7 implies
	// 2*K1-K2-K0 = 7 + 2*(-3-2) - (5-2) = -6.
	const auto mapped = zhangComponentGaugeToProductRow(
		{2, -1}, {0, 2, 3}, {2, -3, 5}, 4, 7);
	BOOST_REQUIRE(mapped.valid);
	const ZhangExactVector expected = {-1, 0, 2, -1};
	BOOST_CHECK(mapped.row == expected);
	BOOST_CHECK_EQUAL(mapped.value, -6);
}

BOOST_AUTO_TEST_CASE(component_gauge_product_row_supports_implicit_reference)
{
	// Datum anchor is the implicit canonical reference at index dimension.
	const auto mapped = zhangComponentGaugeToProductRow(
		{1}, {3, 1}, {0, 4}, 3, -2);
	BOOST_REQUIRE(mapped.valid);
	const ZhangExactVector expected = {0, 1, 0};
	BOOST_CHECK(mapped.row == expected);
	BOOST_CHECK_EQUAL(mapped.value, 2);
}

BOOST_AUTO_TEST_CASE(dual_component_gauge_wl_l1_is_exact)
{
	// [WL,L1] -> [L1,L2] is unimodular: L2=L1-WL, with no rounding.
	const ZhangExactInteger wl = -7;
	const ZhangExactInteger l1 = 19;
	const ZhangExactInteger l2 = l1 - wl;
	BOOST_CHECK_EQUAL(l2, 26);
	BOOST_CHECK_EQUAL(l1 - l2, wl);
}

BOOST_AUTO_TEST_CASE(component_gauge_rows_map_exactly_to_product_lattice)
{
	const auto mapped = zhangComponentGaugeToProductRow(
		{1, -1}, {0, 2, 4}, {3, -2, 5}, 5, 11);
	BOOST_REQUIRE(mapped.valid);
	BOOST_CHECK(zhangExactPrimitiveRowLattice({mapped.row}, 5));
	BOOST_CHECK(zhangIntegerRowLatticeContains(
		zhangExactIdentityMatrix(5), mapped.row).contained);
}

BOOST_AUTO_TEST_CASE(component_gauge_union_increases_exact_rank)
{
	const ZhangExactMatrix target = zhangExactIdentityMatrix(3);
	const auto unionAudit = zhangExactCertifiedUnionAudit(
		target, {{1, 0, 0}}, {4}, {{0, 1, 0}}, {-2});
	BOOST_REQUIRE(unionAudit.consistent);
	BOOST_CHECK_EQUAL(unionAudit.heldRank, 1);
	BOOST_CHECK_EQUAL(unionAudit.newlyFixedRank, 1);
	BOOST_CHECK_EQUAL(unionAudit.combinedCertifiedRank, 2);
}

BOOST_AUTO_TEST_CASE(redundant_component_edge_does_not_increase_rank)
{
	const ZhangExactMatrix target = zhangExactIdentityMatrix(2);
	const auto unionAudit = zhangExactCertifiedUnionAudit(
		target, {{1, -1}}, {3}, {{2, -2}}, {6});
	BOOST_REQUIRE(unionAudit.consistent);
	BOOST_CHECK_EQUAL(unionAudit.heldRank, 1);
	BOOST_CHECK_EQUAL(unionAudit.newlyFixedRank, 0);
	BOOST_CHECK_EQUAL(unionAudit.combinedCertifiedRank, 1);
}

BOOST_AUTO_TEST_CASE(component_gauge_joint_covariance_uses_l1_l2_cross_terms)
{
	MatrixXd design = MatrixXd::Identity(2, 2);
	VectorXd measurement(2); measurement << 3, -1;
	MatrixXd correlated(2, 2); correlated << 1, 0.4, 0.4, 1;
	MatrixXd independent = MatrixXd::Identity(2, 2);
	const auto joint = zhangComponentGaugeGls(measurement, correlated, design);
	const auto diagonal = zhangComponentGaugeGls(measurement, independent, design);
	BOOST_REQUIRE(joint.valid);
	BOOST_REQUIRE(diagonal.valid);
	BOOST_CHECK_SMALL((joint.covariance - correlated).norm(), 1e-12);
	BOOST_CHECK_SMALL((diagonal.covariance - independent).norm(), 1e-12);
	BOOST_CHECK_SMALL((joint.covariance - diagonal.covariance)(0, 1) - 0.4, 1e-12);
}

BOOST_AUTO_TEST_CASE(component_gauge_receiver_group_leave_one_out_refits_gls)
{
	MatrixXd design(4, 2);
	design << 1, 0,
		1, 0,
		0, 1,
		0, 1;
	VectorXd measurement(4); measurement << 1, 3, 4, 4;
	MatrixXd covariance = MatrixXd::Identity(4, 4);
	VectorXd functional(2); functional << 1, -1;
	const auto full = zhangComponentGaugeGls(
		measurement, covariance, design);
	BOOST_REQUIRE(full.valid);
	BOOST_CHECK_SMALL(functional.dot(full.mean) + 2, 1e-12);
	const auto leaveOut = zhangComponentGaugeLeaveOneGroupOut(
		measurement, covariance, design, {0}, functional);
	BOOST_REQUIRE(leaveOut.valid);
	BOOST_CHECK(leaveOut.estimable);
	BOOST_CHECK_EQUAL(leaveOut.removedCount, 1);
	BOOST_CHECK_EQUAL(leaveOut.retainedCount, 3);
	BOOST_CHECK_SMALL(leaveOut.mean + 1, 1e-12);
	BOOST_CHECK_SMALL(leaveOut.variance - 1.5, 1e-12);
	const auto unobservable = zhangComponentGaugeLeaveOneGroupOut(
		measurement, covariance, design, {0, 1}, functional);
	BOOST_CHECK(!unobservable.valid);
	BOOST_CHECK(!unobservable.estimable);
}

BOOST_AUTO_TEST_CASE(
	component_gauge_rejected_candidate_keeps_p1_physical_diagnostics)
{
	// Frozen negative-path block: the dual gauge is estimable and its integer
	// direction is exact, but the deliberately poor fractional mean must fail
	// the 1e-3 reliability gate.  Every P1 diagnostic remains defined before
	// that gate is applied.
	MatrixXd design(4, 2);
	design << 1, 0,
		1, 0,
		0, 1,
		0, 1;
	VectorXd measurement(4); measurement << 0.20, 0.20, 0.70, 0.70;
	MatrixXd covariance = MatrixXd::Identity(4, 4);
	const auto gauge = zhangComponentGaugeGls(
		measurement, covariance, design);
	BOOST_REQUIRE_MESSAGE(gauge.valid, gauge.failureReason);

	VectorXd selectedDirection(2); selectedDirection << 1, -1;
	const double selectedMean = selectedDirection.dot(gauge.mean);
	const double selectedVariance = (selectedDirection.transpose() *
		gauge.covariance * selectedDirection)(0, 0);
	auto frozenRoundPerr = [](double fractional, double variance)
	{
		double competingMass = 0;
		const double factor = -0.25 / variance;
		for (int integer = 1; integer < 10; integer++)
		{
			competingMass += std::exp(
				(integer + 2 * fractional) * integer * factor);
			competingMass += std::exp(
				(integer - 2 * fractional) * integer * factor);
		}
		return competingMass / (1 + competingMass);
	};
	const double selectedPerr = frozenRoundPerr(
		selectedMean - std::round(selectedMean), selectedVariance);
	BOOST_CHECK_GT(selectedPerr, 1e-3);

	const auto product = zhangComponentGaugeToProductRow(
		{1, -1}, {0, 2, 4}, {3, -2, 5}, 5, 0);
	BOOST_REQUIRE(product.valid);
	BOOST_CHECK(zhangExactPrimitiveRowLattice({product.row}, 5));
	const ZhangExactMatrix physicalBasis = zhangExactIdentityMatrix(5);
	const auto physical = zhangExactRowCombination(
		product.row, physicalBasis);
	BOOST_CHECK(physical == product.row);

	Eigen::SelfAdjointEigenSolver<MatrixXd> modes(gauge.covariance);
	BOOST_REQUIRE(modes.info() == Eigen::Success);
	double maximumAlignment = 0;
	for (int column = 0; column < modes.eigenvectors().cols(); column++)
	{
		maximumAlignment = std::max(maximumAlignment,
			std::abs(selectedDirection.dot(modes.eigenvectors().col(column))) /
			selectedDirection.norm());
	}
	BOOST_CHECK_GT(maximumAlignment, 0);

	const auto receiverLeaveOut = zhangComponentGaugeLeaveOneGroupOut(
		measurement, covariance, design, {0}, selectedDirection);
	const auto arcLeaveOut = zhangComponentGaugeLeaveOneGroupOut(
		measurement, covariance, design, {2}, selectedDirection);
	BOOST_REQUIRE(receiverLeaveOut.valid);
	BOOST_REQUIRE(arcLeaveOut.valid);

	const auto alternateProduct = zhangComponentGaugeToProductRow(
		{1, -1}, {0, 2, 4}, {3, -2, 5}, 5, 0);
	BOOST_REQUIRE(alternateProduct.valid);
	BOOST_CHECK(alternateProduct.row == product.row);
	BOOST_CHECK_EQUAL(alternateProduct.value, product.value);
}

BOOST_AUTO_TEST_CASE(component_gauge_residual_conflict_fails_closed)
{
	// Equal deterministic observations with unequal values contain a true null
	// space contradiction; no tolerance may convert it into an integer row.
	MatrixXd design(2, 1); design << 1, 1;
	VectorXd measurement(2); measurement << 2, 3;
	MatrixXd covariance(2, 2); covariance << 1, 1, 1, 1;
	const auto gls = zhangComponentGaugeGls(measurement, covariance, design);
	BOOST_CHECK(!gls.valid);
	BOOST_CHECK_EQUAL(
		gls.failureReason,
		"COMPONENT_GAUGE_COVARIANCE_NULLSPACE_CONFLICT");
	BOOST_CHECK(!gls.nullModes.empty());
	const std::vector<ZhangComponentEdgeId> edges = {{0, 1, 0, 1}, {0, 1, 0, 2}};
	const auto conflicts = zhangLocalizeNullConflicts(gls, edges,
		{E_ObsCode::L1C, E_ObsCode::L1C});
	BOOST_CHECK(!conflicts.empty());
}

BOOST_AUTO_TEST_CASE(
	component_gauge_uses_same_conditioned_snapshot_for_mean_and_covariance)
{
	VectorXd primitiveMean(3); primitiveMean << 1.2, 0.1, 3.0;
	MatrixXd primitiveCovariance = MatrixXd::Identity(3, 3);
	MatrixXd certifiedRow(1, 3); certifiedRow << 1, -1, 0;
	VectorXd certifiedInteger(1); certifiedInteger << 1;
	const auto snapshot = zhangConditionExactProductRows(
		primitiveMean, primitiveCovariance, certifiedRow, certifiedInteger);
	BOOST_REQUIRE(snapshot.valid);
	BOOST_CHECK_SMALL(
		(certifiedRow * snapshot.mean - certifiedInteger).norm(), 1e-12);
	BOOST_CHECK_SMALL((certifiedRow * snapshot.covariance).norm(), 1e-12);

	// Both observations estimate the same component gauge only after their
	// certified within-component affine offsets are propagated from that exact
	// snapshot with the same primitive row matrix B.
	MatrixXd observationRows(2, 3);
	observationRows << 1, 0, -1,
		0, 1, -1;
	VectorXd affineOffsets(2); affineOffsets << 1, 0;
	const VectorXd measurements =
		observationRows * snapshot.mean - affineOffsets;
	const MatrixXd covariance = observationRows * snapshot.covariance *
		observationRows.transpose();
	MatrixXd design = MatrixXd::Ones(2, 1);
	const auto gls = zhangComponentGaugeGls(measurements, covariance, design);
	BOOST_REQUIRE_MESSAGE(gls.valid, gls.failureReason);
	BOOST_CHECK_SMALL(gls.maximumNullResidual, 1e-12);
	BOOST_CHECK_SMALL(measurements(0) - measurements(1), 1e-12);
}

BOOST_AUTO_TEST_CASE(component_gauge_rank_deficiency_keeps_maximal_integer_forest)
{
	// Gauge 0 has both frequency coordinates; gauge 1 has only L1 support.
	// The primitive forest must retain gauge 0 instead of returning rank zero.
	MatrixXd information = MatrixXd::Zero(4, 4);
	information(0, 0) = 1;
	information(2, 2) = 1;
	information(1, 1) = 1;
	const auto forest = zhangMaxEstimableDualGaugeForest(information, 2);
	BOOST_REQUIRE_EQUAL(forest.size(), 1);
	BOOST_CHECK_EQUAL(forest.front(), 0);
}

BOOST_AUTO_TEST_CASE(component_gauge_block_forest_uses_estimable_incidence_edges)
{
	// Four components have two independent dual-frequency bridges: 0-1 and
	// 2-3.  The latter is not a datum-star coordinate, so a global rank test or
	// a star-only fallback would wrongly lose it.  The result must retain both
	// primitive e_a-e_b rows and leave the two blocks independent.
	constexpr int componentCount = 4;
	constexpr int gauges = componentCount - 1;
	MatrixXd information = MatrixXd::Zero(2 * gauges, 2 * gauges);
	auto addDualIncidence = [&information](const VectorXd& incidence)
	{
		information.topLeftCorner(3, 3) += incidence * incidence.transpose();
		information.bottomRightCorner(3, 3) += incidence * incidence.transpose();
	};
	VectorXd zeroOne = VectorXd::Zero(gauges); zeroOne(0) = 1;
	VectorXd twoThree = VectorXd::Zero(gauges); twoThree(1) = 1; twoThree(2) = -1;
	addDualIncidence(zeroOne);
	addDualIncidence(twoThree);
	const std::vector<ZhangComponentEdgeId> candidates = {
		{0, 1, 0, 1}, {0, 2, 0, 2}, {2, 3, 2, 3}};
	const auto forest = zhangMaximumEstimableDualComponentForest(
		information, componentCount, candidates);
	BOOST_REQUIRE_EQUAL(forest.size(), 2);
	std::set<std::pair<int, int>> pairs;
	for (const auto& edge : forest)
		pairs.insert(std::minmax(edge.firstComponent, edge.secondComponent));
	BOOST_CHECK(pairs.contains({0, 1}));
	BOOST_CHECK(pairs.contains({2, 3}));
	BOOST_CHECK(!pairs.contains({0, 2}));
	const auto blocks = zhangDualComponentSupportBlocks(componentCount, forest);
	BOOST_REQUIRE_EQUAL(blocks.size(), 2);
	BOOST_CHECK((blocks[0] == std::vector<int>{0, 1}));
	BOOST_CHECK((blocks[1] == std::vector<int>{2, 3}));
}

BOOST_AUTO_TEST_CASE(component_gauge_forest_uses_cross_frequency_wl_covariance)
{
	// The two frequencies have identical marginal variances.  Only Q12 makes
	// component 0-1 the precise WL bridge, so a marginal-variance score cannot
	// distinguish the candidates while the full WL score must retain it.
	MatrixXd covariance = MatrixXd::Identity(4, 4);
	covariance(0, 2) = covariance(2, 0) = 0.9;
	const MatrixXd information = covariance.inverse();
	const std::vector<ZhangComponentEdgeId> candidates = {
		{0, 1, 0, 1}, {0, 2, 0, 2}, {1, 2, 1, 2}};
	const auto forest = zhangMaximumEstimableDualComponentForest(
		information, 3, candidates);
	BOOST_REQUIRE_EQUAL(forest.size(), 2);
	BOOST_CHECK(std::any_of(forest.begin(), forest.end(), [](const auto& edge)
	{
		return std::min(edge.firstComponent, edge.secondComponent) == 0 &&
			std::max(edge.firstComponent, edge.secondComponent) == 1;
	}));
}

BOOST_AUTO_TEST_CASE(component_gauge_forest_prefers_receiver_and_arc_disjoint_bridge)
{
	MatrixXd covariance = MatrixXd::Identity(4, 4);
	covariance(0, 2) = covariance(2, 0) = 0.9;
	const MatrixXd information = covariance.inverse();
	ZhangComponentEdgeId first{0, 1, 0, 1};
	first.receiverSupport = {"R1"};
	first.arcSupport = {"R1|G01|L1|A"};
	ZhangComponentEdgeId repeated{0, 2, 0, 2};
	repeated.receiverSupport = {"R1"};
	repeated.arcSupport = {"R1|G02|L1|A"};
	ZhangComponentEdgeId independent = repeated;
	independent.firstNode = 7;
	independent.secondNode = 8;
	independent.receiverSupport = {"R2"};
	independent.arcSupport = {"R2|G02|L1|B"};
	const auto forest = zhangMaximumEstimableDualComponentForest(
		information, 3, {first, repeated, independent});
	BOOST_REQUIRE_EQUAL(forest.size(), 2);
	BOOST_CHECK(std::any_of(forest.begin(), forest.end(), [](const auto& edge)
	{
		return edge.firstNode == 7 && edge.secondNode == 8;
	}));
}

BOOST_AUTO_TEST_CASE(component_gauge_real_null_direction_never_becomes_integer_row)
{
	// Only e1-e2 is estimable.  The individual e1 coordinate has a real
	// projection onto that line but is not itself estimable and cannot become
	// a product integer constraint.
	MatrixXd information = MatrixXd::Zero(4, 4);
	VectorXd incidence(2); incidence << 1, -1;
	information.topLeftCorner(2, 2) = incidence * incidence.transpose();
	information.bottomRightCorner(2, 2) = incidence * incidence.transpose();
	VectorXd star = VectorXd::Zero(4); star(0) = 1;
	VectorXd edge = VectorXd::Zero(4); edge(0) = 1; edge(1) = -1;
	BOOST_CHECK(!zhangComponentGaugeFunctionalEstimable(information, star));
	BOOST_CHECK(zhangComponentGaugeFunctionalEstimable(information, edge));
}

BOOST_AUTO_TEST_CASE(
	component_gauge_primitive_forest_recovers_single_reliable_bridge_in_large_block)
{
	// Eight components form a seven-dimensional block.  Only datum--component
	// 1 has a sufficiently precise integer marginal.  The primitive catalogue
	// and rank-interleaved forest logic must retain that rank-one exit instead
	// of treating failure of the complete seven-edge tree as rank zero.
	constexpr int componentCount = 8;
	VectorXd mean = VectorXd::Zero(componentCount - 1);
	MatrixXd covariance = MatrixXd::Identity(componentCount - 1,
		componentCount - 1);
	covariance(0, 0) = 1e-4;
	const auto generated = zhangPrimitiveComponentGaugeCandidates(
		mean, covariance, componentCount);
	std::vector<ZhangNamedPairBeamCandidate> reliable;
	for (const auto& candidate : generated)
	{
		const double perr = zhangTestRoundPerr(
			candidate.fractional, candidate.variance);
		if (perr <= 1e-3)
			reliable.push_back({candidate.row, perr, 0, candidate.variance,
				{candidate.firstComponent, candidate.secondComponent}, 1, 1});
	}
	BOOST_REQUIRE_EQUAL(reliable.size(), 1);
	BOOST_CHECK((reliable.front().nodes == std::vector<int>{0, 1}));
	const auto levels = zhangNamedPairForestBeamLevels(
		reliable, componentCount - 1, 8);
	BOOST_REQUIRE_EQUAL(levels.size(), 1);
	BOOST_CHECK_EQUAL(levels.front().front().selected.size(), 1);
}

BOOST_AUTO_TEST_CASE(
	component_gauge_primitive_forest_recovers_two_independent_bridges)
{
	constexpr int componentCount = 6;
	VectorXd mean = VectorXd::Zero(componentCount - 1);
	MatrixXd covariance = MatrixXd::Identity(componentCount - 1,
		componentCount - 1);
	covariance(0, 0) = 1e-4;
	covariance(1, 1) = 1e-4;
	const auto generated = zhangPrimitiveComponentGaugeCandidates(
		mean, covariance, componentCount);
	std::vector<ZhangNamedPairBeamCandidate> reliable;
	for (const auto& candidate : generated)
	{
		const double perr = zhangTestRoundPerr(
			candidate.fractional, candidate.variance);
		if (perr <= 1e-3)
			reliable.push_back({candidate.row, perr, 0, candidate.variance,
				{candidate.firstComponent, candidate.secondComponent}, 1, 1});
	}
	const auto levels = zhangNamedPairForestBeamLevels(
		reliable, componentCount - 1, 8);
	BOOST_REQUIRE_GE(levels.size(), 2);
	BOOST_CHECK_EQUAL(levels[1].front().selected.size(), 2);
	BOOST_CHECK_EQUAL(levels.size(), 2);
}

BOOST_AUTO_TEST_CASE(
	component_gauge_dense_only_direction_remains_conditioning_only)
{
	// [1,1] is extremely precise, while every primitive component edge is
	// noisy.  A decorrelated dense direction may be retained for posterior
	// conditioning, but it must not create a named graph certificate.
	Vector2d mean = Vector2d::Zero();
	Matrix2d covariance;
	covariance << 1, -0.9999, -0.9999, 1;
	const Vector2d dense = (Vector2d() << 1, 1).finished();
	const double denseVariance =
		(dense.transpose() * covariance * dense)(0, 0);
	BOOST_CHECK_LT(zhangTestRoundPerr(0, denseVariance), 1e-3);
	const auto generated = zhangPrimitiveComponentGaugeCandidates(
		mean, covariance, 3);
	int reliablePrimitiveEdges = 0;
	for (const auto& candidate : generated)
		reliablePrimitiveEdges += zhangTestRoundPerr(
			candidate.fractional, candidate.variance) <= 1e-3;
	BOOST_CHECK_EQUAL(reliablePrimitiveEdges, 0);
}

BOOST_AUTO_TEST_CASE(
	component_gauge_l1_candidates_require_exact_wl_edge_membership)
{
	Vector2d mean = Vector2d::Zero();
	Matrix2d covariance = 1e-4 * Matrix2d::Identity();
	const ZhangExactMatrix acceptedWideLane{{1, 0}};
	const auto candidates = zhangPrimitiveComponentGaugeCandidates(
		mean, covariance, 3, acceptedWideLane);
	BOOST_REQUIRE_EQUAL(candidates.size(), 1);
	BOOST_CHECK((candidates.front().row == ZhangExactVector{-1, 0}));
}

BOOST_AUTO_TEST_CASE(
	component_gauge_square_bridge_reaches_provisional_ledger_while_unreliable)
{
	ZhangComponentGaugeBlockResult block;
	block.targetRank = 5;
	block.newDualRank = 1;
	block.valid = true;
	block.reliable = false;
	block.confirmationRequired = true;
	block.wlRows = {{1, 0}};
	block.wlIntegers = {3};
	block.l1Rows = {{1, 0}};
	block.l1Integers = {7};
	BOOST_CHECK(zhangComponentGaugeBlockProvidesProvisionalEvidence(block));
	block.valid = false;
	BOOST_CHECK(!zhangComponentGaugeBlockProvidesProvisionalEvidence(block));
	block.valid = true;
	block.l1Integers.clear();
	BOOST_CHECK(!zhangComponentGaugeBlockProvidesProvisionalEvidence(block));
}

BOOST_AUTO_TEST_CASE(
	component_gauge_block_family_budget_separates_authoritative_and_provisional)
{
	auto block = [](int rank, double failureProbability, bool authoritative)
	{
		ZhangComponentGaugeBlockResult result;
		result.targetRank = rank;
		result.newDualRank = rank;
		result.valid = true;
		result.reliable = authoritative;
		result.confirmationRequired = !authoritative;
		result.wlFailureProbability = failureProbability / 2;
		result.l1FailureProbability = failureProbability / 2;
		for (int row = 0; row < rank; row++)
		{
			ZhangExactVector exact(2);
			exact[row % 2] = 1;
			result.wlRows.push_back(exact);
			result.wlIntegers.push_back(row + 1);
			result.l1Rows.push_back(exact);
			result.l1Integers.push_back(row + 3);
		}
		return result;
	};
	std::vector<ZhangComponentGaugeBlockResult> blocks;
	blocks.push_back(block(1, 4e-4, true));
	blocks.push_back(block(1, 4e-4, false));
	blocks.push_back(block(2, 9e-4, false));
	auto invalid = block(2, 1e-5, true);
	invalid.valid = false;
	blocks.push_back(std::move(invalid));

	BOOST_CHECK(zhangComponentGaugeBlockProvidesAuthoritativeEvidence(blocks[0]));
	BOOST_CHECK(zhangComponentGaugeBlockProvidesProvisionalEvidence(blocks[1]));
	BOOST_CHECK(!zhangComponentGaugeBlockProvidesAuthoritativeEvidence(blocks[1]));
	BOOST_CHECK(!zhangComponentGaugeBlockProvidesProvisionalEvidence(blocks[3]));
	const auto selected = zhangSelectComponentGaugeEvidenceBlocks(blocks, 1e-3);
	BOOST_CHECK_EQUAL(selected.evidenceBlocks, 3);
	BOOST_CHECK_EQUAL(selected.budgetEligibleBlocks, 3);
	BOOST_CHECK((selected.selectedBlockIndices == std::vector<int>{0, 1}));
	BOOST_CHECK_EQUAL(selected.selectedTotalRank, 2);
	BOOST_CHECK_EQUAL(selected.selectedAuthoritativeRank, 1);
	BOOST_CHECK_EQUAL(selected.selectedAuthoritativeBlocks, 1);
	BOOST_CHECK_EQUAL(selected.selectedProvisionalBlocks, 1);
	BOOST_CHECK_CLOSE(selected.familyFailureProbability, 8e-4, 1e-9);

	const auto tighter = zhangSelectComponentGaugeEvidenceBlocks(blocks, 7e-4);
	BOOST_CHECK((tighter.selectedBlockIndices == std::vector<int>{0}));
	BOOST_CHECK_EQUAL(tighter.selectedTotalRank, 1);
	BOOST_CHECK_EQUAL(tighter.selectedAuthoritativeRank, 1);
}

BOOST_AUTO_TEST_CASE(product_closure_final_union_equals_single_conditioning)
{
	VectorXd mean(2); mean << 2.3, -1.7;
	MatrixXd covariance(2, 2); covariance << 2, 0.3, 0.3, 1;
	MatrixXd first(1, 2); first << 1, 0;
	MatrixXd second(1, 2); second << 0, 1;
	VectorXd firstInteger(1); firstInteger << 2;
	VectorXd secondInteger(1); secondInteger << -2;
	const auto once = zhangConditionExactProductRows(mean, covariance,
		(MatrixXd(2, 2) << 1, 0, 0, 1).finished(),
		(VectorXd(2) << 2, -2).finished());
	const auto stage1 = zhangConditionExactProductRows(mean, covariance, first, firstInteger);
	BOOST_REQUIRE(once.valid);
	BOOST_REQUIRE(stage1.valid);
	const auto stage2 = zhangConditionExactProductRows(stage1.mean, stage1.covariance,
		second, secondInteger);
	BOOST_REQUIRE(stage2.valid);
	BOOST_CHECK_SMALL((once.mean - stage2.mean).norm(), 1e-12);
	BOOST_CHECK_SMALL((once.covariance - stage2.covariance).norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(product_closure_is_not_gated_on_fresh_wide_lane)
{
	BOOST_CHECK(!zhangShouldRunProductClosure(false, false, false, false, false));
	BOOST_CHECK( zhangShouldRunProductClosure(true,  false, false, false, false));
	BOOST_CHECK( zhangShouldRunProductClosure(false, true,  false, false, false));
	BOOST_CHECK( zhangShouldRunProductClosure(false, false, true,  false, false));
	BOOST_CHECK( zhangShouldRunProductClosure(false, false, false, true,  false));
	BOOST_CHECK( zhangShouldRunProductClosure(false, false, false, false, true));
	BOOST_CHECK( zhangShouldRunProductClosure(
		false, false, false, false, false, true));
}

BOOST_AUTO_TEST_CASE(component_null_conflict_quarantine_uses_evidence_order)
{
	BOOST_CHECK_GT(zhangComponentEvidenceTrust("PRODUCT_GAUGE_CERTIFICATE"),
		zhangComponentEvidenceTrust("CURRENT_GENERATION_EXACT_PAIR"));
	BOOST_CHECK_GT(zhangComponentEvidenceTrust("CURRENT_GENERATION_EXACT_PAIR"),
		zhangComponentEvidenceTrust("CROSS_GENERATION_PROJECTED_LEDGER"));
	BOOST_CHECK_GT(zhangComponentEvidenceTrust("CROSS_GENERATION_PROJECTED_LEDGER"),
		zhangComponentEvidenceTrust("TEMPORAL_RECERTIFIED"));
	BOOST_CHECK_GT(zhangComponentEvidenceTrust("TEMPORAL_RECERTIFIED"),
		zhangComponentEvidenceTrust("BESD_CANDIDATE"));
}

BOOST_AUTO_TEST_CASE(product_closure_runs_without_new_network_wl)
{
	BOOST_CHECK(zhangShouldRunProductClosure(false, true, false, false, false));
}

BOOST_AUTO_TEST_CASE(product_gauge_ledger_alone_requires_stable_runtime_id)
{
	BOOST_CHECK(!zhangPersistentProductEvidenceNeedsRuntime(false, false));
	BOOST_CHECK( zhangPersistentProductEvidenceNeedsRuntime(true,  false));
	BOOST_CHECK( zhangPersistentProductEvidenceNeedsRuntime(false, true));
	BOOST_CHECK( zhangPersistentProductEvidenceNeedsRuntime(true,  true));
}

BOOST_AUTO_TEST_CASE(ledger_selected_subset_is_exactly_independent)
{
	const auto hnf = zhangExactRowHermiteNormalForm(
		{{1, 0, 0}, {2, 0, 0}, {0, 1, 0}}, {4, 8, -3});
	BOOST_REQUIRE(hnf.consistent);
	BOOST_CHECK_EQUAL(hnf.basis.size(), 2);
}

BOOST_AUTO_TEST_CASE(ledger_hnf_reduction_preserves_integer_values)
{
	const auto before = zhangExactRowHermiteNormalForm({{1, 1}, {1, -1}}, {5, 1});
	const auto after  = zhangExactRowHermiteNormalForm(before.basis, before.values);
	BOOST_REQUIRE(before.consistent && after.consistent);
	BOOST_CHECK(before.basis == after.basis);
	BOOST_CHECK(before.values == after.values);
}

BOOST_AUTO_TEST_CASE(ledger_exact_rank_equals_numeric_rank)
{
	const auto hnf = zhangExactRowHermiteNormalForm({{1, 0, 1}, {0, 1, -1}}, {2, -4});
	BOOST_REQUIRE(hnf.consistent);
	MatrixXd rows(2, 3); rows << 1,0,1, 0,1,-1;
	BOOST_CHECK_EQUAL(rows.fullPivHouseholderQr().rank(), hnf.basis.size());
}

BOOST_AUTO_TEST_CASE(partial_component_block_survives_other_block_failure)
{
	MatrixXd information = MatrixXd::Zero(6, 6);
	information(0,0) = information(3,3) = 1; // first dual block estimable
	const auto forest = zhangMaxEstimableDualGaugeForest(information, 3);
	BOOST_REQUIRE_EQUAL(forest.size(), 1);
	BOOST_CHECK_EQUAL(forest.front(), 0);
}

BOOST_AUTO_TEST_CASE(null_conflict_localizes_bad_relation_group)
{
	MatrixXd design(2, 1); design << 1, 1;
	VectorXd measurement(2); measurement << 2, 3;
	MatrixXd covariance(2, 2); covariance << 1, 1, 1, 1;
	const auto gls = zhangComponentGaugeGls(measurement, covariance, design);
	const auto conflicts = zhangLocalizeNullConflicts(gls,
		{{0,1,0,1}, {0,1,0,2}}, {E_ObsCode::L1C, E_ObsCode::L1C},
		{101, 202}, {"SEG-A", "SEG-B"},
		{"CURRENT_GENERATION_EXACT_PAIR", "BESD_PROVISIONAL"}, {4, -9});
	BOOST_REQUIRE(!conflicts.empty());
	BOOST_CHECK(!conflicts.front().dominantEdges.empty());
	BOOST_CHECK(std::isfinite(conflicts.front().nullResidual));
	BOOST_REQUIRE_EQUAL(conflicts.front().contributions.size(),
		conflicts.front().dominantEdges.size());
	BOOST_REQUIRE_EQUAL(conflicts.front().affineOffsets.size(),
		conflicts.front().dominantEdges.size());
	BOOST_CHECK_EQUAL(conflicts.front().backendGenerations.front(), 101U);
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_fails_closed_without_joint_evidence)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	certificate.jointNisPassed = false;
	ProductGaugeCertificateLedger ledger;
	const auto update = ledger.observe(100, {certificate}, 1);
	BOOST_REQUIRE(update.valid);
	BOOST_CHECK_EQUAL(update.rejectedCandidates, 1);
	BOOST_CHECK_EQUAL(ledger.activeRank(), 0);
	BOOST_CHECK(ledger.certificates().empty());
}

BOOST_AUTO_TEST_CASE(product_gauge_square_bridge_requires_multi_epoch_confirmation)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate, true);
	ProductGaugeCertificateLedger ledger;
	const auto rejected = ledger.observe(100, {certificate}, 1);
	BOOST_REQUIRE(rejected.valid);
	BOOST_CHECK_EQUAL(rejected.rejectedCandidates, 1);
	BOOST_CHECK(ledger.certificates().empty());
	BOOST_REQUIRE(ledger.observe(130, {certificate}, 2).valid);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	BOOST_CHECK_EQUAL(ledger.certificates().front().state,
		ProductGaugeCertificateState::CONFIRMING);
	BOOST_REQUIRE(ledger.observe(160, {certificate}, 2).valid);
	BOOST_CHECK_EQUAL(ledger.certificates().front().state,
		ProductGaugeCertificateState::ACTIVE);
}

BOOST_AUTO_TEST_CASE(conditioning_only_row_never_becomes_pair_certificate)
{
	ZhangProductIntegerConstraintSet constraints;
	constraints.conditioningOnlyRows = {{1, 1, -2}};
	constraints.conditioningOnlyIntegers = {0};
	BOOST_CHECK(constraints.dualFrequencyCertifiedPairs.empty());
	BOOST_CHECK_EQUAL(constraints.certifiedPairRank, 0);
}

BOOST_AUTO_TEST_CASE(dense_product_gauge_certificate_never_masquerades_as_writer_pair)
{
	ProductGaugeCertificate dense;
	dense.system = E_Sys::GPS;
	dense.firstObservable = E_ObsCode::L1C;
	dense.secondObservable = E_ObsCode::L2W;
	dense.satellite = SatSys(E_Sys::GPS, 2);
	dense.reference = SatSys(E_Sys::GPS, 3);
	dense.productRelationQ[SatSys(E_Sys::GPS, 2)] = 1;
	dense.productRelationQ[SatSys(E_Sys::GPS, 3)] = -2;
	dense.productRelationQ[SatSys(E_Sys::GPS, 4)] = 1;
	BOOST_CHECK(!zhangProductGaugeCertificateIsPrimitivePair(dense));

	ProductGaugeCertificate pair = dense;
	pair.productRelationQ.clear();
	pair.productRelationQ[pair.satellite] = 1;
	pair.productRelationQ[pair.reference] = -1;
	BOOST_CHECK(zhangProductGaugeCertificateIsPrimitivePair(pair));
}

BOOST_AUTO_TEST_CASE(strict_dual_graph_drives_pppar_usable)
{
	ZhangProductIntegerConstraintSet constraints;
	constraints.reliable = true;
	constraints.conditioningRank = 4;
	constraints.certifiedPairRank = 0;
	BOOST_CHECK(!(constraints.reliable && constraints.certifiedPairRank > 0));
	constraints.certifiedPairRank = 1;
	BOOST_CHECK(constraints.reliable && constraints.certifiedPairRank > 0);
}

BOOST_AUTO_TEST_CASE(final_product_constraints_are_conditioned_once)
{
	VectorXd mean(2); mean << 1.3, -2.2;
	MatrixXd covariance(2,2); covariance << 2, .2, .2, 1;
	MatrixXd rows(2,2); rows << 1,0, 0,1;
	VectorXd integers(2); integers << 1, -2;
	const auto fixed = zhangConditionExactProductRows(mean, covariance, rows, integers);
	BOOST_REQUIRE(fixed.valid);
	BOOST_CHECK_SMALL((rows * fixed.mean - integers).lpNorm<Eigen::Infinity>(), 1e-9);
	BOOST_CHECK_SMALL((rows * fixed.covariance * rows.transpose()).lpNorm<Eigen::Infinity>(), 1e-12);
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_survives_backend_generation)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	certificate.componentVersion = 99; // backend generation/version is not identity
	BOOST_REQUIRE(ledger.observe(130, {certificate}, 1).valid);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	BOOST_CHECK(ledger.certificates().front().active);
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_requires_exact_current_segments)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	const auto& active = ledger.certificates().front();
	BOOST_REQUIRE(active.state == ProductGaugeCertificateState::ACTIVE);
	BOOST_REQUIRE(active.active);
	BOOST_REQUIRE(active.currentAlignmentValid);
	BOOST_REQUIRE(active.firstObservable == E_ObsCode::L1C);
	BOOST_REQUIRE(active.secondObservable == E_ObsCode::L2W);
	BOOST_CHECK(zhangProductGaugeCertificateMatchesCurrentSegments(active,
		E_ObsCode::L1C, E_ObsCode::L2W,
		active.satellitePhaseSegments, active.referencePhaseSegments));
	BOOST_CHECK(!zhangProductGaugeCertificateMatchesCurrentSegments(active,
		E_ObsCode::L1C, E_ObsCode::L2W,
		{active.satellitePhaseSegments[0], "CHANGED_SEGMENT"},
		active.referencePhaseSegments));
}

BOOST_AUTO_TEST_CASE(product_gauge_certificate_dies_on_phase_segment_change)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	certificate.satellitePhaseSegments[1] = "G03-L2-SEG1";
	const auto update = ledger.observe(130, {certificate}, 1);
	BOOST_REQUIRE(update.valid);
	BOOST_CHECK_EQUAL(update.retiredSegmentCertificates, 1);
	BOOST_CHECK(!ledger.certificates().front().active);
	BOOST_CHECK(ledger.certificates().front().inactiveReason ==
		ProductGaugeCertificateInactiveReason::PHASE_SEGMENT_RETIRED);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_quiet_epoch_retains_and_schedules_historical_certificate)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.wideLaneInteger = 7;
	certificate.firstSignalInteger = -4;
	certificate.productRelationQ[certificate.satellite] = 1;
	certificate.productRelationQ[certificate.reference] = -1;
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	certificate.componentUuid = "GPS-COMPONENT-A";
	certificate.sourceCertificateIds = {"CERT-001"};
	certificate.failureProbabilityBudget = 1e-6;
	makeProductGaugeEvidenceComplete(certificate);

	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	// No new WL/current certificate at the next epoch must not retire history.
	BOOST_REQUIRE(ledger.observe(130, {}, 1).valid);
	BOOST_REQUIRE_EQUAL(ledger.rawActiveCount(), 1);
	BOOST_REQUIRE_EQUAL(ledger.activeRank(), 1);
	const auto& active = ledger.certificates().front();
	BOOST_CHECK_EQUAL(active.validFrom, 100);
	BOOST_CHECK_EQUAL(active.validTo, 0);
	BOOST_CHECK(active.componentUuid.rfind(
		"SATELLITE_LATTICE_COMPONENT:", 0) == 0);
	BOOST_REQUIRE_EQUAL(active.sourceCertificateIds.size(), 1);
	BOOST_CHECK_EQUAL(active.sourceCertificateIds.front(), "CERT-001");
	BOOST_CHECK_CLOSE(active.failureProbabilityBudget, 1e-6, 1e-9);
	BOOST_CHECK(zhangShouldRunProductClosure(
		false, false, false, ledger.activeRank() > 0, false));
	BOOST_CHECK(zhangProductGaugeCertificateMatchesCurrentSegments(
		active, E_ObsCode::L1C, E_ObsCode::L2W,
		active.satellitePhaseSegments, active.referencePhaseSegments));
	// Canonicalisation reverses G03-G02 and must reverse q and integers too.
	BOOST_CHECK(active.satellite == SatSys(E_Sys::GPS, 2));
	BOOST_CHECK(active.reference == SatSys(E_Sys::GPS, 3));
	BOOST_CHECK_EQUAL(active.wideLaneInteger, -7);
	BOOST_CHECK_EQUAL(active.firstSignalInteger, 4);
	BOOST_CHECK_EQUAL(active.productRelationQ.at(SatSys(E_Sys::GPS, 2)), 1);
	BOOST_CHECK_EQUAL(active.productRelationQ.at(SatSys(E_Sys::GPS, 3)), -1);
}

BOOST_AUTO_TEST_CASE(product_gauge_conflict_has_explicit_inactive_reason)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 2);
	certificate.reference = SatSys(E_Sys::GPS, 3);
	certificate.satellitePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	certificate.referencePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	certificate.firstSignalInteger = 1;
	const auto update = ledger.observe(130, {certificate}, 1);
	BOOST_REQUIRE(update.valid);
	BOOST_CHECK_EQUAL(update.conflicts, 1);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 2);
	BOOST_CHECK(ledger.certificates().front().inactiveReason ==
		ProductGaugeCertificateInactiveReason::CONFLICT);
	BOOST_CHECK(ledger.certificates().back().active);
}

BOOST_AUTO_TEST_CASE(product_components_require_exact_pair_lattice_membership)
{
	const SatSys g01(E_Sys::GPS, 1);
	const SatSys g02(E_Sys::GPS, 2);
	const SatSys g03(E_Sys::GPS, 3);
	const SatSys g04(E_Sys::GPS, 4);
	ProductGaugeCertificate relation;
	relation.system = E_Sys::GPS;
	relation.firstObservable = E_ObsCode::L1C;
	relation.secondObservable = E_ObsCode::L2W;
	relation.satellite = g01;
	relation.reference = g04;
	relation.productRelationQ = {{g01, 1}, {g02, 1}, {g03, -1}, {g04, -1}};
	relation.satellitePhaseSegments = {"G01-L1-SEG0", "G01-L2-SEG0"};
	relation.referencePhaseSegments = {"G04-L1-SEG0", "G04-L2-SEG0"};
	relation.phaseSegmentsBySatellite = {
		{g01, relation.satellitePhaseSegments},
		{g02, {"G02-L1-SEG0", "G02-L2-SEG0"}},
		{g03, {"G03-L1-SEG0", "G03-L2-SEG0"}},
		{g04, relation.referencePhaseSegments}};
	makeProductGaugeEvidenceComplete(relation);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {relation}, 1).valid);
	BOOST_CHECK_EQUAL(ledger.activeRank(), 1);
	// q=e1+e2-e3-e4 has all four satellites in its support, but it does not
	// imply any primitive pair difference.  Co-occurrence must not merge them.
	BOOST_CHECK(ledger.componentUuidFor(g01) != ledger.componentUuidFor(g02));
	BOOST_CHECK(ledger.componentUuidFor(g01) != ledger.componentUuidFor(g04));
	BOOST_CHECK(ledger.certificates().front().componentUuid.find('+') !=
		std::string::npos);
}

BOOST_AUTO_TEST_CASE(
	temporal_component_integer_basis_removes_common_gauge_and_is_reference_invariant)
{
	const auto g02 = zhangComponentRelativeGaugeBasis(5, 0);
	const auto g03 = zhangComponentRelativeGaugeBasis(5, 1);
	BOOST_REQUIRE_EQUAL(g02.size(), 4);
	BOOST_REQUIRE_EQUAL(g03.size(), 4);
	const ZhangExactVector common(5, 1);
	BOOST_CHECK(zhangExactMatrixTimesColumn(g02, common) ==
		ZhangExactVector(4));
	BOOST_CHECK(zhangExactMatrixTimesColumn(g03, common) ==
		ZhangExactVector(4));
	const auto hnf02 = zhangExactRowHermiteNormalForm(g02);
	const auto hnf03 = zhangExactRowHermiteNormalForm(g03);
	BOOST_REQUIRE(hnf02.consistent);
	BOOST_REQUIRE(hnf03.consistent);
	BOOST_CHECK(hnf02.basis == hnf03.basis);
	const auto primitive = zhangIntegerRowLatticeContains(
		g02, ZhangExactVector(5));
	BOOST_REQUIRE_EQUAL(primitive.smithInvariants.size(), 4);
	for (const auto& invariant : primitive.smithInvariants)
	{
		BOOST_CHECK(zhangExactAbs(invariant) == 1);
	}
}

BOOST_AUTO_TEST_CASE(product_relation_named_backward_search_reaches_low_rank)
{
	const std::vector<int> full = {0, 1, 2, 3, 4, 5};
	const auto withSeed = zhangProductNamedBackwardChildren(
		full, {1, 4}, 5, 1);
	BOOST_REQUIRE_EQUAL(withSeed.size(), 1);
	const std::vector<int> expectedSeed = {1, 4};
	BOOST_CHECK_EQUAL_COLLECTIONS(
		withSeed.front().begin(), withSeed.front().end(),
		expectedSeed.begin(), expectedSeed.end());

	std::vector<int> path = full;
	int evaluations = 1;
	while (path.size() > 1)
	{
		const auto children = zhangProductNamedBackwardChildren(
			path, {}, static_cast<int>(path.size()) - 1, 1);
		BOOST_REQUIRE_EQUAL(children.size(), 1);
		path = children.front();
		evaluations++;
	}
	BOOST_CHECK_EQUAL(path.size(), 1);
	BOOST_CHECK_EQUAL(evaluations, 6);
}

BOOST_AUTO_TEST_CASE(product_named_pair_beam_expands_alternate_forests)
{
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{1, 0, 0}, 1e-5, 1.0, 0.2, {0, 3}},
		{{0, 1, 0}, 1e-5, 100.0, 0.2, {1, 3}},
		{{0, 0, 1}, 1e-6, 0.1, 0.1, {2, 3}},
		{{1, -1, 0}, 1e-5, 0.5, 0.2, {0, 1}}};
	const auto levels = zhangNamedPairForestBeamLevels(candidates, 3, 2);
	BOOST_REQUIRE_EQUAL(levels.size(), 3);
	BOOST_REQUIRE_EQUAL(levels[0].size(), 2);
	// Candidate 0 is outside the rank-1 beam.  At rank two, candidate 3 joins
	// the retained most-reliable branch and wins on covered satellites before
	// the much larger gain of a lower-coverage alternative.
	BOOST_CHECK(std::find(
		levels[0][0].selected.begin(), levels[0][0].selected.end(), 0) ==
		levels[0][0].selected.end());
	BOOST_CHECK(std::find(
		levels[0][1].selected.begin(), levels[0][1].selected.end(), 0) ==
		levels[0][1].selected.end());
	const std::vector<int> expectedBest = {2, 3};
	BOOST_CHECK_EQUAL_COLLECTIONS(
		levels[1][0].selected.begin(), levels[1][0].selected.end(),
		expectedBest.begin(), expectedBest.end());
	BOOST_CHECK_EQUAL(levels[1][0].coveredNodes, 4);
	BOOST_CHECK_CLOSE_FRACTION(levels[1][0].summedGain, 0.6, 1e-14);
}

BOOST_AUTO_TEST_CASE(product_named_pair_beam_uses_selected_set_conditional_gain)
{
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{1, 0}, 1e-5, 100.0, 0.2, {0, 2}},
		{{0, 1}, 1e-5,   1.0, 0.2, {1, 2}}};
	const auto levels = zhangNamedPairForestBeamLevels(
		candidates, 2, 1,
		[](const std::vector<int>& selected)
		{
			// Candidate 0 has the larger independent gain, but the product
			// objective says the selected set containing candidate 1 is better.
			return std::find(selected.begin(), selected.end(), 1) != selected.end()
				? 10.0 : 0.0;
		});
	BOOST_REQUIRE(!levels.empty());
	BOOST_REQUIRE(!levels.front().empty());
	BOOST_REQUIRE_EQUAL(levels.front().front().selected.size(), 1U);
	BOOST_CHECK_EQUAL(levels.front().front().selected.front(), 1);
}

BOOST_AUTO_TEST_CASE(product_named_pair_beam_prioritises_dual_graph_gain)
{
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{1, 0}, 1e-5, 100.0, 0.2, {0, 2}, 0},
		{{0, 1}, 1e-5,   1.0, 0.2, {1, 2}, 1}};
	const auto levels = zhangNamedPairForestBeamLevels(
		candidates, 2, 1);

	BOOST_REQUIRE(!levels.empty());
	BOOST_REQUIRE(!levels.front().empty());
	BOOST_REQUIRE_EQUAL(levels.front().front().selected.size(), 1U);
	BOOST_CHECK_EQUAL(levels.front().front().selected.front(), 1);
	BOOST_CHECK_EQUAL(levels.front().front().dualGraphRankGain, 1);
	BOOST_CHECK_CLOSE_FRACTION(
		levels.front().front().summedGain, 1.0, 1e-14);
}

BOOST_AUTO_TEST_CASE(
	product_named_pair_beam_treats_admitted_perr_as_tie_breaker)
{
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{1, 0}, 1e-8, 100.0, 0.1, {0, 2}, 0, 0},
		{{0, 1}, 9e-4,   1.0, 0.2, {1, 2}, 1, 1}};
	const auto levels = zhangNamedPairForestBeamLevels(
		candidates, 2, 1);

	BOOST_REQUIRE(!levels.empty());
	BOOST_REQUIRE(!levels.front().empty());
	BOOST_REQUIRE_EQUAL(levels.front().front().selected.size(), 1U);
	BOOST_CHECK_EQUAL(levels.front().front().selected.front(), 1);
	BOOST_CHECK_EQUAL(levels.front().front().dualGraphRankGain, 1);
	BOOST_CHECK_EQUAL(levels.front().front().productLatticeRankGain, 1);
}

BOOST_AUTO_TEST_CASE(
	product_named_pair_beam_counts_component_cut_only_once)
{
	// Nodes 0, 1 and the implicit reference 3 already share one certified
	// dual-frequency component.  Both candidate edges reach the isolated node
	// 2, so their standalone flags are one but their joint graph-rank gain is
	// only one.
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{0, 1, -1}, 1e-5, 2.0, 0.2, {1, 2}, 1, 1},
		{{1, 0, -1}, 1e-5, 1.0, 0.2, {0, 2}, 1, 1}};
	const ZhangExactMatrix existingSignalRows = {{1, -1, 0}};
	const std::vector<int> existingDualComponents = {0, 0, 2, 0};
	const auto levels = zhangNamedPairForestBeamLevels(
		candidates, 3, 2, {}, std::numeric_limits<int>::max(), nullptr,
		existingSignalRows, existingDualComponents);

	BOOST_REQUIRE_EQUAL(levels.size(), 2U);
	BOOST_REQUIRE_EQUAL(levels[1].size(), 1U);
	BOOST_CHECK_EQUAL(levels[1].front().selected.size(), 2U);
	BOOST_CHECK_EQUAL(levels[1].front().dualGraphRankGain, 1);
	BOOST_CHECK_EQUAL(levels[1].front().productLatticeRankGain, 1);
}

BOOST_AUTO_TEST_CASE(product_named_pair_beam_prioritises_exact_quotient_gain)
{
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{1, 0}, 1e-5, 100.0, 0.2, {0, 2}},
		{{0, 1}, 1e-5,   1.0, 0.2, {1, 2}}};
	const ZhangExactMatrix existingSignalLattice = {{1, 0}};
	const auto levels = zhangNamedPairForestBeamLevels(
		candidates, 2, 1, {}, std::numeric_limits<int>::max(), nullptr,
		existingSignalLattice);

	BOOST_REQUIRE(!levels.empty());
	BOOST_REQUIRE(!levels.front().empty());
	BOOST_REQUIRE_EQUAL(levels.front().front().selected.size(), 1U);
	BOOST_CHECK_EQUAL(levels.front().front().selected.front(), 1);
	BOOST_CHECK_EQUAL(levels.front().front().productLatticeRankGain, 1);
	BOOST_CHECK_CLOSE_FRACTION(
		levels.front().front().summedGain, 1.0, 1e-14);
}

BOOST_AUTO_TEST_CASE(product_named_pair_beam_honours_construction_expansion_cap)
{
	const std::vector<ZhangNamedPairBeamCandidate> candidates = {
		{{1, 0, 0}, 1e-5, 4.0, 0.2, {0, 3}},
		{{0, 1, 0}, 1e-5, 3.0, 0.2, {1, 3}},
		{{0, 0, 1}, 1e-5, 2.0, 0.2, {2, 3}},
		{{1, -1, 0}, 1e-5, 1.0, 0.2, {0, 1}}};
	bool capped = false;
	const auto levels = zhangNamedPairForestBeamLevels(
		candidates, 3, 2, {}, 1, &capped);
	BOOST_CHECK(capped);
	BOOST_REQUIRE(!levels.empty());
	BOOST_REQUIRE_EQUAL(levels.front().size(), 1U);
	BOOST_CHECK_EQUAL(levels.front().front().selected.front(), 0);
}

BOOST_AUTO_TEST_CASE(product_relation_constraints_pull_back_affine_wl_and_l1)
{
	ZhangProductRelationBasis first;
	first.mappableTargetRank = 2;
	first.transform.resize(2, 3);
	first.transform << 1, -1, 0,
		0, 1, -1;
	first.affineOffsets = {2, -1};
	ZhangProductRelationBasis second;
	second.mappableTargetRank = 2;
	second.transform.resize(2, 3);
	second.transform << 1, 0, -1,
		1, -1, 0;
	second.affineOffsets = {1, 3};

	const ZhangExactMatrix wideLaneRows = {{1, -1}};
	const ZhangExactVector wideLaneIntegers = {7};
	const ZhangExactMatrix firstRows = {{2, 1}};
	const ZhangExactVector firstIntegers = {-4};
	ZhangExactMatrix networkRows;
	ZhangExactVector networkIntegers;
	std::string failure;
	BOOST_REQUIRE(zhangPullBackProductIntegerConstraints(
		first, second,
		wideLaneRows, wideLaneIntegers,
		firstRows, firstIntegers,
		networkRows, networkIntegers, failure));
	BOOST_CHECK_EQUAL(failure, "NONE");
	BOOST_REQUIRE_EQUAL(networkRows.size(), 2);
	BOOST_REQUIRE_EQUAL(networkIntegers.size(), 2);
	BOOST_CHECK(networkRows[0] == ZhangExactVector({1, -3, 2}));
	BOOST_CHECK_EQUAL(networkIntegers[0], 2);
	BOOST_CHECK(networkRows[1] == ZhangExactVector({2, -1, -1}));
	BOOST_CHECK_EQUAL(networkIntegers[1], -7);

	// M0: prove the affine pullback, rather than merely its hand-computed
	// coefficients.  For this integer ambiguity vector the product-space WL
	// and L1 equations and the two returned network equations have identical
	// residuals, exactly in Z (no floating-point comparison involved).
	const ZhangExactVector ambiguity = {3, 1, -2};
	auto exactDot = [](const ZhangExactVector& row, const ZhangExactVector& value)
	{
		ZhangExactInteger result = 0;
		for (size_t index = 0; index < row.size(); index++)
		{
			result += row[index] * value[index];
		}
		return result;
	};
	const ZhangExactInteger z1_0 = exactDot({1, -1, 0}, ambiguity) + 2;
	const ZhangExactInteger z1_1 = exactDot({0, 1, -1}, ambiguity) - 1;
	const ZhangExactInteger z2_0 = exactDot({1, 0, -1}, ambiguity) + 1;
	const ZhangExactInteger z2_1 = exactDot({1, -1, 0}, ambiguity) + 3;
	BOOST_CHECK_EQUAL(exactDot(networkRows[0], ambiguity) - networkIntegers[0],
		((z1_0 - z1_1) - (z2_0 - z2_1)) - wideLaneIntegers[0]);
	BOOST_CHECK_EQUAL(exactDot(networkRows[1], ambiguity) - networkIntegers[1],
		2 * z1_0 + z1_1 - firstIntegers[0]);
}

BOOST_AUTO_TEST_CASE(product_relation_constraint_pullback_rejects_noninteger_basis)
{
	ZhangProductRelationBasis first;
	first.mappableTargetRank = 1;
	first.transform = MatrixXd::Constant(1, 1, 0.5);
	first.affineOffsets = {0};
	ZhangProductRelationBasis second = first;
	second.transform(0, 0) = 0;
	ZhangExactMatrix networkRows;
	ZhangExactVector networkIntegers;
	std::string failure;
	BOOST_CHECK(!zhangPullBackProductIntegerConstraints(
		first, second, {{1}}, {0}, {}, {},
		networkRows, networkIntegers, failure));
	BOOST_CHECK_EQUAL(failure, "WL_PRODUCT_TO_NETWORK_MAPPING_FAILED");
	BOOST_CHECK(networkRows.empty());
}

BOOST_AUTO_TEST_CASE(product_integer_ledger_accumulates_exact_physical_rank)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow first;
	first.system = E_Sys::GPS;
	first.firstObservable = E_ObsCode::L1C;
	first.secondObservable = E_ObsCode::L2W;
	first.productRow = {1, -1, 0, 0};
	first.integerValue = 7;
	first.physicalExpansion = {{"L1C|ABCD|G03|V4", 1},
		{"L1C|ABCD|G02|V2", -1}};
	first.canonicalProductExpansion = {
		{"L1C|G03", 1}, {"L1C|G02", -1}};
	first.phaseSegmentFingerprint = "G02|L1C|SEG1;G03|L1C|SEG1;";
	first.backendBasisGeneration = 12;

	auto firstEpoch = ledger.observe(100, {first}, 2);
	BOOST_REQUIRE(firstEpoch.valid);
	BOOST_CHECK_EQUAL(firstEpoch.activeRankAfter, 0);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(!ledger.rows().front().certified);

	// Re-observing in the same epoch cannot manufacture a confirmation.
	auto duplicateEpoch = ledger.observe(100, {first}, 2);
	BOOST_REQUIRE(duplicateEpoch.valid);
	BOOST_CHECK_EQUAL(duplicateEpoch.activeRankAfter, 0);
	BOOST_CHECK_EQUAL(ledger.rows().front().confirmationEpochs, 1);

	auto confirmed = ledger.observe(130, {first}, 2);
	BOOST_REQUIRE(confirmed.valid);
	BOOST_CHECK_EQUAL(confirmed.activeRankAfter, 1);
	BOOST_CHECK(ledger.rows().front().certified);

	auto second = first;
	second.productRow = {0, 1, -1, 0};
	second.integerValue = -3;
	second.physicalExpansion = {{"L1C|EFGH|G08|V1", 1},
		{"L1C|EFGH|G02|V5", -1}};
	second.canonicalProductExpansion = {
		{"L1C|G08", 1}, {"L1C|G02", -1}};
	second.phaseSegmentFingerprint = "G02|L1C|SEG1;G08|L1C|SEG1;";
	auto secondFirst = ledger.observe(160, {second}, 2);
	BOOST_REQUIRE(secondFirst.valid);
	BOOST_CHECK_EQUAL(secondFirst.activeRankAfter, 1);
	auto secondConfirmed = ledger.observe(190, {second}, 2);
	BOOST_REQUIRE(secondConfirmed.valid);
	BOOST_CHECK_EQUAL(secondConfirmed.activeRankAfter, 2);
}

BOOST_AUTO_TEST_CASE(product_integer_ledger_rejects_nonprimitive_transactionally)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow good;
	good.system = E_Sys::GPS;
	good.productRow = {1};
	good.integerValue = 4;
	good.physicalExpansion = {{"L1C|ABCD|G03|V1", 1}};
	good.canonicalProductExpansion = {{"L1C|G03", 1}};
	good.phaseSegmentFingerprint = "G03|L1C|SEG1;";
	BOOST_REQUIRE(ledger.observe(100, {good}, 1).valid);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);

	auto bad = good;
	bad.integerValue = 8;
	bad.physicalExpansion = {{"L1C|EFGH|G05|V1", 2}};
	bad.canonicalProductExpansion = {{"L1C|G05", 2}};
	const auto rejected = ledger.observe(130, {good, bad}, 1);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.failureReason,
		"PRODUCT_LEDGER_ROW_NOT_PRIMITIVE");
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK_EQUAL(ledger.rows().front().integerValue, 4);
}

BOOST_AUTO_TEST_CASE(product_integer_ledger_mature_conflict_is_quarantined)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.firstObservable = E_ObsCode::L1C;
	row.secondObservable = E_ObsCode::L2W;
	row.productRow = {1, -1};
	row.integerValue = 7;
	row.physicalExpansion = {
		{"L1C|ABCD|G03|V1", 1}, {"L1C|ABCD|G02|V1", -1}};
	row.canonicalProductExpansion = {
		{"L1C|G03", 1}, {"L1C|G02", -1}};
	row.phaseSegmentFingerprint = "G02|L1C|SEG1;G03|L1C|SEG1;";
	row.backendBasisGeneration = 12;
	BOOST_REQUIRE(ledger.observe(100, {row}, 2).valid);
	const auto accepted = ledger.observe(130, {row}, 2);
	BOOST_REQUIRE(accepted.valid);
	BOOST_CHECK(zhangProductLedgerWriterCommitAuthorized(accepted));
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(ledger.rows().front().certified);

	auto fresh = row;
	fresh.productRow = {2, -1};
	fresh.integerValue = 3;
	fresh.physicalExpansion = {
		{"L1C|ABCD|G04|V1", 1}, {"L1C|ABCD|G02|V1", -1}};
	fresh.canonicalProductExpansion = {
		{"L1C|G04", 1}, {"L1C|G02", -1}};
	fresh.phaseSegmentFingerprint = "G02|L1C|SEG1;G04|L1C|SEG1;";
	auto replacement = row;
	replacement.integerValue = 8;
	const auto conflict = ledger.observe(160, {fresh, replacement}, 2);
	BOOST_CHECK(!conflict.valid);
	BOOST_CHECK_EQUAL(conflict.conflictingRows, 1);
	BOOST_CHECK_EQUAL(conflict.failureReason,
		"PRODUCT_LEDGER_MATURE_INTEGER_CONFLICT");
	BOOST_CHECK_EQUAL(conflict.activeRankAfter, conflict.activeRankBefore);
	BOOST_CHECK(!zhangProductLedgerWriterCommitAuthorized(conflict));
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK_EQUAL(ledger.rows().front().integerValue, -7);
	BOOST_CHECK_EQUAL(ledger.rows().front().confirmationEpochs, 2);
	BOOST_CHECK(ledger.rows().front().certified);

	const auto stillQuarantined = ledger.observe(190, {replacement}, 2);
	BOOST_CHECK(!stillQuarantined.valid);
	BOOST_CHECK_EQUAL(stillQuarantined.conflictingRows, 1);
	BOOST_CHECK_EQUAL(stillQuarantined.failureReason,
		"PRODUCT_LEDGER_MATURE_INTEGER_CONFLICT");
	BOOST_CHECK(!zhangProductLedgerWriterCommitAuthorized(stillQuarantined));
	BOOST_CHECK(ledger.rows().front().certified);
	BOOST_CHECK_EQUAL(ledger.rows().front().confirmationEpochs, 2);
	BOOST_CHECK_EQUAL(ledger.rows().front().integerValue, -7);
}

BOOST_AUTO_TEST_CASE(product_integer_ledger_canonicalises_negated_relation)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.firstObservable = E_ObsCode::L1C;
	row.secondObservable = E_ObsCode::L2W;
	row.productRow = {1, -1};
	row.integerValue = 7;
	row.physicalExpansion = {
		{"L1C|ABCD|G02|V1", -1}, {"L1C|ABCD|G03|V1", 1}};
	row.canonicalProductExpansion = {
		{"L1C|G02", -1}, {"L1C|G03", 1}};
	row.phaseSegmentFingerprint = "segments-v1";
	row.backendBasisGeneration = 12;
	row.pairCertificate = true;
	row.conditioningOnly = false;
	row.coordinate = "WL";
	row.firstSatellite = "G03";
	row.secondSatellite = "G02";
	BOOST_REQUIRE(ledger.observe(100, {row}, 2).valid);

	auto negated = row;
	for (auto& coefficient : negated.productRow) coefficient = -coefficient;
	for (auto& [identity, coefficient] : negated.physicalExpansion)
		coefficient = -coefficient;
	for (auto& [identity, coefficient] : negated.canonicalProductExpansion)
		coefficient = -coefficient;
	negated.integerValue = -negated.integerValue;
	std::swap(negated.firstSatellite, negated.secondSatellite);
	const auto update = ledger.observe(130, {negated}, 2);
	BOOST_REQUIRE(update.valid);
	BOOST_CHECK_EQUAL(update.freshRows, 0);
	BOOST_CHECK_EQUAL(update.conflictingRows, 0);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(ledger.rows().front().certified);
	BOOST_CHECK_EQUAL(ledger.rows().front().confirmationEpochs, 2);
	BOOST_CHECK(
		ledger.rows().front().physicalExpansion.begin()->second > 0);
}

BOOST_AUTO_TEST_CASE(product_ledger_physical_fingerprint_is_row_local)
{
	const std::map<std::string, ZhangExactInteger> row = {
		{"L2W|EFGH|G02|V4", -1},
		{"L1C|ABCD|G03|V1", 1},
		{"L1C|IJKL|G03|V9", -2}};
	const auto fingerprint = zhangProductPhysicalRowFingerprint(row);
	BOOST_CHECK(!fingerprint.empty());
	BOOST_CHECK(fingerprint.find("L2W|EFGH|G02|V4=-1") != std::string::npos);
	BOOST_CHECK(fingerprint.find("L1C|ABCD|G03|V1=1") != std::string::npos);
	BOOST_CHECK(fingerprint.find("G04") == std::string::npos);
	BOOST_CHECK_EQUAL(std::count(
		fingerprint.begin(), fingerprint.end(), ';'), 3);
}

BOOST_AUTO_TEST_CASE(product_integer_ledger_upgrades_conditioner_to_pair_certificate)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow conditioner;
	conditioner.system = E_Sys::GPS;
	conditioner.firstObservable = E_ObsCode::L1C;
	conditioner.secondObservable = E_ObsCode::L2W;
	conditioner.productRow = {1, 0, -1, 0};
	conditioner.integerValue = 8;
	conditioner.physicalExpansion = {{"L1C|A|G02|V0", 1},
		{"L1C|A|G03|V0", -1}};
	conditioner.canonicalProductExpansion = {
		{"L1C|G02", 1}, {"L1C|G03", -1}};
	conditioner.phaseSegmentFingerprint = "segments-v1";
	const auto first = ledger.observe(100, {conditioner}, 2);
	BOOST_REQUIRE(first.valid);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(ledger.rows().front().conditioningOnly);
	BOOST_CHECK(!ledger.rows().front().pairCertificate);

	auto pair = conditioner;
	pair.source = ZhangProductIntegerLedgerSource::DERIVED_PAIR;
	pair.conditioningOnly = false;
	pair.pairCertificate = true;
	pair.coordinate = "WL";
	pair.firstSatellite = "G02";
	pair.secondSatellite = "G03";
	const auto second = ledger.observe(130, {pair}, 2);
	BOOST_REQUIRE(second.valid);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(ledger.rows().front().certified);
	BOOST_CHECK(ledger.rows().front().pairCertificate);
	BOOST_CHECK(!ledger.rows().front().conditioningOnly);
	BOOST_CHECK_EQUAL(zhangProductIntegerLedgerSourceName(
		ledger.rows().front().source), "DERIVED_PAIR");
	BOOST_CHECK_EQUAL(ledger.rows().front().coordinate, "WL");
	BOOST_CHECK_EQUAL(ledger.rows().front().firstSatellite, "G02");
	BOOST_CHECK_EQUAL(ledger.rows().front().secondSatellite, "G03");
}

BOOST_AUTO_TEST_CASE(product_ledger_exact_affine_membership_separates_conflict)
{
	const ZhangExactMatrix rows{{1, -1, 0}, {0, 1, -1}};
	const ZhangExactVector values{3, 4};
	const auto redundant = zhangExactAffineMembership(
		rows, values, ZhangExactVector{1, 0, -1}, ZhangExactInteger(7));
	BOOST_CHECK(redundant.status ==
		ZhangExactAffineMembershipStatus::EXACT_REDUNDANT);
	BOOST_CHECK_EQUAL(redundant.impliedValue, 7);
	const auto conflict = zhangExactAffineMembership(
		rows, values, ZhangExactVector{1, 0, -1}, ZhangExactInteger(8));
	BOOST_CHECK(conflict.status ==
		ZhangExactAffineMembershipStatus::AFFINE_CONFLICT);
	BOOST_CHECK_EQUAL(conflict.impliedValue, 7);
	const auto fresh = zhangExactAffineMembership(
		rows, values, ZhangExactVector{1, 1, 1}, ZhangExactInteger(0));
	BOOST_CHECK(fresh.status == ZhangExactAffineMembershipStatus::EXACT_NEW);
}

BOOST_AUTO_TEST_CASE(integer_lattice_batch_membership_matches_scalar_queries)
{
	const ZhangExactMatrix rows{{2, 0, 2}, {0, 3, 3}, {1, 1, 1}};
	const ZhangExactMatrix targets{
		{2, 3, 5},
		{1, 1, 1},
		{0, 0, 1},
		{4, -3, 1}};
	const auto batch = zhangIntegerRowLatticeContainsBatch(rows, targets);
	BOOST_REQUIRE_EQUAL(batch.size(), targets.size());
	for (std::size_t target = 0; target < targets.size(); target++)
	{
		const auto scalar = zhangIntegerRowLatticeContains(rows, targets[target]);
		BOOST_CHECK_EQUAL(batch[target].contained, scalar.contained);
		BOOST_CHECK_EQUAL(batch[target].rank, scalar.rank);
		BOOST_CHECK(batch[target].smithInvariants == scalar.smithInvariants);
		if (!scalar.contained) continue;
		BOOST_REQUIRE_EQUAL(batch[target].combination.size(), rows.size());
		BOOST_CHECK(zhangExactRowCombination(
			batch[target].combination, rows) == targets[target]);
	}
}

BOOST_AUTO_TEST_CASE(posterior_effective_basis_excludes_deterministic_integer_row)
{
	MatrixXd covariance = MatrixXd::Zero(3, 3);
	covariance(0, 0) = 4;
	covariance(1, 1) = 1;
	MatrixXd rows(3, 3);
	rows << 1, 0, 0,
		0, 0, 1,
		2, 0, 0;
	const auto basis = zhangPosteriorEffectiveIntegerRowBasis(covariance, rows);
	BOOST_REQUIRE(basis.valid);
	BOOST_CHECK_EQUAL(basis.posteriorRank, 1);
	BOOST_REQUIRE_EQUAL(basis.retainedRows.size(), 1);
	BOOST_CHECK_EQUAL(basis.retainedRows.front(), 0);
	BOOST_REQUIRE_EQUAL(basis.deterministicRows.size(), 2);
	BOOST_CHECK_EQUAL(basis.deterministicRows[0], 1);
	BOOST_CHECK_EQUAL(basis.deterministicRows[1], 2);
}

BOOST_AUTO_TEST_CASE(posterior_conditionally_deterministic_affine_row_is_consistent)
{
	MatrixXd covariance = MatrixXd::Zero(2, 2);
	covariance(0, 0) = 1;
	MatrixXd rows(2, 2);
	rows << 1, 0,
		1, 1;
	VectorXd innovation(2); innovation << 1, 1;
	const auto basis = zhangPosteriorEffectiveIntegerRowBasis(covariance, rows);
	BOOST_REQUIRE(basis.valid);
	BOOST_CHECK_EQUAL(basis.posteriorRank, 1);
	BOOST_REQUIRE_EQUAL(basis.deterministicRows.size(), 1);
	BOOST_CHECK_EQUAL(basis.deterministicRows.front(), 1);
	const auto audit = zhangPosteriorConditionalInnovationResiduals(
		covariance, rows, innovation, basis);
	BOOST_REQUIRE(audit.valid);
	BOOST_CHECK_SMALL(audit.residuals(1), 1e-12);
	BOOST_CHECK_SMALL(audit.maximumAbsoluteResidual, 1e-12);
}

BOOST_AUTO_TEST_CASE(posterior_fully_determined_certificate_needs_no_conditioner)
{
	const MatrixXd covariance = MatrixXd::Zero(2, 2);
	MatrixXd rows(2, 2);
	rows << 1, 0,
		0, 1;
	const VectorXd consistent = VectorXd::Zero(2);
	const auto basis = zhangPosteriorEffectiveIntegerRowBasis(covariance, rows);
	BOOST_REQUIRE(basis.valid);
	BOOST_CHECK_EQUAL(basis.posteriorRank, 0);
	BOOST_CHECK(basis.retainedRows.empty());
	BOOST_REQUIRE_EQUAL(basis.deterministicRows.size(), 2);
	const auto accepted = zhangPosteriorConditionalInnovationResiduals(
		covariance, rows, consistent, basis);
	BOOST_REQUIRE(accepted.valid);
	BOOST_CHECK_SMALL(accepted.maximumAbsoluteResidual, 1e-12);
	VectorXd inconsistent(2); inconsistent << 0, 1;
	const auto rejected = zhangPosteriorConditionalInnovationResiduals(
		covariance, rows, inconsistent, basis);
	BOOST_REQUIRE(rejected.valid);
	BOOST_CHECK_CLOSE(rejected.maximumAbsoluteResidual, 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(
	posterior_effective_conditioner_keeps_consistent_deterministic_rows)
{
	VectorXd mean(3); mean << 0, 2, 7;
	MatrixXd covariance = MatrixXd::Zero(3, 3);
	covariance(0, 0) = 4;
	covariance(2, 2) = 1;
	MatrixXd rows(3, 3);
	rows << 1, 0, 0,
		0, 1, 0,
		1, 1, 0;
	VectorXd integers(3); integers << 3, 2, 5;

	const auto conditioned = zhangConditionPosteriorEffectiveIntegers(
		mean, covariance, rows, integers);
	BOOST_REQUIRE(conditioned.valid);
	BOOST_CHECK(conditioned.conditioned);
	BOOST_CHECK_EQUAL(conditioned.conditionerRank, 1);
	BOOST_CHECK_EQUAL(conditioned.deterministicRows, 2);
	BOOST_CHECK_SMALL(conditioned.maximumDeterministicResidual, 1e-12);
	BOOST_CHECK_CLOSE(conditioned.mean(0), 3.0, 1e-12);
	BOOST_CHECK_CLOSE(conditioned.mean(1), 2.0, 1e-12);
	BOOST_CHECK_CLOSE(conditioned.mean(2), 7.0, 1e-12);
	BOOST_CHECK_SMALL(conditioned.covariance(0, 0), 1e-12);
	BOOST_CHECK_CLOSE(conditioned.covariance(2, 2), 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(
	posterior_effective_conditioner_rejects_deterministic_affine_conflict)
{
	VectorXd mean(2); mean << 0, 2;
	MatrixXd covariance = MatrixXd::Zero(2, 2);
	covariance(0, 0) = 1;
	MatrixXd rows(2, 2);
	rows << 1, 0,
		1, 1;
	VectorXd integers(2); integers << 3, 6;

	const auto rejected = zhangConditionPosteriorEffectiveIntegers(
		mean, covariance, rows, integers);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.conditionerRank, 1);
	BOOST_CHECK_EQUAL(rejected.deterministicRows, 1);
	BOOST_CHECK_CLOSE(rejected.maximumDeterministicResidual, 1.0, 1e-12);
	BOOST_CHECK_EQUAL(rejected.failureReason,
		"POSTERIOR_DETERMINISTIC_AFFINE_CONFLICT");
}

BOOST_AUTO_TEST_CASE(
	posterior_effective_conditioner_handles_canonical_56d_rank44_psd_case)
{
	constexpr int dimension = 56;
	constexpr int stochasticRank = 44;
	VectorXd mean = VectorXd::Zero(dimension);
	MatrixXd covariance = MatrixXd::Zero(dimension, dimension);
	for (int index = 0; index < stochasticRank; index++)
	{
		mean(index) = 0.25;
		covariance(index, index) = 0.5 + 0.01 * index;
	}
	const MatrixXd rows = MatrixXd::Identity(dimension, dimension);
	const VectorXd integers = VectorXd::Zero(dimension);

	const auto basis = zhangPosteriorEffectiveIntegerRowBasis(
		covariance, rows);
	BOOST_REQUIRE(basis.valid);
	BOOST_CHECK_EQUAL(basis.posteriorRank, stochasticRank);
	BOOST_CHECK_EQUAL(basis.retainedRows.size(), stochasticRank);
	BOOST_CHECK_EQUAL(basis.deterministicRows.size(),
		dimension - stochasticRank);

	const auto conditioned = zhangConditionPosteriorEffectiveIntegers(
		mean, covariance, rows, integers);
	BOOST_REQUIRE(conditioned.valid);
	BOOST_CHECK(conditioned.conditioned);
	BOOST_CHECK_EQUAL(conditioned.conditionerRank, stochasticRank);
	BOOST_CHECK_EQUAL(conditioned.deterministicRows,
		dimension - stochasticRank);
	BOOST_CHECK_SMALL(conditioned.maximumDeterministicResidual, 1e-12);
	BOOST_CHECK_SMALL(conditioned.mean.head(stochasticRank).norm(), 1e-12);
	BOOST_CHECK_SMALL(conditioned.covariance.norm(), 1e-10);

	VectorXd conflicting = integers;
	conflicting(dimension - 1) = 1;
	const auto rejected = zhangConditionPosteriorEffectiveIntegers(
		mean, covariance, rows, conflicting);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.failureReason,
		"POSTERIOR_DETERMINISTIC_AFFINE_CONFLICT");
}

BOOST_AUTO_TEST_CASE(product_integer_ledger_transports_backend_generations_canonically)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.firstObservable = E_ObsCode::L1C;
	row.secondObservable = E_ObsCode::L2W;
	row.productRow = {1, -1};
	row.integerValue = 7;
	row.physicalExpansion = {
		{"L1C|ABCD|G03|V4", 1}, {"L1C|ABCD|G02|V2", -1}};
	row.canonicalProductExpansion = {
		{"L1C|G03", 1}, {"L1C|G02", -1}};
	row.phaseSegmentFingerprint = "G02|L1C|SEG1;G03|L1C|SEG1;";
	row.backendBasisGeneration = 12;
	BOOST_REQUIRE(ledger.observe(100, {row}, 1).valid);

	auto nextGeneration = row;
	nextGeneration.backendBasisGeneration = 13;
	nextGeneration.integerValue = -112;
	const auto update = ledger.observe(130, {nextGeneration}, 1);
	BOOST_CHECK(!update.valid);
	BOOST_CHECK_EQUAL(update.freshRows, 0);
	BOOST_CHECK_EQUAL(update.conflictingRows, 1);
	BOOST_CHECK_EQUAL(update.failureReason,
		"PRODUCT_LEDGER_MATURE_INTEGER_CONFLICT");
	BOOST_CHECK(!zhangProductLedgerWriterCommitAuthorized(update));

	auto nextSegment = nextGeneration;
	nextSegment.phaseSegmentFingerprint =
		"G02|L1C|SEG2;G03|L1C|SEG2;";
	nextSegment.integerValue = 32;
	const auto segmentUpdate = ledger.observe(160, {nextSegment}, 1);
	BOOST_REQUIRE(segmentUpdate.valid);
	BOOST_CHECK_EQUAL(segmentUpdate.freshRows, 1);
	BOOST_CHECK_EQUAL(segmentUpdate.conflictingRows, 0);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 2);

	const auto generation12 = ledger.rowsForGeneration(12);
	const auto generation13 = ledger.rowsForGeneration(13);
	BOOST_REQUIRE_EQUAL(generation12.size(), 1);
	BOOST_REQUIRE_EQUAL(generation13.size(), 1);
	BOOST_CHECK_EQUAL(generation12[0].integerValue, -7);
	BOOST_CHECK_EQUAL(generation13[0].integerValue, -32);

	// The current physical chart is retained after the canonical identity is
	// transported across a backend generation.
	ZhangExactVector projected;
	const std::map<std::string, int> currentColumns = {
		{"L1C|ABCD|G02|V2", 0}, {"L1C|ABCD|G03|V4", 1}};
	BOOST_REQUIRE(zhangProjectProductLedgerPhysicalRow(
		generation13.front(), currentColumns, 2, projected));
	BOOST_REQUIRE_EQUAL(projected.size(), 2);
	BOOST_CHECK_EQUAL(projected[0], 1);
	BOOST_CHECK_EQUAL(projected[1], -1);
	const std::map<std::string, int> missingArc = {
		{"L1C|ABCD|G03|V4", 0}};
	BOOST_CHECK(!zhangProjectProductLedgerPhysicalRow(
		generation13.front(), missingArc, 1, projected));
}

BOOST_AUTO_TEST_CASE(r33_product_ledger_rebinds_same_canonical_row_to_current_chart)
{
	// Physical coordinates may change while the exact product statement does not.
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow first;
	first.system = E_Sys::GPS;
	first.firstObservable = E_ObsCode::L1C;
	first.secondObservable = E_ObsCode::L2W;
	first.productRow = {1, -1};
	first.integerValue = 9;
	first.physicalExpansion = {{"L1C|A|G01|V1", 1}, {"L1C|B|G02|V1", -1}};
	first.canonicalProductExpansion = {{"L1C|G01", 1}, {"L1C|G02", -1}};
	first.phaseSegmentFingerprint = "G01|L1C|SEG2;G02|L1C|SEG3;";
	first.backendBasisGeneration = 7;
	BOOST_REQUIRE(ledger.observe(100, {first}, 2).valid);

	auto rebased = first;
	rebased.productRow = {-1, 0, 1, 0};
	rebased.physicalExpansion = {{"L1C|C|G01|V4", 1}, {"L1C|D|G02|V6", -1}};
	rebased.backendBasisGeneration = 8;
	const auto update = ledger.observe(130, {rebased}, 2);
	BOOST_REQUIRE(update.valid);
	BOOST_CHECK_EQUAL(update.freshRows, 0);
	BOOST_CHECK_EQUAL(update.conflictingRows, 0);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	BOOST_CHECK(ledger.rows().front().certified);
	BOOST_CHECK_EQUAL(ledger.rows().front().backendBasisGeneration, 8);
	BOOST_CHECK(ledger.rows().front().physicalExpansion == rebased.physicalExpansion);
}

BOOST_AUTO_TEST_CASE(product_lattice_failure_probability_respects_decimal_budget)
{
	const double configuredBudget = 1e-3;
	const double configuredSuccess = 0.999;
	const double bound = zhangProductFailureProbabilityBound(
		configuredSuccess, configuredBudget);
	BOOST_CHECK_LE(bound, configuredBudget);
	BOOST_CHECK(zhangProductFailureProbabilityPassed(
		bound, configuredBudget));

	const double firstStage = zhangProductFailureProbabilityBound(
		0.9996, configuredBudget);
	const double remaining = configuredBudget - firstStage;
	const double secondStage = zhangProductFailureProbabilityBound(
		1 - remaining, remaining);
	BOOST_CHECK(zhangProductFailureProbabilityPassed(
		firstStage + secondStage, configuredBudget));

	const double unsafe = zhangProductFailureProbabilityBound(
		0.9, configuredBudget);
	BOOST_CHECK_CLOSE_FRACTION(unsafe, 0.1, 1e-14);
	BOOST_CHECK(!zhangProductFailureProbabilityPassed(
		unsafe, configuredBudget));
	BOOST_CHECK_EQUAL(
		zhangProductFailureProbabilityBound(1.01, configuredBudget), 1);
}

BOOST_AUTO_TEST_CASE(product_candidate_pair_rank_requires_exact_named_edge)
{
	BOOST_CHECK(zhangProductCandidateIsNamedPairRow({1, 0, 0}));
	BOOST_CHECK(zhangProductCandidateIsNamedPairRow({1, -1, 0}));
	BOOST_CHECK(zhangProductCandidateIsNamedPairRow({-1, 0, 0}));
	BOOST_CHECK(!zhangProductCandidateIsNamedPairRow({2, -1, 0}));
	BOOST_CHECK(!zhangProductCandidateIsNamedPairRow({1, 1, 0}));
	BOOST_CHECK(!zhangProductCandidateIsNamedPairRow({1, -1, 1}));
}

BOOST_AUTO_TEST_CASE(lambda_reports_selected_suffix_bootstrap_success)
{
	VectorXd conditionalVariances(3);
	conditionalVariances << 4, 0.04, 0.01;
	const double expected =
		std::erf(std::sqrt(1 / (8 * 0.04))) *
		std::erf(std::sqrt(1 / (8 * 0.01)));
	const double selected = lambdaSelectedSuffixBootstrapSuccess(
		conditionalVariances, 2);
	BOOST_CHECK_CLOSE_FRACTION(
		selected, expected, 1e-14);
	BOOST_CHECK_GT(selected,
		lambdaSelectedSuffixBootstrapSuccess(conditionalVariances, 3));
	BOOST_CHECK_EQUAL(
		lambdaSelectedSuffixBootstrapSuccess(conditionalVariances, 0), 0);
}

BOOST_AUTO_TEST_CASE(full_product_lattice_oracle_parser_requires_admissible_complete_rank)
{
	const std::string valid = R"json({
		"status":"FULL_ORACLE_READY",
		"hard_gate_passed":true,
		"oracle":{
			"schema":"ZHANG_FULL_PRODUCT_LATTICE_ORACLE_V1",
			"system":"GPS",
			"reference_satellite":"G02",
			"satellites":["G02","G03","G04"],
			"dual_frequency_rank":2,
			"relations":[
				{"satellite":"G03","reference":"G02",
				 "wl_satellite_minus_reference":3,
				 "l1_satellite_minus_reference":8,
				 "l2_satellite_minus_reference":5},
				{"satellite":"G04","reference":"G02",
				 "wl_satellite_minus_reference":1,
				 "l1_satellite_minus_reference":12,
				 "l2_satellite_minus_reference":11}
			]
		}
	})json";
	const auto oracle = parseZhangFullProductLatticeOracle(valid, 2);
	BOOST_REQUIRE(oracle.valid);
	BOOST_CHECK_EQUAL(oracle.rank, 2);
	BOOST_CHECK_EQUAL(oracle.referenceSatellite, "G02");
	BOOST_REQUIRE_EQUAL(oracle.potentials.size(), 3);
	BOOST_CHECK(oracle.potentials.at("G03").wideLane == 3);
	BOOST_CHECK(oracle.potentials.at("G04").secondSignal == 11);

	std::string inadmissible = valid;
	const auto position = inadmissible.find(
		"\"l2_satellite_minus_reference\":5");
	BOOST_REQUIRE(position != std::string::npos);
	inadmissible.replace(position,
		std::string("\"l2_satellite_minus_reference\":5").size(),
		"\"l2_satellite_minus_reference\":6");
	const auto rejected = parseZhangFullProductLatticeOracle(inadmissible, 2);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.failureReason,
		"ORACLE_WL_L1_L2_NOT_ADMISSIBLE");
	const auto incomplete = parseZhangFullProductLatticeOracle(valid);
	BOOST_CHECK(!incomplete.valid);
	BOOST_CHECK_EQUAL(incomplete.failureReason,
		"ORACLE_EXPECTED_FULL_RANK_MISMATCH");
}

BOOST_AUTO_TEST_CASE(integer_support_quality_fails_closed_without_residuals)
{
	ZhangIntegerArcQuality quality;
	quality.ageEpochs = 100;
	quality.observations = 100;
	const auto missing = zhangEvaluateIntegerSupportQuality(
		quality, ZhangIntegerSupportQualityGates{});
	BOOST_CHECK(!missing.eligibleForIntegerSupport);
	BOOST_CHECK_EQUAL(missing.failureReason, "PHASE_RMS_GATE_FAILED");

	quality.phaseResidualRms = 0.01;
	quality.codeResidualRms = 1;
	quality.phaseResidualMad = 0.01;
	quality.codeResidualMad = 1;
	quality.elevationScore = 0.5;
	quality.whitenedResidualScore = 1;
	const auto accepted = zhangEvaluateIntegerSupportQuality(
		quality, ZhangIntegerSupportQualityGates{});
	BOOST_CHECK(accepted.eligibleForIntegerSupport);
	BOOST_CHECK_EQUAL(accepted.failureReason, "ELIGIBLE");
}

BOOST_AUTO_TEST_CASE(integer_support_residual_audit_uses_postfit_phase_and_code_rows)
{
	int ownerToken = 0;
	KFMeas measurements;
	measurements.obsKeys.resize(2);
	measurements.obsKeys[0].type = KF::PHAS_MEAS;
	measurements.obsKeys[0].str = "R0";
	measurements.obsKeys[0].Sat = SatSys(E_Sys::GPS, 3);
	measurements.obsKeys[1].type = KF::CODE_MEAS;
	measurements.obsKeys[1].str = "R0";
	measurements.obsKeys[1].Sat = SatSys(E_Sys::GPS, 3);
	measurements.VV.resize(2);
	measurements.VV << 0.01, 0.8;
	measurements.postfitRatios.resize(2);
	measurements.postfitRatios << 0.5, 1.5;
	zhangRecordIntegerSupportPostfitResidual(&ownerToken, measurements, 0);
	zhangRecordIntegerSupportPostfitResidual(&ownerToken, measurements, 1);
	const auto summary = zhangIntegerSupportResidualSummary(
		&ownerToken, {"R0", SatSys(E_Sys::GPS, 3)});
	BOOST_CHECK_EQUAL(summary.phaseSamples, 1);
	BOOST_CHECK_EQUAL(summary.codeSamples, 1);
	BOOST_CHECK_CLOSE(summary.phaseRms, 0.01, 1e-10);
	BOOST_CHECK_CLOSE(summary.codeRms, 0.8, 1e-10);
	BOOST_CHECK_CLOSE(summary.maximumWhitenedResidualScore, 1.5, 1e-10);
}

BOOST_AUTO_TEST_CASE(exact_held_quotient_separates_certified_and_unresolved_rank)
{
	// Target is the full four-dimensional integer lattice.  Held evidence
	// certifies two primitive target directions plus one unrelated ambient row.
	const ZhangExactMatrix target = {
		{1, 0, 0, 0, 0},
		{0, 1, 0, 0, 0},
		{0, 0, 1, 0, 0},
		{0, 0, 0, 1, 0}};
	const ZhangExactMatrix held = {
		{1, 0, 0, 0, 0},
		{0, 1, 0, 0, 0},
		{0, 0, 0, 0, 1}};
	const ZhangExactVector values = {7, -3, 99};

	const auto audit = zhangExactHeldQuotientAudit(target, held, values);
	BOOST_REQUIRE(audit.valid);
	BOOST_CHECK_EQUAL(audit.targetRank, 4);
	BOOST_CHECK_EQUAL(audit.heldIntersectionRank, 2);
	BOOST_CHECK_EQUAL(audit.quotientRank, 2);
	BOOST_CHECK(audit.heldIntersectionPrimitiveInTarget);
	BOOST_CHECK(audit.exactClosure);
	BOOST_REQUIRE_EQUAL(audit.heldIntersectionValues.size(), 2);
	BOOST_CHECK_EQUAL(audit.heldIntersectionValues[0], 7);
	BOOST_CHECK_EQUAL(audit.heldIntersectionValues[1], -3);
}

BOOST_AUTO_TEST_CASE(exact_held_quotient_rejects_nonprimitive_intersection)
{
	const ZhangExactMatrix target = {{1, 0}, {0, 1}};
	const ZhangExactMatrix held = {{2, 0}};
	const auto audit = zhangExactHeldQuotientAudit(target, held, {4});
	BOOST_CHECK(!audit.valid);
	BOOST_CHECK_EQUAL(
		audit.failureReason, "HELD_INTERSECTION_NOT_PRIMITIVE_IN_TARGET");
}

BOOST_AUTO_TEST_CASE(exact_dual_product_intersection_finds_maximum_common_lattice)
{
	const ZhangExactMatrix network = {
		{1, 0, 0},
		{0, 1, 0}};
	const ZhangExactVector values = {5, 7};
	const ZhangExactMatrix first = {
		{1, 0, 0},
		{0, 1, 0}};
	const ZhangExactMatrix second = {
		{1, 0, 0},
		{0, 0, 1}};

	const auto audit = zhangExactDualProductIntersection(
		network, values, first, second);
	BOOST_TEST(audit.valid);
	BOOST_TEST(audit.exactClosure);
	BOOST_TEST(audit.networkRank == 2);
	BOOST_TEST(audit.fullProductRank == 2);
	BOOST_TEST(audit.firstSignalIntersectionRank == 2);
	BOOST_TEST(audit.wideLaneIntersectionRank == 1);
	BOOST_TEST(audit.dualIntersectionRank == 1);
	BOOST_REQUIRE(audit.productRows.size() == 1);
	BOOST_TEST(audit.productRows.front()[0] == 1);
	BOOST_TEST(audit.productRows.front()[1] == 0);
	BOOST_TEST(audit.firstSignalIntegers.front() == 5);
	BOOST_TEST(audit.wideLaneIntegers.front() == 0);
}

BOOST_AUTO_TEST_CASE(exact_dual_product_intersection_preserves_full_rank)
{
	const ZhangExactMatrix network = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}};
	const ZhangExactMatrix first = {
		{1, 0, 0},
		{0, 1, 0}};
	const ZhangExactMatrix second = {
		{0, 0, 1},
		{1, 0, 0}};
	const auto audit = zhangExactDualProductIntersection(
		network, {3, 4, 5}, first, second);
	BOOST_TEST(audit.valid);
	BOOST_TEST(audit.dualIntersectionRank == 2);
	BOOST_TEST(audit.productRows.size() == 2);
	BOOST_TEST(audit.firstSignalIntegers.size() == 2);
	BOOST_TEST(audit.wideLaneIntegers.size() == 2);
	BOOST_TEST(audit.productLatticePrimitive);
	BOOST_REQUIRE_EQUAL(audit.productSmithInvariants.size(), 2);
	BOOST_TEST(zhangExactAbs(audit.productSmithInvariants[0]) == 1);
	BOOST_TEST(zhangExactAbs(audit.productSmithInvariants[1]) == 1);
}

BOOST_AUTO_TEST_CASE(exact_affine_row_lattice_intersection_preserves_both_integers)
{
	const ZhangExactMatrix wideLaneRows = {
		{1, 0, 0},
		{0, 1, 0}};
	const ZhangExactVector wideLaneIntegers = {5, 7};
	const ZhangExactMatrix firstSignalRows = {
		{0, 1, 0},
		{0, 0, 1}};
	const ZhangExactVector firstSignalIntegers = {11, 13};

	const auto intersection = zhangExactAffineRowLatticeIntersection(
		wideLaneRows, wideLaneIntegers,
		firstSignalRows, firstSignalIntegers);
	BOOST_REQUIRE(intersection.valid);
	BOOST_CHECK(intersection.exactClosure);
	BOOST_CHECK(intersection.primitive);
	BOOST_CHECK_EQUAL(intersection.firstRank, 2);
	BOOST_CHECK_EQUAL(intersection.secondRank, 2);
	BOOST_REQUIRE_EQUAL(intersection.intersectionRank, 1);
	BOOST_REQUIRE_EQUAL(intersection.rows.size(), 1);
	BOOST_CHECK(intersection.rows[0] == ZhangExactVector({0, 1, 0}));
	BOOST_REQUIRE_EQUAL(intersection.firstValues.size(), 1);
	BOOST_REQUIRE_EQUAL(intersection.secondValues.size(), 1);
	BOOST_CHECK_EQUAL(intersection.firstValues[0], 7);
	BOOST_CHECK_EQUAL(intersection.secondValues[0], 11);
}

BOOST_AUTO_TEST_CASE(exact_affine_row_lattice_intersection_rejects_hidden_conflict)
{
	const ZhangExactMatrix conflictingRows = {
		{1, 0},
		{1, 0}};
	const auto intersection = zhangExactAffineRowLatticeIntersection(
		conflictingRows, {3, 4}, ZhangExactMatrix{{1, 0}}, {3});
	BOOST_CHECK(!intersection.valid);
	BOOST_CHECK_EQUAL(intersection.failureReason,
		"AFFINE_LATTICE_INTERSECTION_INPUT_AFFINE_CONFLICT");
}

BOOST_AUTO_TEST_CASE(
	exact_affine_agreement_keeps_combination_when_basis_rows_disagree)
{
	const ZhangExactMatrix rows = {{1, 0}, {0, 1}};
	const auto agreement = zhangExactAffineAgreementLattice(
		rows, {10, 20}, rows, {11, 19});
	BOOST_REQUIRE(agreement.valid);
	BOOST_CHECK(agreement.exactClosure);
	BOOST_CHECK_EQUAL(agreement.commonRank, 2);
	BOOST_CHECK_EQUAL(agreement.agreementRank, 1);
	BOOST_CHECK_EQUAL(agreement.conflictRank, 1);
	BOOST_REQUIRE_EQUAL(agreement.rows.size(), 1);
	BOOST_REQUIRE_EQUAL(agreement.values.size(), 1);
	const auto firstMembership = zhangIntegerRowLatticeContains(
		rows, agreement.rows.front());
	BOOST_REQUIRE(firstMembership.contained);
	ZhangExactInteger firstValue = 0;
	ZhangExactInteger secondValue = 0;
	for (std::size_t i = 0; i < firstMembership.combination.size(); i++)
	{
		firstValue += firstMembership.combination[i] *
			ZhangExactVector{10, 20}[i];
		secondValue += firstMembership.combination[i] *
			ZhangExactVector{11, 19}[i];
	}
	BOOST_CHECK_EQUAL(firstValue, secondValue);
	BOOST_CHECK_EQUAL(firstValue, agreement.values.front());
}

BOOST_AUTO_TEST_CASE(
	minimal_affine_conflict_isolates_one_edge_and_preserves_complement)
{
	// Candidate 0 contradicts the baseline integer on the first component.
	// Candidate 1 is an independent, conflict-free component and must remain a
	// valid writer constraint rather than being cleared with the whole batch.
	const ZhangExactMatrix baselineRows = {{1, 0, 0}, {0, 1, 0}};
	const ZhangExactVector baselineValues = {4, 5};
	const ZhangExactMatrix candidateRows = {{1, 0, 0}, {0, 0, 1}};
	const ZhangExactVector candidateValues = {6, 7};
	const auto conflict = zhangMinimalAffineConflictCertificate(
		baselineRows, baselineValues, candidateRows, candidateValues);
	BOOST_REQUIRE(conflict.valid);
	BOOST_CHECK(conflict.baselineConsistent);
	BOOST_CHECK(!conflict.combinedConsistent);
	BOOST_CHECK(conflict.minimal);
	BOOST_REQUIRE_EQUAL(conflict.candidateIndices.size(), 1);
	BOOST_CHECK_EQUAL(conflict.candidateIndices.front(), 0);

	auto survivingRows = baselineRows;
	auto survivingValues = baselineValues;
	survivingRows.push_back(candidateRows[1]);
	survivingValues.push_back(candidateValues[1]);
	const auto surviving = zhangExactRowHermiteNormalForm(
		survivingRows, survivingValues);
	BOOST_REQUIRE(surviving.consistent);
	BOOST_CHECK_EQUAL(surviving.basis.size(), 3);
	const auto retained = zhangExactAffineMembership(
		baselineRows, baselineValues, candidateRows[1], candidateValues[1]);
	BOOST_CHECK(retained.status == ZhangExactAffineMembershipStatus::EXACT_NEW);
	ZhangProductIntegerLedgerUpdate isolatedWriterUpdate;
	isolatedWriterUpdate.valid = true;
	isolatedWriterUpdate.inputRows = 1;
	isolatedWriterUpdate.confirmedRows = 1;
	isolatedWriterUpdate.activeRankAfter = 1;
	isolatedWriterUpdate.conflictingRows = 0;
	BOOST_CHECK(zhangProductLedgerWriterCommitAuthorized(isolatedWriterUpdate));
	BOOST_CHECK_GT(isolatedWriterUpdate.activeRankAfter, 0);
}

BOOST_AUTO_TEST_CASE(exact_dual_intersection_recovers_implied_pair_certificates)
{
	const ZhangExactMatrix rows = {
		{1, 1, 0},
		{0, 1, 1}};
	const auto intersection = zhangExactAffineRowLatticeIntersection(
		rows, {5, 8}, rows, {11, 13});
	BOOST_REQUIRE(intersection.valid);
	const auto wideLanePairs = zhangRecoverCertifiedPairRelations(
		intersection.rows, intersection.firstValues, 3, true);
	const auto firstSignalPairs = zhangRecoverCertifiedPairRelations(
		intersection.rows, intersection.secondValues, 3, true);
	auto pairValue = [](const auto& pairs, int first, int second)
		-> std::optional<ZhangExactInteger>
	{
		for (const auto& pair : pairs)
			if (pair.firstNode == first && pair.secondNode == second)
				return pair.value;
		return std::nullopt;
	};
	const auto wideLane02 = pairValue(wideLanePairs, 0, 2);
	const auto firstSignal02 = pairValue(firstSignalPairs, 0, 2);
	BOOST_REQUIRE(wideLane02.has_value());
	BOOST_REQUIRE(firstSignal02.has_value());
	BOOST_CHECK_EQUAL(*wideLane02, -3);
	BOOST_CHECK_EQUAL(*firstSignal02, -2);
}

BOOST_AUTO_TEST_CASE(exact_dual_full_rank_does_not_hide_nonprimitive_index)
{
	const ZhangExactMatrix network = {{2, 0}, {0, 2}};
	const ZhangExactMatrix first = {{1, 0}, {0, 1}};
	const ZhangExactMatrix second = {{0, 0}, {0, 0}};
	const auto audit = zhangExactDualProductIntersection(
		network, {6, 10}, first, second);
	BOOST_REQUIRE(audit.valid);
	BOOST_CHECK_EQUAL(audit.dualIntersectionRank, 2);
	BOOST_CHECK(!audit.productLatticePrimitive);
	BOOST_REQUIRE_EQUAL(audit.productSmithInvariants.size(), 2);
	BOOST_CHECK_EQUAL(zhangExactAbs(audit.productSmithInvariants[0]), 2);
	BOOST_CHECK_EQUAL(zhangExactAbs(audit.productSmithInvariants[1]), 2);
}

BOOST_AUTO_TEST_CASE(partial_physical_intersection_embeds_without_certifying_missing_direction)
{
	ZhangExactMatrix rows{{1, -1}, {0, 1}};
	BOOST_REQUIRE(zhangEmbedCanonicalSubsetRows(rows, {0, 2}, 4));
	BOOST_REQUIRE_EQUAL(rows.size(), 2);
	BOOST_CHECK(rows[0] == ZhangExactVector({1, 0, -1, 0}));
	BOOST_CHECK(rows[1] == ZhangExactVector({0, 0, 1, 0}));
	for (const auto& row : rows)
	{
		BOOST_CHECK_EQUAL(row[1], 0);
		BOOST_CHECK_EQUAL(row[3], 0);
	}
	ZhangExactMatrix invalid{{1, 0}};
	BOOST_CHECK(!zhangEmbedCanonicalSubsetRows(invalid, {1, 1}, 4));
}

BOOST_AUTO_TEST_CASE(exact_certified_union_recomputes_rank_and_target_equality)
{
	const ZhangExactMatrix target = {
		{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	const ZhangExactMatrix held = {{1, 0, 0}};
	const ZhangExactVector heldValues = {5};
	const ZhangExactMatrix partialFixed = {{0, 1, 0}};
	const ZhangExactVector partialValues = {-2};
	const auto partial = zhangExactCertifiedUnionAudit(
		target, held, heldValues, partialFixed, partialValues);
	BOOST_CHECK(partial.consistent);
	BOOST_CHECK_EQUAL(partial.targetRank, 3);
	BOOST_CHECK_EQUAL(partial.heldRank, 1);
	BOOST_CHECK_EQUAL(partial.newlyFixedRank, 1);
	BOOST_CHECK_EQUAL(partial.combinedCertifiedRank, 2);
	BOOST_CHECK(!partial.exactTargetEquality);
	BOOST_CHECK_EQUAL(
		partial.failureReason, "TARGET_LATTICE_NOT_FULLY_CERTIFIED");

	const ZhangExactMatrix completeFixed = {{0, 1, 0}, {0, 0, 1}};
	const ZhangExactVector completeValues = {-2, 8};
	const auto complete = zhangExactCertifiedUnionAudit(
		target, held, heldValues, completeFixed, completeValues);
	BOOST_CHECK(complete.consistent);
	BOOST_CHECK_EQUAL(complete.combinedCertifiedRank, 3);
	BOOST_CHECK(complete.exactTargetEquality);
}

BOOST_AUTO_TEST_CASE(deterministic_quotient_distinguishes_integer_inconsistency)
{
	VectorXd integerMean(2); integerMean << 4, 1.25;
	MatrixXd covariance = MatrixXd::Zero(2, 2);
	covariance(1, 1) = 0.04;
	const auto consistent = zhangAuditDeterministicQuotientModes(
		integerMean, covariance);
	BOOST_REQUIRE(consistent.covarianceValid);
	BOOST_CHECK_EQUAL(consistent.covarianceRank, 1);
	BOOST_CHECK_EQUAL(consistent.nullity, 1);
	BOOST_CHECK(consistent.integerConsistent);
	BOOST_CHECK_EQUAL(consistent.status, "UNTRACKED_DETERMINISTIC_RELATION");

	VectorXd fractionalMean(2); fractionalMean << 4.25, 1.25;
	const auto inconsistent = zhangAuditDeterministicQuotientModes(
		fractionalMean, covariance);
	BOOST_REQUIRE(inconsistent.covarianceValid);
	BOOST_CHECK(!inconsistent.integerConsistent);
	BOOST_CHECK_EQUAL(
		inconsistent.status, "DETERMINISTIC_INTEGER_INCONSISTENCY");
}

BOOST_AUTO_TEST_CASE(integer_gain_frontier_is_reliability_first_and_exact_at_rank_one)
{
	VectorXd mean(2); mean << 3.01, -1.02;
	MatrixXd covariance = MatrixXd::Identity(2, 2) * 1e-4;
	MatrixXd products = MatrixXd::Identity(2, 2);
	MatrixXd productCross = products * covariance;
	const auto frontier = zhangBoundedIntegerProductGainFrontier(
		mean, covariance, productCross, 1, 1e-3, 1e-6, 2, 64,
		(products * covariance * products.transpose()).trace());
	BOOST_REQUIRE(frontier.valid);
	BOOST_CHECK_EQUAL(frontier.status, "COMPLETE");
	BOOST_CHECK(frontier.enumeratedPrimitiveRows > 0);
	BOOST_CHECK(frontier.reliablePrimitiveRows > 0);
	BOOST_REQUIRE(!frontier.points.empty());
	BOOST_CHECK_EQUAL(frontier.points.front().rank, 1);
	BOOST_CHECK(frontier.points.front().reliable);
	BOOST_CHECK(frontier.points.front().exactBoundedOptimum);
	BOOST_CHECK(frontier.points.front().failureProbabilityBound <= 1e-3);
}

BOOST_AUTO_TEST_CASE(integer_gain_frontier_rejects_fractional_low_variance_rows)
{
	VectorXd mean(1); mean << 0.25;
	MatrixXd covariance(1, 1); covariance << 1e-8;
	MatrixXd products(1, 1); products << 1;
	MatrixXd productCross = products * covariance;
	const auto frontier = zhangBoundedIntegerProductGainFrontier(
		mean, covariance, productCross, 2, 1e-3, 1e-6, 1, 128,
		(products * covariance * products.transpose()).trace());
	BOOST_REQUIRE(frontier.valid);
	BOOST_CHECK_EQUAL(frontier.status, "NO_RELIABLE_PRIMITIVE_ROW");
	BOOST_CHECK_EQUAL(frontier.reliablePrimitiveRows, 0);
	BOOST_CHECK(frontier.points.empty());
}

BOOST_AUTO_TEST_CASE(integer_gain_frontier_keeps_reliable_explicit_seed_beyond_bound)
{
	VectorXd mean(2); mean << 0.5, 0;
	MatrixXd covariance = MatrixXd::Identity(2, 2) * 1e-8;
	MatrixXd productCross(1, 2); productCross << 1e-8, 0;
	ZhangExactMatrix seeds = {{2, 1}};
	const auto frontier = zhangBoundedIntegerProductGainFrontier(
		mean, covariance, productCross, 1, 1e-3, 1e-6, 1, 64,
		1e-8, seeds);
	BOOST_REQUIRE(frontier.valid);
	BOOST_CHECK_EQUAL(frontier.explicitSeedRowsAdded, 1);
	BOOST_CHECK_EQUAL(frontier.reliableExplicitSeedRows, 1);
	BOOST_REQUIRE(!frontier.points.empty());
	BOOST_REQUIRE_EQUAL(frontier.points.front().rows.size(), 1);
	BOOST_CHECK_EQUAL(frontier.points.front().rows.front().at(0), 2);
	BOOST_CHECK_EQUAL(frontier.points.front().rows.front().at(1), 1);
}

BOOST_AUTO_TEST_CASE(integer_gain_frontier_sparse_support_reaches_high_dimension)
{
	VectorXd mean = VectorXd::Zero(12);
	MatrixXd covariance = MatrixXd::Identity(12, 12) * 1e-4;
	MatrixXd productCross = MatrixXd::Identity(12, 12) * 1e-4;
	const auto frontier = zhangBoundedIntegerProductGainFrontier(
		mean, covariance, productCross, 2, 1e-3, 1e-6, 2, 64,
		12e-4, {}, 2);
	BOOST_REQUIRE(frontier.valid);
	BOOST_CHECK_EQUAL(frontier.dimension, 12);
	BOOST_CHECK_EQUAL(frontier.maximumEnumerationSupport, 2);
	BOOST_CHECK(frontier.enumeratedPrimitiveRows > 0);
	BOOST_CHECK(frontier.enumeratedPrimitiveRows < 10000);
	BOOST_REQUIRE(!frontier.points.empty());
	BOOST_CHECK_EQUAL(frontier.points.front().rank, 1);
}

BOOST_AUTO_TEST_CASE(product_candidate_generator_adds_dense_real_mode_integer_rows)
{
	constexpr int dimension = 6;
	VectorXd mean = VectorXd::Zero(dimension);
	MatrixXd covariance = MatrixXd::Identity(dimension, dimension) * 1e-4;
	MatrixXd productCross(1, dimension);
	productCross.setOnes();
	productCross *= 1e-2;
	const auto generated = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6, {}, 1, 8, 256);
	BOOST_REQUIRE(generated.valid);
	BOOST_CHECK_EQUAL(generated.failureReason, "NONE");
	BOOST_CHECK(generated.realModeApproximations > 0);
	bool foundDense = false;
	for (const auto& candidate : generated.candidates)
	{
		if (candidate.source != "PRODUCT_GAIN_REAL_MODE_APPROXIMATION") continue;
		int support = 0;
		for (const auto& coefficient : candidate.row)
			support += coefficient != 0;
		if (support > 2 && candidate.reliabilityPassed)
		{
			foundDense = true;
			BOOST_CHECK(candidate.variance > 0);
			BOOST_CHECK(candidate.perr <= 1e-3);
			BOOST_CHECK(candidate.incrementalProductGain > 0);
			break;
		}
	}
	BOOST_CHECK_MESSAGE(foundDense,
		"real product-gain modes must generate legal dense primitive rows");
}

BOOST_AUTO_TEST_CASE(product_candidate_generator_orders_reliability_then_signal_graph_then_gain)
{
	VectorXd mean(3);
	mean << 0, 0, 0.25;
	MatrixXd covariance = MatrixXd::Identity(3, 3) * 1e-6;
	MatrixXd productCross(1, 3);
	productCross << 0.01, 0.01, 100;
	const auto generated = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6, {}, 0, 1, 128);
	BOOST_REQUIRE(generated.valid);
	BOOST_REQUIRE(!generated.candidates.empty());
	BOOST_CHECK(generated.candidates.front().reliabilityPassed);
	// The generator has no paired L1/WL certificate at this stage.  Its
	// graph contribution is therefore signal-specific; dual-frequency rank is
	// assigned only after the WL and conditional-L1 certificates are united.
	BOOST_CHECK_EQUAL(generated.candidates.front().dualGraphRankGain, 0);
	BOOST_CHECK_EQUAL(generated.candidates.front().signalGraphRankGain, 1);
	bool reachedUnreliable = false;
	for (const auto& candidate : generated.candidates)
	{
		if (!candidate.reliabilityPassed) reachedUnreliable = true;
		else BOOST_CHECK(!reachedUnreliable);
	}
}

BOOST_AUTO_TEST_CASE(product_candidate_dual_priority_precedes_dictionary_truncation)
{
	VectorXd mean = VectorXd::Zero(4);
	MatrixXd covariance = MatrixXd::Identity(4, 4) * 1e-6;
	MatrixXd productCross(1, 4);
	productCross << 1e-6, 1, 2, 3;
	const ZhangExactVector lowGainComplementaryPair = {1, 0, 0, 0};
	const ZhangExactMatrix complementary = {lowGainComplementaryPair};

	const auto unprioritised = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6,
		{}, 0, 1, 1);
	const auto prioritised = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6,
		{}, 0, 1, 1, complementary);

	BOOST_REQUIRE(unprioritised.valid);
	BOOST_REQUIRE(prioritised.valid);
	BOOST_REQUIRE_EQUAL(unprioritised.candidates.size(), 1U);
	BOOST_REQUIRE_EQUAL(prioritised.candidates.size(), 1U);
	BOOST_CHECK(unprioritised.candidates.front().row !=
		lowGainComplementaryPair);
	BOOST_CHECK(prioritised.candidates.front().row ==
		lowGainComplementaryPair);
	BOOST_CHECK_EQUAL(
		prioritised.candidates.front().dualGraphRankGain, 1);
	BOOST_CHECK_EQUAL(
		prioritised.dualSupportedReliableRowsBeforeTruncation, 1);
	BOOST_CHECK_EQUAL(
		prioritised.dualSupportedReliableRowsRetained, 1);
	BOOST_CHECK_EQUAL(prioritised.dualGraphRankGainBeforeTruncation, 1);
	BOOST_CHECK_EQUAL(prioritised.dualGraphRankGainRetained, 1);
	BOOST_CHECK_EQUAL(prioritised.dualMembershipBatchTargets, 10);
	BOOST_CHECK_EQUAL(prioritised.dualScalarDecompositionsAvoided, 9);
}

BOOST_AUTO_TEST_CASE(product_candidate_exact_quotient_gain_precedes_refixed_gain)
{
	VectorXd mean = VectorXd::Zero(2);
	MatrixXd covariance = MatrixXd::Identity(2, 2) * 1e-6;
	MatrixXd productCross(1, 2);
	productCross << 1, 0;
	const ZhangExactMatrix existingSignalLattice = {{1, 0}};

	const auto unprioritised = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6,
		{}, 0, 1, 1);
	const auto prioritised = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6,
		{}, 0, 1, 1, {}, existingSignalLattice);

	BOOST_REQUIRE(unprioritised.valid);
	BOOST_REQUIRE(prioritised.valid);
	BOOST_REQUIRE_EQUAL(prioritised.candidates.size(), 1U);
	BOOST_CHECK(zhangIntegerRowLatticeContains(
		existingSignalLattice, unprioritised.candidates.front().row).contained);
	BOOST_CHECK(!zhangIntegerRowLatticeContains(
		existingSignalLattice, prioritised.candidates.front().row).contained);
	BOOST_CHECK_EQUAL(
		prioritised.candidates.front().productLatticeRankGain, 1);
	BOOST_CHECK_EQUAL(
		prioritised.quotientGainingReliableRowsBeforeTruncation, 2);
	BOOST_CHECK_EQUAL(
		prioritised.quotientGainingReliableRowsRetained, 1);
	// Two individually novel rows remain available before truncation, but the
	// exact quotient has only one missing dimension.  Retaining either correct
	// generator therefore loses no quotient rank.
	BOOST_CHECK_EQUAL(prioritised.quotientRankGainBeforeTruncation, 1);
	BOOST_CHECK_EQUAL(prioritised.quotientRankGainRetained, 1);
	BOOST_CHECK_EQUAL(prioritised.quotientMembershipBatchTargets, 3);
	BOOST_CHECK_EQUAL(prioritised.quotientScalarDecompositionsAvoided, 2);
}

BOOST_AUTO_TEST_CASE(product_candidate_priority_completes_complementary_signal_lattice)
{
	ProductIntegerCandidate unsupportedPair;
	unsupportedPair.row = {0, 0, 1};
	unsupportedPair.reliabilityPassed = true;
	unsupportedPair.signalGraphRankGain = 1;
	unsupportedPair.incrementalProductGain = 100;
	unsupportedPair.variance = 1e-4;

	ProductIntegerCandidate denseContained;
	denseContained.row = {1, 1, 0};
	denseContained.reliabilityPassed = true;
	denseContained.signalGraphRankGain = 0;
	denseContained.incrementalProductGain = 1000;
	denseContained.variance = 1e-5;

	ProductIntegerCandidate supportedPair;
	supportedPair.row = {1, -1, 0};
	supportedPair.reliabilityPassed = true;
	supportedPair.signalGraphRankGain = 1;
	supportedPair.incrementalProductGain = 1;
	supportedPair.variance = 1e-3;

	ProductIntegerCandidate unreliableSupportedPair = supportedPair;
	unreliableSupportedPair.row = {1, 0, 0};
	unreliableSupportedPair.reliabilityPassed = false;
	unreliableSupportedPair.incrementalProductGain = 1e6;

	std::vector<ProductIntegerCandidate> candidates = {
		unsupportedPair, denseContained, supportedPair,
		unreliableSupportedPair};
	const ZhangExactMatrix complementarySignalRows = {
		{1, 0, 0}, {0, 1, 0}};
	int batchTargets = 0;
	int scalarDecompositionsAvoided = 0;
	const int supported = zhangPrioritiseProductCandidatesForDualLattice(
		candidates, complementarySignalRows, &batchTargets,
		&scalarDecompositionsAvoided);

	// Only a reliable named pair is allowed to receive dual-graph priority.
	// A dense contained integer row remains conditioning-only, and an
	// unreliable supported pair remains behind every reliable candidate.
	BOOST_CHECK_EQUAL(supported, 1);
	BOOST_REQUIRE_EQUAL(candidates.size(), 4U);
	BOOST_CHECK_EQUAL(candidates.front().row.at(0), 1);
	BOOST_CHECK_EQUAL(candidates.front().row.at(1), -1);
	BOOST_CHECK_EQUAL(candidates.front().dualGraphRankGain, 1);
	BOOST_CHECK_EQUAL(candidates.at(1).dualGraphRankGain, 0);
	BOOST_CHECK_EQUAL(candidates.back().reliabilityPassed, false);
	BOOST_CHECK_EQUAL(candidates.back().dualGraphRankGain, 0);
	BOOST_CHECK_EQUAL(batchTargets, 3);
	BOOST_CHECK_EQUAL(scalarDecompositionsAvoided, 2);
}

BOOST_AUTO_TEST_CASE(
	product_candidate_dual_priority_ignores_already_connected_edges)
{
	ProductIntegerCandidate oldReferenceEdge;
	oldReferenceEdge.row = {1, 0, 0};
	oldReferenceEdge.reliabilityPassed = true;
	oldReferenceEdge.productLatticeRankGain = 0;

	ProductIntegerCandidate oldCycleEdge = oldReferenceEdge;
	oldCycleEdge.row = {1, -1, 0};

	ProductIntegerCandidate connectingEdge = oldReferenceEdge;
	connectingEdge.row = {0, 0, 1};
	connectingEdge.productLatticeRankGain = 1;

	std::vector<ProductIntegerCandidate> candidates = {
		oldReferenceEdge, oldCycleEdge, connectingEdge};
	const ZhangExactMatrix complementarySignalRows = {
		{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	const ZhangExactMatrix existingSignalRows = {
		{1, 0, 0}, {0, 1, 0}};
	int existingDualRank = -1;
	std::vector<int> componentLabels;
	const int gaining = zhangPrioritiseProductCandidatesForDualLattice(
		candidates, complementarySignalRows, nullptr, nullptr,
		existingSignalRows, &existingDualRank, &componentLabels);

	BOOST_CHECK_EQUAL(existingDualRank, 2);
	BOOST_CHECK_EQUAL(gaining, 1);
	BOOST_REQUIRE_EQUAL(componentLabels.size(), 4U);
	BOOST_CHECK_EQUAL(componentLabels[0], componentLabels[1]);
	BOOST_CHECK_EQUAL(componentLabels[0], componentLabels[3]);
	BOOST_CHECK_NE(componentLabels[0], componentLabels[2]);
	BOOST_REQUIRE_EQUAL(candidates.size(), 3U);
	BOOST_CHECK(candidates.front().row == ZhangExactVector({0, 0, 1}));
	BOOST_CHECK_EQUAL(candidates.front().dualGraphRankGain, 1);
	BOOST_CHECK_EQUAL(candidates.at(1).dualGraphRankGain, 0);
	BOOST_CHECK_EQUAL(candidates.at(2).dualGraphRankGain, 0);
}

BOOST_AUTO_TEST_CASE(
	product_final_selection_scores_incremental_completion_not_absolute_forest)
{
	// Nodes 0, 1 and the implicit reference node 3 are already connected in
	// both frequencies.  Re-fixing 0-1 has a non-empty standalone pair forest
	// but contributes no persistent completion rank.  Fixing node 2 to the
	// reference contributes one dual-graph rank and two joint WL/L1 ranks.
	const ZhangExactMatrix existingRows = {
		{1, 0, 0},
		{0, 1, 0}};
	const ZhangExactVector existingValues = {1, 2};
	const ZhangExactMatrix redundantRows = {{1, -1, 0}};
	const ZhangExactVector redundantValues = {-1};
	const ZhangExactMatrix completingRows = {{0, 0, 1}};
	const ZhangExactVector completingValues = {3};

	const auto redundant = zhangScoreProductConstraintCompletion(
		existingRows, existingValues,
		existingRows, existingValues,
		redundantRows, redundantValues,
		redundantRows, redundantValues,
		3);
	const auto completing = zhangScoreProductConstraintCompletion(
		existingRows, existingValues,
		existingRows, existingValues,
		completingRows, completingValues,
		completingRows, completingValues,
		3);
	const auto baseline = zhangScoreProductConstraintCompletion(
		existingRows, existingValues,
		existingRows, existingValues,
		{}, {}, {}, {}, 3);
	const auto completingCached = zhangScoreProductConstraintCompletion(
		existingRows, existingValues,
		existingRows, existingValues,
		completingRows, completingValues,
		completingRows, completingValues,
		3, baseline.existingDualGraphRank);

	BOOST_REQUIRE_MESSAGE(redundant.valid, redundant.failureReason);
	BOOST_REQUIRE_MESSAGE(completing.valid, completing.failureReason);
	BOOST_REQUIRE_MESSAGE(baseline.valid, baseline.failureReason);
	BOOST_REQUIRE_MESSAGE(completingCached.valid,
		completingCached.failureReason);
	BOOST_CHECK_EQUAL(baseline.dualGraphRankGain, 0);
	BOOST_CHECK_EQUAL(completingCached.dualGraphRankGain,
		completing.dualGraphRankGain);
	BOOST_CHECK_EQUAL(completingCached.quotientRankGain,
		completing.quotientRankGain);
	BOOST_CHECK_EQUAL(redundant.existingDualGraphRank, 2);
	BOOST_CHECK_EQUAL(redundant.completedDualGraphRank, 2);
	BOOST_CHECK_EQUAL(redundant.dualGraphRankGain, 0);
	BOOST_CHECK_EQUAL(redundant.quotientRankGain, 0);
	BOOST_CHECK_EQUAL(completing.existingDualGraphRank, 2);
	BOOST_CHECK_EQUAL(completing.completedDualGraphRank, 3);
	BOOST_CHECK_EQUAL(completing.dualGraphRankGain, 1);
	BOOST_CHECK_EQUAL(completing.existingJointLatticeRank, 4);
	BOOST_CHECK_EQUAL(completing.completedJointLatticeRank, 6);
	BOOST_CHECK_EQUAL(completing.quotientRankGain, 2);
	BOOST_CHECK(zhangProductCompletionScoreBetter(completing, redundant));
	BOOST_CHECK(!zhangProductCompletionScoreBetter(redundant, completing));

	const auto conflict = zhangScoreProductConstraintCompletion(
		existingRows, existingValues,
		existingRows, existingValues,
		ZhangExactMatrix{{1, 0, 0}}, ZhangExactVector{9},
		{}, {}, 3);
	BOOST_CHECK(!conflict.valid);
	BOOST_CHECK_EQUAL(conflict.failureReason,
		"PRODUCT_COMPLETION_SCORE_AFFINE_UNION_FAILED");
}

BOOST_AUTO_TEST_CASE(product_candidate_rejects_mismatched_priority_lattice_dimensions)
{
	VectorXd mean = VectorXd::Zero(2);
	MatrixXd covariance = MatrixXd::Identity(2, 2) * 1e-6;
	MatrixXd productCross = MatrixXd::Zero(1, 2);
	const ZhangExactMatrix wrongDimension = {{1, 0, 0}};

	const auto complementaryMismatch = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6,
		{}, 0, 1, 8, wrongDimension);
	BOOST_CHECK(!complementaryMismatch.valid);
	BOOST_CHECK_EQUAL(complementaryMismatch.failureReason,
		"PRODUCT_CANDIDATE_COMPLEMENTARY_LATTICE_DIMENSION_MISMATCH");

	const auto existingMismatch = generateProductIntegerCandidates(
		mean, covariance, productCross, 1e-3, 1e-6,
		{}, 0, 1, 8, {}, wrongDimension);
	BOOST_CHECK(!existingMismatch.valid);
	BOOST_CHECK_EQUAL(existingMismatch.failureReason,
		"PRODUCT_CANDIDATE_EXISTING_LATTICE_DIMENSION_MISMATCH");
}

BOOST_AUTO_TEST_CASE(product_ledger_joint_support_splits_wl_l1_and_mixed_rows)
{
	const ZhangExactMatrix jointRows = {
		{1, 0, -1, 0},
		{0, 1, 0, 0},
		{1, 1, 0, 1}};
	const ZhangExactVector jointValues = {3, 4, 5};
	const auto support = zhangSplitJointProductIntegerSupport(
		jointRows, jointValues, 2);

	BOOST_REQUIRE(support.valid);
	BOOST_CHECK_EQUAL(support.failureReason, "NONE");
	BOOST_CHECK_EQUAL(support.inputRows, 3);
	BOOST_CHECK_EQUAL(support.invalidRows, 0);
	BOOST_CHECK_EQUAL(support.mixedConditioningOnlyRows, 1);
	BOOST_REQUIRE_EQUAL(support.wideLaneRows.size(), 1U);
	BOOST_REQUIRE_EQUAL(support.firstSignalRows.size(), 1U);
	BOOST_CHECK(support.wideLaneRows.front() == ZhangExactVector({1, 0}));
	BOOST_CHECK(support.firstSignalRows.front() == ZhangExactVector({0, 1}));
	BOOST_CHECK_EQUAL(support.wideLaneValues.front(), 3);
	BOOST_CHECK_EQUAL(support.firstSignalValues.front(), 4);
}

BOOST_AUTO_TEST_CASE(product_ledger_joint_support_rejects_affine_conflict)
{
	const ZhangExactMatrix jointRows = {
		{1, 0, -1, 0},
		{1, 0, -1, 0}};
	const ZhangExactVector jointValues = {3, 4};
	const auto support = zhangSplitJointProductIntegerSupport(
		jointRows, jointValues, 2);

	BOOST_CHECK(!support.valid);
	BOOST_CHECK_EQUAL(
		support.failureReason, "JOINT_PRODUCT_SUPPORT_AFFINE_CONFLICT");
}

BOOST_AUTO_TEST_CASE(product_ledger_joint_support_recovers_pure_integer_combination)
{
	// Neither [0,-1] nor the current HNF basis row replacing it is itself the
	// required WL relation.  The exact lattice nevertheless contains
	// [1,-1] = [1,0] + [0,-1], with integer value 2.
	const ZhangExactMatrix jointRows = {{1, 0}, {0, -1}};
	const ZhangExactVector jointValues = {3, -1};
	const auto support = zhangSplitJointProductIntegerSupport(
		jointRows, jointValues, 1);

	BOOST_REQUIRE(support.valid);
	BOOST_REQUIRE_EQUAL(support.wideLaneRows.size(), 1U);
	BOOST_REQUIRE_EQUAL(support.firstSignalRows.size(), 1U);
	BOOST_CHECK_EQUAL(support.wideLaneRows.front().front(), 1);
	BOOST_CHECK_EQUAL(support.wideLaneValues.front(), 2);
	BOOST_CHECK_EQUAL(support.firstSignalRows.front().front(), 1);
	BOOST_CHECK_EQUAL(support.firstSignalValues.front(), 3);
	BOOST_CHECK_EQUAL(support.mixedConditioningOnlyRows, 0);
}

BOOST_AUTO_TEST_CASE(joint_wl_l1_rows_pull_back_exactly_to_l1_l2)
{
	// r_y=[2,-1 | 3,4] in y=[WL,L1].  Since WL=L1-L2,
	// r_x=[2+3,-1+4 | -2,1]=[5,3|-2,1].
	bool dimensionsValid = false;
	const auto rows = zhangWideLaneFirstRowsToSignalRows(
		{{2, -1, 3, 4}}, 2, &dimensionsValid);
	BOOST_REQUIRE(dimensionsValid);
	BOOST_REQUIRE_EQUAL(rows.size(), 1U);
	BOOST_CHECK(rows.front() == ZhangExactVector({5, 3, -2, 1}));

	dimensionsValid = true;
	const auto rejected = zhangWideLaneFirstRowsToSignalRows(
		{{1, 0, 0}}, 2, &dimensionsValid);
	BOOST_CHECK(rejected.empty());
	BOOST_CHECK(!dimensionsValid);
}

BOOST_AUTO_TEST_CASE(
	persistent_joint_lattice_builds_exact_affine_stochastic_quotient)
{
	const ZhangExactMatrix deterministicRows = {
		{1, -1, 0, 0},
		{0, 0, 1, -1}};
	const ZhangExactVector deterministicValues = {3, -2};
	const auto quotient = zhangExactAffineIntegerQuotient(
		deterministicRows, deterministicValues, 4);
	BOOST_REQUIRE_MESSAGE(quotient.valid, quotient.failureReason);
	BOOST_CHECK_EQUAL(quotient.deterministicRank, 2);
	BOOST_CHECK_EQUAL(quotient.quotientRank, 2);
	BOOST_REQUIRE_EQUAL(quotient.smithInvariants.size(), 2);
	for (const auto& invariant : quotient.smithInvariants)
		BOOST_CHECK(zhangExactAbs(invariant) == 1);
	BOOST_CHECK(zhangExactMatrixTimesColumn(
		quotient.deterministicBasis, quotient.particularSolution) ==
		quotient.deterministicValues);

	ZhangExactMatrix kernelTranspose(4, ZhangExactVector(2));
	for (int row = 0; row < 2; row++)
	for (int column = 0; column < 4; column++)
		kernelTranspose[column][row] = quotient.kernelBasis[row][column];
	BOOST_CHECK(zhangExactMultiply(
		quotient.deterministicBasis, kernelTranspose) ==
		zhangExactZeroMatrix(2, 2));
	BOOST_CHECK(zhangExactMultiply(
		quotient.quotientProjector, kernelTranspose) ==
		zhangExactIdentityMatrix(2));

	const ZhangExactVector quotientInteger = {4, -5};
	ZhangExactVector reconstructed = quotient.particularSolution;
	const ZhangExactVector displacement = zhangExactMatrixTimesColumn(
		kernelTranspose, quotientInteger);
	for (int coordinate = 0; coordinate < 4; coordinate++)
		reconstructed[coordinate] += displacement[coordinate];
	BOOST_CHECK(zhangExactMatrixTimesColumn(
		quotient.deterministicBasis, reconstructed) ==
		quotient.deterministicValues);
	ZhangExactVector relative = reconstructed;
	for (int coordinate = 0; coordinate < 4; coordinate++)
		relative[coordinate] -= quotient.particularSolution[coordinate];
	BOOST_CHECK(zhangExactMatrixTimesColumn(
		quotient.quotientProjector, relative) == quotientInteger);
}

BOOST_AUTO_TEST_CASE(persistent_joint_lattice_enforces_nonprimitive_divisibility)
{
    const auto quotient = zhangExactAffineIntegerQuotient({{2, 0}}, {2}, 2);
    BOOST_REQUIRE(quotient.valid);
    BOOST_CHECK_EQUAL(quotient.particularSolution[0],1);
    BOOST_CHECK(!zhangExactAffineIntegerQuotient({{2,0}},{1},2).valid);
}

BOOST_AUTO_TEST_CASE(
	persistent_joint_44d_rank24_exposes_full_rank20_stochastic_quotient)
{
	constexpr int ambient = 44;
	constexpr int deterministicRank = 24;
	ZhangExactMatrix rows;
	ZhangExactVector values;
	for (int row = 0; row < deterministicRank; row++)
	{
		ZhangExactVector unit(ambient);
		unit[row] = 1;
		rows.push_back(std::move(unit));
		values.push_back(row - 7);
	}
	const auto quotient = zhangExactAffineIntegerQuotient(rows, values, ambient);
	BOOST_REQUIRE_MESSAGE(quotient.valid, quotient.failureReason);
	BOOST_CHECK_EQUAL(quotient.deterministicRank, deterministicRank);
	BOOST_CHECK_EQUAL(quotient.quotientRank, 20);

	MatrixXd projector(quotient.quotientRank, ambient);
	for (int row = 0; row < projector.rows(); row++)
		projector.row(row) = zhangExactRowToDouble(
			quotient.quotientProjector[row]).transpose();
	MatrixXd covariance = MatrixXd::Zero(ambient, ambient);
	for (int coordinate = deterministicRank; coordinate < ambient; coordinate++)
		covariance(coordinate, coordinate) = 0.5 + coordinate;
	const MatrixXd quotientCovariance = projector * covariance *
		projector.transpose();
	const auto effective = zhangPosteriorEffectiveIntegerRowBasis(
		quotientCovariance, MatrixXd::Identity(20, 20));
	BOOST_REQUIRE(effective.valid);
	BOOST_CHECK_EQUAL(effective.posteriorRank, 20);
}

BOOST_AUTO_TEST_CASE(no_new_joint_claim_adds_zero_incremental_risk)
{
	BOOST_CHECK_CLOSE(zhangPersistentIncrementalFailureProbability(
		4e-4, false, 1.0), 4e-4, 1e-9);
	BOOST_CHECK_CLOSE(zhangPersistentIncrementalFailureProbability(
		4e-4, true, 2e-4), 6e-4, 1e-9);
	BOOST_CHECK_EQUAL(zhangPersistentIncrementalFailureProbability(
		4e-4, true, std::numeric_limits<double>::quiet_NaN()), 1.0);
}

BOOST_AUTO_TEST_CASE(incremental_risk_budget_never_injects_phantom_floor)
{
	const auto infeasible = zhangAllocateIncrementalRiskBudget(0.0012, 0.001);
	BOOST_CHECK(!infeasible.parentFeasible);
	BOOST_CHECK(!infeasible.searchAuthorized);
	BOOST_CHECK_EQUAL(infeasible.remainingBudget, 0.0);
	BOOST_CHECK(infeasible.status ==
		ZhangIncrementalRiskBudgetStatus::PARENT_INFEASIBLE);

	const auto exhausted = zhangAllocateIncrementalRiskBudget(0.001, 0.001);
	BOOST_CHECK(exhausted.parentFeasible);
	BOOST_CHECK(!exhausted.searchAuthorized);
	BOOST_CHECK_EQUAL(exhausted.remainingBudget, 0.0);
	BOOST_CHECK(exhausted.status ==
		ZhangIncrementalRiskBudgetStatus::BUDGET_EXHAUSTED);

	const auto positive = zhangAllocateIncrementalRiskBudget(0.0004, 0.001);
	BOOST_CHECK(positive.parentFeasible);
	BOOST_CHECK(positive.searchAuthorized);
	BOOST_CHECK_CLOSE(positive.remainingBudget, 0.0006, 1e-9);
	BOOST_CHECK(positive.status ==
		ZhangIncrementalRiskBudgetStatus::SEARCH_WITH_POSITIVE_BUDGET);
}

BOOST_AUTO_TEST_CASE(conditioned_posterior_and_receipt_are_atomic)
{
	const auto rollback = zhangSelectPosteriorReceiptTransaction(true, false);
	BOOST_CHECK(rollback.rollbackRequired);
	BOOST_CHECK(!rollback.useConditionedPosterior);

	const auto committed = zhangSelectPosteriorReceiptTransaction(true, true);
	BOOST_CHECK(!committed.rollbackRequired);
	BOOST_CHECK(committed.useConditionedPosterior);

	const auto untouched = zhangSelectPosteriorReceiptTransaction(false, false);
	BOOST_CHECK(!untouched.rollbackRequired);
	BOOST_CHECK(!untouched.useConditionedPosterior);
}

BOOST_AUTO_TEST_CASE(optional_ledger_cannot_poison_feasible_baseline)
{
	const auto rejected = zhangSelectOptionalLedgerWithoutPoisoningBaseline(
		true, 0.00095706, true, 0.696969, 0.001);
	BOOST_CHECK(rejected.baselineFeasible);
	BOOST_CHECK(rejected.baselinePreserved);
	BOOST_CHECK(!rejected.optionalAuthorized);
	BOOST_CHECK_CLOSE(rejected.selectedFailureProbability, 0.00095706, 1e-9);
	BOOST_CHECK_EQUAL(rejected.status,
		"BASELINE_PRESERVED_OPTIONAL_EXCEEDS_RISK_BUDGET");

	const auto accepted = zhangSelectOptionalLedgerWithoutPoisoningBaseline(
		true, 4e-4, true, 2e-4, 1e-3);
	BOOST_CHECK(accepted.baselinePreserved);
	BOOST_CHECK(accepted.optionalAuthorized);
	BOOST_CHECK_CLOSE(accepted.proposedUnionBound, 6e-4, 1e-9);
	BOOST_CHECK_CLOSE(accepted.selectedFailureProbability, 6e-4, 1e-9);

	const auto malformed = zhangSelectOptionalLedgerWithoutPoisoningBaseline(
		true, 4e-4, true,
		std::numeric_limits<double>::quiet_NaN(), 1e-3);
	BOOST_CHECK(malformed.baselinePreserved);
	BOOST_CHECK(!malformed.optionalAuthorized);
	BOOST_CHECK_CLOSE(malformed.selectedFailureProbability, 4e-4, 1e-9);

	const auto noOptional = zhangSelectOptionalLedgerWithoutPoisoningBaseline(
		true, 4e-4, false, 1.0, 1e-3);
	BOOST_CHECK(noOptional.baselinePreserved);
	BOOST_CHECK(!noOptional.optionalAuthorized);
	BOOST_CHECK_CLOSE(noOptional.selectedFailureProbability, 4e-4, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	persistent_certificate_snapshot_distinguishes_external_posterior_nullity)
{
	const ZhangExactMatrix rows{{1, 0, 0}, {0, 1, 0}};
	const ZhangExactVector values{3, -2};
	const auto snapshot = zhangBuildPersistentCertificateSnapshot(
		rows, values, 3, 3, {"FAMILY-B", "FAMILY-A", "FAMILY-A"},
		4e-4, 8e-4, "POSTERIOR-100|HASH-A");
	BOOST_REQUIRE_MESSAGE(snapshot.valid, snapshot.failureReason);
	BOOST_CHECK_EQUAL(snapshot.deterministicRank, 2);
	BOOST_CHECK_EQUAL(snapshot.expectedPosteriorNullity, 3);
	BOOST_CHECK_EQUAL(snapshot.externalPosteriorNullity, 1);
	BOOST_CHECK_EQUAL(snapshot.certificateFamilyIds.size(), 2);
	BOOST_CHECK_EQUAL(snapshot.posteriorMomentId, "POSTERIOR-100|HASH-A");
	BOOST_CHECK(!snapshot.snapshotId.empty());

	const auto divergence = zhangBuildPersistentCertificateSnapshot(
		rows, values, 3, 1, {"FAMILY-A"}, 4e-4, 8e-4,
		"POSTERIOR-100|HASH-A");
	BOOST_CHECK(!divergence.valid);
	BOOST_CHECK_EQUAL(divergence.failureReason,
		"INTERNAL_PERSISTENT_SNAPSHOT_DIVERGENCE");
}

BOOST_AUTO_TEST_CASE(
	certificate_risk_uses_current_scope_and_charges_family_once)
{
	const auto audit = zhangAuditCertificateFamilyRisk({
		{"FAMILY-A", 2e-4, 6e-4, true},
		{"FAMILY-A", 3e-4, 8e-4, true},
		{"FAMILY-B", 4e-4, 4e-4, false}});
	BOOST_REQUIRE_MESSAGE(audit.valid, audit.failureReason);
	BOOST_CHECK_EQUAL(audit.currentFamilies, 1);
	BOOST_CHECK_EQUAL(audit.lifetimeFamilies, 2);
	BOOST_CHECK_CLOSE(audit.currentCertificateRisk, 3e-4, 1e-9);
	BOOST_CHECK_CLOSE(audit.lifetimeBroadcastRisk, 1.2e-3, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	product_gauge_family_selection_retains_feasible_fallback)
{
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>{4e-4, 4e-4, 4e-4}, 8e-4,
		[](const std::vector<int>& families)
		{
			ZhangProductGaugeFamilySubsetScore score;
			// Family 2 conflicts with every other family.  The established family
			// zero remains a safe fallback while family one can extend it.
			const bool hasTwo = std::find(families.begin(), families.end(), 2) !=
				families.end();
			score.valid = !hasTwo;
			score.dualGraphRank = score.valid ? families.size() : 0;
			score.certifiedSatelliteCount = score.valid
				? score.dualGraphRank + 1 : 0;
			score.largestCertifiedComponentSize =
				score.certifiedSatelliteCount;
			return score;
		}, 8, 4, {0});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(selected.fallbackAvailable);
	BOOST_CHECK(selected.fallbackPreserved);
	BOOST_CHECK((selected.selectedFamilyIndices == std::vector<int>{0, 1}));
}

BOOST_AUTO_TEST_CASE(r26_applied_lattice_image_recovers_combinations_and_affine_offset)
{
	// Neither physical row is itself a target row; their difference cancels r.
	const auto image = zhangAppliedLatticeProductImage(
		{{1, 0, 1}, {0, 1, 1}}, {7, 4}, {{1, -1, 0}}, {2});
	BOOST_REQUIRE(image.consistent);
	BOOST_REQUIRE_EQUAL(image.basis.size(), 1);
	BOOST_CHECK(image.basis[0] == ZhangExactVector{1});
	BOOST_CHECK(image.values[0] == 5);
	// A nonprimitive parent must not silently authorize a unit integer.
	const auto even = zhangAppliedLatticeProductImage({{2, 0}}, {6}, {{1, 0}}, {0});
	BOOST_REQUIRE(even.consistent);
	BOOST_CHECK(even.basis[0] == ZhangExactVector{2});
	BOOST_CHECK(!zhangIntegerRowLatticeContains(even.basis, {1}).contained);
	BOOST_CHECK(!zhangAppliedLatticeProductImage({{1}}, {1}, {{1, 0}}, {0}).consistent);
}

BOOST_AUTO_TEST_CASE(r26_incremental_nis_matches_reference_including_singular_trials)
{
	std::mt19937 gen(260);
	std::normal_distribution<double> normal;
	for (int repeat=0; repeat<12; ++repeat)
	{
		MatrixXd l(24,24);
		for (int i=0;i<l.size();++i) l.data()[i]=normal(gen)*0.1;
		MatrixXd p=l*l.transpose()+MatrixXd::Identity(24,24);
		VectorXd v(24); for(int i=0;i<24;++i) v(i)=normal(gen);
		if(repeat%2) { p.row(23)=p.row(0); p.col(23)=p.col(0); v(23)=v(0)+(repeat%3==0?0.1:0); }
		ZhangIncrementalSubsetNis fast(p,v), reference(p,v);
		for(int i=0;i<24;++i)
		{
			auto a=fast.trial(i), b=reference.trial(i,true);
			BOOST_CHECK_EQUAL(a.valid,b.valid);
			BOOST_CHECK_EQUAL(a.rank,b.rank);
			if(a.valid && b.valid)
			{
				BOOST_CHECK_SMALL(a.nis-b.nis,1e-8);
				if (b.nis < 2*(i+1)) { fast.commit(i,std::move(a)); reference.commit(i,std::move(b)); }
			}
		}
		BOOST_CHECK_GT(fast.fastTrials,0);
	}
}

BOOST_AUTO_TEST_CASE(r27_missing_search_column_keeps_exact_held_consequence)
{
	// x+z=7, y+z=4: neither row survives individually, x-y=3 does.
	const auto kept = zhangExactSurvivingLattice({{1,0,1},{0,1,1}}, {7,4}, {true,true,false});
	BOOST_REQUIRE(kept.consistent);
	BOOST_REQUIRE_EQUAL(kept.basis.size(),1);
	BOOST_CHECK(kept.basis[0] == (ZhangExactVector{1,-1}));
	BOOST_CHECK(kept.values[0] == 3);
	// A nonprimitive consequence stays nonprimitive; no saturation.
	const auto even = zhangExactSurvivingLattice({{2,0,1},{0,2,1}}, {7,3}, {true,true,false});
	BOOST_REQUIRE(even.consistent);
	BOOST_CHECK(even.basis[0] == (ZhangExactVector{2,-2}));
	BOOST_CHECK(even.values[0] == 4);
	BOOST_CHECK(!zhangExactSurvivingLattice({{1,0,1},{1,0,1}}, {7,8}, {true,true,false}).consistent);
}

BOOST_AUTO_TEST_CASE(r51_current_domain_projection_retains_integer_congruence)
{
	// x+z=7, y+z=4 projects to x-y=3 after z retires. A current
	// coordinate that did not exist in the old receipt remains unconstrained.
	const auto projected=zhangExactCurrentDomain(
		{{1,0,1},{0,1,1}}, {7,4}, {0,2,-1}, 4);
	BOOST_REQUIRE(projected.valid);
	BOOST_CHECK_EQUAL(projected.removedColumns,1);
	BOOST_CHECK_EQUAL(projected.touchedRows,2);
	BOOST_REQUIRE_EQUAL(projected.rows.size(),1);
	BOOST_CHECK(projected.rows[0]==(ZhangExactVector{1,0,-1,0}));
	BOOST_CHECK(projected.values[0]==3);

	const auto even=zhangExactCurrentDomain(
		{{2,0,1},{0,2,1}}, {7,3}, {0,2,-1}, 4);
	BOOST_REQUIRE(even.valid);
	BOOST_CHECK(even.rows[0]==(ZhangExactVector{2,0,-2,0}));
	BOOST_CHECK(even.values[0]==4);

	const auto collision=zhangExactCurrentDomain(
		{{1,0},{0,1}}, {1,2}, {0,0}, 1);
	BOOST_CHECK(!collision.valid);
	BOOST_CHECK_EQUAL(collision.failureReason,"CURRENT_COLUMN_COLLISION_OR_RANGE");
}

BOOST_AUTO_TEST_CASE(r27_complete_domain_conditioning_precedes_marginalization)
{
	VectorXd mean=VectorXd::Zero(3), rhs(2); rhs << 7,4;
	MatrixXd covariance(3,3); covariance << 3,0.2,0.8, 0.2,2,0.4, 0.8,0.4,4;
	MatrixXd rows(2,3); rows << 1,0,1, 0,1,1;
	const auto full=zhangConditionPosteriorEffectiveIntegers(mean,covariance,rows,rhs);
	BOOST_REQUIRE(full.valid);
	BOOST_CHECK_SMALL(full.mean(0)-full.mean(1)-3,1e-8);
	// Conditioning only the surviving difference is NOT the same marginal:
	// nuisance correlations still contain information and cannot be discarded.
	MatrixXd difference(1,2); difference << 1,-1;
	VectorXd value(1); value << 3;
	const auto naive=zhangConditionPosteriorEffectiveIntegers(
		VectorXd(mean.head(2)),MatrixXd(covariance.topLeftCorner(2,2)),difference,value);
	BOOST_REQUIRE(naive.valid);
	BOOST_CHECK_GT((full.covariance.topLeftCorner(2,2)-naive.covariance).norm(),0.1);
	const auto repeated=zhangConditionPosteriorEffectiveIntegers(full.mean,full.covariance,rows,rhs);
	BOOST_REQUIRE(repeated.valid);
	BOOST_CHECK_SMALL((repeated.mean-full.mean).norm(),1e-7);
	BOOST_CHECK_SMALL((repeated.covariance-full.covariance).norm(),1e-7);
}

BOOST_AUTO_TEST_CASE(r27_whole_product_image_recovers_combinations_not_rows)
{
	const ZhangExactMatrix raw{{1,0,1},{0,1,1}};
	const ZhangExactMatrix target{{1,-1,0}};
	const auto membership=zhangIntegerRowLatticeContainsBatch(target,raw);
	BOOST_REQUIRE_EQUAL(membership.size(),2);
	BOOST_CHECK(!membership[0].contained && !membership[1].contained);
	const auto image=zhangAppliedLatticeProductImage(raw,{7,4},target,{5});
	BOOST_REQUIRE(image.consistent);
	BOOST_REQUIRE_EQUAL(image.basis.size(),1);
	BOOST_CHECK(image.basis[0] == ZhangExactVector{1});
	BOOST_CHECK(image.values[0] == 8); // (x-y)+5
}

BOOST_AUTO_TEST_CASE(r27_final_baseline_must_pass_same_posterior_not_old_one)
{
	VectorXd joint(2); joint << 0.2,20;
	const auto rejected=assessZhangIntegerCandidateNis(joint,MatrixXd::Identity(2,2),1e-6);
	const auto baseline=assessZhangIntegerCandidateNis(VectorXd(joint.head(1)),MatrixXd::Identity(1,1),1e-6);
	BOOST_REQUIRE(rejected.valid && baseline.valid);
	BOOST_CHECK_GT(rejected.nis,rejected.threshold);
	BOOST_CHECK_LT(baseline.nis,baseline.threshold);
	VectorXd changed(1); changed << 20;
	const auto stale=assessZhangIntegerCandidateNis(changed,MatrixXd::Identity(1,1),1e-6);
	BOOST_CHECK_GT(stale.nis,stale.threshold);
	const auto deterministic=assessZhangIntegerCandidateNis(changed,MatrixXd::Zero(1,1),1e-6);
	BOOST_CHECK(!deterministic.valid);
}

BOOST_AUTO_TEST_CASE(r26_workspace_matches_uncached_psd_and_spd_decisions)
{
	std::mt19937 gen(26);
	std::normal_distribution<double> normal;
	for (int rank : {0, 3, 8})
	{
		MatrixXd factor(8, rank), rows(12, 8);
		for (int i=0; i<factor.size(); ++i) factor.data()[i] = normal(gen);
		for (int i=0; i<rows.size(); ++i) rows.data()[i] = std::round(normal(gen));
		rows.row(11) = rows.row(0) + rows.row(1);
		MatrixXd p = factor * factor.transpose();
		const auto workspace = zhangBuildPosteriorEffectiveWorkspace(p);
		BOOST_REQUIRE(workspace.valid);
		for (int count=1; count<=12; ++count)
		{
			const auto plain = zhangPosteriorEffectiveIntegerRowBasis(p, MatrixXd(rows.topRows(count)));
			const auto cached = zhangPosteriorEffectiveIntegerRowBasis(workspace, MatrixXd(rows.topRows(count)));
			BOOST_CHECK(plain.valid == cached.valid);
			BOOST_CHECK(plain.retainedRows == cached.retainedRows);
			BOOST_CHECK(plain.deterministicRows == cached.deterministicRows);
			BOOST_CHECK_EQUAL(plain.varianceTolerance, cached.varianceTolerance);
		}
	}
	BOOST_CHECK(!zhangBuildPosteriorEffectiveWorkspace(-MatrixXd::Identity(3,3)).valid);
}

BOOST_AUTO_TEST_CASE(r28_snapshot_audit_reuses_ambient_posterior_nullspace)
{
	// The ambient decomposition removes the 17 numerical-floor directions at
	// about 6.5e-13.  A second DPD' decomposition with the legacy 1e-14 floor
	// would incorrectly classify them as stochastic.
	MatrixXd covariance = MatrixXd::Zero(46, 46);
	covariance.diagonal().head(29).setOnes();
	covariance.diagonal().tail(17).setConstant(1e-13);
	const auto workspace = zhangBuildPosteriorEffectiveWorkspace(covariance);
	BOOST_REQUIRE(workspace.valid);
	BOOST_CHECK_EQUAL(workspace.covarianceRank, 29);
	BOOST_CHECK_GT(workspace.eigenTolerance, 1e-13);
	MatrixXd rows = MatrixXd::Zero(17, 46);
	rows.rightCols(17).setIdentity();
	VectorXd innovation = VectorXd::Zero(17);
	const auto shared = zhangAuditPosteriorDeterministicCertificate(
		workspace, rows, innovation);
	BOOST_REQUIRE(shared.valid);
	BOOST_CHECK(shared.deterministic);
	BOOST_CHECK(shared.affineConsistent);
	BOOST_CHECK_EQUAL(shared.projectedStochasticRank, 0);
	BOOST_CHECK_EQUAL(shared.failureReason, "NONE");

	const auto legacy = assessZhangIntegerCandidateNis(
		innovation, rows * covariance * rows.transpose(), 1e-6);
	BOOST_REQUIRE(legacy.valid);
	BOOST_CHECK(!legacy.deterministic);
	BOOST_CHECK_EQUAL(legacy.rank, 17);
}

BOOST_AUTO_TEST_CASE(r28_snapshot_audit_fails_closed_by_specific_cause)
{
	MatrixXd covariance = MatrixXd::Zero(3, 3);
	covariance(0, 0) = 1;
	const auto workspace = zhangBuildPosteriorEffectiveWorkspace(covariance);
	BOOST_REQUIRE(workspace.valid);

	MatrixXd stochasticRow(1, 3); stochasticRow << 1, 0, 0;
	const auto stochastic = zhangAuditPosteriorDeterministicCertificate(
		workspace, stochasticRow, VectorXd::Zero(1));
	BOOST_REQUIRE(stochastic.valid);
	BOOST_CHECK(!stochastic.deterministic);
	BOOST_CHECK_EQUAL(stochastic.projectedStochasticRank, 1);
	BOOST_CHECK_EQUAL(stochastic.failureReason,
		"SNAPSHOT_HAS_STOCHASTIC_COMPONENT");

	MatrixXd deterministicRow(1, 3); deterministicRow << 0, 1, 0;
	VectorXd inconsistent(1); inconsistent << 2e-7;
	const auto affine = zhangAuditPosteriorDeterministicCertificate(
		workspace, deterministicRow, inconsistent);
	BOOST_REQUIRE(affine.valid);
	BOOST_CHECK(affine.deterministic);
	BOOST_CHECK(!affine.affineConsistent);
	BOOST_CHECK_EQUAL(affine.failureReason, "SNAPSHOT_AFFINE_RESIDUAL");
}

BOOST_AUTO_TEST_CASE(r26_historical_subblocks_keep_parent_risk_and_reject_bad_block)
{
	const auto selected = zhangSelectProductGaugeAdmissionFamilies(
		{6e-4,6e-4,6e-4,5e-4}, 1e-3, [](const std::vector<int>& rows)
		{
			ZhangProductGaugeFamilySubsetScore s;
			s.valid = !rows.empty() && std::find(rows.begin(),rows.end(),1)==rows.end();
			s.dualGraphRank = rows.size(); s.largestCertifiedComponentSize = rows.size()+1;
			return s;
		}, 8, 64, {0}, {"parent","parent","parent","different"});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK(selected.selectedFamilyIndices == (std::vector<int>{0,2}));
	BOOST_CHECK_CLOSE(selected.failureProbability,6e-4,1e-9);
	BOOST_CHECK(selected.fallbackPreserved);
}

BOOST_AUTO_TEST_CASE(r26_exact_graph_rejects_rank_plus_one_island_estimate)
{
	// Four named nodes plus reference; edges (0,1) and (2,reference).
	const ZhangExactMatrix rows{{1,-1,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0},
		{0,0,0,0,1,-1,0,0}, {0,0,0,0,0,0,1,0}};
	const auto graph = zhangExactDualGraphScore(rows, {0,0,0,0}, 4);
	BOOST_REQUIRE(graph.valid);
	BOOST_CHECK_EQUAL(graph.graphRank,2);
	BOOST_CHECK_EQUAL(graph.satellites,4);
	BOOST_CHECK_EQUAL(graph.components,2);
	BOOST_CHECK_EQUAL(graph.largest,2);
}

BOOST_AUTO_TEST_CASE(r26_larger_component_cannot_destroy_safe_baseline_rank)
{
	const auto choice=zhangSelectProductGaugeAdmissionFamilies({0.0004,0.0004},0.001,
		[](const std::vector<int>& ids)
		{
			ZhangProductGaugeFamilySubsetScore s;
			s.valid=ids.size()==1;
			if(ids==std::vector<int>{0}) { s.dualGraphRank=6; s.certifiedSatelliteCount=9; s.largestCertifiedComponentSize=4; }
			else { s.dualGraphRank=4; s.certifiedSatelliteCount=5; s.largestCertifiedComponentSize=5; }
			return s;
		},8,64,{0});
	BOOST_REQUIRE(choice.valid);
	BOOST_CHECK(choice.selectedFamilyIndices==std::vector<int>{0});
}

BOOST_AUTO_TEST_CASE(
	product_gauge_risk_receipt_is_immutable_and_lifetime_is_monotone)
{
	ProductGaugeCertificate certificate;
	certificate.system = E_Sys::GPS;
	certificate.firstObservable = E_ObsCode::L1C;
	certificate.secondObservable = E_ObsCode::L2W;
	certificate.satellite = SatSys(E_Sys::GPS, 3);
	certificate.reference = SatSys(E_Sys::GPS, 2);
	certificate.satellitePhaseSegments = {"G03-L1-SEG0", "G03-L2-SEG0"};
	certificate.referencePhaseSegments = {"G02-L1-SEG0", "G02-L2-SEG0"};
	certificate.admissionFamilyId = "IMMUTABLE-FAMILY";
	certificate.admissionFamilyFailureProbabilityBudget = 3e-4;
	certificate.failureProbabilityBudget = 2e-4;
	makeProductGaugeEvidenceComplete(certificate);
	ProductGaugeCertificateLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {certificate}, 1).valid);
	BOOST_REQUIRE_EQUAL(ledger.certificates().size(), 1);
	const auto original = ledger.certificates().front().riskReceipt;
	BOOST_CHECK_EQUAL(original.certificateFamilyId, "IMMUTABLE-FAMILY");
	BOOST_CHECK_CLOSE(ledger.lifetimeBroadcastRisk(), 3e-4, 1e-9);

	auto repeatedView = certificate;
	repeatedView.admissionFamilyFailureProbabilityBudget = 1e-4;
	repeatedView.factorWindowIdentity = "LATER-VIEW";
	BOOST_REQUIRE(ledger.observe(130, {repeatedView}, 1).valid);
	const auto& stored = ledger.certificates().front().riskReceipt;
	BOOST_CHECK_EQUAL(stored.exactLatticeHash, original.exactLatticeHash);
	BOOST_CHECK_EQUAL(stored.factorProvenanceHash,
		original.factorProvenanceHash);
	BOOST_CHECK(!stored.isNewIntegerDecision);
	BOOST_CHECK_CLOSE(ledger.lifetimeBroadcastRisk(), 3e-4, 1e-9);
}

BOOST_AUTO_TEST_CASE(r29_product_tree_is_transported_before_target_build)
{
	const SatSys g01(E_Sys::GPS, 1);
	const SatSys g02(E_Sys::GPS, 2);
	const std::set<ZhangGraphEdge> oldEdges = {
		{"R0", g01}, {"R0", g02}, {"R1", g01}, {"R1", g02}};
	const ZhangGraphBasis oldProduct = zhangBuildSpanningTree(
		oldEdges, "R0", {{"R0", g01}, {"R1", g01}, {"R1", g02}});
	const std::set<ZhangGraphEdge> currentEdges = {
		{"R0", g01}, {"R0", g02}, {"R1", g01}};
	const ZhangGraphBasis current = zhangBuildSpanningTree(
		currentEdges, "R0", {{"R0", g01}, {"R0", g02}, {"R1", g01}});
	BOOST_REQUIRE(oldProduct.connected && current.connected);
	const auto target = zhangBuildSatelliteProductTarget(
		current, oldProduct, g01);
	BOOST_REQUIRE_MESSAGE(target.valid, target.failureReason);
	BOOST_CHECK(target.productTreeTransported);
	BOOST_REQUIRE_EQUAL(target.staleProductTreeEdges.size(), 1);
	const ZhangGraphEdge staleExpected{"R1", g02};
	BOOST_CHECK(target.staleProductTreeEdges.front() == staleExpected);
}

BOOST_AUTO_TEST_CASE(r29_current_family_risk_is_deduplicated_without_floor)
{
	const auto audit = zhangAuditCertificateFamilyRisk({
		{"CURRENT-A", 2e-5, 3e-5, true},
		{"CURRENT-A", 4e-5, 3e-5, true},
		{"CURRENT-B", 1e-5, 2e-5, true}});
	BOOST_REQUIRE(audit.valid);
	BOOST_CHECK_EQUAL(audit.currentFamilies, 2);
	BOOST_CHECK_CLOSE(audit.currentCertificateRisk, 5e-5, 1e-9);
	BOOST_CHECK_CLOSE(audit.lifetimeBroadcastRisk, 5e-5, 1e-9);
	BOOST_CHECK_LT(audit.currentCertificateRisk, 1e-3);
}

BOOST_AUTO_TEST_CASE(r29_explicit_dual_bridge_closes_fragmented_components)
{
	// Four named nodes plus reference. Initially: (0,1) and (2,ref), with
	// node 3 a singleton. The explicit dual-frequency (1,2) target merges the
	// two active components without manufacturing an edge for node 3.
	ZhangExactMatrix rows{
		{1,-1,0,0, 0,0,0,0},
		{0,0,1,0,  0,0,0,0},
		{0,0,0,0,  1,-1,0,0},
		{0,0,0,0,  0,0,1,0}};
	ZhangExactVector values{0,0,0,0};
	const auto before = zhangExactDualGraphProfile(rows, values, 4);
	BOOST_REQUIRE(before.score.valid);
	BOOST_CHECK_NE(before.componentByNode[1], before.componentByNode[2]);
	rows.push_back({0,1,-1,0, 0,0,0,0});
	rows.push_back({0,0,0,0, 0,1,-1,0});
	values.push_back(0);
	values.push_back(0);
	const auto after = zhangExactDualGraphProfile(rows, values, 4);
	BOOST_REQUIRE(after.score.valid);
	BOOST_CHECK_EQUAL(after.score.graphRank, 3);
	BOOST_CHECK_EQUAL(after.score.satellites, 4);
	BOOST_CHECK_EQUAL(after.score.components, 1);
	BOOST_CHECK_EQUAL(after.score.largest, 4);
}

BOOST_AUTO_TEST_CASE(r31_original_decision_dependencies_survive_nested_views)
{
	auto parent=std::make_shared<ZhangIntegerDecisionProof>();
	parent->id="D1"; parent->originalStatement="2a-b=3";
	parent->observationProvenance="epoch1:factor1"; parent->conditionalFailureBound=2e-4;
	auto child=std::make_shared<ZhangIntegerDecisionProof>();
	child->id="D2"; child->originalStatement="c=4";
	child->observationProvenance="epoch2:factor2"; child->conditionalFailureBound=3e-4;
	child->parents={parent};
	const auto closure=zhangDecisionRiskClosure({parent,child,child});
	BOOST_REQUIRE(closure.valid);
	BOOST_CHECK_EQUAL(closure.atoms.size(),2);
	BOOST_CHECK_CLOSE(closure.bound,5e-4,1e-9);
	const auto audit=zhangAuditCertificateFamilyRisk({
		{"NETWORK_VIEW",.001,0,true,{parent}},
		{"GAUGE_HNF_VIEW",.001,0,true,{child}},
		{"TRANSPORTED_VIEW",.001,0,true,{child}}});
	BOOST_REQUIRE(audit.valid);
	BOOST_CHECK_CLOSE(audit.currentCertificateRisk,5e-4,1e-9);
	BOOST_CHECK_EQUAL(audit.currentFamilies,2);
	auto collision=std::make_shared<ZhangIntegerDecisionProof>(*parent);
	BOOST_CHECK(!zhangDecisionRiskClosure({parent,collision}).valid);
	BOOST_CHECK(!zhangAuditCertificateFamilyRisk({
		{"A",0,0,true,{parent}},{"B",0,0,true,{collision}}}).valid);
	BOOST_CHECK(!zhangDecisionRiskClosure({{}}).valid);
	auto unknown=std::make_shared<ZhangIntegerDecisionProof>();
	BOOST_CHECK(!zhangDecisionRiskClosure({unknown}).valid);
	auto other=std::make_shared<ZhangIntegerDecisionProof>(*parent);
	other->id="DISTINCT_DECISION_SAME_INTEGER";
	BOOST_CHECK_CLOSE(zhangDecisionRiskClosure({parent,other}).bound,4e-4,1e-9);
}

BOOST_AUTO_TEST_CASE(r31_incremental_target_image_preserves_congruences_and_rank_one)
{
	const auto even=zhangPrimitiveImageCoordinates({{2,0},{0,4}});
	BOOST_REQUIRE(even.valid);
	BOOST_CHECK(zhangExactMultiply(even.imageGenerators,even.primitiveRows)==ZhangExactMatrix({{2,0},{0,4}}));
	BOOST_CHECK(zhangExactPrimitiveRowLattice(even.primitiveRows,2));
	const auto baseline=zhangExactAffineIntegerQuotient({{1,0}},{7},2);
	BOOST_REQUIRE(baseline.valid);
	BOOST_CHECK_EQUAL(baseline.quotientRank,1);
	ZhangExactMatrix kernel{{baseline.kernelBasis[0][0]},{baseline.kernelBasis[0][1]}};
	const auto missing=zhangPrimitiveImageCoordinates(zhangExactMultiply({{1,0},{0,2}},kernel));
	BOOST_REQUIRE(missing.valid);
	BOOST_CHECK_EQUAL(missing.primitiveRows.size(),1);
	BOOST_CHECK(zhangExactMultiply(missing.imageGenerators,missing.primitiveRows)==zhangExactMultiply({{1,0},{0,2}},kernel));
	const auto zero=zhangPrimitiveImageCoordinates({{0,0},{0,0}});
	BOOST_REQUIRE(zero.valid);
	BOOST_CHECK(zero.primitiveRows.empty());
	BOOST_CHECK(!zhangPrimitiveImageCoordinates({{1,0},{1}}).valid);
}

BOOST_AUTO_TEST_CASE(r31_family_subsets_deduplicate_nested_original_decisions)
{
	auto a=std::make_shared<ZhangIntegerDecisionProof>();
	a->id="A"; a->originalStatement="a=1"; a->observationProvenance="raw1";
	a->conditionalFailureBound=4e-4;
	auto b=std::make_shared<ZhangIntegerDecisionProof>();
	b->id="B"; b->originalStatement="b=2|a=1"; b->observationProvenance="raw2";
	b->conditionalFailureBound=4e-4; b->parents={a};
	const auto selected=zhangSelectProductGaugeAdmissionFamilies(
		std::vector<double>{4e-4,8e-4},1e-3,[](const std::vector<int>& ids)
		{
			ZhangProductGaugeFamilySubsetScore score;
			score.valid=!ids.empty(); score.dualGraphRank=ids.size();
			score.largestCertifiedComponentSize=ids.size()+1;
			score.certifiedSatelliteCount=ids.size()+1;
			return score;
		},8,64,{}, {},{{a},{b}});
	BOOST_REQUIRE(selected.valid);
	BOOST_CHECK_EQUAL(selected.selectedFamilyIndices.size(),2);
	BOOST_CHECK_CLOSE(selected.failureProbability,8e-4,1e-9);
}

BOOST_AUTO_TEST_CASE(r32_unpublishable_posterior_still_charges_parent_risk)
{
	// R31 00:15:30: publishability=false did NOT mean a FLOAT-only posterior.
	const double parent=0.0008849583511821857, extra=0.00014087723960320986;
	for(bool publishable : {false,true})
	{
		const auto risk=zhangSelectOptionalLedgerWithoutPoisoningBaseline(publishable,parent,true,extra,.001);
		BOOST_CHECK(!risk.optionalAuthorized);
		BOOST_CHECK_CLOSE(risk.proposedUnionBound,parent+extra,1e-9);
		BOOST_CHECK_CLOSE(risk.selectedFailureProbability,parent,1e-9);
		const auto transaction=zhangSelectPosteriorReceiptTransaction(true,risk.optionalAuthorized);
		BOOST_CHECK(transaction.rollbackRequired);
		BOOST_CHECK(!transaction.useConditionedPosterior);
	}
}

BOOST_AUTO_TEST_CASE(r32_unpublishable_parent_valid_boundary_and_invalid_inputs)
{
	const auto safe=zhangSelectOptionalLedgerWithoutPoisoningBaseline(false,4e-4,true,2e-4,.001);
	BOOST_CHECK(safe.optionalAuthorized);
	BOOST_CHECK_CLOSE(safe.selectedFailureProbability,6e-4,1e-9);
	BOOST_CHECK(!safe.baselineFeasible); // No claim that a product exists.
	const auto none=zhangSelectOptionalLedgerWithoutPoisoningBaseline(false,4e-4,false,1,.001);
	BOOST_CHECK_CLOSE(none.selectedFailureProbability,4e-4,1e-9);
	BOOST_CHECK(!none.optionalAuthorized);
	BOOST_CHECK(zhangSelectOptionalLedgerWithoutPoisoningBaseline(false,.001,true,0,.001).optionalAuthorized);
	BOOST_CHECK(!zhangSelectOptionalLedgerWithoutPoisoningBaseline(false,.001,true,1e-15,.001).optionalAuthorized);
	for(double invalid : {-1.0,.0011,std::numeric_limits<double>::quiet_NaN()})
		BOOST_CHECK(!zhangSelectOptionalLedgerWithoutPoisoningBaseline(false,invalid,true,0,.001).optionalAuthorized);
}

BOOST_AUTO_TEST_CASE(r32_full_dependency_union_deduplicates_but_retains_conditional_parents)
{
	auto p=std::make_shared<ZhangIntegerDecisionProof>();
	p->id="R32-P"; p->originalStatement="a=1"; p->observationProvenance="raw1"; p->conditionalFailureBound=4e-4;
	auto c=std::make_shared<ZhangIntegerDecisionProof>();
	c->id="R32-C"; c->originalStatement="b=2|a=1"; c->observationProvenance="raw2"; c->conditionalFailureBound=3e-4; c->parents={p};
	const auto total=zhangDecisionRiskClosure({p,c,p});
	BOOST_REQUIRE(total.valid);
	BOOST_CHECK_EQUAL(total.atoms.size(),2);
	const auto risk=zhangSelectOptionalLedgerWithoutPoisoningBaseline(false,p->conditionalFailureBound,true,total.bound-p->conditionalFailureBound,.001);
	BOOST_CHECK(risk.optionalAuthorized);
	BOOST_CHECK_CLOSE(risk.selectedFailureProbability,7e-4,1e-9);
}

BOOST_AUTO_TEST_CASE(r32_late_receipt_failure_requires_atomic_restore_before_retry)
{
	VectorXd before(2); before<<3,5;
	MatrixXd prior=MatrixXd::Identity(2,2);
	VectorXd working=before; working(0)=4;
	MatrixXd posterior=prior; posterior(0,0)=0;
	const auto lateRisk=zhangAllocateIncrementalRiskBudget(.0011,.001);
	const auto late=zhangSelectPosteriorReceiptTransaction(true,lateRisk.parentFeasible);
	BOOST_REQUIRE(late.rollbackRequired);
	if(late.rollbackRequired) { working=before; posterior=prior; }
	BOOST_CHECK((working-before).isZero());
	BOOST_CHECK((posterior-prior).isZero());
	const auto retry=zhangSelectPosteriorReceiptTransaction(false,false);
	BOOST_CHECK(!retry.rollbackRequired && !retry.useConditionedPosterior);
	BOOST_CHECK(zhangSelectPosteriorReceiptTransaction(true,true).useConditionedPosterior);
}

BOOST_AUTO_TEST_CASE(r32_row_proof_hnf_does_not_spread_independent_risk)
{
	auto a=std::make_shared<ZhangIntegerDecisionProof>();
	a->id="R32-A"; a->originalStatement="a=2"; a->observationProvenance="rawA"; a->conditionalFailureBound=2e-4;
	auto b=std::make_shared<ZhangIntegerDecisionProof>();
	b->id="R32-B"; b->originalStatement="b=3"; b->observationProvenance="rawB"; b->conditionalFailureBound=8e-4;
	ZhangExactMatrix rows{{1,0},{0,1}}; ZhangExactVector values{2,3};
	const auto hnf=zhangExactRowHermiteNormalForm(rows,values,true);
	const auto binding=zhangBindExactRowDecisionProofs(rows,values,{{a},{b}},hnf.basis,hnf.values,&hnf.rowTransform);
	BOOST_REQUIRE(binding.valid);
	BOOST_REQUIRE_EQUAL(binding.rowProofs.size(),2);
	for(size_t i=0;i<hnf.basis.size();++i)
	{
		const auto closure=zhangDecisionRiskClosure(binding.rowProofs[i]);
		BOOST_REQUIRE(closure.valid);
		BOOST_CHECK_EQUAL(closure.atoms.size(),1);
		BOOST_CHECK_CLOSE(closure.bound,hnf.basis[i][0]!=0 ? 2e-4:8e-4,1e-9);
	}
	// Survival of only a must not retain independent decision B, even repeatedly.
	auto surviving=zhangBindExactRowDecisionProofs(rows,values,{{a},{b}},{{1,0}},{2});
	BOOST_REQUIRE(surviving.valid);
	for(int step=0;step<10;++step)
	{
		surviving=zhangBindExactRowDecisionProofs({{1,0}},{2},surviving.rowProofs,{{1,0}},{2});
		BOOST_REQUIRE(surviving.valid);
		BOOST_CHECK_EQUAL(zhangDecisionRiskClosure(surviving.rowProofs[0]).atoms.size(),1);
	}
}

BOOST_AUTO_TEST_CASE(r32_survival_cancellation_keeps_used_proofs_and_real_ancestors)
{
	auto a=std::make_shared<ZhangIntegerDecisionProof>();
	a->id="R32-X"; a->originalStatement="x+y=5"; a->observationProvenance="rawX"; a->conditionalFailureBound=2e-4;
	auto b=std::make_shared<ZhangIntegerDecisionProof>();
	b->id="R32-Y"; b->originalStatement="y=3|X"; b->observationProvenance="rawY"; b->conditionalFailureBound=3e-4; b->parents={a};
	const auto cancellation=zhangBindExactRowDecisionProofs({{1,1},{0,1}},{5,3},{{a},{b}},{{1,0}},{2});
	BOOST_REQUIRE(cancellation.valid);
	BOOST_CHECK_CLOSE(zhangDecisionRiskClosure(cancellation.rowProofs[0]).bound,5e-4,1e-9);
	const auto childOnly=zhangBindExactRowDecisionProofs({{1,1},{0,1}},{5,3},{{a},{b}},{{0,1}},{3});
	BOOST_REQUIRE(childOnly.valid);
	BOOST_CHECK_EQUAL(childOnly.rowProofs[0].size(),1);
	BOOST_CHECK_CLOSE(zhangDecisionRiskClosure(childOnly.rowProofs[0]).bound,5e-4,1e-9);
	const auto survival=zhangExactSurvivingLattice({{1,1},{0,1}},{5,3},{true,false},true);
	BOOST_REQUIRE(survival.consistent);
	BOOST_REQUIRE_EQUAL(survival.basis.size(),1);
	const auto native=zhangBindExactRowDecisionProofs({{1,1},{0,1}},{5,3},{{a},{b}},
		{{survival.basis[0][0],0}},survival.values,&survival.rowTransform);
	BOOST_REQUIRE(native.valid);
	BOOST_CHECK_CLOSE(zhangDecisionRiskClosure(native.rowProofs[0]).bound,5e-4,1e-9);
	const auto noLoss=zhangExactSurvivingLattice({{1,0},{0,1}},{2,3},{true,true},true);
	BOOST_CHECK(zhangBindExactRowDecisionProofs({{1,0},{0,1}},{2,3},{{a},{b}},
		noLoss.basis,noLoss.values,&noLoss.rowTransform).valid);
}

BOOST_AUTO_TEST_CASE(r32_row_proofs_fail_closed_on_affine_lattice_or_provenance_mismatch)
{
	auto a=std::make_shared<ZhangIntegerDecisionProof>();
	a->id="R32-VALID"; a->originalStatement="2x=4"; a->observationProvenance="raw"; a->conditionalFailureBound=2e-4;
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{2}},{4},{{a}},{{1}},{2}).valid);
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{2}},{4},{{a}},{{2}},{5}).valid);
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{2}},{4},{{}},{{2}},{4}).valid);
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{2}},{4},{},{{2}},{4}).valid);
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{2}},{4},{{a}},{{2,0}},{4}).valid);
	const ZhangExactMatrix invalidTransform{{2}};
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{2}},{4},{{a}},{{2}},{4},&invalidTransform).valid);
	auto collision=std::make_shared<ZhangIntegerDecisionProof>(*a);
	BOOST_CHECK(!zhangBindExactRowDecisionProofs({{1,0},{0,1}},{2,3},{{a},{collision}},{{1,1}},{5}).valid);
}

BOOST_AUTO_TEST_CASE(r33_canonical_product_expansion_is_reference_invariant)
{
	std::map<std::string, ZhangExactInteger> referenceG03;
	BOOST_REQUIRE(zhangBuildCanonicalProductExpansion(
		{1, -1, 0, 0}, {"G01", "G02"}, "G03",
		E_ObsCode::L1C, E_ObsCode::L2W, referenceG03));

	std::map<std::string, ZhangExactInteger> referenceG01;
	BOOST_REQUIRE(zhangBuildCanonicalProductExpansion(
		{0, -1, 0, 0}, {"G03", "G02"}, "G01",
		E_ObsCode::L1C, E_ObsCode::L2W, referenceG01));

	BOOST_CHECK(referenceG03 == referenceG01);
	BOOST_REQUIRE_EQUAL(referenceG03.size(), 2);
	BOOST_CHECK_EQUAL(referenceG03.at(
		zhangProductCanonicalCoordinateKey(E_ObsCode::L1C, "G01")), 1);
	BOOST_CHECK_EQUAL(referenceG03.at(
		zhangProductCanonicalCoordinateKey(E_ObsCode::L1C, "G02")), -1);
}

BOOST_AUTO_TEST_CASE(r33_incremental_decision_risk_charges_only_new_proof_ids)
{
	auto parent = std::make_shared<ZhangIntegerDecisionProof>();
	parent->id = "R33-PARENT";
	parent->originalStatement = "parent integer";
	parent->observationProvenance = "epoch-1";
	parent->conditionalFailureBound = 4e-4;

	auto child = std::make_shared<ZhangIntegerDecisionProof>();
	child->id = "R33-CHILD";
	child->originalStatement = "child conditional integer";
	child->observationProvenance = "epoch-2";
	child->conditionalFailureBound = 2e-4;
	child->parents = {parent};

	const auto incremental = zhangIncrementalDecisionRisk({parent}, {child});
	BOOST_REQUIRE(incremental.valid);
	BOOST_CHECK_EQUAL(incremental.baselineAtoms, 1);
	BOOST_CHECK_EQUAL(incremental.unionAtoms, 2);
	BOOST_CHECK_EQUAL(incremental.addedAtoms, 1);
	BOOST_CHECK_CLOSE(incremental.bound, 2e-4, 1e-9);

	auto collision = std::make_shared<ZhangIntegerDecisionProof>(*parent);
	const auto rejected = zhangIncrementalDecisionRisk({parent}, {collision});
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.reason, "DECISION_ID_COLLISION");
}

BOOST_AUTO_TEST_CASE(r35_conflict_hypothesis_requires_later_independent_epoch)
{
	ZhangConflictAwareHypothesisSet hypotheses;
	const auto first = hypotheses.observe("H1", 100, "E1", true);
	BOOST_CHECK(first.inserted);
	BOOST_CHECK(!first.confirmed);
	BOOST_CHECK_EQUAL(first.independentEpochs, 1);

	const auto sameEpoch = hypotheses.observe("H1", 100, "E2", true);
	BOOST_CHECK(sameEpoch.sameEpoch);
	BOOST_CHECK(!sameEpoch.confirmed);

	const auto dependent = hypotheses.observe("H1", 130, "E1", true);
	BOOST_CHECK(dependent.status ==
		ZhangConflictHypothesisStatus::DEPENDENT_EVIDENCE_PENDING);
	BOOST_CHECK(!dependent.confirmed);

	const auto confirmed = hypotheses.observe("H1", 160, "E3", true);
	BOOST_CHECK(confirmed.status ==
		ZhangConflictHypothesisStatus::CONFIRMED);
	BOOST_CHECK(confirmed.confirmed);
	BOOST_CHECK_EQUAL(confirmed.independentEpochs, 2);
}

BOOST_AUTO_TEST_CASE(r35_conflict_hypothesis_capacity_is_fail_closed)
{
	ZhangConflictAwareHypothesisSet hypotheses;
	for (int index = 0; index < 4; index++)
		BOOST_CHECK(hypotheses.observe(
			"H" + std::to_string(index), 100,
			"E" + std::to_string(index), true).inserted);
	const auto rejected = hypotheses.observe("H4", 100, "E4", true);
	BOOST_CHECK(rejected.status ==
		ZhangConflictHypothesisStatus::CAPACITY_REJECTED);
	BOOST_CHECK_EQUAL(rejected.candidateCount, 4);
}

BOOST_AUTO_TEST_CASE(r35_component_forest_generates_bounded_joint_families)
{
	const std::vector<ZhangConflictAwareForestEdge> edges = {
		{"A", 0, 1, 0, 1, 8, 0.2},
		{"B", 1, 2, 1, 2, 9, 0.3},
		{"C", 2, 3, 2, 3, 10, 0.4},
		{"D", 0, 3, 0, 3, 7, 0.1},
		{"SAME", 4, 5, 4, 4, 99, 0.01}};
	const auto families = zhangGenerateConflictAwareForestFamilies(edges, 8, 64);
	BOOST_REQUIRE(!families.empty());
	BOOST_CHECK_LE(families.size(), 8);
	bool foundTwo = false;
	bool foundThree = false;
	for (const auto& family : families)
	{
		BOOST_CHECK(family.valid);
		BOOST_CHECK_GE(family.inputEdgeIndices.size(), 2);
		BOOST_CHECK_LE(family.inputEdgeIndices.size(), 3);
		BOOST_CHECK_EQUAL(
			family.trueMergeCount, family.inputEdgeIndices.size());
		for (const auto index : family.inputEdgeIndices)
			BOOST_CHECK_LT(index, 4);
		foundTwo |= family.inputEdgeIndices.size() == 2;
		foundThree |= family.inputEdgeIndices.size() == 3;
	}
	BOOST_CHECK(foundTwo);
	BOOST_CHECK(foundThree);
}

BOOST_AUTO_TEST_CASE(r35_component_forest_selection_forbids_graph_regression)
{
	ZhangConflictAwareProductBranchScore baseline;
	baseline.valid = baseline.safe = true;
	baseline.exactLargestComponentSize = 10;
	baseline.exactSatelliteCount = 16;
	baseline.exactGraphRank = 12;
	baseline.exactLatticeRank = 20;
	baseline.exactComponentCount = 3;
	baseline.temporalOverlap = 6;
	baseline.failureProbability = 4e-4;
	baseline.fingerprint = "BASE";

	auto regressive = baseline;
	regressive.exactLargestComponentSize = 9;
	regressive.exactGraphRank = 15;
	regressive.fingerprint = "REGRESSIVE";
	auto improved = baseline;
	improved.exactLargestComponentSize = 11;
	improved.exactGraphRank = 13;
	improved.exactLatticeRank = 21;
	improved.fingerprint = "IMPROVED";
	const auto selection = zhangSelectConflictAwareProductBranch(
		baseline, {regressive, improved});
	BOOST_REQUIRE(selection.valid);
	BOOST_CHECK_EQUAL(selection.selectedCandidateIndex, 1);
	BOOST_CHECK_EQUAL(selection.selected.fingerprint, "IMPROVED");
}

BOOST_AUTO_TEST_CASE(r35_product_ledger_rejects_epoch_regression_transactionally)
{
	ProductIntegerLedger ledger;
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.firstObservable = E_ObsCode::L1C;
	row.secondObservable = E_ObsCode::L2W;
	row.productRow = {1, -1};
	row.integerValue = 9;
	row.physicalExpansion = {{"L1C|A|G01|V1", 1}, {"L1C|B|G02|V1", -1}};
	row.canonicalProductExpansion = {{"L1C|G01", 1}, {"L1C|G02", -1}};
	row.phaseSegmentFingerprint = "G01|L1C|S1;G02|L1C|S1;";
	row.admissionFailureProbabilityBound = 2e-5;
	BOOST_REQUIRE(ledger.observe(100, {row}, 2).valid);
	BOOST_REQUIRE(ledger.observe(130, {row}, 2).valid);
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	const auto before = ledger.rows().front();

	const auto rejected = ledger.observe(120, {row}, 2);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.failureReason,
		"PRODUCT_LEDGER_EPOCH_REGRESSION");
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), 1);
	const auto& after = ledger.rows().front();
	BOOST_CHECK_EQUAL(after.integerValue, before.integerValue);
	BOOST_CHECK_EQUAL(after.lastConfirmed, before.lastConfirmed);
	BOOST_CHECK_EQUAL(after.confirmationEpochs, before.confirmationEpochs);
	BOOST_CHECK_EQUAL(after.certified, before.certified);
	BOOST_CHECK(after.physicalExpansion == before.physicalExpansion);
	BOOST_CHECK(after.canonicalProductExpansion ==
		before.canonicalProductExpansion);

	const auto sameEpoch = ledger.observe(130, {row}, 2);
	BOOST_REQUIRE(sameEpoch.valid);
	BOOST_CHECK_EQUAL(
		ledger.rows().front().confirmationEpochs, before.confirmationEpochs);
}

BOOST_AUTO_TEST_CASE(r35_product_ledger_rejects_unknown_admission_risk)
{
	auto makeRow = []()
	{
		ProductIntegerLedgerRow row;
		row.system = E_Sys::GPS;
		row.firstObservable = E_ObsCode::L1C;
		row.secondObservable = E_ObsCode::L2W;
		row.productRow = {1, -1};
		row.integerValue = 4;
		row.physicalExpansion = {
			{"L1C|A|G01|V1", 1}, {"L1C|B|G02|V1", -1}};
		row.canonicalProductExpansion = {
			{"L1C|G01", 1}, {"L1C|G02", -1}};
		row.phaseSegmentFingerprint = "G01|L1C|S1;G02|L1C|S1;";
		return row;
	};
	const std::vector<double> invalid = {
		std::numeric_limits<double>::quiet_NaN(),
		std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity(), -1e-6, 1.000001};
	for (const double bound : invalid)
	{
		ProductIntegerLedger ledger;
		auto row = makeRow();
		row.admissionFailureProbabilityBound = bound;
		const auto rejected = ledger.observe(100, {row}, 1);
		BOOST_CHECK(!rejected.valid);
		BOOST_CHECK_EQUAL(rejected.failureReason,
			"PRODUCT_LEDGER_ADMISSION_FAILURE_BOUND_INVALID");
		BOOST_CHECK(ledger.rows().empty());
	}
}

BOOST_AUTO_TEST_CASE(
	r35_product_ledger_invalid_batch_rolls_back_all_existing_fields)
{
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.firstObservable = E_ObsCode::L1C;
	row.secondObservable = E_ObsCode::L2W;
	row.productRow = {1, -1};
	row.integerValue = 4;
	row.physicalExpansion = {
		{"L1C|A|G01|V1", 1}, {"L1C|B|G02|V1", -1}};
	row.canonicalProductExpansion = {
		{"L1C|G01", 1}, {"L1C|G02", -1}};
	row.phaseSegmentFingerprint = "G01|L1C|S1;G02|L1C|S1;";
	row.admissionFailureProbabilityBound = 2e-5;
	ProductIntegerLedger ledger;
	BOOST_REQUIRE(ledger.observe(100, {row}, 1).valid);
	const auto before = ledger.rows();

	auto invalid = row;
	invalid.canonicalProductExpansion = {
		{"L1C|G03", 1}, {"L1C|G04", -1}};
	invalid.physicalExpansion = {
		{"L1C|A|G03|V1", 1}, {"L1C|B|G04|V1", -1}};
	invalid.phaseSegmentFingerprint = "G03|L1C|S1;G04|L1C|S1;";
	invalid.admissionFailureProbabilityBound =
		std::numeric_limits<double>::quiet_NaN();
	const auto rejected = ledger.observe(130, {row, invalid}, 1);
	BOOST_CHECK(!rejected.valid);
	BOOST_CHECK_EQUAL(rejected.failureReason,
		"PRODUCT_LEDGER_ADMISSION_FAILURE_BOUND_INVALID");
	BOOST_REQUIRE_EQUAL(ledger.rows().size(), before.size());
	BOOST_CHECK(ledger.rows()[0].decisionProofs == before[0].decisionProofs);
	BOOST_CHECK(ledger.rows()[0].system == before[0].system);
	BOOST_CHECK(ledger.rows()[0].firstObservable == before[0].firstObservable);
	BOOST_CHECK(ledger.rows()[0].secondObservable == before[0].secondObservable);
	BOOST_CHECK(ledger.rows()[0].productRow == before[0].productRow);
	BOOST_CHECK_EQUAL(ledger.rows()[0].lastConfirmed,
		before[0].lastConfirmed);
	BOOST_CHECK_EQUAL(ledger.rows()[0].confirmationEpochs,
		before[0].confirmationEpochs);
	BOOST_CHECK(ledger.rows()[0].physicalExpansion ==
		before[0].physicalExpansion);
	BOOST_CHECK(ledger.rows()[0].canonicalProductExpansion ==
		before[0].canonicalProductExpansion);
	BOOST_CHECK_EQUAL(ledger.rows()[0].integerValue,
		before[0].integerValue);
	BOOST_CHECK_EQUAL(ledger.rows()[0].coordinate, before[0].coordinate);
	BOOST_CHECK_EQUAL(ledger.rows()[0].firstSatellite,
		before[0].firstSatellite);
	BOOST_CHECK_EQUAL(ledger.rows()[0].secondSatellite,
		before[0].secondSatellite);
	BOOST_CHECK_EQUAL(ledger.rows()[0].phaseSegmentFingerprint,
		before[0].phaseSegmentFingerprint);
	BOOST_CHECK_EQUAL(ledger.rows()[0].backendBasisGeneration,
		before[0].backendBasisGeneration);
	BOOST_CHECK_EQUAL(ledger.rows()[0].admissionFailureProbabilityBound,
		before[0].admissionFailureProbabilityBound);
	BOOST_CHECK_EQUAL(ledger.rows()[0].firstCertified,
		before[0].firstCertified);
	BOOST_CHECK(ledger.rows()[0].source == before[0].source);
	BOOST_CHECK_EQUAL(ledger.rows()[0].conditioningOnly,
		before[0].conditioningOnly);
	BOOST_CHECK_EQUAL(ledger.rows()[0].pairCertificate,
		before[0].pairCertificate);
	BOOST_CHECK_EQUAL(ledger.rows()[0].certified, before[0].certified);
}

BOOST_AUTO_TEST_CASE(
	r35_product_ledger_unrelated_history_does_not_poison_fresh_baseline)
{
	BOOST_CHECK(zhangProductLedgerPresearchAllowsFreshWithoutHistory(
		"NO_CERTIFIED_LEDGER_ROWS"));
	BOOST_CHECK(zhangProductLedgerPresearchAllowsFreshWithoutHistory(
		"NO_CURRENT_SEGMENT_PROJECTABLE_ROWS"));
	BOOST_CHECK(zhangProductLedgerPresearchAllowsFreshWithoutHistory(
		"NO_LEDGER_ROWS_IN_CURRENT_PRODUCT_IMAGE"));
	BOOST_CHECK(!zhangProductLedgerPresearchAllowsFreshWithoutHistory(
		"LEDGER_DETERMINISTIC_AFFINE_CONFLICT"));
	BOOST_CHECK(!zhangProductLedgerPresearchAllowsFreshWithoutHistory(
		"NO_NIS_COMPATIBLE_LEDGER_ROWS"));
	BOOST_CHECK(!zhangProductLedgerPresearchAllowsFreshWithoutHistory(
		"CANONICAL_PRODUCT_TRANSPORT_COMMUTATIVE_SQUARE_FAILED"));
}

BOOST_AUTO_TEST_CASE(r35_product_ledger_sign_authority_is_canonical)
{
	ProductIntegerLedgerRow row;
	row.system = E_Sys::GPS;
	row.productRow = {1, -1};
	row.integerValue = 7;
	// Backend chart begins positive, but the immutable canonical image begins
	// negative.  Canonical orientation must win.
	row.physicalExpansion = {{"A_BACKEND", 1}, {"B_BACKEND", -1}};
	row.canonicalProductExpansion = {{"A_CANONICAL", -1}, {"B_CANONICAL", 1}};
	row.pairCertificate = true;
	row.firstSatellite = "G01";
	row.secondSatellite = "G02";
	auto negated = row;
	for (auto& coefficient : negated.productRow) coefficient = -coefficient;
	for (auto& [identity, coefficient] : negated.physicalExpansion)
		coefficient = -coefficient;
	for (auto& [identity, coefficient] : negated.canonicalProductExpansion)
		coefficient = -coefficient;
	negated.integerValue = -negated.integerValue;
	std::swap(negated.firstSatellite, negated.secondSatellite);

	zhangCanonicaliseProductLedgerRow(row);
	zhangCanonicaliseProductLedgerRow(negated);
	BOOST_CHECK(row.canonicalProductExpansion.begin()->second > 0);
	BOOST_CHECK_EQUAL(row.physicalExpansion.at("A_BACKEND"), -1);
	BOOST_CHECK_EQUAL(row.productRow.front(), -1);
	BOOST_CHECK_EQUAL(row.integerValue, -7);
	BOOST_CHECK(row.physicalExpansion == negated.physicalExpansion);
	BOOST_CHECK(row.canonicalProductExpansion ==
		negated.canonicalProductExpansion);
	BOOST_CHECK(row.productRow == negated.productRow);
	BOOST_CHECK_EQUAL(row.integerValue, negated.integerValue);
	BOOST_CHECK_EQUAL(
		zhangProductLedgerIdentityFingerprint(row),
		zhangProductLedgerIdentityFingerprint(negated));
}

BOOST_AUTO_TEST_CASE(r35_conflict_hypothesis_ttl_uses_independent_evidence)
{
	ZhangConflictAwareHypothesisSet hypotheses;
	BOOST_REQUIRE(hypotheses.observe("H1", 100, "E1", true).inserted);
	BOOST_CHECK(hypotheses.observe("H1", 100, "E2", true).sameEpoch);
	const auto dependent = hypotheses.observe("H1", 130, "E1", true);
	BOOST_CHECK(dependent.dependentEvidence);
	BOOST_CHECK_EQUAL(dependent.lastIndependentEvidenceEpoch, 100);
	BOOST_CHECK_EQUAL(hypotheses.eraseExpiredBefore(101), 1);
	BOOST_CHECK(hypotheses.candidates().empty());

	BOOST_REQUIRE(hypotheses.observe("H2", 200, "F1", true).inserted);
	const auto confirmed = hypotheses.observe("H2", 230, "F2", true);
	BOOST_REQUIRE(confirmed.confirmed);
	BOOST_CHECK_EQUAL(confirmed.confirmationEpoch, 230);
	BOOST_CHECK(!hypotheses.consumeConfirmed("UNKNOWN"));
	BOOST_CHECK(hypotheses.consumeConfirmed("H2"));
	BOOST_CHECK(hypotheses.candidates().empty());

	BOOST_REQUIRE(hypotheses.observe("H3", 300, "G1", true).inserted);
	BOOST_REQUIRE(hypotheses.observe("H3", 330, "G2", true).confirmed);
	BOOST_CHECK_EQUAL(hypotheses.eraseExpiredBefore(331), 1);
	BOOST_CHECK(hypotheses.candidates().empty());
}

BOOST_AUTO_TEST_CASE(r35_component_branch_rejects_degenerate_candidate)
{
	ZhangConflictAwareProductBranchScore emptyBaseline;
	emptyBaseline.valid = emptyBaseline.safe = true;
	emptyBaseline.failureProbability = 0;
	emptyBaseline.fingerprint = "EMPTY_BASELINE";

	auto zeroGraph = emptyBaseline;
	zeroGraph.exactLargestComponentSize = 2;
	zeroGraph.exactSatelliteCount = 2;
	zeroGraph.exactLatticeRank = 1;
	zeroGraph.fingerprint = "ZERO_GRAPH";
	const auto zeroGraphSelection = zhangSelectConflictAwareProductBranch(
		emptyBaseline, {zeroGraph});
	BOOST_REQUIRE(zeroGraphSelection.valid);
	BOOST_CHECK_EQUAL(zeroGraphSelection.selectedCandidateIndex, -1);
	BOOST_CHECK_EQUAL(zeroGraphSelection.selected.fingerprint, "EMPTY_BASELINE");

	auto zeroLattice = zeroGraph;
	zeroLattice.exactGraphRank = 1;
	zeroLattice.exactLatticeRank = 0;
	zeroLattice.fingerprint = "ZERO_LATTICE";
	const auto zeroLatticeSelection = zhangSelectConflictAwareProductBranch(
		emptyBaseline, {zeroLattice});
	BOOST_CHECK_EQUAL(zeroLatticeSelection.selectedCandidateIndex, -1);

	auto oneSatellite = zeroLattice;
	oneSatellite.exactLatticeRank = 1;
	oneSatellite.exactSatelliteCount = 1;
	oneSatellite.fingerprint = "ONE_SATELLITE";
	const auto oneSatelliteSelection = zhangSelectConflictAwareProductBranch(
		emptyBaseline, {oneSatellite});
	BOOST_CHECK_EQUAL(oneSatelliteSelection.selectedCandidateIndex, -1);

	auto valid = zeroLattice;
	valid.exactLatticeRank = 1;
	valid.exactSatelliteCount = 2;
	valid.fingerprint = "VALID";
	const auto validSelection = zhangSelectConflictAwareProductBranch(
		emptyBaseline, {valid});
	BOOST_REQUIRE(validSelection.valid);
	BOOST_CHECK_EQUAL(validSelection.selectedCandidateIndex, 0);
	BOOST_CHECK_EQUAL(validSelection.selected.fingerprint, "VALID");
}

BOOST_AUTO_TEST_CASE(
	r43_historical_whole_lattice_pair_authorization_survives_catalogue_rebuild)
{
	auto proof = std::make_shared<ZhangIntegerDecisionProof>();
	proof->id = "R43-HISTORICAL-FAMILY";
	proof->originalStatement = "mixed WL/L1 affine lattice";
	proof->observationProvenance = "epoch-previous:raw-factors";
	proof->conditionalFailureBound = 2e-4;

	ZhangProductIntegerConstraintSet constraints;
	constraints.productCoordinateDimension = 1;
	constraints.phaseSegmentFingerprint = "G01|L1C|SEG7;G02|L1C|SEG4;";
	// Neither catalogue provenance nor row identity is used as authority.  The
	// first historical row is mixed; only the whole-lattice combination
	// [1,1]-[0,1]=[1,0] exposes the WL pair.
	constraints.currentReauthorizedHistoricalDualRows = {{1, 1}, {0, 1}};
	constraints.currentReauthorizedHistoricalDualValues = {10, 7};
	constraints.currentReauthorizedHistoricalDecisionProofs = {proof};
	constraints.currentReauthorizedHistoricalFamilyIds = {"FAMILY-7"};
	constraints.currentReauthorizedHistoricalSnapshotId = "SNAPSHOT-7";
	constraints.currentReauthorizedSourceMomentId = "SOURCE-MOMENT-7";
	constraints.currentReauthorizedHistoricalFailureProbability = 2e-4;
	constraints.currentReauthorizedHistoricalComponentVersion = 7;
	constraints.currentReauthorizedHistoricalSegmentsValidated = true;
	constraints.currentReauthorizedHistoricalExactTransport = true;

	ZhangCertifiedPairRelation wide;
	wide.firstNode = 0;
	wide.secondNode = 1;
	wide.value = 3;
	wide.coordinate = "WL";
	constraints.certifiedPairs = {wide};
	ZhangCertifiedPairRelation dual;
	dual.firstNode = 0;
	dual.secondNode = 1;
	dual.value = 7;
	dual.coordinate = "L1_AND_WL";
	// Model the final catalogue rebuild erasing transient provenance.
	dual.fromProductGaugeLedger = false;
	constraints.dualFrequencyCertifiedPairs = {dual};

	BOOST_CHECK_EQUAL(zhangFinalizeHistoricalProductPairAuthorization(
		constraints, "DELIVERY-MOMENT-7"), 1);
	BOOST_REQUIRE_EQUAL(constraints.dualFrequencyCertifiedPairs.size(), 1);
	const auto& authorised = constraints.dualFrequencyCertifiedPairs.front();
	BOOST_CHECK(authorised.fromProductGaugeLedger);
	BOOST_CHECK(authorised.currentPosteriorReauthorized);
	BOOST_CHECK(authorised.exactHistoricalTransportWitness);
	BOOST_CHECK_EQUAL(authorised.authorizedBackendBasisGeneration,
		constraints.backendBasisGeneration);
	BOOST_CHECK_EQUAL(authorised.authorizedFamilyIds.size(), 1);
	BOOST_CHECK(!authorised.authorizedDecisionIds.empty());
	BOOST_CHECK_EQUAL(authorised.productGaugeComponentVersion, 7);
	BOOST_CHECK_EQUAL(authorised.productGaugeSnapshotId, "SNAPSHOT-7");
	BOOST_CHECK_EQUAL(authorised.validatedDeliveryMomentId,
		"DELIVERY-MOMENT-7");
	BOOST_CHECK_CLOSE(authorised.authorizedFailureProbability, 2e-4, 1e-9);
}

BOOST_AUTO_TEST_CASE(
	r43_historical_pair_authorization_fails_closed_on_value_segment_or_risk)
{
	auto proof = std::make_shared<ZhangIntegerDecisionProof>();
	proof->id = "R43-RISK";
	proof->originalStatement = "pair";
	proof->observationProvenance = "raw";
	proof->conditionalFailureBound = 2e-4;
	auto makeConstraints = [&]()
	{
		ZhangProductIntegerConstraintSet constraints;
		constraints.productCoordinateDimension = 1;
		constraints.phaseSegmentFingerprint = "CURRENT-SEGMENTS";
		constraints.currentReauthorizedHistoricalDualRows = {{1, 0}, {0, 1}};
		constraints.currentReauthorizedHistoricalDualValues = {3, 7};
		constraints.currentReauthorizedHistoricalDecisionProofs = {proof};
		constraints.currentReauthorizedHistoricalFamilyIds = {"F"};
		constraints.currentReauthorizedHistoricalSnapshotId = "S";
		constraints.currentReauthorizedSourceMomentId = "M";
		constraints.currentReauthorizedHistoricalFailureProbability = 2e-4;
		constraints.currentReauthorizedHistoricalSegmentsValidated = true;
		constraints.currentReauthorizedHistoricalExactTransport = true;
		ZhangCertifiedPairRelation wide;
		wide.firstNode = 0; wide.secondNode = 1;
		wide.value = 3; wide.coordinate = "WL";
		constraints.certifiedPairs = {wide};
		ZhangCertifiedPairRelation dual;
		dual.firstNode = 0; dual.secondNode = 1;
		dual.value = 7; dual.coordinate = "L1_AND_WL";
		dual.currentPosteriorReauthorized = true;
		constraints.dualFrequencyCertifiedPairs = {dual};
		return constraints;
	};

	auto wrongValue = makeConstraints();
	wrongValue.certifiedPairs.front().value = 4;
	BOOST_CHECK_EQUAL(zhangFinalizeHistoricalProductPairAuthorization(
		wrongValue, "DELIVERY"), 0);
	BOOST_CHECK(!wrongValue.dualFrequencyCertifiedPairs.front()
		.currentPosteriorReauthorized);

	auto missingSegments = makeConstraints();
	missingSegments.phaseSegmentFingerprint.clear();
	BOOST_CHECK_EQUAL(zhangFinalizeHistoricalProductPairAuthorization(
		missingSegments, "DELIVERY"), 0);

	auto excessiveRisk = makeConstraints();
	excessiveRisk.currentReauthorizedHistoricalFailureProbability = 2e-3;
	BOOST_CHECK_EQUAL(zhangFinalizeHistoricalProductPairAuthorization(
		excessiveRisk, "DELIVERY"), 0);

	auto missingMoment = makeConstraints();
	BOOST_CHECK_EQUAL(zhangFinalizeHistoricalProductPairAuthorization(
		missingMoment, ""), 0);
}

BOOST_AUTO_TEST_CASE(r43_candidate_first_failure_taxonomy_is_stable)
{
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::EXACT_AFFINE_CONFLICT),
		"EXACT_AFFINE_CONFLICT");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::CURRENT_MAPPING_FAIL),
		"CURRENT_MAPPING_FAIL");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::NONPRIMITIVE_OR_CONGRUENCE_FAIL),
		"NONPRIMITIVE_OR_CONGRUENCE_FAIL");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::RISK_PARENT_INCOMPLETE),
		"RISK_PARENT_INCOMPLETE");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::RISK_BUDGET_EXCEEDED),
		"RISK_BUDGET_EXCEEDED");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::DETERMINISTIC_RESIDUAL_FAIL),
		"DETERMINISTIC_RESIDUAL_FAIL");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::CURRENT_NIS_FAIL),
		"CURRENT_NIS_FAIL");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::ZERO_PAIR_IMAGE),
		"ZERO_PAIR_IMAGE");
	BOOST_CHECK_EQUAL(zhangProductCandidateFirstFailureName(
		ZhangProductCandidateFirstFailure::SEARCH_PRUNED_NOT_TESTED),
		"SEARCH_PRUNED_NOT_TESTED");
}

BOOST_AUTO_TEST_CASE(kfilter_block_guard_rejects_zero_state_before_dgemv)
{
	FilterChunk fullExtent;
	fullExtent.begX = 2;
	fullExtent.numX = -1;
	fullExtent.begH = 3;
	fullExtent.numH = -1;
	resolveKalmanFilterChunkExtents(fullExtent, 12, 15);
	BOOST_CHECK_EQUAL(fullExtent.numX, 10);
	BOOST_CHECK_EQUAL(fullExtent.numH, 12);

	BOOST_CHECK(!activeKalmanFilterStateChunk(0, 1));
	BOOST_CHECK(!activeKalmanFilterStateChunk(1, 0));
	BOOST_CHECK(!activeKalmanFilterStateChunk(-1, 1));
	BOOST_CHECK(activeKalmanFilterStateChunk(1, 1));

	BOOST_CHECK(!validKalmanFilterBlasBlock(
		0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1));
	BOOST_CHECK(!validKalmanFilterBlasBlock(
		0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0));
	BOOST_CHECK(!validKalmanFilterBlasBlock(
		0, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1));
	BOOST_CHECK(validKalmanFilterBlasBlock(
		0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1));
}

BOOST_AUTO_TEST_CASE(r45_pair_ledger_rhs_requires_exact_affine_parent_pullback)
{
    const ZhangExactMatrix rows = {{1, 0, -1, 0}, {0, 1, 0, -1}};
    const ZhangExactVector physical = {4, 7};
    const ZhangExactVector offsets = {10, 3, 2, 1};
    const ZhangExactVector pair = {1, -1, -1, 1};
    ZhangExactInteger rhs = 999;
    // The physical RHS is 4-7=-3, while the product pair is -3+(8-2)=3.
    BOOST_REQUIRE(zhangCertifiedPairPhysicalRhs(pair, 3, {1, -1}, 0,
        rows, physical, offsets, rhs));
    BOOST_CHECK(rhs == -3);
    BOOST_CHECK(!zhangCertifiedPairPhysicalRhs(pair, -3, {1, -1}, 0,
        rows, physical, offsets, rhs));
    BOOST_CHECK(rhs == -3); // failure never publishes a partial result
    BOOST_CHECK(!zhangCertifiedPairPhysicalRhs(pair, 3, {1, 1}, 0,
        rows, physical, offsets, rhs));
    BOOST_CHECK(!zhangCertifiedPairPhysicalRhs(pair, 3, {1, -1}, 1,
        rows, physical, offsets, rhs));
    BOOST_CHECK(!zhangCertifiedPairPhysicalRhs(pair, 3, {1, -1}, 0,
        rows, physical, {}, rhs));
    BOOST_REQUIRE(zhangCertifiedPairPhysicalRhs(pair, -3, {1, -1}, 0,
        rows, physical, {0, 0, 0, 0}, rhs));
    BOOST_CHECK(rhs == -3);
}

BOOST_AUTO_TEST_CASE(r45_bridge_coverage_preserves_total_and_per_attempt_risk_caps)
{
    constexpr double budget = 1e-3;
    for (std::size_t count : {1u, 4u, 5u, 29u, 206u, 400u})
    {
        const double allocation = zhangConnectivityTargetFailureAllocation(budget, count);
        BOOST_CHECK_GT(allocation, 0);
        BOOST_CHECK_LE(allocation, budget / 32);
        BOOST_CHECK_LE(allocation * count, budget / 8 * (1 + 1e-14));
    }
    BOOST_CHECK_EQUAL(zhangConnectivityTargetFailureAllocation(budget, 0), 0);
    BOOST_CHECK_EQUAL(zhangConnectivityTargetFailureAllocation(-budget, 20), 0);
}

BOOST_AUTO_TEST_CASE(r45_mature_conflict_records_exact_scope_without_committing_fresh_rows)
{
    ProductIntegerLedger ledger;
    ProductIntegerLedgerRow row;
    row.system = E_Sys::GPS;
    row.firstObservable = E_ObsCode::L1C;
    row.secondObservable = E_ObsCode::L2W;
    row.productRow = {1, -1};
    row.physicalExpansion = {{"L1C|A|G01|V1", 1}, {"L1C|A|G02|V1", -1}};
    row.canonicalProductExpansion = {{"L1C|G01", 1}, {"L1C|G02", -1}};
    row.phaseSegmentFingerprint = "G01|L1C|SEG0;G02|L1C|SEG0;";
    row.integerValue = 4;
    row.backendBasisGeneration = 7;
    BOOST_REQUIRE(ledger.observe(100, {row}, 1).valid);
    auto conflict = row;
    conflict.integerValue = 5;
    conflict.backendBasisGeneration = 8;
    auto fresh = row;
    fresh.canonicalProductExpansion = {{"L1C|G03", 1}, {"L1C|G04", -1}};
    fresh.phaseSegmentFingerprint = "G03|L1C|SEG0;G04|L1C|SEG0;";
    const auto rejected = ledger.observe(130, {fresh, conflict}, 1);
    BOOST_CHECK(!zhangProductLedgerWriterCommitAuthorized(rejected));
    BOOST_CHECK_EQUAL(rejected.conflictingCandidateIndex, 1);
    BOOST_CHECK(rejected.conflictPhysicalRowsEqual);
    BOOST_CHECK(rejected.incumbentInteger == 4);
    BOOST_CHECK(rejected.proposedInteger == 5);
    BOOST_CHECK_EQUAL(rejected.incumbentBackendGeneration, 7);
    BOOST_CHECK_EQUAL(rejected.proposedBackendGeneration, 8);
    BOOST_CHECK(!rejected.conflictIdentity.empty());
    BOOST_CHECK_EQUAL(ledger.rows().size(), 1);
    BOOST_CHECK(ledger.rows().front().integerValue == 4);
}
