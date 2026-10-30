// Link this driver against production PEA objects (rename main.cpp.o's main,
// as in run_r51_algebra_transport.py). No KFState stubs or replayed core oracle.
#include "common/algebra.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
void require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

int rejectionCalls = 0;
int committedEffects = 0;

bool rejectWithoutRepair(RejectCallbackDetails detail)
{
    ++rejectionCalls;
    require(detail.kfState.activeMeasurementTransaction != nullptr,
        "rejection escaped the authoritative transaction");
    detail.kfState.activeMeasurementTransaction->commitSideEffects.push_back(
        [] { ++committedEffects; });
    return true;
}

bool repairLastMeasurement(RejectCallbackDetails detail)
{
    rejectWithoutRepair(detail);
    detail.kfMeas.R(detail.measIndex, detail.measIndex) = 1e6;
    return true;
}

bool rejectAndRelaxPrior(RejectCallbackDetails detail)
{
    rejectWithoutRepair(detail);
    auto& state = detail.kfState;
    auto& transaction = *state.activeMeasurementTransaction;
    KFStagedFactor transition;
    transition.kind = KFStagedFactorKind::STATE_TRANSITION;
    transition.label = "TEST_ADAPTIVE_PRIOR";
    transition.adaptiveQcTransition = true;
    transition.sourceIndexMap = transition.destinationIndexMap = state.kfIndexMap;
    transition.sourceMean = transition.destinationMean = state.x;
    transition.sourceCovariance = state.P;
    transition.transition = SparseMatrix<double>(1, 1);
    transition.transition.insert(0, 0) = 1;
    transition.processCovariance = SparseMatrix<double>(1, 1);
    transition.processCovariance.insert(0, 0) = 1;
    state.P(0, 0) += 1;
    transition.destinationCovariance = state.P;
    transition.beforeCommitSequence = transaction.provisionalCommitSequence;
    transition.afterCommitSequence = ++transaction.provisionalCommitSequence;
    transaction.stagedFactors.push_back(std::move(transition));
    return true;
}

bool repairToNonfinite(RejectCallbackDetails detail)
{
    rejectWithoutRepair(detail);
    detail.kfMeas.R(detail.measIndex, detail.measIndex) =
        std::numeric_limits<double>::infinity();
    return true;
}

KFState scalarState()
{
    KFKey key; key.type = KF::IONO_STEC; key.str = "QC_SCALAR";
    KFState state;
    state.kfIndexMap = {{key, 0}};
    state.stateTransitionMap.clear();
    state.stateTransitionMap[key][key][0] = 1;
    state.initNoiseMap.clear();
    state.x = VectorXd::Zero(1);
    state.dx = VectorXd::Zero(1);
    state.P = MatrixXd::Identity(1, 1);
    state.prefitOpts.sigma_check = false;
    state.prefitOpts.omega_test = false;
    state.postfitOpts.sigma_check = false;
    state.postfitOpts.omega_test = true;
    state.postfitOpts.state_sigma_threshold = 8;
    state.postfitOpts.meas_sigma_threshold = 6;
    state.postfitOpts.max_iterations = 20;
    state.factorCommitSequence = 7;
    state.metaDataMap["QC_SENTINEL"] = "current-epoch-entry";
    return state;
}

KFMeas scalarMeasurement(const KFState& state, double observed)
{
    KFMeas row;
    row.Y = VectorXd::Constant(1, observed);
    row.V = row.Y;
    row.VV = row.V;
    row.H = MatrixXd::Identity(1, state.x.size());
    row.R = MatrixXd::Identity(1, 1);
    row.H_star = MatrixXd::Zero(1, 0);
    row.uncorrelatedNoise = VectorXd::Zero(0);
    row.prefitRatios = VectorXd::Zero(1);
    row.postfitRatios = VectorXd::Zero(1);
    row.obsKeys = {state.kfIndexMap.begin()->first};
    row.metaDataMaps.resize(1);
    row.componentsMaps.resize(1);
    return row;
}

void checkUnrepairedRollback(bool relaxPrior)
{
    auto state = scalarState();
    auto row = scalarMeasurement(state, 100);
    const auto entryState = state;
    const auto entryRow = row;
    rejectionCalls = committedEffects = 0;
    state.measRejectCallbacks.push_back(relaxPrior ? rejectAndRelaxPrior : rejectWithoutRepair);
    bool factorCalled = false;
    state.authoritativeFactorTransactionCallback =
        [&](const KFState&, const KFMeasurementTransaction&, std::string&)
        { factorCalled = true; return true; };
    std::ostringstream trace;
    const auto result = state.filterKalman(trace, row, "/PPP", true);
    require(result == KFFilterResult::QC_NOT_CONVERGED,
        "20 failed robust tests committed an uncertified posterior");
    require(rejectionCalls == 20, "terminal QC invoked a 21st rejection callback");
    require(committedEffects == 0 && !factorCalled,
        "failed QC published deferred effects or the factor packet");
    require((state.x - entryState.x).norm() == 0 &&
        (state.P - entryState.P).norm() == 0 &&
        state.kfIndexMap == entryState.kfIndexMap &&
        state.stateTransitionMap == entryState.stateTransitionMap,
        "QC failure reset or changed the current transaction entry");
    require((row.R - entryRow.R).norm() == 0 && (row.VV - entryRow.VV).norm() == 0,
        "QC rollback left a rejected measurement model behind");
    require(state.factorCommitSequence == 7 && !state.lastPppTransactionCommitted &&
        !state.lastPppQcConverged && state.activeMeasurementTransaction == nullptr,
        "QC failure changed sequence or retained transaction authority");
    require(state.metaDataMap.at("QC_SENTINEL") == "current-epoch-entry",
        "QC rejection reset state metadata");
    require(state.lastFactorTransactionFailureReason == "FINAL_POSTFIT_QC_NOT_CONVERGED",
        "QC failure lost its reason");
    require(trace.str().find("status=COMMITTED") == std::string::npos &&
        trace.str().find("converged=0 reject_callbacks_enabled=0") != std::string::npos,
        "failed final QC was logged as a committed measurement");
}

void checkLastRepair(bool rejectFactor)
{
    auto state = scalarState();
    auto row = scalarMeasurement(state, 100);
    // A one-iteration budget forces the last callback to change R after the
    // last tested solve. The production 20-iteration policy is not increased.
    state.postfitOpts.max_iterations = 1;
    state.measRejectCallbacks.push_back(repairLastMeasurement);
    rejectionCalls = committedEffects = 0;
    int packets = 0;
    state.authoritativeFactorTransactionCallback =
        [&](const KFState&, const KFMeasurementTransaction& packet, std::string& reason)
        {
            ++packets;
            require(packet.qcConverged && packet.maxIterationsReached &&
                packet.terminalReconciliationPerformed &&
                packet.modelGeneration == packet.solveGeneration,
                "final passing QC did not certify the final model");
            if (rejectFactor) reason = "TEST_FACTOR_REJECTED";
            return !rejectFactor;
        };
    std::ostringstream trace;
    const auto result = state.filterKalman(trace, row, "/PPP", true);
    require(rejectionCalls == 1 && packets == 1, "terminal QC repeated side effects");
    if (rejectFactor)
    {
        require(result == KFFilterResult::FACTOR_TRANSACTION_REJECTED &&
            state.x(0) == 0 && state.P(0, 0) == 1 && row.R(0, 0) == 1 &&
            state.factorCommitSequence == 7 && committedEffects == 0,
            "factor rejection after final QC did not atomically roll back");
        return;
    }
    // Independent scalar Kalman equations: K=P/(P+R)=1/1000001.
    require(result == KFFilterResult::COMMITTED && state.lastPppQcConverged,
        "a genuine passing final QC was rejected at the iteration cap");
    require(std::abs(state.x(0) - 100.0 / 1000001.0) < 1e-12 &&
        std::abs(state.P(0, 0) - 1000000.0 / 1000001.0) < 1e-12,
        "committed scalar posterior does not use the final noise model");
    require(state.factorCommitSequence == 8 && committedEffects == 1,
        "passing transaction published sequence/effects more than once");
    require(trace.str().find("statistics_valid=1 converged=1 reject_callbacks_enabled=0") !=
        std::string::npos, "final passing posterior was not audited");
}

void checkNormalEpochs()
{
    KFKey observed; observed.type = KF::IONO_STEC; observed.str = "OBSERVED";
    KFKey unseen = observed; unseen.str = "UNOBSERVED";
    auto state = scalarState();
    state.kfIndexMap = {{KFState::oneKey, 0}, {observed, 1}, {unseen, 2}};
    state.stateTransitionMap.clear();
    for (const auto& [key, index] : state.kfIndexMap)
        state.stateTransitionMap[key][key][0] = 1;
    state.x = VectorXd::Zero(3); state.x(0) = 1;
    state.dx = VectorXd::Zero(3);
    state.P = MatrixXd::Zero(3, 3); state.P(1, 1) = state.P(2, 2) = 1;
    for (int epoch = 0; epoch < 3; ++epoch)
    {
        auto row = scalarMeasurement(state, 1);
        row.H = MatrixXd::Zero(1, 3); row.H(0, 1) = 1;
        row.obsKeys = {observed};
        row.V = row.Y - row.H * state.x; row.VV = row.V;
        std::ostringstream trace;
        require(state.filterKalman(trace, row, "/PPP", true) == KFFilterResult::COMMITTED,
            "normal epoch rejected fixed ONE or an unobserved state 0/0");
        // n unit-variance observations of 1, prior x=0,P=1: x=n/(n+1),P=1/(n+1).
        const double n = epoch + 1;
        require(std::abs(state.x(1) - n / (n + 1)) < 1e-12 &&
            std::abs(state.P(1, 1) - 1 / (n + 1)) < 1e-12 &&
            state.x(0) == 1 && state.P(0, 0) == 0 &&
            state.x(2) == 0 && state.P(2, 2) == 1 && state.lastPppQcConverged,
            "normal epoch posterior differs from analytic scalar sequence");
    }
}

void checkInvalidStatistics()
{
    auto state = scalarState();
    auto row = scalarMeasurement(state, 1);
    state.postfitRatios = VectorXd::Constant(1, 123);
    row.postfitRatios(0) = 456;
    VectorXd dx = VectorXd::Constant(1, 0.5);
    for (double diagonal : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::quiet_NaN()})
    {
        MatrixXd inverse = MatrixXd::Constant(1, 1, diagonal);
        MatrixXd inverseH = MatrixXd::Constant(1, 1, 0.5);
        KFStatistics stats;
        std::ostringstream trace;
        RejectCallbackDetails detail(trace, state, row);
        require(!state.postFitSigmaChecks(detail, dx, inverse, inverseH, stats, 0, 1, 0, 1, false),
            "invalid omega denominator was silently replaced by a passing zero");
        require(state.postfitRatios(0) == 123 && row.postfitRatios(0) == 456,
            "read-only final QC published diagnostic side effects");
    }
    MatrixXd inverse = MatrixXd::Constant(1, 1, 0.5);
    MatrixXd inverseH = MatrixXd::Zero(1, 1);
    KFStatistics stats;
    std::ostringstream trace;
    RejectCallbackDetails detail(trace, state, row);
    require(!state.postFitSigmaChecks(detail, dx, inverse, inverseH, stats, 0, 1, 0, 1, false),
        "an observed stochastic state with degenerate omega denominator passed");
    require(!state.postFitSigmaChecks(detail, dx, inverse, inverseH, stats, 0, 1, 0, 0, false),
        "an empty QC block passed");
    inverseH(0, 0) = 0.5;
    require(state.postFitSigmaChecks(detail, dx, inverse, inverseH, stats, 0, 1, 0, 1, false) &&
        std::abs(stats.sumOfSquares - 0.5) < 1e-12 &&
        state.postfitRatios(0) == 123 && row.postfitRatios(0) == 456,
        "valid read-only QC changed arrays or disagreed with the analytic omega statistic");
    row.V(0) = std::numeric_limits<double>::quiet_NaN();
    require(!state.postFitSigmaChecks(detail, dx, inverse, inverseH, stats, 0, 1, 0, 1, false),
        "nonfinite innovations passed QC");
}

void checkNonfiniteTransaction()
{
    auto state = scalarState();
    auto row = scalarMeasurement(state, std::numeric_limits<double>::infinity());
    std::ostringstream trace;
    require(state.filterKalman(trace, row, "/PPP", true) == KFFilterResult::NUMERICAL_FAILURE,
        "nonfinite transaction did not fail closed");
    require(state.x(0) == 0 && state.P(0, 0) == 1 && state.factorCommitSequence == 7 &&
        !state.lastPppTransactionCommitted && state.activeMeasurementTransaction == nullptr,
        "nonfinite transaction changed the prior or escaped rollback");

    state = scalarState();
    row = scalarMeasurement(state, 100);
    state.postfitOpts.max_iterations = 1;
    state.measRejectCallbacks.push_back(repairToNonfinite);
    rejectionCalls = committedEffects = 0;
    require(state.filterKalman(trace, row, "/PPP", true) == KFFilterResult::NUMERICAL_FAILURE,
        "nonfinite final model was accepted after a valid first solve");
    require(state.x(0) == 0 && state.P(0, 0) == 1 && row.R(0, 0) == 1 &&
        state.factorCommitSequence == 7 && committedEffects == 0 && rejectionCalls == 1,
        "nonfinite terminal solve leaked the final callback's model/effects");
}
}

int main()
{
    checkUnrepairedRollback(false);
    checkUnrepairedRollback(true);
    checkLastRepair(false);
    checkLastRepair(true);
    checkNormalEpochs();
    checkInvalidStatistics();
    checkNonfiniteTransaction();
    std::cout << "PASS production FLOAT QC: terminal certification, transactional rollback, "
        "deferred effects, final-model analytic posterior, fixed/unobserved states, "
        "three normal epochs, invalid statistics and nonfinite inputs\n";
}
