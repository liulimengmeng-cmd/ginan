// Link this driver to the genuine PEA objects, with only main.cpp.o's main
// renamed in a COPY. Run the identical driver against frozen v2 and candidate
// objects and compare its complete bit-pattern output, not merely its exit.
// The reference below never calls stateTransition()/another predictor: F/Q
// come from closed equations, with a separate scalar-loop covariance oracle.
#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/receiver.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <set>
#include <stdexcept>
#include <tuple>

namespace
{
std::size_t checks = 0, cases = 0;
void require(bool condition, const string& reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}

bool sameBits(const MatrixXd& left, const MatrixXd& right)
{
    return left.rows() == right.rows() && left.cols() == right.cols() &&
        (left.size() == 0 || std::memcmp(left.data(), right.data(), left.size()*sizeof(double)) == 0);
}

auto options(const KFState& state)
{
    return std::make_tuple(state.rts_filename, state.rts_directory, state.rts_lag,
        state.rts_interval, state.rts_smoothed_suffix, state.queue_rts_outputs,
        state.rts_regularisation, state.simulate_filter_only, state.assume_linearity,
        state.advanced_postfits, state.joseph_stabilisation, state.lsq_inverter,
        state.lsqOpts.max_iterations, state.lsqOpts.sigma_check, state.lsqOpts.omega_test,
        state.lsqOpts.meas_sigma_threshold, state.prefitOpts.max_iterations,
        state.prefitOpts.sigma_check, state.prefitOpts.omega_test,
        state.prefitOpts.state_sigma_threshold, state.prefitOpts.meas_sigma_threshold,
        state.postfitOpts.max_iterations, state.postfitOpts.sigma_check,
        state.postfitOpts.omega_test, state.postfitOpts.state_sigma_threshold,
        state.postfitOpts.meas_sigma_threshold, state.chiSquareTest.enable,
        state.chiSquareTest.mode, state.chiSquareTest.sigma_threshold);
}

void checkMetadataAndMoments(const KFState& actual, const KFState& expected)
{
    const auto actualCore=captureZhangCheckpointKfCore(actual);
    const auto expectedCore=captureZhangCheckpointKfCore(expected);
    const bool completeEqual=zhangCheckpointKfCoreBitwiseEqual(actualCore,expectedCore);
    if (!completeEqual)
    {
        const auto archive=[](const auto& value)
        {
            std::ostringstream output(std::ios::binary|std::ios::out);
            boost::archive::binary_oarchive ar(output,boost::archive::no_header);
            ar << value;
            return output.str();
        };
        if (archive(actualCore.time)!=archive(expectedCore.time)) std::cerr<<"CORE_FIELD_DIFF time\n";
        if (archive(actualCore.x)!=archive(expectedCore.x)) std::cerr<<"CORE_FIELD_DIFF x\n";
        if (archive(actualCore.P)!=archive(expectedCore.P)) std::cerr<<"CORE_FIELD_DIFF P\n";
        if (archive(actualCore.dx)!=archive(expectedCore.dx)) std::cerr<<"CORE_FIELD_DIFF dx\n";
        if (archive(actualCore.prefitRatios)!=archive(expectedCore.prefitRatios)) std::cerr<<"CORE_FIELD_DIFF prefitRatios\n";
        if (archive(actualCore.postfitRatios)!=archive(expectedCore.postfitRatios)) std::cerr<<"CORE_FIELD_DIFF postfitRatios\n";
        if (archive(actualCore.kfIndexMap)!=archive(expectedCore.kfIndexMap)) std::cerr<<"CORE_FIELD_DIFF kfIndexMap\n";
        if (archive(actualCore.stateTransitionMap)!=archive(expectedCore.stateTransitionMap)) std::cerr<<"CORE_FIELD_DIFF stateTransitionMap\n";
        if (archive(actualCore.gaussMarkovTauMap)!=archive(expectedCore.gaussMarkovTauMap)) std::cerr<<"CORE_FIELD_DIFF gaussMarkovTauMap\n";
        if (archive(actualCore.gaussMarkovMuMap)!=archive(expectedCore.gaussMarkovMuMap)) std::cerr<<"CORE_FIELD_DIFF gaussMarkovMuMap\n";
        if (archive(actualCore.procNoiseMap)!=archive(expectedCore.procNoiseMap)) std::cerr<<"CORE_FIELD_DIFF procNoiseMap\n";
        if (archive(actualCore.initNoiseMap)!=archive(expectedCore.initNoiseMap)) std::cerr<<"CORE_FIELD_DIFF initNoiseMap\n";
        if (archive(actualCore.sigmaMaxMap)!=archive(expectedCore.sigmaMaxMap)) std::cerr<<"CORE_FIELD_DIFF sigmaMaxMap\n";
        if (archive(actualCore.outageLimitMap)!=archive(expectedCore.outageLimitMap)) std::cerr<<"CORE_FIELD_DIFF outageLimitMap\n";
        if (archive(actualCore.exponentialNoiseMap)!=archive(expectedCore.exponentialNoiseMap)) std::cerr<<"CORE_FIELD_DIFF exponentialNoiseMap\n";
        if (archive(actualCore.pseudoStateMap)!=archive(expectedCore.pseudoStateMap)) std::cerr<<"CORE_FIELD_DIFF pseudoStateMap\n";
        if (archive(actualCore.pseudoParentMap)!=archive(expectedCore.pseudoParentMap)) std::cerr<<"CORE_FIELD_DIFF pseudoParentMap\n";
        if (archive(actualCore.errorCountMap)!=archive(expectedCore.errorCountMap)) std::cerr<<"CORE_FIELD_DIFF errorCountMap\n";
        if (archive(actualCore.filterChunkMap)!=archive(expectedCore.filterChunkMap)) std::cerr<<"CORE_FIELD_DIFF filterChunkMap\n";
        if (archive(actualCore.metaDataMap)!=archive(expectedCore.metaDataMap)) std::cerr<<"CORE_FIELD_DIFF metaDataMap\n";
        if (archive(actualCore.lsqRequired)!=archive(expectedCore.lsqRequired)) std::cerr<<"CORE_FIELD_DIFF lsqRequired\n";
        if (archive(actualCore.sigmaPass)!=archive(expectedCore.sigmaPass)) std::cerr<<"CORE_FIELD_DIFF sigmaPass\n";
        if (archive(actualCore.chiQCPass)!=archive(expectedCore.chiQCPass)) std::cerr<<"CORE_FIELD_DIFF chiQCPass\n";
        if (archive(actualCore.chi2)!=archive(expectedCore.chi2)) std::cerr<<"CORE_FIELD_DIFF chi2\n";
        if (archive(actualCore.dof)!=archive(expectedCore.dof)) std::cerr<<"CORE_FIELD_DIFF dof\n";
        if (archive(actualCore.chi2PerDof)!=archive(expectedCore.chi2PerDof)) std::cerr<<"CORE_FIELD_DIFF chi2PerDof\n";
        if (archive(actualCore.qc)!=archive(expectedCore.qc)) std::cerr<<"CORE_FIELD_DIFF qc\n";
        if (archive(actualCore.id)!=archive(expectedCore.id)) std::cerr<<"CORE_FIELD_DIFF id\n";
        if (archive(actualCore.rtsBasename)!=archive(expectedCore.rtsBasename)) std::cerr<<"CORE_FIELD_DIFF rtsBasename\n";
        if (archive(actualCore.outputResiduals)!=archive(expectedCore.outputResiduals)) std::cerr<<"CORE_FIELD_DIFF outputResiduals\n";
        if (archive(actualCore.outputMongoMeasurements)!=archive(expectedCore.outputMongoMeasurements)) std::cerr<<"CORE_FIELD_DIFF outputMongoMeasurements\n";
        if (archive(actualCore.statisticsMap)!=archive(expectedCore.statisticsMap)) std::cerr<<"CORE_FIELD_DIFF statisticsMap\n";
        if (archive(actualCore.statisticsMapSum)!=archive(expectedCore.statisticsMapSum)) std::cerr<<"CORE_FIELD_DIFF statisticsMapSum\n";
    }
    require(completeEqual, "Complete checkpoint core/metadata changed");
    require(options(actual) == options(expected), "Filter/RTS option metadata changed");
    require(actual.factorCommitSequence == expected.factorCommitSequence &&
        actual.activeMeasurementTransaction == expected.activeMeasurementTransaction &&
        actual.factorCallbackOwnerOverride == expected.factorCallbackOwnerOverride &&
        actual.alternate_ptr == expected.alternate_ptr &&
        actual.lastPppTransactionCommitted == expected.lastPppTransactionCommitted &&
        actual.lastPppQcConverged == expected.lastPppQcConverged,
        "Runtime/transaction identity or commit metadata changed");
    require(actual.stateRejectCallbacks.size() == expected.stateRejectCallbacks.size() &&
        actual.measRejectCallbacks.size() == expected.measRejectCallbacks.size() &&
        bool(actual.acceptedMeasurementFactorCallback) == bool(expected.acceptedMeasurementFactorCallback) &&
        bool(actual.stateTransitionFactorCallback) == bool(expected.stateTransitionFactorCallback) &&
        bool(actual.exactStateTransformCallback) == bool(expected.exactStateTransformCallback) &&
        bool(actual.authoritativeFactorTransactionCallback) == bool(expected.authoritativeFactorTransactionCallback),
        "Callback metadata changed");
    require(actual.filterChunkMap.size() == expected.filterChunkMap.size(), "Chunk count changed");
    for (const auto& [name, chunk] : expected.filterChunkMap)
        require(actual.filterChunkMap.at(name).trace_ptr == chunk.trace_ptr, "Chunk trace pointer changed");
    require(actual.time.cacheN == expected.time.cacheN &&
        actual.time.cacheTime == expected.time.cacheTime &&
        actual.time.cacheString == expected.time.cacheString, "GTime cache was changed by profiling/rollback");
}

struct Fixture
{
    KFState state;
    KFKey one = KFState::oneKey;
    KFKey position{.type=KF::REC_POS, .str="A", .num=0, .comment="position"};
    KFKey rate{.type=KF::REC_POS, .str="A", .num=1, .comment="rate"};
    KFKey clock{.type=KF::REC_CLOCK, .str="B", .comment="clock"};
    KFKey old{.type=KF::IONO_STEC, .str="Z", .num=0, .comment="retired"};
    KFKey added{.type=KF::IONO_STEC, .str="Z", .num=1, .comment="added"};
    MatrixXd F, Q;
    SparseMatrix<double> sparseF;
    VectorXd mean, adjustment, oldEigenMean;
    MatrixXd oldEigenCovariance, scalarCovariance;
    map<KFKey, int> destination;
    GTime target;
    double elapsed = 0;
};

Fixture fixture(double dt, bool grow, bool retire, bool gaussMarkovRate,
                bool simulate, bool zeroInitialRate, bool globalState)
{
    Fixture out;
    auto& state = out.state;
    if (globalState)
    {
        out.clock.str.clear();
        out.clock.Sat = SatSys(E_Sys::GPS, 12);
    }
    state.time.bigTime = 1405213200.L;
    state.time.cacheN = 3;
    state.time.cacheTime = -77;
    state.time.cacheString = "UNCHANGED_TIME_CACHE";
    // This prediction-equivalence fixture is explicitly NOT an outage event.
    // Set before inserting every copy of this key into the model/state maps.
    out.clock.estimatedTime = state.time;
    out.target.bigTime = state.time.bigTime + dt;
    out.elapsed = std::abs(dt);
    state.kfIndexMap = {{out.one,0}, {out.position,1}, {out.rate,2}, {out.clock,3}, {out.old,4}};
    state.x.resize(5);
    state.x << 1., 2., 3., -1., .5;
    state.dx = VectorXd::Constant(5, -7.);
    state.P = MatrixXd::Zero(5,5);
    for (int row=1; row<5; ++row)
        for (int column=1; column<5; ++column)
            state.P(row,column) = row == column ? 5.+row : .125*(row+column);
    state.prefitRatios = VectorXd::Constant(3, .25);
    state.postfitRatios = VectorXd::Constant(4, .75);
    state.stateTransitionMap.clear();
    state.stateTransitionMap[out.one][out.one][0] = 1;
    state.stateTransitionMap[out.position][out.position][0] = 1;
    state.stateTransitionMap[out.position][out.rate][1] = 1;
    state.stateTransitionMap[out.rate][out.rate][0] = 1;
    state.stateTransitionMap[out.clock][out.clock][0] = 1;
    if (!retire) state.stateTransitionMap[out.old][out.old][0] = 1;
    if (grow) state.stateTransitionMap[out.added][out.one][0] = 9;
    state.initNoiseMap = {{out.position,.25}, {out.clock,.5}};
    if (grow) state.initNoiseMap[out.added] = 9;
    if (zeroInitialRate) state.initNoiseMap[out.rate] = 0;
    state.procNoiseMap = {{out.position,.1}, {out.rate,.4}, {out.clock,.2}};
    state.gaussMarkovTauMap = {{out.position,-1}, {out.rate,gaussMarkovRate?60.:-1.}, {out.clock,60}};
    state.gaussMarkovMuMap = {{out.clock,5}};
    if (gaussMarkovRate) state.gaussMarkovMuMap[out.rate] = 2;
    state.exponentialNoiseMap[out.clock] = {.value=.3, .tau=15};
    state.sigmaMaxMap[out.position] = 1e20;
    state.outageLimitMap[out.clock] = 60; // Above the largest tested +30 second gap.
    state.errorCountMap[out.rate] = 7;
    state.pseudoStateMap[out.old][out.added] = 2;
    state.pseudoParentMap[out.added] = out.old;
    state.metaDataMap = {{"runtime","prediction-owner"}, {"sentinel","metadata-is-not-rebuilt"}};
    state.filterChunkMap["sentinel"] = {.id="sentinel", .begX=1, .numX=2, .begH=3, .numH=4};
    state.statisticsMap["epoch"] = 17;
    state.statisticsMapSum["epoch"] = 29;
    state.factorCommitSequence = 42;
    state.lastPppTransactionCommitted = true;
    state.lastPppQcConverged = true;
    state.lsqRequired = true;
    state.sigmaPass = true;
    state.chiQCPass = false;
    state.chi2 = 3.5;
    state.dof = 7;
    state.chi2PerDof = .5;
    state.qc = 2.25;
    state.id = "TRANSITION_OWNER";
    state.output_residuals = true;
    state.outputMongoMeasurements = true;
    state.rts_filename = "DO_NOT_PUBLISH_RTS";
    state.rts_directory = "OPTION_METADATA";
    state.rts_lag = 11;
    state.rts_interval = 3;
    state.rts_regularisation = 2e-11;
    state.simulate_filter_only = simulate;
    state.advanced_postfits = true;
    state.assume_linearity = true;
    state.joseph_stabilisation = true;
    state.prefitOpts.max_iterations = 7;
    state.postfitOpts.max_iterations = 13;
    state.chiSquareTest.enable = false; // No diagnostic NIS is introduced.
    state.lastFactorTransactionFailureReason = "ENTRY_SENTINEL";

    int row=0;
    for (const auto& [key, unused] : state.stateTransitionMap) out.destination[key] = row++;
    out.F = MatrixXd::Zero(row,5);
    out.Q = MatrixXd::Zero(row,row);
    out.sparseF.resize(row,5);
    auto coefficient = [&](const KFKey& key, int column, double value)
    {
        const int index = out.destination.at(key);
        out.F(index,column) += value;
        // Keep explicit zeros too: the legacy sparse construction inserts them.
        out.sparseF.coeffRef(index,column) += value;
    };
    coefficient(out.one,0,1);
    coefficient(out.position,1,1);
    if (gaussMarkovRate)
    {
        const double alpha=std::exp(-dt/60.);
        const double integral=60.*(1.-alpha);
        coefficient(out.position,2,integral);
        coefficient(out.position,0,2.*(dt-integral));
        coefficient(out.rate,2,alpha);
        coefficient(out.rate,0,2.*(1.-alpha));
    }
    else
    {
        coefficient(out.position,2,dt);
        coefficient(out.rate,2,1);
    }
    const double clockAlpha=std::exp(-dt/60.);
    coefficient(out.clock,3,clockAlpha);
    coefficient(out.clock,0,5.*(1.-clockAlpha));
    if (!retire) coefficient(out.old,4,1);
    if (grow) coefficient(out.added,0,9);
    const double t=out.elapsed;
    auto& qpos=out.Q(out.destination.at(out.position),out.destination.at(out.position));
    auto& qrate=out.Q(out.destination.at(out.rate),out.destination.at(out.rate));
    auto& qclock=out.Q(out.destination.at(out.clock),out.destination.at(out.clock));
    qpos=.25;
    qclock=.5;
    if (grow) out.Q(out.destination.at(out.added),out.destination.at(out.added))=9;
    if (t)
    {
        qclock += .3*t;
        qpos += .1/1*t;
        if (gaussMarkovRate)
        {
            qpos += .4/2*60.*60.*(2*t-4*60.*(1-std::exp(-t/60.))+60.*(1-std::exp(-2*t/60.)));
            if (!zeroInitialRate) qrate += .4/2*60.*(1-std::exp(-2*t/60.));
        }
        else if (!zeroInitialRate) qrate += .4/1*t;
        qclock += .2/2*60.*(1-std::exp(-2*t/60.));
    }
    if (simulate) out.Q.setZero();
    out.mean=VectorXd::Zero(row);
    MatrixXd intermediate=MatrixXd::Zero(row,5);
    for (int i=0; i<row; ++i)
        for (int k=0; k<5; ++k)
        {
            out.mean(i) += out.F(i,k)*state.x(k);
            for (int j=0; j<5; ++j) intermediate(i,j) += out.F(i,k)*state.P(k,j);
        }
    out.scalarCovariance=MatrixXd::Zero(row,row);
    for (int i=0; i<row; ++i)
        for (int j=0; j<row; ++j)
        {
            for (int k=0; k<5; ++k) out.scalarCovariance(i,j) += intermediate(i,k)*out.F(j,k);
            out.scalarCovariance(i,j) += out.Q(i,j);
        }
    // Mandatory bitwise old Eigen expression, separate from production code.
    out.oldEigenCovariance=(out.sparseF*state.P*out.sparseF.transpose()+out.Q).eval();
    out.oldEigenMean=out.sparseF*state.x;
    out.adjustment=VectorXd::Zero(row);
    if (row==5) out.adjustment=out.oldEigenMean-state.x;
    return out;
}

// Only for regression of the already-existing receiver_chunking option. This
// is NOT enabled in the joint-network experiment and is not a new reduction
// of cross-station covariance. It covers that branch and new-source-only rows.
void legacyChunkReference(Fixture& test)
{
    struct Chunk { vector<int> rows; std::set<int> columns; };
    map<string,Chunk> chunks;
    for (const auto& [key, transitions] : test.state.stateTransitionMap)
    {
        if (key==test.one) continue;
        auto& chunk=chunks[key.str];
        chunk.rows.push_back(test.destination.at(key));
        for (const auto& [source, coefficients] : transitions)
        {
            if (source==test.one) continue;
            const auto old=test.state.kfIndexMap.find(source);
            if (old!=test.state.kfIndexMap.end()) chunk.columns.insert(old->second);
        }
    }
    test.oldEigenCovariance.setZero();
    for (const auto& [name, chunk] : chunks)
    {
        const vector<int> columns(chunk.columns.begin(),chunk.columns.end());
        if (columns.empty()) test.oldEigenCovariance(chunk.rows,chunk.rows)=test.Q(chunk.rows,chunk.rows);
        else
        {
            const MatrixXd Fblock=test.F(chunk.rows,columns);
            test.oldEigenCovariance(chunk.rows,chunk.rows)=
                (Fblock*test.state.P(columns,columns)*Fblock.transpose()+test.Q(chunk.rows,chunk.rows)).eval();
        }
    }
    for (const auto& [a, i] : test.destination)
        for (const auto& [b, j] : test.destination)
            if (a==test.one || b==test.one || a.str!=b.str) test.scalarCovariance(i,j)=0;
}

void outputBits(const string& name, const KFState& state)
{
    std::cout << "PREDICTION_BITS " << name << " x=";
    const auto emit = [](double value)
    {
        std::uint64_t bits;
        std::memcpy(&bits,&value,sizeof(bits));
        std::cout << std::hex << std::setw(16) << std::setfill('0') << bits;
    };
    for (int i=0; i<state.x.size(); ++i) emit(state.x(i));
    std::cout << " dx=";
    for (int i=0; i<state.dx.size(); ++i) emit(state.dx(i));
    std::cout << " P=";
    for (int j=0; j<state.P.cols(); ++j)
        for (int i=0; i<state.P.rows(); ++i) emit(state.P(i,j));
    std::cout << std::dec << " core_sha="
        << zhangCheckpointKfCoreSha256Streaming(captureZhangCheckpointKfCore(state)) << "\n";
}

void runCase(const string& name, double dt, bool grow, bool retire, bool gmRate,
             bool simulate, bool zeroInitialRate, bool globalState, int mode, bool requestQ)
{
    auto test=fixture(dt,grow,retire,gmRate,simulate,zeroInitialRate,globalState);
    const bool legacyChunked=acsConfig.pppOpts.receiver_chunking && !globalState;
    if (legacyChunked) legacyChunkReference(test);
    auto& state=test.state;
    KFMeasurementTransaction transaction;
    KFState alternate;
    state.activeMeasurementTransaction=&transaction;
    state.factorCallbackOwnerOverride=&state;
    state.alternate_ptr=&alternate;
    std::ostringstream trace;
    state.filterChunkMap.at("sentinel").trace_ptr=&trace;
    const bool expectedSuccess=mode<=1;
    int callbackCalls=0;
    if (mode)
        state.stateTransitionFactorCallback=[&](const KFState& source, GTime time,
            const map<KFKey,int>& oldKeys,const map<KFKey,int>& newKeys,
            const SparseMatrix<double>& F,const MatrixXd& Q,const string& label,
            const VectorXd& mean,const MatrixXd& covariance,std::uint64_t before,std::uint64_t after)
        {
            ++callbackCalls;
            require(&source==&state && source.P.rows()==5 && source.kfIndexMap==oldKeys,
                    "Factor callback source ownership/moments changed");
            require(source.x(1)==2 && source.P(1,3)==.5, "Callback did not retain full source cross-station P");
            require(time==test.target && newKeys==test.destination && label=="KF_STATE_TRANSITION" && before==42 && after==43,
                    "Factor callback metadata/order changed");
            require(sameBits(MatrixXd(F),test.F) && sameBits(Q,test.Q), "Factor tap F/Q changed");
            require(sameBits(covariance,test.oldEigenCovariance), "Tap covariance differs bitwise from old Eigen expression");
            require(sameBits(mean,test.oldEigenMean) && (mean-test.mean).norm()<1e-11,
                    "Tap mean differs from independent/old literal reference");
            if (mode==3) throw std::runtime_error("INJECTED_CALLBACK_EXCEPTION");
            if (mode==4) throw 17;
            if (mode==5) state.lastFactorTransactionFailureReason="INJECTED_REJECTION_REASON";
            return mode==1;
        };
    KFState entry;
    static_cast<KFState_&>(entry)=static_cast<const KFState_&>(state);
    MatrixXd transitionOutput=MatrixXd::Constant(2,3,-99), noiseOutput=MatrixXd::Constant(2,2,-77);
    const MatrixXd transitionSentinel=transitionOutput, noiseSentinel=noiseOutput;
    const bool accepted=state.stateTransition(trace,test.target,&transitionOutput,requestQ?&noiseOutput:nullptr);
    require(accepted==expectedSuccess,name+": unexpected transition acceptance: "+state.lastFactorTransactionFailureReason);
    require(callbackCalls==(mode?1:0), "Factor tap called wrong number of times");
    if (!expectedSuccess)
    {
        checkMetadataAndMoments(state,entry);
        require(sameBits(transitionOutput,transitionSentinel) && sameBits(noiseOutput,noiseSentinel),
                "Rejected prediction leaked external F/Q output");
        if (mode==3) require(state.lastFactorTransactionFailureReason=="STATE_TRANSITION_CALLBACK_EXCEPTION:INJECTED_CALLBACK_EXCEPTION", "Exception reason lost");
        if (mode==4) require(state.lastFactorTransactionFailureReason=="STATE_TRANSITION_CALLBACK_EXCEPTION:UNKNOWN", "Unknown exception reason lost");
        if (mode==5) require(state.lastFactorTransactionFailureReason=="INJECTED_REJECTION_REASON", "Callback rejection reason lost");
    }
    else
    {
        require(sameBits(transitionOutput,test.F), "Published F differs from closed-form reference");
        if (requestQ) require(sameBits(noiseOutput,test.Q), "Published dense Q differs from independent equations");
        else require(sameBits(noiseOutput,noiseSentinel), "Unrequested process output changed");
        require(sameBits(state.P,test.oldEigenCovariance), "Full P is not bitwise equal to original Eigen F*P*F^T+Q");
        require(sameBits(state.x,test.oldEigenMean) && sameBits(state.dx,test.adjustment),
                "x/dx are not bitwise equal to original sparse Eigen prediction");
        require((state.P-test.scalarCovariance).norm()<1e-10 && (state.x-test.mean).norm()<1e-11,
                "Independent scalar prediction/covariance oracle disagrees");
        require((state.dx-test.adjustment).norm()<1e-11, "dx or dimension-change behavior changed");
        if (!simulate && !legacyChunked) require(std::abs(state.P(test.destination.at(test.position),test.destination.at(test.clock)))>0,
                               "Cross-station covariance was deleted");
        // Independently apply only the documented successful metadata changes.
        entry.time=test.target;
        entry.x=test.oldEigenMean;
        entry.P=test.oldEigenCovariance;
        entry.dx=test.adjustment;
        entry.kfIndexMap=test.destination;
        entry.factorCommitSequence=43;
        if (test.elapsed) entry.exponentialNoiseMap.at(test.clock).value*=std::exp(-test.elapsed/15.);
        entry.initNoiseMap.clear();
        for (auto& [key, rows] : entry.stateTransitionMap)
            if (key!=entry.oneKey) rows.erase(entry.oneKey);
        entry.stateTransitionMap[entry.oneKey][entry.oneKey][0]=1;
        if (retire)
        {
            entry.pseudoStateMap.erase(test.old);
            entry.pseudoParentMap.erase(test.added);
        }
        checkMetadataAndMoments(state,entry);
        require(state.lastFactorTransactionFailureReason.empty(), "Successful prediction retained old failure reason");
        outputBits(name,state);
    }
    ++cases;
}

void emptyTransitionRollback()
{
    KFState state;
    state.time.bigTime=1234;
    state.stateTransitionMap.clear();
    state.metaDataMap["sentinel"]="empty-input";
    state.factorCommitSequence=7;
    KFState entry;
    static_cast<KFState_&>(entry)=static_cast<const KFState_&>(state);
    GTime target;
    target.bigTime=1240;
    MatrixXd F=MatrixXd::Constant(1,2,-4), Q=MatrixXd::Constant(3,3,-5);
    const auto beforeF=F, beforeQ=Q;
    std::ostringstream trace;
    require(!state.stateTransition(trace,target,&F,&Q), "Empty transition was accepted");
    checkMetadataAndMoments(state,entry);
    require(state.lastFactorTransactionFailureReason=="STATE_TRANSITION_HAS_NO_STATES" &&
            sameBits(F,beforeF) && sameBits(Q,beforeQ), "Empty transition rollback/publication failed");
    ++cases;
}

struct MutableTransitionTap
{
    int calls=0;
    int action=0;
    bool operator()(const KFState&,GTime,const map<KFKey,int>&,const map<KFKey,int>&,
        const SparseMatrix<double>&,const MatrixXd&,const string&,const VectorXd&,
        const MatrixXd&,std::uint64_t,std::uint64_t)
    {
        ++calls; // Internal value state, not an external/global audit counter.
        if (action==2) throw std::runtime_error("MUTABLE_COUNTER_STD_EXCEPTION");
        if (action==3) throw 19;
        return action==0;
    }
};

void mutableFunctorRollback()
{
    for (int action=0;action<4;++action)
        for (bool requestQ : {false,true})
        {
            KFState state;
            state.time.bigTime=1000;
            state.dx=VectorXd::Constant(1,-9);
            state.stateTransitionFactorCallback=MutableTransitionTap{0,action};
            KFState before;
            static_cast<KFState_&>(before)=static_cast<const KFState_&>(state);
            GTime target;
            target.bigTime=1001;
            MatrixXd Q=MatrixXd::Constant(2,2,-7);
            const auto outputBefore=Q;
            std::ostringstream trace;
            require(state.stateTransition(trace,target,nullptr,requestQ?&Q:nullptr)==(action==0),
                    "Mutable functor changed transition acceptance");
            const auto* targetFunctor=state.stateTransitionFactorCallback.target<MutableTransitionTap>();
            require(targetFunctor && targetFunctor->action==action && targetFunctor->calls==(action==0?1:0),
                    "Mutable callback value counter not restored on reject/throw, or success counter was rolled back");
            if (action!=0)
            {
                checkMetadataAndMoments(state,before);
                require(sameBits(Q,outputBefore),"Mutable rejected tap published Q");
            }
            ++cases;
        }
}

template<class Plain>
void writeMatrix(const std::filesystem::path& path, const Plain& value)
{
    require(!std::filesystem::exists(path), "Do not overwrite prediction evidence");
    std::ofstream stream(path,std::ios::binary);
    require(bool(stream), "Cannot create prediction matrix output");
    // Plain contiguous column-major storage; no dense copies or large strings.
    const char* bytes=reinterpret_cast<const char*>(value.data());
    std::size_t remaining=std::size_t(value.size())*sizeof(double);
    while (remaining)
    {
        const auto count=std::min<std::size_t>(remaining,1024*1024);
        stream.write(bytes,count);
        require(bool(stream), "Prediction matrix write failed");
        bytes+=count;
        remaining-=count;
    }
}

void realCheckpointPrediction(const string& inputPath,const string& outputPath,bool noTap)
{
    namespace fs=std::filesystem;
    const string inputSha="1310d0f6543e552ac318208b79491139d93886de87d8331ee14cb14dc3bf283b";
    require(zhangCheckpointFileSha256(inputPath)==inputSha,"Only the frozen real e121 checkpoint is allowed");
    require(!fs::exists(outputPath),"Real prediction requires a fresh exclusive output directory");
    ZhangCheckpointBundle bundle;
    const auto io=readZhangCheckpointBundle(inputPath,{},bundle);
    require(io.valid,"Real checkpoint read failed: "+io.failureReason);
    require(bundle.manifest.epoch=="2024-07-17 01:00:00" &&
            bundle.manifest.binarySha256=="86f69f9aca5c6e5ba238ba193e1127bb57b8ba5b092f99bb8a13ff102316b288",
            "Real checkpoint epoch/source producer differs");
    const auto inputStates=bundle.kfCore.x.size();
    require(inputStates>0 && inputStates<=20000,"Real core dimension outside controlled bounds");
    const auto sourceTransitions=bundle.kfCore.stateTransitionMap.size();
    const auto sourceProcessModels=bundle.kfCore.procNoiseMap.size();
    require(sourceTransitions>0,"Real checkpoint transition map is absent; never synthesize identity");
    KFState state;
    map<string,Receiver> receiverBindings;
    const auto resolve=[&](const string& id)->Receiver*
    {
        auto& receiver=receiverBindings[id];
        receiver.id=id;
        return &receiver;
    };
    string failure;
    require(restoreZhangCheckpointKfCoreWithReceiverResolver(bundle.kfCore,state,resolve,&failure),
            "Real checkpoint restore failed: "+failure);
    require(zhangCheckpointKfCoreMatchesState(bundle.kfCore,state),"Restored complete core differs from checkpoint");
    const string beforeCoreSha=zhangCheckpointKfCoreSha256Streaming(bundle.kfCore);
    require(state.rts_basename.empty(),"Controlled core prediction will not write inherited RTS output paths");
    GTime target=state.time;
    target.bigTime+=30;
    require(target.to_string(0)=="2024-07-17 01:00:30","Checkpoint+30s is not the fixed target epoch");
    // All transition/process/initialization maps above come from the actual
    // checkpoint. Full receiver event preparation is deliberately NOT replayed.
    // The test controls these existing non-core options equally on old/new.
    state.simulate_filter_only=false;
    acsConfig.pppOpts.receiver_chunking=false;
    bundle=ZhangCheckpointBundle{}; // Release the decoded input P/opaque payloads before prediction.
    require(fs::create_directory(outputPath),"Could not create exclusive real prediction directory");
    const fs::path output(outputPath);
    std::ofstream trace(output/"prediction.trace");
    require(bool(trace),"Cannot create controlled prediction trace");
    int taps=0;
    Eigen::Index transitionRows=0,transitionColumns=0,transitionNonzeros=0;
    if (!noTap) state.stateTransitionFactorCallback=[&](const KFState& source,GTime time,
        const map<KFKey,int>& oldKeys,const map<KFKey,int>& newKeys,
        const SparseMatrix<double>& F,const MatrixXd& Q,const string& label,
        const VectorXd& mean,const MatrixXd& covariance,std::uint64_t before,std::uint64_t after)
    {
        ++taps;
        require(&source==&state && time==target && label=="KF_STATE_TRANSITION" &&
                before==state.factorCommitSequence && after==before+1 &&
                oldKeys==state.kfIndexMap && F.rows()==Eigen::Index(newKeys.size()) &&
                F.cols()==state.x.size() && Q.rows()==F.rows() && Q.cols()==F.rows(),
                "Real prediction factor tap identity/dimension/sequence changed");
        require(mean.allFinite() && covariance.allFinite() && Q.allFinite(),
                "Real prediction factor tap contains nonfinite moments");
        transitionRows=F.rows(); transitionColumns=F.cols(); transitionNonzeros=F.nonZeros();
        writeMatrix(output/"Q.f64",Q);
        // Export the COMPLETE dense representation of sparse F one column at
        // a time: no new n*n dense matrix is retained merely for the test.
        std::ofstream coefficients(output/"F.f64",std::ios::binary);
        require(bool(coefficients),"Cannot create full transition output");
        VectorXd column=VectorXd::Zero(F.rows());
        for (Eigen::Index j=0;j<F.cols();++j)
        {
            column.setZero();
            for (SparseMatrix<double>::InnerIterator element(F,j);element;++element)
                column(element.row())=element.value();
            coefficients.write(reinterpret_cast<const char*>(column.data()),column.size()*sizeof(double));
            require(bool(coefficients),"Full transition output write failed");
        }
        return true;
    };
    require(state.stateTransition(trace,target),"Real production stateTransition rejected: "+state.lastFactorTransactionFailureReason);
    require(taps==(noTap?0:1) && state.time==target && state.x.allFinite() && state.P.allFinite(),
            "Real prediction did not commit one finite target");
    writeMatrix(output/"x.f64",state.x);
    writeMatrix(output/"P.f64",state.P);
    writeMatrix(output/"dx.f64",state.dx);
    const string afterCoreSha=zhangCheckpointKfStateCoreSha256Streaming(state);
    {
        const auto core=captureZhangCheckpointKfCore(state);
        std::ofstream stream(output/"core_after.bin",std::ios::binary);
        require(bool(stream),"Cannot create complete core output");
        boost::archive::binary_oarchive archive(stream,boost::archive::no_header);
        archive << core;
        require(bool(stream),"Complete core output write failed");
    }
    std::ofstream receipt(output/"receipt.json");
    require(bool(receipt),"Cannot create prediction receipt");
    receipt << "{\n  \"schema\": \"REAL_E121_STATE_TRANSITION_CONTROLLED_PREDICTION_V1\",\n"
        << "  \"passed\": true,\n  \"native_exit_code\": 0,\n  \"resource_stop\": null,\n"
        << "  \"source_checkpoint_bundle_sha256\": \""<<inputSha<<"\",\n"
        << "  \"before_core_sha256\": \""<<beforeCoreSha<<"\",\n"
        << "  \"after_core_sha256\": \""<<afterCoreSha<<"\",\n"
        << "  \"input_states\": "<<inputStates<<",\n  \"output_states\": "<<state.x.size()<<",\n"
        << "  \"checkpoint_transition_rows\": "<<sourceTransitions<<",\n"
        << "  \"checkpoint_process_models\": "<<sourceProcessModels<<",\n"
        << "  \"transition_rows\": "<<(noTap?state.x.size():transitionRows)<<",\n  \"transition_columns\": "<<(noTap?inputStates:transitionColumns)<<",\n"
        << "  \"transition_nonzeros\": "<<(noTap?string("null"):std::to_string(transitionNonzeros))<<",\n  \"factor_taps\": "<<taps<<",\n"
        << "  \"all_transition_and_process_maps_restored_not_synthesized\": true,\n"
        << "  \"same_controlled_options\": {\"simulate_filter_only\": false, \"receiver_chunking\": false},\n"
        << "  \"dense_Q_required_by_read_only_tap\": "<<(noTap?"false":"true")<<",\n"
        << "  \"F_and_Q_outputs_requested\": "<<(noTap?"false":"true")<<",\n"
        << "  \"memory_scope\": \""<<(noTap?"CONTROLLED_NO_TAP_PREDICTION_MEMORY":"CONTROLLED_DENSE_Q_ACCEPT_TAP_AUDIT_NOT_NO_CALLBACK_BENCHMARK")<<"\",\n"
        << "  \"full_output_files\": "<<(noTap?"[\"x.f64\",\"P.f64\",\"dx.f64\",\"core_after.bin\"]":"[\"x.f64\",\"P.f64\",\"dx.f64\",\"F.f64\",\"Q.f64\",\"core_after.bin\"]")<<",\n"
        << "  \"matrix_encoding\": \"native IEEE754 float64 column-major; shapes above\",\n"
        << "  \"all_core_fields_byte_exact\": null,\n"
        << "  \"old_new_equivalence_requires_external_complete_file_comparison\": true,\n"
        << "  \"full_epoch_or_physical_model_validation\": false,\n"
        << "  \"scope\": \"Actual checkpoint core prediction only; no observation/receiver-event preparation or measurement/AR update\"\n}\n";
    require(bool(receipt),"Prediction receipt write failed");
    std::cout << "REAL_E121_STATE_TRANSITION_PREDICTION_PASS states="<<inputStates
              <<" output_states="<<state.x.size()<<" taps="<<taps<<" full_epoch=0\n";
}
}

int main(int argc,char** argv)
{
    try
    {
        if ((argc==5 || (argc==6 && string(argv[5])=="--no-tap")) &&
            string(argv[1])=="--checkpoint-core" && string(argv[3])=="--output")
        {
            realCheckpointPrediction(argv[2],argv[4],argc==6);
            return 0;
        }
        require(argc==1,"Usage: driver [--checkpoint-core frozen_checkpoint.bundle --output fresh_directory [--no-tap]]");
        // A true global state forces the original dense fallback even if a
        // caller accidentally requested receiver chunking. Full cross-P remains.
        for (bool chunking : {false,true})
        {
            acsConfig.pppOpts.receiver_chunking=chunking;
            for (int scenario=0; scenario<8; ++scenario)
                for (int mode=0; mode<6; ++mode)
                    for (bool requestQ : {false,true})
                    {
                        const double dt=scenario==0?0.:scenario==1?-7.25:scenario==2?2.5:30.;
                        const bool grow=scenario==3||scenario==5;
                        const bool retire=scenario==4||scenario==5;
                        runCase(std::to_string(chunking)+"_"+std::to_string(scenario)+"_"+
                            std::to_string(mode)+"_"+std::to_string(requestQ),dt,grow,retire,
                            scenario==6,scenario==7,scenario==2,true,mode,requestQ);
                    }
        }
        acsConfig.pppOpts.receiver_chunking=true;
        for (int mode=0; mode<6; ++mode)
            for (bool requestQ : {false,true})
                runCase("EXISTING_CHUNK_BRANCH_"+std::to_string(mode)+"_"+std::to_string(requestQ),
                        2.5,true,true,true,false,false,false,mode,requestQ);
        emptyTransitionRollback();
        mutableFunctorRollback();
        std::cout<<"FLOAT_STATE_TRANSITION_MEMORY_CONTRACT_PASS cases="<<cases<<" checks="<<checks
                 <<" oracle=OLD_LITERAL_EIGEN_BITWISE_AND_INDEPENDENT_SCALAR full_cross_station_P=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr<<"FLOAT_STATE_TRANSITION_MEMORY_CONTRACT_FAIL "<<error.what()<<"\n";
        return 1;
    }
}


