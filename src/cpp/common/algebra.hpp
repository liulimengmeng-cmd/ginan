#pragma once

#include <boost/algorithm/string.hpp>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <functional>
#include <limits>
#include <map>
#include <math.h>
#include <mutex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "common/acsConfig.hpp"
#include "common/eigenIncluder.hpp"
#include "common/enums.h"
#include "common/gTime.hpp"
#include "common/satSys.hpp"
#include "common/trace.hpp"

using boost::algorithm::to_lower;
using boost::algorithm::to_upper;
using std::hash;
using std::lock_guard;
using std::map;
using std::pair;
using std::recursive_mutex;
using std::string;
using std::tuple;
using std::vector;

// forward declaration
struct KFMeasEntryList;
struct Receiver;
struct KFState;
struct KFMeasurementTransaction;

/** Keys used to interface with Kalman filter objects.
 * These have parameters to separate states of different 'type', for different 'Sat's, with
 * different receiver id 'str's and may have a different 'num' (eg xyz->0,1,2)
 *
 * Keys should be used rather than indices for accessing kalman filter state parameters.
 */
struct KFKey
{
    KF        type = KF::NONE;  ///< Key type (From enum)
    SatSys    Sat  = {};        ///< Satellite
    string    str;              ///< String (receiver ID)
    int       num = 0;          ///< Subkey number (eg xyz => 0,1,2)
    string    comment;          ///< Optional comment
    Receiver* rec_ptr = 0;      ///< Pointer to station object for dereferencing

    mutable GTime estimatedTime;

    bool operator!=(const KFKey& b) const;
    bool operator==(const KFKey& b) const;
    bool operator<(const KFKey& b) const;

    string code() const
    {
        string code;

        // Measurements or per-measurement states
        if (type == KF::CODE_MEAS || type == KF::PHAS_MEAS || type == KF::AMBIGUITY ||
            type == KF::Z_AMB || type == KF::CODE_BIAS || type == KF::PHASE_BIAS)
        {
            // IFLC combined
            if (num > 100)
            {
                int num1 = num / 100;
                int num2 = num % 100;

                string code1 = enum_to_string(static_cast<E_ObsCode>(num1));
                string code2 = enum_to_string(static_cast<E_ObsCode>(num2));
                code         = code1 + "-" + code2;
            }

            // Uncombined
            else
            {
                code = enum_to_string(static_cast<E_ObsCode>(num));
            }

            return code;
        }

        // STEC
        if (type == KF::IONO_STEC && acsConfig.pppOpts.ionoOpts.common_ionosphere == false)
        {
            code = (num == CODE) ? "CODE" : "PHASE";

            return code;
        }

        int component = 0;

        // Empirical force coefficients
        if (type >= KF::EMP_D_0 && type <= KF::EMP_Q_4 &&
            (static_cast<int>(type) - static_cast<int>(KF::EMP_D_0)) % 5 != 0)
        {
            // component = E_TrigType::COS + num;

            code = enum_to_string(enum_add(E_TrigType::COS, num));

            return code;
        }

        // Cartesian coordinates
        if ((type >= KF::REC_POS && type <= KF::ACC) ||
            (type > KF::BEGIN_INERTIAL_STATES && type < KF::END_INERTIAL_STATES) ||
            type == KF::ORBIT || type == KF::ORBIT_MEAS || type == KF::XFORM_XLATE ||
            type == KF::XFORM_RTATE || type == KF::XFORM_XLATE_RATE || type == KF::XFORM_RTATE_RATE)
        {
            component = enum_to_int(enum_add(E_StateComponent::X, num));
        }

        // Quaternions
        else if (type == KF::ORIENTATION)
        {
            component = enum_to_int(enum_add(E_StateComponent::W, num));
        }

        // Local tangental coordinates
        else if (type == KF::TROP_GRAD || type == KF::ANT_DELTA)
        {
            component = enum_to_int(enum_add(E_StateComponent::E, num));
        }

        // EOP parameters
        else if (type == KF::EOP || type == KF::EOP_RATE)
        {
            component = enum_to_int(enum_add(E_StateComponent::XP, num));
        }

        code = string(magic_enum::enum_name(static_cast<E_StateComponent>(component)));
        std::replace(code.begin(), code.end(), '_', '-');

        return code;
    }

    /** Create a string with the same spacing as ordinary outputs
     */
    static string emptyString()
    {
        KFKey  key;
        string keyStr = key;
        for (auto& c : keyStr)
        {
            if (c != '\t')
                c = ' ';
        }

        return keyStr;
    }

    operator string() const
    {
        char buff[100];

        snprintf(
            buff,
            sizeof(buff),
            "%10s\t%4s\t%4s\t%7s",
            enum_to_string(type).c_str(),
            Sat.id().c_str(),
            str.c_str(),
            this->code().c_str()
        );

        return string(buff);
    }

    string commaString() const
    {
        char buff[100];
        snprintf(
            buff,
            sizeof(buff),
            "%s,%s,%s,%s",
            enum_to_string(type).c_str(),
            Sat.id().c_str(),
            str.c_str(),
            this->code().c_str()
        );
        string keyStr = buff;

        return keyStr;
    }

    friend ostream& operator<<(ostream& os, const KFKey& kfKey)
    {
        string keyStr = kfKey;
        os << keyStr;

        return os;
    }

    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int& version)
    {
        ar & Sat;
        ar & str;
        ar & num;
        ar & type;
        ar & comment;
    }
};

struct FilterChunk
{
    string id;
    Trace* trace_ptr = nullptr;
    int    begX      = 0;
    int    numX      = 0;
    int    begH      = 0;
    int    numH      = -1;

    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int& version)
    {
        ar & id;
        ar & begX;
        ar & numX;
    }
};

struct ComponentsDetails
{
    double value = 0;
    string eq;  // not valid after combinations
    double var = 0;

    ComponentsDetails& operator+=(const ComponentsDetails& rhs)
    {
        value += rhs.value;
        var += rhs.var;
        return *this;
    }

    ComponentsDetails operator*(double rhs)
    {
        ComponentsDetails newDetails = *this;
        newDetails.value *= rhs;
        newDetails.var *= rhs;
        newDetails.var *= rhs;
        return newDetails;
    }
};

/** Object to hold measurements, design matrices, and residuals for multiple observations
 */
struct KFMeas
{
    GTime    time = GTime::noTime();  ///< Epoch these measurements were recorded
    VectorXd Y;                       ///< Value of the observations (for linear systems)
    VectorXd V;                  ///< Prefit Residual of the observations (for non-linear systems)
    VectorXd VV;                 ///< Postfit Residual of the observations (for non-linear systems)
    VectorXd W;                  ///< Weight (inverse of noise) used in least squares
    MatrixXd R;                  ///< Measurement noise for these observations
    MatrixXd H;                  ///< Design matrix between measurements and state
    MatrixXd H_star;             ///< Design matrix between measurements and noise states
    VectorXd uncorrelatedNoise;  ///< Uncorellated noise for measurements
    VectorXd prefitRatios;       ///< Prefit sigma check or omega test ratios of measurements
    VectorXd postfitRatios;      ///< Postfit sigma check or omega test ratios of measurements

    map<KFKey, int> noiseIndexMap;  ///< Map from key to indexes of parameters in the noise vector
    vector<KFKey>
        obsKeys;  ///< Vector of optional labels for reporting when measurements are removed etc.
    vector<map<string, void*>>                  metaDataMaps;
    vector<map<E_Component, ComponentsDetails>> componentsMaps;

    KFMeas() {

    };

    KFMeas(
        KFMeas&                      kfMeas,       ///< Measurement to form linear combination from
        vector<Triplet<double>>&&    triplets,     ///< Linear combination triplets
        vector<KFKey>&&              obsKeys,      ///< New obs key vector
        vector<map<string, void*>>&& metaDataMaps  ///< Optional new metadata vector
    )
        : obsKeys{obsKeys}, metaDataMaps{metaDataMaps}
    {
        auto F = SparseMatrix<double>(obsKeys.size(), kfMeas.obsKeys.size());

        F.setFromTriplets(triplets.begin(), triplets.end());

        time = kfMeas.time;
        //		Y					= F *	kfMeas.Y;
        V  = F * kfMeas.V;
        VV = V;
        //		W					= F *	kfMeas.W;
        R                 = F * kfMeas.R * F.transpose();
        H                 = F * kfMeas.H;
        H_star            = F * kfMeas.H_star;
        uncorrelatedNoise = kfMeas.uncorrelatedNoise;

        componentsMaps.resize(obsKeys.size());
        for (auto& triplet : triplets)
        {
            auto   newIndex = triplet.row();
            auto   oldIndex = triplet.col();
            double scalar   = triplet.value();

            for (auto& [component, details] : kfMeas.componentsMaps[oldIndex])
            {
                componentsMaps[newIndex][component] += details * scalar;
            }
        }
    }

    KFMeas(
        KFState&         kfState,
        KFMeasEntryList& kfEntryList,
        GTime            measTime        = GTime::noTime(),
        MatrixXd*        noiseMatrix_ptr = nullptr
    );

    int getNoiseIndex(const KFKey& key) const;

    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int& version)
    {
        int rows = H.rows();
        int cols = H.cols();
        ar & rows;
        ar & cols;

        if (ARCHIVE::is_saving::value)
        {
            // just wrote this, we are writing
            map<pair<int, int>, double> H2;

            ar & obsKeys;
            ar & time;
            ar & VV;

            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                {
                    double value = H(i, j);
                    if (value)
                    {
                        H2[{i, j}] = value;
                    }
                }

            ar & H2;
        }
        else
        {
            // we're reading
            map<pair<int, int>, double> H2;

            ar & obsKeys;
            ar & time;
            ar & VV;
            ar & H2;

            H = MatrixXd::Zero(rows, cols);
            R = MatrixXd::Zero(rows, rows);
            V = VectorXd::Zero(rows);

            for (auto& [index, value] : H2)
            {
                H(index.first, index.second) = value;
            }
        }
    }
};

/** Object to hold the values used to initialise new states when adding to the kalman filter object
 */
struct InitialState
{
    bool   estimate         = false;
    bool   use_remote_sigma = false;
    double x                = 0;   ///< State value
    double P                = -1;  ///< State Covariance
    double sigmaMax         = 0;   ///< Sigma limit
    double outageLimit      = 0;   ///< Maxiumum time without state estimation
    double Q                = 0;   ///< Process Noise, -ve indicates infinite (throw away state)
    double tau              = -1;  ///< Correlation Time, default to -1 (inf) (Random Walk)
    double mu               = 0;   ///< Desired Mean Value
    string comment;
};

struct KFStatistics
{
    double averageRatio = 0;
    double sumOfSquares = 0;
};

struct KFState;
struct KalmanModel;
struct KFMeasEntry;

InitialState initialStateFromConfig(const KalmanModel& kalmanModel, int index = 0);

typedef std::ostream Trace;

/** Validate every dense block extent before a Kalman update reaches BLAS.
 * In particular, DGEMV requires a positive leading dimension even when an
 * upstream chunk contains measurements but no state columns. */
inline bool validKalmanFilterBlasBlock(
    int begX,
    int numX,
    int begH,
    int numH,
    int hRows,
    int hCols,
    int pRows,
    int pCols,
    int rRows,
    int rCols,
    int vRows)
{
    return begX >= 0 && begH >= 0 && numX > 0 && numH > 0 &&
        begX + numX <= hCols && begH + numH <= hRows &&
        begX + numX <= pRows && begX + numX <= pCols &&
        begH + numH <= rRows && begH + numH <= rCols &&
        begH + numH <= vRows;
}

/** A filter chunk can update the state only when it contains both state
 * columns and measurement rows.  Measurement-only bookkeeping chunks are
 * valid upstream objects, but they must never be sent to the dense BLAS
 * state-update path. */
inline bool activeKalmanFilterStateChunk(int numX, int numH)
{
    return numX > 0 && numH > 0;
}

struct KFMeasList : vector<KFMeas>
{
};

struct KFMeasEntryList : vector<KFMeasEntry>
{
};

struct RejectCallbackDetails
{
    Trace&   trace;   ///< Trace to output to
    KFState& kfState;
    KFMeas&  kfMeas;  ///< Measurements, noise, and design matrix

    RejectCallbackDetails(Trace& trace, KFState& kfState, KFMeas& kfMeas)
        : trace{trace}, kfState{kfState}, kfMeas{kfMeas}
    {
    }

    KFKey         kfKey;            ///< Key to the state that is flagged as an error
    int           stateIndex = -1;  ///< Index of the state that is flagged as an error
    int           measIndex  = -1;  ///< Index of the measurement that is flagged as an outlier
    E_FilterStage stage      = E_FilterStage::LSQ;  ///< prefit, postfit, least squares
};

typedef bool (*StateRejectCallback)(RejectCallbackDetails rejectCallbackDetails);
typedef bool (*MeasRejectCallback)(RejectCallbackDetails rejectCallbackDetails);

/** Build the stochastic reset used when least squares gives a fresh state its
 * first finite covariance.  Rows selected for initialisation do not inherit
 * the random variable previously occupying the same semantic KF key. */
inline bool makeKFStateInitialisationTransition(
    const MatrixXd&       sourceCovariance,
    const MatrixXd&       destinationCovariance,
    const vector<int>&    initialisedStateIndicies,
    SparseMatrix<double>& transition,
    MatrixXd&             processCovariance
)
{
    if (sourceCovariance.rows() == 0
        || sourceCovariance.rows() != sourceCovariance.cols()
        || destinationCovariance.rows() != sourceCovariance.rows()
        || destinationCovariance.cols() != sourceCovariance.cols()
        || initialisedStateIndicies.empty()
        || !sourceCovariance.allFinite()
        || !destinationCovariance.allFinite())
    {
        return false;
    }

    transition.resize(sourceCovariance.rows(), sourceCovariance.cols());
    transition.setIdentity();
    for (const int index : initialisedStateIndicies)
    {
        if (index < 0 || index >= transition.rows())
        {
            return false;
        }
        transition.coeffRef(index, index) = 0;
    }
    transition.prune([](int, int, double value) { return value != 0; });

    processCovariance = destinationCovariance
        - transition * sourceCovariance * transition.transpose();
    processCovariance = 0.5
        * (processCovariance + processCovariance.transpose());
    return processCovariance.allFinite();
}

struct Exponential
{
    double value = 0;
    double tau   = 0;

    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int& version)
    {
        ar & value;
        ar & tau;
    }
};

/** Resolve the negative "to the end" sentinels before any chunk is screened
 * or passed to BLAS. */
inline void resolveKalmanFilterChunkExtents(FilterChunk& chunk, int stateRows, int measurementRows)
{
    if (chunk.numX < 0)
    {
        chunk.numX = stateRows - chunk.begX;
    }
    if (chunk.numH < 0)
    {
        chunk.numH = measurementRows - chunk.begH;
    }
}

/** Result of an authoritative Kalman measurement transaction. */
enum class KFFilterResult
{
    COMMITTED,
    NO_MEASUREMENTS,
    NUMERICAL_FAILURE,
    FACTOR_TRANSACTION_REJECTED
};

enum class KFStagedFactorKind
{
    STATE_TRANSITION,
    MEASUREMENT
};

/** One factor in an ordered, provisional measurement transaction.
 *
 * Adaptive-QC state relaxations reproduce the numerical prior, but are not
 * independent measurement evidence.
 */
struct KFStagedFactor
{
    KFStagedFactorKind kind = KFStagedFactorKind::STATE_TRANSITION;
    GTime              time = GTime::noTime();
    string             label;
    bool               adaptiveQcTransition = false;

    map<KFKey, int> sourceIndexMap;
    map<KFKey, int> destinationIndexMap;
    VectorXd        sourceMean;
    MatrixXd        sourceCovariance;
    SparseMatrix<double> transition;
    // Robust-QC relaxations normally touch only one state.  Keeping their Q
    // matrices dense made six staged 4k-state relaxations consume more than a
    // gigabyte before the packet could be validated.
    SparseMatrix<double> processCovariance;
    VectorXd             destinationMean;
    MatrixXd             destinationCovariance;

    KFMeas measurement;

    // Production /PPP packets are validated synchronously inside
    // filterKalman().  These non-owning views avoid a second dense copy of
    // H/R and of the posterior covariance.  Tests and offline packets may
    // continue to use the owned fields above.
    const KFMeas*  measurement_ptr = nullptr;
    const VectorXd* destinationMean_ptr = nullptr;
    const MatrixXd* destinationCovariance_ptr = nullptr;

    std::uint64_t beforeCommitSequence = 0;
    std::uint64_t afterCommitSequence  = 0;
};

/** Complete provisional packet for one robust measurement update. */
struct KFMeasurementTransaction
{
    GTime  time = GTime::noTime();
    string suffix;
    string transactionId;
    Trace* trace_ptr = nullptr;

    // One authoritative entry snapshot is sufficient to replay every staged
    // sparse transition.  Per-transition dense covariance snapshots would be
    // prohibitive for the 4k-5k state global PPP filter.
    map<KFKey, int> entryIndexMap;
    VectorXd        entryMean;
    MatrixXd        entryCovariance;

    vector<KFStagedFactor> stagedFactors;

    std::uint64_t entryCommitSequence = 0;
    std::uint64_t provisionalCommitSequence = 0;
    std::size_t   modelGeneration = 0;
    std::size_t   solveGeneration = 0;

    bool qcConverged = true;
    bool maxIterationsReached = false;
    bool terminalReconciliationPerformed = false;

    // Reject callbacks may request mutations of external observation/noise
    // objects.  Defer them until the estimator and all factor consumers have
    // accepted the complete packet, so rollback cannot leave shallow-pointer
    // side effects behind.
    vector<std::function<void()>> commitSideEffects;
    string primaryError;
};

/** Deterministic compact identity for a state mean/covariance pair.
 *
 * This is an FNV-1a audit identity, not a cryptographic digest. Dimensions
 * and IEEE-754 bytes are included; -0 is canonicalised to +0 and non-finite
 * inputs are rejected by returning "NONFINITE".
 */
inline string kfMomentHash(const VectorXd& mean, const MatrixXd& covariance)
{
    if (!mean.allFinite() || !covariance.allFinite())
    {
        return "NONFINITE";
    }

    std::uint64_t hash = 1469598103934665603ULL;
    auto append = [&](const void* data, std::size_t bytes)
    {
        const auto* values = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < bytes; i++)
        {
            hash ^= values[i];
            hash *= 1099511628211ULL;
        }
    };
    auto appendInteger = [&](std::int64_t value)
    {
        append(&value, sizeof(value));
    };
    auto appendDouble = [&](double value)
    {
        if (value == 0)
        {
            value = 0;
        }
        std::uint64_t bits = 0;
        std::memcpy(&bits, &value, sizeof(bits));
        append(&bits, sizeof(bits));
    };

    appendInteger(mean.size());
    appendInteger(covariance.rows());
    appendInteger(covariance.cols());
    for (int i = 0; i < mean.size(); i++)
    {
        appendDouble(mean(i));
    }
    for (int column = 0; column < covariance.cols(); column++)
    for (int row = 0; row < covariance.rows(); row++)
    {
        appendDouble(covariance(row, column));
    }

    std::ostringstream identity;
    identity << std::hex << std::setw(16) << std::setfill('0') << hash;
    return identity.str();
}

/** Production moment operation shared by transaction staging and tests. */
inline bool kfApplyStochasticTransitionMoments(
    const VectorXd& sourceMean,
    const MatrixXd& sourceCovariance,
    const MatrixXd& transition,
    const MatrixXd& processCovariance,
    VectorXd& destinationMean,
    MatrixXd& destinationCovariance)
{
    const int dimension = sourceMean.size();
    if (dimension <= 0
        || sourceCovariance.rows() != dimension
        || sourceCovariance.cols() != dimension
        || transition.rows() != dimension
        || transition.cols() != dimension
        || processCovariance.rows() != dimension
        || processCovariance.cols() != dimension
        || !sourceMean.allFinite()
        || !sourceCovariance.allFinite()
        || !transition.allFinite()
        || !processCovariance.allFinite())
    {
        return false;
    }

    destinationMean = transition * sourceMean;
    destinationCovariance = transition * sourceCovariance
        * transition.transpose() + processCovariance;
    destinationCovariance =
        (0.5 * (destinationCovariance + destinationCovariance.transpose())).eval();
    if (!destinationMean.allFinite() || !destinationCovariance.allFinite())
    {
        return false;
    }
    const double scale = std::max(
        1.0, destinationCovariance.diagonal().cwiseAbs().maxCoeff());
    return destinationCovariance.diagonal().minCoeff() >= -1e-12 * scale;
}

/** Kalman filter object.
 *
 * Contains most persistent parameters and values of state. Includes state vector, covariance, and
 * process noise.
 *
 * This object performs all operations on the kalman filter to ensure that edge cases are included
 * and state kept in a valid configuration.
 */
struct KFState_ : FilterOptions
{
    bool lsqRequired = false;  ///< Uninitialised parameters require least squares calculation

    GTime    time = {};
    VectorXd x;                  ///< State
    MatrixXd P;                  ///< State Covariance
    VectorXd dx;                 ///< Last filter update
    VectorXd prefitRatios;       ///< Prefit sigma check or omega test ratios of states
    VectorXd postfitRatios;      ///< Postfit sigma check or omega test ratios of states

    map<KFKey, int> kfIndexMap;  ///< Map from key to indexes of parameters in the state vector

    map<KFKey, map<KFKey, map<int, double>>> stateTransitionMap;
    map<KFKey, double>                       gaussMarkovTauMap;
    map<KFKey, double>                       gaussMarkovMuMap;
    map<KFKey, double>                       procNoiseMap;
    map<KFKey, double>                       initNoiseMap;
    map<KFKey, double>                       sigmaMaxMap;
    map<KFKey, double>                       outageLimitMap;
    map<KFKey, Exponential>                  exponentialNoiseMap;

    map<KFKey, map<KFKey, double>>
        pseudoStateMap;   ///< Map of pseudo states, and a further map of their coefficients
    map<KFKey, KFKey>
        pseudoParentMap;  ///< Map from ordinary states to their combined pseudo state parent.

    map<KFKey, int> errorCountMap;

    vector<StateRejectCallback> stateRejectCallbacks;
    vector<MeasRejectCallback>  measRejectCallbacks;

    /** Monotone identity for authoritative KF commits observed by read-only
     * factor consumers.  It is deliberately copied with KFState branches but
     * is not part of the estimator state vector. */
    std::uint64_t factorCommitSequence = 0;

    /** Non-owning pointer valid only while filterKalman() executes an atomic
     * packet. It is explicitly cleared by KFState copy operations. */
    KFMeasurementTransaction* activeMeasurementTransaction = nullptr;

    /** Transaction replay may use an immutable state copy while retaining the
     * stable runtime owner guard. */
    const KFState* factorCallbackOwnerOverride = nullptr;

    string lastFactorTransactionFailureReason;
    // These flags describe the last authoritative /PPP packet only.  A later
    // non-PPP filter (for example /AR) must never overwrite this gate.
    bool lastPppTransactionCommitted = false;
    bool lastPppQcConverged = false;

    /** Optional read-only factor taps.  They are invoked only after the
     * corresponding operation has passed its numerical/QC checks. */
    std::function<bool(
        const KFState&,
        const KFMeas&,
        const string&,
        const VectorXd&,
        const MatrixXd&,
        std::uint64_t,
        std::uint64_t
    )> acceptedMeasurementFactorCallback;
    std::function<bool(
        const KFState&,
        GTime,
        const map<KFKey, int>&,
        const map<KFKey, int>&,
        const SparseMatrix<double>&,
        const MatrixXd&,
        const string&,
        const VectorXd&,
        const MatrixXd&,
        std::uint64_t,
        std::uint64_t
    )> stateTransitionFactorCallback;
    std::function<bool(
        const KFState&,
        GTime,
        const map<KFKey, int>&,
        const map<KFKey, int>&,
        const SparseMatrix<double>&,
        const string&,
        const VectorXd&,
        const MatrixXd&,
        std::uint64_t,
        std::uint64_t
    )> exactStateTransformCallback;

    /** Validate and atomically commit a complete ordered factor packet. */
    std::function<bool(
        const KFState&,
        const KFMeasurementTransaction&,
        string&
    )> authoritativeFactorTransactionCallback;

    map<string, FilterChunk> filterChunkMap;

    map<string, string> metaDataMap;

    bool sigmaPass = false;

    bool   chiQCPass  = false;
    double chi2       = 0;
    int    dof        = 0;
    double chi2PerDof = INFINITY;
    double qc         = 0;

    string id = "KFState";

    string rts_basename = "";

    bool output_residuals        = false;
    bool outputMongoMeasurements = false;

    KFState* alternate_ptr = nullptr;

    map<string, int> statisticsMap;
    map<string, int> statisticsMapSum;
};

/** Wrapper to simplify copying with default copy but overriding slightly.
 */
struct KFState : KFState_
{
    mutable recursive_mutex kfStateMutex;

    static const KFKey oneKey;  ///< KFStates generally contain a ONE state as the first element,
                                ///< used for converting matrix additions to matrix multiplications.

    KFState(const KFState& kfState) : KFState_(kfState), kfStateMutex()
    {
        // dont use same rts file unless explicitly copied
        rts_basename.clear();
        activeMeasurementTransaction = nullptr;
        factorCallbackOwnerOverride  = nullptr;
    }

    KFState()
    {
        // initialise all filter state objects with a ONE element for later use.
        x  = VectorXd::Ones(1);
        P  = MatrixXd::Zero(1, 1);
        dx = VectorXd::Zero(1);

        kfIndexMap[oneKey] = 0;

        initFilterEpoch(nullStream);
    }

    KFState& operator=(const KFState& kfState)
    {
        KFState_* thisKfState_ = (KFState_*)this;
        KFState_* thatKfState_ = (KFState_*)&kfState;

        *thisKfState_ = *thatKfState_;

        // dont use same rts file unless explicitly copied
        rts_basename.clear();
        activeMeasurementTransaction = nullptr;
        factorCallbackOwnerOverride  = nullptr;

        return *this;
    }

    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, const unsigned int& version)
    {
        ar & kfIndexMap;
        ar & time;
        ar & x;
        ar & dx;
        ar & filterChunkMap;

        double num;
        int    rows = P.rows();
        ar & rows;

        if (ARCHIVE::is_saving::value)
        {
            for (int i = 0; i < P.rows(); i++)
                for (int j = 0; j <= i; j++)
                {
                    num = P(i, j);

                    ar & num;
                }
        }
        else
        {
            P = MatrixXd(rows, rows);

            for (int i = 0; i < P.rows(); i++)
                for (int j = 0; j <= i; j++)
                {
                    ar & num;

                    P(i, j) = num;
                    P(j, i) = num;
                }
        }
    }

    void initFilterEpoch(Trace& trace);

    int getKFIndex(const KFKey& key) const;

    E_Source getKFValue(
        const KFKey& key,
        double&      value,
        double*      variance       = nullptr,
        double*      adjustment_ptr = nullptr,
        bool         allowAlternate = true
    ) const;

    E_Source getPseudoValue(
        const KFKey& key,
        double&      value,
        double*      variance       = nullptr,
        double*      adjustment_ptr = nullptr
    ) const;

    bool getKFSigma(const KFKey& key, double& sigma);

    bool addKFState(const KFKey& kfKey, const InitialState& initialState = {});

    bool addPseudoState(const KFKey& kfKey, const map<KFKey, double>& coeffMap);

    void setExponentialNoise(const KFKey& kfKey, const Exponential exponential);

    void setAccelerator(
        const KFKey&        element,
        const KFKey&        dotElement,
        const KFKey&        dotDotElement,
        const double        value,
        const InitialState& initialState = {}
    );

    void setKFTrans(
        const KFKey&        dest,
        const KFKey&        source,
        const double        value,
        const InitialState& initialState = {}
    );

    void setKFTransRate(
        const KFKey&        integral,
        const KFKey&        rate,
        const double        value,
        const InitialState& initialRateState     = {},
        const InitialState& initialIntegralState = {}
    );

    void addNoiseElement(const KFKey& obsKey, const double variance);

    void removeState(const KFKey& kfKey, bool allowDeleteParent = true);

    bool stateTransition(
        Trace&    trace,
        GTime     newTime,
        MatrixXd* stm_ptr = nullptr,
        MatrixXd* processNoise_ptr = nullptr
    );

    bool manualStateTransition(
        Trace& trace,
        GTime newTime,
        MatrixXd& stm,
        MatrixXd& procNoise,
        const string& label = "KF_MANUAL_STATE_TRANSITION"
    );

    /** Apply an exact linear change of state coordinates.
     *
     * Each destination state is expressed as a linear combination of states in the current
     * coordinate system.  The full covariance, including all cross-covariances, is transformed
     * with P' = T P T^T.  This is intended for datum/S-basis changes where no stochastic
     * information is added or removed.
     */
    bool applyStateTransform(
        Trace&                                  trace,
        const map<KFKey, map<KFKey, double>>&   transformMap,
        const string&                           label = ""
    );

    void leastSquareSigmaChecks(
        RejectCallbackDetails& callbackDetails,
        MatrixXd&              Pp,
        KFStatistics&          statistics
    );

    void preFitSigmaChecks(
        RejectCallbackDetails& callbackDetails,
        KFStatistics&          statistics,
        int                    begX,
        int                    numX,
        int                    begH,
        int                    numH
    );

    void postFitSigmaChecks(
        RejectCallbackDetails& callbackDetails,
        VectorXd&              dx,
        MatrixXd&              Qinv,
        MatrixXd&              QinvH,
        KFStatistics&          statistics,
        int                    begX,
        int                    numX,
        int                    begH,
        int                    numH
    );

    double stateChiSquare(
        Trace&    trace,
        MatrixXd& Pp,
        VectorXd& dx,
        int       begX,
        int       numX,
        int       begH,
        int       numH
    );

    double measChiSquare(
        Trace&    trace,
        KFMeas&   kfMeas,
        VectorXd& dx,
        int       begX,
        int       numX,
        int       begH,
        int       numH
    );

    double innovChiSquare(Trace& trace, KFMeas& kfMeas, int begX, int numX, int begH, int numH);

    bool kFilter(
        Trace&    trace,
        KFMeas&   kfMeas,
        VectorXd& xp,
        MatrixXd& Pp,
        VectorXd& dx,
        MatrixXd& Qinv,
        MatrixXd& QinvH,
        int       begX           = 0,
        int       numX           = -1,
        int       begH           = 0,
        int       numH           = -1,
        bool      resetOnFailure = true
    );

    bool leastSquare(Trace& trace, KFMeas& kfMeas, VectorXd& xp, MatrixXd& Pp);

    void chiQC(Trace& trace, KFMeas& kfMeas);

    void
    outputStates(Trace& trace, string suffix = "", int iteration = -1, int begX = 0, int numX = -1);

    void outputConditionNumber(Trace& trace);

    void outputCorrelations(Trace& trace);

    void outputMeasurements(Trace& trace, KFMeas& meas);

    bool doStateRejectCallbacks(RejectCallbackDetails rejectDetails);

    bool doMeasRejectCallbacks(RejectCallbackDetails rejectDetails);

    KFFilterResult filterKalman(
        Trace&                    trace,
        KFMeas&                   kfMeas,
        const string&             suffix             = "",
        bool                      innovReady         = false,
        map<string, FilterChunk>* filterChunkMap_ptr = nullptr,
        KFState*                  finalMeasurementPrior_ptr = nullptr
    );

    bool leastSquareInitStates(
        Trace&        trace,
        KFMeas&       kfMeas,
        const string& suffix,
        bool          initCovars   = false,
        bool          innovReady   = false,
        bool          skipLsqCheck = false
    );

    VectorXd getSubState(
        map<KFKey, int>& kfKeyMap,
        MatrixXd*        covarMat_ptr  = nullptr,
        VectorXd*        adjustVec_ptr = nullptr
    ) const;

    void getSubState(map<KFKey, int>& kfKeyMap, KFState& kfState) const;

    KFState getSubState(vector<KF> types, KFMeas* meas_ptr = nullptr) const;

    vector<KFKey> decomposedStateKeys(const KFKey& composedKey) const;
};

/** Object to hold an individual measurement.
 * Includes the measurement itself, (or its innovation) and design matrix entries
 * Adding design matrix entries for states that do not yet exist will create and add new states to
 * the measurement's kalman filter object.
 */
struct KFMeasEntry
{
    KFState* kfState_ptr = nullptr;  ///< Pointer to filter object that measurements are referencing

    double valid = true;  ///< Optional parameter to invalidate a measurement (to avoid needing to
                          ///< delete it and reshuffle a vector)
    double value  = 0;    ///< Value of measurement (for linear systems)
    double noise  = 0;    ///< Noise of measurement
    double innov  = 0;    ///< Innovation of measurement (for non-linear systems)
    KFKey  obsKey = {};   ///< Optional labels to be used in output traces

    map<E_Component, ComponentsDetails> componentsMap;

    map<KFKey, double> noiseElementMap;
    map<KFKey, double> designEntryMap;
    map<KFKey, double> usedValueMap;
    map<KFKey, double> noiseEntryMap;
    map<string, void*> metaDataMap;

    KFMeasEntry(KFState* kfState_ptr, KFKey obsKey = {}) : kfState_ptr(kfState_ptr), obsKey(obsKey)
    {
    }

    KFMeasEntry() {}

    /** Adds a noise element for this measurement
     */
    void addNoiseEntry(
        const KFKey  kfKey,    ///< Key to determine the origin of the noise
        const double value,    ///< Noise entry matrix entry value
        const double variance  ///< Variance of noise element
    )
    {
        if (value == 0 || variance <= 0)
        {
            return;
        }

        noiseElementMap[kfKey] = variance;
        noiseEntryMap[kfKey] += value;
    }

    /** Adds a design matrix entry for this measurement
     */
    void addDsgnEntry(
        const KFKey&        kfKey,  ///< Key to determine which state parameter is affected
        const double        value,  ///< Design matrix entry value
        const InitialState& initialState = {}  ///< Initial conditions for new states
    )
    {
        if (value ==
            0)  // Eugene: Design entry value can be 0 in theory but still needs intial state?
        {
            return;
        }

        if (initialState.Q < 0)
        {
            addNoiseEntry(kfKey, value, initialState.P);
            return;
        }

        if (kfState_ptr)
        {
            auto& kfState = *kfState_ptr;

            kfState.addKFState(kfKey, initialState);

            auto it = kfState.outageLimitMap.find(kfKey);
            if (it != kfState.outageLimitMap.end())
            {
                auto& [editKey, dummy] = *it;

                editKey.estimatedTime = kfState.time;
            }
        }

        usedValueMap[kfKey] = initialState.x;
        designEntryMap[kfKey] += value;
    }

    /** Adds the measurement noise entry for this measurement
     */
    void setNoise(const double value)  ///< Measurement noise matrix entry value
    {
        if (value == 0)
        {
            std::cout << "Zero noise encountered" << "\n";
            // 			return;
        }
        if (std::isinf(value))
        {
            std::cout << "Inf noise encountered" << "\n";
            return;
        }
        else if (std::isnan(value))
        {
            std::cout << "Nan noise encountered" << "\n";
            return;
        }

        this->noise = value;
    }

    /** Adds the actual measurement value for this measurement
     */
    void setValue(const double value)  ///< Actual measurement entry value
    {
        this->value = value;
    }

    /** Adds the innovation value for this measurement
     */
    void setInnov(const double value)  ///< Innovation entry value
    {
        this->innov = value;
    }
};

KFState mergeFilters(const vector<KFState*>& kfStatePointerList, const vector<KF>& stateList);

MatrixXi correlationMatrix(MatrixXd& P);

void outputResiduals(
    Trace&  trace,
    KFMeas& kfMeas,
    string  suffix    = "",
    int     iteration = -1,
    int     begH      = 0,
    int     numH      = -1
);

bool isPositiveSemiDefinite(MatrixXd& mat);

int filter_(
    const double* x,
    const double* P,
    const double* H,
    const double* v,
    const double* R,
    int           n,
    int           m,
    double*       xp,
    double*       Pp
);

// matrix and vector functions
double* mat(int n, int m);
int*    imat(int n, int m);
double* zeros(int n, int m);
double* eye(int n);
double  dot(const double* a, const double* b, int n);
double  norm(const double* a, int n);
void    matcpy(double* A, const double* B, int n, int m);
void    matmul(
    const char*   tr,
    int           n,
    int           k,
    int           m,
    double        alpha,
    const double* A,
    const double* B,
    double        beta,
    double*       C
);
int matinv(double* A, int n);
int solve(const char* tr, const double* A, const double* Y, int n, int m, double* X);
