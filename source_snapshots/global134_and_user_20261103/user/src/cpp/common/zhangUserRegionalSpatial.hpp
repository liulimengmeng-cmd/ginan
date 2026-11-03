#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Experimental spatial model, not a server-supplied user-location correction.
// The same-station satellite difference removes the receiver GF code datum.
// sourceMean/sourceCovariance retain the original complete CLOCK/PHASE/IONO
// ordering. Coordinates below are local east/north kilometres FROM THE USER;
// no evaluation truth coordinate is required or accepted by this interface.
namespace zhangregional
{
inline const char* spatialModelLabel()
{
    return "LOCAL_SD_IONO_PLANE_PLUS_NEAREST4_IDW2_RESIDUAL_NOT_MPFG";
}

struct StationSupport
{
    std::string station;
    double eastKm = 0;
    double northKm = 0;
    int satelliteSourceIndex = -1;
    int referenceSourceIndex = -1;
};

struct PairSupport
{
    std::string satellite;
    std::string reference;
    std::vector<StationSupport> stations;
};

struct Options
{
    int minReferences = 4;
    double maxCondition = 1e6;
    double maxDistanceKm = 150;
    double maxWeightSquaredNorm = 4;
    double perSatelliteFloorMetres = 0.05;
};

struct PairDiagnostic
{
    std::string model = spatialModelLabel();
    std::string satellite;
    std::string reference;
    bool accepted = false;
    std::string reason;
    std::vector<std::string> stations;
    Eigen::VectorXd weights;
    Eigen::VectorXd inverseDistanceWeights;
    Eigen::VectorXd leaveOneOutResiduals;       // normalised by sqrt(m)*sqrt(1+||w_loo||^2)
    Eigen::VectorXd rawLeaveOneOutResiduals;
    Eigen::VectorXd leaveOneOutWeightSquaredNorms;
    Eigen::MatrixXd leaveOneOutWeights;         // omitted reference has coefficient zero
    Eigen::MatrixXd leaveOneOutInverseDistanceWeights;
    double geometryCondition = 0;
    double maximumLeaveOneOutCondition = 0;
    double maxDistanceKm = 0;
    double weightSquaredNorm = 0;
    bool insideConvexHull = false;
};

struct Result
{
    std::string model = spatialModelLabel();
    bool valid = false;
    std::string reason;
    Eigen::MatrixXd prediction;                 // W: target SD ionosphere metres = W * sourceMean
    Eigen::MatrixXd spatialCovariance;          // final experimental spatial mismatch covariance
    Eigen::MatrixXd stationWeights;             // K x globalStations.size(), zero where absent
    std::vector<std::string> globalStations;
    std::vector<PairDiagnostic> acceptedPairs;   // same order as prediction rows
    std::vector<PairDiagnostic> diagnostics;    // original input order, including rejected pairs
    Eigen::MatrixXd residualMatrix;             // normalised LOO residuals, K x global references
    Eigen::MatrixXd formalResidualCovariance;   // sum_r A_r Q A_r^T
    Eigen::MatrixXd rawSpatialCovariance;        // D D^T - formalResidualCovariance
    Eigen::MatrixXd spatialModelCovariance;      // positive spectral part + common-reference floor
};

namespace detail
{
inline bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0;
}

inline bool optionsValid(const Options& options)
{
    return options.minReferences >= 4 && finitePositive(options.maxCondition) &&
        finitePositive(options.maxDistanceKm) && finitePositive(options.maxWeightSquaredNorm) &&
        finitePositive(options.perSatelliteFloorMetres);
}

// Low-order OLS trend plus nearest-four inverse-distance-squared interpolation
// of its residuals, evaluated DIRECTLY at the user. This is NOT a gridded MPFG
// reproduction. The final row is w_plane + (I-H_plane)^T h_IDW, so a plane is
// reproduced and receiver GF offsets still cancel in each same-station SD.
// A three-station LOO design has zero OLS residual space: using all three for
// h is the exact plane-only limit. Ordinary prediction still requires m >= 4.
inline bool spatialWeights(const Eigen::MatrixXd& design, const Eigen::Vector3d& target,
    double maximumCondition, Eigen::VectorXd& weights, double& condition,
    Eigen::VectorXd* inverseDistanceWeights = nullptr)
{
    if (design.rows() < 3 || design.cols() != 3 || !design.allFinite() || !target.allFinite())
        return false;
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(design, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& singular = svd.singularValues();
    if (singular.size() != 3 || !singular.allFinite() || !(singular(0) > 0) ||
        singular(2) <= singular(0) * 1e-12)
        return false;
    condition = singular(0) / singular(2);
    if (!std::isfinite(condition) || condition > maximumCondition)
        return false;
    const Eigen::VectorXd plane = svd.matrixU() * singular.cwiseInverse().asDiagonal() *
        svd.matrixV().transpose() * target;
    std::vector<std::pair<double, int>> distances;
    for (int i = 0; i < design.rows(); ++i)
    {
        const double east = design(i, 1) - target(1);
        const double north = design(i, 2) - target(2);
        distances.emplace_back(east * east + north * north, i);
    }
    std::sort(distances.begin(), distances.end()); // deterministic index tie-break
    const int neighbours = std::min(4, static_cast<int>(distances.size()));
    Eigen::VectorXd h = Eigen::VectorXd::Zero(design.rows());
    if (distances.front().first == 0)
    {
        int coincident = 0;
        while (coincident < neighbours && distances[coincident].first == 0) ++coincident;
        for (int i = 0; i < coincident; ++i) h(distances[i].second) = 1.0 / coincident;
    }
    else
    {
        // Scaling by the nearest squared distance avoids overflow in 1/d^2.
        double normalization = 0;
        for (int i = 0; i < neighbours; ++i)
        {
            const double value = distances.front().first / distances[i].first;
            h(distances[i].second) = value;
            normalization += value;
        }
        h /= normalization;
    }
    weights = plane + h - svd.matrixU() * (svd.matrixU().transpose() * h);
    if (inverseDistanceWeights) *inverseDistanceWeights = h;
    if (!weights.allFinite()) return false;
    const double scale = std::max(1.0, target.cwiseAbs().maxCoeff());
    return (design.transpose() * weights - target).cwiseAbs().maxCoeff() <= 1e-8 * scale;
}

inline double cross(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
    const auto u = b - a;
    const auto v = c - a;
    return u.x() * v.y() - u.y() * v.x();
}

inline bool originInsideHull(std::vector<Eigen::Vector2d> points)
{
    std::sort(points.begin(), points.end(), [](const auto& left, const auto& right)
    {
        return left.x() < right.x() || (left.x() == right.x() && left.y() < right.y());
    });
    points.erase(std::unique(points.begin(), points.end(), [](const auto& a, const auto& b)
    { return a.x() == b.x() && a.y() == b.y(); }), points.end());
    if (points.size() < 3) return false;
    std::vector<Eigen::Vector2d> hull(2 * points.size());
    int count = 0;
    for (const auto& point : points)
    {
        while (count >= 2 && cross(hull[count - 2], hull[count - 1], point) <= 0) --count;
        hull[count++] = point;
    }
    const int lower = count;
    for (int i = static_cast<int>(points.size()) - 2; i >= 0; --i)
    {
        while (count > lower && cross(hull[count - 2], hull[count - 1], points[i]) <= 0) --count;
        hull[count++] = points[i];
    }
    if (count <= 3) return false;
    hull.resize(count - 1);
    double scale = 1;
    for (const auto& point : hull) scale = std::max(scale, point.squaredNorm());
    const double tolerance = 1e-10 * scale;
    for (std::size_t i = 0; i < hull.size(); ++i)
        if (cross(hull[i], hull[(i + 1) % hull.size()], Eigen::Vector2d::Zero()) < -tolerance)
            return false;
    return true;
}
} // namespace detail

inline Result buildSpatialPrediction(const Eigen::VectorXd& sourceMean,
    const Eigen::MatrixXd& sourceCovariance, const std::vector<PairSupport>& pairs,
    const Options& options = {})
{
    Result out;
    const int n = static_cast<int>(sourceMean.size());
    out.prediction.resize(0, n);
    out.spatialCovariance.resize(0, 0);
    if (!detail::optionsValid(options)) { out.reason = "INVALID_SPATIAL_OPTIONS"; return out; }
    if (n <= 0 || sourceCovariance.rows() != n || sourceCovariance.cols() != n ||
        !sourceMean.allFinite() || !sourceCovariance.allFinite())
    { out.reason = "INVALID_JOINT_SOURCE_DIMENSION_OR_VALUES"; return out; }
    const double covarianceScale = std::max(1.0, sourceCovariance.cwiseAbs().maxCoeff());
    if ((sourceCovariance - sourceCovariance.transpose()).cwiseAbs().maxCoeff() > 1e-10 * covarianceScale)
    { out.reason = "ASYMMETRIC_JOINT_SOURCE_COVARIANCE"; return out; }
    const Eigen::MatrixXd q = 0.5 * (sourceCovariance + sourceCovariance.transpose());
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> sourceEigen(q, Eigen::EigenvaluesOnly);
    if (sourceEigen.info() != Eigen::Success || !sourceEigen.eigenvalues().allFinite() ||
        sourceEigen.eigenvalues().minCoeff() < -1e-8 * covarianceScale)
    { out.reason = "INDEFINITE_JOINT_SOURCE_COVARIANCE"; return out; }

    struct Accepted
    {
        PairDiagnostic diagnostic;
        Eigen::RowVectorXd prediction;
        std::map<std::string, Eigen::RowVectorXd> residualRows;
    };
    std::vector<Accepted> accepted;
    std::map<std::string, Eigen::Vector2d> stationPositions;
    std::set<std::string> targets;
    std::string commonReference;
    for (const auto& pair : pairs)
    {
        PairDiagnostic d;
        d.satellite = pair.satellite;
        d.reference = pair.reference;
        auto reject = [&](const std::string& reason)
        { d.reason = reason; out.diagnostics.push_back(d); };
        if (pair.satellite.empty() || pair.reference.empty() || pair.satellite == pair.reference)
        { reject("INVALID_SATELLITE_PAIR"); continue; }
        if (!commonReference.empty() && commonReference != pair.reference)
        { reject("MIXED_COMMON_REFERENCE_NOT_SUPPORTED"); continue; }
        if (targets.count(pair.satellite)) { reject("DUPLICATE_TARGET_PAIR"); continue; }
        const int m = static_cast<int>(pair.stations.size());
        if (m < options.minReferences) { reject("INSUFFICIENT_COMMON_REFERENCES"); continue; }
        Eigen::MatrixXd design(m, 3), differences = Eigen::MatrixXd::Zero(m, n);
        std::set<std::string> names;
        std::vector<Eigen::Vector2d> points;
        std::string inputFailure;
        for (int i = 0; i < m; ++i)
        {
            const auto& station = pair.stations[i];
            if (station.station.empty() || !names.insert(station.station).second)
            { inputFailure = "DUPLICATE_OR_EMPTY_REFERENCE_STATION"; break; }
            if (!std::isfinite(station.eastKm) || !std::isfinite(station.northKm) ||
                station.satelliteSourceIndex < 0 || station.satelliteSourceIndex >= n ||
                station.referenceSourceIndex < 0 || station.referenceSourceIndex >= n ||
                station.referenceSourceIndex == station.satelliteSourceIndex)
            { inputFailure = "INVALID_REFERENCE_COORDINATE_OR_SOURCE_INDEX"; break; }
            const Eigen::Vector2d position(station.eastKm, station.northKm);
            const auto previous = stationPositions.find(station.station);
            if (previous != stationPositions.end() &&
                (previous->second - position).cwiseAbs().maxCoeff() > 1e-8)
            { inputFailure = "INCONSISTENT_SHARED_REFERENCE_COORDINATES"; break; }
            d.stations.push_back(station.station);
            points.push_back(position);
            d.maxDistanceKm = std::max(d.maxDistanceKm, position.norm());
            design.row(i) << 1, station.eastKm, station.northKm;
            differences(i, station.satelliteSourceIndex) = 1;
            differences(i, station.referenceSourceIndex) = -1;
        }
        if (!inputFailure.empty()) { reject(inputFailure); continue; }
        if (d.maxDistanceKm > options.maxDistanceKm) { reject("REFERENCE_DISTANCE_LIMIT"); continue; }
        d.insideConvexHull = detail::originInsideHull(points);
        if (!detail::spatialWeights(design, Eigen::Vector3d(1, 0, 0), options.maxCondition,
                d.weights, d.geometryCondition, &d.inverseDistanceWeights))
        { reject("RANK_DEFICIENT_OR_ILL_CONDITIONED_GEOMETRY"); continue; }
        d.weightSquaredNorm = d.weights.squaredNorm();
        if (d.weightSquaredNorm > options.maxWeightSquaredNorm)
        { reject("USER_PREDICTION_LEVERAGE_LIMIT"); continue; }

        Accepted current;
        d.leaveOneOutResiduals.resize(m);
        d.rawLeaveOneOutResiduals.resize(m);
        d.leaveOneOutWeightSquaredNorms.resize(m);
        d.leaveOneOutWeights = Eigen::MatrixXd::Zero(m, m);
        d.leaveOneOutInverseDistanceWeights = Eigen::MatrixXd::Zero(m, m);
        bool leaveOneOutValid = true;
        for (int omitted = 0; omitted < m; ++omitted)
        {
            Eigen::MatrixXd reduced(m - 1, 3);
            std::vector<int> retained;
            for (int i = 0; i < m; ++i) if (i != omitted)
            { reduced.row(retained.size()) = design.row(i); retained.push_back(i); }
            Eigen::VectorXd weights;
            Eigen::VectorXd inverseDistanceWeights;
            double condition = 0;
            if (!detail::spatialWeights(reduced, design.row(omitted).transpose(), options.maxCondition,
                    weights, condition, &inverseDistanceWeights))
            { leaveOneOutValid = false; break; }
            d.maximumLeaveOneOutCondition = std::max(d.maximumLeaveOneOutCondition, condition);
            Eigen::RowVectorXd residual = differences.row(omitted);
            for (int j = 0; j < m - 1; ++j)
            {
                residual -= weights(j) * differences.row(retained[j]);
                d.leaveOneOutWeights(omitted, retained[j]) = weights(j);
                d.leaveOneOutInverseDistanceWeights(omitted, retained[j]) = inverseDistanceWeights(j);
            }
            d.rawLeaveOneOutResiduals(omitted) = residual.dot(sourceMean);
            d.leaveOneOutWeightSquaredNorms(omitted) = weights.squaredNorm();
            const double normalization = std::sqrt(double(m) * (1 + weights.squaredNorm()));
            residual /= normalization;
            d.leaveOneOutResiduals(omitted) = residual.dot(sourceMean);
            current.residualRows.emplace(pair.stations[omitted].station, residual);
        }
        if (!leaveOneOutValid) { reject("LEAVE_ONE_OUT_GEOMETRY_UNSUPPORTED"); continue; }
        current.prediction = d.weights.transpose() * differences;
        d.accepted = true;
        d.reason = "ACCEPTED_EXPERIMENTAL_SD_IONO_SPATIAL_MODEL";
        current.diagnostic = d;
        accepted.push_back(std::move(current));
        out.diagnostics.push_back(d);
        targets.insert(pair.satellite);
        commonReference = pair.reference;
        for (const auto& station : pair.stations)
            stationPositions[station.station] = {station.eastKm, station.northKm};
    }
    const int k = static_cast<int>(accepted.size());
    if (k == 0) { out.reason = "NO_SUPPORTED_SPATIAL_PAIRS"; return out; }
    for (const auto& item : stationPositions) out.globalStations.push_back(item.first);
    const int g = static_cast<int>(out.globalStations.size());
    out.prediction.resize(k, n);
    out.stationWeights = Eigen::MatrixXd::Zero(k, g);
    out.residualMatrix = Eigen::MatrixXd::Zero(k, g);
    out.formalResidualCovariance = Eigen::MatrixXd::Zero(k, k);
    for (int i = 0; i < k; ++i)
    {
        out.prediction.row(i) = accepted[i].prediction;
        out.acceptedPairs.push_back(accepted[i].diagnostic);
        for (int j = 0; j < static_cast<int>(accepted[i].diagnostic.stations.size()); ++j)
        {
            const auto where = std::lower_bound(out.globalStations.begin(), out.globalStations.end(),
                accepted[i].diagnostic.stations[j]);
            out.stationWeights(i, where - out.globalStations.begin()) = accepted[i].diagnostic.weights(j);
        }
    }
    for (int station = 0; station < g; ++station)
    {
        Eigen::MatrixXd a = Eigen::MatrixXd::Zero(k, n);
        for (int i = 0; i < k; ++i)
        {
            const auto row = accepted[i].residualRows.find(out.globalStations[station]);
            if (row == accepted[i].residualRows.end()) continue;
            a.row(i) = row->second;
            out.residualMatrix(i, station) = a.row(i).dot(sourceMean);
        }
        out.formalResidualCovariance.noalias() += a * q * a.transpose();
    }
    out.rawSpatialCovariance = out.residualMatrix * out.residualMatrix.transpose() -
        out.formalResidualCovariance;
    out.rawSpatialCovariance = 0.5 * (out.rawSpatialCovariance + out.rawSpatialCovariance.transpose()).eval();
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> residualEigen(out.rawSpatialCovariance);
    if (residualEigen.info() != Eigen::Success || !residualEigen.eigenvalues().allFinite())
    { out.reason = "SPATIAL_RESIDUAL_EIGEN_FAILURE"; return out; }
    out.spatialModelCovariance = residualEigen.eigenvectors() *
        residualEigen.eigenvalues().cwiseMax(0.0).asDiagonal() * residualEigen.eigenvectors().transpose();
    // A separate 5 cm L1 mismatch per satellite induces shared-reference noise
    // sigma^2 (I + 11^T) on SDs. It is NEVER a measured product covariance.
    out.spatialModelCovariance += std::pow(options.perSatelliteFloorMetres, 2) *
        (Eigen::MatrixXd::Identity(k, k) + Eigen::MatrixXd::Ones(k, k));
    // Explicit IID spatial-field discrepancy assumption: add prediction-weight
    // leverage, retaining cross-pair terms and positive semidefiniteness.
    const Eigen::MatrixXd propagation = Eigen::MatrixXd::Ones(k, k) +
        out.stationWeights * out.stationWeights.transpose();
    out.spatialCovariance = out.spatialModelCovariance.cwiseProduct(propagation);
    out.spatialCovariance = 0.5 * (out.spatialCovariance + out.spatialCovariance.transpose()).eval();
    if (!out.prediction.allFinite() || !out.spatialCovariance.allFinite())
    { out.reason = "NONFINITE_SPATIAL_PREDICTION"; return out; }
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spatialEigen(out.spatialCovariance, Eigen::EigenvaluesOnly);
    if (spatialEigen.info() != Eigen::Success || spatialEigen.eigenvalues().minCoeff() <= 0)
    { out.reason = "SPATIAL_COVARIANCE_NOT_POSITIVE_DEFINITE"; return out; }
    out.valid = true;
    out.reason = "EXPERIMENTAL_IID_SPATIAL_MISMATCH_NOT_CALIBRATED_USER_ACCURACY";
    return out;
}
} // namespace zhangregional
