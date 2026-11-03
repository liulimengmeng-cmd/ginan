#include "common/zhangUserRegionalSpatial.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
void require(bool condition, const std::string& reason)
{
    if (!condition) throw std::runtime_error(reason);
}

void near(const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected, double tolerance,
    const std::string& reason)
{
    require(actual.rows() == expected.rows() && actual.cols() == expected.cols(), reason + ": dimensions");
    require((actual - expected).cwiseAbs().maxCoeff() <= tolerance, reason);
}

std::vector<zhangregional::PairSupport> supports()
{
    // Source order: one clock, two phases, then station-major [I_A,I_B,I_Q].
    const std::vector<Eigen::Vector2d> positions{{-10, -10}, {10, -10}, {10, 10}, {-10, 10}, {0, 0}};
    std::vector<zhangregional::PairSupport> result;
    for (int target = 0; target < 2; ++target)
    {
        zhangregional::PairSupport pair;
        pair.satellite = target == 0 ? "G12" : "G25";
        pair.reference = "G31";
        for (int station = 0; station < 5; ++station)
            pair.stations.push_back({"R" + std::to_string(station), positions[station].x(),
                positions[station].y(), 3 + 3 * station + target, 5 + 3 * station});
        result.push_back(pair);
    }
    return result;
}

Eigen::VectorXd planeMeans()
{
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(18);
    mean.head(3) << 12, -3, 7;
    const auto input = supports();
    for (const auto& station : input.front().stations)
    {
        const int ordinal = (station.satelliteSourceIndex - 3) / 3;
        const double receiverGf = 10 + 3 * ordinal;
        const double reference = receiverGf + 2 + 0.02 * station.eastKm - 0.01 * station.northKm;
        mean(5 + 3 * ordinal) = reference;
        mean(3 + 3 * ordinal) = reference + 4 + 0.03 * station.eastKm + 0.05 * station.northKm;
        mean(4 + 3 * ordinal) = reference - 2 - 0.04 * station.eastKm + 0.02 * station.northKm;
    }
    return mean;
}

Eigen::MatrixXd mixedCovariance()
{
    Eigen::MatrixXd root = Eigen::MatrixXd::Identity(18, 18) * 0.02;
    for (int row = 0; row < 18; ++row)
        for (int col = 0; col <= row; ++col)
            root(row, col) += 0.001 * (1 + ((row + 2 * col) % 5));
    return root * root.transpose();
}

// Independent implementation for tests: explicit normal equations and hat
// matrix, rather than production SVD or its returned prediction/LOO weights.
Eigen::VectorXd independentTrendIdwWeights(const Eigen::MatrixXd& x,
    const Eigen::Vector3d& target)
{
    const Eigen::Matrix3d normalInverse = (x.transpose() * x).inverse();
    const Eigen::MatrixXd hat = x * normalInverse * x.transpose();
    const Eigen::VectorXd trend = x * normalInverse * target;
    std::vector<std::pair<double, int>> order;
    for (int i = 0; i < x.rows(); ++i)
        order.emplace_back(std::pow(x(i, 1) - target(1), 2) +
            std::pow(x(i, 2) - target(2), 2), i);
    std::sort(order.begin(), order.end());
    const int count = std::min(4, static_cast<int>(order.size()));
    Eigen::VectorXd idw = Eigen::VectorXd::Zero(x.rows());
    if (order.front().first == 0)
    {
        int coincident = 0;
        while (coincident < count && order[coincident].first == 0) ++coincident;
        for (int i = 0; i < coincident; ++i) idw(order[i].second) = 1.0 / coincident;
    }
    else
    {
        for (int i = 0; i < count; ++i) idw(order[i].second) = 1.0 / order[i].first;
        idw /= idw.sum();
    }
    return trend + (Eigen::MatrixXd::Identity(x.rows(), x.rows()) - hat).transpose() * idw;
}

void testPlaneAndGf()
{
    const auto input = supports();
    const Eigen::VectorXd mean = planeMeans();
    const auto result = zhangregional::buildSpatialPrediction(mean, mixedCovariance(), input);
    require(result.valid && result.acceptedPairs.size() == 2, "plane support accepted");
    Eigen::Vector2d expected(4, -2);
    near(result.prediction * mean, expected, 1e-11, "affine plane reproduced at user origin");
    for (const auto& pair : result.acceptedPairs)
    {
        require(pair.insideConvexHull, "user inside square hull");
        require(std::abs(pair.weights.sum() - 1) < 1e-12, "partition of unity");
        require(pair.rawLeaveOneOutResiduals.cwiseAbs().maxCoeff() < 1e-11, "plane LOO residual zero");
    }
    require(result.model.find("NOT_MPFG") != std::string::npos, "experimental non-MPFG model labelled");
    // Analytic trend-plus-IDW weights at an observed centre are interpolatory:
    // w_trend + (I-H)^T e_centre = e_centre, not the old OLS-only 1/5 weights.
    Eigen::MatrixXd expectedW = Eigen::MatrixXd::Zero(2, 18);
    expectedW(0, 15) = 1;
    expectedW(1, 16) = 1;
    expectedW(0, 17) = expectedW(1, 17) = -1;
    near(result.prediction, expectedW, 1e-12, "independent analytic prediction weights");
    Eigen::VectorXd receiverDatum = Eigen::VectorXd::Zero(18);
    for (int station = 0; station < 5; ++station)
        receiverDatum.segment(3 + 3 * station, 3).setConstant(100 * (station + 1));
    for (int target = 0; target < 2; ++target)
        for (int station = 0; station < 5; ++station)
            require(result.prediction(target, 3 + 3 * station + target) +
                result.prediction(target, 5 + 3 * station) == 0,
                "receiver GF coefficients cancel exactly in W");
    near(result.prediction * receiverDatum, Eigen::Vector2d::Zero(), 1e-12, "receiver GF floating product cancels");
    const auto shifted = zhangregional::buildSpatialPrediction(mean + receiverDatum, mixedCovariance(), input);
    require(shifted.valid, "GF shifted support accepted");
    near(shifted.prediction * (mean + receiverDatum), expected, 1e-10, "receiver GF invariance");
    near(shifted.spatialCovariance, result.spatialCovariance, 1e-12, "spatial uncertainty GF invariant");
}

void testFullCovarianceAndFloor()
{
    const Eigen::MatrixXd q = mixedCovariance();
    const auto result = zhangregional::buildSpatialPrediction(planeMeans(), q, supports());
    require(result.valid, "mixed full covariance support accepted");
    // Clock/phase product errors correlate with the same IONO source variables.
    Eigen::MatrixXd observation = Eigen::MatrixXd::Zero(2, 18);
    observation(0, 0) = 1; observation(0, 1) = -1;
    observation(1, 0) = 1; observation(1, 2) = -1;
    const Eigen::MatrixXd cross = observation * q * result.prediction.transpose();
    require(cross.cwiseAbs().maxCoeff() > 1e-7, "nonzero clock-phase/iono cross covariance retained");
    Eigen::MatrixXd combined(4, 18);
    combined.topRows(2) = observation;
    combined.bottomRows(2) = result.prediction;
    Eigen::MatrixXd joint = combined * q * combined.transpose();
    joint.bottomRightCorner(2, 2) += result.spatialCovariance;
    require(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(joint).eigenvalues().minCoeff() > -1e-12,
        "joint clock phase ionosphere covariance PSD");
    near(joint.topRightCorner(2, 2), cross, 1e-14, "joint cross block exact");
    // Plane residuals vanish. Removing formal residual variance has no positive
    // eigenvalues, so the independent result is solely the prescribed floor.
    Eigen::Matrix2d expectedFloor;
    expectedFloor << 0.01, 0.005, 0.005, 0.01; // .05^2*(I+11T)*(1+sum_r e_centre^2)
    near(result.spatialCovariance, expectedFloor, 1e-12, "shared reference floor and weight propagation");
    require(result.spatialCovariance(0, 1) > 0, "shared reference spatial off-diagonal not dropped");
}

void testResidualMathIndependently()
{
    Eigen::VectorXd mean = planeMeans();
    mean(3) += 1.2;
    mean(10) -= 0.7;
    const auto input = supports();
    const Eigen::MatrixXd q = mixedCovariance();
    const auto result = zhangregional::buildSpatialPrediction(mean, q, input);
    require(result.valid, "non-plane data accepted as uncertain spatial field");
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(2, 5);
    Eigen::MatrixXd e = Eigen::MatrixXd::Zero(2, 2);
    for (int omitted = 0; omitted < 5; ++omitted)
    {
        Eigen::MatrixXd a = Eigen::MatrixXd::Zero(2, 18);
        for (int target = 0; target < 2; ++target)
        {
            Eigen::MatrixXd x(4, 3);
            std::vector<int> retained;
            for (int station = 0; station < 5; ++station) if (station != omitted)
            {
                const auto& point = input[target].stations[station];
                x.row(retained.size()) << 1, point.eastKm, point.northKm;
                retained.push_back(station);
            }
            const auto& point = input[target].stations[omitted];
            const Eigen::Vector3d evaluation(1, point.eastKm, point.northKm);
            const Eigen::VectorXd w = independentTrendIdwWeights(x, evaluation);
            a(target, point.satelliteSourceIndex) = 1;
            a(target, point.referenceSourceIndex) = -1;
            for (int j = 0; j < 4; ++j)
            {
                const auto& ref = input[target].stations[retained[j]];
                a(target, ref.satelliteSourceIndex) -= w(j);
                a(target, ref.referenceSourceIndex) += w(j);
            }
            a.row(target) /= std::sqrt(5.0 * (1 + w.squaredNorm()));
            d(target, omitted) = a.row(target).dot(mean);
        }
        e += a * q * a.transpose();
    }
    near(result.residualMatrix, d, 1e-12, "LOO residual map independent trend plus IDW normal equations");
    near(result.formalResidualCovariance, e, 1e-12, "all formal cross-pair residual terms");
    near(result.rawSpatialCovariance, d * d.transpose() - e, 1e-12, "formal variance subtracted once");
    require(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(result.spatialCovariance)
        .eigenvalues().minCoeff() > 0, "positive spatial covariance after spectral projection");
}

void testHandCalculatedNonPlaneIdw()
{
    auto input = supports();
    input.resize(1);
    input.front().stations.resize(4);
    // Four corners of [-1,+1]^2, evaluated at (1/2,1/2). Coordinates supplied
    // to production are relative to that user, so squared distances are
    // [9/2,5/2,1/2,5/2] and h_IDW = [5,9,45,9]/68. Checkerboard data
    // [1,-1,1,-1] is orthogonal to the OLS plane. Thus the prediction is
    // (5-9+45-9)/68 = 8/17, not zero (the OLS-only prediction).
    const std::vector<Eigen::Vector2d> position{{-1.5, -1.5}, {0.5, -1.5},
        {0.5, 0.5}, {-1.5, 0.5}};
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(18);
    for (int i = 0; i < 4; ++i)
    {
        auto& station = input.front().stations[i];
        station.eastKm = position[i].x(); station.northKm = position[i].y();
        mean(station.satelliteSourceIndex) = i % 2 == 0 ? 1 : -1;
    }
    const auto four = zhangregional::buildSpatialPrediction(mean, mixedCovariance(), input);
    require(four.valid, "hand calculated four-reference non-plane case accepted");
    Eigen::Vector4d expectedIdw; expectedIdw << 5, 9, 45, 9; expectedIdw /= 68;
    Eigen::Vector4d expectedWeights; expectedWeights << 8, 9, 42, 9; expectedWeights /= 68;
    near(four.acceptedPairs.front().inverseDistanceWeights, expectedIdw, 1e-12,
        "four-reference inverse squared distances hand calculation");
    near(four.acceptedPairs.front().weights, expectedWeights, 1e-12,
        "four-reference trend residual weights hand calculation");
    require(std::abs((four.prediction * mean)(0) - 8.0 / 17) < 1e-12,
        "four-reference non-plane residual contribution not discarded");

    // Five references: a non-plane unit spike at the centre. The centre is
    // coincident with the user, so IDW exactly reproduces 1, whereas pure OLS
    // returns 1/5. This also exercises the zero-distance branch without 1/0.
    input = supports(); input.resize(1);
    mean.setZero(); mean(15) = 1;
    const auto five = zhangregional::buildSpatialPrediction(mean, mixedCovariance(), input);
    require(five.valid && std::abs((five.prediction * mean)(0) - 1) < 1e-12,
        "five-reference centre spike interpolated exactly");
    require(five.acceptedPairs.front().inverseDistanceWeights(4) == 1,
        "coincident user chooses the coincident reference without division by zero");
    const auto& diagnostic = five.acceptedPairs.front();
    // Independent LOO reference omitting the centre: the remaining corner
    // values are zero, so prediction is zero and raw leave-one-out error is 1.
    require(std::abs(diagnostic.rawLeaveOneOutResiduals(4) - 1) < 1e-12,
        "centre spike independent hand-calculated LOO residual");
}

void testRejections()
{
    const Eigen::VectorXd mean = planeMeans();
    const Eigen::MatrixXd q = mixedCovariance();
    auto check = [&](std::vector<zhangregional::PairSupport> pairs, const std::string& reason)
    {
        const auto result = zhangregional::buildSpatialPrediction(mean, q, pairs);
        require(!result.valid && result.acceptedPairs.empty(), reason);
    };
    auto pairs = supports(); pairs.resize(1); pairs[0].stations.resize(3);
    check(pairs, "fewer than four references rejected");
    pairs = supports(); pairs.resize(1);
    for (int i = 0; i < 5; ++i) {pairs[0].stations[i].eastKm = i; pairs[0].stations[i].northKm = 0;}
    check(pairs, "collinear geometry rejected");
    for (int i = 0; i < 5; ++i) pairs[0].stations[i].northKm = (i % 2) * 1e-8;
    check(pairs, "ill conditioned geometry rejected");
    pairs = supports(); pairs.resize(1); pairs[0].stations.resize(4);
    for (int i = 0; i < 4; ++i) {pairs[0].stations[i].eastKm = i - 1; pairs[0].stations[i].northKm = i == 3 ? 2 : 0;}
    const auto loo = zhangregional::buildSpatialPrediction(mean, q, pairs);
    require(!loo.valid && loo.diagnostics.front().reason == "LEAVE_ONE_OUT_GEOMETRY_UNSUPPORTED",
        "full rank design with singular LOO rejected");
    pairs = supports(); pairs.resize(1); pairs[0].stations[0].eastKm = 151;
    check(pairs, "distance limit rejected");
    pairs = supports(); pairs.resize(1);
    for (auto& point : pairs[0].stations) point.eastKm += 100;
    check(pairs, "large extrapolation leverage rejected");
    pairs = supports(); pairs.push_back(pairs.front());
    const auto duplicate = zhangregional::buildSpatialPrediction(mean, q, pairs);
    require(duplicate.valid && duplicate.acceptedPairs.size() == 2 &&
        duplicate.diagnostics.back().reason == "DUPLICATE_TARGET_PAIR", "redundant target not double counted");
    pairs = supports(); pairs[1].reference = "G19";
    const auto mixed = zhangregional::buildSpatialPrediction(mean, q, pairs);
    require(mixed.valid && mixed.acceptedPairs.size() == 1 &&
        mixed.diagnostics.back().reason == "MIXED_COMMON_REFERENCE_NOT_SUPPORTED", "common reference contract enforced");
    zhangregional::Options options; options.perSatelliteFloorMetres = 0;
    require(!zhangregional::buildSpatialPrediction(mean, q, supports(), options).valid,
        "zero spatial mismatch floor forbidden");
    Eigen::MatrixXd indefinite = q; indefinite(0, 0) = -1;
    require(!zhangregional::buildSpatialPrediction(mean, indefinite, supports()).valid,
        "invalid joint covariance rejected");
}

void testDifferentReferenceSetsAndHull()
{
    auto pairs = supports();
    pairs[1].stations.erase(pairs[1].stations.begin() + 4);
    const auto different = zhangregional::buildSpatialPrediction(planeMeans(), mixedCovariance(), pairs);
    require(different.valid && different.acceptedPairs.size() == 2, "pair-specific common station sets allowed");
    near(different.prediction * planeMeans(), Eigen::Vector2d(4, -2), 1e-11,
        "different station sets preserve plane and GF cancellation");
    require(different.stationWeights(1, 4) == 0, "absent station weight zero");
    pairs = supports(); pairs.resize(1);
    for (auto& station : pairs[0].stations) station.eastKm += 12;
    const auto exterior = zhangregional::buildSpatialPrediction(planeMeans(), mixedCovariance(), pairs);
    require(exterior.valid && !exterior.acceptedPairs.front().insideConvexHull &&
        exterior.acceptedPairs.front().weightSquaredNorm <= 4,
        "bounded outside-hull prediction explicitly marked");
}
} // namespace

int main()
{
    try
    {
        testPlaneAndGf();
        testFullCovarianceAndFloor();
        testResidualMathIndependently();
        testHandCalculatedNonPlaneIdw();
        testRejections();
        testDifferentReferenceSetsAndHull();
        std::cout << "PASS user_regional_spatial_contract: plane-plus-4NN-IDW2-residual-not-MPFG, GF, full-Q, independent-LOO, hand-calculated-4/5-reference, PSD, shared-reference, safety-gates\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL user_regional_spatial_contract: " << error.what() << '\n';
        return 1;
    }
}
