// Read-only authenticated checkpoint inventory; performs no estimation.
#include "common/zhangCheckpoint.hpp"
#include <boost/json.hpp>
#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

static void require(bool pass, const std::string& why)
{ if (!pass) throw std::runtime_error(why); }
static boost::json::value finiteJson(double value)
{ return std::isfinite(value) ? boost::json::value(value) : boost::json::value(nullptr); }

int main(int argc, char** argv)
try {
    require(argc == 5, "Usage: checkpoint.bundle producer-sha256 expected-mode output.json");
    const std::filesystem::path path(argv[1]);
    std::ifstream metadata(path.parent_path() / "checkpoint_manifest.json");
    std::ostringstream text; text << metadata.rdbuf();
    const auto document = boost::json::parse(text.str());
    const auto& doc = document.as_object();
    auto field = [&](const char* name) { return std::string(doc.at(name).as_string().c_str()); };
    ZhangCheckpointExpectations expected;
    expected.experimentMode = argv[3]; expected.binarySha256 = argv[2];
    expected.configSha256 = zhangCheckpointSha256(field("config_text"));
    expected.inputManifestSha256 = zhangCheckpointSha256(field("input_manifest_text"));
    expected.platformFingerprint = field("platform");
    expected.compilerFingerprint = field("compiler");
    expected.linearAlgebraFingerprint = field("linear_algebra");
    expected.endianness = field("endianness");
    require(field("experiment_mode") == argv[3] && field("binary_sha256") == argv[2],
        "WRONG_CHECKPOINT_PRODUCER_OR_MODE");
    require(expected.configSha256 == field("config_sha256") &&
        expected.inputManifestSha256 == field("input_manifest_sha256"), "MANIFEST_HASH_MISMATCH");
    const auto envelope = verifyZhangCheckpointBundleEnvelope(path.string());
    require(envelope.valid, envelope.failureReason);
    ZhangCheckpointBundle bundle;
    const auto decoded = readZhangCheckpointBundle(path.string(), expected, bundle);
    require(decoded.valid, decoded.failureReason);
    const auto& core = bundle.kfCore;
    require(core.P.rows() == core.x.size() && core.P.cols() == core.x.size(), "P_SHAPE_MISMATCH");
    std::map<std::string, std::array<int, 3>> coordinates;
    boost::json::object stateCounts, positions;
    unsigned rateCount = 0;
    for (const auto& [key, index] : core.kfIndexMap) {
        require(index >= 0 && index < core.x.size(), "INVALID_STATE_INDEX");
        const std::string system = key.satellite.sys == E_Sys::NONE ? "SHARED_NONE" :
            enum_to_string(key.satellite.sys);
        const std::string name = system + "/" + enum_to_string(key.type);
        auto* count = stateCounts.if_contains(name);
        stateCounts[name] = count ? count->as_int64() + 1 : 1;
        if (key.type == KF::REC_POS_RATE) ++rateCount;
        if (key.type != KF::REC_POS) continue;
        require(key.number >= 0 && key.number < 3, "INVALID_REC_POS_AXIS");
        auto [it, inserted] = coordinates.try_emplace(key.receiver, std::array<int, 3>{-1, -1, -1});
        require(it->second[key.number] == -1, "DUPLICATE_REC_POS_AXIS");
        it->second[key.number] = index;
    }
    for (const auto& [receiver, indices] : coordinates) {
        boost::json::array xyz, variance, sigma, covariance;
        bool complete = true, valid = true;
        Matrix3d q = Matrix3d::Zero();
        for (int axis = 0; axis < 3; ++axis) {
            const int i = indices[axis];
            complete &= i >= 0;
            const double x = i >= 0 ? core.x(i) : std::numeric_limits<double>::quiet_NaN();
            const double v = i >= 0 ? core.P(i, i) : std::numeric_limits<double>::quiet_NaN();
            xyz.push_back(finiteJson(x)); variance.push_back(finiteJson(v));
            sigma.push_back(finiteJson(v >= 0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN()));
            valid &= std::isfinite(x) && std::isfinite(v) && v > 0;
            for (int other = 0; other < 3; ++other) {
                const double value = i >= 0 && indices[other] >= 0 ? core.P(i, indices[other]) :
                    std::numeric_limits<double>::quiet_NaN();
                covariance.push_back(finiteJson(value)); q(axis, other) = value;
            }
        }
        double minimum = std::numeric_limits<double>::quiet_NaN();
        double radius = std::numeric_limits<double>::quiet_NaN();
        if (complete && q.allFinite()) {
            Eigen::SelfAdjointEigenSolver<Matrix3d> solver((q + q.transpose()) / 2);
            if (solver.info() == Eigen::Success) minimum = solver.eigenvalues().minCoeff();
            radius = std::hypot(core.x(indices[0]), core.x(indices[1]), core.x(indices[2]));
        }
        positions[receiver] = boost::json::object{{"xyz_m", xyz}, {"variance_xyz_m2", variance},
            {"sigma_xyz_m", sigma}, {"covariance_xyz_row_major_m2", covariance},
            {"complete_xyz", complete}, {"finite_xyz_positive_variance", valid},
            {"absolute_ecef_radius_m", finiteJson(radius)},
            {"absolute_ecef_plausible", std::isfinite(radius) && radius > 3e6 && radius < 8e6},
            {"coordinate_block_min_eigenvalue_m2", finiteJson(minimum)},
            {"coordinate_block_positive_definite", std::isfinite(minimum) && minimum > 0}};
    }
    boost::json::object report{{"schema", "CHECKPOINT_COORDINATE_INVENTORY_V1"},
        {"bundle_sha256", zhangCheckpointFileSha256(path.string())},
        {"producer_sha256", argv[2]}, {"experiment_mode", argv[3]},
        {"checkpoint_epoch", bundle.manifest.epoch}, {"runtime_id", bundle.manifest.runtimeId},
        {"state_dimension", core.x.size()}, {"P_entries", core.P.size()},
        {"full_x_finite", core.x.allFinite()}, {"full_P_finite", core.P.allFinite()},
        {"full_P_maximum_asymmetry_mixed_native_units", core.P.size() ?
            finiteJson((core.P - core.P.transpose()).cwiseAbs().maxCoeff()) : boost::json::value(nullptr)},
        {"full_P_minimum_diagonal_mixed_native_units", core.P.rows() ?
            finiteJson(core.P.diagonal().minCoeff()) : boost::json::value(nullptr)},
        {"P_validation_scope", "Full matrix finite, symmetry and diagonal metrics; per-station XYZ block positive-definite. No full-P PSD claim"},
        {"coordinate_state_semantics", "Production ppp_obs.cpp seeds absolute ECEF rRec and uses the existing REC_POS directly in geometry; radius independently audited"},
        {"state_counts_by_system_and_type", stateCounts}, {"receiver_position_rate_count", rateCount},
        {"coordinates", positions}, {"scientific_accuracy_accepted", false}};
    std::ofstream output(argv[4]); output << boost::json::serialize(report) << '\n';
    require(bool(output), "REPORT_WRITE_FAILED");
    std::cout << "CHECKPOINT_COORDINATE_INVENTORY receivers=" << coordinates.size()
              << " dimension=" << core.x.size() << " no_estimation=1\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "CHECKPOINT_COORDINATE_AUDIT_FAILED " << error.what() << '\n'; return 1;
}
