#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"
#include "pea/zhangPppAr.hpp"
#include "pea/zhangUserRegional.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

// Linked to actual PEA objects with --wrap=main. This tests the public regional
// runtime and actual KFMeas covariance assembly, not copied adapter internals.
// It does not claim to exercise pppSatClocks/the complete observation builder;
// that integration boundary requires the separately executed real-data pilot.
namespace
{
namespace fs = std::filesystem;
constexpr double alpha1 = 40.3e16 / (1575.42e6 * 1575.42e6);
const std::string user = "REGIONAL_UNIT_USER";
const std::vector<std::string> satellites{"G01", "G02", "G03"};

void require(bool condition, const std::string& reason)
{
    if (!condition) throw std::runtime_error(reason);
}

void near(const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected,
    double tolerance, const std::string& reason)
{
    require(actual.rows() == expected.rows() && actual.cols() == expected.cols(), reason + ":dimensions");
    require((actual - expected).cwiseAbs().maxCoeff() < tolerance, reason);
}

void write(const fs::path& path, const std::string& contents)
{
    std::ofstream output(path, std::ios::binary);
    output << contents;
    output.close();
    require(bool(output), "fixture write failed: " + path.string());
}

void environment(const char* name, const std::string& value)
{
#ifdef _WIN32
    require(_putenv_s(name, value.c_str()) == 0, "environment setup failed");
#else
    require(setenv(name, value.c_str(), 1) == 0, "environment setup failed");
#endif
}

struct Fixture
{
    fs::path directory;
    std::vector<GTime> times;
    Eigen::MatrixXd covariance;
    Eigen::VectorXd mean;
    std::vector<Eigen::Vector2d> offsets{{-12, -8}, {11, -7}, {9, 14}, {-8, 11}, {2, 3}};

    Fixture()
    {
        directory = fs::current_path() / ("user_regional_joint_fixture_" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        require(fs::create_directory(directory), "fresh isolated fixture directory required");
        fs::create_directory(directory / "bundles");
        GTime start;
        require(str2time("2024 07 17 02 00 00", 0, 19, start, E_TimeSys::GPST) == 0, "fixture GPST parse");
        for (int i = 0; i < 9; ++i) { auto time = start; time.bigTime += 30 * i; times.push_back(time); }
        mean = Eigen::VectorXd::Zero(24);
        for (int sat = 0; sat < 3; ++sat)
        {
            mean(3 * sat) = 10 + 2 * sat;
            mean(3 * sat + 1) = 1 + 0.1 * sat;
            mean(3 * sat + 2) = 2 + 0.2 * sat;
        }
        for (int station = 0; station < 5; ++station)
            for (int sat = 0; sat < 3; ++sat)
                mean(9 + 3 * station + sat) = 20 * station + 0.5 * sat +
                    (0.01 + sat * 0.003) * offsets[station].x() +
                    (0.02 - sat * 0.002) * offsets[station].y();
        Eigen::MatrixXd root = Eigen::MatrixXd::Identity(24, 24) * 0.15;
        for (int row = 0; row < 24; ++row)
            for (int col = 0; col <= row; ++col)
                root(row, col) += 0.002 * (1 + ((3 * row + col) % 7));
        covariance = root * root.transpose(); // independent known SPD source Q

        std::ostringstream geometry;
        geometry << std::setprecision(17) << "receiver,role,x_m,y_m,z_m\n";
        for (int i = 0; i < 5; ++i)
            geometry << "R" << i << ",REFERENCE,6378137," << offsets[i].x() * 1000 << ',' << offsets[i].y() * 1000 << '\n';
        // Intentionally wrong USER values: runtime MUST NOT consume them.
        geometry << user << ",USER,1234,5678,9012\n";
        write(directory / "geometry.csv", geometry.str());
        std::ostringstream products;
        products << std::setprecision(17) << "synthetic legacy33 current products for actual production reader\n";
        for (int epoch : {0, 1, 2, 3, 4, 5, 7, 8})
        {
            makeBundle(epoch);
            for (int sat = 0; sat < 3; ++sat)
                for (int band = 1; band <= 2; ++band)
                {
                    const auto seconds = static_cast<long long>(times[epoch].bigTime);
                    const double clock = mean(3 * sat), phase = mean(3 * sat + band);
                    // Legacy 33-column products intentionally carry no invented
                    // continuity metadata. Current means/PPP authority suffice
                    // for this covariance-only runtime test.
                    products << seconds << ",PRODUCT_FIXED," << satellites[sat] << ',' << (band == 1 ? "L1C" : "L2W")
                        << ',' << clock << ",0.3," << phase << ",0.3,0," << clock - phase
                        << ",0.4,0,0,0,0," << seconds << ",0,initial,1,CURRENT_ALIGNMENT_VALID,1,1,1,1,UNIT_COMPONENT,UNIT_DATUM,"
                        << seconds << ',' << seconds << ",1,1,1,1,NONE\n";
                }
        }
        write(directory / "products.csv", products.str());
    }

    void makeBundle(int epoch)
    {
        const auto bundle = directory / "bundles" / ("epoch_" + std::to_string(epoch));
        fs::create_directory(bundle);
        const auto seconds = static_cast<long long>(times[epoch].bigTime);
        const std::string solution = epoch == 4 ? "FLOAT" : "PRODUCT_FIXED";
        std::ostringstream values, cov, contract;
        values << std::setprecision(17)
            << "gpst_seconds,solution,index,parameter,receiver,satellite,observable,unit,value,source_last_design_gpst,current_epoch_design_support,ppp_authority_valid,ppp_usable,pppar_usable,public_phase_segment,datum_version\n";
        for (int i = 0; i < 24; ++i)
        {
            const bool ion = i >= 9;
            const int sat = ion ? (i - 9) % 3 : i / 3;
            const int part = ion ? 0 : i % 3;
            const auto parameter = ion ? "IONO_L1_M" : (part == 0 ? "CLOCK" : "PHASE");
            const auto code = ion || part == 1 ? "L1C" : (part == 0 ? "NONE" : "L2W");
            values << seconds << ',' << solution << ',' << i << ',' << parameter << ',';
            if (ion) values << "R" << (i - 9) / 3;
            values << ',' << satellites[sat] << ',' << code << ",m," << mean(i) << ',';
            if (ion) values << seconds; else values << "nan";
            values << ',' << (ion ? 1 : 0) << ",1,1,1,UNIT_PUBLIC_SEGMENT,0\n";
        }
        cov << std::setprecision(17) << "gpst_seconds,solution,row_index,column_index,covariance_m2\n";
        for (int i = 0; i < 24; ++i)
            for (int j = i; j < 24; ++j)
            {
                if (epoch == 2 && i == 0 && j == 23) continue; // not a zero: absent cross block
                cov << seconds << ',' << solution << ',' << i << ',' << j << ',' << covariance(i, j) << '\n';
            }
        contract << "schema=REGIONAL_JOINT_V1\nepoch_gpst_seconds=" << seconds << "\nsolution=" << solution
            << "\nsource_numeric_root=INDEPENDENT_SYNTHETIC_ROOT_" << epoch
            << "\nparameters=24\nall_parameters_unit=m\nfull_same_epoch_cross_covariance=1\n"
            << "ionosphere=NETWORK_ESTIMABLE_GF_FIRST_SIGNAL_METRES_NOT_ABSOLUTE_TEC\n"
            << "spatial_prediction_error=NOT_PROVIDED\ntemporal_cross_covariance=NOT_PROVIDED\n";
        // Actual 02:00 bundles contain a repeated optional audit-file list.
        // These are never consumed as extra measurements or source covariance.
        for (const auto* name : {"regional_source_mean.bin", "regional_source_covariance.bin", "regional_source_keys.csv"})
        {
            const std::string contents = "SYNTHETIC_AUDIT_METADATA_NOT_CONSUMED_" + std::to_string(epoch) + "\n";
            write(bundle / name, contents);
            std::string failure;
            contract << "source_audit_file=" << name << " sha256="
                << zhangCheckpointFileSha256((bundle / name).string(), &failure)
                << " bytes=" << contents.size() << '\n';
        }
        if (epoch == 8) contract << "schema=REGIONAL_JOINT_V1\n"; // required singleton remains strict
        write(bundle / "regional_parameters.csv", values.str());
        write(bundle / "regional_covariance.csv", cov.str());
        write(bundle / "regional_contract.txt", contract.str());
        std::ostringstream manifest;
        manifest << "{\"epoch\":\"" << times[epoch].to_string(0) << "\",\"ar\":" << (epoch == 4 ? "false" : "true")
            << ",\"root\":\"SYNTHETIC_BUNDLE_ROOT_" << epoch << "\",\"regional_joint\":{\"schema\":\"REGIONAL_JOINT_V1\",\"files\":{";
        bool first = true;
        for (const auto* name : {"regional_contract.txt", "regional_parameters.csv", "regional_covariance.csv"})
        {
            if (!first) manifest << ',';
            first = false;
            std::string failure;
            auto hash = zhangCheckpointFileSha256((bundle / name).string(), &failure);
            require(hash.size() == 64, "fixture SHA256 computed by public file hash");
            if (epoch == 3 && std::string(name) == "regional_covariance.csv") hash.assign(64, '0');
            manifest << '"' << name << "\":{\"bytes\":" << fs::file_size(bundle / name)
                << ",\"sha256\":\"" << hash << "\"}";
        }
        manifest << "}}}\n";
        write(bundle / "COMMIT.json", manifest.str());
        write(bundle / "SHA256SUMS", zhangCheckpointSha256(manifest.str()) + "  COMMIT.json\n");
    }
};

KFKey ionosphereKey(int satellite)
{
    KFKey key; key.type = KF::IONO_STEC; key.str = user; key.Sat = SatSys(satellites[satellite].c_str());
    key.num = 0;
    return key;
}

void sharedFactors(KFMeasEntry& row, const std::string& identity, const std::vector<double>& factors)
{
    for (int i = 0; i < static_cast<int>(factors.size()); ++i)
    {
        KFKey key; key.type = KF::SAT_CLOCK; key.str = identity; key.num = i;
        row.addNoiseEntry(key, factors[i], 1);
    }
}

KFMeasEntryList physical(KFState& state, bool includeIonosphere, bool includePosition,
    std::string& identity)
{
    KFMeasEntryList entries;
    for (int sat = 0; sat < 3; ++sat)
        for (int part = 0; part < 3; ++part)
        {
            const auto code = part == 2 ? E_ObsCode::L2W : E_ObsCode::L1C;
            std::vector<double> factors;
            std::string current, failure;
            int rank = 0;
            require(zhangUserRegionalNoiseFactors(state.time, user, SatSys(satellites[sat].c_str()),
                code, part != 0, factors, current, &rank, &failure), "public code/phase regional factors: " + failure);
            require(rank == 24 && factors.size() == 24, "one complete 24-column source noise root");
            if (identity.empty()) identity = current;
            require(identity == current, "all code/phase rows share one current-epoch regional source");
            KFMeasEntry row(&state);
            row.obsKey.type = part == 0 ? KF::CODE_MEAS : KF::PHAS_MEAS;
            row.obsKey.str = user; row.obsKey.Sat = SatSys(satellites[sat].c_str()); row.obsKey.num = part;
            if (includeIonosphere)
            {
                InitialState init; init.estimate = true; init.x = 0; init.P = 100; init.Q = 0;
                row.addDsgnEntry(ionosphereKey(sat), part == 0 ? alpha1 : -alpha1, init);
            }
            if (includePosition)
                for (int axis = 0; axis < 3; ++axis)
                {
                    KFKey key; key.type = KF::REC_POS; key.str = user; key.num = axis;
                    row.usedValueMap[key] = axis == 0 ? 6378137 : 0;
                }
            sharedFactors(row, current, factors);
            entries.push_back(std::move(row));
        }
    return entries;
}

Eigen::VectorXd independentWeights(const Fixture& fixture)
{
    Eigen::MatrixXd x(5, 3);
    std::vector<std::pair<double, int>> distance;
    for (int i = 0; i < 5; ++i)
    { x.row(i) << 1, fixture.offsets[i].x(), fixture.offsets[i].y(); distance.emplace_back(fixture.offsets[i].squaredNorm(), i); }
    const Eigen::Matrix3d normalInverse = (x.transpose() * x).inverse();
    const Eigen::MatrixXd hat = x * normalInverse * x.transpose();
    const Eigen::VectorXd plane = x * normalInverse * Eigen::Vector3d(1, 0, 0);
    std::sort(distance.begin(), distance.end());
    Eigen::VectorXd idw = Eigen::VectorXd::Zero(5);
    for (int i = 0; i < 4; ++i) idw(distance[i].second) = 1 / distance[i].first;
    idw /= idw.sum();
    return plane + (Eigen::MatrixXd::Identity(5, 5) - hat).transpose() * idw;
}

std::string checkRuntime(const Fixture& fixture, int epoch)
{
    KFState state; state.time = fixture.times[epoch];
    std::string identity;
    auto entries = physical(state, true, true, identity);
    std::ostringstream trace;
    zhangUserRegionalIonoPseudoObs(trace, state, entries);
    require(entries.size() == 11, "actual runtime adds two independent SD ionosphere rows");
    require(trace.str().find("CURRENT_PHYSICAL_USED_REC_POS") != std::string::npos &&
        trace.str().find("user_truth_geometry_used=0") != std::string::npos,
        "physical operating position used instead of deliberately wrong USER geometry");
    const auto weights = independentWeights(fixture);
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(11, 24);
    for (int sat = 0; sat < 3; ++sat)
        for (int part = 0; part < 3; ++part)
        { d(3 * sat + part, 3 * sat) = 1; if (part) d(3 * sat + part, 3 * sat + part) = -1; }
    for (int sat = 1; sat < 3; ++sat)
        for (int station = 0; station < 5; ++station)
        { d(8 + sat, 9 + 3 * station + sat) = weights(station); d(8 + sat, 9 + 3 * station) = -weights(station); }
    std::map<KFKey, int> spatialColumns;
    int rowIndex = 0;
    for (const auto& entry : entries)
    {
        for (const auto& [key, value] : entry.noiseEntryMap)
        {
            require(entry.noiseElementMap.at(key) == 1, "factor variance one");
            if (key.str.find(":SPACE:") != std::string::npos)
                if (!spatialColumns.count(key)) spatialColumns[key] = spatialColumns.size();
        }
        if (rowIndex >= 9)
        {
            require(entry.obsKey.type == KF::IONOSPHERIC && entry.obsKey.num == 1,
                "independent SD rows use common G01 reference");
            double commonGf = 0;
            for (const auto& [key, value] : entry.designEntryMap)
                if (key.type == KF::IONO_STEC) commonGf += value;
            require(commonGf == 0, "user common GF direction has exactly zero pseudo-design coefficient");
            require(std::abs(entry.innov - (d.row(rowIndex) * fixture.mean)(0)) < 1e-10,
                "pseudo mean equals independent trend-IDW SD prediction");
            bool sharesRoot = false;
            for (const auto& [key, value] : entry.noiseEntryMap) if (key.str == identity) sharesRoot = true;
            require(sharesRoot, "pseudo rows share SAME source factors as code/phase");
        }
        ++rowIndex;
    }
    Eigen::MatrixXd space = Eigen::MatrixXd::Zero(11, spatialColumns.size());
    rowIndex = 0;
    for (const auto& entry : entries)
    {
        for (const auto& [key, value] : entry.noiseEntryMap)
            if (spatialColumns.count(key)) space(rowIndex, spatialColumns.at(key)) = value;
        ++rowIndex;
    }
    const Eigen::MatrixXd spatialGram = space * space.transpose();
    require(spatialGram(9, 10) > 0, "shared-reference spatial off-diagonal retained");
    const Eigen::MatrixXd formal = d * fixture.covariance * d.transpose();
    require(formal.topRightCorner(9, 2).cwiseAbs().maxCoeff() > 1e-5,
        "independent source has genuine CODE/PHASE to IONO cross covariance");
    state.stateTransition(trace, state.time);
    KFMeas actual(state, entries, state.time);
    near(actual.R, formal + spatialGram, 1e-10, "ACTUAL KFMeas R equals independent D Q D^T plus separate spatial Gram");
    near(actual.R.topRightCorner(9, 2), formal.topRightCorner(9, 2), 1e-10,
        "ACTUAL KFMeas cross block not zeroed or double counted");
    require(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(actual.R).eigenvalues().minCoeff() > 0,
        "actual assembled joint R positive definite");
    return identity;
}

void expectRejected(const Fixture& fixture, int epoch, const std::string& expected)
{
    std::vector<double> factors; std::string id, reason;
    bool rejected = false;
    try { zhangUserRegionalNoiseFactors(fixture.times[epoch], user, SatSys("G01"), E_ObsCode::L1C,
        false, factors, id, nullptr, &reason); }
    catch (const std::exception& error) { rejected = std::string(error.what()).find(expected) != std::string::npos; }
    require(rejected, "actual runtime must reject " + expected);
}
} // namespace

extern "C" int __wrap_main()
{
    try
    {
        Fixture fixture;
        environment("ZHANG_USER_REGIONAL_MODE", "SD_PLANE_JOINT");
        environment("ZHANG_USER_REGIONAL_BUNDLE_ROOT", (fixture.directory / "bundles").string());
        environment("ZHANG_USER_REGIONAL_GEOMETRY", (fixture.directory / "geometry.csv").string());
        environment("ZHANG_USER_PRODUCT_ERROR_POLICY", "WHITE");
        acsConfig.zhangPppAr.user_adapter = true;
        acsConfig.zhangPppAr.user_use_full_product_covariance = true;
        acsConfig.zhangPppAr.output_products = false;
        acsConfig.zhangPppAr.product_solution = "PRODUCT_FIXED";
        acsConfig.zhangPppAr.product_filename = (fixture.directory / "products.csv").string();
        acsConfig.zhangPppAr.product_covariance_filename = (fixture.directory / "DELIBERATELY_ABSENT_LEGACY_COV.csv").string();
        acsConfig.zhangFullRank.enable = false;
        acsConfig.pppOpts.ionoOpts.common_ionosphere = true;
        acsConfig.pppOpts.equate_ionospheres = false;
        acsConfig.pppOpts.ionoOpts.use_if_combo = acsConfig.pppOpts.ionoOpts.use_gf_combo = false;
        auto& receiver = acsConfig.getRecOpts(user);
        receiver.ion_stec.estimate = {true}; receiver.ion_stec.process_noise = {0}; receiver.ion_stec.sigma = {100};
        require(zhangUserRegionalEnabled(), "public runtime explicitly enabled");
        const auto first = checkRuntime(fixture, 0);
        const auto second = checkRuntime(fixture, 1);
        require(first != second, "successive white epochs have distinct shared source identities");
        // Missing physical support never creates an ionosphere or position state.
        for (const auto flags : {std::make_pair(false, true), std::make_pair(true, false)})
        {
            KFState state; state.time = fixture.times[1];
            std::string identity;
            auto entries = physical(state, flags.first, flags.second, identity);
            std::ostringstream trace;
            zhangUserRegionalIonoPseudoObs(trace, state, entries);
            require(entries.size() == 9, "missing physical IONO/REC_POS locally skips pseudo rows");
        }
        expectRejected(fixture, 2, "COVARIANCE_CROSS_BLOCKS_UNSPECIFIED");
        expectRejected(fixture, 3, "FILE_SHA256_MISMATCH");
        std::vector<double> factors; std::string identity, reason;
        require(!zhangUserRegionalNoiseFactors(fixture.times[4], user, SatSys("G01"), E_ObsCode::L1C,
            false, factors, identity, nullptr, &reason) && reason == "REQUESTED_SOLUTION_NOT_PUBLISHED",
            "FLOAT-only product gap is a safe no-update, not a fabricated fixed source");
        const auto recovered = checkRuntime(fixture, 5);
        require(recovered != second, "source cache recovers after malformed/gap epochs");
        require(!zhangUserRegionalNoiseFactors(fixture.times[6], user, SatSys("G01"), E_ObsCode::L1C,
            false, factors, identity, nullptr, &reason) && reason == "COMMITTED_EPOCH_UNAVAILABLE",
            "missing committed bundle safely unavailable");
        checkRuntime(fixture, 7);
        expectRejected(fixture, 8, "CONTRACT_DUPLICATE_FIELD");
        require(!fs::exists(fixture.directory / "DELIBERATELY_ABSENT_LEGACY_COV.csv"),
            "regional shared covariance never consults or creates legacy covariance root");
        std::cout << "USER_REGIONAL_JOINT_CONTRACT PASS production_runtime=1 actual_KFMeas_R=1 independent_full_Q_oracle=1 shared_CP_IONO_cross=1 physical_geometry_only=1 malformed_hash_Q_fail_closed=1 optional_audit_list_supported=1 duplicate_singleton_rejected=1 pilot_boundary_NOT_covered=pppSatClocks fixture="
            << fixture.directory << '\n';
        return 0;
    }
    catch (const std::exception& error)
    { std::cerr << "USER_REGIONAL_JOINT_CONTRACT FAIL " << error.what() << '\n'; return 1; }
}
