#include "common/zhangRealtimeFeed.hpp"
#include "pea/zhangUserRegional.hpp"
#include "common/zhangRegionalTropBackground.hpp"
#include "common/zhangUserMultiGnss.hpp"

#include <algorithm>
#include <boost/json.hpp>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "common/acsConfig.hpp"
#include "common/constants.hpp"
#include "common/enumHelpers.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserProductNoiseContract.hpp"
#include "common/zhangUserRegionalSpatial.hpp"
#include "pea/zhangPppAr.hpp"

namespace
{
namespace fs = std::filesystem;
using Eigen::MatrixXd;
using Eigen::VectorXd;

[[noreturn]] void reject(const std::string& reason)
{
    throw std::runtime_error("USER_REGIONAL:" + reason);
}

std::string trim(std::string text)
{
    const auto begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const auto end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

std::string readText(const fs::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) reject("FILE_UNAVAILABLE:" + path.string());
    std::ostringstream content;
    content << stream.rdbuf();
    if (stream.bad()) reject("FILE_READ_FAILED:" + path.string());
    return content.str();
}

std::vector<std::string> csv(const std::string& line)
{
    std::vector<std::string> result;
    std::string value;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
        const char ch = line[i];
        if (ch == '"')
        {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"')
            { value += '"'; ++i; }
            else quoted = !quoted;
        }
        else if (ch == ',' && !quoted)
        { result.push_back(trim(value)); value.clear(); }
        else value += ch;
    }
    if (quoted) reject("CSV_UNTERMINATED_QUOTE");
    result.push_back(trim(value));
    return result;
}

struct CsvTable
{
    std::map<std::string, std::size_t> columns;
    std::vector<std::vector<std::string>> rows;
    explicit CsvTable(const fs::path& path)
    {
        std::ifstream input(path);
        std::string line;
        if (!input || !std::getline(input, line)) reject("CSV_HEADER_MISSING:" + path.string());
        const auto header = csv(line);
        for (std::size_t i = 0; i < header.size(); ++i)
            if (!columns.emplace(header[i], i).second) reject("CSV_DUPLICATE_COLUMN:" + header[i]);
        while (std::getline(input, line))
        {
            if (trim(line).empty()) continue;
            auto fields = csv(line);
            if (fields.size() != header.size()) reject("CSV_ROW_WIDTH:" + path.string());
            rows.push_back(std::move(fields));
        }
        if (input.bad()) reject("CSV_READ_FAILED:" + path.string());
    }
    const std::string& at(const std::vector<std::string>& row, const std::string& key) const
    {
        const auto found = columns.find(key);
        if (found == columns.end()) reject("CSV_COLUMN_MISSING:" + key);
        return row.at(found->second);
    }
};

double number(const std::string& text, bool allowNan = false)
{
    std::size_t end = 0;
    double result = 0;
    try { result = std::stod(text, &end); }
    catch (...) { reject("INVALID_NUMBER:" + text); }
    if (end != text.size() || (!std::isfinite(result) && !(allowNan && std::isnan(result))))
        reject("NONFINITE_OR_TRAILING_NUMBER:" + text);
    return result;
}

int integer(const std::string& text)
{
    const double value = number(text);
    if (value < 0 || value > 100000000 || value != std::floor(value)) reject("INVALID_INTEGER:" + text);
    return static_cast<int>(value);
}

bool flag(const std::string& text)
{
    if (text != "0" && text != "1") reject("INVALID_BOOLEAN:" + text);
    return text == "1";
}

long long epochSeconds(const GTime& time)
{
    const auto value = static_cast<long double>(time.bigTime);
    if (!std::isfinite(value) || std::abs(value - std::round(value)) > 1e-6L)
        reject("NON_INTEGER_PRODUCT_EPOCH");
    return std::llround(value);
}

long long epochText(const std::string& text)
{
    if (text.size() != 19) reject("INVALID_COMMIT_EPOCH:" + text);
    // str2time() accepts six SPACE-separated fields, not ISO punctuation.
    // Keep the manifest spelling strict and verify a round trip so that
    // normalization of an invalid calendar date cannot authorize a bundle.
    std::string fields = text;
    for (std::size_t i = 0; i < fields.size(); ++i)
    {
        const char separator = (i == 4 || i == 7) ? '-' :
            (i == 10 ? ' ' : ((i == 13 || i == 16) ? ':' : '\0'));
        if (separator)
        {
            if (fields[i] != separator) reject("INVALID_COMMIT_EPOCH:" + text);
            fields[i] = ' ';
        }
        else if (fields[i] < '0' || fields[i] > '9') reject("INVALID_COMMIT_EPOCH:" + text);
    }
    GTime time;
    if (str2time(fields.c_str(), 0, fields.size(), time, E_TimeSys::GPST) != 0 || time.to_string(0) != text)
        reject("INVALID_COMMIT_EPOCH:" + text);
    return epochSeconds(time);
}

const boost::json::object object(const std::string& text)
{
    boost::system::error_code error;
    const auto parsed = boost::json::parse(text, error);
    if (error || !parsed.is_object()) reject("INVALID_COMMIT_JSON");
    return parsed.as_object();
}

std::string jsonString(const boost::json::object& obj, const char* key)
{
    const auto* value = obj.if_contains(key);
    if (!value || !value->is_string() || value->as_string().empty()) reject(std::string("COMMIT_FIELD:") + key);
    return std::string(value->as_string());
}

std::map<std::string, std::string> contract(const fs::path& path)
{
    std::map<std::string, std::string> result;
    std::istringstream stream(readText(path));
    std::string line;
    while (std::getline(stream, line))
    {
        const auto equal = line.find('=');
        if (equal == std::string::npos) continue;
        auto key = trim(line.substr(0, equal));
        const auto value = trim(line.substr(equal + 1));
        // REGIONAL_JOINT_V1 optionally lists several source-audit files under
        // this repeated key. They are metadata, not additional observations or
        // covariance. Required singleton fields must STILL reject duplicates.
        if (key == "source_audit_file")
        {
            std::istringstream descriptor(value);
            std::string filename;
            descriptor >> filename;
            if (filename != "regional_source_mean.bin" &&
                filename != "regional_source_covariance.bin" &&
                filename != "regional_source_keys.csv") reject("CONTRACT_AUDIT_FILE_IDENTITY");
            key += ":" + filename;
        }
        if (!result.emplace(key, value).second)
            reject("CONTRACT_DUPLICATE_FIELD");
    }
    return result;
}

struct Geometry
{
    Eigen::Vector3d xyz;
    std::string role;
};

struct EpochLocation
{
    fs::path directory;
    std::string commitHash;
};

bool regionalSwitch(const char* name, bool defaultValue)
{
    const char* value=std::getenv(name);
    const std::string mode=value?value:(defaultValue?"ON":"OFF");
    if (mode!="ON" && mode!="OFF") reject(std::string("UNKNOWN_SWITCH:")+name);
    return mode=="ON";
}

struct Settings
{
    bool enabled = false;
    fs::path root;
    std::map<std::string, Geometry> geometry;
    std::map<long long, EpochLocation> epochs;
    Settings()
    {
        const char* value = std::getenv("ZHANG_USER_REGIONAL_MODE");
        const std::string mode = value ? value : "OFF";
        if (mode == "OFF") return;
        if (mode != "SD_PLANE_JOINT") reject("UNKNOWN_MODE:" + mode);
        enabled = true;
        const char* folder = std::getenv("ZHANG_USER_REGIONAL_BUNDLE_ROOT");
        const char* coordinates = std::getenv("ZHANG_USER_REGIONAL_GEOMETRY");
        if (!folder || !*folder || !coordinates || !*coordinates) reject("ENABLED_PATHS_REQUIRED");
        root = fs::canonical(folder);
        if (!fs::is_directory(root)) reject("BUNDLE_ROOT_NOT_DIRECTORY");
        CsvTable table(coordinates);
        for (const auto& row : table.rows)
        {
            Geometry site;
            site.role = table.at(row, "role");
            if (site.role == "REFERENCE")
            {
                site.xyz << number(table.at(row, "x_m")), number(table.at(row, "y_m")), number(table.at(row, "z_m"));
                if (site.xyz.norm() < 6e6 || site.xyz.norm() > 7e6) reject("REFERENCE_GEOMETRY_ECEF_INVALID");
            }
            else if (site.role == "USER")
                site.xyz.setZero(); // USER geometry values are never read, even for fallback.
            else reject("GEOMETRY_ROLE_INVALID");
            const auto& name = table.at(row, "receiver");
            if (name.empty() || !geometry.emplace(name, std::move(site)).second) reject("GEOMETRY_DUPLICATE_RECEIVER");
        }
        if (!zhangRealtime::enabled())
        for (const auto& entry : fs::directory_iterator(root))
        {
            const auto name = entry.path().filename().string();
            if (!entry.is_directory() || name.ends_with(".pending")) continue;
            const auto manifest = entry.path() / "COMMIT.json";
            if (!fs::is_regular_file(manifest)) continue;
            const auto bytes = readText(manifest);
            const auto doc = object(bytes);
            const long long epoch = epochText(jsonString(doc, "epoch"));
            if (!epochs.emplace(epoch, EpochLocation{fs::canonical(entry.path()), zhangCheckpointSha256(bytes)}).second)
                reject("DUPLICATE_COMMITTED_EPOCH:" + std::to_string(epoch));
        }
        if (zhangRealtime::enabled())
        {
            if (root != zhangRealtime::root()) reject("REALTIME_REGIONAL_ROOT_MISMATCH");
        }
        else if (epochs.empty()) reject("NO_COMMITTED_EPOCHS");
    }
};

Settings& settings()
{
    static Settings value;
    return value;
}

void validateConfiguration()
{
    if (!acsConfig.zhangPppAr.user_adapter || !acsConfig.zhangPppAr.user_use_full_product_covariance ||
        acsConfig.zhangFullRank.enable || acsConfig.zhangPppAr.output_products ||
        !acsConfig.pppOpts.ionoOpts.common_ionosphere || acsConfig.pppOpts.equate_ionospheres ||
        acsConfig.pppOpts.ionoOpts.use_if_combo || acsConfig.pppOpts.ionoOpts.use_gf_combo)
        reject("INCOMPATIBLE_USER_MODEL_CONFIGURATION");
    const char* temporal = std::getenv("ZHANG_USER_PRODUCT_ERROR_POLICY");
    if (!temporal || std::string(temporal) != "WHITE")
        reject("REGIONAL_CURRENT_EPOCH_REQUIRES_EXPLICIT_WHITE_POLICY");
}

struct Parameter
{
    std::string parameter, receiver, satellite, observable;
    int fileIndex = -1;
    double mean = 0;
    bool eligible = false;
};

struct EpochCache
{
    long long epoch = std::numeric_limits<long long>::min();
    std::string solution, noiseId, reason = "NOT_LOADED", numericRoot;
    std::vector<Parameter> parameters;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, int> index;
    VectorXd mean;
    MatrixXd covariance, root;
    int sourceDimension = 0;
    std::uint64_t loadSequence = 0;
    double minimumScaledEigenvalue = 0;
    bool valid = false;
    std::map<std::string,std::string> properties;
};

EpochCache cache;
std::uint64_t regionalLoadSequence = 0;
std::recursive_mutex regionalCacheMutex;

int lookup(const EpochCache& value, const std::string& parameter, const std::string& receiver,
           const std::string& satellite, const std::string& observable)
{
    const auto found = value.index.find({parameter, receiver, satellite, observable});
    return found == value.index.end() ? -1 : found->second;
}

MatrixXd scaledRoot(const MatrixXd& q, double& minimumEigenvalue)
{
    if (q.rows() == 0 || q.rows() != q.cols() || !q.allFinite()) reject("COVARIANCE_DIMENSION_OR_FINITE");
    VectorXd scale(q.rows());
    for (int i = 0; i < q.rows(); ++i)
    {
        if (q(i, i) < 0) reject("COVARIANCE_NEGATIVE_DIAGONAL");
        scale(i) = std::sqrt(std::max(q(i, i), 1e-24));
    }
    const MatrixXd normalized = scale.cwiseInverse().asDiagonal() * q * scale.cwiseInverse().asDiagonal();
    Eigen::SelfAdjointEigenSolver<MatrixXd> eig(0.5 * (normalized + normalized.transpose()));
    if (eig.info() != Eigen::Success || !eig.eigenvalues().allFinite()) reject("SCALED_COVARIANCE_EIGEN_FAILURE");
    minimumEigenvalue = eig.eigenvalues().minCoeff();
    const double tolerance = std::max(1e-12, eig.eigenvalues().cwiseAbs().maxCoeff() * 1e-9);
    if (minimumEigenvalue < -tolerance) reject("COVARIANCE_NOT_PSD");
    std::vector<int> retained;
    for (int i = 0; i < eig.eigenvalues().size(); ++i)
        if (eig.eigenvalues()(i) > 0) retained.push_back(i);
    MatrixXd root(q.rows(), retained.size());
    for (int i = 0; i < static_cast<int>(retained.size()); ++i)
        root.col(i) = scale.asDiagonal() * eig.eigenvectors().col(retained[i]) * std::sqrt(eig.eigenvalues()(retained[i]));
    const MatrixXd rootError = scale.cwiseInverse().asDiagonal() *
        (root * root.transpose() - q) * scale.cwiseInverse().asDiagonal();
    if (!root.allFinite() || !rootError.allFinite() || rootError.norm() > 2e-8 * std::max(1.0, normalized.norm()))
        reject("COVARIANCE_ROOT_RECONSTRUCTION");
    return root;
}

void verifyHashes(const EpochLocation& location, const boost::json::object& doc)
{
    const auto* regional = doc.if_contains("regional_joint");
    if (!regional || !regional->is_object()) reject("COMMIT_REGIONAL_MANIFEST_MISSING");
    const auto& packet = regional->as_object();
    if (jsonString(packet, "schema") != "REGIONAL_JOINT_V1" && jsonString(packet, "schema") != "REGIONAL_JOINT_V2") reject("COMMIT_REGIONAL_SCHEMA");
    const auto* files = packet.if_contains("files");
    if (!files || !files->is_object()) reject("COMMIT_REGIONAL_FILES_MISSING");
    for (const char* name : {"regional_contract.txt", "regional_parameters.csv", "regional_covariance.csv"})
    {
        const auto* descriptor = files->as_object().if_contains(name);
        if (!descriptor || !descriptor->is_object()) reject(std::string("COMMIT_REGIONAL_FILE:") + name);
        const auto& entry = descriptor->as_object();
        const auto* bytes = entry.if_contains("bytes");
        if (!bytes || (!bytes->is_int64() && !bytes->is_uint64())) reject("COMMIT_FILE_BYTES");
        const auto size = bytes->is_uint64() ? bytes->as_uint64() : static_cast<std::uint64_t>(bytes->as_int64());
        const auto path = location.directory / name;
        if (!fs::is_regular_file(path) || fs::file_size(path) != size) reject(std::string("FILE_SIZE_MISMATCH:") + name);
        std::string failure;
        if (zhangCheckpointFileSha256(path.string(), &failure) != jsonString(entry, "sha256"))
            reject(std::string("FILE_SHA256_MISMATCH:") + name);
    }
    const auto sums = readText(location.directory / "SHA256SUMS");
    std::istringstream stream(sums);
    std::string hash, filename;
    bool commitFound = false;
    while (stream >> hash >> filename)
        if (filename == "COMMIT.json")
        {
            if (commitFound || hash != location.commitHash) reject("COMMIT_HASH_MISMATCH");
            commitFound = true;
        }
    if (!commitFound) reject("COMMIT_HASH_MISSING");
}

bool loadEpoch(const GTime& time)
{
    auto& config = settings();
    const long long epoch = epochSeconds(time);
    const auto& solution = acsConfig.zhangPppAr.product_solution;
    if (cache.epoch == epoch && cache.solution == solution) return cache.valid;
    cache = {};
    cache.epoch = epoch;
    cache.solution = solution;
    cache.loadSequence = ++regionalLoadSequence;
    if (zhangRealtime::enabled())
    {
        // The main-epoch barrier owns discovery and its wait deadline. Retain
        // exactly one location, and never cache a live miss across epochs.
        const auto pin = zhangRealtime::current(time);
        config.epochs.clear();
        if (pin.bundle.valid)
            config.epochs.emplace(epoch, EpochLocation{pin.bundle.directory, pin.commitHash});
    }
    const auto location = config.epochs.find(epoch);
    if (location == config.epochs.end()) { cache.reason = "COMMITTED_EPOCH_UNAVAILABLE"; return false; }
    const auto bytes = readText(location->second.directory / "COMMIT.json");
    if (zhangCheckpointSha256(bytes) != location->second.commitHash) reject("COMMIT_CHANGED_AFTER_INDEX");
    const auto doc = object(bytes);
    if (epochText(jsonString(doc, "epoch")) != epoch) reject("COMMIT_EPOCH_MISMATCH");
    verifyHashes(location->second, doc);
    const auto properties = contract(location->second.directory / "regional_contract.txt");
    const auto require = [&](const std::string& key, const std::string& expected)
    {
        const auto found = properties.find(key);
        if (found == properties.end() || found->second != expected) reject("CONTRACT:" + key);
    };
    const auto manifestSchema=jsonString(doc.at("regional_joint").as_object(),"schema");
    require("schema",manifestSchema);
    cache.properties=properties;
    require("all_parameters_unit", "m");
    require("full_same_epoch_cross_covariance", "1");
    require("ionosphere", "NETWORK_ESTIMABLE_GF_FIRST_SIGNAL_METRES_NOT_ABSOLUTE_TEC");
    const auto sourceEpoch = properties.find("epoch_gpst_seconds");
    if (sourceEpoch == properties.end() || number(sourceEpoch->second) != static_cast<double>(epoch)) reject("CONTRACT_EPOCH");
    // Ordinary FLOAT-only publication gaps are not malformed AR epochs.
    const auto sourceSolution = properties.find("solution");
    if (sourceSolution == properties.end()) reject("CONTRACT_SOLUTION_MISSING");
    if (sourceSolution->second != solution) { cache.reason = "REQUESTED_SOLUTION_NOT_PUBLISHED"; return false; }
    if (solution == "PRODUCT_FIXED")
    {
        const auto* ar = doc.if_contains("ar");
        if (!ar || !ar->is_bool() || !ar->as_bool()) reject("FIXED_SOLUTION_WITHOUT_COMMIT_AR_AUTHORITY");
    }
    cache.numericRoot = properties.at("source_numeric_root");
    CsvTable values(location->second.directory / "regional_parameters.csv");
    std::vector<Parameter> parameters(values.rows.size());
    if (parameters.empty() || parameters.size() > 4096) reject("PARAMETER_DIMENSION");
    const auto declared = properties.find("parameters");
    if (declared == properties.end() || integer(declared->second) != static_cast<int>(parameters.size())) reject("PARAMETER_COUNT");
    std::set<int> seen;
    std::set<std::tuple<std::string, std::string, std::string, std::string>> identities;
    for (const auto& row : values.rows)
    {
        if (number(values.at(row, "gpst_seconds")) != static_cast<double>(epoch) || values.at(row, "solution") != solution || values.at(row, "unit") != "m")
            reject("PARAMETER_EPOCH_SOLUTION_UNIT");
        Parameter p;
        p.fileIndex = integer(values.at(row, "index"));
        if (p.fileIndex >= static_cast<int>(parameters.size()) || !seen.insert(p.fileIndex).second) reject("PARAMETER_INDEX");
        p.parameter = values.at(row, "parameter"); p.receiver = values.at(row, "receiver");
        p.satellite = values.at(row, "satellite"); p.observable = values.at(row, "observable");
        p.mean = number(values.at(row, "value"));
        const SatSys sat(p.satellite.c_str());
        if (p.parameter!="TROP_ZTD_M" && p.parameter!="TROP_COEFF" && (!zhangUserGecSystem(sat.sys) || sat.prn <= 0 || sat.id() != p.satellite)) reject("UNSUPPORTED_SATELLITE_ID");
        if (!identities.emplace(p.parameter, p.receiver, p.satellite, p.observable).second) reject("DUPLICATE_PARAMETER_IDENTITY");
        const bool authority = flag(values.at(row, "ppp_authority_valid"));
        const bool usable = flag(values.at(row, "ppp_usable"));
        const bool fresh = flag(values.at(row, "current_epoch_design_support"));
        const double lastDesign = number(values.at(row, "source_last_design_gpst"), true);
        if (p.parameter == "IONO_L1_M")
        {
            const auto station = config.geometry.find(p.receiver);
            if (p.receiver.empty() || p.observable != enum_to_string(zhangUserGecCodes(sat.sys).first) || station == config.geometry.end() || station->second.role != "REFERENCE")
                reject("IONOSPHERE_PARAMETER_IDENTITY");
            if (fresh && (!std::isfinite(lastDesign) || std::abs(lastDesign - epoch) > 1e-6)) reject("IONOSPHERE_FRESHNESS_CONFLICT");
            p.eligible = authority && usable && fresh && std::isfinite(lastDesign) && std::abs(lastDesign - epoch) <= 1e-6;
        }
        else if (p.parameter == "CLOCK" || p.parameter == "PHASE")
        {
            if (!p.receiver.empty() || (p.parameter == "CLOCK" && p.observable != "NONE") ||
                (p.parameter == "PHASE" && p.observable != enum_to_string(zhangUserGecCodes(sat.sys).first) &&
                    p.observable != enum_to_string(zhangUserGecCodes(sat.sys).second))) reject("CLOCK_PHASE_PARAMETER_IDENTITY");
            p.eligible = authority && usable;
        }
        else if (p.parameter=="TROP_ZTD_M" || p.parameter=="TROP_COEFF")
        {
            if (manifestSchema!="REGIONAL_JOINT_V2" || !p.satellite.empty()) reject("TROP_SCHEMA_OR_SATELLITE");
            if (p.parameter=="TROP_ZTD_M") {
                const auto site=config.geometry.find(p.receiver);
                if (site==config.geometry.end() || site->second.role!="REFERENCE" || p.observable!="NONE") reject("TROP_REFERENCE_IDENTITY");
                p.eligible=false; // audit source row; consumer uses the projection
            }
            else {
                if (p.receiver!="REGIONAL" || (p.observable!="0" && p.observable!="1" && p.observable!="2" && p.observable!="3")) reject("TROP_COEFFICIENT_IDENTITY");
                p.eligible=authority && usable && fresh && std::isfinite(lastDesign) && lastDesign==epoch;
            }
        }
        else reject("UNKNOWN_REGIONAL_PARAMETER:" + p.parameter);
        parameters[p.fileIndex] = std::move(p);
    }
    const int n = parameters.size();
    cache.sourceDimension = n;
    MatrixXd q = MatrixXd::Zero(n, n);
    CsvTable covariance(location->second.directory / "regional_covariance.csv");
    std::set<std::pair<int, int>> entries;
    for (const auto& row : covariance.rows)
    {
        if (number(covariance.at(row, "gpst_seconds")) != static_cast<double>(epoch) || covariance.at(row, "solution") != solution)
            reject("COVARIANCE_EPOCH_SOLUTION");
        int a = integer(covariance.at(row, "row_index")), b = integer(covariance.at(row, "column_index"));
        if (a >= n || b >= n) reject("COVARIANCE_INDEX");
        if (a > b) std::swap(a, b);
        if (!entries.emplace(a, b).second) reject("DUPLICATE_COVARIANCE_ENTRY");
        const double value = number(covariance.at(row, "covariance_m2"));
        q(a, b) = q(b, a) = value;
    }
    if (entries.size() != std::size_t(n) * (n + 1) / 2) reject("COVARIANCE_CROSS_BLOCKS_UNSPECIFIED");
    // Validate and factor the COMPLETE supplied source Q once, including rows
    // which are not currently eligible. Inactive rows must not hide a corrupt
    // covariance. Marginal selection below only selects rows of this root.
    const MatrixXd fullRoot = scaledRoot(q, cache.minimumScaledEigenvalue);
    // Marginal selection, NEVER a conditional Schur complement or independent
    // block approximation. Every retained cross block remains supplied by Q.
    std::vector<int> retained;
    for (const auto& p : parameters) if (p.eligible) retained.push_back(p.fileIndex);
    if (retained.empty()) { cache.reason = "NO_CURRENT_PPP_COORDINATES"; return false; }
    cache.mean.resize(retained.size()); cache.covariance.resize(retained.size(), retained.size());
    cache.root.resize(retained.size(), fullRoot.cols());
    for (int i = 0; i < static_cast<int>(retained.size()); ++i)
    {
        const auto& p = parameters[retained[i]];
        cache.parameters.push_back(p); cache.mean(i) = p.mean;
        cache.root.row(i) = fullRoot.row(retained[i]);
        cache.index[{p.parameter, p.receiver, p.satellite, p.observable}] = i;
        for (int j = 0; j < static_cast<int>(retained.size()); ++j) cache.covariance(i, j) = q(retained[i], retained[j]);
    }
    // Verify the same nominal publication is used by the existing user chart
    // and by this covariance source; no product data are rewritten here.
    // A CLOCK can be explicitly unavailable while one PHASE row remains
    // available (for example, first-signal suspension). Verify nominal means
    // against the complete publication, not the eligible marginal. The noise
    // accessor still rejects observations without an eligible clock; this
    // lookup must never promote a source coordinate's usability flag.
    std::map<std::string, double> publishedClockMeans;
    for (const auto& parameter : parameters)
        if (parameter.parameter == "CLOCK")
            publishedClockMeans.emplace(parameter.satellite, parameter.mean);
    for (const auto& p : cache.parameters)
        if (p.parameter == "PHASE")
        {
            ZhangInternalProduct product;
            const auto code = string_to_enum<E_ObsCode>(p.observable);
            if (!queryZhangInternalProduct(time, SatSys(p.satellite.c_str()), code, product) ||
                product.solution != solution || !product.ppp_usable || std::abs(product.phase_m - p.mean) > 1e-8)
                reject("PHASE_MEAN_PUBLICATION_MISMATCH:" + p.satellite + ":" + p.observable);
            const auto clock = publishedClockMeans.find(p.satellite);
            if (clock == publishedClockMeans.end())
                reject("CLOCK_SOURCE_PUBLICATION_MISSING:" + p.satellite);
            if (std::abs(product.clock_m - clock->second) > 1e-8)
                reject("CLOCK_MEAN_PUBLICATION_MISMATCH:" + p.satellite);
        }
    cache.noiseId = "REGIONAL_JOINT:" + std::to_string(epoch) + ":" + solution + ":" +
        jsonString(doc, "root") + ":" + location->second.commitHash;
    cache.valid = true; cache.reason = "NONE";
    return true;
}

Eigen::Vector2d horizontalKm(const Eigen::Vector3d& user, const Eigen::Vector3d& station)
{
    const double longitude = std::atan2(user(1), user(0));
    const double radius = std::hypot(user(0), user(1));
    constexpr double a = 6378137, eccentricitySquared = 0.0066943799901413165;
    double latitude = std::atan2(user(2), radius * (1 - eccentricitySquared));
    for (int i = 0; i < 12; ++i)
    {
        const double normal = a / std::sqrt(1 - eccentricitySquared * std::pow(std::sin(latitude), 2));
        latitude = std::atan2(user(2) + eccentricitySquared * normal * std::sin(latitude), radius);
    }
    const auto delta = station - user;
    return Eigen::Vector2d((-std::sin(longitude) * delta(0) + std::cos(longitude) * delta(1)) / 1000,
        (-std::sin(latitude) * std::cos(longitude) * delta(0) - std::sin(latitude) * std::sin(longitude) * delta(1) + std::cos(latitude) * delta(2)) / 1000);
}

void addFactors(KFMeasEntry& entry, const std::string& identity, const Eigen::RowVectorXd& factors)
{
    for (int i = 0; i < factors.size(); ++i)
    {
        KFKey key;
        key.type = KF::SAT_CLOCK; key.Sat = SatSys{}; key.str = identity; key.num = i;
        entry.addNoiseEntry(key, factors(i), 1);
    }
}

std::vector<zhangregional::PairSupport> pairSupports(
    const std::string& reference, const std::map<std::string, KFKey>& keys, const Eigen::Vector3d& userOperatingPosition)
{
    std::vector<zhangregional::PairSupport> result;
    for (const auto& [satellite, key] : keys)
    {
        if (satellite == reference) continue;
        if (key.Sat.sys != SatSys(reference.c_str()).sys) reject("CROSS_SYSTEM_IONOSPHERE_SD_FORBIDDEN");
        const auto firstCode = enum_to_string(zhangUserGecCodes(key.Sat.sys).first);
        zhangregional::PairSupport pair;
        pair.satellite = satellite; pair.reference = reference;
        for (const auto& [station, position] : settings().geometry)
        {
            if (position.role != "REFERENCE") continue;
            const int left = lookup(cache, "IONO_L1_M", station, satellite, firstCode);
            const int right = lookup(cache, "IONO_L1_M", station, reference, firstCode);
            if (left < 0 || right < 0) continue;
            const auto offset = horizontalKm(userOperatingPosition, position.xyz);
            pair.stations.push_back({station, offset(0), offset(1), left, right});
        }
        result.push_back(std::move(pair));
    }
    return result;
}

bool physicalOperatingPosition(const std::string& receiver, const KFMeasEntryList& entries,
                               const KFState& state, Eigen::Vector3d& position, std::string& origin)
{
    // The USER coordinate in the geometry file is deliberately NOT consumed.
    // Use the actual REC_POS operating point already used to construct a
    // current physical measurement, not evaluation truth or an AR solution.
    bool existing = true;
    for (int component = 0; component < 3; ++component)
    {
        KFKey key; key.type = KF::REC_POS; key.str = receiver; key.num = component;
        if (state.getKFValue(key, position(component), nullptr, nullptr, false) == E_Source::NONE ||
            !std::isfinite(position(component))) { existing = false; break; }
    }
    if (existing && position.norm() > 6e6 && position.norm() < 7e6)
    { origin = "CURRENT_FILTER_REC_POS"; return true; }
    for (const auto& row : entries)
    {
        if (!row.valid || row.obsKey.str != receiver ||
            (row.obsKey.type != KF::CODE_MEAS && row.obsKey.type != KF::PHAS_MEAS)) continue;
        std::set<int> components;
        for (const auto& [key, value] : row.usedValueMap)
            if (key.type == KF::REC_POS && key.str == receiver && key.num >= 0 && key.num < 3 && std::isfinite(value))
            { position(key.num) = value; components.insert(key.num); }
        if (components.size() == 3 && position.allFinite() && position.norm() > 6e6 && position.norm() < 7e6)
        { origin = "CURRENT_PHYSICAL_USED_REC_POS"; return true; }
    }
    origin = "UNAVAILABLE";
    return false;
}
} // namespace

bool zhangUserRegionalEnabled()
{
    const auto& configuration = settings();
    if (configuration.enabled) validateConfiguration();
    return configuration.enabled;
}

bool zhangUserRegionalNoiseFactors(
    const GTime& time, const std::string& receiver, const SatSys& satellite,
    E_ObsCode observable, bool phase, std::vector<double>& factors,
    std::string& noiseId, int* rank, std::string* reason)
{
    factors.clear(); noiseId.clear(); if (rank) *rank = 0;
    if (!zhangUserRegionalEnabled()) { if (reason) *reason = "REGIONAL_MODE_OFF"; return false; }
    const std::lock_guard<std::recursive_mutex> lock(regionalCacheMutex);
    const auto user = settings().geometry.find(receiver);
    if (user == settings().geometry.end() || user->second.role != "USER") reject("USER_GEOMETRY_MISSING:" + receiver);
    const auto ionosphere = initialStateFromConfig(acsConfig.getRecOpts(receiver).ion_stec);
    if (!ionosphere.estimate || ionosphere.Q < 0) reject("USER_STEC_MUST_BE_ESTIMATED");
    if (!loadEpoch(time)) { if (reason) *reason = cache.reason; return false; }
    const int clock = lookup(cache, "CLOCK", "", satellite.id(), "NONE");
    const int bias = lookup(cache, "PHASE", "", satellite.id(), enum_to_string(observable));
    if (clock < 0 || (phase && bias < 0))
    { if (reason) *reason = "CURRENT_PPP_CLOCK_PHASE_PARAMETER_UNAVAILABLE"; return false; }
    Eigen::RowVectorXd row = cache.root.row(clock);
    if (phase) row -= cache.root.row(bias);
    factors.assign(row.data(), row.data() + row.size()); noiseId = cache.noiseId;
    if (rank) *rank = cache.root.cols(); if (reason) *reason = "NONE";
    return true;
}

void zhangUserRegionalIonoPseudoObs(Trace& trace, KFState& state, KFMeasEntryList& entries)
{
    if (!zhangUserRegionalEnabled() || !regionalSwitch("ZHANG_USER_REGIONAL_IONO",true)) return;
    const std::lock_guard<std::recursive_mutex> lock(regionalCacheMutex);
    struct PrecisionGuard
    {
        Trace& stream;
        std::streamsize precision;
        ~PrecisionGuard() { stream.precision(precision); }
    } precisionGuard{trace, trace.precision()};
    std::map<std::pair<std::string, E_Sys>, std::map<std::string, KFKey>> userKeys;
    for (const auto& entry : entries)
        if (entry.valid && (entry.obsKey.type == KF::CODE_MEAS || entry.obsKey.type == KF::PHAS_MEAS))
            for (const auto& [key, coefficient] : entry.designEntryMap)
                if (coefficient != 0 && key.type == KF::IONO_STEC && zhangUserGecSystem(key.Sat.sys) && key.num == 0 && !key.str.empty())
                    userKeys[{key.str, key.Sat.sys}][key.Sat.id()] = key;
    if (userKeys.empty()) return;
    if (!loadEpoch(state.time))
    {
        trace << "\nZHANG_USER_REGIONAL_EPOCH time=" << state.time.to_string(0)
              << " status=SKIPPED reason=" << cache.reason << " physical_rows_unchanged=1";
        return;
    }
    // Spatial prediction is an independent policy module. This adapter owns
    // neither a plane/LOO implementation nor a scientific accuracy claim.
    for (const auto& [group, keys] : userKeys)
    {
        const auto& [receiver, system] = group;
        const auto user = settings().geometry.find(receiver);
        if (user == settings().geometry.end() || user->second.role != "USER") reject("USER_GEOMETRY_MISSING:" + receiver);
        Eigen::Vector3d operatingPosition = Eigen::Vector3d::Zero();
        std::string positionOrigin;
        if (!physicalOperatingPosition(receiver, entries, state, operatingPosition, positionOrigin))
        {
            trace << "\nZHANG_USER_REGIONAL_EPOCH time=" << state.time.to_string(0)
                  << " receiver=" << receiver << " status=SKIPPED reason=CURRENT_REC_POS_OPERATING_POINT_UNAVAILABLE"
                  << " pseudo_rows=0 user_truth_geometry_used=0 physical_rows_unchanged=1";
            continue;
        }
        const auto initTemplate = initialStateFromConfig(acsConfig.getRecOpts(receiver).ion_stec);
        if (!initTemplate.estimate || initTemplate.Q < 0) reject("USER_STEC_MUST_BE_ESTIMATED");
        zhangregional::Result best;
        zhangregional::Result failed;
        std::string reference;
        std::string failedReference;
        std::size_t bestSupport = 0;
        for (const auto& [candidate, key] : keys)
        {
            const auto pairs = pairSupports(candidate, keys, operatingPosition);
            auto trial = zhangregional::buildSpatialPrediction(cache.mean, cache.covariance, pairs);
            std::size_t support = 0;
            for (const auto& pair : trial.acceptedPairs) support += pair.stations.size();
            if (trial.valid && (trial.acceptedPairs.size() > best.acceptedPairs.size() ||
                (trial.acceptedPairs.size() == best.acceptedPairs.size() && support > bestSupport)))
            { best = std::move(trial); reference = candidate; bestSupport = support; }
            else if (!trial.valid && (failedReference.empty() || trial.diagnostics.size() > failed.diagnostics.size()))
            { failed = std::move(trial); failedReference = candidate; }
        }
        if (reference.empty() || !best.valid || best.prediction.rows() == 0)
        {
            trace << "\nZHANG_USER_REGIONAL_EPOCH time=" << state.time.to_string(0)
                  << " receiver=" << receiver << " status=SKIPPED reason=NO_GEOMETRICALLY_SUPPORTED_SD"
                  << " current_ionosphere_states=" << keys.size() << " pseudo_rows=0 physical_rows_unchanged=1"
                  << " spatial_model=LOCAL_SD_IONO_PLANE_PLUS_NEAREST4_IDW2_RESIDUAL_NOT_MPFG"
                  << " spatial_model_calibrated=0 temporal_error_policy=WHITE_APPROXIMATION cross_epoch_Q=UNAVAILABLE"
                  << " geometry_origin=" << positionOrigin << " user_truth_geometry_used=0 cache_load_sequence=" << cache.loadSequence;
            for (const auto& pair : failed.diagnostics)
                trace << "\nZHANG_USER_REGIONAL_SD_SUPPORT time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " satellite=" << pair.satellite << " reference=" << failedReference
                      << " accepted=0 reason=" << pair.reason << " stations=" << pair.stations.size();
            continue;
        }
        for (const auto& pair : best.diagnostics)
            if (!pair.accepted)
                trace << "\nZHANG_USER_REGIONAL_SD_SUPPORT time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " satellite=" << pair.satellite << " reference=" << reference
                      << " accepted=0 reason=" << pair.reason << " stations=" << pair.stations.size();
        const int pairCount = best.acceptedPairs.size();
        if (best.prediction.rows() != pairCount || best.prediction.cols() != cache.mean.size() ||
            best.spatialCovariance.rows() != pairCount || best.spatialCovariance.cols() != pairCount ||
            !best.prediction.allFinite() || !best.spatialCovariance.allFinite()) reject("SPATIAL_RESULT_DIMENSION_OR_FINITE");
        std::set<std::string> independentTargets;
        for (int row = 0; row < pairCount; ++row)
        {
            const auto& pair = best.acceptedPairs[row];
            if (pair.reference != reference || pair.satellite == reference || !keys.contains(pair.satellite) ||
                !independentTargets.insert(pair.satellite).second) reject("SPATIAL_RESULT_NOT_INDEPENDENT_COMMON_REFERENCE_SD");
            std::map<std::string, double> receiverDatum;
            double targetTotal = 0, referenceTotal = 0;
            for (int column = 0; column < best.prediction.cols(); ++column)
            {
                const double coefficient = best.prediction(row, column);
                if (coefficient == 0) continue;
                const auto& parameter = cache.parameters[column];
                if (parameter.parameter != "IONO_L1_M" ||
                    (parameter.satellite != pair.satellite && parameter.satellite != reference))
                    reject("SPATIAL_PREDICTION_NOT_IONOSPHERE_SD");
                receiverDatum[parameter.receiver] += coefficient;
                if (parameter.satellite == pair.satellite) targetTotal += coefficient;
                else referenceTotal += coefficient;
            }
            if (std::abs(targetTotal - 1) > 1e-9 || std::abs(referenceTotal + 1) > 1e-9 ||
                std::any_of(receiverDatum.begin(), receiverDatum.end(), [](const auto& entry) { return std::abs(entry.second) > 1e-9; }))
                reject("SPATIAL_PREDICTION_RECEIVER_GF_DATUM_NOT_ELIMINATED");
        }
        const VectorXd prediction = best.prediction * cache.mean;
        const MatrixXd formal = best.prediction * cache.covariance * best.prediction.transpose();
        const MatrixXd jointFactors = best.prediction * cache.root;
        double spatialMinimum = 0;
        const MatrixXd spaceFactors = scaledRoot(best.spatialCovariance, spatialMinimum);
        const std::string spaceIdentity = cache.noiseId + ":SPACE:" + receiver + ":REF:" + reference;
        // GPS L1, Galileo E1 and BDS-3 B1C all have 1575.42 MHz.
        // loadEpoch enforces those exact first-signal conventions.
        constexpr double alpha1 = 40.3e16 / (1575.42e6 * 1575.42e6);
        const auto& referenceKey = keys.at(reference);
        int added = 0;
        for (int i = 0; i < static_cast<int>(best.acceptedPairs.size()); ++i)
        {
            const auto& pair = best.acceptedPairs[i];
            const auto& targetKey = keys.at(pair.satellite);
            InitialState left = initTemplate, right = initTemplate;
            if (state.getKFValue(targetKey, left.x) == E_Source::NONE || state.getKFValue(referenceKey, right.x) == E_Source::NONE)
            {
                // New states have not yet passed stateTransition. Preserve
                // exactly the operating values used by the physical rows.
                bool haveLeft = false, haveRight = false;
                for (const auto& row : entries)
                {
                    const auto a = row.usedValueMap.find(targetKey), b = row.usedValueMap.find(referenceKey);
                    if (a != row.usedValueMap.end()) { left.x = a->second; haveLeft = true; }
                    if (b != row.usedValueMap.end()) { right.x = b->second; haveRight = true; }
                }
                if (!haveLeft || !haveRight) reject("CURRENT_IONOSPHERE_OPERATING_POINT_MISSING");
            }
            KFMeasEntry row(&state);
            row.obsKey.type = KF::IONOSPHERIC; row.obsKey.str = receiver;
            row.obsKey.Sat = targetKey.Sat; row.obsKey.num = referenceKey.Sat.prn;
            row.obsKey.comment = "REGIONAL_SD_L1_METRES:" + reference;
            row.addDsgnEntry(targetKey, alpha1, left); row.addDsgnEntry(referenceKey, -alpha1, right);
            row.setInnov(prediction(i) - alpha1 * (left.x - right.x));
            row.metaDataMap["pseudoObs"] = reinterpret_cast<void*>(true);
            row.metaDataMap[ZHANG_USER_REGIONAL_SD_RESIDUAL_MARKER] = reinterpret_cast<void*>(true);
            addFactors(row, cache.noiseId, jointFactors.row(i));
            addFactors(row, spaceIdentity, spaceFactors.row(i));
            entries.push_back(std::move(row)); ++added;
            trace << std::setprecision(17) << "\nZHANG_USER_REGIONAL_SD time=" << state.time.to_string(0)
                  << " receiver=" << receiver << " satellite=" << pair.satellite << " reference=" << reference
                  << " stations=" << pair.stations.size() << " predicted_difference_m=" << prediction(i)
                  << " formal_sigma_m=" << std::sqrt(std::max(0.0, formal(i, i)))
                  << " space_sigma_m=" << std::sqrt(std::max(0.0, best.spatialCovariance(i, i)))
                  << " geometry_condition=" << pair.geometryCondition << " max_distance_km=" << pair.maxDistanceKm
                  << " inside_convex_hull=" << pair.insideConvexHull << " weight_squared_norm=" << pair.weightSquaredNorm
                  << " spatial_model=LOCAL_SD_IONO_PLANE_PLUS_NEAREST4_IDW2_RESIDUAL_NOT_MPFG spatial_model_calibrated=0"
                  << " source_freshness=LAST_DESIGN_NOT_PER_LINK_QC absolute_TEC=0 integer_authority_added=0";
        }
        double clockIonoCross = 0, phaseIonoCross = 0;
        for (int i = 0; i < static_cast<int>(cache.parameters.size()); ++i)
            for (int j = 0; j < static_cast<int>(cache.parameters.size()); ++j)
                if (cache.parameters[j].parameter == "IONO_L1_M")
                {
                    const double square = cache.covariance(i, j) * cache.covariance(i, j);
                    if (cache.parameters[i].parameter == "CLOCK") clockIonoCross += square;
                    if (cache.parameters[i].parameter == "PHASE") phaseIonoCross += square;
                }
        trace << "\nZHANG_USER_REGIONAL_EPOCH time=" << state.time.to_string(0)
              << " receiver=" << receiver << " system=" << enum_to_string(system) << " status=APPLIED reference=" << reference << " pairs=" << added
              << " exact_source_Q_dimensions=" << cache.sourceDimension << " retained_Q_dimensions=" << cache.mean.size()
              << " joint_factor_columns=" << cache.root.cols() << " clock_iono_cross_norm_m2=" << std::sqrt(clockIonoCross)
              << " phase_iono_cross_norm_m2=" << std::sqrt(phaseIonoCross)
              << " min_scaled_Q_eigenvalue=" << cache.minimumScaledEigenvalue
              << " min_scaled_spatial_eigenvalue=" << spatialMinimum << " full_joint_R_shared=1"
              << " spatial_model=LOCAL_SD_IONO_PLANE_PLUS_NEAREST4_IDW2_RESIDUAL_NOT_MPFG spatial_model_calibrated=0"
              << " geometry_origin=" << positionOrigin << " user_truth_geometry_used=0"
              << " operating_x_m=" << operatingPosition(0) << " operating_y_m=" << operatingPosition(1)
              << " operating_z_m=" << operatingPosition(2) << " cache_load_sequence=" << cache.loadSequence
              << " receiver_common_GF_direction_preserved=1 temporal_error_policy=WHITE_APPROXIMATION"
              << " temporal_cross_covariance=UNAVAILABLE temporal_independence_certified=0 nis_evaluations=0";
    }
}

#include "pea/zhangUserRegionalTrop.inc"

#include "pea/zhangRegionalCheckpoint.inc"

bool prepareZhangRealtimeRegional(const GTime& time)
{
    std::lock_guard<std::recursive_mutex> lock(regionalCacheMutex);
    return !settings().enabled || loadEpoch(time);
}
