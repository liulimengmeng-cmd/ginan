#include "common/zhangMultiGnssFloatCheckpoint.hpp"
#include "common/navigation.hpp"
#include "pea/zhangCheckpointRuntime.hpp"

#include <algorithm>
#include <iostream>
#include <initializer_list>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>

namespace
{
std::size_t assertions = 0;
void require(bool condition, const std::string& message)
{
    ++assertions;
    if (!condition) throw std::runtime_error(message);
}

struct Environment
{
    std::map<std::string, std::optional<std::string>> previous;
    void set(const std::string& name, const char* value)
    {
        if (!previous.count(name))
        {
            const auto* old = std::getenv(name.c_str());
            previous[name] = old ? std::optional<std::string>(old) : std::nullopt;
        }
        if (value) setenv(name.c_str(), value, 1);
        else unsetenv(name.c_str());
    }
    ~Environment()
    {
        for (const auto& [name, value] : previous)
            if (value) setenv(name.c_str(), value->c_str(), 1);
            else unsetenv(name.c_str());
    }
};

void validEnvironment(Environment& env)
{
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "1");
    env.set("ZHANG_R51_ENABLE", "1");
    env.set("ZHANG_R51_RATIO_ONLY", "1");
    env.set("ZHANG_R51_EPOCH_ONLY_AR", "1");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1405216830");
}

void addSystem(ACSConfig& config, E_Sys system)
{
    std::vector<E_ObsCode> signals;
    switch (system)
    {
        case E_Sys::GPS: signals = {E_ObsCode::L1C, E_ObsCode::L2W}; break;
        case E_Sys::GAL: signals = {E_ObsCode::L1C, E_ObsCode::L5Q}; break;
        case E_Sys::BDS: signals = {E_ObsCode::L1P, E_ObsCode::L5P}; break;
        default: signals = {E_ObsCode::L1C, E_ObsCode::L2W}; break;
    }
    config.process_sys[system] = true;
    config.zhangFullRank.sysOpts[system].baseline_observables = signals;
    config.code_priorities[system] = signals;
    config.zhangPppAr.baseline_observables[system] = signals;
}

// ACSConfig owns a mutex and is intentionally non-copyable. This test fixture
// copies only the public inputs consumed by the policy, never the mutex or any
// producer runtime state. The production function still sees actual ACSConfig.
void copyPolicyFields(ACSConfig& target, const ACSConfig& source)
{
    target.process_sys = source.process_sys;
    target.zhangFullRank = source.zhangFullRank;
    target.code_priorities = source.code_priorities;
    target.zhangPppAr = source.zhangPppAr;
    target.epoch_interval = source.epoch_interval;
    target.start_epoch = source.start_epoch;
    target.end_epoch = source.end_epoch;
    target.max_epochs = source.max_epochs;
    target.outputs_root = source.outputs_root;
    target.ambrOpts = source.ambrOpts;
}

struct FixtureConfig : ACSConfig
{
    FixtureConfig() = default;
    FixtureConfig(const FixtureConfig& other) { *this = other; }
    FixtureConfig& operator=(const FixtureConfig& other)
    {
        copyPolicyFields(*this, other);
        return *this;
    }
};

FixtureConfig validConfig()
{
    FixtureConfig config;
    config.process_sys.clear();
    config.zhangFullRank.sysOpts.clear();
    config.code_priorities.clear();
    config.zhangPppAr.baseline_observables.clear();
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS}) addSystem(config, system);
    config.epoch_interval = 30;
    config.zhangFullRank.enable = true;
    config.zhangPppAr.output_products = true;
    config.zhangPppAr.user_adapter = false;
    config.zhangPppAr.transactional_integer_fixing = true;
    config.zhangPppAr.product_relation_feedback = false;
    config.start_epoch = boost::posix_time::ptime(boost::gregorian::date(2024, 7, 17));
    config.end_epoch = config.start_epoch + boost::posix_time::hours(2);
    return config;
}

void accepted(const ACSConfig& config, const std::string& mode, const std::string& label)
{
    std::string failure = "STALE_FAILURE_MUST_CLEAR";
    require(zhangCheckpointModeForConfig(config, failure) == mode, label + ": mode changed");
    require(failure.empty(), label + ": stale failure retained");
}

void rejected(const ACSConfig& config, const std::string& reason, const std::string& label)
{
    std::string failure;
    require(zhangCheckpointModeForConfig(config, failure).empty(), label + ": incorrectly accepted");
    require(!failure.empty(), label + ": rejection lacks reason");
    if (!reason.empty()) require(failure == reason, label + ": unexpected reason " + failure);
}

void legacyGpsContract()
{
    Environment env;
    for (const auto name : {"ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "ZHANG_R51_ENABLE",
                           "ZHANG_R51_RATIO_ONLY", "ZHANG_R51_EPOCH_ONLY_AR",
                           "ZHANG_R51_AR_START_GPST_SECONDS"}) env.set(name, nullptr);
    auto config = validConfig();
    config.process_sys[E_Sys::GAL] = false;
    config.process_sys[E_Sys::BDS] = false;
    require(!zhangMultiGnssFloatCheckpointActive(config), "missing opt-in activates legacy GPS");
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "1");
    require(!zhangMultiGnssFloatCheckpointActive(config), "disabled non-GPS systems activate multi mode");
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", nullptr);
    const std::string legacy = "E29_GPS_L1C_L2W_ZHANG_FULL_RANK";
    accepted(config, legacy, "legacy GPS byte-exact mode without opt-in");
    // The new strict FLOAT contract must not retroactively tighten legacy GPS
    // schedule, output, authority or search-priority settings.
    config.zhangFullRank.enable = false;
    config.zhangPppAr.output_products = false;
    config.zhangPppAr.user_adapter = true;
    config.zhangPppAr.transactional_integer_fixing = false;
    config.zhangPppAr.product_relation_feedback = true;
    config.code_priorities[E_Sys::GPS] = {E_ObsCode::L1X, E_ObsCode::L2X};
    config.zhangPppAr.baseline_observables.clear();
    config.start_epoch = boost::posix_time::ptime(boost::posix_time::not_a_date_time);
    config.end_epoch = boost::posix_time::ptime(boost::posix_time::not_a_date_time);
    accepted(config, legacy, "legacy unrelated settings preserved");
    config.epoch_interval = 60;
    accepted(config, legacy, "legacy 60 second mode");
    config.epoch_interval = std::numeric_limits<double>::quiet_NaN();
    accepted(config, legacy, "new finite-interval check must remain multi-only");
    config.epoch_interval = 60;
    config.zhangFullRank.sysOpts[E_Sys::GPS].baseline_observables = {E_ObsCode::L2W, E_ObsCode::L1C};
    rejected(config, "CHECKPOINT_ORDERED_BASELINES_MISMATCH", "legacy original ordered-baseline gate");
    std::cout << "PASS legacy GPS byte-name compatibility and unrelated-setting boundaries\n";
}

void multiSystemContract()
{
    Environment env;
    validEnvironment(env);
    const std::string mode = "E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS";
    auto config = validConfig();
    require(zhangMultiGnssFloatCheckpointActive(config), "explicit GEC opt-in not active");
    accepted(config, mode, "joint GEC exact mode");
    config.epoch_interval = 60;
    accepted(config, mode, "joint 60 second mode");
    config.epoch_interval = 30;
    // Equal mode for different insertion order, independent of enum-map order.
    config.process_sys.clear();
    for (const auto system : {E_Sys::BDS, E_Sys::GPS, E_Sys::GAL}) config.process_sys[system] = true;
    accepted(config, mode, "joint deterministic system ordering");
    for (const char* value : std::initializer_list<const char*>{nullptr, "0", "true", "01", "1 "})
    {
        env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", value);
        require(!zhangMultiGnssFloatCheckpointActive(config), "non-exact opt-in activates multi mode");
        rejected(config, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_EXPLICIT_OPT_IN", "explicit opt-in required");
    }
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "1");
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        auto wrong = config;
        std::reverse(wrong.zhangFullRank.sysOpts[system].baseline_observables.begin(),
                     wrong.zhangFullRank.sysOpts[system].baseline_observables.end());
        rejected(wrong, "CHECKPOINT_ORDERED_BASELINES_MISMATCH", "ordered physical baseline");
        wrong = config;
        wrong.code_priorities[system][0] = E_ObsCode::L1X;
        rejected(wrong, "FLOAT_CHECKPOINT_EXACT_SIGNAL_PRIORITIES_MISMATCH", "no X signal substitution");
        wrong = config;
        wrong.zhangPppAr.baseline_observables[system].clear();
        rejected(wrong, "FLOAT_CHECKPOINT_EXACT_SIGNAL_PRIORITIES_MISMATCH", "missing product baseline");
    }
    auto wrong = config;
    addSystem(wrong, E_Sys::GLO);
    rejected(wrong, "FLOAT_CHECKPOINT_UNSUPPORTED_SYSTEM", "GLO unsupported");
    wrong = config;
    wrong.process_sys.clear();
    rejected(wrong, "CHECKPOINT_NO_ACTIVE_SYSTEM", "empty system set");
    wrong = config;
    wrong.epoch_interval = 15;
    rejected(wrong, "CHECKPOINT_REQUIRES_30_OR_60_SECOND_INTERVAL", "unsupported interval");
    for (const double interval : {std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::infinity(),
                                  -std::numeric_limits<double>::infinity()})
    {
        wrong = config;
        wrong.epoch_interval = interval;
        rejected(wrong, "CHECKPOINT_REQUIRES_30_OR_60_SECOND_INTERVAL", "nonfinite multi interval");
    }
    std::cout << "PASS explicit opt-in, ordered GEC signals, unsupported-system and interval gates\n";
}

void authorityAndScheduleContract()
{
    Environment env;
    validEnvironment(env);
    auto config = validConfig();
    const auto rejectAuthority = [&](const ACSConfig& wrong, const char* label)
    { rejected(wrong, "MULTI_GNSS_CHECKPOINT_REQUIRES_ISOLATED_FLOAT_AUTHORITY", label); };
    auto wrong = config;
    wrong.zhangPppAr.product_relation_feedback = true;
    rejectAuthority(wrong, "feedback enabled");
    wrong = config;
    wrong.zhangFullRank.enable = false;
    rejectAuthority(wrong, "full-rank disabled");
    wrong = config;
    wrong.zhangPppAr.output_products = false;
    rejectAuthority(wrong, "products disabled");
    wrong = config;
    wrong.zhangPppAr.user_adapter = true;
    rejectAuthority(wrong, "user adapter");
    wrong = config;
    wrong.zhangPppAr.transactional_integer_fixing = false;
    rejectAuthority(wrong, "no isolated transaction authority");
    for (const auto name : {"ZHANG_R51_ENABLE", "ZHANG_R51_RATIO_ONLY", "ZHANG_R51_EPOCH_ONLY_AR"})
    {
        for (const char* value : std::initializer_list<const char*>{nullptr, "0", "true"})
        {
            env.set(name, value);
            rejectAuthority(config, "R51 isolation flag not exactly one");
        }
        env.set(name, "1");
    }

    // Independent calendar constants, not a second call to production GTime or
    // its epoch_seconds helper. A Unix-origin regression would reject 1405216830.
    constexpr long long gpstEnd = 1405216800;
    constexpr long long unixEnd = 1721181600;
    require(unixEnd - gpstEnd == 315964800, "calendar-origin fixture changed");
    accepted(config, "E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS", "GPST 1980-origin end plus 30 seconds");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1405216800");
    rejected(config, "MULTI_GNSS_FLOAT_CHECKPOINT_END_MUST_PRECEDE_AR", "end equals AR start");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1405216770");
    rejected(config, "MULTI_GNSS_FLOAT_CHECKPOINT_END_MUST_PRECEDE_AR", "AR begins before end");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", nullptr);
    rejected(config, "", "missing AR schedule must not authorise multi-GNSS checkpoint");
    for (const char* value : {"", "nan", "inf", "-1", "0", "1405216830.5", "1405216830x"})
    {
        env.set("ZHANG_R51_AR_START_GPST_SECONDS", value);
        rejected(config, "MULTI_GNSS_FLOAT_CHECKPOINT_INVALID_AR_SCHEDULE", "invalid AR schedule");
    }
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1405216830");
    wrong = config;
    wrong.start_epoch = boost::posix_time::ptime(boost::posix_time::not_a_date_time);
    rejected(wrong, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_BOUNDED_EPOCHS", "missing start epoch");
    wrong = config;
    wrong.end_epoch = boost::posix_time::ptime(boost::posix_time::not_a_date_time);
    rejected(wrong, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_BOUNDED_EPOCHS", "missing end epoch");
    wrong = config;
    wrong.end_epoch = wrong.start_epoch - boost::posix_time::seconds(30);
    rejected(wrong, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_BOUNDED_EPOCHS", "reversed epochs");
    for (const auto special : {boost::posix_time::pos_infin, boost::posix_time::neg_infin})
    {
        wrong = config;
        wrong.start_epoch = boost::posix_time::ptime(special);
        rejected(wrong, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_BOUNDED_EPOCHS", "infinite start epoch");
        wrong = config;
        wrong.end_epoch = boost::posix_time::ptime(special);
        rejected(wrong, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_BOUNDED_EPOCHS", "infinite end epoch");
    }
    std::cout << "PASS isolated FLOAT authority, GPST-1980 origin and strict pre-AR bounded schedule\n";
}

void publicRuntimeModeGateContract()
{
    Environment env;
    validEnvironment(env);
    struct GlobalPolicyRestore
    {
        FixtureConfig saved;
        GlobalPolicyRestore() { copyPolicyFields(saved, acsConfig); }
        ~GlobalPolicyRestore() { copyPolicyFields(acsConfig, saved); }
    } restore;
    auto config = validConfig();
    copyPolicyFields(acsConfig, config);
    KFState state;
    const double epoch[6] = {2024, 7, 17, 0, 0, 0};
    state.time = epoch2time(epoch, E_TimeSys::GPST);
    state.lastPppTransactionCommitted = true;
    state.lastPppQcConverged = true;
    require(bindZhangCheckpointRuntimeId(state, "POLICY-ONLY-RUNTIME"), "public mode fixture runtime binding");
    ZhangPeaControllerCheckpointState controller;
    require(makeZhangPeaPostEpochCheckpointState(1, captureZhangCheckpointTime(state.time), 30, controller).valid,
            "public mode fixture controller");
    const ReceiverMap receivers;
    const Navigation navigation{};
    const std::multimap<std::string, StreamParserPtr> streams;
    const std::map<std::string, bool> done;
    ZhangE29CheckpointIdentity identity;
    identity.runtimeId = "POLICY-ONLY-RUNTIME";
    identity.checkpointId = "POLICY-ONLY-CHECKPOINT";
    const auto check = [&](const char* mode, bool modeMustPass, const char* label)
    {
        ZhangE29CheckpointProvenance provenance;
        provenance.experimentMode = mode;
        // Deliberately absent executable/config/input provenance is an
        // independent downstream sentinel: it proves the real public mode gate
        // was passed, without pretending to capture a scientifically valid KF.
        ZhangCheckpointBundle bundle;
        const auto capture = captureZhangE29CheckpointBundle(state, receivers, navigation,
            streams, done, controller, identity, provenance, bundle);
        ZhangE29CheckpointRestorePlan plan;
        const auto preflight = preflightZhangE29CheckpointBundle(bundle, provenance, state,
            receivers, navigation, streams, done, 30, plan);
        for (const auto* result : {&capture, &preflight})
        {
            require(!result->valid, std::string(label) + ": incomplete fixture incorrectly accepted");
            if (modeMustPass)
                require(result->failureReason == "E29_CHECKPOINT_PROVENANCE_INPUT_MISSING",
                        std::string(label) + ": real public mode gate did not pass: " + result->failureReason);
            else
                require(result->failureReason.rfind("E29_CHECKPOINT_EXPERIMENT_MODE_UNSUPPORTED", 0) == 0,
                        std::string(label) + ": wrong public rejection: " + result->failureReason);
        }
        require(bundle.sections.empty(), std::string(label) + ": early gate mutated output bundle");
        require(!plan.valid, std::string(label) + ": early gate created restore plan");
    };
    const char* joint = "E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS";
    check(joint, true, "actual capture/preflight exact live GEC mode");
    check("E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS_FAKE", false, "no prefix-based mode authorisation");
    check("E29_MULTI_GNSS_FLOAT_V2_GPS_GAL", false, "mode must equal live configuration");
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", nullptr);
    check(joint, false, "actual capture/preflight requires opt-in");
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "1");
    acsConfig.process_sys[E_Sys::BDS] = false;
    check(joint, false, "actual active system set mismatch");
    acsConfig.process_sys[E_Sys::GAL] = false;
    env.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", nullptr);
    env.set("ZHANG_R51_EPOCH_ONLY_AR", nullptr);
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", nullptr);
    check("E29_GPS_L1C_L2W_ZHANG_FULL_RANK", true, "legacy public mode unchanged");
    std::cout << "PASS actual capture/preflight public mode gates: exact live mode, no prefix, opt-in, legacy\n";
}
void us20ColdClassifierContract()
{
    Environment env;
    for (const auto* key : {"ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", "ZHANG_MULTI_GNSS_AR_CHECKPOINT",
            "ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", "ZHANG_GEC134_AR_CONTINUATION_12_6",
            "ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6", "ZHANG_GEC134_LOCAL_ONE_EPOCH_MIGRATION",
            "ZHANG_R51_REUSE_FLOAT_20240717", "ZHANG_R51_REUSE_RATIO_20240717"}) env.set(key, nullptr);
    for (const auto* key : {ZHANG_US20_GLOBAL134_CHECKPOINT_FLAG, "ZHANG_R51_ENABLE",
            "ZHANG_R51_RATIO_ONLY", "ZHANG_R51_EPOCH_ONLY_AR"}) env.set(key, "1");
    env.set(ZHANG_US20_GLOBAL134_STAGE_FLAG, "WARMUP");
    env.set("ZHANG_R51_RATIO_THRESHOLD", "2");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1405216800");
    auto config = validConfig();
    config.outputs_root = string(ZHANG_US20_GLOBAL134_ROOT) + "/outputs";
    config.end_epoch = config.start_epoch + boost::posix_time::hours(7) - boost::posix_time::seconds(30);
    config.max_epochs = 840;
    config.zhangPppAr.deterministic_checkpoint = true;
    config.zhangPppAr.checkpoint_runtime_id = ZHANG_US20_GLOBAL134_RUNTIME_ID;
    config.zhangPppAr.checkpoint_restore_path.clear();
    config.zhangPppAr.checkpoint_output_directory = string(ZHANG_US20_GLOBAL134_ROOT) + "/checkpoints";
    config.zhangPppAr.product_filename = config.outputs_root + "/zhang_internal_products.csv";
    config.zhangPppAr.product_covariance_filename = config.outputs_root + "/zhang_internal_product_covariance.csv";
    config.zhangPppAr.checkpoint_capture_epochs = {"2024-07-17 01:59:30", "2024-07-17 02:59:30",
        "2024-07-17 03:59:30", "2024-07-17 04:59:30", "2024-07-17 05:59:30", "2024-07-17 06:59:30"};
    config.zhangFullRank.product_basis_mode = "SINGLE_S_BASIS_NATIVE";
    require(zhangMultiGnssArCheckpointActive(config), "new cold flag did not activate nine-section container");
    accepted(config, ZHANG_GEC134_AR_MODE, "cold single-S literal main profile");
    auto wrong = config; wrong.max_epochs++;
    rejected(wrong, "US20_GLOBAL134_COLD_EXACT_PROFILE_REQUIRED", "new profile extra epoch");
    wrong = config; wrong.end_epoch += boost::posix_time::seconds(30);
    rejected(wrong, "US20_GLOBAL134_COLD_EXACT_PROFILE_REQUIRED", "new profile extra end time");
    wrong = config; wrong.zhangFullRank.product_basis_mode = "LEGACY_PRODUCT_TREE";
    rejected(wrong, "US20_GLOBAL134_COLD_EXACT_PROFILE_REQUIRED", "new profile legacy product basis");
    wrong = config; wrong.zhangPppAr.checkpoint_restore_path = ZHANG_GEC134_AR_E122_PATH;
    rejected(wrong, "US20_GLOBAL134_COLD_EXACT_PROFILE_REQUIRED", "new profile frozen YAML old checkpoint");
    wrong = config; wrong.zhangPppAr.checkpoint_capture_epochs.erase(wrong.zhangPppAr.checkpoint_capture_epochs.begin());
    rejected(wrong, "US20_GLOBAL134_COLD_EXACT_PROFILE_REQUIRED", "new profile missing warmup checkpoint");
    env.set("ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", "1");
    rejected(config, "US20_GLOBAL134_LEGACY_CAPABILITY_FORBIDDEN", "new profile must not grant old migration");
    env.set("ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", nullptr);
    env.set(ZHANG_US20_GLOBAL134_STAGE_FLAG, "AUTO");
    rejected(config, "US20_GLOBAL134_STAGE_REJECTED", "unknown lifecycle stage");
    env.set(ZHANG_US20_GLOBAL134_STAGE_FLAG, "WARMUP");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1e100");
    rejected(config, "US20_GLOBAL134_COLD_EXACT_PROFILE_REQUIRED", "large finite AR schedule cannot be narrowed");
    env.set("ZHANG_R51_AR_START_GPST_SECONDS", "1405216800");
    config.outputs_root = string(ZHANG_US20_GLOBAL134_PROBE_ROOT) + "/outputs";
    config.end_epoch = config.start_epoch;
    config.max_epochs = 1;
    config.zhangPppAr.checkpoint_output_directory = string(ZHANG_US20_GLOBAL134_PROBE_ROOT) + "/checkpoints";
    config.zhangPppAr.product_filename = config.outputs_root + "/zhang_internal_products.csv";
    config.zhangPppAr.product_covariance_filename = config.outputs_root + "/zhang_internal_product_covariance.csv";
    config.zhangPppAr.checkpoint_capture_epochs = {"2024-07-17 00:00:00"};
    accepted(config, ZHANG_GEC134_AR_MODE, "literal one-epoch cold probe");
    env.set(ZHANG_US20_GLOBAL134_CHECKPOINT_FLAG, nullptr);
    rejected(config, "MULTI_GNSS_FLOAT_CHECKPOINT_REQUIRES_EXPLICIT_OPT_IN", "cold capability is explicit");
    std::cout << "PASS US20 cold main/probe production classifier and preserved capability separation\n";
}
}

int main()
{
    try
    {
        legacyGpsContract();
        multiSystemContract();
        authorityAndScheduleContract();
        publicRuntimeModeGateContract();
        us20ColdClassifierContract();
        std::cout << "PASS multi-GNSS FLOAT checkpoint policy assertions=" << assertions << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL checkpoint policy: " << error.what() << '\n';
        return 1;
    }
}
