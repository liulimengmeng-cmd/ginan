// Small production-classifier contract, no parse, KF or PEA main invocation.
// Link as an independent driver against the frozen complete production objects.
#include "common/zhangMultiGnssFloatCheckpoint.hpp"
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>

namespace
{
std::size_t checks = 0;
void require(bool value, const std::string& reason)
{ ++checks; if (!value) throw std::runtime_error(reason); }
struct Environment
{
    std::map<std::string, std::optional<std::string>> saved;
    void set(const std::string& key, const char* value)
    {
        if (!saved.count(key))
        {
            const auto* old = std::getenv(key.c_str());
            saved[key] = old ? std::optional<std::string>(old) : std::nullopt;
        }
        if (value) setenv(key.c_str(), value, 1); else unsetenv(key.c_str());
    }
    ~Environment()
    {
        for (const auto& [key, value] : saved)
            if (value) setenv(key.c_str(), value->c_str(), 1); else unsetenv(key.c_str());
    }
};
std::unique_ptr<ACSConfig> configuration()
{
    auto config = std::make_unique<ACSConfig>();
    config->process_sys.clear();
    config->zhangFullRank.sysOpts.clear();
    config->code_priorities.clear();
    config->zhangPppAr.baseline_observables.clear();
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        std::vector<E_ObsCode> signals;
        switch (system)
        {
            case E_Sys::GPS: signals = {E_ObsCode::L1C, E_ObsCode::L2W}; break;
            case E_Sys::GAL: signals = {E_ObsCode::L1C, E_ObsCode::L5Q}; break;
            case E_Sys::BDS: signals = {E_ObsCode::L1P, E_ObsCode::L5P}; break;
            default: throw std::runtime_error("Unexpected test constellation");
        }
        config->process_sys[system] = true;
        config->zhangFullRank.sysOpts[system].baseline_observables = signals;
        config->code_priorities[system] = signals;
        config->zhangPppAr.baseline_observables[system] = signals;
    }
    config->epoch_interval = 30;
    config->zhangFullRank.enable = true;
    config->zhangPppAr.output_products = true;
    config->zhangPppAr.user_adapter = false;
    config->zhangPppAr.transactional_integer_fixing = true;
    config->zhangPppAr.product_relation_feedback = false;
    config->start_epoch = boost::posix_time::ptime(boost::gregorian::date(2024, 7, 17));
    config->end_epoch = config->start_epoch + boost::posix_time::hours(2);
    config->zhangPppAr.checkpoint_runtime_id = ZHANG_GEC134_SOURCE_RUNTIME_ID;
    return config;
}
}
int main()
{
    try
    {
        Environment environment;
        for (const auto key : {"ZHANG_MULTI_GNSS_AR_CHECKPOINT", "ZHANG_R51_ENABLE",
                "ZHANG_R51_RATIO_ONLY", "ZHANG_R51_EPOCH_ONLY_AR", "ZHANG_GEC134_FLOAT_TO_AR_MIGRATION",
                ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG}) environment.set(key, "1");
        environment.set("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT", nullptr);
        environment.set("ZHANG_R51_RATIO_THRESHOLD", "2");
        environment.set("ZHANG_R51_AR_START_GPST_SECONDS", ZHANG_GEC134_AR_START);
        environment.set(ZHANG_GEC134_AR_CONTINUATION_FLAG, nullptr);
        auto config = configuration(); std::string reason;
        require(zhangCheckpointModeForConfig(*config, reason) == ZHANG_GEC134_AR_MODE,
            "Original 02:00:00 policy changed");
        config->end_epoch += boost::posix_time::seconds(30);
        require(zhangCheckpointModeForConfig(*config, reason).empty(),
            "Old profile silently widened to 02:00:30");
        environment.set(ZHANG_GEC134_AR_CONTINUATION_FLAG, "1");
        config->max_epochs = 242;
        config->outputs_root = std::string(ZHANG_GEC134_AR_CONTINUATION_ROOT) + "/outputs";
        config->zhangPppAr.checkpoint_restore_path = ZHANG_GEC134_AR_E122_PATH;
        config->zhangPppAr.checkpoint_capture_epochs = {"2024-07-17 01:30:30", "2024-07-17 02:00:30"};
        require(zhangCheckpointModeForConfig(*config, reason) == ZHANG_GEC134_AR_MODE && reason.empty(),
            "Explicit bounded local hour rejected: " + reason);
        config->max_epochs = 243;
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Extra epoch admitted");
        config->max_epochs = 242;
        config->end_epoch += boost::posix_time::seconds(30);
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Extra time admitted");
        config->end_epoch -= boost::posix_time::seconds(60);
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Partial explicit hour admitted");
        config->end_epoch += boost::posix_time::seconds(30);
        config->zhangPppAr.checkpoint_capture_epochs[0] = "2024-07-17 01:30:00";
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Wrong capture admitted");
        config->zhangPppAr.checkpoint_capture_epochs[0] = "2024-07-17 01:30:30";
        config->outputs_root += "_other";
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Wrong root admitted");
        config->outputs_root = std::string(ZHANG_GEC134_AR_CONTINUATION_ROOT) + "/outputs";
        config->zhangPppAr.checkpoint_restore_path = ZHANG_GEC134_SOURCE_CHECKPOINT_PATH;
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "FLOAT source profile admitted");
        config->zhangPppAr.checkpoint_restore_path = ZHANG_GEC134_AR_E122_PATH;
        environment.set("ZHANG_R51_RATIO_THRESHOLD", "3");
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Ratio change admitted");
        environment.set("ZHANG_R51_RATIO_THRESHOLD", "2");
        environment.set("ZHANG_R51_EPOCH_ONLY_AR", "0");
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Temporal AR admitted");
        environment.set("ZHANG_R51_EPOCH_ONLY_AR", "1");
        for (const bool pilot : {true, false})
        {
            const std::string root = pilot ?
                "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01" :
                "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_hour_20261025_01";
            config->end_epoch = config->start_epoch + (pilot ? boost::posix_time::seconds(4020) : boost::posix_time::seconds(7230));
            config->max_epochs = pilot ? 135 : 242;
            config->outputs_root = root + "/outputs";
            config->zhangPppAr.checkpoint_capture_epochs = pilot ? std::vector<std::string>{"2024-07-17 01:07:00"} :
                std::vector<std::string>{"2024-07-17 01:30:30", "2024-07-17 02:00:30"};
            require(zhangCheckpointModeForConfig(*config, reason) == ZHANG_GEC134_AR_MODE,
                "Exact isolated V6 profile rejected: " + reason);
            ++config->max_epochs;
            require(zhangCheckpointModeForConfig(*config, reason).empty(), "V6 extra epoch accepted");
            --config->max_epochs;
            config->end_epoch += boost::posix_time::seconds(30);
            require(zhangCheckpointModeForConfig(*config, reason).empty(), "V6 extra time accepted");
            config->end_epoch -= boost::posix_time::seconds(30);
            config->outputs_root += "_other";
            require(zhangCheckpointModeForConfig(*config, reason).empty(), "V6 unbounded root accepted");
            config->outputs_root = root + "/outputs";
            config->zhangPppAr.checkpoint_capture_epochs[0] = "2024-07-17 01:06:30";
            require(zhangCheckpointModeForConfig(*config, reason).empty(), "V6 wrong capture accepted");
        }
        // Independent literals, not values read back from the profile table.
        // This classifier grants no old-v6-e135 cross-binary restore token.
        for (const bool pilot : {true,false})
        {
            const std::string root = pilot ?
                "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v8_pilot_20261025_01" :
                "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v8_hour_20261025_01";
            config->end_epoch = config->start_epoch + (pilot ? boost::posix_time::seconds(4020) : boost::posix_time::seconds(7230));
            config->max_epochs = pilot ? 135 : 242;
            config->outputs_root = root + "/outputs";
            config->zhangPppAr.checkpoint_restore_path = ZHANG_GEC134_AR_E122_PATH;
            const std::vector<std::string> captures = pilot ? std::vector<std::string>{"2024-07-17 01:07:00"} :
                std::vector<std::string>{"2024-07-17 01:30:30","2024-07-17 02:00:30"};
            config->zhangPppAr.checkpoint_capture_epochs = captures;
            require(zhangCheckpointModeForConfig(*config,reason) == ZHANG_GEC134_AR_MODE && reason.empty(),
                "Exact isolated V8 profile rejected: " + reason);
            ++config->max_epochs;
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 extra epoch accepted");
            config->max_epochs -= 2;
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 partial schedule accepted");
            ++config->max_epochs;
            for (const int offset : {-30,30})
            {
                config->end_epoch += boost::posix_time::seconds(offset);
                require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 changed end epoch accepted");
                config->end_epoch -= boost::posix_time::seconds(offset);
            }
            config->outputs_root += "_other";
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 nonliteral root accepted");
            config->outputs_root = root + "/outputs";
            config->zhangPppAr.checkpoint_capture_epochs[0] = "2024-07-17 01:06:30";
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 wrong capture epoch accepted");
            config->zhangPppAr.checkpoint_capture_epochs = captures;
            config->zhangPppAr.checkpoint_capture_epochs.push_back("2024-07-17 02:01:00");
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 additional capture accepted");
            config->zhangPppAr.checkpoint_capture_epochs = captures;
            for (const auto* source : {ZHANG_GEC134_SOURCE_CHECKPOINT_PATH,
                "/home/rx/GINAN/multi_gnss_float_profile_20261022/server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01/checkpoints/E29-20240717010700-e135-r32f626de9463"})
            {
                config->zhangPppAr.checkpoint_restore_path = source;
                require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 non-e122 configured source accepted");
            }
            config->zhangPppAr.checkpoint_restore_path = ZHANG_GEC134_AR_E122_PATH;
            for (const auto* key : {"ZHANG_R51_RATIO_THRESHOLD","ZHANG_R51_RATIO_ONLY","ZHANG_R51_EPOCH_ONLY_AR"})
            {
                environment.set(key,"0");
                require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 changed science policy accepted");
                environment.set(key,std::string(key) == "ZHANG_R51_RATIO_THRESHOLD" ? "2" : "1");
            }
            config->epoch_interval = 60;
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 changed epoch interval accepted");
            config->epoch_interval = 30;
            config->zhangPppAr.checkpoint_runtime_id += "!";
            require(zhangCheckpointModeForConfig(*config,reason).empty(),"V8 changed runtime owner accepted");
            config->zhangPppAr.checkpoint_runtime_id = ZHANG_GEC134_SOURCE_RUNTIME_ID;
            require(zhangCheckpointModeForConfig(*config,reason) == ZHANG_GEC134_AR_MODE,
                "V8 exact positive failed after rejection matrix");
        }
        environment.set(ZHANG_GEC134_LOCAL_ONE_EPOCH_FLAG, nullptr);
        require(zhangCheckpointModeForConfig(*config, reason).empty(), "Nonlocal profile widened");
        std::cout << "GEC134_AR_CONTINUATION_CLASSIFIER_PASSED checks=" << checks
            << " old_upper_bound_preserved=1 no_PEA_or_KF=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GEC134_AR_CONTINUATION_CLASSIFIER_FAIL " << error.what() << '\n';
        return 1;
    }
}
