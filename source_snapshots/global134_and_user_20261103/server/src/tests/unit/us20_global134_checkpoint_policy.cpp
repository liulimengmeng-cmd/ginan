// Independent literal fixtures: no production profile table generates the
// expected schedules. No PEA/KF/checkpoint payload is executed by this test.
#include "common/zhangUs20Global134Checkpoint.hpp"
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
int checks = 0;
void require(bool condition, const std::string& message)
{
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
ZhangUs20Global134CheckpointPolicy mainPolicy()
{
    ZhangUs20Global134CheckpointPolicy p;
    p.outputsRoot = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/outputs";
    p.binaryPath = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/candidate/pea";
    p.startEpoch = "2024-07-17 00:00:00"; p.endEpoch = "2024-07-17 06:59:30";
    p.productFilename = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/outputs/zhang_internal_products.csv";
    p.covarianceFilename = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/outputs/zhang_internal_product_covariance.csv";
    p.checkpointOutputDirectory = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/checkpoints";
    p.checkpointRuntimeId = "E29-GEC134-FIXEDCOORD-20261023-one_hour";
    p.productBasisMode = "SINGLE_S_BASIS_NATIVE";
    p.systems = {"GPS", "GAL", "BDS"};
    p.captureEpochs = {"2024-07-17 01:59:30", "2024-07-17 02:59:30",
        "2024-07-17 03:59:30", "2024-07-17 04:59:30", "2024-07-17 05:59:30", "2024-07-17 06:59:30"};
    p.receiverCount = 134; p.epochIntervalSeconds = 30; p.maxEpochs = 840;
    p.arStartGpstSeconds = 1405216800; p.ratio = 2; p.jointNisEnabled = false;
    p.independentEpochAr = true; p.fixedCoordinates = true;
    p.physicalIonosphereEnabled = true; p.physicalTroposphereEnabled = true;
    p.stecEstimated = true; p.ztdEstimated = true;
    p.regionalJointEnabled = false; p.regionalTropEnabled = false;
    p.environment = {{"ZHANG_US20_GLOBAL134_CHECKPOINT", "1"}, {"ZHANG_US20_GLOBAL134_STAGE", "WARMUP"},
        {"ZHANG_R51_ENABLE", "1"}, {"ZHANG_R51_RATIO_ONLY", "1"}, {"ZHANG_R51_RATIO_THRESHOLD", "2"},
        {"ZHANG_R51_EPOCH_ONLY_AR", "1"}, {"ZHANG_R51_AR_START_GPST_SECONDS", "1405216800"}};
    return p;
}
void accept(const ZhangUs20Global134CheckpointPolicy& p)
{
    std::string reason = "uncleared";
    require(validateZhangUs20Global134CheckpointPolicy(p, reason), "valid policy rejected: " + reason);
    require(reason.empty(), "successful validator retained failure reason");
}
void reject(const std::function<void(ZhangUs20Global134CheckpointPolicy&)>& change)
{
    auto p = mainPolicy(); change(p); std::string reason;
    require(!validateZhangUs20Global134CheckpointPolicy(p, reason), "unsafe policy accepted");
    require(!reason.empty(), "rejected policy did not explain reason");
}
ZhangUs20Global134CheckpointPolicy arPolicy()
{
    auto p = mainPolicy(); p.environment["ZHANG_US20_GLOBAL134_STAGE"] = "AR";
    p.checkpointRestorePath = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run/checkpoints/E29-20240717015930-e240-r32f626de9463";
    p.environment["ZHANG_E29_RESTORE_BUNDLE"] = p.checkpointRestorePath;
    return p;
}
}

int main()
{
    try {
        auto p = mainPolicy(); accept(p); accept(arPolicy());
        for (const auto& change : std::vector<std::function<void(ZhangUs20Global134CheckpointPolicy&)>>{
            [](auto& x){x.outputsRoot += "/";}, [](auto& x){x.outputsRoot += "/../outputs";},
            [](auto& x){x.outputsRoot = "/home/rx/GINAN/us20_multignss_server_20261027/outputs";},
            [](auto& x){x.binaryPath += ".other";}, [](auto& x){x.startEpoch = "2024-07-17 02:00:00";},
            [](auto& x){x.endEpoch = "2024-07-17 07:00:00";}, [](auto& x){x.maxEpochs = 841;},
            [](auto& x){x.epochIntervalSeconds = 15;}, [](auto& x){x.arStartGpstSeconds += 30;},
            [](auto& x){x.captureEpochs.pop_back();}, [](auto& x){x.captureEpochs[0] = "2024-07-17 02:00:00";},
            [](auto& x){std::swap(x.captureEpochs[0],x.captureEpochs[1]);},
            [](auto& x){x.productFilename += ".other";}, [](auto& x){x.covarianceFilename += ".other";},
            [](auto& x){x.checkpointOutputDirectory += "/";}, [](auto& x){x.checkpointRuntimeId += "other";},
            [](auto& x){x.receiverCount = 20;}, [](auto& x){x.receiverCount = 133;},
            [](auto& x){x.systems.erase("GAL");}, [](auto& x){x.systems.insert("QZS");},
            [](auto& x){x.fixedCoordinates = false;}, [](auto& x){x.productBasisMode = "LEGACY_PRODUCT_TREE";},
            [](auto& x){x.ratio = 3;}, [](auto& x){x.ratio = std::numeric_limits<double>::quiet_NaN();},
            [](auto& x){x.jointNisEnabled = true;}, [](auto& x){x.independentEpochAr = false;},
            [](auto& x){x.physicalIonosphereEnabled = false;}, [](auto& x){x.physicalTroposphereEnabled = false;},
            [](auto& x){x.stecEstimated = false;}, [](auto& x){x.ztdEstimated = false;},
            [](auto& x){x.regionalJointEnabled = true;}, [](auto& x){x.regionalTropEnabled = true;},
            [](auto& x){x.environment["ZHANG_US20_GLOBAL134_CHECKPOINT"] = "0";},
            [](auto& x){x.environment["ZHANG_US20_GLOBAL134_STAGE"] = "ar";},
            [](auto& x){x.environment.erase("ZHANG_US20_GLOBAL134_STAGE");},
            [](auto& x){x.environment["ZHANG_US20_GLOBAL134_STAGE"] = "AR";},
            [](auto& x){x.environment["ZHANG_R51_ENABLE"] = "0";},
            [](auto& x){x.environment["ZHANG_R51_RATIO_ONLY"] = "0";},
            [](auto& x){x.environment["ZHANG_R51_RATIO_THRESHOLD"] = "3";},
            [](auto& x){x.environment["ZHANG_R51_EPOCH_ONLY_AR"] = "0";},
            [](auto& x){x.environment["ZHANG_R51_AR_START_GPST_SECONDS"] = "1405216830";},
            [](auto& x){x.environment["ZHANG_REGIONAL_JOINT_EXPORT"] = "1";},
            [](auto& x){x.environment["ZHANG_REGIONAL_TROPOSPHERE"] = "GPT2_RESIDUAL_ENH";},
            [](auto& x){x.checkpointRestorePath = "/old/e122";},
            [](auto& x){x.environment["ZHANG_E29_RESTORE_BUNDLE"] = "/old/e122";}
        }) reject(change);
        for (const auto* old : {"ZHANG_GEC134_FLOAT_TO_AR_MIGRATION", "ZHANG_GEC134_AR_CONTINUATION_12_6",
                "ZHANG_GEC134_LOCAL_ONE_EPOCH_12_6", "ZHANG_GEC134_LOCAL_ONE_EPOCH_MIGRATION", "ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT",
                "ZHANG_MULTI_GNSS_AR_CHECKPOINT", "ZHANG_R51_REUSE_FLOAT_20240717", "ZHANG_R51_REUSE_RATIO_20240717"}) {
            reject([&](auto& x){x.environment[old] = "1";});
            reject([&](auto& x){x.environment[old] = "0";});
        }
        p.environment["ZHANG_REGIONAL_JOINT_EXPORT"] = "0";
        p.environment["ZHANG_REGIONAL_TROPOSPHERE"] = "OFF"; accept(p);
        require(zhangUs20Global134WarmupTerminal(p,"2024-07-17 01:59:30",240),"warmup exit missing");
        require(!zhangUs20Global134WarmupTerminal(p,"2024-07-17 01:59:30",239),"wrong warmup epoch accepted");
        require(!zhangUs20Global134WarmupTerminal(arPolicy(),"2024-07-17 01:59:30",240),"AR mistaken for warmup exit");
        auto ar = arPolicy(); std::string reason;
        ar.environment["ZHANG_E29_RESTORE_BUNDLE"] += ".other";
        require(!validateZhangUs20Global134CheckpointPolicy(ar,reason),"conflicting restore path accepted");
        auto probe = mainPolicy();
        const std::string oldRoot = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/run";
        const std::string newRoot = "/home/rx/GINAN/us20_global134_osb5h_20261031_v1/probe/run";
        for (auto* field : {&probe.outputsRoot,&probe.productFilename,&probe.covarianceFilename,&probe.checkpointOutputDirectory})
            field->replace(0,oldRoot.size(),newRoot);
        probe.endEpoch = "2024-07-17 00:00:00"; probe.maxEpochs = 1;
        probe.captureEpochs = {"2024-07-17 00:00:00"}; accept(probe);
        require(zhangUs20Global134WarmupTerminal(probe,"2024-07-17 00:00:00",1),"probe terminal missing");
        probe.checkpointRestorePath = newRoot + "/checkpoints/E29-20240717000000-e1-r32f626de9463";
        probe.environment["ZHANG_E29_RESTORE_BUNDLE"] = probe.checkpointRestorePath; accept(probe);
        ZhangUs20Global134CheckpointIdentity source;
        source.binarySha256 = std::string(64,'a'); source.configSha256 = std::string(64,'b');
        source.inputSha256 = std::string(64,'c'); source.binaryPath = ZHANG_US20_GLOBAL134_BINARY;
        source.experimentMode = "E29_MULTI_GNSS_AR_V1_GPS_GAL_BDS";
        source.runtimeId = "E29-GEC134-FIXEDCOORD-20261023-one_hour";
        source.epoch = "2024-07-17 01:59:30"; source.completedEpoch = 240; source.nextEpoch = 241;
        source.graphSchema = 5; source.sectionCount = 9;
        const auto expected = source;
        require(validateZhangUs20Global134StrictRestoreIdentity(arPolicy(),source,expected,reason),"own e240 rejected");
        for (const auto& mutate : std::vector<std::function<void(ZhangUs20Global134CheckpointIdentity&)>>{
            [](auto& x){x.binarySha256[0]='d';}, [](auto& x){x.configSha256[0]='d';},
            [](auto& x){x.inputSha256[0]='d';}, [](auto& x){x.binarySha256.clear();},
            [](auto& x){x.binarySha256 = std::string(64,'G');},
            [](auto& x){x.binaryPath += ".old";}, [](auto& x){x.experimentMode = "E29_MULTI_GNSS_FLOAT_V2";},
            [](auto& x){x.runtimeId += "other";}, [](auto& x){x.epoch = "2024-07-17 01:00:30";},
            [](auto& x){x.completedEpoch = 122;}, [](auto& x){x.nextEpoch = 240;},
            [](auto& x){x.graphSchema = 4;}, [](auto& x){x.sectionCount = 8;}
        }) {
            auto bad = source; mutate(bad);
            require(!validateZhangUs20Global134StrictRestoreIdentity(arPolicy(),bad,expected,reason),"foreign/incomplete source accepted");
        }
        auto probeSource = source; probeSource.epoch = "2024-07-17 00:00:00";
        probeSource.completedEpoch = 1; probeSource.nextEpoch = 2;
        require(validateZhangUs20Global134StrictRestoreIdentity(probe,probeSource,expected,reason),"own probe e1 rejected");
        require(!validateZhangUs20Global134StrictRestoreIdentity(mainPolicy(),source,expected,reason),"cold stage granted restore");
        const auto initialEnvironment = mainPolicy().environment;
        const auto originalScience = zhangUs20Global134ScienceEnvironment(initialEnvironment);
        // These were not present in the original handwritten canonical list.
        // Both adding/removing a setting and changing its value must alter the
        // science fields; future ZHANG switches receive the same treatment.
        for (const auto* key : {"ZHANG_R51_EVIDENCE_SLIP_QC", "ZHANG_SERVER_REACQ_FIX",
                "ZHANG_R51_HANDOFF_FLOAT_RETENTION_SECONDS", "ZHANG_R51_HANDOFF_EVENT_BUDGET_MS",
                "ZHANG_FUTURE_UNKNOWN_ALGORITHM_SWITCH", "ZHANG_E29_RESTORE_BUNDLE_SUFFIX"}) {
            auto actual = initialEnvironment;
            actual[key] = "1";
            const auto first = zhangUs20Global134ScienceEnvironment(actual);
            require(first != originalScience, "new science setting omitted from identity");
            actual[key] = "2";
            require(zhangUs20Global134ScienceEnvironment(actual) != first,
                "changed science setting omitted from identity");
            actual[key] = "";
            require(zhangUs20Global134ScienceEnvironment(actual) != originalScience,
                "present empty science setting mistaken for absent setting");
            actual[key] = "UNSET";
            require(zhangUs20Global134ScienceEnvironment(actual) != originalScience,
                "literal UNSET science setting mistaken for absent setting");
            actual.erase(key);
            require(zhangUs20Global134ScienceEnvironment(actual) == originalScience,
                "science setting removal did not restore exact fields");
        }
        for (const auto* key : {"OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS",
                "OPENBLAS_CORETYPE", "MKL_NUM_THREADS", "MKL_DYNAMIC", "GOTO_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"}) {
            auto actual = initialEnvironment; actual[key] = "1";
            const auto first = zhangUs20Global134ScienceEnvironment(actual);
            require(first != originalScience, "thread/backend environment omitted from identity");
            actual[key] = "2";
            require(zhangUs20Global134ScienceEnvironment(actual) != first,
                "thread/backend environment value change omitted from identity");
        }
        auto lifecycle = initialEnvironment;
        lifecycle["ZHANG_US20_GLOBAL134_STAGE"] = "AR";
        lifecycle["ZHANG_E29_RESTORE_BUNDLE"] = "/own/e240";
        lifecycle["ZHANG_GEC134_RESTORE_ONLY"] = "1";
        require(zhangUs20Global134ScienceEnvironment(lifecycle) == originalScience,
            "three explicit lifecycle controls changed the science fields");
        lifecycle["PATH"] = "/other/path";
        lifecycle["TEMP"] = "/other/temp";
        require(zhangUs20Global134ScienceEnvironment(lifecycle) == originalScience,
            "unrelated shell environment leaked into the scoped science fields");
        std::cout << "US20_GLOBAL134_CHECKPOINT_POLICY_PASS checks=" << checks
            << " cold_start=1 legacy_migration=0 payload_or_native_epoch_test=0\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "US20_GLOBAL134_CHECKPOINT_POLICY_FAIL " << error.what() << '\n';
        return 1;
    }
}
