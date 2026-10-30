#pragma once

#include "common/zhangP0ResourceProbe.hpp"
#include <iostream>
#include <mutex>

// Point diagnostics, deliberately not another elapsed-time scope. No cache,
// estimator state, allocator policy or acceptance decision is changed here.
struct ZhangFloatInputInventory
{
    bool measured = false;
    long long receivers = -1;
    long long receiverObservations = -1;
    long long streams = -1;
    long long queuedEpochs = -1;
    long long queuedObservations = -1;
    long long states = -1;
    long long covarianceRows = -1;
    long long covarianceColumns = -1;
    long long preciseOrbitEpochs = -1;
    long long preciseClockEpochs = -1;
    long long fileStreams = -1;
    long long filePositionTotal = -1;
    long long filePositionMinimum = -1;
    long long filePositionMaximum = -1;
    long long selectedFilePosition = -1;
    long long configReceiverOptions = -1;
    long long configSatelliteOptions = -1;
    long long configCustomAliases = -1;
    long long configFoundOptions = -1;
    long long configDocumentationEntries = -1;
    long long configuredYamlFiles = -1;
    long long configYamlDefaults = -1;
    long long configAvailableOptions = -1;
};

inline bool zhangFloatInputMemoryEnabled()
{
    static const bool enabled = []
    {
        const char* setting = std::getenv("ZHANG_FLOAT_PROFILE");
        return setting && setting[0] == '1' && setting[1] == '\0';
    }();
    return enabled;
}

inline std::string zhangFloatInputMemoryLine(
    const char* phase,
    const char* boundary,
    const std::string& epoch,
    const std::string& receiver,
    long long monotonicNs,
    const ZhangP0Resources& resources,
    const ZhangFloatInputInventory& inventory)
{
    std::ostringstream line;
    line << std::setprecision(17)
         << "ZHANG_FLOAT_INPUT_MEMORY phase=" << phase
         << " boundary=" << boundary << " epoch=" << std::quoted(epoch)
         << " receiver=" << std::quoted(receiver)
         << " monotonic_ns=" << monotonicNs
         << " rss_kib=" << resources.rssKiB << " hwm_kib=" << resources.hwmKiB
         << " anon_kib=" << resources.anonKiB << " swap_kib=" << resources.swapKiB
         << " vmsize_kib=" << resources.vmKiB
         << " allocator_used_bytes=" << resources.allocatorUsedBytes
         << " allocator_mmap_bytes=" << resources.allocatorMmapBytes
         << " minor_faults=" << resources.minorFaults
         << " major_faults=" << resources.majorFaults
         << " read_bytes=" << resources.readBytes << " write_bytes=" << resources.writeBytes
         << " user_s=" << resources.userSeconds << " system_s=" << resources.systemSeconds
         << " inventory_measured=" << inventory.measured << " unknown_inventory_value=-1"
         << " receivers=" << inventory.receivers
         << " receiver_observations=" << inventory.receiverObservations
         << " streams=" << inventory.streams
         << " queued_epochs=" << inventory.queuedEpochs
         << " queued_observations=" << inventory.queuedObservations
         << " states=" << inventory.states
         << " covariance_rows=" << inventory.covarianceRows
         << " covariance_columns=" << inventory.covarianceColumns
         << " precise_orbit_epochs=" << inventory.preciseOrbitEpochs
         << " precise_clock_epochs=" << inventory.preciseClockEpochs
         << " file_streams=" << inventory.fileStreams
         << " file_position_total=" << inventory.filePositionTotal
         << " file_position_min=" << inventory.filePositionMinimum
         << " file_position_max=" << inventory.filePositionMaximum
         << " selected_file_position=" << inventory.selectedFilePosition
         << " config_receiver_options=" << inventory.configReceiverOptions
         << " config_satellite_options=" << inventory.configSatelliteOptions
         << " config_custom_aliases=" << inventory.configCustomAliases
         << " config_found_options=" << inventory.configFoundOptions
         << " config_documentation_entries=" << inventory.configDocumentationEntries
         << " configured_yaml_files=" << inventory.configuredYamlFiles
         << " config_yaml_defaults=" << inventory.configYamlDefaults
         << " config_available_options=" << inventory.configAvailableOptions
         << " point_sample_not_phase_peak=1 no_elapsed_timer=1";
    return line.str();
}

inline void zhangFloatInputMemoryBoundary(
    const char* phase,
    const char* boundary,
    const std::string& epoch,
    const std::string& receiver,
    const ZhangFloatInputInventory& inventory) noexcept
{
    if (!zhangFloatInputMemoryEnabled()) return;
    try
    {
        const auto resources = zhangP0ReadResources();
        const auto monotonicNs = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
        const auto message = zhangFloatInputMemoryLine(
            phase,boundary,epoch,receiver,monotonicNs,resources,inventory);
        static std::mutex outputMutex;
        std::lock_guard<std::mutex> lock(outputMutex);
        // Flush independently of pooled TRACE buffering so a resource stop
        // still leaves the last completed input-lifecycle boundary in the log.
        std::clog << message << std::endl;
    }
    catch (...)
    {
        // Diagnostics must not change the numerical or publication path.
    }
}
