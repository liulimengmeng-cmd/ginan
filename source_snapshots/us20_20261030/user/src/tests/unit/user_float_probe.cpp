#include "common/algebra.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

struct Network;
struct ReceiverMap;

// Diagnostic-only link interposition; the production function is invoked
// unchanged. The authoritative moments are read BEFORE any epoch AR output.
extern "C" void realUserOutput(
    Trace&, GTime, Network&, ReceiverMap&, KFState&, bool, bool, bool)
    asm("__real__Z32perEpochPostProcessingAndOutputsRSo5GTimeR7NetworkR11ReceiverMapR7KFStatebbb");
extern "C" void probeUserOutput(
    Trace& trace, GTime time, Network& net, ReceiverMap& receivers,
    KFState& state, bool empty, bool rts, bool firstRts)
    asm("__wrap__Z32perEpochPostProcessingAndOutputsRSo5GTimeR7NetworkR11ReceiverMapR7KFStatebbb");

extern "C" void probeUserOutput(
    Trace& trace, GTime time, Network& net, ReceiverMap& receivers,
    KFState& state, bool empty, bool rts, bool firstRts)
{
    const char* directory = std::getenv("ZHANG_USER_FLOAT_PROBE_DIR");
    if (directory && !rts)
    {
        const std::filesystem::path path = std::filesystem::path(directory) /
            (std::to_string(static_cast<long long>(state.time.bigTime)) + ".moments");
        if (std::filesystem::exists(path)) throw std::runtime_error("FLOAT_PROBE_WILL_NOT_OVERWRITE");
        std::ofstream file(path, std::ios::binary);
        const std::int64_t dimension = state.x.size();
        file.write(reinterpret_cast<const char*>(&dimension), sizeof(dimension));
        file.write(reinterpret_cast<const char*>(state.x.data()), state.x.size() * sizeof(double));
        file.write(reinterpret_cast<const char*>(state.P.data()), state.P.size() * sizeof(double));
        if (!file) throw std::runtime_error("FLOAT_PROBE_WRITE_FAILED");
        const auto precision = trace.precision();
        trace << std::setprecision(17);
        for (const auto& [key, rate] : state.procNoiseMap)
            if (key.type == KF::TROP)
                trace << "\nZHANG_USER_FLOAT_PROBE_TROP time=" << state.time.to_string(0)
                      << " receiver=" << key.str << " variance_rate_m2_per_second=" << rate;
        trace.precision(precision);
    }
    realUserOutput(trace, time, net, receivers, state, empty, rts, firstRts);
}
