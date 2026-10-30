#pragma once
#include <atomic>
#include <csignal>
#include <cstdlib>

// Run control is separate from the authenticated numerical configuration.
// A signal must not mutate max_epochs and invalidate checkpoint provenance.
class ZhangEpochStopRequest
{
    std::atomic<std::sig_atomic_t> requested_{0};
    static_assert(decltype(requested_)::is_always_lock_free,
                  "Signal stop requests require a lock-free atomic");
public:
    bool request() noexcept
    { return requested_.exchange(1,std::memory_order_relaxed)==0; }
    bool pending() const noexcept
    { return requested_.load(std::memory_order_relaxed)!=0; }
    bool interruptedCaptureSchedule(bool captureMode,bool hasPending,bool fatal) const noexcept
    { return pending() && captureMode && hasPending && !fatal; }
    int exitCode(bool fatal) const noexcept
    { return fatal?EXIT_FAILURE:pending()?128+SIGINT:EXIT_SUCCESS; }
};
