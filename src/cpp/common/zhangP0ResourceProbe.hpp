#pragma once

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef __linux__
#include <sys/resource.h>
#include <unistd.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#endif

// Opt-in diagnostics only: no gates, state changes, allocator trimming or
// signal handlers. RSS probes are boundary samples, not sub-stage peak claims.
struct ZhangP0Resources
{
	long long rssKiB = -1;
	long long hwmKiB = -1;
	long long vmKiB = -1;
	long long swapKiB = -1;
	long long anonKiB = -1;
	long long majorFaults = -1;
	long long minorFaults = -1;
	long long readBytes = -1;
	long long writeBytes = -1;
	long long allocatorUsedBytes = -1;
	long long allocatorMmapBytes = -1;
	double userSeconds = -1;
	double systemSeconds = -1;
};

inline bool zhangP0ResourceEnabled()
{
	static const bool enabled = []
	{
		const char* value = std::getenv("ZHANG_P0_RESOURCE_PROBE");
		return value && std::string(value) == "1";
	}();
	return enabled;
}

inline ZhangP0Resources zhangP0ReadResources()
{
	ZhangP0Resources resources;
#ifdef __linux__
	std::ifstream status("/proc/self/status");
	std::string line;
	while (std::getline(status, line))
	{
		std::istringstream input(line);
		std::string key;
		long long value;
		if (!(input >> key >> value))
		{
			continue;
		}
		if      (key == "VmRSS:")   resources.rssKiB  = value;
		else if (key == "VmHWM:")   resources.hwmKiB  = value;
		else if (key == "VmSize:")  resources.vmKiB   = value;
		else if (key == "VmSwap:")  resources.swapKiB = value;
		else if (key == "RssAnon:") resources.anonKiB = value;
	}

	std::ifstream io("/proc/self/io");
	while (std::getline(io, line))
	{
		std::istringstream input(line);
		std::string key;
		long long value;
		if (!(input >> key >> value))
		{
			continue;
		}
		if (key == "read_bytes:")  resources.readBytes  = value;
		if (key == "write_bytes:") resources.writeBytes = value;
	}

	struct rusage usage {};
	if (getrusage(RUSAGE_SELF, &usage) == 0)
	{
		resources.majorFaults  = usage.ru_majflt;
		resources.minorFaults  = usage.ru_minflt;
		resources.userSeconds  = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1e-6;
		resources.systemSeconds = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1e-6;
	}
#if defined(__GLIBC__) && defined(__GLIBC_PREREQ)
#if __GLIBC_PREREQ(2, 33)
	const auto allocator = mallinfo2();
	resources.allocatorUsedBytes = allocator.uordblks;
	resources.allocatorMmapBytes = allocator.hblkhd;
#endif
#endif
#endif
	return resources;
}

class ZhangP0ResourceScope
{
	using Clock = std::chrono::steady_clock;
	inline static thread_local std::vector<ZhangP0ResourceScope*> stack_;

	std::ostream& out_;
	std::string name_;
	std::string epoch_;
	std::uint64_t bytes_;
	bool enabled_;
	Clock::time_point begin_;
	double children_ = 0;

	void record(const char* boundary, double inclusive = 0, double exclusive = 0)
	{
		const auto resources = zhangP0ReadResources();
		std::ostringstream line;
		line << std::setprecision(17)
			 << "\nZHANG_P0_RESOURCE phase=" << name_
			 << " epoch=" << std::quoted(epoch_)
			 << " boundary=" << boundary
			 << " monotonic_ns="
			 << std::chrono::duration_cast<std::chrono::nanoseconds>(
					Clock::now().time_since_epoch()).count()
			 << " dense_payload_bytes=" << bytes_
			 << " rss_kib=" << resources.rssKiB
			 << " hwm_kib=" << resources.hwmKiB
			 << " vmsize_kib=" << resources.vmKiB
			 << " swap_kib=" << resources.swapKiB
			 << " anon_kib=" << resources.anonKiB
			 << " major_faults=" << resources.majorFaults
			 << " minor_faults=" << resources.minorFaults
			 << " read_bytes=" << resources.readBytes
			 << " write_bytes=" << resources.writeBytes
			 << " allocator_used_bytes=" << resources.allocatorUsedBytes
			 << " allocator_mmap_bytes=" << resources.allocatorMmapBytes
			 << " user_s=" << resources.userSeconds
			 << " system_s=" << resources.systemSeconds
			 << " inclusive_s=" << inclusive
			 << " exclusive_s=" << exclusive
			 << " boundary_sample_not_phase_peak=1";
		out_ << line.str();
	}

public:
	ZhangP0ResourceScope(
		std::ostream& out,
		std::string name,
		std::string epoch = {},
		std::uint64_t bytes = 0)
		: out_(out)
		, name_(std::move(name))
		, epoch_(std::move(epoch))
		, bytes_(bytes)
		, enabled_(zhangP0ResourceEnabled())
	{
		if (!enabled_)
		{
			return;
		}
		begin_ = Clock::now();
		try
		{
			stack_.push_back(this);
			record("BEGIN");
		}
		catch (...)
		{
			if (!stack_.empty() && stack_.back() == this)
			{
				stack_.pop_back();
			}
			enabled_ = false;
		}
	}

	~ZhangP0ResourceScope() noexcept
	{
		if (!enabled_)
		{
			return;
		}
		try
		{
			const double elapsed = std::chrono::duration<double>(Clock::now() - begin_).count();
			record("END", elapsed, elapsed - children_);
			if (!stack_.empty() && stack_.back() == this)
			{
				stack_.pop_back();
				if (!stack_.empty())
				{
					stack_.back()->children_ += elapsed;
				}
			}
		}
		catch (...)
		{
			// A diagnostic sink failure must not change AR acceptance.
			if (!stack_.empty() && stack_.back() == this)
			{
				stack_.pop_back();
			}
		}
	}

	ZhangP0ResourceScope(const ZhangP0ResourceScope&) = delete;
	ZhangP0ResourceScope& operator=(const ZhangP0ResourceScope&) = delete;
};
