#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <utility>

struct ZhangP0GraphStamp
{
	std::string runtimeIdentity;
	std::string epochIdentity;
	int system = 0;
	int event = 0;
	int representation = 0;
	int productDatum = 0;
	int floatGauge = 0;
	int integerComponent = 0;
	int runtimeEpoch = -1;
	bool available = false;

	auto fields() const
	{
		return std::tie(
			runtimeIdentity,
			epochIdentity,
			system,
			event,
			representation,
			productDatum,
			floatGauge,
			integerComponent,
			runtimeEpoch,
			available);
	}

	bool operator<(const ZhangP0GraphStamp& other) const
	{
		return fields() < other.fields();
	}

	bool operator==(const ZhangP0GraphStamp& other) const
	{
		return fields() == other.fields();
	}
};

// This cache is deliberately invocation-local, not static or thread-local.
// std::map preserves references to owned values. A failed loader cannot leave
// behind a partially initialized cache hit.
template<typename Key, typename Value>
class ZhangP0InvocationCache
{
	std::map<Key, Value> entries_;

public:
	std::size_t requests = 0;
	std::size_t builds = 0;

	ZhangP0InvocationCache() = default;
	ZhangP0InvocationCache(const ZhangP0InvocationCache&) = delete;
	ZhangP0InvocationCache& operator=(const ZhangP0InvocationCache&) = delete;

	template<typename Loader>
	const Value& get(const Key& key, Loader&& loader)
	{
		++requests;
		auto found = entries_.find(key);
		if (found != entries_.end())
		{
			return found->second;
		}
		Value value = std::forward<Loader>(loader)();
		auto inserted = entries_.emplace(key, std::move(value));
		++builds;
		return inserted.first->second;
	}

	std::size_t size() const noexcept
	{
		return entries_.size();
	}
};
