#pragma once

#include "common/zhangIntegerAudit.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>

// Exact append-only affine row-HNF for one block-search invocation. This is
// not a cross-epoch cache and does not replace integer-affine compatibility.
class ZhangP1JointHnf
{
	static bool rectangular(const ZhangExactMatrix& rows, std::size_t dimension)
	{
		for (const auto& row : rows)
		{
			if (row.size() != dimension)
			{
				return false;
			}
		}
		return true;
	}

	ZhangExactRowHnf committed_;
	std::size_t dimension_;
	std::uint64_t generation_ = 0;

public:
	struct Trial
	{
		const ZhangP1JointHnf* owner = nullptr;
		std::uint64_t generation = 0;
		ZhangExactRowHnf hnf;
	};

	ZhangP1JointHnf(
		const ZhangExactMatrix& rows,
		const ZhangExactVector& values,
		std::size_t dimension)
		: dimension_(dimension)
	{
		if (rows.size() != values.size() || !rectangular(rows, dimension))
		{
			throw std::invalid_argument("P1_HNF_INITIAL_DIMENSION");
		}
		committed_ = zhangExactRowHermiteNormalForm(rows, values);
	}

	ZhangP1JointHnf(const ZhangP1JointHnf&) = delete;
	ZhangP1JointHnf& operator=(const ZhangP1JointHnf&) = delete;

	Trial assess(
		const ZhangExactMatrix& rows,
		const ZhangExactVector& values) const
	{
		if (rows.size() != values.size() || !rectangular(rows, dimension_))
		{
			throw std::invalid_argument("P1_HNF_CANDIDATE_DIMENSION");
		}
		Trial trial;
		trial.owner = this;
		trial.generation = generation_;
		if (!committed_.consistent)
		{
			trial.hnf = committed_;
			return trial;
		}
		auto all = committed_.basis;
		auto rhs = committed_.values;
		all.insert(all.end(), rows.begin(), rows.end());
		rhs.insert(rhs.end(), values.begin(), values.end());
		trial.hnf = zhangExactRowHermiteNormalForm(all, rhs);
		return trial;
	}

	bool commit(Trial&& trial)
	{
		if (trial.owner != this
			|| trial.generation != generation_
			|| !trial.hnf.consistent)
		{
			return false;
		}
		committed_ = std::move(trial.hnf);
		++generation_;
		trial.owner = nullptr;
		return true;
	}

	const ZhangExactRowHnf& current() const noexcept
	{
		return committed_;
	}
};
