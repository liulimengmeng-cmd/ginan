#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "common/zhangIntegerConditioner.hpp"

struct ZhangProductBasePosteriorIdentity
{
    std::string owner, epoch, root;
    std::uint64_t factorSequence = 0;
    bool operator==(const ZhangProductBasePosteriorIdentity& other) const
    {
        return owner == other.owner && epoch == other.epoch && root == other.root &&
            factorSequence == other.factorSequence;
    }
};

// KFKey::operator== intentionally omits comment. It is part of the serialized
// coordinate descriptor, so the product-base factor binds it explicitly.
// rec_ptr/estimatedTime are runtime annotations, not coordinate fields.
struct ZhangProductBaseColumnDescriptor
{
    int column=0, type=0, satelliteSystem=0, satellitePrn=0, observable=0;
    std::string receiver, comment;
    bool operator==(const ZhangProductBaseColumnDescriptor& other) const
    {
        return column==other.column && type==other.type &&
            satelliteSystem==other.satelliteSystem && satellitePrn==other.satellitePrn &&
            observable==other.observable && receiver==other.receiver && comment==other.comment;
    }
};

/** One layered invocation's original FLOAT decomposition, not its conditioned P.
 *
 * The production caller owns a shared_ptr<const GinAR_mtx> with no mutable
 * alias. An aliasing shared_ptr retains that existing P without making another
 * covariance copy. Identity alone never authorizes reuse: complete ordered keys
 * and every covariance coefficient must also match. Expanded evidence/nuisance
 * columns use the original uncached conditioner. No process-global state.
 */
template<class Key> class ZhangProductBasePosteriorWorkspace
{
    const ZhangProductBasePosteriorIdentity identity_;
    const std::vector<Key> keys_;
    const std::shared_ptr<const MatrixXd> covariance_;
    ZhangPosteriorEffectiveWorkspace workspace_;
    bool attempted_ = false;
public:
    std::uint64_t requests = 0, hits = 0, builds = 0, fallbacks = 0;
    std::string lastReason = "NOT_REQUESTED";

    ZhangProductBasePosteriorWorkspace(ZhangProductBasePosteriorIdentity identity,
        std::vector<Key> keys, std::shared_ptr<const MatrixXd> covariance)
        : identity_(std::move(identity)), keys_(std::move(keys)),
          covariance_(std::move(covariance)) {}
    ZhangProductBasePosteriorWorkspace(const ZhangProductBasePosteriorWorkspace&) = delete;
    ZhangProductBasePosteriorWorkspace& operator=(const ZhangProductBasePosteriorWorkspace&) = delete;

    std::size_t sourceCovarianceBytes() const
    { return covariance_ ? std::size_t(covariance_->size()) * sizeof(double) : 0; }
    std::size_t effectiveRootBytes() const
    { return std::size_t(workspace_.squareRoot.size()) * sizeof(double); }
    const std::string& frozenSourceRoot() const { return identity_.root; }

    const ZhangPosteriorEffectiveWorkspace* acquire(
        const ZhangProductBasePosteriorIdentity& identity, const std::vector<Key>& keys,
        const MatrixXd& covariance, bool requirePinnedImmutableSource = false)
    {
        ZhangPhaseTimer timer("PRODUCT_BASE_WORKSPACE_IDENTITY_AND_ACQUIRE");
        ++requests;
        auto fallback = [&](const char* reason) -> const ZhangPosteriorEffectiveWorkspace*
        { ++fallbacks; lastReason = reason; return nullptr; };
        if (identity_.owner.empty() || identity_.epoch.empty() || identity_.root.empty() ||
            !(identity == identity_)) return fallback("SOURCE_OWNER_EPOCH_SEQUENCE_OR_ROOT_CHANGED");
        // In production this proves the reused root token belongs to the same
        // still-owned immutable snapshot, before all value checks below. The
        // object pointer alone never grants reuse. A fresh/expanded P falls back
        // without building a huge scientific-format fingerprint string.
        if (requirePinnedImmutableSource && covariance_.get()!=&covariance)
            return fallback("SOURCE_IS_NOT_THE_PINNED_IMMUTABLE_FLOAT_MATRIX");
        if (keys != keys_) return fallback("SOURCE_ORDER_CHANGED_OR_EVIDENCE_EXPANDED");
        if (!covariance_ || covariance.rows() <= 0 || covariance.rows() != covariance.cols() ||
            keys.size() != std::size_t(covariance.rows()) ||
            covariance.rows() != covariance_->rows() || covariance.cols() != covariance_->cols() ||
            !(covariance.array() == covariance_->array()).all())
            return fallback("SOURCE_COVARIANCE_NOT_EXACTLY_EQUAL");
        bool builtNow = false;
        if (!attempted_)
        {
            attempted_ = true; builtNow = true; ++builds;
            workspace_ = zhangBuildPosteriorEffectiveWorkspace(*covariance_);
        }
        if (!workspace_.valid) return fallback("SOURCE_EFFECTIVE_FACTOR_INVALID_USE_ORIGINAL_CHECKS");
        if (!builtNow) ++hits;
        lastReason = "EXACT_IMMUTABLE_FLOAT_FACTOR";
        return &workspace_;
    }
};
