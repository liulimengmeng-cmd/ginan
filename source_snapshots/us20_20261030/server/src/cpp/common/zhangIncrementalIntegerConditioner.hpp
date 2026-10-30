#pragma once

#include "common/zhangIntegerConditioner.hpp"
#include "common/zhangParallelProjection.hpp"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

/** One search-stage posterior, with transactional square-root updates.
 *
 * This is not a cache of unrelated P matrices. The mean, covariance and root
 * are owned together, and are changed only by committing a matching ticket.
 * Failed integer proofs leave all three unchanged. Numerical rank/affine
 * checks still use the production conditioner. Rebase on ambiguous rank,
 * accumulated roundoff, and periodically; never add jitter to exact rows.
 */
class ZhangIncrementalIntegerConditioner
{
    VectorXd mean_;
    MatrixXd covariance_;
    ZhangPosteriorEffectiveWorkspace workspace_;
    std::size_t generation_ = 0;
    bool strictAudit_ = false;
    const bool legacy_ = [] {
        const char* value = std::getenv("ZHANG_R51_BLOCK_CONDITIONER");
        if (!value || std::string(value) == "incremental") return false;
        if (std::string(value) == "legacy") return true;
        throw std::invalid_argument("ZHANG_R51_BLOCK_CONDITIONER must be incremental or legacy");
    }();
    std::ofstream replay_;

    template<class T> void write(const T& value)
    { replay_.write(reinterpret_cast<const char*>(&value), sizeof(value)); }
    template<class T> void writeDense(const T& value)
    { replay_.write(reinterpret_cast<const char*>(value.data()), value.size()*sizeof(double)); }

    void build()
    {
        const auto start = std::chrono::steady_clock::now();
        workspace_ = zhangBuildPosteriorEffectiveWorkspace(covariance_);
        ++decompositions;
        if(workspace_.valid && generation_==0)
        {
            // A pre-existing numerical null space may change under the old
            // eigenvalue cutoff when the spectral scale falls. Do not replace
            // that policy with a transported-root rank decision. The fast path
            // starts only from a well-conditioned full-rank quotient.
            strictAudit_ = workspace_.covarianceRank != covariance_.rows();
            if(strictAudit_) fallbackReason="INITIAL_NUMERICAL_NULL_SPACE";
            if(workspace_.squareRoot.cols()>0)
            {
                const VectorXd eigenvalues=workspace_.squareRoot.colwise().squaredNorm();
                if(eigenvalues.minCoeff()<eigenvalues.maxCoeff()*1e-8)
                { strictAudit_=true; fallbackReason="INITIAL_ILL_CONDITIONED_QUOTIENT"; }
            }
        }
        buildSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }

public:
    struct Trial
    {
        ZhangPosteriorEffectiveConditioningResult result;
    private:
        friend class ZhangIncrementalIntegerConditioner;
        const ZhangIncrementalIntegerConditioner* owner = nullptr;
        std::size_t generation = 0;
        ZhangPosteriorEffectiveWorkspace next;
    };
    std::size_t decompositions = 0, prepares = 0, commits = 0, rebases = 0;
    double buildSeconds = 0, conditionSeconds = 0, rootUpdateSeconds = 0;
    std::string fallbackReason="NONE";

    ZhangIncrementalIntegerConditioner(VectorXd mean, MatrixXd covariance,
                                      const std::string& replayPath = {})
        : mean_(std::move(mean)), covariance_(std::move(covariance))
    {
        if (!replayPath.empty())
        {
            // The caller supplies a new, writable diagnostic directory.
            if (std::ifstream(replayPath).good())
                throw std::runtime_error("Refuse to overwrite integer-conditioner replay");
            replay_.open(replayPath, std::ios::binary);
            replay_.exceptions(std::ios::badbit | std::ios::failbit);
            const int n = mean_.size();
            write(n); writeDense(mean_); writeDense(covariance_);
        }
    }
    const VectorXd& mean() const { return mean_; }
    const MatrixXd& covariance() const { return covariance_; }

    Trial prepare(const MatrixXd& rows, const VectorXd& integers)
    {
        ++prepares;
        if (replay_.is_open())
        {
            const char event = 'P'; const int count = rows.rows();
            write(event); write(count); writeDense(rows); writeDense(integers);
        }
        Trial trial; trial.owner = this; trial.generation = generation_;
        const auto legacyTrial = [&]()
        {
            const auto start = std::chrono::steady_clock::now();
            trial.result = zhangConditionPosteriorEffectiveIntegers(mean_, covariance_, rows, integers);
            ++decompositions;
            conditionSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            return trial;
        };
        if (legacy_ || strictAudit_) return legacyTrial();
        // The legacy cutoff depends on the current spectral scale. Re-audit
        // when that scale collapses, including a fully constrained quotient;
        // otherwise old numerical-null directions could silently resurrect.
        const bool scaleCollapsed=workspace_.valid && workspace_.squareRoot.cols()>0 &&
            std::max(1.,workspace_.maximumEigenvalue)>4.*std::max(1.,
                workspace_.squareRoot.colwise().squaredNorm().maxCoeff());
        if (!workspace_.valid || workspace_.covarianceRank==0 || scaleCollapsed ||
            (generation_ > 0 && generation_ % 16 == 0)) build();
        if (!workspace_.valid)
        {
            trial.result.failureReason = workspace_.failureReason;
            return trial;
        }
        if(strictAudit_) return legacyTrial();
        // Root transport is numeric work, never an integer certificate.
        {
            const auto conditionStart = std::chrono::steady_clock::now();
            ZhangPosteriorEffectiveRowBasis effective;
            trial.result = zhangConditionPosteriorEffectiveIntegers(
                mean_, covariance_, rows, integers, 1e-7, &workspace_, &effective);
            conditionSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-conditionStart).count();
            if (!trial.result.valid) return trial;
            const auto rootStart = std::chrono::steady_clock::now();
            trial.next = workspace_;
            if (!trial.result.conditioned) return trial;
            const ZhangParallelRowProjection projection(rows);
            const MatrixXd whitened = projection.multiply(workspace_.squareRoot);
            MatrixXd q = MatrixXd::Zero(whitened.cols(), effective.posteriorRank);
            int column = 0;
            bool ambiguous = false;
            for (int row = 0; row < rows.rows(); ++row)
            {
                const double variance = effective.conditionalVariances[row];
                ambiguous |= variance > effective.varianceTolerance / 256 &&
                             variance < effective.varianceTolerance * 256;
            }
            for (const int row : effective.retainedRows)
            {
                VectorXd direction = whitened.row(row).transpose();
                for (int pass = 0; pass < 2; ++pass)
                    for (int j = 0; j < column; ++j)
                        direction -= q.col(j).dot(direction) * q.col(j);
                const double norm = direction.norm();
                if (!(norm > 0) || !std::isfinite(norm)) { ambiguous = true; break; }
                q.col(column++) = direction / norm;
            }
            // Compare in the small constraint space before trusting transport.
            const MatrixXd cross = projection.multiply(covariance_);
            const MatrixXd expected = projection.covariance(cross);
            const MatrixXd actual = whitened * whitened.transpose();
            const double scale = std::max(1.0, expected.cwiseAbs().maxCoeff());
            const double errorLimit = scale * std::numeric_limits<double>::epsilon() *
                covariance_.rows() * 512;
            ambiguous |= column != effective.posteriorRank ||
                (actual-expected).cwiseAbs().maxCoeff() > errorLimit;
            // Ill-conditioned small solves can amplify subtraction error far
            // beyond a transported root's roundoff. Keep their old full-PSD
            // audit rather than treating mathematical PSD as a numeric proof.
            Eigen::LDLT<MatrixXd> conditionCheck(expected);
            const bool illConditioned=conditionCheck.info()!=Eigen::Success ||
                !conditionCheck.isPositive() || conditionCheck.rcond()<1e-10;
            if(illConditioned)
            {
                strictAudit_=true;
                fallbackReason="ILL_CONDITIONED_CONSTRAINT";
                return legacyTrial();
            }
            if (ambiguous)
            {
                ++rebases;
                strictAudit_=true;
                fallbackReason="RANK_OR_ROOT_GRAM_AMBIGUITY";
                return legacyTrial();
            }
            const MatrixXd removed = workspace_.squareRoot * q;
            trial.next.squareRoot.noalias() -= removed * q.transpose();
            trial.next.covarianceRank = std::max(0, workspace_.covarianceRank-effective.posteriorRank);
            if (!trial.next.squareRoot.allFinite())
            { strictAudit_=true; fallbackReason="ROOT_UPDATE_NONFINITE"; ++rebases; return legacyTrial(); }
            rootUpdateSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-rootStart).count();
            return trial;
        }
    }

    bool commit(Trial&& trial)
    {
        if (trial.owner != this || trial.generation != generation_ || !trial.result.valid) return false;
        if (replay_.is_open()) { const char event = 'C'; write(event); replay_.flush(); }
        mean_ = std::move(trial.result.mean);
        covariance_ = std::move(trial.result.covariance);
        if (!legacy_ && !strictAudit_) workspace_ = std::move(trial.next);
        trial.owner = nullptr;
        ++generation_; ++commits;
        return true;
    }
};
