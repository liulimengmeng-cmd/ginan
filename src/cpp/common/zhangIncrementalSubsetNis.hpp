#pragma once
#include "common/eigenIncluder.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

struct ZhangSubsetNisTrial
{
    bool valid = false, fast = false;
    int rank = 0;
    double nis = std::numeric_limits<double>::quiet_NaN();
    double nullInnovation = 0;
    MatrixXd inverse;
};

// Local to ONE immutable (S,v) and ONE ordered subset admission. The SPD
// Schur path is allowed only far from the legacy eigen cutoff. Rank-deficient
// trials and near-threshold decisions use the original eigensolver, including
// its deterministic residual gate. No jitter or inverse of a singular matrix.
class ZhangIncrementalSubsetNis
{
    const MatrixXd& covariance;
    const VectorXd& innovation;
    std::vector<int> selected;
    MatrixXd inverse;
    int fastSinceReference = 0;
public:
    int fastTrials = 0, referenceTrials = 0;
    ZhangIncrementalSubsetNis(const MatrixXd& s, const VectorXd& v) : covariance(s), innovation(v) {}
    ZhangSubsetNisTrial trial(int candidate, bool forceReference = false)
    {
        ZhangSubsetNisTrial r;
        auto indices = selected; indices.push_back(candidate);
        const int n = indices.size(), k = selected.size();
        const VectorXd v = innovation(indices);
        const MatrixXd s = covariance(indices,indices);
        // The empty base has a 0x0 inverse (lda=0).  Eigen's BLAS path
        // validates lda before its zero-size quick return.  Seed with the
        // unchanged scalar reference test, then use Schur updates for k>0.
        if (k > 0 && !forceReference && fastSinceReference < 16 && inverse.rows()==k)
        {
            VectorXd b = s.topRightCorner(k,1);
            VectorXd u = inverse*b;
            const double schur = s(k,k)-b.dot(u);
            if (std::isfinite(schur) && schur>0)
            {
                r.inverse = MatrixXd::Zero(n,n);
                r.inverse.topLeftCorner(k,k) = inverse + u*u.transpose()/schur;
                r.inverse.topRightCorner(k,1) = -u/schur;
                r.inverse.bottomLeftCorner(1,k) = -u.transpose()/schur;
                r.inverse(k,k) = 1/schur;
                const double lower = 1/r.inverse.norm();
                const double upper = s.trace(); // SPD eigenvalue upper bound
                const double conditionBound = upper/lower;
                // Strict separation from 1e-12*lambda_max and from numerical
                // inverse-update uncertainty. This is deliberately conservative.
                if (r.inverse.allFinite() && lower > std::max(1e-14,1e-12*upper)*128 &&
                    conditionBound*n*std::numeric_limits<double>::epsilon() < 1e-9)
                {
                    r.nis = v.dot(r.inverse*v);
                    r.rank=n; r.valid=std::isfinite(r.nis) && r.nis>=0; r.fast=r.valid;
                    if (r.valid) { ++fastTrials; return r; }
                }
            }
        }
        ++referenceTrials;
        r = {};
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(s);
        if (eigen.info()!=Eigen::Success || !eigen.eigenvalues().allFinite()) return r;
        const double tolerance=std::max(1e-14,1e-12*eigen.eigenvalues().maxCoeff());
        const VectorXd coordinates=eigen.eigenvectors().transpose()*v;
        VectorXd inv=VectorXd::Zero(n);
        for(int i=0;i<n;++i)
            if(eigen.eigenvalues()(i)>tolerance) { inv(i)=1/eigen.eigenvalues()(i); ++r.rank; }
            else r.nullInnovation=std::max(r.nullInnovation,std::abs(coordinates(i)));
        r.nis=coordinates.dot(inv.asDiagonal()*coordinates);
        r.valid=r.rank>0 && r.nullInnovation<=1e-7 && std::isfinite(r.nis);
        if(r.rank==n) r.inverse=eigen.eigenvectors()*inv.asDiagonal()*eigen.eigenvectors().transpose();
        return r;
    }
    void commit(int candidate, ZhangSubsetNisTrial&& r)
    {
        selected.push_back(candidate);
        inverse=std::move(r.inverse);
        fastSinceReference = r.fast ? fastSinceReference+1 : 0;
    }
};
