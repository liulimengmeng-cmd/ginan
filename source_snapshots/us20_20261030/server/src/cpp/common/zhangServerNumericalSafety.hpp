#pragma once
#include "common/eigenIncluder.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

// Explicitly scoped server rollout; do not silently change user PPP policies.
inline bool zhangServerNumericalSafetyEnabled()
{
    const char* value=std::getenv("ZHANG_SERVER_REACQ_FIX");
    if (!value || std::string(value)=="0") return false;
    if (std::string(value)=="1") return true;
    throw std::invalid_argument("ZHANG_SERVER_REACQ_FIX must be 0 or 1");
}

struct ZhangCovarianceHealth
{
    bool valid=false;
    double minimumEigenvalue=0,scale=1;
    std::string reason;
};
inline ZhangCovarianceHealth zhangCheckFullCovariance(const MatrixXd& p)
{
    ZhangCovarianceHealth out;
    if (!p.rows() || p.rows()!=p.cols() || !p.allFinite())
    { out.reason="NONFINITE_OR_DIMENSION"; return out; }
    out.scale=std::max(1.,p.cwiseAbs().maxCoeff());
    if ((p-p.transpose()).cwiseAbs().maxCoeff()>1e-11*out.scale)
    { out.reason="ASYMMETRIC"; return out; }
    // A positive diagonal is necessary but is never a substitute for full PSD.
    if ((p.diagonal().array()<0).any())
    { out.reason="NEGATIVE_VARIANCE"; return out; }
    Eigen::SelfAdjointEigenSolver<MatrixXd> eig((p+p.transpose())*.5,Eigen::EigenvaluesOnly);
    if (eig.info()!=Eigen::Success || !eig.eigenvalues().allFinite())
    { out.reason="EIGENSOLVER_FAILED"; return out; }
    out.minimumEigenvalue=eig.eigenvalues().minCoeff();
    out.scale=std::max(1.,eig.eigenvalues().cwiseAbs().maxCoeff());
    // Stricter than the downstream AR root's 1e-9 relative rejection gate.
    // No eigenvalue clipping and no mutation of the submitted covariance.
    out.valid=out.minimumEigenvalue>=-1e-10*out.scale;
    out.reason=out.valid?"VALID":"NOT_PSD";
    return out;
}

struct ZhangIndependentStateReset { int index; double mean,variance; };
inline bool zhangBuildIndependentStateReset(int size,int one,
    const std::vector<ZhangIndependentStateReset>& resets,MatrixXd& f,MatrixXd& q)
{
    if (size<1 || one<0 || one>=size || resets.empty()) return false;
    std::vector<bool> seen(size,false);
    f=MatrixXd::Identity(size,size); q=MatrixXd::Zero(size,size);
    for (const auto& reset:resets)
    {
        if (reset.index<0 || reset.index>=size || reset.index==one || seen[reset.index] ||
            !std::isfinite(reset.mean) || !std::isfinite(reset.variance) || reset.variance<=0)
            return false;
        seen[reset.index]=true;
        f.row(reset.index).setZero(); f(reset.index,one)=reset.mean;
        q(reset.index,reset.index)=reset.variance;
    }
    return true;
}
