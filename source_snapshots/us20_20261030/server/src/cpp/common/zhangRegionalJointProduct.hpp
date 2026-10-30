#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Target-space projection only. No factorisation/conditioning of the estimator,
// no NIS, and no invented independence across products or epochs.
using ZhangRegionalRow = std::vector<std::pair<int,double>>;
struct ZhangRegionalParameter {
    std::string parameter,receiver,satellite,observable;
    double sourceUpdated=0;
};
struct ZhangRegionalJointSnapshot {
    bool valid=false,pppCommitted=false,pppQcValid=false;
    std::string solution,root;
    double epoch=0,minimumEigenvalue=0;
    std::vector<ZhangRegionalParameter> parameters;
    std::vector<ZhangRegionalRow> jacobian;
    Eigen::VectorXd mean;
    Eigen::MatrixXd covariance;
};
inline bool zhangRegionalJointEnabled() {
    const char* p=std::getenv("ZHANG_REGIONAL_JOINT_EXPORT");
    return p && std::string(p)=="1";
}
inline bool zhangRegionalTropEnabled() {
    const char* value=std::getenv("ZHANG_REGIONAL_TROPOSPHERE");
    const std::string mode=value?value:"OFF";
    if (mode=="OFF") return false;
    if (mode!="GPT2_RESIDUAL_ENH" || !zhangRegionalJointEnabled())
        throw std::runtime_error("REGIONAL_TROP_MODE_OR_JOINT_EXPORT_INVALID");
    return true;
}
inline void zhangRegionalValidate(const Eigen::VectorXd& mean,const Eigen::MatrixXd& q) {
    if (q.rows()!=mean.size() || q.cols()!=mean.size() || !mean.allFinite() || !q.allFinite())
        throw std::runtime_error("REGIONAL_JOINT_NONFINITE_OR_DIMENSION");
    const double scale=std::max(1.0,q.cwiseAbs().maxCoeff());
    if ((q-q.transpose()).cwiseAbs().maxCoeff()>1e-10*scale)
        throw std::runtime_error("REGIONAL_JOINT_ASYMMETRIC");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(q,Eigen::EigenvaluesOnly);
    if (eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff() < -1e-8*scale)
        throw std::runtime_error("REGIONAL_JOINT_INDEFINITE");
}
inline std::pair<Eigen::VectorXd,Eigen::MatrixXd> zhangRegionalProject(
    const Eigen::VectorXd& x,const Eigen::MatrixXd& p,const std::vector<ZhangRegionalRow>& rows) {
    if (p.rows()!=x.size() || p.cols()!=x.size()) throw std::runtime_error("REGIONAL_SOURCE_DIMENSION");
    Eigen::VectorXd mean=Eigen::VectorXd::Zero(rows.size());
    Eigen::MatrixXd q=Eigen::MatrixXd::Zero(rows.size(),rows.size());
    for (int r=0;r<int(rows.size());++r) {
        for (auto [i,a]:rows[r]) {
            if (i<0 || i>=x.size() || !std::isfinite(a)) throw std::runtime_error("REGIONAL_BAD_SOURCE_ROW");
            mean(r)+=a*x(i);
        }
        for (int c=0;c<=r;++c) {
            for (auto [i,a]:rows[r]) for (auto [j,b]:rows[c]) {
                if (j<0 || j>=x.size() || !std::isfinite(b)) throw std::runtime_error("REGIONAL_BAD_SOURCE_ROW");
                q(r,c)+=a*b*p(i,j);
            }
            q(c,r)=q(r,c);
        }
    }
    return {std::move(mean),std::move(q)};
}
inline std::vector<ZhangRegionalRow> zhangRegionalTransformRows(
    const Eigen::MatrixXd& t,const std::vector<ZhangRegionalRow>& raw) {
    if (t.cols()!=int(raw.size())) throw std::runtime_error("REGIONAL_TRANSFORM_DIMENSION");
    std::vector<ZhangRegionalRow> out(t.rows());
    for(int r=0;r<t.rows();++r) {
        std::map<int,double> terms;
        for(int c=0;c<t.cols();++c) if(t(r,c)!=0)
            for(auto [i,a]:raw[c]) terms[i]+=t(r,c)*a;
        for(auto [i,a]:terms) if(a!=0) out[r].push_back({i,a});
    }
    return out;
}
