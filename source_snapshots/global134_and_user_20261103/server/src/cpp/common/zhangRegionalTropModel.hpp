#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// Regional total-ZTD = the SAME GPT2 background + a low-order residual.
// This is an engineering adaptation of regional atmospheric modelling, not
// the seven-parameter MOFC/ERA5 vertical-profile model. Every coefficient is
// in metres: horizontal coordinates are divided by 20 km, height by 1 km.
// The formal covariance is conditional on fixed reference geometry and GPT2.
namespace zhangtrop
{
inline constexpr const char* modelName = "GPT2_ZTD_PLUS_REGIONAL_ENH_RESIDUAL_V1";
inline constexpr double horizontalScale = 20000;
inline constexpr double heightScale = 1000;
inline constexpr double modelVarianceFloor = 0.001; // m^2, NOT sigma or 1 mm

inline Eigen::Vector3d geodetic(const Eigen::Vector3d& xyz)
{
    constexpr double a = 6378137, e2 = 0.0066943799901413165;
    const double p = std::hypot(xyz.x(), xyz.y());
    double lat = std::atan2(xyz.z(), p * (1-e2));
    for (int i=0; i<12; ++i)
    {
        const double n = a / std::sqrt(1-e2*std::pow(std::sin(lat),2));
        lat = std::atan2(xyz.z()+e2*n*std::sin(lat), p);
    }
    const double n = a / std::sqrt(1-e2*std::pow(std::sin(lat),2));
    return {lat, std::atan2(xyz.y(),xyz.x()), p/std::cos(lat)-n};
}

inline Eigen::Vector4d basis(const Eigen::Vector3d& origin, const Eigen::Vector3d& xyz)
{
    const auto llh = geodetic(origin);
    const Eigen::Vector3d d = xyz-origin;
    const double e = -std::sin(llh.y())*d.x()+std::cos(llh.y())*d.y();
    const double n = -std::sin(llh.x())*std::cos(llh.y())*d.x()
        -std::sin(llh.x())*std::sin(llh.y())*d.y()+std::cos(llh.x())*d.z();
    return {1,e/horizontalScale,n/horizontalScale,(geodetic(xyz).z()-llh.z())/heightScale};
}

struct Fit
{
    bool valid = false;
    std::string reason = "NOT_FITTED";
    Eigen::Vector4d coefficients = Eigen::Vector4d::Zero();
    Eigen::MatrixXd jacobian; // coefficients = jacobian * reference residuals
    Eigen::Matrix4d leverage = Eigen::Matrix4d::Zero(); // J J^T, not formal Q
    Eigen::VectorXd leaveOneOutResiduals, leaveOneOutFormalVariances;
    int order = 0;
    double condition = 0, looRms = 0;
    double spatialVariance = modelVarianceFloor;
};

inline bool inverseDesign(const Eigen::MatrixXd& design, Eigen::MatrixXd& inverse,
                          double& condition)
{
    const int k=design.cols();
    if (k<1 || k>4 || design.rows()<k || !design.allFinite()) return false;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(design,Eigen::ComputeThinU|Eigen::ComputeThinV);
    const auto s=svd.singularValues();
    if (s.size()!=k || !s.allFinite() || s(0)<=0 || s(k-1)<=s(0)*1e-10) return false;
    condition=s(0)/s(k-1);
    if (condition>1e4) return false;
    inverse=svd.matrixV()*s.cwiseInverse().asDiagonal()*svd.matrixU().transpose();
    return inverse.allFinite() &&
        (inverse*design-Eigen::MatrixXd::Identity(k,k)).cwiseAbs().maxCoeff()<1e-8;
}

inline Fit fit(const Eigen::MatrixXd& design, const Eigen::VectorXd& residual,
               const Eigen::MatrixXd& covariance)
{
    Fit out;
    const int n=design.rows();
    if (n<8 || design.cols()!=4 || residual.size()!=n || covariance.rows()!=n || covariance.cols()!=n)
    { out.reason="INSUFFICIENT_REFERENCES_OR_DIMENSION"; return out; }
    if (!design.allFinite() || !residual.allFinite() || !covariance.allFinite())
    { out.reason="NONFINITE_SOURCE"; return out; }
    const double scale=std::max(1e-16,covariance.cwiseAbs().maxCoeff());
    if ((covariance-covariance.transpose()).cwiseAbs().maxCoeff()>1e-10*scale)
    { out.reason="ASYMMETRIC_SOURCE_Q"; return out; }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(covariance/scale,Eigen::EigenvaluesOnly);
    if (eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff() < -1e-9)
    { out.reason="INDEFINITE_SOURCE_Q"; return out; }
    // Geometry-only, predeclared fallback; never choose order from user truth.
    // Require all leave-one-out designs to support the SAME chosen model.
    for (int k : {4,3,1})
    {
        Eigen::MatrixXd j;
        double condition=0;
        if (!inverseDesign(design.leftCols(k),j,condition)) continue;
        bool supported=true;
        for (int omitted=0; omitted<n; ++omitted)
        {
            Eigen::MatrixXd reduced(n-1,k), ji;
            for (int i=0,l=0; i<n; ++i) if(i!=omitted) reduced.row(l++)=design.row(i).head(k);
            double c=0;
            if (!inverseDesign(reduced,ji,c) || (design.row(omitted).head(k)*ji).squaredNorm()>16)
            { supported=false; break; }
        }
        if (!supported) continue;
        out.order=k; out.condition=condition;
        out.jacobian=Eigen::MatrixXd::Zero(4,n); out.jacobian.topRows(k)=j;
        break;
    }
    if (!out.order) { out.reason="RANK_OR_CONDITION"; return out; }
    out.coefficients=out.jacobian*residual;
    out.leverage=out.jacobian*out.jacobian.transpose();
    out.leaveOneOutResiduals.resize(n); out.leaveOneOutFormalVariances.resize(n);
    double excess=0;
    for (int omitted=0; omitted<n; ++omitted)
    {
        Eigen::MatrixXd reduced(n-1,out.order), inverse;
        for (int i=0,k=0; i<n; ++i) if (i!=omitted) reduced.row(k++)=design.row(i).head(out.order);
        double condition=0;
        if (!inverseDesign(reduced,inverse,condition))
        { out.reason="LEAVE_ONE_OUT_RANK_OR_CONDITION"; return out; }
        const Eigen::RowVectorXd w=design.row(omitted).head(out.order)*inverse;
        if (!w.allFinite() || w.squaredNorm()>16)
        { out.reason="LEAVE_ONE_OUT_EXCESSIVE_LEVERAGE"; return out; }
        Eigen::VectorXd contrast=Eigen::VectorXd::Zero(n); contrast(omitted)=1;
        for (int i=0,k=0; i<n; ++i) if (i!=omitted) contrast(i)=-w(k++);
        const double error=contrast.dot(residual);
        const double variance=contrast.dot(covariance*contrast);
        if (!std::isfinite(variance) || variance < -1e-9*scale)
        { out.reason="LEAVE_ONE_OUT_VARIANCE_INVALID"; return out; }
        out.leaveOneOutResiduals(omitted)=error;
        out.leaveOneOutFormalVariances(omitted)=std::max(0.0,variance);
        // Estimate EXTRA local model variance, after source estimation noise.
        // The floor remains experimental, not a certified integrity bound.
        excess+=(error*error-std::max(0.0,variance))/(1+w.squaredNorm());
    }
    out.looRms=std::sqrt(out.leaveOneOutResiduals.squaredNorm()/n);
    if (out.leaveOneOutResiduals.cwiseAbs().maxCoeff()>0.20)
    { out.reason="LEAVE_ONE_OUT_GROSS_MISMATCH"; return out; }
    // Fixed extra USER prediction variance budget. LOO is diagnostic/gating
    // only: no online tuning on this experiment's observations or user truth.
    // LOO cannot detect a network-wide common ZTD bias. This is not an
    // integrity bound and does not represent temporal error correlations.
    out.spatialVariance=modelVarianceFloor;
    out.valid=true; out.reason="NONE";
    return out;
}
} // namespace zhangtrop
