#pragma once
#include "common/eigenIncluder.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <stdexcept>
#include <string>

struct ZhangUserJointChartResult
{
    Eigen::VectorXd mean;
    Eigen::MatrixXd covariance;
    Eigen::MatrixXd integerMap;
    Eigen::VectorXd integerOffset;
};

inline void zhangUserRequireJointPsd(const Eigen::MatrixXd& covariance)
{
    if (!covariance.rows() || covariance.rows()!=covariance.cols() || !covariance.allFinite() ||
        (covariance-covariance.transpose()).cwiseAbs().maxCoeff()>
            1e-11*std::max(1.0,covariance.cwiseAbs().maxCoeff()))
        throw std::runtime_error("INVALID_FULL_JOINT_COVARIANCE");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
    if (solver.info()!=Eigen::Success || solver.eigenvalues().minCoeff()<
        -1e-11*std::max(1.0,solver.eigenvalues().cwiseAbs().maxCoeff()))
        throw std::runtime_error("FULL_JOINT_COVARIANCE_NOT_PSD");
}

// Exact affine coordinate change on an augmented [x,delta] posterior. The
// caller MUST provide the full joint covariance and a same-instant chart
// displacement contract. Qt from a product file is not accepted as Qdelta.
// delta is deliberately retained, so n's dependence on it is not forgotten.
// This primitive does not invent or install an inter-epoch stochastic model.
inline ZhangUserJointChartResult zhangUserTransportJointChart(
    const Eigen::VectorXd& mean,const Eigen::MatrixXd& jointCovariance,
    const Eigen::MatrixXd& T,const Eigen::MatrixXd& E,const Eigen::VectorXd& shift,
    const Eigen::MatrixXd& K,const Eigen::VectorXd& integerOffset,
    const std::string& covarianceTarget,bool fullJointCrossBlocksProvided)
{
    const int p=T.rows(), q=E.cols();
    if (covarianceTarget!="SAME_INSTANT_CHART_DISPLACEMENT" || !fullJointCrossBlocksProvided)
        throw std::runtime_error("JOINT_CHART_DISPLACEMENT_STATISTICS_UNAVAILABLE");
    if (p<=0 || T.cols()!=p || E.rows()!=p || shift.size()!=p || mean.size()!=p+q ||
        jointCovariance.rows()!=p+q || K.cols()!=p+q || K.rows()!=integerOffset.size() ||
        !T.allFinite() || !E.allFinite() || !mean.allFinite() || !shift.allFinite() ||
        !K.allFinite() || !integerOffset.allFinite()) throw std::runtime_error("INVALID_JOINT_CHART_DIMENSIONS");
    zhangUserRequireJointPsd(jointCovariance);
    Eigen::FullPivLU<Eigen::MatrixXd> factor(T);
    if (!factor.isInvertible()) throw std::runtime_error("NONINVERTIBLE_REAL_CHART_MAP");
    Eigen::MatrixXd A=Eigen::MatrixXd::Identity(p+q,p+q);
    A.topLeftCorner(p,p)=T; A.topRightCorner(p,q)=E;
    Eigen::MatrixXd inverse=Eigen::MatrixXd::Identity(p+q,p+q);
    inverse.topLeftCorner(p,p)=factor.inverse();
    inverse.topRightCorner(p,q)=-inverse.topLeftCorner(p,p)*E;
    Eigen::VectorXd d=Eigen::VectorXd::Zero(p+q); d.head(p)=shift;
    ZhangUserJointChartResult out;
    out.mean=A*mean+d;
    out.covariance=A*jointCovariance*A.transpose();
    out.integerMap=K*inverse;
    out.integerOffset=integerOffset-out.integerMap*d;
    return out;
}

// R is positive definite observation-only noise. B carries unconstrained real
// chart nuisance parameters with no finite prior. The returned map eliminates
// their whitened column space; a full-row-rank B returns zero observations.
inline Eigen::MatrixXd zhangUserUnknownChartProjection(
    const Eigen::MatrixXd& R,const Eigen::MatrixXd& B)
{
    if (!R.rows() || R.rows()!=R.cols() || B.rows()!=R.rows() || !R.allFinite() || !B.allFinite() ||
        (R-R.transpose()).cwiseAbs().maxCoeff()>1e-11*std::max(1.0,R.cwiseAbs().maxCoeff()))
        throw std::runtime_error("INVALID_UNKNOWN_CHART_PROJECTION_INPUT");
    Eigen::LLT<Eigen::MatrixXd> factor(R);
    if (factor.info()!=Eigen::Success) throw std::runtime_error("OBSERVATION_NOISE_NOT_POSITIVE_DEFINITE");
    const Eigen::MatrixXd whiten=factor.matrixL().solve(Eigen::MatrixXd::Identity(R.rows(),R.rows()));
    if (B.cols()==0) return whiten;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(whiten*B,Eigen::ComputeFullU);
    svd.setThreshold(1e-12);
    const int residual=R.rows()-svd.rank();
    return svd.matrixU().rightCols(residual).transpose()*whiten;
}
