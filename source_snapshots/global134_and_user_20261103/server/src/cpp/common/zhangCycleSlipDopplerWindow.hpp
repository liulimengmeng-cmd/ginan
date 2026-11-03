#pragma once

// Causal event-window GF/Doppler contrasts. This helper does not
// decide slips: the caller must also require independent phase support.
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <vector>
#include "common/zhangCycleSlipEvidenceState.hpp"

struct ZhangSlipDopplerWindow
{
    bool valid=false;
    Eigen::VectorXd residual;
    Eigen::MatrixXd covariance,phaseMap;
    double rateVarianceFactor=0;
};

inline ZhangSlipDopplerWindow zhangSlipDopplerWindow(
    const std::vector<ZhangSlipSample>& history,int prior)
{
    ZhangSlipDopplerWindow out;
    const int count=history.size(),past=prior-1,post=count-prior;
    if(past<8 || post<1)return out;
    std::vector<double> errors,scales;
    Eigen::VectorXd intervals(count-1),data(count-1),diagonal(count-1),adjacent(count-1);
    auto median=[](std::vector<double> values)
    {std::sort(values.begin(),values.end());const int n=values.size();
     return n%2 ? values[n/2] : .5*(values[n/2-1]+values[n/2]);};
    for(int k=0;k<count;k++)
    {
        if(!history[k].gfRateValid || !std::isfinite(history[k].gfRate) ||
            history[k].covariance[5]<=0 || history[k].signature!=history.back().signature)return out;
        if(!k)continue;
        const double dt=history[k].time-history[k-1].time;
        if(dt<=0 || dt>45)return out;
        intervals[k-1]=dt;
        data[k-1]=history[k].value[1]-history[k-1].value[1]-
            .5*dt*(history[k].gfRate+history[k-1].gfRate);
        if(k<prior)
        {errors.push_back(data[k-1]);scales.push_back(std::sqrt(history[k].covariance[5]+
            history[k-1].covariance[5]+.008*.008*std::pow(dt/30,2)));}
    }
    const double bias=median(errors);
    std::vector<double> normalized;
    for(int k=0;k<past;k++)normalized.push_back(std::abs(errors[k]-bias)/scales[k]);
    const double scale=std::max(1.,1.25*1.4826*median(normalized));
    const double dt=intervals[past],rawVariance=history[prior-1].covariance[5]+history[prior].covariance[5];
    const double q=.005*.005*std::pow(dt/30,6);
    // Calibration uses only pre-event intervals. The unknown constant bias is
    // fitted below with its FULL predictive uncertainty, not independent
    // per-interval copies of an estimated median's error.
    out.rateVarianceFactor=std::max(0.,(rawVariance+.008*.008*std::pow(dt/30,2))*scale*scale-
        rawVariance-q)/(.25*dt*dt*rawVariance);
    for(int k=0;k<count-1;k++)
    {
        const double interval=intervals[k];
        diagonal[k]=(history[k].covariance[5]+history[k+1].covariance[5])*
            (1+.25*interval*interval*out.rateVarianceFactor)+.005*.005*std::pow(interval/30,6);
        adjacent[k]=k ? history[k].covariance[5]*(-1+
            .25*intervals[k-1]*interval*out.rateVarianceFactor) : 0;
    }
    // Two tridiagonal solves: C_past^{-1}1 and C_past^{-1}e_last.
    // Only the first future interval shares a raw endpoint with the past.
    Eigen::VectorXd pivots(past),multipliers=Eigen::VectorXd::Zero(past);
    pivots[0]=diagonal[0];
    for(int k=1;k<past;k++)
    {multipliers[k]=adjacent[k]/pivots[k-1];pivots[k]=diagonal[k]-multipliers[k]*adjacent[k];}
    if(!pivots.allFinite() || pivots.minCoeff()<=0)return out;
    auto solve=[&](Eigen::VectorXd value)
    {for(int k=1;k<past;k++)value[k]-=multipliers[k]*value[k-1];
     value.array()/=pivots.array();
     for(int k=past-2;k>=0;k--)value[k]-=multipliers[k+1]*value[k+1];return value;};
    const Eigen::VectorXd inverseOne=solve(Eigen::VectorXd::Ones(past));
    Eigen::VectorXd endpoint=Eigen::VectorXd::Zero(past);endpoint[past-1]=1;
    const Eigen::VectorXd inverseLast=solve(endpoint);
    const double information=inverseOne.sum(),boundary=adjacent[past];
    if(!std::isfinite(information) || information<=0)return out;
    Eigen::VectorXd biasDirection=Eigen::VectorXd::Ones(post);
    biasDirection[0]-=boundary*inverseOne[past-1];
    Eigen::MatrixXd prediction=biasDirection*(inverseOne/information).transpose();
    prediction.row(0)+=boundary*inverseLast.transpose();
    out.residual=data.tail(post)-prediction*data.head(past);
    out.covariance=biasDirection*biasDirection.transpose()/information;
    for(int k=0;k<post;k++)
    {out.covariance(k,k)+=diagonal[past+k];
     if(k)out.covariance(k,k-1)+=adjacent[past+k],out.covariance(k-1,k)+=adjacent[past+k];}
    out.covariance(0,0)-=boundary*boundary*inverseLast[past-1];
    // Coefficients of all raw GF phase samples, needed for exact cross-
    // covariance with the IF and third-frequency phase contrasts.
    out.phaseMap=Eigen::MatrixXd::Zero(post,count);
    for(int k=0;k<post;k++)
    {out.phaseMap(k,prior+k)=1;out.phaseMap(k,prior+k-1)=-1;
     for(int j=0;j<past;j++)
     {out.phaseMap(k,j)+=prediction(k,j);out.phaseMap(k,j+1)-=prediction(k,j);}}
    out.valid=out.residual.allFinite() && out.covariance.allFinite() && out.covariance.ldlt().isPositive();
    return out;
}
