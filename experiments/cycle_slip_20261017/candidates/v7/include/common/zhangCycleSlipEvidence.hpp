#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <functional>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "common/zhangCycleSlipEvidenceState.hpp"

// Preserve the nominal two-sided four-sigma false-alarm budget of the old
// scalar GF detector, but account for the dimension of the correlated test.
inline double zhangSlipPhaseGate(int dimension,int multiplicity=1)
{
    if(multiplicity>1)
    {
        static const std::array<std::array<double,3>,8> family=[]
        {
            std::array<std::array<double,3>,8> out{};
            for(int count=1;count<=8;count++)
            for(int dimension=1;dimension<=3;dimension++)
            {
                const double alpha=std::erfc(4/std::sqrt(2.))/count;
                double low=0,high=128;
                for(int i=0;i<64;i++)
                {
                    const double x=.5*(low+high);
                    const double tail=dimension==1 ? std::erfc(std::sqrt(x/2)) : dimension==2 ? std::exp(-x/2) :
                        std::erfc(std::sqrt(x/2))+std::sqrt(2*x/3.141592653589793)*std::exp(-x/2);
                    if(tail>alpha)low=x;else high=x;
                }
                out[count-1][dimension-1]=.5*(low+high);
            }
            return out;
        }();
        return family[std::max(1,std::min(8,multiplicity))-1][std::max(1,std::min(3,dimension))-1];
    }
    const double alpha=std::erfc(4/std::sqrt(2.));
    if(dimension<=1)return 16;
    if(dimension==2)return -2*std::log(alpha);
    static const double three=[]
    {
        const double probability=std::erfc(4/std::sqrt(2.));double low=0,high=128;
        for(int i=0;i<64;i++)
        {
            const double value=.5*(low+high);
            const double tail=std::erfc(std::sqrt(value/2))+
                std::sqrt(2*value/3.141592653589793)*std::exp(-value/2);
            if(tail>probability)low=value;else high=value;
        }
        return .5*(low+high);
    }();
    return three;
}

// Explicit raw input adapter shared by production and replay. Pseudorange is
// needed only for MW, never for inventing an ionospheric phase correction.
inline ZhangSlipSample zhangMakeSlipSample(double time, int signature,
    const std::array<double,3>& lambda, const std::array<double,3>& phaseMetres,
    const std::array<double,3>& phaseVariance, const std::array<double,2>& codeMetres,
    const std::array<double,2>& codeVariance, bool third, bool code,
    double rangeMinusSatelliteClock, bool geometry)
{
    ZhangSlipSample out;
    out.time=time; out.signature=signature; out.third=third;
    out.code=code; out.geometry=geometry; out.wavelength=lambda;
    const double a=lambda[1]*lambda[1]/(lambda[1]*lambda[1]-lambda[0]*lambda[0]);
    const double b=1-a;
    Eigen::Matrix<double,4,5> transform=Eigen::Matrix<double,4,5>::Zero();
    transform(0,0)=a; transform(0,1)=b;
    transform(1,0)=1; transform(1,1)=-1;
    const double wide=lambda[0]*lambda[1]/(lambda[1]-lambda[0]);
    transform(2,0)=1/lambda[0]; transform(2,1)=-1/lambda[1];
    transform(2,3)=-lambda[1]/((lambda[0]+lambda[1])*wide);
    transform(2,4)=-lambda[0]/((lambda[0]+lambda[1])*wide);
    if (third)
    {
        const double r=(lambda[1]*lambda[1]-lambda[0]*lambda[0])/
            (lambda[2]*lambda[2]-lambda[0]*lambda[0]);
        transform(3,0)=1-r; transform(3,1)=-1; transform(3,2)=r;
    }
    Eigen::Matrix<double,5,1> raw, variance;
    raw << phaseMetres[0],phaseMetres[1],phaseMetres[2],codeMetres[0],codeMetres[1];
    variance << phaseVariance[0],phaseVariance[1],phaseVariance[2],codeVariance[0],codeVariance[1];
    const Eigen::Vector4d value=transform*raw;
    const Eigen::Matrix4d covariance=transform*variance.asDiagonal()*transform.transpose();
    for (int i=0;i<4;i++)
    {
        out.value[i]=value[i];
        for (int j=0;j<4;j++) out.covariance[4*i+j]=covariance(i,j);
    }
    out.rawIf=out.value[0]-rangeMinusSatelliteClock;
    out.value[0]=out.rawIf;
    return out;
}

// Polynomial extrapolation weights; centred/scaled time avoids epoch-sized
// normal equations. Exactly the same weights propagate measurement covariance.
inline Eigen::VectorXd zhangSlipPredictionWeights(const std::vector<ZhangSlipSample>& history,
    double time, int degree)
{
    const int count=history.size();
    degree=std::min(degree,count-1);
    Eigen::MatrixXd design(count,degree+1);
    for (int i=0;i<count;i++)
    {
        const double dt=(history[i].time-time)/30;
        design(i,0)=1;
        for (int j=1;j<=degree;j++) design(i,j)=design(i,j-1)*dt;
    }
    Eigen::VectorXd target=Eigen::VectorXd::Zero(degree+1); target[0]=1;
    return design*(design.transpose()*design).ldlt().solve(target);
}

// Average constant/linear/quadratic PAST-only phase forecasts, including their
// between-model uncertainty. Always fitting a quadratic needlessly amplifies
// endpoint noise on a stable IF/auxiliary phase series. No current candidate
// value participates in the model scores or in these forecast weights.
struct ZhangSlipPhaseForecast
{Eigen::MatrixXd weights,extraCovariance;};

inline ZhangSlipPhaseForecast zhangSlipPhaseForecast(const std::vector<ZhangSlipSample>& history,
    const std::vector<double>& times,int component)
{
    const int count=history.size(),future=times.size();
    Eigen::VectorXd values(count),precision(count);
    for(int k=0;k<count;k++)
    {
        values[k]=history[k].value[component];
        precision[k]=1/std::max(1e-12,history[k].covariance[5*component]);
    }
    std::array<Eigen::MatrixXd,3> weights,variances;
    std::array<Eigen::VectorXd,3> means;
    std::array<double,3> cost{},probability{};
    double minimum=std::numeric_limits<double>::infinity();
    for(int degree=0;degree<=2;degree++)
    {
        Eigen::MatrixXd design(count,degree+1);
        for(int k=0;k<count;k++)
        {
            const double t=(history[k].time-history.back().time)/30;design(k,0)=1;
            for(int j=1;j<=degree;j++)design(k,j)=design(k,j-1)*t;
        }
        const Eigen::MatrixXd normal=design.transpose()*precision.asDiagonal()*design;
        auto decomposition=normal.ldlt();
        const auto coefficients=decomposition.solve(design.transpose()*precision.asDiagonal()*values).eval();
        const Eigen::VectorXd error=values-design*coefficients;
        cost[degree]=(error.array().square()*precision.array()).sum()+(degree+1)*std::log(count);
        minimum=std::min(minimum,cost[degree]);
        Eigen::MatrixXd prediction(future,degree+1);
        for(int k=0;k<future;k++)
        {const double t=(times[k]-history.back().time)/30;prediction(k,0)=1;
         for(int j=1;j<=degree;j++)prediction(k,j)=prediction(k,j-1)*t;}
        weights[degree]=prediction*decomposition.solve(design.transpose()*precision.asDiagonal());
        means[degree]=weights[degree]*values;
        variances[degree]=weights[degree]*precision.cwiseInverse().asDiagonal()*weights[degree].transpose();
    }
    double normalization=0;
    for(int degree=0;degree<=2;degree++)normalization+=probability[degree]=std::exp(-.5*(cost[degree]-minimum));
    ZhangSlipPhaseForecast out;out.weights=Eigen::MatrixXd::Zero(future,count);
    Eigen::VectorXd mean=Eigen::VectorXd::Zero(future);
    for(int degree=0;degree<=2;degree++)
    {probability[degree]/=normalization;out.weights+=probability[degree]*weights[degree];mean+=probability[degree]*means[degree];}
    Eigen::MatrixXd totalVariance=Eigen::MatrixXd::Zero(future,future);
    for(int degree=0;degree<=2;degree++)
    {const Eigen::VectorXd difference=means[degree]-mean;
     totalVariance+=probability[degree]*(variances[degree]+difference*difference.transpose());}
    out.extraCovariance=2*(totalVariance-out.weights*precision.cwiseInverse().asDiagonal()*out.weights.transpose());
    out.extraCovariance=.5*(out.extraCovariance+out.extraCovariance.transpose()).eval();
    // Two independent component-model choices may be correlated; twice their
    // diagonal uncertainty dominates any admissible 2x2 cross-model spread.
    return out;
}

inline Eigen::VectorXd zhangSlipPhasePredictionWeights(const std::vector<ZhangSlipSample>& history,
    double time,int component,double& extraVariance)
{
    extraVariance=0;
    if(history.size()<5)return zhangSlipPredictionWeights(history,time,history.size()>=2 ? 1 : 0);
    const auto out=zhangSlipPhaseForecast(history,{time},component);
    extraVariance=std::max(0.,out.extraCovariance(0,0));return out.weights.row(0).transpose();
}

struct ZhangSlipTrendDiagnostic
{
    bool qualified=false,assisted=false;
    int window=0,stage=0;
    double consistency=0,score=0,support=0,fit=0;
};

struct ZhangSlipGfPrediction
{
    bool valid=false;
    int window=0;
    double modelVariance=0,consistency=0;
    Eigen::VectorXd weights;
};

// Qualify the trend BEFORE inspecting the candidate. Nested predictions share
// raw phase, so their disagreement covariance must include that correlation.
inline ZhangSlipGfPrediction zhangSlipQualifiedGfPrediction(
    const std::vector<ZhangSlipSample>& history,double time)
{
    ZhangSlipGfPrediction out;
    if(history.size()<16 || time-history.back().time<=0 || time-history.back().time>45)return out;
    for(int count:{32,16})
    {
        if(history.size()<(size_t)count)continue;
        std::array<Eigen::VectorXd,2> weights;
        bool valid=true;double modelVariance=0;
        for(int half=0;half<2;half++)
        {
            const int n=count/(half+1);
            std::vector<ZhangSlipSample> past(history.end()-n,history.end());
            Eigen::MatrixXd design(n,3);Eigen::VectorXd values(n);double noise=0;
            for(int k=0;k<n;k++)
            {const double t=(past[k].time-time)/30;design.row(k)<<1,t,t*t;
             values[k]=past[k].value[1];noise+=past[k].covariance[5]/n;}
            const auto fit=(design.transpose()*design).ldlt().solve(design.transpose()*values).eval();
            const Eigen::VectorXd error=values-design*fit;
            double score=0;for(int k=0;k<n;k++)score+=error[k]*error[k]/std::max(1e-12,past[k].covariance[5]);
            const double dof=n-3,x=-std::log(.001);
            valid &= score<=dof+2*std::sqrt(dof*x)+2*x;
            const auto local=zhangSlipPredictionWeights(past,time,2);
            weights[half]=Eigen::VectorXd::Zero(history.size());weights[half].tail(n)=local;
            if(!half)modelVariance=std::max(1e-5*(time-history.back().time),
                std::max(0.,error.squaredNorm()/dof-noise)*(1+local.squaredNorm()));
        }
        if(!valid)continue;
        const Eigen::VectorXd difference=weights[0]-weights[1];
        double level=0,variance=0;
        for(int k=0;k<(int)history.size();k++)
        {level+=difference[k]*history[k].value[1];variance+=difference[k]*difference[k]*history[k].covariance[5];}
        if(variance<=0 || level*level>4*variance)continue;
        out.valid=true;out.window=count;out.weights=weights[0];
        out.modelVariance=modelVariance;out.consistency=level*level/variance;return out;
    }
    return out;
}

inline void zhangAcceptSlipSample(ZhangSlipEvidenceState& state, const ZhangSlipSample& sample)
{
    state.accepted.push_back(sample);
    if (state.accepted.size()>32) state.accepted.erase(state.accepted.begin());
}

// peers are the OTHER links at this receiver and same constellation. The caller
// supplies a snapshot taken before any current-epoch decisions are committed.
// This prevents both candidate self-validation and observation-order dependence.
struct ZhangSlipClockPeer
{
    double previousTime=0, delta=0, variance=0;
    std::vector<ZhangSlipSample> history;
    ZhangSlipSample current;
    bool historyAvailable=false;
};

// Apply ONE fixed peer contrast to the whole temporal window. Integrating
// epoch-wise medians creates a random walk with unrepresented time covariance;
// changing the peer contrast inside a window also manufactures discontinuities.
// Here the common receiver clock cancels algebraically, including clock jumps.
// Peer selection is robust (and not an independent integrity certificate).
inline bool zhangSlipClockFreeWindow(std::vector<ZhangSlipSample>& window,
    const std::vector<ZhangSlipClockPeer>& peers, int& partners)
{
    if(window.size()<2)return false;
    std::vector<std::vector<const ZhangSlipSample*>> available;
    for(const auto& peer:peers)
    {
        if(!peer.historyAvailable)continue;
        std::vector<const ZhangSlipSample*> samples;
        for(const auto& target:window)
        {
            const ZhangSlipSample* found=nullptr;
            if(peer.current.time==target.time)found=&peer.current;
            else for(const auto& old:peer.history)
                if(old.time==target.time){found=&old;break;}
            if(!found || !found->geometry || found->signature!=peer.current.signature)break;
            samples.push_back(found);
        }
        if(samples.size()==window.size())available.push_back(std::move(samples));
    }
    if(available.size()<3)return false;
    std::vector<bool> usable(available.size(),true);
    for(size_t k=1;k<window.size();k++)
    {
        std::vector<double> increments;
        for(const auto& samples:available)increments.push_back(samples[k]->rawIf-samples[k-1]->rawIf);
        auto sorted=increments;std::sort(sorted.begin(),sorted.end());
        const auto middle=sorted.size()/2;
        const double median=sorted.size()%2 ? sorted[middle] : .5*(sorted[middle-1]+sorted[middle]);
        for(size_t p=0;p<available.size();p++)
        {
            const double variance=available[p][k]->covariance[0]+available[p][k-1]->covariance[0];
            if(std::abs(increments[p]-median)>.01+6*std::sqrt(2*variance))usable[p]=false;
        }
    }
    partners=std::count(usable.begin(),usable.end(),true);
    if(partners<3)return false;
    for(size_t k=0;k<window.size();k++)
    {
        double clock=0,variance=0;
        for(size_t p=0;p<available.size();p++)if(usable[p])
        {
            clock+=(available[p][k]->rawIf-available[p][0]->rawIf)/partners;
            variance+=available[p][k]->covariance[0]/(partners*partners);
        }
        window[k].value[0]=window[k].rawIf-clock;
        // The common reference-epoch error cancels because forecast/step
        // contrasts sum to zero. Other raw epochs are independent in this
        // stated model; no integrated-clock white-noise fiction is needed.
        window[k].covariance[0]+=variance;
    }
    return true;
}

inline const ZhangSlipSample* zhangSlipClockReference(const ZhangSlipEvidenceState& state)
{
    // Use the latest common time. Changing a long interval's median peer set
    // is not a gauge-invariant clock transport and can manufacture phase steps.
    for(auto it=state.screening.rbegin();it!=state.screening.rend();++it)
        if(it->geometry)return &*it;
    return nullptr;
}

struct ZhangSlipDopplerEvidence
{
    bool valid=false;
    double innovation=0,variance=0,previousTime=0;
    std::array<double,4> currentCross{},previousCross{};
    std::vector<std::pair<double,std::array<double,4>>> pastCross;
};

inline double zhangSlipMedian(std::vector<double> values)
{
    std::sort(values.begin(),values.end());const auto middle=values.size()/2;
    return values.size()%2 ? values[middle] : .5*(values[middle-1]+values[middle]);
}

inline ZhangSlipDopplerEvidence zhangSlipDopplerEvidence(
    const ZhangSlipEvidenceState& state,const ZhangSlipSample& sample)
{
    ZhangSlipDopplerEvidence out;
    if(!sample.gfRateValid || !std::isfinite(sample.gfRate) || state.screening.empty())return out;
    const auto& previous=state.screening.back();const double dt=sample.time-previous.time;
    if(!previous.gfRateValid || previous.signature!=sample.signature || dt<=0 || dt>45)return out;
    std::vector<double> errors,scales;
    for(size_t k=1;k<state.accepted.size();k++)
    {
        const auto& a=state.accepted[k-1];const auto& b=state.accepted[k];const double interval=b.time-a.time;
        if(!a.gfRateValid || !b.gfRateValid || a.signature!=b.signature || interval<=0 || interval>45)continue;
        const double error=b.value[1]-a.value[1]-.5*interval*(a.gfRate+b.gfRate);
        errors.push_back(error);
        scales.push_back(std::sqrt(a.covariance[5]+b.covariance[5]+.008*.008*std::pow(interval/30,2)));
    }
    const double rawVariance=sample.covariance[5]+previous.covariance[5];
    double variance=rawVariance+.03*.03*std::pow(dt/30,2),bias=0;
    if(errors.size()>=8)
    {
        bias=zhangSlipMedian(errors);std::vector<double> normalized;
        for(size_t i=0;i<errors.size();i++)normalized.push_back(std::abs(errors[i]-bias)/scales[i]);
        const double scale=std::max(1.,1.25*1.4826*zhangSlipMedian(normalized));
        variance=(rawVariance+.008*.008*std::pow(dt/30,2))*scale*scale;
        // Include finite-past median calibration uncertainty. Never estimate
        // a noise scale from the current candidate or its injected phase step.
        variance*=1+1.570796326794897/errors.size();
    }
    out.valid=true;out.previousTime=previous.time;out.variance=variance;
    out.innovation=sample.value[1]-previous.value[1]-.5*dt*(sample.gfRate+previous.gfRate)-bias;
    for(int i=0;i<4;i++)
    {out.currentCross[i]=sample.covariance[4*i+1];out.previousCross[i]=previous.covariance[4*i+1];}
    if(state.accepted.size()>=5 && state.accepted.back().time==previous.time)
    {
        // Whiten correlated phase/Doppler differences in O(n), without a GF
        // trend model. Adjacent differences share phase with opposite signs
        // and Doppler with equal signs: treating them as independent loses
        // information and lets the first noisy endpoint hide a small step.
        std::vector<const ZhangSlipSample*> recent;
        double nextTime=sample.time;
        for(auto it=state.accepted.rbegin();it!=state.accepted.rend();++it)
        {
            if(!it->gfRateValid || it->signature!=sample.signature || nextTime-it->time>45)break;
            recent.push_back(&*it);nextTime=it->time;
        }
        std::reverse(recent.begin(),recent.end());recent.push_back(&sample);
        const double q=.005*.005*std::pow(dt/30,6);
        const double factor=rawVariance>0 ? std::max(0.,variance-rawVariance-q)/(.25*dt*dt*rawVariance) : 0;
        std::vector<double> multipliers(recent.size()-1),beta(recent.size()-1);
        double innovation=0,innovationVariance=0,priorInterval=0;
        bool valid=recent.size()>=6;
        for(size_t k=1;k<recent.size() && valid;k++)
        {
            const auto& a=*recent[k-1];const auto& b=*recent[k];const double interval=b.time-a.time;
            const double diagonal=(a.covariance[5]+b.covariance[5])*(1+.25*interval*interval*factor)+
                .005*.005*std::pow(interval/30,6);
            const double adjacent=k==1 ? 0 : a.covariance[5]*(-1+.25*priorInterval*interval*factor);
            const double multiplier=k==1 ? 0 : adjacent/innovationVariance;
            const double error=b.value[1]-a.value[1]-.5*interval*(a.gfRate+b.gfRate)-bias;
            innovation=error-multiplier*innovation;
            innovationVariance=diagonal-multiplier*adjacent;
            multipliers[k-1]=multiplier;priorInterval=interval;
            valid=innovationVariance>0 && std::isfinite(innovationVariance) && std::isfinite(innovation);
        }
        if(valid)
        {
            beta.back()=1;
            for(int k=(int)beta.size()-2;k>=0;k--)beta[k]=-multipliers[k+1]*beta[k+1];
            out.innovation=innovation;out.variance=innovationVariance;
            for(size_t k=0;k+1<recent.size();k++)
            {
                const double phaseCoefficient=(k ? beta[k-1] : 0)-beta[k];
                std::array<double,4> cross{};
                for(int i=0;i<4;i++)cross[i]=-phaseCoefficient*recent[k]->covariance[4*i+1];
                out.pastCross.push_back({recent[k]->time,cross});
            }
        }
    }
    return out;
}

// Bounded retrospective step-vs-smooth fit. The step column is residualised
// against a quadratic (constant for MW); a smooth GF trend is never a step
// certificate. A radial score of 32 bounds the at-most-58-way two-dimensional
// Gaussian scan below 1e-5 per window (union bound). This is a NOMINAL local error budget,
// not a measured network false alarm guarantee; replay/calibration remain gates.
inline bool zhangSlipScreenStep(const std::vector<ZhangSlipSample>& history,
    Eigen::Vector4d& step, Eigen::Matrix4d& stepCov, double& score, double& eventTime)
{
    const int count=history.size();
    if(count<10 || !history.back().geometry) return false;
    bool allThird=history.back().third;
    for(const auto& sample:history)
    {
        if(!sample.code)return false;
        allThird &= sample.third && sample.thirdSignature==history.back().thirdSignature;
    }
    double maximum=0;
    Eigen::MatrixXd smooth(count,3);
    for(int k=0;k<count;k++)
    {
        if(history[k].signature!=history.back().signature || !history[k].geometry)return false;
        const double t=(history[k].time-history.back().time)/30;
        smooth.row(k)<<1,t,t*t;
    }
    const Eigen::MatrixXd smoothInverse=(smooth.transpose()*smooth).ldlt().solve(smooth.transpose());
    for(int split=5;split<=count-2;split++)
    {
        std::array<Eigen::VectorXd,4> weights;
        Eigen::VectorXd indicator=Eigen::VectorXd::Zero(count);indicator.tail(count-split).setOnes();
        const Eigen::VectorXd direction=indicator-smooth*(smoothInverse*indicator);
        const double information=direction.squaredNorm();
        if(information<=1e-12)continue;
        weights[0]=direction/information;weights[1]=weights[0];weights[3]=weights[0];
        const double fraction=double(count-split)/count;
        weights[2]=(indicator.array()-fraction).matrix()/((count-split)*(1-fraction));
        Eigen::Vector4d amplitude=Eigen::Vector4d::Zero(); Eigen::Matrix4d covariance=Eigen::Matrix4d::Zero();
        for(int i=0;i<4;i++)
            for(int k=0;k<count;k++)
            {
                amplitude[i]+=weights[i][k]*history[k].value[i];
                for(int j=0;j<4;j++)
                    covariance(i,j)+=weights[i][k]*weights[j][k]*history[k].covariance[4*i+j];
            }
        // Model floors prevent precise raw phase from treating tiny unmodelled
        // antenna/windup/clock interpolation effects as a certified step.
        covariance(0,0)+=.01*.01;
        covariance(2,2)+=.05*.05;
        covariance(3,3)+=.002*.002;
        std::vector<int> components{0}; if(allThird)components.push_back(3);
        Eigen::VectorXd residual(components.size());Eigen::MatrixXd cov(components.size(),components.size());
        for(int i=0;i<(int)components.size();i++)
        {
            residual[i]=amplitude[components[i]];
            for(int j=0;j<(int)components.size();j++)cov(i,j)=covariance(components[i],components[j]);
        }
        auto decomposition=cov.ldlt();if(!decomposition.isPositive())continue;
        const double candidate=residual.dot(decomposition.solve(residual));
        if(candidate>maximum)
        {maximum=candidate;step=amplitude;stepCov=covariance;eventTime=history[split].time;}
    }
    score=maximum;
    // chi-square(1) tail at 28 is < 1.22e-7; at most 58 splits.
    return maximum>(allThird ? 32 : 28);
}

// Causal, recent common-cycle matched tests: freeze the reference BEFORE
// each candidate event. Post-event samples may not train a smooth baseline
// through their own jump. Bonferroni allocates 1e-5/window across eight
// scalar tests, including both signs. It is an early quarantine detector,
// never a claim that the fitted common integer has been identified.
inline bool zhangSlipScreenCommonMode(const std::vector<ZhangSlipSample>& raw,
    const std::vector<ZhangSlipClockPeer>& peers,bool rawClock,
    Eigen::Vector4d& step,Eigen::Matrix4d& stepCov,double& score,double& eventTime,
    std::array<double,4>* diagnostic=nullptr,ZhangSlipTrendDiagnostic* trendDiagnostic=nullptr)
{
    static const double gate=[]
    {
        double low=0,high=128;
        for(int i=0;i<64;i++)
        {const double x=.5*(low+high);if(std::erfc(std::sqrt(x/2))>1e-5/16)low=x;else high=x;}
        return .5*(low+high);
    }();
    double best=gate;
    bool winnerGf=false;int winnerWindow=0;double winnerSupport=0,winnerFit=0;
    for(int post=1;post<=8;post++)
    {
        if(raw.size()<(size_t)(post+8))continue;
        const int prior=std::min(64,(int)raw.size()-post),count=prior+post;
        std::vector<ZhangSlipSample> history(raw.end()-count,raw.end());int partners=0;
        if(rawClock && !zhangSlipClockFreeWindow(history,peers,partners))continue;
        bool valid=true,third=history.back().third;
        for(int k=0;k<count;k++)
        {
            const auto& old=history[k];
            valid &= old.geometry && old.code && old.signature==history.back().signature;
            third &= old.third && old.thirdSignature==history.back().thirdSignature;
        }
        if(!valid)continue;
        const std::vector<ZhangSlipSample> before(history.begin(),history.begin()+prior);
        std::vector<double> times;for(int k=prior;k<count;k++)times.push_back(history[k].time);
        const int dim=third ? 2 : 1;
        std::array<ZhangSlipPhaseForecast,2> forecast;
        std::array<double,2> discrepancy{};
        for(int c=0;c<dim;c++)
        {
            const int component=c ? 3 : 0;
            forecast[c]=zhangSlipPhaseForecast(before,times,component);
            Eigen::MatrixXd design(prior,3);Eigen::VectorXd values(prior);double noise=0;
            for(int k=0;k<prior;k++)
            {
                const double t=(before[k].time-before.back().time)/30;
                design.row(k)<<1,t,t*t;values[k]=before[k].value[component];noise+=before[k].covariance[5*component]/prior;
            }
            const auto coefficients=(design.transpose()*design).ldlt().solve(design.transpose()*values).eval();
            discrepancy[c]=std::max(0.,(values-design*coefficients).squaredNorm()/(prior-3)-noise);
        }
        const auto& l=history.back().wavelength;const double a=l[1]*l[1]/(l[1]*l[1]-l[0]*l[0]);
        Eigen::VectorXd common(dim);common[0]=a*l[0]+(1-a)*l[1];
        if(third)
        {const double r=(l[1]*l[1]-l[0]*l[0])/(l[2]*l[2]-l[0]*l[0]);common[1]=(1-r)*l[0]-l[1]+r*l[2];}
        Eigen::VectorXd residual(post*dim),target(post*dim);
        Eigen::MatrixXd cov=Eigen::MatrixXd::Zero(post*dim,post*dim);
        for(int i=0;i<post;i++)for(int c=0;c<dim;c++)
        {
            const int component=c ? 3 : 0,row=i*dim+c;
            residual[row]=history[prior+i].value[component];target[row]=common[c];
            for(int k=0;k<prior;k++)residual[row]-=forecast[c].weights(i,k)*before[k].value[component];
            for(int j=0;j<post;j++)for(int d=0;d<dim;d++)
            {
                const int other=d ? 3 : 0,col=j*dim+d;
                if(i==j)cov(row,col)+=history[prior+i].covariance[4*component+other];
                for(int k=0;k<prior;k++)cov(row,col)+=forecast[c].weights(i,k)*forecast[d].weights(j,k)*before[k].covariance[4*component+other];
                if(c==d)
                {
                    const double rate=c ? .00005 : .0003;
                    cov(row,col)+=forecast[c].extraCovariance(i,j)+discrepancy[c]+
                        rate*rate*(times[i]-before.back().time)*(times[j]-before.back().time);
                }
            }
        }
        auto decomposition=cov.ldlt();if(!decomposition.isPositive())continue;
        const auto inverseTarget=decomposition.solve(target).eval();const double precision=target.dot(inverseTarget);
        if(precision<=0)continue;
        const double phaseMatched=residual.dot(inverseTarget)/std::sqrt(precision);
        double matched=phaseMatched;bool gfAssisted=false;double gfFit=0;int gfWindow=0;
        const auto gf=zhangSlipQualifiedGfPrediction(before,times.front());
        if(gf.valid && std::abs(phaseMatched)>=2)
        {
            std::vector<ZhangSlipSample> gfPast(before.end()-gf.window,before.end());
            Eigen::MatrixXd gfWeights(post,prior);gfWeights.setZero();
            for(int i=0;i<post;i++)gfWeights.row(i).tail(gf.window)=
                zhangSlipPredictionWeights(gfPast,times[i],2).transpose();
            Eigen::VectorXd joinedResidual(post*(dim+1)),joinedTarget(post*(dim+1));
            joinedResidual.head(post*dim)=residual;joinedTarget.head(post*dim)=target;
            joinedTarget.tail(post).setConstant(l[0]-l[1]);
            Eigen::MatrixXd joined=Eigen::MatrixXd::Zero(post*(dim+1),post*(dim+1));
            joined.topLeftCorner(post*dim,post*dim)=cov;
            for(int i=0;i<post;i++)
            {
                joinedResidual[post*dim+i]=history[prior+i].value[1];
                for(int k=0;k<prior;k++)joinedResidual[post*dim+i]-=gfWeights(i,k)*before[k].value[1];
                for(int j=0;j<post;j++)
                {
                    double variance=i==j ? history[prior+i].covariance[5] : 0;
                    for(int k=0;k<prior;k++)variance+=gfWeights(i,k)*gfWeights(j,k)*before[k].covariance[5];
                    // Same discrepancy-rate convention as the one-step GF
                    // predictor, but with shared Brownian time covariance.
                    variance+=gf.modelVariance*std::min(times[i]-before.back().time,times[j]-before.back().time)/
                        (times.front()-before.back().time);
                    joined(post*dim+i,post*dim+j)=variance;
                }
                for(int j=0;j<post;j++)for(int c=0;c<dim;c++)
                {
                    const int component=c ? 3 : 0;
                    double cross=i==j ? history[prior+i].covariance[4*component+1] : 0;
                    for(int k=0;k<prior;k++)cross+=gfWeights(i,k)*forecast[c].weights(j,k)*before[k].covariance[4*component+1];
                    joined(post*dim+i,j*dim+c)=joined(j*dim+c,post*dim+i)=cross;
                }
            }
            auto jointFactor=joined.ldlt();
            if(jointFactor.isPositive())
            {
                const auto inverse=jointFactor.solve(joinedTarget).eval();const double information=joinedTarget.dot(inverse);
                if(information>0)
                {
                    const double jointMatched=joinedResidual.dot(inverse)/std::sqrt(information);
                    const Eigen::VectorXd error=joinedResidual-joinedTarget*(jointMatched>=0 ? 1 : -1);
                    const double fit=error.dot(jointFactor.solve(error)),dof=joinedResidual.size(),x=-std::log(.001);
                    if(jointMatched*phaseMatched>0 && std::abs(jointMatched)>std::abs(phaseMatched) &&
                        fit<=dof+2*std::sqrt(dof*x)+2*x)
                    {matched=jointMatched;gfAssisted=true;gfFit=fit;gfWindow=gf.window;}
                }
            }
        }
        if(diagnostic && matched*matched>(*diagnostic)[0])
            *diagnostic={matched*matched,(double)post,precision,residual[0]};
        if(matched*matched<=best)continue;
        Eigen::MatrixXd repeat=Eigen::MatrixXd::Zero(post*dim,dim);
        for(int k=0;k<post;k++)repeat.block(k*dim,0,dim,dim).setIdentity();
        const Eigen::MatrixXd inverseRepeat=decomposition.solve(repeat);
        const Eigen::MatrixXd information=repeat.transpose()*inverseRepeat;
        const Eigen::MatrixXd phaseCov=information.ldlt().solve(Eigen::MatrixXd::Identity(dim,dim));
        const Eigen::MatrixXd phaseMap=phaseCov*inverseRepeat.transpose();
        const Eigen::VectorXd phaseStep=phaseMap*residual;
        const Eigen::VectorXd error=phaseStep-common*(matched>=0 ? 1 : -1);
        if(error.dot(information*error)>zhangSlipPhaseGate(dim))continue;
        // Recover the complete phase/MW target covariance from raw-variable
        // coefficients; do not drop their shared phase cross-covariance.
        Eigen::MatrixXd map=Eigen::MatrixXd::Zero(4,4*count);
        for(int out=0;out<dim;out++)for(int i=0;i<post;i++)for(int c=0;c<dim;c++)
        {
            const int component=c ? 3 : 0,output=out ? 3 : 0;
            const double weight=phaseMap(out,i*dim+c);map(output,4*(prior+i)+component)+=weight;
            for(int k=0;k<prior;k++)map(output,4*k+component)-=weight*forecast[c].weights(i,k);
        }
        for(int component:{1,2})
        {
            double past=0,now=0;
            for(int k=0;k<count;k++)(k<prior ? past : now)+=1/std::max(1e-12,history[k].covariance[5*component]);
            for(int k=0;k<count;k++)map(component,4*k+component)=
                (k<prior ? -1/past : 1/now)/std::max(1e-12,history[k].covariance[5*component]);
        }
        Eigen::Vector4d amplitude=Eigen::Vector4d::Zero();Eigen::Matrix4d covariance=Eigen::Matrix4d::Zero();
        for(int k=0;k<count;k++)
        {
            Eigen::Map<const Eigen::Vector4d> values(history[k].value.data());
            Eigen::Matrix4d rawCov;for(int i=0;i<4;i++)for(int j=0;j<4;j++)rawCov(i,j)=history[k].covariance[4*i+j];
            const auto block=map.middleCols(4*k,4).eval();amplitude+=block*values;covariance+=block*rawCov*block.transpose();
        }
        for(int c=0;c<dim;c++)for(int d=0;d<dim;d++)covariance(c ? 3 : 0,d ? 3 : 0)=phaseCov(c,d);
        best=matched*matched;step=amplitude;stepCov=covariance;eventTime=history[prior].time;
        winnerGf=gfAssisted;winnerWindow=gfWindow;winnerSupport=phaseMatched;winnerFit=gfFit;
    }
    if(best<=gate)return false;
    if(trendDiagnostic && winnerGf)
    {trendDiagnostic->qualified=true;trendDiagnostic->assisted=true;trendDiagnostic->stage=2;
     trendDiagnostic->window=winnerWindow;trendDiagnostic->score=best;
     trendDiagnostic->support=winnerSupport;trendDiagnostic->fit=winnerFit;}
    score=best;return true;
}

inline ZhangSlipEvidenceDecision zhangInspectSlipEvidence(ZhangSlipEvidenceState& state,
    ZhangSlipSample sample, std::vector<ZhangSlipClockPeer> peers, bool hardBreak=false,
    const ZhangSlipDopplerEvidence& doppler={},ZhangSlipTrendDiagnostic* trendDiagnostic=nullptr)
{
    ZhangSlipEvidenceDecision out;
    if (!std::isfinite(sample.time) || !std::isfinite(sample.rawIf)) return out;
    for (int i=0;i<4;i++)
        if (!std::isfinite(sample.value[i]) || !std::isfinite(sample.covariance[4*i+i]) ||
            sample.covariance[4*i+i]<0) return out;
    const bool observationGap=!state.screening.empty() &&
        (sample.time<=state.screening.back().time || sample.time-state.screening.back().time>90);
    if (hardBreak || observationGap || (!state.accepted.empty() &&
        (sample.time<=state.accepted.back().time ||
         sample.signature!=state.accepted.back().signature)))
        { state.accepted.clear();state.screening.clear();state.phaseSuspect=false; }
    if (state.accepted.empty())
    {
        if(!sample.geometry) return out; // no independent clock/geometry baseline
        zhangAcceptSlipSample(state,sample);state.screening.push_back(sample);out.action=ZhangSlipAction::WARMUP;
        out.usePhase=true; return out;
    }
    const auto& last=state.accepted.back();
    auto rateEvidence=doppler.valid ? doppler : zhangSlipDopplerEvidence(state,sample);
    out.dopplerAvailable=rateEvidence.valid;out.dopplerInnovation=rateEvidence.innovation;
    out.dopplerSigma=std::sqrt(rateEvidence.variance);
    const auto* reference=zhangSlipClockReference(state);
    if(!peers.empty())
        for(const auto& old:state.screening)
            if(old.geometry && old.time==peers.front().previousTime){reference=&old;break;}
    const auto& clockReference=reference ? *reference : last;
    peers.erase(std::remove_if(peers.begin(),peers.end(),[&](const auto& p)
        { return p.previousTime!=clockReference.time || !std::isfinite(p.delta) || p.variance<=0; }),peers.end());
    const bool rawClock=std::any_of(peers.begin(),peers.end(),[](const auto& p){return p.historyAvailable;});
    bool geometry=sample.geometry && clockReference.geometry && peers.size()>=3;
    double clockVariance=0;
    if(geometry && rawClock)
    {
        sample.clock=0;sample.value[0]=sample.rawIf;
    }
    else if (geometry)
    {
        std::sort(peers.begin(),peers.end(),[](auto a,auto b){return a.delta<b.delta;});
        const auto middle=peers.size()/2;
        const double clock=(peers.size()%2) ? peers[middle].delta :
            .5*(peers[middle-1].delta+peers[middle].delta);
        // Conservative heteroscedastic median uncertainty, not independent votes.
        for (const auto& p:peers) clockVariance=std::max(clockVariance,p.variance);
        clockVariance*=1.570796326794897/peers.size();
        sample.clock=clockReference.clock+clock;
        sample.value[0]-=sample.clock;
        out.clockPartners=peers.size();
    }
    else sample.value[0]=last.value[0]; // disabled below, not a fake measurement
    sample.geometry=geometry;
    state.screening.push_back(sample);
    if(state.screening.size()>96)state.screening.erase(state.screening.begin());

    auto forecast=state.accepted;forecast.push_back(sample);
    if(rawClock)
    {
        geometry &= zhangSlipClockFreeWindow(forecast,peers,out.clockPartners);
        if(!geometry){out.action=ZhangSlipAction::NO_EVIDENCE;return out;}
    }

    std::array<Eigen::VectorXd,4> weights;
    std::array<double,4> forecastModelVariance{};
    bool useThird=sample.third, useMw=sample.code;
    for (int i=0;i<4;i++)
    {
        std::vector<ZhangSlipSample> usable;std::vector<int> positions;
        for(int k=0;k<(int)state.accepted.size();k++)
        {
            const auto& old=forecast[k];
            if((i==0 && !old.geometry) || (i==2 && !old.code) ||
               (i==3 && (!old.third || old.thirdSignature!=sample.thirdSignature)))continue;
            usable.push_back(old);positions.push_back(k);
        }
        weights[i]=Eigen::VectorXd::Zero(state.accepted.size());
        if(usable.empty()) {if(i==0)geometry=false;if(i==2)useMw=false;if(i==3)useThird=false;continue;}
        if(i==3 && usable.size()<2)useThird=false;
        const auto local=(i==0 || i==3) ?
            zhangSlipPhasePredictionWeights(usable,sample.time,i,forecastModelVariance[i]) :
            zhangSlipPredictionWeights(usable,sample.time,
                i==2 ? 0 : (usable.size()>=5 ? 2 : (usable.size()>=2 ? 1 : 0)));
        for(int j=0;j<(int)positions.size();j++)weights[i][positions[j]]=local[j];
    }
    Eigen::Vector4d residual;
    Eigen::Matrix4d covariance;
    for (int i=0;i<4;i++)
    {
        double prediction=0;
        for (int k=0;k<(int)state.accepted.size();k++)
            prediction+=weights[i][k]*forecast[k].value[i];
        residual[i]=forecast.back().value[i]-prediction;
        for (int j=0;j<4;j++)
        {
            covariance(i,j)=forecast.back().covariance[4*i+j];
            for (int k=0;k<(int)state.accepted.size();k++)
                covariance(i,j)+=weights[i][k]*weights[j][k]*forecast[k].covariance[4*i+j];
        }
    }
    const double dt=sample.time-last.time;
    // Explicit model discrepancy floors: IF includes unmodelled troposphere,
    // windup/antenna variations; triple removes first-order ionosphere only.
    covariance(0,0)+=clockVariance+std::pow(.0003*dt,2)+forecastModelVariance[0];
    // Forecast ionospheric-model error is distinct from raw phase noise. Use
    // accepted-past lack of fit only; the current suspect never trains Q.
    double gfModelVariance=1e-5*dt;
    if(state.accepted.size()>=5)
    {
        const int count=state.accepted.size();Eigen::MatrixXd design(count,3);Eigen::VectorXd values(count);
        double measurementVariance=0;
        for(int k=0;k<count;k++)
        {
            const double t=(state.accepted[k].time-sample.time)/30;
            design.row(k)<<1,t,t*t;values[k]=state.accepted[k].value[1];
            measurementVariance+=state.accepted[k].covariance[5]/count;
        }
        const auto coefficients=(design.transpose()*design).ldlt().solve(design.transpose()*values).eval();
        const double lackOfFit=(values-design*coefficients).squaredNorm()/(count-3);
        gfModelVariance=std::max(gfModelVariance,
            std::max(0.,lackOfFit-measurementVariance)*(1+weights[1].squaredNorm()));
    }
    covariance(1,1)+=gfModelVariance;
    covariance(2,2)+=1e-4*dt;
    covariance(3,3)+=std::pow(.00005*dt,2)+forecastModelVariance[3];
    // A longer accepted-past window reduces extrapolation noise, but must not
    // claim that real IF/troposphere or auxiliary phase curvature is exact.
    // Add only PAST model lack of fit, never the current innovation.
    for(int component:{0,3})
    {
        std::vector<int> valid;
        for(size_t k=0;k<state.accepted.size();k++)
            if(component==0 || (forecast[k].third && forecast[k].thirdSignature==sample.thirdSignature))valid.push_back(k);
        if(valid.size()<5)continue;
        Eigen::MatrixXd design(valid.size(),3);Eigen::VectorXd values(valid.size());double noise=0;
        for(size_t k=0;k<valid.size();k++)
        {
            const auto& old=forecast[valid[k]];const double t=(old.time-sample.time)/30;
            design.row(k)<<1,t,t*t;values[k]=old.value[component];noise+=old.covariance[5*component]/valid.size();
        }
        const auto coefficients=(design.transpose()*design).ldlt().solve(design.transpose()*values).eval();
        const double lack=(values-design*coefficients).squaredNorm()/(valid.size()-3);
        covariance(component,component)+=std::max(0.,lack-noise)*(1+weights[component].squaredNorm());
    }
    for (int i=0;i<4;i++) { out.innovation[i]=residual[i]; out.sigma[i]=std::sqrt(covariance(i,i)); }
    const bool gfAlarm=std::abs(residual[1])>4*out.sigma[1];
    const bool mwAlarm=useMw && std::abs(residual[2])>4*out.sigma[2];
    std::vector<int> indices;
    if (geometry) indices.push_back(0);
    // GF forecast mismatch alone is not a phase discontinuity certificate:
    // unconstrained ionospheric acceleration can produce precisely that signal.
    // Keep GF in integer identification, not in this phase-integrity trigger.
    if (useThird) indices.push_back(3);
    Eigen::VectorXd phase(indices.size()); Eigen::MatrixXd phaseCov(indices.size(),indices.size());
    for (int i=0;i<(int)indices.size();i++)
    {
        phase[i]=residual[indices[i]];
        for (int j=0;j<(int)indices.size();j++) phaseCov(i,j)=covariance(indices[i],indices[j]);
    }
    if (!indices.empty())
    {
        auto decomposition=phaseCov.ldlt();
        if (!decomposition.isPositive()) return out;
        out.phaseScore=phase.dot(decomposition.solve(phase));
    }
    if(!geometry)
    {
        // A third-frequency GF check alone has a common phase null direction.
        // Passing it is not evidence that all primary ambiguities survived.
        out.action=ZhangSlipAction::NO_EVIDENCE;return out;
    }
    const int family=rateEvidence.valid ? 5 : 2;
    bool phaseAlarm=out.phaseScore>zhangSlipPhaseGate(indices.size(),family) || (gfAlarm && mwAlarm);
    const bool originalPhaseAlarm=phaseAlarm;
    if(!phaseAlarm && geometry && state.accepted.size()>=16)
    {
        const std::vector<ZhangSlipSample> past(forecast.begin(),forecast.end()-1);
        const auto gf=zhangSlipQualifiedGfPrediction(past,sample.time);
        if(trendDiagnostic)
        {trendDiagnostic->qualified=gf.valid;trendDiagnostic->window=gf.window;trendDiagnostic->consistency=gf.consistency;}
        if(gf.valid)
        {
            const int size=indices.size();
            Eigen::VectorXd joinedResidual(size+1),common(size+1);
            Eigen::MatrixXd joined=Eigen::MatrixXd::Zero(size+1,size+1);
            joined.topLeftCorner(size,size)=phaseCov;joinedResidual.head(size)=phase;
            joinedResidual[size]=sample.value[1];joined(size,size)=sample.covariance[5]+gf.modelVariance;
            for(int k=0;k<(int)state.accepted.size();k++)
            {joinedResidual[size]-=gf.weights[k]*forecast[k].value[1];
             joined(size,size)+=gf.weights[k]*gf.weights[k]*forecast[k].covariance[5];}
            for(int i=0;i<size;i++)
            {
                const int component=indices[i];double cross=sample.covariance[4*component+1];
                for(int k=0;k<(int)state.accepted.size();k++)
                    cross+=weights[component][k]*gf.weights[k]*forecast[k].covariance[4*component+1];
                joined(i,size)=joined(size,i)=cross;
            }
            const auto& l=sample.wavelength;const double a=l[1]*l[1]/(l[1]*l[1]-l[0]*l[0]);
            common[0]=a*l[0]+(1-a)*l[1];common[size]=l[0]-l[1];
            if(useThird){const double r=(l[1]*l[1]-l[0]*l[0])/(l[2]*l[2]-l[0]*l[0]);common[1]=(1-r)*l[0]-l[1]+r*l[2];}
            auto decomposition=joined.ldlt();
            if(decomposition.isPositive())
            {
                const auto inverse=decomposition.solve(common).eval();const double information=common.dot(inverse);
                const auto direction=common.head(size).eval();const auto phaseInverse=phaseCov.ldlt().solve(direction).eval();
                const double phaseInformation=direction.dot(phaseInverse);
                if(information>0 && phaseInformation>0)
                {
                    const double matched=joinedResidual.dot(inverse)/std::sqrt(information);
                    const double support=phase.dot(phaseInverse)/std::sqrt(phaseInformation);
                    const Eigen::VectorXd error=joinedResidual-common*(matched>=0 ? 1 : -1);
                    const double fit=error.dot(decomposition.solve(error));
                    const bool assisted=matched*matched>zhangSlipPhaseGate(1,family) &&
                        matched*support>0 && std::abs(support)>=2 && fit<=zhangSlipPhaseGate(size+1,family);
                    if(trendDiagnostic)
                    {trendDiagnostic->score=matched*matched;trendDiagnostic->support=support;
                     trendDiagnostic->fit=fit;trendDiagnostic->assisted=assisted;trendDiagnostic->stage=1;}
                    if(assisted)
                    {phaseAlarm=true;out.changePointScore=matched*matched;out.changePointTime=sample.time;}
                }
            }
        }
    }
    if(rateEvidence.valid && sample.gfRateValid && state.accepted.size()>=5)
    {
        // Check LEVEL consensus, not just endpoint derivative agreement.
        // The current GF cancels from (forecast residual - Doppler residual),
        // so the candidate cannot validate its own predictor. Do not add the
        // GF extrapolator's model-error allowance to this model-mismatch test.
        double consensusVariance=rateEvidence.variance-sample.covariance[5];
        for(size_t k=0;k<state.accepted.size();k++)
        {
            const auto& old=state.accepted[k];double endpointWeight=0;
            if(rateEvidence.pastCross.empty() && old.time==rateEvidence.previousTime)endpointWeight=1;
            else for(const auto& [time,cross]:rateEvidence.pastCross)
                if(time==old.time && old.covariance[5]>0)endpointWeight=cross[1]/old.covariance[5];
            consensusVariance+=weights[1][k]*weights[1][k]*old.covariance[5]-
                2*weights[1][k]*endpointWeight*old.covariance[5];
        }
        const double discrepancy=residual[1]-rateEvidence.innovation;
        if(consensusVariance>0 && discrepancy*discrepancy<=4*consensusVariance)
        {
            const int size=indices.size();Eigen::MatrixXd joined(size+1,size+1);Eigen::VectorXd value(size+1);
            joined.topLeftCorner(size,size)=phaseCov;value.head(size)=phase;
            joined(size,size)=covariance(1,1);value[size]=residual[1];
            for(int i=0;i<size;i++)joined(i,size)=joined(size,i)=covariance(indices[i],1);
            auto decomposition=joined.ldlt();
            if(decomposition.isPositive())
            {
                const double score=value.dot(decomposition.solve(value));
                const double bounded=out.phaseScore+std::min(20.,std::max(0.,score-out.phaseScore));
                if(bounded>zhangSlipPhaseGate(size+1,family))
                {out.dopplerAssisted=true;out.dopplerScore=bounded;phaseAlarm=true;}
            }
        }
    }
    if(rateEvidence.valid && std::isfinite(rateEvidence.innovation) && rateEvidence.variance>0)
    {
        const int size=indices.size();
        Eigen::MatrixXd joined=Eigen::MatrixXd::Zero(size+1,size+1);
        Eigen::VectorXd joinedResidual(size+1);
        joined.topLeftCorner(size,size)=phaseCov;joinedResidual.head(size)=phase;
        joined(size,size)=rateEvidence.variance;joinedResidual[size]=rateEvidence.innovation;
        for(int i=0;i<size;i++)
        {
            const int component=indices[i];double cross=rateEvidence.currentCross[component];
            for(size_t k=0;k<state.accepted.size();k++)
                if(rateEvidence.pastCross.empty() && state.accepted[k].time==rateEvidence.previousTime)
                    cross+=weights[component][k]*rateEvidence.previousCross[component];
                else for(const auto& [time,values]:rateEvidence.pastCross)
                    if(state.accepted[k].time==time)cross+=weights[component][k]*values[component];
            joined(i,size)=joined(size,i)=cross;
        }
        auto decomposition=joined.ldlt();
        if(decomposition.isPositive())
        {
            // DGF observes phase *rate*, so a smooth ionospheric trend is not
            // a jump. Bound its influence: an isolated bad Doppler must not
            // override otherwise clean phase evidence. This is only an early
            // quarantine trigger; the exact-integer retirement gate is intact.
            const double score=joinedResidual.dot(decomposition.solve(joinedResidual));
            const double bounded=out.phaseScore+std::min(20.,std::max(0.,score-out.phaseScore));
            out.dopplerScore=std::max(out.dopplerScore,bounded);
            double integerFit=std::numeric_limits<double>::infinity();
            const auto& l=sample.wavelength;const double a=l[1]*l[1]/(l[1]*l[1]-l[0]*l[0]);
            const double r=useThird ? (l[1]*l[1]-l[0]*l[0])/(l[2]*l[2]-l[0]*l[0]) : 0;
            if(bounded>zhangSlipPhaseGate(size+1,family))
            for(int n1=-2;n1<=2;n1++)for(int n2=-2;n2<=2;n2++)
                for(int n3=(useThird ? -2 : 0);n3<=(useThird ? 2 : 0);n3++)
                {
                    if(!n1 && !n2)continue;
                    Eigen::VectorXd error=joinedResidual;
                    error[0]-=a*l[0]*n1+(1-a)*l[1]*n2;
                    if(useThird)error[1]-=(1-r)*l[0]*n1-l[1]*n2+r*l[2]*n3;
                    error[size]-=l[0]*n1-l[1]*n2;
                    integerFit=std::min(integerFit,error.dot(decomposition.solve(error)));
                }
            out.dopplerIntegerFit=integerFit;
            out.dopplerAssisted |= bounded>zhangSlipPhaseGate(size+1,family) && integerFit<=zhangSlipPhaseGate(size+1,family);
            // A predeclared weak common-cycle direction is one scalar test,
            // not a search over 125 integer hypotheses. This improves power
            // for the observed near-null mode without lowering the family's
            // nominal false-alarm budget. It is NOT an integer certificate.
            Eigen::VectorXd common(size+1);
            common[0]=a*l[0]+(1-a)*l[1];
            if(useThird)common[1]=(1-r)*l[0]-l[1]+r*l[2];
            common[size]=l[0]-l[1];
            const auto inverseCommon=decomposition.solve(common).eval();
            const double information=common.dot(inverseCommon);
            const auto phaseDirection=common.head(size).eval();
            const auto phaseInverse=phaseCov.ldlt().solve(phaseDirection).eval();
            const double phaseInformation=phaseDirection.dot(phaseInverse);
            if(information>0 && phaseInformation>0)
            {
                const double matched=joinedResidual.dot(inverseCommon)/std::sqrt(information);
                const double support=phase.dot(phaseInverse)/std::sqrt(phaseInformation);
                const Eigen::VectorXd error=joinedResidual-common*(matched>=0 ? 1 : -1);
                const double fit=error.dot(decomposition.solve(error));
                out.dopplerAssisted |= matched*matched>zhangSlipPhaseGate(1,family) &&
                    matched*support>0 && std::abs(support)>=2 && fit<=zhangSlipPhaseGate(size+1,family);
            }
            phaseAlarm |= out.dopplerAssisted;
        }
    }
    bool screenedStep=false;
    ZhangSlipTrendDiagnostic windowTrend;
    Eigen::Vector4d step;Eigen::Matrix4d stepCov;
    if(!phaseAlarm && zhangSlipScreenCommonMode(state.screening,peers,rawClock,
        step,stepCov,out.changePointScore,out.changePointTime,nullptr,&windowTrend))
    {
        screenedStep=!windowTrend.assisted;phaseAlarm=true;residual=step;covariance=stepCov;
        if(windowTrend.assisted && trendDiagnostic)*trendDiagnostic=windowTrend;
        for(int i=0;i<4;i++){out.innovation[i]=residual[i];out.sigma[i]=std::sqrt(covariance(i,i));}
    }
    auto screening=state.screening;
    if(screening.size()>64)screening.erase(screening.begin(),screening.end()-64);
    int screenPartners=0;
    const bool screenClock=!rawClock || zhangSlipClockFreeWindow(screening,peers,screenPartners);
    if(!phaseAlarm && screenClock && zhangSlipScreenStep(screening,step,stepCov,
        out.changePointScore,out.changePointTime))
    {
        screenedStep=true;
        phaseAlarm=true;residual=step;covariance=stepCov;
        for(int i=0;i<4;i++){out.innovation[i]=residual[i];out.sigma[i]=std::sqrt(covariance(i,i));}
    }
    if (!state.phaseSuspect && !phaseAlarm && !(gfAlarm && mwAlarm) && !(gfAlarm && !geometry && !useThird))
    {
        out.usePhase=true;
        out.useCode=!mwAlarm;
        out.action=mwAlarm ? ZhangSlipAction::CODE_OUTLIER :
            (gfAlarm ? ZhangSlipAction::IONO_TREND : ZhangSlipAction::ACCEPT);
        // A code outlier must not poison the MW reference or reset a phase arc.
        if (!mwAlarm) zhangAcceptSlipSample(state,sample);
        return out;
    }
    out.action=ZhangSlipAction::PHASE_SUSPECT;
    state.phaseSuspect=true; // a growing extrapolation variance is NOT evidence of recovery
    if(!originalPhaseAlarm && phaseAlarm && !screenedStep)return out;
    // Absence of geometry leaves a null direction. Never declare continuity
    // from a passing GF/MW test or manufacture an integer certificate in it.
    if (!geometry || !useMw) return out;
    // A retrospective fit may have lacked a continuous auxiliary frequency.
    if(screenedStep)
        for(const auto& old:state.screening)
            useThird &= old.third && old.thirdSignature==sample.thirdSignature;
    indices=screenedStep ? std::vector<int>{0,2} : std::vector<int>{0,1,2};
    if(useThird) indices.push_back(3);
    const int dimension=useThird ? 3 : 2, measurements=indices.size();
    Eigen::MatrixXd transform=Eigen::MatrixXd::Zero(measurements,dimension), cov(measurements,measurements);
    Eigen::VectorXd innovation(measurements);
    const auto& l=sample.wavelength;
    const double a=l[1]*l[1]/(l[1]*l[1]-l[0]*l[0]);
    transform(0,0)=a*l[0]; transform(0,1)=(1-a)*l[1];
    const int mwRow=screenedStep ? 1 : 2;
    transform(mwRow,0)=1; transform(mwRow,1)=-1;
    if(!screenedStep){transform(1,0)=l[0];transform(1,1)=-l[1];}
    if (useThird)
    {
        const double r=(l[1]*l[1]-l[0]*l[0])/(l[2]*l[2]-l[0]*l[0]);
        const int tripleRow=measurements-1;
        transform(tripleRow,0)=(1-r)*l[0];transform(tripleRow,1)=-l[1];transform(tripleRow,2)=r*l[2];
    }
    for(int i=0;i<measurements;i++)
    {
        innovation[i]=residual[indices[i]];
        for(int j=0;j<measurements;j++)cov(i,j)=covariance(indices[i],indices[j]);
    }
    auto decomposition=cov.ldlt();
    if (!decomposition.isPositive()) return out;
    const Eigen::MatrixXd normal=transform.transpose()*decomposition.solve(transform);
    const Eigen::VectorXd centre=normal.ldlt().solve(transform.transpose()*decomposition.solve(innovation));
    if (!centre.allFinite() || centre.cwiseAbs().maxCoeff()>100000) return out;
    Eigen::LLT<Eigen::MatrixXd> cholesky(normal);
    if(cholesky.info()!=Eigen::Success)return out;
    const Eigen::MatrixXd upper=cholesky.matrixU();
    // Complete local 2/3-D confidence ellipsoid, not a box around a float
    // solution. A node-budget stop is UNKNOWN, never a uniqueness certificate.
    // Extra coverage margin for receiver-clock interpolation/model error. The
    // nominal Gaussian tail is not a measured real-world failure probability;
    // null/injected holdouts, not ratio alone, gate deployment.
    constexpr double aperture=64;
    const Eigen::VectorXd continuousError=innovation-transform*centre;
    const double continuousCost=std::max(0.,continuousError.dot(decomposition.solve(continuousError)));
    if(!std::isfinite(continuousCost) || continuousCost>aperture)return out;
    const double zeroCost=innovation.dot(decomposition.solve(innovation));
    const bool zeroCompatible=std::isfinite(zeroCost) && zeroCost<=aperture+1e-8;
    std::size_t nodes=0,points=0;bool exhausted=false,zeroSeen=false,ambiguousZero=false;
    double best=std::numeric_limits<double>::infinity();
    std::array<int,3> current{},candidate{};std::array<bool,3> possibleSlip{};
    std::function<void(int,double)> enumerate=[&](int axis,double cost)
    {
        if(ambiguousZero)return;
        if(exhausted || ++nodes>20000){exhausted=true;return;}
        if(axis<0)
        {
            Eigen::VectorXd integer(dimension);
            for(int i=0;i<dimension;i++)integer[i]=current[i];
            const Eigen::VectorXd error=innovation-transform*integer;
            const double score=error.dot(decomposition.solve(error));
            if(!std::isfinite(score) || score>aperture+1e-8)return;
            points++;zeroSeen |= current==std::array<int,3>{0,0,0};
            for(int i=0;i<dimension;i++)possibleSlip[i]=possibleSlip[i] || current[i]!=0;
            if(score<best){best=score;candidate=current;}
            // Once zero and a nonzero integer are both possible, neither
            // retirement nor release can be authorized. No full enumeration
            // is needed to return UNKNOWN (especially for a stale predictor).
            ambiguousZero=zeroCompatible && current!=std::array<int,3>{0,0,0};
            return;
        }
        double shift=0;
        for(int j=axis+1;j<dimension;j++)shift+=upper(axis,j)*(current[j]-centre[j]);
        const double midpoint=centre[axis]-shift/upper(axis,axis);
        const double radius=std::sqrt(std::max(0.,aperture+1e-8-cost))/std::abs(upper(axis,axis));
        if(!std::isfinite(radius) || radius>1e6 || std::abs(midpoint)>1e6){exhausted=true;return;}
        const int low=std::ceil(midpoint-radius-1e-9),high=std::floor(midpoint+radius+1e-9);
        for(int integer=low;integer<=high && !exhausted && !ambiguousZero;integer++)
        {
            current[axis]=integer;
            const double component=upper(axis,axis)*(integer-midpoint);
            const double next=cost+component*component;
            if(next<=aperture+1e-8)enumerate(axis-1,next);
        }
    };
    enumerate(dimension-1,continuousCost);
    if(exhausted || ambiguousZero || points==0)return out;
    out.bestScore=best;
    if(points==1 && best<=9 && zeroSeen)
    {
        // Release a prior quarantine only with uniquely supported zero-step
        // evidence, never just because the stale predictor became uncertain.
        state.phaseSuspect=false;out.action=ZhangSlipAction::ACCEPT;out.usePhase=true;
        zhangAcceptSlipSample(state,sample);return out;
    }
    // Excluding zero is insufficient when no integer step fits the data well:
    // that is model/outlier uncertainty, not permission to retire a phase arc.
    if(zeroSeen || !phaseAlarm || best>9)return out;
    out.breakArc=true;
    if(points==1 && best<=9)
    {
        out.action=ZhangSlipAction::INTEGER_STEP;out.usePhase=true;out.integer=candidate;
        out.ratio=aperture/std::max(best,1e-12); // lower bound, NOT the network AR ratio
        for(int i=0;i<3;i++)out.retire[i]=candidate[i]!=0;
    }
    else
    {
        // Several integers remain, but zero-step is excluded. Retire the union
        // conservatively; do not pretend to identify a particular integer.
        out.action=ZhangSlipAction::UNRESOLVED_STEP;out.retire=possibleSlip;
    }
    state.phaseSuspect=false;
    state.accepted.clear();state.screening.clear();sample.clock=0;sample.value[0]=sample.rawIf;
    zhangAcceptSlipSample(state,sample);
    state.screening.push_back(sample);
    return out;
}
