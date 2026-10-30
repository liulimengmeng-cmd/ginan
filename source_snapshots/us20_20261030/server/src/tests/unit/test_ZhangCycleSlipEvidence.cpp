#include "common/zhangCycleSlipEvidence.hpp"
#include "common/zhangPhaseArcDecision.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>

static void check(bool condition,const char* text) { if(!condition)throw std::runtime_error(text); }
static const std::array<double,3> wavelength{299792458./1575420000,299792458./1227600000,299792458./1176450000};
static ZhangSlipSample synthetic(int epoch,std::array<int,3> step={},double clock=0,double codeStep=0)
{
    const double iono=2+.02*epoch+.001*epoch*epoch;
    std::array<double,3> phase;
    for(int i=0;i<3;i++)phase[i]=1000+clock-std::pow(wavelength[i]/wavelength[0],2)*iono+wavelength[i]*step[i];
    return zhangMakeSlipSample(epoch*30,123,wavelength,phase,{4e-6,4e-6,4e-6},
        {1000+clock+iono+codeStep,1000+clock+std::pow(wavelength[1]/wavelength[0],2)*iono},
        {.16,.16},true,true,1000,true);
}
static std::vector<ZhangSlipClockPeer> peers(double time,double delta=0)
{return {{time,delta,.00008},{time,delta,.00008},{time,delta,.00008},{time,delta,.00008}};}

// Independent archive writer matching the previously shipped sample layout.
// No call to the new sample serializer: old version-0 bytes must migrate with
// missing Doppler explicitly unavailable, not silently invented as zero rate.
struct LegacySlipSample
{
    double time=0,rawIf=0,clock=0;int signature=0,thirdSignature=0;
    bool geometry=false,third=false,code=false;
    std::array<double,4> value{};std::array<double,16> covariance{};std::array<double,3> wavelength{};
    template<class A>void serialize(A& a,const unsigned int&)
    {a & time & rawIf & clock & signature & thirdSignature & geometry & third & code & value & covariance & wavelength;}
};
struct LegacySlipState
{
    std::vector<LegacySlipSample> accepted,screening;bool phaseSuspect=false;
    template<class A>void serialize(A& a,const unsigned int&){a & accepted & screening & phaseSuspect;}
};
static void extendedEvidenceUnitTests()
{
    const double alpha=std::erfc(4/std::sqrt(2.));
    for(int family:{1,2,4,6})for(int dim=1;dim<=3;dim++)
    {
        const double q=zhangSlipPhaseGate(dim,family);
        const double tail=dim==1 ? std::erfc(std::sqrt(q/2)) : dim==2 ? std::exp(-q/2) :
            std::erfc(std::sqrt(q/2))+std::sqrt(2*q/3.141592653589793)*std::exp(-q/2);
        check(std::abs(tail*family-alpha)<1e-13,"instantaneous family did not preserve its nominal tail budget");
    }
    {
        constexpr int prior=32,post=3,count=prior+post,dim=2;
        std::vector<ZhangSlipSample> before,after;
        for(int k=0;k<count;k++)(k<prior ? before : after).push_back(synthetic(k));
        const auto gf=zhangSlipQualifiedGfPrediction(before,after.front().time);
        check(gf.valid && gf.window==32,"quadratic past GF failed qualification");
        Eigen::MatrixXd design(prior,3),future(post,3);
        for(int k=0;k<prior;k++){const double t=k-prior+1.;design.row(k)<<1,t,t*t;}
        for(int k=0;k<post;k++){const double t=k+1.;future.row(k)<<1,t,t*t;}
        // Independent QR inverse, not the predictor's normal equations.
        const Eigen::MatrixXd gfWeights=future*design.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(prior,prior));
        check((gf.weights-gfWeights.row(0).transpose()).norm()<1e-10,"GF predictor disagrees with independent QR");
        auto nonPolynomial=before;
        for(int k=0;k<prior;k++)nonPolynomial[k].value[1]+=.25*std::sin(.7*k);
        check(!zhangSlipQualifiedGfPrediction(nonPolynomial,after.front().time).valid,
            "non-polynomial past GF was accepted as an exact trend");
        std::array<ZhangSlipPhaseForecast,2> forecast;
        for(int c=0;c<dim;c++)forecast[c].weights=Eigen::MatrixXd::Ones(post,prior)/prior;
        Eigen::VectorXd values(4*count);
        Eigen::MatrixXd rawCov=Eigen::MatrixXd::Zero(4*count,4*count);
        for(int k=0;k<count;k++)
        {
            const auto& sample=k<prior ? before[k] : after[k-prior];
            for(int c=0;c<4;c++)
            {values[4*k+c]=sample.value[c];for(int d=0;d<4;d++)rawCov(4*k+c,4*k+d)=sample.covariance[4*c+d];}
        }
        Eigen::MatrixXd contrast=Eigen::MatrixXd::Zero(post*(dim+1),4*count);
        for(int i=0;i<post;i++)
        {
            for(int c=0;c<dim;c++)
            {contrast(i*dim+c,4*(prior+i)+(c ? 3 : 0))=1;
             for(int k=0;k<prior;k++)contrast(i*dim+c,4*k+(c ? 3 : 0))=-1./prior;}
            contrast(post*dim+i,4*(prior+i)+1)=1;
            for(int k=0;k<prior;k++)contrast(post*dim+i,4*k+1)=-gfWeights(i,k);
        }
        Eigen::MatrixXd process=Eigen::MatrixXd::Zero(post*(dim+1),post*(dim+1));
        for(int i=0;i<post;i++)for(int j=0;j<post;j++)
        {
            for(int c=0;c<dim;c++)process(i*dim+c,j*dim+c)=std::pow(c ? .00005 : .0003,2)*30*30*(i+1)*(j+1);
            process(post*dim+i,post*dim+j)=gf.modelVariance*std::min(i+1,j+1);
        }
        const Eigen::VectorXd expectedResidual=contrast*values;
        const Eigen::MatrixXd expectedCov=contrast*rawCov*contrast.transpose()+process;
        const Eigen::VectorXd target=Eigen::VectorXd::Ones(post*dim);
        const auto moments=zhangSlipAugmentPhaseWithGfWindow(before,after,forecast,
            expectedResidual.head(post*dim),expectedCov.topLeftCorner(post*dim,post*dim),target,gf);
        check(moments.valid && (moments.residual-expectedResidual).norm()<5e-10 &&
            (moments.covariance-expectedCov).norm()<1e-12,
            "GF temporal/cross covariance disagrees with independent raw-variable matrix construction");
        check(std::abs(moments.covariance(0,post*dim))>1e-12 &&
            moments.covariance(post*dim,post*dim+1)>=gf.modelVariance,
            "GF window dropped shared-phase or future-time covariance");
    }
    {
        constexpr int count=20,prior=18,past=prior-1,post=count-prior,m=count-1;
        std::vector<ZhangSlipSample> history;
        for(int k=0;k<count;k++)
        {
            auto sample=synthetic(k);sample.gfRateValid=true;
            sample.gfRate=(std::pow(wavelength[1]/wavelength[0],2)-1)*(.02+.002*k)/30;
            sample.value[1]+=.001*std::cos(k);
            for(auto& value:sample.covariance)value*=1+.02*k;
            history.push_back(sample);
        }
        const auto window=zhangSlipDopplerWindow(history,prior);
        const double rawPair=history[prior-1].covariance[5]+history[prior].covariance[5];
        const double rateFactor=(.008*.008-.005*.005)/(225*rawPair);
        check(window.valid && std::abs(window.rateVarianceFactor-rateFactor)<1e-12,
            "Doppler window floor calibration changed unexpectedly");
        Eigen::MatrixXd transform=Eigen::MatrixXd::Zero(m,2*count+m);
        Eigen::VectorXd data=Eigen::VectorXd::Zero(2*count+m),variance(2*count+m);
        for(int k=0;k<count;k++)
        {data[k]=history[k].value[1];data[count+k]=history[k].gfRate;
         variance[k]=history[k].covariance[5];variance[count+k]=rateFactor*history[k].covariance[5];}
        for(int k=0;k<m;k++)
        {transform(k,k)=-1;transform(k,k+1)=1;transform(k,count+k)=transform(k,count+k+1)=-15;
         transform(k,2*count+k)=1;variance[2*count+k]=.005*.005;}
        const Eigen::MatrixXd covariance=transform*variance.asDiagonal()*transform.transpose();
        const Eigen::MatrixXd inverse=covariance.topLeftCorner(past,past).colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(past,past));
        const Eigen::VectorXd inverseOne=inverse*Eigen::VectorXd::Ones(past);
        const Eigen::MatrixXd cross=covariance.bottomLeftCorner(post,past);
        const Eigen::VectorXd bias=Eigen::VectorXd::Ones(post)-cross*inverseOne;
        const Eigen::MatrixXd prediction=cross*inverse+bias*(inverseOne/inverseOne.sum()).transpose();
        Eigen::MatrixXd contrast(post,m);contrast.leftCols(past)=-prediction;contrast.rightCols(post).setIdentity();
        const Eigen::VectorXd expected=contrast*transform*data;
        const Eigen::MatrixXd expectedCov=contrast*covariance*contrast.transpose();
        check((window.residual-expected).norm()<1e-10 && (window.covariance-expectedCov).norm()<1e-12 &&
            (window.phaseMap-contrast*transform.leftCols(count)).norm()<1e-10,
            "Doppler window disagrees with independent dense raw-noise/GLS-bias reference");
        history.back().gfRateValid=false;
        check(!zhangSlipDopplerWindow(history,prior).valid,"missing Doppler became an invented rate");
    }
}

static void unitTests()
{
    extendedEvidenceUnitTests();
    ZhangSlipEvidenceState smooth;
    for(int epoch=0;epoch<30;epoch++)
    {
        const auto out=zhangInspectSlipEvidence(smooth,synthetic(epoch),peers(30*(epoch-1)));
        check(out.usePhase && !out.breakArc,"smooth ionosphere rejected");
        check(smooth.accepted.size()<=32,"unbounded history");
        check(smooth.screening.size()<=96,"unbounded screening history");
    }
    for(auto slip:std::vector<std::array<int,3>>{{1,0,0},{0,1,0},{1,1,0},{1,1,1},
        {-1,-1,-1},{9,7,0},{77,60,0},{60,77,0},{0,0,1}})
    {
        auto state=smooth;
        const auto out=zhangInspectSlipEvidence(state,synthetic(30,slip),peers(870));
        check(!out.usePhase || out.breakArc,"injected integer silently used");
        if(!out.breakArc) check(state.accepted.size()==smooth.accepted.size(),"suspect learned into history");
        if(out.breakArc) check(out.integer==slip,"wrong integer certificate");
    }
    auto state=smooth;
    auto gfTrend=synthetic(30);gfTrend.value[1]+=.5;
    auto trend=zhangInspectSlipEvidence(state,gfTrend,peers(870));
    check(trend.usePhase && !trend.breakArc && trend.action==ZhangSlipAction::IONO_TREND,
        "GF-only ionospheric model mismatch retired or quarantined phase");
    state=smooth;
    auto out=zhangInspectSlipEvidence(state,synthetic(30,{},300000),peers(870,300000));
    check(out.usePhase && !out.breakArc,"common receiver clock jump became a slip");
    state=smooth;
    out=zhangInspectSlipEvidence(state,synthetic(30,{},0,100),peers(870));
    check(out.usePhase && !out.breakArc && !out.useCode,"code outlier retired phase");
    state=smooth;
    out=zhangInspectSlipEvidence(state,synthetic(30,{1,0,0}),{});
    check(!out.usePhase && !out.breakArc,"missing geometry evidence certified a slip");
    state=smooth;
    out=zhangInspectSlipEvidence(state,synthetic(30),{});
    check(!out.usePhase && !out.breakArc,"missing clock evidence did not quarantine current phase");
    out=zhangInspectSlipEvidence(state,synthetic(31),peers(870));
    check(out.usePhase && !out.breakArc,"temporary external-clock gap permanently poisoned clean tracking");
    state=smooth;
    for(int epoch=30;epoch<55;epoch++)
    {
        const auto* reference=zhangSlipClockReference(state);
        out=zhangInspectSlipEvidence(state,synthetic(epoch,{1,0,0}),peers(reference ? reference->time : 870));
        check(out.breakArc || !out.usePhase || epoch>30,
            "first strong injection silently accepted");
        if(out.breakArc)break;
        check(!out.usePhase,"unresolved slip silently released as extrapolation variance grew");
    }
    state=smooth;
    auto missingAuxiliary=synthetic(30);missingAuxiliary.third=false;
    out=zhangInspectSlipEvidence(state,missingAuxiliary,peers(870));
    check(out.usePhase && !out.breakArc && state.accepted.size()==31,
        "auxiliary loss cleared primary history or retired a primary arc");
    auto reacquired=synthetic(31);reacquired.thirdSignature=99;
    out=zhangInspectSlipEvidence(state,reacquired,peers(900));
    check(out.usePhase && !out.breakArc,"new auxiliary code inherited a false third-frequency step");
    state=smooth;
    auto codeLoss=synthetic(30);codeLoss.code=false;
    out=zhangInspectSlipEvidence(state,codeLoss,peers(870));
    check(out.usePhase && !out.breakArc && state.accepted.size()==31,"code absence cleared phase history");
    state=smooth;
    auto invalid=synthetic(30);invalid.covariance[0]=std::numeric_limits<double>::quiet_NaN();
    out=zhangInspectSlipEvidence(state,invalid,peers(870));
    check(!out.usePhase && !out.breakArc,"nonfinite covariance accepted");
    const auto sample=synthetic(30);
    auto reversePeers=peers(870);std::reverse(reversePeers.begin(),reversePeers.end());
    auto ordered=smooth,reversed=smooth;
    const auto forward=zhangInspectSlipEvidence(ordered,sample,peers(870));
    const auto backward=zhangInspectSlipEvidence(reversed,sample,reversePeers);
    check(forward.phaseScore==backward.phaseScore && forward.action==backward.action,"peer order changed decision");
    std::stringstream archive;
    {boost::archive::binary_oarchive writer(archive);writer << smooth;}
    ZhangSlipEvidenceState restored;
    {boost::archive::binary_iarchive reader(archive);reader >> restored;}
    auto a=zhangInspectSlipEvidence(smooth,synthetic(30,{1,1,0}),peers(870));
    auto b=zhangInspectSlipEvidence(restored,synthetic(30,{1,1,0}),peers(870));
    check(a.action==b.action && a.integer==b.integer && a.phaseScore==b.phaseScore,"restart changed decision");
    ZhangSlipEvidenceState bounded;
    for(int epoch=0;epoch<100;epoch++)
        zhangInspectSlipEvidence(bounded,synthetic(epoch),peers(30*(epoch-1)));
    check(bounded.accepted.size()==32 && bounded.screening.size()==96,"long history exceeded or failed to fill bound");
    std::stringstream fullArchive;
    {boost::archive::binary_oarchive writer(fullArchive);writer<<bounded;}
    ZhangSlipEvidenceState fullRestored;
    {boost::archive::binary_iarchive reader(fullArchive);reader>>fullRestored;}
    check(fullRestored.screening.size()==96 && fullRestored.accepted.size()==32,"full-window restart lost history");
    LegacySlipState legacy;LegacySlipSample legacySample;legacySample.time=30;legacySample.rawIf=12.5;
    legacy.accepted.push_back(legacySample);legacy.screening=legacy.accepted;
    std::stringstream legacyArchive;
    {boost::archive::binary_oarchive writer(legacyArchive);writer<<legacy;}
    ZhangSlipEvidenceState migrated;
    {boost::archive::binary_iarchive reader(legacyArchive);reader>>migrated;}
    check(migrated.accepted.size()==1 && migrated.accepted[0].rawIf==12.5 &&
        !migrated.accepted[0].gfRateValid,"old sample archive invented Doppler or failed migration");
    ZhangSlipEvidenceState dopplerState;
    {
        std::vector<ZhangSlipSample> past;
        for(int k=0;k<32;k++){auto sample=synthetic(k);sample.value[0]=2;sample.rawIf=2;past.push_back(sample);}
        const std::vector<double> times{960,990,1020,1050};
        const auto block=zhangSlipPhaseForecast(past,times,0);
        const double variance=past.front().covariance[0];
        std::array<Eigen::MatrixXd,3> reference;
        std::array<double,3> probability{};double normalization=0;
        for(int degree=0;degree<=2;degree++)
        {
            Eigen::MatrixXd design(32,degree+1),future(4,degree+1);
            for(int k=0;k<32;k++)
            {const double t=k-31.;design(k,0)=1;for(int j=1;j<=degree;j++)design(k,j)=design(k,j-1)*t;}
            for(int k=0;k<4;k++)
            {const double t=(times[k]-930)/30;future(k,0)=1;for(int j=1;j<=degree;j++)future(k,j)=future(k,j-1)*t;}
            // Independent QR pseudoinverse, not the forecast's normal solve.
            reference[degree]=future*(design/std::sqrt(variance)).colPivHouseholderQr().solve(
                Eigen::MatrixXd::Identity(32,32))/std::sqrt(variance);
            normalization+=probability[degree]=std::pow(32.,-.5*(degree+1));
        }
        Eigen::MatrixXd mean=Eigen::MatrixXd::Zero(4,32),cov=Eigen::MatrixXd::Zero(4,4);
        for(int d=0;d<3;d++){probability[d]/=normalization;mean+=probability[d]*reference[d];cov+=probability[d]*variance*reference[d]*reference[d].transpose();}
        cov=2*(cov-variance*mean*mean.transpose());
        check((block.weights-mean).norm()<1e-10 && (block.extraCovariance-cov).norm()<1e-12,
            "joint future forecast disagrees with independent QR covariance reference");
        check(std::abs(block.extraCovariance(0,1))>1e-12,"future forecast silently discarded time covariance");
    }
    auto dopplerSample=[](int epoch)
    {
        auto value=synthetic(epoch);value.gfRateValid=true;
        value.gfRate=(std::pow(wavelength[1]/wavelength[0],2)-1)*(.02+.002*epoch)/30;
        return value;
    };
    for(int epoch=0;epoch<40;epoch++)
        check(zhangInspectSlipEvidence(dopplerState,dopplerSample(epoch),peers(30*(epoch-1))).usePhase,
            "smooth ionosphere with consistent Doppler rejected");
    {
        auto referenceState=dopplerState;referenceState.accepted[30].value[1]+=.003;
        auto current=dopplerSample(40);current.value[1]+=.005;
        const auto evaluated=zhangSlipDopplerEvidence(referenceState,current);
        auto observations=referenceState.accepted;observations.push_back(current);
        const int n=observations.size(),m=n-1;
        Eigen::MatrixXd transform=Eigen::MatrixXd::Zero(m,2*n);
        Eigen::VectorXd data(2*n),rawVariance(2*n);
        const double rawPair=observations[0].covariance[5]*2;
        const double calibrated=(rawPair+.008*.008)*(1+1.570796326794897/(n-2));
        const double rateFactor=(calibrated-rawPair-.005*.005)/(225*rawPair);
        for(int k=0;k<n;k++)
        {
            data[k]=observations[k].value[1];data[n+k]=observations[k].gfRate;
            rawVariance[k]=observations[k].covariance[5];rawVariance[n+k]=rateFactor*rawVariance[k];
            if(k)
            {transform(k-1,k)=1;transform(k-1,k-1)=-1;transform(k-1,n+k)=-15;transform(k-1,n+k-1)=-15;}
        }
        // Independent dense construction from raw variables. It deliberately
        // does NOT call the tridiagonal recursion or its covariance builder.
        const Eigen::VectorXd error=transform*data;
        const Eigen::MatrixXd covariance=transform*rawVariance.asDiagonal()*transform.transpose()+
            .005*.005*Eigen::MatrixXd::Identity(m,m);
        const auto prior=covariance.topLeftCorner(m-1,m-1).ldlt();
        const Eigen::VectorXd lastCross=covariance.col(m-1).head(m-1);
        const Eigen::VectorXd regression=prior.solve(lastCross);
        const double expected=error[m-1]-regression.dot(error.head(m-1));
        const double expectedVariance=covariance(m-1,m-1)-regression.dot(lastCross);
        check(std::abs(evaluated.innovation-expected)<1e-10 &&
            std::abs(evaluated.variance-expectedVariance)<1e-12,
            "Doppler whitening disagrees with independent raw-variable covariance construction");
        Eigen::VectorXd contrast=Eigen::VectorXd::Zero(m);contrast[m-1]=1;contrast.head(m-1)=-regression;
        const Eigen::VectorXd phaseCoefficients=transform.leftCols(n).transpose()*contrast;
        double actualCross=current.covariance[1],expectedCross=current.covariance[1];
        for(int k=0;k<n-1;k++)
        {
            expectedCross-=phaseCoefficients[k]*observations[k].covariance[1]/(n-1);
            for(const auto& [time,cross]:evaluated.pastCross)
                if(time==observations[k].time)actualCross+=cross[0]/(n-1);
        }
        check(std::abs(actualCross-expectedCross)<1e-12,"Doppler/IF shared-phase covariance sign wrong");
    }
    auto blunder=dopplerSample(40);blunder.gfRate+=.05;
    const auto blunderResult=zhangInspectSlipEvidence(dopplerState,blunder,peers(1170));
    check(blunderResult.usePhase && !blunderResult.breakArc,"isolated Doppler blunder retired/quarantined clean phase");
    auto lostDoppler=dopplerSample(41);lostDoppler.gfRateValid=false;
    check(zhangInspectSlipEvidence(dopplerState,lostDoppler,peers(1200)).usePhase,"Doppler loss falsely broke phase continuity");
    {
        std::stringstream archive;
        {boost::archive::binary_oarchive writer(archive);writer<<dopplerState;}
        bool refused=false;
        try{LegacySlipState incompatible;boost::archive::binary_iarchive reader(archive);reader>>incompatible;}
        catch(const boost::archive::archive_exception&){refused=true;}
        check(refused,"old sample reader silently accepted a new Doppler layout");
    }
    ZhangSignalTracking tracking;
    check(!zhangDecidePhaseTracking(tracking,0,true,1,0,45).breakArc,"first tracking sample broke arc");
    check(zhangDecidePhaseTracking(tracking,1,true,1,30,45).breakArc,"LLI missed");
    check(zhangDecidePhaseTracking(tracking,2,true,1,60,45).breakArc,"half cycle missed");
    check(!zhangDecidePhaseTracking(tracking,2,true,1,90,45).useCurrentPhase,"half cycle used");
    check(zhangDecidePhaseTracking(tracking,0,true,1,120,45).breakArc,"half cycle recovery inherited arc");
    check(zhangDecidePhaseTracking(tracking,0,true,1,210,45).breakArc,"gap missed");
    std::cout << "cycle_slip_unit=PASS\n";
}

static int monteCarlo(int seeds,int offset=0,bool inject=false)
{
    size_t tested=0,quarantine=0,confirmed=0,detected=0,misidentified=0,newQuarantines=0,wrongMasks=0,retiredCases=0,identifiedCases=0;
    std::vector<int> detectionDelays;
    std::vector<int> missedSeeds;
    std::vector<int> eventExcludedSeeds;
    std::map<std::string,size_t> firstActions;
    std::array<size_t,4> patternCases{},patternDetections{};
    size_t instantCases=0,within30=0,within60=0,preExcludedCases=0;
    for(int seed=0;seed<seeds;seed++)
    {
        std::mt19937 random(seed+offset+7319);std::normal_distribution<double> normal;
        std::array<ZhangSlipEvidenceState,6> states;
        const std::array<std::array<int,3>,4> patterns{{{1,0,0},{1,1,0},{-1,-1,-1},{0,0,1}}};
        const int pattern=(seed+offset)%patterns.size();
        const auto step=patterns[pattern];bool seen=false,certified=false,retired=false;
        for(int epoch=0;epoch<120;epoch++)
        {
            std::array<ZhangSlipSample,6> samples;
            std::array<double,6> gfRates{};
            for(int sat=0;sat<6;sat++)
            {
                double iono=3+sat+.01*epoch+.0003*epoch*epoch;
                if(std::getenv("ZHANG_SLIP_TEST_IONO_STRESS"))
                    iono+=(1+.2*sat)*(.5*std::sin(.11*epoch)+.7*std::tanh((epoch-60.)/3));
                const double clock=epoch>=60 ? 300000 : 0;
                std::array<double,3> phase;
                const double sigma=.002+.001*sat;
                const bool phaseTrend=std::getenv("ZHANG_SLIP_TEST_PHASE_TREND");
                const double geometryError=phaseTrend ? .03*std::sin(.07*epoch+sat)+
                    .05*std::tanh((epoch-60.)/8)*(sat+1)/6 : 0;
                const double windup=phaseTrend ? .1*std::sin(.05*epoch+sat) : 0;
                for(int f=0;f<3;f++)phase[f]=1000+clock+sat*10+geometryError+wavelength[f]*windup-
                    std::pow(wavelength[f]/wavelength[0],2)*iono+sigma*normal(random)+
                    (inject && sat==5 && epoch>=60 ? wavelength[f]*step[f] : 0);
                samples[sat]=zhangMakeSlipSample(epoch*30,123,wavelength,phase,
                    {sigma*sigma,sigma*sigma,sigma*sigma},
                    {1000+clock+geometryError+iono+200*sigma*normal(random),
                     1000+clock+geometryError+std::pow(wavelength[1]/wavelength[0],2)*iono+200*sigma*normal(random)},
                    {40000*sigma*sigma,40000*sigma*sigma},true,true,1000,true);
            }
            const bool useDoppler=std::getenv("ZHANG_SLIP_TEST_DOPPLER");
            if(useDoppler)
                for(int sat=0;sat<6;sat++)
                {
                    double derivative=.01+.0006*epoch;
                    if(std::getenv("ZHANG_SLIP_TEST_IONO_STRESS"))
                        derivative+=(1+.2*sat)*(.055*std::cos(.11*epoch)+(.7/3)/std::pow(std::cosh((epoch-60.)/3),2));
                    gfRates[sat]=(std::pow(wavelength[1]/wavelength[0],2)-1)*derivative/30+
                        .0005*(.002+.001*sat)/.007*normal(random);
                    if(std::getenv("ZHANG_SLIP_TEST_PHASE_TREND"))
                        gfRates[sat]+=(wavelength[0]-wavelength[1])*.005*std::cos(.05*epoch+sat)/30;
                    if(std::getenv("ZHANG_SLIP_TEST_DOPPLER_BLUNDER") && sat==5 && epoch==60)gfRates[sat]+=.05;
                    samples[sat].gfRate=gfRates[sat];samples[sat].gfRateValid=true;
                }
            const auto before=states;
            for(int sat=0;sat<6;sat++)
            {
                std::vector<ZhangSlipClockPeer> clock;
                const auto* reference=zhangSlipClockReference(before[sat]);
                for(int peer=0;peer<6;peer++)
                    if(reference && peer!=sat && !before[peer].screening.empty())
                    {
                        const auto& history=before[peer].screening;
                        const auto previous=std::find_if(history.begin(),history.end(),[&](const auto& old)
                            {return old.time==reference->time && old.geometry;});
                        if(previous==history.end())continue;
                        clock.push_back({previous->time,samples[peer].rawIf-previous->rawIf,
                            samples[peer].covariance[0]+previous->covariance[0],history,samples[peer],true});
                    }
                ZhangSlipTrendDiagnostic trend;
                auto out=zhangInspectSlipEvidence(states[sat],samples[sat],clock,false,{},&trend);
                if(sat==5 && epoch==60 && (!out.usePhase || out.breakArc))eventExcludedSeeds.push_back(seed+offset);
                if(inject && sat==5 && epoch==59 && !out.usePhase)preExcludedCases++;
                if(std::getenv("ZHANG_SLIP_TEST_TRACE") && sat==std::atoi(std::getenv("ZHANG_SLIP_TEST_TRACE")))
                {
                    auto raw=before[sat].screening;raw.push_back(samples[sat]);
                    Eigen::Vector4d step;Eigen::Matrix4d covariance;double score=0,time=0;
                    std::array<double,4> diagnostic{};
                    zhangSlipScreenCommonMode(raw,clock,true,step,covariance,score,time,&diagnostic);
                    std::cerr<<"QC_TEST epoch="<<epoch<<" accepted="<<before[sat].accepted.size()
                        <<" last="<<(before[sat].accepted.empty() ? -1 : before[sat].accepted.back().time)
                        <<" action="<<zhangSlipActionName(out.action)<<" phase_score="<<out.phaseScore
                        <<" change_score="<<out.changePointScore<<" best="<<out.bestScore
                        <<" doppler="<<out.dopplerAvailable<<" de="<<out.dopplerInnovation
                        <<" ds="<<out.dopplerSigma<<" dscore="<<out.dopplerScore<<" dfit="<<out.dopplerIntegerFit;
                    std::cerr<<" early_max="<<diagnostic[0]<<" early_post="<<diagnostic[1]
                        <<" early_info="<<diagnostic[2]<<" early_first_if="<<diagnostic[3];
                    std::cerr<<" early_stage="<<trend.stage<<" trend_qualified="<<trend.qualified
                        <<" trend_window="<<trend.window<<" trend_score="<<trend.score
                        <<" trend_support="<<trend.support<<" trend_fit="<<trend.fit<<" trend_assisted="<<trend.assisted;
                    for(int component=0;component<4;component++)
                        std::cerr<<" e"<<component<<'='<<out.innovation[component]<<" s"<<component<<'='<<out.sigma[component];
                    std::cerr<<'\n';
                }
                if(epoch>=8){tested++;quarantine+=!out.usePhase;confirmed+=out.breakArc;}
                // A certificate on a no-slip link is itself a wrong integer,
                // not merely a false-retirement counter. Previously this
                // field only covered misidentification at the injected link.
                if(epoch>=8 && out.action==ZhangSlipAction::INTEGER_STEP &&
                    (!inject || sat!=5 || epoch<60))misidentified++;
                if(!inject && epoch>=8 && out.breakArc)
                    std::cerr<<"FALSE_RETIREMENT seed_offset="<<seed+offset<<" station_link="<<sat
                        <<" epoch="<<epoch<<" action="<<zhangSlipActionName(out.action)
                        <<" phase_score="<<out.phaseScore<<" change_score="<<out.changePointScore
                        <<" best_score="<<out.bestScore<<'\n';
                if(epoch>=8 && !before[sat].phaseSuspect && states[sat].phaseSuspect)newQuarantines++;
                if(inject && sat==5 && epoch>=60)
                {
                    if(!seen && (!out.usePhase || out.breakArc))
                    {
                        detectionDelays.push_back((epoch-60)*30);firstActions[zhangSlipActionName(out.action)]++;
                        instantCases+=epoch==60;within30+=epoch<=61;within60+=epoch<=62;
                        if(epoch>62)std::cerr<<"DELAYED_DETECTION seed_offset="<<seed+offset
                            <<" delay_s="<<(epoch-60)*30<<" phase_score="<<out.phaseScore
                            <<" change_score="<<out.changePointScore<<'\n';
                    }
                    seen |= !out.usePhase || out.breakArc;
                    if(!certified && out.action==ZhangSlipAction::INTEGER_STEP)
                    {misidentified+=out.integer!=step;certified=true;identifiedCases++;}
                    if(!retired && out.breakArc)
                    {
                        bool missing=false;
                        for(int f=0;f<3;f++)missing |= step[f]!=0 && !out.retire[f];
                        wrongMasks+=missing;retired=true;retiredCases++;
                    }
                }
            }
        }
        detected+=seen;
        if(inject){patternCases[pattern]++;patternDetections[pattern]+=seen;if(!seen)missedSeeds.push_back(seed+offset);}
    }
    std::cout<<"{\"seeds\":"<<seeds<<",\"offset\":"<<offset<<",\"tested_link_epochs\":"<<tested
        <<",\"quarantines\":"<<quarantine<<",\"retirements\":"<<confirmed
        <<",\"new_quarantine_events\":"<<newQuarantines
        <<",\"injected_cases\":"<<(inject ? seeds : 0)<<",\"detected_cases\":"<<detected
        <<",\"wrong_integer_certificates\":"<<misidentified<<",\"wrong_retirement_masks\":"<<wrongMasks
        <<",\"retired_cases\":"<<retiredCases<<",\"identified_cases\":"<<identifiedCases;
    std::cout<<",\"current_epoch_detected_cases\":"<<instantCases<<",\"within_30s_cases\":"<<within30
        <<",\"within_60s_cases\":"<<within60<<",\"already_excluded_before_injection\":"<<preExcludedCases;
    std::sort(detectionDelays.begin(),detectionDelays.end());
    std::cout<<",\"max_detection_delay_s\":"<<(detectionDelays.empty() ? -1 : detectionDelays.back())
        <<",\"median_detection_delay_s\":"<<(detectionDelays.empty() ? -1 : detectionDelays[detectionDelays.size()/2])
        <<",\"p95_detection_delay_s\":"<<(detectionDelays.empty() ? -1 : detectionDelays[(detectionDelays.size()-1)*95/100])
        <<",\"first_detection_actions\":{";
    bool firstAction=true;
    for(const auto& [action,count]:firstActions)
    {if(!firstAction)std::cout<<',';firstAction=false;std::cout<<'\"'<<action<<"\":"<<count;}
    std::cout<<'}'
        <<",\"patterns\":[";
    for(int i=0;i<4;i++)
    {if(i)std::cout<<',';std::cout<<"{\"cases\":"<<patternCases[i]<<",\"detected\":"<<patternDetections[i]<<'}';}
    std::cout<<"],\"missed_seed_offsets\":[";
    for(size_t i=0;i<missedSeeds.size();i++){if(i)std::cout<<',';std::cout<<missedSeeds[i];}
    std::cout<<"],\"event_epoch_excluded_seed_offsets\":[";
    for(size_t i=0;i<eventExcludedSeeds.size();i++){if(i)std::cout<<',';std::cout<<eventExcludedSeeds[i];}
    std::cout<<"]}\n";
    return inject ? (misidentified || wrongMasks || detected!=(size_t)seeds) :
        ((misidentified || confirmed || (tested && double(newQuarantines)/tested>1e-4)) ? 1 : 0);
}

// Fixture columns are raw observations, geometry and variances, not precomputed
// expected decisions. Every scenario starts from its own empty detector states.
static int replay(const std::string& filename,int n1,int n2,int n3,double injection)
{
    struct Row {int satellite;ZhangSlipSample sample;bool hard;};
    std::ifstream input(filename);check(bool(input),"fixture missing");
    std::map<double,std::vector<Row>> epochs;
    std::string line;
    while(std::getline(input,line))
    {
        if(line.empty() || line[0]=='#')continue;
        std::istringstream row(line);
        double time,range,el,p1,p2,l1,l2,l5,phaseVar,codeVar;int satellite,lli1,lli2,lli5;
        check(bool(row>>time>>satellite>>range>>el>>l1>>l2>>l5>>p1>>p2>>phaseVar>>codeVar>>lli1>>lli2>>lli5),"bad fixture row");
        if(el<10)continue;
        if(satellite==25 && time>=injection) {l1+=n1;l2+=n2;if(l5!=0)l5+=n3;}
        const bool third=l5!=0 && !(lli5&3);
        auto sample=zhangMakeSlipSample(time,123,wavelength,
            {l1*wavelength[0],l2*wavelength[1],l5*wavelength[2]},
            {phaseVar,phaseVar,third ? phaseVar : 0},{p1,p2},{codeVar,codeVar},third,true,range,true);
        sample.thirdSignature=third ? 1 : 0;
        double d1=0,d2=0;
        if(row>>d1>>d2)
        {sample.gfRateValid=d1!=0 && d2!=0 && !std::getenv("ZHANG_SLIP_TEST_NO_DOPPLER");
         sample.gfRate=-wavelength[0]*d1+wavelength[1]*d2;}
        epochs[time].push_back({satellite,sample,bool((lli1|lli2)&3)});
    }
    std::map<int,ZhangSlipEvidenceState> states;
    std::cout << "time,satellite,action,use_phase,break_arc,clock_partners,if,gf,mw,triple,sigma_if,sigma_gf,sigma_mw,sigma_triple,phase_score,change_score,change_time,ratio,n1,n2,n3,doppler_available,doppler_assisted,doppler_error,doppler_sigma,doppler_score\n";
    for(const auto& [time,rows]:epochs)
    {
        const auto before=states;
        for(const auto& row:rows)
        {
            std::vector<ZhangSlipClockPeer> clock;
            for(const auto& peer:rows)
            {
                if(peer.satellite==row.satellite || peer.hard || !before.count(peer.satellite))continue;
                const auto* reference=before.count(row.satellite) ? zhangSlipClockReference(before.at(row.satellite)) : nullptr;
                if(!reference)continue;
                const auto& history=before.at(peer.satellite).screening;
                const auto previous=std::find_if(history.begin(),history.end(),[&](const auto& old)
                    {return old.time==reference->time && old.geometry && old.signature==peer.sample.signature;});
                if(previous==history.end())continue;
                clock.push_back({previous->time,peer.sample.rawIf-previous->rawIf,
                    peer.sample.covariance[0]+previous->covariance[0],history,peer.sample,true});
            }
            auto out=zhangInspectSlipEvidence(states[row.satellite],row.sample,clock,row.hard);
            std::cout << std::setprecision(15) << time << ',' << row.satellite << ',' << zhangSlipActionName(out.action)
                << ',' << out.usePhase << ',' << out.breakArc << ',' << out.clockPartners;
            for(auto value:out.innovation)std::cout << ',' << value;
            for(auto value:out.sigma)std::cout << ',' << value;
            std::cout << ',' << out.phaseScore << ',' << out.changePointScore << ',' << out.changePointTime << ',' << out.ratio;
            for(auto value:out.integer)std::cout << ',' << value;
            std::cout<<','<<out.dopplerAvailable<<','<<out.dopplerAssisted<<','<<out.dopplerInnovation
                <<','<<out.dopplerSigma<<','<<out.dopplerScore;
            std::cout << '\n';
        }
    }
    return 0;
}
int main(int argc,char** argv)
{
    try
    {
        if(argc==3 && std::string(argv[1])=="--monte-carlo") return monteCarlo(std::stoi(argv[2]));
        if(argc==4 && std::string(argv[1])=="--monte-carlo") return monteCarlo(std::stoi(argv[2]),std::stoi(argv[3]));
        if(argc==4 && std::string(argv[1])=="--monte-carlo-slip") return monteCarlo(std::stoi(argv[2]),std::stoi(argv[3]),true);
        if(argc==7 && std::string(argv[1])=="--replay")
            return replay(argv[2],std::stoi(argv[3]),std::stoi(argv[4]),std::stoi(argv[5]),std::stod(argv[6]));
        unitTests();return 0;
    }
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
