// Real production kFilter, with an independent information-form oracle.
// This tests the numerical kernel without physical graph augmentation.
#include "common/algebra.hpp"
#include "common/zhangFloatThreading.hpp"
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

static void check(bool ok,const char* why)
{ if(!ok) throw std::runtime_error(why); }

static void fixture(KFState& state,KFMeas& meas,int n=7,int m=11)
{
    state.x=VectorXd::LinSpaced(n,-0.3,0.4);
    state.dx=VectorXd::Zero(n);
    VectorXd u=VectorXd::LinSpaced(n,0.01,0.09);
    state.P=MatrixXd::Identity(n,n)+u*u.transpose();
    meas.H.resize(m,n);
    for(int j=0;j<n;++j) for(int i=0;i<m;++i)
        meas.H(i,j)=0.2*std::sin(0.31*(i+1)*(j+2));
    VectorXd r=VectorXd::LinSpaced(m,0.01,0.08);
    meas.R=MatrixXd::Identity(m,m)*0.7+r*r.transpose();
    meas.V=VectorXd::LinSpaced(m,-0.03,0.02);
    meas.Y=meas.V;meas.VV=meas.V;
    meas.H_star=MatrixXd::Zero(m,0);meas.uncorrelatedNoise.resize(0);
    state.advanced_postfits=false;state.postfitOpts.omega_test=false;
}

static void independentPosterior()
{
    for(const bool joseph:{false,true}) for(const bool symmetric:{false,true})
    {
        setenv("ZHANG_R51_FLOAT_MATRIX_MODE",symmetric?"symmetric":"legacy",1);
        ZhangFloatThreading::ScopedKernelContext floatContext(true);
        KFState state;KFMeas meas;fixture(state,meas);
        state.joseph_stabilisation=joseph;
        const MatrixXd prior=state.P;const VectorXd priorMean=state.x;
        const MatrixXd information=prior.inverse()+
            meas.H.transpose()*meas.R.ldlt().solve(meas.H);
        const MatrixXd expectedP=information.ldlt().solve(MatrixXd::Identity(7,7));
        const VectorXd expectedX=priorMean+expectedP*meas.H.transpose()*meas.R.ldlt().solve(meas.V);
        VectorXd xp=state.x,dx=state.dx;MatrixXd pp=state.P,qi,qh;
        std::ostringstream directTrace;
        check(state.kFilter(directTrace,meas,xp,pp,dx,qi,qh,0,7,0,11),"direct full-state update failed");
        check(directTrace.str().find("direct_full=1")!=std::string::npos,"complete nonalias destination was not used");
        check(directTrace.str().find(std::string("ZHANG_FLOAT_COVARIANCE_MODE symmetric=")+
            (symmetric && !joseph ? "1" : "0"))!=std::string::npos,
            "requested symmetric/Joseph kernel coverage was not actually executed");
        check((xp-expectedX).cwiseAbs().maxCoeff()<3e-12 &&
            (pp-expectedP).cwiseAbs().maxCoeff()<3e-12,"independent information-form posterior differs");
        check((state.P-prior).isZero() && (state.x-priorMean).isZero(),"nonalias update changed authoritative prior");
        // A separate aliasing call is forced through the original staged
        // implementation. It does not share the direct output's matrices.
        KFState alias;KFMeas aliasMeas;fixture(alias,aliasMeas);alias.joseph_stabilisation=joseph;
        VectorXd ax=alias.x,ad=alias.dx;MatrixXd ai,ah;std::ostringstream aliasTrace;
        check(alias.kFilter(aliasTrace,aliasMeas,ax,alias.P,ad,ai,ah,0,7,0,11),"staged alias update failed");
        check(aliasTrace.str().find("direct_full=0")!=std::string::npos,"aliasing prior entered direct path");
        check((ax-xp).isZero() && (alias.P-pp).isZero() && (ad-dx).isZero(),"direct/staged arithmetic changed");
    }
}

static void partialAndFailure()
{
    setenv("ZHANG_R51_FLOAT_MATRIX_MODE","legacy",1);
    ZhangFloatThreading::ScopedKernelContext floatContext(true);
    KFState state;KFMeas meas;fixture(state,meas,7,9);
    VectorXd xp=state.x,dx=state.dx;MatrixXd pp=state.P,qi,qh;
    const MatrixXd prior=state.P;
    const MatrixXd p=prior.block(1,1,4,4),h=meas.H.block(2,1,4,4),r=meas.R.block(2,2,4,4);
    const MatrixXd gain=p*h.transpose()*(h*p*h.transpose()+r).ldlt().solve(MatrixXd::Identity(4,4));
    MatrixXd expected=prior;expected.block(1,1,4,4)=p-gain*h*p;
    std::ostringstream trace;
    check(state.kFilter(trace,meas,xp,pp,dx,qi,qh,1,4,2,4),"partial update failed");
    check(trace.str().find("direct_full=0")!=std::string::npos,"partial chunk entered full destination path");
    check((pp-expected).cwiseAbs().maxCoeff()<3e-12,"partial covariance/cross-blocks changed");
    fixture(state,meas,1,1);state.x(0)=1.7e308;state.P.setIdentity();
    meas.H.setIdentity();meas.R.setIdentity();meas.V(0)=1.7e308;
    xp=state.x;dx=state.dx;pp=state.P;std::ostringstream rejected;
    check(!state.kFilter(rejected,meas,xp,pp,dx,qi,qh,0,1,0,1,true),"nonfinite posterior accepted");
    check((xp-state.x).isZero() && (pp-state.P).isZero() && dx.isZero(),"direct destination did not preserve failure reset");
    check(rejected.str().find("NONFINITE KALMAN POSTERIOR")!=std::string::npos,"failure gate did not execute");
}

static void symmetricMemoryOracle()
{
    const Eigen::Index n=513,leading=n+7;
    for(int threads:{1,4})
    {
        omp_set_num_threads(threads);
        MatrixXd storage=MatrixXd::Constant(leading,n,-42.0);
        for(Eigen::Index j=0;j<n;++j) for(Eigen::Index i=0;i<n;++i)
            storage(i,j)=std::sin(0.01*(i+2)*(j+1))+0.001*(i-j);
        const MatrixXd before=storage.topRows(n);
        const MatrixXd expected=(before+before.transpose()).eval()/2;
        ZhangFloatThreading::ScopedKernelContext floatContext(true);
        ZhangFloatThreading::averageSymmetricInPlace(storage.data(),leading,n);
        check((storage.topRows(n)-expected).isZero(),"parallel pairwise symmetry differs from independent expression");
        check((storage.bottomRows(7).array()==-42.0).all(),"symmetry helper touched padded storage");
    }
    double large=std::numeric_limits<double>::max();
    ZhangFloatThreading::averageSymmetricInPlace(&large,1,1);
    check(!std::isfinite(large),"diagonal overflow behavior was silently changed");
}

int main()
{
    independentPosterior();partialAndFailure();symmetricMemoryOracle();
    std::cout<<"FLOAT_POSTERIOR_DESTINATION_CONTRACT_PASS independent_information=1 "
        "joseph_standard_symmetric=1 alias_staging_exact=1 partial_crossQ=1 "
        "nonfinite_reset=1 parallel_padded_symmetry=1\n";
}
