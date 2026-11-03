// Paired old/new production-kernel benchmark, not a real network epoch.
// Compile the same driver with each producer's own complete headers/objects.
#include "common/algebra.hpp"
#include "common/zhangFloatThreading.hpp"
#include <openssl/sha.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void check(bool ok,const char* why)
{ if(!ok) throw std::runtime_error(why); }
static std::string digest(const double* data,std::size_t count)
{
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(data),count*sizeof(double),hash);
    std::ostringstream out;
    for(unsigned char value:hash) out<<std::hex<<std::setw(2)<<std::setfill('0')<<int(value);
    return out.str();
}
int main()
{
    constexpr int n=1536,m=768,repeats=3;
    setenv("ZHANG_R51_FLOAT_THREADS","12",1);
    setenv("ZHANG_FLOAT_MATRIX_THREADS","6",1);
    setenv("ZHANG_FLOAT_MEMORY_MODE","parallel",1);
    setenv("ZHANG_R51_FLOAT_MATRIX_MODE","legacy",1);
    ZhangFloatThreading::ScopedFloatThreads threads;
    KFState state;KFMeas measurement;
    state.x=VectorXd::LinSpaced(n,-.3,.4);state.dx=VectorXd::Zero(n);
    const VectorXd u=VectorXd::LinSpaced(n,.001,.009);
    state.P=MatrixXd::Identity(n,n)+u*u.transpose();
    state.advanced_postfits=false;state.postfitOpts.omega_test=false;
    measurement.H.resize(m,n);
    for(int j=0;j<n;++j) for(int i=0;i<m;++i)
        measurement.H(i,j)=.003*std::sin(.001*(i+1)*(j+2));
    measurement.R=MatrixXd::Identity(m,m)*.7;
    measurement.V=VectorXd::LinSpaced(m,-.03,.02);
    measurement.Y=measurement.V;measurement.VV=measurement.V;
    measurement.H_star=MatrixXd::Zero(m,0);measurement.uncorrelatedNoise.resize(0);
    std::string xHash,pHash,dxHash;
    bool observedDirect=false;
    for(int repeat=0;repeat<repeats;++repeat)
    {
        VectorXd x=state.x,dx=state.dx;MatrixXd p=state.P,qi,qh;
        std::ostringstream trace;
        const auto start=std::chrono::steady_clock::now();
        check(state.kFilter(trace,measurement,x,p,dx,qi,qh,0,n,0,m),"benchmark kernel failed");
        const auto end=std::chrono::steady_clock::now();
        const auto nextX=digest(x.data(),x.size()),nextP=digest(p.data(),p.size()),
            nextDx=digest(dx.data(),dx.size());
        if(repeat) check(xHash==nextX && pHash==nextP && dxHash==nextDx,
            "identical frozen input gave different raw posterior bytes");
        xHash=nextX;pHash=nextP;dxHash=nextDx;
        observedDirect=trace.str().find("direct_full=1")!=std::string::npos;
        std::cout<<"FLOAT_BUFFER_BENCH_SAMPLE repeat="<<repeat
            <<" seconds="<<std::setprecision(17)<<std::chrono::duration<double>(end-start).count()
            <<" states="<<n<<" measurements="<<m<<" float_threads=12 matrix_threads=6\n";
    }
    std::cout<<"FLOAT_BUFFER_BENCH_PASS x_raw_sha256="<<xHash
        <<" P_raw_sha256="<<pHash<<" dx_raw_sha256="<<dxHash
        <<" observed_direct_destination="<<observedDirect
        <<" network_epoch_test=0 no_accuracy_claim=1\n";
}
