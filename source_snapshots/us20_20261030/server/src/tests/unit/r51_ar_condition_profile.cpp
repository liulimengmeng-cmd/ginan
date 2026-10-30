#include "common/zhangIntegerConditioner.hpp"
#include "common/blasThreading.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>
#include <omp.h>

int main(int argc,char** argv)
{
    const int threads=argc>1?std::stoi(argv[1]):1;
    const int n=argc>2?std::stoi(argv[2]):1175;
    const int steps=3,block=32;
    ZhangArThreading::active=argc>3 && std::stoi(argv[3])!=0;
    omp_set_num_threads(12);Eigen::setNbThreads(1);
    BlasThreading::openblas_set_num_threads(threads);
    MatrixXd l(n,17);
    VectorXd initial(n);
    for(int i=0;i<n;++i) {
        initial(i)=.02*std::sin(.13*i);
        for(int j=0;j<17;++j)l(i,j)=.03*std::sin(.11*i+.17*j);
    }
    const MatrixXd prior=l*l.transpose()+MatrixXd::Identity(n,n);
    MatrixXd covariance=prior;VectorXd mean=initial;
    const auto started=std::chrono::steady_clock::now();
    for(int step=0;step<steps;++step) {
        MatrixXd a=MatrixXd::Zero(block,n);
        for(int i=0;i<block;++i)a(i,step*block+i)=1;
        const auto conditioned=zhangConditionPosteriorEffectiveIntegers(mean,covariance,a,VectorXd::Zero(block));
        if(!conditioned.valid)throw std::runtime_error(conditioned.failureReason);
        mean=conditioned.mean;covariance=conditioned.covariance;
    }
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    double meanError,covarianceError;
    {
        BlasThreading::ScopedOpenBlasThreadLimit serial(1);
        MatrixXd a=MatrixXd::Zero(steps*block,n);
        a.leftCols(steps*block).setIdentity();
        MatrixXd cross=prior*a.transpose();
        MatrixXd s=a*cross;
        VectorXd expectedMean=initial-cross*s.ldlt().solve(a*initial);
        MatrixXd expectedCovariance=prior-cross*s.ldlt().solve(cross.transpose());
        meanError=(mean-expectedMean).cwiseAbs().maxCoeff();
        covarianceError=(covariance-expectedCovariance).cwiseAbs().maxCoeff();
    }
    if(meanError>2e-11 || covarianceError>2e-11)
        throw std::runtime_error("full mean/covariance independent oracle mismatch");
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(15)<<"CONDITION threads="<<threads
        <<" adaptive="<<ZhangArThreading::active<<" dimension="<<n
        <<" steps="<<steps<<" block="<<block<<" seconds="<<elapsed
        <<" mean_error="<<meanError<<" covariance_error="<<covarianceError
        <<" peak_rss_kib="<<usage.ru_maxrss<<'\n';
}
