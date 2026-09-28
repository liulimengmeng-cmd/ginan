#pragma once

#include "common/eigenIncluder.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// Bound these memory-bandwidth-heavy, row-independent operations.  The
// environment override allows a controlled server comparison without making
// the number of OpenMP workers equal to the number of ambiguity rows.
inline int zhangR51ParallelThreads()
{
#ifdef _OPENMP
    const int maximum=std::max(1,omp_get_max_threads());
    int threads=std::min(maximum,8);
    if(const char* value=std::getenv("ZHANG_R51_AR_THREADS")) {
        char* end=nullptr;
        const long requested=std::strtol(value,&end,10);
        if(end!=value && *end=='\0' && requested>0)
            threads=static_cast<int>(std::min<long>(requested,maximum));
    }
    return threads;
#else
    return 1;
#endif
}

// A is immutable for the lifetime of this object.  Dense BLAS remains the
// fallback: traversing sparse supports is beneficial only for genuinely sparse
// AR equality rows.  Every worker owns complete output rows or distinct Gram
// entries; neither the posterior nor the input rows are modified.
struct ZhangParallelRowProjection
{
    const Eigen::MatrixXd& rows;
    std::vector<std::vector<std::pair<int,double>>> supports;
    std::size_t nonzeros=0;
    bool sparse=false;

    explicit ZhangParallelRowProjection(const Eigen::MatrixXd& a):rows(a)
    {
        supports.resize(rows.rows());
        for(int r=0;r<rows.rows();++r)
        for(int c=0;c<rows.cols();++c)if(rows(r,c)!=0) {
            supports[r].emplace_back(c,rows(r,c));++nonzeros;
        }
        sparse=rows.rows()>0 && rows.cols()>0 &&
            nonzeros*5<=static_cast<std::size_t>(rows.rows())*rows.cols();
    }

    Eigen::MatrixXd multiply(const Eigen::MatrixXd& right) const
    {
        Eigen::MatrixXd result(rows.rows(),right.cols());
        if(!sparse) {result.noalias()=rows*right;return result;}
        const int threads=zhangR51ParallelThreads();
        const bool parallel=threads>1 && rows.rows()>=32 && right.cols()>=128;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
#endif
        for(int r=0;r<rows.rows();++r) {
            result.row(r).setZero();
            for(const auto& [c,value]:supports[r])
                result.row(r)+=value*right.row(c);
        }
        return result;
    }

    Eigen::VectorXd multiply(const Eigen::VectorXd& right) const
    {
        if(!sparse)return rows*right;
        Eigen::VectorXd result(rows.rows());
        for(int r=0;r<rows.rows();++r) {
            double sum=0;
            for(const auto& [c,value]:supports[r])sum+=value*right(c);
            result(r)=sum;
        }
        return result;
    }

    // cross = A*P.  The exact upper triangle has one writer per (i,j).
    Eigen::MatrixXd covariance(const Eigen::MatrixXd& cross) const
    {
        Eigen::MatrixXd result(rows.rows(),rows.rows());
        if(!sparse) {
            result.noalias()=cross*rows.transpose();
            return (0.5*(result+result.transpose())).eval();
        }
        const int threads=zhangR51ParallelThreads();
        const bool parallel=threads>1 && rows.rows()>=64;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
#endif
        for(int i=0;i<rows.rows();++i)
        for(int j=i;j<rows.rows();++j) {
            double forward=0;
            for(const auto& [c,value]:supports[j])forward+=cross(i,c)*value;
            double sum=forward;
            if(i!=j) {
                double reverse=0;
                for(const auto& [c,value]:supports[i])reverse+=cross(j,c)*value;
                sum=0.5*(forward+reverse);
            }
            result(i,j)=sum;result(j,i)=sum;
        }
        return result;
    }
};

// Q=A*L and S=Q*Q'.  Batch the independent eigenvector products into a GEMM
// instead of issuing one matrix-vector product for every active equality.
inline Eigen::MatrixXd zhangR51ConstraintRightBasis(
    const Eigen::MatrixXd& whitened,const Eigen::MatrixXd& eigenvectors,
    const Eigen::VectorXd& eigenvalues,double tolerance)
{
    int rank=0;
    for(int i=0;i<eigenvalues.size();++i)rank+=eigenvalues(i)>tolerance;
    Eigen::MatrixXd active(eigenvectors.rows(),rank);
    Eigen::VectorXd inverseRoots(rank);
    int column=0;
    for(int i=0;i<eigenvalues.size();++i)if(eigenvalues(i)>tolerance) {
        active.col(column)=eigenvectors.col(i);
        inverseRoots(column)=1/std::sqrt(eigenvalues(i));
        ++column;
    }
    Eigen::MatrixXd right=whitened.transpose()*active;
    const int threads=zhangR51ParallelThreads();
    const bool parallel=threads>1 && rank>=64 && whitened.cols()>=512;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
#endif
    for(int c=0;c<rank;++c)right.col(c)*=inverseRoots(c);
    return right;
}
