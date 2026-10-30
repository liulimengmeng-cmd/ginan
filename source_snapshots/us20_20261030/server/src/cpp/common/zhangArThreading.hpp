#pragma once
#include "common/zhangFloatThreading.hpp"

namespace ZhangArThreading
{
inline thread_local bool active=false;

// Independent from FLOAT: keep the old inherited settings outside R51 AR.
class ScopedThreads
{
    const bool previousActive=active;
    const int previousEigen=Eigen::nbThreads();
#ifdef _OPENMP
    const int previousOmp=omp_get_max_threads();
    const int previousDynamic=omp_get_dynamic();
    const int previousLevels=omp_get_max_active_levels();
#endif
    static int parse(bool enabled)
    {
        if(!enabled)return 0;
        try {return ZhangFloatThreading::requestedThreads(
            std::getenv("ZHANG_R51_AR_THREADS"),ZhangFloatThreading::availableThreads());}
        catch(const std::invalid_argument&) {
            throw std::invalid_argument("ZHANG_R51_AR_THREADS must be all, legacy, or 1..available CPUs");
        }
    }
    const int requested;
    BlasThreading::ScopedOpenBlasThreadLimit blas{requested};
public:
    explicit ScopedThreads(bool enabled) : requested(parse(enabled)),blas(requested)
    {
        active=enabled && requested>0;
        if(!requested)return;
#ifdef _OPENMP
        omp_set_dynamic(0);omp_set_max_active_levels(1);omp_set_num_threads(requested);
#endif
        Eigen::setNbThreads(1);
    }
    ~ScopedThreads()
    {
        active=previousActive;
        if(!requested)return;
#ifdef _OPENMP
        omp_set_num_threads(previousOmp);omp_set_dynamic(previousDynamic);omp_set_max_active_levels(previousLevels);
#endif
        Eigen::setNbThreads(previousEigen);
    }
    int threads() const{return requested;}
    ScopedThreads(const ScopedThreads&)=delete;
    ScopedThreads& operator=(const ScopedThreads&)=delete;
};

inline int eigensolverThreads()
{
#ifdef GINAN_USE_OPENBLAS
    return std::min(BlasThreading::openblas_get_num_threads(),ZhangFloatThreading::physicalCores());
#else
    return 0;
#endif
}

// Eigen's tridiagonalisation issues many small BLAS operations. A 1175-state
// production-conditioner profile was ~18x slower at BLAS=12 than BLAS=1.
// Bound only this backend; surrounding GEMMs and LAPACK DSYEVD keep their
// own parallel budgets. Non-R51 callers retain their inherited policy.
inline Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSelfAdjoint(const MatrixXd& covariance)
{
    BlasThreading::ScopedOpenBlasThreadLimit serialEigen(active?1:0);
    return Eigen::SelfAdjointEigenSolver<MatrixXd>(covariance);
}
}
