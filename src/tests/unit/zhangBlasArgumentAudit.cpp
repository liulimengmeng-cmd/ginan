#include <algorithm>
#include <sstream>
#include <stdexcept>

extern "C" void __real_dgemv_(const char*, const int*, const int*, const double*,
    const double*, const int*, const double*, const int*, const double*, double*, const int*);

extern "C" void __wrap_dgemv_(const char* trans, const int* m, const int* n,
    const double* alpha, const double* a, const int* lda, const double* x,
    const int* incx, const double* beta, double* y, const int* incy)
{
    if (*m < 0 || *n < 0 || *lda < std::max(1, *m) || *incx == 0 || *incy == 0)
    {
        std::ostringstream message;
        message << "INVALID_DGEMV_ARGUMENT trans=" << *trans << " m=" << *m
                << " n=" << *n << " lda=" << *lda << " incx=" << *incx
                << " incy=" << *incy;
        throw std::runtime_error(message.str());
    }
    __real_dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
