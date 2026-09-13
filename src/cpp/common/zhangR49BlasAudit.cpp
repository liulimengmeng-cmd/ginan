#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <stdexcept>
#include <cctype>
// The configured OpenBLAS backend is LP64 (32-bit Fortran INTEGER).
extern "C" void __real_dgemv_(const char*,const int*,const int*,const double*,const double*,const int*,const double*,const int*,const double*,double*,const int*);
extern "C" void __wrap_dgemv_(const char* t,const int* m,const int* n,const double* a,const double* p,const int* lda,const double* x,const int* ix,const double* b,double* y,const int* iy) {
    const char trans=std::toupper(*t);
    if((trans!='N' && trans!='T' && trans!='C') || *m<0 || *n<0 || *lda<std::max(1,*m) || *ix==0 || *iy==0) {
        std::fprintf(stderr,"R49_INVALID_DGEMV_ARGUMENT TRANS=%c M=%d N=%d LDA=%d INCX=%d INCY=%d interface=FORTRAN_LP64\n",*t,*m,*n,*lda,*ix,*iy);
        void* frames[40];int size=backtrace(frames,40);backtrace_symbols_fd(frames,size,2);
        throw std::runtime_error("R49_INVALID_DGEMV_ARGUMENT");
    }
    __real_dgemv_(t,m,n,a,p,lda,x,ix,b,y,iy);
}
namespace {
struct BlasInjection {
    BlasInjection() {
        if(!std::getenv("ZHANG_R49_BLAS_INJECT_INVALID")) return;
        char t='N';int m=2,n=1,lda=1,inc=1;double a=1,b=0,p[2]={},x[1]={},y[2]={};
        try {__wrap_dgemv_(&t,&m,&n,&a,p,&lda,x,&inc,&b,y,&inc);}
        catch(const std::runtime_error&) {std::fprintf(stderr,"R49_BLAS_INJECTION_CAUGHT\n");std::_Exit(86);}
        std::_Exit(87);
    }
} injection;
}
