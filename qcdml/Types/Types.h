#ifndef __types_h__
#define __types_h__

#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#include <stdbool.h>
#endif

typedef void*           voidPtr;
typedef char*           CstringRef;
#ifdef __cplusplus
typedef std::complex<double>  doubleComplex;
typedef std::complex<float>   floatComplex;
#else
// C99
typedef double complex  doubleComplex;
typedef float  complex  floatComplex;
#endif

#endif
