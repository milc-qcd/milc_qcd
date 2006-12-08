#ifndef _MILC_COMPLEX_H
#define _MILC_COMPLEX_H

/*============================================================================*/
/*									      */
/* Complex Numbers							      */
/*									      */
/* Typedefs are included for type complex (single-precision) and type         */
/* double_complex (double_precision) complex numbers.  At this time, the      */
/* functions cannot be overloaded, so there are separate routines for the     */
/* single and double precision types.  All of the macros, however, will work  */
/* with both types and mix types freely.				      */
/*									      */
/* The following functions are provided in single and double precision:       */
/*									      */
/*   complex cmplx(Real r, Real i);           (r,i)   		      */
/*   complex cadd(complex *a, complex *b);      *a + *b   		      */
/*   complex cmul(complex *a, complex *b);      *a * *b   		      */
/*   complex csub(complex *a, complex *b);      *a - *b   		      */
/*   complex cdiv(complex *a, complex *b);      *a / *b   		      */
/*   complex conjg(complex *a);	       conjugate of *a   		      */
/*   complex cexp(complex *a);                  exp(*a)   		      */
/*   complex clog(complex *a);                  ln(*a)    		      */
/*   complex csqrt(complex *a);      sqrt(a)		  		      */
/*   complex ce_itheta(Real theta);      exp(i*theta)   		      */
/*									      */
/* The following macros are provided, which work for BOTH single and double   */
/* precision and for mixtures:						      */
/*									      */
/* 1) Macros which appear to return values (Real or double, as appropriate): */
/*    cabs(*a)        magnitude of the complex number *a		      */
/*    cabs_sq(*a)     square of the magnitude (faster than cabs)              */
/*    carg(*a)        phase of the complex number *a			      */
/*									      */
/* 2) Macro to convert from single to double or double to single:             */
/*    set_complex_equal(*a,*b)    do *b=*a by components to convert           */
/*									      */
/* 3) Macros for fast in-line operations:				      */
/*    CONJG(a,b)        b = conjg(a)					      */
/*    CADD(a,b,c)       c = a + b					      */
/*    CSUM(a,b)         a += b						      */
/*    CSUB(a,b,c)       c = a - b					      */
/*    CMUL(a,b,c)       c = a * b					      */
/*    CDIV(a,b,c)       c = a / b					      */
/*    CMUL_J(a,b,c)     c = a * conjg(b)				      */
/*    CMULJ_(a,b,c)     c = conjg(a) * b				      */
/*    CMULJJ(a,b,c)     c = conjg(a*b)					      */
/*    CNEGATE(a,b)      b = -a						      */
/*    CMUL_I(a,b)       b = ia						      */
/*    CMUL_MINUS_I(a,b) b = -ia						      */
/*    CMULREAL(a,b,c)   c = ba with b real and a complex                      */
/*    CDIVREAL(a,b,c)   c = a/b with a complex and b real		      */
/*    CSUM_TPI(a,b)     a += i*b with a and b complex    		      */
/*    CSUM_TMI(a,b)     a += -i*b with a and b complex    		      */
/*    									      */
/*============================================================================*/

/* On the paragon, under OSF, complex is defined in math.h, but not 
  quite the way we did it, so redefine it:
*/
#include "../include/precision.h"

#ifdef HPUX
#define complex complexx
#endif

/* Rename to prevent conflict with gcc 3.xx standard functions */
#define cexp  cexp_single
#define clog  clog_single
#define csqrt csqrt_single

/* generic precision complex number definition */
/* specific for float complex */
typedef struct {   
  float real;	   
  float imag; 
} fcomplex;  

/* specific for double complex */
typedef struct {
   double real;
   double imag;
} double_complex;

#if (PRECISION==1)
#define complex fcomplex
#else
#define complex dcomplex
#endif

/* Alternative name for double complex */
typedef double_complex dcomplex;

/* define complex as a union to ensure alignment to doubleword boundary */
/*typedef union {    ** standard complex number declaration for single- **
  Real f[2];             ** precision complex numbers **
   double dummy;
} complex;
typedef struct {           ** standard complex number declaration for double- **
   double f[2];		   ** precision complex numbers			      **
} double_complex; */
/*#define real f[0] */
/*#define imag f[1] */


/* Generic Precision Function Prototypes for Complex Numbers */
complex cmplx(  Real x, Real y );
complex cadd( complex *a, complex *b );
complex cmul( complex *a, complex *b );
complex csub( complex *a, complex *b );
complex cdiv( complex *a, complex *b );
complex conjg( complex *a );
complex cexp( complex *a );       
complex clog( complex *a );        
complex csqrt( complex *z );        
complex ce_itheta( Real theta );    

double_complex dcmplx( double x, double y );
double_complex dcadd( double_complex *a, double_complex *b );
double_complex dcmul( double_complex *a, double_complex *b );
double_complex dcsub( double_complex *a, double_complex *b );
double_complex dcdiv( double_complex *a, double_complex *b );
double_complex dconjg(  double_complex *a );
double_complex dcexp(  double_complex *a ); 
double_complex dclog(  double_complex *a );  
double_complex dcsqrt( double_complex *z );  
double_complex dce_itheta( double theta );   

/* Macros for Complex Numbers */
								/* *b = *a    */
#define set_complex_equal(a,b) { (*b).real=(*a).real; (*b).imag=(*a).imag; }
								/*    |*a|    */
#define cabs(a) (sqrt( (*a).real*(*a).real + (*a).imag*(*a).imag ) )
								/*  *a * *a*  */
#define dcabs cabs
#define cabs_sq(a) ( (*a).real*(*a).real + (*a).imag*(*a).imag )
								/* phase(*a)  */
#define carg(a) (atan2((double)(*a).imag, (double)(*a).real ) )
								/*   b = a*   */
#define dcarg carg
#define CONJG(a,b) { (b).real = (a).real; (b).imag = -(a).imag; }
								/*  c = a + b */
#define CADD(a,b,c) { (c).real = (a).real + (b).real;  \
		      (c).imag = (a).imag + (b).imag; }
								/*  a += b    */
#define CSUM(a,b) { (a).real += (b).real; (a).imag += (b).imag; }
								/*  c = a - b */
#define CSUB(a,b,c) { (c).real = (a).real - (b).real;  \
		      (c).imag = (a).imag - (b).imag; }
								/*  c = a * b */
#define CMUL(a,b,c) { (c).real = (a).real*(b).real - (a).imag*(b).imag; \
		      (c).imag = (a).real*(b).imag + (a).imag*(b).real; }
								/* c = a / b  */
#define CDIV(a,b,c) { double t = (b).real*(b).real + (b).imag*(b).imag; \
		      (c).real = ((a).real*(b).real + (a).imag*(b).imag)/t; \
		      (c).imag = ((a).imag*(b).real - (a).real*(b).imag)/t; }
								/* c = a * b* */
#define CMUL_J(a,b,c) { (c).real = (a).real*(b).real + (a).imag*(b).imag; \
	  	        (c).imag = (a).imag*(b).real - (a).real*(b).imag; }
								/* c = a* * b */
#define CMULJ_(a,b,c) { (c).real = (a).real*(b).real + (a).imag*(b).imag; \
		        (c).imag = (a).real*(b).imag - (a).imag*(b).real; }
								/* c = (a*b)* */
#define CMULJJ(a,b,c) { (c).real =  (a).real*(b).real - (a).imag*(b).imag; \
		        (c).imag = -(a).real*(b).imag - (a).imag*(b).real; }
								/* b = - a    */
#define CNEGATE(a,b) { (b).real = -(a).real; (b).imag = -(a).imag; }
								/* b =  ia    */
#define CMUL_I(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
								/* b = -ia    */
#define CMUL_MINUS_I(a,b) { (b).real = (a).imag; (b).imag = -(a).real; }
								/* c = ba     */
#define CMULREAL(a,b,c) { (c).real = (b) * (a).real; (c).imag = (b)*(a).imag; }
								/* c = a/b    */
#define CDIVREAL(a,b,c) { (c).real = (a).real/(b); (c).imag = (a).imag/(b); }

                                                               /* a += i*b */
#define CSUM_TPI(a,b) { (a).real -= (b).imag; (a).imag +=  (b).real; } 

                                                                /* a += -i*b */
#define CSUM_TMI(a,b) { (a).real += (b).imag; (a).imag -=  (b).real; }

#endif	/* _MILC_COMPLEX_H */
