#ifndef _MILC_DATATYPES_H
#define _MILC_DATATYPES_H

/************************************************************/
/* From old complex.h */

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

/* Alternative name for double complex */
typedef double_complex dcomplex;

/* specific for long double complex */
typedef struct {
   long double real;
   long double imag;
} long_double_complex;


#if (MILC_PRECISION==1)
#define complex fcomplex
#else
#define complex dcomplex
#endif

#ifdef HAVE_QUDA
// When using QUDA, we need to set the site struct member arrays alignment to a multiple of 16 bytes
#define ALIGNAS(n) __attribute__((aligned(n)))
#define ALIGNMENT ALIGNAS(16)
#else
#define ALIGNMENT
#endif


/************************************************************/
/* From old su3.h */

typedef struct { fcomplex e[3][3]; } fsu3_matrix;
typedef struct { fcomplex c[3]; } fsu3_vector;
typedef struct { 
  fcomplex m01,m02,m12; 
  float m00im,m11im,m22im; 
  float space; } fanti_hermitmat;

typedef struct { dcomplex e[3][3]; } dsu3_matrix;
typedef struct { dcomplex c[3]; } dsu3_vector;
typedef struct { 
  dcomplex m01,m02,m12; 
  double m00im,m11im,m22im; 
  double space; } danti_hermitmat;

#if (MILC_PRECISION==1)

#define su3_matrix      fsu3_matrix
#define su3_vector      fsu3_vector
#define anti_hermitmat  fanti_hermitmat

#else

#define su3_matrix      dsu3_matrix
#define su3_vector      dsu3_vector
#define anti_hermitmat  danti_hermitmat

#endif

/* For KS spectroscopy */
typedef struct {
  su3_vector ** v;  /* For the nc vector fields */
  int nc;           /* number of vectors */
} ks_prop_field;

/* Used in HISQ codes */
/* Rank 4 tensor for storing derivatives */
typedef struct { fcomplex t4[3][3][3][3]; } fsu3_tensor4;
typedef struct { dcomplex t4[3][3][3][3]; } dsu3_tensor4;

#if (MILC_PRECISION==1)
#define su3_tensor4 fsu3_tensor4
#else
#define su3_tensor4 dsu3_tensor4
#endif

/* SU(2) */
typedef struct { complex e[2][2]; } su2_matrix;

/* Wilson vectors */

/* e.g. 					     */
/* wilson_propagator prop;                           */
/* prop.c[ci].d[si].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          su3_vector            */
/* ----------->                wilson_vector         */
/* ----->                      spin_wilson_vector    */
/* e.g. 					     */
/* wilson_matrix matr;                               */
/* matr.d[si].c[ci].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          su3_vector            */
/* ----------->                wilson_vector         */
/* ----->                      color_wilson_vector   */

/* Object with two Dirac and two color indices. A given element
   of a "wilson_propagator" is accessed by
   object.c[color1].d[spin1].d[spin2].c[color2].real , etc.
   As alway, "d" denotes a Dirac index and "c" a color index.
   "1" refers to the source, "2" to the sink.
*/

typedef struct { fsu3_vector d[4]; } fwilson_vector;
typedef struct { fsu3_vector h[2]; } fhalf_wilson_vector;
typedef struct { fwilson_vector c[3]; } fcolor_wilson_vector;
typedef struct { fwilson_vector d[4]; } fspin_wilson_vector;
typedef struct { fcolor_wilson_vector d[4]; } fwilson_matrix;
typedef struct { fspin_wilson_vector c[3]; } fwilson_propagator;

typedef struct { dsu3_vector d[4]; } dwilson_vector;
typedef struct { dsu3_vector h[2]; } dhalf_wilson_vector;
typedef struct { dwilson_vector c[3]; } dcolor_wilson_vector;
typedef struct { dwilson_vector d[4]; } dspin_wilson_vector;
typedef struct { dcolor_wilson_vector d[4]; } dwilson_matrix;
typedef struct { dspin_wilson_vector c[3]; } dwilson_propagator;

#if (MILC_PRECISION==1)

#define wilson_vector       fwilson_vector
#define half_wilson_vector  fhalf_wilson_vector
#define color_wilson_vector fcolor_wilson_vector
#define spin_wilson_vector  fspin_wilson_vector
#define wilson_matrix       fwilson_matrix
#define wilson_propagator   fwilson_propagator

#else

#define wilson_vector       dwilson_vector
#define half_wilson_vector  dhalf_wilson_vector
#define color_wilson_vector dcolor_wilson_vector
#define spin_wilson_vector  dspin_wilson_vector
#define wilson_matrix       dwilson_matrix
#define wilson_propagator   dwilson_propagator

#endif

/* For Wilson spectroscopy */
typedef struct {
  spin_wilson_vector ** swv;
  int nc;
} wilson_prop_field;

#endif /* _MILC_DATATYPES_H */
