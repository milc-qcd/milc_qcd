/******************  su3_dot.c  (in su3.a) ******************************
*									*
* complex su3_dot( su3_vector *a, su3_vector *b )			*
* return dot product of two su3_vectors: a^dagger b			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex su3_dot( su3_vector *a, su3_vector *b ){

#ifndef FAST
complex temp1,temp2;
    CMULJ_(a->c[0],b->c[0],temp1) 
    CMULJ_(a->c[1],b->c[1],temp2)
    CSUM(temp1,temp2);
    CMULJ_(a->c[2],b->c[2],temp2)
    CSUM(temp1,temp2);
    return(temp1);

#else /* RS6000 version */

#ifdef NATIVEDOUBLE
  register double ar,ai,br,bi,cr,ci;
#else
  register Real ar,ai,br,bi,cr,ci;
#endif
  register complex cc;

  ar=a->c[0].real;  ai=a->c[0].imag;
  br=b->c[0].real;  bi=b->c[0].imag;
  cr = ar*br + ai*bi;
  ci = ar*bi - ai*br;

  ar=a->c[1].real;  ai=a->c[1].imag;
  br=b->c[1].real;  bi=b->c[1].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;

  ar=a->c[2].real;  ai=a->c[2].imag;
  br=b->c[2].real;  bi=b->c[2].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;

  cc.real = cr;
  cc.imag = ci;
  return(cc);

#endif
}
