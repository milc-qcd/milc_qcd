/******************  wvec_dot.c  (in su3.a) ****************************/
/* MIMD version 7 */
/*									*
* complex wvec_dot(a,b) wilson_vector *a,*b;				*
* return dot product of two wilson_vectors					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex wvec_dot( wilson_vector *a, wilson_vector *b ){

#ifndef FAST
  complex temp1,temp2;
  register int i;
  temp1.real = temp1.imag = 0.0;
  for(i=0;i<4;i++){
    CMULJ_(a->d[i].c[0],b->d[i].c[0],temp2); CSUM(temp1,temp2);
    CMULJ_(a->d[i].c[1],b->d[i].c[1],temp2); CSUM(temp1,temp2);
    CMULJ_(a->d[i].c[2],b->d[i].c[2],temp2); CSUM(temp1,temp2);
  }
  return(temp1);

#else

#ifdef NATIVEDOUBLE
  register double ar,ai,br,bi,cr,ci;
#else
  register Real ar,ai,br,bi,cr,ci;
#endif
  register complex cc;


  ar=a->d[0].c[0].real;  ai=a->d[0].c[0].imag;
  br=b->d[0].c[0].real;  bi=b->d[0].c[0].imag;
  cr = ar*br + ai*bi;
  ci = ar*bi - ai*br;
  ar=a->d[0].c[1].real;  ai=a->d[0].c[1].imag;
  br=b->d[0].c[1].real;  bi=b->d[0].c[1].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[0].c[2].real;  ai=a->d[0].c[2].imag;
  br=b->d[0].c[2].real;  bi=b->d[0].c[2].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;

  ar=a->d[1].c[0].real;  ai=a->d[1].c[0].imag;
  br=b->d[1].c[0].real;  bi=b->d[1].c[0].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[1].c[1].real;  ai=a->d[1].c[1].imag;
  br=b->d[1].c[1].real;  bi=b->d[1].c[1].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[1].c[2].real;  ai=a->d[1].c[2].imag;
  br=b->d[1].c[2].real;  bi=b->d[1].c[2].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;

  ar=a->d[2].c[0].real;  ai=a->d[2].c[0].imag;
  br=b->d[2].c[0].real;  bi=b->d[2].c[0].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[2].c[1].real;  ai=a->d[2].c[1].imag;
  br=b->d[2].c[1].real;  bi=b->d[2].c[1].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[2].c[2].real;  ai=a->d[2].c[2].imag;
  br=b->d[2].c[2].real;  bi=b->d[2].c[2].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;

  ar=a->d[3].c[0].real;  ai=a->d[3].c[0].imag;
  br=b->d[3].c[0].real;  bi=b->d[3].c[0].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[3].c[1].real;  ai=a->d[3].c[1].imag;
  br=b->d[3].c[1].real;  bi=b->d[3].c[1].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;
  ar=a->d[3].c[2].real;  ai=a->d[3].c[2].imag;
  br=b->d[3].c[2].real;  bi=b->d[3].c[2].imag;
  cr += ar*br + ai*bi;
  ci += ar*bi - ai*br;

  cc.real = cr;
  cc.imag = ci;
  return(cc);

#endif
}
