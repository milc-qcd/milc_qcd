/********************  s_m_wvec.c  (in su3.a) ********************
*
*void scalar_mult_wvec(wilson_vector *src, Real s, wilson_vector *dest)
*  Multiply a Wilson vector by a scalar
* dest  <-  s*src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_wvec( wilson_vector *src, Real s, wilson_vector *dest){

#ifndef NATIVEDOUBLE
register int i;
    for(i=0;i<4;i++)scalar_mult_su3_vector( &(src->d[i]), s, &(dest->d[i]));

#else /* RS6000 version */

  register double ss;
  ss = s;

  dest->d[0].c[0].real = ss*src->d[0].c[0].real;
  dest->d[0].c[0].imag = ss*src->d[0].c[0].imag;
  dest->d[0].c[1].real = ss*src->d[0].c[1].real;
  dest->d[0].c[1].imag = ss*src->d[0].c[1].imag;
  dest->d[0].c[2].real = ss*src->d[0].c[2].real;
  dest->d[0].c[2].imag = ss*src->d[0].c[2].imag;

  dest->d[1].c[0].real = ss*src->d[1].c[0].real;
  dest->d[1].c[0].imag = ss*src->d[1].c[0].imag;
  dest->d[1].c[1].real = ss*src->d[1].c[1].real;
  dest->d[1].c[1].imag = ss*src->d[1].c[1].imag;
  dest->d[1].c[2].real = ss*src->d[1].c[2].real;
  dest->d[1].c[2].imag = ss*src->d[1].c[2].imag;

  dest->d[2].c[0].real = ss*src->d[2].c[0].real;
  dest->d[2].c[0].imag = ss*src->d[2].c[0].imag;
  dest->d[2].c[1].real = ss*src->d[2].c[1].real;
  dest->d[2].c[1].imag = ss*src->d[2].c[1].imag;
  dest->d[2].c[2].real = ss*src->d[2].c[2].real;
  dest->d[2].c[2].imag = ss*src->d[2].c[2].imag;

  dest->d[3].c[0].real = ss*src->d[3].c[0].real;
  dest->d[3].c[0].imag = ss*src->d[3].c[0].imag;
  dest->d[3].c[1].real = ss*src->d[3].c[1].real;
  dest->d[3].c[1].imag = ss*src->d[3].c[1].imag;
  dest->d[3].c[2].real = ss*src->d[3].c[2].real;
  dest->d[3].c[2].imag = ss*src->d[3].c[2].imag;

#endif
}
