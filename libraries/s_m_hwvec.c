/********************  s_m_hwvec.c  (in su3.a) ********************
*
*void scalar_mult_hwvec(half_wilson_vector *src, Real s,
	half_wilson_vector *dest)
*  Multiply a half Wilson vector by a scalar
* dest  <-  s*src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_hwvec( half_wilson_vector *src, Real s,
    half_wilson_vector *dest ){

#ifndef NATIVEDOUBLE
  register int i;
  for(i=0;i<2;i++)scalar_mult_su3_vector( &(src->h[i]), s, &(dest->h[i]));

#else /* RS6000 version */

  register double ss;
  ss = s;

  dest->h[0].c[0].real = ss*src->h[0].c[0].real;
  dest->h[0].c[0].imag = ss*src->h[0].c[0].imag;
  dest->h[0].c[1].real = ss*src->h[0].c[1].real;
  dest->h[0].c[1].imag = ss*src->h[0].c[1].imag;
  dest->h[0].c[2].real = ss*src->h[0].c[2].real;
  dest->h[0].c[2].imag = ss*src->h[0].c[2].imag;

  dest->h[1].c[0].real = ss*src->h[1].c[0].real;
  dest->h[1].c[0].imag = ss*src->h[1].c[0].imag;
  dest->h[1].c[1].real = ss*src->h[1].c[1].real;
  dest->h[1].c[1].imag = ss*src->h[1].c[1].imag;
  dest->h[1].c[2].real = ss*src->h[1].c[2].real;
  dest->h[1].c[2].imag = ss*src->h[1].c[2].imag;

#endif
}
