/********************  s_m_a_wvec.c  (in su3.a) ********************
*
*void scalar_mult_add_wvec(wilson_vector *src1, wilson_vector *src2,
	Real s, wilson_vector *dest)
*  Multiply a Wilson vector by a scalar and add to another vector
* dest  <-  src1 + s*src2
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_add_wvec(wilson_vector *src1,wilson_vector *src2,Real
	s, wilson_vector *dest) {

#ifndef NATIVEDOUBLE
  register int i;
  for(i=0;i<4;i++)scalar_mult_add_su3_vector( &(src1->d[i]), &(src2->d[i]),
					     s, &(dest->d[i]));
  
#else /* RS6000 version */
  
  register double ss;
  ss = s;
  
  dest->d[0].c[0].real = src1->d[0].c[0].real + ss*src2->d[0].c[0].real;
  dest->d[0].c[0].imag = src1->d[0].c[0].imag + ss*src2->d[0].c[0].imag;
  dest->d[0].c[1].real = src1->d[0].c[1].real + ss*src2->d[0].c[1].real;
  dest->d[0].c[1].imag = src1->d[0].c[1].imag + ss*src2->d[0].c[1].imag;
  dest->d[0].c[2].real = src1->d[0].c[2].real + ss*src2->d[0].c[2].real;
  dest->d[0].c[2].imag = src1->d[0].c[2].imag + ss*src2->d[0].c[2].imag;
  
  dest->d[1].c[0].real = src1->d[1].c[0].real + ss*src2->d[1].c[0].real;
  dest->d[1].c[0].imag = src1->d[1].c[0].imag + ss*src2->d[1].c[0].imag;
  dest->d[1].c[1].real = src1->d[1].c[1].real + ss*src2->d[1].c[1].real;
  dest->d[1].c[1].imag = src1->d[1].c[1].imag + ss*src2->d[1].c[1].imag;
  dest->d[1].c[2].real = src1->d[1].c[2].real + ss*src2->d[1].c[2].real;
  dest->d[1].c[2].imag = src1->d[1].c[2].imag + ss*src2->d[1].c[2].imag;
  
  dest->d[2].c[0].real = src1->d[2].c[0].real + ss*src2->d[2].c[0].real;
  dest->d[2].c[0].imag = src1->d[2].c[0].imag + ss*src2->d[2].c[0].imag;
  dest->d[2].c[1].real = src1->d[2].c[1].real + ss*src2->d[2].c[1].real;
  dest->d[2].c[1].imag = src1->d[2].c[1].imag + ss*src2->d[2].c[1].imag;
  dest->d[2].c[2].real = src1->d[2].c[2].real + ss*src2->d[2].c[2].real;
  dest->d[2].c[2].imag = src1->d[2].c[2].imag + ss*src2->d[2].c[2].imag;
  
  dest->d[3].c[0].real = src1->d[3].c[0].real + ss*src2->d[3].c[0].real;
  dest->d[3].c[0].imag = src1->d[3].c[0].imag + ss*src2->d[3].c[0].imag;
  dest->d[3].c[1].real = src1->d[3].c[1].real + ss*src2->d[3].c[1].real;
  dest->d[3].c[1].imag = src1->d[3].c[1].imag + ss*src2->d[3].c[1].imag;
  dest->d[3].c[2].real = src1->d[3].c[2].real + ss*src2->d[3].c[2].real;
  dest->d[3].c[2].imag = src1->d[3].c[2].imag + ss*src2->d[3].c[2].imag;

#endif
}
