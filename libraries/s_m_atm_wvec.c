/*****************  s_m_atm_wvec.c  (in su3.a) ********************
*
*void scalar_mult_addtm_wvec(wilson_vector *src1, wilson_vector *src2,
	Real s, wilson_vector *dest)
*  Multiply a Wilson vector by a scalar and add to minus one times
*   another vector
* dest  <-  (-1)*src1 + s*src2
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_addtm_wvec(wilson_vector *src1,wilson_vector *src2,
	Real s,wilson_vector *dest){

#ifndef NATIVEDOUBLE
  register int i,j;
  for(i=0;i<4;i++){	/*spins*/
    for(j=0;j<3;j++){  /*colors*/
      dest->d[i].c[j].real = -src1->d[i].c[j].real +
	s*src2->d[i].c[j].real;
      dest->d[i].c[j].imag = -src1->d[i].c[j].imag +
	s*src2->d[i].c[j].imag;
    }
  }

#else /* RS6000 version */

  register double ss;
  ss = s;

  dest->d[0].c[0].real = -src1->d[0].c[0].real + ss*src2->d[0].c[0].real;
  dest->d[0].c[0].imag = -src1->d[0].c[0].imag + ss*src2->d[0].c[0].imag;
  dest->d[0].c[1].real = -src1->d[0].c[1].real + ss*src2->d[0].c[1].real;
  dest->d[0].c[1].imag = -src1->d[0].c[1].imag + ss*src2->d[0].c[1].imag;
  dest->d[0].c[2].real = -src1->d[0].c[2].real + ss*src2->d[0].c[2].real;
  dest->d[0].c[2].imag = -src1->d[0].c[2].imag + ss*src2->d[0].c[2].imag;

  dest->d[1].c[0].real = -src1->d[1].c[0].real + ss*src2->d[1].c[0].real;
  dest->d[1].c[0].imag = -src1->d[1].c[0].imag + ss*src2->d[1].c[0].imag;
  dest->d[1].c[1].real = -src1->d[1].c[1].real + ss*src2->d[1].c[1].real;
  dest->d[1].c[1].imag = -src1->d[1].c[1].imag + ss*src2->d[1].c[1].imag;
  dest->d[1].c[2].real = -src1->d[1].c[2].real + ss*src2->d[1].c[2].real;
  dest->d[1].c[2].imag = -src1->d[1].c[2].imag + ss*src2->d[1].c[2].imag;

  dest->d[2].c[0].real = -src1->d[2].c[0].real + ss*src2->d[2].c[0].real;
  dest->d[2].c[0].imag = -src1->d[2].c[0].imag + ss*src2->d[2].c[0].imag;
  dest->d[2].c[1].real = -src1->d[2].c[1].real + ss*src2->d[2].c[1].real;
  dest->d[2].c[1].imag = -src1->d[2].c[1].imag + ss*src2->d[2].c[1].imag;
  dest->d[2].c[2].real = -src1->d[2].c[2].real + ss*src2->d[2].c[2].real;
  dest->d[2].c[2].imag = -src1->d[2].c[2].imag + ss*src2->d[2].c[2].imag;

  dest->d[3].c[0].real = -src1->d[3].c[0].real + ss*src2->d[3].c[0].real;
  dest->d[3].c[0].imag = -src1->d[3].c[0].imag + ss*src2->d[3].c[0].imag;
  dest->d[3].c[1].real = -src1->d[3].c[1].real + ss*src2->d[3].c[1].real;
  dest->d[3].c[1].imag = -src1->d[3].c[1].imag + ss*src2->d[3].c[1].imag;
  dest->d[3].c[2].real = -src1->d[3].c[2].real + ss*src2->d[3].c[2].real;
  dest->d[3].c[2].imag = -src1->d[3].c[2].imag + ss*src2->d[3].c[2].imag;

#endif
}

