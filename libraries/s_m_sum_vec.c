/****************  s_m_sum_vec.c  (in su3.a) ****************************
*									*
* void scalar_mult_sum_su3_vector( su3_vector *a, su3_vector *b, Real s )*
* A <- A + s*B,   A and B vectors 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* a <- a + s*b, vectors */
void scalar_mult_sum_su3_vector(su3_vector *a, su3_vector *b, Real s){

#ifndef NATIVEDOUBLE
register int i;
    for(i=0;i<3;i++){
	a->c[i].real += s*b->c[i].real;
	a->c[i].imag += s*b->c[i].imag;
    }

#else /* RS6000 version */

  register double ss;

  ss = s;

  a->c[0].real += ss*b->c[0].real;
  a->c[0].imag += ss*b->c[0].imag;
  a->c[1].real += ss*b->c[1].real;
  a->c[1].imag += ss*b->c[1].imag;
  a->c[2].real += ss*b->c[2].real;
  a->c[2].imag += ss*b->c[2].imag;

#endif
}
