/******************  s_m_vec.c  (in su3.a) ******************************
*									*
* void scalar_mult_su3_vector( su3_vector *a, Real s, su3_vector *c)	*
* C <- s*A,  A and C vectors 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* c <- s*a, vectors */
void scalar_mult_su3_vector( su3_vector *a, Real s, su3_vector *c){

#ifndef NATIVEDOUBLE
register int i;
    for(i=0;i<3;i++){
	c->c[i].real = s*a->c[i].real;
	c->c[i].imag = s*a->c[i].imag;
    }

#else /* RS6000 version */

  register double ss;

  ss = s;

  c->c[0].real = ss*a->c[0].real;
  c->c[0].imag = ss*a->c[0].imag;
  c->c[1].real = ss*a->c[1].real;
  c->c[1].imag = ss*a->c[1].imag;
  c->c[2].real = ss*a->c[2].real;
  c->c[2].imag = ss*a->c[2].imag;

#endif
}
