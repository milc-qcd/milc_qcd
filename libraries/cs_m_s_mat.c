/****************  cs_m_s_mat.c  (in su3.a) *****************************
*									*
* void c_scalar_mult_sub_su3mat( su3_matrix *a, su3_matrix *b,		*
*	complex *s, su3_matrix *c)					*
* C <- A - s*B,   A,B and C matrices 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* c <- a - s*b, matrices */
void c_scalar_mult_sub_su3mat( su3_matrix *a, su3_matrix *b, complex *s,
	su3_matrix *c){

#ifndef NATIVEDOUBLE
register int i,j;
complex t;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	t = cmul(&b->e[i][j], s);
	c->e[i][j] = csub(&a->e[i][j], &t);
    }

#else
register int i,j;
register double sr,si,br,bi,cr,ci;

    sr = (*s).real; si = (*s).imag;

    for(i=0;i<3;i++)for(j=0;j<3;j++){
	br=b->e[i][j].real; bi=b->e[i][j].imag;

	cr = sr*br - si*bi;
	ci = sr*bi + si*br;

	c->e[i][j].real = a->e[i][j].real - cr;
	c->e[i][j].imag = a->e[i][j].imag - ci;
    }
#endif
}
