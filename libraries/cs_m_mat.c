/****************  cs_m_mat.c  (in su3.a) *******************************
*									*
* void c_scalar_mult_su3mat( su3_matrix *b, complex *s, su3_matrix *c)	*
* C <- s*B,   B and C matrices 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* c <- s*b, matrices */
void c_scalar_mult_su3mat( su3_matrix *b, complex *s, su3_matrix *c ){

#ifndef NATIVEDOUBLE
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	c->e[i][j] = cmul(&b->e[i][j], s);
	/* old: c->e[i][j].real = s.real*b->e[i][j].real-s.imag*b->e[i][j].imag;
	c->e[i][j].imag = s.real*b->e[i][j].imag + s.imag*b->e[i][j].real; */
    }

#else
register int i,j;
register double sr,si,br,bi,cr,ci;

    sr = (*s).real; si = (*s).imag;

    for(i=0;i<3;i++)for(j=0;j<3;j++){
	br=b->e[i][j].real; bi=b->e[i][j].imag;

	cr = sr*br - si*bi;
	ci = sr*bi + si*br;

	c->e[i][j].real = cr;
	c->e[i][j].imag = ci;
    }
#endif
}
