/******************  cs_m_a_mat.c  (in su3.a) ***************************
*									*
*  c_scalar_mult_add_su3mat( su3_matrix *ma, su3_matrix *m2,		*
*	complex *phase, su3_matrix *m3)					*
*  multiply an su3 matrix by a complex scalar and add it to another	*
*  matrix:   m3 <- m1 + number*m2 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_add_su3mat( su3_matrix *m1, su3_matrix *m2,
	complex *phase, su3_matrix *m3){

#ifndef NATIVEDOUBLE
register int i,j;
complex t;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	t = cmul(&m2->e[i][j],phase);
	m3->e[i][j] = cadd(&m1->e[i][j],&t);
    }

#else
register int i,j;
register double sr,si,br,bi,cr,ci;

    sr = (*phase).real; si = (*phase).imag;

    for(i=0;i<3;i++)for(j=0;j<3;j++){
	br=m2->e[i][j].real; bi=m2->e[i][j].imag;

	cr = sr*br - si*bi;
	ci = sr*bi + si*br;

	m3->e[i][j].real = m1->e[i][j].real + cr;
	m3->e[i][j].imag = m1->e[i][j].imag + ci;
    }
#endif
}
