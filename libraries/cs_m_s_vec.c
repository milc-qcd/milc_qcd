/*******************  cs_m_s_vec.c  (in su3.a) **************************
*									*
*  c_scalar_mult_sub_su3vec()						*
*  multiply an su3 vector by a complex scalar and subtract it from	*
*  another vector:  v1 <- v1 - number*v2 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_sub_su3vec( su3_vector *v1, complex *phase, su3_vector *v2 ){

#ifndef NATIVEDOUBLE
register int i;
complex t;
    for(i=0;i<3;i++){
	t = cmul(&v2->c[i],phase);
	v1->c[i] = csub(&v1->c[i],&t);
    }
#else
register int i;
register double sr,si,br,bi,cr,ci;

    sr = (*phase).real; si = (*phase).imag;

    for(i=0;i<3;i++){
	br=v2->c[i].real; bi=v2->c[i].imag;

	cr = sr*br - si*bi;
	ci = sr*bi + si*br;

	v1->c[i].real -= cr;
	v1->c[i].imag -= ci;
    }
#endif
}
