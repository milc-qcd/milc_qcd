/*******************  cs_m_vec.c  (in su3.a) ****************************
*									*
*  c_scalar_mult_su3vec():						*
*  multiply an su3 vector by a complex scalar				*
*  dest <- number*src 							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_su3vec( su3_vector *src, complex *phase, su3_vector *dest ){

#ifndef NATIVEDOUBLE
register int i;
    for(i=0;i<3;i++){
	dest->c[i] = cmul(&src->c[i],phase);
    }

#else
register int i;
register double sr,si,br,bi,cr,ci;

    sr = (*phase).real; si = (*phase).imag;

    for(i=0;i<3;i++){
	br=src->c[i].real; bi=src->c[i].imag;

	cr = sr*br - si*bi;
	ci = sr*bi + si*br;

	dest->c[i].real = cr;
	dest->c[i].imag = ci;
    }
#endif
}
