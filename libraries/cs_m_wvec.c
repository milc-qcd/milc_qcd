/********************  	cs_m_wvec.c  (in su3.a) ********************
*
*void c_scalar_mult_wvec(wilson_vector *src, complex *s, wilson_vector *dest)
*  Multiply a Wilson vector by a complex scalar and add to another vector
* dest  <-  s * src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_wvec(wilson_vector *src, complex *phase,
			wilson_vector *dest) {

#ifndef FAST
register int i,j;
complex t;
    for(i=0;i<4;i++){
           for(j=0;j<3;j++){
		CMUL( src->d[i].c[j], *phase, dest->d[i].c[j] );
           }
    }

#else
register int i,j;
#ifdef NATIVEDOUBLE
register double sr,si,br,bi;
#else
register Real sr,si,br,bi;
#endif

    sr = (*phase).real; si = (*phase).imag;

    for(i=0;i<4;i++){
	for(j=0;j<3;j++){
	    br=src->d[i].c[j].real; bi=src->d[i].c[j].imag;

	    dest->d[i].c[j].real = sr*br - si*bi;
	    dest->d[i].c[j].imag = sr*bi + si*br;
	}
    }
#endif
}
