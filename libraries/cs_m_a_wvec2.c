/********************  cs_m_a_wvec2.c  (in su3.a) ********************
*
*void c_scalar_mult_add_wvec2(wilson_vector *src1, wilson_vector *src2,
	complex s, wilson_vector *dest)
*  Multiply a Wilson vector by a  complex scalar and add to another vector
* dest  <-  src1 + s*src2
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_add_wvec2( wilson_vector *src1,wilson_vector *src2,
	complex s, wilson_vector *dest ){
	wilson_vector src3;
	register int i,j;

	scalar_mult_add_wvec( src1, src2, (s.real), dest );

	for(i=0;i<4;i++) {
	for(j=0;j<3;j++) {
		src3.d[i].c[j].real = -(src2->d[i].c[j].imag);
		src3.d[i].c[j].imag = src2->d[i].c[j].real;
	}
	}

	scalar_mult_add_wvec( dest, &src3, (s.imag), dest);

}
