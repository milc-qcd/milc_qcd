/*****************  su3_s_m_a_hwvec_proj.c  (in su3.a) *******************
*									 *
* void scalar_mult_add hwvec_proj( su3_matrix *a, half_wilson_vector *b, * 
*    half_wilson_vector *c, Real *s, su3_matrix *d )	                 *
*  D  <- A + B[0] * C[0]_adj * s[0] + B[1] *C[1]_adj * s[1]              *
*/

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_add_hwvec_proj( su3_matrix * const a, 
				 half_wilson_vector * const b, 
				 half_wilson_vector * const c, 
				 Real * const s, su3_matrix *d )
{
    int i,j;
    Real tmp0,tmp1;

#ifdef FAST
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	tmp1 = b->h[0].c[i].real * c->h[0].c[j].real;
	tmp0 = b->h[0].c[i].imag * c->h[0].c[j].imag;
	d->e[i][j].real = a->e[i][j].real + (tmp0 + tmp1)*s[0];

	tmp1 = b->h[0].c[i].real * c->h[0].c[j].imag;
	tmp0 = b->h[0].c[i].imag * c->h[0].c[j].real;
	d->e[i][j].imag += (tmp0 - tmp1)*s[0];

	tmp1 = b->h[1].c[i].real * c->h[1].c[j].real;
	tmp0 = b->h[1].c[i].imag * c->h[1].c[j].imag;
	d->e[i][j].real += (tmp0 + tmp1)*s[1];

	tmp1 = b->h[1].c[i].real * c->h[1].c[j].imag;
	tmp0 = b->h[1].c[i].imag * c->h[1].c[j].real;
	d->e[i][j].imag += (tmp0 - tmp1)*s[1];
    }
#else
    su3_matrix tmat;

    su3_projector(&(b->h[0]), &(c->h[0]), &tmat);
    scalar_mult_add_su3_matrix(d, &tmat,  s[0], d );
    su3_projector(&(b->h[1]), &(c->h[1]), &tmat);
    scalar_mult_add_su3_matrix(d, &tmat,  s[1], d );
#endif
}
