/***************** su3_proj_w.c  (in su3.a) ****************************/
/* MIMD version 7 */
/*									*
* void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c )
* C  <- sum over spins of outer product of A.d[i] and B.d[i]		*
*  C_ij = sum( A_i * B_adjoint_j )					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c ){
register int i,j,k;
register complex cc;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	c->e[i][j] = cmplx(0.0,0.0);
	for(k=0;k<4;k++){
	    CMUL_J( a->d[k].c[i], b->d[k].c[j], cc ); CSUM( c->e[i][j], cc );
	}
    }
}

#else
#ifdef NATIVEDOUBLE   /* RS6000 version */

void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c ){
  register int i,j;
  register double ar,ai,br,bi,cr,ci;

    for(i=0;i<3;i++)for(j=0;j<3;j++){
	ar=a->d[0].c[i].real;  ai=a->d[0].c[i].imag;
	br=b->d[0].c[j].real;  bi=b->d[0].c[j].imag;
	cr = ar*br + ai*bi;
	ci = ai*br - ar*bi;

	ar=a->d[1].c[i].real;  ai=a->d[1].c[i].imag;
	br=b->d[1].c[j].real;  bi=b->d[1].c[j].imag;
	cr += ar*br + ai*bi;
	ci += ai*br - ar*bi;

	ar=a->d[2].c[i].real;  ai=a->d[2].c[i].imag;
	br=b->d[2].c[j].real;  bi=b->d[2].c[j].imag;
	cr += ar*br + ai*bi;
	ci += ai*br - ar*bi;

	ar=a->d[3].c[i].real;  ai=a->d[3].c[i].imag;
	br=b->d[3].c[j].real;  bi=b->d[3].c[j].imag;
	cr += ar*br + ai*bi;
	ci += ai*br - ar*bi;

	c->e[i][j].real = cr;
	c->e[i][j].imag = ci;
    }
}
#else
void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c ){
register int i,j,k;
register Real tmp_r,tmp_i,tmp2;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	tmp_r = tmp_i = 0.0;
	for(k=0;k<4;k++){
	    tmp2 = a->d[k].c[i].real * b->d[k].c[j].real; tmp_r = tmp_r + tmp2;
	    tmp2 = a->d[k].c[i].imag * b->d[k].c[j].imag; tmp_r = tmp_r + tmp2;
	    tmp2 = a->d[k].c[i].imag * b->d[k].c[j].real; tmp_i = tmp_i + tmp2;
	    tmp2 = a->d[k].c[i].real * b->d[k].c[j].imag; tmp_i = tmp_i - tmp2;
	}

	c->e[i][j].real = tmp_r;
	c->e[i][j].imag = tmp_i;
    }
}
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif /* end ifdef FAST */
