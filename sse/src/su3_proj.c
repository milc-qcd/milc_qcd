/*****************  su3_proj.c  (in su3.a) ******************************
*									*
* void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c )	*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	CMUL_J( a->c[i], b->c[j], c->e[i][j] );
    }
}

#else
#ifdef NATIVEDOUBLE   /* RS6000 version */

void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c ){

  register int i,j;
  register double ar,ai,br,bi;

    for(i=0;i<3;i++){
	ar=a->c[i].real;  ai=a->c[i].imag;
	for(j=0;j<3;j++){
	    br=b->c[j].real;  bi=b->c[j].imag;
	    c->e[i][j].real = ar*br + ai*bi;
	    c->e[i][j].imag = ai*br - ar*bi;
	}
    }
}
#else

void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c ){
register int i,j;
register Real tmp,tmp2;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	tmp2 = a->c[i].real * b->c[j].real;
	tmp = a->c[i].imag * b->c[j].imag;
	c->e[i][j].real = tmp + tmp2;
	tmp2 = a->c[i].real * b->c[j].imag;
	tmp = a->c[i].imag * b->c[j].real;
	c->e[i][j].imag = tmp - tmp2;
    }
}
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif /* end ifdef FAST */
