/****************  m_mat_na.c  (in su3.a) *******************************
*									*
* void mult_su3_na( su3_matrix *a,*b,*c )				*
* matrix multiply, second matrix is adjoint 				*
* C  <-  A*B_adjoint							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_na(  su3_matrix *a, su3_matrix *b, su3_matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	x.real=x.imag=0.0;
	for(k=0;k<3;k++){
	    CMUL_J( a->e[i][k] , b->e[j][k] , y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}

/* "Hand coded" routines, clearer coding is up above */
#else

void mult_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
int i,j;
register Real t,ar,ai,br,bi,cr,ci;
    for(i=0;i<3;i++)for(j=0;j<3;j++){

	ar=a->e[i][0].real; ai=a->e[i][0].imag;
	br=b->e[j][0].real; bi=b->e[j][0].imag;
	cr=ar*br; t=ai*bi; cr += t;
	ci=ai*br; t=ar*bi; ci -= t;

	ar=a->e[i][1].real; ai=a->e[i][1].imag;
	br=b->e[j][1].real; bi=b->e[j][1].imag;
	t=ar*br; cr += t; t=ai*bi; cr += t;
	t=ar*bi; ci -= t; t=ai*br; ci += t;

	ar=a->e[i][2].real; ai=a->e[i][2].imag;
	br=b->e[j][2].real; bi=b->e[j][2].imag;
	t=ar*br; cr += t; t=ai*bi; cr += t;
	t=ar*bi; ci -= t; t=ai*br; ci += t;

	c->e[i][j].real=cr;
	c->e[i][j].imag=ci;
    }
}

#endif	/* End of "#ifdef FAST" */
