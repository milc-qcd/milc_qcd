/******************  m_mat_an.c  (in su3.a) *****************************
*									*
* void mult_su3_an( su3_matrix *a,*b,*c )				*
* matrix multiply, first matrix is adjoint 				*
* C <-  A_adjoint*B							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	x.real=x.imag=0.0;
	for(k=0;k<3;k++){
	    CMULJ_( a->e[k][i] , b->e[k][j], y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}

/* "Hand coded" routines, clearer coding is up above */
#else

void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
  int j;

#ifdef NATIVEDOUBLE
  register double a0r,a0i,a1r,a1i,a2r,a2i;
  register double b0r,b0i,b1r,b1i,b2r,b2i;
#else
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
#endif

  for(j=0;j<3;j++){
    
    a0r=a->e[0][0].real; a0i=a->e[0][0].imag;
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->e[1][0].real; a1i=a->e[1][0].imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->e[2][0].real; a2i=a->e[2][0].imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;
    
    c->e[0][j].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
    c->e[0][j].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
    a0r=a->e[0][1].real; a0i=a->e[0][1].imag;
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->e[1][1].real; a1i=a->e[1][1].imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->e[2][1].real; a2i=a->e[2][1].imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;
    
    c->e[1][j].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
    c->e[1][j].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
    a0r=a->e[0][2].real; a0i=a->e[0][2].imag;
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->e[1][2].real; a1i=a->e[1][2].imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->e[2][2].real; a2i=a->e[2][2].imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;
    
    c->e[2][j].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
    c->e[2][j].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
    }
}

#endif	/* End of "#ifdef FAST" */
