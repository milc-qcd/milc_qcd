/****************  m_matvec.c  (in su3.a) *******************************
*									*
* void mult_su3_mat_vec( su3_matrix *a, su3_vector *b,*c )		*
* matrix times vector multiply, no adjoints 				*
*  C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c  ){
register int i,j;
register complex x,y;
    for(i=0;i<3;i++){
	x.real=x.imag=0.0;
	for(j=0;j<3;j++){
	    CMUL( a->e[i][j] , b->c[j] , y )
	    CSUM( x , y );
	}
	c->c[i] = x;
    }
}
#else
#ifdef NATIVEDOUBLE   /* RS6000 version */
void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c ){

  register double a0r,a0i,a1r,a1i,a2r,a2i;
  register double b0r,b0i,b1r,b1i,b2r,b2i;
  
  a0r=a->e[0][0].real; a0i=a->e[0][0].imag;
  b0r=b->c[0].real;    b0i=b->c[0].imag;
  a1r=a->e[0][1].real; a1i=a->e[0][1].imag;
  b1r=b->c[1].real;    b1i=b->c[1].imag;
  a2r=a->e[0][2].real; a2i=a->e[0][2].imag;
  b2r=b->c[2].real;    b2i=b->c[2].imag;

  c->c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  c->c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=a->e[1][0].real; a0i=a->e[1][0].imag;
  b0r=b->c[0].real;    b0i=b->c[0].imag;
  a1r=a->e[1][1].real; a1i=a->e[1][1].imag;
  b1r=b->c[1].real;    b1i=b->c[1].imag;
  a2r=a->e[1][2].real; a2i=a->e[1][2].imag;
  b2r=b->c[2].real;    b2i=b->c[2].imag;

  c->c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  c->c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=a->e[2][0].real; a0i=a->e[2][0].imag;
  b0r=b->c[0].real;    b0i=b->c[0].imag;
  a1r=a->e[2][1].real; a1i=a->e[2][1].imag;
  b1r=b->c[1].real;    b1i=b->c[1].imag;
  a2r=a->e[2][2].real; a2i=a->e[2][2].imag;
  b2r=b->c[2].real;    b2i=b->c[2].imag;

  c->c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  c->c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

}

#else
void mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c ){
int i;
register Real t,ar,ai,br,bi,cr,ci;
    for(i=0;i<3;i++){

	ar=a->e[i][0].real; ai=a->e[i][0].imag;
	br=b->c[0].real; bi=b->c[0].imag;
	cr=ar*br; t=ai*bi; cr -= t;
	ci=ar*bi; t=ai*br; ci += t;

	ar=a->e[i][1].real; ai=a->e[i][1].imag;
	br=b->c[1].real; bi=b->c[1].imag;
	t=ar*br; cr += t; t=ai*bi; cr -= t;
	t=ar*bi; ci += t; t=ai*br; ci += t;

	ar=a->e[i][2].real; ai=a->e[i][2].imag;
	br=b->c[2].real; bi=b->c[2].imag;
	t=ar*br; cr += t; t=ai*bi; cr -= t;
	t=ar*bi; ci += t; t=ai*br; ci += t;

	c->c[i].real=cr;
	c->c[i].imag=ci;
    }
}
#endif	/* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#infdef FAST" */
