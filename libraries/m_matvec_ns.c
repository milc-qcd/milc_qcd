/*****************  m_matvec_ns.c  (in su3.a) ***************************
*									*
* void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b,*c )		*
* su3_matrix times su3_vector multiply and subtract from another	*
*  su3_vector 								*
*  C  <-  C - A*B 							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
/* su3_matrix times su3_vector multiply and subtract from another su3_vector */
/* c  <-  A*b-c */
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c ){
register int i,j;
register complex x,y;
    for(i=0;i<3;i++){
	x.real=x.imag=0.0;
	for(j=0;j<3;j++){
	    CMUL( a->e[i][j] , b->c[j] , y )
	    CSUM( x , y );
	}
	c->c[i].real -= x.real;
	c->c[i].imag -= x.imag;
    }
}

#else
#ifdef NATIVEDOUBLE
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c ){

  register double c0r,c0i,c1r,c1i,c2r,c2i;
  register double br,bi,a0,a1,a2;

  c0r = c->c[0].real;
  c0i = c->c[0].imag;
  c1r = c->c[1].real;
  c1i = c->c[1].imag;
  c2r = c->c[2].real;
  c2i = c->c[2].imag;

  br=b->c[0].real;    bi=b->c[0].imag;
  a0=a->e[0][0].real;
  a1=a->e[1][0].real;
  a2=a->e[2][0].real;

  c0r -= a0*br;
  c1r -= a1*br;
  c2r -= a2*br;
  c0i -= a0*bi;
  c1i -= a1*bi;
  c2i -= a2*bi;

  a0=a->e[0][0].imag;
  a1=a->e[1][0].imag;
  a2=a->e[2][0].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  br=b->c[1].real;    bi=b->c[1].imag;
  a0=a->e[0][1].real;
  a1=a->e[1][1].real;
  a2=a->e[2][1].real;

  c0r -= a0*br;
  c1r -= a1*br;
  c2r -= a2*br;
  c0i -= a0*bi;
  c1i -= a1*bi;
  c2i -= a2*bi;

  a0=a->e[0][1].imag;
  a1=a->e[1][1].imag;
  a2=a->e[2][1].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  br=b->c[2].real;    bi=b->c[2].imag;
  a0=a->e[0][2].real;
  a1=a->e[1][2].real;
  a2=a->e[2][2].real;

  c0r -= a0*br;
  c1r -= a1*br;
  c2r -= a2*br;
  c0i -= a0*bi;
  c1i -= a1*bi;
  c2i -= a2*bi;

  a0=a->e[0][2].imag;
  a1=a->e[1][2].imag;
  a2=a->e[2][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  c->c[0].real = c0r;
  c->c[0].imag = c0i;
  c->c[1].real = c1r;
  c->c[1].imag = c1i;
  c->c[2].real = c2r;
  c->c[2].imag = c2i;

}

#else
void mult_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b, su3_vector *c ){
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

	c->c[i].real -= cr;
	c->c[i].imag -= ci;
    }
}
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#ifdef FAST" */
