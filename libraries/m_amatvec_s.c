/*******************  m_amatvec_s.c  (in su3.a) *************************
*									*
*  void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b,*c )	*
*  adjoint matrix times vector multiply and add to another vector 	*
*  C  <-  C + A_adjoint*B						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c ){
register int i,j;
register complex x,y,z;
    for(i=0;i<3;i++){
	x.real=x.imag=0.0;
	for(j=0;j<3;j++){
	    CONJG( a->e[j][i], z );
	    CMUL( z , b->c[j], y )
	    CSUM( x , y );
	}
	c->c[i].real += x.real;
	c->c[i].imag += x.imag;
    }
}

#else
void mult_adj_su3_mat_vec_sum( su3_matrix *a, su3_vector *b, su3_vector *c ){

#ifdef NATIVEDOUBLE
  register double c0r,c0i,c1r,c1i,c2r,c2i;
  register double br,bi,a0,a1,a2;
#else
  register Real c0r,c0i,c1r,c1i,c2r,c2i;
  register Real br,bi,a0,a1,a2;
#endif

  br=b->c[0].real;    bi=b->c[0].imag;
  a0=a->e[0][0].real;
  a1=a->e[0][1].real;
  a2=a->e[0][2].real;

  c0r = a0*br;
  c1r = a1*br;
  c2r = a2*br;
  c0i = a0*bi;
  c1i = a1*bi;
  c2i = a2*bi;

  a0=a->e[0][0].imag;
  a1=a->e[0][1].imag;
  a2=a->e[0][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  br=b->c[1].real;    bi=b->c[1].imag;
  a0=a->e[1][0].real;
  a1=a->e[1][1].real;
  a2=a->e[1][2].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=a->e[1][0].imag;
  a1=a->e[1][1].imag;
  a2=a->e[1][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  br=b->c[2].real;    bi=b->c[2].imag;
  a0=a->e[2][0].real;
  a1=a->e[2][1].real;
  a2=a->e[2][2].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=a->e[2][0].imag;
  a1=a->e[2][1].imag;
  a2=a->e[2][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  c->c[0].real += c0r;
  c->c[0].imag += c0i;
  c->c[1].real += c1r;
  c->c[1].imag += c1i;
  c->c[2].real += c2r;
  c->c[2].imag += c2i;
}

#endif	/* End of "#ifdef FAST" */
