/***************** m_amv_4vec.c in su3.a    *****************************
*									*
*  void mult_adj_su3_mat_4vec( su3_matrix *mat,			        *
*  su3_vector *src, su3_vector *dest0, *dest1, *dest2, *dest3 )		*
*  Multiply an su3_vector by an array of four adjoint su3_matrices,	*
*  result in four SEPARATE su3_vectors.		           		*
*  desti  <-  A_adjoint[i] * src					*
*  See also m_amv_4dir.c for the case desti = dest[i]                   *
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST
void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3 ) {
    mult_adj_su3_mat_vec( mat+0, src, dest0 );
    mult_adj_su3_mat_vec( mat+1, src, dest1 );
    mult_adj_su3_mat_vec( mat+2, src, dest2 );
    mult_adj_su3_mat_vec( mat+3, src, dest3 );
}

#else
/* Fast code, with subroutines inlined */

void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3  ){
  register int n;
#ifdef NATIVEDOUBLE
  register double c0r,c0i,c1r,c1i,c2r,c2i;
  register double br,bi,a0,a1,a2;
#else
  register Real c0r,c0i,c1r,c1i,c2r,c2i;
  register Real br,bi,a0,a1,a2;
#endif
  register su3_matrix *a;
  register su3_vector *b,*c;
  su3_vector *cc[4] ;
  
  cc[0] = dest0 ; cc[1] = dest1 ;
  cc[2] = dest2 ; cc[3] = dest3 ;

  a = mat ; c=dest0 ; b = src;
  for(n=0;n<4;n++,a++,c=cc[n]){

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

  c->c[0].real = c0r;
  c->c[0].imag = c0i;
  c->c[1].real = c1r;
  c->c[1].imag = c1i;
  c->c[2].real = c2r;
  c->c[2].imag = c2i;
  }
}
#endif	/* End of "#ifndef FAST" */
