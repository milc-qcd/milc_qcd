/**************  m_mat_hwvec.c  (in su3.a) ***********************
*									*
* void mult_su3_mat_hwvec(su3_matrix *mat,				*
*	half_wilson_vector *src,*dest)					*
*  multiply a Wilson half-vector by a matrix				*
*  dest  <-  mat*src							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST

void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest ){
    mult_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) );
    mult_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) );
}

#else  /* Fast version */


void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest ){

#ifdef NATIVEDOUBLE
  register double a0r,a0i,a1r,a1i,a2r,a2i;
  register double b0r,b0i,b1r,b1i,b2r,b2i;
#else
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
#endif
  
/*    mult_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) ); */

  a0r=mat->e[0][0].real;    a0i=mat->e[0][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[0][1].real;    a1i=mat->e[0][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[0][2].real;    a2i=mat->e[0][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=mat->e[1][0].real;    a0i=mat->e[1][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][1].real;    a1i=mat->e[1][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[1][2].real;    a2i=mat->e[1][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=mat->e[2][0].real;    a0i=mat->e[2][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[2][1].real;    a1i=mat->e[2][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][2].real;    a2i=mat->e[2][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

/*    mult_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) ); */

  a0r=mat->e[0][0].real;    a0i=mat->e[0][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[0][1].real;    a1i=mat->e[0][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[0][2].real;    a2i=mat->e[0][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=mat->e[1][0].real;    a0i=mat->e[1][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][1].real;    a1i=mat->e[1][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[1][2].real;    a2i=mat->e[1][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=mat->e[2][0].real;    a0i=mat->e[2][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[2][1].real;    a1i=mat->e[2][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][2].real;    a2i=mat->e[2][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

}

#endif /* "ifndef FAST" */
