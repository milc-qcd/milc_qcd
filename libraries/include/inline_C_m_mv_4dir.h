/****************  m_mv_s_4dir.c  (in su3.a) *****************************
*									*
* void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b[0123],*c )*
* Multiply the elements of an array of four su3_matrices by the		*
* four su3_vectors, and add the results to				*
* produce a single su3_vector.						*
* C  <-  A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3]			*
*/
#define _inline_C_mult_su3_mat_vec_sum_4dir(  aa , bb0, bb1, bb2, bb3, cc ) \
{ \
  register int n; \
  register Real c0r,c0i,c1r,c1i,c2r,c2i; \
  register Real br,bi,a0,a1,a2; \
  register su3_matrix *mat; \
  register su3_vector *b; \
  c0r = c0i = c1r = c1i = c2r = c2i = 0.0; \
  mat = (aa); \
  for(n=0;n<4;n++,mat++){ \
  switch(n){ \
    case(0): b=(bb0); break; \
    case(1): b=(bb1); break; \
    case(2): b=(bb2); break; \
    case(3): b=(bb3); break; \
    default: b=NULL; \
  } \
  br=b->c[0].real;    bi=b->c[0].imag; \
  a0=mat->e[0][0].real; \
  a1=mat->e[1][0].real; \
  a2=mat->e[2][0].real; \
  c0r += a0*br; \
  c1r += a1*br; \
  c2r += a2*br; \
  c0i += a0*bi; \
  c1i += a1*bi; \
  c2i += a2*bi; \
  a0=mat->e[0][0].imag; \
  a1=mat->e[1][0].imag; \
  a2=mat->e[2][0].imag; \
  c0r -= a0*bi; \
  c1r -= a1*bi; \
  c2r -= a2*bi; \
  c0i += a0*br; \
  c1i += a1*br; \
  c2i += a2*br; \
  br=b->c[1].real;    bi=b->c[1].imag; \
  a0=mat->e[0][1].real; \
  a1=mat->e[1][1].real; \
  a2=mat->e[2][1].real; \
  c0r += a0*br; \
  c1r += a1*br; \
  c2r += a2*br; \
  c0i += a0*bi; \
  c1i += a1*bi; \
  c2i += a2*bi; \
  a0=mat->e[0][1].imag; \
  a1=mat->e[1][1].imag; \
  a2=mat->e[2][1].imag; \
  c0r -= a0*bi; \
  c1r -= a1*bi; \
  c2r -= a2*bi; \
  c0i += a0*br; \
  c1i += a1*br; \
  c2i += a2*br; \
  br=b->c[2].real;    bi=b->c[2].imag; \
  a0=mat->e[0][2].real; \
  a1=mat->e[1][2].real; \
  a2=mat->e[2][2].real; \
  c0r += a0*br; \
  c1r += a1*br; \
  c2r += a2*br; \
  c0i += a0*bi; \
  c1i += a1*bi; \
  c2i += a2*bi; \
  a0=mat->e[0][2].imag; \
  a1=mat->e[1][2].imag; \
  a2=mat->e[2][2].imag; \
  c0r -= a0*bi; \
  c1r -= a1*bi; \
  c2r -= a2*bi; \
  c0i += a0*br; \
  c1i += a1*br; \
  c2i += a2*br; \
  } \
  (cc)->c[0].real = c0r; \
  (cc)->c[0].imag = c0i; \
  (cc)->c[1].real = c1r; \
  (cc)->c[1].imag = c1i; \
  (cc)->c[2].real = c2r; \
  (cc)->c[2].imag = c2i; \
}
