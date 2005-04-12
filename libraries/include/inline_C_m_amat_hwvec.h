/**************  m_amat_hwvec.c  (in su3.a) **********************
*									*
*  void mult_adj_su3_mat_hwvec( su3_matrix *mat,			*
*	half_wilson_vector *src,*dest )					*
*  multiply a Wilson half-vector by the adjoint of a matrix		*
*/


#define _inline_C_mult_adj_su3_mat_hwvec( mat, src, dest ) {\
\
  Real _a0r,_a0i,_a1r,_a1i,_a2r,_a2i;\
  Real _b0r,_b0i,_b1r,_b1i,_b2r,_b2i;\
\
\
  _a0r=(mat)->e[0][0].real;   _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[0].c[0].real; _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[1][0].real;   _a1i=(mat)->e[1][0].imag;\
  _b1r=(src)->h[0].c[1].real; _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[2][0].real;   _a2i=(mat)->e[2][0].imag;\
  _b2r=(src)->h[0].c[2].real; _b2i=(src)->h[0].c[2].imag;\
  \
  (dest)->h[0].c[0].real = _a0r*_b0r + _a0i*_b0i + _a1r*_b1r + _a1i*_b1i + _a2r*_b2r + _a2i*_b2i;\
  (dest)->h[0].c[0].imag = _a0r*_b0i - _a0i*_b0r + _a1r*_b1i - _a1i*_b1r + _a2r*_b2i - _a2i*_b2r;\
  \
  _a0r=(mat)->e[0][1].real;   _a0i=(mat)->e[0][1].imag;\
  _b0r=(src)->h[0].c[0].real; _b0i=(src)->h[0].c[0].imag;  \
  _a1r=(mat)->e[1][1].real;   _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[0].c[1].real; _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[2][1].real;   _a2i=(mat)->e[2][1].imag;\
  _b2r=(src)->h[0].c[2].real; _b2i=(src)->h[0].c[2].imag;\
  \
  (dest)->h[0].c[1].real = _a0r*_b0r + _a0i*_b0i + _a1r*_b1r + _a1i*_b1i + _a2r*_b2r + _a2i*_b2i;\
  (dest)->h[0].c[1].imag = _a0r*_b0i - _a0i*_b0r + _a1r*_b1i - _a1i*_b1r + _a2r*_b2i - _a2i*_b2r;\
  \
  _a0r=(mat)->e[0][2].real;   _a0i=(mat)->e[0][2].imag;\
  _b0r=(src)->h[0].c[0].real; _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[1][2].real;   _a1i=(mat)->e[1][2].imag;\
  _b1r=(src)->h[0].c[1].real; _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[2][2].real;   _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[0].c[2].real; _b2i=(src)->h[0].c[2].imag;\
  \
  (dest)->h[0].c[2].real = _a0r*_b0r + _a0i*_b0i + _a1r*_b1r + _a1i*_b1i + _a2r*_b2r + _a2i*_b2i;\
  (dest)->h[0].c[2].imag = _a0r*_b0i - _a0i*_b0r + _a1r*_b1i - _a1i*_b1r + _a2r*_b2i - _a2i*_b2r;\
\
\
/*    mult_adj_su3_(mat)_vec((mat), &((src)->h[1]), &((dest)->h[1]) ); */\
\
  _a0r=(mat)->e[0][0].real;   _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[1].c[0].real; _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[1][0].real;   _a1i=(mat)->e[1][0].imag;\
  _b1r=(src)->h[1].c[1].real; _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[2][0].real;   _a2i=(mat)->e[2][0].imag;\
  _b2r=(src)->h[1].c[2].real; _b2i=(src)->h[1].c[2].imag;\
  \
  (dest)->h[1].c[0].real = _a0r*_b0r + _a0i*_b0i + _a1r*_b1r + _a1i*_b1i + _a2r*_b2r + _a2i*_b2i;\
  (dest)->h[1].c[0].imag = _a0r*_b0i - _a0i*_b0r + _a1r*_b1i - _a1i*_b1r + _a2r*_b2i - _a2i*_b2r;\
  \
  _a0r=(mat)->e[0][1].real;   _a0i=(mat)->e[0][1].imag;\
  _b0r=(src)->h[1].c[0].real; _b0i=(src)->h[1].c[0].imag;  \
  _a1r=(mat)->e[1][1].real;   _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[1].c[1].real; _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[2][1].real;   _a2i=(mat)->e[2][1].imag;\
  _b2r=(src)->h[1].c[2].real; _b2i=(src)->h[1].c[2].imag;\
  \
  (dest)->h[1].c[1].real = _a0r*_b0r + _a0i*_b0i + _a1r*_b1r + _a1i*_b1i + _a2r*_b2r + _a2i*_b2i;\
  (dest)->h[1].c[1].imag = _a0r*_b0i - _a0i*_b0r + _a1r*_b1i - _a1i*_b1r + _a2r*_b2i - _a2i*_b2r;\
  \
  _a0r=(mat)->e[0][2].real;   _a0i=(mat)->e[0][2].imag;\
  _b0r=(src)->h[1].c[0].real; _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[1][2].real;   _a1i=(mat)->e[1][2].imag;\
  _b1r=(src)->h[1].c[1].real; _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[2][2].real;   _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[1].c[2].real; _b2i=(src)->h[1].c[2].imag;\
  \
  (dest)->h[1].c[2].real = _a0r*_b0r + _a0i*_b0i + _a1r*_b1r + _a1i*_b1i + _a2r*_b2r + _a2i*_b2i;\
  (dest)->h[1].c[2].imag = _a0r*_b0i - _a0i*_b0r + _a1r*_b1i - _a1i*_b1r + _a2r*_b2i - _a2i*_b2r;\
\
}
