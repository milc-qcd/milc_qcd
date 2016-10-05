/* void su3mat_copy( su3_matrix *a, su3_matrix *b )			*
* Copy an su3 matrix:  B <- A   						*
*/
#define _inline_C_su3mat_copy(a,b) \
  (b)->e[0][0].real = (a)->e[0][0].real; \
  (b)->e[0][0].imag = (a)->e[0][0].imag; \
  (b)->e[0][1].real = (a)->e[0][1].real; \
  (b)->e[0][1].imag = (a)->e[0][1].imag; \
  (b)->e[0][2].real = (a)->e[0][2].real; \
  (b)->e[0][2].imag = (a)->e[0][2].imag; \
\
  (b)->e[1][0].real = (a)->e[1][0].real; \
  (b)->e[1][0].imag = (a)->e[1][0].imag; \
  (b)->e[1][1].real = (a)->e[1][1].real; \
  (b)->e[1][1].imag = (a)->e[1][1].imag; \
  (b)->e[1][2].real = (a)->e[1][2].real; \
  (b)->e[1][2].imag = (a)->e[1][2].imag; \
\
  (b)->e[2][0].real = (a)->e[2][0].real; \
  (b)->e[2][0].imag = (a)->e[2][0].imag; \
  (b)->e[2][1].real = (a)->e[2][1].real; \
  (b)->e[2][1].imag = (a)->e[2][1].imag; \
  (b)->e[2][2].real = (a)->e[2][2].real; \
  (b)->e[2][2].imag = (a)->e[2][2].imag;


