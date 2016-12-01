/* void sub_su3_matrix(a,b,c) su3_matrix *a,*b,*c;			*
* subtract su3 matrices:  C  <- A - B 					*
*/
#define _inline_C_sub_su3_matrix(a,b,c) \
  (c)->e[0][0].real = (a)->e[0][0].real - (b)->e[0][0].real; \
  (c)->e[0][0].imag = (a)->e[0][0].imag - (b)->e[0][0].imag; \
  (c)->e[0][1].real = (a)->e[0][1].real - (b)->e[0][1].real; \
  (c)->e[0][1].imag = (a)->e[0][1].imag - (b)->e[0][1].imag; \
  (c)->e[0][2].real = (a)->e[0][2].real - (b)->e[0][2].real; \
  (c)->e[0][2].imag = (a)->e[0][2].imag - (b)->e[0][2].imag; \
\
  (c)->e[1][0].real = (a)->e[1][0].real - (b)->e[1][0].real; \
  (c)->e[1][0].imag = (a)->e[1][0].imag - (b)->e[1][0].imag; \
  (c)->e[1][1].real = (a)->e[1][1].real - (b)->e[1][1].real; \
  (c)->e[1][1].imag = (a)->e[1][1].imag - (b)->e[1][1].imag; \
  (c)->e[1][2].real = (a)->e[1][2].real - (b)->e[1][2].real; \
  (c)->e[1][2].imag = (a)->e[1][2].imag - (b)->e[1][2].imag; \
\
  (c)->e[2][0].real = (a)->e[2][0].real - (b)->e[2][0].real; \
  (c)->e[2][0].imag = (a)->e[2][0].imag - (b)->e[2][0].imag; \
  (c)->e[2][1].real = (a)->e[2][1].real - (b)->e[2][1].real; \
  (c)->e[2][1].imag = (a)->e[2][1].imag - (b)->e[2][1].imag; \
  (c)->e[2][2].real = (a)->e[2][2].real - (b)->e[2][2].real; \
  (c)->e[2][2].imag = (a)->e[2][2].imag - (b)->e[2][2].imag;

