/* void add_su3_matrix( su3_matrix *a, su3_matrix *b,	su3_matrix *c) * 
* C <- A + B,   A,B and C matrices 					*
*/
#define _inline_C_add_su3_matrix(aa,bb,cc) \
  { register su3_matrix *aaa,*bbb,*ccc; register int iii; \
  aaa=(aa); bbb=(bb); ccc=(cc); \
    for( iii=0; iii<3; iii++ ){ \
  (ccc)->e[iii][0].real = (aaa)->e[iii][0].real + (bbb)->e[iii][0].real; \
  (ccc)->e[iii][0].imag = (aaa)->e[iii][0].imag + (bbb)->e[iii][0].imag; \
  (ccc)->e[iii][1].real = (aaa)->e[iii][1].real + (bbb)->e[iii][1].real; \
  (ccc)->e[iii][1].imag = (aaa)->e[iii][1].imag + (bbb)->e[iii][1].imag; \
  (ccc)->e[iii][2].real = (aaa)->e[iii][2].real + (bbb)->e[iii][2].real; \
  (ccc)->e[iii][2].imag = (aaa)->e[iii][2].imag + (bbb)->e[iii][2].imag; \
    } \
  }

#if 0
// For mysterious reasons, the following version doesn't work in all use cases */

/* void add_su3_matrix( su3_matrix *a, su3_matrix *b,	su3_matrix *c) * 
* C <- A + B,   A,B and C matrices 					*
*/
#define _inline_C_add_su3_matrix(aa,bb,cc) \
  (cc)->e[0][0].real = (aa)->e[0][0].real + (bb)->e[0][0].real; \
  (cc)->e[0][0].imag = (aa)->e[0][0].imag + (bb)->e[0][0].imag; \
  (cc)->e[0][1].real = (aa)->e[0][1].real + (bb)->e[0][1].real; \
  (cc)->e[0][1].imag = (aa)->e[0][1].imag + (bb)->e[0][1].imag; \
  (cc)->e[0][2].real = (aa)->e[0][2].real + (bb)->e[0][2].real; \
  (cc)->e[0][2].imag = (aa)->e[0][2].imag + (bb)->e[0][2].imag; \
\
  (cc)->e[1][0].real = (aa)->e[1][0].real + (bb)->e[1][0].real; \
  (cc)->e[1][0].imag = (aa)->e[1][0].imag + (bb)->e[1][0].imag; \
  (cc)->e[1][1].real = (aa)->e[1][1].real + (bb)->e[1][1].real; \
  (cc)->e[1][1].imag = (aa)->e[1][1].imag + (bb)->e[1][1].imag; \
  (cc)->e[1][2].real = (aa)->e[1][2].real + (bb)->e[1][2].real; \
  (cc)->e[1][2].imag = (aa)->e[1][2].imag + (bb)->e[1][2].imag; \
\
  (cc)->e[2][0].real = (aa)->e[2][0].real + (bb)->e[2][0].real; \
  (cc)->e[2][0].imag = (aa)->e[2][0].imag + (bb)->e[2][0].imag; \
  (cc)->e[2][1].real = (aa)->e[2][1].real + (bb)->e[2][1].real; \
  (cc)->e[2][1].imag = (aa)->e[2][1].imag + (bb)->e[2][1].imag; \
  (cc)->e[2][2].real = (aa)->e[2][2].real + (bb)->e[2][2].real; \
  (cc)->e[2][2].imag = (aa)->e[2][2].imag + (bb)->e[2][2].imag;


#endif
