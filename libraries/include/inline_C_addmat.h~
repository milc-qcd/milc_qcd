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

