/* void scalar_mult_add_su3_vector( su3_vector *a, su3_vector *b,	*
*	Real s, su3_vector *c)						*
* C <- A + s*B,   A,B and C vectors 					*
*/
#define _inline_C_scalar_mult_add_su3_vector(aa,bb,ss,cc) \
  { register su3_vector *aaa,*bbb,*ccc; register Real sss; \
  aaa=(aa); bbb=(bb); ccc=(cc); sss=(ss); \
  (ccc)->c[0].real = (aaa)->c[0].real + (sss)*(bbb)->c[0].real; \
  (ccc)->c[0].imag = (aaa)->c[0].imag + (sss)*(bbb)->c[0].imag; \
  (ccc)->c[1].real = (aaa)->c[1].real + (sss)*(bbb)->c[1].real; \
  (ccc)->c[1].imag = (aaa)->c[1].imag + (sss)*(bbb)->c[1].imag; \
  (ccc)->c[2].real = (aaa)->c[2].real + (sss)*(bbb)->c[2].real; \
  (ccc)->c[2].imag = (aaa)->c[2].imag + (sss)*(bbb)->c[2].imag; \
  }
