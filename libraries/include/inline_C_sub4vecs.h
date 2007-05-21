/* void sub_four_su3_vecs( su3_vector *a,*b1,*b2,*b3,*b4) 		*
 * A  <-  A - B1 - B2 - B3 - B4						*/
#define _inline_C_sub_four_su3vecs( a, b1, b2, b3, b4 ) \
  { register su3_vector *aaa,*bbb1,*bbb2,*bbb3,*bbb4;	\
    aaa=(a); bbb1=(b1); bbb2=(b2); bbb3=(b3); bbb4=(b4); \
    (aaa)->c[0].real -= (bbb1)->c[0].real; \
    (aaa)->c[0].imag -= (bbb1)->c[0].imag; \
    (aaa)->c[1].real -= (bbb1)->c[1].real; \
    (aaa)->c[1].imag -= (bbb1)->c[1].imag; \
    (aaa)->c[2].real -= (bbb1)->c[2].real; \
    (aaa)->c[2].imag -= (bbb1)->c[2].imag; \
    (aaa)->c[0].real -= (bbb2)->c[0].real; \
    (aaa)->c[0].imag -= (bbb2)->c[0].imag; \
    (aaa)->c[1].real -= (bbb2)->c[1].real; \
    (aaa)->c[1].imag -= (bbb2)->c[1].imag; \
    (aaa)->c[2].real -= (bbb2)->c[2].real; \
    (aaa)->c[2].imag -= (bbb2)->c[2].imag; \
    (aaa)->c[0].real -= (bbb3)->c[0].real; \
    (aaa)->c[0].imag -= (bbb3)->c[0].imag; \
    (aaa)->c[1].real -= (bbb3)->c[1].real; \
    (aaa)->c[1].imag -= (bbb3)->c[1].imag; \
    (aaa)->c[2].real -= (bbb3)->c[2].real; \
    (aaa)->c[2].imag -= (bbb3)->c[2].imag; \
    (aaa)->c[0].real -= (bbb4)->c[0].real; \
    (aaa)->c[0].imag -= (bbb4)->c[0].imag; \
    (aaa)->c[1].real -= (bbb4)->c[1].real; \
    (aaa)->c[1].imag -= (bbb4)->c[1].imag; \
    (aaa)->c[2].real -= (bbb4)->c[2].real; \
    (aaa)->c[2].imag -= (bbb4)->c[2].imag; \
}
