/* void add_su3_vector(a,b,c) su3_vector *a,*b,*c;                       *
* add su3 vectors:  C <-  A - B                                    */
#define _inline_C_add_su3_vector( aa,bb,cc ) \
    { register su3_vector *aaa,*bbb,*ccc; \
    aaa=(aa); bbb=(bb); ccc=(cc); \
    (ccc)->c[0].real = (aaa)->c[0].real + (bbb)->c[0].real; \
    (ccc)->c[0].imag = (aaa)->c[0].imag + (bbb)->c[0].imag; \
    (ccc)->c[1].real = (aaa)->c[1].real + (bbb)->c[1].real; \
    (ccc)->c[1].imag = (aaa)->c[1].imag + (bbb)->c[1].imag; \
    (ccc)->c[2].real = (aaa)->c[2].real + (bbb)->c[2].real; \
    (ccc)->c[2].imag = (aaa)->c[2].imag + (bbb)->c[2].imag; \
    }
