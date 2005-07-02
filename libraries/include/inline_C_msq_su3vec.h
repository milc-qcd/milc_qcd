
/* Real magsq_su3vec( su3_vector *a )                                   *
* return squared magnitude of an SU3 vector */

#define _inline_C_magsq_su3vec( aa ) \
  ( (aa)->c[0].real * (aa)->c[0].real + \
    (aa)->c[0].imag * (aa)->c[0].imag + \
    (aa)->c[1].real * (aa)->c[1].real + \
    (aa)->c[1].imag * (aa)->c[1].imag + \
    (aa)->c[2].real * (aa)->c[2].real + \
    (aa)->c[2].imag * (aa)->c[2].imag )

