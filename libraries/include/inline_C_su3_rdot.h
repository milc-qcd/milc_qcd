/* Real su3_rdot( su3_vector *a, su3_vector *b )			*
* return real part of dot product of two su3_vectors			*
*/

#define _inline_C_su3_rdot( aa, bb ) \
  ( (aa)->c[0].real * (bb)->c[0].real + \
    (aa)->c[0].imag * (bb)->c[0].imag + \
    (aa)->c[1].real * (bb)->c[1].real + \
    (aa)->c[1].imag * (bb)->c[1].imag + \
    (aa)->c[2].real * (bb)->c[2].real + \
    (aa)->c[2].imag * (bb)->c[2].imag )
