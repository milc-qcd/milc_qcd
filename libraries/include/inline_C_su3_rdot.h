/* float su3_rdot( su3_vector *a, su3_vector *b )			*
* return real part of dot product of two su3_vectors			*
*/

#define _inline_C_su3_rdot( aa, bb ) \
  ( ILvecpt = (aa), ILvecpt2 = (bb), \
    (ILvecpt)->c[0].real * (ILvecpt2)->c[0].real + \
    (ILvecpt)->c[0].imag * (ILvecpt2)->c[0].imag + \
    (ILvecpt)->c[1].real * (ILvecpt2)->c[1].real + \
    (ILvecpt)->c[1].imag * (ILvecpt2)->c[1].imag + \
    (ILvecpt)->c[2].real * (ILvecpt2)->c[2].real + \
    (ILvecpt)->c[2].imag * (ILvecpt2)->c[2].imag )
