
/* float magsq_su3vec( su3_vector *a )                                   *
* return squared magnitude of an SU3 vector */

#define _inline_C_magsq_su3vec( aa ) \
  ( ILvecpt = (aa), \
    (ILvecpt)->c[0].real * (ILvecpt)->c[0].real + \
    (ILvecpt)->c[0].imag * (ILvecpt)->c[0].imag + \
    (ILvecpt)->c[1].real * (ILvecpt)->c[1].real + \
    (ILvecpt)->c[1].imag * (ILvecpt)->c[1].imag + \
    (ILvecpt)->c[2].real * (ILvecpt)->c[2].real + \
    (ILvecpt)->c[2].imag * (ILvecpt)->c[2].imag )

