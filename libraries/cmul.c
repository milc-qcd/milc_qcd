/********************** cmul.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* multiply two complex numbers */
#include "../include/config.h"
#include "../include/complex.h"

complex cmul( complex *a, complex *b ) {
    complex c;
    c.real = (*a).real * (*b).real - (*a).imag * (*b).imag;
    c.imag = (*a).imag * (*b).real + (*a).real * (*b).imag;
    return(c);
}
