/********************** cdiv.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* Divide two complex numbers */
#include "../include/config.h"
#include "../include/complex.h"

complex cdiv( complex *a, complex *b ) {
    complex c;
    Real scale;
    scale = 1.0/((*b).real*(*b).real+(*b).imag*(*b).imag);
    c.real = scale*((*a).real*(*b).real + (*a).imag*(*b).imag);
    c.imag = scale*((*a).imag*(*b).real - (*a).real*(*b).imag);
    return(c);
}
