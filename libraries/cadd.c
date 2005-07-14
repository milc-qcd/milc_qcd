/********************** cadd.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* add two complex numbers */
#include "../include/config.h"
#include "../include/complex.h"

complex cadd( complex *a, complex *b ) {
    complex c;
    c.real = (*a).real + (*b).real;
    c.imag = (*a).imag + (*b).imag;
    return(c);
}
