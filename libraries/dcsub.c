/********************** dcsub.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* double complex subtract */
#include "../include/config.h"
#include "../include/complex.h"

double_complex dcsub( double_complex *a, double_complex *b ){
    double_complex c;
    c.real = (*a).real - (*b).real;
    c.imag = (*a).imag - (*b).imag;
    return(c);
}
