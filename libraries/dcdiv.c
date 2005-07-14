/********************** dcdiv.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* double complex divide */
#include "../include/config.h"
#include "../include/complex.h"

double_complex dcdiv( double_complex *a, double_complex *b ){
    double_complex c;
    double scale;
    scale = 1.0/((*b).real*(*b).real+(*b).imag*(*b).imag);
    c.real = scale*((*a).real*(*b).real + (*a).imag*(*b).imag);
    c.imag = scale*((*a).imag*(*b).real - (*a).real*(*b).imag);
    return(c);
}
