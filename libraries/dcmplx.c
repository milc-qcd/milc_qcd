/********************** dcmplx.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* make a double complex number from two double precision reals */
#include "../include/config.h"
#include "../include/complex.h"

double_complex dcmplx( double x, double y ){
    double_complex c;
    c.real = x; c.imag = y;
    return(c);
}
