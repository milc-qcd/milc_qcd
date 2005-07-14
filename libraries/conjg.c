/********************** conjg.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* complex conjugate */
#include "../include/config.h"
#include "../include/complex.h"

complex conjg( complex *a ){
    complex c;
    c.real = (*a).real;
    c.imag = -(*a).imag;
    return(c);
}
