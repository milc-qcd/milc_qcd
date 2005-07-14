/********************** dconjg.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* double precision complex conjugate */
#include "../include/config.h"
#include "../include/complex.h"

double_complex dconjg(  double_complex *a ){
    double_complex c;
    c.real = (*a).real;
    c.imag = -(*a).imag;
    return(c);
}
