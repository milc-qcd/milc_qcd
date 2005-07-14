/********************** dcsqrt.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* double precision complex square root */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dcsqrt( double_complex *z ){
double_complex c;
double theta,r;
    r = sqrt(hypot(z->real,z->imag));
    theta = 0.5*atan2(z->imag,z->real);
    c = dce_itheta(theta);
    c.real *=r; c.imag *= r;
    return(c);
}

