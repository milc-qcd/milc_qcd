/********************** csqrt.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* complex square root */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex csqrt( complex *z ){
complex c;
Real theta,r;
    r = sqrt(hypot(z->real,z->imag));
    theta = 0.5*atan2(z->imag,z->real);
    c = ce_itheta(theta);
    c.real *=r; c.imag *= r;
    return(c);
}
