/********************** cexp.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* complex exponential */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex cexp( complex *a ){
    complex c;
    Real mag;
    mag = (Real)exp( (double)(*a).real );
    c.real = mag*(Real)cos( (double)(*a).imag );
    c.imag = mag*(Real)sin( (double)(*a).imag );
    return(c);
}
