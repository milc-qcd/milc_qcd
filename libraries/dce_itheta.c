/********************** dce_itheta.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* double complex exp( i*theta ) */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dce_itheta( double theta ){
    double_complex c;
    c.real = (double)cos( (double)theta );
    c.imag = (double)sin( (double)theta );
    /* there must be a more efficient way */
    return( c );
}
