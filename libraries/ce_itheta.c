/********************** ce_itheta.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* exp( i*theta ) */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex ce_itheta( Real theta ){
    complex c;
    c.real = (Real)cos( (double)theta );
    c.imag = (Real)sin( (double)theta );
    /* there must be a more efficient way */
    return( c );
}
