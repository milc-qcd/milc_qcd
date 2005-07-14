/********************** dclog.c (in complex.a) **********************/
/* MIMD version 7 */
/* Subroutines for operations on complex numbers */
/* double complex logarithm */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dclog(  double_complex *a ){
    double_complex c;
    c.real = 0.5*(double)log((double)((*a).real*(*a).real+(*a).imag*(*a).imag));
    c.imag = (double)atan2( (double)(*a).imag, (double)(*a).real );
    return(c);
}
