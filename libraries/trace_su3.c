/*******************  trace_su3.c  (in su3.a) ***************************
*									*
* complex trace_su3(a) su3_matrix *a;					*
* return complex trace of an SU3 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Complex trace of an SU3 matrix */
complex trace_su3( su3_matrix *a ) {
register complex t1,t2;
    CADD(a->e[0][0],a->e[1][1],t1);
    CADD(t1,a->e[2][2],t2);
    return(t2);
}
