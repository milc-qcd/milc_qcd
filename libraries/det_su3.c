/******************  det_su3.c  (in su3.a) ******************************
*									*
* complex det_su3( su3_matrix *a )					*
* Complex determinant of an SU3 matrix 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* FIX THIS - more efficient to take cross product of first two
   rows, dot with third. */
complex det_su3( su3_matrix *a ) {
register complex cc,dd,sum;
    CMUL(a->e[0][0],a->e[1][1],cc);
    CMUL(cc,a->e[2][2],sum);
    CMUL(a->e[0][0],a->e[1][2],cc);
    CMUL(cc,a->e[2][1],dd);
    CSUB(sum,dd,sum);
    CMUL(a->e[0][1],a->e[1][2],cc);
    CMUL(cc,a->e[2][0],dd);
    CADD(sum,dd,sum);
    CMUL(a->e[0][1],a->e[1][0],cc);
    CMUL(cc,a->e[2][2],dd);
    CSUB(sum,dd,sum);
    CMUL(a->e[0][2],a->e[1][0],cc);
    CMUL(cc,a->e[2][1],dd);
    CADD(sum,dd,sum);
    CMUL(a->e[0][2],a->e[1][1],cc);
    CMUL(cc,a->e[2][0],dd);
    CSUB(sum,dd,sum);
    return(sum);
}
