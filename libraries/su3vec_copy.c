/*****************  su3vec_copy.c  (in su3.a) ***************************
*									*
* void su3vec_copy( su3_vector *a, su3_vector *b )			*
* Copy an su3 vector:  B <- A   					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Copy a su3 vector:  b <- a   */
void su3vec_copy( su3_vector *a, su3_vector *b ){
register int i;
    for(i=0;i<3;i++){
	b->c[i].real = a->c[i].real;
	b->c[i].imag = a->c[i].imag;
    }
}
