/*****************  su3mat_copy.c  (in su3.a) ***************************
*									*
* void su3mat_copy( su3_matrix *a, su3_matrix *b )			*
* Copy an su3 matrix:  B <- A   						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Copy a su3 matrix:  b <- a   */
void su3mat_copy( su3_matrix *a, su3_matrix *b ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	b->e[i][j].real = a->e[i][j].real;
	b->e[i][j].imag = a->e[i][j].imag;
    }
}
