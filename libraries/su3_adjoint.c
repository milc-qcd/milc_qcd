/******************  su3_adjoint.c  (in su3.a) **************************
*									*
* void su3_adjoint( su3_matrix *a, su3_matrix *b )			*
* B  <- A_adjoint,  adjoint of an SU3 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* adjoint of an SU3 matrix */
void su3_adjoint( su3_matrix *a, su3_matrix *b ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	CONJG( a->e[j][i], b->e[i][j] );
    }
}
