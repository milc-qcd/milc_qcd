/*******************  submat.c  (in su3.a) ******************************
*									*
* void sub_su3_matrix(a,b,c) su3_matrix *a,*b,*c;			*
* subtract su3 matrices:  C  <- A - B 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* subtract su3 matrices */
void sub_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c ) {
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	CSUB( a->e[i][j], b->e[i][j], c->e[i][j] );
    }
}
