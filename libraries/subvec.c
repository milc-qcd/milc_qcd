/*********************  subvec.c  (in su3.a) ****************************
*									*
* void sub_su3_vector(a,b,c) su3_vector *a,*b,*c;			*
* subtract su3 vectors:  C <-  A - B 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* subtract su3 vectors */
void sub_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c ){
register int i;
    for(i=0;i<3;i++){
	CSUB( a->c[i], b->c[i], c->c[i] );
    }
}
