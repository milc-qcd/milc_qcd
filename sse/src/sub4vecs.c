/*****************  sub4vecs.c  (in su3.a) ******************************
*									*
*  Subtract four su3_vectors from an su3_vector				*
* void sub_four_su3_vecs( su3_vector *a,*b1,*b2,*b3,*b4) 		*
* A  <-  A - B1 - B2 - B3 - B4						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* subtract four su3 vectors */
#ifndef FAST
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
    su3_vector *b3, su3_vector *b4 ){
register int i;
    for(i=0;i<3;i++){
	CSUB( a->c[i], b1->c[i], a->c[i] );
	CSUB( a->c[i], b2->c[i], a->c[i] );
	CSUB( a->c[i], b3->c[i], a->c[i] );
	CSUB( a->c[i], b4->c[i], a->c[i] );
    }
}
#else
void sub_four_su3_vecs( su3_vector *a, su3_vector *b1, su3_vector *b2,
    su3_vector *b3, su3_vector *b4 ){
	CSUB( a->c[0], b1->c[0], a->c[0] );
	CSUB( a->c[1], b1->c[1], a->c[1] );
	CSUB( a->c[2], b1->c[2], a->c[2] );
	CSUB( a->c[0], b2->c[0], a->c[0] );
	CSUB( a->c[1], b2->c[1], a->c[1] );
	CSUB( a->c[2], b2->c[2], a->c[2] );
	CSUB( a->c[0], b3->c[0], a->c[0] );
	CSUB( a->c[1], b3->c[1], a->c[1] );
	CSUB( a->c[2], b3->c[2], a->c[2] );
	CSUB( a->c[0], b4->c[0], a->c[0] );
	CSUB( a->c[1], b4->c[1], a->c[1] );
	CSUB( a->c[2], b4->c[2], a->c[2] );
}
#endif
