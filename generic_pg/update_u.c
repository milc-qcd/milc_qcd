/*************** update_u.c ************************************/
/* MIMD version 6 */

/* UMH: Combined with Schroedinger functional version, Jan 2000 */

/* update the link matrices					*
*  								*
*  Go to sixth order in the exponential of the momentum		*
*  matrices, since unitarity seems important.			*
*  Evaluation is done as:					*
*	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
*	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
*								*
*/

#include "generic_pg_includes.h"

void update_u(Real eps) {

register int i,dir;
register site *s;
su3_matrix *link,temp1,temp2,htemp;
register Real t2,t3,t4,t5,t6;

/* Temporary by-hand optimization until pgcc compiler bug is fixed */
t2 = eps/2.0;
t3 = eps/3.0;
t4 = eps/4.0;
t5 = eps/5.0;
t6 = eps/6.0;
    FORALLSITES(i,s){
#ifdef SCHROED_FUN
	for(dir=XUP; dir <=TUP; dir++) if(dir==TUP || s->t>0){
#else
	for(dir=XUP; dir <=TUP; dir++){
#endif
	    uncompress_anti_hermitian( &(s->mom[dir]) , &htemp );
	    link = &(s->link[dir]);
	    mult_su3_nn(&htemp,link,&temp1);
            /**scalar_mult_add_su3_matrix(link,&temp1,eps/6.0,&temp2);**/
scalar_mult_add_su3_matrix(link,&temp1,t6,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
            /**scalar_mult_add_su3_matrix(link,&temp1,eps/5.0,&temp2);**/
scalar_mult_add_su3_matrix(link,&temp1,t5,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
            /**scalar_mult_add_su3_matrix(link,&temp1,eps/4.0,&temp2);**/
scalar_mult_add_su3_matrix(link,&temp1,t4,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    /**scalar_mult_add_su3_matrix(link,&temp1,eps/3.0,&temp2);**/
scalar_mult_add_su3_matrix(link,&temp1,t3,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    /**scalar_mult_add_su3_matrix(link,&temp1,eps/2.0,&temp2);**/
scalar_mult_add_su3_matrix(link,&temp1,t2,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,eps    ,&temp2); 
	    su3mat_copy(&temp2,link);
	}
    }
} /* update_u */
