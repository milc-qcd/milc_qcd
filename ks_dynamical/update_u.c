/*************** update_u.c ************************************/
/* MIMD version 6 */

/* update the link matrices					*
*  								*
*  Go to sixth order in the exponential of the momentum		*
*  matrices, since unitarity seems important.			*
*  Evaluation is done as:					*
*	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
*	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
*								*
*/

#include "ks_dyn_includes.h"

void update_u( Real eps ) {

register int i,dir;
register site *s;
su3_matrix *link,temp1,temp2,htemp;
register Real t2,t3,t4,t5,t6;
/**TEMP**
Real gf_x,gf_av,gf_max;
int gf_i,gf_j;
**END TEMP **/

/**double dtime,dtime2,dclock();
dtime = -dclock();**/

/* Temporary by-hand optimization until pgcc compiler bug is fixed */
t2 = eps/2.0;
t3 = eps/3.0;
t4 = eps/4.0;
t5 = eps/5.0;
t6 = eps/6.0;

/** TEMP **
gf_av=gf_max=0.0;
**END TEMP**/
    FORALLSITES(i,s){
	for(dir=XUP; dir <=TUP; dir++){
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
/** FIND AVERAGE AND MAXIMUM MOMENTUM **/
/**TEMP **
for(gf_x=0,gf_i=0;gf_i<3;gf_i++)for(gf_j=0;gf_j<3;gf_j++)
gf_x+=cabs_sq(&(htemp.e[gf_i][gf_j]));
gf_av += gf_x;
if(gf_x > gf_max) gf_max = gf_x;
** END TEMP **/
	}
    }
/**TEMP**
g_floatsum( &gf_av );
g_floatmax( &gf_max );
gf_av /= (4*volume);
if(this_node==0)printf("MOM: %e %e\n",gf_av,gf_max);
**END TEMP **/
/**dtime += dclock();
if(this_node==0)printf("LINK_UPDATE: time = %e  mflops = %e\n",
dtime, (double)(5616.0*volume/(1.0e6*dtime*numnodes())) );**/
} /* update_u */
