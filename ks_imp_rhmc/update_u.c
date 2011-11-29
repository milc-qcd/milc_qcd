/*************** update_u.c ************************************/
/* MIMD version 7 */

/* update the link matrices					*
*  								*
*  Go to eight order in the exponential of the momentum		*
*  matrices, since higher order integrators use large step	*
*  Evaluation is done as:					*
*	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
*	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
*								*
*/

#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_u( Real eps ){

register int i,dir;
register site *s;
su3_matrix *link,temp1,temp2,htemp;
register Real t2,t3,t4,t5,t6,t7,t8;
/**TEMP**
Real gf_x,gf_av,gf_max;
int gf_i,gf_j;
**END TEMP **/

/**double dtime,dtime2,dclock();**/
/**dtime = -dclock();**/

/* Take divisions out of site loop (can't be done by compiler) */
t2 = eps/2.0;
t3 = eps/3.0;
t4 = eps/4.0;
t5 = eps/5.0;
t6 = eps/6.0;
t7 = eps/7.0;
t8 = eps/8.0;

/** TEMP **
gf_av=gf_max=0.0;
**END TEMP**/
#ifdef FN
 invalidate_fermion_links(fn_links);
 //  free_fn_links(&fn_links);
 //  free_fn_links(&fn_links_dmdu0);
#endif

    FORALLSITES(i,s){
	for(dir=XUP; dir <=TUP; dir++){
	    uncompress_anti_hermitian( &(s->mom[dir]) , &htemp );
	    link = &(s->link[dir]);
	    mult_su3_nn(&htemp,link,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t8,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t7,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t6,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t5,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t4,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t3,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,t2,&temp2);
	    mult_su3_nn(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix(link,&temp1,eps    ,&temp2); 
	    su3mat_copy(&temp2,link);
	}
    }
/**dtime += dclock();
node0_printf("LINK_UPDATE: time = %e  mflops = %e\n",
dtime, (double)(5616.0*volume/(1.0e6*dtime*numnodes())) );**/
} /* update_u */
