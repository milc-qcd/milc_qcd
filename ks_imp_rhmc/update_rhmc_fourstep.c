/********** update_fourstep.c ****************************************************/
/* MIMD version 7 */

/*
 See Takaishi and de Forcrand hep-lat/-0505020
 Update lattice.
	update U by (epsilon/2)(1+2a) 
	    compute PHI for light and strange quarks
	    compute X for light and strange quarks for each term in rational 
		function
	update H, by epsilon*(1+g)
	update U by epsilon*(1-a-b)
	    compute PHI and X
	update H, by epsilon*(1-g)
	update U by epsilon*(1+2b)
	    compute PHI and X
	update H, by epsilon*(1-g)
	update U by epsilon*(1-a-b)
	    compute PHI and X
	update H, by epsilon*(1-g)
	update U by epsilon*(1-a-b)
	    compute PHI and X
	update H, by epsilon*(1+g)
	update U by (epsilon/2)(1+2a) 

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
  a,b and g are adjustable parameters.  At a=b=c=0 this is just four
  leapfrog steps of length epsilon.
  Omelyan algorithm is a = b = -.1136, g=0
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */

#define mat_invert mat_invert_uml
/**#define mat_invert mat_invert_cg**/

#include "rationals.h"

int update()  {
int step, iters=0;
Real final_rsq;
double startaction,endaction,d_action();
Real xrandom;
int i,j; site *s;
su3_vector *multi_x[MAX_RAT_ORDER];
su3_vector *sumvec;
int alg_flag; 
Real a,b,g;

    a = -0.11; b= -0.09; g = -0.00;
    //a = b = -0.1; g = 0.0; //TEST

    node0_printf("Fourstep integration, steps= %d eps= %e a= .3f b= %.3f  g= %.3f\n",steps,epsilon,a,b,g);
    if (steps %4 != 0 ){
	node0_printf("BONEHEAD! need steps%4 == 0\n");
	exit(0);
    }
    alg_flag = 0; //default - fermion force with multiplier = 1

   /* allocate space for multimass solution vectors */
   for(i=0;i<MAX_RAT_ORDER;i++) multi_x[i]=(su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
    sumvec = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );

    /* refresh the momenta */
    ranmom();

    /* generate a pseudofermion configuration only at start*/
    grsource_imp_rhmc( F_OFFSET(phi1), mass1, A_GR_1, B_GR_1, GRSOURCE_ORDER_1, EVEN,
	multi_x,sumvec);
    grsource_imp_rhmc( F_OFFSET(phi2), mass2, A_GR_2, B_GR_2, GRSOURCE_ORDER_2, EVEN,
	multi_x,sumvec);

    /* find action */
    startaction=d_action_rhmc(multi_x,sumvec);
    /* copy link field to old_link */
    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));

    /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
    for(step=4; step <= steps; step+=4){
 
	// alg_flag= -N skips some force terms, alg_flag= +N does them with 3X weight

	/* update U's and H.   See header comment */
     	update_u( 0.5*epsilon*(1.0+2.0*a) );
	update_h_rhmc( alg_flag, epsilon*(1.0+g), multi_x);
     	update_u( epsilon*(1.0-a-b) );
	update_h_rhmc( alg_flag, epsilon*(1.0-g), multi_x);
    	update_u( epsilon*(1.0+2*b) );
	update_h_rhmc( alg_flag, epsilon*(1.0-g), multi_x);
     	update_u( epsilon*(1.0-a-b) );
	update_h_rhmc( alg_flag, epsilon*(1.0+g), multi_x);
     	update_u( 0.5*epsilon*(1.0+2.0*a) );

        /* reunitarize the gauge field */
	rephase( OFF );
        reunitarize();
	rephase( ON );
/*TEMP - monitor action*/if(step%6==0)d_action_rhmc(multi_x,sumvec);

    }	/* end loop over microcanonical steps */

    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    endaction=d_action_rhmc(multi_x,sumvec);
    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
#ifdef HMC
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( (double)(startaction-endaction) ) < xrandom ){
	if(steps > 0)
	    gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
#ifdef FN
	valid_longlinks=valid_fatlinks=0;
#endif
	node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
    }
    else {
	node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
    }
#else // not HMC
	node0_printf("CHECK: delta S = %e\n", (double)(endaction-startaction));
#endif // HMC

    /* free multimass solution vector storage */
    for(i=0;i<MAX_RAT_ORDER;i++)free(multi_x[i]);
    free(sumvec);

    if(steps > 0)return (iters/steps);
    else return(-99);
}

