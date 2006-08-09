/********** update_omelyan.c ****************************************************/
/* MIMD version 7 */

/*
 See Takaishi and de Forcrand hep-lat/-0505020
 Update lattice.
    compute PHI for both factors in fermion determinant at beginning
	update U by (epsilon/2)*lambda
	compute X for light and strange quarks for each term in rational 
		function (or "both factors in det")
	update H, by epsilon
	update U by epsilon * (2-lambda)
	compute X for light and strange quarks for each term in rational 
		function
	update H, by epsilon
	update U by (epsilon/2)*lambda

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
  "lambda" is adjustable parameter.  At lambda=1 this is just two
  leapfrog steps of length epsilon. Note my lambda is 4X the lambda in
  Takaishi and de Forcrand and my epsilon is 1/2 theirs.  This makes
  epsilon the same as for the leapfrog algorithm.
  Omelyan "optimum lambda" is 4*0.1932
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
Real lambda;

    lambda = 0.8;
    node0_printf("Omelyan integration, steps= %d eps= %e lambda= %e\n",steps,epsilon,lambda);
    if (steps %2 != 0 ){
	node0_printf("BONEHEAD! need even number of steps\n");
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
    for(step=2; step <= steps; step+=2){
 
	// alg_flag= -N skips some force terms, alg_flag= +N does them with 3X weight
    	//if(step%6==4) alg_flag= +2;
	//else	      alg_flag= -2;

	/* update U's and H's - see header comment */
     	update_u(0.5*epsilon*lambda);
	update_h_rhmc( alg_flag, epsilon, multi_x);
     	update_u(epsilon*(2.0-lambda));
	update_h_rhmc( alg_flag, epsilon, multi_x);
    	update_u(0.5*epsilon*lambda);

        /* reunitarize the gauge field */
	rephase( OFF );
        reunitarize();
	rephase( ON );
/*TEMP - monitor action*/ //if(step%6==0)d_action_rhmc(multi_x,sumvec);

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

