/********** update_rhmc.c ****************************************************/
/* MIMD version 7 */

/*  D.T. first try at RHMC version 12/05
 *  D.T. 03/07 Added 2G1F and generalized to gang together multiple psfermions

 Update lattice by a molecular dynamics trajectory.
 Contains a selection of integration algorithms

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.

 Integrators:
    LEAPFROG - traditional 
    OMELYAN -
       See Takaishi and de Forcrand hep-lat/0505020
      "lambda" is adjustable parameter.  At lambda=1 this is just two
      leapfrog steps of length epsilon. Note my lambda is 4X the lambda in
      Takaishi and de Forcrand and my epsilon is 1/2 theirs.  This makes
      epsilon the same as for the leapfrog algorithm.
      Omelyan "optimum lambda" 4*0.1932

   2G1F -
     Two omelyan steps for gauge force per one Omelyan step for fermion
     ("epsilon" is time for one fermion step)
        update U to epsilon*(1/4-alpha/2)
        Update H by epsilon*1/2*gauge_force
        update U to epsilon*(1/2-beta)
        Update H by epsilon*fermion_force
        update U to epsilon*(3/4+alpha/2)
        Update H by epsilon*1/2*gauge_force

        update U to epsilon*(5/4-alpha/2)
        Update H by epsilon*1/2*gauge_force
        update U to epsilon*(3/2+beta)
        Update H by epsilon*fermion_force
        update U to epsilon*(7/4+alpha/2)
        Update H by epsilon*1/2*gauge_force
        update U to epsilon*(2)
   3G1F -
     Three omelyan steps for gauge force per one Omelyan step for fermion
     ("epsilon" is time for one fermion step)
        update U to epsilon*(1/6-alpha/3)
        Update H by epsilon*1/3*gauge_force
        update U to epsilon*(1/2-beta)
        Update H by epsilon*fermion_force
        update U to epsilon*(3/6+alpha/3)
        Update H by epsilon*1/3*gauge_force

        update U to epsilon*(5/6-alpha/3)
        Update H by epsilon*1/3*gauge_force
        update U to epsilon*(7/6+alpha/3)
        Update H by epsilon*1/3*gauge_force

        update U to epsilon*(9/6-alpha/3)
        Update H by epsilon*1/3*gauge_force
        update U to epsilon*(3/2+beta)
        Update H by epsilon*fermion_force
        update U to epsilon*(11/6+alpha/3)
        Update H by epsilon*1/3*gauge_force
        update U to epsilon*(2)
   2EPS_3TO1
        Trial version using different step sizes for the two factors in the determinant.
        eg 3*eps for light/strange ratio, eps for strange^(3/4)
        Three Omelyan steps for gauge and factor_two force, one leapfrog step for factor 
        one force (should upgrade to Omelyan for each)
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */

#define mat_invert mat_invert_uml
/**#define mat_invert mat_invert_cg**/

int update()  {
  int step, iters=0;
  double startaction,endaction;
#ifdef HMC
  Real xrandom;
#endif
  int i,j;
  su3_vector **multi_x;
  int n_multi_x;	// number of vectors in multi_x
  su3_vector *sumvec;
  int iphi, int_alg;
  Real lambda, alpha, beta; // parameters in integration algorithms

  int_alg = INT_ALG;
  switch(int_alg){
    case INT_LEAPFROG:
      node0_printf("Leapfrog integration, steps= %d eps= %e\n",steps,epsilon);
      n_multi_x = max_rat_order;
      for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
      if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
    break;
    case INT_OMELYAN:
      lambda = 0.8;
      node0_printf("Omelyan integration, steps= %d eps= %e lambda= %e\n",steps,epsilon,lambda);
      if (steps %2 != 0 ){
        node0_printf("BONEHEAD! need even number of steps\n");
        terminate(1);
      }
      n_multi_x = max_rat_order;
      for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
      if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
    break;
    case INT_2G1F:
      alpha = 0.1; beta = 0.1;
      node0_printf("Omelyan integration, 2 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
          steps,epsilon,alpha,beta);
      if (steps %2 != 0 ){
          node0_printf("BONEHEAD! need even number of steps\n"); terminate(0);
      }
      n_multi_x = max_rat_order;
      for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
      if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
    break;
    case INT_3G1F:
      alpha = 0.1; beta = 0.1;
      node0_printf("Omelyan integration, 3 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
          steps,epsilon,alpha,beta);
      if (steps %2 != 0 ){
          node0_printf("BONEHEAD! need even number of steps\n"); terminate(0);
      }
      n_multi_x = max_rat_order;
      for(j=0,i=0; i<n_pseudo; i++){j+=rparam[i].MD.order;}
      if(j>n_multi_x)n_multi_x=j; // Fermion force needs all multi_x at once in this algorithm
    break;
    case INT_2EPS_3TO1:
      lambda = 0.8;
      node0_printf("Omelyan integration, steps= %d eps= %e lambda= %e\n",steps,epsilon,lambda);
      node0_printf("3X step size leapfrog for factor 1 in determinant\n");
      if ( n_pseudo != 2 ){
          node0_printf("BONEHEAD! need exactly two pseudofermions\n"); terminate(0);
      }
      if (steps %6 != 0 ){
          node0_printf("BONEHEAD! need 6N number of steps\n"); terminate(0);
      }
      n_multi_x = max_rat_order;
    break;
    default:
      node0_printf("No integration algorithm, or unknown one\n");
      terminate(1);
    break;
  }
  
  /* allocate space for multimass solution vectors */

  multi_x = (su3_vector **)malloc(n_multi_x*sizeof(su3_vector *));
  if(multi_x == NULL){
    printf("update: No room for multi_x\n");
    terminate(1);
  }
  for(i=0;i<n_multi_x;i++){
    multi_x[i]=(su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
    if(multi_x[i] == NULL){
      printf("update: No room for multi_x\n");
      terminate(1);
    }
  }

  sumvec = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
  if( sumvec==NULL ){
    printf("update: No room for sumvec\n"); 
      terminate(1);
  }
  
  /* refresh the momenta */
  ranmom();
  
  /* generate a pseudofermion configuration only at start*/
  // NOTE used to clear xxx here.  May want to clear all solutions for reversibility
  load_ferm_links(&fn_links, &ks_act_paths);
  for(iphi = 0; iphi < n_pseudo; iphi++){
    grsource_imp_rhmc( F_OFFSET(phi[iphi]), &(rparam[iphi].GR), EVEN,
		       multi_x,sumvec, rsqmin_gr[iphi], niter_gr[iphi],
		       prec_gr[iphi], &fn_links);
  }

  /* find action */
  startaction=d_action_rhmc(multi_x,sumvec);
#ifdef HMC
  /* copy link field to old_link */
  gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
#endif
  
  switch(int_alg){
    case INT_LEAPFROG:
      /* do "steps" microcanonical steps"  */
      for(step=1; step <= steps; step++){
        /* update U's to middle of interval */
        update_u(0.5*epsilon);
        /* now update H by full time interval */
        update_h_rhmc( epsilon, multi_x);
        /* update U's by half time step to get to even time */
        update_u(epsilon*0.5);
        /* reunitarize the gauge field */
        rephase( OFF ); reunitarize(); rephase( ON );
        /*TEMP - monitor action*/if(step%4==0)d_action_rhmc(multi_x,sumvec);
      }	/* end loop over microcanonical steps */
    break;
    case INT_OMELYAN:
      /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
      for(step=2; step <= steps; step+=2){
        /* update U's and H's - see header comment */
        update_u(0.5*epsilon*lambda);
        update_h_rhmc( epsilon, multi_x);
        update_u(epsilon*(2.0-lambda));
        update_h_rhmc( epsilon, multi_x);
        update_u(0.5*epsilon*lambda);
        /* reunitarize the gauge field */
        rephase( OFF ); reunitarize(); rephase( ON );
        /*TEMP - monitor action*/ //if(step%4==0)d_action_rhmc(multi_x,sumvec);
      }	/* end loop over microcanonical steps */
    break;
    case INT_2G1F:
        /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
        for(step=2; step <= steps; step+=2){
	    /* update U's and H's - see header comment */
     	    update_u( epsilon*( (0.25-0.5*alpha) ) );
	    update_h_gauge( 0.5*epsilon);
     	    update_u( epsilon*( (0.5-beta)-(0.25-0.5*alpha) ) );
	    update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (0.75+0.5*alpha)-(0.5-beta) ) );
	    update_h_gauge( 0.5*epsilon);

     	    update_u( epsilon*( (1.25-0.5*alpha)-(0.75+0.5*alpha) ) );
	    update_h_gauge( 0.5*epsilon);
     	    update_u( epsilon*( (1.5+beta)-(1.25-0.5*alpha) ) );
	    update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (1.75+0.5*alpha)-(1.5+beta) ) );
	    update_h_gauge( 0.5*epsilon);
     	    update_u( epsilon*( (2.0)-(1.75+0.5*alpha) ) );

            /* reunitarize the gauge field */
	    rephase( OFF );
            reunitarize();
	    rephase( ON );
            /*TEMP - monitor action*/ //if(step%6==0)d_action_rhmc(multi_x,sumvec);
        }	/* end loop over microcanonical steps */
    break;
    case INT_3G1F:
        /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
        for(step=2; step <= steps; step+=2){
	    /* update U's and H's - see header comment */
     	    update_u( epsilon*( (1.0/6.0-alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (0.5-beta)-(1.0/6.0-alpha/3.0) ) );
	    update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (3.0/6.0+alpha/3.0)-(0.5-beta) ) );
	    update_h_gauge( epsilon/3.0);

     	    update_u( epsilon*( (5.0/6.0-alpha/3.0)-(3.0/6.0+alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (7.0/6.0+alpha/3.0)-(5.0/6.0-alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);

     	    update_u( epsilon*( (9.0/6.0-alpha/3.0)-(7.0/6.0+alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (1.5+beta)-(9.0/6.0-alpha/3.0) ) );
	    update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (11.0/6.0+alpha/3.0)-(1.5+beta) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (2.0)-(11.0/6.0+alpha/3.0) ) );

            /* reunitarize the gauge field */
	    rephase( OFF );
            reunitarize();
	    rephase( ON );
            /*TEMP - monitor action*/ //if(step%6==0)d_action_rhmc(multi_x,sumvec);
        }	/* end loop over microcanonical steps */
    break;
    case INT_2EPS_3TO1:
        /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
        for(step=6; step <= steps; step+=6){
	    /* update U's and H's - first Omelyan step */
     	    update_u(0.5*epsilon*lambda);
	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff );
     	    update_u(epsilon*( 1.0 + 0.5*(1-lambda) )); // to time = (3/2)*epsilon
            eo_fermion_force_rhmc( 3.0*epsilon,  &rparam[0].MD,
				   multi_x, F_OFFSET(phi[0]), rsqmin_md[0], 
				   niter_md[0], prec_md[0], prec_ff );

     	    update_u(epsilon*( 0.5*(1.0-lambda) ));

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff );

    	    update_u(0.5*epsilon*lambda);

	    /* update U's and H's - second Omelyan step */
     	    update_u(0.5*epsilon*lambda);

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff );

     	    update_u(epsilon*(2.0-lambda));

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff );

    	    update_u(0.5*epsilon*lambda);

	    /* update U's and H's - third Omelyan step */
     	    update_u(0.5*epsilon*lambda);

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff );

     	    update_u(epsilon*(  0.5*(1.0-lambda) )); // to time 2*epsilon + epsilon/2

            eo_fermion_force_rhmc( 3.0*epsilon,  &rparam[0].MD,
				   multi_x, F_OFFSET(phi[0]), rsqmin_md[0], 
				   niter_md[0], prec_md[0], prec_ff );

     	    update_u(epsilon*( 1.0 + 0.5*(1.0-lambda) ));

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff );

    	    update_u(0.5*epsilon*lambda);

            /* reunitarize the gauge field */
	    rephase( OFF ); reunitarize(); rephase( ON );
            /*TEMP - monitor action*/ //if(step%6==0)d_action_rhmc(multi_x,sumvec);

        }	/* end loop over microcanonical steps */
    break;
    default:
      node0_printf("No integration algorithm, or unknown one\n");
      terminate(1);
    break;
  }

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
    free_fn_links(&fn_links);
    free_fn_links(&fn_links_dmdu0);
#endif
    node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
  }
  else {
    node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
  }
#else  // not HMC
  node0_printf("CHECK: delta S = %e\n", (double)(endaction-startaction));
#endif // HMC
  
  /* free multimass solution vector storage */
  for(i=0;i<n_multi_x;i++)free(multi_x[i]);
  free(sumvec);
  
  if(steps > 0)return (iters/steps);
  else return(-99);
}

/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *ks_int_alg_opt_chr( void )
{
  switch(INT_ALG){
  case INT_LEAPFROG:
    return "INT_LEAPFROG";
    break;
  case INT_OMELYAN:
    return "INT_OMELYAN";
    break;
  case INT_2G1F:
    return "INT_2G1F";
    break;
  case INT_3G1F:
    return "INT_3G1F";
    break;
  case INT_2EPS_3TO1:
    return "INT_2EPS_3TO1";
    break;
  case INT_4MN4FP:
    return "INT_4MN4FP";
    break;
  case INT_4MN5FV:
    return "INT_4MN5FV";
    break;
  case INT_FOURSTEP:
    return "INT_FOURSTEP";
    break;
  case INT_PLAY:
    return "INT_PLAY";
    break;
  default:
    return "UNKNOWN";
  }
  return NULL;
}
