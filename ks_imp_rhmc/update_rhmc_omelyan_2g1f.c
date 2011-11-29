OBSOLETE!! multi_x is no longer sized correctly for more than one pseudofermion
/********** update_omelyan.c ****************************************************/
/* MIMD version 7 */

/*
 See Takaishi and de Forcrand hep-lat/-0505020
 Update lattice.
 Two gauge steps for one fermion force step ("epsilon" is time for one fermion step)
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

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
  "alpha" and "beta" are adjustable parameters.  At alpha=beta=0 this
   is just leapfrog integration.
  Omelyan "optimum alpha" is 2*(0.25-0.1932) ~ 0.1
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */

int update()  {
  int step, iters=0;
  Real final_rsq;
  double startaction,endaction,d_action();
  Real xrandom;
  int i,j; site *s;
  su3_vector *multi_x[MAX_RAT_ORDER];
  su3_vector *sumvec;
  Real alpha,beta;
  int iphi;
  
  alpha = 0.1;
  beta = 0.1;
  node0_printf("Omelyan integration, 2 gauge for one 1 fermion step, steps= %d eps= %e alpha= %e beta= %e\n",
	       steps,epsilon,alpha,beta);
  if (steps %2 != 0 ){
    node0_printf("BONEHEAD! need even number of steps\n");
    exit(0);
  }
  
  /* allocate space for multimass solution vectors */
  for(i=0;i<MAX_RAT_ORDER;i++) multi_x[i]=(su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
  sumvec = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
  
  /* refresh the momenta */
  ranmom();
  
  /* generate a pseudofermion configuration only at start*/
  for(iphi = 0; iphi < nphi; iphi++){
    grsource_imp_rhmc( F_OFFSET(phi[iphi]), &(rparam[iphi].GR), EVEN,
		       multi_x,sumvec, rsqmin_gr[iphi], niter_gr[iphi]);
  }

  /* find action */
  startaction=d_action_rhmc(multi_x,sumvec);
  /* copy link field to old_link */
  gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
  
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
    invalidate_fermion_links(fn_links);
    //  invalidate_all_ferm_links(&fn_links);
    //  invalidate_all_ferm_links(&fn_links_dmdu0);
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

