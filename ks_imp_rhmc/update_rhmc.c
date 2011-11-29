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
#ifdef MILC_GLOBAL_DEBUG
#include "../include/su3_mat_op.h"
#endif

#define mat_invert mat_invert_uml
/**#define mat_invert mat_invert_cg**/

int update()  {
  int step, iters=0;
  double startaction,endaction;
#ifdef MILC_GLOBAL_DEBUG
  double tempaction;
#endif
#ifdef HMC
  Real xrandom;
#endif
  int i,j;
  su3_vector **multi_x;
  int n_multi_x;	// number of vectors in multi_x
  su3_vector *sumvec;
  int iphi, int_alg, inaik, jphi, n;
  Real lambda, alpha, beta; // parameters in integration algorithms
  imp_ferm_links_t** fn;

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
  iphi=0;
#if FERM_ACTION == HISQ
  n = fermion_links_get_n_naiks(fn_links);
#else
  n = 1;
#endif
  for( inaik=0; inaik<n; inaik++ ) {
    for( jphi=0; jphi<n_pseudo_naik[inaik]; jphi++ ) {
      restore_fermion_links_from_site(fn_links, prec_gr[iphi]);
      fn = get_fm_links(fn_links);
      grsource_imp_rhmc( F_OFFSET(phi[iphi]), &(rparam[iphi].GR), EVEN,
			 multi_x, sumvec, rsqmin_gr[iphi], niter_gr[iphi],
			 prec_gr[iphi], fn[inaik], inaik, 
			 rparam[iphi].naik_term_epsilon);
      iphi++;
    }
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
        iters += update_h_rhmc( epsilon, multi_x);
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
        iters += update_h_rhmc( epsilon, multi_x);
        update_u(epsilon*(2.0-lambda));
        iters += update_h_rhmc( epsilon, multi_x);
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
	    iters += update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (0.75+0.5*alpha)-(0.5-beta) ) );
	    update_h_gauge( 0.5*epsilon);

     	    update_u( epsilon*( (1.25-0.5*alpha)-(0.75+0.5*alpha) ) );
	    update_h_gauge( 0.5*epsilon);
     	    update_u( epsilon*( (1.5+beta)-(1.25-0.5*alpha) ) );
	    iters += update_h_fermion( epsilon, multi_x);
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
#ifdef MILC_GLOBAL_DEBUG
            global_current_time_step = step-1;
            node0_printf( "Current time step: %d\n", global_current_time_step );
#endif /* MILC_GLOBAL_DEBUG */
	    /* update U's and H's - see header comment */
     	    update_u( epsilon*( (1.0/6.0-alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (0.5-beta)-(1.0/6.0-alpha/3.0) ) );
	    iters += update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (3.0/6.0+alpha/3.0)-(0.5-beta) ) );
	    update_h_gauge( epsilon/3.0);

     	    update_u( epsilon*( (5.0/6.0-alpha/3.0)-(3.0/6.0+alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (7.0/6.0+alpha/3.0)-(5.0/6.0-alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
            {
            double max_delta_phase = 0.0;
            double max_delta_norm = 0.0;
            site *s;
            int idir;
            double delta_phase,delta_norm;
            su3_matrix Wdiff;
            double min_det_V, max_det_V;
            double min_eigen,max_eigen,min_denom;
            double Wm1unit; /* norm of (W^+W - I) */
            double flag_detV=1;
            FORALLSITES(i,s) {
              for( idir=XUP;idir<=TUP;idir++ ) {
                /* phase deviation */
                delta_phase = fabs( lattice[i].phase_Y[idir] -
                                    lattice[i].phase_Y_previous[idir] );
                if( delta_phase > max_delta_phase ) max_delta_phase = delta_phase;
                /* Wlink norm deviation */
                sub_su3_matrix( &(lattice[i].Wlink[idir]),
                                &(lattice[i].Wlink_previous[idir]), &Wdiff );
                delta_norm = su3_norm_frob( &Wdiff );
                if( delta_norm > max_delta_norm ) max_delta_norm = delta_norm;

                if( 1==flag_detV ) {
                  min_det_V = lattice[i].Vdet[idir];
                  max_det_V = min_det_V;
                  min_eigen = lattice[i].gmin[idir];
                  max_eigen = lattice[i].gmax[idir];
                  min_denom = lattice[i].denom[idir];
                  Wm1unit = lattice[i].unitW1[idir];
                  flag_detV = 0;
                }
                else {
                  if( min_det_V > lattice[i].Vdet[idir] )
                    min_det_V = lattice[i].Vdet[idir];
                  if( max_det_V < lattice[i].Vdet[idir] )
                    max_det_V = lattice[i].Vdet[idir];
                  if( min_eigen > lattice[i].gmin[idir] )
                    min_eigen = lattice[i].gmin[idir];
                  if( max_eigen < lattice[i].gmax[idir] )
                    max_eigen = lattice[i].gmax[idir];
                  if( min_denom > lattice[i].denom[idir] )
                    min_denom = lattice[i].denom[idir];
                  if( Wm1unit < lattice[i].unitW1[idir] )
                    Wm1unit = lattice[i].unitW1[idir];
                }
              }
            }
            g_doublemax( &max_delta_phase );
            g_doublemax( &max_delta_norm );
            min_det_V = -min_det_V;
            g_doublemax( &min_det_V );
            g_doublemax( &max_det_V );
            min_det_V = -min_det_V;
            min_eigen = -min_eigen;
            g_doublemax( &min_eigen );
            g_doublemax( &max_eigen );
            min_eigen = -min_eigen;
            min_denom = -min_denom;
            g_doublemax( &min_denom );
            min_denom = -min_denom;
            node0_printf("PHASE_Y maximum jump: %28.14g\n", max_delta_phase );
            node0_printf("NORM_W maximum jump: %28.14g\n", max_delta_norm );
            node0_printf("DET_V minimum: %28.14g\n", min_det_V);
            node0_printf("DET_V maximum: %28.14g\n", max_det_V);
            node0_printf("(V^+V) eigenvalue minimum: %28.14g\n", min_eigen);
            node0_printf("(V^+V) eigenvalue maximum: %28.14g\n", max_eigen);
            node0_printf("denom=ws*(us*vs-ws)  minimum: %28.14g\n", min_denom);
            node0_printf("Deviation from unitary |W^+W-1|  maximum: %28.14g\n", Wm1unit);
            }
#endif /* HISQ_REUNITARIZATION_DEBUG */
            global_current_time_step = step;
            node0_printf( "Current time step: %d\n", global_current_time_step );
#endif /* MILC_GLOBAL_DEBUG */

     	    update_u( epsilon*( (9.0/6.0-alpha/3.0)-(7.0/6.0+alpha/3.0) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (1.5+beta)-(9.0/6.0-alpha/3.0) ) );
	    iters += update_h_fermion( epsilon, multi_x);
     	    update_u( epsilon*( (11.0/6.0+alpha/3.0)-(1.5+beta) ) );
	    update_h_gauge( epsilon/3.0);
     	    update_u( epsilon*( (2.0)-(11.0/6.0+alpha/3.0) ) );

            /* reunitarize the gauge field */
	    rephase( OFF );
            reunitarize();
	    rephase( ON );
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
            {
            double max_delta_phase = 0.0;
            double max_delta_norm = 0.0;
            site *s;
            int idir;
            double delta_phase,delta_norm;
            su3_matrix Wdiff;
            double min_det_V, max_det_V;
            double min_eigen,max_eigen,min_denom;
            double Wm1unit; /* norm of (W^+W - I) */
            double flag_detV=1;
            FORALLSITES(i,s) {
              for( idir=XUP;idir<=TUP;idir++ ) {
                /* phase deviation */
                delta_phase = fabs( lattice[i].phase_Y[idir] -
                                    lattice[i].phase_Y_previous[idir] );
                if( delta_phase > max_delta_phase ) max_delta_phase = delta_phase;
                /* Wlink norm deviation */
                sub_su3_matrix( &(lattice[i].Wlink[idir]),
                                &(lattice[i].Wlink_previous[idir]), &Wdiff );
                delta_norm = su3_norm_frob( &Wdiff );
                if( delta_norm > max_delta_norm ) max_delta_norm = delta_norm;

                if( 1==flag_detV ) {
                  min_det_V = lattice[i].Vdet[idir];
                  max_det_V = min_det_V;
                  min_eigen = lattice[i].gmin[idir];
                  max_eigen = lattice[i].gmax[idir];
                  min_denom = lattice[i].denom[idir];
                  Wm1unit = lattice[i].unitW1[idir];
                  flag_detV = 0;
                }
                else {
                  if( min_det_V > lattice[i].Vdet[idir] )
                    min_det_V = lattice[i].Vdet[idir];
                  if( max_det_V < lattice[i].Vdet[idir] )
                    max_det_V = lattice[i].Vdet[idir];
                  if( min_eigen > lattice[i].gmin[idir] )
                    min_eigen = lattice[i].gmin[idir];
                  if( max_eigen < lattice[i].gmax[idir] )
                    max_eigen = lattice[i].gmax[idir];
                  if( min_denom > lattice[i].denom[idir] )
                    min_denom = lattice[i].denom[idir];
                  if( Wm1unit < lattice[i].unitW1[idir] )
                    Wm1unit = lattice[i].unitW1[idir];
                }
              }
            }
            g_doublemax( &max_delta_phase );
            g_doublemax( &max_delta_norm );
            min_det_V = -min_det_V;
            g_doublemax( &min_det_V );
            g_doublemax( &max_det_V );
            min_det_V = -min_det_V;
            min_eigen = -min_eigen;
            g_doublemax( &min_eigen );
            g_doublemax( &max_eigen );
            min_eigen = -min_eigen;
            min_denom = -min_denom;
            g_doublemax( &min_denom );
            min_denom = -min_denom;
            node0_printf("PHASE_Y maximum jump: %28.14g\n", max_delta_phase );
            node0_printf("NORM_W maximum jump: %28.14g\n", max_delta_norm );
            node0_printf("DET_V minimum: %28.14g\n", min_det_V);
            node0_printf("DET_V maximum: %28.14g\n", max_det_V);
            node0_printf("(V^+V) eigenvalue minimum: %28.14g\n", min_eigen);
            node0_printf("(V^+V) eigenvalue maximum: %28.14g\n", max_eigen);
            node0_printf("denom=ws*(us*vs-ws)  minimum: %28.14g\n", min_denom);
            node0_printf("Deviation from unitary |W^+W-1|  maximum: %28.14g\n", Wm1unit);
            }
#endif /* HISQ_REUNITARIZATION_DEBUG */
            /*TEMP - monitor action*/ tempaction = d_action_rhmc(multi_x,sumvec);
            /*TEMP - monitor action*/ node0_printf("DELTA_SO_FAR %e\n",tempaction-startaction);
#endif /* MILC_GLOBAL_DEBUG */
        }	/* end loop over microcanonical steps */
    break;
    case INT_2EPS_3TO1:
#if FERM_ACTION == HISQ
      printf("update(%d): INT_2EPS_3TO1 is not supported for HISQ\n",
	     this_node);
      terminate(1);
#endif
        /* do "steps" microcanonical steps (one "step" = one force evaluation)"  */
        for(step=6; step <= steps; step+=6){
	    /* update U's and H's - first Omelyan step */
     	    update_u(0.5*epsilon*lambda);
	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff,
				   fn_links );
     	    update_u(epsilon*( 1.0 + 0.5*(1-lambda) )); // to time = (3/2)*epsilon
            eo_fermion_force_rhmc( 3.0*epsilon,  &rparam[0].MD,
				   multi_x, F_OFFSET(phi[0]), rsqmin_md[0], 
				   niter_md[0], prec_md[0], prec_ff,
				   fn_links );

     	    update_u(epsilon*( 0.5*(1.0-lambda) ));

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff,
				   fn_links );

    	    update_u(0.5*epsilon*lambda);

	    /* update U's and H's - second Omelyan step */
     	    update_u(0.5*epsilon*lambda);

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff,
				   fn_links );

     	    update_u(epsilon*(2.0-lambda));

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff,
				   fn_links );

    	    update_u(0.5*epsilon*lambda);

	    /* update U's and H's - third Omelyan step */
     	    update_u(0.5*epsilon*lambda);

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff,
				   fn_links );

     	    update_u(epsilon*(  0.5*(1.0-lambda) )); // to time 2*epsilon + epsilon/2

            eo_fermion_force_rhmc( 3.0*epsilon,  &rparam[0].MD,
				   multi_x, F_OFFSET(phi[0]), rsqmin_md[0], 
				   niter_md[0], prec_md[0], prec_ff,
				   fn_links );

     	    update_u(epsilon*( 1.0 + 0.5*(1.0-lambda) ));

	    rephase(OFF); imp_gauge_force(epsilon,F_OFFSET(mom)); rephase(ON);
            eo_fermion_force_rhmc( epsilon,  &rparam[1].MD,
				   multi_x, F_OFFSET(phi[1]), rsqmin_md[1], 
				   niter_md[1], prec_md[1], prec_ff,
				   fn_links );

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
    invalidate_fermion_links(fn_links);
    //    free_fn_links(&fn_links);
    //    free_fn_links(&fn_links_dmdu0);
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
