/************************* integrate3d.c **********************/
/* Integrators and associated utilities for the Wilson Flow */

#include "wilson_flow_includes.h"

#ifdef SPHALERON

static int region = FULLVOL;
static Real this_stepsize = 0.;
static Real this_stoptime = 0.;

/* Highest level function that integrates the flow, called from control.c */
void
run_gradient_flow_region( int region_flag ) {

  /* RK integration variables */
  int i;
  double flowtime;

  /* Wilson flow output variables */
  double Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge;
  double old_value=0, new_value=0;
  double der_value=0;

  su3_matrix ***link_last_flow = NULL;

  char RTAG[7] = "UNDEF", RSPACE[7] = "      ";
  region = region_flag;
  switch ( region ) {
    case FULLVOL:
      this_stepsize = stepsize;
      this_stoptime = stoptime;
      sprintf(RSPACE," ");
      sprintf(RTAG,":");
      break;
    case BULK:
      this_stepsize = stepsize_bulk;
      this_stoptime = stoptime_bulk;
      sprintf(RTAG,"_BULK:");
      break;
  case BOUNDARY:
      this_stepsize = stepsize_bdry;
      this_stoptime = stoptime_bdry;
      sprintf(RTAG,"_BDRY:");
      break;
    default:
      node0_printf("Invalid region specification: %d\n",region);
      fflush(stdout); terminate(1);
  }
  if ( this_stepsize == 0. && this_stepsize * this_stoptime < 0. ) {
    node0_printf("Invalid stepsize/stoptime combination: %g %g\n",this_stepsize,this_stoptime);
    fflush(stdout); terminate(1);
  }

  /* Print flow output column labels */
  node0_printf("#LABEL%stime Clover_t Clover_s Plaq_t Plaq_s Rect_t Rect_s charge\n",RSPACE);
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  node0_printf("#LABEL2 time stepsize distance local_tol/distance\n");
#endif
  fflush(stdout);

  /* Calculate and print initial flow output */
  if ( region == FULLVOL ) {
    fmunu_fmunu(&Et_C, &Es_C, &charge);
    gauge_action_w_s( &Et_W, &Es_W, &Et_S, &Es_S );
  }
  if ( region == BULK ) {
    fmunu_fmunu_bulk(&Et_C, &Es_C, &charge);
    gauge_action_w_s_bulk( &Et_W, &Es_W, &Et_S, &Es_S );
  }
  if ( region == BOUNDARY ) {
    link_last_flow = new_last_flow_links();
    update_last_flow_links( link_last_flow );

    fmunu_fmunu_bdry( link_last_flow, &Et_C, &Es_C, &charge);
    gauge_action_w_s_bdry( link_last_flow, &Et_W, &Es_W, &Et_S, &Es_S );
  }
#if (MILC_PRECISION==1)
  node0_printf("GFLOW%s %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", RTAG, 0.0, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#else
  node0_printf("GFLOW%s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", RTAG, 0.0, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#if (MILC_PRECISION==1)
    node0_printf("ADAPT%s %.6g %.6g %.6g %.6g\n", RTAG, 0.0, this_stepsize, 0.0, 0.0 );
#else
    node0_printf("ADAPT%s %.16g %.16g %.16g %.16g\n", RTAG, 0.0, this_stepsize, 0.0, 0.0 );
#endif
#endif
  fflush(stdout);

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  steps_rejected = 0; // count rejected steps in adaptive schemes
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  is_first_step = 1; // need to know the first step for FSAL
  // set the permutation array for FSAL, this saves copying
  // K[3] to K[0] after each step
  indK[0] = 0; indK[1] = 1; indK[2] = 2; indK[3] = 3;
#endif
  is_final_step = 0;
  flowtime = 0;
  i = 0;
  /* Loop over the flow time */
  while( this_stoptime==AUTO_STOPTIME || ( flowtime<this_stoptime && is_final_step==0 ) ) {
    /* Adjust last time step to fit exactly stoptime_bulk */
    if( this_stepsize>this_stoptime-flowtime && this_stoptime!=AUTO_STOPTIME ) {
      this_stepsize = this_stoptime-flowtime;
      is_final_step = 1;
    }
//      printf("%g\n", this_stepsize);


    if ( i > 0 && region == BOUNDARY ) 
      update_last_flow_links( link_last_flow );
    /* Perform one flow step (most of the computation is here) */
    flow_step_region();

    flowtime += this_stepsize;
    i++;

    /* Calculate and print current flow output */
    if ( region == FULLVOL ) {
      fmunu_fmunu(&Et_C, &Es_C, &charge);
      gauge_action_w_s( &Et_W, &Es_W, &Et_S, &Es_S );
    }
    if ( region == BULK ) {
      fmunu_fmunu_bulk(&Et_C, &Es_C, &charge);
      gauge_action_w_s_bulk( &Et_W, &Es_W, &Et_S, &Es_S );
    }
    if ( region == BOUNDARY ) {
      fmunu_fmunu_bdry( link_last_flow, &Et_C, &Es_C, &charge);
      gauge_action_w_s_bdry( link_last_flow, &Et_W, &Es_W, &Et_S, &Es_S );
    }
#if (MILC_PRECISION==1)
    node0_printf("GFLOW%s %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", RTAG, flowtime, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#else
    node0_printf("GFLOW%s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", RTAG, flowtime, Et_C, Es_C, Et_W, Es_W, Et_S, Es_S, charge);
#endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
#if (MILC_PRECISION==1)
    node0_printf("ADAPT%s %.6g %.6g %.6g %.6g\n", RTAG, flowtime, this_stepsize, dist, local_tol/dist );
#else
    node0_printf("ADAPT%s %.16g %.16g %.16g %.16g\n", RTAG, flowtime, this_stepsize, dist, local_tol/dist );
#endif
#endif
    fflush(stdout);

    /* Automatic determination of stoptime_bulk:                         */
    /*  t^2 E > 0.45 and d/dt { t^2 E } > 0.35                      */
    /*  Bounds need to be adjusted with scale determination cutoff  */
    if( this_stoptime==AUTO_STOPTIME ) {

      old_value = new_value;
      new_value = flowtime*flowtime*(Et_C+Es_C);
      der_value = flowtime*(new_value-old_value)/this_stepsize;

      if( new_value > 0.45 && der_value > 0.35 )
        break;
    } /* end: auto stoptime_bulk */

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    if( is_final_step==0 ) {
      // adjust step size for the next step except if it is final
      this_stepsize = this_stepsize * SAFETY * pow( local_tol/dist, 1/3. );
    }
#endif

  } /* end: flowtime loop */

  /* Save and print the number of steps */
  total_steps = i;
  node0_printf("Number of steps = %i\n", total_steps);
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  node0_printf("Number of rejected steps = %i\n", steps_rejected);
#endif
  fflush(stdout);
  if ( region == BOUNDARY ) {
    destroy_last_flow_links( link_last_flow );
  }
}

#if GF_INTEGRATOR==INTEGRATOR_LUSCHER || GF_INTEGRATOR==INTEGRATOR_CK \
 || GF_INTEGRATOR==INTEGRATOR_BBB || GF_INTEGRATOR==INTEGRATOR_CF3
/* A single step for a 2N-storage Runge-Kutta scheme
 * where the right hand side of the flow equation is evaluated
 * and the fields are updated
 *  A: accumulating matrix over all smear steps in a single time step
 *  S: staple (action dependent); recalculated before each smear step
 *  U: gauge links
 *  cA,cB: constants
 *    Calculating a single smear step is done by
 *    (this follows usual convention on 2N-storage RK schemes):
 *    A = cA*A + proj(S*U) -> update accumulation matrix
 *    U = exp(cB*A)*U   -> update gauge links
 */
void
integrate_RK_2N_one_stage_region( Real cA, Real cB )
{
  register int dir, i;
  register site *s;

  /* Temporary matrix holders */
  anti_hermitmat *Acur, tempA1;
  su3_matrix *U, tempS1, tempS2;

  /* Calculate the new staple */
  staple_region( region );

  // dumpmat( &(lattice->staple[XUP]) );
  // terminate(1);
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region) 
    FORALLDIRSUP(dir, region) {
      /* Retrieve the current link and accumulation matrix */
      U = &(s->link[dir]);
      Acur = &(s->accumulate[dir]);

      /* Update the accumulation matrix A = cA*A + proj(U*S) */
      mult_su3_na( U, &(s->staple[dir]), &tempS1 );
      anti_hermitian_traceless_proj( &tempS1, &tempA1 );
      scalar_mult_add_ah( &tempA1, Acur, cA, Acur );

      /* Update the links U = exp(cB*A)*U */
      scalar_mult_ah( Acur, cB, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
  }
}

/* one step of a low-storage Runge-Kutta scheme,
   this includes Luscher, arXiv:1006.4518 [hep-lat]
   or any other 2N-storage scheme */
void
integrate_RK_2N_region()
{
  register int dir, i;
  register site *s;

  /* Clear the accumulation matrix */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region)
      clear_anti_hermitian(&(s->accumulate[dir]));

  /* Infinitesimal stout smearing */
  for( i=0; i<N_stages; i++ ) {
    /* be careful with stepsize: Nathan's convention on the staple
       is such that stepsize should be taken negative */
    integrate_RK_2N_one_stage_region( A_2N[i], -B_2N[i]*this_stepsize );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4 || GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8
/* generic integrator of Runge-Kutta-Munthe-Kaas type,
   requires dexpinv evaluation at each step (nested commutators) */
void
integrate_RKMK_generic_region() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, *Atemp;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region)
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );

  // loop over RK stages
  for( i_rk=0; i_rk<N_stages; i_rk++ ) {
    if( i_rk!=0 ) {
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      IF_REGION(s, region) 
        FORALLDIRSUP(dir, region) {
          Atemp = &(s->accumulate[dir]);
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            scalar_mult_add_ah( Atemp, &(s->K[j_rk][dir]), a_RK[i_rk][j_rk], Atemp );
          }
          // update the link
          scalar_mult_ah( Atemp, -this_stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, exp_order );
          mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
      }
    }

    // get the right hand side of the flow equation from the staple
    // and store in s->K[i_rk]
    staple_region( region );

    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    IF_REGION(s, region)
      FORALLDIRSUP(dir, region) {
        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[i_rk][dir]) );
        if( i_rk!=0 )
          dexpinv( &(s->accumulate[dir]), &(s->K[i_rk][dir]), p_order, &(s->K[i_rk][dir]));
    }
  }
  // final RK stage
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region) {
      clear_anti_hermitian( &tempA1 );
      for( i_rk=0; i_rk<N_stages; i_rk++ ) {
        // accumulate b_1*K1 + b_2*K2 + ...
        scalar_mult_add_ah( &tempA1, &(s->K[i_rk][dir]), b_RK[i_rk], &tempA1 );
      }
      // update the link
      scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3
/* Third order Runge-Kutta-Munthe-Kaas type,
   requires a single commutator at the last stage */
void
integrate_RKMK3_region() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, tempA2, *Atemp;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region)
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );

  // loop over RK stages
  for( i_rk=0; i_rk<N_stages; i_rk++ ) {
    if( i_rk!=0 ) {  
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      IF_REGION(s, region)
        FORALLDIRSUP(dir, region) {
          Atemp = &(s->accumulate[dir]);
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            scalar_mult_add_ah( Atemp, &(s->K[j_rk][dir]), a_RK[i_rk][j_rk], Atemp );
          }
          // update the link
          scalar_mult_ah( Atemp, -this_stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, exp_order );
          mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
      }
    }

    // get the right hand side of the flow equation from the staple
    // and store in s->K[i_rk]
    staple_region( region );

    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
      IF_REGION(s, region)
      FORALLDIRSUP(dir, region) {
        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[i_rk][dir]) );
    }
  }
  // final RK stage
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region) {
      clear_anti_hermitian( &tempA1 );
      for( i_rk=0; i_rk<N_stages; i_rk++ ) {
        // accumulate b_1*K1 + b_2*K2 + ...
        scalar_mult_add_ah( &tempA1, &(s->K[i_rk][dir]), b_RK[i_rk], &tempA1 );
      }
      // the only commutator in this scheme is at the last stage
      commutator_ah( &tempA1, &(s->K[0][dir]), &tempA2 );
      scalar_mult_add_ah( &tempA1, &tempA2, -this_stepsize/6, &tempA1 );
      // update the link
      scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
void
integrate_adapt_RK_2N_one_stage_region( Real cA, Real cB, int istep )
{
  register int dir, i;
  register site *s;

  /* Temporary matrix holders */
  anti_hermitmat *Acur, tempA1;
  su3_matrix *U, tempS1, tempS2;

  /* Calculate the new staple */
  staple_region( region );

  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region) {
      /* Retrieve the current link and accumulation matrix */
      U = &(s->link[dir]);
      Acur = &(s->accumulate[dir]);

      /* Update the accumulation matrix A = cA*A + proj(U*S) */
      mult_su3_na( U, &(s->staple[dir]), &tempS1 );
      anti_hermitian_traceless_proj( &tempS1, &(s->K[istep][dir]) );
      scalar_mult_add_ah( &(s->K[istep][dir]), Acur, cA, Acur );

      /* Update the links U = exp(cB*A)*U */
      scalar_mult_ah( Acur, cB, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
  }
}

/* Adaptive scheme based on Luscher's in 2N-storage format */
void
integrate_adapt_RK_2N_region()
{
  register int dir, i;
  register site *s;
  int is_repeat = 1;
  anti_hermitmat tempA1;
  su3_matrix tempS1, tempS2;
  Real temp;

  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region) {
      /* Clear the accumulation matrix */
      clear_anti_hermitian(&(s->accumulate[dir]));
      /* Store the initial state of the gauge field */
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );
  }

  do {
    /* Make one RK step */
    for( i=0; i<N_stages; i++ ) {
      /* be careful with stepsize: Nathan's convention on the staple
         is such that stepsize should be taken negative */
      integrate_adapt_RK_2N_one_stage_region( A_2N[i], -B_2N[i]*this_stepsize, i );
    }

    dist = 0;
    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    IF_REGION(s, region)
      FORALLDIRSUP(dir, region) {
        /* Construct lower order approximation */
        scalar_mult_ah( &(s->K[0][dir]), Lambda[0], &tempA1 );
        scalar_mult_add_ah( &tempA1, &(s->K[1][dir]), Lambda[1], &tempA1 );
        scalar_mult_add_ah( &tempA1, &(s->K[2][dir]), Lambda[2], &tempA1 );
        // NOTE: Nathan Brown's convention: stepsize is negative
        scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, exp_order );
        mult_su3_nn( &tempS1, &(s->link0[dir]), &tempS2 );
        /* Calculate distance between the two approximations */
        temp = su3mat_distance( &(s->link[dir]), &tempS2 );
        /* Find the maximum over the local volume */
        if( dist<temp ) dist = temp;
    }
    /* Get the global maximum distance */
    g_floatmax( &dist );
    /* Check if tolerance is exceeded, redo the step
       except if it is final */
    if( dist>local_tol && is_final_step==0 ) {
      // adjust step size
      this_stepsize = this_stepsize * SAFETY * pow( local_tol/dist, 1/3. );
      // record failed step
      steps_rejected++;
      // copy over the original state of the gauge field
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)
      IF_REGION(s, region)
        FORALLDIRSUP(dir, region)
          su3mat_copy( &(s->link0[dir]), &(s->link[dir]) );
    }
    else {
      is_repeat = 0;
    }
  } while( is_repeat==1 );
}
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
/* Bogacki-Shampine adaptive 3(2) embedded pair */
void
integrate_adapt_bs_region() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, tempA2, *Atemp;
  Real temp;
  int is_repeat = 1;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  // and if the step gets rejected
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    FORALLDIRSUP(dir, region)
      su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );

  // get the first force evaluation on the very first step of integration
  if( is_first_step==1 ) {
    staple_region( region );
    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    IF_REGION(s, region)
      FORALLDIRSUP(dir, region) {
        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[0]][dir]) );
    }
    is_first_step = 0;
  }

  do {
    // loop over RK stages, skip 0 due to FSAL
    for( i_rk=1; i_rk<N_stages; i_rk++ ) {
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      IF_REGION(s, region)
        FORALLDIRSUP(dir, region) {

            Atemp = &(s->accumulate[dir]);
            clear_anti_hermitian( Atemp );
            for( j_rk=0; j_rk<i_rk; j_rk++ ) {
              // accumulate a_i1*K1 + a_i2*K2 + ...
              scalar_mult_add_ah( Atemp, &(s->K[indK[j_rk]][dir]), a_RK[i_rk][j_rk], Atemp );
            }
            // update the link
            scalar_mult_ah( Atemp, -this_stepsize, Atemp );
            exp_anti_hermitian( Atemp, &tempS1, exp_order );
            mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
        }
        // get the right hand side of the flow equation from the staple
        // and store in s->K[i_rk]
        // NOTE: here FSAL property is used, so the force in K[0] is
        //       already filled from the previous step
        staple_region( region );

        FORALLSITES(i, s) 
        IF_BLOCKED(s, block_stride)
        IF_REGION(s, region)
          FORALLDIRSUP(dir, region) {
            mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
            anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[i_rk]][dir]) );
        }
    }
    // final RK stage that gives fourth-order local approximation
    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    IF_REGION(s, region)
      FORALLDIRSUP(dir, region) {
        clear_anti_hermitian( &tempA1 );
        for( i_rk=0; i_rk<N_stages; i_rk++ ) {
          // accumulate b_1*K1 + b_2*K2 + ...
          scalar_mult_add_ah( &tempA1, &(s->K[indK[i_rk]][dir]), b_RK[i_rk], &tempA1 );
        }
        // the only commutator in this scheme is at the last stage
        commutator_ah( &tempA1, &(s->K[indK[0]][dir]), &tempA2 );
        scalar_mult_add_ah( &tempA1, &tempA2, -this_stepsize/6, &tempA1 );
        // update the link
        scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, exp_order );
        mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
    }
    // additional stage
    staple_region( region );
    dist = 0;
    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    IF_REGION(s, region)
      FORALLDIRSUP(dir, region) {
        mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[3]][dir]) );
        clear_anti_hermitian( &tempA1 );
        for( i_rk=0; i_rk<4; i_rk++ ) {
          // accumulate b'_1*K1 + b'_2*K2 + b'_3*K3 + b'_4*K4
          // NOTE: b' coefficients are stored as a_RK[3][0], a_RK[3][1], etc.
          scalar_mult_add_ah( &tempA1, &(s->K[indK[i_rk]][dir]), a_RK[3][i_rk], &tempA1 );
        }
        // get the lower (third) order estimate
        scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, exp_order );
        mult_su3_nn( &tempS1, &(s->link0[dir]), &tempS2 );
        /* Calculate distance between the two approximations */
        temp = su3mat_distance( &(s->link[dir]), &tempS2 );
        /* Find the maximum over the local volume */
        if( dist<temp ) dist = temp;
    }
    /* Get the global maximum distance */
    g_floatmax( &dist );
    /* Check if tolerance is exceeded, redo the step
       except if it is final */
    if( dist>local_tol && is_final_step==0 ) {
      // adjust step size
      this_stepsize = this_stepsize * SAFETY * pow( local_tol/dist, 1/3. );
      // record failed step
      steps_rejected++;
      // copy over the original state of the gauge field
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)
      IF_REGION(s, region)
        FORALLDIRSUP(dir, region)
          su3mat_copy( &(s->link0[dir]), &(s->link[dir]) );
    }
    else {
      is_repeat = 0;
    }
  } while( is_repeat==1 );
  // permute indices to read the force on the next step
  i = indK[0];
  indK[0] = indK[3];
  indK[3] = i;
}
#endif

#endif