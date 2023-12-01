/************************* integrate3d.c **********************/
/* Integrators and associated utilities for the Wilson Flow */

#include "wilson_flow_includes.h"

#ifdef REGIONS

static int region = FULLVOL;
static int this_exp_order = 0;
static Real this_stepsize = 0.;
static Real this_stoptime = 0.;

static su3_matrix **stap = NULL, **link0 = NULL;
static anti_hermitmat ***ahK = NULL, **accu = NULL;

#if ( GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER \
  || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 \
  || GF_INTEGRATOR==INTEGRATOR_ADAPT_BS ) 
static Real this_local_tol = 0.;
#endif

static
void print_adapt( char *TAG, double flowtime, double stepsize, double dist, double ratio ) {

  #if (MILC_PRECISION==1)
      node0_printf("%s %.6g %.6g %.6g %.6g\n", TAG, flowtime, stepsize, dist, stepsize );
  #else
      node0_printf("%s %.16g %.16g %.16g %.16g\n", TAG, flowtime, stepsize, dist, stepsize );
  #endif
}

/* Highest level function that integrates the flow, called from control.c */
void
run_gradient_flow( int region_flag ) {

  /* RK integration variables */
  int i;
  double flowtime;

  /* Wilson flow output variables */
  double Et_WS[2]={0.0,0.0}, Es_WS[2]={0.0,0.0};
  double Et_C[2] ={0.0,0.0}, Es_C[2] ={0.0,0.0}, charge[4]={0.0,0.0,0.0,0.0};
  double old_value=0, new_value=0;
  double der_value=0;
  char RTAG[3][12];


#ifdef SPHALERON
  // su3_matrix ***link_last_flow = NULL;
  su3_matrix **link_last_flow = NULL;
  #ifdef HALF_LATTICE_TEST
    double Et_WS_lwr[2]={0.0,0.0}, Es_WS_lwr[2]={0.0,0.0};
    double Et_C_lwr[2] ={0.0,0.0}, Es_C_lwr[2] ={0.0,0.0}, charge_lwr[4]={0.0,0.0,0.0,0.0};
    double Et_WS_upr[2]={0.0,0.0}, Es_WS_upr[2]={0.0,0.0};
    double Et_C_upr[2] ={0.0,0.0}, Es_C_upr[2] ={0.0,0.0}, charge_upr[4]={0.0,0.0,0.0,0.0};
  #endif
#endif

  region = region_flag;
  this_stepsize = stepsize;
  this_stoptime = stoptime;
  #ifdef REGIONS
    switch ( region ) {
      case FULLVOL:
        this_exp_order = exp_order;
        this_stepsize = stepsize;
        this_stoptime = stoptime;
      #if ( GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER \
        || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 \
        || GF_INTEGRATOR==INTEGRATOR_ADAPT_BS ) 
        this_local_tol = local_tol;
      #endif
        sprintf(RTAG[0],"GFLOW:");
        sprintf(RTAG[1],"ADAPT:");
        sprintf(RTAG[2]," ");
        break;
      case BULK:
        this_exp_order = exp_order_bulk;
        this_stepsize = stepsize_bulk;
        this_stoptime = stoptime_bulk;
      #if ( GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER \
        || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 \
        || GF_INTEGRATOR==INTEGRATOR_ADAPT_BS ) 
        this_local_tol = local_tol_bulk;
      #endif
        sprintf(RTAG[0],"GFLOW_BULK:");
        sprintf(RTAG[1],"ADAPT_BULK:");
        sprintf(RTAG[2],"      ");
        break;
    case BOUNDARY:
        this_exp_order = exp_order_bdry;
        this_stepsize = stepsize_bdry;
        this_stoptime = stoptime_bdry;
      #if ( GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER \
        || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 \
        || GF_INTEGRATOR==INTEGRATOR_ADAPT_BS ) 
        this_local_tol = local_tol_bdry;
      #endif
        sprintf(RTAG[0],"GFLOW_BDRY:");
        sprintf(RTAG[1],"ADAPT_BDRY:");
        sprintf(RTAG[2],"      ");
        break;
      default:
        node0_printf("Invalid region specification: %d\n",region);
        fflush(stdout); terminate(1);
    }
  #endif
  if ( this_stepsize == 0. || this_stepsize * this_stoptime < 0. || this_exp_order == 0 ) {
    node0_printf("Invalid stepsize/stoptime/exp_order combination: %g %g %d\n",
      this_stepsize,this_stoptime,this_exp_order);
    fflush(stdout); terminate(1);
  }

  /* Print flow output column labels */
  #if (REPORT == OLDREPORT)
    node0_printf("#LABEL%stime Clover_t Clover_s Plaq_t Plaq_s Rect_t Rect_s charge\n",RTAG[2]);
  #else
    node0_printf("#LABEL%stime Clover_t Clover_s iClover_t iClover_s Plaq_t Plaq_s Rect_t Rect_s charge icharge\n",RTAG[2]);
  #endif
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  node0_printf("#LABEL2 time stepsize distance local_tol/distance\n");
#endif
  fflush(stdout);

  /* Calculate and print initial flow output */
  if ( region == FULLVOL ) {
    fmunu_fmunu_full( Et_C, Es_C, charge );
    gauge_action_w_s_full( Et_WS, Es_WS );
  }
  #ifdef SPHALERON
    if ( region == BULK ) {
      fmunu_fmunu_bulk( Et_C, Es_C, charge );
      gauge_action_w_s_bulk( Et_WS, Es_WS );
    }
    if ( region == BOUNDARY ) {
      // link_last_flow = new_last_flow_links();
      // update_last_flow_links( link_last_flow );

      link_last_flow = new_field( N_LAST_FLOW );
      update_last_flow_links( link_last_flow );
      
      fmunu_fmunu_bdry( link_last_flow, Et_C, Es_C, charge );
      gauge_action_w_s_bdry( link_last_flow, Et_WS, Es_WS );
      
      #ifdef HALF_LATTICE_TEST
        fmunu_fmunu_lwr_bdry( link_last_flow, Et_C_lwr, Es_C_lwr, charge_lwr );
        gauge_action_w_s_lwr_bdry( link_last_flow, &Et_WS_lwr, &Es_WS_lwr );
        fmunu_fmunu_upr_bdry( link_last_flow, Et_C_upr, Es_C_upr, charge_upr );
        gauge_action_w_s_upr_bdry( link_last_flow, &Et_WS_upr, &Es_WS_upr );
        print_observables( "_LWRB", 0.0, Et_WS_lwr, Es_WS_lwr, Et_C_lwr, Es_C_lwr, charge_lwr );
        print_observables( "_UPRB", 0.0, Et_WS_upr, Es_WS_upr, Et_C_upr, Es_C_upr, charge_upr );
      #endif
    }
  #endif
  print_observables( RTAG[0], 0.0, Et_WS, Es_WS, Et_C, Es_C, charge );

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  print_adapt( RTAG[1], 0.0, this_stepsize, 0.0, 0.0 );
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

  #ifdef USE_FIELD
    stap = new_field( 4 );
    accu = new_anti_hermitian_field( 4 );
    #if GF_INTEGRATOR==INTEGRATOR_RKMK3
      ahK = new_anti_hermitian_twodim_field( 3, 4 ); /* right-hand-side in RK method for all stages */
    #elif GF_INTEGRATOR==INTEGRATOR_RKMK4
      ahK = new_anti_hermitian_twodim_field( 4, 4 ); /* right-hand-side in RK method for all stages */
    #elif GF_INTEGRATOR==INTEGRATOR_RKMK5
      ahK = new_anti_hermitian_twodim_field( 6, 4 ); /* right-hand-side in RK method for all stages */
    #elif GF_INTEGRATOR==INTEGRATOR_RKMK8
      ahK = new_anti_hermitian_twodim_field( 13, 4 ); /* right-hand-side in RK method for all stages */
    #elif GF_INTEGRATOR==INTEGRATOR_ADAPT_EULER || \
          GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
      ahK = new_anti_hermitian_twodim_field( 3, 4 ); /* right-hand-side in RK method for all stages */
    #elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
      ahK = new_anti_hermitian_twodim_field( 4, 4 ); /* right-hand-side in RK method for all stages */
    #endif
  #endif

  /* Loop over the flow time */
  while( this_stoptime==AUTO_STOPTIME || ( flowtime<this_stoptime && is_final_step==0 ) ) {
    /* Adjust last time step to fit exactly stoptime_bulk */
    if( this_stepsize>this_stoptime-flowtime && this_stoptime!=AUTO_STOPTIME ) {
      this_stepsize = this_stoptime-flowtime;
      is_final_step = 1;
    }
    // printf("%g\n", this_stepsize);

    #ifdef SPHALERON
      if ( i > 0 && region == BOUNDARY ) 
        update_last_flow_links( link_last_flow );
    #endif

    /* Perform one flow step (most of the computation is here) */
    flow_step();
    flowtime += this_stepsize;
    i++;

    /* Calculate and print current flow output */
    if ( region == FULLVOL ) {
      fmunu_fmunu_full( Et_C, Es_C, charge );
      gauge_action_w_s_full( Et_WS, Es_WS );
    }
    #ifdef SPHALERON
      if ( region == BULK ) {
        fmunu_fmunu_bulk( Et_C, Es_C, charge );
        gauge_action_w_s_bulk( Et_WS, Es_WS );
      }
      if ( region == BOUNDARY ) {
        fmunu_fmunu_bdry( link_last_flow, Et_C, Es_C, charge );
        gauge_action_w_s_bdry( link_last_flow, Et_WS, Es_WS );
        #ifdef HALF_LATTICE_TEST
          fmunu_fmunu_lwr_bdry( link_last_flow, Et_C_lwr, Es_C_lwr, charge_lwr );
          gauge_action_w_s_lwr_bdry( link_last_flow, Et_WS_lwr, Es_WS_lwr );
          fmunu_fmunu_upr_bdry( link_last_flow, Et_C_upr, Es_C_upr, charge_upr );
          gauge_action_w_s_upr_bdry( link_last_flow, Et_WS_upr, Es_WS_upr );
        #endif
      }
    #endif
    #if ( REPORT != NO_REPORT )
      print_observables( RTAG[0], flowtime, Et_WS, Es_WS, Et_C, Es_C, charge );
    #endif

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    print_adapt( RTAG[1], flowtime, this_stepsize, dist, this_local_tol/dist );
#endif
    fflush(stdout);

    /* Automatic determination of stoptime_bulk:                         */
    /*  t^2 E > 0.45 and d/dt { t^2 E } > 0.35                      */
    /*  Bounds need to be adjusted with scale determination cutoff  */
    if( this_stoptime==AUTO_STOPTIME ) {

      old_value = new_value;
      new_value = flowtime*flowtime*(Et_C[0]+Es_C[0]);
      der_value = flowtime*(new_value-old_value)/this_stepsize;

      if( new_value > 0.45 && der_value > 0.35 )
        break;
    } /* end: auto stoptime_bulk */

#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    if( is_final_step==0 ) {
      // adjust step size for the next step except if it is final
      this_stepsize = this_stepsize * SAFETY * pow( this_local_tol/dist, 1/3. );
    }
#endif

  } /* end: flowtime loop */

  #if REPORT == NO_REPORT
    print_observables( RTAG[0], flowtime, Et_WS, Es_WS, Et_C, Es_C, charge );
  #endif
  #ifdef SPHALERON
    if ( region == BOUNDARY ) {
      print_observables( "ACCUM_BDRY", flowtime, Et_WS, Es_WS, Et_C, Es_C, charge );
      q_acc[0] += charge[2];
      q_acc[1] += charge[3];
      // destroy_last_flow_links( &link_last_flow );
      destroy_field( &link_last_flow );
    }
  #endif
  /* Save and print the number of steps */
  total_steps = i;
  node0_printf("Number of steps = %i\n", total_steps);
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
  GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  node0_printf("Number of rejected steps = %i\n", steps_rejected);
#endif
  fflush(stdout);
  #ifdef USE_FIELD
    if (ahK != NULL )
      destroy_anti_hermitian_twodim_field( &ahK );
    destroy_anti_hermitian_field( &accu );
    destroy_field( &stap );
  #endif
}

#if GF_INTEGRATOR==INTEGRATOR_EULER || \
    GF_INTEGRATOR==INTEGRATOR_LUSCHER || GF_INTEGRATOR==INTEGRATOR_CK \
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
integrate_RK_2N_one_stage( Real cA, Real cB )
{
  register int dir, i;
  register site *s;

  /* Temporary matrix holders */
  anti_hermitmat *Acur, tempA1;
  su3_matrix *U, tempS1, tempS2;


  /* Calculate the new staple */
  #ifdef USE_FIELD
    su3_matrix *STAP;
    staple( region, stap );
  #else
    staple( region );
  #endif

  // dumpmat( &(lattice->staple[XUP]) );
  // terminate(1);
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
    // FORALLDIRSUP(dir, region) 
    FORALLUPDIR(dir)
    IF_LINK_IN_REGION(s, dir, region) 
    {
      // IF_BOUNDARY(s) { node0_printf("dir %d  i %d coord %d %d %d %d\n", dir, i, s->x, s->y, s->z, s->t); }
      /* Retrieve the current link and accumulation matrix */
      U = &(s->link[dir]);
      #ifdef USE_FIELD
        Acur = &(accu[dir][i]);
      #else
        Acur = &(s->accumulate[dir]);
      #endif

      /* Update the accumulation matrix A = cA*A + proj(U*S) */
      #ifdef USE_FIELD
        STAP = &(stap[dir][i]);
        mult_su3_na( U, STAP, &tempS1 );
      #else
        mult_su3_na( U, &(s->staple[dir]), &tempS1 );
      #endif
      anti_hermitian_traceless_proj( &tempS1, &tempA1 );
      scalar_mult_add_ah( &tempA1, Acur, cA, Acur );

      /* Update the links U = exp(cB*A)*U */
      scalar_mult_ah( Acur, cB, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
    }
    // if (region == BOUNDARY) { terminate(1); }
}

/* one step of a low-storage Runge-Kutta scheme,
   this includes Luscher, arXiv:1006.4518 [hep-lat]
   or any other 2N-storage scheme */
void
integrate_RK_2N()
{
  register int dir, i;
  register site *s;

  /* Clear the accumulation matrix */
  #ifdef USE_FIELD
    clear_anti_hermitian_field( accu, 4);
  #else
    FORALLSITES(i, s)
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
      // FORALLDIRSUP(dir, region)
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
        clear_anti_hermitian(&(s->accumulate[dir]));
  #endif

  /* Infinitesimal stout smearing */
  for( i=0; i<N_stages; i++ ) {
    /* be careful with stepsize: Nathan's convention on the staple
       is such that stepsize should be taken negative */
    integrate_RK_2N_one_stage( A_2N[i], -B_2N[i]*this_stepsize );
  }
}
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4 || GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8
/* generic integrator of Runge-Kutta-Munthe-Kaas type,
   requires dexpinv evaluation at each step (nested commutators) */
void
integrate_RKMK_generic() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, *Atemp;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  #ifdef USE_FIELD
    su3_matrix *LINK0, *STAP;
    link0 = new_links_from_site( region, 4 );
  #else
    FORALLSITES(i, s)
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
      // FORALLDIRSUP(dir, region)
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
        su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );
  #endif

  // loop over RK stages
  for( i_rk=0; i_rk<N_stages; i_rk++ ) {
    if( i_rk!=0 ) {
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region) 
        // FORALLDIRSUP(dir, region)
        FORALLUPDIR(dir)
        IF_LINK_IN_REGION(s, dir, region) 
        {
          #ifdef USE_FIELD
            Atemp = &(accu[dir][i]);
          #else
            Atemp = &(s->accumulate[dir]);
          #endif
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            #ifdef USE_FIELD
              scalar_mult_add_ah( Atemp, &(ahK[j_rk][dir][i]), a_RK[i_rk][j_rk], Atemp );
            #else
              scalar_mult_add_ah( Atemp, &(s->K[j_rk][dir]), a_RK[i_rk][j_rk], Atemp );
            #endif
          }
          // update the link
          scalar_mult_ah( Atemp, -this_stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, this_exp_order );
          #ifdef USE_FIELD
            LINK0 = &(link0[dir][i]);
            mult_su3_nn( &tempS1, LINK0, &(s->link[dir]) );
          #else
            mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
          #endif
      }
    }

    // get the right hand side of the flow equation from the staple
    // and store in s->K[i_rk]
    #ifdef USE_FIELD
      staple( region, stap );
    #else
      staple( region );
    #endif

    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
      // FORALLDIRSUP(dir, region) 
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
      {
        #ifdef USE_FIELD
          STAP = &(stap[dir][i]);
          mult_su3_na( &(s->link[dir]), STAP, &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(ahK[i_rk][dir][i]) );
          if( i_rk!=0 )
            dexpinv( &(accu[dir][i]), &(ahK[i_rk][dir][i]), p_order, &(ahK[i_rk][dir][i]) );
        #else
          mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(s->K[i_rk][dir]) );
          if( i_rk!=0 )
            dexpinv( &(s->accumulate[dir]), &(s->K[i_rk][dir]), p_order, &(s->K[i_rk][dir]) );
        #endif
    }
  }
  // final RK stage
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region)
    // FORALLDIRSUP(dir, region) 
    FORALLUPDIR(dir)
    IF_LINK_IN_REGION(s, dir, region) 
    {
      clear_anti_hermitian( &tempA1 );
      for( i_rk=0; i_rk<N_stages; i_rk++ ) {
        // accumulate b_1*K1 + b_2*K2 + ...
        #ifdef USE_FIELD
          scalar_mult_add_ah( &tempA1, &(ahK[i_rk][dir][i]), b_RK[i_rk], &tempA1 );
        #else
          scalar_mult_add_ah( &tempA1, &(s->K[i_rk][dir]), b_RK[i_rk], &tempA1 );
        #endif
      }
      // update the link
      scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
      #ifdef USE_FIELD
        LINK0 = &(link0[dir][i]);
        mult_su3_nn( &tempS1, LINK0, &(s->link[dir]) );
      #else
        mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
      #endif
  }
  #ifdef USE_FIELD
    destroy_field( &link0 );
  #endif
}
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3
/* Third order Runge-Kutta-Munthe-Kaas type,
   requires a single commutator at the last stage */
void
integrate_RKMK3() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, tempA2, *Atemp;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  #ifdef USE_FIELD
    su3_matrix *LINK0, *STAP;
    link0 = new_links_from_site( region, 4 );
  #else
    FORALLSITES(i, s)
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
      // FORALLDIRSUP(dir, region)
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
        su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );
  #endif

  // loop over RK stages
  for( i_rk=0; i_rk<N_stages; i_rk++ ) {
    if( i_rk!=0 ) {  
      // FORALLDIRSUP(dir, region) 
      FORALLUPDIR(dir)
        FORALLSITES(i, s) 
        IF_BLOCKED(s, block_stride)
        // IF_REGION(s, region)
        IF_LINK_IN_REGION(s, dir, region) 
        {
          #ifdef USE_FIELD
            Atemp = &(accu[dir][i]);
          #else
            Atemp = &(s->accumulate[dir]);
          #endif
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            #ifdef USE_FIELD
              scalar_mult_add_ah( Atemp, &(ahK[j_rk][dir][i]), a_RK[i_rk][j_rk], Atemp );
            #else
              scalar_mult_add_ah( Atemp, &(s->K[j_rk][dir]), a_RK[i_rk][j_rk], Atemp );
            #endif
          }
          // update the link
          scalar_mult_ah( Atemp, -this_stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, this_exp_order );
          #ifdef USE_FIELD
            LINK0 = &(link0[dir][i]);
            mult_su3_nn( &tempS1, LINK0, &(s->link[dir]) );
          #else
            mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
          #endif
      }
    }

    // get the right hand side of the flow equation from the staple
    // and store in s->K[i_rk]
    #ifdef USE_FIELD
      staple( region, stap );
    #else
      staple( region );
    #endif

    // FORALLDIRSUP(dir, region)
    FORALLUPDIR(dir) 
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region)
      IF_LINK_IN_REGION(s, dir, region) 
      {
        #ifdef USE_FIELD
          STAP = &(stap[dir][i]);
          mult_su3_na( &(s->link[dir]), STAP, &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(ahK[i_rk][dir][i]) );
        #else
          mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(s->K[i_rk][dir]) );
        #endif

    }
  }

  // final RK stage
  // FORALLDIRSUP(dir, region) 
  FORALLUPDIR(dir) 
  {

    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
    IF_LINK_IN_REGION(s, dir, region) 
    {
      clear_anti_hermitian( &tempA1 );
      for( i_rk=0; i_rk<N_stages; i_rk++ ) {
        // accumulate b_1*K1 + b_2*K2 + ...
        #ifdef USE_FIELD
          scalar_mult_add_ah( &tempA1, &(ahK[i_rk][dir][i]), b_RK[i_rk], &tempA1 );
        #else
          scalar_mult_add_ah( &tempA1, &(s->K[i_rk][dir]), b_RK[i_rk], &tempA1 );
        #endif
      }
      // the only commutator in this scheme is at the last stage
      #ifdef USE_FIELD
        commutator_ah( &tempA1, &(ahK[0][dir][i]), &tempA2 );
      #else
        commutator_ah( &tempA1, &(s->K[0][dir]), &tempA2 );
      #endif
      scalar_mult_add_ah( &tempA1, &tempA2, -this_stepsize/6, &tempA1 );
      // update the link
      scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
      #ifdef USE_FIELD
        LINK0 = &(link0[dir][i]);
        mult_su3_nn( &tempS1, LINK0, &(s->link[dir]) );
      #else
        mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
      #endif
    }
  }
  #ifdef USE_FIELD
    destroy_field( &link0 );
  #endif
}
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
void
integrate_adapt_RK_2N_one_stage( Real cA, Real cB, int istep )
{
  register int dir, i;
  register site *s;

  /* Temporary matrix holders */
  anti_hermitmat *Acur, tempA1;
  su3_matrix *U, tempS1, tempS2;

  /* Calculate the new staple */
  #ifdef USE_FIELD
    su3_matrix *STAP;
    staple( region, stap );
  #else
    staple( region );
  #endif

  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region)
    // FORALLDIRSUP(dir, region) 
    FORALLUPDIR(dir) 
    IF_LINK_IN_REGION(s, dir, region) 
    {
      /* Retrieve the current link and accumulation matrix */
      U = &(s->link[dir]);
      #ifdef USE_FIELD
        Acur = &(accu[dir][i]);
      #else
        Acur = &(s->accumulate[dir]);
      #endif

      /* Update the accumulation matrix A = cA*A + proj(U*S) */
      #ifdef USE_FIELD
        STAP = &(stap[dir][i]);
        mult_su3_na( U, STAP, &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(ahK[istep][dir][i]) );
        scalar_mult_add_ah( &(ahK[istep][dir][i]), Acur, cA, Acur );
      #else
        mult_su3_na( U, &(s->staple[dir]), &tempS1 );
        anti_hermitian_traceless_proj( &tempS1, &(s->K[istep][dir]) );
        scalar_mult_add_ah( &(s->K[istep][dir]), Acur, cA, Acur );
      #endif

      /* Update the links U = exp(cB*A)*U */
      scalar_mult_ah( Acur, cB, &tempA1 );
      exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
      mult_su3_nn( &tempS1, U, &tempS2 );
      su3mat_copy( &tempS2, U );
    }
}

/* Adaptive scheme based on Luscher's in 2N-storage format */
void
integrate_adapt_RK_2N()
{
  register int dir, i;
  register site *s;
  int is_repeat = 1;
  anti_hermitmat tempA1;
  su3_matrix tempS1, tempS2;
  Real temp;

  #ifdef USE_FIELD
    /* Clear the accumulation matrix */
    clear_anti_hermitian_field( accu, 4 );
    /* Store the initial state of the gauge field */
    su3_matrix *LINK0;
    link0 = new_links_from_site( region, 4 );
  #else
    FORALLSITES(i, s)
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
      // FORALLDIRSUP(dir, region) 
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
      {
        /* Clear the accumulation matrix */
        clear_anti_hermitian(&(s->accumulate[dir]));
        /* Store the initial state of the gauge field */
        su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );
      }
  #endif

  do {
    /* Make one RK step */
    for( i=0; i<N_stages; i++ ) {
      /* be careful with stepsize: Nathan's convention on the staple
         is such that stepsize should be taken negative */
      integrate_adapt_RK_2N_one_stage( A_2N[i], -B_2N[i]*this_stepsize, i );
    }

    dist = 0;
    // FORALLDIRSUP(dir, region) 
    FORALLUPDIR(dir)
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region)
      IF_LINK_IN_REGION(s, dir, region) 
      {
        /* Construct lower order approximation */
        #ifdef USE_FIELD
          scalar_mult_ah( &(ahK[0][dir][i]), Lambda[0], &tempA1 );
          scalar_mult_add_ah( &tempA1, &(ahK[1][dir][i]), Lambda[1], &tempA1 );
          scalar_mult_add_ah( &tempA1, &(ahK[2][dir][i]), Lambda[2], &tempA1 );
        #else
          scalar_mult_ah( &(s->K[0][dir]), Lambda[0], &tempA1 );
          scalar_mult_add_ah( &tempA1, &(s->K[1][dir]), Lambda[1], &tempA1 );
          scalar_mult_add_ah( &tempA1, &(s->K[2][dir]), Lambda[2], &tempA1 );
        #endif
        // NOTE: Nathan Brown's convention: stepsize is negative
        scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
        #ifdef USE_FIELD
          LINK0 = &(link0[dir][i]);
          mult_su3_nn( &tempS1, LINK0, &tempS2 );
        #else
          mult_su3_nn( &tempS1, &(s->link0[dir]), &tempS2 );
        #endif
        /* Calculate distance between the two approximations */
        temp = su3mat_distance( &(s->link[dir]), &tempS2 );
        /* Find the maximum over the local volume */
        if( dist<temp ) dist = temp;
      }
    /* Get the global maximum distance */
    g_floatmax( &dist );
    /* Check if tolerance is exceeded, redo the step
       except if it is final */
    if( dist>this_local_tol && is_final_step==0 ) {
      // adjust step size
      this_stepsize = this_stepsize * SAFETY * pow( this_local_tol/dist, 1/3. );
      // record failed step
      steps_rejected++;
      // copy over the original state of the gauge field
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region)
        // FORALLDIRSUP(dir, region) 
        FORALLUPDIR(dir) 
        IF_LINK_IN_REGION(s, dir, region) 
        {
          #ifdef USE_FIELD
            LINK0 = &(link0[dir][i]);
            su3mat_copy( LINK0, &(s->link[dir]) );
          #else
            su3mat_copy( &(s->link0[dir]), &(s->link[dir]) );
          #endif
        }
    }
    else {
      is_repeat = 0;
    }
  } while( is_repeat==1 );
  #ifdef USE_FIELD
    destroy_field( &link0 );
  #endif
}
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
/* Bogacki-Shampine adaptive 3(2) embedded pair */
void
integrate_adapt_bs() {

  register int dir, i, i_rk, j_rk;
  register site *s;
  su3_matrix tempS1, tempS2;
  anti_hermitmat tempA1, tempA2, *Atemp;
  Real temp;
  int is_repeat = 1;

  // store the initial state of the gauge field,
  // it is used at every stage in this integrator format
  // and if the step gets rejected
  #ifdef USE_FIELD
    su3_matrix *LINK0, *STAP;
    link0 = new_links_from_site( region, 4 );
  #else
    FORALLSITES(i, s)
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
    //   FORALLDIRSUP(dir, region)
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
        su3mat_copy( &(s->link[dir]), &(s->link0[dir]) );
  #endif

  // get the first force evaluation on the very first step of integration
  if( is_first_step==1 ) {
    #ifdef USE_FIELD
      staple( region, stap );
    #else
      staple( region );
    #endif
    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
    //   FORALLDIRSUP(dir, region) 
      FORALLUPDIR(dir) 
      IF_LINK_IN_REGION(s, dir, region) 
      {
        #ifdef USE_FIELD
          STAP = &(stap[dir][i]);
          mult_su3_na( &(s->link[dir]), STAP, &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(ahK[indK[0]][dir][i]) );
        #else
          mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[0]][dir]) );
        #endif
    }
    is_first_step = 0;
  }

  do {
    // loop over RK stages, skip 0 due to FSAL
    for( i_rk=1; i_rk<N_stages; i_rk++ ) {
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region)
      //   FORALLDIRSUP(dir, region) 
        FORALLUPDIR(dir) 
        IF_LINK_IN_REGION(s, dir, region) 
        {

          #ifdef USE_FIELD
            Atemp = &(accu[dir][i]);
          #else
            Atemp = &(s->accumulate[dir]);
          #endif
          clear_anti_hermitian( Atemp );
          for( j_rk=0; j_rk<i_rk; j_rk++ ) {
            // accumulate a_i1*K1 + a_i2*K2 + ...
            #ifdef USE_FIELD
              scalar_mult_add_ah( Atemp, &(ahK[indK[j_rk]][dir][i]), a_RK[i_rk][j_rk], Atemp );
            #else
              scalar_mult_add_ah( Atemp, &(s->K[indK[j_rk]][dir]), a_RK[i_rk][j_rk], Atemp );
            #endif
          }
          // update the link
          scalar_mult_ah( Atemp, -this_stepsize, Atemp );
          exp_anti_hermitian( Atemp, &tempS1, this_exp_order );
          #ifdef USE_FIELD
            LINK0 = &(link0[dir][i]);
            mult_su3_nn( &tempS1, LINK0, &(s->link[dir]) );
          #else
            mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
          #endif
        }
        // get the right hand side of the flow equation from the staple
        // and store in s->K[i_rk]
        // NOTE: here FSAL property is used, so the force in K[0] is
        //       already filled from the previous step
        #ifdef USE_FIELD
          staple( region, stap );
        #else
          staple( region );
        #endif

        FORALLSITES(i, s) 
        IF_BLOCKED(s, block_stride)
        // IF_REGION(s, region)
        //   FORALLDIRSUP(dir, region) 
          FORALLUPDIR(dir)
          IF_LINK_IN_REGION(s, dir, region) 
          {
            #ifdef USE_FIELD
              STAP = &(stap[dir][i]);
              mult_su3_na( &(s->link[dir]), STAP, &tempS1 );
              anti_hermitian_traceless_proj( &tempS1, &(ahK[indK[i_rk]][dir][i]) );
            #else
              mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
              anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[i_rk]][dir]) );
            #endif
          }
    }
    // final RK stage that gives fourth-order local approximation
    FORALLSITES(i, s) 
    IF_BLOCKED(s, block_stride)
    // IF_REGION(s, region)
    //   FORALLDIRSUP(dir, region) 
      FORALLUPDIR(dir)
      IF_LINK_IN_REGION(s, dir, region) 
      {
        clear_anti_hermitian( &tempA1 );
        for( i_rk=0; i_rk<N_stages; i_rk++ ) {
          // accumulate b_1*K1 + b_2*K2 + ...
          #ifdef USE_FIELD
            scalar_mult_add_ah( &tempA1, &(ahK[indK[i_rk]][dir][i]), b_RK[i_rk], &tempA1 );
          #else
            scalar_mult_add_ah( &tempA1, &(s->K[indK[i_rk]][dir]), b_RK[i_rk], &tempA1 );
          #endif
        }
        // the only commutator in this scheme is at the last stage
        #ifdef USE_FIELD
          commutator_ah( &tempA1, &(ahK[indK[0]][dir][i]), &tempA2 );
        #else
          commutator_ah( &tempA1, &(s->K[indK[0]][dir]), &tempA2 );
        #endif
        scalar_mult_add_ah( &tempA1, &tempA2, -this_stepsize/6, &tempA1 );
        // update the link
        scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
        #ifdef USE_FIELD
          LINK0 = &(link0[dir][i]);
          mult_su3_nn( &tempS1, LINK0, &(s->link[dir]) );
        #else
          mult_su3_nn( &tempS1, &(s->link0[dir]), &(s->link[dir]) );
        #endif
    }
    // additional stage
    #ifdef USE_FIELD
      staple( region, stap );
    #else
      staple( region );
    #endif
    dist = 0;
    // FORALLDIRSUP(dir, region) 
    FORALLUPDIR(dir)
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region)
      IF_LINK_IN_REGION(s, dir, region) 
      {
        #ifdef USE_FIELD
          STAP = &(stap[dir][i]);
          mult_su3_na( &(s->link[dir]), STAP, &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(ahK[indK[3]][dir][i]) );
        #else
          mult_su3_na( &(s->link[dir]), &(s->staple[dir]), &tempS1 );
          anti_hermitian_traceless_proj( &tempS1, &(s->K[indK[3]][dir]) );
        #endif  
        clear_anti_hermitian( &tempA1 );
        for( i_rk=0; i_rk<4; i_rk++ ) {
          // accumulate b'_1*K1 + b'_2*K2 + b'_3*K3 + b'_4*K4
          // NOTE: b' coefficients are stored as a_RK[3][0], a_RK[3][1], etc.
          #ifdef USE_FIELD
            scalar_mult_add_ah( &tempA1, &(ahK[indK[i_rk]][dir][i]), a_RK[3][i_rk], &tempA1 );
          #else
            scalar_mult_add_ah( &tempA1, &(s->K[indK[i_rk]][dir]), a_RK[3][i_rk], &tempA1 );
          #endif
        }
        // get the lower (third) order estimate
        scalar_mult_ah( &tempA1, -this_stepsize, &tempA1 );
        exp_anti_hermitian( &tempA1, &tempS1, this_exp_order );
        #ifdef USE_FIELD
          LINK0 = &(link0[dir][i]);
          mult_su3_nn( &tempS1, LINK0, &tempS2 );
        #else
          mult_su3_nn( &tempS1, &(s->link0[dir]), &tempS2 );
        #endif
        /* Calculate distance between the two approximations */
        temp = su3mat_distance( &(s->link[dir]), &tempS2 );
        /* Find the maximum over the local volume */
        if( dist<temp ) dist = temp;
    }
    /* Get the global maximum distance */
    g_floatmax( &dist );
    /* Check if tolerance is exceeded, redo the step
       except if it is final */
    if( dist>this_local_tol && is_final_step==0 ) {
      // adjust step size
      this_stepsize = this_stepsize * SAFETY * pow( this_local_tol/dist, 1/3. );
      // record failed step
      steps_rejected++;
      // copy over the original state of the gauge field
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)
      // IF_REGION(s, region)
      //   FORALLDIRSUP(dir, region)
        FORALLUPDIR(dir)
        IF_LINK_IN_REGION(s, dir, region) 
        {
          #ifdef USE_FIELD
            LINK0 = &(link0[dir][i]);
            su3mat_copy( LINK0, &(s->link[dir]) );
          #else
            su3mat_copy( &(s->link0[dir]), &(s->link[dir]) );
          #endif
        }
    }
    else {
      is_repeat = 0;
    }
  } while( is_repeat==1 );
  // permute indices to read the force on the next step
  i = indK[0];
  indK[0] = indK[3];
  indK[3] = i;
  #ifdef USE_FIELD
    destroy_field( &link0 );
  #endif
}
#endif

#endif