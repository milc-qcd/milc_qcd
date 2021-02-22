/*********************  control.c ***********************/
/* MIMD version 7 */
/* Main procedure for pure gauge SU3 */

/* original code by UMH */
/* 2/19/98 Version 5 port CD */

/* This version does some required number of smearing iterations,
   then transforms to axial gauge and computes simple, i.e. one-
   plaquette, glueball operators and (time-like) Wilson loops
   for the computation of the heavy quark potential. */

#define CONTROL
#include "hvy_qpot_includes.h"

int main(int argc, char *argv[])  {
  int todo,lsmeared=0;
  int sm_lev;
  su3_matrix /* *links,*/ *links;
  int prompt;
  double dtime;
  double etime;
  char myname[] = "main [control.c]";
  
  initialize_machine(&argc,&argv);
  
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  
  g_sync();
  /* set up */
  prompt = setup();
  
  /* loop over input sets */
  while( readin(prompt) == 0){

    if(prompt == 2)continue;
    dtime = -dclock();
    tot_smear = 0;
    
#ifndef WLOOP_MEAS 
#ifdef COULOMB
#ifdef GFIXONLY
    if (fixflag==NO_GAUGE_FIX) {
     node0_printf("Should actually fix a gauge when doing only gauge fixing\n");
     fflush(stderr);
     terminate(1);
    }
    if (saveflag == FORGET) {
     node0_printf("Should save the gauge-fixed lattice when doing only gauge fixing\n");
     fflush(stderr);
     terminate(1);
    }
#endif // END #ifdef GFIXONLY
    if (fixflag==COULOMB_GAUGE_FIX) {
     node0_printf("Fixing to Coulomb gauge, %d steps, %.1e tol\n",GAUGE_FIX_STEPS,GAUGE_FIX_TOL);
     fixflag = COULOMB_GAUGE_FIX;
     etime = -dclock();
     gaugefix(cor_dir,(Real)1.8,GAUGE_FIX_STEPS,GAUGE_FIX_TOL);
     etime += dclock();
     node0_printf("Gauge fixing: %e seconds\n",etime); fflush(stdout);
    } else {
     node0_printf("Skip gauge fixing by manual override -- make sure this is deliberate and sane!\n");
    }
#else
    etime=0;
#endif // END #ifdef COULOMB

#ifndef GFIXONLY
    /* Loop over the different smearing levels */
    for(sm_lev=0; sm_lev < no_smear_level; sm_lev++ ){

      etime = -dclock();
      /* Do the smearing iterations */
      for(todo=smear_num[sm_lev]; todo > 0; --todo ){
#ifdef SMEARING
        smearing();
        lsmeared=1;
#endif
      }
      if(lsmeared!=0)
#ifdef HYP_4D_SMEARING 
        node0_printf("HYP 4D smearing completed\n");
#endif
#ifdef APE_4D_SMEARING 
        node0_printf("APE 4D smearing completed\n");
#endif 
#ifdef APE_1D2_SMEARING 
        node0_printf("APE 1D2 smearing completed\n");
#endif 
      etime += dclock();
      tot_smear += smear_num[sm_lev];
#ifdef SMEARING
      node0_printf("%d Smearing steps (total %d): %e seconds\n",smear_num[sm_lev],tot_smear,etime);
#endif // END #ifndef SMEARING

#ifdef ENLARGE_MAX_XYZ_AFTER_SMEARING
      if (sm_lev > 0)
        for (int mu=XUP; mu<TUP;mu++) 
          maxc[mu] = (2*maxc[mu]<nc[mu]/2?2*maxc[mu]:nc[mu]/2);
#elif (defined ENLARGE_MAX_X_AFTER_SMEARING)
      if (sm_lev > 0)
        maxc[XUP] = (2*maxc[XUP]<nc[XUP]/2?2*maxc[XUP]:nc[XUP]/2);
#endif 

      links = create_G_from_site();

#ifndef PLANE_AVERAGED_PLC
#ifndef NEW_HVY_POT
      /* keep old function interface available, too. At least for the time being. */
      hvy_pot( links, max_t, max_x );
#else 
/* previous code distinguished at compile time between new and old hvy_pot versions, to be phased out and replaced */
/* future code to distinguish at run time between new and old hvy_pot algorithms, to be phased out and replaced */
/* hqp_switch function distinguishes and switches between both algorithms */
/* hvy_pot needs to be revised to provide interface and functionality of new_hvy_pot */
      hqp_cycle_spatial( links, hqp_alg );
//      hqp_switch( links, hqp_alg );
#endif // END #ifndef NEW_HVY_POT
#else
      plane_averaged_plc( links );
#endif // END #ifndef PLANE_AVERAGED_PLC

      free(links);
    } // END for(sm_lev=0; sm_lev < no_smear_level; sm_lev++ )
#endif // END #ifndef GFIXONLY

#else // ELSE OF #ifndef WLOOP_MEAS 
    
    /* fix to axial gauge */
    if( startflag != CONTINUE){
      ax_gauge();
      fixflag = AXIAL_GAUGE_FIX;
      tot_smear = 0;
    }
    
    /* Compute unsmeared simple, i.e one-plaquette, glueball operators */
    if( no_smear_level > 0 ){
      gball_simp(tot_smear);
    }
    
    /* Loop over the different smearing levels */
    for(sm_lev=0; sm_lev < no_smear_level; sm_lev++ ){
      
      /* Do the smearing iterations */
      for(todo=smear_num[sm_lev]; todo > 0; --todo ){
	smearing();
        lsmeared=1;
      }

      if(this_node==0 && lsmeared!=0)
#ifdef HYP_3D_SMEARING
        printf("HYP 3D SMEARING COMPLETED\n"); 
#endif
#ifdef APE_3D_SMEARING
        printf("APE 3D SMEARING COMPLETED\n"); 
#endif
#ifdef APE_2D_SMEARING
        printf("APE 2D SMEARING COMPLETED\n"); 
#endif
#ifdef APE_1D_SMEARING
        printf("APE 1D SMEARING COMPLETED\n"); 
#endif
#ifdef APE_1D2_SMEARING
        printf("APE 1D2 SMEARING COMPLETED\n"); 
#endif
      tot_smear += smear_num[sm_lev];
      
      /* Compute simple, i.e one-plaquette, glueball operators */
      gball_simp(tot_smear);
      
#ifndef HYBRIDS_MEASURE
      /* Compute on-axis time-like Wilson loops */
      w_loop1(tot_smear);
#else
      /* Compute on-axis time-like hybrid potential loops 
	 and on-axis time-like Wilson loops */
      hybrid_loop1(tot_smear);
#endif
      
      /* Compute off-axis time-like Wilson loops, if desired */
      if( off_axis_flag == 1 ){
	w_loop2(tot_smear);
      }
    }
#endif
    if(this_node==0)printf("RUNNING COMPLETED\n");
    
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
    }
    fflush(stdout);
    dtime = -dclock();
    
#ifndef SMEARING
    /* save lattice if requested */
    if( saveflag != FORGET ){
      save_lattice( saveflag, savefile, stringLFN );
    }
#endif
  }
  return 0;
}
