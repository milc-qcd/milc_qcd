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
#ifndef COULOMB
  int todo,sm_lev;
#else
  su3_matrix *links, *ape_links;
#endif
  int prompt;
  double dtime;
  
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
    
#ifdef COULOMB
    if(this_node == 0) 
      printf("Fixing to Coulomb gauge\n");
    fixflag = COULOMB_GAUGE_FIX;
    
    gaugefix(TUP,(Real)1.8,500,GAUGE_FIX_TOL);
    
    links = create_G_from_site();
    ape_links = create_G();
    ape_smear_field_dir( links, TUP, ape_links, staple_weight, u0, 0, ape_iter, 0.0 ); 
    free(links);
    hvy_pot( ape_links, max_t, max_x );
    free(ape_links);
#else
    
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
      }
      if(this_node==0)printf("SMEARING COMPLETED\n"); 
      tot_smear += smear_num[sm_lev];
      
      /* Compute simple, i.e one-plaquette, glueball operators */
      gball_simp(tot_smear);
      
      /* Compute on-axis time-like Wilson loops */
      /** w_loop1(tot_smear); **/
      
      /* Compute on-axis time-like hybrid potential loops 
	 and on-axis time-like Wilson loops */
      hybrid_loop1(tot_smear);
      
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
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      save_lattice( saveflag, savefile, stringLFN );
    }
  }
  return 0;
}
