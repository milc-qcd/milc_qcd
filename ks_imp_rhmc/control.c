/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This file is for lattice generation with the RHMC algorithm */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "rationals.h"
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int i;
  int meascount,traj_done;
  int prompt;
  Real ssplaq,stplaq;
  int s_iters, avs_iters, avspect_iters, avbcorr_iters;
  double dtime, dclock();
  
  initialize_machine(argc,argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();
  node0_printf("quarks1 orders: (md,gr,fa) = %d, %d, %d\n",MOLDYN_ORDER_1,GRSOURCE_ORDER_1,ACTION_ORDER_1);
  node0_printf("quarks2 orders: (md,gr,fa) = %d, %d, %d\n",MOLDYN_ORDER_2,GRSOURCE_ORDER_2,ACTION_ORDER_2);

  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    /* perform warmup trajectories */
    dtime = -dclock();
    for( traj_done=0; traj_done < warms; traj_done++ ){
      update();
    }
    node0_printf("WARMUPS COMPLETED\n"); fflush(stdout);
    
    /* perform measuring trajectories, reunitarizing and measuring 	*/
    meascount=0;		/* number of measurements 		*/
    avspect_iters = avs_iters = avbcorr_iters = 0;
    for( traj_done=0; traj_done < trajecs; traj_done++ ){ 
      
      /* do the trajectories */
      s_iters=update();
      
      /* measure every "propinterval" trajectories */
      if( (traj_done%propinterval)==(propinterval-1) ){
	
	/* call gauge_variable fermion_variable measuring routines */
	/* results are printed in output file */
	rephase(OFF);
	g_measure( );
	rephase(ON);
	f_meas_imp( F_OFFSET(phi1), F_OFFSET(xxx1), mass1 );
	f_meas_imp( F_OFFSET(phi2), F_OFFSET(xxx2), mass2 );
#ifdef FN
	if( !valid_fatlinks )load_fatlinks();
#endif
	avs_iters += s_iters;
	++meascount;
	fflush(stdout);
      }
    }	/* end loop over trajectories */
    
    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
    if(meascount>0)  {
      node0_printf("average cg iters for step= %e\n",
		   (double)avs_iters/meascount);
    }
    
    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
      printf("total_iters = %d\n",total_iters);
    }
    fflush(stdout);
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      rephase( OFF );
      save_lattice( saveflag, savefile, stringLFN );
      rephase( ON );
    }
  }
#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  normal_exit(0);
  return 0;
}
