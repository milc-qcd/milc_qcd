/************************* control_2spec.c *******************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE */
/* Main procedure for SU3 with dynamical fermions 			*/
/* Naik plus fat link fermions, general gauge action */
/* Truncated version just for spectrum for two masses */

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

int main( int argc, char **argv ){
    int i;
    int prompt;
    Real ssplaq,stplaq,rpbpp,rpbpm;
    Real f_energy,f_pressure;	/* fermionic energy and pressure */
    Real f_action;			/* fermionic action */
    Real rsq;
    int avs_iters,avspect_iters,
        avbcorr_iters;
    double dtime, dclock();

    initialize_machine(argc,argv);
    g_sync();
    /* set up */
    prompt = setup();
    /* loop over input sets */
    while( readin(prompt) == 0){
	dtime = -dclock();

	if( warms != 0 || trajecs != 0 ){
	    node0_printf("Can't do updates\n"); terminate(0);
	}

	avspect_iters = avs_iters = avbcorr_iters = 0;
        rephase( OFF );
        gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL,
           F_OFFSET(tempmat1),F_OFFSET(tempvec[0]),
		 0,NULL,NULL,0,NULL,NULL);
        rephase( ON );
#ifdef FN
	valid_fatlinks = valid_longlinks = 0;
#endif
	avspect_iters += spectrum2( mass1, F_OFFSET(phi1),
				    F_OFFSET(xxx1) );
	/* Fix TUP Coulomb gauge - gauge links only*/
	rephase( OFF );
	gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL,
		 F_OFFSET(tempmat1),F_OFFSET(tempvec[0]),
		 0,NULL,NULL,0,NULL,NULL);
	rephase( ON );
#ifdef FN
	valid_fatlinks = valid_longlinks = 0;
#endif
        avspect_iters += nl_spectrum( mass1, F_OFFSET(phi1),F_OFFSET(xxx1),
					F_OFFSET(tempmat1),F_OFFSET(staple));
	avspect_iters += spectrum_mom( mass1, mass1, 
				       F_OFFSET(phi1), 1e-1 );
	avspect_iters += spectrum_nlpi2( mass1, mass1, 
					 F_OFFSET(phi1), 1e-1 );
	avspect_iters += spectrum2( mass2, F_OFFSET(phi1),
				    F_OFFSET(xxx1) );
        avspect_iters += nl_spectrum( mass2, F_OFFSET(phi1),F_OFFSET(xxx1),
					F_OFFSET(tempmat1),F_OFFSET(staple));
	avspect_iters += spectrum_mom( mass2, mass2, 
				       F_OFFSET(phi1), 1e-1 );
	avspect_iters += spectrum_nlpi2( mass2, mass2, 
					 F_OFFSET(phi1), 1e-1 );
	fflush(stdout);
	node0_printf("RUNNING COMPLETED\n"); fflush(stdout);
	    node0_printf("average cg iters for spectrum = %e\n",
		(double)avspect_iters);

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);

	/* save lattice if requested */
        if( saveflag != FORGET ){
          rephase( OFF );
          save_lattice( saveflag, savefile );
          rephase( ON );
        }
   }
    return 0;
}
