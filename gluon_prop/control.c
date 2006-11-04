/************************ control.c ******************************/
/* MIMD version 7 */
/* Main procedure for pure gauge SU3 */

/* This version combines code for the refreshed molecular dynamics
   algorithm with the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM 
   should be defined.  (Actually, the changes to control.c are minimal
   and the real differences will appear in update.c */

#define CONTROL
#include "gluon_prop_includes.h"
#define NULL_FP -1

int main(int argc, char *argv[])  {
    int prompt;
    int cg_iter;
    double dssplaq,dstplaq;
    double dtime;
    int key[4];
#define restrict rstrict /* C-90 T3D cludge */
    int restrict[4];
    int first_set = 1;

    initialize_machine(argc,argv);
#ifdef HAVE_QDP
    QDP_initialize(&argc, &argv);
#endif
    /* Remap standard I/O */
    if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

    g_sync();
    /* set up */
    prompt = setup();

    /* Set up Fourier transform */
    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 1;
    key[TUP] = 1;
/*    setup_restrict_fourier(key, restrict); */

    /* loop over input sets */
    while( readin(prompt) == 0){

	dtime = -dclock();
 
	/* gaugefix if requested */
	if( fixflag == COULOMB_GAUGE_FIX){
	    if( first_set == 1){
		key[TUP] = 0;
		setup_restrict_fourier(key, restrict);
		first_set = 0;
	    }
#ifdef GFIX
	    gaugefix(TUP, (Real)1.8, 500, (Real)1.e-7,
		     F_OFFSET(tempmat1), F_OFFSET(tempvec[0]),
		     0, NULL, NULL, 0, NULL, NULL);
	    gaugefixfft(TUP, (Real)(-0.10), 1750, (Real)GAUGE_FIX_TOL,
			F_OFFSET(staple), F_OFFSET(tempmat1),
			F_OFFSET(tempmat2), 0, NULL, NULL, 0, NULL, NULL);
	    if(this_node==0)printf("FIXED TO COULOMB GAUGE\n");
	    /* call plaquette measuring process */
	    d_plaquette(&dssplaq,&dstplaq);
	    if(this_node==0)printf("Gauge fixed PLAQ: %e %e\n",dssplaq,dstplaq);
	    if(this_node==0){
		printf("WARNING: minimal Coulomb gauge not implemented!\n");
		printf("WARNING: Gluon propagator NOT adopded to Coulomb gauge!\n");
	    }
	    fflush(stdout);
#endif
	}
	else if( fixflag == LANDAU_GAUGE_FIX){
	    if( first_set == 1){
		key[TUP] = 1;
		setup_restrict_fourier(key, restrict);
		first_set = 0;
	    }
#ifdef GFIX
	    gaugefix(8, (Real)1.8, 500, (Real)1.e-7,
		     F_OFFSET(tempmat1), F_OFFSET(tempvec[0]),
		     0, NULL, NULL, 0, NULL, NULL);
	    gaugefixfft(8, (Real)(-0.08), 1750, (Real)GAUGE_FIX_TOL,
			F_OFFSET(staple), F_OFFSET(tempmat1),
			F_OFFSET(tempmat2), 0, NULL, NULL, 0, NULL, NULL);
	    if(this_node==0)printf("FIXED TO LANDAU GAUGE\n");
	    /* call plaquette measuring process */
	    d_plaquette(&dssplaq,&dstplaq);
	    if(this_node==0)printf("Gauge fixed PLAQ: %e %e\n",dssplaq,dstplaq);
	    fflush(stdout);
#endif
	}
	else if( first_set == 1){
	    key[TUP] = 1;
	    setup_restrict_fourier(key, restrict);
	    first_set = 0;
	}

#ifdef GLUON_PROP
	/* Now compute the gluon propagator */
	gluon_prop();
#endif

#ifdef QUARK_PROP
	/* Now compute the quark propagator */
#ifdef FN
	valid_fn_links = 0;
	valid_fn_links_dmdu0 = 0;
#endif
#ifdef QUARK_RENORM
	cg_iter = quark_renorm();
#else
	cg_iter = quark_prop();
#endif
#endif

	if(this_node==0)printf("RUNNING COMPLETED\n");
#ifdef QUARK_PROP
	if(this_node==0)printf("CG iters for quark propagator = %d\n", cg_iter);
#endif

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
#ifdef HAVE_QDP
    QDP_finalize();
#endif
    normal_exit(0);
    return 0;
}
