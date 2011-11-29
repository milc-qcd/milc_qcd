/************************ control.c ******************************/
/* MIMD version 7 */
/* Main procedure for pure gauge SU3 */

/* This version contains code for the hybrid Monte Carlo algorithm,
   for which HMC_ALGORITHM should be defined, and for the overrelaxed/
   quasi heat bath algorithm, for which ORA_ALGORITHM should be defined.
   (Actually, the changes to control.c are minimal
   and the real differences will appear in update.c */
/* This version includes measurements of smeared Wilson loops and
   glueball operators. */

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

#define CONTROL
#include "symanzik_sl32_includes.h"

int main(int argc, char *argv[]){
int meascount,todo;
int prompt;
double dssplaq,dstplaq;
complex plp;
double dtime;

#ifdef PLCOR
int key[4];
#endif

 initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

 g_sync();
    /* set up */
    prompt = setup();

/* set up polyakof loop correlator */
#ifdef PLCOR
    setup_ploop_corr(key);
#endif

    /* loop over input sets */
    while( readin(prompt) == 0){

	/* set up loop tables */
	make_loop_table();

	/* perform warmup trajectories */
	dtime = -dclock();

	for(todo=warms; todo > 0; --todo ){
	    update();
	}
	if(this_node==0)printf("WARMUPS COMPLETED\n");

	/* perform measuring trajectories, reunitarizing and measuring 	*/
	meascount=0;		/* number of measurements 		*/
	plp = cmplx(99.9,99.9);
	for(todo=trajecs; todo > 0; --todo ){ 
	    /* do the trajectories */
	    update();

	    /* measure plaquettes and Polyakov loop every trajectory */
#ifndef G_MEASURE
	    /* call plaquette measuring process */
	    d_plaquette(&dssplaq,&dstplaq);

	    /* call the Polyakov loop measuring program */
	    plp = ploop();

	    if(this_node==0)printf("GMES %e %e %e %e %e\n",
		(double)plp.real,(double)plp.imag,99.9,
		dssplaq,dstplaq);
	    /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */
#else
	    /* do "extensive local" measurements */
	    g_measure();
#endif

	    /* measure other stuff every "propinterval" trajectories */
	    if(((todo-1)%propinterval) == 0){
	    
#ifdef PLCOR
                ploop_corr();
                monte_block(steps_rg);
                ploop_corr_m();
#endif


#ifdef BLOCK
/*
		torelon();
*/
		monte_block_ape_a(steps_rg);

#endif

	        ++meascount;

		fflush(stdout);
	    }
	}	/* end loop over trajectories */

#ifdef PLCOR
                print_ploop_corr_m(trajecs);
                print_ploop_corr(trajecs);
#endif

	if(this_node==0)printf("RUNNING COMPLETED\n");

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	}
	fflush(stdout);

	/* save lattice if requested */
	if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
	}
    }

    normal_exit(0);
    return 0;
}
