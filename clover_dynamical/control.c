/******** control.c ********/
/* Main procedure for SU3 with dynamical clover fermions 
   and Symanzik/Tadpole improved gauge action	*/
/* MIMD version 7 */

/* Modifications
   1996 Created by Matt Wingate
   6/6/98 Version 5
*/

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL

#include "cl_dyn_includes.h"

int main(int argc, char *argv[])  {
int meascount,todo;
int prompt;
double dssplaq,dstplaq;
int m_iters,s_iters,avm_iters,avs_iters,avspect_iters;
#ifdef SPECTRUM
int spect_iters;
#endif
complex plp;
double dtime;

 initialize_machine(&argc,&argv);
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

/* set up loop tables */
        make_loop_table2();

	/* perform warmup trajectories */
	dtime = -dclock();

	/* call plaquette measuring process */
/* Check: compute initial plaquette (T. D.) */
                d_plaquette(&dssplaq,&dstplaq);
                if(this_node==0)printf("START %e %e %e\n",
                    dssplaq,dstplaq, 
                        dssplaq+dstplaq);

	for(todo=warms; todo > 0; --todo ){
	    update();
	}
	if(this_node==0)printf("WARMUPS COMPLETED\n");

	/* perform measuring trajectories, reunitarizing and measuring 	*/
	meascount=0;		/* number of measurements 		*/
	plp = cmplx(99.9,99.9);
	avspect_iters = avm_iters = avs_iters = 0;
	for(todo=trajecs; todo > 0; --todo ){ 

	    /* do the trajectories */
	    s_iters=update();

	    /* Do "local" measurements every trajectory! */
            /* The action from the RG trans */
            gauge_action(&dssplaq);
            if(this_node==0)printf("ACTION_V  %e  %e\n",
                dssplaq,(dssplaq)/(double)(volume*6));

	    /* call plaquette measuring process */
	    d_plaquette(&dssplaq,&dstplaq);

	    /* call the Polyakov loop measuring program */
	    plp = ploop();

            /* generate a pseudofermion configuration */
	    boundary_flip(MINUS);
	    m_iters = f_measure_cl();
	    boundary_flip(PLUS);

	    ++meascount;
	    avm_iters += m_iters;
	    avs_iters += s_iters;
	           
	    if(this_node==0)printf("GMES %e %e %e %e %e\n",
		(double)plp.real,(double)plp.imag,(double)m_iters,
		dssplaq,dstplaq);
	    /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */

	    /* measure other stuff every "propinterval" trajectories */
	    if(((todo-1)%propinterval) == 0){

	      fixflag = NO_GAUGE_FIX;
#ifdef SPECTRUM 
#ifdef SCREEN 
	      gaugefix(ZUP,(Real)1.5,500,(Real)GAUGE_FIX_TOL);
	      invalidate_this_clov(gen_clov);
		fixflag = COULOMB_GAUGE_FIX;
		boundary_flip(MINUS);
		spect_iters = s_props_cl();
		avspect_iters += spect_iters;
		boundary_flip(PLUS);
#else	/* spectrum in time direction */
		gaugefix(TUP,(Real)1.5,500,(Real)GAUGE_FIX_TOL);
		invalidate_this_clov(gen_clov);
		fixflag = COULOMB_GAUGE_FIX;
/* commented 15 OCT 95, MBW....we'll use periodic b.c. for spect. */
/*		boundary_flip(MINUS); */
/* Don't need t_props if we are doing w_spectrum  - C. DeTar */
/*		spect_iters = t_props_cl();
		avspect_iters += spect_iters; */
		spect_iters = w_spectrum_cl();
		avspect_iters += spect_iters;
/*		boundary_flip(PLUS); */
#endif	/* end ifndef SCREEN */
#endif	/* end ifdef SPECTRUM */

	    }
	    fflush(stdout);

	}	/* end loop over trajectories */

	if(this_node==0)printf("RUNNING COMPLETED\n");
/* Check: compute final plaquette (T. D.) */
                d_plaquette(&dssplaq,&dstplaq);
                if(this_node==0)printf("STOP %e %e %e\n",
                    dssplaq,dstplaq,
                        dssplaq+dstplaq);

	if(meascount>0)  {
	    if(this_node==0)printf("average cg iters for step= %e\n",
		(double)avs_iters/meascount);
	    if(this_node==0)printf("average cg iters for measurement= %e\n",
		(double)avm_iters/meascount);
#ifdef SPECTRUM
	    if(this_node==0)printf("average cg iters for spectrum = %e\n",
		(double)avspect_iters/meascount);
#endif
	}

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
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
