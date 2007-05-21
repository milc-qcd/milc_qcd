/************************ control.c ******************************/
/* MIMD version 6 */
/* Main procedure for SU3 with dynamical Wilson fermions */

/* This version combines code for the PHI algorithm (approriate for 2
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL

#include "wi_dyn_includes_qdp.h"

int main(int argc, char *argv[]){
int meascount,todo;
int prompt;
double ssplaq,stplaq;
int m_iters,s_iters,avm_iters,avs_iters,avspect_iters;
#ifdef SPECTRUM
int spect_iters;
#endif
#ifdef QUARK
int disp_iters;
#endif
complex plp;
double dtime;

  initialize_machine(&argc,&argv);
  QDP_initialize(&argc, &argv);
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

	/* perform warmup trajectories */
	dtime = -dclock();
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
	    if(steps>0)s_iters=update(); else s_iters=-99;

	    /* measure every "propinterval" trajectories */
	    if((todo%propinterval) == 0){
	    
                /* generate a pseudofermion configuration */
		boundary_flip(MINUS);
		m_iters = f_measure2();

#ifdef SPECTRUM 
#ifdef SCREEN 
		boundary_flip(PLUS);
		gaugefix(ZUP,(Real)1.5,100,(Real)GAUGE_FIX_TOL);
		boundary_flip(MINUS);
		spect_iters = s_props();
		avspect_iters += spect_iters;
		spect_iters = w_spectrum_z();
		avspect_iters += spect_iters;
#ifdef QUARK
		boundary_flip(PLUS);
		/* Lorentz gauge*/
		gaugefix(8,(Real)1.5,100,(Real)GAUGE_FIX_TOL);
		boundary_flip(MINUS);
                disp_iters = quark();
                avspect_iters += disp_iters;
#endif /* ifdef QUARK */
#else	/* spectrum in time direction */
		boundary_flip(PLUS);
		gaugefix(TUP,(Real)1.5,100,(Real)GAUGE_FIX_TOL);
		boundary_flip(MINUS);
		spect_iters = t_props();
		avspect_iters += spect_iters;
		spect_iters = w_spectrum();
		avspect_iters += spect_iters;
#endif	/* end ifndef SCREEN */
#endif	/* end ifdef SPECTRUM */
		boundary_flip(PLUS);

	        /* call plaquette measuring process */
		d_plaquette(&ssplaq,&stplaq);

		/* call the Polyakov loop measuring program */
		plp = ploop();

		avm_iters += m_iters;
		avs_iters += s_iters;
	           
	        ++meascount;
	        if(this_node==0)printf("GMES %e %e %e %e %e\n",
		    (double)plp.real,(double)plp.imag,(double)m_iters,
		    (double)ssplaq,(double)stplaq);
		/* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */

		fflush(stdout);
	    }
	}	/* end loop over trajectories */

	if(this_node==0)printf("RUNNING COMPLETED\n");
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

	/* save lattice if requested */
        if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
        }
    }
    QDP_finalize();
    return 0;
}
