/************************ control.c ******************************/
/* MIMD version 6 */
/* Main procedure for pure gauge SU3 */

/* This version combines code for the refreshed molecular dynamics
   algorithm with the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM 
   should be defined.  (Actually, the changes to control.c are minimal
   and the real differences will appear in update.c */

/* This is the version for the Schroedinger functional simulation */
/* 1/23/97 and 1/5/98 UMH */

#define CONTROL
#include "schroed_pg_includes.h"

int main(int argc, char *argv[])  {
int meascount,todo;
int prompt;
double dssplaq,dstplaq,ds_deta,bd_plaq;
double dtime;

initialize_machine(&argc,&argv);

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

	/* perform measuring trajectories, reunitarizing and measuring  */
	meascount=0;            /* number of measurements               */
	ds_deta = 0.0;
	bd_plaq = 0.0;
	for(todo=trajecs; todo > 0; --todo ){ 

	    /* do the trajectories */
	    update();

	    /* measure every "propinterval" trajectories */
	    if((todo%propinterval) == 0){
            
		/* call plaquette measuring process */
		d_plaquette(&dssplaq,&dstplaq);

		/* call the coupling measuring process */
		if(bc_flag > 0){
		    coupling(&ds_deta,&bd_plaq);
		}
		++meascount;
		if(this_node==0)printf("GMES %e %e %e %e %e\n",
		    ds_deta,bd_plaq,99.9,dssplaq,dstplaq);
		/* dS/deta bd_plaq dummy ss_plaq st_plaq */

		fflush(stdout);
	    }
	}       /* end loop over trajectories */

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
