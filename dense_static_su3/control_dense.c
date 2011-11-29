/**************************** control_dense.c ******************************/
/* Main procedure for pure gauge SU2 */
/* Almost quenched at nonzero density */
/* MIMD version 4 */
/* NEEDS UPDATING */

/* We use over-relaxed + Metropolis updatings in this code. */

#define CONTROL
#include "su3_dense_includes.h"

#ifdef LIGHT_PBP
int niter;
Real rsqmin;
#endif

int main( int argc, char **argv )  {
   int readin();
   int meascount,traj_done;
   int prompt;
   Real ssplaq,stplaq;
   Real rsq;
   double dtime;
   int setup();
   int i;
   
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
 
      for(traj_done=0; traj_done < warms; traj_done++ ){
	 update_dense();
      }
      if(this_node==0)printf("THERMALIZATION COMPLETED\n");
      
      /* perform measuring trajectories, reunitarizing and measuring  */
      meascount=0;            /* number of measurements               */
      for(traj_done=0; traj_done <trajecs; traj_done++ ){
	 
	 /* do the trajectories */
	 update_dense();

	 /* measure every "propinterval" trajectories */
	 if( (traj_done%propinterval) == (propinterval-1) ){
            
	    /* measure density and P-loop */
	    measure();

	    ++meascount;
	    plaquette(&ssplaq,&stplaq);
	    if(this_node==0)printf("GMES\t%e\t%e\n",
			     (double)ssplaq, (double)stplaq);

#ifdef LIGHT_PBP
	    /* measure pbp and other light quark observables */
	    rephase( ON );
	    mass=0.01;
	    niter=300;
	    rsqmin=0.0001*0.0001;
	    nflavors=1;
            grsource(EVEN);
	    load_ferm_links(&fn_links);
	    ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
		       niter, rsqmin, PRECISION, EVEN, &rsq, &fn_links);
	    f_measure();
	    rephase( OFF );
#endif
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
}
