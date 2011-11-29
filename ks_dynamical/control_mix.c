/* Main procedure for SU3 with dynamical staggered fermions */
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL
#include "ks_dyn_includes.h"

int main(int argc, char *argv[]){
int meascount,traj_done;
int prompt;
double ssplaq,stplaq;
Real rpbpp,rpbpm;
Real f_energy,f_pressure;	/* fermionic energy and pressure */
Real f_action;			/* fermionic action */
Real rsq;
int m_iters,s_iters;
int avm_iters,avs_iters;
complex plp;

double dtime;

 initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
 g_sync();
    /* set up */
    prompt = setup();
    setup_analyze();

    /* loop over input sets */
    while( readin(prompt) == 0){
	/* perform warmup trajectories */
	dtime = -dclock();
	for(traj_done=0; traj_done < warms; traj_done++ ){
            update();
	}
	if(this_node==0)printf("WARMUPS COMPLETED\n");
if(warms==0 && trajecs==0 && startflag==CONTINUE)cool_half();

	/* perform measuring trajectories, reunitarizing and measuring 	*/
	meascount=0;		/* number of measurements 		*/
	plp = cmplx(99.9,99.9);
	init_analyze();
	for(traj_done=1; traj_done <=trajecs; traj_done++ ){ 

	    /* do the trajectories */
	    s_iters=update();

	    /* measure every "propinterval" trajectories */
	    if( (traj_done%propinterval) == (propinterval-1) ){
	    
                /* generate a pseudofermion configuration */
		grsource(EVEN);
                /* do conjugate gradient to get (Madj M)inverse * phi  */
		load_ferm_links(&fn_links);
		m_iters=ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
				   niter, rsqmin, PRECISION, EVEN, &rsq,
				   &fn_links);

	        /* call plaquette measuring process */
		plaquette(&ssplaq,&stplaq);
		ssplaq *= -1.0; stplaq *= -1.0;	/* KS phases change sign */

	        /* call Psi-bar-Psi and fermion energy/pressure measurement */
	        /* it also measures the pseudofermionic action 		*/
		f_measure(&rpbpp,&rpbpm,&f_energy,&f_pressure,&f_action);
		/* call the Polyakov loop measuring program */
		plp = ploop();
		plp.real *= -1.0; plp.imag *= -1.0;	/* KS phases! */


		/* Perform project-specific measurements */
		analyze(meascount);

		avm_iters += m_iters;
		avs_iters += s_iters;
	        ++meascount;
	        if(this_node==0)printf("GMES %e %e %e %e %e\n",
		    (double)plp.real,(double)plp.imag,(double)m_iters,
		    (double)ssplaq,(double)stplaq);
		/* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */
	        if(this_node==0)printf("FMES %e %e %e %e %e\n",(double)rpbpp,
		    (double)rpbpm,(double)f_energy,(double)f_pressure,
		    (double)f_action);

		fflush(stdout);
	    }
	}	/* end loop over trajectories */

	if(this_node==0)printf("RUNNING COMPLETED\n");
	if(meascount>0)  {
	    if(this_node==0)printf("average cg iters for step= %e\n",
		(double)avs_iters/meascount);
	    if(this_node==0)printf("average cg iters for measurement= %e\n",
		(double)avm_iters/meascount);
	}
	end_analyze(meascount);

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
}


cool_half(){
register int i,j,k,dir;
register site *s;
    if(this_node==0)printf("COOLING HALF LATTICE\n");
    rephase( OFF );
    FORALLSITES(i,s){
	if(s->z < nz/2) continue;
	for(dir=XUP;dir<=TUP;dir++)for(j=0;j<3;j++)for(k=0;k<3;k++){
	    if(j==k)s->link[dir].e[j][k] = cmplx(1.0,0.0);
	    else    s->link[dir].e[j][k] = cmplx(0.0,0.0);
	}
    }
    rephase( ON );
}
