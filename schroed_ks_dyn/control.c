/************************ control.c ******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical fermions 			*/
/* NEEDS UPGRADING TO ASQTAD */
/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

/* This is the version for the Schroedinger functional simulation */
/* 1/29/97 UMH */

#define CONTROL
#include "schroed_ks_includes.h"

int main(int argc, char *argv[])  {
int meascount,todo;
int prompt;
double dssplaq,dstplaq,ds_deta,bd_plaq;
Real rsq,rpbpp,rpbpm;
Real f_energy,f_pressure;	/* fermionic energy and pressure */
Real f_action;			/* fermionic action */
int m_iters,s_iters,avm_iters,avs_iters;
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

	/* perform measuring trajectories, reunitarizing and measuring 	*/
	meascount=0;		/* number of measurements 		*/
	ds_deta = 0.0;
	bd_plaq = 0.0;
	avm_iters = avs_iters = 0;
	for(todo=trajecs; todo > 0; --todo ){ 

	    /* do the trajectories */
	    s_iters=update();

	    /* call plaquette measuring process */
	    d_plaquette(&dssplaq,&dstplaq);
	    dssplaq *= -1.0; dstplaq *= -1.0;	/* KS phases change sign */

	    /* call the coupling measuring process */
	    if(bc_flag > 0){
		coupling(&ds_deta,&bd_plaq);
		ds_deta *= -1.0;	/* KS phases change sign */
		bd_plaq *= -1.0;	/* KS phases change sign */
	    }

	    /* generate a pseudofermion configuration */
	    grsource(EVEN);
	    /* do conjugate gradient to get (Madj M)inverse * phi  */
	    load_ferm_links(&fn_links);
	    m_iters=ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			       niter, nrestart, rsqmin, PRECISION, EVEN,
			       &rsq,  &fn_links);

	    /* call Psi-bar-Psi and fermion energy/pressure measurement */
	    /* it also measures the pseudofermionic action 		*/
	    f_measure(&rpbpp,&rpbpm,&f_energy,&f_pressure,&f_action);

	    avm_iters += m_iters;
	    avs_iters += s_iters;
	    ++meascount;

	    if(this_node==0)printf("GMES %e %e %e %e %e\n",
		ds_deta,bd_plaq,(double)m_iters,dssplaq,dstplaq);
	    /* dS/deta bd_plaq cg_iters ss_plaq st_plaq */

	    if(this_node==0)printf("FMES %e %e %e %e %e\n",(double)rpbpp,
		(double)rpbpm,(double)f_energy,(double)f_pressure,
		(double)f_action);

	    fflush(stdout);
	}	/* end loop over trajectories */

	if(this_node==0)printf("RUNNING COMPLETED\n");
	if(meascount>0)  {
	    if(this_node==0)printf("average cg iters for step= %e\n",
		(double)avs_iters/meascount);
	    if(this_node==0)printf("average cg iters for measurement= %e\n",
		(double)avm_iters/meascount);
	}

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);

	/* save lattice if requested */
	if( saveflag != FORGET ){
	    rephase_sf( OFF );
	    save_lattice( saveflag, savefile, stringLFN );
	}
    }
    return 0;
}
