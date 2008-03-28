/***************** control.c *****************************************/

/* Main procedure for quenched SU3 clover fermions 			*/
/* MIMD version 7 */

/* This version computes propagators for clover fermions on a supplied
   background field config and ties them together according to the
   input parameters. */

/* Modifications ...
   
   8/8/98  Rearranged loop so hadron propagators are written as soon
           as they are calculated. C.D.
   8/10/96 Installed new propagator IO and added timing C.D.
   5/29/07 Generalized the algorithm C.D.
 */

/* Comment out if you want to suppress detailed timing */

#define CONTROL
#include "cl_inv_includes.h"
#include <string.h>
#ifdef HAVE_QDP
#include "lattice_qdp.h"
#include <qdp.h>
#endif

int main(int argc, char *argv[])
{
  int prompt;
  int i, iq0, iq1;
  double starttime, endtime, dtime;
  wilson_prop_field quark_prop[MAX_QK];
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();
  /* loop over input sets */

  while( readin(prompt) == 0){
    
    starttime=dclock();
    
    total_iters=0;
    
    /**************************************************************/
    /* Set up gauge field */
    
    if( param.fixflag == COULOMB_GAUGE_FIX)
      {
	if(this_node == 0) 
	  printf("Fixing to Coulomb gauge\n");

	STARTTIME;
	gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);
	ENDTIME("gauge fix");
      }
    else
      if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");
    
    /* save lattice if requested */
    if( param.saveflag != FORGET ){
      savelat_p = save_lattice( param.saveflag, param.savefile, 
				param.stringLFN );
    }

    if(this_node==0)printf("END OF HEADER\n");
    
    /* Loop over quarks */

    for(i=0; i<param.num_qk; i++){
      
      /**************************************************************/
      /* Read and/or generate quark propagator */

      quark_prop[i] = create_wp_field();

      if(quark_prop[i] == NULL){
	printf("main(%d): No room for prop\n",this_node);
	terminate(1);
      }
      
      
      if(param.qk_type[i] == CLOVER_TYPE)
	{
	  
	  if(this_node==0)printf("Kappa= %g source %s residue= %g rel= %g\n",
				 (double)param.dcp[i].Kappa,
				 param.src_wqs[i].descrp,
				 (double)param.qic.resid,
				 (double)param.qic.relresid);
	  
	  total_iters += get_wprop_to_wp_field(param.startflag_w[i], 
					       param.startfile_w[i], 
					       param.saveflag_w[i], 
					       param.savefile_w[i],
					       quark_prop[i], 
					       &param.src_wqs[i], &param.qic, 
					       &param.dcp[i],
					       param.check[i]);
	}
	
      else /* KS_PROP */
	{

	  if(this_node==0)printf("Mass= %g source %s residue= %g rel= %g\n",
				 (double)param.ksp[i].mass,
				 param.src_ksqs[i].descrp,
				 (double)param.qic.resid,
				 (double)param.qic.relresid);
	  
	  total_iters += get_ksprop_to_wp_field(param.startflag_ks[i], 
						param.startfile_ks[i], 
						param.saveflag_ks[i], 
						param.savefile_ks[i],
						quark_prop[i], 
						&param.src_ksqs[i], &param.qic, 
						&param.ksp[i],
						param.check[i]);
	}
      
    } /* quarks */

    /* Compute the meson propagators */

    for(i = 0; i < param.num_pair; i++){

      /* Index for the quarks making up this meson */
      iq0 = param.qkpair[i][0];
      iq1 = param.qkpair[i][1];

      node0_printf("Mesons for quarks %d and %d\n",iq0,iq1);

      spectrum_cl(quark_prop[iq0], quark_prop[iq1], i);
    }

    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
    fflush(stdout);

    for(i = 0; i < param.num_qk; i++)
      destroy_wp_field(quark_prop[i]);
  } /* readin(prompt) */

  return 0;
}
