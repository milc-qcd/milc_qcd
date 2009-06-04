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
  int i, j, iq0, iq1, oldiq0, oldiq1, oldip0;
  double starttime, endtime, dtime;
  wilson_prop_field prop[MAX_PROP];
  wilson_prop_field quark[MAX_QK];
  
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
	invalidate_this_clov(gen_clov);
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

    STARTTIME;
    for(i=0; i<param.num_prop; i++){
      
      /**************************************************************/
      /* Read and/or generate quark propagator */

      prop[i] = create_wp_field();

      if(prop[i] == NULL){
	printf("main(%d): No room for prop\n",this_node);
	terminate(1);
      }
      
      
      if(param.prop_type[i] == CLOVER_TYPE)
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
					       prop[i], 
					       &param.src_wqs[i], &param.qic, 
					       &param.dcp[i],
					       param.check[i]);
#ifdef CLOV_LEAN
	  /* Free clover prop memory if we have saved the prop to disk */
	  if(param.saveflag_w[i] != FORGET){
	    destroy_wp_field(prop[i]); prop[i] = NULL;
	    clear_wqs(&param.src_wqs[i]);
	    node0_printf("destroy prop[%d]\n",i);
	  }
#endif
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
						prop[i], 
						&param.src_ksqs[i], &param.qic, 
						&param.ksp[i],
						param.check[i]);
#ifdef CLOV_LEAN
	  /* (We don't free naive prop memory, since we don't save it
	     as a clover prop in get_ksprop_to_wp_field) */
#endif
	}
      
    } /* propagators */
    ENDTIME("compute propagators");

    /*****************************************************************/
    /* Complete the quark propagators by applying the sink operators
       to either the raw propagator or by building on an existing quark
       propagator */
    
    oldip0 = -1;
    oldiq0 = -1;
    oldiq1 = -1;
    for(j=0; j<param.num_qk; j++){
      STARTTIME;
      i = param.prop_for_qk[j];

      if(param.parent_type[j] == PROP_TYPE){
#ifdef CLOV_LEAN
	/* Restore clover prop[i] from file. */
	/* But first destroy the old one, unless we still need it */
	if(oldip0 >= 0 && oldip0 != i)
	  if(param.prop_type[oldip0] == CLOVER_TYPE &&
	     param.saveflag_w[oldip0] != FORGET){
	    destroy_wp_field(prop[oldip0]);  prop[oldip0] = NULL;
	    node0_printf("destroy prop[%d]\n",oldip0);
	  }
	
	/* In this case we won't need any old quarks */
	if(oldiq0 >= 0)
	  if(param.saveflag_q[oldiq0] != FORGET){
	    destroy_wp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq0);
	  }

	if(oldiq1 >= 0)
	  if(param.saveflag_q[oldiq1] != FORGET){
	    destroy_wp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq1);
	  }

	if(prop[i] == NULL)
	  prop[i] = reread_wprop_to_wp_field(param.saveflag_w[i], 
					     param.savefile_w[i]);
#endif
	/* Apply sink operator quark[j] <- Op[j] prop[i] */
	quark[j] = w_sink_op(&param.snk_wqs[j], prop[i]);
	oldip0 = i;
	oldiq0 = -1;
      }
      else { /* QUARK_TYPE */
#ifdef CLOV_LEAN
	/* Restore quark[i] from file */
	/* But first destroy the old ones, unless we still need one of them */

	/* In this case we won't need the old prop */
	if(oldip0 >= 0)
	   if(param.prop_type[oldip0] == CLOVER_TYPE &&
	      param.saveflag_w[oldip0] != FORGET){
	     destroy_wp_field(prop[oldip0]); prop[oldip0] = NULL;
	     node0_printf("destroy prop[%d]\n",oldip0);
	   }

	if(oldiq0 >= 0 && oldiq0 != i)
	  if(param.saveflag_q[oldiq0] != FORGET){
	    destroy_wp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq0);
	  }

	if(oldiq1 >= 0 && oldiq1 != i)
	  if(param.saveflag_q[oldiq1] != FORGET){
	    destroy_wp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	    node0_printf("destroy quark[%d]\n",oldiq1);
	  }

	if(quark[i] == NULL)
	  quark[i] = reread_wprop_to_wp_field(param.saveflag_q[i], 
					      param.savefile_q[i]);
	
#endif
	/* Apply sink operator quark[j] <- Op[j] quark[i] */
	quark[j] = w_sink_op(&param.snk_wqs[j], quark[i]);
	oldip0 = -1;
	oldiq0 = i;
      }
	
      /* Save the resulting quark[j] if requested */
      dump_wprop_from_wp_field( param.saveflag_q[j], 
				param.savefile_q[j], quark[j]);
      oldiq1 = j;
      ENDTIME("generate sink operator");
    }
#ifdef CLOV_LEAN
    /* Free remainig memory */
    if(oldip0 >= 0)
       if(param.prop_type[oldip0] == CLOVER_TYPE &&
	  param.saveflag_w[oldip0] != FORGET){
	 destroy_wp_field(prop[oldip0]); prop[oldip0] = NULL;
	 node0_printf("destroy prop[%d]\n",oldip0);
       }
    
    if(oldiq0 >= 0)
      if(param.saveflag_q[oldiq0] != FORGET){
	destroy_wp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    
    if(oldiq1 >= 0)
      if(param.saveflag_q[oldiq1] != FORGET){
	destroy_wp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
#endif

    /* Now destroy all remaining propagator fields */

    for(i = 0; i < param.num_prop; i++){
      if(prop[i] != NULL)node0_printf("destroy prop[%d]\n",i);
      destroy_wp_field(prop[i]);
      prop[i] = NULL;
    }
    
    /****************************************************************/
    /* Compute the meson propagators */

    STARTTIME;
    for(i = 0; i < param.num_pair; i++){

      /* Index for the quarks making up this meson */
      iq0 = param.qkpair[i][0];
      iq1 = param.qkpair[i][1];

      node0_printf("Mesons for quarks %d and %d\n",iq0,iq1);

#ifdef CLOV_LEAN
      /* Restore quarks from file and free old memory */
      /* We try to reuse props that are already in memory, so we don't
         destroy them immediately, but wait to see if we need
         them again for the next pair. */
      if(i > 0 && oldiq0 != iq0 && oldiq0 != iq1)
	if(param.saveflag_q[oldiq0] != FORGET){
	  destroy_wp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq0);
	}
      
      if(i > 0 && oldiq1 != iq0 && oldiq1 != iq1)
	if(param.saveflag_q[oldiq1] != FORGET){
	  destroy_wp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	  node0_printf("destroy quark[%d]\n",oldiq1);
	}

      if(quark[iq0] == NULL)
	quark[iq0] = 
	  reread_wprop_to_wp_field(param.saveflag_q[iq0], 
				   param.savefile_q[iq0]);

      if(quark[iq1] == NULL){
	quark[iq1] = 
	  reread_wprop_to_wp_field(param.saveflag_q[iq1], 
				   param.savefile_q[iq1]);
      }
#endif

      /* Tie together to generate hadron spectrum */
      spectrum_cl(quark[iq0], quark[iq1], i);

      /* Remember, in case we need to free memory */
      oldiq0 = iq0;
      oldiq1 = iq1;
    }
#ifdef CLOV_LEAN
    /* Free any remaining quark prop memory */
    if(quark[oldiq0] != NULL)
      if(param.saveflag_q[oldiq0] != FORGET){
	destroy_wp_field(quark[oldiq0]); quark[oldiq0] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    if(quark[oldiq1] != NULL)
      if(param.saveflag_q[oldiq1] != FORGET){
	destroy_wp_field(quark[oldiq1]); quark[oldiq1] = NULL;
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
#endif
    ENDTIME("tie hadron correlators");

    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
    fflush(stdout);

    for(i = 0; i < param.num_qk; i++){
      if(quark[i] != NULL)node0_printf("destroy quark[%d]\n",i);
      destroy_wp_field(quark[i]); quark[i] = NULL;
    }

    destroy_ape_links_3D();
  } /* readin(prompt) */

  return 0;
}
