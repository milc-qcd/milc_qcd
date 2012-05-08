/***************** control.c *****************************************/

/* Main procedure for clover spectrosopy 			*/
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
#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif

int main(int argc, char *argv[])
{
  int prompt;
  int i, j, k, iq0, iq1;
#ifdef CLOV_LEAN
  int oldiq0, oldiq1, oldip0;
#endif
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif
  wilson_prop_field *prop[MAX_PROP];
  wilson_prop_field *quark[MAX_QK];
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  
  starttime=dclock();
    
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");

  /* loop over input sets */

  while( readin(prompt) == 0){

    if(prompt == 2)continue;
    
    total_iters=0;
    
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif
    
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
    } else {
      savelat_p = NULL;
    }

    if(this_node==0)printf("END OF HEADER\n");
    
    /**************************************************************/


    /* Loop over the propagators */

    STARTTIME;
    for(i=0; i<param.num_prop; i++){
      
      /**************************************************************/
      /* Read and/or generate quark propagator */

      if(param.prop_type[i] == CLOVER_TYPE)
	{
	  int ncolor = convert_ksource_to_color(param.src_qs[i].nsource);
	  
	  prop[i] = create_wp_field(ncolor);
      
	  if(this_node==0)printf("Kappa= %g source %s residue= %g rel= %g\n",
				 (double)param.dcp[i].Kappa,
				 param.src_qs[i].descrp,
				 (double)param.qic[i].resid,
				 (double)param.qic[i].relresid);
	  
	  total_iters += get_wprop_to_wp_field(param.startflag_w[i], 
					       param.startfile_w[i], 
					       param.saveflag_w[i], 
					       param.savefile_w[i],
					       prop[i], 
					       &param.src_qs[i],
					       &param.qic[i], 
					       &param.dcp[i],
					       param.bdry_phase[i],
					       param.coord_origin,
					       param.check[i]);
#ifdef CLOV_LEAN
	  /* Free clover prop memory if we have saved the prop to disk */
	  if(param.saveflag_w[i] != FORGET){
	    free_wp_field(prop[i]);
	    clear_qs(&param.src_qs[i]);
	    node0_printf("destroy prop[%d]\n",i);
	  }
#endif
	}
	
      else if(param.prop_type[i] == KS_TYPE) /* KS_TYPE */
	{
	  prop[i] = create_wp_field(param.src_qs[i].ncolor);

	  if(this_node==0)printf("Mass= %g source %s residue= %g rel= %g\n",
				 (double)param.ksp[i].mass,
				 param.src_qs[i].descrp,
				 (double)param.qic[i].resid,
				 (double)param.qic[i].relresid);
	  
	  total_iters += get_ksprop_to_wp_field(param.startflag_ks[i], 
						param.startfile_ks[i], 
						param.saveflag_ks[i], 
						param.savefile_ks[i],
						prop[i], 
						&param.src_qs[i],
						&param.qic[i], 
						&param.ksp[i],
						param.bdry_phase[i],
						param.coord_origin,
						param.check[i]);
#ifdef CLOV_LEAN
	  /* (We don't free naive prop memory, since we don't save it
	     as a clover prop in get_ksprop_to_wp_field) */
#endif
	}
      else /* KS4_TYPE */
	{
	  prop[i] = create_wp_field(param.src_qs[i].ncolor);

	  if(this_node==0)printf("Mass= %g source %s residue= %g rel= %g\n",
				 (double)param.ksp[i].mass,
				 param.src_qs[i].descrp,
				 (double)param.qic[i].resid,
				 (double)param.qic[i].relresid);
	  
	  total_iters += get_ksprop4_to_wp_field(param.startflag_w[i], 
						 param.startfile_w[i], 
						 param.saveflag_w[i], 
						 param.savefile_w[i],
						 prop[i], 
						 &param.src_qs[i],
						 &param.qic[i], 
						 &param.ksp[i],
						 param.bdry_phase[i],
						 param.coord_origin,
						 param.check[i]);
	}
      
    } /* propagators */
    ENDTIME("compute propagators");

    /*****************************************************************/
    /* Complete the quark propagators by applying the sink operators
       to either the raw propagator or by building on an existing quark
       propagator */
    
    STARTTIME;
#ifdef CLOV_LEAN
    oldip0 = -1;
    oldiq0 = -1;
    oldiq1 = -1;
#endif
    for(j=0; j<param.num_qk; j++){
      i = param.prop_for_qk[j];

      if(param.parent_type[j] == PROP_TYPE){
#ifdef CLOV_LEAN
	/* Restore clover prop[i] from file. */
	/* But first destroy the old one, unless we still need it */
	if(oldip0 >= 0 && oldip0 != i)
	  if(param.prop_type[oldip0] == CLOVER_TYPE &&
	     param.saveflag_w[oldip0] != FORGET){
	    free_wp_field(prop[oldip0]);
	    node0_printf("destroy prop[%d]\n",oldip0);
	  }
	
	/* In this case we won't need any old quarks */
	if(oldiq0 >= 0)
	  if(param.saveflag_q[oldiq0] != FORGET){
	    free_wp_field(quark[oldiq0]);
	    node0_printf("destroy quark[%d]\n",oldiq0);
	  }

	if(oldiq1 >= 0)
	  if(param.saveflag_q[oldiq1] != FORGET){
	    free_wp_field(quark[oldiq1]);
	    node0_printf("destroy quark[%d]\n",oldiq1);
	  }

	if(prop[i]->swv[0] == NULL)
	  reread_wprop_to_wp_field(param.saveflag_w[i], param.savefile_w[i], prop[i]);
#endif
	/* Apply sink operator quark[j] <- Op[j] prop[i] */
	quark[j] = create_wp_field_copy(prop[i]);
	wp_sink_op(&param.snk_qs_op[j], quark[j]);
#ifdef CLOV_LEAN
	oldip0 = i;
	oldiq0 = -1;
#endif

	/* Can we delete any props now? */
	/* For each prop, scan ahead to see if it is no longer needed. */
	for(i = 0; i < param.num_prop; i++)
	  if(prop[i]->swv[0] != NULL){
	    int wont_need = 1;
	    for(k = j + 1; k < param.num_qk; k++)
	      if(param.parent_type[k] == PROP_TYPE &&
		 param.prop_for_qk[k] == i)
		wont_need = 0;  /* Still need this one */
	    if(wont_need){
	      free_wp_field(prop[i]);
	      node0_printf("destroy prop[%d]\n",i);
	    }
	  }
      }
      else { /* QUARK_TYPE */
#ifdef CLOV_LEAN
	/* Restore quark[i] from file */
	/* But first destroy the old ones, unless we still need one of them */

	/* In this case we won't need the old prop */
	if(oldip0 >= 0)
	   if(param.prop_type[oldip0] == CLOVER_TYPE &&
	      param.saveflag_w[oldip0] != FORGET){
	     free_wp_field(prop[oldip0]);
	     node0_printf("destroy prop[%d]\n",oldip0);
	   }

	if(oldiq0 >= 0 && oldiq0 != i)
	  if(param.saveflag_q[oldiq0] != FORGET){
	    free_wp_field(quark[oldiq0]);
	    node0_printf("destroy quark[%d]\n",oldiq0);
	  }

	if(oldiq1 >= 0 && oldiq1 != i)
	  if(param.saveflag_q[oldiq1] != FORGET){
	    free_wp_field(quark[oldiq1]);
	    node0_printf("destroy quark[%d]\n",oldiq1);
	  }

	if(quark[i]->swv[0] == NULL)
	  reread_wprop_to_wp_field(param.saveflag_q[i], param.savefile_q[i], quark[i]);
	
#endif
	/* Apply sink operator quark[j] <- Op[j] quark[i] */
	quark[j] = create_wp_field_copy(quark[i]);
	wp_sink_op(&param.snk_qs_op[j], quark[j]);
#ifdef CLOV_LEAN
	oldip0 = -1;
	oldiq0 = i;
#endif
      }
	
      /* Save the resulting quark[j] if requested */
      dump_wprop_from_wp_field( param.saveflag_q[j], 
				param.savefile_q[j], quark[j]);
#ifdef CLOV_LEAN
      oldiq1 = j;
#endif
    }

#ifdef CLOV_LEAN
    /* Free remaining memory */
    if(oldip0 >= 0)
       if(param.prop_type[oldip0] == CLOVER_TYPE &&
	  param.saveflag_w[oldip0] != FORGET){
	 free_wp_field(prop[oldip0]);
	 node0_printf("destroy prop[%d]\n",oldip0);
       }
    
    if(oldiq0 >= 0)
      if(param.saveflag_q[oldiq0] != FORGET){
	free_wp_field(quark[oldiq0]);
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    
    if(oldiq1 >= 0)
      if(param.saveflag_q[oldiq1] != FORGET){
	free_wp_field(quark[oldiq1]);
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
#endif

    /* Now destroy all remaining propagator fields */

    for(i = 0; i < param.num_prop; i++){
      if(prop[i] != NULL)node0_printf("destroy prop[%d]\n",i);
      destroy_wp_field(prop[i]);
      prop[i] = NULL;
    }
    ENDTIME("generate quarks");
    
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
	  free_wp_field(quark[oldiq0]);
	  node0_printf("destroy quark[%d]\n",oldiq0);
	}
      
      if(i > 0 && oldiq1 != iq0 && oldiq1 != iq1)
	if(param.saveflag_q[oldiq1] != FORGET){
	  free_wp_field(quark[oldiq1]);
	  node0_printf("destroy quark[%d]\n",oldiq1);
	}

      if(quark[iq0]->swv[0] == NULL)
	reread_wprop_to_wp_field(param.saveflag_q[iq0], param.savefile_q[iq0], quark[iq0]);

      if(quark[iq1]->swv[0] == NULL){
	reread_wprop_to_wp_field(param.saveflag_q[iq1], param.savefile_q[iq1], quark[iq1]);
      }
#endif

      /* Tie together to generate hadron spectrum */
      spectrum_cl(quark[iq0], quark[iq1], i);

      /* Remember, in case we need to free memory */
#ifdef CLOV_LEAN
      oldiq0 = iq0;
      oldiq1 = iq1;
#endif
    }
#ifdef CLOV_LEAN
    /* Free any remaining quark prop memory */
    if(quark[oldiq0]->swv[0] != NULL)
      if(param.saveflag_q[oldiq0] != FORGET){
	free_wp_field(quark[oldiq0]);
	node0_printf("destroy quark[%d]\n",oldiq0);
      }
    if(quark[oldiq1]->swv[0] != NULL)
      if(param.saveflag_q[oldiq1] != FORGET){
	free_wp_field(quark[oldiq1]);
	node0_printf("destroy quark[%d]\n",oldiq1);
      }
#endif
    ENDTIME("tie hadron correlators");

    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
    printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
    fflush(stdout);

    for(i = 0; i < param.num_qk; i++){
      if(quark[i] != NULL)node0_printf("destroy quark[%d]\n",i);
      destroy_wp_field(quark[i]);
    }

    destroy_ape_links_3D(ape_links);

    /* Destroy fermion links (possibly created in make_prop()) */

#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;

  } /* readin(prompt) */

#ifdef HAVE_QUDA
  qudaFinalize();
#endif

  return 0;
}
