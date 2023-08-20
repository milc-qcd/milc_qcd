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

/* Apply momentum twist of parent to ape_links */
void momentum_twist_ape_links(int i, int sign){
  int ip = get_eve(i);
  Real *bp = param.bdry_phase[ip];
  Real mybdry_phase[4] = {bp[0], bp[1], bp[2], 0.};
  momentum_twist_links(mybdry_phase, sign, ape_links);
}

int main(int argc, char *argv[])
{
  char myname[] = "main";
  int prompt;
  int i, is, j, k, iq0, iq1;
#ifdef CLOV_LEAN
  int oldiq0, oldiq1, oldip0;
#endif
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif
  ks_prop_field *v_source[MAX_SOURCE];
  wilson_prop_field *wv_source[MAX_SOURCE];
  wilson_prop_field *prop[MAX_PROP];
  wilson_prop_field *quark[MAX_QK];
  
  initialize_machine(&argc,&argv);

  for(i = 0; i < MAX_SOURCE; i++)v_source[i] = NULL;
  for(i = 0; i < MAX_SOURCE; i++)wv_source[i] = NULL;
  for(i = 0; i < MAX_PROP; i++)prop[i] = NULL;
  for(i = 0; i < MAX_QK; i++)quark[i] = NULL;

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

	/* (Re)construct APE smeared links after gauge fixing.  
	   No KS phases here! */
	destroy_ape_links_4D(ape_links);
	ape_links = ape_smear_4D( param.staple_weight, param.ape_iter );
	if(param.time_bc == 0)apply_apbc( ape_links, param.coord_origin[3] );
	refresh_ape_links = 1;
	ape_links_ks_phases = OFF;
	
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

    /* Create sources */

    STARTTIME;

    /* Base sources */
    for(k=0; k<param.num_base_source; k++){
      quark_source *qs = &param.base_src_qs[k];

      if(qs->field_type == WILSON_FIELD){
	
	wv_source[k] = create_wp_field(qs->ncolor);
	
	if(qs->saveflag != FORGET){
	  char *fileinfo = create_w_QCDML();
	  w_source_open_dirac(qs, fileinfo);
	  free(fileinfo);
	}
	
	/* Make the source */
	wilson_vector *src = create_wv_field();
	
	for(int ksource = 0; ksource < qs->nsource; ksource++){
	  int spin = convert_ksource_to_spin(ksource);
	  int color = convert_ksource_to_color(ksource);
	  
	  
	  /* Create a base source */
	  if(wv_source_field(src, qs)){
	    printf("%s(%d): error getting source\n",myname,this_node);
	    terminate(1);
	  }

	  /* Copy to source[k] */
	  insert_swv_from_wv(wv_source[k]->swv[color], spin, src);
	    
	  /* Write the source, if requested */
	  if(qs->saveflag != FORGET){
	    if(w_source_dirac( src, qs ) != 0)
	      node0_printf("%s: Error writing source\n",myname);
	  }
	} /* ksource */
  
	destroy_wv_field(src);
	  
	if(qs->saveflag != FORGET) w_source_close(qs);

      }
      else if(qs->field_type == KS_FIELD){
	v_source[k] = create_ksp_field(qs->ncolor);
	
	if(qs->saveflag != FORGET){
	  char *fileinfo = create_w_QCDML();
	  w_source_open_ks(qs, fileinfo);
	  free(fileinfo);
	}
	
	for(int color = 0; color < qs->ncolor; color++){
	  
	  /* Create a base source */
	  if(v_source_field(v_source[k]->v[color], qs)){
	    printf("%s(%d): error getting source\n",myname,this_node);
	    terminate(1);
	  }
	  
	  /* Write the source, if requested */
	  if(qs->saveflag != FORGET){
	    qs->color = color;
	    if(w_source_ks( v_source[k]->v[color], qs ) != 0)
	      node0_printf("Error writing source\n");
	  }
	} /* color */
	
	if(qs->saveflag != FORGET) w_source_close(qs);
      } /* qs->field_type == KS_FIELD */
      else {
	node0_printf("%s: Unrecognized source field type %d\n", myname, qs->field_type);
	terminate(1);
      }
    } /* k */


    /* Modified sources */
    for(is=param.num_base_source; is<param.num_base_source+param.num_modified_source; is++){

      quark_source *qs = &param.base_src_qs[is];

      if(qs->field_type == WILSON_FIELD){
	wv_source[is] = create_wp_field(qs->ncolor);

	if(qs->saveflag != FORGET){
	  char *fileinfo = create_w_QCDML();
	  w_source_open_dirac(qs, fileinfo);
	  free(fileinfo);
	}
	
	/* Copy parent source */
	int p = param.parent_source[is];
	if(wv_source[p] == NULL){
	  node0_printf("Parent source %d is missing\n", p);
	  terminate(1);
	}
	copy_wp_field(wv_source[is],  wv_source[p]);
	
	for(int ksource = 0; ksource < qs->nsource; ksource++){
	  int spin = convert_ksource_to_spin(ksource);
	  int color = convert_ksource_to_color(ksource);
	  
	  /* Apply operator*/
	  wilson_vector *src = create_wv_field();
	  extract_wv_from_swv(src, wv_source[is]->swv[color], spin);
	  wv_field_op(src, qs->op, qs->subset, qs->t0);
	  insert_swv_from_wv(wv_source[is]->swv[color], spin, src);
	  
	  destroy_wv_field(src);
	} /* ksource */
	
	if(qs->saveflag != FORGET) w_source_close(qs);
      }

      else if(qs->field_type == KS_FIELD){
	v_source[k] = create_ksp_field(qs->ncolor);

	if(qs->saveflag != FORGET){
	  char *fileinfo = create_ks_XML();
	  w_source_open_ks(qs, fileinfo);
	  free(fileinfo);
	}
	
	/* Copy parent source */
	int p = param.parent_source[is];
	if(v_source[p] == NULL){
	  node0_printf("Parent source %d is missing\n", p);
	  terminate(1);
	}
	copy_ksp_field(v_source[is],  v_source[p]);
	
	for(int color = 0; color < qs->ncolor; color++){
	  
	  /* Apply operator*/
	  v_field_op(v_source[is]->v[color], qs->op, qs->subset, qs->t0);
	  
	  /* Write the source, if requested */
	  if(qs->saveflag != FORGET){
	    qs->color = color;
	    if(w_source_ks( v_source[is]->v[color], qs ) != 0)
	      node0_printf("Error writing source\n");
	  }
	} /* color */
	
	if(qs->saveflag != FORGET) w_source_close(qs);
      }
      else{
	node0_printf("%s: Unrecognized source field type %d\n", myname, qs->field_type);
	terminate(1);
      }

    } /* is */

    ENDTIME("Create sources");

    /**************************************************************/


    /* Loop over the propagators */

    STARTTIME;
    for(i=0; i<param.num_prop; i++){
      node0_printf("******* Creating propagator %d ********\n",i);fflush(stdout);
      
      /**************************************************************/
      /* Read and/or generate quark propagator */

      if(param.prop_type[i] == CLOVER_TYPE)
	{
	  int ncolor = convert_ksource_to_color(param.src_qs[i].nsource);

	  /* Convert KS source to Dirac source if need be */
	  is = param.source[i];
	  wilson_prop_field *tmp_src = wv_source[is];
	  if (param.base_src_qs[is].field_type == KS_FIELD){
	    tmp_src = create_wp_field(ncolor);
	    for(int color = 0; color < ncolor; color++)
	      for(int spin = 0; spin < 4; spin++)
		insert_swv_from_v(tmp_src->swv[color], spin, spin,
				  v_source[is]->v[color]);
	  }

	  prop[i] = create_wp_field(ncolor);
      
	  node0_printf("Generate Dirac propagator\n");
	  node0_printf("Kappa= %g source %s residue= %g rel= %g\n",
		       (double)param.dcp[i].Kappa,
		       param.src_qs[i].descrp,
		       (double)param.qic[i].resid,
		       (double)param.qic[i].relresid);
	  
	  /* For clover_info */
	  wqstmp = param.src_qs[i];
	  dcptmp = param.dcp[i];

	  total_iters += get_wprop_to_wp_field(param.prop_type[i],
					       param.startflag_w[i], 
					       param.startfile_w[i], 
					       param.saveflag_w[i], 
					       param.savefile_w[i],
					       prop[i],
					       tmp_src,
					       &param.src_qs[i],
					       &param.qic[i], 
					       (void *)&param.dcp[i],
					       param.bdry_phase[i],
					       param.coord_origin,
					       param.check[i]);


	  if (param.base_src_qs[is].field_type == KS_FIELD)
	    destroy_wp_field(tmp_src);

#ifdef CLOV_LEAN
	  /* Free clover prop memory if we have saved the prop to disk */
	  if(param.saveflag_w[i] != FORGET){
	    free_wp_field(prop[i]);
	    clear_qs(&param.src_qs[i]);
	    node0_printf("destroy prop[%d]\n",i);
	  }
#endif
	}
	
      /* ------------------------------------------- */
      else if(param.prop_type[i] == IFLA_TYPE)
	{
	  int ncolor = convert_ksource_to_color(param.src_qs[i].nsource);
	  

	  prop[i] = create_wp_field(ncolor);
      
	  node0_printf("Generate Dirac IFLA propagator\n");
	  if(this_node==0)printf("Kappa= %g source %s residue= %g rel= %g\n",
				 (double)param.nap[i].kapifla,
				 param.src_qs[i].descrp,
				 (double)param.qic[i].resid,
				 (double)param.qic[i].relresid);
	 
	  /* For clover_info */
	  wqstmp = param.src_qs[i];
	  naptmp = param.nap[i];
	  
	  total_iters += get_wprop_to_wp_field(param.prop_type[i],
					       param.startflag_w[i], 
					       param.startfile_w[i], 
					       param.saveflag_w[i], 
					       param.savefile_w[i],
					       prop[i], 
					       wv_source[param.source[i]],
					       &param.src_qs[i], 
					       &param.qic[i], 
					       (void *)&param.nap[i],
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
      else if(param.prop_type[i] == KS_TYPE || param.prop_type[i] == KS0_TYPE ) /* KS_TYPE */
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
						v_source[param.source[i]],
						&param.src_qs[i],
						&param.qic[i], 
						&param.ksp[i],
						param.bdry_phase[i],
						param.coord_origin,
						param.check[i],
						param.prop_type[i] == KS0_TYPE);
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
						 v_source[param.source[i]],
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
      node0_printf("******* Creating quark %d ********\n",j); fflush(stdout);
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
	/* Before applying operator, apply momentum twist to
	   ape_links.  Use the momentum twist of the parent
	   propagator */
	momentum_twist_ape_links(i, +1);
	/* Apply sink operator quark[j] <- Op[j] prop[i] */
	quark[j] = create_wp_field_copy(prop[i]);
	wp_sink_op(&param.snk_qs_op[j], quark[j]);
	/* Remove twist */
	momentum_twist_ape_links(i, -1);
#ifdef CLOV_LEAN
	oldip0 = i;
	oldiq0 = -1;
#endif
      }
      else if(param.parent_type[j] == QUARK_TYPE) { /* QUARK_TYPE */
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
	momentum_twist_ape_links(i, +1);
	quark[j] = create_wp_field_copy(quark[i]);
	wp_sink_op(&param.snk_qs_op[j], quark[j]);
	momentum_twist_ape_links(i, -1);
#ifdef CLOV_LEAN
	oldip0 = -1;
	oldiq0 = i;
#endif
      } else { /* COMBO_TYPE */
	int k;
	int nc = quark[param.combo_qk_index[j][0]]->nc;
	/* Create a zero field */
	quark[j] = create_wp_field(nc);
	/* Compute the requested linear combination */
	for(k = 0; k < param.num_combo[j]; k++){
	  wilson_prop_field *q = quark[param.combo_qk_index[j][k]];
	  if(nc != q->nc){
	    printf("Error: Attempting to combine an inconsistent number of colors: %d != %d\n",nc, q->nc);
	    terminate(1);
	  }
	  scalar_mult_add_wprop_field(quark[j], q, param.combo_coeff[j][k], quark[j]);
	}
      }
	
      /* Save the resulting quark[j] if requested */
      dump_wprop_from_wp_field( param.saveflag_q[j], param.savetype_q[j],
				param.savefile_q[j], quark[j]);

      /* Can we delete any props and quarks now? */
      /* If nothing later depends on a prop or quark, free it up. */
      for(i = 0; i < param.num_prop; i++)
	if( prop[i]->swv[0] != NULL  &&  param.prop_dep_qkno[i] < j ){
	  free_wp_field(prop[i]);
	  node0_printf("free prop[%d]\n",i);
	}
      
      for(i = 0; i < j; i++)
	if( quark[i]->swv[0] != NULL  &&  param.quark_dep_qkno[i] < j ){
	  free_wp_field(quark[i]);
	  node0_printf("free quark[%d]\n",i);
	}
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
      quark[i] = NULL;
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

  free_lattice();

#ifdef HAVE_QUDA
  finalize_quda();
#endif
  normal_exit(0);

  return 0;
}
