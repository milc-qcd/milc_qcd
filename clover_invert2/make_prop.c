/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "cl_inv_includes.h"

int get_wprop_to_wp_field(int startflag, char startfile[], 
			  int saveflag, char savefile[],
			  wilson_prop_field wp,
			  wilson_quark_source *my_wqs,
			  quark_invert_control *my_qic,
			  dirac_clover_param *my_dcp,
			  int check)
{
  int color, spin;
  int cl_cg = CL_CG;
  int avs_iters;
  int tot_iters = 0;
  int status;
  wilson_vector *dst;
  w_prop_file *fp_in, *fp_out; 

  dst = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(dst == NULL){
    printf("get_wprop_to_wp_field(%d): No room for dst\n",this_node);
    terminate(1);
  }

  /* For clover_info */
  wqstmp = *my_wqs;
  dcptmp = *my_dcp;

  /* Open files for Wilson propagators, if requested */
  fp_in  = r_open_wprop(startflag, startfile);
  fp_out = w_open_wprop(saveflag,  savefile, my_wqs->type);

  /* Loop over source colors and spins */
  for(spin=0;spin<4;spin++)
    for(color=0;color<3;color++){
      
	
      status = reload_wprop_sc_to_field(startflag, fp_in, my_wqs, 
					spin, color, dst, 1);
      
      /* (Re)construct propagator */
      
      if(startflag == FRESH)my_qic->start_flag = START_ZERO_GUESS;
      else                  my_qic->start_flag = START_NONZERO_GUESS;      
      
      /* Complete the source structure */
      my_wqs->color = color;
      my_wqs->spin = spin;
      
      if(check || startflag == FRESH){  /* Provision for suppressing checking */
	
	/* solve for dst */

	switch (cl_cg) {
	case BICG:
	  avs_iters =
	    (Real)wilson_invert_field_wqs(my_wqs, w_source_field, dst,
					  bicgilu_cl_field,
					  my_qic,(void *)my_dcp);
	  break;
	case HOP:
	  avs_iters = 
	    (Real)wilson_invert_field_wqs(my_wqs, w_source_field, dst,
					  hopilu_cl_field,
					  my_qic,(void *)my_dcp);
	  break;
	case MR:
	  avs_iters = 
	    (Real)wilson_invert_field_wqs(my_wqs, w_source_field, dst,
					  mrilu_cl_field,
					  my_qic,(void *)my_dcp);
	  break;
	case CG:
	  avs_iters = 
	    (Real)wilson_invert_field_wqs(my_wqs, w_source_field, dst,
					  cgilu_cl_field,
					  my_qic,(void *)my_dcp);
	  break;
	default:
	  node0_printf("main(%d): Inverter choice %d not supported\n",
		       this_node,cl_cg);
	}
      }
      
      /* Copy solution vector to wilson_prop_field wp */
      copy_wp_from_wv( wp, dst, color, spin);
      
      /* save solution if requested */
      save_wprop_sc_from_field( saveflag, fp_out, my_wqs, spin, color, dst, 
				"", 1);
      
      tot_iters += avs_iters;
    } /* source color, spin */
  
  /* close files for wilson propagators */
  r_close_wprop(startflag, fp_in);
  w_close_wprop(saveflag,  fp_out);

  free(dst);

  return tot_iters;
}

int get_ksprop_to_wp_field(int startflag, char startfile[], 
			   int saveflag, char savefile[],
			   wilson_prop_field wp,
			   ks_quark_source *my_ksqs,
			   quark_invert_control *my_qic,
			   ks_param *my_ksp,
			   int check)
{
  int color;
  int avs_iters;
  int tot_iters = 0;
  int status;
  int ks_source_r[4];
  su3_vector *dst;
  ks_prop_file *fp_in, *fp_out; 

  dst = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(dst == NULL){
    printf("get_ksprop_to_wp_field(%d): No room for dst\n",this_node);
    terminate(1);
  }

  /* For ksprop_info */
  ksqstmp = *my_ksqs;
  ksptmp  = *my_ksp;

  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);
  fp_out = w_open_ksprop(saveflag,  savefile, my_ksqs->type);

  /* Loop over source colors */
  rephase( ON);
  for(color=0;color<3;color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, my_ksqs, 
				      color, dst, 1);
    /* (Re)construct propagator */
    
    if(startflag == FRESH) my_qic->start_flag = START_ZERO_GUESS;
    else                   my_qic->start_flag = START_NONZERO_GUESS;      
    
    /* Complete the source structure */
    my_ksqs->color = color;
    
    /* Check the solution */
    if(check || startflag == FRESH){
      if(make_path_table(&ks_act_paths, NULL, my_ksp->mass))
	invalidate_all_ferm_links(&fn_links);
      load_ferm_links(&fn_links, &ks_act_paths);

      /* In most use cases we will be reading a precomputed staggered
	 propagator, so we use the less optimized mat_invert_cg_field
	 algorithm, instead of mat_invert_uml_field here to avoid
	 "reconstructing", and so overwriting the odd-site solution.
	 This would be a degradation if the propagator was precomputed
	 in double precision, and we were doing single precision
	 here. If we start using this code to compute the propagator
	 from a fresh start, we could make the selection of the
	 inverter algorithm dependent on the propagator start flag.
	 FRESH -> uml and otherwise cg. */
      avs_iters =
	(Real)ks_invert_ksqs(my_ksqs, ks_source_field, dst,
			     mat_invert_cg_field,
			     my_qic, my_ksp->mass, &fn_links);
    }
    
    /* save solution if requested */
    save_ksprop_c_from_field( saveflag, fp_out, my_ksqs, color, dst, "", 1);
    
    /* Convert KS prop to naive prop (su3_vector maps to
       spin_wilson_vector) for a given source color */
    ks_source_r[0] = my_ksqs->x0;  ks_source_r[1] = my_ksqs->y0;
    ks_source_r[2] = my_ksqs->z0;  ks_source_r[3] = my_ksqs->t0;
    
    convert_ksprop_to_wprop_swv(wp[color], dst, ks_source_r);
    
    tot_iters += avs_iters;
  } /* source color */
  rephase( OFF);
  
  /* close files for staggered propagators */
  r_close_ksprop(startflag, fp_in);
  w_close_ksprop(saveflag,  fp_out);

  free(dst);

  return tot_iters;
}

