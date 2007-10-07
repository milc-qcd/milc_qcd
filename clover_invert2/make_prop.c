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
    printf("get_wprop_to_field(%d): No room for dst\n",this_node);
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
      
      if(startflag == FRESH){
	my_qic->start_flag = START_ZERO_GUESS;
      } else {
	my_qic->start_flag = START_NONZERO_GUESS;      
      }
      
      /* Complete the source structure */
      my_wqs->color = color;
      my_wqs->spin = spin;
      
      if(check){  /* Provision for suppressing checking */
	
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

