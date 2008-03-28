/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "pw_nr_meson_includes.h"

static void convert_wprop_to_block_pauli(block_pauli_propagator *bp,
					 wilson_propagator *wp)
{
  int i;
  site *s;
  int c0,s0,c1,s1;

  /* Convert from MILC to FNAL basis */
  convert_wprop_milc_to_fnal_field(wp);

  /* Then put the diagonal spin blocks in the up and dn members of
     the Pauli propagator */
  
  FORALLSITES(i,s){
    for(c0 = 0; c0 < 3; c0++)
      for(s0 = 0; s0 < 2; s0++)
	for(s1 = 0; s1 < 2; s1++)
	  for(c1 = 0; c1 < 3; c1++){
	    bp[i].up.c[c0].d[s0].d[s1].c[c1] = 
	      wp[i].c[c0].d[s0].d[s1].c[c1];
	    bp[i].dn.c[c0].d[s0].d[s1].c[c1] = 
	      wp[i].c[c0].d[s0+2].d[s1+2].c[c1];
	  }
  }
}

int get_wprop_to_field(int startflag, char startfile[], 
		      int saveflag, char savefile[],
		      block_pauli_propagator bp[],
		      wilson_quark_source *my_wqs,
		      quark_invert_control *my_qic,
		      dirac_clover_param *my_dcp)
{
  int i;
  site *s;
  int color, spin;
  int cl_cg = CL_CG;
  int avs_iters;
  int tot_iters = 0;
  wilson_propagator *wp;
  wilson_vector *dst;

  /* Make space for propagator field */

  wp = (wilson_propagator *)malloc(sites_on_node*sizeof(wilson_propagator));
  if(wp == NULL){
    printf("get_wprop_to_field(%d): No room for prop\n",this_node);
    terminate(1);
  }

  /* Load quark prop in MILC basis */

  reload_wprop_to_field(startflag, startfile, my_wqs, wp, 1);

  /* (Re)construct propagator */

  dst = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(dst == NULL){
    printf("get_wprop_to_field(%d): No room for dst\n",this_node);
    terminate(1);
  }
  
  for(color = 0; color < 3; color++)
    for(spin = 0; spin < 4; spin++)
      {
	if(startflag == FRESH){
	  my_qic->start_flag = START_ZERO_GUESS;
	} else {
	  my_qic->start_flag = START_NONZERO_GUESS;      
	
	  /* Copy prop to make initial guess */
	  FORALLSITES(i,s){
	    dst[i] = wp[i].c[color].d[spin];
	  }
	}
	  
	/* Complete the source structure */
	my_wqs->color = color;
	my_wqs->spin = spin;
	
	/* compute the propagator.  Result in wp. */
	
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

	/* Store result in wp */
	FORALLSITES(i,s){
	  wp[i].c[color].d[spin] = dst[i];
	}
	tot_iters += avs_iters;
      } /* source color, spin */

  free(dst);

  /* Save if requested */

  save_wprop_from_field(saveflag, savefile, my_wqs, wp, "dummy recXML", 1);

  /* Convert to block Pauli basis */

  convert_wprop_to_block_pauli(bp, wp);

  free(wp);

  return tot_iters;
}

