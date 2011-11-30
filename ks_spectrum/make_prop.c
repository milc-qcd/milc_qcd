/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "ks_spectrum_includes.h"
#include "../include/dslash_ks_redefine.h"

// void read_ksprop_to_ksp_field(int startflag, char startfile[], 
// 			      quark_source *my_ksqs, ks_prop_field *ksprop)
// {
//   int color;
//   int status;
//   ks_prop_file *fp_in; 
// 
//   /* Open files for KS propagators, if requested */
//   /* We need to modify our code so we can get the number of colors here,
//      if it is not already 3. */
//   fp_in  = r_open_ksprop(startflag, startfile);
// 
//   /* Loop over source colors */
//   for(color=0;color<ksprop->nc;color++){
//     
//     /* Read color vector (and source as appropriate) from file */
//     status = reload_ksprop_c_to_field(startflag, fp_in, my_ksqs, 
// 				      color, ksprop->v[color], 1);
//   } /* source color */
//   
//   /* clean up */
//   r_close_ksprop(startflag, fp_in);
// 
//   if(startflag != FRESH)
//     node0_printf("Restored propagator from %s\n",startfile);
// 
// }

/* Solve for the propagator (if requested) for all members of the set */

int solve_ksprop(int num_prop, int startflag[], char startfile[][MAXFILENAME],
		 int saveflag[], char savefile[][MAXFILENAME],
		 ks_prop_field *ksprop[],
		 quark_source *my_ksqs,
		 quark_invert_control my_qic[],
		 ks_param my_ksp[],
		 Real bdry_phase[],
		 int r0[4],
		 int check)
{

  int color;
  int i,j;
  int status = 0;
  char *fileinfo;
  int tot_iters = 0;
  ks_prop_file *fp_in[MAX_PROP], *fp_out[MAX_PROP]; 
  su3_vector **dst, *src;
  imp_ferm_links_t **fn = NULL;
  Real mybdry_phase[4];
  imp_ferm_links_t **fn_multi = NULL;
  int n_naiks = fermion_links_get_n_naiks(fn_links);
  char myname[] = "solve_ksprop";

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

  /* Offset for setting the staggered phases in FN links.
     They can be different for each propagator
     NO, they must be the same! */

//  r0[0] = my_ksqs->x0;
//  r0[1] = my_ksqs->y0;
//  r0[2] = my_ksqs->z0;
//  r0[3] = my_ksqs->t0;
  
  ksqstmp = *my_ksqs;     /* For ksprop_info. Source is common to the set */
  ksptmp  = my_ksp[0];      /* For action parameters */

  /* Open files for KS propagators, if requested */

  for(j = 0; j < num_prop; j++){
    fp_in[j]  = r_open_ksprop(startflag[j], startfile[j]);
    fp_out[j] = w_open_ksprop(saveflag[j],  savefile[j], my_ksqs->type);
  }

  /* Provision for writing the source to a file */
  
  if(my_ksqs->saveflag != FORGET){
    fileinfo = create_ks_XML();
    w_source_open_ks(my_ksqs, fileinfo);
    free(fileinfo);
  }
  
  /* Construct fermion links if we will need them */

  if(check != CHECK_NO || startflag[0] == FRESH ){
    
    restore_fermion_links_from_site(fn_links, my_qic[0].prec);
    fn = get_fm_links(fn_links);

    /* Apply twisted boundary conditions and move KS phases and
       time boundary to new origin r0, if requested */
    /* This operation applies the phase to the boundary FN links */
    for(j = 0; j < n_naiks; j++){
      set_boundary_twist_fn(fn[j], mybdry_phase, r0);
      boundary_twist_fn(fn[j], ON);
    }

    /* Copy pointers for fermion links, based on Naik epsilon indices */
    fn_multi = (imp_ferm_links_t **)
      malloc(sizeof(imp_ferm_links_t *)*num_prop);
    for(j = 0; j < num_prop; j++)
      fn_multi[j] = fn[my_ksp[j].naik_term_epsilon_index];
    
  }

  /* Check (or produce) the solution if requested */
  for(color = 0; color < my_ksqs->ncolor; color++){
    
    node0_printf("%s: color = %d\n",myname, color);

    /* Load the propagators */

    for(j = 0; j < num_prop; j++){
      status = reload_ksprop_c_to_field(startflag[j], fp_in[j], my_ksqs, 
					color, ksprop[j]->v[color], 1);
    }
      
    /* Solve for the propagator if the starting guess is zero
       or we didn't say not to solve. */
    if(check != CHECK_NO || startflag[0] == FRESH){
      /* Create the source common to this inversion */
      src = create_v_field();

      if(v_source_field(src, my_ksqs)){
	printf("%s(%d): error getting source\n",myname,this_node);
	terminate(1);
      };
    
      /* Write the source, if requested */
      if(my_ksqs->saveflag != FORGET){
	if(w_source_ks( src, my_ksqs ) != 0){
	  node0_printf("Error writing source\n");
	}
      }
      
      if(check != CHECK_SOURCE_ONLY){

	/* Apply the momentum twist to the source.  This U(1) gauge
	   transformation converts the boundary twist on the gauge field
	   above into the desired volume twist. We do it this way to
	   make our coding compatible with QOP, which does only a
	   surface twist. */
	
	/* The time phase is special.  It is applied only on the
	   boundary, so we don't gauge-transform it to the volume here.
	   It is used only for switching between periodic and
	   antiperiodic bc's */
	
	mybdry_phase[3] = 0; 
	rephase_v_field(src, mybdry_phase, r0, 1);
	mybdry_phase[3] = bdry_phase[3]; 
	
	/* Make a list of num_prop solution pointers for the current color */
	dst = (su3_vector **)malloc(num_prop*sizeof(su3_vector *));
	for(j = 0; j < num_prop; j++){
	  dst[j] = ksprop[j]->v[color];
	}
	
	if(num_prop == 1){
	  
	  /* Single mass inversion */
	  
	  if(startflag[0] == FRESH) my_qic[0].start_flag = START_ZERO_GUESS;
	  else                      my_qic[0].start_flag = START_NONZERO_GUESS;
	  
	  /* When we start from a preloaded solution we use the less
	     optimized mat_invert_cg_field algorithm, instead of
	     mat_invert_uml_field here to avoid "reconstructing", and so
	     overwriting the odd-site solution.  This would be a
	     degradation if the propagator was precomputed in double
	     precision, and we were doing single precision here. When we
	     are computing the propagator from a fresh start, we use
	     the optimized mat_invert_uml_field algorithm. */
	  
	  if(startflag[0] == FRESH){
	    tot_iters += mat_invert_uml_field(src, dst[0], 
	      my_qic+0, my_ksp[0].mass, fn_multi[0]);
	  } else {
	    tot_iters += mat_invert_cg_field(src, dst[0],
	     my_qic+0, my_ksp[0].mass, fn_multi[0]);
	  }
	} else {
	  
	  /* Multimass inversion */
	  tot_iters = 
	    mat_invert_multi(src, dst, my_ksp, num_prop, my_qic, fn_multi);
	}
	
	
	/* Transform solution, completing the U(1) gauge transformation */
	mybdry_phase[3] = 0; 
	for(j = 0; j < num_prop; j++){
	  rephase_v_field(dst[j], mybdry_phase, r0, -1);
	}
	mybdry_phase[3] = bdry_phase[3]; 
	
	/* save solutions if requested */
	for(i = 0; i < num_prop; i++){
	  save_ksprop_c_from_field( saveflag[i], fp_out[i], my_ksqs, 
				    color, dst[i], "", 1);
	}
	
#ifdef DEBUG_NAIVE
	rephase( OFF );
	for(i = 0; i < num_prop; i++)
	  {
	    spin_wilson_vector *swv = create_swv_field();
	    wilson_vector *wv  = create_wv_field();
	    wilson_vector *wvsrc = create_wv_field();
	    su3_vector *v = create_v_field();
	    int ks_source_r[4] = {0,0,0,0};   /* Hypercube corners */
	    int spin;
	    
	    convert_ksprop_to_wprop_swv(swv, dst[i], ks_source_r, r0);
	    for(spin = 0; spin < 4; spin++){
	      copy_wv_from_swv(wv, swv, spin);
	      clear_wv_field(wvsrc);
	      copy_wv_from_v(wvsrc, src, spin);
	      check_naive(wv, wvsrc, my_ksp[i].mass, 1e-5);
	    }
	    
	    destroy_v_field(v);
	    destroy_wv_field(wv);
	    destroy_wv_field(wvsrc);
	    destroy_swv_field(swv);
	  }
	rephase( ON );
#endif
	
	/* Clean up */
	free(dst);
      } /* if(check != CHECK_SOURCE_ONLY) */

      destroy_v_field(src);

    } /* if(check != CHECK_NO || startflag[0] == FRESH)} */

  } /* color */

  
  /* Unapply twisted boundary conditions on the fermion links and
     restore conventional KS phases and antiperiodic BC, if
     changed. */
  if(check != CHECK_NO || startflag[0] == FRESH ){
    for(j = 0; j < n_naiks; j++)
      boundary_twist_fn(fn[j], OFF);
  }
    
  for(i = 0; i < num_prop; i++){
    w_close_ksprop(saveflag[i],  fp_out[i]);
    if(saveflag[i] != FORGET)
      node0_printf("Saved propagator to %s\n",savefile[i]);
  }

  return tot_iters;

}
/* Dump wilson propagator field to file */

void dump_ksprop_from_ksp_field(int saveflag, char savefile[], 
				ks_prop_field *ksprop){
  quark_source dummy_ksqs;
  
  /* When we dump a propagator, we don't keep the source information.
     Normally we want the source information for checking consistency
     with the Dirac operator we are using, since if we know the
     source, we can run the propagator through the solver to check it.
     Here the propagator we are dumping is usually one that has
     already had sink operators applied to it, so it won't satisfy the
     Dirac equation, anyway.
     
     So we take a default source type "UNKOWN".  A minimal source
     record is still written to the file, since we don't have any
     propagator file formats without sources records.  The minimal
     source record is a null complex field on time slice zero. */

  /* For clover_info.c */

  init_qs(&dummy_ksqs);
  ksqstmp = dummy_ksqs;   /* For ksprop_info.c */
  save_ksprop_from_ksp_field(saveflag, savefile, "", &dummy_ksqs, ksprop, 1);
  clear_qs(&dummy_ksqs); /* Free any allocations */
}

/* Create a ks_prop_field and restore it from a dump file */

ks_prop_field *reread_ksprop_to_ksp_field(int saveflag, char savefile[], 
					  int nc){
  char myname[] = "reread_ksprop_to_ksp_field";
  quark_source dummy_ksqs;
  ks_prop_field *ksprop;
  int rereadflag = convert_outflag_to_inflag_ksprop(saveflag);

  if(rereadflag == FRESH){
    printf("%s(%d) Can't reread file %s when saveflag is %d\n",
	   myname, this_node, savefile, saveflag);
    terminate(1);
  }

  ksprop = create_ksp_field(nc);

  init_qs(&dummy_ksqs);
  reload_ksprop_to_ksp_field(rereadflag, savefile, &dummy_ksqs, ksprop, 1);
  clear_qs(&dummy_ksqs); /* Free any allocations */
  return ksprop;
}
