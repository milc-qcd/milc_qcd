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

int solve_ksprop(int set_type,
		 int num_prop, int startflag[], char startfile[][MAXFILENAME],
		 int saveflag[], char savefile[][MAXFILENAME],
		 ks_prop_field *ksprop[],
		 ks_prop_field *source[],
                 quark_source *my_ksqs[],
		 quark_invert_control my_qic[],
		 ks_param my_ksp[],
		 Real charge,
		 Real bdry_phase[],
		 int r0[4],
		 int check)
{

  int color;
  int i,j;
  int status = 0;
  char *fileinfo;
  int tot_iters = 0;
  su3_vector **dst;
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

  ksqstmp = *my_ksqs[0];     /* For ksprop_info. Source is common to the set */
  ksptmp  = my_ksp[0];      /* For action parameters */

  /* Construct fermion links if we will need them */

  if(check != CHECK_NO || startflag[0] == FRESH ){

#ifdef U1_FIELD
    /* Apply U(1) phases if we are using it */
    u1phase_on(charge, u1_A);
    invalidate_fermion_links(fn_links);
#endif

    restore_fermion_links_from_site(fn_links, my_qic[0].prec);
    fn = get_fm_links(fn_links);

    /* Apply twisted boundary conditions and move KS phases, if
       requested */
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
    
  } /* check != CHECK_NO || startflag[0] == FRESH */

  /* Check (or produce) the solution if requested */
  
  /* Load the propagators */
  for(j = 0; j < num_prop; j++){
    if(startflag[j] != FRESH){
      status = reload_ksprop_to_ksp_field(startflag[j], startfile[j],
					  my_ksqs[j], source[j], ksprop[j], 1);
      if(status != 0){
	node0_printf("Failed to reload propagator\n");
	terminate(1);
      } else {
	node0_printf("Restored propagator from %s\n",startfile[j]);
      }
    }
  }

  /* Loop over source colors.  They should be the same for all sources. */
  for(color = 0; color < source[0]->nc; color++){
    
    node0_printf("%s: color = %d\n",myname, color);

    /* List pointers to sources for this color */
    su3_vector **src = (su3_vector **)malloc(num_prop*sizeof(su3_vector *));
    for(j = 0; j < num_prop; j++) src[j] = source[j]->v[color];

    /* Solve for the propagator if the starting guess is zero
       or we didn't say not to solve. */
    if(check != CHECK_NO || startflag[0] == FRESH){

      /* List pointers to solutions for the current color */
      dst = (su3_vector **)malloc(num_prop*sizeof(su3_vector *));
      for(j = 0; j < num_prop; j++) dst[j] = ksprop[j]->v[color];
	
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
	for(j = 0; j < num_prop; j++)
	  rephase_v_field(src[j], mybdry_phase, r0, 1);
	mybdry_phase[3] = bdry_phase[3]; 
	
	if(startflag[0] != FRESH){

	  /* Apply the momentum twist to the initial guess */
	  mybdry_phase[3] = 0; 
	  for(j = 0; j < num_prop; j++){
	    rephase_v_field(dst[j], mybdry_phase, r0, 1);
	  }
	  mybdry_phase[3] = bdry_phase[3]; 
	} /* startflag[0] != FRESH */

	if(num_prop == 1){
	  
	  /* Single mass inversion */
	  
	  /* When we start from a preloaded solution we use the less
	     optimized mat_invert_cg_field algorithm, instead of
	     mat_invert_uml_field here to avoid "reconstructing", and so
	     overwriting the odd-site solution.  This would be a
	     degradation if the propagator was precomputed in double
	     precision, and we were doing single precision here. When we
	     are computing the propagator from a fresh start, we use
	     the optimized mat_invert_uml_field algorithm. */
	  
	  if(startflag[0] == FRESH){
	    mat_invert_uml_field(src[0], dst[0], my_qic+0, my_ksp[0].mass, 
				 fn_multi[0]);
	  } else {
	    mat_invert_cg_field(src[0], dst[0], my_qic+0, my_ksp[0].mass, 
				fn_multi[0]);
	  }
	} else {
	  /* If we have restored any propagator, we use the single-mass inverter */
	  /* In most use cases they are either all restored, or all fresh */

	  if(startflag[0] != FRESH){
	    for(j = 0; j < num_prop; j++){
	      if(set_type == MULTIMASS_SET)
		/* Multimass inversion */
		mat_invert_cg_field(src[0], dst[j], my_qic+j, my_ksp[j].mass, 
				    fn_multi[j]);
	      else
		/* Multisource inversion */
		mat_invert_cg_field(src[j], dst[j], my_qic+j, my_ksp[0].mass, 
				    fn_multi[j]);
	    }
	  } else {

	    if(set_type == MULTIMASS_SET)
	      /* Multimass inversion */
	      mat_invert_multi(src[0], dst, my_ksp, num_prop, my_qic, fn_multi);
	    else {
	      /* Multisource inversion */
	      int num_src = num_prop;  /* Should change to num_prop * ncolors */
	      mat_invert_block_uml(src, dst, my_ksp[0].mass, num_src, my_qic, fn_multi[0]);
	    }
	  }
	}
	
	/* Transform solutions, completing the U(1) gauge transformation */
	mybdry_phase[3] = 0; 
	for(j = 0; j < num_prop; j++){
	  rephase_v_field(dst[j], mybdry_phase, r0, -1);
	}
	mybdry_phase[3] = bdry_phase[3]; 
	
      } else {
	
	/* Copy source to solution(s) so we can use it there */
	for(j = 0; j < num_prop; j++)
	  copy_v_field(dst[j], src[j]);
      }  /* if(check != CHECK_SOURCE_ONLY) */
      
      /* Clean up */
      free(dst);
      free(src);
      
    } /* if(check != CHECK_NO || startflag[0] == FRESH)} */

  } /* color */

  /* save solutions if requested */
  for(j = 0; j < num_prop; j++){
    status = save_ksprop_from_ksp_field( saveflag[j], savefile[j], "",
					 my_ksqs[j], source[j], ksprop[j], 1);
    if(status != 0){
      node0_printf("Failed to write propagator\n");
      terminate(1);
    }
    if(saveflag[j] != FORGET)
      node0_printf("Saved propagator to %s\n",savefile[j]);
  }

  if(check != CHECK_NO || startflag[0] == FRESH ){

  /* Unapply twisted boundary conditions on the fermion links and
     restore conventional KS phases and antiperiodic BC, if
     changed. */
    for(j = 0; j < n_naiks; j++)
      boundary_twist_fn(fn[j], OFF);
  }
    
#ifdef U1_FIELD
  /* Unapply the U(1) field phases */
  u1phase_off();
  invalidate_fermion_links(fn_links);
#endif

  if(fn_multi != NULL)free(fn_multi);

  return tot_iters;

}
/* Dump KS propagator field to file */

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
  save_ksprop_from_ksp_field(saveflag, savefile, "", &dummy_ksqs, NULL, ksprop, 1);
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
    printf("%s(%d) Can't reread file %s when saveflag is FRESH\n",
	   myname, this_node, savefile);
    terminate(1);
  }

  ksprop = create_ksp_field(nc);

  init_qs(&dummy_ksqs);
  reload_ksprop_to_ksp_field(rereadflag, savefile, &dummy_ksqs, NULL, ksprop, 1);
  clear_qs(&dummy_ksqs); /* Free any allocations */
  return ksprop;
}
