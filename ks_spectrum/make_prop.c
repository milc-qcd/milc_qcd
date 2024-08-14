/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "ks_spectrum_includes.h"
#include "../include/dslash_ks_redefine.h"
#include <assert.h>

/************************************************************/

static int
load_ksprops(int num_prop, int startflag[], char startfile[][MAXFILENAME],
	     quark_source *my_ksqs[], ks_prop_field *source[],
	     ks_prop_field *ksprop[])
{
  int have_initial_guess = 0;
  for(int j = 0; j < num_prop; j++){
    if(startflag[j] != FRESH){
      int status = reload_ksprop_to_ksp_field(startflag[j], startfile[j],
					  my_ksqs[j], source[j], ksprop[j], 1);
      if(status != 0){
	node0_printf("Failed to reload propagator\n");
	terminate(1);
      } else {
	node0_printf("Restored propagator from %s\n",startfile[j]);
	have_initial_guess = 1;
      }
    }
  }
  return have_initial_guess;

}

/************************************************************/

static void
save_ksprops(int num_prop, int saveflag[], char savefile[][MAXFILENAME],
	     quark_source *my_ksqs[], ks_prop_field *source[],
	     ks_prop_field *ksprop[])
{
  /* Save solutions if requested */
  for(int j = 0; j < num_prop; j++){
    int status = save_ksprop_from_ksp_field( saveflag[j], savefile[j], "",
					     my_ksqs[j], source[j], ksprop[j], 1);
    if(status != 0){
      node0_printf("Failed to write propagator\n");
      terminate(1);
    }
    if(saveflag[j] != FORGET)
      node0_printf("Saved propagator to %s\n",savefile[j]);
  }
}

/************************************************************/

static void
copy_ksprops(int num_prop, ks_prop_field *ksprop[], ks_prop_field *source[])
{
  for(int j = 0; j < num_prop; j++){
    assert(ksprop[j]->nc == source[j]->nc);
    int nc = source[j]->nc;
    for(int color = 0; color < nc; color++)
      copy_v_field(ksprop[j]->v[color], source[j]->v[color]);
  }
}
  
/************************************************************/
/* Solve for the propagator (if requested) for all members of the set */

int solve_ksprop(enum set_type set_type, enum inv_type inv_type,
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

  int tot_iters = 0;
  Real mybdry_phase[4];
  int nc = source[0]->nc;   /* Must be the same for all in this set */
  int num_src = num_prop*nc;
  char myname[] = "solve_ksprop";

  /* Consistency checks */

  /* All members of the set must have matching inversion parameters,
     naik epsilons, and numbers of colors */

  for(int j = 1; j < num_prop; j++){
    if(my_qic[j].resid != my_qic[0].resid){
      if(set_type == MULTIMASS_SET){
	node0_printf("WARNING: %s: inversion error parameters do not match within the set\n", myname);
	node0_printf("WARNING: %s: will correct this in the refinement step.\n", myname);
      } else {
	node0_printf("ERROR: %s: found a nonmatching inversion error parameter in the set\n", myname);
	terminate(1);
      }
    }
    if(my_ksp[j].naik_term_epsilon_index != my_ksp[0].naik_term_epsilon_index){
      if(set_type == MULTIMASS_SET){
	node0_printf("WARNING: %s: Naik parameters do not match within the set\n", myname);
	node0_printf("WARNING: %s: will correct for this in the refinement step.\n", myname);
      } else {
	node0_printf("ERROR: %s: Naik epsilon  mismatch within set\n", myname);
	terminate(1);
      }
    }
    if(source[j]->nc != source[0]->nc){
      node0_printf("ERROR: %s: Source color mismatch within set\n", myname);
	terminate(1);
    }
  }

  /* All members of the set must have matching inv_types */
  for(int j = 0; j < num_prop; j++){
    if(my_qic[j].inv_type != inv_type){
      node0_printf("ERROR: %s: inversion type mismatch within set %d %d\n", myname,my_qic[j].inv_type,inv_type);
      terminate(1);
    }
  }

  switch(set_type){

  case SINGLES_SET:
    break;

  case MULTIMASS_SET:

    /* Sources and gauge links must be the same for all masses */
    for(int j = 1; j < num_prop; j++){
      if(source[j] != source[0]){
	node0_printf("ERROR: %s: source mismatch within set\n", myname);
	terminate(1);
      }
    }
    break;

  case MULTICOLORSOURCE_SET:
  case MULTISOURCE_SET:
    
    /* Masses must be the same for all sources */
    
    for(int j = 1; j < num_prop; j++){
      if(my_ksp[j].mass != my_ksp[0].mass){
	node0_printf("ERROR: %s: mass mismatch within set\n", myname);
	terminate(1);
      }
    }
    break;
    
  default:
    
    node0_printf("Can't handle set type %d\n", set_type);
    terminate(1);
  }
   
  /* Local copy of bdry_phase */
  for(int i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

  /* Load any requested propagators from files */
  int have_initial_guess =
    load_ksprops(num_prop, startflag, startfile, my_ksqs, source, ksprop);

  /* If we are reusing the input propagators and not recomputing, we are done */
  if(check == CHECK_NO){
    save_ksprops(num_prop, saveflag, savefile, my_ksqs, source, ksprop);
    return 0;
  }

  /* "CHECK_SOURCE_ONLY" means no inversion: copy source to destination */
  if(check == CHECK_SOURCE_ONLY){
    copy_ksprops(num_prop, ksprop, source);
    node0_printf("Copied %d sources to propagators\n", num_prop);
    return 0;
  }

  /* Construct fermion links for the solver */

#ifdef U1_FIELD
  /* Apply U(1) phases if we are using it */
  u1phase_on(charge, u1_A);
  invalidate_fermion_links(fn_links);
#endif
  
  restore_fermion_links_from_site(fn_links, my_qic[0].prec);
  imp_ferm_links_t **fn = get_fm_links(fn_links);
  
  /* Apply twisted boundary conditions and move KS phases, if
     requested */
  /* This operation applies the phase to the boundary FN links */
  int n_naiks = fermion_links_get_n_naiks(fn_links);
  for(int j = 0; j < n_naiks; j++){
    set_boundary_twist_fn(fn[j], mybdry_phase, r0);
    boundary_twist_fn(fn[j], ON);
  }
  
  /* Copy pointers for fermion links, based on Naik epsilon indices */
  imp_ferm_links_t **fn_multi = (imp_ferm_links_t **)
    malloc(sizeof(imp_ferm_links_t *)*num_prop);
  for(int j = 0; j < num_prop; j++)
      fn_multi[j] = fn[my_ksp[j].naik_term_epsilon_index];
  
  /* Apply the momentum twist to the sources. */
  
  /* This U(1) gauge transformation converts the boundary twist on
     the gauge field above into the desired volume twist. We do it
     this way to make our coding compatible with QOP, which does
     only a surface twist. */
  
  /* The time phase is special.  It is applied only on the
     boundary, so we don't gauge-transform it to the volume here.
     It is used only for switching between periodic and
     antiperiodic bc's */
  
  mybdry_phase[3] = 0;
  for(int j = 0; j < num_prop; j++){
    /* Because source[] is a list of pointers that could have
       duplicates, we need to make sure the rephasing is done
       only once on each source */
    int rephase_already_done = 0;
    for(int k = 0; k < j; k++){
      if(source[j] == source[k]){
	rephase_already_done = 1;
	break;
      }
      if(!rephase_already_done)
	for(int color = 0; color < nc; color++)
	  rephase_v_field(source[j]->v[color], mybdry_phase, r0, 1);
    }
  }
  mybdry_phase[3] = bdry_phase[3];
  
  /* Apply the momentum twist to the initial guess(es) if any */
  
  if(have_initial_guess){
    mybdry_phase[3] = 0; 
    for(int j = 0; j < num_prop; j++){
      if(startflag[j] != FRESH)
	for(int color = 0; color < nc; color++)
	  rephase_v_field(ksprop[j]->v[color], mybdry_phase, r0, 1);
    }
    mybdry_phase[3] = bdry_phase[3]; 
  }
  
  /* Make a list of source and destination pointers */

  su3_vector **src = (su3_vector **)malloc(num_src*sizeof(su3_vector *));
  su3_vector **dst = (su3_vector **)malloc(num_src*sizeof(su3_vector *));

  for(int j = 0; j < num_prop; j++)
    for(int color = 0; color < nc; color++){
      src[num_prop*color+j] = source[j]->v[color];
      dst[num_prop*color+j] = ksprop[j]->v[color];
    }

    /* When we start from a preloaded solution we use the
	 unpreconditioned mat_invert_cg_field algorithm, instead of
	 the preconditioned mat_invert_uml_field here to avoid
	 "reconstructing", and so overwriting the odd-site solution.
	 This would be a degradation if the propagator were
	 precomputed in double precision, and we were doing single
	 precision here. When we are computing the propagator from a
	 fresh start, we use the preconditioned algorithm. So we
	 always use CGTYPE if we are startin from an initial guess
      */

  int it = inv_type;
  int st = set_type;
  if(have_initial_guess){
    if(set_type != SINGLES_SET || my_qic[0].inv_type == CGTYPE)
      node0_printf("Because an initial guess is given, treating as SINGLES and CG\n");
    set_type = SINGLES_SET;
    my_qic[0].inv_type = CGTYPE;
  }
    
  switch(set_type){

  case(SINGLES_SET):

    for(int k = 0; k < num_src; k++){
      int color = k/num_prop;
      int j = k % num_prop;
      node0_printf("%s: color index = %d; mass = %f\n", myname,
		   j/num_prop, my_ksp[j].mass);
      mat_invert_field(src[k], dst[k], my_qic+j, my_ksp[j].mass, fn_multi[j]);
    }
    break;
    
  case(MULTIMASS_SET):
    
    for(int color = 0; color < nc; color++){
      node0_printf("%s: color index = %d; all masses\n", myname, color);
      mat_invert_multi(src[num_prop*color], &dst[num_prop*color], my_ksp,
		       num_prop, my_qic, fn_multi);
    }
    break;
    
  case(MULTISOURCE_SET):

    for(int color = 0; color < nc; color++){
      node0_printf("%s: color index = %d; mass = %f\n", myname, color, my_ksp[0].mass);
      mat_invert_block(&src[num_prop*color], &dst[num_prop*color], my_ksp[0].mass,
		       num_prop, my_qic, fn_multi[0]);
    }
    break;

  case(MULTICOLORSOURCE_SET):

    node0_printf("%s: all colors; mass = %f\n", myname, my_ksp[0].mass);
    mat_invert_block(src, dst, my_ksp[0].mass, num_src, my_qic, fn_multi[0]);
    break;

    default:
      node0_printf("Can't handle set type %d\n", set_type);
      terminate(1);

  } /* switch(set_type) */

  my_qic[0].inv_type = it;
  set_type = st;   /* In case we use it again */
  
  /* Transform solution and restore source */

  mybdry_phase[3] = 0;
  for(int j = 0; j < num_prop; j++){
    for(int color = 0; color < nc; color++)
      rephase_v_field(ksprop[j]->v[color], mybdry_phase, r0, -1);
    int rephase_already_done = 0;
    for(int k = 0; k < j; k++){
      if(source[j] == source[k]){
	rephase_already_done = 1;
	break;
      }
      if(!rephase_already_done)
	for(int color = 0; color < nc; color++)
	  rephase_v_field(source[j]->v[color], mybdry_phase, r0, -1);
    }
  }
  mybdry_phase[3] = bdry_phase[3];

  /* If we are reusing the input propagators and not recomputing, we are done */
  save_ksprops(num_prop, saveflag, savefile, my_ksqs, source, ksprop);

  /* Unapply twisted boundary conditions on the fermion links and
     restore conventional KS phases and antiperiodic BC, if
     changed. */
  for(int j = 0; j < n_naiks; j++)
    boundary_twist_fn(fn[j], OFF);

#ifdef U1_FIELD
  /* Unapply the U(1) field phases */
  u1phase_off();
  invalidate_fermion_links(fn_links);
#endif
  
  /* Clean up */
  free(dst);
  free(src);
  if(fn_multi != NULL)free(fn_multi);
  
  return tot_iters;
}

/************************************************************/
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
     
     So we save with all zero source fields */

  /* Set up an empty source */
  ks_prop_field *dummy_src = create_ksp_field(ksprop->nc);

  init_qs(&dummy_ksqs);
  ksqstmp = dummy_ksqs;   /* For ksprop_info.c */
  save_ksprop_from_ksp_field(saveflag, savefile, "", &dummy_ksqs, dummy_src, ksprop, 1);
  clear_qs(&dummy_ksqs); /* Free any allocations */

  destroy_ksp_field(dummy_src);
}

/************************************************************/
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
