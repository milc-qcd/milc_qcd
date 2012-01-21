/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "cl_inv_includes.h"
#include "../include/fermion_links.h"

// DEBUG
// Forces an unsophisticated static propagator -- NOT FOR MPP!
void static_prop_wv(wilson_vector *dst, wilson_vector *src, quark_source *my_qs){
  int i, jsrc;
  int t0 = my_qs->t0;
  site *s;
  

  FORALLSITES(i,s){
    /* src index at the dst spatial coordinate */
    jsrc = node_index(s->x,s->y,s->z,t0);
    dst[i] = src[jsrc];
  }
}


int get_wprop_to_wp_field(int startflag, char startfile[], 
			  int saveflag, char savefile[],
			  wilson_prop_field *wp,
			  quark_source *my_wqs,
			  quark_invert_control *my_qic,
			  dirac_clover_param *my_dcp,
			  Real bdry_phase[],
			  int r0[4],
			  int check)
{
  int i, color, spin, ksource;
  int cl_cg = CL_CG;
  int avs_iters = 0;
  int tot_iters = 0;
  int status;
  char *fileinfo;
  wilson_vector *dst;
  w_prop_file *fp_in, *fp_out; 
  char myname[] = "get_wprop_to_wp_field";
  Real mybdry_phase[4];
#ifdef IOTIME
  int io_timing = 1;
#else
  int io_timing = 0;
#endif

  node0_printf("%s: Generate Dirac propagator for kappa %g\n", 
	       myname, my_dcp->Kappa);

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

//  r0[0] = my_wqs->x0;
//  r0[1] = my_wqs->y0;
//  r0[2] = my_wqs->z0;
//  r0[3] = my_wqs->t0;

  dst = create_wv_field();

  /* For clover_info */
  wqstmp = *my_wqs;
  dcptmp = *my_dcp;

  /* Open files for Wilson propagators, if requested */
  fp_in  = r_open_wprop(startflag, startfile);
  fp_out = w_open_wprop(saveflag,  savefile, my_wqs->type);

  /* Provision for writing the source to a file */
  if(my_wqs->saveflag != FORGET){
    fileinfo = create_ws_XML(startfile, my_wqs);
    w_source_open_dirac(my_wqs, fileinfo);
    free(fileinfo);
  }
  
  /* Apply twist to boundary links of the gauge field in site structure */
  boundary_twist_site(mybdry_phase, r0, +1);

  /* Loop over source colors and spins */
  for(ksource = 0; ksource < my_wqs->nsource; ksource++){
    spin = convert_ksource_to_spin(ksource);
    color = convert_ksource_to_color(ksource);

    status = reload_wprop_sc_to_field(startflag, fp_in, my_wqs, 
				      spin, color, dst, 1);
    if(status != 0){
      printf("%s(%d): Error reloading propagator\n", myname, this_node);
      terminate(1);
    }
      
    /* (Re)construct propagator */
    
    if(startflag == FRESH)my_qic->start_flag = START_ZERO_GUESS;
    else                  my_qic->start_flag = START_NONZERO_GUESS;      
    
    /* Solve for the propagator if the starting guess is zero
       or we didn't say not to solve. */
    if(check != CHECK_NO || startflag == FRESH){

      wilson_vector *src;

      /* Make the source */
      src = create_wv_field();
      
      /* Create the source */
      if(wv_source_field(src, my_wqs)){
	printf("%s(%d): error getting source\n",myname, this_node);
	terminate(1);
      };
      
      /* Write the source, if requested */
      if(my_wqs->saveflag != FORGET){
	if(w_source_dirac( src, my_wqs ) != 0){
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
	rephase_wv_field(src, mybdry_phase, r0, 1);
	mybdry_phase[3] = bdry_phase[3]; 
	
	/* solve for dst */
	
	switch (cl_cg) {
	case BICG:
	  avs_iters = bicgilu_cl_field(src, dst, my_qic,(void *)my_dcp);
	  break;
	case HOP:
	  avs_iters = hopilu_cl_field(src, dst, my_qic,(void *)my_dcp);
	  break;
	case MR:
	  avs_iters = mrilu_cl_field(src, dst, my_qic,(void *)my_dcp);
	  break;
	case CG:
	  avs_iters = cgilu_cl_field(src, dst, my_qic,(void *)my_dcp);
	  break;
	default:
	  node0_printf("main(%d): Inverter choice %d not supported\n",
		       this_node,cl_cg);
	}

	report_status(my_qic);

	// DEBUG
	//static_prop_wv(dst, src, my_wqs);

      }

      destroy_wv_field(src);
    }
    
    /* Transform solution, completing the U(1) gauge transformation */
    mybdry_phase[3] = 0; 
    rephase_wv_field(dst, mybdry_phase, r0, -1);
    mybdry_phase[3] = bdry_phase[3]; 
    
    /* Copy solution vector to wilson_prop_field wp */
    copy_wp_from_wv( wp, dst, color, spin);
    
    /* save solution if requested */
    save_wprop_sc_from_field( saveflag, fp_out, my_wqs, spin, color, dst, 
			      "", io_timing);
    
    tot_iters += avs_iters;
  } /* ksource */
  
  /* Unapply twist to boundary links in gauge field */
  boundary_twist_site(mybdry_phase, r0, -1);
  
  /* clean up */
  if(my_wqs->saveflag != FORGET)
    w_source_close(my_wqs);

  clear_qs(my_wqs);
  r_close_wprop(startflag, fp_in);
  if(startflag != FRESH)
    node0_printf("Restored propagator from %s\n",startfile);

  w_close_wprop(saveflag,  fp_out);
  if(saveflag != FORGET)
    node0_printf("Saved propagator to %s\n",savefile);

  destroy_wv_field(dst);

  return tot_iters;
}

// DEBUG
// Forces an unsophisticated static propagator -- NOT MPP!
void static_prop_v(su3_vector *dst, su3_vector *src, quark_source *my_ksqs){
  int i, jsrc;
  int t0 = my_ksqs->t0;
  site *s;
  

  FORALLSITES(i,s){
    /* src index at the dst spatial coordinate */
    jsrc = node_index(s->x,s->y,s->z,t0);
    dst[i] = src[jsrc];
  }
}

int get_ksprop_to_wp_field(int startflag, char startfile[], 
			   int saveflag, char savefile[],
			   wilson_prop_field *wp,
			   quark_source *my_ksqs,
			   quark_invert_control *my_qic,
			   ks_param *my_ksp,
			   Real bdry_phase[],
			   int r0[4],
			   int check)
{
  int i, color;
  int avs_iters = 0;
  int tot_iters = 0;
  int status;
  int ks_source_r[4];
  int naik_epsilon_index = my_ksp->naik_term_epsilon_index;
  su3_vector *dst, *src;
  ks_prop_file *fp_in, *fp_out; 
  char *fileinfo;
  imp_ferm_links_t **fn = NULL;
  Real mybdry_phase[4];
#ifdef IOTIME
  int io_timing = 1;
#else
  int io_timing = 0;
#endif
  char myname[] = "get_ksprop_to_wp_field";

  node0_printf("%s: Generate naive propagator of mass %g\n",myname, my_ksp->mass);

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

  dst = create_v_field();

  /* For ksprop_info */
  ksqstmp = *my_ksqs;
  ksptmp  = *my_ksp;

//  r0[0] = my_ksqs->x0;
//  r0[1] = my_ksqs->y0;
//  r0[2] = my_ksqs->z0;
//  r0[3] = my_ksqs->t0;

  /* Phases must be in before constructing the fermion links */
  rephase( ON );
  invalidate_fermion_links(fn_links);

  /* Need FN links if we will call the inverter */
  if(check || startflag == FRESH){
    if(fn_links == NULL)
      fn_links = create_fermion_links_from_site(my_qic->prec, n_naiks, eps_naik);
    else
      restore_fermion_links_from_site(fn_links, my_qic->prec);

    fn = get_fm_links(fn_links);
    
    /* Apply twisted boundary conditions if requested */
    /* This operation applies the phase to the boundary FN links */
    set_boundary_twist_fn(fn[naik_epsilon_index], mybdry_phase, r0);
    boundary_twist_fn(fn[naik_epsilon_index], ON);
  }


  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);
  fp_out = w_open_ksprop(saveflag, savefile, my_ksqs->savetype);

  /* Provision for writing the source to a file */
  if(my_ksqs->saveflag != FORGET){
    fileinfo = create_ks_XML();
    w_source_open_ks(my_ksqs, fileinfo);
    free(fileinfo);
  }
  
  /* Check (or produce) the solution if requested */
  /* Loop over source colors */
  for(color=0;color<my_ksqs->ncolor;color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, my_ksqs, 
				      color, dst, 1);
    if(status != 0){
      printf("%s(%d): Error reloading propagator\n", myname, this_node);
      terminate(1);
    }
      
    /* (Re)construct propagator */
    
    if(startflag == FRESH) my_qic->start_flag = START_ZERO_GUESS;
    else                   my_qic->start_flag = START_NONZERO_GUESS;      
    
    /* Solve for the propagator if the starting guess is zero
       or we didn't say not to solve. */
    if(check != CHECK_NO || startflag == FRESH){

      /* Create the source common to this inversion */
      src = create_v_field();
      
      if(v_source_field(src, my_ksqs)){
	printf("%s(%d): error getting source\n",myname,this_node);
	terminate(1);
      };

	/* Write the source, if requested */
      if(my_ksqs->saveflag != FORGET){
	if(w_source_ks( src, my_ksqs ) != 0)
	  node0_printf("Error writing source\n");
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
	
	if(startflag == FRESH){
	  avs_iters += mat_invert_uml_field(src, dst, 
	    my_qic, my_ksp->mass, fn[naik_epsilon_index]);
	} else {
	  avs_iters += mat_invert_cg_field(src, dst,
 	    my_qic, my_ksp->mass, fn[naik_epsilon_index]);
	}

	// DEBUG
	//static_prop_v(dst, src, my_ksqs);
      }

      destroy_v_field(src);

      /* Transform solution, completing the U(1) gauge transformation */
      mybdry_phase[3] = 0; 
      rephase_v_field(dst, mybdry_phase, r0, -1);
      mybdry_phase[3] = bdry_phase[3]; 

    }
    
    /* save solution if requested */
    save_ksprop_c_from_field( saveflag, fp_out, my_ksqs, color, dst, "", io_timing);
    
    /* Convert KS prop to naive prop (su3_vector maps to
       spin_wilson_vector) for a given source color */
    ks_source_r[0] = my_ksqs->x0;  ks_source_r[1] = my_ksqs->y0;
    ks_source_r[2] = my_ksqs->z0;  ks_source_r[3] = my_ksqs->t0;
    
    convert_ksprop_to_wprop_swv(wp->swv[color], dst, ks_source_r, r0);
    
    tot_iters += avs_iters;
  } /* source color */

  /* Unapply twisted boundary conditions on the fermion links */
  if(check || startflag == FRESH)
      boundary_twist_fn(fn[naik_epsilon_index], OFF);


  /* Turn off staggered phases in the gauge links */
  rephase( OFF );
  invalidate_fermion_links(fn_links);

  /* clean up */
  if(my_ksqs->saveflag != FORGET)
    w_source_close(my_ksqs);

  clear_qs(my_ksqs);
  r_close_ksprop(startflag, fp_in);
  w_close_ksprop(saveflag,  fp_out);

  if(startflag != FRESH)
    node0_printf("Restored propagator from %s\n",startfile);

  if(saveflag != FORGET)
    node0_printf("Saved propagator to %s\n",savefile);

  destroy_v_field(dst);

  return tot_iters;
}

/* Sum (accumulate) a wilson vector (a += b) */
static void sum_wvec(wilson_vector *a, wilson_vector *b){
  int c,s;
  for(s = 0; s < 4; s++)
    for(c = 0; c < 3; c++)
      CSUM(a->d[s].c[c],b->d[s].c[c]);
}

/* Accumulate specific components of a spin wilson vector */
static void add_ksprop4_to_swv(spin_wilson_vector *swv_sum, 
			       spin_wilson_vector *swv, 
			       int spin_snk, int spin_src){
  int i;

  FORALLFIELDSITES(i){
    sum_wvec(&swv_sum[i].d[spin_src], &swv[i].d[spin_snk]);
  }
}

/* In this case we generate a naive propagator from a Dirac source */
/* It is assumed that the Dirac source has been prepared in the
   "staggered" basis.  See ext_src/make_ext_src.c */

int get_ksprop4_to_wp_field(int startflag, char startfile[], 
			    int saveflag, char savefile[],
			    wilson_prop_field *wp,
			    quark_source *my_qs,
			    quark_invert_control *my_qic,
			    ks_param *my_ksp,
			    Real bdry_phase[],
			    int r0[4],
			    int check)
{
  int i, color_src, spin_src, spin_snk, ksource;
  int avs_iters = 0;
  int tot_iters = 0;
  int ks_source_r[4] = {0,0,0,0};   /* Hypercube corners */
  int naik_epsilon_index = my_ksp->naik_term_epsilon_index;
  su3_vector *dst, *src;
  spin_wilson_vector *swv;
  w_prop_file *fp_out; 
  imp_ferm_links_t **fn = NULL;
  Real mybdry_phase[4];
#ifdef IOTIME
  int io_timing = 1;
#else
  int io_timing = 0;
#endif
  char myname[] = "get_ksprop4_to_wp_field";

  node0_printf("%s: Generate KS4 propagator of mass %g\n",myname, my_ksp->mass);

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

  dst = create_v_field();
  swv  = create_swv_field();
  clear_wp_field(wp);

  /* For ksprop_info */
  ksqstmp = *my_qs;
  ksptmp  = *my_ksp;

//  r0[0] = my_qs->x0;
//  r0[1] = my_qs->y0;
//  r0[2] = my_qs->z0;
//  r0[3] = my_qs->t0;

  /* Phases must be in before constructing the fermion links */
  rephase( ON );
  invalidate_fermion_links(fn_links);

  /* If we will call the inverter, we need the FN links */
  if(check == CHECK_YES || startflag == FRESH){
    if(fn_links == NULL)
      fn_links = create_fermion_links_from_site(my_qic->prec, n_naiks, eps_naik);
    else
      restore_fermion_links_from_site(fn_links, my_qic->prec);
  
    fn = get_fm_links(fn_links);
    /* Apply twisted boundary conditions if requested */
    /* This operation applies the phase to the boundary FN links */
    set_boundary_twist_fn(fn[naik_epsilon_index], mybdry_phase, r0);
    boundary_twist_fn(fn[naik_epsilon_index], ON);
  }
	
  /* Open file for Wilson propagator, if requested */
  fp_out = w_open_wprop(saveflag, savefile, my_qs->type);

  /* Check (or produce) the solution if requested */
  /* Loop over source colors and spins */
  for(ksource = 0; ksource < my_qs->nsource; ksource++){
    /* Complete the source structure */
    spin_src = convert_ksource_to_spin(ksource);
    color_src = convert_ksource_to_color(ksource);
    for(spin_snk=0;spin_snk<4;spin_snk++){
      my_qs->spin_snk = spin_snk;

      /* (Re)construct propagator */
      
      if(startflag == FRESH) my_qic->start_flag = START_ZERO_GUESS;
      else                   my_qic->start_flag = START_NONZERO_GUESS;      
    
      /* Solve for the propagator if the starting guess is zero
	 or we said to solve. */
      if(check == CHECK_YES || startflag == FRESH){

	/* Create the source common to this inversion */
	src = create_v_field();
      
	if(v_source_field(src, my_qs)){
	  printf("%s(%d): error getting source\n",myname,this_node);
	  terminate(1);
	};


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
	
	if(startflag == FRESH){
	  avs_iters += mat_invert_uml_field(src, dst, 
	    my_qic, my_ksp->mass, fn[naik_epsilon_index]);
	} else {
	  avs_iters += mat_invert_cg_field(src, dst,
	   my_qic, my_ksp->mass, fn[naik_epsilon_index]);
	}

	// DEBUG
	//static_prop_v(dst, src, my_qs);
	destroy_v_field(src);
	
	/* Transform solution, completing the U(1) gauge transformation */
	mybdry_phase[3] = 0; 
	rephase_v_field(dst, mybdry_phase, r0, -1);
	mybdry_phase[3] = bdry_phase[3]; 
	
      }
    
      /* Convert KS prop to naive prop (su3_vector maps to
	 spin_wilson_vector) for a given source color */
      
      convert_ksprop_to_wprop_swv(swv, dst, ks_source_r, r0);

      /* Add a component of the naive propagator to the result */
      node0_printf("Adding spin_snk %d spin_src %d for color_src %d ksource %d\n",
		   spin_snk, spin_src, color_src, ksource);
      add_ksprop4_to_swv(wp->swv[color_src], swv, spin_snk, spin_src);

      tot_iters += avs_iters;
    } /* spin_snk */

    /* save solution if requested */
    if(saveflag != FORGET){
      wilson_vector *wv = create_wv_field();
      copy_wv_from_swv(wv, wp->swv[color_src], spin_src);
      save_wprop_sc_from_field( saveflag, fp_out, my_qs, 
				spin_src, color_src, wv, "", io_timing);
      destroy_wv_field(wv);
    }
  } /* ksource -> color_src, spin_src */

  /* Unapply twisted boundary conditions on the fermion links */
  if(check == CHECK_YES || startflag == FRESH)
    boundary_twist_fn(fn[naik_epsilon_index], OFF);

  /* Turn off staggered phases in the gauge links */
  rephase( OFF );
  invalidate_fermion_links(fn_links);

  /* clean up */
  clear_qs(my_qs);
  
  w_close_wprop(saveflag,  fp_out);
  if(saveflag != FORGET)
    node0_printf("Saved propagator to %s\n",savefile);

  destroy_swv_field(swv);
  destroy_v_field(dst);
  
  return tot_iters;
}

/* Dump wilson propagator field to file */

void dump_wprop_from_wp_field(int saveflag, char savefile[], 
			      wilson_prop_field *wp){
  quark_source dummy_wqs;
#ifdef IOTIME
  int io_timing = 1;
#else
  int io_timing = 0;
#endif
  
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

  init_qs(&dummy_wqs);
  wqstmp = dummy_wqs;   /* For clover_info.c */
  save_wprop_from_wp_field(saveflag, savefile, &dummy_wqs, wp, "", io_timing);
  clear_qs(&dummy_wqs); /* Free any allocations */
}

/* Create a wilson_prop_field and restore it from a dump file */

void reread_wprop_to_wp_field(int saveflag, char savefile[], wilson_prop_field *wp){
  char myname[] = "reread_wprop_to_wp_field";
  quark_source dummy_wqs;
  int rereadflag = convert_outflag_to_inflag_wprop(saveflag);

  if(rereadflag == FRESH){
    printf("%s(%d) Can't reread file %s when saveflag is %d\n",
	   myname, this_node, savefile, saveflag);
    terminate(1);
  }

  rebuild_wp_field(wp);

  init_qs(&dummy_wqs);
  reload_wprop_to_wp_field(rereadflag, savefile, &dummy_wqs, wp, 1);
  clear_qs(&dummy_wqs); /* Free any allocations */
}
