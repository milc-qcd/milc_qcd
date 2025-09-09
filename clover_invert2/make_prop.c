/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "cl_inv_includes.h"
#include "../include/fermion_links.h"
#include <string.h>

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


int get_wprop_to_wp_field(int prop_type, int startflag, char startfile[], 
			  int saveflag, char savefile[],
			  wilson_prop_field *wp,
			  wilson_prop_field *source,
			  quark_source *my_wqs,
			  quark_invert_control *my_qic,
			  void *my_dcp,
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
  char myname[] = "get_wprop_to_wp_field";
  Real mybdry_phase[4];
#ifdef IOTIME
  int io_timing = 1;
#else
  int io_timing = 0;
#endif

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

  /* Apply twist to boundary links of the gauge field in site structure */
  boundary_twist_site(mybdry_phase, r0, +1);

  /* Load the propagator if needed */
  status = reload_wprop_to_wp_field(startflag, startfile, my_wqs,
				    source, wp, 1);

  if(status != 0){
    printf("%s(%d): Error reloading propagator\n", myname, this_node);
    terminate(1);
  }

  /* Loop over source colors and spins */
  for(ksource = 0; ksource < my_wqs->nsource; ksource++){
    spin = convert_ksource_to_spin(ksource);
    color = convert_ksource_to_color(ksource);

    node0_printf("%s: spin = %d color = %d\n",myname, spin, color);

    /* Source for this color */
    wilson_vector *src = create_wv_field();
    extract_wv_from_swv(src, source->swv[color], spin);
      
    /* Solve for the propagator if the starting guess is zero
       or we didn't say not to solve. */
    if(check != CHECK_NO || startflag == FRESH){

      wilson_vector *dst = create_wv_field();
      if(startflag != FRESH)
	extract_wv_from_swv( dst, wp->swv[color], spin);

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
	
	if(startflag != FRESH){

	  /* Apply the momentum twist to the initial guess */
	  mybdry_phase[3] = 0; 
	  rephase_wv_field(dst, mybdry_phase, r0, 1);
	  mybdry_phase[3] = bdry_phase[3]; 
	} /* startflag[0] != FRESH */

	/* solve for dst */

	if(prop_type == CLOVER_TYPE){

	  switch (cl_cg) {
	  case BICG:
	    avs_iters = bicgilu_cl_field(src, dst, my_qic, my_dcp);
	    break;
	  case HOP:
	    avs_iters = hopilu_cl_field(src, dst, my_qic, my_dcp);
	    break;
	  case MR:
	    avs_iters = mrilu_cl_field(src, dst, my_qic, my_dcp);
	    break;
	  case CG:
	    avs_iters = cgilu_cl_field(src, dst, my_qic, my_dcp);
	    break;
	  default:
	    node0_printf("%s(%d): Inverter choice %d not supported\n",
			 myname, this_node,cl_cg);
	  }

	} else if(prop_type == IFLA_TYPE){

#ifdef HAVE_QOP
	  switch (cl_cg) {
	  case BICG:
	    avs_iters = bicgilu_cl_field_ifla(src, dst, my_qic, my_dcp);
	    break;
	  default:
	    node0_printf("%s(%d): Inverter choice %d not supported\n",
			 myname, this_node,cl_cg);
	  }
#else
	  node0_printf("%s: QOP compilation required for IFLA\n", myname);
	  terminate(1);
#endif
	} else {
	  node0_printf("%s: Unsupported propagator type\n", myname);
	  terminate(1);
	}

	report_status(my_qic);
	
	/* Transform solution, completing the U(1) gauge transformation */
	mybdry_phase[3] = 0; 
	rephase_wv_field(dst, mybdry_phase, r0, -1);
	mybdry_phase[3] = bdry_phase[3];

      } else {

	/* Copy source to solution(s) so we can use it there */
	copy_wv_field(dst, src);
      }  /* check != CHECK_SOURCE_ONLY */

      /* Copy solution to wp */
      insert_swv_from_wv(wp->swv[color], spin, dst);
      
      /* Clean up */
      destroy_wv_field(dst);

    } /* if(check != CHECK_NO || startflag[0] == FRESH)} */
    
    destroy_wv_field(src);
      
  } /* ksource */
    
    
  /* save solution if requested */
  save_wprop_from_wp_field( saveflag, savefile,"", my_wqs,
			    source, wp, io_timing);
    
  tot_iters += avs_iters;
  
  /* Unapply twist to boundary links in gauge field */
  boundary_twist_site(mybdry_phase, r0, -1);
  
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

/* In this case we generate a naive propagator from a KS vector
   source */

int get_ksprop_to_wp_field(int startflag, char startfile[], 
			   int saveflag, char savefile[],
			   wilson_prop_field *wp,
			   ks_prop_field *source,
			   quark_source *my_ksqs,
			   quark_invert_control *my_qic,
			   ks_param *my_ksp,
			   Real bdry_phase[],
			   int r0[4],
			   int check,
			   int is_ks0_type)
{
  int i, color;
  int avs_iters = 0;
  int tot_iters = 0;
  int status;
  char *fileinfo;
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

  /* For ksprop_info */
  ksqstmp = *my_ksqs;
  ksptmp  = *my_ksp;

  /* Need FN links if we will call the inverter */
  int naik_epsilon_index;
  imp_ferm_links_t *fn_naik = NULL;

  if(check == CHECK_YES || startflag == FRESH){

    rephase( ON );
    invalidate_fermion_links(fn_links);

    if(fn_links == NULL)
      fn_links = create_fermion_links_from_site(my_qic->prec, n_naiks, eps_naik);
    else
      restore_fermion_links_from_site(fn_links, my_qic->prec);

    naik_epsilon_index = my_ksp->naik_term_epsilon_index;
    fn_naik = get_fm_links(fn_links, naik_epsilon_index);
    
    /* Apply twisted boundary conditions and move KS phases, if
       requested */
    /* This operation applies the phase to the boundary FN links */
    set_boundary_twist_fn(fn_naik, mybdry_phase, r0);
    boundary_twist_fn(fn_naik, ON);
  } /* check != CHECK_NO || startflag[0] == FRESH */

  /* Check (or produce) the solution if requested */

  /* Load the propagator if needed */
  ks_prop_field *ksprop = NULL;
  if(startflag != FRESH){
    ksprop = create_ksp_field(my_ksqs->ncolor);
    status = reload_ksprop_to_ksp_field(startflag, startfile,
					my_ksqs, source, ksprop, io_timing);
    if(status != 0){
      node0_printf("Failed to reload propagator\n");
      terminate(1);
    } else {
      node0_printf("Restored propagator from %s\n",startfile);
    }
  }

  /* Loop over source colors */
  for(color=0;color<my_ksqs->ncolor;color++){
    
    node0_printf("%s: color = %d\n",myname, color);

    /* Source for this color */
    su3_vector *src = source->v[color];

    /* Solve for the propagator if the starting guess is zero
       or we didn't say not to solve. */
    if(check != CHECK_NO || startflag == FRESH){

      su3_vector *dst;
      if(ksprop == NULL)
	dst = create_v_field();
      else
	dst = ksprop->v[color];

      
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
	
	if(startflag != FRESH){

	  /* Apply the momentum twist to the initial guess */
	  mybdry_phase[3] = 0; 
	  rephase_v_field(dst, mybdry_phase, r0, 1);
	  mybdry_phase[3] = bdry_phase[3]; 
	} /* startflag != FRESH */

	/* In most use cases we will be reading a precomputed staggered
	   propagator, so we use the less optimized mat_invert_field
	   algorithm, instead of mat_invert_uml_field here to avoid
	   "reconstructing", and so overwriting the odd-site solution.
	   This would be a degradation if the propagator was precomputed
	   in double precision, and we were doing single precision
	   here. If we start using this code to compute the propagator
	   from a fresh start, we could make the selection of the
	   inverter algorithm dependent on the propagator start flag.
	   FRESH -> uml and otherwise cg. */
	
	enum inv_type it = my_qic[0].inv_type;
	if(startflag != FRESH)my_qic[0].inv_type = CGTYPE;
	avs_iters += mat_invert_field(src, dst, my_qic, my_ksp->mass, fn_naik);
	my_qic[0].inv_type = it;

	/* Transform solution, completing the U(1) gauge transformation */
	mybdry_phase[3] = 0; 
	rephase_v_field(dst, mybdry_phase, r0, -1);
	mybdry_phase[3] = bdry_phase[3]; 

      } else {
	
	/* Copy source to solution(s) so we can use it there */
	copy_v_field(dst, src);
      } /* check != CHECK_SOURCE_ONLY */

      /* Convert KS prop to naive prop (su3_vector maps to
	 spin_wilson_vector) for a given source color,
	 The conversion depends on the source and sink locations.
	 
	 Since the source location is buried by the inversion,
	 the KS source must have support on only one site in
	 the source hypercube, and that source location
	 should be in the same hypercube location as
	 the my_ksqs coordinate, x0, y0, z0, t0.
	 
	 SHOULD CHECK SOURCE TO MAKE SURE IT SATISFIES THESE CRITERIA
	 
	 For the KS0_TYPE, we promise that the source is at the corner
	 of the hypercube. */
      
      int ks_source_r[4] = {0,0,0,0};
      if(!is_ks0_type){
	ks_source_r[0] = my_ksqs->x0;  ks_source_r[1] = my_ksqs->y0;
	ks_source_r[2] = my_ksqs->z0;  ks_source_r[3] = my_ksqs->t0;
      }
      
      convert_ksprop_to_wprop_swv(wp->swv[color], dst, ks_source_r, r0);

      if(ksprop == NULL)
	destroy_v_field(dst);
    } /* check != CHECK_NO || startflag == FRESH */
    
    tot_iters += avs_iters;

  } /* color */

  if(startflag != FRESH)
    destroy_ksp_field(ksprop);
  
  /* save solutions if requested */
  /* We don't save the source, because the current file format 
     is not designed for a KS source and a naive sink propagator */
  status = save_wprop_from_wp_field( saveflag, savefile, "", my_ksqs,
				     NULL, wp, io_timing);

  if(status != 0){
    node0_printf("Failed to write propagator\n");
    terminate(1);
  }
  if(saveflag != FORGET)
    node0_printf("Saved naive propagator to %s\n",savefile);

  if(check != CHECK_NO || startflag == FRESH){

    /* Unapply twisted boundary conditions on the fermion links and
       restore conventional KS phases and BC, if changed. */
    boundary_twist_fn(fn_naik, OFF);
    destroy_fn_links(fn_naik);

    /* Turn off staggered phases in the gauge links */
    rephase( OFF );
    invalidate_fermion_links(fn_links);
  }

  /* Why do this ?? */
  clear_qs(my_ksqs);

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

/* In this case we generate a naive propagator from a Dirac propagator
   -- i.e. a Dirac source that carries two spin and two color indices.
   Such a case arises when we calculate an extended naive propagator
   from a Dirac propagator and we need to keep track of the spin and
   color of the Dirac propagator source. */

/* It is assumed that the Dirac source has been prepared in the
   "staggered" basis.  See ext_src/make_ext_src.c */

int get_ksprop4_to_wp_field(int startflag, char startfile[], 
			    int saveflag, char savefile[],
			    wilson_prop_field *wp,
			    ks_prop_field *source,
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
  su3_vector *dst = NULL;
  spin_wilson_vector *swv;
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

  swv  = create_swv_field();
  clear_wp_field(wp);  /* Is this necessary? */

  /* For ksprop_info */
  ksqstmp = *my_qs;
  ksptmp  = *my_ksp;

  /* If we will call the inverter, we need the FN links */
  int naik_epsilon_index;
  imp_ferm_links_t *fn_naik = NULL;
  if(check == CHECK_YES || startflag == FRESH){
    /* Phases must be in before constructing the fermion links */
    rephase( ON );
    invalidate_fermion_links(fn_links);

    if(fn_links == NULL)
      fn_links = create_fermion_links_from_site(my_qic->prec, n_naiks, eps_naik);
    else
      restore_fermion_links_from_site(fn_links, my_qic->prec);
  
    naik_epsilon_index = my_ksp->naik_term_epsilon_index;
    fn_naik = get_fm_links(fn_links, naik_epsilon_index);

    /* Apply twisted boundary conditions and move KS phases, if
       requested */
    /* This operation applies the phase to the boundary FN links */
    set_boundary_twist_fn(fn_naik, mybdry_phase, r0);
    boundary_twist_fn(fn_naik, ON);
  }
	
  /* We do not load a precomputed propagator here.  We don't have an
     appropriate file format for a ksprop4 type */
  
  /* Loop over source colors and spins */
  for(ksource = 0; ksource < my_qs->nsource; ksource++){
    spin_src = convert_ksource_to_spin(ksource);
    color_src = convert_ksource_to_color(ksource);
    for(spin_snk=0;spin_snk<4;spin_snk++){

      node0_printf("%s: spin_src = %d color_src = %d spin_snk %d\n",
		   myname, spin_src, color_src, spin_snk);

      /* Complete the source structure */
      my_qs->spin_snk = spin_snk;

      /* Source for this color */
      su3_vector *src = source->v[color_src];

      /* (Re)construct propagator */
      
      /* Solve for the propagator if the starting guess is zero
	 or we said to solve. */
      if(check == CHECK_YES || startflag == FRESH){

	dst = create_v_field();

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
	   propagator, so we use the less optimized mat_invert_field
	   algorithm, instead of mat_invert_uml_field here to avoid
	   "reconstructing", and so overwriting the odd-site solution.
	   This would be a degradation if the propagator was precomputed
	   in double precision, and we were doing single precision
	   here. If we start using this code to compute the propagator
	   from a fresh start, we could make the selection of the
	   inverter algorithm dependent on the propagator start flag.
	   FRESH -> uml and otherwise cg. */
	
	enum inv_type it = my_qic[0].inv_type;
	if(startflag != FRESH)my_qic[0].inv_type = CGTYPE;
	avs_iters += mat_invert_field(src, dst, my_qic, my_ksp->mass, fn_naik);
	my_qic[0].inv_type = it;

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

      destroy_v_field(dst);

      tot_iters += avs_iters;
    } /* spin_snk */

  } /* ksource -> color_src, spin_src */

  /* save solution if requested */
  /* We don't save the source, because the current file format 
     is not designed for a KS4 source and a naive sink propagator */
  save_wprop_from_wp_field( saveflag, savefile, "", my_qs, NULL, wp, 
			    io_timing);

  /* Free memory */
  if(check == CHECK_YES || startflag == FRESH){
    destroy_fn_links(fn_naik);
    
    /* Turn off staggered phases in the gauge links */
    rephase( OFF );
    invalidate_fermion_links(fn_links);
  }

  /* Why do this ?? */
  clear_qs(my_qs);
  
  destroy_swv_field(swv);
  
  return tot_iters;
}

/* Dump wilson propagator field to file */

void dump_wprop_from_wp_field(int saveflag, int savetype, char savefile[], 
			      wilson_prop_field *wp){
  quark_source dummy_wqs;
#ifdef IOTIME
  int io_timing = 1;
#else
  int io_timing = 0;
#endif

  /* Two output formats are supported.  We can write a file in
     propagator format or source format.  The latter would be suitable
     as an extended source. */

  if(saveflag == FORGET)return;
  
  init_qs(&dummy_wqs);
  wqstmp = dummy_wqs;   /* For clover_info.c */

  if(savetype == DIRAC_PROPAGATOR_FILE){
    /* Propagator format */

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
    
    save_wprop_from_wp_field(saveflag, savefile, "", &dummy_wqs, NULL, wp, io_timing);
  } else {
    /* Source file format */
    int ksource, ncolor, color, spin;
    wilson_vector *wv = create_wv_field();
    char fileinfo[] = "";

    /* Fill relevant source structure parameters */
    dummy_wqs.savetype = savetype;
    /* For the moment we do this only for sources on all time slices. */  
    dummy_wqs.t0 = ALL_T_SLICES;
    dummy_wqs.saveflag = saveflag;
    strcpy(dummy_wqs.save_file,savefile);

    if(w_source_open_dirac(&dummy_wqs, fileinfo) != 0){
      node0_printf("dump_wprop_from_wp_field: Quitting: Can't open source file for writing\n");
      terminate(1);
    }
    ncolor = wp->nc;
    for(ksource = 0; ksource < 4*ncolor; ksource++){
      dummy_wqs.ksource = ksource;
      color = convert_ksource_to_color(ksource);
      spin = convert_ksource_to_spin(ksource);
      copy_wv_from_wp(wv, wp, color, spin);
      w_source_dirac(wv, &dummy_wqs);
    }
    destroy_wv_field(wv);
    w_source_close(&dummy_wqs);
  }
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
  reload_wprop_to_wp_field(rereadflag, savefile, &dummy_wqs, NULL, wp, 1);
  clear_qs(&dummy_wqs); /* Free any allocations */
}
