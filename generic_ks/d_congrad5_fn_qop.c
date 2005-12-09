/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* 2/2005 D. Renner and C. Jung */
/* 12/2005 C. DeTar upgrade to new Level 3 API */

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/d_congrad5_fn_qop.c,v 1.6 2005/12/09 16:59:02 detar Exp $";

/*
 * $Log: d_congrad5_fn_qop.c,v $
 * Revision 1.6  2005/12/09 16:59:02  detar
 * Support new version of qop.h with parity-dependent create_from_raw
 *
 * Revision 1.5  2005/12/04 18:19:57  detar
 * Add Log header
 *
 */

#include "generic_ks_includes.h"
#include <qop.h>

/* Load QOP_FermionLinksAsqtad object from MILC fat and long links */
static void load_fermion_links_asqtad( QOP_FermionLinksAsqtad** qop_links )
{
  su3_matrix **raw_fat_links, **raw_long_links;

  /* Map fat and long links to raw format */

  raw_fat_links  = create_raw_G_from_field_links(t_fatlink,EVENANDODD);
  if(raw_fat_links == NULL)terminate(1);
  raw_long_links = create_raw_G_from_field_links(t_longlink,EVENANDODD);
  if(raw_long_links == NULL)terminate(1);

#if 0
  // Release for memory savings.  Links are recomputed later.
  free_fatlinks();
  free_longlinks();
  valid_fatlinks = 0;
  valid_longlinks = 0;
#endif

  /* Map raw to QOP format */

  *qop_links = QOP_asqtad_create_L_from_raw((Real **)raw_fat_links, 
					    (Real **)raw_long_links,
					    QOP_EVENODD);
  destroy_raw_G(raw_fat_links);
  destroy_raw_G(raw_long_links);
  return;
}

/* Map color vector field from site structure to QOP field */

static void load_V_from_site( QOP_ColorVector** qop, 
			      field_offset milc, 
			      int parity)
{
  su3_vector *raw;

  raw = create_raw_V_from_site(milc, parity);
  if(raw == NULL)terminate(1);
  *qop = QOP_create_V_from_raw((Real *)raw, milc2qop_parity(parity));
  destroy_raw_V(raw);

  return;
}

/* Map color vector from QOP field to site */

static void unload_V_to_site( field_offset milc, QOP_ColorVector *qop,
			      int parity){
  su3_vector *raw;

  raw = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  if(raw == NULL){
    printf("unload_V_to_site: No room for raw vector\n");
    terminate(1);
  }

  QOP_extract_V_to_raw((Real *)raw, qop, milc2qop_parity(parity));
  unload_raw_V_to_site(milc, raw, parity);

  destroy_raw_V(raw);
}

/* Load inversion args for Level 3 inverter */

static void set_qop_invert_arg( QOP_invert_arg_t* qop_invert_arg, 
					   int max_iters, 
					   int max_restart, int milc_parity )
{
  qop_invert_arg->max_iter   = max_iters;
  qop_invert_arg->restart    = max_restart;
  qop_invert_arg->evenodd    = milc2qop_parity(milc_parity);
  qop_invert_arg->final_iter = 0;
  qop_invert_arg->final_flop = 0;
  qop_invert_arg->final_sec  = 0;
}

static void set_qop_resid_arg( QOP_resid_arg_t *qop_resid_arg,
					  Real min_resid_sq )
{
  qop_resid_arg->rsqmin   = min_resid_sq;
}

/* General MILC interface for Level 3 inverter */

static int ks_congrad_qop_generic( QOP_FermionLinksAsqtad* qop_links, 
			    QOP_invert_arg_t *qop_invert_arg, 
			    QOP_resid_arg_t **qop_resid_arg[], 
			    Real *masses[], int nmass[], 
			    QOP_ColorVector **qop_sol[], 
			    QOP_ColorVector* qop_src[], 
			    int nsrc,		    
			    Real* final_rsq_ptr )
{
  int isrc, imass;

  if(nsrc == 1 && nmass[0] == 1)
    QOP_asqtad_invert( qop_links, qop_invert_arg, qop_resid_arg[0][0],
		       masses[0][0], qop_sol[0][0], qop_src[0] );
  else
    QOP_asqtad_invert_multi( qop_links, qop_invert_arg, qop_resid_arg,
			     masses, nmass, qop_sol, qop_src, nsrc );

  /* For now we return the largest value */
  *final_rsq_ptr = 0;
  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++){
      if(*final_rsq_ptr < qop_resid_arg[isrc][imass]->final_rsq)
	*final_rsq_ptr = qop_resid_arg[isrc][imass]->final_rsq;
    }

  return qop_invert_arg->final_iter;
}

/* Map MILC fields to QOP format and call generic QOP driver */

#define MAXSRC 20

int ks_congrad_qop(int niter, Real rsqmin, 
		   Real *masses[], int nmass[], 
		   field_offset milc_srcs[], field_offset *milc_sols[],
		   int nsrc, Real* final_rsq_ptr, int milc_parity )
{
  int isrc, imass;
  QOP_FermionLinksAsqtad *qop_links = NULL;
  QOP_ColorVector **qop_sol[MAXSRC], *qop_src[MAXSRC];
  QOP_invert_arg_t qop_invert_arg;
  QOP_resid_arg_t  **qop_resid_arg[MAXSRC];
  int iterations_used = 0;
  int max_restart = 5;              /* Hard wired at the moment */

#ifdef CGTIME
  
  double dtimec;
  double nflop = 1187;
  if( milc_parity == EVENANDODD ) nflop *= 2;
  
#endif

  if(nsrc > MAXSRC){
    printf("ks_congrad_qop: too many sources\n");
    terminate(1);
  }

  /* Create fat and long links if necessary */

  if( valid_fatlinks  != 1 ) load_fatlinks();
  if( valid_longlinks != 1 ) load_longlinks();

#ifdef CGTIME
  dtimec = -dclock(); 
#endif

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("ks_congrad: Error initializing QOP\n");
    terminate(1);
  }

  /* Set qop_invert_arg and qop_resid_arg */
  set_qop_invert_arg( & qop_invert_arg, niter, 
				 max_restart, milc_parity );
  /* Pointers for residual errors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_resid_arg[isrc] = 
      (QOP_resid_arg_t **)malloc(sizeof(QOP_resid_arg_t *)*nmass[isrc]);
    if(qop_resid_arg[isrc] == NULL){
      printf("ks_congrad_qop: Can't allocate qop_resid_arg_t*\n");
      terminate(1);
    }
    /* For now the residuals the same for all sources and masses */
    for(imass = 0; imass < nmass[isrc]; imass++){
      qop_resid_arg[isrc][imass] = (QOP_resid_arg_t *)malloc(sizeof(QOP_resid_arg_t ));
      if(qop_resid_arg[isrc][imass] == NULL){
	printf("ks_congrad_qop: Can't allocate qop_resid_arg_t\n");
	terminate(1);
      }
      set_qop_resid_arg( qop_resid_arg[isrc][imass], rsqmin );
    }
  }
  
  /* Map MILC fat and long links to QOP links object */

  load_fermion_links_asqtad( &qop_links );

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nmass[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("ks_congrad_qop: Can't allocate qop_sol\n");
      terminate(1);
    }
  }

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    load_V_from_site( &qop_src[isrc], milc_srcs[isrc], milc_parity);
    for(imass = 0; imass < nmass[isrc]; imass++){
      load_V_from_site( &qop_sol[isrc][imass], 
			milc_sols[isrc][imass], milc_parity);
    }
  }
  
  /* Call QOP inverter via restart driver */

  iterations_used = ks_congrad_qop_generic( qop_links, & qop_invert_arg, 
    qop_resid_arg, masses, nmass, qop_sol, qop_src, nsrc, final_rsq_ptr );
  
  /* Map qop solutions to MILC site structure   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      unload_V_to_site( milc_sols[isrc][imass], 
			  qop_sol[isrc][imass], milc_parity );

  /* Free QOP fields  */

  QOP_destroy_L(qop_links);          
  qop_links  = NULL;
  for(isrc = 0; isrc < nsrc; isrc++){
    QOP_destroy_V(qop_src[isrc]);    
    qop_src[isrc] = NULL;
    for(imass = 0; imass < nmass[isrc]; imass++){
      QOP_destroy_V(qop_sol[isrc][imass]);     
      free(qop_sol[isrc]);
      qop_sol[isrc] = NULL;
      free(qop_resid_arg[isrc][imass]);
    }
    free(qop_resid_arg[isrc]);
  }

#ifdef CGTIME
  {
    dtimec += dclock();
    node0_printf("CONGRAD5(total): time = %e iters = %d mflops = %e\n",
		 dtimec,qop_invert_arg.final_iter,
		 qop_invert_arg.final_flop/(1.0e6*dtimec) );
    fflush(stdout);
  }
#endif
  return iterations_used;
}

/* Standard MILC interface for the Asqtad inverter */

int ks_congrad( field_offset milc_src, field_offset milc_sol, Real mass,
	        int niter, Real rsqmin, int milc_parity, Real* final_rsq_ptr )
{
  int iterations_used;
  static Real t_mass;
  Real *masses[1];
  int nmass[1], nsrc;
  field_offset milc_srcs[1], milc_sols0[1], *milc_sols[1];

  /* Set up general source and solution pointers for one mass, one source */
  nsrc = 1;
  milc_srcs[0] = milc_src;

  nmass[0] = 1;
  t_mass = mass;
  masses[0] = &t_mass;

  milc_sols0[0] = milc_sol;
  milc_sols[0] =  milc_sols0;

  iterations_used = ks_congrad_qop( niter, rsqmin, 
				    masses, nmass, milc_srcs,
				    milc_sols, nsrc, final_rsq_ptr,
				    milc_parity );

  total_iters += iterations_used;
  return( iterations_used );
}

