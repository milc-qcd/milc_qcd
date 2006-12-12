/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* NOTE: This code is actually an include file for d_congrad5_fn_qop_F.c
   and d_congrad5_fn_qop_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_CONGRAD_QOP_SITE2SITE
   KS_CONGRAD_QOP_SITE2FIELD   
   KS_CONGRAD_MILC2QOP

*/

/* Redefinitions according to requested precision */

#if ( QOP_Precision == 1 )

#define KS_CONGRAD_QOP_SITE2SITE  ks_congrad_qop_F_site2site
#define KS_CONGRAD_QOP_SITE2FIELD ks_congrad_qop_F_site2field
#define KS_CONGRAD_MILC2QOP       ks_congrad_milc2qop_F
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_F_asqtad_fermion_links
#define MASSREAL float

#else

#define KS_CONGRAD_QOP_SITE2SITE  ks_congrad_qop_D_site2site
#define KS_CONGRAD_QOP_SITE2FIELD ks_congrad_qop_D_site2field
#define KS_CONGRAD_MILC2QOP       ks_congrad_milc2qop_D
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_D_asqtad_fermion_links
#define MASSREAL double

#endif

/* 2/2005 D. Renner and C. Jung */
/* 12/2005 C. DeTar upgrade to new Level 3 API */

/*
 * $Log: d_congrad5_fn_qop_P.c,v $
 * Revision 1.2  2006/12/12 18:07:15  detar
 * Correct mixed precision features.  Add 1sum variant of the QDP inverter.
 *
 * Revision 1.1  2006/12/09 13:52:38  detar
 * Add mixed precision capability for KS inverter in QOP and QDP
 *
 * Revision 1.15  2006/11/13 03:05:26  detar
 * Add timing for remapping and make separate from timing for computation.
 *
 * Revision 1.14  2006/11/04 23:41:14  detar
 * Add QOP and QDP support for FN fermion links
 * Create QDP version of fermion_links_fn_multi
 * Add nrestart parameter for ks_congrad
 *
 * Revision 1.13  2006/10/12 03:43:58  detar
 * Move load_fermion_links_asqtad to (new) fermion_links_asqtad_qop.c
 * to prepare for level 3 link fattening
 *
 * Revision 1.12  2006/09/09 20:12:50  detar
 * Fix qop_invert_arg and split out fermion_links_fn.c from quark_stuff.c
 *
 * Revision 1.11  2006/08/13 15:07:24  detar
 * Adjust entry points for RHMC code and Level 3 multicg wrappers
 *
 * Revision 1.10  2006/03/11 04:24:51  detar
 * Set pointers to null after freeing them
 *
 * Revision 1.9  2006/02/25 16:35:29  detar
 * Fix printf error message
 *
 * Revision 1.8  2005/12/12 23:18:18  detar
 * Correct the name of QOP_asqtad_destroy_L and remove an unused declaration.
 *
 * Revision 1.7  2005/12/09 17:07:33  detar
 * Move cvsheader def
 *
 * Revision 1.6  2005/12/09 16:59:02  detar
 * Support new version of qop.h with parity-dependent create_from_raw
 *
 * Revision 1.5  2005/12/04 18:19:57  detar
 * Add Log header
 *
 */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

#ifdef CGTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

/*#define CGDEBUG*/

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/d_congrad5_fn_qop_P.c,v 1.2 2006/12/12 18:07:15 detar Exp $";


/* Load inversion args for Level 3 inverter */

static void 
set_qop_invert_arg( QOP_invert_arg_t* qop_invert_arg, 
		    int max_iters, 
		    int max_restart, int milc_parity )
{
  qop_invert_arg->max_iter   = max_restart*max_iters;
  qop_invert_arg->restart    = max_iters;
  qop_invert_arg->evenodd    = milc2qop_parity(milc_parity);
}


static QOP_resid_arg_t ***
create_qop_resid_arg( int nsrc, int nmass[], Real min_resid_sq )
{
  QOP_resid_arg_t ***res_arg;
  char myname[] = "create_qop_resid_arg";
  int isrc,imass;

  /* Pointers for residual errors */
  res_arg = (QOP_resid_arg_t ***)malloc(sizeof(QOP_resid_arg_t **)*nsrc);
  if(res_arg == NULL){
    printf("%s(%d): Can't allocate res_arg*\n",myname,this_node);
    terminate(1);
  }
  for(isrc = 0; isrc < nsrc; isrc++){
    res_arg[isrc] = 
      (QOP_resid_arg_t **)malloc(sizeof(QOP_resid_arg_t *)*nmass[isrc]);
    if(res_arg[isrc] == NULL){
      printf("%s(%d): Can't allocate res_arg*\n",myname,this_node);
      terminate(1);
    }
    for(imass = 0; imass < nmass[isrc]; imass++){
      res_arg[isrc][imass] = 
	(QOP_resid_arg_t *)malloc(sizeof(QOP_resid_arg_t ));
      if(res_arg[isrc][imass] == NULL){
	printf("%s(%d): Can't allocate res_arg\n",myname,this_node);
	terminate(1);
      }
      /* For now the residuals are the same for all sources and masses */
      res_arg[isrc][imass]->rsqmin = min_resid_sq;
    }
  }
  return res_arg;
}

static void
destroy_qop_resid_arg(QOP_resid_arg_t ***res_arg, int nsrc, int nmass[])
{
  int isrc, imass;

  for(isrc = 0; isrc < nsrc; isrc++){
    for(imass = 0; imass < nmass[isrc]; imass++){
      free(res_arg[isrc][imass]);
    }
    free(res_arg[isrc]);
  }
  free(res_arg);
}

/* General MILC wrapper for Level 3 inverter */

static int 
ks_congrad_qop_generic( QOP_FermionLinksAsqtad* qop_links, 
			QOP_invert_arg_t *qop_invert_arg,
			QOP_resid_arg_t  ***qop_resid_arg,
			MASSREAL *masses[], int nmass[], 
			QOP_ColorVector **qop_sol[], 
			QOP_ColorVector* qop_src[], 
			int nsrc,		    
			Real* final_rsq_ptr )
{
  int isrc, imass;
  int iters;
  QOP_info_t info;
  
  if(nsrc == 1 && nmass[0] == 1)
    QOP_asqtad_invert( &info, qop_links, qop_invert_arg, qop_resid_arg[0][0],
		       masses[0][0], qop_sol[0][0], qop_src[0] );
  else
    QOP_asqtad_invert_multi( &info, qop_links, qop_invert_arg, qop_resid_arg,
			     masses, nmass, qop_sol, qop_src, nsrc );

  /* For now we return the largest value and total iterations */
  *final_rsq_ptr = 0;
  iters = 0;
  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++){
      if(*final_rsq_ptr < qop_resid_arg[isrc][imass]->final_rsq)
	*final_rsq_ptr = qop_resid_arg[isrc][imass]->final_rsq;
      iters += qop_resid_arg[isrc][imass]->final_iter;
#ifdef CGDEBUG
      if(nsrc > 1 || nmass[isrc] > 1)
	node0_printf("CONGRAD5(src %d,mass %d): iters = %d resid = %e\n",
	       isrc, imass,
	       qop_resid_arg[isrc][imass]->final_iter,
	       qop_resid_arg[isrc][imass]->final_rsq);
#endif
    }

#ifdef CGTIME
  node0_printf("CONGRAD5: time = %e (fn_qop %s) ",
	       info.final_sec,qop_prec[QOP_Precision-1]);
  for(isrc = 0; isrc < nsrc; isrc++)
    node0_printf("nmass[%d] = %d iters = %d ",
		 isrc,nmass[isrc],qop_resid_arg[isrc][0]->final_iter);
  node0_printf("mflops = %e\n", info.final_flop/(1.0e6*info.final_sec) );
  fflush(stdout);
#endif

  return iters;
}

#define MAXSRC 20

/* Map MILC fields to QOP format and call generic QOP driver */
/* This version is for site sources and site sinks */

int
KS_CONGRAD_QOP_SITE2SITE(int niter, int nrestart, Real rsqmin, 
			 MASSREAL *masses[], int nmass[], 
			 field_offset milc_srcs[], 
			 field_offset *milc_sols[],
			 int nsrc, Real* final_rsq_ptr, int milc_parity )
{
  int isrc, imass;
  QOP_FermionLinksAsqtad *qop_links;
  QOP_ColorVector **qop_sol[MAXSRC], *qop_src[MAXSRC];
  int iterations_used = 0;
  QOP_invert_arg_t qop_invert_arg;
  QOP_resid_arg_t  ***qop_resid_arg;
  double remaptime;

  if(nsrc > MAXSRC){
    printf("ks_congrad_qop_site2site: too many sources\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("ks_congrad_qop_site2site: Error initializing QOP\n");
    terminate(1);
  }

  /* Map MILC fat and long links to QOP links object */

  qop_links = CREATE_QOP_ASQTAD_FERMION_LINKS( );

  /* Set qop_invert_arg */
  set_qop_invert_arg( & qop_invert_arg, niter, 
		      nrestart, milc_parity );
  
  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nmass, rsqmin );

  remaptime = -dclock(); 

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nmass[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("ks_congrad_qop_site2site: Can't allocate qop_sol\n");
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
  
  /* Call QOP inverter */

  remaptime += dclock();
  iterations_used = ks_congrad_qop_generic( qop_links, &qop_invert_arg,
    qop_resid_arg, masses, nmass, qop_sol, qop_src, nsrc, final_rsq_ptr );
  remaptime -= dclock();
  
  /* Map qop solutions to MILC site structure   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      unload_V_to_site( milc_sols[isrc][imass], 
			  qop_sol[isrc][imass], milc_parity );

  /* Free QOP fields  */

  for(isrc = 0; isrc < nsrc; isrc++){
    QOP_destroy_V(qop_src[isrc]);    
    qop_src[isrc] = NULL;
    for(imass = 0; imass < nmass[isrc]; imass++){
      QOP_destroy_V(qop_sol[isrc][imass]);     
      free(qop_sol[isrc]);
      qop_sol[isrc] = NULL;
    }
  }

  remaptime += dclock();
#ifdef CGTIME
#ifdef REMAP
    node0_printf("CGREMAP:  time = %e\n",remaptime);
#endif
#endif
  destroy_qop_resid_arg(qop_resid_arg, nsrc, nmass);
  qop_resid_arg = NULL;

  return iterations_used;
}

/* Map MILC fields to QOP format and call generic QOP driver */
/* This version is for site sources and field sinks */

int KS_CONGRAD_QOP_SITE2FIELD(int niter, int nrestart, Real rsqmin, 
			      MASSREAL *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[],
			      int nsrc, Real* final_rsq_ptr, int milc_parity )
{
  int isrc, imass;
  QOP_FermionLinksAsqtad *qop_links;
  QOP_ColorVector **qop_sol[MAXSRC], *qop_src[MAXSRC];
  int iterations_used = 0;
  double remaptime;
  QOP_resid_arg_t  ***qop_resid_arg;
  QOP_invert_arg_t qop_invert_arg;

  if(nsrc > MAXSRC){
    printf("ks_congrad_qop_site2field: too many sources\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("ks_congrad_qop_site2field: Error initializing QOP\n");
    terminate(1);
  }

  /* Map MILC fat and long links to QOP links object */

  qop_links = CREATE_QOP_ASQTAD_FERMION_LINKS();

  /* Set qop_invert_arg */
  set_qop_invert_arg( & qop_invert_arg, niter, 
		      nrestart, milc_parity );

  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nmass, rsqmin );

  remaptime = -dclock(); 

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nmass[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("ks_congrad_qop_site2field: Can't allocate qop_sol\n");
      terminate(1);
    }
  }

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    load_V_from_site( &qop_src[isrc], milc_srcs[isrc], milc_parity);
    for(imass = 0; imass < nmass[isrc]; imass++){
      load_V_from_field( &qop_sol[isrc][imass], 
			 milc_sols[isrc][imass], milc_parity);
    }
  }
  
  /* Call QOP inverter */

  remaptime += dclock();
  iterations_used = ks_congrad_qop_generic( qop_links, &qop_invert_arg,
     qop_resid_arg, masses, nmass, qop_sol, qop_src, nsrc, final_rsq_ptr );
  remaptime -= dclock();
  
  /* Map qop solutions to MILC field   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      unload_V_to_field( milc_sols[isrc][imass], 
			 qop_sol[isrc][imass], milc_parity );

  /* Free QOP fields  */

  for(isrc = 0; isrc < nsrc; isrc++){
    QOP_destroy_V(qop_src[isrc]);    
    qop_src[isrc] = NULL;
    for(imass = 0; imass < nmass[isrc]; imass++){
      QOP_destroy_V(qop_sol[isrc][imass]);     
    }
    free(qop_sol[isrc]);
    qop_sol[isrc] = NULL;
  }

  remaptime += dclock();

#ifdef CGTIME
#ifdef REMAP
    node0_printf("CGREMAP:  time = %e\n",remaptime);
#endif
#endif

  destroy_qop_resid_arg(qop_resid_arg, nsrc, nmass);
  qop_resid_arg = NULL;

  return iterations_used;
}


int 
KS_CONGRAD_MILC2QOP( field_offset milc_src, field_offset milc_sol, Real mass,
		     int niter, int nrestart, Real rsqmin, 
		     int milc_parity, Real* final_rsq_ptr )
{
  int iterations_used;
  static MASSREAL t_mass;
  MASSREAL *masses[1];
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

  iterations_used = 
    KS_CONGRAD_QOP_SITE2SITE( niter, nrestart, rsqmin, 
			      masses, nmass, milc_srcs,
			      milc_sols, nsrc, final_rsq_ptr,
			      milc_parity );
  
  return  iterations_used;
}
