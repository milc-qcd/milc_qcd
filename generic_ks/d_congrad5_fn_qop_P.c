/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* NOTE: This code is actually an include file for d_congrad5_fn_qop_F.c
   and d_congrad5_fn_qop_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names)

   KS_CONGRAD_QOP_SITE2SITE
   KS_CONGRAD_QOP_SITE2FIELD   
   KS_CONGRAD_QOP_FIELD2FIELD   
   KS_CONGRAD_MILC2QOP

*/

/* Redefinitions according to requested precision */

#if ( QOP_Precision == 1 )

#define KS_CONGRAD_QOP_SITE2SITE  ks_congrad_qop_F_site2site
#define KS_CONGRAD_QOP_SITE2FIELD ks_congrad_qop_F_site2field
#define KS_CONGRAD_QOP_FIELD2FIELD ks_congrad_qop_F_field2field
#define KS_CONGRAD_MILCFIELD2QOP  ks_congrad_milcfield2qop_F
#define KS_CONGRAD_MILC2QOP       ks_congrad_milc2qop_F
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_F_asqtad_fermion_links
#define MASSREAL float
#define QOP_L qop_F_l

#else

#define KS_CONGRAD_QOP_SITE2SITE  ks_congrad_qop_D_site2site
#define KS_CONGRAD_QOP_SITE2FIELD ks_congrad_qop_D_site2field
#define KS_CONGRAD_QOP_FIELD2FIELD ks_congrad_qop_D_field2field
#define KS_CONGRAD_MILCFIELD2QOP  ks_congrad_milcfield2qop_D
#define KS_CONGRAD_MILC2QOP       ks_congrad_milc2qop_D
#define CREATE_QOP_ASQTAD_FERMION_LINKS create_qop_D_asqtad_fermion_links
#define MASSREAL double
#define QOP_L qop_D_l

#endif

/* 2/2005 D. Renner and C. Jung */
/* 12/2005 C. DeTar upgrade to new Level 3 API */

/*
 * $Log: d_congrad5_fn_qop_P.c,v $
 * Revision 1.5  2007/11/09 16:42:41  detar
 * Pull FN link calculation out of inverters
 *
 * Revision 1.4  2007/10/09 20:10:14  detar
 * Add ferm_links_t and ks_action_paths structures and pass them as params
 *
 * Revision 1.3  2007/05/21 05:06:48  detar
 * Change stopping condition to true residual.
 *
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

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_ks/d_congrad5_fn_qop_P.c,v 1.5 2007/11/09 16:42:41 detar Exp $";


/* Load inversion args for Level 3 inverter */

static void 
set_qop_invert_arg( QOP_invert_arg_t* qop_invert_arg, 
		    quark_invert_control *qic )
{
  qic->nrestart = 10;
  qop_invert_arg->max_iter     = qic->nrestart * qic->max;
  qop_invert_arg->restart      = qic->max;
  qop_invert_arg->max_restarts = qic->nrestart;
  qop_invert_arg->evenodd      = milc2qop_parity(qic->parity);
}


static QOP_resid_arg_t ***
create_qop_resid_arg( int nsrc, int nmass[], 
		      quark_invert_control *qic )
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
      res_arg[isrc][imass]->rsqmin = qic->resid;
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
			int nsrc, quark_invert_control *qic )
{
  int isrc, imass;
  int iters;
  QOP_info_t info;

#ifndef OLD_QOPQDP_NORM  
  /* Since version 0.9.0 the conventions for QOP_asqtad_invert* do not
     give results consistent with other calls to ks_congrad and
     ks_multicg... when the parity is EVENODD.  You should change the
     call and use mat_invert_uml instead.  */

  if(qop_invert_arg->evenodd == QOP_EVENODD){
    printf("ks_congrad_qop_generic: EVENANDODD is currently not supported.\n");
    terminate(1);
  }
#endif

  if(nsrc == 1 && nmass[0] == 1)
    QOP_asqtad_invert( &info, qop_links, qop_invert_arg, qop_resid_arg[0][0],
		       masses[0][0], qop_sol[0][0], qop_src[0] );
  else
    QOP_asqtad_invert_multi( &info, qop_links, qop_invert_arg, qop_resid_arg,
			     masses, nmass, qop_sol, qop_src, nsrc );

  /* For now we return the largest value and total iterations */
  qic->final_rsq = 0;
  qic->converged = 1;
  iters = 0;
  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++){
      if(this_node == 0){
	if(qic->resid < qop_resid_arg[isrc][imass]->final_rsq){
	  /** suppressed until we clear up the normalization conventions
	  qic->converged = 0;
	  printf("ks_congrad_qop_generic: NOT CONVERGED (src %d, mass %d),\n",
		 isrc, imass);
	  printf("after %d iters (max iters %d max %d restarts) rsq. = %e wanted %e\n",
		 qop_resid_arg[isrc][imass]->final_iter,
		 qic->max, qic->nrestart,
		 qop_resid_arg[isrc][imass]->final_rsq,qic->resid);
	  **/
	}
      }
      if(qic->final_rsq < qop_resid_arg[isrc][imass]->final_rsq)
	qic->final_rsq = qop_resid_arg[isrc][imass]->final_rsq;
      iters += qop_resid_arg[isrc][imass]->final_iter;
#ifdef CGDEBUG
      if(nsrc > 1 || nmass[isrc] > 1)
	node0_printf("CONGRAD5(src %d, mass %d): iters = %d resid = %e\n",
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

#ifndef OLD_QOPQDP_NORM
/* Starting with qopqdp-0.9.0 the normalization of the solution is changed */

#define NORMFACT(a) 4.*(a)*(a)

static
void qop_to_milc_normalization_site(field_offset sol, MASSREAL mass,
				    int parity)
{
  site *s;
  int i;
  Real x = 1/(NORMFACT(mass));

  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector((su3_vector *)F_PT(s,sol), x,
			   (su3_vector *)F_PT(s,sol));
  }
}

static
void qop_to_milc_normalization_field(su3_vector *sol, MASSREAL mass,
				     int parity)
{
  site *s;
  int i;
  Real x = 1/(NORMFACT(mass));

  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector(sol+i, x, sol+i);
  }
}

static
void milc_to_qop_normalization_site(field_offset sol, MASSREAL mass,
				    int parity)
{
  site *s;
  int i;
  Real x = NORMFACT(mass);

  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector((su3_vector *)F_PT(s,sol), x,
			   (su3_vector *)F_PT(s,sol));
  }
}

static
void milc_to_qop_normalization_field(su3_vector *sol, MASSREAL mass,
				     int parity)
{
  site *s;
  int i;
  Real x = NORMFACT(mass);

  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector(sol+i, x, sol+i);
  }
}

#endif

#define MAXSRC 20

/* Map MILC fields to QOP format and call generic QOP driver */
/* This version is for site sources and site sinks */

int
KS_CONGRAD_QOP_SITE2SITE(quark_invert_control *qic,
			 MASSREAL *masses[], int nmass[], 
			 field_offset milc_srcs[], 
			 field_offset *milc_sols[], int nsrc,
			 ferm_links_t *fn)
{
  int isrc, imass;
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

  CREATE_QOP_ASQTAD_FERMION_LINKS(fn);

  /* Set qop_invert_arg */
  set_qop_invert_arg( &qop_invert_arg, qic );
  
  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nmass, qic );

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

#ifndef OLD_QOPQDP_NORM
  /* Convert proposed solutions from MILC to QOP normalization convention */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      milc_to_qop_normalization_site(milc_sols[isrc][imass],
				     masses[isrc][imass],qic->parity);
#endif

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_src[isrc] = create_V_from_site( milc_srcs[isrc], qic->parity);
    for(imass = 0; imass < nmass[isrc]; imass++){
      qop_sol[isrc][imass] = 
	create_V_from_site( milc_sols[isrc][imass], qic->parity);
    }
  }
  
  /* Call QOP inverter */

  remaptime += dclock();
  iterations_used = ks_congrad_qop_generic( fn->QOP_L, &qop_invert_arg,
    qop_resid_arg, masses, nmass, qop_sol, qop_src, nsrc, qic );
  remaptime -= dclock();
  
  /* Map qop solutions to MILC site structure   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      unload_V_to_site( milc_sols[isrc][imass], 
			  qop_sol[isrc][imass], qic->parity );

#ifndef OLD_QOPQDP_NORM
  /* Convert solutions to MILC ks_congrad normalization convention */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      qop_to_milc_normalization_site(milc_sols[isrc][imass],
				     masses[isrc][imass],qic->parity);
#endif

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

int KS_CONGRAD_QOP_SITE2FIELD(quark_invert_control *qic,
			      MASSREAL *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[], int nsrc,
			      ferm_links_t *fn )
{
  int isrc, imass;
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

  CREATE_QOP_ASQTAD_FERMION_LINKS(fn);

  /* Set qop_invert_arg */
  set_qop_invert_arg( & qop_invert_arg, qic );

  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nmass, qic );

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

#ifndef OLD_QOPQDP_NORM
  /* Convert proposed solution from MILC to QOP normalization convention */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      milc_to_qop_normalization_field(milc_sols[isrc][imass],
				      masses[isrc][imass],qic->parity);
#endif

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_src[isrc] = create_V_from_site( milc_srcs[isrc], qic->parity);
    for(imass = 0; imass < nmass[isrc]; imass++){
      qop_sol[isrc][imass] = 
	create_V_from_field( milc_sols[isrc][imass], qic->parity);
    }
  }
  
  /* Call QOP inverter */

  remaptime += dclock();
  iterations_used = ks_congrad_qop_generic( fn->QOP_L, &qop_invert_arg,
     qop_resid_arg, masses, nmass, qop_sol, qop_src, nsrc, qic );
  remaptime -= dclock();
  
  /* Map qop solutions to MILC field   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      unload_V_to_field( milc_sols[isrc][imass], 
			 qop_sol[isrc][imass], qic->parity );

#ifndef OLD_QOPQDP_NORM
  /* Convert solutions to MILC ks_congrad normalization convention */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      qop_to_milc_normalization_field(milc_sols[isrc][imass],
				      masses[isrc][imass],qic->parity);
#endif

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


/* Map MILC fields to QOP format and call generic QOP driver */
/* This version is for field sources and sinks */

int KS_CONGRAD_QOP_FIELD2FIELD(quark_invert_control *qic,
			       MASSREAL *masses[], int nmass[], 
			       su3_vector *milc_srcs[], 
			       su3_vector **milc_sols[], int nsrc,
			       ferm_links_t *fn )
{
  int isrc, imass;
  QOP_ColorVector **qop_sol[MAXSRC], *qop_src[MAXSRC];
  int iterations_used = 0;
  double remaptime;
  QOP_resid_arg_t  ***qop_resid_arg;
  QOP_invert_arg_t qop_invert_arg;

  if(nsrc > MAXSRC){
    printf("ks_congrad_qop_field2field: too many sources\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("ks_congrad_qop_field2field: Error initializing QOP\n");
    terminate(1);
  }

  /* Map MILC fat and long links to QOP links object */

  CREATE_QOP_ASQTAD_FERMION_LINKS(fn);

  /* Set qop_invert_arg */
  set_qop_invert_arg( & qop_invert_arg, qic );

  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nmass, qic );

  remaptime = -dclock(); 

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_ColorVector **)malloc(sizeof(QOP_ColorVector *)*nmass[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("ks_congrad_qop_field2field: Can't allocate qop_sol\n");
      terminate(1);
    }
  }

#ifndef OLD_QOPQDP_NORM
  /* Convert proposed solutions from MILC to QOP normalization convention */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      milc_to_qop_normalization_field(milc_sols[isrc][imass],
				      masses[isrc][imass],qic->parity);
#endif

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_src[isrc] = create_V_from_field( milc_srcs[isrc], qic->parity);
    for(imass = 0; imass < nmass[isrc]; imass++){
      qop_sol[isrc][imass] = 
	create_V_from_field( milc_sols[isrc][imass], qic->parity);
    }
  }
  
  /* Call QOP inverter */

  remaptime += dclock();
  iterations_used = ks_congrad_qop_generic( fn->QOP_L, &qop_invert_arg,
     qop_resid_arg, masses, nmass, qop_sol, qop_src, nsrc, qic );
  remaptime -= dclock();
  
  /* Map qop solutions to MILC field   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      unload_V_to_field( milc_sols[isrc][imass], 
			 qop_sol[isrc][imass], qic->parity );

#ifndef OLD_QOPQDP_NORM
  /* Convert solutions to MILC ks_congrad normalization convention */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(imass = 0; imass < nmass[isrc]; imass++)
      qop_to_milc_normalization_field(milc_sols[isrc][imass],
				      masses[isrc][imass],qic->parity);
#endif

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
KS_CONGRAD_MILCFIELD2QOP( su3_vector *milc_src, su3_vector *milc_sol, 
			  quark_invert_control *qic, Real mass,
			  ferm_links_t *fn )
{
  int iterations_used;
  static MASSREAL t_mass;
  MASSREAL *masses[1];
  int nmass[1], nsrc;
  su3_vector *milc_srcs[1], *milc_sols0[1], **milc_sols[1];

  /* Set up general source and solution pointers for one mass, one source */
  nsrc = 1;
  milc_srcs[0] = milc_src;

  nmass[0] = 1;
  t_mass = mass;
  masses[0] = &t_mass;

  milc_sols0[0] = milc_sol;
  milc_sols[0] =  milc_sols0;

  iterations_used = 
    KS_CONGRAD_QOP_FIELD2FIELD( qic, masses, nmass, milc_srcs,
				milc_sols, nsrc, fn );
  
  return  iterations_used;
}

int 
KS_CONGRAD_MILC2QOP( field_offset milc_src, field_offset milc_sol, 
		     quark_invert_control *qic, Real mass,
		     ferm_links_t *fn )
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
    KS_CONGRAD_QOP_SITE2SITE( qic, masses, nmass, milc_srcs,
			      milc_sols, nsrc, fn );
  
  return  iterations_used;
}
