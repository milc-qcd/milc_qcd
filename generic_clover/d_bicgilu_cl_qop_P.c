/******* d_bicgilu_cl_qop_P.c - BiCG ILU ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* NOTE: This code is actually an include file for d_bicgilu_cl_F.c
   and d_bicgilu_cl_D.c, so any edits should be consistent with this
   purpose. 

   NOTE: Several local procedures are defined "static".  That makes it
   possible to use the same name without confusion in the "D" and "F"
   variants of the compilation.  If you need to make them global
   externals, they must be given precision-specific names as we did
   for BICGILU_MILC2QOP.  Otherwise, the linker will find that they
   are multiply defined.

*/

/* Precision-dependent types

   NYREAL
   MYSU3_MATRIX
*/

#include "generic_clover_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_clover_qop.h"
#include "../include/generic_qopqdp.h"
#include <string.h>

/*#define CG_DEBUG*/

/* Redefinitions according to requested precision */

#if ( QOP_PrecisionInt == 1 )

#define BICGILU_CL_QOP_GENERIC bicgilu_cl_qop_generic_F
#define BICGILU_MILC2QOP   bicgilu_cl_milc2qop_F
#define MYREAL float
#define MYSU3_MATRIX fsu3_matrix
#define CREATE_RAW4_G_FROM_SITE create_raw4_F_G_from_site
#define DESTROY_RAW4_G destroy_raw4_F_G

#else

#define BICGILU_CL_QOP_GENERIC bicgilu_cl_qop_generic_D
#define BICGILU_MILC2QOP   bicgilu_cl_milc2qop_D
#define MYREAL double
#define MYSU3_MATRIX dsu3_matrix
#define CREATE_RAW4_G_FROM_SITE create_raw4_D_G_from_site
#define DESTROY_RAW4_G destroy_raw4_D_G

#endif

/* 4/29/07 C. DeTar */

/*
 * $Log: d_bicgilu_cl_qop_P.c,v $
 * Revision 1.10  2012/11/23 23:45:47  detar
 * Add OK action support (IFLA)
 *
 * Revision 1.9  2012/02/16 16:49:56  detar
 * Initialize QOP_info
 *
 * Revision 1.8  2011/11/29 20:23:54  detar
 * Set uniform convention for qic residuals.  See include/generic_quark_types.h
 *
 * Revision 1.7  2009/06/04 16:37:09  detar
 * Make clover term persistent. Accommodate changes to generic_clover/make_clov2.c
 *
 * Revision 1.6  2009/04/05 18:12:39  detar
 * Move #if 0
 *
 * Revision 1.5  2008/04/18 23:11:54  detar
 * Remove QOP_verbose setting
 *
 * Revision 1.4  2008/04/18 15:37:03  detar
 * Support qopqdp-0.11.2
 *
 * Revision 1.3  2008/03/28 15:48:09  detar
 * Report absolute and relative residuals
 *
 * Revision 1.2  2007/12/14 04:37:48  detar
 * Add final restart member to qic.
 *
 * Revision 1.1  2007/05/21 04:43:58  detar
 * Support Level 3 inversion
 *
 */

#ifdef CGTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_clover/d_bicgilu_cl_qop_P.c,v 1.10 2012/11/23 23:45:47 detar Exp $";

#if 1

/********************************************************************/
/* Load Wilson clover parameters                                    */
/********************************************************************/

static void 
load_qop_wilson_coeffs(QOP_wilson_coeffs_t *c, Real clov_c)
{
  c->clov_s       = clov_c;
  c->clov_t       = clov_c;
  c->aniso        = 0.;
}

/********************************************************************/
/* Load Wilson clover parameters (IFLA) Bugra 09/02/09              */
/********************************************************************/
static void
load_qop_wilson_ifla_coeffs(QOP_wilson_ifla_coeffs_t *c, 
			    newaction_ifla_param *nap)
{  
  c->kapifla = nap->kapifla;
  c->kappa_s = nap->kappa_s;
  c->kappa_t = nap->kappa_t;
  c->r_s     = nap->r_s;
  c->r_t     = nap->r_t;
  c->zeta    = nap->zeta;
  c->c_E     = nap->c_E;
  c->c_B     = nap->c_B;
  c->c_1     = nap->c_1;
  c->c_2     = nap->c_2;
  c->c_3     = nap->c_3;
  c->c_4     = nap->c_4;
  c->c_5     = nap->c_5;
  c->c_EE    = nap->c_EE;
  c->u0      = nap->u0;
};

/********************************************************************/
/* Create Wilson fermion links using QOP                            */
/********************************************************************/

static QOP_FermionLinksWilson *
create_qop_wilson_fermion_links( Real clov )
{
  QOP_FermionLinksWilson *qop_links = NULL;
  QOP_info_t info = {0., 0., 0, 0, 0};
  QOP_GaugeField *links;
  QOP_wilson_coeffs_t coeffs;
  double remaptime;

  /* Load coeffs structure */
  load_qop_wilson_coeffs(&coeffs, clov);

  /* Map SU(3) gauge field to G type */
  remaptime = -dclock(); 
  links = create_G_from_site4(F_OFFSET(link),EVENANDODD);
  remaptime += dclock();

  /* Create links */
  qop_links = QOP_wilson_create_L_from_G(&info, &coeffs, links);

  QOP_destroy_G(links);

#ifdef FLTIME
#ifdef REMAP
  node0_printf("FLREMAP:  time = %e\n",remaptime);
#endif
  if(info.final_sec > 0){
    double mflops = 0.;
    mflops = (Real)info.final_flop/(1e6*info.final_sec);
    node0_printf("FLTIME:  time = %e (cl_qop) terms = 1 mflops = %e\n",
		 info.final_sec, mflops );
  }
#endif
  return qop_links;
}

#else

/********************************************************************/
/* Use plain MILC code to construct clover term */
/********************************************************************/

static QOP_FermionLinksWilson *
create_qop_wilson_fermion_links( Real clov ){

  clover *milc_clov = gen_clov;
  MYSU3_MATRIX **raw_links;
  MYREAL *raw_clov;
  QOP_FermionLinksWilson *qop_links;

  /* Construct raw QOP clover term from MILC clover term */
  if(clov == 0){
    raw_clov = NULL;
  }
  else{
    raw_clov = (MYREAL *)malloc(72*sites_on_node*sizeof(MYREAL));
    if(raw_clov == NULL){
      printf("create_qop_wilson_fermion_links(%d): no room for raw_clov\n",
	     this_node);
      terminate(1);
    }
    // milc_clov = create_clov(); /* Note Real clov has no kappa factor! */
    if(milc_clov == NULL) terminate(1);
    compute_clov(milc_clov,clov);
    map_milc_clov_to_qop_raw(raw_clov, milc_clov);
    //free_this_clov(milc_clov);
  }

  raw_links = CREATE_RAW4_G_FROM_SITE(F_OFFSET(link), EVENANDODD);
  if(raw_links == NULL)terminate(1);

  /* Map QOP/QDP raw to QOP/QDP structure */

  qop_links = QOP_wilson_create_L_from_raw((MYREAL **)raw_links, raw_clov, 
					   QOP_EVENODD);
  
  DESTROY_RAW4_G(raw_links); raw_links = NULL;
  free(raw_clov);

  return qop_links;
}

#endif

/********************************************************************/
/* Destroy fermion links                                            */
/********************************************************************/

static void
destroy_qop_wilson_fermion_links( QOP_FermionLinksWilson *qop_links )
{
  QOP_wilson_destroy_L(qop_links);
  qop_links = NULL;
}

/********************************************************************/
/* Load inversion args for Level 3 inverter                         */
/********************************************************************/

static void 
set_qop_invert_arg( QOP_invert_arg_t* qop_invert_arg, 
		    quark_invert_control *qic, int milc_parity )
{
  qop_invert_arg->max_iter     = qic->nrestart*qic->max;
  qop_invert_arg->restart      = qic->max;
  qop_invert_arg->max_restarts = qic->nrestart;
  qop_invert_arg->evenodd      = milc2qop_parity(milc_parity);
}

/********************************************************************/
/* Load residual values for stopping                                */
/********************************************************************/

static QOP_resid_arg_t ***
create_qop_resid_arg( int nsrc, int nmass[], Real min_resid_sq,
		      Real min_relresid_sq)
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
      res_arg[isrc][imass]->relmin = min_relresid_sq;
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

/********************************************************************/
/* General MILC wrapper for Level 3 inverter                        */
/********************************************************************/

void
BICGILU_CL_QOP_GENERIC( int prop_type, QOP_info_t *info,
			QOP_FermionLinksWilson *qop_links, 
			QOP_invert_arg_t *qop_invert_arg,
			QOP_resid_arg_t  ***qop_resid_arg,
			void *dmps[], int nkappa[],
			QOP_DiracFermion **qop_sol[], 
			QOP_DiracFermion *qop_src[], 
			int nsrc)
{
  MYREAL **kappas;
  int isrc, ikappa;
  char myname[] = "bicgilu_cl_qop_generic";

  
  kappas = (MYREAL **)malloc(sizeof(MYREAL *)*nsrc);
  for(isrc = 0; isrc < nsrc; isrc++)
    kappas[isrc] = (MYREAL *)malloc(sizeof(MYREAL)*nkappa[isrc]);

#ifdef CG_DEBUG
  QOP_verbose(QOP_VERB_HI);
#endif
  
  if(prop_type == CLOVER_TYPE){
    for(isrc = 0; isrc < nsrc; isrc++)
      for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
	dirac_clover_param *dcp = (dirac_clover_param *)dmps[0];
	kappas[isrc][ikappa] = dcp->Kappa;
      }

    if(nsrc == 1 && nkappa[0] == 1)
      QOP_wilson_invert( info, qop_links, qop_invert_arg, qop_resid_arg[0][0],
			 kappas[0][0], qop_sol[0][0], qop_src[0] );
    else
      QOP_wilson_invert_multi( info, qop_links, qop_invert_arg, qop_resid_arg,
			       kappas, nkappa, qop_sol, qop_src, nsrc );

  } else if(prop_type == IFLA_TYPE){

    QOP_wilson_ifla_coeffs_t dcof;
    /* If there are multiple kappas, we assume that the parameters
       are the same for all kappas in this solve */
    load_qop_wilson_ifla_coeffs(&dcof, (newaction_ifla_param *)dmps[0] );
    kappas[0][0] = dcof.kapifla;

    if(nsrc == 1 && nkappa[0] == 1){
      QOP_wilson_ifla_invert( info, qop_links, qop_invert_arg, qop_resid_arg[0][0],
			      kappas[0][0], &dcof, qop_sol[0][0], qop_src[0] );
    } else {
      node0_printf("%s: Multishift bicgilu with OK action not provided\n", myname);
      terminate(1);
    }

  } else {

    node0_printf("%s: Unrecognized prop type\n", myname);
    terminate(1);

  }

  for(isrc = 0; isrc < nsrc; isrc++)
    free(kappas[isrc]);
  free(kappas);
}

/* temporary hack until we get a more flexible QOP inverter that can
   handle our gamma matrix conventions */

static void
gamma5_flip(wilson_vector *milc, int parity){
  int i,c;
  site *s;
  FORSOMEPARITY(i,s,parity){
    for(c = 0; c < 3; c++){
      CNEGATE(milc[i].d[2].c[c],milc[i].d[2].c[c]);
      CNEGATE(milc[i].d[3].c[c],milc[i].d[3].c[c]);
    }
  }
}

#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MAXSRC 20


/********************************************************************/
/* Map MILC fields to QOP format and call generic QOP driver        */
/********************************************************************/

static int
bicgilu_cl_qop(int prop_type, int nsrc, int nkappa[], 
	       quark_invert_control *qic,
	       void *dmps[], wilson_vector *milc_srcs[], 
	       wilson_vector **milc_sols[],
	       int *final_restart,
               Real *final_rsq_ptr, Real* final_relrsq_ptr, 
	       int milc_parity )
{
  int isrc, ikappa, iters = 0;
  QOP_FermionLinksWilson *qop_links;
  QOP_DiracFermion **qop_sol[MAXSRC], *qop_src[MAXSRC];
  QOP_invert_arg_t qop_invert_arg;
  QOP_resid_arg_t  ***qop_resid_arg;
  QOP_info_t info = {0., 0., 0, 0, 0};
  double remaptime;
  int i;
  site *s;
  static float t_kappa;
  float *kappas[1] = { &t_kappa };
  char myname[] = "bicgilu_cl_qop";

  if(nsrc > MAXSRC){
    printf("%s: too many sources\n",myname);
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s: Error initializing QOP\n",myname);
    terminate(1);
  }

  /* Create QOP links object and extract kappa */

  if(prop_type == CLOVER_TYPE){
    /* If there are multiple kappas, we assume that the Clov_c and u0
       are the same for all kappas in this solve */
    dirac_clover_param *dcp 
      = (dirac_clover_param *)dmps[0]; /* Cast pass-through pointer */
    Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
    Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
    Real clov = Clov_c/(U0*U0*U0); /* Full clover coefficient */
    t_kappa = dcp->Kappa;
    
    qop_links = create_qop_wilson_fermion_links( clov );

  } else { /* IFLA (OK) type */
    newaction_ifla_param *nap 
      = (newaction_ifla_param *)dmps[0]; /* Cast pass-through pointer */
    t_kappa = nap->kapifla;

    qop_links = create_qop_wilson_fermion_links( 0 );

  }

  /* Set qop_invert_arg */
  set_qop_invert_arg( &qop_invert_arg, qic, milc_parity );
  
  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nkappa, 
					(qic->resid)*(qic->resid),
					(qic->relresid)*(qic->relresid));

  remaptime = -dclock(); 

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_DiracFermion **)malloc(sizeof(QOP_DiracFermion *)*nkappa[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("%s: Can't allocate qop_sol\n",myname);
      terminate(1);
    }
  }

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    gamma5_flip(milc_srcs[isrc], milc_parity);  /* compensate for QOP gamma */
    qop_src[isrc] = create_D_from_field( milc_srcs[isrc], milc_parity);
    gamma5_flip(milc_srcs[isrc], milc_parity);  /* restore the source */
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      /* Adjust normalization for MILC conventions */
      gamma5_flip(milc_sols[isrc][ikappa], milc_parity);  /* compensate for QOP gamma */
      FORALLSITES(i,s){
	scalar_mult_wvec( milc_sols[isrc][ikappa]+i, 2.*kappas[isrc][ikappa],
			  milc_sols[isrc][ikappa]+i);
      }
      qop_sol[isrc][ikappa] = 
	create_D_from_field( milc_sols[isrc][ikappa], milc_parity);
    }
  }

  /* Call QOP inverter */

  remaptime += dclock();
  BICGILU_CL_QOP_GENERIC( prop_type, &info, qop_links, 
	&qop_invert_arg, qop_resid_arg, dmps, nkappa, qop_sol, qop_src, nsrc);
  remaptime -= dclock();

#ifdef CGTIME
  node0_printf("CGTIME: time = %e (wilson_qop %s) ",
	       info.final_sec,qop_prec[QOP_PrecisionInt-1]);
  for(isrc = 0; isrc < nsrc; isrc++)
    node0_printf("nkappa[%d] = %d iters = %d ",
		 isrc,nkappa[isrc],qop_resid_arg[isrc][0]->final_iter);
  node0_printf("mflops = %e\n", info.final_flop/(1.0e6*info.final_sec) );
  fflush(stdout);
#endif

  /* Map qop solutions to MILC fields   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      unload_D_to_field( milc_sols[isrc][ikappa], 
			 qop_sol[isrc][ikappa], milc_parity );
      /* Adjust normalization for MILC conventions */
      gamma5_flip(milc_sols[isrc][ikappa], milc_parity);  /* compensate for QOP gamma */
      FORALLSITES(i,s){
	scalar_mult_wvec( milc_sols[isrc][ikappa]+i, 1/(2.*kappas[isrc][ikappa]),
			  milc_sols[isrc][ikappa]+i);
      }
    }

  /* Collect performance statistics */

  /* For now we return the largest value and total iterations */
  *final_rsq_ptr = 0;
  *final_relrsq_ptr = 0;
  *final_restart = 0;
  iters = 0;
  for(isrc = 0; isrc < nsrc; isrc++)
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      /* QOP routines return the ratios of the squared norms */
      *final_rsq_ptr =    MAX(*final_rsq_ptr, qop_resid_arg[isrc][ikappa]->final_rsq);
      *final_relrsq_ptr = MAX(*final_relrsq_ptr, qop_resid_arg[isrc][ikappa]->final_rel);
      *final_restart =    MAX(*final_restart,  qop_resid_arg[isrc][ikappa]->final_restart);
      iters += qop_resid_arg[isrc][ikappa]->final_iter;
#ifdef CG_DEBUG
      if(nsrc > 1 || nkappa[isrc] > 1)
	node0_printf("CONGRAD5(src %d,kappa %d): iters = %d resid = %e relresid = %e\n",
	       isrc, ikappa,
	       qop_resid_arg[isrc][ikappa]->final_iter,
	       sqrt(qop_resid_arg[isrc][ikappa]->final_rsq),
 	       sqrt(qop_resid_arg[isrc][ikappa]->final_rel));
#endif
    }

  /* Free QOP fields  */

  destroy_qop_wilson_fermion_links( qop_links );

  for(isrc = 0; isrc < nsrc; isrc++){
    QOP_destroy_D(qop_src[isrc]);    
    qop_src[isrc] = NULL;
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      QOP_destroy_D(qop_sol[isrc][ikappa]);     
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
  destroy_qop_resid_arg(qop_resid_arg, nsrc, nkappa);
  qop_resid_arg = NULL;

  return iters;
}

/********************************************************************/
/* Inverter interface for specific precision                        */
/********************************************************************/
int 
BICGILU_MILC2QOP( int prop_type, wilson_vector *milc_src, 
		  wilson_vector *milc_sol, 
		  quark_invert_control *qic, void *dmp)
{
  int iterations_used;

  /* Since the wilson QOP anticipates multiple sources with multiple
     masses for each, we humor it by setting up its source and sink
     arguments that way.  But for now we support only one source
     with possibly multiple masses (kappas) for that source */

  int nsrc = 1;
  int nkappa[1] = { 1 };

  /* Provision for multiple sources */
  wilson_vector *milc_srcs[1] = { milc_src };

  /* Provision for a solution for each kappa for each source */
  wilson_vector *milc_sols0[1] = { milc_sol };
  wilson_vector **milc_sols[1] = { milc_sols0 };

  /* Provision for separate propagator parameters (kappas) for each
     mass (one source) */
  void *dmps[1] = { dmp };

  /* Performance results */
  Real final_rsq_val;           
  Real final_relrsq_val;           
  int final_restart;

  iterations_used = 
    bicgilu_cl_qop( prop_type, nsrc, nkappa, qic, dmps, milc_srcs, 
		    milc_sols, &final_restart, &final_rsq_val,
		    &final_relrsq_val, EVENANDODD );


  qic->size_r = 0;  /* We don't see the cumulative resid with QOP */
  qic->size_relr = 0;
  qic->final_rsq = final_rsq_val;
  qic->final_relrsq = final_relrsq_val;
  qic->final_iters = iterations_used;
  qic->final_restart = final_restart;
  return  iterations_used;
}

