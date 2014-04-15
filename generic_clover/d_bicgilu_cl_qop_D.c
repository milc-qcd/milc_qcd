/******************** d_bicgilu_cl_qop_D.c ****************************/
/* MILC Version 7 */

/* For specifically double precision inversions using QOP routines */

/* 4/29/07 C. DeTar */

#undef QOP_PrecisionInt
/* QOP precision in this file is double unless explicitly declared single */
#define QOP_PrecisionInt 2

/* Backward compatibility*/
#ifdef SINGLE_FOR_DOUBLE
#define HALF_MIXED
#endif

/* QOP supports only "single for double" */
#if ! defined(HALF_MIXED) && ! defined(MAX_MIXED)

/********************************************************************/
/* This is the standard double-precision interface */

#include "d_bicgilu_cl_qop_P.c"

/********************************************************************/

#else

/********************************************************************/
/* With this version we do double precision inversion in single precision */
/* Good only for one kappa per source */
/* 6/15/09 C. DeTar */

/* This procedure requires qop_qdp, since we don't yet have generic
   QOP support for the necessary linear algebra */
#include <qop_qdp.h>
#include "generic_clover_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_clover_qop.h"
#include "../include/generic_qopqdp.h"
#include <string.h>

/* 6/15/09 C. DeTar */

/*
 * $Log:
 */

// #ifdef CGTIME
// static const char *qop_prec[2] = {"F", "D"};
// #endif

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_clover/d_bicgilu_cl_qop_D.c,v 1.6 2013/12/26 17:42:48 detar Exp $";

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
  dsu3_matrix **raw_links;
  double *raw_clov;
  QOP_FermionLinksWilson *qop_links;

  /* Construct raw QOP clover term from MILC clover term */
  if(clov == 0){
    raw_clov = NULL;
  }
  else{
    raw_clov = (double *)malloc(72*sites_on_node*sizeof(double));
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

  raw_links = create_raw4_D_G_from_site(F_OFFSET(link), EVENANDODD);
  if(raw_links == NULL)terminate(1);

  /* Map QOP/QDP raw to QOP/QDP structure */

  qop_links = QOP_wilson_create_L_from_raw((double **)raw_links, raw_clov, 
					   QOP_EVENODD);
  
  destroy_raw4_D_G(raw_links); raw_links = NULL;
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
/* Convert MILC parity to QDP subset                                */
/********************************************************************/
static QDP_Subset
milc2qdp_subset(int milc_parity){
  switch(milc_parity){
  case(EVEN):       return QDP_even;
  case(ODD ):       return QDP_odd;
  case(EVENANDODD): return QDP_all;
  default:
    printf("milc2qdp_subset: Bad MILC parity %d\n", milc_parity);
    terminate(1);
  }
  return NULL;
}

/********************************************************************/
/* Load inversion args for Level 3 inverter                         */
/********************************************************************/

static void 
set_qop_invert_arg_norestart( QOP_invert_arg_t* qop_invert_arg, 
			      quark_invert_control *qic, int milc_parity )
{
  qop_invert_arg->max_iter     = qic->max;
  qop_invert_arg->restart      = qic->max;
  qop_invert_arg->max_restarts = 1;
  qop_invert_arg->evenodd      = milc2qop_parity(milc_parity);
}

/********************************************************************/
/* Load residual values for stopping                                */
/********************************************************************/

static QOP_resid_arg_t ***
create_qop_resid_arg( int nsrc, int nmass[], Real min_resid_sq, Real min_rel_sq )
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
      res_arg[isrc][imass]->relmin = min_rel_sq;
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

#define MAXSRC 20

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

/***********************************************************/
/* Create a zero QOP DiracFermion field                    */
/***********************************************************/
static QOP_F3_DiracFermion *
create_qop_DiracFermion_F(){
  QOP_F3_DiracFermion *a;
  QDP_F3_DiracFermion *qdp_a;
  
  qdp_a = QDP_F3_create_D();
  QDP_F3_D_eq_zero(qdp_a, QDP_all);
  a = QOP_F3_convert_D_from_qdp(qdp_a);
  return a;
}

/****************************************************************************/
/* Update double-precision solution vector from single-precision correction */
/****************************************************************************/

static void
update_qop_solution( QOP_DiracFermion **qop_sol_D[], 
		     QLA_D_Real norm_resid[],
		     QOP_F3_DiracFermion **qop_add_F[], 
		     int nsrc, QDP_Subset subset){
  QDP_D3_DiracFermion *qdp_sol_D;
  QDP_F3_DiracFermion *qdp_add_F;
  QDP_D3_DiracFermion *qdp_add_D;
  int i;

  qdp_add_D = QDP_D3_create_D();
  for(i = 0; i < nsrc; i++){
    qdp_sol_D = QOP_D3_convert_D_to_qdp(qop_sol_D[i][0]);
    qdp_add_F = QOP_F3_convert_D_to_qdp(qop_add_F[i][0]);
    QDP_DF3_D_eq_D(qdp_add_D, qdp_add_F, subset);
    /* Rescale the correction */
    QDP_D3_D_eq_r_times_D(qdp_add_D, norm_resid+i, qdp_add_D, subset);
    QDP_D3_D_peq_D(qdp_sol_D, qdp_add_D, subset);
    /* Restore QOP fields */
    qop_add_F[i][0] = QOP_F3_convert_D_from_qdp(qdp_add_F);
    qop_sol_D[i][0] = QOP_D3_convert_D_from_qdp(qdp_sol_D);
  }
  QDP_D3_destroy_D(qdp_add_D);
}


/******************************************/
/* Compute residual vectors from solution */
/******************************************/

static void
compute_qdp_residuals( int prop_type,
		       QDP_D3_DiracFermion *qdp_resid[], 
		       QDP_D3_DiracFermion *qdp_rhs[],  
		       QOP_FermionLinksWilson *qop_links, 
		       QOP_DiracFermion **qop_sol[], 
		       void *dmps[],
		       float *kappas[], int nkappa[],
		       int nsrc, int milc_parity ){
  QOP_info_t info = {0., 0., 0, 0, 0};
  QDP_D3_DiracFermion *qdp_prod = QDP_D3_create_D();
  QOP_DiracFermion *qop_prod;
  QOP_evenodd_t qop_parity = milc2qop_parity(milc_parity);
  QDP_Subset subset = milc2qdp_subset(milc_parity);

  int i;

  for(i = 0; i < nsrc; i++){
    qop_prod = QOP_D3_convert_D_from_qdp(qdp_prod);
    if(prop_type == CLOVER_TYPE)
      QOP_wilson_dslash(&info, qop_links, kappas[i][0], +1,
			qop_prod, qop_sol[i][0], qop_parity, qop_parity);
    else {/* IFLA_TYPE */
      QOP_wilson_ifla_coeffs_t dcof;
      /* If there are multiple kappas, we assume that the parameters
	 are the same for all kappas in this solve */
      load_qop_wilson_ifla_coeffs(&dcof, (newaction_ifla_param *)dmps[i] );
      QOP_wilson_ifla_dslash(&info, qop_links, kappas[i][0], +1,
		     &dcof, qop_prod, qop_sol[i][0], qop_parity, qop_parity);
    }
    qdp_prod = QOP_D3_convert_D_to_qdp( qop_prod );
    
    QDP_D3_D_eq_D_minus_D( qdp_resid[i], qdp_rhs[i], qdp_prod, subset );
  }
  QDP_D3_destroy_D(qdp_prod);
}

/**********************************************************************/
/* FNAL relative residual norm                                        */
/**********************************************************************/

/* This algorithm is rather clumsy. */

static Real 
qdp_relative_residue(QDP_D3_DiracFermion *p, QDP_D3_DiracFermion *q, 
		     QDP_Subset subset)
{
  QDP_D_Real *ratio, *num, *den, *ones;
  QDP_Int *ok;
  QLA_D_Real one = 1., residue;

  num = QDP_D_create_R();
  den = QDP_D_create_R();
  ones = QDP_D_create_R();
  QDP_D_R_eq_r(ones, &one, subset);
  ratio = QDP_D_create_R();
  ok = QDP_create_I();
  QDP_D_R_eq_zero(ratio, subset);
  
  QDP_D3_R_eq_norm2_D(num, p, subset);
  QDP_D3_R_eq_norm2_D(den, q, subset);

  /* Set mask to 1 if the denominator element is zero */
  QDP_D_I_eq_R_eq_R(ok, ratio, den, subset);
  /* Replace any zeros with ones in denominator to prevent division by zero */
  QDP_D_R_eq_R_mask_I(den, ones, ok, subset);
  /* Now take ratios */
  QDP_D_R_eq_R_divide_R(ratio, num, den, subset);
  /* Replace ratios with 1 where original denominator was zero */
  QDP_D_R_eq_R_mask_I(ratio, ones, ok, subset);
  /* Total the ratios */
  QDP_D_r_eq_sum_R(&residue, ratio, subset);

  QDP_D_destroy_R(num);
  QDP_D_destroy_R(den);
  QDP_D_destroy_R(ratio);
  QDP_D_destroy_R(ones);
  QDP_destroy_I(ok);

  /* Normalize.  Only choices are full volume and half-volume  */
  if(subset == QDP_all)
    return residue/QDP_volume();
  else
    return 2*residue/QDP_volume();
}

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

/**********************************************************************/
/* Driver for doing double precision inversion using single precision */
/**********************************************************************/

int 
bicgilu_cl_qop_single_for_double( int prop_type,
				  QOP_FermionLinksWilson *qop_links, 
				  quark_invert_control *qic, int milc_parity,
				  void *dmps[],
				  float *kappas[], int nkappa[], 
				  QOP_DiracFermion **qop_sol[], 
				  QOP_DiracFermion *qop_src[], 
				  int nsrc,		    
				  int *final_restart,
				  Real *final_rsq_ptr )
{
  int i, iters, iters_F = 0;
  int converged;
  int nrestart;
  int max_restarts = qic->nrestart;
  int isrc, ikappa;
  int final_restart_F;
  Real final_rsq_F, final_relrsq_F;
  Real resid_F = 3e-7;   /* The limits of a single precision inversion */
  Real rel_F = 0;   /* The limits of a single precision inversion */
  QOP_invert_arg_t qop_invert_arg;
  QOP_resid_arg_t  ***qop_resid_arg_F;
  QOP_info_t info_F = {0., 0., 0, 0, 0}, info = {0., 0., 0, 0, 0};
  QDP_Subset subset = milc2qdp_subset(milc_parity);
  QOP_F3_FermionLinksWilson *qop_links_F;
  QOP_F3_DiracFermion **qop_sol_F[MAXSRC], *qop_rhs_F[MAXSRC];
  QDP_F3_DiracFermion *qdp_rhs_F[MAXSRC];
  QDP_D3_DiracFermion *qdp_src[MAXSRC], *qdp_resid[MAXSRC];
  QDP_D3_DiracFermion *qdp_sol;
  Real relresid2[MAXSRC];
  Real resid2[MAXSRC];
  QLA_D_Real norm2_src[MAXSRC], norm2_resid[MAXSRC], norm_resid[MAXSRC], scale_resid;
  char myname[] = "bicgilu_cl_qop_single_for_double";
  
  /* Only one kappa allowed per source for this algorithm */
  for(i = 0; i < nsrc; i++){
    if(nkappa[i] > 1){
      printf("%s: nkappa[%d] = %d != 1\n",myname,i,nkappa[i]);
      terminate(1);
    }
  }
  
  /* Set qop_invert_arg */
  /* We don't do restarts for the single precision step */
  /* We interpret "qic->nrestart" to mean the max number of calls to
     the single-precision inverter */
  set_qop_invert_arg_norestart( & qop_invert_arg, qic, milc_parity );
  
  /* Pointers for residual errors */
  /* For now we set the residual to something sensible for single precision */
  qop_resid_arg_F = create_qop_resid_arg( nsrc, nkappa, resid_F*resid_F, rel_F*rel_F);

  /* Create a single precision copy of the links object */
  qop_links_F = QOP_FD3_wilson_create_L_from_L( qop_links );

  /* Take norm of source and create temporaries */

  for(i = 0; i < nsrc; i++){
    qdp_src[i] = QOP_D3_convert_D_to_qdp( qop_src[i] );
    QDP_D3_r_eq_norm2_D( norm2_src+i, qdp_src[i], subset );
    qdp_resid[i] = QDP_D3_create_D();
    qdp_rhs_F[i] = QDP_F3_create_D();
    qop_sol_F[i] = (QOP_F3_DiracFermion **)malloc(sizeof(QOP_F3_DiracFermion *));
  }


  /* Main loop */

  nrestart = 0;
  converged = 0;
  iters = 0;

  info.final_sec = -dclock();
  info.final_flop = 0;
  info.status = QOP_SUCCESS;

  while(1){
    /* Create new residual vectors from the result */
    /* r = src - A sol */
    compute_qdp_residuals( prop_type, qdp_resid, qdp_src, 
			   qop_links, qop_sol, dmps, kappas, 
			   nkappa, nsrc, milc_parity );

    /* Compute two different norms */
    qic->final_rsq = 0;
    qic->final_relrsq = 0;
    for(i = 0; i < nsrc; i++){
      qdp_sol = QOP_convert_D_to_qdp( qop_sol[i][0] );
      relresid2[i] = qdp_relative_residue( qdp_resid[i], qdp_sol, subset );
      qop_sol[i][0] = QOP_convert_D_from_qdp( qdp_sol );
      qic->final_relrsq = (relresid2[i] > qic->final_relrsq) ? relresid2[i] : qic->final_relrsq;

      QDP_D3_r_eq_norm2_D( norm2_resid+i, qdp_resid[i], subset );
      resid2[i] = norm2_resid[i]/norm2_src[i];
      qic->final_rsq = (resid2[i] > qic->final_rsq) ? resid2[i] : qic->final_rsq;
#ifdef CG_DEBUG
      node0_printf("%s: double precision restart %d resid2 = %.2e vs %.2e relresid2 = %.2e vs %.2e\n",
		   myname, nrestart, resid2[i], qic->resid * qic->resid, relresid2[i], 
		   qic->relresid * qic->relresid );
#endif
    }
    *final_rsq_ptr = qic->final_rsq;  /* Use Cartesian norm for now */
    *final_restart = nrestart;

    /* Stop when converged */
    converged = 1;
    for(i = 0; i < nsrc; i++){
      if((qic->resid > 0 && resid2[i] > qic->resid * qic->resid) || 
	 (qic->relresid > 0 && relresid2[i] > qic->relresid * qic->relresid)){
	converged = 0;
	break;
      }
    }

    if(converged || nrestart++>=max_restarts)break;

    for(i = 0; i < nsrc; i++){
      /* Scale the RHS to avoid underflow */
      norm_resid[i] = sqrt(norm2_resid[i]);
      scale_resid = 1./norm_resid[i];
      QDP_D3_D_eq_r_times_D(qdp_resid[i], &scale_resid, qdp_resid[i], subset);
      /* Scaled residual becomes the new source */
      QDP_FD3_D_eq_D( qdp_rhs_F[i], qdp_resid[i],  subset);
      qop_rhs_F[i] = QOP_F3_convert_D_from_qdp( qdp_rhs_F[i]);
      /* Prepare to solve in single precision by creating a single
	 precision copy of the source.  Set the trial solution to zero. */
      qop_sol_F[i][0] = create_qop_DiracFermion_F();
    }


    /* Solve in single precision */
    double dtime = -dclock();
    info_F.final_flop = 0.;
    bicgilu_cl_qop_generic_F( prop_type, &info_F, qop_links_F, 
	  &qop_invert_arg, qop_resid_arg_F, dmps, nkappa, qop_sol_F, 
  	  qop_rhs_F, nsrc);
    dtime += dclock();

    /* Report performance statistics */
    
    /* For now we return the largest value and total iterations */
    final_rsq_F = 0;
    final_relrsq_F = 0;
    final_restart_F = 0;
    iters_F = 0;
    for(isrc = 0; isrc < nsrc; isrc++)
      for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
	/* QOP routines return the ratios of the squared norms */
	final_rsq_F =    MAX(final_rsq_F, qop_resid_arg_F[isrc][ikappa]->final_rsq);
	final_relrsq_F = MAX(final_relrsq_F, qop_resid_arg_F[isrc][ikappa]->final_rel);
	final_restart_F =    MAX(final_restart_F,  qop_resid_arg_F[isrc][ikappa]->final_restart);
	iters_F += qop_resid_arg_F[isrc][ikappa]->final_iter;
	if(nsrc > 1 || nkappa[isrc] > 1)
	  node0_printf("BICG(src %d,kappa %d): iters = %d resid = %e relresid = %e\n",
		       isrc, ikappa,
		       qop_resid_arg_F[isrc][ikappa]->final_iter,
		       sqrt(qop_resid_arg_F[isrc][ikappa]->final_rsq),
		       sqrt(qop_resid_arg_F[isrc][ikappa]->final_rel));
      }
    
#ifdef CGTIME
    node0_printf("%s: single precision iters = %d status %d final_rsq %.2e wanted %2e final_rel %.2e wanted %.2e\n",
		 myname, iters_F, info_F.status, final_rsq_F, resid_F * resid_F, final_relrsq_F, rel_F);
    node0_printf("time = %g flops = %e mflops = %g\n", dtime, info_F.final_flop, 
		 info_F.final_flop/(1.0e6*dtime) );
    fflush(stdout);
#endif

    /* Add single-precision result to double precision solution (with rescaling) */
    update_qop_solution( qop_sol, norm_resid, qop_sol_F, nsrc, subset );

    for(i = 0; i < nsrc; i++){
      QOP_F3_destroy_D(qop_sol_F[i][0]);
      /* Convert back */
      qdp_rhs_F[i] = QOP_F3_convert_D_to_qdp(qop_rhs_F[i]);
    }

    info.final_flop += info_F.final_flop;
    iters += iters_F;
  }

  /* Clean up */

  for(i = 0; i < nsrc; i++){
    QDP_F3_destroy_D( qdp_rhs_F[i] );
    QDP_D3_destroy_D( qdp_resid[i] );
    /* Must restore qop_src in case the caller reuses it */
    qop_src[i] = QOP_D3_convert_D_from_qdp( qdp_src[i] );
    free(qop_sol_F[i]);
  }

  QOP_F3_wilson_destroy_L( qop_links_F );
  destroy_qop_resid_arg(qop_resid_arg_F, nsrc, nkappa);
  qop_resid_arg_F = NULL;

  if(!converged){
    node0_printf("%s: NOT Converged after %d iters and %d restarts\n",
		 myname, iters, nrestart);
  }

  info.final_sec += dclock();
#ifdef CGTIME
  node0_printf("CGTIME: time = %e (wilson_qop FD) ", info.final_sec);
  for(isrc = 0; isrc < nsrc; isrc++)
    node0_printf("nkappa[%d] = %d tot_iters = %d ",
		 isrc,nkappa[isrc],iters);
  node0_printf("mflops = %e\n", info.final_flop/(1.0e6*info.final_sec) );
  fflush(stdout);
#endif

  return iters;
}

#define MAXSRC 20

/******************************************************************************/
/* Map MILC fields to QOP format via QDP.  Call QOP single-precision driver   */
/******************************************************************************/

static int
bicgilu_cl_qop_qdp(int prop_type, int nsrc, int nkappa[], 
		   quark_invert_control *qic,
		   void *dmps[], wilson_vector *milc_srcs[], 
		   wilson_vector **milc_sols[],
		   int *final_restart,
		   Real* final_rsq_ptr, int milc_parity )
{
  int isrc, ikappa;
  QOP_FermionLinksWilson *qop_links;
  QOP_DiracFermion **qop_sol[MAXSRC], *qop_src[MAXSRC];
  int iterations_used = 0;
  double remaptime;
  int i;
  site *s;
  static float t_kappa;
  float *kappas[1] = { &t_kappa };
  char myname[] = "bicgilu_cl_qop_qdp";
  
  /* Check dimension */
  if(nsrc > MAXSRC){
    printf("%s: too many sources\n",myname);
    terminate(1);
  }
  
  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s: Error initializing QOP\n",myname);
    terminate(1);
  }
  
  /* Create QOP links object (double precision) and extract kappa */

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
  remaptime = -dclock(); 

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_DiracFermion **)malloc(sizeof(QOP_DiracFermion *)*nkappa[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("bicgilu_cl_qop_qdp: Can't allocate qop_sol\n");
      terminate(1);
    }
  }

  /* Map MILC source and sink to double-precision QOP fields */
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
  iterations_used = bicgilu_cl_qop_single_for_double( prop_type,
    qop_links, qic,  milc_parity, dmps, kappas, nkappa, qop_sol, qop_src, nsrc, 
    final_restart, final_rsq_ptr );
  remaptime -= dclock();
  
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

  return iterations_used;
}

/********************************************************************/
/* Inverter interface for double precision                        */
/********************************************************************/
int 
bicgilu_cl_milc2qop_D( int prop_type, wilson_vector *milc_src, 
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
  int final_restart;

  iterations_used = 
    bicgilu_cl_qop_qdp( prop_type, nsrc, nkappa, qic, dmps, milc_srcs,
			milc_sols, &final_restart, &final_rsq_val,
			EVENANDODD );

  qic->size_r = 0;  /* We don't see the cumulative resid with QOP*/
  qic->size_relr = 0;
  qic->final_rsq = final_rsq_val;
  qic->final_iters = iterations_used;
  qic->final_restart = final_restart;
  return  iterations_used;
}

#endif
