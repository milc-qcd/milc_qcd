/******* inc_eigcg.c - Incremental eigCG for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  
*/

/* Hiroshi Ohno */
/* About eigCG algorithm, see A. Stathopoulos and K. Orginos [arXiv:0707.0131] */
#include "generic_ks_includes.h"
#include "../include/blas_lapack.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <string.h>

#include <quda_milc_interface.h>

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/* eigCG */
/* Solve a linear equation as well as calculate Nvecs eigenpairs of -Dslash^2 simultaniously. */
/* Lanczos part restarts after first m steps and then every (m - 2*Nvecs) steps */
/* This version looks at the initial vector every "niter" passes. */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations before restarting.
   max_restarts = max number of restarts
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.

   reinitialize after niters iterations and try once more.
*/
int ks_eigCG_parity_gpu(su3_vector *src, su3_vector *dest, double *eigVal, su3_vector **eigVec, int m, int Nvecs, int Nvecs_max, int curr_idx, int last_rhs_flag, quark_invert_control *qic, Real mass, imp_ferm_links_t *fn){

  /*************/
  char myname[] = "ks_eigCG_parity_gpu";

  QudaInvertArgs_t inv_args;
  int i;
  double dtimec;

  /* Unpack structure */
  int niter        = qic->max;                    /* maximum number of iters per restart */
  int max_restarts = qic->nrestart;               /* maximum restarts */
  Real rsqmin      = qic->resid*qic->resid;       /* desired residual - 
						     normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relrsqmin   = qic->relresid*qic->relresid; /* desired relative residual (FNAL)*/
  int parity       = qic->parity;                 /* EVEN, ODD */

  int otherparity = (parity == EVEN) ? ODD : EVEN;

  if(fn == NULL){
    printf("%s(%d): Called with NULL fn\n", myname, this_node);
    terminate(1);
  }

  dtimec = -dclock(); 

  qic->size_r        = 0;
  qic->size_relr     = 1.;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq     = 0.;
  qic->final_relrsq  = 0.;

  /* Source norm */
  source_norm = dzero;
  FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:source_norm)){
    source_norm += (double)magsq_su3vec(src+i);
  } END_LOOP_OMP
  g_doublesum(&source_norm);
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Provision for trivial solution */
  if(source_norm == dzero){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
      memset(dest+i, 0, sizeof(su3_vector));
    } END_LOOP_OMP

    dtimec += dclock();
#ifdef CGTIME
    node0_printf("CONGRAD5_EIGCG: time = %e (fn %s) masses = 1 iters = %d\n",
		 dtimec, prec_label[PRECISION-1], qic->final_iters);
#endif

    return 0;
  }

  initialize_quda();
 
  inv_args.evenodd =  ( parity == EVEN ) ? QUDA_EVEN_PARITY : inv_args.evenodd = QUDA_ODD_PARITY;

  inv_args.max_iter = qic->max*qic->nrestart;
#if defined(MAX_MIXED) || defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif
#ifdef CG_DEBUG
  inv_args.mixed_precision = 0;//enforce full precision
#endif
  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);

  const int quda_precision = qic->prec;
  const int ritz_precision = (inv_args.mixed_precision == 0) ? quda_precision : 1;//1 "means" single precision...

  double residual, relative_residual;
  int iteration; 
  int deflation_grid = Nvecs_max / Nvecs;
  double restart_tol = 5e+3*rsqmin;//hardcoded for a moment

  qudaEigCGInvert(PRECISION,
	          quda_precision, 
	          mass,
	          inv_args,
	          qic->resid,
	          qic->relresid,
	          fatlink, 
	          longlink,
                  u0,//??
	          src, 
	          dest,
                  (void*)eigVec,//array of ritz vectors
                  eigVals,
                  ritz_precision,
                  m, Nvecs,
                  deflation_grid,
                  restart_tol,//e.g.: 5e+3*target_residual
                  curr_idx,//current rhs
                  last_rhs_flag,//is this the last rhs to solve?
	          &residual,
	          &relative_residual, 
	          &iteration);

   qic->size_r        = (Real)residual;
   qic->size_relr     = relative_residual;
   qic->final_iters   = iteration;
   qic->final_restart = 0;
   qic->converged     = 1;

#ifdef CG_DEBUG
   node0_printf("iter=%d, rsq/src= %e, relrsq= %e,\n",
		 iteration, (double)qic->size_r, (double)qic->size_relr);
#endif

  return iteration;
}

/* Incremental eigCG */
/* !!! This computes eigenpairs of 4m^2-Dslash^2 !!! */

int ks_inc_eigCG_parity_gpu( su3_vector *src, su3_vector *dest, double *eigVal,
			 su3_vector **eigVec, eigcg_params *eigcgp, quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn ){

  int i, parity, m, Nvecs, Nvecs_max, iteration, last_rhs;
  static int curr_rhs;

  dtimec = -dclock();

  parity = qic->parity;
  m = eigcgp->m;
  Nvecs = eigcgp->Nvecs;
  //Nvecs_curr = eigcgp->Nvecs_curr;
  Nvecs_max = eigcgp->Nvecs_max;

  last_rhs = eigcgp->last_rhs;//not defined!
  
  /* Solve a linear equation */
  dtimec3 = -dclock();
  iteration = ks_eigCG_parity_gpu(src, dest, eigVal, eigVec, m, Nvecs, Nvecs_max, curr_rhs, last_rhs, qic, mass, fn);

  dtimec += dclock();

#ifdef CGTIME
  node0_printf("INC_EIGCG: total time                 %e\n", dtimec);
#endif

  return iteration;
}
