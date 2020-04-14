/******* d_congrad5_fn_gpu.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* GPU version of d_congrad5_fn_milc.c.  Can be compiled together with it. */
// Note that the restart criterion used by QUDA is different from 
// the restart criterion used by the MILC code

#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#define LOOPEND
#include "../include/loopend.h"
#include <string.h>
#include "../include/generic.h"

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/* Backward compatibility*/
#ifdef SINGLE_FOR_DOUBLE
#define HALF_MIXED
#endif

/********************************************************************/
/* Solution of the normal equations for a single site parity        */
/********************************************************************/

// Note that the restart criterion used by QUDA is different from 
// the restart criterion used by the MILC code
int ks_congrad_parity_gpu(su3_vector *t_src, su3_vector *t_dest, 
			  quark_invert_control *qic, Real mass,
			  imp_ferm_links_t *fn)
{

  char myname[] = "ks_congrad_parity_gpu";
  QudaInvertArgs_t inv_args;
  int i;
  double dtimec = -dclock();
#ifdef CGTIME
  double nflop = 1187;
#endif

//  if(qic->relresid != 0.){
//    printf("%s: GPU code does not yet support a Fermilab-type relative residual\n",myname);
//    terminate(1);
//  }
 
  /* Initialize qic */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq = 0.;
  qic->final_relrsq = 0.;

  /* Compute source norm */
  double source_norm = 0.0;
  FORSOMEFIELDPARITY(i,qic->parity){
    source_norm += (double)magsq_su3vec( &t_src[i] );
  } END_LOOP;
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY(i,qic->parity){
      memset(t_dest + i, 0, sizeof(su3_vector));
    } END_LOOP;

    dtimec += dclock();
#ifdef CGTIME
    if(this_node==0){
      printf("CONGRAD5: time = %e (fn_QUDA %s) masses = 1 iters = %d mflops = %e\n",
	     dtimec, prec_label[qic->prec-1], qic->final_iters, 
	     (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
      fflush(stdout);}
#endif
    
    return 0;
  }

  /* Initialize QUDA parameters */

  initialize_quda();
 
  if(qic->parity == EVEN){
	  inv_args.evenodd = QUDA_EVEN_PARITY;
  }else if(qic->parity == ODD){
	  inv_args.evenodd = QUDA_ODD_PARITY;
  }else{
    printf("%s: Unrecognised parity\n",myname);
    terminate(2);
  }

  inv_args.max_iter = qic->max*qic->nrestart;
#if defined(MAX_MIXED)
  inv_args.mixed_precision = 2;
#elif defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif

  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);
  const int quda_precision = qic->prec;

  double residual, relative_residual;
  int num_iters = 0;

  // for newer versions of QUDA we need to invalidate the gauge field if the links are new
  if ( fn != get_fn_last() || fresh_fn_links(fn) ){
    cancel_quda_notification(fn);
    set_fn_last(fn);
    num_iters = -1;
    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
  }

  inv_args.naik_epsilon = fn->eps_naik;

#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  qudaInvert(MILC_PRECISION,
	     quda_precision, 
	     mass,
	     inv_args,
	     qic->resid,
	     qic->relresid,
	     fatlink, 
	     longlink,
	     t_src, 
	     t_dest,
	     &residual,
	     &relative_residual, 
	     &num_iters);
  

  qic->final_rsq = residual*residual;
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0; 

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn_QUDA %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[quda_precision-1], qic->final_iters, 
	   (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  return num_iters;
}

int ks_congrad_block_parity_gpu(int nsrc, su3_vector **t_src, su3_vector **t_dest, 
				quark_invert_control *qic, Real mass,
				imp_ferm_links_t *fn)
{

  char myname[] = "ks_congrad_block_parity_gpu";
  QudaInvertArgs_t inv_args;
  int i;
  double dtimec = -dclock();
#ifdef CGTIME
  double nflop = 1187;  // FIXME Wrong flops
#endif

  /* Initialize qic */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq = 0.;
  qic->final_relrsq = 0.;

  /* Compute source norm */
  double source_norm = 0.0;
  FORSOMEFIELDPARITY(i,qic->parity){
    source_norm += (double)magsq_su3vec( &t_src[0][i] );
  } END_LOOP;
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY(i,qic->parity){
      memset(t_dest + i, 0, sizeof(su3_vector));
    } END_LOOP;
    
    dtimec += dclock();
#ifdef CGTIME
    if(this_node==0){
      printf("CONGRAD5: time = %e (fn_QUDA %s) masses = 1 iters = %d mflops = %e\n",
	     dtimec, prec_label[MILC_PRECISION-1], qic->final_iters,
	     (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
      fflush(stdout);}
#endif
    
    return 0;
  }
  
  /* Initialize QUDA parameters */

  initialize_quda();

  if(qic->parity == EVEN){
          inv_args.evenodd = QUDA_EVEN_PARITY;
  }else if(qic->parity == ODD){
          inv_args.evenodd = QUDA_ODD_PARITY;
  }else{
    printf("%s: Unrecognised parity\n",myname);
    terminate(2);
  }

  inv_args.max_iter = qic->max*qic->nrestart;
#if defined(MAX_MIXED) || defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif

  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);
  const int quda_precision = qic->prec;

  double residual, relative_residual;
  int num_iters;

  // for newer versions of QUDA we need to invalidate the gauge field if the links are new
  static imp_ferm_links_t *fn_last = NULL;
  if ( fn != fn_last || fresh_fn_links(fn) ){
    cancel_quda_notification(fn);
    fn_last = fn;
    num_iters = -1;
    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
  }

  inv_args.naik_epsilon = fn->eps_naik;

#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  qudaInvertMsrc(MILC_PRECISION,
                 quda_precision,
                 mass,
                 inv_args,
                 qic->resid,
                 qic->relresid,
                 fatlink,
                 longlink,
                 (void**)t_src,
		 (void**)t_dest,
                 &residual,
                 &relative_residual,
                 &num_iters,
		 nsrc);


  qic->final_rsq = residual*residual;
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0;

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn_QUDA %s) masses = 1 iters = %d mflops = %e\n",
           dtimec, prec_label[quda_precision-1], qic->final_iters,
           (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  return num_iters;
}

/********************************************************************/
/* Multigrid solution of the full Dirac equation for both parities  */
/********************************************************************/

int mat_invert_mg_field_gpu(su3_vector *src, su3_vector *dst, 
			    quark_invert_control *qic,
			    Real mass, imp_ferm_links_t *fn){

  char myname[] = "mat_invertt_mg_field_gpu";
  QudaInvertArgs_t inv_args;
  int i;
  double dtimec = -dclock();
#ifdef CGTIME
  double nflop = 1187;
#endif

  /* Initialize qic */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq = 0.;
  qic->final_relrsq = 0.;

  /* Compute source norm */
  double source_norm = 0.0;
  FORSOMEFIELDPARITY(i,qic->parity){
    source_norm += (double)magsq_su3vec( &t_src[i] );
  } END_LOOP;
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("mat_invert_mg_field_gpu: source_norm = %e\n", (double)source_norm);
#endif

  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY(i,qic->parity){
      memset(t_dest + i, 0, sizeof(su3_vector));
    } END_LOOP;

    dtimec += dclock();
#ifdef CGTIME
    if(this_node==0){
      printf("CONGRAD5: time = %e (fn_QUDA_MG %s) masses = 1 iters = %d mflops = %e\n",
	     dtimec, prec_label[qic->prec-1], qic->final_iters, 
	     (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
      fflush(stdout);}
#endif
    
    return 0;
  }

  /* Initialize QUDA parameters */

  initialize_quda();

  inv_args.max_iter = qic->max*qic->nrestart;
#if defined(MAX_MIXED)
  inv_args.mixed_precision = 2;
#elif defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif

  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);
  const int quda_precision = qic->prec;

  // TBD: we should persist this.
  static void *mg_preconditioner = NULL;

  double residual, relative_residual;
  int num_iters = 0;

  inv_args.naik_epsilon = fn->eps_naik;

#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  // For now staggered MG only supports fields of both parities.
  // How do we check this?


  // for newer versions of QUDA we need to invalidate the gauge field if the links are new
  if ( fn != get_fn_last() || fresh_fn_links(fn) ){

    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
    cancel_quda_notification(fn);
    set_fn_last(fn);
    num_iters = -1;

    node0_printf("%s: setting up the MG inverter\n", myname);
    /* Set up the MG inverter when the links change */
    mg_preconditioner = qudaSetupMultigrid(MILC_PRECISION,
                              quda_precision,
                              mass,
                              inv_args,
                              milc_fatlink,
                              milc_longlink,
                              qic->mgparamfile);
  }

  
  // Just BiCGstab for now
  qudaInvertMG(MILC_PRECISION,
	       quda_precision, 
	       mass,
	       inv_args,
	       qic->resid,
	       qic->relresid,
	       fatlink, 
	       longlink,
	       mg_preconditioner,
	       t_src, 
	       t_dest,
	       &residual,
	       &relative_residual, 
	       &num_iters);
  
  qic->final_rsq = residual*residual;
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0; 

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn_QUDA_MG %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[quda_precision-1], qic->final_iters, 
	   (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

    // TBD: this object should be persistent between solves.
  qudaCleanupMultigrid(mg_preconditioner);

  return num_iters;
}

