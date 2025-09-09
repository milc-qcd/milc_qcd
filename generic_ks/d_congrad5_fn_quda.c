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

  // node0_printf("Entered %s\n", myname);

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
    num_iters = -1;  /* This hack signals QUDA */
    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
  }

  inv_args.naik_epsilon = fn->eps_naik;

#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  // Setup for deflation (and eigensolve) on GPU
  int parity = qic->parity;
  int blockSize = param.eigen_param.blockSize;
  static Real previous_mass = -1.0;
  static bool first_solve=true;

  QudaEigensolverArgs_t eig_args;
  eig_args.struct_size = 1192; // Could also use sizeof(QudaEigensolverArgs_t) to automagically update, but using a static number will catch the case when the struct is updated by QUDA but MILC is not updated
  eig_args.block_size = blockSize;
  eig_args.n_conv = (param.eigen_param.Nvecs_in > param.eigen_param.Nvecs) ? param.eigen_param.Nvecs_in : param.eigen_param.Nvecs;
  eig_args.n_ev_deflate = ( parity == EVEN && qic->deflate ) ? param.eigen_param.Nvecs : 0; // Only deflate even solves for now
  eig_args.n_ev = eig_args.n_conv;
  eig_args.n_kr = (param.eigen_param.Nkr < eig_args.n_ev ) ? 2*eig_args.n_ev : param.eigen_param.Nkr;
  eig_args.tol = param.eigen_param.tol;
  eig_args.max_restarts = param.eigen_param.MaxIter;
  eig_args.poly_deg = param.eigen_param.poly.norder;
  eig_args.a_min = param.eigen_param.poly.minE;
  eig_args.a_max = param.eigen_param.poly.maxE;
  strcpy( eig_args.vec_infile, param.ks_eigen_startfile );
  strcpy( eig_args.vec_outfile, param.ks_eigen_savefile );
  eig_args.vec_in_parity = QUDA_EVEN_PARITY; // TODO: Update when we add support for odd parity eigenvector files
  eig_args.preserve_evals = ( first_solve || fabs(mass - previous_mass) < 1e-6 ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.batched_rotate = param.eigen_param.batchedRotate;
  eig_args.save_prec = QUDA_SINGLE_PRECISION; // add to input parameters?
  eig_args.partfile = param.eigen_param.partfile ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.io_parity_inflate = QUDA_BOOLEAN_FALSE;
  eig_args.use_norm_op = ( parity == EVENANDODD ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.use_pc = ( parity != EVENANDODD) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.tol_restart = param.eigen_param.tol_restart;
  eig_args.eig_type = ( eig_args.block_size > 1 ) ? QUDA_EIG_BLK_TR_LANCZOS : QUDA_EIG_TR_LANCZOS;  /* or QUDA_EIG_IR_ARNOLDI, QUDA_EIG_BLK_IR_ARNOLDI */
  eig_args.spectrum = QUDA_SPECTRUM_SR_EIG; /* Smallest Real. Other options: LM, SM, LR, SR, LI, SI */
  eig_args.qr_tol = eig_args.tol;
  eig_args.require_convergence = QUDA_BOOLEAN_TRUE;
  eig_args.check_interval = 10;
  eig_args.use_dagger = QUDA_BOOLEAN_FALSE;
  eig_args.compute_gamma5 = QUDA_BOOLEAN_FALSE;
  eig_args.compute_svd = QUDA_BOOLEAN_FALSE;
  eig_args.use_eigen_qr = QUDA_BOOLEAN_TRUE;
  eig_args.use_poly_acc = QUDA_BOOLEAN_TRUE;
  eig_args.arpack_check = QUDA_BOOLEAN_FALSE;
  eig_args.compute_evals_batch_size = 16;
  eig_args.preserve_deflation = QUDA_BOOLEAN_TRUE;
  
  if(param.eigen_param.eigPrec == 2) {
    eig_args.prec_eigensolver = QUDA_DOUBLE_PRECISION;
  } else if(param.eigen_param.eigPrec == 1) {
    eig_args.prec_eigensolver = QUDA_SINGLE_PRECISION;
  } else if(param.eigen_param.eigPrec == 0) {
    eig_args.prec_eigensolver = QUDA_HALF_PRECISION;
  } else {
    printf("%s: Unrecognized eigensolver precision\n",myname);
    terminate(2);
  }

  previous_mass = mass;
  first_solve = false;

  qudaInvertDeflatable(MILC_PRECISION,
	     quda_precision, 
	     mass,
	     inv_args,
	     eig_args,
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

    node0_printf("Calling check_invert_field2\n"); fflush(stdout);
    check_invert_field2(t_src, t_dest, mass, 2e-5, fn, qic->parity);
  
  return num_iters;
}

/********************************************************************/
/* Solution of the normal equations for a single site parity with   */
/* multiple right sides                                             */
/********************************************************************/

int ks_congrad_block_parity_gpu(int nsrc, su3_vector **t_src, su3_vector **t_dest, 
				quark_invert_control *qic, Real mass,
				imp_ferm_links_t *fn)
{

  char myname[] = "ks_congrad_block_parity_gpu";

#if 0
  /* Debug: Solve separately, rather than batch */
  int num_iters = 0;
  for(int i = 0; i < nsrc; i++){
    num_iters += ks_congrad_parity_gpu(t_src[i], t_dest[i], qic, mass, fn);
  }
  return num_iters;
#else
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
  for(int j = 0; j < nsrc; j++){
    FORSOMEFIELDPARITY(i,qic->parity){
      source_norm += (double)magsq_su3vec( &t_src[j][i] );
    } END_LOOP;
  }
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solution and return zero iterations */
    for(int j = 0; j< nsrc; j++){
      FORSOMEFIELDPARITY(i,qic->parity){
	memset(t_dest[j] + i, 0, sizeof(su3_vector));
      } END_LOOP;
    }
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
    node0_printf("%s: Using QUDA's block solver with EVEN parity %x\n", myname);
  }else if(qic->parity == ODD){
    inv_args.evenodd = QUDA_ODD_PARITY;
    node0_printf("%s: Using QUDA's block solver with ODD parity %x\n", myname);
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

  // Setup for deflation (and eigensolve) on GPU
  int parity = qic->parity;
  int blockSize = param.eigen_param.blockSize;
  static Real previous_mass = -1.0;
  static bool first_solve=true;

  QudaEigensolverArgs_t eig_args;
  eig_args.struct_size = 1192; // Could also use sizeof(QudaEigensolverArgs_t) to automagically update, but using a static number will catch the case when the struct is updated by QUDA but MILC is not updated
  eig_args.block_size = blockSize;
  eig_args.n_conv = (param.eigen_param.Nvecs_in > param.eigen_param.Nvecs) ? param.eigen_param.Nvecs_in : param.eigen_param.Nvecs;
  eig_args.n_ev_deflate = ( parity == EVEN && qic->deflate ) ? param.eigen_param.Nvecs : 0; // Only deflate even solves for now
  eig_args.n_ev = eig_args.n_conv;
  eig_args.n_kr = (param.eigen_param.Nkr < eig_args.n_ev ) ? 2*eig_args.n_ev : param.eigen_param.Nkr;
  eig_args.tol = param.eigen_param.tol;
  eig_args.max_restarts = param.eigen_param.MaxIter;
  eig_args.poly_deg = param.eigen_param.poly.norder;
  eig_args.a_min = param.eigen_param.poly.minE;
  eig_args.a_max = param.eigen_param.poly.maxE;
  strcpy( eig_args.vec_infile, param.ks_eigen_startfile );
  strcpy( eig_args.vec_outfile, param.ks_eigen_savefile );
  eig_args.vec_in_parity = QUDA_EVEN_PARITY; // TODO: Update when we add support for odd parity eigenvector files
  eig_args.preserve_evals = ( first_solve || fabs(mass - previous_mass) < 1e-6 ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.batched_rotate = param.eigen_param.batchedRotate;
  eig_args.save_prec = QUDA_SINGLE_PRECISION; // add to input parameters?
  eig_args.partfile = param.eigen_param.partfile ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.io_parity_inflate = QUDA_BOOLEAN_FALSE;
  eig_args.use_norm_op = ( parity == EVENANDODD ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.use_pc = ( parity != EVENANDODD) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.tol_restart = param.eigen_param.tol_restart;
  eig_args.eig_type = ( eig_args.block_size > 1 ) ? QUDA_EIG_BLK_TR_LANCZOS : QUDA_EIG_TR_LANCZOS;  /* or QUDA_EIG_IR_ARNOLDI, QUDA_EIG_BLK_IR_ARNOLDI */
  eig_args.spectrum = QUDA_SPECTRUM_SR_EIG; /* Smallest Real. Other options: LM, SM, LR, SR, LI, SI */
  eig_args.qr_tol = eig_args.tol;
  eig_args.require_convergence = QUDA_BOOLEAN_TRUE;
  eig_args.check_interval = 10;
  eig_args.use_dagger = QUDA_BOOLEAN_FALSE;
  eig_args.compute_gamma5 = QUDA_BOOLEAN_FALSE;
  eig_args.compute_svd = QUDA_BOOLEAN_FALSE;
  eig_args.use_eigen_qr = QUDA_BOOLEAN_TRUE;
  eig_args.use_poly_acc = QUDA_BOOLEAN_TRUE;
  eig_args.arpack_check = QUDA_BOOLEAN_FALSE;
  eig_args.compute_evals_batch_size = 16;
  eig_args.preserve_deflation = QUDA_BOOLEAN_TRUE;
  
  if(param.eigen_param.eigPrec == 2) {
    eig_args.prec_eigensolver = QUDA_DOUBLE_PRECISION;
  } else if(param.eigen_param.eigPrec == 1) {
    eig_args.prec_eigensolver = QUDA_SINGLE_PRECISION;
  } else if(param.eigen_param.eigPrec == 0) {
    eig_args.prec_eigensolver = QUDA_HALF_PRECISION;
  } else {
    printf("%s: Unrecognized eigensolver precision\n",myname);
    terminate(2);
  }

  previous_mass = mass;
  first_solve = false;

  qudaInvertMsrcDeflatable(MILC_PRECISION,
     quda_precision,
     mass,
     inv_args,
     eig_args,
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


  // MILC's convention impled from d_congrad5_fn_milc.c is that final_rsq, final_relrsq, and final_iters
  // are based on the values from the last solve, which qudaInvertMsrc respects.
  qic->final_rsq = residual * residual;
  qic->final_relrsq = relative_residual * relative_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0;

  dtimec += dclock();

#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn_QUDA %s) masses = 1 srcs = %d iters = %d mflops = %e\n",
           dtimec, prec_label[quda_precision-1], nsrc,qic->final_iters,
           (double)(nflop*nsrc*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  for(int j = 0; j < nsrc; j++){
    node0_printf("Calling check_invert_field2 for case %d\n", j); fflush(stdout);
    check_invert_field2(t_src[j], t_dest[j], mass, 2e-5, fn, qic->parity);
  }

  // On the other hand, MILC expects the returned value to be the aggregate number of iterations
  // performed by each solve if it was performed sequentially. This can be approximated by
  // the number of iterations for a single solve times the number of sources.
  return num_iters * nsrc;

#endif /* if 1 */
  return num_iters;
}


