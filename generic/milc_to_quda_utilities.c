/******************* milc_to_quda_utilities.c ************************/
/* For the QUDA/GPU interface */

//#include <cuda.h>  // DEBUG-JNS
//#include <cuda_runtime.h> // DEBUG-JNS
#include "generic_includes.h"
#include "../include/generic_quda.h"
#include <string.h>

static int is_quda_initialized = 0;


int initialize_quda(void){

  QudaInitArgs_t init_args;

#if defined(SET_QUDA_SILENT)
  init_args.verbosity = QUDA_SILENT;
#elif defined(SET_QUDA_VERBOSE)
  init_args.verbosity = QUDA_VERBOSE;
#elif defined(SET_QUDA_DEBUG_VERBOSE)
  init_args.verbosity = QUDA_DEBUG_VERBOSE;
#else
  init_args.verbosity = QUDA_SUMMARIZE; /* default */
#endif

  const int dim[4] = {nx, ny, nz, nt};
  int status = 0;

  if(is_quda_initialized)return status;

  init_args.layout.device = 0; 								// only valid for single-gpu build
  init_args.layout.latsize = dim;
  init_args.layout.machsize = get_logical_dimensions();

  /* Tell QUDA which communicator we are using, in case we have split it */
  qudaSetMPICommHandle(mycomm());
  qudaInit(init_args);

  //cudaDeviceSetLimit(cudaLimitPrintfFifoSize,128*1024*1024); // DEBUG-JNS


  if(status == 0)
    is_quda_initialized = 1;

  return status;

}

void finalize_quda(void){
  /* Free any QUDA-resident deflation space.  A space can be created by the
     deflated CG, the exact-current path, or the QUDA eigensolver, so clean up
     for any of those consumers.  qudaCleanUpDeflationSpace() is a safe no-op
     when nothing was allocated. */
#if defined(USE_CG_GPU) || defined(USE_EIG_GPU) || defined(USE_CURRENT_GPU)
  qudaCleanUpDeflationSpace();
#endif
#if defined(USE_CG_GPU) && defined(MULTIGRID)
  mat_invert_mg_cleanup();
#endif
  qudaFinalize();
}

/** Load QUDA eig_args from the MILC eigen_param structure
*
*   The eig_args struct is needed whenever QUDA handles eigenvectors,
*   e.g., for an eigensolve, or for a deflated CG solve.
**/
void load_quda_default_eig_args(QudaEigensolverArgs_t *eig_args, int quda_does_eigensolve){

  char myname[] = "load_quda_default_eig_args";
  
  // Set default values
  
  // The param.eigen_param structure does not necessarily have the
  // QUDA parameters needed above.  So these are mostly dummy eig_args.  
    
  // QUDA requires the eig_args structure in order to set up the QUDA
  // deflation space and run the QUDA deflated solver.  So to be safe,
  // we set dummy values where they should be irrelevant.
  eig_args->struct_size = 1192;
  eig_args->block_size = 1;
  // When doing a fresh eigensolve, n_conv is the number of eigenvectors requested
  // when loading from file, n_conv is the total number of eigenvectors in the file
  eig_args->n_conv = param.eigen_param.Nvecs_in;
  // Number of eigenvectors to use for deflation must be set with each inversion
  eig_args->n_ev_deflate = param.eigen_param.Nvecs_in;
  eig_args->n_ev = param.eigen_param.Nvecs_in;
  eig_args->n_kr = param.eigen_param.Nvecs_in + 20;
  eig_args->tol = 1e-8;
  eig_args->max_restarts = 0;
  eig_args->poly_deg = 0;
  eig_args->a_min = 0.;
  eig_args->a_max = 0.;
  eig_args->preserve_evals = QUDA_BOOLEAN_TRUE; // Default to preserving the eigenvalues
  eig_args->batched_rotate = 0;
  /** Precision at which eigenvectors are written to disk.  This follows the
     * precision at which they were computed/held (eigensolver_prec), not the
     * compiled MILC precision: eigensolver_prec 2 -> double, otherwise single.
     * (eigensolver_prec 0 is half, which cannot be saved, so it maps to single;
     * single storage is lossless for single- or half-precision eigenvectors.)
     * This lets, e.g., a double-precision build with eigensolver_prec 1 store
     * a single-precision deflation space in single-precision files.
  **/
  eig_args->save_prec = (param.eigen_param.eigPrec == 2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  eig_args->partfile = QUDA_BOOLEAN_FALSE;
  eig_args->io_parity_inflate = QUDA_BOOLEAN_FALSE;
  eig_args->use_norm_op = QUDA_BOOLEAN_FALSE;
  eig_args->use_pc = QUDA_BOOLEAN_TRUE;
  eig_args->tol_restart = 1e-2;
  eig_args->eig_type = QUDA_EIG_TR_LANCZOS;
  eig_args->spectrum = QUDA_SPECTRUM_SR_EIG; /* Smallest Real. Other options: LM, SM, LR, SR, LI, SI */
  eig_args->qr_tol = 1e-8;
  eig_args->require_convergence = QUDA_BOOLEAN_FALSE;
  eig_args->check_interval = 10;
  eig_args->use_dagger = QUDA_BOOLEAN_FALSE;
  eig_args->compute_gamma5 = QUDA_BOOLEAN_FALSE;
  eig_args->compute_svd = QUDA_BOOLEAN_FALSE;
  eig_args->use_eigen_qr = QUDA_BOOLEAN_TRUE;
  eig_args->use_poly_acc = QUDA_BOOLEAN_TRUE;
  eig_args->arpack_check = QUDA_BOOLEAN_FALSE;
  eig_args->compute_evals_batch_size = 16;
  eig_args->preserve_deflation = QUDA_BOOLEAN_TRUE;
  strcpy( eig_args->vec_infile, "" );
  strcpy( eig_args->vec_outfile, param.ks_eigen_savefile );
  /* Precision at which QUDA holds/applies the deflation space, from the
     input-file parameter eigensolver_prec.  Applied whether QUDA computes
     the eigenvectors or loads them from file. */
  if(param.eigen_param.eigPrec == 2) {
    eig_args->prec_eigensolver = QUDA_DOUBLE_PRECISION;
  } else if(param.eigen_param.eigPrec == 1) {
    eig_args->prec_eigensolver = QUDA_SINGLE_PRECISION;
  } else if(param.eigen_param.eigPrec == 0) {
    eig_args->prec_eigensolver = QUDA_HALF_PRECISION;
  } else {
    printf("%s: Unrecognized eigensolver precision\n",myname);
    terminate(2);
  }

  if(quda_does_eigensolve){
    
    // In this case QUDA does the eigensolve and we need the proper eig_args
    
    int blockSize = param.eigen_param.blockSize;
    
    // Block size for block Lanczos
    eig_args->block_size = blockSize;
    // The number of converged eigenvectors requested
    eig_args->n_conv = (param.eigen_param.Nvecs_in > param.eigen_param.Nvecs) ? param.eigen_param.Nvecs_in : param.eigen_param.Nvecs;
    // The size of the eigenvector search space
    eig_args->n_ev = eig_args->n_conv;
    // Size of the Krylov subspace to use in the eigensolver, should be 1.5-2 times n_ev
    eig_args->n_kr = (param.eigen_param.Nkr < eig_args->n_ev ) ? 2*eig_args->n_ev : param.eigen_param.Nkr;
    // Ritz residual tolerance
    eig_args->tol = param.eigen_param.tol;
    // Lanczos restarts
    eig_args->max_restarts = param.eigen_param.MaxIter;
    // Polynomial acceleration
    eig_args->poly_deg = param.eigen_param.poly.norder;
    eig_args->a_min = param.eigen_param.poly.minE;
    eig_args->a_max = param.eigen_param.poly.maxE;
    // The max number of extra eigenvectors that solver may allocate to perform a Ritz rotation
    eig_args->batched_rotate = param.eigen_param.batchedRotate;
    eig_args->partfile = param.eigen_param.partfile ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
    eig_args->tol_restart = param.eigen_param.tol_restart;
    // Lanczos or block Lanczos
    eig_args->eig_type = ( eig_args->block_size > 1 ) ? QUDA_EIG_BLK_TR_LANCZOS : QUDA_EIG_TR_LANCZOS;
    eig_args->qr_tol = eig_args->tol;
    eig_args->require_convergence = QUDA_BOOLEAN_TRUE;
    strcpy( eig_args->vec_infile, param.ks_eigen_startfile );
  }
} // load_quda_default_eig_args

void print_quda_eig_args(QudaEigensolverArgs_t *eig_args){

  char myname[] = "print_quda_eig_args";

  node0_printf("%s:", myname);
  node0_printf(" struct_size = %lu\n", eig_args->struct_size);
  node0_printf(" n_conv = %d\n", eig_args->n_conv);
  node0_printf(" n_ev_deflate = %d\n", eig_args->n_ev_deflate);
  node0_printf(" n_ev = %d\n", eig_args->n_ev);
  node0_printf(" n_kr = %d\n", eig_args->n_kr);
  node0_printf(" eig_type = %d\n", eig_args->eig_type);
  node0_printf(" spectrum = %d\n", eig_args->spectrum);
  node0_printf(" block_size = %d\n", eig_args->block_size);
  node0_printf(" max_restarts = %d\n", eig_args->max_restarts);
  node0_printf(" batched_rotate = %d\n", eig_args->batched_rotate);
  node0_printf(" tol_restart = %e\n", eig_args->tol_restart);
  node0_printf(" use_poly_acc = %d\n", eig_args->use_poly_acc);
  node0_printf(" poly_deg = %d\n", eig_args->poly_deg);
  node0_printf(" a_min = %e\n", eig_args->a_min);
  node0_printf(" a_max = %e\n", eig_args->a_max);
  node0_printf(" compute_evals_batch_size = %d\n", eig_args->compute_evals_batch_size);
  node0_printf(" tol = %e\n", eig_args->tol);
  node0_printf(" qr_tol = %e\n", eig_args->qr_tol);
  node0_printf(" prec_eigensolver = %d\n", eig_args->prec_eigensolver);
  node0_printf(" require_convergence = %d\n", eig_args->require_convergence);
  node0_printf(" check_interval = %d\n", eig_args->check_interval);
  node0_printf(" preserve_evals = %d\n", eig_args->preserve_evals);
  node0_printf(" preserve_deflation = %d\n", eig_args->preserve_deflation);
  node0_printf(" use_norm_op = %d\n", eig_args->use_norm_op);
  node0_printf(" use_pc = %d\n", eig_args->use_pc);
  node0_printf(" use_dagger = %d\n", eig_args->use_dagger);
  node0_printf(" use_eigen_qr = %d\n", eig_args->use_eigen_qr);
  node0_printf(" compute_gamma5 = %d\n", eig_args->compute_gamma5);
  node0_printf(" compute_svd = %d\n", eig_args->compute_svd);
  node0_printf(" arpack_check = %d\n", eig_args->arpack_check);
  node0_printf(" vec_infile = %s\n", eig_args->vec_infile);
  node0_printf(" vec_outfile = %s\n", eig_args->vec_outfile);
  node0_printf(" save_prec = %d\n", eig_args->save_prec);
  node0_printf(" partfile = %d\n", eig_args->partfile);
  node0_printf(" io_parity_inflate = %d\n", eig_args->io_parity_inflate);

} // print_quda_eig_args
/* milc_to_quda_utilities */














