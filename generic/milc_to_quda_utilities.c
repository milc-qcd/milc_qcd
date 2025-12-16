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
#ifdef USE_CG_GPU
  qudaCleanUpDeflationSpace();
#ifdef MULTIGRID
  mat_invert_mg_cleanup();
#endif
#endif
  qudaFinalize();
}

/* Load QUDA eigen_args from the MILC eigen_param structure */
void load_quda_default_eig_args(QudaEigensolverArgs_t *eig_args){

  char myname[] = "load_quda_default_eig_args";

#ifdef USE_EIG_GPU

  // QUDA eig_args
  
  int blockSize = param.eigen_param.blockSize;

  eig_args->struct_size = 1192;
  eig_args->block_size = blockSize;
  eig_args->n_conv = param.eigen_param.Nvecs;
  eig_args->n_ev_deflate = param.eigen_param.Nvecs;
  eig_args->n_ev = param.eigen_param.Nvecs;
  eig_args->n_kr = param.eigen_param.Nkr;
  eig_args->tol = param.eigen_param.tol;
  eig_args->max_restarts = param.eigen_param.MaxIter;
  eig_args->poly_deg = param.eigen_param.poly.norder;
  eig_args->a_min = param.eigen_param.poly.minE;
  eig_args->a_max = param.eigen_param.poly.maxE;
  eig_args->preserve_evals = QUDA_BOOLEAN_TRUE; // Default to preserving the eigenvalues
  eig_args->batched_rotate = param.eigen_param.batchedRotate;
  /** With save_prec, we are currently saving the eigenvectors in double
   * precision if running MILC in double precision. At some point, we may
   * want to be able to control this separately via the parameters input file
   **/
  eig_args->save_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  eig_args->partfile = param.eigen_param.partfile ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args->io_parity_inflate = QUDA_BOOLEAN_FALSE;
  eig_args->use_norm_op = QUDA_BOOLEAN_FALSE;
  eig_args->use_pc = QUDA_BOOLEAN_TRUE;
  eig_args->tol_restart = param.eigen_param.tol_restart;
  eig_args->eig_type = ( eig_args->block_size > 1 ) ? QUDA_EIG_BLK_TR_LANCZOS : QUDA_EIG_TR_LANCZOS;  /* or QUDA_EIG_IR_ARNOLDI, QUDA_EIG_BLK_IR_ARNOLDI */
  eig_args->spectrum = QUDA_SPECTRUM_SR_EIG; /* Smallest Real. Other options: LM, SM, LR, SR, LI, SI */
  eig_args->qr_tol = eig_args->tol;
  eig_args->require_convergence = QUDA_BOOLEAN_TRUE;
  eig_args->check_interval = 10;
  eig_args->use_dagger = QUDA_BOOLEAN_FALSE;
  eig_args->compute_gamma5 = QUDA_BOOLEAN_FALSE;
  eig_args->compute_svd = QUDA_BOOLEAN_FALSE;
  eig_args->use_eigen_qr = QUDA_BOOLEAN_TRUE;
  eig_args->use_poly_acc = QUDA_BOOLEAN_TRUE;
  eig_args->arpack_check = QUDA_BOOLEAN_FALSE;
  eig_args->compute_evals_batch_size = 16;
  eig_args->preserve_deflation = QUDA_BOOLEAN_TRUE;
  strcpy( eig_args->vec_infile, param.ks_eigen_startfile );
  strcpy( eig_args->vec_outfile, param.ks_eigen_savefile );
  
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

#else

  // In this case we are not using the QUDA eigensolver and we are not
  // asking QUDA to read its own eigenvector file

  // The param.eigen_param structure does not necessarily have the
  // QUDA parameters needed above.  So these are mostly dummy eig_args.  

  // QUDA requires the eig_args structure in order to set up the QUDA
  // deflation space and run the QUDA deflated solver.  So to be safe,
  // we set dummy values where they should be irrelevant.

  eig_args->struct_size = 1192;
  eig_args->block_size = 1;
  eig_args->n_conv = param.eigen_param.Nvecs;
  eig_args->n_ev_deflate = param.eigen_param.Nvecs;
  eig_args->n_ev = param.eigen_param.Nvecs;
  eig_args->n_kr = param.eigen_param.Nvecs + 10;
  eig_args->tol = 1e-8;
  eig_args->max_restarts = 0;
  eig_args->poly_deg = 0;
  eig_args->a_min = 0.;
  eig_args->a_max = 0.;
  eig_args->preserve_evals = QUDA_BOOLEAN_TRUE;
  eig_args->batched_rotate = 0;
  eig_args->save_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  eig_args->partfile = QUDA_BOOLEAN_FALSE;
  eig_args->io_parity_inflate = QUDA_BOOLEAN_FALSE;
  eig_args->use_norm_op = QUDA_BOOLEAN_FALSE;
  eig_args->use_pc = QUDA_BOOLEAN_TRUE;
  eig_args->tol_restart = 1e-2;
  eig_args->eig_type = QUDA_EIG_TR_LANCZOS;
  eig_args->spectrum = QUDA_SPECTRUM_SR_EIG;
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
  eig_args->prec_eigensolver = QUDA_DOUBLE_PRECISION;
#endif

}
/* milc_to_quda_utilities */
