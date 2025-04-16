/******** read_eigen_param.c *********/
/* MIMD version 7 */

/* Standardize eigenparameter input */

#include "generic_ks_includes.h"
#define IF_OK if(status==0)

int read_ks_eigen_param(ks_eigen_param *eigen_param, int status, int prompt){

#if defined(HAVE_PRIMME)

  IF_OK status += get_i(stdin, prompt,"Max_Rayleigh_iters", &eigen_param->MaxIter );
  IF_OK status += get_i(stdin, prompt,"Restart_Rayleigh", &eigen_param->Restart );
  IF_OK status += get_f(stdin, prompt,"eigenval_tolerance", &eigen_param->tol );
  
#elif defined(HAVE_ARPACK)
  
  IF_OK status += get_i(stdin, prompt,"Max_Rayleigh_iters", &eigen_param->MaxIter );
  IF_OK status += get_i(stdin, prompt,"nArnoldi", &eigen_param->nArnoldi );
  IF_OK status += get_f(stdin, prompt,"eigenval_tolerance", &eigen_param->tol );
  
#elif defined(HAVE_GRID) && defined(USE_EIG_GPU)
  
  IF_OK status += get_i(stdin, prompt, "Max_Lanczos_restart_iters", &eigen_param->MaxIter );
  IF_OK status += get_f(stdin, prompt, "eigenval_tolerance", &eigen_param->tol );
  IF_OK status += get_i(stdin, prompt, "Lanczos_max", &eigen_param->Nmax );
  IF_OK status += get_i(stdin, prompt, "Lanczos_restart", &eigen_param->Nrestart );
  IF_OK status += get_i(stdin, prompt, "Lanczos_reorth_period", &eigen_param->reorth_period );
  IF_OK status += get_f(stdin, prompt, "Chebyshev_alpha", &eigen_param->poly.minE );
  IF_OK status += get_f(stdin, prompt, "Chebyshev_beta", &eigen_param->poly.maxE );
  IF_OK status += get_i(stdin, prompt, "Chebyshev_order", &eigen_param->poly.norder );
  IF_OK status += get_s(stdin, prompt, "diag_algorithm", param.eigen_param.diagAlg );

#elif defined(HAVE_QUDA) && defined(USE_CG_GPU) && !defined(USE_EIG_GPU)
  node0_printf("ERROR: When using QUDA for CG and wanting FRESH eigenvectors, only the QUDA eigensolver is allowed!\n");
  node0_printf("ERROR: Recompile with WANT_QUDA=true, WANT_CG_GPU=true, and WANT_EIG_GPU=true\n");
  terminate(1);

#elif defined(HAVE_QUDA) && defined(USE_EIG_GPU)
  
  IF_OK status += get_i(stdin, prompt, "Max_Lanczos_restart_iters", &eigen_param->MaxIter );    
  IF_OK status += get_f(stdin, prompt, "eigenval_tolerance", &eigen_param->tol );
  IF_OK status += get_i(stdin, prompt, "Lanczos_max", &eigen_param->Nkr );
  IF_OK status += get_i(stdin, prompt, "Lanczos_restart", &eigen_param->Nrestart );
  IF_OK status += get_i(stdin, prompt, "eigensolver_prec", &eigen_param->eigPrec );
  IF_OK status += get_i(stdin, prompt, "batched_rotate", &eigen_param->batchedRotate );
  IF_OK status += get_f(stdin, prompt, "Chebyshev_alpha", &eigen_param->poly.minE );
  IF_OK status += get_f(stdin, prompt, "Chebyshev_beta", &eigen_param->poly.maxE );
  IF_OK status += get_i(stdin, prompt, "Chebyshev_order", &eigen_param->poly.norder );
  IF_OK status += get_i(stdin, prompt, "block_size", &param.eigen_param.blockSize );

#elif defined(HAVE_QDP)

  /* Incomplete Needs updating */
  IF_OK status += get_i(stdin, prompt,"Max_Rayleigh_iters", &eigen_param->MaxIter );
  IF_OK status += get_i(stdin, prompt,"Restart_Rayleigh", &eigen_param->Restart );
  IF_OK status += get_i(stdin, prompt,"Kalkreuter_iters", &eigen_param->Kiters );
  IF_OK status += get_f(stdin, prompt,"eigenval_tolerance", &eigen_param->tol );
 
#else
  
  /* Kalkreuter_Ritz */
  IF_OK status += get_i(stdin, prompt,"Max_Rayleigh_iters", &eigen_param->MaxIter );
  IF_OK status += get_i(stdin, prompt,"Restart_Rayleigh", &eigen_param->Restart );
  IF_OK status += get_i(stdin, prompt,"Kalkreuter_iters", &eigen_param->Kiters );
  IF_OK status += get_f(stdin, prompt,"eigenval_tolerance", &eigen_param->tol );
  IF_OK status += get_f(stdin, prompt,"error_decrease", &eigen_param->error_decr);
  
#endif
  
#ifdef POLY_EIGEN

  /* Chebyshev preconditioner */

#ifdef HAVE_ARPACK
  IF_OK status += get_i(stdin, prompt,"which_poly", &eigen_param->poly.which_poly );
#endif
  IF_OK status += get_i(stdin, prompt,"norder", &eigen_param->poly.norder);
  IF_OK status += get_f(stdin, prompt,"eig_start", &eigen_param->poly.minE);
  IF_OK status += get_f(stdin, prompt,"eig_end", &eigen_param->poly.maxE);
  
#ifdef HAVE_ARPACK
  IF_OK status += get_f(stdin, prompt,"poly_param_1", &eigen_param->poly.poly_param_1  );
  IF_OK status += get_f(stdin, prompt,"poly_param_2", &eigen_param->poly.poly_param_2  );
  IF_OK status += get_i(stdin, prompt,"eigmax", &eigen_param->poly.eigmax );
#endif

  /* TODO: Support QUDA and GRID Cheby cases */

#endif
  
  return status;
}
