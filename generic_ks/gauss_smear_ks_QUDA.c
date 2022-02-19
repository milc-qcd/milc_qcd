/************************** gauss_smear_ks_QUDA.c ********************************/
/* MIMD version 7 */
/*
 * Create a Gauss-smeared source using QUDA
 */

#include "generic_ks_includes.h"

#ifdef USE_GSMEAR_QUDA

#include <string.h>
#include <assert.h>

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

#define GS_TIME
/* #define GS_DEBUG */

/* Compute two-link if it is non-zero. */
static int compute_2link = 1;

void
gauss_smear_v_field_QUDA(su3_vector *src, su3_matrix *t_links,
                         Real width, int iters, int t0)
{
  char myname[] = "gauss_smear_v_field_QUDA";

#ifdef GS_DEBUG
  node0_printf( "%s: Start\n", myname );
#endif

#ifdef GS_TIME
  double dtimec;
#endif
  
  /* initialize QUDA */
  if( initialize_quda() )
  {
    node0_printf( "%s: FATAL. QUDA has not been initialized.\n", myname );
    terminate(1);
  }
  
  /* input parameters ***************************/
  int laplaceDim = 3;
  /**********************************************/

  if( laplaceDim > 3 && t0 != ALL_T_SLICES )
  {
    node0_printf( "%s: [Warning] t0 is ignored for d>3 dimensional Laplacian.\n", myname );
    t0 = ALL_T_SLICES;
  }
  
#ifdef GS_DEBUG
  setVerbosityQuda( QUDA_DEBUG_VERBOSE, "", stdout );
#endif
  
  /* QUDA inverter setup ************************/
  QudaInvertParam qip = newQudaInvertParam();

  qip.verbosity = QUDA_SUMMARIZE;

  qip.dslash_type = QUDA_ASQTAD_DSLASH;
  qip.laplace3D = laplaceDim;
  qip.Ls = 1;

  qip.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  qip.mass = 0.0;
  qip.kappa = 1.0;

  qip.dagger = QUDA_DAG_NO;
  qip.gauge_smear = QUDA_BOOLEAN_FALSE;

  qip.cpu_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  qip.cuda_prec = qip.cpu_prec;
  qip.cuda_prec_sloppy = qip.cpu_prec;
  qip.cuda_prec_refinement_sloppy = qip.cuda_prec;
  qip.dirac_order = QUDA_DIRAC_ORDER;
  qip.input_location = QUDA_CPU_FIELD_LOCATION;
  qip.output_location = QUDA_CPU_FIELD_LOCATION;
  /* Removed from QUDA */
  /* qip.sp_pad = 0; */ 
  /* qip.cl_pad = 0; */

  /* Not used, but need to be set. */
  qip.inv_type = QUDA_CG_INVERTER;
  qip.solution_type = QUDA_MAT_SOLUTION;
  qip.solve_type = QUDA_DIRECT_SOLVE;
  qip.matpc_type = QUDA_MATPC_EVEN_EVEN;
  qip.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  qip.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  qip.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
  qip.tol = 1e-14;
  qip.tol_restart = 1e-14;
  qip.maxiter = 1;
  qip.reliable_delta = 0.1;
  qip.use_alternative_reliable = QUDA_BOOLEAN_FALSE;
  qip.use_sloppy_partial_accumulator = QUDA_BOOLEAN_FALSE;
  qip.solution_accumulator_pipeline = 0;
  qip.pipeline = 0;
  qip.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  qip.heavy_quark_check = 0;
  qip.tol_hq = 0;
  qip.Nsteps = 0;
  qip.inv_type_precondition = QUDA_INVALID_INVERTER;
  qip.tol_precondition = 1e-14;
  qip.maxiter_precondition = 1;
  qip.verbosity_precondition = QUDA_SILENT;
  qip.cuda_prec_precondition = QUDA_INVALID_PRECISION;
  qip.cuda_prec_eigensolver = QUDA_INVALID_PRECISION;
  qip.gcrNkrylov = 1;
  qip.ca_basis = QUDA_POWER_BASIS;
  qip.ca_lambda_min = 0.0;
  qip.ca_lambda_max = 0.0;
  qip.native_blas_lapack = QUDA_BOOLEAN_TRUE;
  /* End of *************** QUDA inverter setup */

  /* QUDA gauge setup ***************************/
  QudaGaugeParam qgp = newQudaGaugeParam();

  qgp.type = QUDA_SU3_LINKS;

  const int * nsquares = get_logical_dimensions();
  qgp.X[0] = nx / nsquares[0];
  qgp.X[1] = ny / nsquares[1];
  qgp.X[2] = nz / nsquares[2];
  qgp.X[3] = nt / nsquares[3];  

  qgp.cpu_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  qgp.cuda_prec = qgp.cpu_prec;
  qgp.cuda_prec_sloppy = qgp.cuda_prec;
  qgp.cuda_prec_precondition= qgp.cuda_prec;
  qgp.cuda_prec_eigensolver = qgp.cuda_prec;
  qgp.cuda_prec_refinement_sloppy = qgp.cuda_prec;

  qgp.reconstruct = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_eigensolver = QUDA_RECONSTRUCT_NO;
  qgp.reconstruct_refinement_sloppy = QUDA_RECONSTRUCT_NO;

  qgp.gauge_order = QUDA_MILC_GAUGE_ORDER;
  qgp.anisotropy = 1.0;
  qgp.t_boundary = QUDA_PERIODIC_T;
  qgp.tadpole_coeff = u0;
  qgp.gauge_fix = QUDA_GAUGE_FIXED_NO;
  qgp.staggered_phase_type = QUDA_STAGGERED_PHASE_NO;

  int pad_size = 0;
  int x_face_size = qgp.X[1] * qgp.X[2] * qgp.X[3] / 2;
  int y_face_size = qgp.X[0] * qgp.X[2] * qgp.X[3] / 2;
  int z_face_size = qgp.X[0] * qgp.X[1] * qgp.X[3] / 2;
  int t_face_size = qgp.X[0] * qgp.X[1] * qgp.X[2] / 2;
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
  pad_size = MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);

  qgp.ga_pad = pad_size;
  qgp.mom_ga_pad = 0;
  /* End of ****************** QUDA gauge setup */
  
  /* Load gauge field */
  loadGaugeQuda( (void*) t_links, &qgp );

#ifdef GS_TIME
  dtimec = -dclock();
#endif
#ifdef GS_DEBUG
  node0_printf( "%s: Gaussian smearing starts.\n", myname );
#endif
  
  /* Run gaussian smearing */
  performTwoLinkGaussianSmearNStep( (void*) src, &qip, iters, width, compute_2link, t0 );

#ifdef GS_DEBUG
  node0_printf( "%s: Gaussian smearing ends.\n", myname );
#endif  
#ifdef GS_TIME
  dtimec += dclock();
  node0_printf( "[GS_TIME] QUDA two-link Gaussian smearing: time = %g s, iters = %d\n", dtimec, iters );
#endif
  
  /* Disable two-link calculation in subsequent calls */
  if( compute_2link != 0 ) compute_2link = 0;  
  
#ifdef GS_DEBUG
  node0_printf( "%s: End\n", myname );
#endif  

  return ;
}


#else /* #ifdef USE_GSMEAR_QUDA */

void
gauss_smear_v_field_QUDA(su3_vector *src, su3_matrix *t_links,
        Real width, int iters, int t0)
{
  char myname[] = "gauss_smear_v_field_QUDA";
  node0_printf( "%s: Requires compilation with the QUDA library\n", myname );
  terminate(1);

  return ;  
}

#endif
