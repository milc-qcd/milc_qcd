/************ eigen_stuff_QUDA.c **************/
/* Eigenvalue and Eigevector computation routines.
 * This version uses QUDA (https://github.com/lattice/quda)
 * MIMD version 7 */
/**********************************************/

#include "generic_ks_includes.h"

#ifdef USE_EIG_GPU

#include <string.h>
#include <assert.h>

/* #include <lattice.h> */

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

/* #define EIG_DEBUG */
/* #define EIGTIME */

/* Compute eigenvalues and eigenvectors of the Kogut-Susskind
 * dslash^2. */
int ks_eigensolve_QUDA( su3_vector ** eigVec,
                        Real * eigVal,
                        ks_eigen_param * eigen_param,
                        int init )
{
  char myname[] = "ks_eigensolve_QUDA";
  int ii; /* loop index */

#ifdef EIG_DEBUG
  node0_printf( "%s: Start\n", myname );
#endif

  /* initialize QUDA */
  if( initialize_quda() )
  {
    node0_printf( "%s: FATAL. QUDA has not been initialized.\n", myname );
    terminate(1);
  }
  
#ifdef EIGTIME
  double dtimec;
#endif

  /* input parameters for QUDA eigensolver **********************/
  int Nvecs = eigen_param->Nvecs;
  int Nvecs_in = eigen_param->Nvecs_in;
  double tol = eigen_param->tol;
  int maxIter = eigen_param->MaxIter;
  int Nkr = eigen_param->Nkr;
  double aMin = eigen_param->poly.minE;
  double aMax = eigen_param->poly.maxE;
  int polyDeg = eigen_param->poly.norder;
  int blockSize = eigen_param->blockSize;
  int parity = eigen_param->parity;
  int batchedRotate = eigen_param->batchedRotate;
  int precEigensolver = eigen_param->eigPrec;
  /**************************************************************/
  
  /* QUDA inverter setup *************************/  
  QudaInvertParam qip = newQudaInvertParam();

  qip.verbosity = QUDA_VERBOSE; /* SILENT, SUMMARIZE, VERBOSE, DEBUG_VERBOSE */

  qip.dslash_type = QUDA_ASQTAD_DSLASH;

  if( parity == EVEN )
  {
    qip.solution_type = QUDA_MATPC_SOLUTION;
    qip.solve_type = QUDA_NORMOP_PC_SOLVE;
    qip.matpc_type = QUDA_MATPC_EVEN_EVEN;
  }
  else if( parity == ODD )
  {
    qip.solution_type = QUDA_MATPC_SOLUTION;
    qip.solve_type = QUDA_NORMOP_PC_SOLVE;
    qip.matpc_type = QUDA_MATPC_ODD_ODD;
  }
  else /*if ( parity == EVENANDODD ) */
  {
    qip.solution_type = QUDA_MAT_SOLUTION;
    qip.solve_type = QUDA_DIRECT_SOLVE;
    qip.matpc_type = QUDA_MATPC_EVEN_EVEN; /* meaningless */
  }
  qip.dagger = QUDA_DAG_NO;
  qip.mass_normalization = QUDA_MASS_NORMALIZATION;

  qip.mass = 0; /* set mass=0 */
  qip.kappa = 1.0 / ( 8.0 + qip.mass ); /* for Laplace operator */
  qip.laplace3D = 4; /* for Laplace operator */
  
  qip.cpu_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  qip.cuda_prec = qip.cpu_prec;
  qip.cuda_prec_sloppy = qip.cuda_prec;
  qip.cuda_prec_refinement_sloppy = qip.cuda_prec;
  qip.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  qip.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  qip.dirac_order = QUDA_DIRAC_ORDER;

  qip.input_location = QUDA_CPU_FIELD_LOCATION;
  qip.output_location = QUDA_CPU_FIELD_LOCATION;

  qip.native_blas_lapack = QUDA_BOOLEAN_TRUE;
  
  /* Below are meaningless in eigensovler, but need to be set. */
  qip.inv_type = QUDA_CG_INVERTER;
  qip.tol = 1e-7;
  qip.tol_restart = 5e+3 * qip.tol;
  qip.maxiter = 1000;
  qip.reliable_delta = 0.1;
  qip.use_alternative_reliable = QUDA_BOOLEAN_FALSE;
  qip.use_sloppy_partial_accumulator = QUDA_BOOLEAN_FALSE;
  qip.solution_accumulator_pipeline = 0;
  qip.pipeline = 0;
  qip.Ls = 1;
  qip.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  qip.heavy_quark_check = 0;
  qip.tol_hq = 0;
  qip.Nsteps = 2;
  qip.inv_type_precondition = QUDA_INVALID_INVERTER;
  qip.tol_precondition = 1e-1;
  qip.maxiter_precondition = 10;
  qip.verbosity_precondition = QUDA_SILENT;
  qip.cuda_prec_precondition = QUDA_INVALID_PRECISION;
  qip.cuda_prec_eigensolver = QUDA_INVALID_PRECISION;
  qip.gcrNkrylov = 10;
  qip.ca_basis = QUDA_POWER_BASIS;
  qip.ca_lambda_min = 0.0;
  qip.ca_lambda_max = -1.0;
  /***********************************************/  

  /* QUDA gauge setup **********************************/
  QudaGaugeParam qgp = newQudaGaugeParam();

  qgp.type = ( qip.dslash_type == QUDA_STAGGERED_DSLASH || qip.dslash_type == QUDA_LAPLACE_DSLASH ) ? QUDA_SU3_LINKS : QUDA_ASQTAD_FAT_LINKS;

  const int * nsquares = get_logical_dimensions();
  
  qgp.X[0] = nx / nsquares[0];
  qgp.X[1] = ny / nsquares[1];
  qgp.X[2] = nz / nsquares[2];
  qgp.X[3] = nt / nsquares[3];

  qgp.cpu_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  qgp.cuda_prec = qgp.cpu_prec;
  qgp.cuda_prec_sloppy = qgp.cuda_prec;
  qgp.cuda_prec_precondition = qgp.cuda_prec;
  //qgp.cuda_prec_eigensolver = qgp.cuda_prec;

  if(precEigensolver == 2) {
    qgp.cuda_prec_eigensolver = QUDA_DOUBLE_PRECISION;
  } else if(precEigensolver == 1) {
    qgp.cuda_prec_eigensolver = QUDA_SINGLE_PRECISION;
  } else if(precEigensolver == 0) {
    qgp.cuda_prec_eigensolver = QUDA_HALF_PRECISION;
  } else {
    printf("%s: Unrecognized eigensolver precision\n",myname);
    terminate(2);
  }

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
  qgp.staggered_phase_type = QUDA_STAGGERED_PHASE_MILC;

  int pad_size = 0;
  int x_face_size = qgp.X[1] * qgp.X[2] * qgp.X[3] / 2;
  int y_face_size = qgp.X[0] * qgp.X[2] * qgp.X[3] / 2;
  int z_face_size = qgp.X[0] * qgp.X[1] * qgp.X[3] / 2;
  int t_face_size = qgp.X[0] * qgp.X[1] * qgp.X[2] / 2;
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )
  pad_size = MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);
  int fat_pad = pad_size;
  int link_pad = 3*pad_size;
  
  qgp.ga_pad = fat_pad; 
  qgp.mom_ga_pad = 0;
  /*****************************************************/

  /* load gauge field ********************************/
  imp_ferm_links_t * fn = get_fm_links( fn_links, 0);
  su3_matrix * fatLinks = get_fatlinks(fn);
  su3_matrix * longLinks = get_lnglinks(fn);

  /* qgp.type = QUDA_ASQTAD_FAT_LINKS; */
  loadGaugeQuda( (void*) fatLinks, &qgp );

  if( qip.dslash_type == QUDA_ASQTAD_DSLASH )
  {
    qgp.type = QUDA_ASQTAD_LONG_LINKS;      
    qgp.ga_pad = link_pad;
    qgp.staggered_phase_type = QUDA_STAGGERED_PHASE_NO;
    qgp.scale = - ( 1.0 + fn->eps_naik ) / ( 24.0 * qgp.tadpole_coeff * qgp.tadpole_coeff );
  
    loadGaugeQuda( (void*) longLinks, &qgp );
  }
  destroy_fn_links(fn);
  /***************************************************/
  
  /* quda eigensolver setup *************************/  
  QudaEigParam qep = newQudaEigParam();
  qep.invert_param = &qip;

  qep.eig_type = ( blockSize > 1 ) ? QUDA_EIG_BLK_TR_LANCZOS : QUDA_EIG_TR_LANCZOS;  /* or QUDA_EIG_IR_ARNOLDI, QUDA_EIG_BLK_IR_ARNOLDI */
    
  qep.spectrum = QUDA_SPECTRUM_SR_EIG; /* Smallest Real. Other options: LM, SM, LR, SR, LI, SI */
  qep.n_conv = Nvecs;
  qep.n_ev_deflate = qep.n_conv;
  qep.n_ev = qep.n_conv;
  qep.n_kr = Nkr;
  qep.block_size = blockSize;

  qep.tol = tol;
  qep.qr_tol = qep.tol;
  qep.max_restarts = maxIter;
  qep.batched_rotate = batchedRotate;
  qep.require_convergence = QUDA_BOOLEAN_TRUE;
  qep.check_interval = 1;

  qep.use_norm_op = ( parity == EVENANDODD ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  qep.use_pc = ( parity != EVENANDODD) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;  

  qep.use_dagger = QUDA_BOOLEAN_FALSE;
  qep.compute_gamma5 = QUDA_BOOLEAN_FALSE;
  qep.compute_svd = QUDA_BOOLEAN_FALSE;
  qep.use_eigen_qr = QUDA_BOOLEAN_TRUE;

  qep.use_poly_acc = QUDA_BOOLEAN_TRUE;
  qep.poly_deg = polyDeg;
  qep.a_min = aMin;
  qep.a_max = aMax;

  qep.arpack_check = QUDA_BOOLEAN_FALSE;
  strcpy( qep.arpack_logfile, "" );

  strcpy( qep.vec_infile, "" );
  strcpy( qep.vec_outfile, "" );
  qep.save_prec = (MILC_PRECISION==2) ? QUDA_DOUBLE_PRECISION : QUDA_SINGLE_PRECISION;
  qep.io_parity_inflate = QUDA_BOOLEAN_FALSE;
  /**************************************************/  

  /* print input parameters for eigensolver ******************/
  node0_printf( "========= parameters for eigensolver =========\n" );
  node0_printf( "Number of wanted eigenvalues: %d\n", qep.n_conv );
  node0_printf( "Krylov subspace size: %d\n", qep.n_kr );
  node0_printf( "Eigenvalue equation tolerance: %e\n", qep.tol );
  node0_printf( "Maximum iterations of Lanczos restarts: %d\n", qep.max_restarts );
  node0_printf( "Chebyshev polynomial - alpha (lower bound for exclusion): %g\n", qep.a_min );
  node0_printf( "Chebyshev polynomial - beta (upper bound for exclusion): %g\n", qep.a_max );
  node0_printf( "Chebyshev polynomial order: %d\n", qep.poly_deg );
  node0_printf( "Block size: %d\n", qep.block_size );
  node0_printf( "Even-odd parity: %d\n", parity );
  node0_printf( "==============================================\n" );
  /* End of ********* print input parameters for eigensolver */
  
  void ** eigVec_QUDA = (void**) malloc( Nvecs * sizeof(void*) );
  for( ii=0; ii<Nvecs; ii++ )
  {
    eigVec_QUDA[ii] = (void*) eigVec[ii];
  }
  double _Complex * eigVal_QUDA = (double _Complex *) malloc( Nvecs * sizeof( double _Complex) );
  
#ifdef EIGTIME
  dtimec = -dclock();
#endif
#ifdef EIG_DEBUG
  node0_printf( "%s: Eigensolver starts.\n", myname );  
#endif

  /* QUDA's eigensolver using Thick Restarted (Block) Lanczos algorithm */
  eigensolveQuda( eigVec_QUDA, eigVal_QUDA, &qep );

#ifdef EIG_DEBUG
  node0_printf( "%s: Eigensolver ended.\n", myname );  
#endif
#ifdef EIGTIME
  dtimec += dclock();
  if( blockSize > 1 )
  {
    node0_printf( "[EIGTIME] QUDA BLKTRL: time = %g s\n", dtimec );
  }
  else
  {
    node0_printf( "[EIGTIME] QUDA TRL: time = %g s\n", dtimec );
  }    
#endif

  /* copy QUDA eigenvalues to MILC */
  for( ii=0; ii<Nvecs; ii++ )
  {
    eigVal[ii] = *( ((double*)eigVal_QUDA) + 2*ii );
  }
  
  /* print eigenvalues */
  node0_printf( "BEGIN RESULTS\n" );
  for( ii=0; ii<Nvecs; ii++ )
  {
#ifdef EIG_DEBUG
    node0_printf( "Eigenvalue(%i) = %.15e \n", ii, eigVal[ii] );
#else
    node0_printf( "Eigenvalue(%i) = %g \n", ii, eigVal[ii] );
#endif
  }

  /* free **************/
  free( eigVal_QUDA );
  /* End of ******* free */
  
#ifdef EIG_DEBUG
  node0_printf( "%s: End\n", myname );
#endif

  return 0;
}

#else /* #ifdef QUDA_EIG */

int ks_eigensolve_QUDA( su3_vector ** eigeVec, Real * eigVal, ks_eigen_param * eigen_param, int init )
{
  char myname[] = "ks_eigensolve_QUDA";
  node0_printf( "%s: Requires compilation with the QUDA library\n", myname );
  terminate(1);

  return 0;
}

#endif
