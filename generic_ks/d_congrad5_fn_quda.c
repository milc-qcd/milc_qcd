/******* d_congrad5_fn_gpu.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* GPU version of d_congrad5_fn_milc.c.  Can be compiled together with it. */
// Note that the restart criterion used by QUDA is different from 
// the restart criterion used by the MILC code

#if ! ( defined(USE_CG_GPU) && defined(HAVE_QUDA))
#error "d_congrad5_fn_quda.c requires USE_CG_GPU and HAVE_QUDA"
#endif

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
/* Solution of the normal equations for a single site parity with   */
/* multiple right sides                                             */
/********************************************************************/

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

  /* Load eig_args structure with default values */
  QudaEigensolverArgs_t eig_args;
  
#ifdef USE_EIG_GPU
  // Here we use QUDA for the eigensolution or for reading its own eigenvector file
  int quda_does_eigensolve = (param.ks_eigen_startflag == FRESH);
  load_quda_default_eig_args(&eig_args, quda_does_eigensolve);
  strcpy( eig_args.vec_infile, param.ks_eigen_startfile );
  eig_args.n_conv = (param.eigen_param.Nvecs_in > param.eigen_param.Nvecs) ? param.eigen_param.Nvecs_in : param.eigen_param.Nvecs;
  eig_args.n_ev = eig_args.n_conv;
#else
  // Here, eigenvectors were loaded from file(s) by MILC or by a non-QUDA eigensolver
  int quda_does_eigensolve = 0;
  load_quda_default_eig_args(&eig_args, quda_does_eigensolve);
#endif
      
  // Adjustments to default eig_args
#ifdef USE_EIG_GPU
  /* tol_restart is read from the input file only when QUDA owns the eigensolve
     (USE_EIG_GPU); otherwise keep the default from load_quda_default_eig_args(). */
  eig_args.tol_restart = param.eigen_param.tol_restart;
#endif
  eig_args.n_ev_deflate = ( qic->deflate ) ? param.eigen_param.Nvecs : 0;
  // QUDA currently doesn't support deflation with the relative residual stopping condition 
  if(qic->relresid > 0.) eig_args.n_ev_deflate = 0;

  int parity = qic->parity;
  eig_args.use_norm_op = ( parity == EVENANDODD ) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_args.use_pc = ( parity != EVENANDODD) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;

  if(eig_args.n_ev_deflate == 0){
    node0_printf("Solving for %d source(s) without deflation for parity %d\n", nsrc, parity);
  } else {
    node0_printf("Solving for %d source(s) with deflation for parity %d\n", nsrc, parity);
  }
 
#ifdef CG_DEBUG 
  print_quda_eig_args(&eig_args);
#endif

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
#ifdef CG_DEBUG
  for(int j = 0; j < nsrc; j++){
    node0_printf("Calling check_invert_field2 for case %d\n", j); fflush(stdout);
    check_invert_field2(t_src[j], t_dest[j], mass, 2e-5, fn, qic->parity);
  }
#endif
  // On the other hand, MILC expects the returned value to be the aggregate number of iterations
  // performed by each solve if it was performed sequentially. This can be approximated by
  // the number of iterations for a single solve times the number of sources.
  return num_iters * nsrc;
}

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

  int nsrc = 1;
  su3_vector *my_src[] = { t_src };
  su3_vector *my_dest[] = { t_dest };

  int num_iters = ks_congrad_block_parity_gpu(nsrc, my_src, my_dest, qic, mass, fn);

  return num_iters;
}


