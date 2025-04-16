/************** ks_multicg_offset_gpu.c **************************/
/* MIMD version 7 */

/* The following headers are supplied with the MILC code */
#include <string.h>
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"

#include <quda.h>
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#define LOOPEND
#include "../include/loopend.h"
#include <string.h>

/* Backward compatibility*/
#ifdef SINGLE_FOR_DOUBLE
#define HALF_MIXED
#endif

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

// these are used to store the most recent fermion link field passed to QUDA
static int naik_term_epsilon_index = -1; 
static imp_ferm_links_t *fn_last = NULL;

// return the most recent fermion link field passed to QUDA
imp_ferm_links_t* get_fn_last() {
  return fn_last;
}

// update the fermion link field passed to QUDA
void set_fn_last(imp_ferm_links_t *fn_last_new) {
  fn_last = fn_last_new;
}

int ks_multicg_offset_field_gpu(
    su3_vector *src,
    su3_vector **psim,
    ks_param *ksp,
    int num_offsets,
    quark_invert_control *qic,
    imp_ferm_links_t *fn
    )
{
  int i,j;
  char myname[] = "ks_multicg_offset_field_gpu";

#ifdef CGTIME
  double dtimec = -dclock();
  double nflop = 1205 + 15*num_offsets;
#endif

  if(qic[0].relresid != 0.){
    printf("%s: WARNING; GPU code does not yet support a Fermilab-type relative residual\n", myname);
    //    terminate(1);
  }

  /* Initialize structure */
  for(j = 0; j < num_offsets; j++){
    qic[j].final_rsq     = 0.;
    qic[j].final_relrsq  = 0.; /* No relative residual in use here */
    qic[j].size_r        = 0.;
    qic[j].size_relr     = 0.;
    qic[j].final_iters   = 0;
    qic[j].final_restart = 0;  /* No restarts with this algorithm */
    qic[j].converged     = 1;
  }

  if( num_offsets==0 )return(0);

  /* Compute source norm */
  double source_norm = 0.0;
  FORSOMEFIELDPARITY(i,qic[0].parity){
    source_norm += (double)magsq_su3vec( &src[i] );
  } END_LOOP;
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("ks_multicg_offset_gpu: source_norm = %e\n", (double)source_norm);
#endif

  /* Provide for trivial solution */
  if(source_norm == 0.0){
    /* Zero the solutions and return zero iterations */
    for(j = 0; j < num_offsets; j++){
      FORSOMEFIELDPARITY(i,qic->parity){
	memset(psim[j] + i, 0, sizeof(su3_vector));
      } END_LOOP;
    }

#ifdef CGTIME
  dtimec += dclock();
  if(this_node==0){
    printf("CONGRAD5: time = %e (multicg_offset_QUDA %s) masses = %d iters = %d mflops = %e\n",
	   dtimec, prec_label[qic[0].prec-1], num_offsets, qic->final_iters, 
	   (double)(nflop*volume*qic->final_iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

    return 0;
  }

  QudaInvertArgs_t inv_args;

  if(qic[0].parity == EVEN){
    inv_args.evenodd = QUDA_EVEN_PARITY;
  }else if(qic[0].parity == ODD){
    inv_args.evenodd = QUDA_ODD_PARITY;
  }else{
    printf("%s: Unrecognised parity\n", myname);
    terminate(2);
  }

  if(qic[0].parity==EVENANDODD){
    node0_printf("%s: EVENANDODD not supported\n", myname);
    terminate(1);
  }

  double* offset = (double*)malloc(num_offsets*sizeof(double));
  double* residue = (double*)malloc(num_offsets*sizeof(double));

  for(i=0; i<num_offsets; ++i){
    offset[i] = ksp[i].offset;
    residue[i] = ksp[i].residue;
#if defined(SET_QUDA_VERBOSE) || defined(SET_QUDA_DEBUG_VERBOSE)
    node0_printf("offset[%d] = %g\n",i,offset[i]);
#endif
  }


  // The following arrays temporarily hold the final residuals
  double* residual = (double*)malloc(num_offsets*sizeof(double));
  double* relative_residual = (double*)malloc(num_offsets*sizeof(double));
  double* final_residual = (double*)malloc(num_offsets*sizeof(double));
  double* final_relative_residual = (double*)malloc(num_offsets*sizeof(double));

  for(i=0; i<num_offsets; ++i){

   residual[i]          = qic[i].resid;
   if (i>0) {
#if defined(MAX_MIXED) || defined(HALF_MIXED)
     if (residue[i] != 0) {
       // scale the shifted residual relative to the residue
       residual[i] = fabs(residue[0] / residue[i]) * residual[0];
       if (residual[i] < 1e-14) residual[i] = 1e-14;
     } else {
       residual[i] = qic[i].resid; // for a mixed-precision solver use residual for higher shifts
     }
#else
     residual[i] = 0; // a unmixed solver should iterate until breakdown to agree with CPU behavior
#endif
   }
   relative_residual[i] = qic[i].relresid;

#if defined(SET_QUDA_VERBOSE) || defined(SET_QUDA_DEBUG_VERBOSE)
   node0_printf("residual[%d] = %g relative %g\n",i, residual[i], relative_residual[i]);
#endif
  }

  inv_args.max_iter = qic[0].max*qic[0].nrestart;
#if defined(MAX_MIXED) || defined(HALF_MIXED) // never do half precision with multi-shift solver
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif

  int num_iters = 0; // number of iterations taken
  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);

  initialize_quda();

  // for newer versions of QUDA we need to invalidate the gauge field if the naik term changes to prevent caching
  if ( fn != get_fn_last() || fresh_fn_links(fn) ){
    cancel_quda_notification(fn);
    set_fn_last(fn);
    num_iters = -1;
    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", myname);
  }

  if ( naik_term_epsilon_index != ksp[0].naik_term_epsilon_index) {
    num_iters = -1; // temporary back door hack to invalidate gauge fields since naik index has changed
    naik_term_epsilon_index = ksp[0].naik_term_epsilon_index;
    node0_printf("%s: naik_epsilon: Signal QUDA to refresh links\n", myname);
  }

  inv_args.naik_epsilon = ksp[0].naik_term_epsilon;
  if (inv_args.naik_epsilon != fn->eps_naik) {
    node0_printf("%s: naik_epsilon in action (%e) does not match value in link (%e)\n",
                 myname, inv_args.naik_epsilon, fn->eps_naik);
    terminate(1);
  }

#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  qudaMultishiftInvert(
		       MILC_PRECISION,
		       qic[0].prec,
		       num_offsets,
		       offset,
		       inv_args,
		       residual,
		       relative_residual,
		       fatlink,
		       longlink,
		       (void *)src,
		       (void **)psim,
		       final_residual,
		       final_relative_residual,
		       &num_iters);

  for(i=0; i<num_offsets; ++i){
    qic[i].final_rsq = final_residual[i]*final_residual[i];
    qic[i].final_relrsq = final_relative_residual[i]*final_relative_residual[i];
    qic[i].final_iters = num_iters;

    // check for convergence
    if(relative_residual[i]){
      qic[i].converged = (final_relative_residual[i] <= qic[i].relresid) ? 1 :0;
    }else{
      relative_residual[i] = 0.;// /* No relative residual in use here */ qic[i].relresid;
    }
    // Cumulative residual. Not used in practice
    qic[i].size_r = 0.0;
    qic[i].size_relr = 0.0;
  }

  free(residue);
  free(offset);
  free(residual);
  free(relative_residual);
  free(final_residual);
  free(final_relative_residual);

#ifdef CGTIME
  dtimec += dclock();
  if(this_node==0){
    printf("CONGRAD5: time = %e (multicg_offset_QUDA %s) masses = %d iters = %d mflops = %e\n",
	   dtimec, prec_label[qic[0].prec-1], num_offsets, num_iters,
	   (double)(nflop)*volume*num_iters/(1.0e6*dtimec*numnodes()));
    fflush(stdout);}
#endif

  return num_iters;
}
