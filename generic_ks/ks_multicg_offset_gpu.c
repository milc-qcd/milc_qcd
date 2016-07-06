/************** ks_multicg_offset_gpu.c **************************/
/* MIMD version 7 */

/* The following headers are supplied with the MILC code */
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"

#include <quda.h>
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#define LOOPEND
#include "../include/loopend.h"

/* Backward compatibility*/
#ifdef SINGLE_FOR_DOUBLE
#define HALF_MIXED
#endif

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

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
    printf("%s: GPU code does not yet support a Fermilab-type relative residual\n", myname);
    terminate(1);
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
    printf("CONGRAD5: time = %e (fn %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], qic->final_iters, 
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
  for(i=0; i<num_offsets; ++i) offset[i] = ksp[i].offset;

  double tmp;
  for(i=0; i<num_offsets; ++i){
    tmp = ksp[i].offset;
    offset[i] = tmp;
    node0_printf("offset[%d] = %g\n",i,offset[i]);
  }


  // The following arrays temporarily hold the final residuals
  double* residual = (double*)malloc(num_offsets*sizeof(double));
  double* relative_residual = (double*)malloc(num_offsets*sizeof(double));
  double* final_residual = (double*)malloc(num_offsets*sizeof(double));
  double* final_relative_residual = (double*)malloc(num_offsets*sizeof(double));

  for(i=0; i<num_offsets; ++i){

   residual[i]          = qic[i].resid;
   if (i>0){
#if defined(MAX_MIXED) || defined(HALF_MIXED)
       residual[i] = qic[i].resid; // for a mixed precision solver use residual for higher shifts
#else
       residual[i] = 0; // a unmixed solver should iterate until breakdown to agreee with CPU behavior
#endif
   }
   relative_residual[i] = qic[i].relresid;
   node0_printf("residual[%d] = %g relative %g\n",i, residual[i], relative_residual[i]);
  }

  inv_args.max_iter  = qic[0].max*qic[0].nrestart;
#if defined(MAX_MIXED) || defined(HALF_MIXED)
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif

  int num_iters = 0; // number of iterations taken
  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);

  initialize_quda();

  // for newer versions of QUDA we need to invalidate the gauge field if the naik term changes to prevent caching
  static int naik_term_epsilon_index = -1; 
  if ( naik_term_epsilon_index != ksp[0].naik_term_epsilon_index) {
    num_iters = -1; // temporary back door hack to invalidate gauge fields since naik index has changed
    naik_term_epsilon_index = ksp[0].naik_term_epsilon_index;
  }

  qudaMultishiftInvert(
		       PRECISION,
		       qic[0].prec,
		       num_offsets,
		       offset,
		       inv_args,
		       residual,
		       relative_residual,
		       fatlink,
		       longlink,
                       u0,
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

  free(offset);
  free(residual);
  free(relative_residual);
  free(final_residual);
  free(final_relative_residual);

#ifdef CGTIME
  dtimec += dclock();
  if(this_node==0){
    printf("CONGRAD5: time = %e (multicg_offset_QUDA %s) masses = %d iters = %d mflops = %e\n",
	   dtimec,prec_label[qic[0].prec-1],num_offsets,num_iters,
	   (double)(nflop)*volume*
	   num_iters/(1.0e6*dtimec*numnodes()));
    fflush(stdout);}
#endif

  return num_iters;
}
