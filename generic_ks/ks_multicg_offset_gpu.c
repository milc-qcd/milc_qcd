/************** ks_multicg_offset_gpu.c **************************/
/* MIMD version 7 */

/* The following headers are supplied with the MILC code */
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/loopend.h"

/* The following headers are supplied with QUDA */
#include <quda.h>
#include <milc_interface.h>


int ks_multicg_offset_field_gpu(
    su3_vector *src,
    su3_vector **psim,
    ks_param *ksp, 
    int num_offsets,
    quark_invert_control *qic,
    imp_ferm_links_t *fn
    )
{
  int i;
  char myname[] = "ks_multicg_offset_field_gpu";
  if(qic[0].relresid != 0.){
    printf("%s: GPU code does not yet support a Fermilab-type relative residual\n", myname);
    terminate(1);
  }

  int quda_parity = 0;
  if(qic[0].parity == EVEN){
    quda_parity = 0;
  }else if(qic[0].parity == ODD){
    quda_parity = 1;
  }else{
    printf("%s: Unrecognised parity\n", myname);
    terminate(2);
  }


  double* offset = (double*)malloc(num_offsets*sizeof(double));
  for(i=0; i<num_offsets; ++i) offset[i] = ksp[i].offset;

  double tmp;
  for(i=0; i<num_offsets; ++i){
    tmp = ksp[i].offset;
    offset[i] = tmp;
    printf("offset[%d] = %lf\n",i,offset[i]);
  }


  for(i=0; i<num_offsets; ++i){
    printf("ksp[%d].offset = %lf\n",i,ksp[i].offset);      
  }


  // The following arrays temporarily hold the final residuals
  double* residual = (double*)malloc(num_offsets*sizeof(double));
  double* fermilab_residual = (double*)malloc(num_offsets*sizeof(double));

  const int maxiter  = qic[0].max*qic[0].nrestart;
  int num_iters; // number of iterations taken
  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);

  // Work out the grid dimensions
  int dim[4] = {nx,ny,nz,nt};
  int const* ng = nodegeom();
  int node_geom[4] = {ng[0], ng[1], ng[2], ng[3]};

  printf("Calling qudaMultishiftInvert\n");


  const int quda_precision = 2;
  qudaMultishiftInvert(quda_parity,
                   dim,
		   node_geom,
                   PRECISION,
                   quda_precision,
                   num_offsets,
                   offset,
                   qic[0].resid,
                   qic[0].relresid,
                   maxiter,
                   fatlink,
                   longlink,
                   src,
		       (void **)psim,
                   residual,
                   fermilab_residual,
                   &num_iters);


  for(i=0; i<num_offsets; ++i){
    qic[i].final_rsq = residual[i]*residual[i];
    qic[i].final_relrsq = fermilab_residual[i]*fermilab_residual[i];
    qic[i].final_iters = num_iters;

    // check for convergence
    qic[i].converged = (residual[i] <= qic[0].resid) ? 1 : 0;

    // Cumulative residual. Not used in practice
    qic[i].size_r = 0.0;
    qic[i].size_relr = 0.0;
  }
                  
  free(offset); 
  free(residual);
  free(fermilab_residual);

  return num_iters;
}

