/******* d_congrad5_fn_gpu.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* GPU version of d_congrad5_fn_milc.c.  Can be compiled together with it. */
// Note that the restart criterion used by QUDA is different from 
// the restart criterion used by the MILC code

/* The following headers are supplied with the MILC code */
#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#define LOOPEND
#include "../include/loopend.h"
#include <string.h>

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/* The following headers are supplied with QUDA */
#include <milc_interface.h>

int ks_congrad_parity_gpu(su3_vector *t_src, su3_vector *t_dest, 
			  quark_invert_control *qic, Real mass,
			  imp_ferm_links_t *fn)
{

  char myname[] = "ks_congrad_parity_gpu";

  if(qic->relresid != 0.){
    printf("%s: GPU code does not yet support a Fermilab-type relative residual\n",myname);
    terminate(1);
  }
  
  int quda_parity = 0;
  if(qic->parity == EVEN){
    quda_parity = 0;
  }else if(qic->parity == ODD){
    quda_parity = 1;
  }else{
    printf("%s: Unrecognised parity\n",myname);
    terminate(2);
  }

  int dim[4] = {nx, ny, nz, nt};
  const int maxiter = qic->max*qic->nrestart;
  int num_iters;

  su3_matrix* fatlink = get_fatlinks(fn);
  su3_matrix* longlink = get_lnglinks(fn);

  int const* ng = nodegeom();
  int node_geom[4] = {ng[0], ng[1], ng[2], ng[3]};

  printf("Calling qudaInvert\n");
  double residual, fermilab_residual;

  const double restart_tolerance = 1e-3;

  qudaInvert(quda_parity, 
					   dim, 
						 node_geom,
					   PRECISION,
				                 qic->prec, 
						 mass,
						 qic->resid,
						 qic->relresid,
						 restart_tolerance,
						 maxiter,
						 fatlink, 
						 longlink,
						 t_src, 
						 t_dest,
             &residual,
			       &fermilab_residual, 
					   &num_iters);

  qic->final_rsq = residual*residual;
  qic->final_relrsq = fermilab_residual*fermilab_residual;
  qic->final_iters = num_iters;

  // check for convergence 
  qic->converged = (residual < qic->resid) ? 1 : 0;

  // Cumulative residual. Not used in practice 
  qic->size_r = 0.0;
  qic->size_relr = 0.0; 
 
  return num_iters;
}

  


