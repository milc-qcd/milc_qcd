/******* d_bicgilu_cl_qphixj.c - BiCGstab-ILU for  clover fermions on Xeon Phi ****/
/* MIMD version 7 */

/* 1/16/17 C DeTar */

#include "generic_clover_includes.h"
#include "../include/generic_qphixj.h"
#include "../include/gammatypes.h"
#include "../include/openmp_defs.h"
#include <assert.h>

#ifdef SCHROED_FUN
# error "QPhiXJ bicgstab does not support SCHROED_FUN"
#endif

/* Compute clover propagator from the source src.  Answer in dest. */
int bicgilu_cl_field_qphixj ( // Return value is number of iterations taken 
			     wilson_vector *src,
			     wilson_vector *dest,
			     quark_invert_control *qic,
			     void *dmp  		   // parameters defining the Dirac matrix
			      )			   
{

  char myname[] = "bicgilu_cl_field_inner";
  int flag = qic->start_flag;  // O: use a zero initial guess;
  // 1: use dest as starting guess
  dirac_clover_param *dcp = (dirac_clover_param *)dmp; 
  Real kappa  = dcp->Kappa;
  Real clov_c = dcp->Clov_c;
  Real u0     = dcp->U0; 
  Real CKU0   = kappa*clov_c/(u0*u0*u0);

  /* We don't do precision conversions yet */
  assert(MILC_PRECISION == QPHIXJ_PrecisionInt);
  
  /* Reset statistics */
  qic->size_r = 0;
  qic->size_relr = 0;
  qic->final_rsq = 0;
  qic->final_relrsq = 0;
  qic->final_iters = 0;
  qic->final_restart = 0;

  /* Handle trivial case */
  if(kappa == 0.){
    copy_wv_field(dest,src);
    return 0;
  }

  /* Clear dest if starting with a zero guess */
  if(flag == 0) {
#ifdef CG_DEBUG
    node0_printf("dest_0=0\n");fflush(stdout);
#endif
    clear_wv_field(dest);
  }
  
  // Clover term and its inverse
  
  clover *milc_clov = gen_clov;
  if(milc_clov == NULL){
    printf("%s(%d): milc_clov == NULL\n", myname, this_node);	
    terminate(1);
  }
  
  double dtime2 = -dclock();
  /* Compute R_e and R_o and put in "clov" and "clov_diag" */
  compute_clov(milc_clov, CKU0); // CKU0 = coefficient of the clover term

  /* Invert R_o only, leaving R_e on even sites and 1/R_o on odd sites 
     in "clov" and "clov_diag" */
  compute_clovinv(milc_clov, ODD);

  dtime2 += dclock();
#ifdef CGTIME
  node0_printf("Time to compute clover term %.2e\n", dtime2);
#endif

#ifdef CGTIME
  double dtime = -dclock();
#endif
  
  dtime2 = -dclock();
  /* LU transformation of the source -- result returned in place */
  /* r = L^(-1)*r. See cl_solver_utilities.c for notation. */
  /* Set up the LU decomposition */
  int is_startedo = 0, is_startede = 0;
  int num_iters = 0;
  wilson_vector *r = create_wv_field();
  wilson_vector *tmp = create_wv_field();
  msg_tag *tago[8],*tage[8];

  /* now we copy src to temporary */
  copy_wv_field(r, src);

  /* src = L^(-1)*src */
  Real size_src = ilu_xfm_source(NULL, r, tmp, kappa, &is_startede, tage);
  Real size_src2 = size_src*size_src;

  dtime2 += dclock();
#ifdef CGTIME
  node0_printf("Time for RB preconditioning %.2e\n", dtime2);
#endif

  /* Provision for trivial solution */
  if(size_src2 == 0.0){
    clear_wv_field(dest);
    destroy_wv_field(r);
    destroy_wv_field(tmp);
    is_startede = is_startedo = 0;
    
    cleanup_tmp_links();
    cleanup_dslash_wtemps();
    return 0;
  }
  
  /* Initial guess */
  /* set dest = src on odd sites.  
     Required for reconstructing the odd checkerboard solution */
  int i; site *s;
  FORODDFIELDSITES_OMP(i, ) {
    copy_wvec( &(r[i]), &(dest[i]) );
  } END_LOOP_OMP;

  /* Solve the preconditioned problem on even sites */
  /* Precision conversion is possible */
  dtime2 = -dclock();
  if(qic->prec == 1)
    bicgilu_cl_qphixj_inner_F(milc_clov, kappa, r, dest, qic);
  else
    bicgilu_cl_qphixj_inner_D(milc_clov, kappa, r, dest, qic);

  dtime2 += dclock();
#ifdef CGTIME
  node0_printf("Time in bicbilu_cl_qphixj_inner_P %.2e\n", dtime2);
#endif

  dtime2 = -dclock();
  /* Reconstruct the solution on odd sites */
  ilu_xfm_dest(dest, tmp, kappa, &is_startedo, tago);
  
  dtime2 += dclock();
#ifdef CGTIME
  node0_printf("Time to reconstruct odd-site solution %.2e\n", dtime2);
#endif

  num_iters = qic->final_iters;
  
  /* Clean up */
  destroy_wv_field(r);
  destroy_wv_field(tmp);

#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    static const char *qphixj_prec[2] = {"F", "D"};
    if(num_iters==0)
      printf("BICGILU: NO iterations taken size_r= %.2e rel %.2e\n",
	     qic->final_rsq, qic->final_relrsq);
#ifdef CGTIME
    else{
      printf("CGTIME: time = %.2e (bicgilu qphixj %s) size_r= %.2e relr= %.2e iters= %d MF = %.1f\n",
	     dtime,qphixj_prec[qic->prec-1],qic->final_rsq,qic->final_relrsq,qic->final_iters,
	     (double)8742*qic->final_iters*even_sites_on_node/(dtime*1e6));
    }
#endif
    fflush(stdout);
  }    
  
  for(int i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;
  cleanup_tmp_links();
  cleanup_dslash_wtemps();


  return num_iters;  

} /* bicgilu_cl_field_qphixj */

