/****** fermion_force_hisq_multi_quda.c  -- ******************/
/* MIMD version 7 */
/* Multisource fermion force wrapper for QUDA */
/* May 5, 2024 C. DeTar: Split from fermion_force_hisq_multi.c */

#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#define LOOPEND
#include "../include/loopend.h"
#include <string.h>
#include "../include/generic.h"

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

static void
write_path_coeffs_to_array(ks_component_paths paths, double array[6])
{
  array[0] = paths.act_path_coeff.one_link; 
  array[1] = paths.act_path_coeff.naik;
  array[2] = paths.act_path_coeff.three_staple;
  array[3] = paths.act_path_coeff.five_staple;
  array[4] = paths.act_path_coeff.seven_staple;
  array[5] = paths.act_path_coeff.lepage;
  return;
}

void 
fermion_force_multi_hisq_quda(info_t* info, int prec, Real eps, Real *residues, 
			      su3_vector **multi_x, int nterms,
			      fermion_links_t *fl)
{
  // Get the terms we need from the fermion links structure */
  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  hisq_auxiliary_t *aux = get_hisq_auxiliary(fl);
  int n_naiks = fermion_links_get_n_naiks(fl);
  double *eps_naik = fermion_links_get_eps_naik(fl);
  su3_matrix *U_link = aux->U_link;
  su3_matrix *V_link = aux->V_link;
  su3_matrix *W_unitlink = aux->W_unitlink;
  int i;
  site *s;
  
  if(prec != MILC_PRECISION){
    node0_printf("fermion_force_multi_hisq_quda: WARNING, precision requests not yet supported. Using %d.\n",
		 MILC_PRECISION);
  }
  initialize_quda();
  
  int n_orders_naik_current = n_order_naik_total;

  Real* one_hop_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));
  Real* three_hop_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));

  for(i=0; i<n_orders_naik_current; ++i){
    one_hop_coeff[i] = 2.0*residues[i];
    three_hop_coeff[i] = ap->p2.act_path_coeff.naik*2.0*residues[i];
  }

  n_orders_naik_current = 0;
  for(int inaik=1; inaik<n_naiks; ++inaik) n_orders_naik_current += n_orders_naik[inaik];

  Real* one_hop_naik_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));
  Real* three_hop_naik_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));

  i = 0;
  int n_naik_shift = n_orders_naik[0];
  for(int inaik=1; inaik<n_naiks; ++inaik){
    for(int j=0; j<n_orders_naik[inaik]; j++){
      one_hop_naik_coeff[i] = ap->p3.act_path_coeff.one_link*eps_naik[inaik]*2.0*residues[n_naik_shift+j];
      three_hop_naik_coeff[i++] = ap->p3.act_path_coeff.naik*eps_naik[inaik]*2.0*residues[n_naik_shift+j];
    }
    n_naik_shift += n_orders_naik[inaik];
  }

  double** coeff = (double**)malloc( (n_order_naik_total + n_orders_naik_current) *sizeof(double*));
  for(int term=0; term<(n_order_naik_total+n_orders_naik_current); ++term) coeff[term] = (double*)malloc(2*sizeof(double));

  for(int term=0; term<n_order_naik_total; ++term){
    coeff[term][0] = one_hop_coeff[term];
    coeff[term][1] = three_hop_coeff[term];
  }
  
  for(int term=0; term<n_orders_naik_current; term++){
    coeff[term+n_order_naik_total][0] = one_hop_naik_coeff[term];
    coeff[term+n_order_naik_total][1] = three_hop_naik_coeff[term];
  }

  double level2_coeff[6];
  write_path_coeffs_to_array(ap->p2, level2_coeff);
  double fat7_coeff[6];
  write_path_coeffs_to_array(ap->p1, fat7_coeff);
  
  /* Set HISQ reunitarization parameters for QUDA */
  /* Compare with default values set in su3_mat_op.c */
  QudaHisqParams_t params;

#ifndef HISQ_REUNIT_ALLOW_SVD
  params.reunit_allow_svd = 1;
#else
  params.reunit_allow_svd = HISQ_REUNIT_ALLOW_SVD;
#endif

#ifndef HISQ_REUNIT_SVD_ONLY
  params.reunit_svd_only = 0;
#else
  params.reunit_svd_only = HISQ_REUNIT_SVD_ONLY;
#endif

#ifndef HISQ_REUNIT_SVD_ABS_ERROR
  params.reunit_svd_abs_error = 1e-8;
#else
  params.reunit_svd_abs_error = HISQ_REUNIT_SVD_ABS_ERROR;
#endif

#ifndef HISQ_REUNIT_SVD_REL_ERROR
  params.reunit_svd_rel_error = 1e-8;
#else
  params.reunit_svd_rel_error = HISQ_REUNIT_SVD_REL_ERROR;
#endif

#ifndef HISQ_FORCE_FILTER
  params.force_filter = 5e-5;
#else
  params.force_filter = HISQ_FORCE_FILTER;
#endif

  Real* momentum = (Real*)qudaAllocatePinned(sites_on_node*4*sizeof(anti_hermitmat));

  qudaHisqParamsInit(params);
  qudaHisqForce(MILC_PRECISION, n_order_naik_total, n_orders_naik_current, eps, coeff, (void**)multi_x,
                level2_coeff, fat7_coeff, W_unitlink, V_link, U_link, momentum);

  // append result
  FORALLSITES_OMP(i,s,){
    for(int dir=0; dir<4; ++dir){
      for(int j=0; j<10; ++j){
	*((Real*)(&(s->mom[dir])) + j) += *(momentum + (4*i+ dir)*10 + j);
      }
    }
  } END_LOOP_OMP

  // free memory
  qudaFreePinned(momentum);

  for (int term=0; term<(n_order_naik_total+n_orders_naik_current); ++term) free(coeff[term]);
  free(coeff);
  free(one_hop_coeff);
  free(three_hop_coeff);
  free(one_hop_naik_coeff);
  free(three_hop_naik_coeff);

  return;
} //fn_fermion_force_multi_hisq_wrapper_mx_gpu

