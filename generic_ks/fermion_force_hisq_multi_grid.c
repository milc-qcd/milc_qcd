/******************* fermion_force_hisq_multi_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

/* PLACEHOLDER */

#include "generic_ks_includes.h"
#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include <string.h>


#define GRID_MAX_NAIK 4
/* HISQ datatypes*/
typedef struct {
  int n_naiks;
  double eps_naik[GRID_MAX_NAIK];
  GRID_hisq_unitarize_group_t ugroup;
  GRID_hisq_unitarize_method_t umethod;
  double fat7_one_link;
  double fat7_three_staple;
  double fat7_five_staple;
  double fat7_seven_staple;
  double fat7_lepage;
  double asqtad_one_link;
  double asqtad_three_staple;
  double asqtad_five_staple;
  double asqtad_seven_staple;
  double asqtad_lepage;
  double asqtad_naik;
  double difference_one_link;
  double difference_naik;
} GRID_hisq_coeffs_t;
#define GRID_HISQ_COEFFS_ZERO \
  ((GRID_hisq_coeffs_t){1, GRID_EPS_NAIK_ZERO, GRID_UNITARIZE_U3, \
      GRID_UNITARIZE_RATIONAL, 0,0,0,0,0, 0,0,0,0,0,0, 0,0})

extern GRID_4Dgrid *grid_full;

static GRID_hisq_coeffs_t
assemble_action_coeffs(ks_action_paths *ap, int n_naiks, double eps_naik[]){
  static GRID_hisq_coeffs_t hc = GRID_HISQ_COEFFS_ZERO;
  hc.n_naiks = n_naiks;
  for(int i = 0; i < n_naiks; i++){
    hc.eps_naik[i] = eps_naik[i];
  }
  hc.ugroup              = UNITARIZATION_GROUP
  hc.umethod             = UNITARIZATION_METHOD
  hc.fat7_one_link       = ap->p1.
  hc.fat7_three_staple   = 
  hc.fat7_five_staple    = 
  hc.fat7_seven_staple   = 
  hc.fat7_lepage         = 
  hc.asqtad_one_link     = 
  hc.asqtad_three_staple = 
  hc.asqtad_five_staple  = 
  hc.asqtad_seven_staple = 
  hc.asqtad_lepage       = 
  hc.asqtad_naik         = 
  hc.difference_one_link = 
  hc.difference_naik     = 
}

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
fermion_force_multi_hisq_grid(info_t* info, int prec, Real eps, Real *residues, 
			      su3_vector **multi_x, int nterms,
			      fermion_links_t *fl)
{
  char myname[] = "fermion_force_multi_hisq_grid";
  GRID_info_t grid_info;

  // Setup copied from fermion_force_hisq_multi_quda.c
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
    node0_printf("%s: WARNING, precision requests not yet supported. Using %d.\n",
		 myname, MILC_PRECISION);
  }

  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  int n_orders_naik_current = n_order_naik_total;

  /* Create tables of uncorrected level 2 one-hop and three-hop
     coefficients, one for each rational function term, multiplied by
     the rational function residue */
  
  Real* one_hop_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));
  Real* three_hop_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));

  for(i=0; i<n_orders_naik_current; ++i){
    one_hop_coeff[i] = 2.0*residues[i];
    three_hop_coeff[i] = ap->p2.act_path_coeff.naik*2.0*residues[i];
  }

  /* Count terms with nonzero Naik epsilons */
  n_orders_naik_current = 0;
  for(int inaik=1; inaik<n_naiks; ++inaik) n_orders_naik_current += n_orders_naik[inaik];

  /* Create tables of the Naik epsilon corrections to the one-hop and
     three-hop coefficients */

  Real* one_hop_naik_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));
  Real* three_hop_naik_coeff = (Real*)malloc(n_orders_naik_current*sizeof(Real));

  i = 0;
  int n_naik_shift = n_orders_naik[0];
  for(int inaik=1; inaik<n_naiks; ++inaik){
    for(int j=0; j<n_orders_naik[inaik]; j++){
      one_hop_naik_coeff[i] = ap->p3.act_path_coeff.one_link *
	eps_naik[inaik] * 2.0*residues[n_naik_shift+j];
      three_hop_naik_coeff[i++] = ap->p3.act_path_coeff.naik *
	eps_naik[inaik] * 2.0*residues[n_naik_shift+j];
    }
    n_naik_shift += n_orders_naik[inaik];
  }

  /* Collect all the one-hop and three-hop coefficients for level 2
     smearing in a table with two parts.  The first contains all the
     unmodified one-term and three-term coefficients.  The second part
     contains the modifications of the one-hop and three-hop
     coefficients coming from nonzero Naik epsilons */

  double** coeff = (double**)malloc( (n_order_naik_total + n_orders_naik_current) *
				     sizeof(double*));
  for(int term=0; term<(n_order_naik_total+n_orders_naik_current); ++term)
    coeff[term] = (double*)malloc(2*sizeof(double));

  for(int term=0; term<n_order_naik_total; ++term){
    coeff[term][0] = one_hop_coeff[term];
    coeff[term][1] = three_hop_coeff[term];
  }
  
  for(int term=0; term<n_orders_naik_current; term++){
    coeff[term+n_order_naik_total][0] = one_hop_naik_coeff[term];
    coeff[term+n_order_naik_total][1] = three_hop_naik_coeff[term];
  }

  /* Collect the six path coefficients for the level 2 smearing
     in a table */
  double level2_coeff[6];
  write_path_coeffs_to_array(ap->p2, level2_coeff);
  /* Collect the six path coefficiens for the level 1 (fat7) smearing */
  double fat7_coeff[6];
  write_path_coeffs_to_array(ap->p1, fat7_coeff);
  
  /* Set HISQ reunitarization parameters for GRID */
  /* Compare with default values set in su3_mat_op.c */
  GRID_HisqParams_t params;

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

  anti_hermitmat* momentum = (anti_hermmat *)malloc(sites_on_node*4*sizeof(anti_hermitmat));

  GRID_HisqParamsInit(params);

  if(prec == 1)
    GRID_F3_hisq_force(&grid_info, MILC_PRECISION, prec, n_order_naik_total,
		       n_orders_naik_current, eps, coeff, multi_x, level2_coeff,
		       fat7_coeff, W_unitlink, V_link, U_link, momentum, grid_full);
  else
    GRID_D3_hisq_force(&grid_info, MILC_PRECISION, prec, n_order_naik_total,
		       n_orders_naik_current, eps, coeff, multi_x, level2_coeff,
		       fat7_coeff, W_unitlink, V_link, U_link, momentum, grid_full);


    GRID_D3_hisq_force(&grid_info, fl, residues, nterms, multi_x, n_orders_naik,
		       grid_full);

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

  for (int term=0; term<(n_order_naik_total+n_orders_naik_current); ++term)
    free(coeff[term]);
  free(coeff);
  free(one_hop_coeff);
  free(three_hop_coeff);
  free(one_hop_naik_coeff);
  free(three_hop_naik_coeff);

  return;
  
}

/* fermion_links_hisq_load_grid.c */
