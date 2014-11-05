/****************** fermion_force_asqtad_gpu.c *******************/
/* MIMD version 7 */
/* fermion force optimized for the Asqtad action 
 * GPU version
 * J. Foley 6/15/2014
 */

/* External entry points in this file
   
    eo_fermion_force_twoterms_site
   
 */

#include "generic_ks_includes.h"        /* definitions files and prototypes */
#include <string.h>

/* This routine is valid only for Asqtad, so requires the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

static void 
fermion_force_oprod_site(Real eps, Real weight1, Real weight2,
			 field_offset x1_off, field_offset x2_off,
			 su3_matrix* one_hop_oprod[4],
			 su3_matrix* three_hop_oprod[4])
{
  su3_vector *v[2];
  site* s;
  int i, j, dir;
  msg_tag* mtag[2];
  
  { // copy the quark-field information to su3_vector fields
    v[0] = (su3_vector*)malloc(sites_on_node*sizeof(su3_vector));
    v[1] = (su3_vector*)malloc(sites_on_node*sizeof(su3_vector));

    if(v[0] == NULL) printf("fermion_force_oprod_site: v[0] not allocated\n");
    if(v[1] == NULL) printf("fermion_force_oprod_site: v[1] not allocated\n");  

    FORALLSITES(i,s){
      v[0][i] = *(su3_vector*)F_PT(s,x1_off);
      v[1][i] = *(su3_vector*)F_PT(s,x2_off);
    }
  }

  for(dir=XUP; dir<=TUP; ++dir){
    memset(one_hop_oprod[dir], 0, sites_on_node*sizeof(su3_matrix));
    memset(three_hop_oprod[dir], 0, sites_on_node*sizeof(su3_matrix));
  }

  double** combined_coeff;
  combined_coeff = (double**)malloc(2*sizeof(double*));
  combined_coeff[0] = (double*)malloc(2*sizeof(double));
  combined_coeff[1] = (double*)malloc(2*sizeof(double));

  combined_coeff[0][0] = 2*eps*weight1;
  combined_coeff[0][1] = 2*eps*weight2;
  combined_coeff[1][0] = 2*eps*weight1;
  combined_coeff[1][1] = 2*eps*weight2;

  void* oprod[2] = {one_hop_oprod, three_hop_oprod};

  qudaComputeOprod(PRECISION, 2, combined_coeff, v, oprod);

  free(combined_coeff[0]);
  free(combined_coeff[1]);
  free(combined_coeff);

  // Cleanup
  free(v[0]);
  free(v[1]);
}     

void 
eo_fermion_force_twoterms_site_gpu(Real eps, Real weight1, Real weight2, 
				   field_offset x1_off, field_offset x2_off,
				   int prec, fermion_links_t *fl){
  int i, j, dir;
  site* s;

  initialize_quda();

  double act_path_coeff[6];
#ifdef FFTIME
  double dtime;
  dtime=-dclock();
#endif  

  ks_action_paths *ap = get_action_paths(fl);

  act_path_coeff[0] = ap->p.act_path_coeff.one_link;
  act_path_coeff[1] = ap->p.act_path_coeff.naik;
  act_path_coeff[2] = ap->p.act_path_coeff.three_staple;
  act_path_coeff[3] = ap->p.act_path_coeff.five_staple;
  act_path_coeff[4] = ap->p.act_path_coeff.seven_staple;
  act_path_coeff[5] = ap->p.act_path_coeff.lepage;


  su3_matrix* one_link_oprod[4];
  su3_matrix* three_link_oprod[4];

  for(dir=XUP; dir<=TUP; ++dir){
    one_link_oprod[dir] = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));
    three_link_oprod[dir] = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));
  }


#ifdef FFTIME
  double otime;
  otime = -dclock();
#endif
  // compute the one-link and three-link outer products
  fermion_force_oprod_site(eps, weight1, weight2,
      x1_off, x2_off, one_link_oprod, three_link_oprod);

#ifdef FFTIME
  otime += dclock();
  node0_printf("OPTIME: time = %e (asqtad3)\n",otime);
#endif

  Real *momentum = (Real*)malloc(sites_on_node*4*sizeof(anti_hermitmat));
  Real *gauge = (Real*)malloc(sites_on_node*4*sizeof(su3_matrix));
  memset(momentum, 0, sites_on_node*4*sizeof(anti_hermitmat));

  // Populate gauge field 
  FORALLSITES(i,s){
    for(dir=0; dir<4; ++dir){
      for(j=0; j<18; ++j){
        gauge[(4*i + dir)*18 + j] = ((Real*)(&(s->link[dir])))[j];
      }
    }
  }


  // Call the quda function
  qudaAsqtadForce(PRECISION, act_path_coeff, 
      (const void* const*)one_link_oprod,
      (const void* const*)three_link_oprod,
      gauge , momentum);

  Real mtemp;
  // Copy result to site structure
  FORALLSITES(i,s){
    for(dir=0; dir<4; ++dir){
      for(j=0; j<10; ++j){
        mtemp = *(momentum + (4*i+ dir)*10 + j);
        *((Real*)(&(s->mom[dir])) + j) += *(momentum + (4*i+ dir)*10 + j);
      }
    }
  }

  for(dir=XUP; dir<=TUP; ++dir){
    free(one_link_oprod[dir]);
    free(three_link_oprod[dir]);
  }
  free(momentum);
  free(gauge);

#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME: time = %e (asqtad3)\n",dtime);
#endif
  return;
}


