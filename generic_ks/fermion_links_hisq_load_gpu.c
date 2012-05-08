/**************** fermion_links_hisq_load_gpu.c **********************/
/* MILC Version 7 */

/* Foley 2012 */

/* Entry points 

   load_hisq_aux_links_gpu

*/

#include "generic_ks_includes.h"
#include "../include/info.h"
#include <string.h>

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

void 
load_hisq_aux_links_gpu(info_t *info, ks_action_paths_hisq *ap, 
			hisq_auxiliary_t *aux, su3_matrix *links)
{
  char myname[] = "load_hisq_aux_links_gpu";
  
  if(ap == NULL){
    printf("%s(%d): KS action paths not initialized\n", myname, this_node);
  }

  // load U links
  memcpy(aux->U_link, links, 4*sizeof(su3_matrix)*sites_on_node);
  

  double path_coeff[6];
  path_coeff[0] = ap->p1.act_path_coeff.one_link;
  path_coeff[1] = ap->p1.act_path_coeff.naik;
  path_coeff[2] = ap->p1.act_path_coeff.three_staple;
  path_coeff[3] = ap->p1.act_path_coeff.five_staple;
  path_coeff[4] = ap->p1.act_path_coeff.seven_staple;
  path_coeff[5] = ap->p1.act_path_coeff.lepage;

  QudaFatLinkArgs_t fatlink_args;
  fatlink_args.su3_source = 1; // Is the incoming field an SU(3) gauge field? 
			       // If so, run SU(3) optimized QUDA code.
  fatlink_args.use_pinned_memory = 0; // Use page-locked memory in QUDA?

  initialize_quda();

  // Right now, if aux->V_link == NULL 
  // the level1 fat link is not copied from the GPU back to the CPU.
  qudaLoadUnitarizedLink(PRECISION, fatlink_args, path_coeff, aux->U_link, aux->V_link, aux->W_unitlink);

  return;
}

/* fermion_links_hisq_load_gpu.c */
