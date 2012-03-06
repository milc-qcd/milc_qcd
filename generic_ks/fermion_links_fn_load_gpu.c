/**************** fermion_links_fn_load_gpu.c **********************/
/* MILC Version 7 */

/* Originally Guochun Shi and Steve Gottlieb 2010 */
/* Revised Foley 2012 */

/* Entry points 

   load_fatlinks_gpu

*/

#include "generic_ks_includes.h"
#include "../include/info.h"

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

void  
load_fatlinks_gpu(info_t *info, su3_matrix *fat, ks_component_paths *p, su3_matrix *links)
{

  double path_coeff[6];
  path_coeff[0] = p->act_path_coeff.one_link;
  path_coeff[1] = p->act_path_coeff.naik;
  path_coeff[2] = p->act_path_coeff.three_staple;
  path_coeff[3] = p->act_path_coeff.five_staple;
  path_coeff[4] = p->act_path_coeff.seven_staple;
  path_coeff[5] = p->act_path_coeff.lepage;

  QudaFatLinkArgs_t fatlink_args;
  fatlink_args.su3_source = 0; // Cannot guarantee that the incoming field is an SU(3) gauge-field 
			       // Need a workaround for this
  fatlink_args.use_pinned_memory = 0;
 
  initialize_quda();

  qudaLoadFatLink(PRECISION, fatlink_args, path_coeff, links, fat);
  return;
}

/* fermion_links_fn_load_gpu.c */
