/******************* fermion_links_hisq_load_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

#include "generic_ks_includes.h"
#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include <string.h>

void
load_fatlinks_gpu(info_t *info, su3_matrix *fat, ks_component_paths *p, su3_matrix *thin_links){
  if(MILC_PRECISION == 1)
    printf("ERROR. load_fatlinks_grid_F\n");
  //    load_fatlinks_grid_F(info, fat, p, thin_links);
  else
    load_fatlinks_grid_D(info, fat, p, thin_links);
}    
  
void
load_hisq_aux_links_gpu(info_t *info, ks_action_paths_hisq *ap,
			hisq_auxiliary_t *aux, su3_matrix *links){
  if(MILC_PRECISION == 1)
    load_hisq_aux_links_grid_F(info, ap, aux, links);
  else
    load_hisq_aux_links_grid_D(info, ap, aux, links);

}

void load_fatlonglinks_gpu(info_t *info, su3_matrix *fatlinks, su3_matrix *longlinks,
			   ks_component_paths *p, su3_matrix *links){
  if(MILC_PRECISION == 1)
    load_fatlonglinks_grid_F(info, fatlinks, longlinks, p, links);
  else
    load_fatlonglinks_grid_D(info, fatlinks, longlinks, p, links);
}

/* fermion_links_hisq_load_grid.c */
