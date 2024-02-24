/**************** fermion_links_fn_load_gpu.c **********************/
/* MILC Version 7 */

/* UNDER CONSTRUCTION */

/* NOTE: This code is actually an include file for fermion_links_fn_load_grid_F.c
   and fermion_links_fn_load_grid_D.c, so any edits should be consistent with this
   purpose. */

/* Entry points (must be redefined to precision-specific names) 

   LOAD_FATLINKS_GPU
   LOAD_FN_LINKS_GPU
   LOAD_HISQ_AUX_LINKS_GPU
*/

/* Redefinitions according to requested precision */

#if ( GRID_PrecisionInt == 1 )

#define LOAD_FATLINKS_GPU        load_fatlinks_grid_F
#define LOAD_FATLONGLINKS_GPU    load_fatlonglinks_grid_F
#define LOAD_HISQ_AUX_LINKS_GPU  load_hisq_aux_links_grid_F
#define MYREAL GRID_F_Real

#else

#define LOAD_FATLINKS_GPU       load_fatlinks_grid_D
#define LOAD_FATLONGLINKS_GPU   load_fatlonglinks_grid_D
#define LOAD_HISQ_AUX_LINKS_GPU load_hisq_aux_links_grid_D
#define MYREAL GRID_D_Real

#endif

#include "generic_ks_includes.h"
#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include <string.h>

extern GRID_4Dgrid *grid_full;

void  
LOAD_FATLINKS_GPU(info_t *info, su3_matrix *fat, ks_component_paths *p, su3_matrix *thin_links)
{
  char myname[] = "load_fatlinks_gpu";

  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  double path_coeff[6];
  path_coeff[0] = p->act_path_coeff.one_link;
  path_coeff[1] = p->act_path_coeff.naik;
  path_coeff[2] = p->act_path_coeff.three_staple;
  path_coeff[3] = p->act_path_coeff.five_staple;
  path_coeff[4] = p->act_path_coeff.seven_staple;
  path_coeff[5] = p->act_path_coeff.lepage;

  GRID_info_t grid_info;
  printf("Calling GRID_hisq_links\n"); fflush(stdout);
  GRID_hisq_links(&grid_info, path_coeff, fat, NULL, thin_links, grid_full);
  dumpmat(fat);

  /* Fatlinks */
  info->final_flop = 61632.*volume/numnodes();
  if( p->act_path_coeff.three_staple == 0.0 &&
      p->act_path_coeff.lepage == 0.0 &&
      p->act_path_coeff.five_staple == 0.0)
    info->final_flop = 72.*volume/numnodes();
  /* Longlinks */
  info->final_flop += 1728.*volume/numnodes();
}

void
LOAD_HISQ_AUX_LINKS_GPU(info_t *info, ks_action_paths_hisq *ap,
			hisq_auxiliary_t *aux, su3_matrix *links)
{
  char myname[] = "load_hisq_aux_links_gpu";

  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  if(ap == NULL){
    printf("%s(%d): KS action paths not initialized\n", myname, this_node);
  }

  // load U links (is this really necessary since we have extracted "links" already?)
  memcpy(aux->U_link, links, 4*sizeof(su3_matrix)*sites_on_node);

  double path_coeff[6];
  path_coeff[0] = ap->p1.act_path_coeff.one_link;
  path_coeff[1] = ap->p1.act_path_coeff.naik;
  path_coeff[2] = ap->p1.act_path_coeff.three_staple;
  path_coeff[3] = ap->p1.act_path_coeff.five_staple;
  path_coeff[4] = ap->p1.act_path_coeff.seven_staple;
  path_coeff[5] = ap->p1.act_path_coeff.lepage;

  GRID_info_t grid_info;
  GRID_hisq_aux_links(&grid_info, path_coeff, aux->U_link,
		      aux->V_link, aux->W_unitlink, grid_full);
  /*
    The above equates to
    - load_V_from_U: 61632 flops
    - load_Y_from_V: (as CPU code: presently not counted)
    - load_W_from_Y: (as CPU code: presently not counted)
  */
  info->final_flop = 61632.*volume/numnodes();
  if( ap->p1.act_path_coeff.three_staple == 0.0 &&
      ap->p1.act_path_coeff.lepage == 0.0 &&
      ap->p1.act_path_coeff.five_staple == 0.0)
    info->final_flop = 72.*volume/numnodes();

  return;
}

void LOAD_FATLONGLINKS_GPU(info_t *info, su3_matrix *fat, su3_matrix *lng, ks_component_paths *p,
			   su3_matrix *thin_links)
{
  char myname[] = "load_fatlonglinks_gpu";

  if(! grid_initialized()){
    node0_printf("%s: FATAL Grid has not been initialized\n", myname);
    terminate(1);
  }

  double path_coeff[6];
  path_coeff[0] = p->act_path_coeff.one_link;
  path_coeff[1] = p->act_path_coeff.naik;
  path_coeff[2] = p->act_path_coeff.three_staple;
  path_coeff[3] = p->act_path_coeff.five_staple;
  path_coeff[4] = p->act_path_coeff.seven_staple;
  path_coeff[5] = p->act_path_coeff.lepage;

  double dtime = -dclock();
  GRID_info_t grid_info;
  printf("Calling GRID_hisq_links2\n"); fflush(stdout);
  GRID_hisq_links(&grid_info, path_coeff, fat, lng, thin_links, grid_full);

#if 0
  if(want_back)
    load_fn_backlinks(fn);
  else
    destroy_fn_backlinks(fn);
#endif

  dtime += dclock();
  info->final_sec = dtime;
}

/* fermion_links_fn_load_grid.c */
