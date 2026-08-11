/**************** fermion_links_fn_load_quda.c **********************/
/* MILC Version 7 */

/* qudaLoadKSLink code originally Guochun Shi and Steve Gottlieb 2010 */
/* Revised Foley 2012 */

/* Entry points 

   load_fatlinks_gpu
   load_fatlonglinks_gpu
   load_hisq_aux_links_gpu
*/

#include "generic_ks_includes.h"
#include "../include/info.h"
#include "../include/generic_quda.h"
#include <string.h>

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
 
  initialize_quda();

  qudaLoadKSLink(MILC_PRECISION, fatlink_args, path_coeff, links, fat, NULL);

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
load_fatlonglinks_gpu(info_t *info, su3_matrix *fatlinks, su3_matrix *longlinks,
		      ks_component_paths *p, su3_matrix *links)
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
 
  initialize_quda();

  // qudaLoadUnitarizedLink(MILC_PRECISION, fatlink_args, path_coeff, links, fatlinks, longlinks, NULL);
  qudaLoadKSLink(MILC_PRECISION, fatlink_args, path_coeff, links, fatlinks, longlinks);

  /* Release QUDA's device memory pool before the solver's loadGaugeQuda runs.

     The HISQ link build leaves large blocks in QUDA's pool (three
     extended gauge fields plus one unextended, freed to the cache rather than
     to the driver). loadGaugeQuda then builds the gauge tower, and its small
     per-dimension ghost arrays take those blocks WHOLE, whic results in some
     device memory stranded within a live allocation where no later flush can
     reach.

     OFF BY DEFAULT. Enable at runtime by setting MILC_FLUSH_QUDA_POOL_AFTER_LINKS=1.
     Note that there is a performance cost, and this should not be used in
     gauge generation where the links are rebuilt on every MD step and
     flushing the pool would discard the solver's cached working set
     hundreds of times per trajectory.
   
     This is the LAST of the two computeKSLinkQuda calls a HISQ build makes (the
     other is the aux U->V->W build in load_hisq_aux_links_gpu below, whose
     blocks this call has already reused), so one flush here catches both. */
  {
    static int flush_quda_pool = -1;
    if(flush_quda_pool < 0){
      char *env = getenv("MILC_FLUSH_QUDA_POOL_AFTER_LINKS");
      flush_quda_pool = (env != NULL && strcmp(env, "0") != 0);
      if(flush_quda_pool)
	node0_printf("load_fatlonglinks_gpu: flushing the QUDA device memory pool after the link build\n");
    }
    if(flush_quda_pool) flushPoolQuda(QUDA_MEMORY_DEVICE);
  }

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
load_hisq_aux_links_gpu(info_t *info, ks_action_paths_hisq *ap,
			hisq_auxiliary_t *aux, su3_matrix *links)
{
  char myname[] = "load_hisq_aux_links_gpu";

  if(ap == NULL){
    printf("%s(%d): KS action paths not initialized\n", myname, this_node);
  }

  // load U links (is this really necessary since we have extracted "links" already?)
  copy_m_array_field(aux->U_link, links, 4);

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

  initialize_quda();

  // Right now, if aux->V_link == NULL
  // the level1 fat link is not copied from the GPU back to the CPU.
  qudaLoadUnitarizedLink(MILC_PRECISION, fatlink_args, path_coeff, aux->U_link, aux->V_link, aux->W_unitlink);

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

void load_fn_links_gpu(info_t *info, fn_links_t *fn, ks_action_paths *ap,
		       su3_matrix *links, int want_back)
{
  ks_component_paths *p = &ap->p;
  double dtime = -dclock();

  load_fatlonglinks_gpu(info, fn->fat, fn->lng, p, links);

  if(want_back)
    load_fn_backlinks(fn);
  else
    destroy_fn_backlinks(fn);

  dtime += dclock();
  info->final_sec = dtime;
}

/* fermion_links_fn_load_quda.c */
