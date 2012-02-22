/******* d_congrad5_fn.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  

   Previously called d_congrad5_fn_fewsums.c and d_congrad5_fn_tmp.c
   Previous functionality now appears in d_congrad5_fn_milc.c.
*/

#include "generic_ks_includes.h"
#include "../include/fermion_links.h"


/* API for field arguments */

int ks_congrad_field( su3_vector *src, su3_vector *dest, 
		      quark_invert_control *qic, Real mass,
		      imp_ferm_links_t *fn)
{
  int iters = 0;
  int parity = qic->parity;

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_parity(src, dest, qic, mass, fn);
    report_status(qic);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_parity(src, dest, qic, mass, fn);
    report_status(qic);
  }

  qic->parity = parity;
  return iters;
}

/* API for field arguments.  This one never uses the GPU. */

int ks_congrad_field_cpu( su3_vector *src, su3_vector *dest, 
			  quark_invert_control *qic, Real mass,
			  imp_ferm_links_t *fn)
{
  int iters = 0;
  int parity = qic->parity;

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_parity_cpu(src, dest, qic, mass, fn);
    report_status(qic);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_parity_cpu(src, dest, qic, mass, fn);
    report_status(qic);
  }

  qic->parity = parity;
  return iters;
}

/* API for site arguments */

int ks_congrad_site( field_offset src, field_offset dest, 
		     quark_invert_control *qic, Real mass,
		     imp_ferm_links_t *fn)
{
  int i;
  site *s;
  int iters = 0;
  su3_vector *t_src, *t_dest;
  int parity = qic->parity;

  /* Map src and dest from site to field of correct precision */
  
  t_src  = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  t_dest = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(t_src == NULL || t_dest == NULL){
    printf("ks_congrad_site(%d): No room for temporaries\n",this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    t_src[i]  = *((su3_vector *)F_PT(s,src) );
    t_dest[i] = *((su3_vector *)F_PT(s,dest));
  }

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_parity(t_src, t_dest, qic, mass, fn );
    report_status(qic);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_parity(t_src, t_dest, qic, mass, fn );
    report_status(qic);
  }

  /* Map solution to site structure */

  FORALLSITES(i,s){
    *((su3_vector *)F_PT(s,dest)) = t_dest[i];
  }

  qic->parity = parity;

  free(t_src); free(t_dest);

  return iters;
}

/* Traditional MILC API for site arguments and no relative residual test */

int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int prec,
		int parity, Real *final_rsq,
		imp_ferm_links_t *fn){
  int iters;
  quark_invert_control qic;

  /* Pack structure */
  qic.prec      = prec;  /* Currently ignored */
  qic.min       = 0;
  qic.max       = niter;
  qic.nrestart  = nrestart;
  qic.parity    = parity;
  qic.start_flag = 0;
  qic.nsrc      = 1;
  qic.resid     = sqrt(rsqmin);
  qic.relresid  = 0;     /* Suppresses this test */

  /* Solve the system */
  iters = ks_congrad_site( src, dest, &qic, mass, fn );

  /* Unpack the results */
  *final_rsq    = qic.final_rsq;
  total_iters += iters;
  return iters;
}
