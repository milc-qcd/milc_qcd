#ifndef _FERMION_LINKS_MILC_H
#define _FERMION_LINKS_MILC_H
/******************** fermion_links_milc.h **********************************
*									*
*  Definitions for fermion links
*  MIMD version 7 							*
*									*
*/

#include "../include/ks_action_paths.h"
#include "../include/imp_ferm_links.h"
#include "../include/su3.h"
#include "../include/info.h"

/*********************************************************************/
/* MILC fermion links structures                                     */
/*********************************************************************/

/* 

   Structure nesting for both FN and EO actions

   milc_fm_links_t
     fm_ap_links_t *fm_ap
       ks_action_paths *ap
       fm_links_t *fm
     fm_ap_links_t *fm_ap_du0
       ks_action_paths *ap
       imp_ferm_links_t *fm

       (Note that in imp_ferm_links.h imp_ferm_links_t is equated
        to fn_links_t or eo_links_t, depending on the action)

   Structure nesting for HISQ actions

   milc_hisq_links_t
     hisq_links_t *hisq
       ks_action_paths_hisq *ap
       hisq_auxiliary_t *aux
         su3_matrix *U_link;
         su3_matrix *V_link;
	 su3_matrix *Y_unitlink; 
	 su3_matrix *W_unitlink; 
       fn_links_t *fn[MAX_NAIK]
       fn_links_t *fn_deps

       (Here fn_links_t is the only option.)
 */

/*********************************************************************/
typedef struct {
  ks_action_paths *ap;     // Paths defining the action
  imp_ferm_links_t *fm;    // The usual links
} fm_ap_links_t;

typedef struct {
  fm_ap_links_t *fm_ap;
  fm_ap_links_t *fm_ap_du0; /* Derivative of links wrto u0 */
} milc_fm_links_t;

/*********************************************************************/
typedef struct {
  int phases_in;          // track KS phases in the V and Y links
  int WeqY;               // true if W = Y
  su3_matrix *U_link;     // original gauge matrices, stored as four fields
  su3_matrix *V_link;     // first iteration of fattening
  su3_matrix *Y_unitlink; // unitary projection of V_link, U(3)
  su3_matrix *W_unitlink; // special unitary projection of Y_link, SU(3)
} hisq_auxiliary_t;

typedef struct {
  ks_action_paths_hisq *ap;
  hisq_auxiliary_t *aux;          // Intermediate links needed for fermion force
  imp_ferm_links_t *fn[MAX_NAIK]; // Table of links depending on epsilon
  imp_ferm_links_t *fn_deps;      // Derivative of links wrto epsilon.
} hisq_links_t;

typedef struct {
  hisq_links_t *hisq;
} milc_hisq_links_t;


/********************************************************************/
/* Fermion links routines */
/********************************************************************/

/* fermion_links_eo_load_milc.c */
void load_eo_links(info_t *info, eo_links_t *eo, ks_action_paths *ap, 
		   su3_matrix *links, int want_back);

/* fermion_links_fn_load_milc.c */
void load_fn_links(info_t *info, fn_links_t *fm, ks_action_paths *ap, 
		   su3_matrix *links, int want_back);

void
load_lnglinks(info_t *info, su3_matrix *lng, ks_component_paths *p,
	      su3_matrix *links );

void
load_fatlinks_cpu(info_t *info, su3_matrix *fat, ks_component_paths *p, 
		  su3_matrix *links);

void
load_fatlinks_gpu(info_t *info, su3_matrix *fat, ks_component_paths *p, 
		  su3_matrix *links);

void 
load_hisq_aux_links_gpu(info_t *info, ks_action_paths_hisq *ap, 
                        hisq_auxiliary_t *aux, su3_matrix *links);


#ifdef USE_FL_GPU
#define load_fatlinks load_fatlinks_gpu
#define load_hisq_aux_links load_hisq_aux_links_gpu
#else
#define load_fatlinks load_fatlinks_cpu
#define load_hisq_aux_links load_hisq_aux_links_cpu
#endif

/* fermion_links_hisq_load_milc.c */
void create_hisq_links_milc(info_t *info, fn_links_t **fn, fn_links_t **fn_deps,
			  hisq_auxiliary_t **aux, ks_action_paths_hisq *ap, 
			  su3_matrix *links, int want_deps, int want_back);

void destroy_hisq_links_milc(ks_action_paths_hisq *ap, hisq_auxiliary_t *aux, 
			     fn_links_t **fn, fn_links_t *fn_deps);


//hisq_auxiliary_t *create_hisq_auxiliary_t(ks_action_paths_hisq *ap,
//					  su3_matrix *links);
void destroy_hisq_auxiliary_t(hisq_auxiliary_t *aux);

/* fermion_links_fn_twist_milc.c */
link_phase_info_t *create_link_phase_info(void);
void destroy_link_phase_info(link_phase_info_t *lp);
void set_boundary_twist_fn(fn_links_t *fn_links, Real bdry_phase[4], int r0[4]);
void boundary_twist_fn(fn_links_t *fn_links, int flag);
void custom_rephase( su3_matrix **internal_links, int flag, int *status_now );

/* ff_opt.c */
void mult_adj_su3_fieldlink_lathwvec( su3_matrix *link,
				      half_wilson_vector **src_pt, 
				      half_wilson_vector *dest);
void mult_su3_sitelink_lathwvec( int dir, 
				 half_wilson_vector **src_pt, 
				 half_wilson_vector *dest);
void scalar_mult_add_lathwvec_proj(anti_hermitmat *mom, 
				   half_wilson_vector *back, 
				   half_wilson_vector *forw, Real coeff[2]);
void scalar_mult_add_lathwvec(half_wilson_vector *dest, 
			      half_wilson_vector *src, Real s[2]);
void mult_su3_fieldlink_lathwvec( su3_matrix *link,
				  half_wilson_vector **src_pt, 
				  half_wilson_vector *dest);

#endif /* _FERMION_LINKS_MILC_H */

