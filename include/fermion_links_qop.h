#ifndef _FERMION_LINKS_QOP_H
#define _FERMION_LINKS_QOP_H
/************************ fermion_links_qop.h ***************************
*									*
*  Definitions for fermion links
*  MIMD version 7 							*
*									*
*/

#include <qop.h>
#include "../include/su3.h"
#include "../include/imp_ferm_links.h"
#include "../include/ks_action_coeffs_qop.h"
#include "../include/hisq_links_qop.h"

/*
  
  Structure nesting for asqtad actions

  qop_asqtad_links_t
    qop_fm_ac_links_t *fm_ac
       QOP_asqtad_coeffs_t *ac
       fn_links_qop_t *fm
         link_phase_info_t *phase
         QOP_F3_FermionLinksAsqtad al_F;
         QOP_D3_FermionLinksAsqtad al_D;
    qop_fm_ac_links_t *fm_ac_du0
       QOP_asqtad_coeffs_t *ac;
       fn_links_qop_t *fm
         link_phase_info_t *phase
	 int al_F_allocated;
         QOP_F3_FermionLinksAsqtad al_F;
         QOP_D3_FermionLinksAsqtad al_D;

  Structure nesting for HISQ actions

   qop_hisq_links_t
     qop_hisq_ac_links_t *hisq_ac
       QOP_hisq_coeffs_t *ac
       hisq_links_qop_t *hl
         link_phase_info_t *phase;
         QOP_F3_FermionLinksHisq hl_F;
         QOP_D3_FermionLinksHisq hl_D;
       fn_links_qop_t *fn[MAX_NAIK]
         link_phase_info_t *phase;
         QOP_F3_FermionLinksAsqtad al_F;
         QOP_D3_FermionLinksAsqtad al_D;
       fn_links_qop_t *fn_deps

   where (in qopqdp-0.17.0)

   QOP_hisq_coeffs_t
     int n_naiks;
     double eps_naik[MAX_NAIK];
     int umethod;
     int ugroup;
     double fat7_one_link;
     double fat7_three_staple;
     double fat7_five_staple;
     double fat7_seven_staple;
     double asqtad_one_link;
     double asqtad_three_staple;
     double asqtad_five_staple;
     double asqtad_seven_staple;
     double asqtad_lepage;
     double asqtad_naik;
     double difference_one_link;
     double difference_naik;

   QOP_FermionLinksHisq
     int n_naiks, WeqY;
     QDP_ColorMatrix *U_link[4];
     QDP_ColorMatrix *V_link[4];
     QDP_ColorMatrix *Y_unitlink[4];
     QDP_ColorMatrix *W_unitlink[4];
     QOP_FermionLinksAsqtad **fn;
     QOP_FermionLinksAsqtad fn_deps;

*/

/*********************************************************************/
/* QOP fermion links structures                                      */
/*********************************************************************/

/* asqtad */

typedef struct {
  QOP_asqtad_coeffs_t *ac;
  fn_links_qop_t *fm;
} qop_fm_ac_links_t; 

typedef struct {
  qop_fm_ac_links_t *fm_ac;
  qop_fm_ac_links_t *fm_ac_du0;
} qop_asqtad_links_t;

/*********************************************************************/
/* hisq */

typedef struct {
  QOP_hisq_coeffs_t *ac;
  hisq_links_qop_t *hl;
  fn_links_qop_t *fn[MAX_NAIK];
  fn_links_qop_t *fn_deps;
} qop_hisq_ac_links_t;

typedef struct {
  qop_hisq_ac_links_t *hisq_ac;
} qop_hisq_links_t;


/* PLACEHOLDER FOR NOW */
typedef struct {
  int phases_in;         // track KS phases in the V and Y links
  int WeqY;              // true if W = Y
  su3_matrix *U_link[4]; // original gauge matrices, stored as four fields
  su3_matrix *V_link[4]; // first iteration of fattening
  su3_matrix *Y_unitlink[4]; // unitary projection of V_link, U(3)
  su3_matrix *W_unitlink[4]; // special unitary projection of Y_link, SU(3)
} hisq_auxiliary_t;


#include "../include/macros.h"

/* fermion_links_generic_qop_F.c */

//void create_qop_F_links_from_milc_fn(fn_links_qop_t *ql, fn_links_t *fl);
//void create_qop_F_fermion_links( fn_links_qop_t *ql, fn_links_t *fl );
//void invalidate_qop_ferm_links_F( fn_links_qop_t *fn );

/* fermion_links_asqtad_qop_F.c */
/* fermion_links_hisq_qop_F.c */

//void load_ferm_links_F(fermion_links_t *fl);
//void load_ferm_links_dmdu0_F(fermion_links_t *fl);
//void look_at_link_F(fermion_links_t *fl, int *x, int dir);
//void invalidate_all_ferm_links_F( fermion_links_t *fl );

/* fermion_links_generic_qop_D.c */

//void create_qop_D_fermion_links( fn_links_qop_t *ql, fn_links_t *fl );
//void invalidate_qop_ferm_links_D( fn_links_qop_t *ql );

/* fermion_links_asqtad_qop_D.c */
/* fermion_links_hisq_qop_D.c */

//void load_ferm_links_D(fermion_links_t *fl);
//void load_ferm_links_dmdu0_D(fermion_links_t *fl);
//void look_at_link_D(fermion_links_t *fl, int *x, int dir);
//void invalidate_all_ferm_links_D( fermion_links_t *fl );

/* fermion_links_fn_twist_qop.c */

link_phase_info_t *create_link_phase_info(void);
void destroy_link_phase_info(link_phase_info_t *lp);
void set_boundary_twist_fn(fn_links_qop_t *fn_links, Real bdry_phase[4], int r0[4]);
void boundary_twist_fn(fn_links_qop_t *fn_links, int flag);

/********************************************************************/
/* Action coefficient conversion */
/********************************************************************/

#include "../include/ks_action_paths.h"

/* load_qop_asqtad_coeffs_F.c */

void load_qop_F_asqtad_coeffs(QOP_asqtad_coeffs_t *c, Real weight,
			      ks_action_paths *ap);

/* load_qop_asqtad_coeffs_D.c */

void load_qop_D_asqtad_coeffs(QOP_asqtad_coeffs_t *c, Real weight,
			      ks_action_paths *ap);

/* load_qop_hisq_coeffs_F.c */

void load_qop_F_hisq_coeffs(QOP_hisq_coeffs_t *c, Real weight,
			    ks_action_paths_hisq *ap);

/* load_qop_hisq_coeffs_D.c */

void load_qop_D_hisq_coeffs(QOP_hisq_coeffs_t *c, Real weight,
			    ks_action_paths_hisq *ap);


#endif /* _FERMION_LINKS_QOP_H */
