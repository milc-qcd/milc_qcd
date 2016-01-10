#ifndef _GENERIC_KS_QPHIX_H
#define _GENERIC_KS_QPHIX_H
/******************** generic_ks_qphix.h *********************************
*  MIMD version 7 		  				       *
*/

#include <qphix.h>
#include "../include/generic_ks.h"


/* d_congrad5_fn_qphix_F.c */

fsu3_matrix *
create_backlinks_without_adjoint_F(su3_matrix *t);

void destroy_backlinks_F(fsu3_matrix *);

QPHIX_F3_FermionLinksAsqtad *
create_qphix_F_L_from_field(fn_links_t *fn, int parity);

int ks_congrad_parity_qphix_F( su3_vector *src,
			       su3_vector *sol,
			       quark_invert_control *qic,
			       Real mass,
			       fn_links_t *fn);
/* d_congrad5_fn_qphix_D.c */


dsu3_matrix *
create_backlinks_without_adjoint_D(su3_matrix *t);

void destroy_backlinks_D(dsu3_matrix *);

QPHIX_D3_FermionLinksAsqtad *
create_qphix_D_L_from_field(fn_links_t *fn, int parity);

int ks_congrad_parity_qphix_D( su3_vector *src,
			       su3_vector *sol,
			       quark_invert_control *qic,
			       Real mass,
			       fn_links_t *fn);

/* ks_multicg_offset_qphix_F.c */

int ks_multicg_offset_field_qphix_F(       /* Return value is number of iterations taken */
    su3_vector *src,	       /* source vector (type su3_vector) */
    su3_vector **psim,	       /* solution vectors */
    ks_param *ksp,	       /* the offsets */
    int num_offsets,	       /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    imp_ferm_links_t *fn             /* Storage for fat and Naik links */
			);

/* ks_multicg_offset_qphix_D.c */

int ks_multicg_offset_field_qphix_D(       /* Return value is number of iterations taken */
    su3_vector *src,	       /* source vector (type su3_vector) */
    su3_vector **psim,	       /* solution vectors */
    ks_param *ksp,	       /* the offsets */
    int num_offsets,	       /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    imp_ferm_links_t *fn             /* Storage for fat and Naik links */
			);

#endif /* _GENERIC_KS_QPHIX_H */
