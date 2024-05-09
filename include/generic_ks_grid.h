#ifndef _GENERIC_KS_GRID_H
#define _GENERIC_KS_GRID_H
/******************** generic_ks_grid.h *********************************
*  MIMD version 7 		  				       *
*/

#include "../include/mGrid/mGrid.h"
#include "../include/generic_ks.h"


/* d_congrad5_fn_grid_F.c */

int ks_congrad_parity_grid_F( su3_vector *src,
			      su3_vector *sol,
			      quark_invert_control *qic,
			      Real mass,
			      fn_links_t *fn);

int ks_congrad_block_parity_grid_F( int nrhs,
				    su3_vector *src[],
				    su3_vector *sol[],
				    quark_invert_control *qic,
				    Real mass,
				    fn_links_t *fn);
/* d_congrad5_fn_grid_D.c */

int ks_congrad_parity_grid_D( su3_vector *src,
			      su3_vector *sol,
			      quark_invert_control *qic,
			      Real mass,
			      fn_links_t *fn);

int ks_congrad_mixed_parity_grid_D( su3_vector *src,
				    su3_vector *sol,
				    quark_invert_control *qic,
				    Real mass,
				    fn_links_t *fn);

int ks_congrad_block_parity_grid_D( int nrhs,
				    su3_vector *src[],
				    su3_vector *sol[],
				    quark_invert_control *qic,
				    Real mass,
				    fn_links_t *fn);

int ks_congrad_mixed_block_parity_grid_D( int nrhs,
					  su3_vector *src[],
					  su3_vector *sol[],
					  quark_invert_control *qic,
					  Real mass,
					  fn_links_t *fn);

/* ks_multicg_offset_grid_F.c */

int ks_multicg_offset_field_grid_F(       /* Return value is number of iterations taken */
    su3_vector *src,	       /* source vector (type su3_vector) */
    su3_vector **psim,	       /* solution vectors */
    ks_param *ksp,	       /* the offsets */
    int num_offsets,	       /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    imp_ferm_links_t *fn             /* Storage for fat and Naik links */
			);

/* ks_multicg_offset_grid_D.c */

int ks_multicg_offset_field_grid_D(       /* Return value is number of iterations taken */
    su3_vector *src,	       /* source vector (type su3_vector) */
    su3_vector **psim,	       /* solution vectors */
    ks_param *ksp,	       /* the offsets */
    int num_offsets,	       /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    imp_ferm_links_t *fn             /* Storage for fat and Naik links */
			);

#endif /* _GENERIC_KS_GRID_H */
