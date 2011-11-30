#ifndef _GENERIC_KS_QOP_H
#define _GENERIC_KS_QOP_H
/******************** generic_ks_qop.h *********************************
*  MIMD version 7 		  				       *
*/

#include <qop.h>
#include "../include/generic_ks.h"

/* d_congrad5_fn_qopqdp.c */

void initialize_congrad( void );
void finalize_congrad( void );

/* d_congrad5_fn_qop_F.c */

int ks_congrad_qop_F_site2site( quark_invert_control *qic,
				float *masses[], int nmass[], 
				field_offset milc_srcs[], 
				field_offset *milc_sols[], int nsrc,
				imp_ferm_links_t *fn );

int ks_congrad_qop_F_site2field( quark_invert_control *qic,
				 float *masses[], int nmass[], 
				 field_offset milc_srcs[], 
				 su3_vector **milc_sols[], int nsrc,
				 imp_ferm_links_t *fn );

int ks_congrad_qop_F_field2field( quark_invert_control *qic,
				  float *masses[], int nmass[], 
				  su3_vector *milc_srcs[], 
				  su3_vector **milc_sols[], int nsrc,
				  imp_ferm_links_t *fn );

int ks_congrad_milcfield2qop_F( su3_vector *milc_src, su3_vector *milc_sol, 
				quark_invert_control *qic, Real mass,
				imp_ferm_links_t *fn );

int ks_congrad_milc2qop_F( field_offset milc_src, field_offset milc_sol, 
			   quark_invert_control *qic, Real mass,
			   imp_ferm_links_t *fn );

/* d_congrad5_fn_qop_D.c */

int ks_congrad_qop_D_site2site( quark_invert_control *qic,
				double *masses[], int nmass[], 
				field_offset milc_srcs[], 
				field_offset *milc_sols[], int nsrc,
				imp_ferm_links_t *fn );

int ks_congrad_qop_D_site2field( quark_invert_control *qic,
			      double *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[], int nsrc,
				 imp_ferm_links_t *fn );

int ks_congrad_qop_D_field2field( quark_invert_control *qic,
				  double *masses[], int nmass[], 
				  su3_vector *milc_srcs[], 
				  su3_vector **milc_sols[], int nsrc,
				  imp_ferm_links_t *fn );

int ks_congrad_milcfield2qop_D( su3_vector *milc_src, su3_vector *milc_sol, 
				quark_invert_control *qic, Real mass,
				imp_ferm_links_t *fn );

int ks_congrad_milc2qop_D( field_offset milc_src, field_offset milc_sol, 
			   quark_invert_control *qic, Real mass,
			   imp_ferm_links_t *fn );

/* dslash_fn_qop_milc_F.c */

void cleanup_gathers_qop_milc_F(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_qop_milc_temps_F(void);
void dslash_fn_qop_milc_F( su3_matrix *fatlinks, su3_matrix *longlinks,
			 su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_qop_milc_field_special_F(su3_matrix *fatlinks, 
				      su3_matrix *longlinks,
				      su3_vector *src, su3_vector *dest,
				      int parity, msg_tag **tag, int start );

/* dslash_fn_qop_milc_F.c */

void cleanup_gathers_qop_milc_D(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_qop_milc_temps_D(void);
void dslash_fn_qop_milc_D( su3_matrix *fatlinks, su3_matrix *longlinks,
			 su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_qop_milc_field_special_D(su3_matrix *fatlinks, 
				      su3_matrix *longlinks,
				      su3_vector *src, su3_vector *dest,
				      int parity, msg_tag **tag, int start );

/* ks_multicg_offset_qop.c */

int ks_multicg_offset_field_F(       /* Return value is number of iterations taken */
    su3_vector *src,	       /* source vector (type su3_vector) */
    su3_vector **psim,	       /* solution vectors */
    ks_param *ksp,	       /* the offsets */
    int num_offsets,	       /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    imp_ferm_links_t *fn             /* Storage for fat and Naik links */
			);

int ks_multicg_offset_field_D(       /* Return value is number of iterations taken */
    su3_vector *src,	       /* source vector (type su3_vector) */
    su3_vector **psim,	       /* solution vectors */
    ks_param *ksp,	       /* the offsets */
    int num_offsets,	       /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    imp_ferm_links_t *fn             /* Storage for fat and Naik links */
			);

int ks_multicg_mass_field_F(	      /* Return value is number of iterations taken */
    su3_vector *src,	      /* source vector (type su3_vector) */
    su3_vector **psim,	      /* solution vectors (preallocated) */
    Real *masses,	      /* the masses */
    int num_masses,	      /* number of masses */
    quark_invert_control *qic,/* inversion parameters */
    imp_ferm_links_t *fn            /* Storage for fat and Naik links */
			);

int ks_multicg_mass_field_D(	      /* Return value is number of iterations taken */
    su3_vector *src,	      /* source vector (type su3_vector) */
    su3_vector **psim,	      /* solution vectors (preallocated) */
    Real *masses,	      /* the masses */
    int num_masses,	      /* number of masses */
    quark_invert_control *qic,/* inversion parameters */
    imp_ferm_links_t *fn            /* Storage for fat and Naik links */
			);

#endif /* _GENERIC_KS_QOP_H */
