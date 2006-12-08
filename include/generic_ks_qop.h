#ifndef _GENERIC_KS_QOP_H
#define _GENERIC_KS_QOP_H
/******************** generic_ks_qop.h *********************************
*  MIMD version 7 		  				       *
*/

#include <qop.h>

/* d_congrad5_fn_qopqdp.c */

void initialize_congrad( void );
void finalize_congrad( void );

/* d_congrad5_fn_qop_F.c */

int ks_congrad_qop_F_site2site(int niter, int nrestart, Real rsqmin, 
			     float *masses[], int nmass[], 
			     field_offset milc_srcs[], 
			     field_offset *milc_sols[],
			     int nsrc, Real* final_rsq_ptr, int milc_parity );

int ks_congrad_qop_F_site2field(int niter, int nrestart, Real rsqmin, 
			      float *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[],
			      int nsrc, Real* final_rsq_ptr, int milc_parity );

int 
ks_congrad_milc2qop_F( field_offset milc_src, field_offset milc_sol, Real mass,
		       int niter, int nrestart, Real rsqmin, 
		       int milc_parity, Real* final_rsq_ptr );

/* d_congrad5_fn_qop_D.c */

int ks_congrad_qop_D_site2site(int niter, int nrestart, Real rsqmin, 
			     double *masses[], int nmass[], 
			     field_offset milc_srcs[], 
			     field_offset *milc_sols[],
			     int nsrc, Real* final_rsq_ptr, int milc_parity );

int ks_congrad_qop_D_site2field(int niter, int nrestart, Real rsqmin, 
			      double *masses[], int nmass[], 
			      field_offset milc_srcs[], 
			      su3_vector **milc_sols[],
			      int nsrc, Real* final_rsq_ptr, int milc_parity );

int 
ks_congrad_milc2qop_D( field_offset milc_src, field_offset milc_sol, Real mass,
		       int niter, int nrestart, Real rsqmin,
		       int milc_parity, Real* final_rsq_ptr );

/* dslash_fn_qop_milc_F.c */

void cleanup_gathers_qop_milc_F(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_qop_milc_temps_F();
void dslash_fn_qop_milc_F( su3_matrix *fatlinks, su3_matrix *longlinks,
			 su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_qop_milc_field_special_F(su3_matrix *fatlinks, 
				      su3_matrix *longlinks,
				      su3_vector *src, su3_vector *dest,
				      int parity, msg_tag **tag, int start );

/* dslash_fn_qop_milc_F.c */

void cleanup_gathers_qop_milc_D(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_qop_milc_temps_D();
void dslash_fn_qop_milc_D( su3_matrix *fatlinks, su3_matrix *longlinks,
			 su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_qop_milc_field_special_D(su3_matrix *fatlinks, 
				      su3_matrix *longlinks,
				      su3_vector *src, su3_vector *dest,
				      int parity, msg_tag **tag, int start );

/* fermion_links_asqtad_qop_F.c */

QOP_FermionLinksAsqtad *create_qop_F_asqtad_fermion_links( void );
void destroy_qop_F_asqtad_fermion_links( void );
void load_fn_links_F(void);
void invalidate_fn_links_F(void);

/* fermion_links_asqtad_qop_D.c */

QOP_FermionLinksAsqtad *create_qop_D_asqtad_fermion_links( void );
void destroy_qop_D_asqtad_fermion_links( void );
void load_fn_links_D(void);
void invalidate_fn_links_D(void);

/* ks_multicg_offset_qop.c */

int ks_multicg_offset_F(/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
			);

int ks_multicg_offset_D(/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
			);

int ks_multicg_mass_F(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors (preallocated) */
    float *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
			);

int ks_multicg_mass_D(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors (preallocated) */
    double *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr	/* final residue squared */
			);

/* load_qop_asqtad_coeffs_F.c */

void load_qop_F_asqtad_coeffs(QOP_asqtad_coeffs_t *c, Real weight,
			    Real *act_path_coeff);

/* load_qop_asqtad_coeffs_D.c */

void load_qop_D_asqtad_coeffs(QOP_asqtad_coeffs_t *c, Real weight,
			    Real *act_path_coeff);

#endif /* _GENERIC_KS_QOP_H */
