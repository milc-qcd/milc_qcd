#ifndef _GENERIC_KS_QDP_H
#define _GENERIC_KS_QDP_H
/******************** generic_ks_qdp.h *********************************
*  MIMD version 7 							*
*/

#include <qdp.h>

typedef struct {
  QLA_Real one_link     ; 
  QLA_Real naik         ;
  QLA_Real three_staple ;
  QLA_Real five_staple  ;
  QLA_Real seven_staple ;
  QLA_Real lepage       ;
} asqtad_path_coeff;

/* dslash_fn_qdp_F.c */

void setup_dslash_F(void);

void dslash_qdp_F_fn_special2(QDP_ColorVector *src, QDP_ColorVector *dest,
			      QDP_Subset parity, QDP_ColorVector *temp[]);
void dslash_qdp_fn_F(QDP_ColorVector *src, QDP_ColorVector *dest, 
		     QDP_Subset parity);

/* dslash_fn_qdp_D.c */

void setup_dslash_D(void);

void dslash_qdp_D_fn_special2(QDP_ColorVector *src, QDP_ColorVector *dest,
			      QDP_Subset parity, QDP_ColorVector *temp[]);
void dslash_qdp_fn_D(QDP_ColorVector *src, QDP_ColorVector *dest, 
		     QDP_Subset parity);

void fn_fermion_force_qdp( QDP_ColorMatrix *force[], QDP_ColorMatrix *gf[], 
		      asqtad_path_coeff *coeffs, QDP_ColorVector *in_pt[], 
		      Real eps[], int nsrc );

int ks_congrad_milc2qdp_F(field_offset f_src, field_offset f_dest, Real mass,
			  int niter, int nrestart, Real rsqmin, int parity, 
			  Real *final_rsq_ptr);

int ks_congrad_milc2qdp_D(field_offset f_src, field_offset f_dest, Real mass,
			  int niter, int nrestart, Real rsqmin, int parity, 
			  Real *final_rsq_ptr);

int ks_congrad_qdp_F(QDP_ColorVector *src, QDP_ColorVector *dest, 
		     QLA_Real mass, int niter, int nrestart, QLA_Real rsqmin, 
		     QDP_Subset parity, QLA_Real *final_rsq_ptr);


/* ks_multicg_offset_qdp.c */

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



#endif /* _GENERIC_KS_QDP_H */
