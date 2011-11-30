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
			   Real eps[], int nsrc,
			   ferm_links_t *fn );

int ks_congrad_milcfield2qdp_F(su3_vector *f_src, su3_vector *f_dest, 
			       quark_invert_control *qic, Real mass,
			       ferm_links_t *fn );

int ks_congrad_milcfield2qdp_D(su3_vector *f_src, su3_vector *f_dest, 
			       quark_invert_control *qic, Real mass,
			       ferm_links_t *fn );

int ks_congrad_milc2qdp_F(field_offset f_src, field_offset f_dest, 
			  quark_invert_control *qic, Real mass,
			  ferm_links_t *fn );

int ks_congrad_milc2qdp_D(field_offset f_src, field_offset f_dest, 
			  quark_invert_control *qic, Real mass,
			  ferm_links_t *fn );

int ks_congrad_qdp_F(QDP_ColorVector *src, QDP_ColorVector *dest, 
		     quark_invert_control *qic, QLA_Real mass,
		     ferm_links_t *fn );

int ks_congrad_qdp_D(QDP_ColorVector *src, QDP_ColorVector *dest, 
		     quark_invert_control *qic, QLA_Real mass,
		     ferm_links_t *fn );

/* fermion_force_asqtad_qdp_D.c */

void fermion_force_asqtad_qdp_D( QDP_D3_ColorMatrix *force[], 
				 QDP_D3_ColorMatrix *gf[], 
				 asqtad_path_coeff *coeffs, 
				 QDP_D3_ColorVector *in_pt[], 
				 double eps[], int nsrc,
				 ferm_links_t *fn );
void eo_fermion_force_oneterm_D( Real eps, Real weight, su3_vector *x_off,
				 ferm_links_t *fn );
void eo_fermion_force_twoterms_D( Real eps, Real weight1, Real weight2, 
				  su3_vector *x1_off, su3_vector *x2_off,
				  ferm_links_t *fn );
void fermion_force_asqtad_multi_D( Real eps, Real *residues,
				   su3_vector *xxx[], int nsrc,
				   ferm_links_t *fn);
void fermion_force_asqtad_block_D( Real eps, Real *residues, 
				   su3_vector **xxx, int nterms, 
				   int veclength,
				   ferm_links_t *fn );


/* fermion_force_asqtad_qdp_F.c */

void fermion_force_asqtad_qdp_F( QDP_F3_ColorMatrix *force[], 
				 QDP_F3_ColorMatrix *gf[], 
				 asqtad_path_coeff *coeffs, 
				 QDP_F3_ColorVector *in_pt[], 
				 float eps[], int nsrc,
				  ferm_links_t *fn );
void eo_fermion_force_oneterm_F( Real eps, Real weight, su3_vector *x_off,
				  ferm_links_t *fn );
void eo_fermion_force_twoterms_F( Real eps, Real weight1, Real weight2, 
				  su3_vector *x1_off, su3_vector *x2_off,
				  ferm_links_t *fn );
void fermion_force_asqtad_multi_F( Real eps, Real *residues,
				   su3_vector *xxx[], int nsrc,
				  ferm_links_t *fn);
void fermion_force_asqtad_block_F( Real eps, Real *residues, 
				   su3_vector **xxx, int nterms, int veclength,
				   ferm_links_t *fn );


/* ks_multicg_offset_qdp.c */

int ks_multicg_offset_field_F( /* Return value is number of iterations taken */
    su3_vector *src,	 /* source vector (type su3_vector) */
    su3_vector **psim,	 /* solution vectors */
    Real *offsets,	 /* the offsets */
    int num_offsets,	 /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    ferm_links_t *fn       /* Storage for fat and Naik links */
			);

int ks_multicg_offset_field_D( /* Return value is number of iterations taken */
    su3_vector *src,	 /* source vector (type su3_vector) */
    su3_vector **psim,	 /* solution vectors */
    Real *offsets,	 /* the offsets */
    int num_offsets,	 /* number of offsets */
    quark_invert_control *qic,  /* inversion parameters */
    ferm_links_t *fn       /* Storage for fat and Naik links */
			 );

int ks_multicg_mass_field_F(	 /* Return value is number of iterations taken */
    su3_vector *src,	 /* source vector (type su3_vector) */
    su3_vector **psim,	 /* solution vectors (preallocated) */
    Real *masses,	 /* the masses */
    int num_masses,	 /* number of masses */
    quark_invert_control *qic,  /* inversion parameters */
    ferm_links_t *fn       /* Storage for fat and Naik links */
			);

int ks_multicg_mass_field_D(	 /* Return value is number of iterations taken */
    su3_vector *src,	 /* source vector (type su3_vector) */
    su3_vector **psim,	 /* solution vectors (preallocated) */
    Real *masses,	 /* the masses */
    int num_masses,	 /* number of masses */
    quark_invert_control *qic,  /* inversion parameters */
    ferm_links_t *fn       /* Storage for fat and Naik links */
			);



#endif /* _GENERIC_KS_QDP_H */
