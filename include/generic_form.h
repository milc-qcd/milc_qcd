#ifndef _GENERIC_FORM_H
#define _GENERIC_FORM_H
/************************ generic_form.h *************************************
*									*
*  Macros and declarations for miscellaneous generic routines           *
*  This header is for codes that call generic_form routines             *
*  MIMD version 7 							*
*									*
*/

#include "../include/macros.h"
#include "../include/su3.h"
#include "../generic_form/gammatypes.h"

void c_scale_wilson_vector2(wilson_vector *m , complex *scale);
void copy_site_wilson_vector(field_offset src, field_offset dest) ;
void flip_source_re(field_offset quark_prop);
int load_momentum_from_disk(int mom_in[][3], char filename[], int max_mom);
int load_momentum(int prompt, char *label, int *no_mom, int mom_in[][3], int max_mom);
void load_smearing(field_offset where_smear, char filename[]);
void mult_gamma(int phase, gamma_matrix *g1, gamma_matrix *g2, gamma_matrix *g3);
void make_gammas(gamma_matrix *gamma);
void mult_sw_by_gamma_l(spin_wilson_vector * src,
			spin_wilson_vector * dest, int dir);
void mult_sw_by_gamma_r(spin_wilson_vector * src,
			spin_wilson_vector * dest, int dir);
void meson_cont_mom(complex prop[],
		    field_offset src1,field_offset src2,
		    int base_pt, int q_stride, int op_stride,
		    gamma_corr gamma_table[], int no_gamma_corr);
void meson_cont_mom_lean(complex prop[],
		    field_offset src1,field_offset src2,
		    int base_pt, int q_stride, int op_stride,
		    int w_meson_store_t[],
		    gamma_corr gamma_table[], int no_gamma_corr);

void meson_cont_mom_lean2(
  complex prop[],        /* where result is stored */
  field_offset src1,     /* quark propagator (spin_wilson_vector) */
  field_offset src2,     /* quark propagator (spin_wilson_vector) */
  int base_pt,           /* used for indexing prop (see prop_pt below) */
  int q_stride,          /* used for indexing prop (see prop_pt below) */
  int op_stride,         /* used for indexing prop (see prop_pt below) */
  int w_meson_store_t[], /* for compact storage */
  int w_meson_my_t[],    /* for compact storage */
  int w_meson_nstore,    /* for compact storage */
  int no_q_momenta,       /* number of q values */
  int q_momstore[][3],   /* q values themselves */
  int no_gamma_corr,     /* number of pairs */
  gamma_corr gamma_table[], /* table of gamma matrix pairs */
  field_offset tmp,      /* vector of complex numbers */
  int dimtmp             /* maximum number of complex values in tmp */
  );

void load_wilson_source(field_offset src, field_offset dest,int color,int spin);

void load_wvec(wilson_vector *dest, complex *z, int spin, int colour) ;




#endif /* _GENERIC_FORM_H */
