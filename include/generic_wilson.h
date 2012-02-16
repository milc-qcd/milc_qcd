#ifndef _GENERIC_WILSON_H
#define _GENERIC_WILSON_H
/************************ generic_wilson.h ******************************
*									*
*  Macros and declarations for generic_wilson routines                  *
*  This header is for codes that call generic_wilson routines           *
*  MIMD version 7 							*
*									*
*/

#include "../include/comdefs.h"
#include "../include/gammatypes.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"
#include "../include/su3.h"
#ifdef HAVE_QIO
#include <qio.h>
#endif

/* For various inversion routines.  Not used sytematically yet. CD */
enum guess_params { START_ZERO_GUESS = 0 ,  START_NONZERO_GUESS } ;  

/* baryon_cont.c */

void baryon_cont1(wilson_prop_field * src1, wilson_prop_field * src2, 
		  wilson_prop_field * src3, 
		  int chi_b[4][4], int eps[3][3][3], complex *prop);

void baryon_cont2(wilson_prop_field * src1, wilson_prop_field * src2, 
		  wilson_prop_field * src3, 
		  int chi_b[4][4], int eps[3][3][3], complex *prop);

/* boundary_flip.c */
void boundary_flip(int sign );

/* canopy2weyl_rot.c */
void weyl2canopy_site(field_offset src, field_offset dest);
void canopy2weyl_site(field_offset src, field_offset dest);
void weyl2canopy_field(wilson_propagator *src, wilson_propagator *dest);
void canopy2weyl_field(wilson_propagator *src, wilson_propagator *dest);
void convert_wprop_fnal_to_milc_site(field_offset wprop);
void convert_wprop_fnal_to_milc_field(wilson_propagator *wprop);
void convert_wprop_milc_to_fnal_site(field_offset wprop);
void convert_wprop_milc_to_fnal_field(wilson_propagator *wprop);

/* gauss_smear_w.c */
void gauss_smear_wv_field(wilson_vector *src, su3_matrix *t_links,
			  Real width, int iters, int t0);
void gauss_smear_wv_site(field_offset src, su3_matrix *t_links,
			 Real width, int iters, int t0);

/* meson_cont.c */
void meson_cont_site(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop);
void meson_cont_field(spin_wilson_vector *src1, spin_wilson_vector *src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop);

/* staggered2naive.c */

void convert_ksprop_to_wprop_swv(spin_wilson_vector *swv, 
				 su3_vector *ksp, int r[], int r0[]);
void convert_naive_to_staggered_wv(wilson_vector *wv, int r[], int r0[]);
void convert_staggered_to_naive_wv(wilson_vector *wv, int r[], int r0[]);
void check_naive(wilson_vector *dst, wilson_vector *src, Real mass, Real tol);
void dslash_naive(wilson_vector *src, wilson_vector *dst);

/* w_source_h.c (also has an ask_quark_source) */
Real *make_template(Real gamma, int cutoff);
void w_source_h(field_offset src,quark_source *wqs);
void free_source_template(void);

void bj_to_weyl( wilson_vector *src, wilson_vector *dest);
void dslash_w_site(field_offset src,field_offset dest,
	    int isign,int parity);
void dslash_w_site_special(field_offset src,field_offset dest,
		    int isign,int parity,msg_tag **tag,int is_started);
void dslash_w_field( wilson_vector *src, wilson_vector *dest, 
		     int isign, int parity);
void dslash_w_field_special(wilson_vector *src, wilson_vector *dest,
			    int isign,int parity,msg_tag **tag,int is_started);
void dslash_w_3D_site( field_offset src, field_offset dest, 
		       int isign, int parity);
void dslash_w_3D_field( wilson_vector *src, wilson_vector *dest, 
			int isign, int parity);
void hop_w_field( wilson_vector *src, wilson_vector *dest, 
		  int isign, int iphase, int parity, int dir);

void cleanup_dslash_w_3D_temps(void);
void cleanup_dslash_wtemps(void);
void cleanup_tmp_links(void);

/* wilson_invert.c */

int wilson_invert_site(    /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (source already created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp            /* Passthrough Dirac matrix parameters */
    );

int wilson_invert_field( /* Return value is number of iterations taken */
    wilson_vector *src, /* type wilson_vector (where source has been created)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess) */
    int (*invert_func_field)(wilson_vector *src, wilson_vector *dest,
			     quark_invert_control *qic,
			     void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp                 /* Passthrough Dirac matrix parameters */
    );

int wilson_invert_site_wqs( /* Return value is number of iterations taken */
    field_offset src, 	/* where source is to be created */
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    void (*source_func)(field_offset src, 
			quark_source *wqs),  /* source function */
    quark_source *wqs, /* source parameters */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp           /* Passthrough Dirac matrix parameters */
    );

int wilson_invert_field_wqs( /* Return value is number of iterations taken */
    quark_source *wqs, /* source parameters */
    int (*source_func_field)(wilson_vector *src, 
			     quark_source *wqs),  /* source function */
    wilson_vector *dest,  /* type wilson_vector (answer and initial guess) */
    int (*invert_func_field)(wilson_vector *src, wilson_vector *dest,
			     quark_invert_control *qic, void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp           /* Passthrough Dirac matrix parameters */
    );

/* w_meson.c */
void w_meson_site(field_offset src1,field_offset src2,complex *prop[10]);
void w_meson_field(spin_wilson_vector *src1, spin_wilson_vector *src2,
		   complex *prop[10]);

/* w_meson_mom.c */
void meson_cont_mom(
  complex **prop,           /* prop[m][t] is where result is accumulated */
  spin_wilson_vector *src1, /* quark propagator */
  spin_wilson_vector *src2, /* quark propagator */
  int no_q_momenta,         /* no of unique mom/parity values (gt p) */
  int **q_momstore,         /* q_momstore[p] are the momentum components */
  char **q_parity,          /* q_parity[p] the parity of each mom component */
  int no_gamma_corr,        /* # of gamma src/snk combinations (gt g) */
  int num_corr_mom[],       /* number of momentum/parity values for each corr */
  int **corr_table,         /* c = corr_table[g][k] correlator index */
  int p_index[],            /* p = p_index[c] is the momentum index */
  int gout[],               /* gout[c] is the sink gamma */
  int gin[],                /* gin[c] is the source gamma */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_fact[],        /* meson_fact[c] scales the correlator */
  int corr_index[],         /* m = corr_index[c] is the correlator index */
  int r0[]                  /* spatial origin for defining FT phases */
		    );


/* w_meson_open_mom.c */

void meson_open_mom(
  wilson_propagator **prop, /* prop[m][t] is where result is accumulated */
  spin_wilson_vector *src1, /* quark propagator (to become antiquark) */
  spin_wilson_vector *src2, /* quark propagator */
  int no_q_momenta,         /* no of unique mom/parity values (gt p) */
  int **q_momstore,         /* q_momstore[p] are the momentum components */
  char **q_parity,          /* q_parity[p] the parity of each mom component */
  int no_gamma_corr,        /* # of gamma src/snk combinations (gt g) */
  int num_corr_mom[],       /* number of momentum/parity values for each corr */
  int **corr_table,         /* c = corr_table[g][k] correlator index */
  int p_index[],            /* p = p_index[c] is the momentum index */
  int gout[],               /* gout[c] is the sink gamma */
  int gin[],                /* gin[c] is the source gamma */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_factor[],      /* meson_factor[c] scales the correlator */
  int corr_index[],         /* m = corr_index[c] is the correlator index */
  int r0[]                  /* spatial origin for defining FT phases */
		    );


/* w_baryon.c */
void w_baryon(wilson_prop_field * src1,wilson_prop_field * src2,
	      wilson_prop_field * src3, complex *prop[4]);
void w_baryon_hl(wilson_prop_field * src1,wilson_prop_field * src2,
		 wilson_prop_field * src3, complex *prop[6]);

#endif /* _GENERIC_WILSON_H */

