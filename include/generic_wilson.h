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

/* For various inversion routines.  Not used sytematically yet. CD */
enum guess_params { START_ZERO_GUESS = 0 ,  START_NONZERO_GUESS } ;  

#define ALL_T_SLICES -1

/* Structure defining Wilson (or clover) quark source */
/* There must be a color and spin member */
/* Add other members to suit the generic_wilson code 
   that builds the source.  Ignore the members you don't need. */
typedef struct {
  int type;           /* source type for most source builders */
  char descrp[30];    /* alpha description for most */
  int color;          /* source color */
  int spin;           /* source spin  */
  int wall_cutoff;    /* half size of box for w_source_h */
  int parity;         /* even or odd sites for w_source_h */
  Real r0;            /* source size for gaussian, width for gauge invt  */
  int iters;          /* iterations for gauge invariant source */
  int x0,y0,z0,t0;    /* source coordinates for most */ 
  char source_file[MAXFILENAME]; /* file name for some sources */
  complex *c_src_ptr; /* Pointer for complex field storage */
  wilson_vector *wv_src_ptr; /* Pointer for wilson vector field storage */
  int src_pointer ;   /* smearing function (for the moment, only
		         clover_finite_p_vary/create_wilson_source.c) */
} wilson_quark_source;

/* For quark source routines */
/* The Weyl representation types are included for w_source_h */
enum source_type { 
  POINT = 1, GAUSSIAN, CUTOFF_GAUSSIAN,
  POINT_WEYL, CUTOFF_GAUSSIAN_WEYL, COVARIANT_GAUSSIAN,
  COMPLEX_FIELD_FILE, DIRAC_FIELD_FILE, COMPLEX_FIELD_STORE,
  DIRAC_FIELD_STORE } ;

/* baryon_cont.c */

void baryon_cont1(wilson_prop_field src1, wilson_prop_field src2, 
		  wilson_prop_field src3, 
		  int chi_b[4][4], int eps[3][3][3], complex *prop);

void baryon_cont2(wilson_prop_field src1, wilson_prop_field src2, 
		  wilson_prop_field src3, 
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

/* gammas.c */
void mult_sw_by_gamma_l(spin_wilson_vector * src,spin_wilson_vector * dest, int dir);

void mult_sw_by_gamma_r(spin_wilson_vector * src,spin_wilson_vector * dest, int dir);
gamma_matrix_t gamma_mat(enum gammatype i);
int gamma_index(char *label);


/* gauss_smear_w.c */
void gauss_smear_field(wilson_vector *src, Real width, int iters, int t0);
void gauss_smear_site(field_offset src, Real width, int iters, int t0);

/* meson_cont.c */
void meson_cont_site(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop);
void meson_cont_field(spin_wilson_vector *src1, spin_wilson_vector *src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop);

/* w_source.c */
int ask_quark_source( FILE *fp, int prompt, int *type, char *descrp );
int get_quark_source( FILE *fp, int prompt, wilson_quark_source *wqs );
int get_quark_sink( FILE *fp, int prompt, wilson_quark_source *wqs );
void w_source_site(field_offset src,wilson_quark_source *wqs);
void w_source_field(wilson_vector *src,wilson_quark_source *wqs);
void w_sink_site(field_offset snk, wilson_quark_source *wqs);
void w_sink_field(wilson_vector *snk, wilson_quark_source *wqs);
void w_sink_scalar(field_offset snk,wilson_quark_source *wqs);

/* w_source_h.c (also has an ask_quark_source) */
Real *make_template(Real gamma, int cutoff);
void w_source_h(field_offset src,wilson_quark_source *wqs);
void free_source_template();

void bj_to_weyl( wilson_vector *src, wilson_vector *dest);
void dslash_w_site(field_offset src,field_offset dest,
	    int isign,int parity);
void dslash_w_site_special(field_offset src,field_offset dest,
		    int isign,int parity,msg_tag **tag,int is_started);
void dslash_w_field( wilson_vector *src, wilson_vector *dest, 
		     int isign, int parity);
void dslash_w_field_special(wilson_vector *src, wilson_vector *dest,
			    int isign,int parity,msg_tag **tag,int is_started);
void dslash_w_3D( field_offset src, field_offset dest, int isign, int parity);

void cleanup_dslash_temps();
void cleanup_tmp_links();

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
			wilson_quark_source *wqs),  /* source function */
    wilson_quark_source *wqs, /* source parameters */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp           /* Passthrough Dirac matrix parameters */
    );

int wilson_invert_field_wqs( /* Return value is number of iterations taken */
    wilson_quark_source *wqs, /* source parameters */
    void (*source_func_field)(wilson_vector *src, 
			      wilson_quark_source *wqs),  /* source function */
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
  complex ***prop,        /* where result is stored */
  spin_wilson_vector *src1, /* quark propagator */
  spin_wilson_vector *src2, /* quark propagator */
  int no_q_momenta,       /* number of q values */
  int q_momstore[][3],   /* q values themselves */
  int no_gamma_corr,     /* number of meson types */
  int meson_index[],     /* where in prop the propagator is accumulated */
  
  int gin[],              /* Gamma matrix type for source */
  int gout[],             /* Gamma matrix type for sink */
  complex meson_phase[]   /* phase factor to apply to correlator */
		    );
complex decode_phase(char *label);

/* w_baryon.c */
void w_baryon(wilson_prop_field src1,wilson_prop_field src2,
	      wilson_prop_field src3, complex *prop[4]);
void w_baryon_hl(wilson_prop_field src1,wilson_prop_field src2,
		 wilson_prop_field src3, complex *prop[6]);

#endif /* _GENERIC_WILSON_H */

