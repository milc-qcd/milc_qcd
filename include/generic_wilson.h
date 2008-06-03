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

#define ALL_T_SLICES -1
#define MAXDESCRP 128
#define MAXSRCLABEL 8

/* Structure defining Wilson (or clover) quark source */
/* There must be a color and spin member */
/* Add other members to suit the generic_wilson code 
   that builds the source.  Ignore the members you don't need. */
typedef struct {
  int type;           /* source type for most source builders */
  char descrp[MAXDESCRP];  /* alpha description for most */
  char label[MAXSRCLABEL]; /* Abbreviation of description */
  int color;          /* source color */
  int spin;           /* source spin  */
  int wall_cutoff;    /* half size of box for w_source_h */
  int parity;         /* even or odd sites for w_source_h */
  Real r0;            /* source size for gaussian, width for gauge invt  */
  Real d1;            /* Fermilab 3D rotation parameter */
  int iters;          /* iterations for gauge invariant source */
  int x0,y0,z0,t0;    /* source coordinates for most */ 
  char source_file[MAXFILENAME]; /* file name for some sources */
  int flag;           /* mode of reading or writing for some sources */
  int mom[3];         /* insertion momentum for some sources */
#ifdef HAVE_QIO
  QIO_Reader *infile;
  QIO_Writer *outfile;
#endif
  int file_initialized;
  complex *c_src;     /* Pointer for complex source field storage */
  wilson_vector *wv_src; /* Pointer for wilson vector source field storage */
  int ksource;        /* Counter for a list of sources */
  int src_pointer ;   /* smearing function (for the moment, only
		         clover_finite_p_vary/create_wilson_source.c) */
} wilson_quark_source;

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

/* dirac_utilities.c */
double start_timing(void);
void print_timing(double dtime, char *str);
void clear_swv_field(spin_wilson_vector *swv);
void clear_wp_field(wilson_prop_field wp);
void clear_wv_field(wilson_vector *wv);
spin_wilson_vector *create_swv_field(void);
wilson_prop_field create_wp_field(void);
wilson_vector *create_wv_field(void);
void copy_wp_from_wv(wilson_prop_field wp, wilson_vector *wv, 
		     int color, int spin);
void copy_wp_field(wilson_prop_field wpcopy, wilson_prop_field wp);
void copy_wv_from_swv(wilson_vector *wv, spin_wilson_vector *swv, int spin);
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field wp, 
		     int color, int spin);
void copy_wv_from_wprop(wilson_vector *wv, wilson_propagator *wprop, 
			int color, int spin);
wilson_prop_field create_wp_field_copy(wilson_prop_field w);
void destroy_swv_field(spin_wilson_vector *swv);
void destroy_wv_field(wilson_vector *wv);
void destroy_wp_field(wilson_prop_field wp);
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field wp, int color);

/* gammas.c */
void mult_w_by_gamma(wilson_vector * src, wilson_vector * dest, int dir);
void mult_sw_by_gamma_l(spin_wilson_vector * src,spin_wilson_vector * dest, int dir);
void mult_sw_by_gamma_mat_l(spin_wilson_vector * src, 
			    spin_wilson_vector * dest, 
			    gamma_matrix_t *gm);
void mult_sw_by_gamma_r(spin_wilson_vector * src,spin_wilson_vector * dest, int dir);
void mult_sw_by_gamma_mat_r(spin_wilson_vector * src,
			    spin_wilson_vector * dest, 
			    gamma_matrix_t *gm);
void mult_gamma_by_gamma(gamma_matrix_t *g1, gamma_matrix_t *g2, 
			 gamma_matrix_t *g3);
gamma_matrix_t gamma_mat(enum gammatype i);
void gamma_adj(gamma_matrix_t *dest, gamma_matrix_t *src);
void gamma_transp(gamma_matrix_t *dest, gamma_matrix_t *src);
void gamma_conj(gamma_matrix_t *dest, gamma_matrix_t *src);
int gamma_index(char *label);
char *gamma_label(int index);

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

/* staggered2naive.c */

void convert_ksprop_to_wprop(wilson_propagator *wp, su3_matrix *ksp,
			     int ks_source_r[]);
void convert_ksprop_to_wprop_swv(spin_wilson_vector *swv, 
				 su3_vector *ksp, int r[]);

/* w_source.c */
void alloc_wqs_wv_src(wilson_quark_source *wqs);
void alloc_wqs_c_src(wilson_quark_source *wqs);
int ask_output_w_quark_source_file( FILE *fp, int prompt, 
				    int *flag, int *source_type,
				    int *t0, char *descrp, char *filename);
int ask_w_quark_source( FILE *fp, int prompt, int *type, char *descrp );
int choose_usqcd_w_file_type(int source_type);
void clear_wqs(wilson_quark_source *wqs);
int get_w_quark_source( FILE *fp, int prompt, wilson_quark_source *wqs );
int get_w_quark_sink(FILE *fp, int prompt, wilson_quark_source *wqs);
void init_wqs(wilson_quark_source *wqs);
void r_close_w_source(wilson_quark_source *wqs);
void r_open_w_source(wilson_quark_source *wqs);
void w_close_w_source(wilson_quark_source *wqs);
void w_open_w_source(wilson_quark_source *wqs, char *fileinfo);
int w_source_site(field_offset src, wilson_quark_source *wqs);
int w_source_field(wilson_vector *src, wilson_quark_source *wqs);
int w_source_write(wilson_vector *src, wilson_quark_source *wqs);
void w_sink_site(field_offset snk, wilson_quark_source *wqs);
void w_sink_field(complex *snk, wilson_quark_source *wqs);

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
void dslash_w_3D_site( field_offset src, field_offset dest, 
		       int isign, int parity);
void dslash_w_3D_field( wilson_vector *src, wilson_vector *dest, 
			int isign, int parity);

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
    int (*source_func_field)(wilson_vector *src, 
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
// void meson_cont_mom(
//   complex ***prop,          /* where result is stored */
//   spin_wilson_vector *src1, /* quark propagator */
//   spin_wilson_vector *src2, /* quark propagator */
//   int no_q_momenta,         /* number of q values */
//   int q_momstore[][3],      /* q values themselves */
//   int mom_index[],          /* momentum hash */
//   int no_gamma_corr,        /* number of meson types */
//   int meson_index[],        /* meson hash */
//   int gin[],                /* Gamma matrix type for source */
//   int gout[],               /* Gamma matrix type for sink */
//   complex meson_phase[],    /* phase factor to apply to correlator */
//   int **do_corr,            /* do_corr[m][p] = 1 to compute corr */
//   char **corr_parity        /* parity of correlator */
//		    );

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
  int gin[],                /* gin[c] is the source gamma */
  int gout[],               /* gout[c] is the sink gamma */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_fact[],        /* meson_fact[c] scales the correlator */
  int corr_index[]          /* m = corr_index[c] is the correlator index */
		    );


int decode_phase(char *label);

/* w_baryon.c */
void w_baryon(wilson_prop_field src1,wilson_prop_field src2,
	      wilson_prop_field src3, complex *prop[4]);
void w_baryon_hl(wilson_prop_field src1,wilson_prop_field src2,
		 wilson_prop_field src3, complex *prop[6]);

#endif /* _GENERIC_WILSON_H */

