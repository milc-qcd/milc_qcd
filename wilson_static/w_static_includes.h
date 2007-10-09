/****************** w_static_includes.h  ******************************/
/*
 *  Include files for the wilson_static application
 */

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "params.h"
#include "../include/comdefs.h"	/* definitions and variables for communications */
#include "../include/io_lat.h"
#include "../include/io_ksprop.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif

/*************************************************** 
   prototypes for functions in high level code 
   ***************************************************/

void wall_sink_h(field_offset src,wilson_vector *chi_out,
		 int parity,int keep_parity,int x0,int y0,int z0)  ; 

void w_meson_hop(double *prop[], field_offset  heavy_quark, field_offset  light_quark, 
		 wilson_vector *heavy_wall, wilson_vector *light_wall, 
		 int parity, int spin, int wallflag)  ;


void hopping(field_offset src, field_offset temp, 
	     field_offset light_quark,
	     int nhop,Real kappa_c,int parity_of_source,
	     int color, int spin, int wallflag, FILE  *fp_m_out, int fb_m_out) ;


void sum(Real kappa_h, Real kappa_c, int nhop, int file_nhop,
	   FILE * fp_m_in, int fb_m_in, FILE * fp_k_out, int writeflag) ;


void sum_light(int file_nhop, FILE *fp_m_in, int fb_m_in,
	       FILE *fp_k_out) ;


void mult_by_gamma_right_vec( wilson_vector *src, wilson_vector *dest, int dir) ;

void light_meson(field_offset light_quark,int color,int spin,
		 int wallflag, FILE  *fp_m_out, int fb_m_out) ;


FILE *(r_ascii_m_i(char *filenam,int i1,int *file_hops))  ;

void r_ascii_m(FILE *fp,int spin,int color,int iter,double *prop[]) ;

void r_ascii_m_f(FILE *fp, char *filenam) ;

int r_binary_m_i(char *filenam,int i1,int *file_hops)  ;
void r_binary_m(int fp, int spin, int color, int iter,double *prop[]) ;
void r_binary_m_f(int fp,char *filenam) ;

FILE *(w_ascii_m_i(char *filenam,int i1)) ;
FILE *(a_ascii_m_i(char *filenam,int i1)) ;
void w_ascii_m(FILE *fp,int spin, int color,int iter,double *prop[] ) ;
void w_ascii_m_f(FILE *fp, char *filenam) ;

int w_binary_m_i(char *filenam, int i1) ;
int a_binary_m_i(char *filenam,int i1) ;
void w_binary_m(int fp, int spin, int color, int iter, double *prop[]) ;
void w_binary_m_f(int fp, char *filenam) ;

int  setup_h()  ;
int get_smearing_funcs_code(char savebuf[], params *par_buf ) ;
int initial_set() ;
void make_lattice() ;
int readin(int prompt)  ;

void kappa_dslash(field_offset phi, field_offset temp, Real kappa, int parity);

void dagger_quark_prop(field_offset conj_quark, field_offset quark);
void contract_a_src(Real *vary_matrix, int block_pt, int block_smear[4] );

void setup_static_prop();
void setup_static_prop_bparam();

void fft_smeared_wline();
void fft_negmom_quark();

void write_vary_matrix(Real *vary_matrix) ;


void setup_vary(complex *meson, int nodata);
void calc_smeared_meson(complex *meson, 
			field_offset quark,field_offset quark_smear,
			int colour, int spin);
void write_smear_mesonx(complex *meson);
void calc_vary_matrix();

void c_scalar_wilsonvec(wilson_vector *m, complex *phase);
void mult_wilson_vec_matdag(wilson_vector *out_vec,
			    wilson_vector *in_vec, su3_matrix *m);
void create_local_oper(int which_smear, char filename[]);
void create_wall_oper(int which_smear, char filename[]);
void create_expon_oper(int which_smear,   Real decay , char filename[]);
void create_2S_polyexpon_oper(int which_smear, Real decay, Real A, char filename[]);
void create_3S_polyexpon_oper(int which_smear, Real decay, Real B, Real D,
char filename[]);
void read_expon_params(Real *decay, char savebuf[]  ) ;
void read_twoS_params(Real *decay, Real *A, char savebuf[]  ) ;
void read_threeS_params(Real *decay, Real *B, Real *D, char savebuf[]  ) ;

void load_smearing();
void buildup_strip(field_offset quark, int colour, int spin);
void   setup_timeslice_fft();
void  zero_strip_quark();
void   static_prop() ; 
void create_smearing_funcs(void) ;


/********** end of the function prototypes ***/

/* end of w_static_includes.h  */

