#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"	/* definitions and variables for communications */
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../generic_form/gammatypes.h"
#include "../include/generic_form.h"
#include "../include/dirs.h"

/* some macros that are useful  **/

/*** memory location for the form factor operators ******/

#define LIGHT_FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) \
         t + w_meson_nstore*(zonk_pt + no_zonked_light*(seq_pt + no_sequential*(spect_pt + no_spectator*(q_pt + no_q_values*(p_pt + no_p_values*oper_pt )))))

#define HEAVY_FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) \
         t + w_meson_nstore*(zonk_pt + no_zonked_heavy*(seq_pt + no_sequential*(spect_pt + no_spectator*(q_pt + no_q_values*(p_pt + no_p_values*oper_pt )))))


#define LL_TWOPT_FORM_WHERE(t,zonk_pt,spect_pt,k_pt,oper_pt) \
         t + w_meson_nstore*(zonk_pt + no_zonked_light*(spect_pt + no_spectator*(k_pt + no_k_values*oper_pt )))

#define HL_TWOPT_FORM_WHERE(t,zonk_pt,spect_pt,k_pt,oper_pt) \
         t + w_meson_nstore*(zonk_pt + no_zonked_heavy*(spect_pt + no_spectator*(k_pt + no_k_values*oper_pt )))




/**** these flags are used for the binary IO of the three point functions  ****/
enum  form_type { HEAVY_TO_HEAVY = 10 , HEAVY_TO_LIGHT  } ; 

/**** these flags are used for the binary IO of the two point functions  ****/
enum  two_form_type { HL_2PT_BAG = 20 , HL_REL_2PT  , LL_2PT , HL_LOCAL_SINK } ; 

enum meson_spec_choices { SETUP_CORR = 101 , CALCULATE_SPECTRUM , WRITE_RESULTS   } ; 

/** debug ouput statement ***/
#define IF_VERBOSE_ON(ff)  if(this_node==0 && verbose_flag >= ff )

/* prototypes for functions in high level code */
void local_wilson_source(field_offset src,int color,int spin) ;

void setup_w_meson_store();
void setup_HL3_corr(complex **corr, int *corr_dim , int *corr_stride) ;
void setup_HH3_corr(complex **corr, int *corr_dim, int *corr_stride );

void sequential_source(
  field_offset ans,          /* size wilson_vector
				on input: 1st propagator in sequence
				on output: result of inversion */
  field_offset work,         /* for building seq source - size wilson_vector */
  int px, int py, int pz,    /* momentum injected at sequential source */
  int t_final,               /* time slice of sequential source */
  int color,int spin,        /* color and spin of 1st source (unused) */
  Real Kappa,               /* for seq source inversion */
  int inverter_type,         /* HOPILU or CGILU */
  int MaxMR,                 /* MR or CG max iterations per restart */
  int nrestart,              /* MR or CG max restarts */
  Real RsdMR ,              /* MR or CG residual tolerance */
  int p_insert               /* (unused) */
  );

void contract_HH3(complex *corr, field_offset q_zonked, field_offset q_sequential,
field_offset q_rot, int p_pt, int seq_pt, int zonked_pt, int spect_pt) ;

void contract_HL3(complex *corr, field_offset q_zonked, field_offset q_sequential,
field_offset q_rot, int p_pt, int seq_pt, int zonked_pt, int spect_pt) ;


void setup_LL2_corr(complex **two_corr, int *two_corr_dim ) ;
void setup_HL2_corr(complex **two_corr, int *two_corr_dim ) ;
void setup_HL2_corr_with_rotations(complex **two_corr, int *two_corr_dim, int *two_corr_stride ) ;

void do_inversion(field_offset src,field_offset dest,
			   int spin, int color, Real kappa, int MaxMR, Real RsdMR )   ;

void contract_HL2(complex *corr, field_offset q_zonked, field_offset q_spectator,
int seq_pt, int spect_pt) ;

void contract_HL2_with_rotations(complex *corr, field_offset q_zonked, field_offset q_spectator,
field_offset q_rot, int zonked_pt, int spect_pt) ;

void contract_LL2(complex *corr, field_offset q_zonked, field_offset q_spectator,
int zonked_pt, int spect_pt) ; 


void contract_seq_twopt(complex *corr, field_offset quark, int spectate_pt, int seq_pt) ;

void finish_HL3_corr(complex *corr, int corr_dim , int copy_stride) ;
void finish_HH3_corr(complex *corr, int corr_dim , int copy_stride) ;

void finish_LL2_GG_corr(complex *corr, int corr_dim ) ;
void finish_HL2_GE_corr(complex *corr, int corr_dim ) ;
void finish_HL2_GG_corr(complex *corr, int corr_dim ) ;
void finish_HL2_GL_corr(complex *corr, int corr_dim , int copy_stride ) ;

void sink_smear_light_quark(field_offset quark, 
			    field_offset fft_work_a, field_offset fft_work_b, 
			    field_offset smear_func ) ; 

void sink_smear_light_pos_quark(field_offset quark, 
			    field_offset fft_work_a, field_offset fft_work_b, 
			    field_offset smear_func_mom ) ;

void sink_smear_light_pos_quark_large(field_offset quark, 
			    field_offset fft_work_a, field_offset fft_work_b, 
			    field_offset smear_func_mom )  ;

int readin(int prompt)  ;

void calc_heavy_light_form(void) ;
void setup_control_form(void) ;
int  setup_h(void)  ;

FILE *write_prop_form_header(complex *corr, char filename[], 
	     int hl_flag, int no_zonked, int numer_of_operators, int no_copies, int dim);
void write_corr_item(FILE *fp,complex *corr,char *filename);
void write_prop_form_close(FILE *fp,char *filename);


FILE *write_prop_twopt_header(complex *corr, char filename[], 
int hl_flag, int no_k_one , int no_k_two, int numer_of_operators, int no_copies,int dim) ;
void write_prop_twopt_close(FILE *fp, char *filename);


void setup_timeslice_fft() ;

int load_quark_read_info(int prompt, char savebuf[], char filename[], char mess[]) ;
void load_in_filename(char filename[] , int prompt, char match_string[], char prompt_message[] ) ;

void c_scalar_self_mult_wvec(wilson_vector *ans, complex s)  ;
void c_scalar_self_mult_spin_wvec(spin_wilson_vector *ans, complex s)  ;

void meson_spectrum(field_offset quark,int t_source, int ikap, int num_kap, int what_to_do,char *filename) ; 
void copy_site_spin_wilson_vector(field_offset src, field_offset dest) ;

void load_in_zonked_light(int color, int k_zonked_light)  ;
void load_in_zonked_light2(int color, int spin, int k_zonked_light,
			   field_offset dest);
void load_in_zonked_light_ssink(int color, int spin, int k_zonked_light,
			   field_offset dest);
void load_in_zonked_heavy_local(int color, int spin, int k_zonked_heavy,
			   field_offset dest);
void load_in_zonked_heavy_smear(int color, int spin, int k_zonked_heavy,
			   field_offset dest);
void generate_heavy_zonked(int color, int spin, 
	  Real Kappa, int inverter_type, field_offset dest);

void set_zonked_save_intermediate() ;
int set_reload_flag(int l_save_flag) ;

void clover_rotate_fermilab(field_offset src, field_offset dest) ;
void dslash_space(field_offset src,field_offset dest,int isign,int parity) ;
void clover_rotate(field_offset src, field_offset dest) ;

char *oper_name(int n )  ;
char *three_oper_name(int n , int copy_pt ) ;
char *oper_name_and_corrections(int n , int copy_pt ) ;

void dump_ll_twopt_text(complex *corr,  char* filename, int no_oper, Real norm) ;
void dump_hl_twopt_text(complex *corr, char* filename, int no_oper, char tag[], Real norm)  ;
void dump_hl_twopt_text_oper_corrections(complex *corr, char* filename, int no_oper, int copy_stride, char tag[], Real norm)  ;

void dump_hh_prop_form_text(complex *corr,  char* filename, int no_oper, int copy_stride, Real norm) ;
void dump_hl_prop_form_text(complex *corr,  char* filename, int no_oper, int copy_stride, Real norm) ;

/* prop_form_includes.h */






