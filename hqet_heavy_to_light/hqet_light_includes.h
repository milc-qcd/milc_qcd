/****************** hqet_light_includes.h ******************************/
/*
*  Include files for HQET-light application
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
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_form.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../generic_form/gammatypes.h"
#include "../include/dirs.h"

/* some macros that are useful  **/

/*** memory location for the form factor operators ******/

#define HQET_FORM_WHERE(t,zonk_pt,spect_pt,q_pt,vel_pt,oper_pt) \
         t + nt*(zonk_pt + no_zonked_light*(spect_pt + no_spectator*(q_pt + no_q_values*(vel_pt + novel*oper_pt ))))


#define TWOPT_FORM_WHERE(t,zonk_pt,spect_pt,q_pt,oper_pt) \
         t + nt*(zonk_pt + no_zonked_light*(spect_pt + no_spectator*(q_pt + no_q_values*oper_pt )))

#define TWOPT_SEQ_WHERE(t,spect_pt,v_pt) \
         t + nt*(spect_pt + no_spectator*( v_pt ))


/** debug ouput statement ***/
#define IF_VERBOSE_ON(ff)  if(this_node==0 && verbose_flag >= ff )

/* prototypes for functions in high level code */
int setup_hqet_form();
int readin(int prompt);
int load_velocity_from_disk(params *par_buf, char filename[], int maxvel);

void sub3_su3_matrix(su3_matrix *a, su3_matrix*b, su3_matrix*c, su3_matrix *d);
void apply_hqet_proj(field_offset ans, int v_pt);

void twopt_sequential_corr(complex *corr, field_offset hqet_prop,
			   int v_pt, int spect_pt, int color) ;
void smeared_sequential_source(field_offset ans,
			       field_offset hqet_prop, 
			       int tsrc, int v_pt) ;
void contract_hqet_to_light(complex *corr, field_offset , field_offset , 
			    field_offset ,int vel_pt,  int zonked_pt, 
			    int spect_pt) ;
void contract_light_twopt(complex *corr, field_offset q_zonked, 
			  field_offset q_sequential,
			  int zonked_pt, int spect_pt) ;

void setup_hlcorr(complex **hl_corr, int *hl_corr_dim ) ;
void finish_hlcorr(complex *hl_corr, int hl_corr_dim ) ;

void setup_twocorr(complex **two_corr, int *two_corr_dim ) ;
void finish_twocorr(complex *two_corr, int two_corr_dim ) ;

void setup_seq_corr(complex **seq_corr, int *seq_corr_dim ) ;
void finish_seqcorr(complex *seq_corr, int seq_corr_dim ) ; 

void clover_rotate(field_offset src, field_offset dest);

void copy_lattice_spin_wilson_vector(field_offset out, field_offset in) ; 

void calc_hqet_light_form(void) ;
void setup_control_hqet_form(void) ;

void write_hqet_form_corr(complex *corr, char filename[], int dim) ; 
void write_twopt(complex *corr, char filename[], int dim) ;
void write_seq_twopt(complex *corr, char filename[], int dim) ;

void generate_hqet_prop(field_offset hqet_prop, int tstart, int tend , 
			int tcurrent, int v1, int v2)   ;
void generate_hqet_prop_back(field_offset hqet_prop, int tstart, int tend , 
			     int tcurrent, int v1, int v2) ;
void evolve_hqet_forwards_and_backwards(field_offset hqet_prop, int tsrc, int v);
void sub3_su3_matrix(su3_matrix *a, su3_matrix*b, su3_matrix*c, su3_matrix *d) ;
void mult_su3_nn_z_inc(complex z,su3_matrix *a,su3_matrix *b,su3_matrix *c) ;
void mult_su3_an_z(complex z,su3_matrix *a,su3_matrix *b,su3_matrix *c );

void apply_hqet_proj(field_offset ans, int v_pt) ;
void setup_timeslice_fft() ;
void smear_hqet_prop(field_offset hqet_prop,int which_smear, int vpt);

/* hqet_light_includes.h */
