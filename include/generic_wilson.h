#ifndef _GENERIC_WILSON_H
#define _GENERIC_WILSON_H
/************************ generic_wilson.h ******************************
*									*
*  Macros and declarations for generic_wilson routines                  *
*  This header is for codes that call generic_wilson routines           *
*  MIMD version 6 							*
*									*
*/

#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/macros.h"
#include "../include/generic_quark_types.h"

/* For various inversion routines.  Not used sytematically yet. CD */
enum guess_params { START_ZERO_GUESS = 0 ,  START_NONZERO_GUESS } ;  

int wilson_invert(      /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (source already created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    field_offset sav,   /* type wilson_vector (for saving source) */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp            /* Passthrough Dirac matrix parameters */
    );

int wilson_invert_lean( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
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

void boundary_flip(int sign );
int congrad_w(int niter,Real rsqmin,Real *final_rsq_ptr);

void copy_site_wilson_vector(field_offset src, field_offset dest);

int cgilu_w(             /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int bicgilu_w(          /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

int mrilu_w_or(          /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    );

/* For quark source routines */
/* The Weyl representation types are included for w_source_h */
enum source_type { 
  POINT = 1, GAUSSIAN, CUTOFF_GAUSSIAN,
  POINT_WEYL, CUTOFF_GAUSSIAN_WEYL } ;

/* w_source.c */
void w_source(field_offset src,wilson_quark_source *wqs);
void w_sink(field_offset snk,wilson_quark_source *wqs);
void w_sink_scalar(field_offset snk,wilson_quark_source *wqs);
int ask_quark_source( int prompt, int *type, char *descrp );

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

void cleanup_tmp_links();
void meson_cont(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		complex *prop);
void w_meson(field_offset src1,field_offset src2,complex *prop[10]);
void d_w_meson(field_offset src1,field_offset src2,double_complex *prop[10]);
void w_baryon(field_offset src1,field_offset src2,field_offset src3,
	      complex *prop[4]);
void w_baryon_hl(field_offset src1,field_offset src2,
		 field_offset src3, complex *prop[6]);
#endif /* _GENERIC_WILSON_H */

