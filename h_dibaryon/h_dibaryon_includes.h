/****************** h_dibaryon_includes.h ******************************/
/*
*  Include files for the clover H-dibaryon application
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
#include "../include/comdefs.h"	/* definitions and variables for communications */
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"

typedef struct {
  int dq[3];
  int wt;
} diqkwf;

typedef struct {
  int dq;   /* diquark color-spin index */
  int cs;   /* color-spin index */
  int wt;   /* weight */
} twoplusoneqkwf;

typedef struct {
  int cs[3]; /* color-spin indices */
  int wt;    /* weight */
} threeqkwf;

typedef struct {
  complex *c;
  char *label;
  int done;
} propagator;

/* prototypes for functions in high level code */
int setup_H_cl();
int readin(int prompt);
void diquarkprop(field_offset qk, field_offset diqk);
void nr_propagator(field_offset qk, field_offset nr, int forw_back);  
void w_hdibaryon(field_offset src_u, field_offset src_s, propagator prop[]);
void w_nrbaryon(field_offset src_1, field_offset src_2, 
		field_offset src_2di, propagator prop[]);
void lam_lam(diqkwf *wf, int *nterm, Real *norm);
void nuc_xi(diqkwf *wf, int *nterm, Real *norm);

#define MAXCHANNEL 3
#define MAXWF 200
#define MAXLABEL 5
void make_channel_wfs(diqkwf channel_wf[MAXCHANNEL][MAXWF],
		      int channel_terms[MAXCHANNEL],
		      char channel_label[MAXCHANNEL][MAXLABEL], 
		      Real channel_norm[MAXCHANNEL],
		      int *nchannel);
void make_nucleon_wf(twoplusoneqkwf wf[],int *terms, Real *norm);
void make_lambda_wf(threeqkwf wf[], int *terms, Real *norm);
void w_nrbaryon(field_offset src_1, 
		field_offset src_2, field_offset src_2di, 
		propagator prop[]) ;
/*  h_dibaryon_includes.h */




