#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice.
*/

#include "defines.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/su3.h"

/* Begin definition of site structure */

/* The lattice is an array of sites.  */
typedef struct {
  /* The first part is standard to all programs */
  /* coordinates of this site */
  short x,y,z,t;
  /* is it even or odd? */
  char parity;
  /* my index in the array */
  int index;

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
  su3_matrix link[4]; /* gauge field */

  /* Temporary matricies for staple, smoothing, and field strength */
  su3_matrix staple[4]; /* staple for each link */
  anti_hermitmat accumulate[4]; /* accumulation matrix for smearing */
  su3_matrix fieldstrength[6]; /* components of fmunu */

} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars */
/* Initialization parameters */
EXTERN	int nx,ny,nz,nt; 
EXTERN  int volume;
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;

/* Flow Parameters */
EXTERN  Real stepsize;
EXTERN  Real stoptime;
EXTERN  int total_steps;
EXTERN  int exp_order;
EXTERN  char flow_description[20];
EXTERN  int stapleflag;

EXTERN int startflag;
EXTERN int saveflag;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME]; 
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/


/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;      /* number of sites on this node */
EXTERN	int even_sites_on_node; /* number of even sites on this node */
EXTERN	int odd_sites_on_node;  /* number of odd sites on this node */
EXTERN	int number_of_nodes;    /* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Temporary su3 matricies for gathers */
#define N_TEMPORARY 7
EXTERN su3_matrix *tempmat[N_TEMPORARY];

/* Generic pointers, for gather routines */
#define N_POINTERS 9
EXTERN char ** gen_pt[N_POINTERS];

#endif /* _LATTICE_H */
