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
#include <stdint.h>

/* Begin definition of site structure */

/* The lattice is an array of sites.  */
typedef struct {
  /* The first part is standard to all programs */
  /* coordinates of this site */
  short x,y,z,t;
  /* is it even or odd? */
  char parity;
  /* my index in the array */
  uint32_t index;
  /* The state information for a random number generator */
  double_prn site_prn;

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
  su3_matrix link[4] ALIGNMENT; /* gauge field */
#if GF_INTEGRATOR==INTEGRATOR_RKMK3 || GF_INTEGRATOR==INTEGRATOR_RKMK4 || \
    GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  /* gauge field at the beginning of RK step */
  su3_matrix link0[4] ALIGNMENT;
#endif
  /* Temporary matrices for staple, smoothing, and field strength */
  su3_matrix staple[4]; /* staple for each link */
#if GF_INTEGRATOR==INTEGRATOR_RKMK3
  anti_hermitmat K[3][4]; /* right-hand-side in RK method for all stages */
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4
  anti_hermitmat K[4][4]; /* right-hand-side in RK method for all stages */
#elif GF_INTEGRATOR==INTEGRATOR_RKMK5
  anti_hermitmat K[6][4]; /* right-hand-side in RK method for all stages */
#elif GF_INTEGRATOR==INTEGRATOR_RKMK8
  anti_hermitmat K[13][4]; /* right-hand-side in RK method for all stages */
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
  anti_hermitmat K[3][4]; /* right-hand-side in RK method for all stages */
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  anti_hermitmat K[4][4]; /* right-hand-side in RK method for all stages */
#endif
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
EXTERN  size_t volume;
#ifdef ANISOTROPY
EXTERN  Real ani;
#endif
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
/* Integrator parameters */
// Maximum number of stages (for storing coefficients)
#define MAX_RK_NS 13
// number of stages
EXTERN int N_stages;
// order of the method
EXTERN int p_order;
/* 2N-storage schemes */
#if GF_INTEGRATOR==INTEGRATOR_LUSCHER || GF_INTEGRATOR==INTEGRATOR_CK \
 || GF_INTEGRATOR==INTEGRATOR_BBB || GF_INTEGRATOR==INTEGRATOR_CF3
// A, B coefficients
EXTERN Real A_2N[MAX_RK_NS];
EXTERN Real B_2N[MAX_RK_NS];
/* RKMK schemes */
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3 || GF_INTEGRATOR==INTEGRATOR_RKMK4 || GF_INTEGRATOR==INTEGRATOR_RKMK5 || GF_INTEGRATOR==INTEGRATOR_RKMK8
// RK coefficients in Butcher table
EXTERN Real a_RK[MAX_RK_NS][MAX_RK_NS];
EXTERN Real b_RK[MAX_RK_NS];
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
// A, B coefficients
EXTERN Real A_2N[MAX_RK_NS];
EXTERN Real B_2N[MAX_RK_NS];
// second-order coefficients
EXTERN Real Lambda[3];
// rejeted steps in adaptive integrators
EXTERN int steps_rejected;
// local tolerance for adaptive integrators
EXTERN Real local_tol;
// distance between two approximations in adaptive schemes
EXTERN Real dist;
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
// RK coefficients in Butcher table
EXTERN Real a_RK[MAX_RK_NS][MAX_RK_NS];
EXTERN Real b_RK[MAX_RK_NS];
// rejeted steps in adaptive integrators
EXTERN int steps_rejected;
// local tolerance for adaptive integrators
EXTERN Real local_tol;
// distance between two approximations in adaptive schemes
EXTERN Real dist;
// to use FSAL property permute indices in the storage array
EXTERN int indK[4];
EXTERN int is_first_step;
#endif
// flag if the integration step is final
EXTERN int is_final_step;

EXTERN int startflag;
EXTERN int saveflag;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/


/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;      /* number of sites on this node */
EXTERN	size_t even_sites_on_node; /* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;  /* number of odd sites on this node */
EXTERN	int number_of_nodes;    /* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;
EXTERN char hostname[128];

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Temporary su3 matrices for gathers */
#define N_TEMPORARY 7
EXTERN su3_matrix *tempmat[N_TEMPORARY];

/* Generic pointers, for gather routines */
#define N_POINTERS 9
EXTERN char ** gen_pt[N_POINTERS];

#endif /* _LATTICE_H */
