#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
typedef struct {
  /* The first part is standard to all programs */
  /* coordinates of this site */
  short x,y,z,t;
  /* is it even or odd? */
  char parity;
  /* my index in the array */
  int index;
#ifdef SITERAND
  /* The state information for a random number generator */
/*
  double_prn site_prn;
*/
  /* align to double word boundary (kludge for Intel compiler) */
  int space1;
#endif
  
  /* Now come the physical fields, program dependent */
  /* gauge field */
  su3_matrix link[4];
  double ch_dens;
#ifdef HYP
  su3_matrix blocked_link[32];
#else
  su3_matrix fatlink[4];
#endif

  /* temporary matrices */
  su3_matrix tempmat1,tempmat2,staple;
} site;


/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars */
EXTERN  int nx,ny,nz,nt;        /* lattice dimensions */
EXTERN  int volume;                     /* volume of lattice = nx*ny*nz*nt */
EXTERN  int iseed;              /* random number seed */
EXTERN  Real ape_weight;       /* weight parameter in APE blocking */
EXTERN  int sweeps,hits,measinterval;
EXTERN  int total_sweeps;
EXTERN  char startfile[MAXFILENAME],savefile[MAXFILENAME],topofile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  int startflag;  /* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN  int saveflag;   /* do with lattice: 1=save; */
EXTERN  int savetopoflag; /* do with topo file */
#ifdef HYP
EXTERN Real alpha;
EXTERN Real alpha2;
EXTERN Real alpha3;
#endif

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN  int sites_on_node;              /* number of sites on this node */
EXTERN  int even_sites_on_node; /* number of even sites on this node */
EXTERN  int odd_sites_on_node;  /* number of odd sites on this node */
EXTERN  int number_of_nodes;    /* number of nodes in use */
EXTERN  int this_node;          /* node number of this node */

EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
/*
EXTERN  double_prn node_prn ;
*/

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8    /* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* for loop_inst.c */
#define nist 2
#define max_inst_length 16
#define inreps 2
#define max_num 400

/* global defns for general action */
EXTERN int inst_ind[nist][max_inst_length],inst_length[nist];
EXTERN int inst_table[nist][max_num][max_inst_length],inst_num[nist];
EXTERN int loop_char[max_num][2],ch;
EXTERN Real loop_coeff_inst[nist][inreps];
EXTERN int eps[nist][max_num];
#endif /* _LATTICE_H */
