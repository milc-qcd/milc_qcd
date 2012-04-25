#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 2
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/random.h"  /* For double_prn */
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"

/* The lattice is an array of sites.  */
struct site {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;

    /* Now come the physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4];

  /*	su3_matrix s_link,s_link_f,t_link_f; */
#ifndef COULOMB
        su3_matrix diag,staple,tempmat1,tempmat2;
#endif
};
typedef struct site site;


/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;			/* volume of lattice = nx*ny*nz*nt */
EXTERN  Real u0;
EXTERN  int max_t, max_x;
EXTERN  Real staple_weight;
EXTERN  int ape_iter;
EXTERN	int no_smear_level,smear_num[5],off_axis_flag;
EXTERN  int tot_smear;  /* running total of smearing steps for lattice */
EXTERN	Real smear_fac;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN  int fixflag;    /* gauge fix flag */
EXTERN	int total_iters;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8
EXTERN char ** gen_pt[N_POINTERS];

#endif /* _LATTICE_H */
