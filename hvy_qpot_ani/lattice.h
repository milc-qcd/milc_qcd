#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 2
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "hvy_qpot_includes.h"    
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
        /* The state information for a random number generator */
        double_prn site_prn;
        /* align to double word boundary (kludge for Intel compiler) */
        int space1;

    /* Now come the physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4] ALIGNMENT;

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
EXTERN  size_t volume;			/* volume of lattice = nx*ny*nz*nt */
EXTERN  Real u0;

/*  details of the correlations to be considered */
EXTERN  short cor_dir;
EXTERN  int min_ct;
EXTERN  int max_x, max_y, max_z, max_t;
EXTERN  int xc[NDIRS];
EXTERN  int nc[4];
EXTERN  int maxc[4];
EXTERN int locx[4];
#ifdef ANISOTROPY
EXTERN  short ani_dir; /* direction of anisotropy */
EXTERN  Real ani_xiq; /* bare quark anisotropy */
#  ifdef ONEDIM_ANISO_TEST
EXTERN  Real iso_xiq; /* bare quark isotropic link factor for debugging */
#  endif
#endif
EXTERN  Real max_r;
EXTERN  int off_axis_flag;
EXTERN  int hqp_alg; /* choosing an algorithm for the Coulomb gauge hvy_pot */

/*  details of the smearing ot be used */
EXTERN	int no_smear_level,smear_num[5];
EXTERN  int tot_smear;  /* running total of smearing steps for lattice */
#if (defined HYP_3D_SMEARING || defined HYP_4D_SMEARING)
EXTERN  Real hyp_alpha1; /* parameters for HYP smearing */
EXTERN  Real hyp_alpha2;
EXTERN  Real hyp_alpha3;
#else  /* APE smearing */
EXTERN  Real staple_weight;
EXTERN  int ape_iter;
#endif
EXTERN  Real smear_fac;
#if (defined APE_1D_SMEARING || defined APE_1D2_SMEARING)
EXTERN  int stap_dir;
#endif

/*  details of starting and saving */
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN  int fixflag;    /* gauge fix flag */
//EXTERN	int total_iters;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
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
