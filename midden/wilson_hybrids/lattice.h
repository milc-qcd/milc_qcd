#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 6
   wilson_hybrids program
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h" /* For quark_invert_control, etc */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/random.h"    /* For double_prn */
#include "../include/io_lat.h"    /* For gauge_file */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
struct site {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	int index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

    /* Now come the physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4];

	/* Use unions so smeared gauge fields overlap wilson vectors.
	   Programmer must use only one at a time!!  In particular, this
	   means that you must do all your smearing, and compute the
	   field strength, before doing any conjugate gradients, since
	   the CG will overwrite the smeared links. */
	union {
	    wilson_vector Uvec[5];	/* cg vectors */
	    su3_matrix Umat[6];	/* smearlinks, tempmat1 and tempmat2 */
	} U;
	union {
	    wilson_vector Vvec[3];	/* */
	    su3_matrix Vmat[4];	/* */
	} V;

	/* Wilson complex vectors */
#define	      g_rand  U.Uvec[0]	/* gaussian random vector */
#define	      psi  U.Uvec[1]	/* solution vector */
#define	      chi  U.Uvec[2]	/* source vector */
#define	      mp   U.Uvec[3]	/* another CG vector */
	wilson_vector ttt;	/* internal conjugate gradient vector */
#ifdef BICONGRAD
#define		      vtmp  U.Uvec[4]	/* more cg vectors for bicongrad */
	wilson_vector sss;	/* more cg vectors for bicongrad */
#endif
#ifndef LEAN
        wilson_vector p;        /* conjugate gradient change vector */
	wilson_vector r; 	/* residue */
#endif
        /* wilson half vector (temporary used in dslash_w) */
#ifndef LEAN
        half_wilson_vector htmp[8];
#else
        half_wilson_vector htmp[2];
#endif

#ifdef SMEAR
	/* smeared links, and temporary matrix */
#define		   smearlink  U.Umat	/* su3_matrix[4] */
#define		   templink   V.Vmat	/* su3_matrix[4] */
#define		   tempmat1 U.Umat[4]	/* su3_matrix */
#define		   tempmat2 U.Umat[5]	/* su3_matrix */
#endif
        /* components of F_mu_nu, use defines to give  components names */
        su3_matrix field_strength[6];
        /* quark and antiquark propagators */
#define		quark_source V.Vvec[0]	/* wilson_vector */
#define		quark_prop V.Vvec[1]	/* wilson_vector */
#define		anti_prop V.Vvec[2]	/* wilson_vector */

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
EXTERN	int iseed;		/* random number seed */
EXTERN	int niter;
EXTERN	int source_start, source_inc, n_sources;
	/* source time, increment for it, and number of source slices */
EXTERN	Real rsqprop,beta,kappa;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN	int startflag;	/* beginning lattice flag */
EXTERN  int saveflag;	/* end lattice flag */
EXTERN	int total_iters;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN quark_invert_control qic;
EXTERN dirac_wilson_param dwp;
EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8
EXTERN char ** gen_pt[N_POINTERS];

#endif /* _LATTICE_H */
