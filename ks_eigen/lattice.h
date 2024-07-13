#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 4
   This file defines global scalars and the fields in the lattice.

   Directory for dynamical improved KS action.  Allow:
	naik term in fermion action
	fat links in fermion action
	arbitrary paths in gauge action (eg Symanzik imp.)

   If "FN" is defined,
     Includes storage for Naik improvement (longlink[4], templongvec[4],
     gen_pt[16], etc.
     Includes storage for "fat links"  (fatlink[4])
*/


#include "defines.h"
#include "params.h"
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For ferm_links_t and ks_action_paths */
#include "../include/random.h"
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/fermion_links.h"

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

/* The lattice is an array of sites.  */
#define MOM_SITE   /* If there is a mom member of the site struct */
typedef struct {
    /* The first part is standard to all programs */
	/* coordinates of this site */
	short x,y,z,t;
	/* is it even or odd? */
	char parity;
	/* my index in the array */
	uint32_t index;
#ifdef SITERAND
	/* The state information for a random number generator */
	double_prn site_prn;
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;
#endif

/* ------------------------------------------------------------ */
/*   Now come the physical fields, program dependent            */
/* ------------------------------------------------------------ */
	/* gauge field */
	su3_matrix link[4] ALIGNMENT;	/* the fundamental field */
#ifdef HMC_ALGORITHM
 	su3_matrix old_link[4];
	/* For accept/reject */
#endif

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4] ALIGNMENT;

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

	/* 3 element complex vectors */
 	su3_vector phi;		/* Gaussian random source vector */
 	su3_vector resid;	/* conjugate gradient residual vector */
 	su3_vector cg_p;	/* conjugate gradient change vector */
 	su3_vector xxx;		/* solution vector = Kinverse * phi */
 	su3_vector ttt;		/* temporary vector, for K*ppp */
 	su3_vector g_rand;	/* Gaussian random vector*/
	/* Use trick of combining xxx=D^adj D)^(-1) on even sites with
	   Dslash time this on odd sites when computing fermion force */
	
#ifdef PHI_ALGORITHM
 	su3_vector old_xxx;	/* For predicting next xxx */
#endif
#ifdef SPECTRUM
	su3_vector propmat[3];	/* For three source colors */
	su3_vector propmat2[3];	/* nl_spectrum() */
	su3_matrix tempmat2;
	/* for spectrum_imp() */
	/**su3_vector quark_source, quark_prop, anti_prop;**/
#define quark_source propmat2[0]
#define quark_prop propmat2[1]
#define anti_prop propmat2[2]
#endif

	/* temporary vectors and matrices */
	su3_vector tempvec[4];	/* One for each direction */
#ifdef FN
	su3_vector templongvec[4];	/* One for each direction */
        su3_vector templongv1;
#endif
	su3_matrix tempmat1,staple;
} site;


/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* The following are global scalars 
   beta is overall gauge coupling factor (6/g^2)
   mass is quark mass
   u0 is tadpole improvement factor, perhaps (plaq/3)^(1/4)
*/
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  size_t volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	uint32_t iseed;		/* random number seed */
#ifdef FIX_NODE_GEOM
EXTERN  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimensions */
#ifdef FIX_IONODE_GEOM
EXTERN int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
EXTERN	int niter,nrestart,nflavors;
EXTERN  Real u0;
EXTERN	Real rsqmin,rsqprop;
EXTERN  params param;           /* user input parameters */
EXTERN	int total_iters;
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
        /* source time, increment for it, and number of source slices */

/*******/

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

EXTERN  gauge_file *savelat_p;
EXTERN  gauge_file *start_u1lat_p;
EXTERN  char utc_date_time[64];
EXTERN  char hostname[128];

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN su3_matrix *ape_links;
EXTERN int refresh_ape_links;
EXTERN int ape_links_ks_phases;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 16
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN fermion_links_t    *fn_links;

EXTERN int n_naiks;	/* Only one (dummy) Naik correction (for HISQ) */
EXTERN double eps_naik[1];	/* Naik correction */

EXTERN Real *u1_A;
EXTERN Real g_splaq,g_tplaq;	/* global U(1) plaquette measures */

/* For eigenpair calculation */
EXTERN int Nvecs_tot;
EXTERN Real *eigVal; /* eigenvalues of D^dag D */
EXTERN su3_vector **eigVec; /* eigenvectors */

#endif /* _LATTICE_H */
