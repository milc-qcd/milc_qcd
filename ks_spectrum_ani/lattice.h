#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

/* Modifications:

   5/30/07 Moved most input parameters into a global param structure. CD

*/

#include "defines.h"
#include "params.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/generic_ks.h" /* For ferm_links_t and ks_action_paths */
#include "../include/generic_clover.h" /* For clover */
#include "../include/fermion_links.h"  /* For fermion_links_t */

/* Begin definition of site structure */

#include "../include/su3.h"

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
	/* align to double word boundary (kludge for Intel compiler) */
	int space1;

    /* Now come the physical fields, program dependent */
	/* gauge field */
	su3_matrix link[4] ALIGNMENT;

	/* The Kogut-Susskind phases, which have been absorbed into
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

  /* WE WANT TO GET RID OF THESE */

} site;

/* End definition of site structure */

/* Definition of globals */


#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  uint32_t iseed;
EXTERN  int niter, nrestart;
EXTERN  size_t volume;		/* volume of lattice = nx*ny*nz*nt */
#ifdef FIX_NODE_GEOM
EXTERN  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimension */
#endif
#ifdef FIX_SUBNODE_GEOM
EXTERN  int subnode_geometry[4];  /* Specifies fixed "nsubsquares" (i.e. 4D
			    hypercubes) for the PE ranks on each node in each
			    coordinate direction.  Must be divisors of
			    the node sublattice dimensions -- that is
			    full lattice dims divided by node_geometry */
#endif
#ifdef FIX_IONODE_GEOM
EXTERN int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
EXTERN  params param;           /* user input parameters */
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN	int total_iters;
EXTERN	int hisq_svd_counter;
EXTERN	int hisq_force_filter_counter;
EXTERN  Real u0, mass;  /* for Asqtad, etc. fermions only! */
EXTERN	Real rsqmin,rsqprop; /* for Asqtad, etc. fermions only! */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	size_t sites_on_node;		/* number of sites on this node */
EXTERN	size_t even_sites_on_node;	/* number of even sites on this node */
EXTERN	size_t odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/* Temporary for clover_info and ksprop_info */
EXTERN quark_source ksqstmp;
EXTERN ks_param ksptmp;
EXTERN quark_source wqstmp;
EXTERN dirac_clover_param dcptmp;
EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;
EXTERN gauge_file *start_u1lat_p;
EXTERN char utc_date_time[64];
EXTERN char hostname[128];

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN  double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN su3_matrix *ape_links;
EXTERN int refresh_ape_links;
EXTERN int ape_links_ks_phases;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define P_OFFSET 8
#define N_POINTERS 26+P_OFFSET	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN fermion_links_t        *fn_links;

#ifdef ANISOTROPY
EXTERN short ani_dir; /* direction of anisotropy */
EXTERN Real ani_xiq; /* bare quark anisotropy */
EXTERN Real ani_u0; /* anisotropic tadpole factor, for Asqtad fermions */
#  ifdef ONEDIM_ANISO_TEST
EXTERN Real iso_xiq; /* bare quark isotropic link factor for debugging */
#  endif
#endif
#ifdef CORDIR
EXTERN short is_permuted; /* permutation level of lattice directions */
EXTERN short cor_dir; /* time slice direction of correlation functions */
#endif


EXTERN Real *u1_A;
EXTERN Real g_splaq,g_tplaq;	/* global U(1) plaquette measures */

/* For eigenpair calculation */
EXTERN int Nvecs_tot;
EXTERN int Nvecs_alloc;
EXTERN Real *eigVal; /* eigenvalues of D^dag D */
EXTERN su3_vector **eigVec; /* eigenvectors */

#endif /* _LATTICE_H */
