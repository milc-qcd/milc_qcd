#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_cl.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

/* Modifications:
   8/13/97 wilson_propagator and spin_wilson_vector to su3.h C.D.
   8/10/96  Changed handling of propagator file names; changed same for param C.D.
*/

#include "defines.h"
#include "../include/generic_wilson.h"
#include "../include/generic_quark_types.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/macros.h"   /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/generic_clover.h" /* For clover */

/* Begin definition of site structure */

#include "../include/su3.h"

typedef struct {
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

        /* wilson half vector (temporary used in dslash_w_site) */
        half_wilson_vector htmp[MAXHTMP];
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
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
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
EXTERN	int niter,nrestart;
#define MAX_KAP 6
EXTERN	Real kappa,source_r0,kap[MAX_KAP],resid[MAX_KAP],relresid[MAX_KAP];
EXTERN	Real clov_c,u0;
EXTERN	int num_kap;		/* max number of kappa's <= MAX_KAP */
EXTERN	char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN  char savefile_w[MAX_KAP][MAXFILENAME], 
  startfile_w[MAX_KAP][MAXFILENAME];
EXTERN  char scratchstem_w[MAXFILENAME];
EXTERN  int scratchflag;        /* scratch file mode: SAVE_SERIAL, SAVE_CHECKPOINT */
EXTERN	int startflag;		/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;		/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int startflag_w[MAX_KAP];  /* beginning wilson: 
			   RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_KAP];   /* save propagator: SAVE_ASCII, SAVE_BINARY
				   SAVE_PARALLEL */
EXTERN	int total_iters;
EXTERN  char spectrum_request[MAX_SPECTRUM_REQUEST]; /* request list for spectral measurements */
EXTERN  Real sink_r0;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN wilson_quark_source wqs;
EXTERN wilson_quark_source wqstmp;  /* Temporary */

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN  double_prn node_prn ;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif /* _LATTICE_H */
