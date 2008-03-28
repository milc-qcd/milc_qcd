#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

/* Modifications:

   5/30/07 Moved most input parameters into a global param structure. CD

*/

#include "params.h"
#include "../include/random.h"   /* For double_prn */
#include "../include/io_lat.h"    /* For gauge_file */

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

	/* The Kogut-Susskind phases, which have been absorbed into 
		the matrices.  Also the antiperiodic boundary conditions.  */
 	Real phase[4];

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
EXTERN  int niter, nrestart;
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN  params param;           /* user input parameters */
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN	int total_iters;
EXTERN  int phases_in; /* 1 if KS and BC phases absorbed into matrices */
EXTERN  Real u0, mass;  /* for Asqtad, etc. fermions only! */
EXTERN	Real rsqmin,rsqprop; /* for Asqtad, etc. fermions only! */

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

/* Temporary for clover_info and ksprop_info */
EXTERN ks_quark_source ksqstmp;
EXTERN ks_param ksptmp;
EXTERN wilson_quark_source wqstmp;
EXTERN dirac_clover_param dcptmp;
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
#define N_POINTERS 16	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* Storage for definition of the quark action */
EXTERN ferm_links_t        fn_links;
EXTERN ks_action_paths ks_act_paths;

#endif /* _LATTICE_H */
