#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_w.h ********************************/

/* include file for MIMD version 7
   This file defines global scalars and the fields in the lattice. */

/** version 4 ==> version 5 port by mcneile ***/

#include "defines.h"
#include "../include/generic_wilson.h" 
#include "../include/generic_quark_types.h" /* For wilson_quark_source */
#include "../include/random.h"  /* For double_prn */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/generic_clover.h" /* For clover */

/* Begin definition of site structure */

#include "../include/su3.h"

/* The lattice is an array of sites.  Within each node the even sites will
   be stored first, then the odd sites.
*/


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

	/* wilson complex vectors */
 	wilson_vector psi;   /* solution vector for mrilu, cgilu or jacobi */
 	wilson_vector chi;	 /* source vector */
 	wilson_vector mp;		/* another CG vector
--- used as a temp vector for mrilu and jacobi */

	/* wilson half vector (temporary used in dslash) */
#define MAXHTMP 8
	half_wilson_vector htmp[MAXHTMP];

        su3_matrix staple;   /* so we can calculate the plaquette */

} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/**** Global parameters as defines ***/

#ifdef CONTROL
int n_spins = 4;                /* Number of spins generated */
int source_loc[4] = {0,0,0,0};  /* Source location */
int spins[4] = {0,1,2,3};       /* List of spins generated */
#else
extern int n_spins;
extern int source_loc[4];
extern int spins[4];
#endif


/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int  volume;	/* volume of lattice = nx*ny*nz*nt */
EXTERN  int nkap;  /****** number of kappa values ******/
#define MAX_NKAP 20 /* maximum number of kappa values */
EXTERN  Real cappa[MAX_NKAP];  /******kappa******/
EXTERN  Real kappa; /*nickname for cappa[] in kappa loop***/
EXTERN  int start_kap,start_spin,start_color; /******starting values of 
                                              kappa, spin and color 
                                              in their loops******/
EXTERN  int end_kap,end_spin,end_color; /******end values of 
                                              kappa, spin and color 
                                              in their loops******/
EXTERN	int niter,nrestart,nhop,flag;
EXTERN	Real rsqmin,rsqprop,beta,kappa_c,width;
EXTERN	char startfile[MAXFILENAME],savefile_w[MAX_NKAP][MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char startfile_w[MAX_NKAP][MAXFILENAME],savefile_m[MAX_NKAP][MAXFILENAME];
EXTERN  char savefile[MAXFILENAME];  /** what to do with the gauge file ***/
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int saveflag ;

EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int startflag_w[MAX_NKAP]; /* beginning wilson: START_ASCII, START_BINARY */
EXTERN	int saveflag_w[MAX_NKAP];	/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int saveflag_m;	/* save meson props: SAVE_ASCII, SAVE_BINARY */
EXTERN	int total_iters;
EXTERN  int wall_separation;
EXTERN  int source_parity;
EXTERN  int nchannels;  /* = NCHANNELS +2 if extra_sink chosen */

EXTERN Real *wall_template; /* pointer to array of wall heights */
EXTERN  Real byterevReal;  /* flag to bytereverse lattice */

EXTERN  int fixflag;    /* whether to gauge fix **/



/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN  int	sites_on_node;		/* number of sites on this node */
EXTERN  int	even_sites_on_node;	/* number of even sites on this node */
EXTERN  int	odd_sites_on_node;	/* number of odd sites on this node */
EXTERN  int	number_of_nodes;	/* number of nodes in use */
EXTERN  int  this_node;		/* node number of this node */

EXTERN wilson_quark_source wqs;
EXTERN wilson_quark_source wqstmp, wqstmp2;  /* Temporary */

EXTERN quark_invert_control qic;
EXTERN dirac_wilson_param dwp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];


/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif /* _LATTICE_H */
