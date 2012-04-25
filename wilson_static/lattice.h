#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice_w.h ********************************/
/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/lattice.h,v 1.12 2012/04/25 03:20:09 detar Exp $   **/
/* include file for MIMD heavy-light, version 4
   This file defines global scalars and the fields in the lattice. */

/** version 4 ==> version 5 port by DeTar ***/

#include "defines.h"
#include "../include/generic_wilson.h" 
#include "../include/generic_quark_types.h"  /* For quark_source */
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/random.h"    /* For double_prn */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/generic_clover.h" /* For clover */

/* Begin definition of site structure */

#include "../include/complex.h"
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

#define MAXHTMP 8
	/* wilson half vector (temporary used in dslash) */
	half_wilson_vector htmp[MAXHTMP];

        su3_matrix staple;   /* so we can calculate the plaquette */

     /****** variables for the static variational code ***/
	/* Quark propagator stripped (1 +/- g4)/2  **/
	su3_matrix strip_quark ; 

        /* Wilson line **/
	su3_matrix w_line ; 

        /* Smeared Wilson line, to save FFT's use a block **/
	su3_matrix smear_w_line[4] ; 
	su3_matrix fftwork[4] ;

	/* The smearing functions in position space ****/
#define MAX_SMEAR 10
	complex smear_func[MAX_SMEAR] ;

	/* Work space for the Fourier transform routines  ***/
	complex fftwk1 ;
	complex fftwk2 ;

	su3_matrix su3fftwk1 ;
	su3_matrix su3fftwk2 ;

} site;


/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

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
EXTERN  char savefile[MAXFILENAME];  /** ending gauge file name ***/
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int saveflag ;

EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN	int startflag_w[MAX_NKAP]; /* beginning wilson: START_ASCII, START_BINARY */
EXTERN	int saveflag_w[MAX_NKAP];	/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int saveflag_m;	/* save meson props: SAVE_ASCII, SAVE_BINARY */
EXTERN  int wall_separation;
EXTERN  int source_parity;
EXTERN  int nchannels;  /* = NCHANNELS +2 if extra_sink chosen */

EXTERN Real *wall_template; /* pointer to array of wall heights */
EXTERN  Real byterevReal;  /* flag to bytereverse lattice */

EXTERN  int fixflag;    /* whether to gauge fix **/


/* global variables required for the static-variational code ***/
EXTERN	char vary_out[MAXFILENAME] ,  smear_meson_out[MAXFILENAME] ; 
EXTERN  int nosmear ;  /* The number of smearing functions **/
EXTERN	char smearfile_in[MAX_SMEAR][MAXFILENAME]; /** the names of the smearing functions ***/

EXTERN	Real smear_code[MAX_SMEAR][5]; /** the code for smearing functions ***/
EXTERN	int total_iters;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN  int	sites_on_node;		/* number of sites on this node */
EXTERN  int	even_sites_on_node;	/* number of even sites on this node */
EXTERN  int	odd_sites_on_node;	/* number of odd sites on this node */
EXTERN  int	number_of_nodes;	/* number of nodes in use */
EXTERN  int  this_node;		/* node number of this node */
EXTERN quark_source wqs;
EXTERN quark_source wqstmp, wqstmp2;  /* Temporary */

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
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

/* Storage for the clover term */
EXTERN clover *gen_clov;

#endif /* _LATTICE_H */
