/****************************** lattice.h ********************************/

/* include file for the MIMD QCD program, version 4
   This file defines global scalars and the fields in the lattice. */

/* Modifications:
   8/10/96  Changed handling of propagator file names; changed same for param C.D.
*/

/* to shrink the size of the code... */
typedef struct { wilson_vector d[4]; } spin_wilson_vector;
typedef struct { spin_wilson_vector c[3]; } wilson_propagator;

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif
#define PI 3.14159265358979323846
#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03
#define VERSION_NUMBER 59354
#define MAX_KAP 6
/* #define SITERAND */	/* Use site-based random number generators */
#define NO_GAUGE_FIX 30    /* constants for gauge fixing */
#define COULOMB_GAUGE_FIX 31

/* Options which control the layout ( for 2d-planes layout )
   With the "ACCORDION" option, the slice numbers are folded, meaning
   that for alternate values of "x" the "y" coordinates run forwards
   and backwards.

   With the "GRAYCODE" option the node numbers are gray coded so that
   adjacent bunches of slices will physically be on adjacent nodes
   in a hypercube architecture

   With the "EVENFIRST" option the even sites are listed contiguously
   in the first part of lattice[], and the odd sites in the last part.
   In this version there must be an equal number of even and odd sites
   in each plane - in other words one of the two shortest dimensions must
   be even.
*/
/*#define ACCORDION */
#ifdef MULTIPLE_MASSES
#define GRAYCODE	/* Use "GRAYCODE" option for Craig's parallel
			   quark propagator i/o routines */
#endif

#define EVENFIRST 

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int niter,wallflag;
EXTERN  char walldescrp[3][30];  /* Describes type of source */
EXTERN	float beta,kappa,width,source_r0,kap[MAX_KAP],resid[MAX_KAP],r0[MAX_KAP];
EXTERN	int num_kap;		/* max number of kappa's <= MAX_KAP */
EXTERN	char startfile[80], savefile[80];
EXTERN  char savefile_w[MAX_KAP][80], startfile_w[MAX_KAP][80];
EXTERN  char scratchstem_w[80];
EXTERN	int startflag;		/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;		/* save lattice: SAVE_ASCII, SAVE_BINARY */
EXTERN	int startflag_w[MAX_KAP];  /* beginning wilson: 
			   RELOAD_ASCII, RELOAD_BINARY, RELOAD_PARALLEL */
EXTERN	int saveflag_w[MAX_KAP];   /* save propagator: SAVE_ASCII, SAVE_BINARY
				   SAVE_PARALLEL */
EXTERN	int total_iters;

/* For Craig McNeile's fast parallel quark propagator io routines */
EXTERN	int layout_flag ;
enum layout_flag_set { HYPER_GRAY_EVENFIRST = 100 } ;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()" */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */


/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* The lattice is an array of sites.  Within each node the even sites will
   be stored first, then the odd sites.
*/


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

	/* wilson complex vectors */
 	wilson_vector psi;	/* solution vector */
 	wilson_vector chi;	/* source vector */
#ifdef BI
	wilson_vector sss;	/* internal biconjugate gradient vector */
	wilson_vector tmp;	/* auxiliary for other internal biCG vectors */
#endif
 /* 	wilson_vector p; conjugate gradient change vector:
				overwrites half the source */
#define IOTMP mp
 	wilson_vector mp;	/* another CG vector */
 /*	wilson_vector r; residue: overwrites half the source */
	/* wilson half vector (temporary used in dslash_w_site) */
	half_wilson_vector htmp[2];

	/* storage for one quark_propagator, for four source spins, three source colors */
	wilson_propagator quark_propagator;
};
typedef struct site site;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING */
EXTERN char ** gen_pt[N_POINTERS];

/* structure for passing simulation parameters to each node */
typedef struct {
	int nx,ny,nz,nt;	/* lattice dimensions */
	int startflag;	/* what to do for beginning lattice */
	int fixflag;    /* whether to gauge fix */
	int saveflag;	/* what to do for saving lattice */
	int startflag_w[MAX_KAP];	/* what to do for beginning wilson vector */
	int saveflag_w[MAX_KAP];	/* what to do for saving wilson vector */
	int wallflag;	/* wall or point source */
	int num_kap;	/* number of kappa's */
	float beta,kappa;	/* gauge coupling, quark hopping parameter */
	float kap[MAX_KAP];	/* kappa values for multiple propagators */
	float resid[MAX_KAP];	/* residue for invertion convergence */
	float r0[MAX_KAP];	/* ``width'' of source */
	int niter; 	/* maximum number of c.g. iterations */
	char startfile[80];
	char startfile_w[MAX_KAP][80];
	char savefile[80];
	char savefile_w[MAX_KAP][80];
	char scratchstem_w[80];
}  params;

/* macros for "field offset" and "field pointer", used when fields
  are arguments to subroutines */
/* Usage:  fo = F_OFFSET( field ), where "field" is the name of a field
  in lattice.
     address = F_PT( &site , fo ), where &site is the address of the
  site and fo is a field_offset.  Usually, the result will have to be
  cast to a pointer to the appropriate type. (It is naturally a char *).
*/
typedef int field_offset;
#define F_OFFSET(a) \
  ((field_offset)(((char *)&(lattice[0]. a ))-((char *)&(lattice[0])) ))
#define F_PT( site , fo )  ((char *)( site ) + (fo)) 

/* macros to loop over sites of a given parity.
   Usage:  
	int i;
	site *s;
	FOREVENSITES(i,s){
	    commands, where s is a pointer to the current site and i is
	    the index of the site on the node
	}
*/

#define FORALLSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)
#ifdef EVENFIRST
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<even_sites_on_node;i++,s++)
#define FORODDSITES(i,s) \
    for(i=even_sites_on_node,s= &(lattice[i]);i<sites_on_node;i++,s++)
#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 ),  \
    s= &(lattice[i]); \
    i< ( (choice)==EVEN ? even_sites_on_node : sites_on_node); \
    i++,s++)
/**#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 ),  \
    s= &(lattice[i]), \
    last_in_loop = ((choice)==EVEN ? even_sites_on_node : sites_on_node); \
    i< last_in_loop; \
    i++,s++)**/
#else
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)if(s->parity==EVEN)
#define FORODDSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)if(s->parity==ODD)
#define FORSOMEPARITY(i,s,choice) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)if( (s->parity & (choice)) != 0)
#endif	/* end ifdef EVENFIRST */

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;
