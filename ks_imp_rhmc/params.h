#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"

#define MAX_MASS_PBP 8
/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt; /* lattice dimensions */
#ifdef FIX_NODE_GEOM
  int node_geometry[4];  /* Specifies fixed "nsquares" (i.e. 4D
			    hypercubes) for the compute nodes in each
			    coordinate direction.  Must be divisors of
			    the lattice dimension */
#ifdef FIX_IONODE_GEOM
  int ionode_geometry[4]; /* Specifies fixed "nsquares" for I/O
			     partitions in each coordinate direction,
			     one I/O node for each square.  The I/O
			     node is at the origin of the square.
			     Must be divisors of the node_geometry. */
#endif
#endif
  int iseed;	          /* for random numbers */
  Real beta;              /* gauge coupling */
  int n_dyn_masses;       /* number of dynamical masses */
  Real dyn_mass[MAX_DYN_MASSES];  /* List of dynamical masses */
  int dyn_flavors[MAX_DYN_MASSES]; /* Numbers of dynamical flavors */
  Real u0;                /* tadpole parameter */
  /*  REPEATING BLOCK */
  int warms;	/* the number of warmup trajectories */
  int trajecs;	/* the number of real trajectories */
  int steps;	/* number of steps for updating */
  int propinterval;     /* number of trajectories between measurements */
  int niter; 	/* maximum number of c.g. iterations */
  int nrestart;   /* maximum number of c.g. restarts */
  int niter_md[MAX_N_PSEUDO], niter_fa[MAX_N_PSEUDO], niter_gr[MAX_N_PSEUDO];
  int prec_md[MAX_N_PSEUDO], prec_fa[MAX_N_PSEUDO], prec_gr[MAX_N_PSEUDO];
  int prec_ff;   /* fermion force precision */
  Real rsqmin_md[MAX_N_PSEUDO], rsqmin_fa[MAX_N_PSEUDO], rsqmin_gr[MAX_N_PSEUDO];
  Real rsqprop;  /* for deciding on convergence */
  Real epsilon;	/* time step */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  int n_pseudo; /* Number of pseudofermion fields */
  char rparamfile[MAXFILENAME];
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  /* PBP and related quantities */
  int num_pbp_masses;   /* Number of masses for pbp calculation */
  quark_invert_control qic_pbp[MAX_MASS_PBP];
  int prec_pbp;         /* Precision of the pbp calculation (1 or 2) */
  int npbp_reps;     /* Number of random sources for pbp calculation */
  ks_param ksp_pbp[MAX_MASS_PBP];
} params;

#endif /* _PARAMS_H */
