#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"
#include "../include/precision.h"

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
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
  int iseed;	/* for random numbers */
  Real beta;      /* gauge coupling */
  int n_dyn_masses;       /* number of dynamical masses */
  Real dyn_mass[MAX_DYN_MASSES];  /* List of dynamical masses */
  int dyn_flavors[MAX_DYN_MASSES]; /* Numbers of dynamical flavors */
  Real u0; /* tadpole parameter */
  /*  REPEATING BLOCK */
  int fixflag;    /* whether to gauge fix */
  Real gauge_fix_tol;
  int rshift[4];
  Real bdry_phase[4];      /* For introducing twisted boundary conditions */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  char gauge_fix_description[MAXFILENAME];
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
