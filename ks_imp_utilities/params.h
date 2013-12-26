#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"
#include "../include/generic_quark_types.h"

#define MAX_MASS 8

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
  /*  REPEATING BLOCK */
  Real u0; /* tadpole parameter */
  int nmass;    /* number of masses */
  ks_param ksp[MAX_MASS];
  quark_invert_control qic[MAX_MASS];
  int niter; 	/* maximum number of c.g. iterations */
  int nrestart;	/* maximum number of c.g. restarts */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  int savelongflag, savefatflag;  /* same for longlinks and fatlinks */
  int srcflag[MAX_MASS]; /* what to do for source lattice */
  int ansflag[MAX_MASS]; /* what to do for answer lattice */
  
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char savelongfile[MAXFILENAME],savefatfile[MAXFILENAME];
  char srcfile[MAX_MASS][MAXFILENAME],ansfile[MAX_MASS][MAXFILENAME];
  int inverttype;
}  params;

#endif /* _PARAMS_H */
