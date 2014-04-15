#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */

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
        int iseed;      /* for random numbers */
   /*  REPEATING BLOCK */
#ifdef HYP
        Real alpha;
        Real alpha2;     /* weight parameters in HYP APE blocking */
        Real alpha3;
#else
        Real ape_weight;     /* weight parameter in APE blocking */
#endif
        int sweeps;     /* the number of real trajectories */
        int hits;       /* number of SU(2) hits for SU(3) projection */
        int measinterval;     /* smoothing steps between measurements */
        int startflag;  /* what to do for beginning lattice */
        int fixflag;    /* whether to gauge fix */
        int saveflag;   /* what to do with lattice at end */
        int savetopoflag;  /* what to do with topo file at end */
        char startfile[MAXFILENAME],savefile[MAXFILENAME],topofile[MAXFILENAME];
	char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
