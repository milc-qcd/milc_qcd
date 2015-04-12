#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/complex.h"

#define MAXFILENAME  256   /* ASCII string length for all file names */

#define MAXBIN 120  /* Maximum number of correlator bins */
#define RMAX 5      /* Start binning for r > RMAX */
#define MAXRAND 64  /* Maximum number of random sources in a file */
#define NMU 4       /* Number of vector current components per site */
#define MAXFLAV 4   /* Max number of flavors */

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;	/* lattice dimensions */
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
  char job_id[MAXFILENAME]; /* Usually encoded by scripts */

  /*  REPEATING BLOCK */
  int nflav;
  int nrand;           
  char fname[MAXFLAV][MAXFILENAME];
  char corrfile[MAXFILENAME];
  Real charges[MAXFLAV];
  Real mass[MAXFLAV];   /* Filled in later only on node 0 */

}  params;

#endif /* _PARAMS_H */
