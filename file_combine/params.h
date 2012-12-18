#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"

#define MAX_FILES 6

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
  int nfile;	                /* number of files */
  int file_type;                /* Type of file */
  int startflag[MAX_FILES];	/* what to do for beginning wilson vector */
  char startfile[MAX_FILES][MAXFILENAME];
  int ncolor, nspin;
  int t0;
  Real coeff[MAX_FILES];
  int saveflag, savetype;
  char savefile[MAXFILENAME];
}  params;

#endif /* _PARAMS_H */
