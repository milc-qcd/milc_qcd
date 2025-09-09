#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/imp_ferm_links.h"
#include <stdint.h>

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
  uint32_t iseed;	/* for random numbers */
  char job_id[MAXFILENAME]; /* Usually encoded by scripts */

  /*  REPEATING BLOCK */
  Real mass; /* gauge coupling, quark mass */
  Real u0; /* tadpole parameter */
  int coord_origin[4];  /* Origin of coordinates for KS phases and time_bc */
  int time_bc;          /* 0 for antiperiodic, 1 for periodic */
  int fixflag;    /* whether to gauge fix */
  Real charge;     /* charge for the Dirac matrix */
  char charge_label[32];  /* for eigenpair label */
  ks_eigen_param eigen_param; /* Parameters for eigensolver */
  int ks_eigen_startflag; /* what to do for beginning eigenvectors */
  char ks_eigen_startfile[MAXFILENAME]; /* KS eigenvector file to be loaded */
  int ks_eigen_saveflag; /* eigenvector file type */
  char ks_eigen_savefile[MAXFILENAME]; /* eigenvector output file name */
  int start_u1flag;	/* what to do for beginning u(1) lattice */
  int save_u1flag;	/* what to do with ending u(1) lattice */
  char start_u1file[MAXFILENAME]; /* U(1) gauge file */
  char save_u1file[MAXFILENAME]; /* U(1) gauge file */
  int startflag;  /* what to do for beginning lattice */
  char startfile[MAXFILENAME]; /* starting lattice file name */
  int saveflag;	/* what to do for saving lattice */
  char savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
}  params;

#endif /* _PARAMS_H */
