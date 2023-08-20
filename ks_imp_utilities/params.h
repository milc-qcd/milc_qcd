#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/imp_ferm_links.h"
#include <stdint.h>

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
  uint32_t iseed;	/* for random numbers */
  /*  REPEATING BLOCK */
  Real u0; /* tadpole parameter */
  int coord_origin[4];  /* Origin of coordinates for KS phases and time_bc */
  int time_bc;          /* 0 for antiperiodic, 1 for periodic */
  int nmass;    /* number of masses */
  ks_param ksp[MAX_MASS];
  quark_invert_control qic[MAX_MASS];
  int niter; 	/* maximum number of c.g. iterations */
  int nrestart;	/* maximum number of c.g. restarts */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  int savelongflag, savefatflag;  /* same for longlinks and fatlinks */
  int withKSphases;  /* T/F include KS phases in output fat/long links */
  int srcflag[MAX_MASS]; /* what to do for source lattice */
  int ansflag[MAX_MASS]; /* what to do for answer lattice */
  
  char ks_eigen_startfile[MAXFILENAME]; /* KS eigenvector file to be loaded */
  char ks_eigen_savefile[MAXFILENAME]; /* KS eigenvector file to be saved */
  int ks_eigen_startflag; /* what to do for beginning eigenvectors */
  int ks_eigen_saveflag; /* what to do for ending eigenvectors */
  ks_eigen_param eigen_param; /* Parameters for eigensolver */

  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char savelongfile[MAXFILENAME],savefatfile[MAXFILENAME];
  char stringLFNlong[MAXFILENAME],stringLFNfat[MAXFILENAME];  /** ILDG LFN if applicable ***/
  char srcfile[MAX_MASS][MAXFILENAME],ansfile[MAX_MASS][MAXFILENAME];
  int inverttype;
}  params;

#endif /* _PARAMS_H */
