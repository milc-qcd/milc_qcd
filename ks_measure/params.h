#ifndef _PARAMS_H
#define _PARAMS_H

#include "defines.h"
#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h"
#include "../include/generic_wilson.h"
#include "../include/gammatypes.h"

#define MAX_MASS_PBP 8
#define MAX_SET 8
#define MAX_PBP_MASSES 32

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
  int iseed;
  char job_id[MAXFILENAME]; /* Usually encoded by scripts */

  /*  REPEATING BLOCK */
  int startflag;	/* what to do for beginning lattice */
  Real u0;
  Real staple_weight;
  int ape_iter;
  int saveflag;	/* what to do for saving lattice */
  char startfile[MAXFILENAME];  /* Gauge file */
  char savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  int num_set;                  /* Number of sets */
  int num_pbp_masses[MAX_SET];   /* Number of masses for pbp calculation */
  int begin_pbp_masses[MAX_SET]; /* index of beginning propagator in this set */
  int end_pbp_masses[MAX_SET]; /* index of ending propagator in this set */
  int prec_pbp[MAX_SET];         /* Precision of the pbp calculation (1 or 2) */
  int npbp_reps[MAX_SET];     /* Number of random sources for pbp calculation */
  quark_invert_control qic_pbp[MAX_MASS_PBP];
  ks_param ksp_pbp[MAX_MASS_PBP];
  int set[MAX_MASS_PBP];  /* The set to which the propagator belongs */
  char mass_label[MAX_MASS_PBP][32];
}  params;


#endif /* _PARAMS_H */
