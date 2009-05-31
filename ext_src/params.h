#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"
#include "../include/generic_ks.h" /* For ks_quark_source */
#include "../include/generic_wilson.h"  /* For wilson_quark_source */

#define MAX_QK 6
#define MAX_SINK_LABEL 32
#define CLOVER_TYPE 0
#define KS_TYPE 1

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
  int num_qk;	/* number of quarks */
  int qk_type[MAX_QK];          /* 0 clover 1 KS */
  int startflag_w[MAX_QK];	/* what to do for beginning wilson vector */
  int startflag_ks[MAX_QK];	/* what to do for beginning wilson vector */
  char startfile_w[MAX_QK][MAXFILENAME];
  char startfile_ks[MAX_QK][MAXFILENAME];
  int dst_type[MAX_QK];          /* Extended source type 0 clover 1 KS */
  wilson_quark_source dst_wqs[MAX_QK];
  wilson_quark_source snk_wqs[MAX_QK];
  ks_quark_source dst_ksqs[MAX_QK];
  ks_quark_source snk_ksqs[MAX_QK];
  int snk_gam[MAX_QK]; 
}  params;

#endif /* _PARAMS_H */
