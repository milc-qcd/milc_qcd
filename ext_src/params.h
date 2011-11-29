#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"

#define MAX_QK 6
#define MAX_T0 32
#define MAX_SINK_LABEL 32
#define CLOVER_TYPE 0
#define KS_TYPE 1
#define KS4_TYPE 2

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
  int iseed;
  int num_qk;	/* number of quarks */
  int qk_type[MAX_QK];          /* 0 clover 1 KS */
  int startflag_w[MAX_QK];	/* what to do for beginning wilson vector */
  int startflag_ks[MAX_QK];	/* what to do for beginning wilson vector */
  char startfile_w[MAX_QK][MAXFILENAME];
  char startfile_ks[MAX_QK][MAXFILENAME];
  int ncolor[MAX_QK];
  int dst_type[MAX_QK];          /* Extended source type 0 clover 1 KS */
  int r_offset[MAX_QK][4];       /* Shift of origin for meson correlator */
  int num_t0[MAX_QK];            /* Number of time slices for each quark */
  quark_source dst_qs[MAX_QK][MAX_T0];
  quark_source_sink_op snk_qs_op[MAX_QK];
  int snk_gam[MAX_QK]; 
}  params;

#endif /* _PARAMS_H */
