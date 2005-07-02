/******************* d_congrad5_two_src.c **********************************/
/* wrapper for two calls to ks_congrad */
/* MIMD version 7 */

#include "generic_ks_includes.h"

int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset src1,    /* source vector (type su3_vector) */
    field_offset src2,
    field_offset dest1,	/* solution vectors */
    field_offset dest2,
    Real mass1,
    Real mass2,
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    Real  *final_rsq_ptr 	/* final residue squared */
    )
{

  int myiters = 0;

  myiters += ks_congrad( src1, dest1, mass1, niter, rsqmin, parity, 
			 final_rsq_ptr );
  myiters += ks_congrad( src2, dest2, mass2, niter, rsqmin, parity, 
			 final_rsq_ptr );

  return myiters;
}
