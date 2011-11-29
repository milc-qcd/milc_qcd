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
    int nrestart,       /* maximum number of restarts */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for the inversion */
    int parity,		/* parity to be worked on */
    Real  *final_rsq_ptr, 	/* final residue squared */
    imp_ferm_links_t *fn       /* Storage for fermion links */
    )
{

  int myiters = 0;

  myiters += ks_congrad( src1, dest1, mass1, niter, nrestart, rsqmin, prec,
			 parity, final_rsq_ptr, fn );
  myiters += ks_congrad( src2, dest2, mass2, niter, nrestart, rsqmin, prec, 
			 parity, final_rsq_ptr, fn );

  return myiters;
}
