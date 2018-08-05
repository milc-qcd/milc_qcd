/******************* d_congrad5_fn_qphix.c ************************/
/* For the QPhiX interface */
/* MIMD version 7 */

#include "../include/generic_qphix.h"
#include "../include/generic_ks_qphix.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*! \brief call to the qphix_ks_congrad_parity.
 *
 * Choose the inversion precision
 */
int
ks_congrad_parity_qphix ( su3_vector *src
			, su3_vector *sol
			, quark_invert_control *qic
			, Real mass
			, fn_links_t *fn)			 
{
  int iterations_used;
  
  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_parity_qphix_F( src, sol, qic, mass, fn );
  else
    iterations_used = 
      ks_congrad_parity_qphix_D( src, sol, qic, mass, fn );
  
  total_iters += iterations_used;
  return iterations_used;
}

int ks_congrad_block_parity_qphix( int nsrc, su3_vector **t_src, su3_vector **t_dest, 
				   quark_invert_control *qic, Real mass,
				   imp_ferm_links_t *fn){
  /* FAKE version for now */
  int iters = 0;
  for(int i = 0; i < nsrc; i++)
    iters += ks_congrad_parity_qphix(t_src[i], t_dest[i], qic, mass, fn);
  return iters;
}


