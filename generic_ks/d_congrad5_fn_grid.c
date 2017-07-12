/******************* d_congrad5_fn_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*! \brief call to the grid_ks_congrad_parity.
 *
 * Choose the inversion precision
 */
int
ks_congrad_parity_grid ( su3_vector *src
			, su3_vector *sol
			, quark_invert_control *qic
			, Real mass
			, fn_links_t *fn)			 
{
  int iterations_used;
  
  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_parity_grid_F( src, sol, qic, mass, fn );
  else
    iterations_used = 
      ks_congrad_parity_grid_D( src, sol, qic, mass, fn );
  
  total_iters += iterations_used;
  return iterations_used;
}
