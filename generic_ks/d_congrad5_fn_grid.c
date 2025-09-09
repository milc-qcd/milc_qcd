/******************* d_congrad5_fn_grid.c ************************/
/* For the Grid interface */
/* MIMD version 7 */

/* TODO: standardize API with at least d_congrad5_fn_quda.c and possibly d_congrad5_fn_milc.c */

#include "../include/generic_grid.h"
#include "../include/generic_ks_grid.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*! \brief call to ks_congrad_parity_grid.
 *
 * Choose the inversion precision
 */
int
ks_congrad_parity_gpu ( su3_vector *src
			, su3_vector *sol
			, quark_invert_control *qic
			, Real mass
			, fn_links_t *fn)			 
{
  int iterations_used;
  
  if(qic->prec == 1){

#if defined(HALF_MIXED) || defined(MAX_MIXED)
    node0_printf("ERROR: Mixed precision in GRID is supported only for double precision\n");
    node0_printf("Solving unmixed in single precision");
#endif
    iterations_used = 
      ks_congrad_parity_grid_F( src, sol, qic, mass, fn );

  } else {

#if defined(HALF_MIXED) || defined(MAX_MIXED)
    node0_printf("Using ks_congrad_mixed_parity_grid_D\n");
    iterations_used = 
      ks_congrad_mixed_parity_grid_D( src, sol, qic, mass, fn );
#else    
    iterations_used = 
      ks_congrad_parity_grid_D( src, sol, qic, mass, fn );
#endif

  }
  
  return iterations_used;
}

/*! \brief call to ks_congrad_block_parity_grid.
 *
 * Choose the inversion precision
 */
int
ks_congrad_block_parity_gpu ( int nrhs
			      , su3_vector *src[]
			      , su3_vector *sol[]
			      , quark_invert_control *qic
			      , Real mass
			      , fn_links_t *fn)			 
{
  int iterations_used;
  
  if(qic->prec == 1){

#if defined(HALF_MIXED) || defined(MAX_MIXED)
    node0_printf("ERROR: Mixed precision in GRID is supported only for double precision\n");
    node0_printf("Solving unmixed in single precision");
#endif
    iterations_used = 
      ks_congrad_block_parity_grid_F( nrhs, src, sol, qic, mass, fn );
  } else {

#if defined(HALF_MIXED) || defined(MAX_MIXED)
    node0_printf("Using ks_congrad_block_mixed_parity_grid_D\n");
    iterations_used = 
      ks_congrad_mixed_block_parity_grid_D( nrhs, src, sol, qic, mass, fn );
#else
    iterations_used = 
      ks_congrad_block_parity_grid_D( nrhs, src, sol, qic, mass, fn );
#endif
  }
  
  return iterations_used;
}
