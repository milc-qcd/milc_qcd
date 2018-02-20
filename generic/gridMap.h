#ifndef GRID_MAP_H
#define GRID_MAP_H

#include <Grid/Grid.h>
//#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"

extern "C" {
#include "generic_includes.h"
#include "../include/openmp_defs.h"
}

class GridMap {

 public:

  GridMap( void ){};
	  
template<typename FT>
  struct GRID_ColorVector_struct<FT> *
  create_V_from_vec( su3_vector *src, GRID_evenodd_t evenodd);

template<typename FT>
  void 
  extract_V_to_vec( su3_vector *dest, struct GRID_ColorVector_struct<FT> *gff, 
		     GRID_evenodd_t evenodd );

template<typename FT>
  void 
  destroy_V( struct GRID_ColorVector_struct<FT> *gff );

template<typename FT>
  struct GRID_FermionLinksAsqtad_struct<FT>  *
  asqtad_create_L_from_MILC( su3_matrix *thnlinks, su3_matrix *fatlinks, su3_matrix *lnglinks, 
			     GRID_evenodd_t evenodd );

template<typename FT>
  void  
  asqtad_destroy_L( struct GRID_FermionLinksAsqtad_struct<FT> *gl );

};

#endif // GRID_MAP_H
