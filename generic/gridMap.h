#ifndef GRID_MAP_H
#define GRID_MAP_H

#include <Grid/Grid.h>
//#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"

extern "C" {
  void get_coords(int coords[], int node, size_t index);

#include "../include/milc_datatypes.h"
#include "../include/macros.h"

extern	size_t sites_on_node;		/* number of sites on this node */
extern	size_t even_sites_on_node;	/* number of even sites on this node */
extern  int    this_node;
}

using namespace Grid;
//using namespace Grid::QCD;
using namespace std;

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

template<typename LatticeGaugeField>
  struct GRID_ColorMatrix_struct<LatticeGaugeField> *
  create_M( GridCartesian *CGrid );

template<typename LatticeGaugeField>
  static void  
  destroy_M( struct GRID_ColorMatrix_struct<LatticeGaugeField> *mat );

template<typename LatticeGaugeField>
  struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>  *
  asqtad_create_L( GridCartesian *CGrid );

template<typename LatticeGaugeField, typename Complex>
  struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>  *
  asqtad_create_L_from_MILC( su3_matrix *thnlinks, su3_matrix *fatlinks, su3_matrix *lnglinks, 
			     GridCartesian *CGrid );

template<typename LatticeGaugeField, typename Complex>
static void
asqtad_extract_MILC_from_L( su3_matrix *fat, su3_matrix *lng,
			    struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *fn,
			    GridCartesian *CGrid );
  
template<typename FT>
  void  
  asqtad_destroy_L( struct GRID_FermionLinksAsqtad_struct<FT> *gl );

};


using namespace Grid;
//using namespace Grid::QCD;
using namespace std;

extern Coordinate squaresize;

static void
indexToCoords(uint64_t idx, Coordinate &x){

  int r[4];

  // Gets the lattice coordinates from the MILC index
  get_coords(r, this_node, idx);
  // For Grid, we need the coordinates within the sublattice hypercube for the current MPI rank
  // NOTE: Would be better to provide a get_subl_coords() in MILC layout_*.c
  for(int i = 0; i < 4; i++)
    x[i] = r[i] % squaresize[i];

  //printf("Converted %d to %d %d %d %d\n", idx, x[0], x[1], x[2], x[3]); fflush(stdout);
}

// Copy MILC su3_matrix to Grid vLorentzColourMatrix
// Precision conversion can happen here

template<typename sobj, typename Complex>
static void milcSU3MatrixToGrid(su3_matrix *in, sobj &out){

  for (int mu=0; mu<4; mu++)
    for (int i=0; i<Nc; i++)
      for (int j=0; j<Nc; j++)
	out._internal[mu]._internal._internal[i][j] = Complex(in[mu].e[i][j].real, in[mu].e[i][j].imag);
}

// Copy Grid vLorentzColourMatrix to MILC su3_matrix to
// Precision conversion can happen here

template<typename sobj, typename Complex>
static void gridToMilcSU3Matrix(sobj &in, su3_matrix *out){
  for (int mu=0; mu<4; mu++)
    for (int i=0; i<Nc; i++)
      for (int j=0; j<Nc; j++){
	out[mu].e[i][j].real = in._internal[mu]._internal._internal[i][j].real();
	out[mu].e[i][j].imag = in._internal[mu]._internal._internal[i][j].imag();
      }
}

// Map a flattened MILC gauge field (4 matrices per site) to a Grid LatticeGaugeField
template<typename LatticeGaugeField, typename Complex>
void milcGaugeFieldToGrid(su3_matrix *in, LatticeGaugeField *out){

  typedef typename LatticeGaugeField::vector_object vobj;
  typedef typename vobj::scalar_object sobj;

  GridBase *grid = out->Grid();
  int lsites = grid->lSites();
  std::vector<sobj> scalardata(lsites);

  #pragma omp parallel for
    for (size_t milc_idx = 0; milc_idx < sites_on_node; milc_idx++){
      Coordinate x(4);
      indexToCoords(milc_idx, x);
      int grid_idx;
      Coordinate lx(4);
      for (int i = 0; i < 4; i++)lx[i] = x[i];
      Lexicographic::IndexFromCoor(lx, grid_idx, grid->_ldimensions);
      milcSU3MatrixToGrid<sobj, Complex>(in + 4*milc_idx, scalardata[grid_idx]);
    }
  
  vectorizeFromLexOrdArray(scalardata, *out);
}

// Map a Grid LatticeGaugeField  to a flattened MILC gauge field (4 matrices per site)
template<typename LatticeGaugeField, typename Complex>
void gridToMilcGaugeField(su3_matrix *out, LatticeGaugeField *in){

  typedef typename LatticeGaugeField::vector_object vobj;
  typedef typename vobj::scalar_object sobj;

  GridBase *grid = in->Grid();
  int lsites = grid->lSites();
  std::vector<sobj> scalardata(lsites);

  unvectorizeToLexOrdArray(scalardata, *in);

#pragma omp parallel for
    for (size_t milc_idx = 0; milc_idx < sites_on_node; milc_idx++){
      Coordinate x(4);
      indexToCoords(milc_idx, x);
      int grid_idx;
      Coordinate lx(4);
      for (int i = 0; i < 4; i++)lx[i] = x[i];
      Lexicographic::IndexFromCoor(lx, grid_idx, grid->_ldimensions);
      gridToMilcSU3Matrix<sobj, Complex>(scalardata[grid_idx], out + 4*milc_idx);
    }
  
}

#endif // GRID_MAP_H
