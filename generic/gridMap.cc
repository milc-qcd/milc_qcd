// Mapping between MILC and Grid types

#include <Grid/Grid.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/mGrid/mGrid_assert.h"

extern "C" {
#include "generic_includes.h"
#include "../include/openmp_defs.h"
}


using namespace Grid;
using namespace Grid::QCD;
using namespace std;

extern GridCartesian         *CGrid;
extern GridRedBlackCartesian *RBGrid;
extern std::vector<int> squaresize;

static void
indexToCoords(uint64_t idx, std::vector<int> &x){

  // Gets the lattice coordinates from the MILC index
  get_coords(x.data(), this_node, idx);
  // For Grid, we need the coordinates within the sublattice hypercube for the current MPI rank
  // NOTE: Would be better to provide a get_subl_coords() in MILC layout_*.c
  for(int i = 0; i < 4; i++)
    x[i] %= squaresize[i];

  //printf("Converted %d to %d %d %d %d\n", idx, x[0], x[1], x[2], x[3]); fflush(stdout);
}


// Create the color vector interface object
// and Map a MILC color vector field from MILC to Grid layout
// Precision conversion takes place in the copies if need be

template<typename ImprovedStaggeredFermion, typename ColourVector, typename Complex>
static struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *
create_V_from_vec( su3_vector *src, int milc_parity){

  node0_printf("Entered create_V_from_vec\n"); fflush(stdout);

  struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out;

  out = (struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *) 
    malloc(sizeof(struct GRID_ColorVector_struct<ImprovedStaggeredFermion>));
  GRID_ASSERT( out != NULL, GRID_MEM_ERROR );
  //  GRID_ASSERT( milc_parity != EVENANDODD, GRID_FAIL );  // We don't support EVENANDODD

  switch (milc_parity)
    {
    case EVEN:
      out->cv = new typename ImprovedStaggeredFermion::FermionField(RBGrid);
      out->cv->checkerboard = Even;
      break;

    case ODD:
      out->cv = new typename ImprovedStaggeredFermion::FermionField(RBGrid);
      out->cv->checkerboard = Odd;
      break;

    case EVENANDODD:
      out->cv = new typename ImprovedStaggeredFermion::FermionField(CGrid);
      break;
      
    default:
      break;
    }

  GRID_ASSERT(out->cv != NULL, GRID_FAIL);

  int loopend= (milc_parity)==EVEN ? even_sites_on_node : sites_on_node ;
  int loopstart=((milc_parity)==ODD ? even_sites_on_node : 0 );

  auto start = std::chrono::system_clock::now();
  PARALLEL_FOR_LOOP
    for( uint64_t idx = loopstart; idx < loopend; idx++){

      std::vector<int> x(4);
      indexToCoords(idx,x);

      ColourVector cVec;
      for(int col=0; col<Nc; col++){
	cVec._internal._internal._internal[col] = 
	  Complex(src[idx].c[col].real, src[idx].c[col].imag);
      }
      
      pokeLocalSite(cVec, *(out->cv), x);
      
    }
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Mapped vector field in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n" << std::flush;
  
  return out;
}

// Map a color vector field from Grid layout to MILC layout
template<typename ImprovedStaggeredFermion, typename ColourVector>
static void extract_V_to_vec( su3_vector *dest, 
			      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *src, 
			      int milc_parity ){
  uint64_t idx;

  FORSOMEFIELDPARITY_OMP(idx, milc_parity, )
    {
      std::vector<int> x(4);
      indexToCoords(idx, x);
      ColourVector cVec;

      peekLocalSite(cVec, *(src->cv), x);

      for(int col = 0; col < Nc; col++)
	{
	  dest[idx].c[col].real = cVec._internal._internal._internal[col].real();
	  dest[idx].c[col].imag = cVec._internal._internal._internal[col].imag();
	}
    } END_LOOP_OMP;

  return;
}

// free color vector
template<typename ImprovedStaggeredFermion>
static void 
destroy_V( struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *V ){

  if (V->cv != NULL) delete V->cv;
  if (V != NULL) free(V);

  return;
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

// Map a flattened MILC gauge field (4 matrices per site) to a Grid LatticeGaugeField
template<typename LatticeGaugeField, typename Complex>
static void milcGaugeFieldToGrid(su3_matrix *in, LatticeGaugeField &out, int milc_parity){

  typedef typename LatticeGaugeField::vector_object vobj;
  typedef typename vobj::scalar_object sobj;

  GridBase *grid = out._grid;
  int lsites = grid->lSites();
  std::vector<sobj> scalardata(lsites);

  int loopend= (milc_parity)==EVEN ? even_sites_on_node : sites_on_node ;
  int loopstart=((milc_parity)==ODD ? even_sites_on_node : 0 );

  PARALLEL_FOR_LOOP
    for (uint64_t milc_idx = loopstart; milc_idx < loopend; milc_idx++){
      std::vector<int> x(4);
      indexToCoords(milc_idx, x);
      int grid_idx;
      Lexicographic::IndexFromCoor(x, grid_idx, grid->LocalDimensions());
      milcSU3MatrixToGrid<sobj, Complex>(in + 4*milc_idx, scalardata[grid_idx]);
    }
  
  vectorizeFromLexOrdArray(scalardata, out);
}

template<typename ColourMatrix>
static void dumpGrid(ColourMatrix out){

  for (int i=0; i<Nc; i++){
    for (int j=0; j<Nc; j++)
      std::cout << "(" << out._internal._internal._internal[i][j].real() << "," << 
	out._internal._internal._internal[i][j].imag() << ") ";
    std::cout << "\n";
  }
  std::cout << "\n";
}

// Create asqtad fermion links object from MILC fields
// Precision conversion takes place in the copies if need be

template<typename LatticeGaugeField, typename LatticeColourMatrix, typename Complex>
static struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>  *
asqtad_create_L_from_MILC( su3_matrix *thn, su3_matrix *fat, 
			   su3_matrix *lng, int milc_parity ){

  // We support only even and odd.
  GRID_ASSERT(milc_parity == EVENANDODD, GRID_MEM_ERROR);

  struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *out;

  out = (struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *) malloc(sizeof(struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>));
  GRID_ASSERT(out != NULL, GRID_MEM_ERROR);


  out->thnlinks = NULL;  // We don't need this one
  out->fatlinks = new LatticeGaugeField(CGrid);
  out->lnglinks = new LatticeGaugeField(CGrid);
  
  auto start = std::chrono::system_clock::now();

  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(fat, *out->fatlinks, milc_parity);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(lng, *out->lnglinks, milc_parity);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Mapped gauge fields in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n" << std::flush;
  return out;
}

  // free aasqtad fermion links
template<typename LatticeGaugeField>
static void  
asqtad_destroy_L( struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *Link ){

  if (Link == NULL) return;
  
  if (Link->thnlinks != NULL) delete Link->thnlinks;
  if (Link->fatlinks != NULL) delete Link->fatlinks;
  if (Link->lnglinks != NULL) delete Link->lnglinks;
  
  delete Link;
  
  Link = NULL;
}

//====================================================================//
// The GRID C API for mapping between MILC and GRID types

// Map a Dirac vector field from MILC layout to GRID layout
GRID_F3_ColorVector *
GRID_F3_create_V_from_vec( su3_vector *src, int milc_parity ){
  return create_V_from_vec<ImprovedStaggeredFermionF, ColourVectorF, ComplexF>( src, milc_parity );
}
  
// Map a Dirac vector field from MILC layout to QPhiX layout
GRID_D3_ColorVector *
GRID_D3_create_V_from_vec( su3_vector *src, int milc_parity ){
  return create_V_from_vec<ImprovedStaggeredFermionD, ColourVectorD, ComplexD>( src, milc_parity );
}
  
// Map a color vector field from GRID layout to MILC layout
void 
GRID_F3_extract_V_to_vec( su3_vector *dest, GRID_F3_ColorVector *gcv, int milc_parity ){
  extract_V_to_vec<ImprovedStaggeredFermionF, ColourVectorF>( dest, gcv, milc_parity );
}

// Map a color vector field from GRID layout to MILC layout
void 
GRID_D3_extract_V_to_vec( su3_vector *dest, GRID_D3_ColorVector *gcv, int milc_parity ){
  extract_V_to_vec<ImprovedStaggeredFermionD, ColourVectorD>( dest, gcv, milc_parity );
}

// free color vector
void  
GRID_F3_destroy_V( GRID_F3_ColorVector *gcv ){
  destroy_V<ImprovedStaggeredFermionF>( gcv );
}

// free Dirac vector
void  
GRID_D3_destroy_V( GRID_D3_ColorVector *gcv ){
  destroy_V<ImprovedStaggeredFermionD>( gcv );
}

// create asqtad fermion links from MILC
GRID_F3_FermionLinksAsqtad  *
GRID_F3_asqtad_create_L_from_MILC( su3_matrix *thnlinks, su3_matrix *fatlinks, su3_matrix *lnglinks, 
				   int milc_parity ){
  return asqtad_create_L_from_MILC<LatticeGaugeFieldF, LatticeColourMatrixF, ComplexF>( thnlinks, fatlinks, lnglinks, milc_parity );
}

// create asqtad fermion links from MILC
GRID_D3_FermionLinksAsqtad  *
GRID_D3_asqtad_create_L_from_MILC( su3_matrix *thnlinks, su3_matrix *fatlinks, su3_matrix *lnglinks, 
				   int milc_parity ){
  return asqtad_create_L_from_MILC<LatticeGaugeFieldD, LatticeColourMatrixD, ComplexD>( thnlinks, fatlinks, lnglinks, milc_parity );
}

// free asqtad fermion links
void  
GRID_F3_asqtad_destroy_L( GRID_F3_FermionLinksAsqtad *gl ){
  asqtad_destroy_L<LatticeGaugeFieldF>( gl );
}

// free asqtad fermion links
void  
GRID_D3_asqtad_destroy_L( GRID_D3_FermionLinksAsqtad *gl ){
  asqtad_destroy_L<LatticeGaugeFieldD>( gl );
}


