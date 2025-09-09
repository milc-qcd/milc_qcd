// Mapping between MILC and Grid types

#include "../include/openmp_defs.h"
#include <Grid/Grid.h>
#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/mGrid/mGrid_assert.h"
#include "gridMap.h"

using namespace Grid;
using namespace std;

// Create the color vector interface object

template<typename ImprovedStaggeredFermion>
static struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *
create_V(int milc_parity, GridCartesian *CGrid, GridRedBlackCartesian *RBGrid){
  struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out;

  out = (struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *) 
    malloc(sizeof(struct GRID_ColorVector_struct<ImprovedStaggeredFermion>));
  GRID_ASSERT( out != NULL, GRID_MEM_ERROR );


  switch (milc_parity)
    {
    case EVEN:
      out->cv = new typename ImprovedStaggeredFermion::FermionField(RBGrid);
      GRID_ASSERT(out->cv != NULL, GRID_MEM_ERROR);
      out->cv->Checkerboard() = Even;
      break;

    case ODD:
      out->cv = new typename ImprovedStaggeredFermion::FermionField(RBGrid);
      GRID_ASSERT(out->cv != NULL, GRID_MEM_ERROR);
      out->cv->Checkerboard() = Odd;
      break;

    case EVENANDODD:
      out->cv = new typename ImprovedStaggeredFermion::FermionField(CGrid);
      GRID_ASSERT(out->cv != NULL, GRID_MEM_ERROR);
      break;
      
    default:
      break;
    }

  return out;
}

// Create the block color vector interface object

template<typename ImprovedStaggeredFermion5D>
static struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *
create_nV(int n, int milc_parity,
	  GridCartesian *FCGrid, GridRedBlackCartesian *FRBGrid,
	  GridCartesian *CGrid, GridRedBlackCartesian *RBGrid ){

  struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *out;
  out = (struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *) 
    malloc(sizeof(struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D>));
  GRID_ASSERT( out != NULL, GRID_MEM_ERROR );

  if(milc_parity == EVEN || milc_parity == ODD){
    std::cout << "Constructing 5D field\n" << std::flush;
    out->cv = new typename ImprovedStaggeredFermion5D::FermionField(FRBGrid);
    std::cout << "Done constructing 5D field\n" << std::flush;
    GRID_ASSERT(out->cv != NULL, GRID_MEM_ERROR);
    out->cv->Checkerboard() = milc_parity == EVEN ? Even : Odd ;
  } else {
    out->cv = new typename ImprovedStaggeredFermion5D::FermionField(FCGrid);
    GRID_ASSERT(out->cv != NULL, GRID_MEM_ERROR);
  }

  return out;
}

// free color vector
template<typename ImprovedStaggeredFermion>
static void 
destroy_V( struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *V ){

  if (V->cv != NULL) delete V->cv;
  if (V != NULL) free(V);

  return;
}

// free block color vector
template<typename ImprovedStaggeredFermion5D>
static void 
destroy_nV( struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *V ){

  if (V->cv != NULL) {
    delete V->cv;
  }
  if (V != NULL) free(V);

  return;
}

// Create the color vector interface object
// and Map a MILC color vector field from MILC to Grid layout
// Precision conversion takes place in the copies if need be

template<typename ImprovedStaggeredFermion, typename ColourVector, typename Complex>
static struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *
create_V_from_vec( su3_vector *src, int milc_parity,
		   GridCartesian *CGrid, GridRedBlackCartesian *RBGrid){

  struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out;

  out = create_V<ImprovedStaggeredFermion>(milc_parity, CGrid, RBGrid);
  autoView(Dst_cv, (*(out->cv)), CpuWrite);

  int loopend= (milc_parity)==EVEN ? even_sites_on_node : sites_on_node ;
  int loopstart=((milc_parity)==ODD ? even_sites_on_node : 0 );

  auto start = std::chrono::system_clock::now();
  #pragma omp parallel for 
    for( size_t idx = loopstart; idx < loopend; idx++){

      Coordinate x(4);
      indexToCoords(idx,x);

      ColourVector cVec;
      for(int col=0; col<Nc; col++){
	cVec._internal._internal._internal[col] = 
	  Complex(src[idx].c[col].real, src[idx].c[col].imag);
      }

      pokeLocalSite(cVec, Dst_cv, x);

    }
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  //  std::cout << "Mapped vector field in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
  //	    << "\n" << std::flush;
  
  return out;
}

// Create the blocked color vector interface object
// and Map a set of MILC color vector field from MILC to Grid layout
// Precision conversion takes place in the copies if need be

template<typename ImprovedStaggeredFermion5D, typename ColourVector, typename Complex>
static struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *
create_nV_from_vecs( su3_vector *src[], int n, int milc_parity,
		     GridCartesian *FCGrid, GridRedBlackCartesian *FRBGrid,
		     GridCartesian *CGrid, GridRedBlackCartesian *RBGrid){

  struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *out;

  out = create_nV<ImprovedStaggeredFermion5D>(n, milc_parity, FCGrid, FRBGrid, CGrid, RBGrid);
  autoView(Dst_cv, (*(out->cv)), CpuWrite);

  int loopend= (milc_parity)==EVEN ? even_sites_on_node : sites_on_node ;
  int loopstart=((milc_parity)==ODD ? even_sites_on_node : 0 );


  std::cout << "create_nv_from_vecs: ColourVector size  = " << sizeof(ColourVector)  
	    << " ColourVectorField size = " << sizeof(*(out->cv)) << "\n" << std::flush;
  auto start = std::chrono::system_clock::now();
  #pragma omp parallel for
    for( size_t idx = loopstart; idx < loopend; idx++){

      Coordinate x(4);
      indexToCoords(idx,x);
      Coordinate x5(5);
      for( int d = 0; d < 4; d++ )
	x5[d+1] = x[d];
      
      for( int j = 0; j < n; j++ ){
	x5[0] = j;
	ColourVector cVec;
	for(int col=0; col<Nc; col++){
	  cVec._internal._internal._internal[col] = 
	    Complex(src[j][idx].c[col].real, src[j][idx].c[col].imag);
	}

        pokeLocalSite(cVec, Dst_cv, x5);

      }
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
  autoView(Src_cv, (*(src->cv)), CpuRead);

  FORSOMEFIELDPARITY_OMP(idx, milc_parity, )
    {
      Coordinate x(4);
      indexToCoords(idx, x);
      ColourVector cVec;

      peekLocalSite(cVec, Src_cv, x);

      for(int col = 0; col < Nc; col++)
	{
	  dest[idx].c[col].real = cVec._internal._internal._internal[col].real();
	  dest[idx].c[col].imag = cVec._internal._internal._internal[col].imag();
	}
    } END_LOOP_OMP;

  return;
}

// Map a color vector field from Grid layout to MILC layout
template<typename ImprovedStaggeredFermion5D, typename ColourVector>
static void extract_nV_to_vecs( su3_vector *dest[], int n,
				struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *src, 
				int milc_parity ){
  uint64_t idx;

  autoView(Src_cv, (*(src->cv)), CpuRead);
  FORSOMEFIELDPARITY_OMP(idx, milc_parity, )
    {
      Coordinate x(4);
      indexToCoords(idx, x);
      Coordinate x5(5);
      for( int d = 0; d < 4; d++ )
	x5[d+1] = x[d];

      for( int j = 0; j < n; j++ ){
	x5[0] = j;

	ColourVector cVec;
        peekLocalSite(cVec, Src_cv, x5);
	
	for(int col = 0; col < Nc; col++)
	  {
	    dest[j][idx].c[col].real = cVec._internal._internal._internal[col].real();
	    dest[j][idx].c[col].imag = cVec._internal._internal._internal[col].imag();
	  }
      }
    } END_LOOP_OMP;

  return;
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

// Create and allocate color matrix structure
template<typename LatticeGaugeField>
struct GRID_ColorMatrix_struct<LatticeGaugeField> *
create_M( GridCartesian *CGrid ){
  struct GRID_ColorMatrix_struct<LatticeGaugeField> *out;

  out = ( struct GRID_ColorMatrix_struct<LatticeGaugeField> * )
    malloc(sizeof(struct GRID_ColorMatrix_struct<LatticeGaugeField>));
  GRID_ASSERT(out != NULL, GRID_MEM_ERROR);

  out->links = new LatticeGaugeField(CGrid);
  GRID_ASSERT(out->links != NULL, GRID_MEM_ERROR);

  return out;
}

// delete color matrix structure
template<typename LatticeGaugeField>
static void  
destroy_M( struct GRID_ColorMatrix_struct<LatticeGaugeField> *mat ){

  if (mat == NULL) return;
  
  if (mat->links != NULL) delete mat->links;

  free(mat);
}

// Create and allocate asqtad fermion links structure
template<typename LatticeGaugeField>
struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>  *
asqtad_create_L( GridCartesian *CGrid ){

  struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *out;

  out = (struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *)
    malloc(sizeof(struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>));
  GRID_ASSERT(out != NULL, GRID_MEM_ERROR);

  out->thnlinks = NULL;  // We don't need this one
  out->fatlinks = new LatticeGaugeField(CGrid);
  out->lnglinks = new LatticeGaugeField(CGrid);
  GRID_ASSERT(out->fatlinks != NULL, GRID_MEM_ERROR);
  GRID_ASSERT(out->lnglinks != NULL, GRID_MEM_ERROR);

  return out;
}  


// Create asqtad fermion links object from MILC fields
// Precision conversion takes place in the copies if need be

template<typename LatticeGaugeField, typename Complex>
static struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField>  *
asqtad_create_L_from_MILC( su3_matrix *thn, su3_matrix *fat, 
			   su3_matrix *lng, GridCartesian *CGrid ){

  struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *out =
    asqtad_create_L<LatticeGaugeField>(CGrid);
  
  auto start = std::chrono::system_clock::now();

  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(fat, out->fatlinks);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(lng, out->lnglinks);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Mapped gauge fields in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n" << std::flush;
  return out;
}

// Extract MILC asqtad fermion links from Grid object
// Precision conversion takes place in the copies if need be

template<typename LatticeGaugeField, typename Complex>
static void
asqtad_extract_MILC_from_L( su3_matrix *fat, su3_matrix *lng,
			    struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *fn,
			    GridCartesian *CGrid ){

  auto start = std::chrono::system_clock::now();

  gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, fn->fatlinks);
  gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, fn->lnglinks);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Extracted gauge fields in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n" << std::flush;
}

// free asqtad fermion links
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

// Create a 4D, full-grid wrapper
GRID_4Dgrid *
GRID_create_grid(void){
  const Coordinate latt_size    = GridDefaultLatt();
  const Coordinate simd_layoutF = GridDefaultSimd(Nd,vComplexF::Nsimd());
  const Coordinate simd_layoutD = GridDefaultSimd(Nd,vComplexD::Nsimd());
  const Coordinate mpi_layout   = GridDefaultMpi();

  GridCartesian *CGridF  = new GridCartesian(latt_size,simd_layoutF,mpi_layout);
  GRID_ASSERT(CGridF != NULL, GRID_MEM_ERROR);
  GridCartesian *CGridD  = new GridCartesian(latt_size,simd_layoutD,mpi_layout);
  GRID_ASSERT(CGridD != NULL, GRID_MEM_ERROR);
  GRID_4Dgrid *g4D = (GRID_4Dgrid *)malloc(sizeof(struct GRID_4Dgrid_struct));
  GRID_ASSERT(g4D != NULL, GRID_MEM_ERROR);
  g4D->gridF = CGridF;
  g4D->gridD = CGridD;
  return g4D;
}

// Create a 4D red-black-grid wrapper
GRID_4DRBgrid *
GRID_create_RBgrid(GRID_4Dgrid *grid_full){
  GridRedBlackCartesian *RBGridF = new GridRedBlackCartesian(grid_full->gridF);
  GRID_ASSERT(RBGridF != NULL, GRID_MEM_ERROR);
  GridRedBlackCartesian *RBGridD = new GridRedBlackCartesian(grid_full->gridD);
  GRID_ASSERT(RBGridD != NULL, GRID_MEM_ERROR);
  GRID_4DRBgrid *g4D = (GRID_4DRBgrid *)malloc(sizeof(struct GRID_4DRBgrid_struct));
  GRID_ASSERT(g4D != NULL, GRID_MEM_ERROR);
  g4D->gridF = RBGridF;
  g4D->gridD = RBGridD;
  return g4D;
}

// Create a 5D full-grid wrapper
GRID_5Dgrid *
GRID_create_5Dgrid(int n, GRID_4Dgrid *grid_full){
  GRID_5Dgrid *g5D = (GRID_5Dgrid *)malloc(sizeof(struct GRID_5Dgrid_struct));
  GRID_ASSERT(g5D != NULL, GRID_MEM_ERROR);
  GridCartesian *FGridF = SpaceTimeGrid::makeFiveDimGrid(n, grid_full->gridF);
  GridCartesian *FGridD = SpaceTimeGrid::makeFiveDimGrid(n, grid_full->gridD);
  g5D->gridF = FGridF;
  g5D->gridD = FGridD;
  return g5D;
}

// Create a 5D RB-grid wrapper
GRID_5DRBgrid *
GRID_create_5DRBgrid(int n, GRID_4Dgrid *grid_full){
  GRID_5DRBgrid *g5D = (GRID_5DRBgrid *)malloc(sizeof(struct GRID_5DRBgrid_struct));
  GRID_ASSERT(g5D != NULL, GRID_MEM_ERROR);
  std::cout << "Constructing 5D grid for " << n << " fields\n" << std::flush;
  auto start = std::chrono::system_clock::now();
  GridRedBlackCartesian *FRBGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(n, grid_full->gridF);
  GridRedBlackCartesian *FRBGridD = SpaceTimeGrid::makeFiveDimRedBlackGrid(n, grid_full->gridD);
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Construct 5D grids in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n" << std::flush;
  g5D->gridF = FRBGridF;
  g5D->gridD = FRBGridD;
  return g5D;
}

void 
GRID_destroy_4Dgrid(GRID_4Dgrid *grid){
  if(grid){
    if(grid->gridF)
      delete grid->gridF;
    if(grid->gridD)
      delete grid->gridD;
    free(grid);
  }
}

void 
GRID_destroy_4DRBgrid(GRID_4DRBgrid *grid){
  if(grid){
    if(grid->gridF)
      delete grid->gridF;
    if(grid->gridD)
      delete grid->gridD;
    free(grid);
  }
}

void 
GRID_destroy_5Dgrid(GRID_5Dgrid *grid){
  if(grid){
    if(grid->gridF)
      delete grid->gridF;
    if(grid->gridD)
      delete grid->gridD;
    free(grid);
  }
}

void 
GRID_destroy_5DRBgrid(GRID_5DRBgrid *grid){
  if(grid){
    if(grid->gridF)
      delete grid->gridF;
    if(grid->gridD)
      delete grid->gridD;
    free(grid);
  }
}

// create color vector
GRID_F3_ColorVector *
GRID_F3_create_V( int milc_parity, GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb ){
  return create_V<ImprovedStaggeredFermionF>( milc_parity, grid_full->gridF, grid_rb->gridF );
}

// create block color vector
GRID_F3_ColorVectorBlock *
GRID_F3_create_nV( int n, int milc_parity, 
		   GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb,
		   GRID_4Dgrid *grid_full,GRID_4DRBgrid *grid_rb ){
  return create_nV<ImprovedStaggeredFermion5DF>( n, milc_parity, grid_5D->gridF, grid_5Drb->gridF,
					  grid_full->gridF, grid_rb->gridF );
}

// create color vector
GRID_D3_ColorVector *
GRID_D3_create_V( int milc_parity, GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb ){
  return create_V<ImprovedStaggeredFermionD>( milc_parity, grid_full->gridD, grid_rb->gridD );
}

// ceate block color vector
GRID_D3_ColorVectorBlock *
GRID_D3_create_nV( int n, int milc_parity,
                   GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
                   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb ){
  return create_nV<ImprovedStaggeredFermion5DD>( n, milc_parity, grid_5D->gridD, grid_5Drb->gridD,
					  grid_full->gridD, grid_rb->gridD );
}

// free color vector
void  
GRID_F3_destroy_V( GRID_F3_ColorVector *gcv ){
  destroy_V<ImprovedStaggeredFermionF>( gcv );
}

// free block color vector
void  
GRID_F3_destroy_nV( GRID_F3_ColorVectorBlock *gcv ){
  destroy_nV<ImprovedStaggeredFermion5DF>( gcv );
}

// free color vector
void  
GRID_D3_destroy_V( GRID_D3_ColorVector *gcv ){
  destroy_V<ImprovedStaggeredFermionD>( gcv );
}

// free block color vector
void  
GRID_D3_destroy_nV( GRID_D3_ColorVectorBlock *gcv ){
  destroy_nV<ImprovedStaggeredFermion5DD>( gcv );
}

// Map a Dirac vector field from MILC layout to GRID layout
GRID_F3_ColorVector  *
GRID_F3_create_V_from_vec( su3_vector *src, int milc_parity,
			   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb){
  return create_V_from_vec<ImprovedStaggeredFermionF, ColourVectorF, ComplexF>( src, milc_parity,
							grid_full->gridF, grid_rb->gridF);
}
  
// Map a Dirac vector field from MILC layout to GRID layout
GRID_F3_ColorVectorBlock *
GRID_F3_create_nV_from_vecs( su3_vector *src[], int n, int milc_parity,
			     GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
			     GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb ){
  return create_nV_from_vecs<ImprovedStaggeredFermion5DF, ColourVectorF,
			     ComplexF>( src, n, milc_parity,
			     grid_5D->gridF, grid_5Drb->gridF, grid_full->gridF,grid_rb->gridF );
}
  
// Map a Dirac vector field from MILC layout to QPhiX layout
GRID_D3_ColorVector *
GRID_D3_create_V_from_vec( su3_vector *src, int milc_parity,
			   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb){
  return create_V_from_vec<ImprovedStaggeredFermionD, ColourVectorD, ComplexD>( src, milc_parity,
							grid_full->gridD, grid_rb->gridD);
}
  
// Map a blocked Dirac vector field from MILC layout to QPhiX layout
GRID_D3_ColorVectorBlock *
GRID_D3_create_nV_from_vecs( su3_vector *src[], int n, int milc_parity,
                             GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb,
                             GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb ){
  return create_nV_from_vecs<ImprovedStaggeredFermion5DD, ColourVectorD,
			     ComplexD>( src, n, milc_parity,
					grid_5D->gridD, grid_5Drb->gridD,
					grid_full->gridD, grid_rb->gridD );
}
  
// Map a color vector field from GRID layout to MILC layout
void 
GRID_F3_extract_V_to_vec( su3_vector *dest, GRID_F3_ColorVector *gcv, int milc_parity ){
  extract_V_to_vec<ImprovedStaggeredFermionF, ColourVectorF>( dest, gcv, milc_parity );
}

// Map a block color vector field from GRID layout to MILC layout
void 
GRID_F3_extract_nV_to_vecs( su3_vector *dest[], int n, GRID_F3_ColorVectorBlock *gcv, int milc_parity ){
  extract_nV_to_vecs<ImprovedStaggeredFermion5DF, ColourVectorF>( dest, n, gcv, milc_parity );
}

// Map a color vector field from GRID layout to MILC layout
void 
GRID_D3_extract_V_to_vec( su3_vector *dest, GRID_D3_ColorVector *gcv, int milc_parity ){
  extract_V_to_vec<ImprovedStaggeredFermionD, ColourVectorD>( dest, gcv, milc_parity );
}

// Map a block color vector field from GRID layout to MILC layout
void 
GRID_D3_extract_nV_to_vecs( su3_vector *dest[], int n, GRID_D3_ColorVectorBlock *gcv, int milc_parity ){
  extract_nV_to_vecs<ImprovedStaggeredFermion5DD, ColourVectorD>( dest, n, gcv, milc_parity );
}

// create asqtad fermion links from MILC
GRID_F3_FermionLinksAsqtad  *
GRID_F3_asqtad_create_L_from_MILC( su3_matrix *thnlinks, su3_matrix *fatlinks, su3_matrix *lnglinks, 
				   GRID_4Dgrid *grid_full ){
  return asqtad_create_L_from_MILC<LatticeGaugeFieldF, ComplexF>( thnlinks, 
                  fatlinks, lnglinks, grid_full->gridF );
}

// create asqtad fermion links from MILC
GRID_D3_FermionLinksAsqtad  *
GRID_D3_asqtad_create_L_from_MILC( su3_matrix *thnlinks, su3_matrix *fatlinks, su3_matrix *lnglinks, 
				   GRID_4Dgrid *grid_full ){
  return asqtad_create_L_from_MILC<LatticeGaugeFieldD, ComplexD>( thnlinks, 
		   fatlinks, lnglinks, grid_full->gridD );
}

// extract MILC links from fermion link structure
void
GRID_F3_extract_MILC_from_L( su3_matrix *fatlinks, su3_matrix *lnglinks, GRID_F3_FermionLinksAsqtad  *fn,
			     GRID_4Dgrid *grid_full ){
  asqtad_extract_MILC_from_L<LatticeGaugeFieldF, ComplexF>( fatlinks,
		   lnglinks, fn, grid_full->gridF );
}

// extract MILC links from the asqtad link structure
void
GRID_D3_extract_MILC_from_L( su3_matrix *fatlinks, su3_matrix *lnglinks, GRID_D3_FermionLinksAsqtad  *fn,
			     GRID_4Dgrid *grid_full ){
  asqtad_extract_MILC_from_L<LatticeGaugeFieldD, ComplexD>( fatlinks,
		    lnglinks, fn, grid_full->gridD );
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

// load lattice gauge field from MILC
void
GRID_F3_load_M_from_mat4( su3_matrix *mat4, GRID_F3_ColorMatrix *out )
{
  milcGaugeFieldToGrid<LatticeGaugeFieldF,ComplexF>( mat4, out->links );
}

// load lattice gauge field from MILC
void
GRID_D3_load_M_from_mat4( su3_matrix *mat4, GRID_D3_ColorMatrix *out )
{
  milcGaugeFieldToGrid<LatticeGaugeFieldD,ComplexD>( mat4, out->links );
}

// Create the color vector array interface object
template< typename ImprovedStaggeredFermion >
static struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > *
create_V_array( int n, int milc_parity, GridCartesian * CGrid, GridRedBlackCartesian * RBGrid )
{
  int i; // loop index
  struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > * out;

  out = ( struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > * )
    malloc( sizeof( struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > ) );
  GRID_ASSERT( out != NULL, GRID_MEM_ERROR );

  switch( milc_parity )
  {
  case EVEN:
    out->N = n;
    out->cv = new std::vector< typename ImprovedStaggeredFermion::FermionField >( n, RBGrid );
    GRID_ASSERT( out->cv != NULL, GRID_MEM_ERROR );
    
    for( i=0; i<n; i++ )
    {
      ( ( *(out->cv) )[i] ).Checkerboard() = Even;
    }
    break;

  case ODD:
    out->N = n;
    out->cv = new std::vector< typename ImprovedStaggeredFermion::FermionField >( n, RBGrid );
    GRID_ASSERT( out->cv != NULL, GRID_MEM_ERROR );
    
    for( i=0; i<n; i++ )
    {
      ( ( *(out->cv) )[i] ).Checkerboard() = Odd;
    }
    break;

  case EVENANDODD:
    out->N = n;
    out->cv = new std::vector< typename ImprovedStaggeredFermion::FermionField >( n, CGrid );
    GRID_ASSERT( out->cv != NULL, GRID_MEM_ERROR );
    break;

  default:
    break;
  }

  return out;
}

// Free color vector array
template< typename ImprovedStaggeredFermion >
static void destroy_V_array( struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > * V )
{
  if( V->cv != NULL ) delete V->cv;
  if( V != NULL ) free(V);
}

// Create the color vector array interface object
// and map a MILC color vector field array from MILC to Grid layout
// Precision conversion takes place in the copies if need be
template< typename ImprovedStaggeredFermion, typename ColourVector, typename Complex >
static struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > *
create_V_array_from_vec_array( su3_vector ** src, int n, int milc_parity, GridCartesian * CGrid,
                               GridRedBlackCartesian * RBGrid )
{
  size_t i; // loop index
  
  struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > * out;

  out = create_V_array< ImprovedStaggeredFermion >( n, milc_parity, CGrid, RBGrid );

  int loopend = (milc_parity)==EVEN ? even_sites_on_node : sites_on_node;
  int loopstart = (milc_parity)==ODD ? even_sites_on_node : 0;

  for( i=0; i<n; i++ )
  {
    autoView( Dst_cv,  ( *(out->cv) )[i], CpuWrite);
#pragma omp parallel for
    for( size_t idx=loopstart; idx<loopend; idx++ )
    {
      Coordinate x(4);
      indexToCoords( idx, x );

      ColourVector cVec;

      for( int col=0; col<Nc; col++ )
      {
        cVec._internal._internal._internal[col] =
          Complex( (*(src+i))[idx].c[col].real, (*(src+i))[idx].c[col].imag );
      }
      pokeLocalSite( cVec, Dst_cv, x );      
    }
  }
  
  return out;
}

// Map a color vector field array from Grid layout to MILC layout
template< typename ImprovedStaggeredFermion, typename ColourVector >
static void extract_V_array_to_vec_array(
  su3_vector ** dest, int n,
  struct GRID_ColorVectorArray_struct< ImprovedStaggeredFermion > * src,
  int milc_parity )
{
  size_t i; // loop index
  
  int loopend = (milc_parity)==EVEN ? even_sites_on_node : sites_on_node;
  int loopstart = (milc_parity)==ODD ? even_sites_on_node : 0;

  for( i=0; i<n; i++ )
  {
    autoView( Src_cv, ( *(src->cv) )[i], CpuRead );
#pragma omp parallel for
    for( size_t idx=loopstart; idx<loopend; idx++ )
    {
      Coordinate x(4);
      indexToCoords( idx, x );

      ColourVector cVec;
      peekLocalSite( cVec, Src_cv, x );

      for( int col=0; col<Nc; col++ )
      {
        (*(dest+i))[idx].c[col].real = cVec._internal._internal._internal[col].real();
        (*(dest+i))[idx].c[col].imag = cVec._internal._internal._internal[col].imag();
      }
    }
  }

  return ;
}

// Create color vector array
GRID_F3_ColorVectorArray * GRID_F3_create_V_array(
  int n,
  int milc_parity,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb )
{
  return create_V_array< ImprovedStaggeredFermionF >( n, milc_parity, grid_full->gridF, grid_rb->gridF );
}

// Create color vector array
GRID_D3_ColorVectorArray * GRID_D3_create_V_array(
  int n,
  int milc_parity,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb )
{
  return create_V_array< ImprovedStaggeredFermionD >( n, milc_parity, grid_full->gridD, grid_rb->gridD );
}

// Free color vector array
void GRID_F3_destroy_V_array( GRID_F3_ColorVectorArray * V )
{
  destroy_V_array< ImprovedStaggeredFermionF >( V );
}

// Free color vector array
void GRID_D3_destroy_V_array( GRID_D3_ColorVectorArray * V )
{
  destroy_V_array< ImprovedStaggeredFermionD >( V );
}

// Create color vector array from MILC type
GRID_F3_ColorVectorArray * GRID_F3_create_V_array_from_vec_array(
  su3_vector ** src,
  int n,
  int milc_parity,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb )
{
  return create_V_array_from_vec_array< ImprovedStaggeredFermionF, ColourVectorF, ComplexF >(
    src, n, milc_parity, grid_full->gridF, grid_rb->gridF );
}

// Create color vector array from MILC type
GRID_D3_ColorVectorArray * GRID_D3_create_V_array_from_vec_array(
  su3_vector ** src,
  int n,
  int milc_parity,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb )
{
  return create_V_array_from_vec_array< ImprovedStaggeredFermionD, ColourVectorD, ComplexD >(
    src, n, milc_parity, grid_full->gridD, grid_rb->gridD );
}

// Copy color vector array from Grid structure to MILC type
void GRID_F3_extract_V_array_to_vec_array(
  su3_vector ** dest,
  int n,
  GRID_F3_ColorVectorArray * src,
  int milc_parity )
{
  extract_V_array_to_vec_array< ImprovedStaggeredFermionF, ColourVectorF >( dest, n, src, milc_parity );
  return ;
}

// Copy color vector array from Grid structure to MILC type
void GRID_D3_extract_V_array_to_vec_array(
  su3_vector ** dest,
  int n,
  GRID_D3_ColorVectorArray * src,
  int milc_parity )
{
  extract_V_array_to_vec_array< ImprovedStaggeredFermionD, ColourVectorD >( dest, n, src, milc_parity );
  return ;
}


