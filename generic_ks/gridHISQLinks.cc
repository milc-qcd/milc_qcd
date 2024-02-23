// Wrappers for Grid Staggered Link Fattening

#include <omp.h>
#include <Grid/Grid.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"
#include "../include/mGrid/mGrid_assert.h"
#include "../include/macros.h"
#include "../generic/gridMap.h"

extern "C" {
  void dumpmat( su3_matrix *mat );
}

using namespace Grid;

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void
hisqLinks (GRID_info_t *info,
	   double path_coeff[],
	   su3_matrix *fat,
	   su3_matrix *lng,
	   su3_matrix *in,
	   GridCartesian *CGrid)
{
  auto start = std::chrono::system_clock::now();

  // Instantiate the Smear_HISQ class
  Smear_HISQ<Gimpl> HL(CGrid, path_coeff);

  // Copy MILC-formatted thin links
  LatticeGaugeField Umu(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(in, &Umu);

  // Allocate space for output fat and long links
  LatticeGaugeField fatlinks(CGrid);
  GRID_ASSERT(&fatlinks != NULL, GRID_MEM_ERROR);
  LatticeGaugeField lnglinks(CGrid);
  GRID_ASSERT(&lnglinks != NULL, GRID_MEM_ERROR);

  if(lng != NULL){
    HL.smear(fatlinks, lnglinks, Umu);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, &lnglinks);
  }
  else{
    HL.smear(fatlinks, lnglinks, Umu);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "generate fat and long links " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
}

	
template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void
hisqAuxLinks (GRID_info_t *info,
	      double path_coeff[],
	      su3_matrix *U, su3_matrix *V, su3_matrix *W,
	      GridCartesian *CGrid)
{
  auto start = std::chrono::system_clock::now();

  // Load U links
  LatticeGaugeField Ugrid(CGrid);

  // Do the first level fattening
  hisqLinks<LatticeGaugeField, Gimpl, Complex>(info, path_coeff, V, NULL, U, CGrid);

  LatticeGaugeField Vgrid(CGrid);
  LatticeGaugeField Wgrid(CGrid);

  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(V, &Vgrid);

  // Do the reunitarization
  Smear_HISQ<Gimpl> HL(CGrid, path_coeff);
  HL.projectU3(Wgrid, Vgrid);
  
  gridToMilcGaugeField<LatticeGaugeField, Complex>(V, &Vgrid);
  gridToMilcGaugeField<LatticeGaugeField, Complex>(W, &Wgrid);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "generate HISQ aux links " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
}

	
//====================================================================//
// The GRID API for link fattening

void GRID_F3_hisq_links(GRID_info_t *info,
			double path_coeff[],
			su3_matrix *fat,
			su3_matrix *lng,
			su3_matrix *in,
			GRID_4Dgrid *grid_full)
{
  hisqLinks<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, path_coeff, fat, lng, in, grid_full->gridF);
}

void GRID_D3_hisq_links(GRID_info_t *info,
			double path_coeff[],
			su3_matrix *fat,
			su3_matrix *lng,
			su3_matrix *in,
			GRID_4Dgrid *grid_full)
{
  hisqLinks<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, path_coeff, fat, lng, in, grid_full->gridD);
}

void GRID_F3_hisq_aux_links(GRID_info_t *info,
			    double path_coeff[],
			    su3_matrix *U, su3_matrix *V, su3_matrix *W,
			    GRID_4Dgrid *grid_full)
{
  hisqAuxLinks<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, path_coeff, U, V, W, grid_full->gridF);
}

void GRID_D3_hisq_aux_links(GRID_info_t *info,
			    double path_coeff[],
			    su3_matrix *U, su3_matrix *V, su3_matrix *W,
			    GRID_4Dgrid *grid_full)
{
  hisqAuxLinks<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, path_coeff, U, V, W, grid_full->gridD);
}

