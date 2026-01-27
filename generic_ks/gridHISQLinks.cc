// Wrappers for Grid Staggered Link Fattening

#include <omp.h>

#undef HMC

#include "../include/macros.h"
extern "C" {
#include "../include/fermion_links.h"
}

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"
//#include "../include/mGrid/mGrid_assert.h"

#include "../generic/gridMap.h"
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/HISQSmearing.h>

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

  // Copy MILC-formatted thin links
  LatticeGaugeField Umu(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(in, &Umu);

  // Allocate space for output fat and long links
  LatticeGaugeField fatlinks(CGrid);
  GRID_ASSERT(&fatlinks != NULL);
  LatticeGaugeField lnglinks(CGrid);
  GRID_ASSERT(&lnglinks != NULL);

  // Set action coefficients
  RealD pc_one_link       = path_coeff[0];
  RealD pc_naik           = path_coeff[1];
  RealD pc_three_staple   = path_coeff[2];
  RealD pc_five_staple    = path_coeff[3];
  RealD pc_seven_staple   = path_coeff[4];
  RealD pc_lepage         = path_coeff[5];
  bool backupSVD = false;
  RealD svdTolerance = 0.;
  RealD eigenCutoff = 0.;
  
  HISFContext ctx(pc_one_link, pc_three_staple, pc_five_staple, pc_seven_staple,
		  pc_lepage, pc_naik, backupSVD, svdTolerance, eigenCutoff);

  // Instantiate the HISQ fermion implementation class
  bool calculateStaggeredPhases = true;
  HighlyImprovedStaggeredFermionImpl<Gimpl> HL(CGrid, calculateStaggeredPhases);

  if(lng != NULL){
    HL.smear(fatlinks, lnglinks, Umu, ctx);
    std::cout << "Done with smear" << std::endl << std::flush;
    gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, &lnglinks);
  }
  else{
    HL.smear(fatlinks, lnglinks, Umu, ctx);
    std::cout << "Done with smear" << std::endl << std::flush;
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

  // Set action coefficients
  RealD pc_one_link       = path_coeff[0];
  RealD pc_naik           = path_coeff[1];
  RealD pc_three_staple   = path_coeff[2];
  RealD pc_five_staple    = path_coeff[3];
  RealD pc_seven_staple   = path_coeff[4];
  RealD pc_lepage         = path_coeff[5];
  bool backupSVD = false;
  RealD svdTolerance = 0.;
  RealD eigenCutoff = 0.;
  
  HISFContext ctx(pc_one_link, pc_three_staple, pc_five_staple, pc_seven_staple,
		  pc_lepage, pc_naik, backupSVD, svdTolerance, eigenCutoff);

  // Instantiate the HISQ fermion implementation class
  bool calculateStaggeredPhases = true;
  HighlyImprovedStaggeredFermionImpl<Gimpl> HL(CGrid, calculateStaggeredPhases);

  // Do the reunitarization
  HL.project(Wgrid, Vgrid);
  
  gridToMilcGaugeField<LatticeGaugeField, Complex>(V, &Vgrid);
  gridToMilcGaugeField<LatticeGaugeField, Complex>(W, &Wgrid);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "generate HISQ aux links " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
}

//====================================================================//
// The GRID C API for link fattening

void GRID_F3_hisq_links(GRID_info_t *info,
			double path_coeff[],
			su3_matrix *fat,
			su3_matrix *lng,
			su3_matrix *in,
			GRID_4Dgrid *grid_full)
{
  //  std::cout << "GRID_F3_hisq_links is not supported yet" << std::endl;
  //  assert(0);
  // hisqLinks<LatticeGaugeFieldF, StaggeredImplF, ComplexF>(info, path_coeff, fat, lng, in, grid_full->gridF);
}

void GRID_D3_hisq_links(GRID_info_t *info,
			double path_coeff[],
			su3_matrix *fat,
			su3_matrix *lng,
			su3_matrix *in,
			GRID_4Dgrid *grid_full)
{
  hisqLinks<LatticeGaugeFieldD, StaggeredImplD, ComplexD>(info, path_coeff, fat, lng, in, grid_full->gridD);
}

void GRID_F3_hisq_aux_links(GRID_info_t *info,
			    double path_coeff[],
			    su3_matrix *U, su3_matrix *V, su3_matrix *W,
			    GRID_4Dgrid *grid_full)
{
  std::cout << "GRID_F3_hisq_aux_links" << std::endl;
  assert(0);
  // hisqAuxLinks<LatticeGaugeFieldF, StaggeredImplF, ComplexF>(info, path_coeff, U, V, W, grid_full->gridF);
}

void GRID_D3_hisq_aux_links(GRID_info_t *info,
			    double path_coeff[],
			    su3_matrix *U, su3_matrix *V, su3_matrix *W,
			    GRID_4Dgrid *grid_full)
{
  hisqAuxLinks<LatticeGaugeFieldD, StaggeredImplD, ComplexD>(info, path_coeff, U, V, W, grid_full->gridD);
}

