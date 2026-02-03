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
#include <Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>

using namespace Grid;

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void hisqLinks(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* fat,
  su3_matrix* lng,
  su3_matrix* in,
  GridCartesian* CGrid,
  bool reunitarize = false,
  bool filter = false
) {
  // start timer
  auto start = std::chrono::system_clock::now();

  // Copy MILC-formatted thin links
  LatticeGaugeField Umu(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(in, &Umu);

  // Allocate space for output fat and long links
  LatticeGaugeField fatlinks(CGrid);
  GRID_ASSERT(&fatlinks != NULL);
  LatticeGaugeField lnglinks(CGrid);
  GRID_ASSERT(&lnglinks != NULL);

  // reunitarization "force filter" & backup SVD
  Real eigenvalue_cutoff = (filter) ? HISQ_FORCE_FILTER : 0.;
  Real rel_svd_tol = HISQ_REUNIT_SVD_REL_ERROR;
  Real abs_svd_tol = HISQ_REUNIT_SVD_ABS_ERROR;
  bool allow_svd = false;
  bool svd_only = false;
  
#ifdef HISQ_REUNIT_ALLOW_SVD
  allow_svd = true;
#endif

#ifdef HISQ_REUNIT_SVD_ONLY
  svd_only = true;
  rel_svd_tol = 0.;
#endif

  // Instantiate context object
  HISFContext ctx(
    path_coeff[0],         // 1-link
    path_coeff[2],         // 3-link
    path_coeff[3],         // 5-link
    path_coeff[4],         // 7-link
    path_coeff[5],         // Lepage
    path_coeff[1],         // Naik
    allow_svd || svd_only, // backup SVD
    rel_svd_tol,           // relative SVD tolerance
    abs_svd_tol,           // absolute SVD tolerance
    eigenvalue_cutoff      // reunit eig cutoff = "force filter"
  );

  // Instantiate the HISQ fermion implementation class
  bool calculateStaggeredPhases = false;
  HighlyImprovedStaggeredFermionImpl<Gimpl> hisq(CGrid, calculateStaggeredPhases);

  // Smear according to context
  if (lng != NULL) {
    hisq.smear(fatlinks, lnglinks, Umu, ctx);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, &lnglinks);
  } else {
    hisq.smear(fatlinks, Umu, ctx);
  }
  std::cout << "Done with smear" << std::endl << std::flush;
  if (reunitarize) hisq.project(fatlinks, fatlinks, ctx); // reunitarize
  gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);

  // finish timing
  if (reunitarize) {
    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    std::cout << "generate reunitarized fat links "
        << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
        << std::endl;
  } else {
    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    std::cout << "generate hisq and long links "
        << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
        << std::endl;
  }
}

	
template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void hisqAuxLinks(
  GRID_info_t* info,
  double path_coeff[],
  su3_matrix* U,
  su3_matrix* V,
  su3_matrix* W,
  GridCartesian* CGrid,
  bool filter = false
) {
  // Do the first level fattening w/ additional reunitarization
  // NOTE: The auxiliary field V is not returned.
  hisqLinks<LatticeGaugeField, Gimpl, Complex>(info, path_coeff, W, NULL, U, CGrid, true, filter);
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
  std::cout << "GRID_F3_hisq_links is not supported yet" << std::endl;
  assert(0);
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

