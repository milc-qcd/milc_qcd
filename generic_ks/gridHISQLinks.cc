// Wrappers for Grid Staggered Link Fattening

#include <omp.h>
#include <Grid/Grid.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"
#include "../include/mGrid/mGrid_assert.h"
#include "../include/fermion_links.h"
#include "../include/macros.h"
#include "../generic/gridMap.h"

using namespace Grid;

#if 0
// residues and multi_x are indexed by the pseudofermion fields
// multi_x[i] points to a color vector field.
//   The fieldss for each Naik mass are grouped together and
//   the set is concatenated in the order of n_orders_naik.
// n_orders_naik gives the number of pseudofermion fields for each
//   Naik mass.
// deriv[] is indexed by the spacetime dimension

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void
hisqForce (GRID_info_t *info,
	   fermion_links_t *fl,
	   Real residues[],
	   su3_vector *multi_x[],
	   int n_orders_naik[],
	   su3_matrix *deriv[],
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
    std::cout << "Done with smear" << std::endl << std::flush;
    gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, &lnglinks);
  }
  else{
    HL.smear(fatlinks, lnglinks, Umu);
    std::cout << "Done with smear" << std::endl << std::flush;
    gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "generate fat and long links " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
}
#endif
	
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
    std::cout << "Done with smear" << std::endl << std::flush;
    gridToMilcGaugeField<LatticeGaugeField, Complex>(fat, &fatlinks);
    gridToMilcGaugeField<LatticeGaugeField, Complex>(lng, &lnglinks);
  }
  else{
    HL.smear(fatlinks, lnglinks, Umu);
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

/* Q contains the force accumulated thus far and V contains the link
   matrix to be unitarized. Then W contains the unitarized link
   matrix.  Calculate the derivatives of W with respect to respect
   to V: dW/dV and d(W^+)/dV (at fixed V^+ !), where W=V(V^+V)^-1/2
   Return dW = Tr[Q dW/dV) + Tr(Q^+ dW+/dV) */

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void
reunitDeriv(GRID_info_t *info,
	    su3_matrix *V, su3_matrix *dW, su3_matrix *Q, GridCartesian *CGrid)
{
  auto start = std::chrono::system_clock::now();


  // V is the nonunitarized matrix
  LatticeGaugeField Vgrid(CGrid);
  LatticeGaugeField dWgrid(CGrid);
  LatticeGaugeField Qgrid(CGrid);

  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(V, &Vgrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Q, &Qgrid);

  // Calculate the derivative
  Smear_HISQ<Gimpl> RD(CGrid, 0, 0, 0, 0, 0, 0);
  RD.ddVprojectU3(dWgrid, Vgrid, Qgrid, HISQ_FORCE_FILTER);

  gridToMilcGaugeField<LatticeGaugeField, Complex>(dW, &dWgrid);
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Time to do reunit deriv "
	    << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) << " ms"
	    << "\n";
  std::chrono::duration<double, milli>(info->final_sec) = elapsed;
  info->final_sec /= 1e3;

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
  // hisqLinks<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, path_coeff, fat, lng, in, grid_full->gridF);
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
  std::cout << "GRID_F3_hisq_aux_links" << std::endl;
  assert(0);
  // hisqAuxLinks<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, path_coeff, U, V, W, grid_full->gridF);
}

void GRID_D3_hisq_aux_links(GRID_info_t *info,
			    double path_coeff[],
			    su3_matrix *U, su3_matrix *V, su3_matrix *W,
			    GRID_4Dgrid *grid_full)
{
  hisqAuxLinks<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, path_coeff, U, V, W, grid_full->gridD);
}

#if 0
//====================================================================//
// The GRID C API for the fermion force

void GRID_F3_hisq_force(GRID_info_t *info,
			fermion_links_t *fl,
			Real residues[],
			int nterms,
			su3_vector *multi_x[],
			int n_orders_naik[],
			su3_matrix *deriv[],
			GRID_4Dgrid *grid_full)
{
  hisqForce<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, fl, residues,
							  multi_x, n_orders_naik,
							  deriv, grid_full->gridF);
}

void GRID_D3_hisq_force(GRID_info_t *info,
			fermion_links_t *fl,
			Real residues[],
			int nterms,
			su3_vector *multi_x[],
			int n_orders_naik[],
			su3_matrix *deriv[],
			GRID_4Dgrid *grid_full)
{
  hisqForce<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, fl, residues,
							  multi_x, n_orders_naik,
							  deriv, grid_full->gridF);
}

#endif

//====================================================================//
// The GRID C API for testing the reunitarization derivative

void GRID_F3_reunit_deriv( GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			   su3_matrix *Q, GRID_4Dgrid * grid_full ){
  std::cout << "GRID_F3_reunit_deriv is not supported yet" << std::endl;
  assert(0);
  // reunitDeriv<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, V, dW, Q, grid_full->gridF);
}

void GRID_D3_reunit_deriv( GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			   su3_matrix *Q, GRID_4Dgrid * grid_full ){
  reunitDeriv<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, V, dW, Q, grid_full->gridD);
}

