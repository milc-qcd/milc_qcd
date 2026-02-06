// Wrappers for Grid Staggered Link Fattening

#include <omp.h>

#undef HMC

#include "../include/macros.h"

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"
#include "../include/mGrid/mGrid_assert.h"

#include "../generic/gridMap.h"
#include <Grid/Grid.h>
#include <Grid/qcd/action/gauge/PeriodicPlaqPlusRectanglePlusParallelogramAction.h>

using namespace Grid;

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static void gaugeForce (
  GRID_info_t* info,
  su3_matrix* Umilc,
  double beta,
  double u0,
  int nf,
  su3_matrix *deriv,
  GridCartesian* CGrid
) {
  
  auto start = std::chrono::system_clock::now();

  LatticeGaugeField Umu(CGrid), UForce(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Umilc, &Umu);

  // Instantiate the gauge action class
  PeriodicSymanzikOneLoopGaugeAction<Gimpl> action(CGrid, beta, u0, nf);
  
  action.deriv(Umu, UForce);

  gridToMilcGaugeField<LatticeGaugeField, Complex>(deriv, &UForce);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Grid gauge force " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << std::endl;
}

template<typename LatticeGaugeField, typename Gimpl, typename Complex>
static double gaugeAction (
  GRID_info_t* info,
  su3_matrix* Umilc,
  double beta,
  double u0,
  int nf,
  GridCartesian* CGrid
) {
  
  auto start = std::chrono::system_clock::now();

  LatticeGaugeField Umu(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Umilc, &Umu);

  PeriodicSymanzikOneLoopGaugeAction<Gimpl> action(CGrid, beta, u0, nf);
  
  auto actionValue = action.S(Umu);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Grid gauge action " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << std::endl;
  return actionValue;
}

//====================================================================//
// The GRID C API for the gauge force

#if 0
void GRID_F3_gauge_force(GRID_info_t *info,
			 su3_matrix *Umilc,
			 Real beta,
			 Real u0,
			 int nf,
			 su3_matrix *deriv,
			 GRID_4Dgrid *grid_full)
{
  gaugeForce<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, Umilc, beta, u0, nf,
							   deriv, grid_full->gridF);
}
#endif

void GRID_D3_gauge_force(GRID_info_t *info,
			 su3_matrix *Umilc,
			 Real beta,
			 Real u0,
			 int nf,
			 su3_matrix *deriv,
			 GRID_4Dgrid *grid_full)
{
  gaugeForce<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, Umilc, beta, u0, nf,
							   deriv, grid_full->gridD);
}

#if 0
double GRID_F3_gauge_action(GRID_info_t *info,
			    su3_matrix *Umilc,
			    Real beta,
			    Real u0,
			    int nf,
			    GRID_4Dgrid *grid_full)
{
  return gaugeAction<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, Umilc, beta, u0, nf,
								   grid_full->gridF);
}
#endif

double GRID_D3_gauge_action(GRID_info_t *info,
			    su3_matrix *Umilc,
			    Real beta,
			    Real u0,
			    int nf,
			    GRID_4Dgrid *grid_full)
{
  return gaugeAction<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, Umilc, beta, u0, nf,
								   grid_full->gridD);
}


