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
  double cp,
  double cr,
  double cpg,
  su3_matrix *deriv,
  GridCartesian* CGrid
) {
  
  auto start = std::chrono::system_clock::now();

  LatticeGaugeField Umu(CGrid), UForce(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Umilc, &Umu);

  // Factors of 2*Nc to match MILC conventions */
  OneLoopGaugeActionContext ctx(beta, 2.*Nc*cp, 2.*Nc*cr, 2.*Nc*cpg);

  // Instantiate the gauge action class
  PeriodicPlaqPlusRectanglePlusParallelogramGaugeAction<Gimpl> action(CGrid, ctx);
  
  action.deriv(Umu, UForce, ctx);

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
  double cp,
  double cr,
  double cpg,
  GridCartesian* CGrid
) {
  
  auto start = std::chrono::system_clock::now();

  LatticeGaugeField Umu(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Umilc, &Umu);

  OneLoopGaugeActionContext ctx(beta, cp, cr, cpg);

  PeriodicPlaqPlusRectanglePlusParallelogramGaugeAction<Gimpl> action(CGrid, ctx);
  
  double actionValue = action.S(Umu);

  // A MILC convention
  actionValue *= 3.;

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
			 float beta,
			 double cp,
			 double cr,
			 double cpg,
			 su3_matrix *deriv,
			 GRID_4Dgrid *grid_full)
{
  gaugeForcePhase<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, Umilc, beta, cp, cr, cpg,
							   deriv, grid_full->gridF);
}
#endif

void GRID_D3_gauge_force(GRID_info_t *info,
			 su3_matrix *Umilc,
			 double beta,
			 double cp,
			 double cr,
			 double cpg,
			 su3_matrix *deriv,
			 GRID_4Dgrid *grid_full)
{
  gaugeForce<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, Umilc, beta, cp, cr, cpg,
							   deriv, grid_full->gridD);
}

#if 0
double GRID_F3_gauge_action(GRID_info_t *info,
			    su3_matrix *Umilc,
			    float beta,
			    double cp,
			    double cr,
			    double cpg,
			    GRID_4Dgrid *grid_full)
{
  return gaugeAction<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, Umilc, beta, cp, cr, cpg,
								   grid_full->gridF);
}
#endif

double GRID_D3_gauge_action(GRID_info_t *info,
			    su3_matrix *Umilc,
			    double beta,
			    double cp,
			    double cr,
			    double cpg,
			    GRID_4Dgrid *grid_full)
{
  return gaugeAction<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, Umilc, beta, cp, cr, cpg,
								   grid_full->gridD);
}


