#include <Grid/Grid.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"

extern "C" {
#include "generic_ks_includes.h"
}

#include "../include/mGrid/mGrid_assert.h"
#include "../include/mGrid/mGrid_internal.h"

using namespace Grid;
using namespace Grid::QCD;

extern Grid::GridCartesian         *CGrid;
extern Grid::GridRedBlackCartesian *RBGrid;

template<typename FT, typename LatticeGaugeField, typename ImprovedStaggeredFermion>
static void
asqtadInvert (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
	      GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg, FT mass, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in )
{
  typedef typename ImprovedStaggeredFermion::FermionField FermionField;

  // Must recognize the parity flag 
  GRID_ASSERT((inv_arg->parity  == GRID_EVENODD) || (inv_arg->parity  == GRID_EVEN) || 
	      (inv_arg->parity  == GRID_ODD), GRID_FAIL);
  // In and out fields must be on the same lattice
  GRID_ASSERT(in->cv->_grid == out->cv->_grid,  GRID_FAIL);

  // Note: the first argument is ignored here
  auto start = std::chrono::system_clock::now();
  ImprovedStaggeredFermion Ds(*(asqtad->lnglinks), *(asqtad->lnglinks), *(asqtad->fatlinks), *CGrid, *RBGrid, 2.*mass);
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Instantiate ImprovedStaggeredFermion Ds " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";

  // Instantiate the inverter. The last arg = false says don't abort if no convergence 
  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max*inv_arg->nrestart, false);

  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {
	GRID_ASSERT((in->cv->_grid == CGrid) && (out->cv->_grid == CGrid), GRID_FAIL);

	std::cout << "WARNING: inversion with EVENODD is untested.\n";
	MdagMLinearOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	auto start = std::chrono::system_clock::now();
	CG(HermOp, *(in->cv), *(out->cv));
	res_arg->final_iter = CG.IterationsToComplete;
	res_arg->final_rsq = CG.TrueResidual*CG.TrueResidual;
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
	break;
      }
      
    case GRID_EVEN:
    case GRID_ODD:
      {
	// For EVEN or ODD parity we support only the case that in and out fields are checkerboards
	// and the checkerboards must be the same
	GRID_ASSERT((in->cv->_grid == RBGrid) && (out->cv->_grid == RBGrid)
		    && (in->cv->checkerboard == out->cv->checkerboard), GRID_FAIL);

	SchurStaggeredOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	
	auto start = std::chrono::system_clock::now();
	CG(HermOp, *(in->cv), *(out->cv));
	res_arg->final_iter = CG.IterationsToComplete;
	res_arg->final_rsq = CG.TrueResidual*CG.TrueResidual;
	std::cout << "iters = " << CG.IterationsToComplete << "\n" << std::flush;
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
#if 0
	FermionField tmp(RBGrid);
	HermOp.Mpc(*(out->cv), tmp);
	Real check_resid = axpy_norm(tmp, -1., *(in->cv), tmp)/norm2(*(in->cv));
	std::cout << "check resid " << sqrt(check_resid) << "\n";
#endif
      }
    }
}

template<typename FT, typename LatticeGaugeField, typename ImprovedStaggeredFermion>
static void
asqtadInvertMulti (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
		   GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg[],
		   FT mass[], int nmass, struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out[], 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in )
{
  // In and out fields must be on the same lattice
  //  std::cout << "Checking that out and in belong to the same grid\n" << std::flush;
  //  std::cout << "in is " << in->cv->_grid << "\n" << std::flush;
  for(int i = 0; i < nmass; i++){
    //    std::cout << "out[" << i << "] is " << out[i]->cv->_grid << "\n" << std::flush;
    GRID_ASSERT(in->cv->_grid == out[i]->cv->_grid,  GRID_FAIL);
  }

  typedef typename ImprovedStaggeredFermion::FermionField FermionField;

  std::vector<RealD> poles(nmass);
  std::vector<RealD> tolerances(nmass);
  for(int i = 0; i < nmass; i++){
    poles[i] = 4.*mass[i]*mass[i];
    tolerances[i] = Grid::sqrt(res_arg[i]->resid);
  }
  
  MultiShiftFunction Shifts(nmass, 0., 0.);
  Shifts.order = nmass;
  Shifts.norm = 1.;
  Shifts.poles = poles;
  Shifts.tolerances = tolerances;

  auto start = std::chrono::system_clock::now();
  ImprovedStaggeredFermion Ds(*(asqtad->lnglinks), *(asqtad->lnglinks), *(asqtad->fatlinks), 
			      *CGrid, *RBGrid, 0.);
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Instantiate ImprovedStaggeredFermion Ds " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
  info->misc_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	
  ConjugateGradientMultiShift<FermionField> MSCG(inv_arg->max*inv_arg->nrestart, Shifts);

  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {
	GRID_ASSERT((in->cv->_grid == CGrid), GRID_FAIL);
	std::vector<FermionField> outvec(nmass, CGrid);
#if 0
	auto start = std::chrono::system_clock::now();
	for(int i = 0; i < nmass; i++){
	  GRID_ASSERT((out[i]->cv->_grid == CGrid), GRID_FAIL);
	  outvec[i] = *(out[i]->cv);
	}
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	std::cout << "Pack initial guess in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
#endif

	MdagMLinearOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);

	auto start = std::chrono::system_clock::now();
	int iters = 0;
	MSCG(HermOp, *(in->cv), outvec);
	for(int i = 0; i < nmass; i++){
	  res_arg[i]->final_iter = MSCG.IterationsToComplete[i];
	  iters += MSCG.IterationsToComplete[i];
	  res_arg[i]->final_rsq = MSCG.TrueResidual[i]*MSCG.TrueResidual[i];
	}
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
	for(int i = 0; i < nmass; i++)
	  outvec[i] = *(out[i]->cv);
	break;
      }
      
    case GRID_EVEN:
    case GRID_ODD:
      {
	// For EVEN or ODD parity we support only the case that in and out fields are checkerboards
	// and the checkerboards is the same
	GRID_ASSERT(in->cv->_grid == RBGrid, GRID_FAIL);
	std::vector<FermionField> outvec(nmass, RBGrid);
#if 0
	auto start = std::chrono::system_clock::now();
	for(int i = 0; i < nmass; i++){
	  GRID_ASSERT((out[i]->cv->_grid == RBGrid) && 
		      (in->cv->checkerboard == out[i]->cv->checkerboard), GRID_FAIL);
	  outvec[i] = *(out[i]->cv);
	}
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	std::cout << "Pack initial guess in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
#endif
	
	SchurStaggeredOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	
	// Do the solve
	start = std::chrono::system_clock::now();
	MSCG(HermOp, *(in->cv), outvec);
	end = std::chrono::system_clock::now();
	elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed)
		  << "\n" << std::flush;
	
	int iters = 0;
	for(int i = 0; i < nmass; i++){
	  res_arg[i]->final_iter = MSCG.IterationsToComplete[i];
	  iters += MSCG.IterationsToComplete[i];
	  res_arg[i]->final_rsq = MSCG.TrueResidual[i]*MSCG.TrueResidual[i];
	}
	std::cout << "iters summed over masses = " << iters << "\n" << std::flush;

	// Unpack the solution
	start = std::chrono::system_clock::now();
	for(int i = 0; i < nmass; i++)
	  *(out[i]->cv) = outvec[i];
	end = std::chrono::system_clock::now();
	elapsed = end - start;
	std::cout << "Unpack outved in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";

	break;
      }
    default:
      {
	node0_printf("asqtadInvertMulti: Unrecognized parity %d\n", inv_arg->parity);
	terminate(1);
      }
    }
  
}

//  BlockConjugateGradient<FermionField> BCG(res_arg->resid, inv_arg->max*inv_arg->nrestart);
//  MultiRHSConjugateGradient<FermionField> mCG(res_arg->resid, inv_arg->max*inv_arg->nrestart);

/*	HAZ WRAPPERS PARA EL MULTIRHS Y EL BLOCKCG, LUEGO SINGLE PRECISION Y VUELVES A MIRAR	*/


//====================================================================//
// The GRID API for the inverter

// Single mass inverter

void GRID_F3_asqtad_invert(GRID_info_t *info,
			   GRID_F3_FermionLinksAsqtad *asqtad,
			   GRID_invert_arg_t *inv_arg,
			   GRID_resid_arg_t *res_arg,
			   float mass,
			   GRID_F3_ColorVector *out,
			   GRID_F3_ColorVector *in)
{
  asqtadInvert<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg, res_arg, mass, out, in);
}

void GRID_D3_asqtad_invert(GRID_info_t *info,
			   GRID_D3_FermionLinksAsqtad *asqtad,
			   GRID_invert_arg_t *inv_arg,
			   GRID_resid_arg_t *res_arg,
			   double mass,
			   GRID_D3_ColorVector *out,
			   GRID_D3_ColorVector *in)
{
  asqtadInvert<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg, res_arg, mass, out, in);
}

// Multimass inverter

void GRID_F3_asqtad_invert_multi(GRID_info_t *info,
				 GRID_F3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg[],
				 float mass[], int nmass,
				 GRID_F3_ColorVector *out[],
				 GRID_F3_ColorVector *in)
{
  asqtadInvertMulti<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg, res_arg, mass, nmass, out, in);
}


void GRID_D3_asqtad_invert_multi(GRID_info_t *info,
				 GRID_D3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg[],
				 double mass[], int nmass,
				 GRID_D3_ColorVector *out[],
				 GRID_D3_ColorVector *in)
{
  asqtadInvertMulti<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg, res_arg, mass, nmass, out, in);
}
