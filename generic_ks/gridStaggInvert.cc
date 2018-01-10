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

// extern Grid::GridCartesian         *CGrid;
// extern Grid::GridRedBlackCartesian *RBGrid;

template<typename FT, typename LatticeGaugeField, typename ImprovedStaggeredFermion>
static void
asqtadInvert (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
	      GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg, FT mass, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in,
	      GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  typedef typename ImprovedStaggeredFermion::FermionField FermionField;

  // Must recognize the parity flag 
  GRID_ASSERT((inv_arg->parity  == GRID_EVENODD) || (inv_arg->parity  == GRID_EVEN) || 
	      (inv_arg->parity  == GRID_ODD), GRID_FAIL);
  // In and out fields must be on the same lattice
  GRID_ASSERT(in->cv->_grid == out->cv->_grid,  GRID_FAIL);

  GridCartesian *CGrid = grid_full->grid;
  GridRedBlackCartesian *RBGrid = grid_rb->grid;

  // Note: the first argument is ignored here
  auto start = std::chrono::system_clock::now();
  ImprovedStaggeredFermion Ds(*(asqtad->lnglinks), *(asqtad->lnglinks), *(asqtad->fatlinks), 
			      *CGrid, *RBGrid, 2.*mass);
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
	auto end = std::chrono::system_clock::now();
	res_arg->final_iter = CG.IterationsToComplete;
	res_arg->final_rsq = CG.TrueResidual*CG.TrueResidual;
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
	auto end = std::chrono::system_clock::now();
	res_arg->final_iter = CG.IterationsToComplete;
	res_arg->final_rsq = CG.TrueResidual*CG.TrueResidual;
	auto elapsed = end - start;
	std::cout << "iters = " << CG.IterationsToComplete << "\n" << std::flush;
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
		   FT mass[], int nmass, 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out[], 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in,
		   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{

  GridCartesian *CGrid = grid_full->grid;
  GridRedBlackCartesian *RBGrid = grid_rb->grid;

  // In and out fields must be on the same lattice

  for(int i = 0; i < nmass; i++){
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

	MdagMLinearOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);

	auto start = std::chrono::system_clock::now();
	MSCG(HermOp, *(in->cv), outvec);
	auto end = std::chrono::system_clock::now();
	int iters = 0;
	for(int i = 0; i < nmass; i++){
	  res_arg[i]->final_iter = MSCG.IterationsToComplete[i];
	  iters += MSCG.IterationsToComplete[i];
	  res_arg[i]->final_rsq = MSCG.TrueResidual[i]*MSCG.TrueResidual[i];
	}
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
	std::cout << "Unpack outvec in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
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

template<typename FT, typename LatticeGaugeField, typename ImprovedStaggeredFermion5D>
static void
asqtadInvertBlock (GRID_info_t *info, 
		   struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
		   GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg,
		   FT mass, int nrhs, 
		   struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *out, 
		   struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5D> *in,
		   GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
		   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  // In and out fields must be on the same lattice
  GRID_ASSERT(in->cv->_grid == out->cv->_grid,  GRID_FAIL);

  GridCartesian *FCGrid = grid_5D->grid;
  GridRedBlackCartesian *FRBGrid = grid_5Drb->grid;
  GridCartesian *CGrid = grid_full->grid;
  GridRedBlackCartesian *RBGrid = grid_rb->grid;

  typedef typename ImprovedStaggeredFermion5D::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermion5D::ComplexField ComplexField; 

  // Call using c1 = c2 = 2. and u0 = 1. to neutralize link rescaling
  ImprovedStaggeredFermion5D Ds(*(asqtad->lnglinks), *(asqtad->fatlinks), *FCGrid, *FRBGrid, 
				*CGrid, *RBGrid, mass, 2., 2., 1.);
  std::cout << "instantiating 5D CG with resid " << res_arg->resid << " and " << inv_arg->max*inv_arg->nrestart << " iters\n" << std::flush;
  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max*inv_arg->nrestart, false);

  int blockDim = 0;

  BlockConjugateGradient<FermionField>    BCGrQ(BlockCGrQ, blockDim, res_arg->resid, inv_arg->max*inv_arg->nrestart);
  BlockConjugateGradient<FermionField>    BCG  (BlockCG, blockDim, res_arg->resid, inv_arg->max*inv_arg->nrestart);
  BlockConjugateGradient<FermionField>    mCG  (CGmultiRHS, blockDim, res_arg->resid, inv_arg->max*inv_arg->nrestart);

  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {

	MdagMLinearOperator<ImprovedStaggeredFermion5D,FermionField> HermOp(Ds);
	CG(HermOp, *(in->cv), *(out->cv));

	break;
      }
      
    case GRID_EVEN:
    case GRID_ODD:
      {

	SchurStaggeredOperator<ImprovedStaggeredFermion5D,FermionField> HermOp(Ds);

	int blockDim = 0;

	std::cout << "Calling 5D CG for " << nrhs << " sources\n" << std::flush;

	// 5D CG
	Ds.ZeroCounters();
	//	*(out->cv) = zero;
	std::cout << "Running CG\n" << std::flush;
	auto start = std::chrono::system_clock::now();
	CG(HermOp, *(in->cv), *(out->cv));
	auto end = std::chrono::system_clock::now();
	res_arg->final_iter = CG.IterationsToComplete;
	res_arg->final_rsq = CG.TrueResidual*CG.TrueResidual;
	auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
	Ds.Report();
	
#if 0
	// multiRHS
	Ds.ZeroCounters();
	*(out->cv) = zero;
	mCG(HermOp, *(in->cv), *(out->cv));
	Ds.Report();
	
	// Block CG
	Ds.ZeroCounters();
	*(out->cv) = zero;
	BCGrQ(HermOp, *(in->cv), *(out->cv));
	Ds.Report();
#endif	
	break;
      }

    default:
      {
	node0_printf("asqtadInvertBlock: Unrecognized parity %d\n", inv_arg->parity);
	terminate(1);
      }
    }

  std::cout << "Leaving asqtadInvertBlock\n" << std::flush;
}
	
//====================================================================//
// The GRID API for the inverter

// Single mass inverter

void GRID_F3_asqtad_invert(GRID_info_t *info,
			   GRID_F3_FermionLinksAsqtad *asqtad,
			   GRID_invert_arg_t *inv_arg,
			   GRID_resid_arg_t *res_arg,
			   float mass,
			   GRID_F3_ColorVector *out,
			   GRID_F3_ColorVector *in,
			   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvert<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg, res_arg, mass, out, in, grid_full, grid_rb);
}

void GRID_D3_asqtad_invert(GRID_info_t *info,
			   GRID_D3_FermionLinksAsqtad *asqtad,
			   GRID_invert_arg_t *inv_arg,
			   GRID_resid_arg_t *res_arg,
			   double mass,
			   GRID_D3_ColorVector *out,
			   GRID_D3_ColorVector *in,
			   GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvert<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg, res_arg, mass, out, in, grid_full, grid_rb);
}

// Multimass inverter

void GRID_F3_asqtad_invert_multi(GRID_info_t *info,
				 GRID_F3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg[],
				 float mass[], int nmass,
				 GRID_F3_ColorVector *out[],
				 GRID_F3_ColorVector *in,
				 GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvertMulti<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg, res_arg, mass, nmass, out, in, grid_full, grid_rb);
}


void GRID_D3_asqtad_invert_multi(GRID_info_t *info,
				 GRID_D3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg[],
				 double mass[], int nmass,
				 GRID_D3_ColorVector *out[],
				 GRID_D3_ColorVector *in,
				 GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvertMulti<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg, res_arg, mass, nmass, out, in, grid_full, grid_rb);
}

// Block CG inverter

void GRID_F3_asqtad_invert_block(GRID_info_t *info,
				 GRID_F3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg,
				 float mass, int nrhs,
				 GRID_F3_ColorVectorBlock *out,
				 GRID_F3_ColorVectorBlock *in,
				 GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
				 GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvertBlock<float, LatticeGaugeFieldF, ImprovedStaggeredFermion5DF>(info, asqtad, inv_arg, res_arg, mass, nrhs, out, in, grid_5D, grid_5Drb, grid_full, grid_rb);
}


void GRID_D3_asqtad_invert_block(GRID_info_t *info,
				 GRID_D3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg,
				 double mass, int nrhs,
				 GRID_D3_ColorVectorBlock *out,
				 GRID_D3_ColorVectorBlock *in,
				 GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
				 GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvertBlock<double, LatticeGaugeFieldD, ImprovedStaggeredFermion5DD>(info, asqtad, inv_arg, res_arg, mass, nrhs, out, in, grid_5D, grid_5Drb, grid_full, grid_rb);
}
