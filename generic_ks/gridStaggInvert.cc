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
  ImprovedStaggeredFermion Ds(*(asqtad->lnglinks), *(asqtad->lnglinks), *(asqtad->fatlinks), *CGrid, *RBGrid, 2.*mass);

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
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::seconds>(elapsed) 
		  << "seconds\n";
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
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::seconds>(elapsed) 
		  << "seconds\n";
#if 0
	FermionField tmp(RBGrid);
	HermOp.Mpc(*(out->cv), tmp);
	Real check_resid = axpy_norm(tmp, -1., *(in->cv), tmp)/norm2(*(in->cv));
	std::cout << "check resid " << sqrt(check_resid) << "\n";
#endif
      }
    }
}

#if 0
template<typename FT>
static int 
asqtadInvert_mrhs (GRID_info_t *info, GRID_FermionLinksAsqtad_struct<ImprovedStaggeredFermion> *asqtad, 
		   GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg,
		   GRID_Real<FT> mass, GRID_ColorVector5D<ImprovedStaggeredFermion> *out, GRID_ColorVector5D<ImprovedStaggeredFermion> *in)
{
  GRID_ASSERT((inv_arg->parity  == GRID_EVENODD) || (inv_arg->parity  == GRID_EVEN) || (inv_arg->parity  == GRID_ODD), GRID_FAIL);
  GRID_ASSERT(in->cv->_grid == out->cv->_grid,  GRID_FAIL);

  ImprovedStaggeredFermion Ds(*(asqtad->thnlinks), *(asqtad->lnglinks), *(asqtad->fatlinks), 
				       *CGrid, *RBGrid, mass);
  ImprovedStaggeredFermion5DR Ds(Umu,Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass);
  
  *(out->cv) = zero;


  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max*inv_arg->nrestart);
  BlockConjugateGradient<FermionField> BCG(res_arg->resid, inv_arg->max*inv_arg->nrestart);
  MultiRHSConjugateGradient<FermionField> mCG(res_arg->resid, inv_arg->max*inv_arg->nrestart);


  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {
	MdagMLinearOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	CG(HermOp, *(in->cv), *(out->cv));
	
	break;
      }
      
    case GRID_EVEN:
    case GRID_ODD:
      {
	if (in->cv->_grid == RBGrid) {
	  SchurDiagMooeeOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	  CG(HermOp, *(in->cv), *(out->cv));
	} else {
	  FermionField *srcHalf = new FermionField(RBGrid);
	  FermionField *outHalf = new FermionField(RBGrid);
	  pickCheckerboard(inv_arg->parity,*srcHalf,*in->cv);
	  SchurDiagMooeeOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	  
	  auto start = std::chrono::system_clock::now();
	  
	  CG(HermOp, *srcHalf, *outHalf);
	  
	  auto end = std::chrono::system_clock::now();
	  auto elapsed = end - start;
	  std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::seconds>(elapsed) 
		    << "seconds\n";
	  
	  setCheckerboard(*out->cv, *outHalf);
	  delete	srcHalf, outHalf;
	}
	//SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);
	//SchurSolver(Ds,*(in->cv),*(out->cv));
	
	break;
      }
    }
}
#endif

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

#if 0
// Multimass inverter

void GRID_F3_asqtad_invert_multi(GRID_info_t *info,
				 GRID_F3_FermionLinksAsqtad *asqtad,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg,
				 float mass[],
				 int nmass,
				 GRID_F3_ColorVector *out[],
				 GRID_F3_ColorVector *in[])
{
  asqtadInvertMulti<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg, res_arg, mass, nmass, out, in);
}


void GRID_F3_asqtad_invert(GRID_info_t *info,
			   GRID_F3_FermionLinksAsqtad *asqtad,
			   GRID_invert_arg_t *inv_arg,
			   GRID_resid_arg_t *res_arg,
			   float mass,
			   GRID_F3_ColorVector *out,
			   GRID_F3_ColorVector *in)
{
  asqtadInvertMulti<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg, res_arg, mass, nmass, out, in);
}
#endif
