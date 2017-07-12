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
static int 
asqtadInvert (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
	      GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg, FT mass, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in )
{
  typedef typename ImprovedStaggeredFermion::FermionField FermionField;

  GRID_ASSERT((inv_arg->parity  == GRID_EVENODD) || (inv_arg->parity  == GRID_EVEN) || 
	      (inv_arg->parity  == GRID_ODD), GRID_FAIL);
  GRID_ASSERT(in->cv->_grid == out->cv->_grid,  GRID_FAIL);

  ImprovedStaggeredFermion Ds(*(asqtad->thnlinks), *(asqtad->lnglinks), *(asqtad->fatlinks), 
			      *CGrid, *RBGrid, mass);

  //*(out->cv) = zero;

  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max);

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


  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max);
  BlockConjugateGradient<FermionField> BCG(1.0e-8,10000);
  MultiRHSConjugateGradient<FermionField> mCG(1.0e-8,10000);


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



