// Wrappers for Grid Staggered Inverters

#include <omp.h>
#include <Grid/Grid.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"

#define GRID_5DCG    0
#define GRID_MRHSCG  1
#define GRID_BLOCKCG 2

#include "../include/mGrid/mGrid_assert.h"
#include "../include/macros.h"

using namespace Grid;
//using namespace Grid::QCD;

typedef typename ImprovedStaggeredFermion5DF::FermionField FermionField5DF;
typedef typename ImprovedStaggeredFermion5DD::FermionField FermionField5DD;

template<typename Field>
static void
restartCG (int nrestart, GRID_resid_arg_t *res_arg, 
	   ConjugateGradient<Field> &CG, LinearOperatorBase<Field> &LinOp,
	   const Field &in, Field &out)
{
  int i;
#ifdef CG_DEBUG
  GridLogIterative.Active(1);   // Turns on iterative logging. Normally off.
#endif
  res_arg->final_iter = 0;
  for(i = 0; i < nrestart; i++){
    CG(LinOp, in, out);
    res_arg->final_iter += CG.IterationsToComplete;
    res_arg->final_rsq = CG.TrueResidual*CG.TrueResidual;
    std::cout << "Comparing true residual " << CG.TrueResidual << " with desired residual" <<  res_arg->resid << std::endl;
    if(CG.TrueResidual < res_arg->resid)break;
    std::cout << "Restart " << i << " after iteration " 
	      << res_arg->final_iter << ", true residual is "
	      << CG.TrueResidual << " target " << res_arg->resid << "\n";
  }
  res_arg->final_restart = i;
}

template<typename FieldD, typename FieldF>
static void
restartMixedCG (int nrestart, GRID_resid_arg_t *res_arg, 
		MixedPrecisionConjugateGradient<FieldD, FieldF> &MPCG, 
		const FieldD &in, FieldD &out)
{
  int i;
#ifdef CG_DEBUG
  GridLogIterative.Active(1);   // Turns on iterative logging. Normally off.
#endif
  res_arg->final_iter = 0;
  for(i = 0; i < nrestart; i++){
    MPCG(in, out);
    res_arg->final_iter += MPCG.TotalInnerIterations + MPCG.TotalFinalStepIterations + MPCG.TotalOuterIterations;
    res_arg->final_rsq = MPCG.TrueResidual*MPCG.TrueResidual;
    if(MPCG.TrueResidual < res_arg->resid)break;
    std::cout << "Restart " << i << " after iteration " 
	      << res_arg->final_iter << ", true residual is "
	      << MPCG.TrueResidual << " target " << res_arg->resid << "\n";
  }
  res_arg->final_restart = i;
}

template<typename FT, typename LatticeGaugeField, typename ImprovedStaggeredFermion>
static void
asqtadInvert (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
	      GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg, FT mass, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out, 
	      struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in,
	      GridCartesian *CGrid, GridRedBlackCartesian *RBGrid)
{
  typedef typename ImprovedStaggeredFermion::FermionField FermionField;

  // Must recognize the parity flag 
  GRID_ASSERT((inv_arg->parity  == GRID_EVENODD) || (inv_arg->parity  == GRID_EVEN) || 
	      (inv_arg->parity  == GRID_ODD), GRID_FAIL);
  // In and out fields must be on the same lattice
  GRID_ASSERT(in->cv->Grid() == out->cv->Grid(),  GRID_FAIL);

  auto start = std::chrono::system_clock::now();
  // Call using c1 = c2 = 2. and u0 = 1. to neutralize link rescaling -- probably ignored anyway.
  ImprovedStaggeredFermion Ds(*CGrid, *RBGrid, 2.*mass, 2., 2., 1.);
  Ds.ImportGaugeSimple(*(asqtad->lnglinks), *(asqtad->fatlinks));
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Instantiate ImprovedStaggeredFermion Ds " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";

  // Instantiate the inverter. The last arg = false says don't abort if no convergence 
  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max, false);

  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {
	GRID_ASSERT((in->cv->Grid() == CGrid) && (out->cv->Grid() == CGrid), GRID_FAIL);

	std::cout << "WARNING: inversion with EVENODD is untested.\n";
	MdagMLinearOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	auto start = std::chrono::system_clock::now();
	restartCG<FermionField>(inv_arg->nrestart, res_arg, CG, HermOp, 
				*(in->cv), *(out->cv));
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
	GRID_ASSERT((in->cv->Grid() == RBGrid) && (out->cv->Grid() == RBGrid)
		    && (in->cv->Checkerboard() == out->cv->Checkerboard()), GRID_FAIL);

	SchurStaggeredOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);
	
	auto start = std::chrono::system_clock::now();
	restartCG<FermionField>(inv_arg->nrestart, res_arg, CG, HermOp, 
				*(in->cv), *(out->cv));
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

static void
asqtadInvertMixed (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldF> *asqtad_f,
		   struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldD> *asqtad_d, 
		   GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg, RealD mass, 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermionD> *out, 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermionD> *in,
		   GridCartesian *CGrid_f, GridRedBlackCartesian *RBGrid_f,
		   GridCartesian *CGrid_d, GridRedBlackCartesian *RBGrid_d)
{

  typedef typename ImprovedStaggeredFermionF::FermionField FermionFieldF;
  typedef typename ImprovedStaggeredFermionD::FermionField FermionFieldD;

// Must recognize the parity flag  We don't support EVENODD
  GRID_ASSERT((inv_arg->parity  == GRID_EVEN) || 
	      (inv_arg->parity  == GRID_ODD), GRID_FAIL);

  // In and out fields must be on the same double-precision lattice 
  GRID_ASSERT((in->cv->Grid() == RBGrid_d) && (out->cv->Grid() == RBGrid_d)
	      && (in->cv->Checkerboard() == out->cv->Checkerboard()), GRID_FAIL);

  auto start = std::chrono::system_clock::now();

  // Set up action using c1 = c2 = 2. and u0 = 1. to neutralize link rescaling -- probably ignored anyway.
  ImprovedStaggeredFermionF Ds_f(*CGrid_f, *RBGrid_f, 2.*mass, 2., 2., 1.);
  Ds_f.ImportGaugeSimple(*(asqtad_f->lnglinks), *(asqtad_f->fatlinks));
  SchurStaggeredOperator<ImprovedStaggeredFermionF,FermionFieldF> HermOp_f(Ds_f);

  ImprovedStaggeredFermionD Ds_d(*CGrid_d, *RBGrid_d, 2.*mass, 2., 2., 1.);
  SchurStaggeredOperator<ImprovedStaggeredFermionD,FermionFieldD> HermOp_d(Ds_d);
  Ds_d.ImportGaugeSimple(*(asqtad_d->lnglinks), *(asqtad_d->fatlinks));

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Instantiate ImprovedStaggeredFermion Ds_d and Ds_f " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";

  // Instantiate and run the inverter
  MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF>
    MPCG(res_arg->resid, inv_arg->maxInner, inv_arg->max, RBGrid_f, HermOp_f, HermOp_d);

  start = std::chrono::system_clock::now();
  restartMixedCG<FermionFieldD, FermionFieldF>(inv_arg->nrestart, res_arg, MPCG, *(in->cv), *(out->cv));
  end = std::chrono::system_clock::now();
  elapsed = end - start;
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

template<typename FT, typename LatticeGaugeField, typename ImprovedStaggeredFermion>
static void
asqtadInvertMulti (GRID_info_t *info, struct GRID_FermionLinksAsqtad_struct<LatticeGaugeField> *asqtad, 
		   GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg[],
		   FT mass[], int nmass, 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *out[], 
		   struct GRID_ColorVector_struct<ImprovedStaggeredFermion> *in,
		   GridCartesian *CGrid, GridRedBlackCartesian *RBGrid)
{
  // In and out fields must be on the same lattice

  for(int i = 0; i < nmass; i++){
    GRID_ASSERT(in->cv->Grid() == out[i]->cv->Grid(),  GRID_FAIL);
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
  ImprovedStaggeredFermion Ds(*CGrid, *RBGrid, 0., 2., 2., 1.);
  Ds.ImportGaugeSimple(*(asqtad->lnglinks), *(asqtad->fatlinks));
  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Instantiate ImprovedStaggeredFermion Ds " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
  info->misc_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	
  ConjugateGradientMultiShift<FermionField> MSCG(inv_arg->max*inv_arg->nrestart, Shifts);

  GRID_ASSERT(inv_arg->parity == GRID_EVENODD || inv_arg->parity == GRID_EVEN || 
	      inv_arg->parity == GRID_ODD, GRID_FAIL);

  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {
	GRID_ASSERT((in->cv->Grid() == CGrid), GRID_FAIL);
	std::vector<FermionField> outvec(nmass, CGrid);

	MdagMLinearOperator<ImprovedStaggeredFermion,FermionField> HermOp(Ds);

	auto start = std::chrono::system_clock::now();
	MSCG(HermOp, *(in->cv), outvec);
	auto end = std::chrono::system_clock::now();
	int iters = 0;
	for(int i = 0; i < nmass; i++){
	  res_arg[i]->final_iter = MSCG.IterationsToCompleteShift[i];
	  iters += MSCG.IterationsToCompleteShift[i];
	  res_arg[i]->final_rsq = MSCG.TrueResidualShift[i]*MSCG.TrueResidualShift[i];
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
	GRID_ASSERT(in->cv->Grid() == RBGrid, GRID_FAIL);
	std::vector<FermionField> outvec(nmass, RBGrid);
	for(int i = 0; i < nmass; i++)
	  outvec[i].Checkerboard() = inv_arg->parity == EVEN ? Even : Odd ;
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
	  res_arg[i]->final_iter = MSCG.IterationsToCompleteShift[i];
	  iters += MSCG.IterationsToCompleteShift[i];
	  res_arg[i]->final_rsq = MSCG.TrueResidualShift[i]*MSCG.TrueResidualShift[i];
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
		   GridCartesian *FCGrid, GridRedBlackCartesian *FRBGrid,
		   GridCartesian *CGrid, GridRedBlackCartesian *RBGrid)
{
#if 0
  std::cout << "Block inverter currently unavailable" << std::endl;
#else
  // In and out fields must be on the same lattice
  GRID_ASSERT(in->cv->Grid() == out->cv->Grid(),  GRID_FAIL);

  typedef typename ImprovedStaggeredFermion5D::FermionField FermionField; 

  // Call using c1 = c2 = 2. and u0 = 1. to neutralize link rescaling -- probably ignored anyway.
  ImprovedStaggeredFermion5D Ds(*FCGrid, *FRBGrid, *CGrid, *RBGrid, 2.*mass, 2., 2., 1.);
  Ds.ImportGaugeSimple(*(asqtad->lnglinks), *(asqtad->fatlinks));
  std::cout << "Instantiating 5D CG with resid " << res_arg->resid << " and " << inv_arg->max << " iters\n" << std::flush;
  ConjugateGradient<FermionField> CG(res_arg->resid, inv_arg->max, false);

  int blockDim = 0;

  BlockConjugateGradient<FermionField>    BCGrQ(BlockCGrQ, blockDim, res_arg->resid, inv_arg->max*inv_arg->nrestart);
  BlockConjugateGradient<FermionField>    BCG  (BlockCG, blockDim, res_arg->resid, inv_arg->max*inv_arg->nrestart);
  BlockConjugateGradient<FermionField>    mCG  (CGmultiRHS, blockDim, res_arg->resid, inv_arg->max*inv_arg->nrestart);

  GRID_ASSERT(inv_arg->parity == GRID_EVENODD || inv_arg->parity == GRID_EVEN || 
	      inv_arg->parity == GRID_ODD, GRID_FAIL);

  switch (inv_arg->parity)
    {
    case GRID_EVENODD:
      {

	MdagMLinearOperator<ImprovedStaggeredFermion5D,FermionField> HermOp(Ds);
	restartCG<FermionField>(inv_arg->nrestart, res_arg, CG, HermOp, 
				*(in->cv), *(out->cv));
	
	break;
      }
      
    case GRID_EVEN:
    case GRID_ODD:
      {

	SchurStaggeredOperator<ImprovedStaggeredFermion5D,FermionField> HermOp(Ds);

#if GRID_MULTI_CG == GRID_5DCG

	// 5D CG

	std::cout << "Running 5D CG for " << nrhs << " sources\n" << std::flush;
	auto start = std::chrono::system_clock::now();
	restartCG<FermionField>(inv_arg->nrestart, res_arg, CG, HermOp, 
				*(in->cv), *(out->cv));
	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";

#elif GRID_MULTI_CG == GRID_MRHSCG

	// multiRHS
	std::cout << "Running multiRHS CG for " << nrhs << " sources\n" << std::flush;
	auto start = std::chrono::system_clock::now();
	mCG(HermOp, *(in->cv), *(out->cv));
	auto end = std::chrono::system_clock::now();
	res_arg->final_iter = mCG.IterationsToComplete;
	res_arg->final_rsq = mCG.TrueResidual*mCG.TrueResidual;
	auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";

#else

	// Block CG
	std::cout << "Running Block CG for " << nrhs << " sources\n" << std::flush;
        auto start = std::chrono::system_clock::now();
	BCGrQ(HermOp, *(in->cv), *(out->cv));
        auto end = std::chrono::system_clock::now();
	res_arg->final_iter = BCGrQ.IterationsToComplete;
        res_arg->final_rsq = BCGrQ.TrueResidual*mCG.TrueResidual;
        auto elapsed = end - start;
	info->final_sec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/1000.;
	std::cout << "Inverted in " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
		  << "\n";
#endif	
	break;
      }
    }
#endif
}
	
static void
asqtadInvertMixedBlock (GRID_info_t *info,
			struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldF> *asqtad_f,
			struct GRID_FermionLinksAsqtad_struct<LatticeGaugeFieldD> *asqtad_d, 
			GRID_invert_arg_t *inv_arg, GRID_resid_arg_t *res_arg, RealD mass, int nrhs,
			struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5DD> *out, 
			struct GRID_ColorVectorBlock_struct<ImprovedStaggeredFermion5DD> *in,
			GridCartesian *FCGrid_f, GridRedBlackCartesian *FRBGrid_f,
			GridCartesian *FCGrid_d, GridRedBlackCartesian *FRBGrid_d,
			GridCartesian *CGrid_f, GridRedBlackCartesian *RBGrid_f,
			GridCartesian *CGrid_d, GridRedBlackCartesian *RBGrid_d)
{

  // Must recognize the parity flag  We don't support EVENODD
  GRID_ASSERT((inv_arg->parity  == GRID_EVEN) || 
	      (inv_arg->parity  == GRID_ODD), GRID_FAIL);

  // In and out fields must be on the same double-precision lattice 
  GRID_ASSERT(in->cv->Grid() == out->cv->Grid(), GRID_FAIL);
  
  //  GRID_ASSERT((in->cv->Grid() == RBGrid_d) && (out->cv->Grid() == RBGrid_d)
  //	      && (in->cv->Checkerboard() == out->cv->Checkerboard()), GRID_FAIL);
  
  auto start = std::chrono::system_clock::now();

  // Set up action using c1 = c2 = 2. and u0 = 1. to neutralize link rescaling -- probably ignored anyway.
  ImprovedStaggeredFermion5DF Ds_f(*FCGrid_f, *FRBGrid_f, *CGrid_f, *RBGrid_f, 2.*mass, 2., 2., 1.);
  Ds_f.ImportGaugeSimple(*(asqtad_f->lnglinks), *(asqtad_f->fatlinks));
  SchurStaggeredOperator<ImprovedStaggeredFermion5DF,FermionField5DF> HermOp_f(Ds_f);

  ImprovedStaggeredFermion5DD Ds_d(*FCGrid_d, *FRBGrid_d, *CGrid_d, *RBGrid_d, 2.*mass, 2., 2., 1.);
  Ds_d.ImportGaugeSimple(*(asqtad_d->lnglinks), *(asqtad_d->fatlinks));
  SchurStaggeredOperator<ImprovedStaggeredFermion5DD,FermionField5DD> HermOp_d(Ds_d);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Instantiate ImprovedStaggeredFermion Ds_d and Ds_f " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";

  // Instantiate and run the inverter
  MixedPrecisionConjugateGradient<FermionField5DD, FermionField5DF>
    MPCG(res_arg->resid, inv_arg->maxInner, inv_arg->max, FRBGrid_f, HermOp_f, HermOp_d);

  start = std::chrono::system_clock::now();
  restartMixedCG<FermionField5DD, FermionField5DF>(inv_arg->nrestart, res_arg, MPCG, *(in->cv), *(out->cv));
  end = std::chrono::system_clock::now();
  elapsed = end - start;
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
  asqtadInvert<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg,
			  res_arg, mass, out, in, grid_full->gridF, grid_rb->gridF);
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
  asqtadInvert<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg,
			   res_arg, mass, out, in, grid_full->gridD, grid_rb->gridD);
}

// Single mass, single source mixed-precision inverter

void GRID_D3_asqtad_invert_mixed(GRID_info_t *info,
				 GRID_F3_FermionLinksAsqtad *asqtad_f,
				 GRID_D3_FermionLinksAsqtad *asqtad_d,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg,
				 double mass,
				 GRID_D3_ColorVector *out,
				 GRID_D3_ColorVector *in,
				 GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvertMixed(info, asqtad_f, asqtad_d, inv_arg,
		    res_arg, mass, out, in, grid_full->gridF, grid_rb->gridF, grid_full->gridD, grid_rb->gridD);
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
  asqtadInvertMulti<float, LatticeGaugeFieldF, ImprovedStaggeredFermionF>(info, asqtad, inv_arg,
			       res_arg, mass, nmass, out, in, grid_full->gridF, grid_rb->gridF);
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
  asqtadInvertMulti<double, LatticeGaugeFieldD, ImprovedStaggeredFermionD>(info, asqtad, inv_arg,
				res_arg, mass, nmass, out, in, grid_full->gridD, grid_rb->gridD);
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
  asqtadInvertBlock<float, LatticeGaugeFieldF, ImprovedStaggeredFermion5DF>(info, asqtad, inv_arg,
			       res_arg, mass, nrhs, out, in, grid_5D->gridF, grid_5Drb->gridF,
			       grid_full->gridF, grid_rb->gridF);
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
  asqtadInvertBlock<double, LatticeGaugeFieldD, ImprovedStaggeredFermion5DD>(info, asqtad, inv_arg,
			        res_arg, mass, nrhs, out, in, grid_5D->gridD, grid_5Drb->gridD,
			        grid_full->gridD, grid_rb->gridD);
}

// Mixed-precision Block CG inverter (double-precision only)

void GRID_D3_asqtad_invert_mixed_block(GRID_info_t *info,
				       GRID_F3_FermionLinksAsqtad *asqtad_f,
				       GRID_D3_FermionLinksAsqtad *asqtad_d,
				       GRID_invert_arg_t *inv_arg,
				       GRID_resid_arg_t *res_arg,
				       double mass, int nrhs,
				       GRID_D3_ColorVectorBlock *out,
				       GRID_D3_ColorVectorBlock *in,
				       GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
				       GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb)
{
  asqtadInvertMixedBlock(info, asqtad_f, asqtad_d, inv_arg,
			 res_arg, mass, nrhs, out, in,
			 grid_5D->gridF, grid_5Drb->gridF, grid_5D->gridD, grid_5Drb->gridD,
			 grid_full->gridF, grid_rb->gridF, grid_full->gridD, grid_rb->gridD);
}

