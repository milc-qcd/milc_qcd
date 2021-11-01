//------------- gridStaggEigen.cc --------------
// MILC-Grid interfaces for eigensolver routines
// 
// impResLanczos() : implicitly restarted Lanczos
//----------------------------------------------

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h>

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"

using namespace std;
using namespace Grid;

// #define EIG_DEBUG
// #define EIG_TIME

// Run implicitly restarted Lanczos routine in Grid.
// Mainly use ImplicitlyRestartedLanczos class.
//
// FT: floating point number type.
// LatticeGaugeField: lattice gauge field class.
// ImprovedStaggeredFermion: improved staggered fermion class.
//
// eigVecs: resulting eigenvectors will be saved
// eigVals: resulting eigenvalues will be saved
// link: improved fermion links
// eig_arg: parameters defined in include/mGrid/mGrid_int.h
// mass: valence quark mass
// FGrid: full grid
// RBGrid: checkerborad (half) grid
//
template< typename FT,
          typename LatticeGaugeField,
          typename ImprovedStaggeredFermion >
static void impResLanczos( GRID_ColorVectorArray_struct<ImprovedStaggeredFermion> * eigVecs,
                           FT * eigVals,
                           GRID_FermionLinksAsqtad_struct<LatticeGaugeField> * link,
                           GRID_eig_arg_t * eig_arg,
                           FT mass,
                           GridCartesian * FGrid,
                           GridRedBlackCartesian * RBGrid )
{
  char myname[] = "impResLanczos";
  int ii; // loop index
  
#ifdef EIG_DEBUG
  cout << myname << ": Start" << endl;
#endif
#ifdef EIGTIME
  RealD dtimec = 0;
#endif
  
  typedef typename ImprovedStaggeredFermion::FermionField FermionField;

  // Coefficients for improved staggered action
  // They are already considered in smearing
  RealD c1 = 1.0;
  RealD c2 = 1.0;
  // Tadpole improvement
  RealD u0 = 1.0; 

  // Dirac operator
  // Factor two comes to match convention
  ImprovedStaggeredFermion Ds( *FGrid, *RBGrid, 2.*mass, 2.*c1, 2.*c2, u0 ); 
  Ds.ImportGaugeSimple( *(link->lnglinks), *(link->fatlinks) );

  // Even-odd preconditioned D^dagger D operator
  SchurStaggeredOperator< ImprovedStaggeredFermion, FermionField > HermOp(Ds);

  // Parameters for Chebyshev polynomial
  ChebyParams chebyParams;
  chebyParams.alpha = eig_arg->chebyParams.alpha;
  chebyParams.beta = eig_arg->chebyParams.beta;
  chebyParams.Npoly = eig_arg->chebyParams.Npoly;

  // Chebyshev polynomial
  Chebyshev< FermionField > chebyshev( chebyParams );
  FunctionHermOp< FermionField > poly( chebyshev, HermOp );

  // Original operator to calculate true eigenvalues by Rayleigh quotient
  PlainHermOp< FermionField > op( HermOp );

  // Parameters for implicitly restarted Lanczos
  int Nstop = eig_arg->Nstop;
  int Nm = eig_arg->Nm;
  int Nk = eig_arg->Nk;
  FT tol = eig_arg->tol;
  int maxIter = eig_arg->maxIter;
  int restartMin = eig_arg->restartMin;
  int reorth_period = eig_arg->reorth_period;
  FT betastp = 0; // not used
  IRLdiagonalisation diag = (IRLdiagonalisation) eig_arg->diag;

  // ImplicitlyRestartedLanczos class defined in
  // Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h
  ImplicitlyRestartedLanczos< FermionField >
    irl( poly, op, Nstop, Nk, Nm, tol, maxIter, betastp, restartMin, reorth_period, diag );

  // eigenvalues and eigenvectors 
  vector< RealD > _eigVals( Nm, 0 );
  eigVecs->cv->resize( Nm, RBGrid );

  // initial vector for Lanczos
  FermionField src(RBGrid);
  src.Checkerboard() = eig_arg->parity;

  // assign initial vector values
  autoView( src_v, src, AcceleratorWrite );
  for( ii=0; ii<RBGrid->oSites(); ii+=1 )
  {
    src_v[ii] = 1.0;
  }

  // Total number of converged eigenvalues
  int Nconv;

#ifdef EIGTIME
  dtimec -= usecond()/1.0e6;
#endif
  
  // Calculate eigenvalues and eigenvectors
  irl.calc( _eigVals, *(eigVecs->cv), src, Nconv );

#ifdef EIGTIME
  dtimec += usecond()/1.0e6;
  cout << "[irl.calc] time = " << dtimec << " s" << endl;
#endif

  // Copy resulting eigenvalues
  for( ii=0; ii<Nstop; ii+=1 )
  {
    *(eigVals+ii) = _eigVals[ii];
  }
  
#ifdef EIG_DEBUG
  cout << myname << ": End" << endl;
#endif
  
  return ;
}

// Double preicision version wrapper
void GRID_D3_implicitly_restarted_lanczos(
  GRID_D3_ColorVectorArray * eigVecs,
  double * eigVals,
  GRID_D3_FermionLinksAsqtad * link,
  GRID_eig_arg_t * eig_arg,
  double mass,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb )
{
  impResLanczos< RealD,
                 LatticeGaugeFieldD,
                 ImprovedStaggeredFermionD >
    ( eigVecs, eigVals, link, eig_arg, mass, grid_full->gridD, grid_rb->gridD );

  return ;
}

// Single preicision version wrapper
void GRID_F3_implicitly_restarted_lanczos(
  GRID_F3_ColorVectorArray * eigVecs,
  float * eigVals,
  GRID_F3_FermionLinksAsqtad * link,
  GRID_eig_arg_t * eig_arg,
  float mass,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb )
{
  impResLanczos< RealF,
                 LatticeGaugeFieldF,
                 ImprovedStaggeredFermionF >
    ( eigVecs, eigVals, link, eig_arg, mass, grid_full->gridD, grid_rb->gridD );
  
  return ;
}
