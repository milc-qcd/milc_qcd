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
#include "../include/mGrid/mGrid_assert.h"

#include "../generic/gridMap.h"
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/HISQSmearing.h>
#include <Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>

using namespace Grid;

#if 0
template<typename T>
static HISQParameters<T> get_hisq_param(int n_naiks,
					std::array<T,GRID_MAX_NAIK> eps_naiks, fermion_links_t *fl){

  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  Real fat7_c1    = ap->p1.act_path_coeff.one_link ;
  Real fat7_c3    = ap->p1.act_path_coeff.three_staple ;
  Real fat7_c5    = ap->p1.act_path_coeff.five_staple ;
  Real fat7_c7    = ap->p1.act_path_coeff.seven_staple ;

  Real asqtad_c1  = ap->p2.act_path_coeff.one_link ;
  Real asqtad_c3  = ap->p2.act_path_coeff.three_staple ;
  Real asqtad_c5  = ap->p2.act_path_coeff.five_staple ;
  Real asqtad_c7  = ap->p2.act_path_coeff.seven_staple ;
  Real asqtad_clp = ap->p2.act_path_coeff.lepage ;
  Real cnaik      = ap->p2.act_path_coeff.naik ;

  Real diff_c1    = ap->p3.act_path_coeff.one_link ;
  Real diff_cnaik = ap->p3.act_path_coeff.naik ;
  int ugroup     = ap->ugroup;
  int umethod    = ap->umethod;
  
  HISQParameters<T> hisq_param(n_naiks  , eps_naiks ,
	  fat7_c1  , fat7_c3  , fat7_c5  , fat7_c7  , 0.,
	  asqtad_c1, asqtad_c3, asqtad_c5, asqtad_c7, asqtad_clp,
	  cnaik    , diff_c1     , diff_cnaik);
  return hisq_param;
}

#endif

// residues and multi_x are indexed by the pseudofermion fields
// multi_x[i] points to a color vector field.
//   The fieldss for each Naik mass are grouped together and
//   the set is concatenated in the order of n_orders_naik.
// n_orders_naik gives the number of pseudofermion fields for each
//   Naik mass.
// deriv[] is indexed by the spacetime dimension

template<typename LatticeGaugeField, typename FermionField, typename Gimpl, typename Complex>
static void hisqForce (
  GRID_info_t* info,
  void* fl_void,
  Real residues[],
  su3_vector* multi_x[],
  int n_orders_naik[],
  su3_matrix* deriv,
  GridCartesian* CGrid
) {
  fermion_links_t* fl = (fermion_links_t*)fl_void;
  
  auto start = std::chrono::system_clock::now();

  hisq_auxiliary_t* aux = get_hisq_auxiliary(fl);
  su3_matrix* Umilc = aux->U_link;
  su3_matrix* Vmilc = aux->V_link;
  su3_matrix* Wmilc = aux->W_unitlink;

  LatticeGaugeField Umu(CGrid), Vmu(CGrid), Wmu(CGrid), UForce(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Umilc, &Umu);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Vmilc, &Vmu);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Wmilc, &Wmu);

  // -- coefficient preparation -- //

  ks_action_paths_hisq* ap = get_action_paths_hisq(fl);
  Real eigenvalue_cutoff = HISQ_FORCE_FILTER;
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
  
  HISFContext fatCtx(
    ap->p1.act_path_coeff.one_link,
    ap->p1.act_path_coeff.three_staple,
    ap->p1.act_path_coeff.five_staple,
    ap->p1.act_path_coeff.seven_staple,
    0.0,
    0.0
  );

  HISFContext asqCtx(
    ap->p2.act_path_coeff.one_link,
    ap->p2.act_path_coeff.three_staple,
    ap->p2.act_path_coeff.five_staple,
    ap->p2.act_path_coeff.seven_staple,
    ap->p2.act_path_coeff.lepage,
    ap->p2.act_path_coeff.naik
  );

  // -- naik preparation -- //

  int n_naiks = fermion_links_get_n_naiks(fl);
  Real* eps_naik = fermion_links_get_eps_naik(fl);
  std::vector<RealD> eps_naiks(n_naiks);
  for(int i = 0; i < n_naiks; i++) eps_naiks[i] = eps_naik[i];

  // Make orders_naik
  std::vector<int> orders_naik(n_naiks);
  int nterms = 0;
  for(int i = 0; i< n_naiks; i++) {
    orders_naik[i] = n_orders_naik[i];
    nterms += n_orders_naik[i];
  }

  // Make vecdt
  std::vector<Real> vecdt(nterms);
  for(int i = 0; i < nterms; i++)
    // Need a factor of 2 to match the MILC-code force.
    vecdt[i] = 2.*residues[i];
  
  // Make vecx
  std::vector<FermionField> vecx(nterms,CGrid);
  for(int i = 0; i < nterms; i++) {
    milcVectorFieldToGrid<FermionField, Complex>(multi_x[i], &vecx[i]);
  }

  // -- force calculation -- //

  // MILC-specific context
  MILCContext milcCtx(fatCtx, asqCtx, vecdt, eps_naiks, orders_naik);

  // Instantiate the HISQ fermion implementation class
  bool calculateStaggeredPhases = false;
  HighlyImprovedStaggeredFermionImpl<Gimpl> hisq(CGrid, calculateStaggeredPhases);
  
  // Calculate derivative
  hisq.milcSmearDerivative(UForce, Wmu, Vmu, Umu, vecx, milcCtx);
  
  gridToMilcGaugeField<LatticeGaugeField, Complex>(deriv, &UForce);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "Grid fermion force " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	          << std::endl;
}

/*
template<typename LatticeGaugeField, typename FermionField, typename Gimpl, typename Complex>
static void
hisqForce (GRID_info_t *info,
	   void *fl_void,
	   Real residues[],
	   su3_vector *multi_x[],
	   int n_orders_naik[],
	   su3_matrix *deriv,
	   GridCartesian *CGrid)
{

  fermion_links_t *fl = (fermion_links_t *)fl_void;
  
  auto start = std::chrono::system_clock::now();

  // Sort out the Gimpl. This handles BCs and part of the precision. 
  INHERIT_GIMPL_TYPES(Gimpl);
  typedef typename Gimpl::FermionField   FF;
  typedef typename Gimpl::GaugeField     GF;
  typedef typename Gimpl::GaugeLinkField LF;
  typedef typename Gimpl::ComplexField   CF;
  typedef typename Gimpl::Scalar ComplexScalar;
  typedef decltype(real(ComplexScalar())) RealScalar;
  typedef iColourMatrix<ComplexScalar> ComplexColourMatrix;

  hisq_auxiliary_t *aux = get_hisq_auxiliary(fl);
  su3_matrix *Umilc = aux->U_link;
  su3_matrix *Vmilc = aux->V_link;
  su3_matrix *Wmilc = aux->W_unitlink;

  LatticeGaugeField Umu(CGrid), Vmu(CGrid), Wmu(CGrid), UForce(CGrid);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Umilc, &Umu);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Vmilc, &Vmu);
  milcGaugeFieldToGrid<LatticeGaugeField, Complex>(Wmilc, &Wmu);

  int n_naiks = fermion_links_get_n_naiks(fl);
  Real *eps_naik = fermion_links_get_eps_naik(fl);
  std::vector<RealD> eps_naiks(n_naiks);
  for(int i = 0; i < n_naiks; i++)
    eps_naiks[i] = eps_naik[i];
  
  //  HISQParameters<Real> hisq_param = get_hisq_param(n_naiks, eps_naiks, fl);

  bool allow_svd = false, svd_only = false;
  Real svd_rel_error = HISQ_REUNIT_SVD_REL_ERROR;
  Real svd_abs_error = HISQ_REUNIT_SVD_ABS_ERROR;
  Real force_filter  = HISQ_FORCE_FILTER;

#ifdef HISQ_REUNIT_ALLOW_SVD
  allow_svd = true;
#endif

#ifdef HISQ_REUNIT_SVD_ONLY
  svd_only = true;
#endif

  // Make orders_naik
  std::vector<int> orders_naik(n_naiks);
  int nterms = 0;
  for(int i = 0; i< n_naiks; i++){
    orders_naik[i] = n_orders_naik[i];
    nterms += n_orders_naik[i];
  }

  // Make vecdt
  std::vector<Real> vecdt(nterms);
  for(int i = 0; i < nterms; i++)
    // Need a factor of 2 to match the MILC-code force.
    vecdt[i] = 2.*residues[i];
  
  // Make vecx
  std::vector<FermionField> vecx(nterms,CGrid);
  for(int i = 0; i < nterms; i++){
    milcVectorFieldToGrid<FermionField, Complex>(multi_x[i], &vecx[i]);
  }
  
  //  HISQReunitSVDParameters<Real> hisq_SVD(allow_svd, svd_only, svd_rel_error,
  //					 svd_abs_error, force_filter);

#if 0
  // Set action coefficients
  RealD pc_one_link       = path_coeff[0];
  RealD pc_naik           = path_coeff[1];
  RealD pc_three_staple   = path_coeff[2];
  RealD pc_five_staple    = path_coeff[3];
  RealD pc_seven_staple   = path_coeff[4];
  RealD pc_lepage         = path_coeff[5];
  bool backupSVD = allow_SVD;
  RealD svdTolerance = 0.;  // IS ABS OR REL ERROR??
  RealD eigenCutoff = 0.;   // WHAT ??
  
  HISFContext ctx(pc_one_link, pc_three_staple, pc_five_staple, pc_seven_staple,
		  pc_lepage, pc_naik, backupSVD, svdTolerance, eigenCutoff);
#endif

  // Instantiate the HISQ fermion implementation class
  bool calculateStaggeredPhases = true;
  HighlyImprovedStaggeredFermionImpl<Gimpl> HL(CGrid, calculateStaggeredPhases);
  
  HL.milcSmearDerivative(UForce, Wmu, Vmu, Umu, vecx, vecdt, orders_naik, eps_naiks);
  
  gridToMilcGaugeField<LatticeGaugeField, Complex>(deriv, &UForce);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "generate fat and long links " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
}
*/
	
#if 0	
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

  bool allow_svd = false, svd_only = false;
  Real svd_rel_error = HISQ_REUNIT_SVD_REL_ERROR;
  Real svd_abs_error = HISQ_REUNIT_SVD_ABS_ERROR;
  Real force_filter  = HISQ_FORCE_FILTER;

#ifdef HISQ_REUNIT_ALLOW_SVD
  allow_svd = true;
#endif

#ifdef HISQ_REUNIT_SVD_ONLY
  svd_only = true;
#endif

  // Calculate the derivative
  // We don't need hisq_param for the derivative
  std::array<Real,GRID_MAX_NAIK> eps_naiks;
  HISQParameters<Real> hisq_param(0., eps_naiks, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
  LatticeGaugeField Umu(CGrid), Vmu(CGrid), Wmu(CGrid), UForce(CGrid);
  HISQReunitSVDParameters<Real> hisq_SVD(allow_svd, svd_only, svd_rel_error,
					 svd_abs_error, force_filter);
  Force_HISQ<Gimpl> RD(CGrid, hisq_param, Wmu, Vmu, Umu, hisq_SVD);
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
#endif	

//====================================================================//
// The GRID C API for the fermion force

// Single precision is not supported
#if 0
void GRID_F3_hisq_force(GRID_info_t *info,
			void *fl,
			Real residues[],
			su3_vector *multi_x[],
			int n_orders_naik[],
			su3_matrix *deriv,
			GRID_4Dgrid *grid_full)
{
  hisqForce<LatticeGaugeFieldF, ImprovedStaggeredFermionF::FermionField, StaggeredImplF, ComplexF>(info, fl, residues,
							  multi_x, n_orders_naik,
							  deriv, grid_full->gridF);
}
#endif

void GRID_D3_hisq_force(GRID_info_t *info,
			void *fl,
			Real residues[],
			su3_vector *multi_x[],
			int n_orders_naik[],
			su3_matrix *deriv,
			GRID_4Dgrid *grid_full)
{
  hisqForce<LatticeGaugeFieldD, ImprovedStaggeredFermionD::FermionField, StaggeredImplD, ComplexD>(info, fl, residues,
							  multi_x, n_orders_naik,
							  deriv, grid_full->gridD);
}

//====================================================================//
// The GRID C API for testing the reunitarization derivative

#if 0
void GRID_F3_reunit_deriv( GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			   su3_matrix *Q, GRID_4Dgrid * grid_full ){
  //  std::cout << "GRID_F3_reunit_deriv is not supported yet" << std::endl;
  //  assert(0);
  reunitDeriv<LatticeGaugeFieldF, StaggeredImplF, ComplexF>(info, V, dW, Q, grid_full->gridF);
}

void GRID_D3_reunit_deriv( GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			   su3_matrix *Q, GRID_4Dgrid * grid_full ){
  reunitDeriv<LatticeGaugeFieldD, StaggeredImplD, ComplexD>(info, V, dW, Q, grid_full->gridD);
}

#endif
