// Wrappers for Grid Staggered Link Fattening

#include <omp.h>

#include "../include/macros.h"
#include "../include/fermion_links.h"

#include "../include/mGrid/mGrid_internal.h"
#include "../include/mGrid/mGrid.h"
#include "../include/milc_datatypes.h"
#include "../include/mGrid/mGrid_assert.h"

#include "../generic/gridMap.h"
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/HISQSmearing.h>

using namespace Grid;

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
	   void *fl_void,
	   Real residues[],
	   su3_vector *multi_x[],
	   int n_orders_naik[],
	   su3_matrix *deriv[],
	   GridCartesian *CGrid)
{
  fermion_links_t *fl = (fermion_links_t *)fl_void;
  
  auto start = std::chrono::system_clock::now();

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
  std::array<Real,GRID_MAX_NAIK> eps_naiks;
  for(int i = 0; i < n_naiks; i++)
    eps_naiks[i] = eps_naik[i];
  
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
  
  HISQParameters<Real> hisq_param(n_naiks  , eps_naiks ,
	  fat7_c1  , fat7_c3  , fat7_c5  , fat7_c7  , 0.,
	  asqtad_c1, asqtad_c3, asqtad_c5, asqtad_c7, asqtad_clp,
	  cnaik    , diff_c1     , diff_cnaik);

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
  
  HISQReunitSVDParameters<Real> hisq_SVD(allow_svd, svd_only, svd_rel_error,
					 svd_abs_error, force_filter);
    
  Force_HISQ<Gimpl> HF(CGrid, hisq_param, Wmu, Vmu, Umu, hisq_SVD);

  HF.ddVprojectU3(UForce, Umu, Umu, 5e-5);

  gridToMilcGaugeField<LatticeGaugeField, Complex>(deriv, &UForce);

  auto end = std::chrono::system_clock::now();
  auto elapsed = end - start;
  std::cout << "generate fat and long links " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed) 
	    << "\n";
}
	
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
  //  std::cout << "GRID_F3_hisq_links is not supported yet" << std::endl;
  //  assert(0);
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

//====================================================================//
// The GRID C API for the fermion force

void GRID_F3_hisq_force(GRID_info_t *info,
			void *fl,
			Real residues[],
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
			void *fl,
			Real residues[],
			su3_vector *multi_x[],
			int n_orders_naik[],
			su3_matrix *deriv[],
			GRID_4Dgrid *grid_full)
{
  hisqForce<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, fl, residues,
							  multi_x, n_orders_naik,
							  deriv, grid_full->gridF);
}

//====================================================================//
// The GRID C API for testing the reunitarization derivative

void GRID_F3_reunit_deriv( GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			   su3_matrix *Q, GRID_4Dgrid * grid_full ){
  //  std::cout << "GRID_F3_reunit_deriv is not supported yet" << std::endl;
  //  assert(0);
  reunitDeriv<LatticeGaugeFieldF, PeriodicGimplF, ComplexF>(info, V, dW, Q, grid_full->gridF);
}

void GRID_D3_reunit_deriv( GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			   su3_matrix *Q, GRID_4Dgrid * grid_full ){
  reunitDeriv<LatticeGaugeFieldD, PeriodicGimplD, ComplexD>(info, V, dW, Q, grid_full->gridD);
}

