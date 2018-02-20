#include "qphixjClovInvert.h"
#include <omp.h>

#include <qphix/clover_dslash_def.h>
#include <qphix/clover_dslash_body.h>
#include <qphix/clover.h>

#include <qphix/invcg.h>
#include <qphix/invbicgstab.h>
#include <qphix/print_utils.h>

#include <cstdlib>
#include <assert.h>

using namespace std;
using namespace QPhiX;

extern QPHIXJ_vec_layout_t vecLayout;

template<typename FT, int V>
void
QPHIXJClovInvert::runClov(QPHIXJ_info_t *info,
			  QPHIXJ_FermionLinksWilson_struct<FT, V> *ql, 
			  QPHIXJ_invert_arg_t *inv_arg, 
			  QPHIXJ_resid_arg_t *res_arg, FT kappa,  // Unused!!
			  QPHIXJ_DiracFermion_struct<FT, V> *qdf_out, 
			  QPHIXJ_DiracFermion_struct<FT, V> *qdf_in) 
{
  typedef typename Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT,V,QPHIX_SOALEN,COMPRESS>::CloverBlock Clover;

  bool verbose = false;

  // Work out local lattice size
  const int *machsize = get_logical_dimensions();
  int subLattSize[4] = { nx/machsize[0], ny/machsize[1], nz/machsize[2], nt/machsize[3] };

  // Work out the size of checkerboarded X-dimension
  int lX1h = subLattSize[0]/2;
  int lY = subLattSize[1];
  int lZ = subLattSize[2];
  int lT = subLattSize[3];
  
  masterPrintf("Initializing Dslash\n");

  Geometry<FT,V,QPHIX_SOALEN,COMPRESS> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

  // Unpack the fermion links structure: gauge field
  Gauge *packed_gauge_cb0 = ql->packed_gauge_cb0;
  Gauge *packed_gauge_cb1 = ql->packed_gauge_cb1;
  
  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  
  // Unpack the fermion links structure: clover term
  Clover *A_cb0 = ql->A_cb0;
  Clover *A_inv_cb1 = ql->A_inv_cb1;

  // Unpack the fermion links structure: Dirac spinor source field
  QPHIXJ_evenodd_t parity_in = qdf_in->parity;
  Spinor *p_even = qdf_in->p_even;
  Spinor *p_odd  = qdf_in->p_odd;
  Spinor *psi_s[2] = { p_even, p_odd };

  // Unpack the fermion links structure: Dirac spinor solution field
  QPHIXJ_evenodd_t parity_out = qdf_out->parity;
  Spinor *c_even = qdf_out->p_even;
  Spinor *c_odd  = qdf_out->p_odd;
  Spinor *chi_s[2] = { c_even, c_odd };

  double t_boundary= 1.;
  double coeff_s = 1.;
  double coeff_t = 1.;

  if ( do_dslash ) { 
    
    ClovDslash<FT,V,QPHIX_SOALEN,COMPRESS> D32(&geom, t_boundary, coeff_s, coeff_t);

    int isign = -1;   // The sign is opposite that of the MILC code
    int cb = 1;
    //int cb = 0;
	
    int source_cb = 1 - cb;
    int target_cb = cb;
    int iters = 1;
	
    masterPrintf("Timing on cb=%d isign=%d\n", cb, isign);
    masterPrintf("=============================\n");
    
    double start = omp_get_wtime();
    
    // Apply Optimized Dslash
    D32.dslash(chi_s[target_cb],	
	       psi_s[source_cb],
	       u_packed[target_cb],
	       A_inv_cb1,
	       isign, 
	       target_cb);
	  
    double end = omp_get_wtime();
    double time = end - start;
    CommsUtils::sumDouble(&time);
    time /= (double)CommsUtils::numNodes();
    
    masterPrintf("\t %d iterations in %e seconds\n", iters, time);
    masterPrintf("\t %e usec/iteration\n", 1.0e6*time/(double)iters);
    double Gflops = 1824.0f*(double)(iters)*(double)((nx/2)*ny*ny*nt)/1.0e9;
    double perf = Gflops/time;
    masterPrintf("\t Performance: %g GFLOPS total\n", perf);
    masterPrintf("\t              %g GFLOPS / node\n", perf/(double)CommsUtils::numNodes());
    
  } // do dslash
  
  masterPrintf("Creating EvenOdd Clover Op\n");

  EvenOddCloverOperator<FT,V,QPHIX_SOALEN,COMPRESS> 
    M(u_packed, A_cb0, A_inv_cb1, &geom, t_boundary, coeff_s, coeff_t);

  if ( do_m ) { 
    
    int target_cb = 0;
    int isign = -1;
      
    masterPrintf("Timing M: isign=%d\n",  isign);
    masterPrintf("=============================\n");
    
    double start = omp_get_wtime();
    
    // Apply Optimized Dirac Matrix
    int iters = 1;
    M(chi_s[0], psi_s[0], isign, target_cb);
    
    double end = omp_get_wtime();
    double time = end - start;
    
    CommsUtils::sumDouble(&time);
    time /= (double)CommsUtils::numNodes();
    
    masterPrintf("\t %d iterations in %e seconds\n", iters, time);
    masterPrintf("\t %e usec/iteration\n", 1.0e6*time/(double)iters);
    double flops_per_iter = 3696.0f;
    double Gflops = flops_per_iter*(double)(iters)*(double)((nx/2)*ny*nz*nt)/1.0e9;
    double perf = Gflops/time;
    masterPrintf("\t Performance: %g GFLOPS total\n", perf);
    masterPrintf("\t              %g GFLOPS / node\n", perf/(double)CommsUtils::numNodes());
    
  }
  
  // Get max iters and residual target
  int max_iters = inv_arg->max;
  double rsd_target = res_arg->resid;

  int niters = 0;
  double rsd_final = 0.;
  
  masterPrintf("Creating BiCGStab Solver\n");
  
  InvBiCGStab<FT,V,QPHIX_SOALEN,COMPRESS> solver2(M, max_iters);

  masterPrintf("max_iters %d rsd_target %e\n", max_iters, rsd_target);

  if( do_bicgstab ) {

    int target_cb;

    // We don't support off-diagonal inversions and we don't support
    // evenodd for the moment
    assert(parity_in == parity_out);
    assert(parity_in == QPHIXJ_EVEN);
    if(parity_in == QPHIXJ_EVEN)target_cb = 0;
    else target_cb = 1;
    
    unsigned long site_flops;
    unsigned long mv_apps;
    int isign=-1;
    
    double start = omp_get_wtime();
    solver2(chi_s[0], psi_s[0], rsd_target, niters, rsd_final, site_flops, mv_apps, 
	    isign, target_cb, verbose);
    double end = omp_get_wtime();
    
    unsigned long num_cb_sites=lX1h*lY*lZ*lT;
    unsigned long total_flops = (site_flops + (3696)*mv_apps)*num_cb_sites;
    
    masterPrintf("Solver iters=%d\n", niters);
    masterPrintf("Solver Time=%g(s)\n", (end-start));
    masterPrintf("BICGSTAB GFLOPS/rank=%g\n", 1.0e-9*(double)(total_flops)/(end -start));
    
    info->final_sec = end-start;
    info->final_flop = total_flops;
    info->status = QPHIXJ_SUCCESS;  // Need to be more honest!!

    // Collect statistics (some of these are not supported)
    res_arg->final_rsq = rsd_final;
    res_arg->final_iter = niters;
    res_arg->final_restart = 0;
    res_arg->final_rel = 0.;
    res_arg->size_r = 0.;
    res_arg->size_relr = 0.;
  }
}


template<typename FT, int V>
void
QPHIXJClovInvert::run(QPHIXJ_info_t *info,
		      QPHIXJ_FermionLinksWilson_struct<FT, V> *ql, QPHIXJ_invert_arg_t *inv_arg,
		      QPHIXJ_resid_arg_t *res_arg, FT kappa,
		      QPHIXJ_DiracFermion_struct<FT, V> *qdf_out, 
		      QPHIXJ_DiracFermion_struct<FT, V> *qdf_in)
{
  
  if ( QPHIX_SOALEN > VECLEN_SP ) { 
    masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,VECLEN_SP);
    abort();
  }
  
  if ( compress12 ) { 
    runClov<FT,V>(info, ql, inv_arg, res_arg, kappa, qdf_out, qdf_in);
  }
  else { 
    runClov<FT,V>(info, ql, inv_arg, res_arg, kappa, qdf_out, qdf_in);
  }
}  

//====================================================================//
// The QHIXJ API for the inverter

// Debugging switches
static bool do_dslash = false;
static bool do_m = false;
static bool do_bicgstab = true;

// NOTE: kappa is unused here.  The kappa dependence is folded into the clover term.

void QPHIXJ_F3_wilson_invert(QPHIXJ_info_t *info,
			     QPHIXJ_F3_FermionLinksWilson *ql,
			     QPHIXJ_invert_arg_t *inv_arg,
			     QPHIXJ_resid_arg_t *res_arg,
			     float kappa,
			     QPHIXJ_F3_DiracFermion *qdf_out,
			     QPHIXJ_F3_DiracFermion *qdf_in)
{
  QPHIXJClovInvert CI(&vecLayout, do_dslash, do_m, do_bicgstab);
  CI.run<float, VECLEN_SP>(info, ql, inv_arg, res_arg, kappa, 
			   qdf_out, qdf_in);
}

void QPHIXJ_D3_wilson_invert(QPHIXJ_info_t *info,
			     QPHIXJ_D3_FermionLinksWilson *ql,
			     QPHIXJ_invert_arg_t *inv_arg,
			     QPHIXJ_resid_arg_t *res_arg,
			     double kappa,
			     QPHIXJ_D3_DiracFermion *qdf_out,
			     QPHIXJ_D3_DiracFermion *qdf_in)
{
  QPHIXJClovInvert CI(&vecLayout, do_dslash, do_m, do_bicgstab);
  CI.run<double, VECLEN_DP>(info, ql, inv_arg, res_arg, kappa, 
			    qdf_out, qdf_in);
}



