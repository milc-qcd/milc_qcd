/******* d_bicgilu_cl_gpu.c - BiCGstab-ILU for  clover fermions on GPU's ****/
/* MIMD version 7 */

/* 3/22/12 J Foley */

#include "generic_clover_includes.h"
#include <quda_milc_interface.h>
#include "../include/generic_quda.h"
#include "../include/gammatypes.h"

/* Backward compatibility*/
#ifdef SINGLE_FOR_DOUBLE
#define HALF_MIXED
#endif

#define make_map_milc_clov_to_quda_raw(P, MILCFLOAT) \
void map_milc_clov_to_quda_raw_##P(MILCFLOAT *raw_clov, clover *milc_clov){\
  int i; \
  MILCFLOAT *r;\
\
  printf("sites_on_node = %d\n", sites_on_node);\
  for(i=0; i<sites_on_node; ++i){\
\
    r = raw_clov + 72*i;\
    r[0] = milc_clov->clov_diag[i].di[0][0]; \
    r[1] = milc_clov->clov_diag[i].di[0][1]; \
    r[2] = milc_clov->clov_diag[i].di[0][2]; \
    r[3] = milc_clov->clov_diag[i].di[0][3]; \
    r[4] = milc_clov->clov_diag[i].di[0][4]; \
    r[5] = milc_clov->clov_diag[i].di[0][5]; \
    r[6] = milc_clov->clov[i].tr[0][0].real; \
    r[7] = milc_clov->clov[i].tr[0][0].imag; \
    r[8] = milc_clov->clov[i].tr[0][1].real; \
    r[9] = milc_clov->clov[i].tr[0][1].imag; \
    r[10] = milc_clov->clov[i].tr[0][3].real; \
    r[11] = milc_clov->clov[i].tr[0][3].imag; \
    r[12] = milc_clov->clov[i].tr[0][6].real; \
    r[13] = milc_clov->clov[i].tr[0][6].imag; \
    r[14] = milc_clov->clov[i].tr[0][10].real; \
    r[15] = milc_clov->clov[i].tr[0][10].imag; \
    r[16] = milc_clov->clov[i].tr[0][2].real; \
    r[17] = milc_clov->clov[i].tr[0][2].imag; \
    r[18] = milc_clov->clov[i].tr[0][4].real; \
    r[19] = milc_clov->clov[i].tr[0][4].imag; \
    r[20] = milc_clov->clov[i].tr[0][7].real; \
    r[21] = milc_clov->clov[i].tr[0][7].imag; \
    r[22] = milc_clov->clov[i].tr[0][11].real; \
    r[23] = milc_clov->clov[i].tr[0][11].imag; \
    r[24] = milc_clov->clov[i].tr[0][5].real; \
    r[25] = milc_clov->clov[i].tr[0][5].imag; \
    r[26] = milc_clov->clov[i].tr[0][8].real; \
    r[27] = milc_clov->clov[i].tr[0][8].imag; \
    r[28] = milc_clov->clov[i].tr[0][12].real; \
    r[29] = milc_clov->clov[i].tr[0][12].imag; \
    r[30] = milc_clov->clov[i].tr[0][9].real; \
    r[31] = milc_clov->clov[i].tr[0][9].imag; \
    r[32] = milc_clov->clov[i].tr[0][13].real; \
    r[33] = milc_clov->clov[i].tr[0][13].imag; \
    r[34] = milc_clov->clov[i].tr[0][14].real; \
    r[35] = milc_clov->clov[i].tr[0][14].imag; \
\
    r[36] = milc_clov->clov_diag[i].di[1][0]; \
    r[37] = milc_clov->clov_diag[i].di[1][1]; \
    r[38] = milc_clov->clov_diag[i].di[1][2]; \
    r[39] = milc_clov->clov_diag[i].di[1][3]; \
    r[40] = milc_clov->clov_diag[i].di[1][4]; \
    r[41] = milc_clov->clov_diag[i].di[1][5]; \
    r[42] = milc_clov->clov[i].tr[1][0].real; \
    r[43] = milc_clov->clov[i].tr[1][0].imag; \
    r[44] = milc_clov->clov[i].tr[1][1].real; \
    r[45] = milc_clov->clov[i].tr[1][1].imag; \
    r[46] = milc_clov->clov[i].tr[1][3].real; \
    r[47] = milc_clov->clov[i].tr[1][3].imag; \
    r[48] = milc_clov->clov[i].tr[1][6].real; \
    r[49] = milc_clov->clov[i].tr[1][6].imag; \
    r[50] = milc_clov->clov[i].tr[1][10].real; \
    r[51] = milc_clov->clov[i].tr[1][10].imag; \
    r[52] = milc_clov->clov[i].tr[1][2].real; \
    r[53] = milc_clov->clov[i].tr[1][2].imag; \
    r[54] = milc_clov->clov[i].tr[1][4].real; \
    r[55] = milc_clov->clov[i].tr[1][4].imag; \
    r[56] = milc_clov->clov[i].tr[1][7].real; \
    r[57] = milc_clov->clov[i].tr[1][7].imag; \
    r[58] = milc_clov->clov[i].tr[1][11].real; \
    r[59] = milc_clov->clov[i].tr[1][11].imag; \
    r[60] = milc_clov->clov[i].tr[1][5].real; \
    r[61] = milc_clov->clov[i].tr[1][5].imag; \
    r[62] = milc_clov->clov[i].tr[1][8].real; \
    r[63] = milc_clov->clov[i].tr[1][8].imag; \
    r[64] = milc_clov->clov[i].tr[1][12].real; \
    r[65] = milc_clov->clov[i].tr[1][12].imag; \
    r[66] = milc_clov->clov[i].tr[1][9].real; \
    r[67] = milc_clov->clov[i].tr[1][9].imag; \
    r[68] = milc_clov->clov[i].tr[1][13].real; \
    r[69] = milc_clov->clov[i].tr[1][13].imag; \
    r[70] = milc_clov->clov[i].tr[1][14].real; \
    r[71] = milc_clov->clov[i].tr[1][14].imag; \
  } \
}
// make_map_clov_to_quda_raw

#if ( PRECISION == 1 )
make_map_milc_clov_to_quda_raw(F, float)
#define map_milc_clov_to_quda_raw map_milc_clov_to_quda_raw_F
#else // PRECISION == 2
make_map_milc_clov_to_quda_raw(D, double)
#define map_milc_clov_to_quda_raw map_milc_clov_to_quda_raw_D
#endif


int bicgilu_cl_field_gpu ( // Return value is number of iterations taken 
  wilson_vector *src,
  wilson_vector *dest,
  quark_invert_control *qic,
  void *dmp  		   // parameters defining the Dirac matrix
)			   
{
   int flag = qic->start_flag;  // O: use a zero initial guess;
			        // 1: use dest

   dirac_clover_param *dcp = (dirac_clover_param *)dmp; 
   Real kappa  = dcp->Kappa;
   Real clov_c = dcp->Clov_c;
   Real u0     = dcp->U0; 
   Real CKU0   = kappa*clov_c/(u0*u0*u0);
#ifdef CGTIME
   double dtime = -dclock();
#endif

  clover *milc_clov = gen_clov;
  if(milc_clov == NULL){
    printf("bicgilu_cl_field_gpu(%d): milc_clov == NULL\n", this_node);	
    terminate(1);
  }
  
  Real* raw_clov; Real* raw_clov_inv;
  raw_clov     = (Real*)malloc(72*sites_on_node*sizeof(Real));
  if(raw_clov == NULL){
    printf("bicgilu_cl_field_gpu(%d): no room for raw_clov\n",this_node);
    terminate(1);
  }
  raw_clov_inv = (Real*)malloc(72*sites_on_node*sizeof(Real));   
  if(raw_clov_inv == NULL){
    printf("bicgilu_cl_field_gpu(%d): no room for raw_clov_inv\n",this_node);
    terminate(2);
  }

  // solution vector from quda is returned to quda_dest 
  wilson_vector* quda_dest = (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  if(quda_dest == NULL){
    printf("bicgilu_cl_field_gpu(%d): no room for quda_dest\n",this_node);
    terminate(3);
  }

  wilson_vector* quda_src = (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
  if(quda_src == NULL){
    printf("bicgilu_cl_field_gpu(%d): no room for quda_src\n",this_node);
    terminate(4);
  }

  compute_clov(milc_clov, CKU0); // CKU0 = coefficient of the clover term
  map_milc_clov_to_quda_raw(raw_clov, milc_clov);
  
  // Note that compute_clovinv overwrites milc_clov
  compute_clovinv(milc_clov, EVENANDODD);
  map_milc_clov_to_quda_raw(raw_clov_inv, milc_clov);
  
  void* links = create_G_from_site();
  // setup QUDA
  initialize_quda();
  
  QudaInvertArgs_t inv_args;
  inv_args.max_iter          = qic->max*qic->nrestart;
  inv_args.restart_tolerance = 1e-3;
#ifdef MAX_MIXED
  inv_args.mixed_precision = 2;
#else
#ifdef HALF_MIXED
  inv_args.mixed_precision = 1;
#else
  inv_args.mixed_precision = 0;
#endif
#endif
  
  int dir;
  for(dir=0; dir<4; ++dir) inv_args.boundary_phase[dir] = boundary_phase[dir]; // boundary_phase is a global array
  const int quda_precision   = qic->prec;
  
  double residual, relative_residual;
  int num_iters = 0;
  
  qic->relresid = 0.0;
  gamma_matrix_t g5 = gamma_mat(G5);

  int i;
  FORALLFIELDSITES(i){
    mult_w_by_gamma_mat_l(&(src[i]), &(quda_src[i]), &g5);
  }

  qudaCloverInvert(PRECISION, 
		   quda_precision,
		   kappa,
		   inv_args,					
		   qic->resid,
		   qic->relresid,
		   links,
		   raw_clov,
		   raw_clov_inv,
		   quda_src, 
		   quda_dest,
		   &residual, 
		   &relative_residual,
		   &num_iters);
  qic->final_rsq    = residual*residual; 
  qic->final_relrsq = relative_residual*relative_residual;
  qic->final_iters  = num_iters; 


  // Apply gamma[5] matrix to quda_dest to get the solution in MILC format
  FORALLFIELDSITES(i){
    mult_w_by_gamma_mat_l(&(quda_dest[i]), &(dest[i]), &g5);
  }
 
   
  // check for convergence 
  qic->converged = (residual <= qic->resid) ? 1 : 0;
  // Cumulative residual. Not used in practice
  qic->size_r = 0.0;
  qic->size_relr = 0.0;
  
  destroy_G(links);
  free(quda_src);
  free(quda_dest);
  if(raw_clov) free(raw_clov);
  if(raw_clov_inv) free(raw_clov_inv);
  
#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(num_iters==0)
      printf("BICGILU: NO iterations taken size_r= %.2e rel %.2e\n",
	     qic->final_rsq, qic->final_relrsq);
#ifdef CGTIME
    else
      printf("CGTIME: time = %.2e (bicgilu GPU) size_r= %.2e relr= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->final_rsq,qic->final_relrsq,qic->final_iters,
	     (double)8742*qic->final_iters*even_sites_on_node/(dtime*1e6));
#endif
  }    
  
  return num_iters;  
} // invert_cl_field_gpu

