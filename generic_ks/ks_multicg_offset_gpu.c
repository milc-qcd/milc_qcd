/************** ks_multicg_offset_gpu.c **************************/
/* MIMD version 7 */
#include <test_util.h>
#include <blas_reference.h>
#include <staggered_dslash_reference.h>
#include <quda.h>
#include <string.h>
#include <gauge_quda.h>

#include "generic_ks_includes.h"        /* definitions files and prototypes */
#include "../include/dslash_ks_redefine.h"

#include "../include/loopend.h"

extern FullGauge cudaFatLinkPrecise;
extern FullGauge cudaFatLinkSloppy;
extern FullGauge cudaLongLinkPrecise;
extern FullGauge cudaLongLinkSloppy;
#define mySpinorSiteSize 6



/* Interface for call with offsets = 4 * mass * mass */
int ks_multicg_offset_gpu(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    Real *_offsets,	/* the offsets */
    int num_offsets,	/* number of offsets */
    int niter,		/* maximal number of CG interations */
    Real rsqmin,	/* desired residue squared */
    int prec,           /* internal precision for inversion (ignored) */
    int parity,		/* parity to be worked on */
    Real *final_rsq_ptr,/* final residue squared */
    imp_ferm_links_t **fn /* Storage for fat and Naik links */
    )
{
    
    //return ks_multicg_offset(src, psim, _offsets, num_offsets, niter, rsqmin, prec, parity, final_rsq_ptr, fn);
    int i,dir;
    QudaGaugeParam gaugeParam;
    QudaInvertParam inv_param;
    
#ifdef CGTIME
    static const char *milc_prec[2] = {"F", "D"};
    double dtimec = -dclock(); 
#endif
    int nflop = 1205 + 15*num_offsets;
    if(parity==EVENANDODD)nflop *=2;

#ifdef GPU_VERBOSE
    inv_param.verbosity = QUDA_VERBOSE;
#else
    inv_param.verbosity = QUDA_SILENT;
#endif

#ifdef GPU_USE_12_RECON
    QudaReconstructType link_recon = QUDA_RECONSTRUCT_12;
#else
    QudaReconstructType link_recon = QUDA_RECONSTRUCT_NO;
#endif

#if (PRECISION==1)
    QudaPrecision spinor_prec = QUDA_SINGLE_PRECISION;
    QudaPrecision  link_prec = QUDA_SINGLE_PRECISION;
    QudaPrecision  cpu_prec = QUDA_SINGLE_PRECISION;
    QudaReconstructType link_recon_sloppy =  link_recon;
    QudaPrecision spinor_prec_sloppy = spinor_prec;
    QudaPrecision  link_prec_sloppy = link_prec;
#else
    QudaPrecision spinor_prec = QUDA_DOUBLE_PRECISION;
    QudaPrecision  link_prec = QUDA_DOUBLE_PRECISION;
    QudaPrecision  cpu_prec = QUDA_DOUBLE_PRECISION;
    QudaReconstructType link_recon_sloppy =  link_recon;
    QudaPrecision spinor_prec_sloppy = spinor_prec;
    QudaPrecision  link_prec_sloppy = link_prec;

#endif

    su3_matrix *fatlink[4];
    su3_matrix *longlink[4];
    
    int V = nx*ny*nz*nt;
    int Vh = V/2;
    
 
    gaugeParam.X[0] = nx;
    gaugeParam.X[1] = ny;
    gaugeParam.X[2] = nz;
    gaugeParam.X[3] = nt;

    gaugeParam.cpu_prec = cpu_prec;

    gaugeParam.cuda_prec = link_prec;
    gaugeParam.reconstruct = link_recon;

    gaugeParam.cuda_prec_sloppy = link_prec_sloppy;
    gaugeParam.reconstruct_sloppy = link_recon_sloppy;

    gaugeParam.anisotropy = u0;

    inv_param.inv_type = QUDA_CG_INVERTER;

    gaugeParam.t_boundary = QUDA_ANTI_PERIODIC_T;
    gaugeParam.gauge_order = QUDA_QDP_GAUGE_ORDER;

    inv_param.tol = sqrt(rsqmin);
    inv_param.maxiter = niter;
    inv_param.reliable_delta = 1e-3;

    inv_param.mass_normalization = QUDA_MASS_NORMALIZATION;
    inv_param.cpu_prec = cpu_prec;
    inv_param.cuda_prec = spinor_prec;
    inv_param.cuda_prec_sloppy = spinor_prec_sloppy;
    inv_param.solution_type = QUDA_MATDAG_MAT_SOLUTION;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
    inv_param.dirac_order = QUDA_DIRAC_ORDER;
    inv_param.dslash_type = QUDA_ASQTAD_DSLASH;
    gaugeParam.ga_pad = nx*ny*nz;
    inv_param.sp_pad = nx*ny*nz;
    inv_param.cl_pad = nx*ny*nz;

    size_t gSize = (gaugeParam.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
    size_t sSize = (inv_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
    
    for (dir = 0; dir < 4; dir++) {
        fatlink[dir] = (su3_matrix*)malloc(V*gaugeSiteSize*gSize);
        longlink[dir] = (su3_matrix*)malloc(V*gaugeSiteSize*gSize);
	if (fatlink[dir] == NULL || longlink[dir]== NULL){
	  fprintf(stderr, "ERROR: malloc failed for fatlink and longlink\n");
	  exit(1);
	}
	
    }

    for(dir = 0; dir < 4; dir++) {
        for(i=0;i < V; i++){
	  fatlink[dir][i] = fn->fat[4*i+dir];
	  longlink[dir][i] = fn->lng[4*i+dir];
        }
    }
    
    
#if 0
    for (i =0;i < 4 ;i++){
      int dir = 2*i;
      link_sanity_check(longlink[i], V, gaugeParam.cpu_prec, dir, &gaugeParam);
    }
#endif
    
    void *spinorIn;
    su3_vector* src_vec;
    
    src_vec= (su3_vector*)malloc(V*mySpinorSiteSize*sSize);
    if (src_vec == NULL){
      fprintf(stderr, "Error: malloc failed for src in multimass CG\n");
      exit(1);
    }

    if (parity == ODD){
      spinorIn = ((char*)src_vec) + Vh*mySpinorSiteSize*sSize;
    }else{
      spinorIn= (void*)src_vec;
    }
    
    site* s;
    int j=0;
    FORSOMEPARITY(i,s,parity){	
	su3vec_copy((su3_vector *)F_PT(s,src), (su3_vector*)(((char*)spinorIn)+j*mySpinorSiteSize*sSize));
	j++;
    }END_LOOP
    

    void *spinorOut[num_offsets];
    for(i=0;i < num_offsets;i++){
	if (parity == ODD){
	    spinorOut[i]= (void*)(psim[i] + Vh);
	}else{
	    spinorOut[i]= (void*)psim[i];
	}
    }

    int device=0;    
    initQuda(device);

    gaugeParam.type = QUDA_ASQTAD_FAT_LINKS;
    gaugeParam.reconstruct = gaugeParam.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
    loadGaugeQuda(fatlink, &gaugeParam);

    gaugeParam.type = QUDA_ASQTAD_LONG_LINKS;
    gaugeParam.reconstruct= link_recon;
    gaugeParam.reconstruct_sloppy = link_recon_sloppy;
    loadGaugeQuda(longlink, &gaugeParam);

    switch(parity){
    case EVEN:
      inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
      inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
      break;
    case ODD:
      inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
      inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
      break;
    default:
      fprintf(stderr, "ERROR: invalid parity\n");
      exit(1);
      break;
    }
    
    double offsets[num_offsets];
    for(i=0;i < num_offsets;i++){
	offsets[i] = _offsets[i];
    }
    double rsd_sq;
    invertMultiShiftQuda(spinorOut, spinorIn, &inv_param, offsets, num_offsets, &rsd_sq);
    int iters = inv_param.iter;
    *final_rsq_ptr =rsd_sq;
    
#ifdef GPU_CG_CHECK
    {   
	printf("parity =%d , final residue square is =%e\n", parity,  *final_rsq_ptr);
	su3_vector* psim_check[num_offsets];
	
	for(i=0;i < num_offsets; i++){
	    psim_check[i] = (su3_vector*)malloc(V*mySpinorSiteSize*sSize);
	    if (psim_check[i] == NULL){
		fprintf(stderr, "ERROR: malloc failed for psim_check\n");
		exit(1);	    
	    }
	}
	
	int cpu_iters = ks_multicg_offset_cpu(src, psim_check, _offsets, num_offsets, niter, rsqmin, prec, parity, final_rsq_ptr, fn);
	printf("GPU iterations=%d, CPU iterations =%d\n", iters, cpu_iters);
#define DELTA 0.0001
	for(i=0;i < num_offsets;i++){
	    printf("matching %dth solution\n", i);
	    FORSOMEPARITY(j,s,parity){	
		if ( fabs(psim[i][j].c[0].real - psim_check[i][j].c[0].real) > DELTA
		     || fabs(psim[i][j].c[0].imag - psim_check[i][j].c[0].imag) > DELTA
		     || fabs(psim[i][j].c[1].real - psim_check[i][j].c[1].real) > DELTA
		     || fabs(psim[i][j].c[1].imag - psim_check[i][j].c[1].imag) > DELTA
		     || fabs(psim[i][j].c[2].real - psim_check[i][j].c[2].real) > DELTA
		     || fabs(psim[i][j].c[2].imag - psim_check[i][j].c[2].imag) > DELTA){
		    fprintf(stderr, "Error: %dth mass solution does not match in j=%d\n", i, j);
		    printf("from gpu: psim[%d][%d]      =(%f,%f), (%f,%f), (%f,%f)\n", 
			   i, j, 
			   psim[i][j].c[0].real, psim[i][j].c[0].imag, 
			   psim[i][j].c[1].real, psim[i][j].c[1].imag, 
			   psim[i][j].c[2].real, psim[i][j].c[2].imag);
		    printf("from cpu: psim_check[%d][%d]=(%f,%f), (%f,%f), (%f,%f)\n",
			   i, j, 
			   psim_check[i][j].c[0].real, psim_check[i][j].c[0].imag, 
			   psim_check[i][j].c[1].real, psim_check[i][j].c[1].imag, 
			   psim_check[i][j].c[2].real, psim_check[i][j].c[2].imag);
			   
		    exit(1);
		}
	    }END_LOOP			
	}
	

	//solution residue checking
	int other_parity = ODD;
	if (parity == ODD){
	  other_parity = EVEN;
	}	
	int len = V;
	
	void* spinorCheck = malloc(len*mySpinorSiteSize*sSize);
	if (spinorCheck == NULL){
	  printf("ERROR: malloc failed for spinorCheck\n");
	  exit(1);
	}
	printf("checking the solution\n");

	if (parity ==EVENANDODD){
	  printf("ERROR: Checking for evenodd parity is not suppported yet\n");
	  exit(1);
	}
	//compute source norm
	double src_rsq = 0;
	site* s;
	FORSOMEPARITYDOMAIN(i,s,parity){
	  src_rsq += (double)magsq_su3vec( &src_vec[i]);
	} END_LOOP

	for(i=0;i < num_offsets;i++){
	    printf("GPU %dth solution: ", i);
	    double gpu_rsq = 0;
	    dslash_fn_field(psim[i], spinorCheck, other_parity, fn);
	    dslash_fn_field(spinorCheck, spinorCheck, parity, fn);
	    su3_vector* x = (su3_vector*)psim[i];
	    su3_vector* DDx = (su3_vector*)spinorCheck;
	    su3_vector mAx; //minus Ax
	    su3_vector resid;
	    FORSOMEPARITYDOMAIN(j,s,parity){
	      scalar_mult_add_su3_vector( &DDx[j], &x[j], -offsets[i], &mAx );
	      add_su3_vector( &mAx, &src_vec[j], &resid );      
	      gpu_rsq += (double)magsq_su3vec( &resid);
	    } END_LOOP
	    printf("GPU relative residual, requested = %g, actual = %g\n", inv_param.tol, sqrt(gpu_rsq/src_rsq));
	}
		
	for(i=0;i < num_offsets;i++){
	    printf("CPU %dth solution: ", i);
	    double cpu_rsq = 0;
	
	    dslash_fn_field(psim_check[i], spinorCheck, other_parity, fn);
	    dslash_fn_field(spinorCheck, spinorCheck, parity, fn);
	    su3_vector* x = (su3_vector*)psim_check[i];
	    su3_vector* DDx = (su3_vector*)spinorCheck;
	    su3_vector mAx; //minus Ax
	    su3_vector resid;
	    FORSOMEPARITYDOMAIN(j,s,parity){
	      scalar_mult_add_su3_vector( &DDx[j], &x[j], -offsets[i], &mAx );
	      add_su3_vector( &mAx, &src_vec[j], &resid );      
	      cpu_rsq += (double)magsq_su3vec( &resid);
	    } END_LOOP
	    printf("CPU relative residual, requested = %g, actual = %g\n", inv_param.tol, sqrt(cpu_rsq/src_rsq));
	}	
	
	free(spinorCheck);	
	for(i=0;i < num_offsets; i++){
	  free(psim_check[i]);
	}
	
	
    }

#endif


#ifdef CGTIME
	    dtimec += dclock();
	    if(this_node==0){
	      printf("CONGRAD5: time = %e (multicg_offset_gpu %s) masses = %d iters = %d mflops = %e\n",
		     dtimec,milc_prec[PRECISION-1],num_offsets,iters,
		     (double)(nflop)*volume*
		     iters/(1.0e6*dtimec*numnodes()));
		fflush(stdout);}
#endif

   total_iters += iters;

    for (dir = 0; dir < 4; dir++) {
        free(fatlink[dir]);
        free(longlink[dir]);
    }
    free(src_vec);

    freeGaugeField(&cudaFatLinkPrecise);
    freeGaugeField(&cudaLongLinkPrecise);
    if (gaugeParam.cuda_prec_sloppy != gaugeParam.cuda_prec){
      freeGaugeField(&cudaFatLinkSloppy);
    }
    if ( (gaugeParam.cuda_prec_sloppy != gaugeParam.cuda_prec)
	 || (link_recon != link_recon_sloppy)){
      freeGaugeField(&cudaLongLinkSloppy);
    }

    return iters;
    
}
