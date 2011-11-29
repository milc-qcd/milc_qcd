/******* d_congrad5_fn_gpu.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* GPU version of d_congrad5_fn_milc.c.  Can be compiled together with it. */

/* The following headers are supplied with QUDA */
#include <test_util.h>
#include <blas_reference.h>
#include <quda.h>
#include <string.h>
#include <gauge_quda.h>
#include <mpicomm.h>

/* The following headers are supplied with the MILC code */
#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1


#define LOOPEND
#include "../include/loopend.h"


extern FullGauge cudaFatLinkPrecise;
extern FullGauge cudaFatLinkSloppy;
extern FullGauge cudaLongLinkPrecise;
extern FullGauge cudaLongLinkSloppy;

#define MAX(a,b) ((a)>(b)?(a):(b))

#define mySpinorSiteSize 6
extern int numnodes(void);
int V;
int Vh;

#if 0

static void
display_link_internal(float* link)
{
    int i, j;
    
    for (i = 0;i < 3; i++){
	for(j=0;j < 3; j++){
	    printf("(%.10f,%.10f) \t", link[i*3*2 + j*2], link[i*3*2 + j*2 + 1]);
	}
	printf("\n");
    }
    printf("\n");
    return;
}

static void 
accumulateConjugateProduct(float *a, float *b, float *c, int sign) 
{
    a[0] += sign * (b[0]*c[0] - b[1]*c[1]);
    a[1] -= sign * (b[0]*c[1] + b[1]*c[0]);
}

static int
link_sanity_check_internal(float* link, int dir, int ga_idx, QudaGaugeParam* gaugeParam, int oddBit)
{
    //printf("link sanity check is called\n");
    
    int ret =0;
    
    float refc_buf[6];
    float* refc = &refc_buf[0];

    memset((void*)refc, 0, sizeof(refc_buf));

    float* a = link;
    float* b = link + 6;
    float* c = link + 12;
    
    accumulateConjugateProduct(refc + 0*2, a + 1*2, b + 2*2, +1);
    accumulateConjugateProduct(refc + 0*2, a + 2*2, b + 1*2, -1);
    accumulateConjugateProduct(refc + 1*2, a + 2*2, b + 0*2, +1);
    accumulateConjugateProduct(refc + 1*2, a + 0*2, b + 2*2, -1);
    accumulateConjugateProduct(refc + 2*2, a + 0*2, b + 1*2, +1);
    accumulateConjugateProduct(refc + 2*2, a + 1*2, b + 0*2, -1);
    
    int X1h=gaugeParam->X[0]/2;
    int X1 =gaugeParam->X[0];    
    int X2 =gaugeParam->X[1];
    int X3 =gaugeParam->X[2];
    int X4 =gaugeParam->X[3];
    double t_boundary = (gaugeParam->t_boundary ==QUDA_ANTI_PERIODIC_T)? -1.0:1.0;

   double u0 = gaugeParam->anisotropy;
   double coff= -u0*u0*24;
   //coff = (dir < 6) ? coff : ( (ga_idx >= (X4-3)*X1h*X2*X3 )? t_boundary : 1); 

   //float u0 = (dir < 6) ? gaugeParam->anisotropy : ( (ga_idx >= (X4-3)*X1h*X2*X3 )? t_boundary : 1); 

#if 1  

   {
       int index = fullLatticeIndex(ga_idx, oddBit);
       int i4 = index /(X3*X2*X1);
       int i3 = (index - i4*(X3*X2*X1))/(X2*X1);
       int i2 = (index - i4*(X3*X2*X1) - i3*(X2*X1))/X1;
       int i1 = index - i4*(X3*X2*X1) - i3*(X2*X1) - i2*X1;
       
       int coords[4];
       get_coords(coords, 0, ga_idx + oddBit*Vh);
       
    

       if (dir == 0) {
	   if (i4 % 2 == 1){
	       coff *= -1;
	   }
       }
       
       if (dir == 2){
	   if ((i4+i1) % 2 == 1){
	       coff *= -1;
	   }
       }
       if (dir == 4){
	   if ( (i4+i1+i2) % 2 == 1){
	       coff *= -1;
	   }
       }       
       if (dir == 6){	   
	   if (ga_idx >= (X4-3)*X1h*X2*X3 ){
	       coff *= -1;
	   }
       }     

       printf("local ga_idx =%d, index=%d, i4,3,2,1 =%d %d %d %d, coords=%d %d %d %d\n", 
	      ga_idx, index, i4, i3, i2,i1, coords[3], coords[2], coords[1], coords[0]);       
       
   }
#endif

   refc[0]*=coff; refc[1]*=coff; refc[2]*=coff; refc[3]*=coff; refc[4]*=coff; refc[5]*=coff;
   
    
    double delta = 0.0001;
    int i;
    for (i =0;i < 6; i++){
	double diff =  refc[i] -  c[i];
	double absdiff = diff > 0? diff: (-diff);
	if (absdiff  > delta){
	    printf("ERROR: sanity check failed for link\n");
	    display_link_internal(link);
	    printf("refc = (%.10f,%.10f) (%.10f,%.10f) (%.10f,%.10f)\n", 
		   refc[0], refc[1], refc[2], refc[3], refc[4], refc[5]);
	    printf("dir=%d, ga_idx=%d, coff=%f, t_boundary=%f\n",dir, ga_idx,coff, t_boundary);
	    printf("X=%d %d %d %d, X1h=%d\n", gaugeParam->X[0], X2, X3, X4, X1h);
	    return -1;
	}
	
    }
    

    return ret;
}


//this len must be V
int
link_sanity_check(void* link, int len, int precision, int dir, QudaGaugeParam* gaugeParam)
{
    int i;
    int rc = 0;
    
    if (precision == QUDA_DOUBLE_PRECISION){
	printf("NOT SUPPORTED! ERROR\n");
	exit(1);	
    }else if (precision == QUDA_SINGLE_PRECISION){
	float* mylink = (float*)link;
	
	//even
	for (i=0;i < len/2 ;i++){
	    rc = link_sanity_check_internal(mylink + gaugeSiteSize*i, dir, i, gaugeParam, 0);
	    if (rc != 0){
		printf("ERROR: even link sanity check failed, i=%d\n",i);
		display_link_internal(mylink+gaugeSiteSize*i);
		exit(1);
	    }		
	}

	//odd
	mylink = mylink + gaugeSiteSize*Vh;
	for (i=0;i < len/2 ;i++){
	    rc = link_sanity_check_internal(mylink + gaugeSiteSize*i, dir, i, gaugeParam,1);
	    if (rc != 0){
		printf("ERROR: odd link sanity check failed, i=%d\n", i);
		display_link_internal(mylink+gaugeSiteSize*i);
		exit(1);
	    }		
	}	

    }
    
    return rc;
}
#endif

extern int
ks_congrad_parity_cpu( su3_vector *t_src, su3_vector *t_dest, 
		       quark_invert_control *qic, Real mass,
		       imp_ferm_links_t *fn);


static void
set_params(QudaGaugeParam* gaugeParam, QudaInvertParam* inv_param,
	   int X1, int  X2, int X3, int X4,
	   QudaPrecision cpu_prec, QudaPrecision prec, QudaPrecision prec_sloppy,
	   QudaReconstructType link_recon, QudaReconstructType link_recon_sloppy,
	   double mass, double tol, int maxiter, double reliable_delta,
	   double tadpole_coeff)
{
  int i;
  
#ifdef GPU_VERBOSE
  inv_param->verbosity = QUDA_VERBOSE;
#elif defined(GPU_SUMMARIZE)
  inv_param->verbosity = QUDA_SUMMARIZE;
#else
  inv_param->verbosity = QUDA_SILENT;
#endif

  gaugeParam->X[0] = X1;
  gaugeParam->X[1] = X2;
  gaugeParam->X[2] = X3;
  gaugeParam->X[3] = X4;

  gaugeParam->cpu_prec = cpu_prec;    
  gaugeParam->cuda_prec = prec;
  gaugeParam->reconstruct = link_recon;  
  gaugeParam->cuda_prec_sloppy = prec_sloppy;
  gaugeParam->reconstruct_sloppy = link_recon_sloppy;
  gaugeParam->gauge_fix = QUDA_GAUGE_FIXED_NO;
  gaugeParam->tadpole_coeff = tadpole_coeff;
  gaugeParam->t_boundary = QUDA_ANTI_PERIODIC_T;
  gaugeParam->gauge_order = QUDA_QDP_GAUGE_ORDER;
  gaugeParam->ga_pad = X1*X2*X3/2;

  inv_param->inv_type = QUDA_CG_INVERTER;    
  inv_param->mass = mass;
  inv_param->tol = tol;
  inv_param->maxiter = maxiter;
  inv_param->reliable_delta = 1e-1;

  inv_param->solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
  inv_param->solve_type = QUDA_NORMEQ_PC_SOLVE;
  inv_param->matpc_type = QUDA_MATPC_EVEN_EVEN;
  inv_param->dagger = QUDA_DAG_NO;
  inv_param->mass_normalization = QUDA_MASS_NORMALIZATION;

  inv_param->cpu_prec = cpu_prec;
  inv_param->cuda_prec = prec; 
  inv_param->cuda_prec_sloppy = prec_sloppy;
  inv_param->preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param->dirac_order = QUDA_DIRAC_ORDER;
  inv_param->dslash_type = QUDA_ASQTAD_DSLASH;
  inv_param->dirac_tune = QUDA_TUNE_NO;
  inv_param->preserve_dirac = QUDA_PRESERVE_DIRAC_NO;
  inv_param->sp_pad = X1*X2*X3/2;
  inv_param->use_init_guess = QUDA_USE_INIT_GUESS_YES;
}

  
  
  
  
int 
ks_congrad_parity_gpu(su3_vector *t_src, su3_vector *t_dest, 
    quark_invert_control *qic, Real mass,
    imp_ferm_links_t *fn)
{

  node0_printf("In ks_congrad_parity_gpu\n");

  switch(qic->parity){
    case EVEN:
    case ODD:
      break;
    default:
      fprintf(stderr, "ERROR: invalid parity\n");
      exit(1);
      break;
  }


  int i, dir;
  su3_matrix *t_fatlink = get_fatlinks(fn); 
  su3_matrix *t_longlink = get_lnglinks(fn);


  QudaGaugeParam gaugeParam = newQudaGaugeParam();
  QudaInvertParam inv_param = newQudaInvertParam();


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
#else
  QudaPrecision spinor_prec = QUDA_DOUBLE_PRECISION;
  QudaPrecision  link_prec = QUDA_DOUBLE_PRECISION;
  QudaPrecision  cpu_prec = QUDA_DOUBLE_PRECISION;    
  QudaReconstructType link_recon_sloppy =  link_recon;
#endif


#ifdef MIXED_PREC_CG_GPU 
  QudaPrecision spinor_prec_sloppy = QUDA_SINGLE_PRECISION;
  QudaPrecision  link_prec_sloppy =  QUDA_SINGLE_PRECISION;
#else
  QudaPrecision spinor_prec_sloppy = spinor_prec;
  QudaPrecision  link_prec_sloppy = link_prec;
#endif 


  su3_matrix *fatlink[4];  // links have to be passed to loadGaugeQuda 
  su3_matrix *longlink[4]; // loadGaugeQuda as two-dimensional arrays
  int X[4]; 

  int even_odd_change = 0;
  //initialize communication 
  {
    // why initialize it here and then immediately reassign below?
    int gridsize[4]={1,1,1,numnodes()};

    //hyper_prime
    const int *nsquares = get_logical_dimensions();

    for(i=0;i < 4; i++){
      gridsize[i]= nsquares[i];
    }   

    int total_num_processes = 1;
    int d;
    for (d=0; d<4; d++) total_num_processes *= gridsize[d];
    int size = -1;
    size = numnodes();
    if (total_num_processes != size)
      errorQuda("Number of processes %d must match requested MPI volume %d",
          size, total_num_processes);

    comm_set_gridsize(gridsize[0], gridsize[1], gridsize[2], gridsize[3]);  
    comm_init();


    X[0] = gaugeParam.X[0] = nx/gridsize[0];
    X[1] = gaugeParam.X[1] = ny/gridsize[1];
    X[2] = gaugeParam.X[2] = nz/gridsize[2];
    X[3] = gaugeParam.X[3] = nt/gridsize[3];

    for(i=0;i < 4;i++){
      if(X[i] % 2 == 1 && get_logical_coordinate()[i] % 2 ==1 ){
        even_odd_change = 1 - even_odd_change;
      } 
    }  
  } // end initialization of communication


  V = X[0]*X[1]*X[2]*X[3];
  Vh = V/2;
  set_params(&gaugeParam, &inv_param,
      X[0], X[1], X[2],X[3],
      cpu_prec, link_prec, link_prec_sloppy,
      link_recon, link_recon_sloppy, mass, sqrt(qic->resid), qic->max*qic->nrestart, 1e-3,
      u0);

  int device = 0; //device is not used in initQuda in MPI_COMMS case
  initQuda(device);

  size_t gSize = (gaugeParam.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
  size_t sSize = (inv_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);

  for (dir = 0; dir < 4; dir++) {
    fatlink[dir] = (su3_matrix*)malloc(V*gaugeSiteSize*gSize);
    longlink[dir] = (su3_matrix*)malloc(V*gaugeSiteSize*gSize);
    memset(fatlink[dir], 0, V*gaugeSiteSize*gSize);
    memset(longlink[dir], 0, V*gaugeSiteSize*gSize);
  }

  for(dir = 0; dir < 4; dir++) {
    if (even_odd_change){
      for(i=0;i < Vh; i++){
        fatlink[dir][Vh+i] = t_fatlink[4*i+dir];
        longlink[dir][Vh+i] = t_longlink[4*i+dir];
      }

      for(i=Vh;i < V; i++){
        fatlink[dir][i-Vh] = t_fatlink[4*i+dir];
        longlink[dir][i-Vh] = t_longlink[4*i+dir];
      }

    }else{
      for(i=0;i < Vh; i++){
        fatlink[dir][i] = t_fatlink[4*i+dir];
        longlink[dir][i] = t_longlink[4*i+dir];
      }

      for(i=Vh;i < V; i++){
        fatlink[dir][i] = t_fatlink[4*i+dir];
        longlink[dir][i] = t_longlink[4*i+dir];
      }
    }
  }


  memset(&cudaFatLinkPrecise, 0, sizeof(FullGauge));         
  memset(&cudaFatLinkSloppy, 0, sizeof(FullGauge));         
  memset(&cudaLongLinkPrecise, 0, sizeof(FullGauge));         
  memset(&cudaLongLinkSloppy, 0, sizeof(FullGauge));      
  int tmp = MAX(X[1]*X[2]*X[3]/2, X[0]*X[2]*X[3]/2);
  tmp = MAX(tmp, X[0]*X[1]*X[3]/2);
  tmp = MAX(tmp, X[0]*X[1]*X[2]/2);
  int fat_pad = tmp;
  int long_pad =  3*tmp;

  gaugeParam.type = QUDA_ASQTAD_FAT_LINKS;
  gaugeParam.ga_pad = fat_pad; 
  gaugeParam.reconstruct = gaugeParam.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
  loadGaugeQuda(fatlink, &gaugeParam);    

  gaugeParam.type = QUDA_ASQTAD_LONG_LINKS;
  gaugeParam.ga_pad = long_pad;;
  gaugeParam.reconstruct= link_recon;
  gaugeParam.reconstruct_sloppy = link_recon_sloppy;
  loadGaugeQuda(longlink, &gaugeParam);



  void *spinorIn;
  void *spinorOut;

  switch(qic->parity){
    case EVEN:
      spinorIn = t_src;
      spinorOut= t_dest;
      break;
    case ODD:
      spinorIn = ((char*)t_src) + Vh*mySpinorSiteSize*sSize;
      spinorOut= ((char*)t_dest)+ Vh*mySpinorSiteSize*sSize;
      break;
    default:
      fprintf(stderr, "ERROR: invalid parity\n");
      exit(1);
  }

  QudaParity local_parity;
  if (even_odd_change){
    local_parity = (qic->parity == EVEN)? QUDA_ODD_PARITY: QUDA_EVEN_PARITY;
  }else{
    local_parity= (qic->parity == EVEN)? QUDA_EVEN_PARITY: QUDA_ODD_PARITY;
  }

  if(local_parity ==  QUDA_EVEN_PARITY){
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;      
  }else{
    inv_param.matpc_type = QUDA_MATPC_ODD_ODD;      
  }

  invertQuda(spinorOut, spinorIn, &inv_param);

  int iters = inv_param.iter;

#ifdef GPU_CG_CHECK
  {

#define DELTA 0.01
    su3_vector* t_dest_check = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
    if (t_dest_check == NULL){
      printf("ERROR: malloc faild for t_dest_check\n");
      exit(1);
    }
    memset(t_dest_check, 0, sizeof(su3_vector)*sites_on_node);
    int cpu_iters=0;
    if(qic->parity == EVEN || qic->parity == ODD){
      cpu_iters += ks_congrad_parity_cpu(t_src, t_dest_check, qic, mass, fn);    
    }else{
      printf("ERROR: inalid parity(%d)\n", qic->parity);
    }
    printf("GPU iterations=%d, CPU iters =%d\n", iters, cpu_iters);


    void *spinorCheck = malloc(V*mySpinorSiteSize*sSize);

    int parity = EVEN;
    int other_parity = ODD;
    if (qic->parity == ODD){
      parity = ODD;
      other_parity = EVEN;
    }


    //solution residue checking

    //compute source norm
    double gpu_rsq = 0;
    double cpu_rsq = 0;
    double src_rsq = 0;
    site* s;
    FORSOMEPARITYDOMAIN(i,s,parity){
      src_rsq += (double)magsq_su3vec( &t_src[i]);
    } END_LOOP
    g_doublesum(&src_rsq);
    //compute residue for GPU solution
    dslash_fn_field(t_dest, spinorCheck, other_parity, fn);
    dslash_fn_field(spinorCheck, spinorCheck, parity, fn);
    FORSOMEPARITYDOMAIN(i,s,parity){
      su3_vector* DDx = (su3_vector*)spinorCheck;
      su3_vector* x = (su3_vector*)t_dest;
      su3_vector mAx; //minus Ax
      scalar_mult_add_su3_vector( &DDx[i], &x[i], -4*mass*mass,
          &mAx );
      su3_vector resid;
      add_su3_vector( &mAx, &t_src[i], &resid );      
      gpu_rsq += (double)magsq_su3vec( &resid);
    } END_LOOP
    g_doublesum(&gpu_rsq);
    printf("GPU: requested relative residual=%g, actual relative residual = %g, parity=%d\n", 
        inv_param.tol, sqrt(gpu_rsq/src_rsq), parity);

    //compute residue for CPU solution
    su3_vector* sol = t_dest_check;
    dslash_fn_field(sol, spinorCheck, other_parity, fn);
    dslash_fn_field(spinorCheck, spinorCheck, parity, fn);
    FORSOMEPARITYDOMAIN(i,s,parity){
      su3_vector* DDx = (su3_vector*)spinorCheck;
      su3_vector* x = (su3_vector*)sol;
      su3_vector mAx; //minus Ax
      scalar_mult_add_su3_vector( &DDx[i], &x[i], -4*mass*mass,
          &mAx );
      su3_vector resid;
      add_su3_vector( &mAx, &t_src[i], &resid );      
      cpu_rsq += (double)magsq_su3vec( &resid);
    } END_LOOP
    g_doublesum(&cpu_rsq);
    printf("CPU: requested relative residue=%g, actual relative residue = %g, parity=%d\n", 
        inv_param.tol, sqrt(cpu_rsq/src_rsq),parity);

    FORSOMEPARITY(i,s, qic->parity){
      if ( fabs( t_dest[i].c[0].real - t_dest_check[i].c[0].real) > DELTA
          || fabs( t_dest[i].c[0].imag - t_dest_check[i].c[0].imag) > DELTA
          || fabs( t_dest[i].c[1].real - t_dest_check[i].c[1].real) > DELTA
          || fabs( t_dest[i].c[1].imag - t_dest_check[i].c[1].imag) > DELTA
          || fabs( t_dest[i].c[2].real - t_dest_check[i].c[2].real) > DELTA
          || fabs( t_dest[i].c[2].imag - t_dest_check[i].c[2].imag) > DELTA){
        fprintf(stderr, "ERROR: in function %s, %dth data not match\n", __FUNCTION__ , i);
        fprintf(stderr, "qic->parity =%d\n", qic->parity);
        printf("CPU: dest[%d]      =(%f,%f), (%f,%f), (%f,%f)\n", 
            i, 
            t_dest[i].c[0].real, t_dest[i].c[0].imag, 
            t_dest[i].c[1].real, t_dest[i].c[1].imag, 
            t_dest[i].c[2].real, t_dest[i].c[2].imag);
        printf("CPU: dest_check[%d]=(%f,%f), (%f,%f), (%f,%f)\n", 
            i, 
            t_dest_check[i].c[0].real, t_dest_check[i].c[0].imag, 
            t_dest_check[i].c[1].real, t_dest_check[i].c[1].imag, 
            t_dest_check[i].c[2].real, t_dest_check[i].c[2].imag);		

        exit(1);
      }
    }END_LOOP

    free(t_dest_check);
    free(spinorCheck);

  }
#endif


  total_iters += iters;

  for (dir = 0; dir < 4; dir++) {
    free(fatlink[dir]);
    free(longlink[dir]);
  }

  freeGaugeField(&cudaFatLinkPrecise);
  freeGaugeField(&cudaLongLinkPrecise);
  if (gaugeParam.cuda_prec_sloppy != gaugeParam.cuda_prec){
    freeGaugeField(&cudaFatLinkSloppy);
  }else{
    memset(&cudaFatLinkSloppy,0, sizeof(cudaFatLinkSloppy));
  }

  if ( (gaugeParam.cuda_prec_sloppy != gaugeParam.cuda_prec)
      || (link_recon != link_recon_sloppy)){
    freeGaugeField(&cudaLongLinkSloppy);
  }else{
    memset(&cudaLongLinkSloppy,0, sizeof(cudaLongLinkSloppy));
  }
  return iters;
}




