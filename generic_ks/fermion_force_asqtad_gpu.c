/****************** fermion_force_asqtad_gpu.c *******************/
/* MIMD version 7 */

#if 0


//headers from quda library
#include <quda.h>
#include <enum_quda.h>
#include <misc.h>
#include <dslash_reference.h>
#include <gauge_quda.h>
#include <hw_quda.h>

#include "generic_ks_includes.h"

#ifdef LOOPEND
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#else
#define END_LOOP        /* define it to be nothing */
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#define MAX_NUM_PATHS 64    

extern void fermion_force_reference(float eps, float weight1, float weight2, 
				    void* act_path_coeff, void* temp_x, void* sitelink, void* mom);

static void
print_anti_hermitmat(anti_hermitmat* a)
{
    printf("(%f, %f) (%f, %f) (%f, %f)\n", a->m01.real, a->m01.imag, 
	   a->m02.real,a->m02.imag, a->m12.real, a->m12.imag);
    printf("(%f, %f) (%f, %f)\n", a->m00im, a->m11im, a->m22im, a->space);
    
}
static void
anti_hermitmat_sanity_check(anti_hermitmat* a)
{
    if ( (a->m01.real  != a->m01.real)
	 ||(a->m01.imag  != a->m01.imag)
	 ||(a->m02.real != a->m02.real)
	 ||(a->m02.imag != a->m02.imag)
	 ||(a->m12.real != a->m12.real)
	 ||(a->m12.imag != a->m12.imag)
	 ||(a->m00im != a->m00im)
	 ||(a->m11im != a->m11im)
	 ||(a->m22im != a->m22im)
	 ||(a->space != a->space) ){
	
	printf("ERROR: nan encounted in anti_hermitmat\n");
	printf("the anti_hermitmat is\n");
	print_anti_hermitmat(a);
	exit(1);
    }
    
    return;
}

static int 
compare_anti_hermitmat(anti_hermitmat* a, anti_hermitmat* b, double tol)
{
    
    anti_hermitmat_sanity_check(a);
    anti_hermitmat_sanity_check(b);

#define ERROUT if (diff > tol){			\
	return 1;				\
    }
    
    double diff;

    diff = fabs(a->m01.real - b->m01.real); ERROUT;        
    diff = fabs(a->m01.imag - b->m01.imag); ERROUT;

    diff = fabs(a->m02.real - b->m02.real); ERROUT;
    diff = fabs(a->m02.imag - b->m02.imag); ERROUT;
    
    diff = fabs(a->m12.real - b->m12.real); ERROUT;
    diff = fabs(a->m12.imag - b->m12.imag); ERROUT;

    diff = fabs(a->m00im - b->m00im); ERROUT;
    diff = fabs(a->m11im - b->m11im); ERROUT;
    diff = fabs(a->m22im - b->m22im); ERROUT;

    diff = fabs(a->space - b->space); ERROUT;    
    
    return 0;
    
}


static QudaGaugeParam gaugeParam;
static FullGauge cudaSiteLink;
static FullMom cudaMom;
static FullHw cudaHw;
void
eo_fermion_force_twoterms_field_gpu( Real eps, Real weight1, Real weight2,
				     half_wilson_vector *temp_x, int prec,
				     fermion_links_t *fn,
				     ks_action_paths *ap)
{
    int i;

    site* st;
    
    quda_set_verbose(0);
#if (PRECISION==1)
    QudaReconstructType link_recon = QUDA_RECONSTRUCT_12;
    QudaPrecision  link_prec = QUDA_SINGLE_PRECISION;
    QudaPrecision cpu_prec = QUDA_SINGLE_PRECISION;    
    QudaPrecision hw_prec = QUDA_SINGLE_PRECISION;
#else    
    QudaReconstructType link_recon = QUDA_RECONSTRUCT_12;
    QudaPrecision  link_prec = QUDA_DOUBLE_PRECISION;
    QudaPrecision  cpu_prec = QUDA_DOUBLE_PRECISION;    
    printf("ERROR: double precision in %s not supported yet\n", __FUNCTION__);
    exit(1);
#endif

#ifdef FFTIME
  int nflop = 433968;
  double dtime;

  dtime=-dclock();
#endif
  
    gaugeParam.X[0] = nx;
    gaugeParam.X[1] = ny;
    gaugeParam.X[2] = nz;
    gaugeParam.X[3] = nt;
    setDims(gaugeParam.X);
    
    gaugeParam.blockDim = 64;
    gaugeParam.cpu_prec = cpu_prec;
    
    gaugeParam.cuda_prec = link_prec;
    gaugeParam.reconstruct = link_recon;
    gauge_param = &gaugeParam;   

    fermion_force_init_cuda(&gaugeParam);
    
    su3_matrix* sitelink;
    sitelink = malloc(sites_on_node * 4*sizeof(su3_matrix));
    if (sitelink == NULL){
	printf("ERROR: malloc failed for sitelink in function %s\n", __FUNCTION__);
	exit(1);
    }
    FORALLSITES(i, st){
	sitelink[4*i] = st->link[0];
	sitelink[4*i + 1] = st->link[1];
	sitelink[4*i + 2] = st->link[2];
	sitelink[4*i + 3] = st->link[3];	
    }

#if 0
    site_link_sanity_check(sitelink, V, gaugeParam.cpu_prec, &gaugeParam);
#endif
    
    createLinkQuda(&cudaSiteLink, &gaugeParam);
    loadLinkToGPU(cudaSiteLink, sitelink, &gaugeParam);

    anti_hermitmat* mom;
    mom = malloc(sites_on_node * 4* sizeof(anti_hermitmat));
    if (mom == NULL){
	printf("ERROR: malloc failed for mom in function %s\n", __FUNCTION__);
	exit(1);
    }
    FORALLSITES(i, st){
	mom[4*i] = st->mom[0];
	mom[4*i + 1] = st->mom[1];
	mom[4*i + 2] = st->mom[2];
	mom[4*i + 3] = st->mom[3];
    }
    createMomQuda(&cudaMom, &gaugeParam);
    loadMomToGPU(cudaMom, mom, &gaugeParam);      
    
    cudaHw = createHwQuda(gaugeParam.X, hw_prec);
    loadHwToGPU(cudaHw, temp_x, hw_prec);

    //fermion_force_reference(eps, weight1, weight2, ap->act_path_coeff, temp_x, sitelink,mom);

    fermion_force_cuda(eps, weight1, weight2, ap->act_path_coeff, cudaHw, cudaSiteLink, cudaMom, &gaugeParam);
    storeMomToCPU(mom, cudaMom, &gaugeParam);
    
#ifdef GPU_FF_CHECK
    
    eo_fermion_force_twoterms_field_cpu(eps, weight1, weight2, temp_x, prec, fn, ap);
    for(i=0, st=lattice;i < sites_on_node; i++, st++){
	anti_hermitmat* cpumom = (anti_hermitmat *)(&st->mom[0]);
	double tol = 0.00001;

	int k;
	for( k=0;k < 4; k++){
	    if (compare_anti_hermitmat(cpumom + k , mom + 4*i + k, tol) != 0){
		printf("ERROR: CPU and GPU result does not match for momemtum(i=%d, k=%d)\n", i, k);
		printf("CPU result:\n");
		print_anti_hermitmat(cpumom+k);
		printf("GPU result:\n");
		print_anti_hermitmat(mom + 4*i + k);
		exit(1);
	    }
	}

    }
#endif

    FORALLSITES(i, st){
	st->mom[0] = mom[4*i];
	st->mom[1] = mom[4*i + 1];
	st->mom[2] = mom[4*i + 2];
	st->mom[3] = mom[4*i + 3];
    }
    
    free(mom);
    free(sitelink);

    freeLinkQuda(&cudaSiteLink);
    freeMomQuda(&cudaMom);
    freeHwQuda(cudaHw);

#ifdef FFTIME
    dtime += dclock();
    node0_printf("FFTIME(GPU):  time = %e (asqtad3) terms = 2 mflops = %e\n",dtime,
		 (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

    
} 


#endif
