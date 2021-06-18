/********************** contraction_cpu.c *********************************/
/* MIMD version 7 */

/* Taken from ks_meson_mom.c to emulate a call to QUDA */
/* CPU version of the GPU KS contraction code */

#include "generic_ks_includes.h"
#include <string.h>
#include "../include/openmp_defs.h"
#ifdef OMP
#include <omp.h>
#endif
#include <assert.h>

#include "../include/static_cast.h"

/* Calculate FT weight factor */

#include <limits.h>
#define QUDA_INVALID_ENUM INT_MIN

// enum_quda.h  describes corr_parity
typedef enum QudaFFTSymmType_t {
  QUDA_FFT_SYMM_ODD  = 1,  // sin(phase)
  QUDA_FFT_SYMM_EVEN = 2,  // cos(phase)
  QUDA_FFT_SYMM_EO   = 3,  // exp(-i phase)
  QUDA_FFT_SYMM_INVALID = QUDA_INVALID_ENUM
} QudaFFTSymmType;

/*******************************************/
typedef struct {
  double_complex* meson_q;     /* cache aligned thread local storage. order meson_q[k*nt+t] */
  void*    alloc_base;  /* base address of this allocation */
} meson_storage_t;

static meson_storage_t*
create_meson_q_thread(int nt, int max_threads, int num_corr_mom){
  char myname[] = "create_meson_q_thread";
  
  meson_storage_t* threadstore = static_cast(meson_storage_t*,malloc(max_threads*sizeof(meson_storage_t)));
    if(threadstore == NULL){
      printf("%s(%d): No room for meson_q_thread array\n",myname,this_node);
    }
  size_t allocsz = nt*num_corr_mom*sizeof(double_complex);
  size_t align = 128; /* bytes; cache line is 64b on x86_64 and 128b on ppc64 */
  allocsz += align; // padding
  for(int mythread=0; mythread<max_threads; mythread++) {
    threadstore[mythread].alloc_base = malloc(allocsz);
    //printf("threadstore[%d].alloc_base = %p [%lu]\n", mythread, threadstore[mythread].alloc_base, allocsz);
    if(threadstore[mythread].alloc_base == NULL){
      printf("%s(%d): No room for meson_q_thread array\n",myname,this_node);
    }
    off_t offset = align - static_cast(size_t,threadstore[mythread].alloc_base) % align;
    threadstore[mythread].meson_q = static_cast(double_complex*,threadstore[mythread].alloc_base + offset);
    assert(static_cast(size_t,threadstore[mythread].meson_q) % align == 0);
    //printf("threadstore[%d].meson_q = %p\n", mythread, threadstore[mythread].meson_q);
  }

  /* first touch initialization by owning thread */
  #pragma omp parallel
  {
    #ifdef OMP
    int mythread = omp_get_thread_num();
    #else
    int mythread = 0;
    #endif
    //printf("thread %d touch %p\n",mythread,threadstore[mythread].meson_q);

    for(int j=0; j<nt*num_corr_mom; ++j)
      {
	threadstore[mythread].meson_q[j].real = 0.;
	threadstore[mythread].meson_q[j].imag = 0.;
      }
  }
  return threadstore;
}

/*******************************************/
static void
destroy_meson_q_thread(meson_storage_t* threadstore, int max_threads){
  if(threadstore == NULL)return;
  for(int mythread=0; mythread<max_threads; mythread++) {
      if(threadstore[mythread].alloc_base != NULL)
	free(threadstore[mythread].alloc_base);
  }
  free(threadstore);
}

/*******************************************/
static Real
sum_meson_q(double_complex *meson_q, meson_storage_t* threadstore, int nonzero[],
	    int max_threads, int nt, int num_corr_mom){

  for(int mythread=0; mythread<max_threads; mythread++) {
    for(int t = 0; t < nt; t++)if(nonzero[t]){
	for(int k=0; k<num_corr_mom; k++)
	  {
	    int idx = nt*k + t;
	    meson_q[idx].real += threadstore[mythread].meson_q[idx].real;
	    meson_q[idx].imag += threadstore[mythread].meson_q[idx].imag;
	    threadstore[mythread].meson_q[idx].real = 0.;
	    threadstore[mythread].meson_q[idx].imag = 0.; // Prevent re-add
	  }
      }
  }

  Real flops = (Real)sites_on_node*8*num_corr_mom;
  return flops;
}
      
/*******************************************/
/* Calculate a single Fourier phase factor */
static complex ff(Real theta, QudaFFTSymmType parity, complex tmp)
{
  complex z; // = {0.,0.};
  
  if(parity == QUDA_FFT_SYMM_EVEN){
    Real costh = cos(theta);
    z.real = tmp.real*costh;
    z.imag = tmp.imag*costh;
  }
  else if(parity == QUDA_FFT_SYMM_ODD){
    Real sinth = sin(theta);
    z.real = -tmp.imag*sinth;
    z.imag =  tmp.real*sinth;
  }
  else if(parity == QUDA_FFT_SYMM_EO){
    Real costh = cos(theta);
    Real sinth = sin(theta);
    z.real = tmp.real*costh-tmp.imag*sinth;
    z.imag = tmp.imag*costh+tmp.real*sinth;
  }
  else{
    printf("ff(%d): bad parity %d\n", this_node, parity);
    terminate(1);
  }
  return z;
} /* ff */

#if 0 //UNUSED
/*******************************************/
/* Create a table of Fourier phases, one for each momentum for each site */

static complex *
create_ftfact(int nx, int ny, int nz, int nt, int num_corr_mom,
	      int **corr_mom, char **corr_parity, int *r0, Real *flops){

  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 

  complex *ftfact = (complex *)malloc(num_corr_mom*sites_on_node*sizeof(complex));
  if(ftfact == NULL)
    {
      printf("(%d): No room for FT phases\n",this_node);
      terminate(1);
    }
  
  /* ftfact contains factors such as cos(kx*x)*sin(ky*y)*exp(ikz*z)
     with factors of cos, sin, and exp selected according to the
     requested component parity */
  
  int i;
  site *s;
  FORALLSITES_OMP(i,s,) {
    for(int k=0; k<num_corr_mom; k++)
      {
	int px = corr_mom[k][0];
	int py = corr_mom[k][1];
	int pz = corr_mom[k][2];
	
	char ex = corr_parity[k][0];
	char ey = corr_parity[k][1];
	char ez = corr_parity[k][2];

	complex tmp;
	
	tmp.real = 1.;
	tmp.imag = 0.;
	
	tmp = ff(factx*(s->x-r0[0])*px, ex, tmp);
	tmp = ff(facty*(s->y-r0[1])*py, ey, tmp);
	tmp = ff(factz*(s->z-r0[2])*pz, ez, tmp);
	
	ftfact[k+num_corr_mom*i] = tmp;
      }
  } END_LOOP_OMP;
  
  *flops += (Real)sites_on_node*18*num_corr_mom;
  return ftfact;
}
#endif //UNUSED
  
/*******************************************/

static void
destroy_ftfact(complex *ftfact ){
  if(ftfact != NULL)
    free(ftfact);
}

/*******************************************/
/* Put this in an appropriate header */

typedef struct {
  int num_corr_mom;  /* Number of sink momenta */
  int *corr_mom;  /* List of four component momenta modes corr_mom[mode,dir=0..3]  */
  QudaFFTSymmType *corr_parity; /* The "parity" of the FT component corr_parity[mode,dir=0..3] */
  int *r0;  /* The coordinate origin for the Fourier phases */
  Real flops; /* Return value */
  Real dtime; /* Return value */
} QudaContractArgs_t;

void qudaContractFT(int milc_precision,
		  QudaContractArgs_t *cont_args,
		  su3_vector *antiquark,  /* Color vector field (propagator) */
		  su3_vector *quark,   /* Color vector field (propagator) */
                  double_complex meson_q[]  /* Resulting hadron correlator indexed by time and momentum: idx=nt*k+t */
		  )
{

  Real dtime = -dclock();
  char myname[] = "qudaContractFT";
  Real flops = 0;

  int num_corr_mom = cont_args->num_corr_mom;
  int *corr_mom = cont_args->corr_mom;
  QudaFFTSymmType *corr_parity = cont_args->corr_parity;
  int *r0 = cont_args->r0;

  node0_printf("CPU contraction code 'qudaContractFT'\n");

#ifdef OMP
  /* max_threads=getenv("OMP_NUM_THREADS"); */
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 1;
#endif

  /* Working space for threaded time-slice reductions */
  meson_storage_t* threadstore = create_meson_q_thread(nt, max_threads, num_corr_mom);

  /* Fourier factors for FT */
  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 

  /* For avoiding adding unecessary zeros */
  int *nonzero = (int *)malloc(nt*sizeof(int));
  if(nonzero == NULL){
    printf("%s(%d): No room for nonzero array\n",myname,this_node);
    terminate(1);
  }
  
  for(int t = 0; t < nt; t++)nonzero[t] = 0;
  
  /* Do FT on "meson" for momentum projection - 
     Result in meson_q.  We use a dumb FT because there 
     are usually very few momenta needed. */

  int i;
  site *s;
  FORALLSITES_OMP(i,s,) {
#ifdef OMP
    int mythread=omp_get_thread_num();
#else
    int mythread=0;
#endif

    /* Color-vector inner product on site i */
    complex meson = su3_dot(antiquark+i, quark+i);

    /* The time coordinate for this site */
    int st = s->t;
    double real = meson.real;
    double imag = meson.imag;
    nonzero[st] = 1;

    double_complex* meson_q_thread = threadstore[mythread].meson_q;

    /* Each thread accumulates its own time-slice values in meson_q_thread
       Each thread works with all of the momenta */

    for(int k=0; k<num_corr_mom; k++)
      {
	/* compute Fourier phase */
	int px = corr_mom[4*k+0];
	int py = corr_mom[4*k+1];
	int pz = corr_mom[4*k+2];
	char ex = corr_parity[4*k+0];
	char ey = corr_parity[4*k+1];
	char ez = corr_parity[4*k+2];
	complex fourier_fact; fourier_fact.real=1.0; fourier_fact.imag=0.0;
	fourier_fact = ff(factx*(s->x-r0[0])*px, ex, fourier_fact);
	fourier_fact = ff(facty*(s->y-r0[1])*py, ey, fourier_fact);
	fourier_fact = ff(factz*(s->z-r0[2])*pz, ez, fourier_fact);

	int idx = k*nt + st;
	meson_q_thread[idx].real += 
	  real*fourier_fact.real - imag*fourier_fact.imag;
	meson_q_thread[idx].imag += 
	  real*fourier_fact.imag + imag*fourier_fact.real;
      }
  } END_LOOP_OMP;
  flops += (Real)sites_on_node*18*num_corr_mom; // Fourier phase; does not count sin(x), cos(x)?
  flops += (Real)num_corr_mom*8*sites_on_node; // contraction

  /* sum meson_q over all the threads */
  flops += sum_meson_q(meson_q, threadstore, nonzero,
		       max_threads, nt, num_corr_mom);
  
  destroy_meson_q_thread(threadstore, max_threads);

  //TODO: do global reduction sum on meson_q. NOTE: must remove global reductions from e.g. ks_spectrum/spectrum_ks.c
  g_veccomplexsum(meson_q, nt*num_corr_mom);

  dtime += dclock();
  cont_args->dtime = dtime;
  cont_args->flops = flops;
}
