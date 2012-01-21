/****** eigen_stuff_PRIMME.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* This version uses PRIMME (UMH Jul 2010e
* MIMD version 7
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of the Kogut-Susskind dslash^2. 
*/


/* If you do not define the USE_DSLASH_SPECIAL then calls to the standard   *
 * dslash are made. If you do define USE_DSLASH_SPECIAL then DSLASH_SPECIAL *
 * is used                                                                  */
#ifdef FN
#define USE_DSLASH_SPECIAL
#endif


/* Include files */
#include "generic_ks_includes.h"

#ifdef PRIMME

#include "../include/primme.h"
#include "../include/dslash_ks_redefine.h"
#include <string.h>

static int mxv_kalk;
static imp_ferm_links_t **fn;
static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme) ;

/************************************************************************/

/* The matrix times vector routine give to PRIMME */
static void ks_mxv(void *x, void *y, int *blockSize, primme_params *primme){
  site* s;
  int i,j,iblock;
  int maxn;
  int parity = active_parity;
  double* xx;
  Real* yy;
  su3_vector tmp1[sites_on_node], tmp2[sites_on_node];

  mxv_kalk++;
  if(parity == EVENANDODD){
    maxn=sites_on_node*3;                       /*local size of matrix*/
  }
  else
    maxn=sites_on_node*3/2;                     /*local size of matrix*/

  /* This routine gets a number of vectors (stored consequtively) which
   * need to be mutliplied by the matrix */
    
  for (iblock=0;iblock<*blockSize;iblock++)
  {
    /* Copy double precsion vector to single precision su3_vectors */
    xx=((double*) x)+2*iblock*maxn;
    FORSOMEPARITY(i,s,parity){
      clearvec(&tmp1[i]);
      yy= &(tmp1[i].c[0].real);
      for(j=0;j<6;j++) *(yy++) = *(xx++);
    }

    Matrix_Vec_mult(tmp1,tmp2,parity,fn[0]);
	
    /* And copy the result back to a complex vector */
    xx=((double*) y)+2*iblock*maxn;
    FORSOMEPARITY(i,s,parity){
      yy= &(tmp2[i].c[0].real);
      for(j=0;j<6;j++) *(xx++) = *(yy++);
    }

  }

}


/*****************************************************************************/
int Kalkreuter_PRIMME(su3_vector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int Kiters ){

  int maxnev=Nvecs;       /* number of eigenvalues to compute*/
  int maxn;
  double * evals , *rnorms;    /*work space*/
  double_complex* evecs;    

  static primme_params primme;

  int i,j,k;
  int ret;
  site* s;
  double *xx;		/* for copying */
  Real *yy;		/* for copying */
  int parity = active_parity;

//  int total_iters=0 ;
//  Matrix Array,V ;
//  register site *s ;
//  register  int i ;
//  su3_vector *vec ;
//  su3_vector **MeigVec ;
//  double max_error = 1.0e+10 ;
//  double *grad, *err ;
//  int iter = 0 ;
//  int *converged ;
//  Real ToleranceG ;
#ifdef EIGTIME
  double dtimec;
#endif

  mxv_kalk = 0;
  fn = get_fm_links(fn_links);

  if(parity == EVENANDODD){
    maxn=sites_on_node*3;			/*local size of matrix*/
  }
  else
    maxn=sites_on_node*3/2;			/*local size of matrix*/

  /*allocate memory, EV finder works in double precision only*/
  evals=malloc(Nvecs*sizeof(double));
  if (evals==NULL) exit(1);
  evecs=malloc(Nvecs*maxn*sizeof(double_complex));
  if (evecs==NULL) exit(1);
  rnorms=malloc(Nvecs*sizeof(double_complex));
  if (rnorms==NULL) exit(1);

  /* Initiallize all the eigenvectors to a random vector and
     convert to double precision temporary fields */
  for(j=0;j<Nvecs;j++) {
    grsource_plain(F_OFFSET(g_rand), parity);  
    xx = (double*)&(evecs[0].real)+2*j*maxn;
    FORSOMEPARITY(i,s,parity){
      yy = &(lattice[i].g_rand.c[0].real);
      for(k=0;k<6;k++) *(xx++) = *(yy++);
    }
  }

  /*set the parameters of the EV finder*/
  primme_initialize(&primme);

  primme.n=maxn*number_of_nodes;		/* global size of matrix */
  primme.nLocal=maxn;				/* local volume */
  primme.maxOuterIterations=MaxIter;

  primme.numProcs=number_of_nodes;          
  primme.procID=this_node;
  primme.globalSumDouble=par_GlobalSumDouble;	/* the wrapper function to do global sums */

  primme.matrixMatvec =ks_mxv;			/* the matrix on vector product */

  ret = primme_set_method(DEFAULT_MIN_MATVECS, &primme);

  primme.printLevel=2;
  primme.target=primme_smallest;
  primme.eps=Tolerance;
  primme.numEvals=maxnev;
  primme.initSize=Nvecs;

  /*
  primme_display_params(primme);
  */

#ifdef EIGTIME
  dtimec = -dclock();
#endif

  /* call the actual EV finder*/
  ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

  if (ret!=0) /*check return value */
  {
    node0_printf("Return Value %i !\n",ret);
    fflush(stdout);
    exit(1);
  }

  /* copy Evectors and Evalues in global arrays 
  * (convert from double to single precision) */
  for(j=0;j<Nvecs;j++) {
    eigVal[j]=evals[j];
    xx = (double*)&(evecs[0].real)+2*j*maxn;
    FORSOMEPARITY(i,s,parity){
      yy= &(eigVec[j][i].c[0].real);
      for(k=0;k<6;k++) *(yy++) = *(xx++);
    }
  }

#ifdef EIGTIME
  dtimec += dclock();
  node0_printf("KAULKRITER: time = %e iters = %d iters/vec = %e\n",
	   dtimec,total_iters, (double)(total_iters)/Nvecs);
#endif

  node0_printf("mxv operations for eigenvecs %d\n",mxv_kalk);

  node0_printf("BEGIN RESULTS\n");
  for(i=0;i<Nvecs;i++){
    node0_printf("Eigenvalue(%i) = %g \n", i,eigVal[i]);
  }

  free(evals);
  free(evecs);
  free(rnorms);
  primme_Free(&primme);
  cleanup_Matrix() ;
  return mxv_kalk;
}

/*****************************************************************************/

static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme) 
{
    int i;
    for (i=0;i<*count;i++) *((double*)recvBuf+i)=*((double*)sendBuf+i);
    g_vecdoublesum((double*)recvBuf,*count);

}

#else

/* Stub to allow compilation (but not execution) in case PRIMME is not available */

int Kalkreuter_PRIMME(su3_vector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int Kiters )
{
  node0_printf("Kalkreuter_PRIMME: Requires compilation with the PRIMME package\n");
  terminate(1);

  return 0;
}

#endif
