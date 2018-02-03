/****** eigen_stuff_PRIMME.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* This version uses PRIMME (UMH Jul 2010e
* MIMD version 7
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of the Kogut-Susskind dslash^2. 
*/

/* Include files */
#include "generic_ks_includes.h"

#ifdef PRIMME

/* NOTE: The PRIMME release version has clashing definitions for complex functions, so we provide
   a stripped-down version.  This must be checked against future releases. 
   Also, watch out for the definition of PRIMME_INT = long? */
#include "../include/primme.h"
#include "../include/dslash_ks_redefine.h"
#include <string.h>

static int mxv_kalk;
static int mxv_precond_kalk;
static ks_eigen_param *my_eigen_param;
static imp_ferm_links_t *my_fn;
static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme, int *ierr) ;

/************************************************************************/
//static void ks_mxv(void *x, void *y, int *blockSize, primme_params *primme){

/* The matrix times vector routine given to PRIMME */
static void ks_mxv(void *x, long *ldx, void *y, long *ldy,
		   int *blockSize, struct primme_params *primme, int *ierr){
  site* s;
  int i,j,iblock;
  int maxn;
  int parity = my_eigen_param->parity;
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

    Matrix_Vec_mult(tmp1,tmp2,my_eigen_param,my_fn);
	
    /* And copy the result back to a complex vector */
    xx=((double*) y)+2*iblock*maxn;
    FORSOMEPARITY(i,s,parity){
      yy= &(tmp2[i].c[0].real);
      for(j=0;j<6;j++) *(xx++) = *(yy++);
    }

  }

  *ierr = 0 ;
}


/* The matrix times vector preconditioning routine given to PRIMME */
static void ks_precond_mxv(void *x, long *ldx, void *y, long *ldy,
			   int *blockSize, struct primme_params *primme, 
			   int *ierr){
  site* s;
  int i,j,iblock;
  int maxn;
  int parity = my_eigen_param->parity;
  double* xx;
  Real* yy;
  su3_vector tmp1[sites_on_node], tmp2[sites_on_node];

  node0_printf("ks_precond_mxv called\n");

  mxv_precond_kalk++;
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

    Precond_Matrix_Vec_mult(tmp1,tmp2,my_eigen_param,my_fn);
	
    /* And copy the result back to a complex vector */
    xx=((double*) y)+2*iblock*maxn;
    FORSOMEPARITY(i,s,parity){
      yy= &(tmp2[i].c[0].real);
      for(j=0;j<6;j++) *(xx++) = *(yy++);
    }

  }

  *ierr = 0 ;
}


/*****************************************************************************/
int ks_eigensolve_PRIMME(su3_vector **eigVec, double *eigVal, 
			 ks_eigen_param *eigen_param, int init){

  my_eigen_param = eigen_param;  /* Save for mxv call-back function */
  int Nvecs   = eigen_param->Nvecs;
  int MaxIter = eigen_param->MaxIter;
  int parity  = eigen_param->parity;
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
  mxv_precond_kalk = 0;
  my_fn = get_fm_links(fn_links)[0];

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
  
  su3_vector *gr0 = create_v_field();
  for(j=0;j<Nvecs;j++) {
    grsource_plain_field( gr0, parity);  
    xx = (double*)&(evecs[0].real)+2*j*maxn;
    FORSOMEFIELDPARITY(i,parity){
      yy = &gr0[i].c[0].real;
      for(k=0;k<6;k++) *(xx++) = *(yy++);
    }
  }
  destroy_v_field(gr0);

  /*set the parameters of the EV finder*/
  primme_initialize(&primme);

  primme.n=maxn*number_of_nodes;		/* global size of matrix */
  primme.nLocal=maxn;				/* local volume */
  primme.maxOuterIterations=MaxIter;

  primme.numProcs=number_of_nodes;          
  primme.procID=this_node;
  primme.globalSumReal=par_GlobalSumDouble;	/* the wrapper function to do global sums */

  primme.matrixMatvec =ks_mxv;			/* the matrix on vector product */

  //  ret = primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
  ret = primme_set_method(PRIMME_DYNAMIC, &primme);

  primme.printLevel=2;
#ifdef MATVEC_PRECOND
  primme.target=primme_largest;
#else
  primme.target=primme_smallest;
#endif
  primme.eps=eigen_param->tol;
  primme.numEvals=maxnev;
  primme.initSize=Nvecs;
#ifdef PRIMME_PRECOND
  primme.applyPreconditioner = ks_precond_mxv;
#endif

  /*
  primme_display_params(primme);
  */

#ifdef EIGTIME
  dtimec = -dclock();
#endif

  /* call the actual EV finder*/
  ret = zprimme(NULL, NULL, NULL, &primme);
  printf("PRIMME workspace int = %d long int = %ld\n", primme.intWorkSize, primme.realWorkSize); fflush(stdout);
  primme_display_params(primme);
  ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

  if (ret!=0){ /*check return value */
    node0_printf("ks_eigensolve_PRIMME: zprimme error\nCall stack:\n");
    /**    primme_PrintStackTrace(primme); NOT SUPPORTED in 2.0**/ 
    fflush(stdout);
    exit(1);
  }

  cleanup_Matrix() ;

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

#ifdef MATVEC_PRECOND
  /* Reset eigenvalues from eigenvectors */
  reset_eigenvalues(eigVec, eigVal, Nvecs, parity, my_fn);
#endif

#ifdef EIGTIME
  dtimec += dclock();
  node0_printf("KAULKREUTER: time = %e iters = %d iters/vec = %e\n",
	   dtimec,total_iters, (double)(total_iters)/Nvecs);
#endif

  node0_printf("mxv operations for eigenvecs %d\n",mxv_kalk);
  node0_printf("mxv precond operations for eigenvecs %d\n",mxv_precond_kalk);
  node0_printf("BEGIN RESULTS\n");
  for(i=0;i<Nvecs;i++){
    node0_printf("Eigenvalue(%i) = %g \n", i,eigVal[i]);
  }

  free(evals);
  free(evecs);
  free(rnorms);
  primme_free(&primme);
  return mxv_kalk;
}

/*****************************************************************************/

static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme, int *ierr) 
{
    int i;
    for (i=0;i<*count;i++) *((double*)recvBuf+i)=*((double*)sendBuf+i);
    g_vecdoublesum((double*)recvBuf,*count);

    *ierr = 0 ;
}

#else

/* Stub to allow compilation (but not execution) in case PRIMME is not available */

 int ks_eigensolve_PRIMME(su3_vector **eigVec, double *eigVal, ks_eigen_param *eigen_param, int init)
{
  node0_printf("ks_eigensolve_PRIMME: Requires compilation with the PRIMME package\n");
  terminate(1);

  return 0;
}

#endif
