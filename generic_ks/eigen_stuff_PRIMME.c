/****** eigen_stuff_PRIMME.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* This version uses PRIMME (2.1 released on Apr. 4 2017)
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

/* NOTE: The PRIMME release version has clashing definitions for complex functions against the MILC version of complex
   type, "complex.h", so we provide a stripped-down version, "Complex_primme.h" in the include directory.
   This must be checked against future releases of MILC code.
   Also, watch out for the definition of PRIMME_INT = long? <- NO, I replaced long to PRIMME_INT */
#include "../include/primme.h"
#include "../include/dslash_ks_redefine.h"
#include <string.h>

static int mxv;
static imp_ferm_links_t **fn;
static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme, int *ierr);

/************************************************************************/

/* The matrix times vector routine given to PRIMME */
static void ks_mxv(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
		   int *blockSize, struct primme_params *primme, int *ierr){
  site* s;
  int i,j,iblock;
  int maxn;
  int parity = active_parity;
  double* xx;
  Real* yy;
  su3_vector tmp1[sites_on_node], tmp2[sites_on_node];

  mxv++;
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
  //  node0_printf("mxv %d\n",mxv);
  *ierr = 0 ;
}

/* Externally preconditioned matrix times vector routine */
static void ks_cheb_mxv(void *x0, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
		    int *blockSize, struct primme_params *primme, int *ierr){
  site* s;
  double el = ((Real *) (*primme).matrix)[0]; // The left end of the Chebyshev interval.
  double er = ((Real *) (*primme).matrix)[1]; // The right end of the Chebyshev interval.
  int chebD = ((double *) (*primme).matrix)[2]; // Degree of the Chebyshev polynomial.
  double c = (er-el)/2; // Half the length of the Chebyshev interval; also indicates the focal distance of the ellipse in (Saad, 1982)
  double d = (er+el)/2; // The center of the ellipse
  double sigma1 = -c/d, sigmaNew, sigmaOld=sigma1;

  int i, j, maxn, parity = active_parity;
  double* xx; //The assumption is that we use either dprimme(..) or zprimme(...)
  Real* yy;
  su3_vector temp1[sites_on_node], temp2[sites_on_node], x[sites_on_node];

  if(parity == EVENANDODD){
    maxn=sites_on_node*3;                       /*local size of input matrix*/
  }
  else
    maxn=sites_on_node*3/2;                     /*local size of input matrix*/

  for (int iblock=0;iblock<*blockSize;iblock++){
    /* Copy double precsion vector to single precision su3_vectors: temp1=x0_iblock in value*/
    xx=((double*) x0)+2*iblock*maxn;
    FORSOMEPARITY(i,s,parity){
      clearvec(&temp1[i]);
      yy= &(temp1[i].c[0].real);
      for(j=0;j<6;j++) *(yy++) = *(xx++);
    }

    /* n=0 */
    mxv++;    
    Matrix_Vec_mult(temp1,x,parity,fn[0]); // x = D^2 x0
    FORSOMEPARITY(i,s,parity){
      scalar_mult_sub_su3_vector(x+i,temp1+i,d,x+i); // x = x-d*temp1 
      scalar_mult_su3_vector(x+i,sigma1/c,x+i); //x = (sigma1/c)*x
    }

    /* n=1, 2, ... , chebD-1 */
    for(int n=1;n<chebD;n++){
      sigmaNew = 1/(2/sigma1-sigmaOld);
      FORSOMEPARITY(i,s,parity){
	clearvec(&temp2[i]);
	su3vec_copy(x+i,temp2+i); // temp2 = x 
      }
      mxv++;
      Matrix_Vec_mult(temp2,x,parity,fn[0]); //x = D^2 temp2
      FORSOMEPARITY(i,s,parity){
        scalar_mult_sub_su3_vector(x+i,temp2+i,d,x+i); // x = x-d*temp2
        scalar_mult_su3_vector(x+i,2*sigmaNew/c,x+i); // x = (2*sigmaNew/c)*x                                                                                        
        scalar_mult_sub_su3_vector(x+i,temp1+i,sigmaOld*sigmaNew,x+i); //x = x-sigmaNew*sigmaOld*temp1
        su3vec_copy(temp2+i,temp1+i); // temp1 = temp2
      }
      sigmaOld = sigmaNew;
    } /* n */

    /* And copy the result back to a complex vector */
    xx=((double*) y)+2*iblock*maxn;
    FORSOMEPARITY(i,s,parity){
      yy= &(x[i].c[0].real);
      for(j=0;j<6;j++) *(xx++) = *(yy++);
    }
  } /* iblock */
  *ierr=0;
}

/* Preconditioning operator given to PRIMME: y = M^{-1}x where M=K^{-1} is a preconditioning operator*/
static void ks_cheb_pre(void *x0, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
			int *blockSize, struct primme_params *primme, int *ierr);

/*****************************************************************************/
int cheb_PRIMME(su3_vector **eigVec, double *eigVal, Real Tolerance, int Nvecs, int MaxIter, Real *chebInfo, int flag){
  /* flag == 0: no preconditioning
     flag == 1: external Chebyshev preconditioning
     flag == 2: internal Chebyshev preconditioning (not impremented yet)
  */
  int maxn;                    /*local size of matrix*/
  double * evals , *rnorms;    /*work space*/
  double_complex* evecs;    
  double pShifts[1] = {0.2};
  static primme_params primme;

  int j, ret;
  site* s;
  double *xx;		/* for copying */
  Real *yy;		/* for copying */
  int parity = active_parity;

#ifdef EIGTIME
  double dtimec;
#endif

  mxv = 0;
  fn = get_fm_links(fn_links);

  if(parity == EVENANDODD){
    maxn=sites_on_node*3;
  }
  else
    maxn=sites_on_node*3/2;

  /*allocate memory, EV finder works in double precision only*/
  evals=malloc(Nvecs*sizeof(double));
  if (evals==NULL) exit(1);
  evecs=malloc(Nvecs*maxn*sizeof(double_complex));
  if (evecs==NULL) exit(1);
  rnorms=malloc(Nvecs*sizeof(double_complex));
  if (rnorms==NULL) exit(1);

#if 0
  /* Initiallize all the eigenvectors to a random vector and
     convert to double precision temporary fields */
  su3_vector *gr0 = create_v_field();
  for(int i=0;i<Nvecs;i++) {
    grsource_plain_field( gr0, parity);  
    xx = (double*)&(evecs[0].real)+2*i*maxn;
    FORSOMEFIELDPARITY(j,parity){
      yy = &gr0[j].c[0].real;
      for(int k=0;k<6;k++) *(xx++) = *(yy++);
    }
  }
  destroy_v_field(gr0);
#endif

  /*set the parameters of the EV finder*/
  primme_initialize(&primme);

  primme.n=maxn*number_of_nodes;		/* global size of matrix */
  primme.nLocal=maxn;				/* local size of a vector */
  primme.numProcs=number_of_nodes;          
  primme.procID=this_node;
  primme.globalSumReal=par_GlobalSumDouble;	/* the wrapper function to do global sums */

  if(flag==0){
    /* No Chebyshev preconditioning is used. */
    primme.matrixMatvec = ks_mxv;		/* the matrix on vector product */
    primme.target=primme_smallest;
  }else if(flag==1){
    /* External Chebyshev preconditioning is applied. */
    primme.matrixMatvec = ks_cheb_mxv;          /* the matrix on vector product */
    primme.matrix = chebInfo;                   /* Chebyshev interval + degree of polynomial */
    primme.target=primme_largest;
  }else if(flag==2){
    /* Internal Chevyshev preconditioning is applied. */
    // under construction
    primme.correctionParams.precondition = 1;   /* Flag indicating the use of preconditioning */
  }
  ret = primme_set_method(PRIMME_DYNAMIC, &primme);

  primme.printLevel=2;
  primme.eps=Tolerance;
  primme.numEvals=Nvecs;
  //primme.initSize=Nvecs; /* # initial guess vectors stored in evecs */

  /* Optimaized Parameter Setting */
  primme.correctionParams.robustShifts = 1; // led to faster convergence with 0 for tol=1e-8  
  primme.locking=1;
#if 1 /* James and Xiao-Yong's optimal setting */
  primme.minRestartSize=120;  /* relevant if locking != 0 */
  primme.maxBasisSize=192;    /* relevant if locking != 0 */
  primme.maxBlockSize=8;
  primme.restartingParams.maxPrevRetain=2;
#endif

#ifdef EIGTIME
  dtimec = -dclock();
#endif

  /* call the actual EV finder*/
  ret = zprimme(NULL, NULL, NULL, &primme); // On return, intWorksize and realWorkSize contains the size in bytes required for parametes set.
  node0_printf("PRIMME workspace int = %d long int = %ld\n", primme.intWorkSize, primme.realWorkSize); fflush(stdout); // the above line is necessary only for this line.
  //primme_display_params(primme); //<- nodes other than the master one also prints the same info.
  ret = zprimme(evals, (PRIMME_COMPLEX_DOUBLE*)evecs, rnorms, &primme);
  node0_printf("max eval: %lf \n",primme.stats.estimateMaxEVal);

  if (ret!=0){ /*check return value */
    node0_printf("PRIMME: zprimme error (Error Code: %d)\nCall stack:\n",ret);
    /**    primme_PrintStackTrace(primme); NOT SUPPORTED in 2.1**/ 
    fflush(stdout);
    exit(1);
  }

  if(flag==0){
    /* copy Evectors and Evalues in global arrays 
       (convert from double to single precision) */
    for(int i=0;i<Nvecs;i++) {
      eigVal[i]=evals[i];
      xx = (double*)&(evecs[0].real)+2*i*maxn;
      FORSOMEPARITY(j,s,parity){
	yy= &(eigVec[i][j].c[0].real);
      for(int k=0;k<6;k++) *(yy++) = *(xx++);
      }
    }
  }/* construct eigenvalues for the obtained eigenvectors if exernal preconditining is used */
  else if(flag==1) construct_eigen_vals(eigVec, eigVal, Nvecs, parity, fn[0]); 

#ifdef EIGTIME
  dtimec += dclock();
  node0_printf("PRIMME: time = %e ",dtimec);//primme.stats.elapsedTime contains time spent by the call to zprimme(...)
#endif

  node0_printf("mxv operations for eigenvecs %d\n",mxv); //primme.stats.numMatvecs should contains the same info.

  node0_printf("BEGIN RESULTS\n");
  for(int i=0;i<Nvecs;i++){
    node0_printf("Eigenvalue(%i) = %g \n", i,eigVal[i]);
  }

  free(evals);
  free(evecs);
  free(rnorms);
  primme_free(&primme);
  cleanup_Matrix() ;
  return mxv;
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

int cheb_PRIMME(su3_vector **eigVec, double *eigVal, Real Tolerance, int Nvecs, int MaxIter, double *intvl)
{
  node0_printf("chev_PRIMME: Requires compilation with the PRIMME package\n");
  terminate(1);

  return 0;
}

#endif
