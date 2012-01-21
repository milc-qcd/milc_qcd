/****************** eigen_stuff_PRIMME.c ****************/
/* MIMD version 7 */

#include "arb_ov_includes.h"

#ifdef PRIMME

#include "../include/primme.h"
static int mxv_kalk;

void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme) ;

/* Matrix on vector routine for which to compute the eigenvalues */

//static const Real jd_shift=1.e-6;

void av_ov (void *x, void *y, int *blockSize, primme_params *primme)
{
    site* s;
    int i, j, iblock;
    double* xx;
    Real* yy;
    wilson_vector tmp1[sites_on_node], tmp2[sites_on_node];

    /*this routine gets a number of vectors (stored consequtively) which
     * need to be mutliplied by the matrix */
    
    for (iblock=0;iblock<*blockSize;iblock++)
    {
	/*copy double precsion vector to single precision wilson_vectors */
	if (kind_of_h0==HOVERLAP)
	{
	    xx=((double*) x)+12*iblock*sites_on_node;
	    FORALLSITES(i,s){
		clear_wvec(&tmp1[i]);
		if (chirality_flag==1)
		{
		    yy= &(tmp1[i].d[0].c[0].real);
		} else {
		    yy= &(tmp1[i].d[2].c[0].real);
		}
		for(j=0;j<12;j++) *(yy++) = *(xx++)  ;
	    }
	} else {
	    xx=((double*) x)+24*iblock*sites_on_node;
	    FORALLSITES(i,s){
		yy= &(tmp1[i].d[0].c[0].real);
		for(j=0;j<24;j++) *(yy++) = *(xx++)  ;
	    }
	}

	Matrix_Vec_mult(tmp1,tmp2);
	
	/* And copy the result back to a complex vector */
	if (kind_of_h0==HOVERLAP)
	{
	    xx=((double*) y)+12*iblock*sites_on_node;
	    FORALLSITES(i,s){
//		scalar_mult_add_wvec(&tmp2[i],&tmp1[i],jd_shift,&tmp2[i]);
		if (chirality_flag==1)
		{
		    yy= &(tmp2[i].d[0].c[0].real);
		} else {
		    yy= &(tmp2[i].d[2].c[0].real);
		}
		for(j=0;j<12;j++) *(xx++) = *(yy++)  ;
	    }
	} else {
	    xx=((double*) y)+24*iblock*sites_on_node;
	    FORALLSITES(i,s){
		yy= &(tmp2[i].d[0].c[0].real);
		for(j=0;j<24;j++) *(xx++) = *(yy++)  ;
	    }
	}


    }
    


}

int Kalkreuter_PRIMME(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
	Real RelTol, int Nvecs, int MaxIter, 
	int Restart, int Kiters, int parity)
{
    int maxnev=Nvecs;       /* number of eigenvalues to compute*/
    int maxn;
    double * evals , *rnorms;    /*work space*/
    double_complex* evecs;    

    static primme_params primme;

    int ivec,j,  ret, i;
    site* s;
    double *xx;             /*for copying*/
    Real *yy;
    long int nn0=ndelta0;

    mxv_kalk=0;
    
    if (kind_of_h0==HOVERLAP)
    {
	maxn=sites_on_node*12/2;                     /*local size of matrix*/
    } else {
	maxn=sites_on_node*12;                       /*local size of matrix*/
    }


    /*allocate memory, EV finder works in double precision only*/
    evals=malloc(Nvecs*sizeof(double));
    if (evals==NULL) exit(1);
    evecs=malloc(Nvecs*maxn*sizeof(double_complex));
    if (evecs==NULL) exit(1);
    rnorms=malloc(Nvecs*sizeof(double_complex));
    if (rnorms==NULL) exit(1);

    /*copy initial guesses into double precision temporary fields*/
    if (kind_of_h0==HOVERLAP){
    for (ivec=0;ivec<Nvecs;ivec++)
    {	
	xx=&evecs[0].real+12*ivec*sites_on_node;
	FORALLSITES(i,s){
	    if (chirality_flag==1)
	    {
		yy= &(eigVec[ivec][i].d[0].c[0].real);
	    } else {
		yy= &(eigVec[ivec][i].d[2].c[0].real);
	    }
	    for(j=0;j<12;j++) *(xx++) = *(yy++)  ;
	}
    }
    } else {
	for (ivec=0;ivec<Nvecs;ivec++)
	{	
	    xx=&evecs[0].real+24*ivec*sites_on_node;
	    FORALLSITES(i,s){
		yy= &(eigVec[ivec][i].d[0].c[0].real);
		for(j=0;j<24;j++) *(xx++) = *(yy++)  ;
	    }
	}
    }


    /*set the parameters of the EV finder*/
    primme_initialize(&primme);

    primme.n=maxn*number_of_nodes;                       /*global size of matrix*/
    primme.nLocal=maxn;           /*local volume*/
    primme.maxOuterIterations=MaxIter;

    primme.numProcs=number_of_nodes;          
    primme.procID=this_node;
    primme.globalSumDouble=par_GlobalSumDouble;/*the wrapper function to do global sums*/

    primme.matrixMatvec =av_ov;                  /*the matrix on vector product*/


    ret = primme_set_method(DEFAULT_MIN_MATVECS, &primme);
    /*
    if (kind_of_h0==HOVERLAP) primme.printLevel=3;
    else */ primme.printLevel=2;
    
    primme.target=primme_smallest;
    
    primme.eps=Tolerance;
 
    primme.numEvals=maxnev;
    primme.initSize=Nvecs;

    /*
    primme_display_params(primme);
    */
    /* call the actual EV finder*/
    ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

    if (ret!=0) /*check return value */
    {
	node0_printf("Return Value %i !\n",ret);
	fflush(stdout);
	exit(1);
    }

 
    if (HOVERLAP==kind_of_h0)
    {
	for (ivec=0;ivec<Nvecs;ivec++)
	{
	    /* copy Evectors and Evalues in global arrays 
	     * (convert from double to single precision) */
//	    evals[ivec]-=jd_shift;
	    eigVal[ivec]=evals[ivec];
	    xx=(double*)&(evecs[0].real)+12*ivec*sites_on_node;
	    FORALLSITES(i,s){
		if (chirality_flag==1)
		{
		    yy= &(eigVec[ivec][i].d[0].c[0].real);
		} else {
		    yy= &(eigVec[ivec][i].d[2].c[0].real);
		}
		for(j=0;j<12;j++) *(yy++) = *(xx++)  ;
	    }
	}
    } else {

	for (ivec=0;ivec<Nvecs;ivec++)
	{
	    eigVal[ivec]=evals[ivec];
	    /* copy Evectors and Evalues in global arrays 
	     * (convert from double to single precision) */
	    xx=(double*)&(evecs[0].real)+24*ivec*sites_on_node;
	    FORALLSITES(i,s){
		yy= &(eigVec[ivec][i].d[0].c[0].real);
		for(j=0;j<24;j++) *(yy++) = *(xx++)  ;
	    }
	}
    }
    node0_printf("mxv operations for eigenvecs %ld\n",ndelta0-nn0);
    
    free(evals);
    free(evecs);
    free(rnorms);
    primme_Free(&primme);
    mxv_kalk=ndelta0-nn0;
    return mxv_kalk;

}

void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme) 
{
    int i;
    for (i=0;i<*count;i++) *((double*)recvBuf+i)=*((double*)sendBuf+i);
    g_vecdoublesum((double*)recvBuf,*count);

}

#else

/* Stub to allow compilation (but not execution) in case PRIMME is not available */

int Kalkreuter_PRIMME(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
	Real RelTol, int Nvecs, int MaxIter, 
	int Restart, int Kiters, int parity)
{
  node0_printf("Kalkreuter_PRIMME: Requires compilation with the PRIMME package\n");
  terminate(1);

  return 0;
}

#endif
