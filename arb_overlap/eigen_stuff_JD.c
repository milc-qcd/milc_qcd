/****************** eigen_stuff_JD.c ****************/
/* MIMD version 7 */

#include "arb_ov_includes.h"
#include "primme.h"
static int mxv_kalk;
void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme) ;

void GramSchmidt(wilson_vector **vector, int Num) ;/* Not used */
void copy_Vector(wilson_vector *src, wilson_vector *res ) ; 
void norm2(wilson_vector *vec, double *norm ); 
void dot_product(wilson_vector *vec1, wilson_vector *vec2, 
		 double_complex *dot) ;
void complex_vec_mult_sub(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2) ;
void complex_vec_mult_add(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2) ;
void double_vec_mult(double *a, wilson_vector *vec1, 
		     wilson_vector *vec2) ;
void double_vec_mult_sub(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2) ;
void double_vec_mult_add(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2) ;
void dax_p_by(double *a, wilson_vector *vec1, double *b, wilson_vector *vec2) ; 
void vec_plus_double_vec_mult(wilson_vector *vec1, double *a, wilson_vector *vec2) ; 

void normalize(wilson_vector *vec) ;
void project_out(wilson_vector *vec, wilson_vector **vector, int Num);
void constructArray(wilson_vector **eigVec,wilson_vector **MeigVec,
 Matrix *A,  int *converged) ;
/* void RotateBasis(wilson_vector **eigVec, Matrix *V) ; */


/* Matrix on vector routine for which to compute the eigenvalues */

/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. */
void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res);


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

int Kalkreuter(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
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

/*  Projects out the *vectors from the  vec. Num is the Number of vectors  *
 * The vectors are assumed to be orthonormal.                              */
void project_out(wilson_vector *vec, wilson_vector **vector, int Num)
{
  register int i ;
  double_complex cc ;
  for(i=Num-1;i>-1;i--)
    {
      dot_product(vector[i], vec, &cc) ;
      complex_vec_mult_sub(&cc, vector[i], vec);
    }
}

/* Copies scr to res */
void copy_Vector(wilson_vector *src, wilson_vector *res)
{
  memcpy((void *)res, (void *)src, sites_on_node*sizeof(wilson_vector)) ;
}

/* Returns the 2-norm of a fermion vector */
void norm2(wilson_vector *vec, double *norm)
{
  double n ;
  register site *s;
  register  int i;
  
  n=0 ; 
  FORALLSITES(i,s){
    n += (double)magsq_wvec(&(vec[i]));
  }
  *norm = n ;
  g_doublesum(norm);
  *norm = sqrt(*norm) ;
}
 
/* Returns the dot product of two fermion vectors */
void dot_product(wilson_vector *vec1, wilson_vector *vec2, 
		   double_complex *dot) 
{
  double re,im ;
  register site *s;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORALLSITES(i,s){
    cc = wvec_dot( &(vec1[i]), &(vec2[i]) );
    re += (double)cc.real ;
    im += (double)cc.imag ;
  }
  dot->real = re ; dot->imag = im ;
  g_dcomplexsum(dot);
}

/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
void complex_vec_mult_sub(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2)
{
  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= -(Real)(cc->real) ; sc.imag= -(Real)(cc->imag) ; 
  FORALLSITES(i,s){
    c_scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]),(&sc),&(vec2[i])) ;
  }
}

/* Returns vec2 = vec2 + cc*vec1   cc is a double complex   */
void complex_vec_mult_add(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2)
{
  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; sc.imag= (Real)(cc->imag) ;
  FORALLSITES(i,s){
    c_scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]), (&sc),&(vec2[i])) ;
  }
}

/* Returns vec2 = vec2 - rr*vec1   rr is a double    */
void double_vec_mult_sub(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2)
{
  register site *s;
  register  int i;
  
  FORALLSITES(i,s){
    scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]), -(Real)*rr, &(vec2[i]));
  }
}

/* Returns vec2 = vec2 + rr*vec1   rr is a double    */
void double_vec_mult_add(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2)
{
  register site *s;
  register  int i;
  
  FORALLSITES(i,s){
    scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]), (Real)*rr, &(vec2[i]));
  }
}

/* Returns vec2 = a*vec1   a is a double vec2 can be vec1*/
void double_vec_mult(double *a, wilson_vector *vec1, 
		     wilson_vector *vec2) 
{
  
  register site *s;
  register  int i;
  
  FORALLSITES(i,s){ 
    scalar_mult_wvec( &(vec1[i]),(Real)*a, &(vec2[i])) ;
  }
}

/* Returns vec2 = vec1 + a*vec2 */
void vec_plus_double_vec_mult(wilson_vector *vec1, double *a, wilson_vector *vec2) 
{
  register site *s;
  register  int i;
  
  FORALLSITES(i,s){ 
    scalar_mult_wvec( &(vec2[i]), (Real)(*a), &(vec2[i]) ) ;
    add_wilson_vector( &(vec1[i]), &(vec2[i]), &(vec2[i]) ) ;
  }
}

/* Returns vec1 = a*vec1 + b*vec2   a,b are double    */
void dax_p_by(double *a, wilson_vector *vec1, double *b, wilson_vector *vec2) 
{
  register site *s;
  register  int i;
  
  FORALLSITES(i,s){
    scalar_mult_wvec(&(vec1[i]), (Real)(*a), &(vec1[i])) ;
    scalar_mult_add_wvec(&(vec1[i]), &(vec2[i]), (Real)(*b), 
			       &(vec1[i])) ;
  }
}

/* normalizes the vecror vec. Work only on parity. */
void normalize(wilson_vector *vec)
{
  double norm ;
  norm2(vec,&norm) ;
  norm = 1.0/norm ;
  double_vec_mult(&norm,vec,vec) ;
}

/* Returns the projected matrix A and the error of each eigenvector */
void constructArray(wilson_vector **eigVec,wilson_vector **MeigVec, 
Matrix *A, int *converged)
{
  int i,j,Nvecs ;
  /*
  wilson_vector *res,*grad ;
  */
  double_complex cc,Aij,Aji ;

  Nvecs = A->N ;
  /*
  res = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  grad = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  */
  for(i=0;i<Nvecs;i++)
    {
      if(converged[i]==0)Matrix_Vec_mult(eigVec[i],MeigVec[i]) ;
      dot_product(MeigVec[i], eigVec[i], &cc) ;
      A->M[i][i].real = cc.real ; A->M[i][i].imag = 0.0 ;
      /*
      copy_Vector(MeigVec[i], grad);
      double_vec_mult_sub(&cc.real, eigVec[i], grad);
      norm2(grad, &err[i]);
      */
      for(j=i+1;j<Nvecs;j++)
	{
	  dot_product(MeigVec[i], eigVec[j], &cc) ;
	  Aij=cc ;
	  CONJG(cc, Aji) ;
	  dot_product(eigVec[j],MeigVec[i], &cc) ;
	  CSUM(Aji,cc) ;
	  CONJG(cc, cc) ;
	  CSUM(Aij,cc) ;
	  CMULREAL(Aij,0.5,A->M[i][j]) ;
	  CMULREAL(Aji,0.5,A->M[j][i]) ;
	}
    }
  /* 
  free(grad) ;
  free(res) ; */
}


void RotateBasis(wilson_vector **eigVec, Matrix *V)
{
  wilson_vector **Tmp ;
  register site *s;
  register  int i,N ;
  int j,k ;
  
  N = V->N ;
  /* Allocate the temporary vectors needed */
  Tmp = (wilson_vector **)malloc(N*sizeof(wilson_vector *));
  for(j=0;j<N;j++)
    Tmp[j]=(wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));

  for(j=0;j<N;j++)
    {
      FORALLSITES(i,s){
	clear_wvec(&Tmp[j][i]) ;
      }
      for(k=0;k<N;k++)
	complex_vec_mult_add(&V->M[k][j],eigVec[k],Tmp[j]) ;
    }
  
  /* Copy rotated basis to the eigVec and free temporaries */
  for(i=0;i<N;i++)
    {
      copy_Vector(Tmp[i], eigVec[i]) ; 
      free(Tmp[i]) ;
    }
  free(Tmp) ;
}

/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. */
void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res)
{  
  register site *s;
  register  int i;
  /* chirality check */
  complex ctmp;
  double_complex cn;
  Real chirality,cd;
  wilson_vector wtmp;

  if(kind_of_h0 == HOVERLAP){
    /* check src for chirality */
    cn=dcmplx(0.0,0.0);
    cd=0.0;
    FORALLSITES(i,s){
      mult_by_gamma(&(src[i]),&wtmp,GAMMAFIVE);
      cd += magsq_wvec(&(src[i]));
      ctmp =  wvec_dot(&(src[i]),&wtmp);
      CSUM(cn,ctmp);
    }
    g_dcomplexsum(&cn);
    g_floatsum(&cd);
    chirality=cn.real/cd;
 if(fabs(chirality) < 0.95) node0_printf("chirality %e\n",chirality);
    if(chirality >0.0) chirality_flag=1;
    if(chirality <0.0) chirality_flag= -1;

    hdelta0_field(src,res,kind_of_h0);
    /*
      node0_printf(" Matrix_Vec_mult: chirality_flag %d\n",chirality_flag);
    */
  }
  else 
  {
      hdelta0_field(src,res,kind_of_h0);
  }

}





