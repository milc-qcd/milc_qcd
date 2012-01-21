/****** eigen_stuff_helpers.c  ****************************************/
/* For the eigenvalue and eigevector computation routines.
* K.O. 8/99 Started. 
* MIMD version 7
*/

#include "arb_ov_includes.h"	/* definitions files and prototypes */

#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif

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


/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. */
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

