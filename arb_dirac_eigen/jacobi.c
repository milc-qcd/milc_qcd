/**************************** jacobi.c ****************************/
/* MIMD version 6 */
#include <malloc.h>
#include <math.h>
#include "../include/complex.h"
#include "jacobi.h"
#include <stdio.h>


/** #define DEBUG **/

#define sign(a) ( a<0 ? -1.0 : 1.0 ) 


/* Allocate space for an NxN Matrix */
Matrix AllocateMatrix(int N)
{
  Matrix A ;
  register int i;

  A.N = N ;
  A.M = (double_complex **)malloc(N*sizeof(double_complex *));
  A.M[0] = (double_complex *)malloc(N*N*sizeof(double_complex));
  for(i=1;i<N;i++) A.M[i] =  A.M[0] + N*i ;
  return (A) ;
}

/* Deallocate the space of Matrix A */
void deAllocate(Matrix *A)
{
  free(A->M[0]) ; 
  free(A->M) ;
}

/* Returns the Frobenius norm of the matrix A                          *
 * FrobeniusNorm(A) = sqrt(Sum_ij |A_ij|^2)                            */
double FrobeniusNorm(Matrix *A)
{
  register int i,N,j ;
  register double tmp ;
  
  N = A->N ;
  tmp = 0.0 ;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      tmp += A->M[i][j].real*A->M[i][j].real+A->M[i][j].imag*A->M[i][j].imag ;
  return sqrt(tmp)  ;
}

/* OffDiag(A) = sqrt(Sum_ij(i!=j) |A_ij|^2)                            */
double OffDiag(Matrix *A)
{
  register int i,N,j ;
  register double tmp ;
  
  N = A->N ;
  tmp = 0.0 ;
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      {
	tmp+=A->M[i][j].real*A->M[i][j].real+A->M[i][j].imag*A->M[i][j].imag ;
	tmp+=A->M[j][i].real*A->M[j][i].real+A->M[j][i].imag*A->M[j][i].imag ;
      }
  return sqrt(tmp)  ;
}

/* Matrix - Matrix Multiply  C = C + AB  */
void MatrixMult(Matrix *A, Matrix *B, Matrix *C)
{
  register int i,j,k,N;
  register double Aik_r, Aik_i, Bkj_r, Bkj_i ;
  
  N=A->N ;
  for(i=0;i<N;i++)
    for(k=0;k<N;k++)
      {
	Aik_r = A->M[i][k].real ;
	Aik_i = A->M[i][k].imag ;
	for(j=0;j<N;j++)
	  {
	    Bkj_r=B->M[k][j].real ;
	    Bkj_i=B->M[k][j].imag ;
	    C->M[i][j].real += Aik_r*Bkj_r - Aik_i*Bkj_i ;
	    C->M[i][j].imag += Aik_r*Bkj_i + Aik_i*Bkj_r ;
	  }
      }
}

/* Copy A to B */
void MatrixCopy(Matrix *A, Matrix *B)
{
  register int i,j,N;
  
  N=A->N ;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      {
	B->M[i][j].real = A->M[i][j].real ;
	B->M[i][j].imag = A->M[i][j].imag ;
      }
}

/* Sets Matrix A to zero */
void ZeroMatrix(Matrix *A)
{
  register int i,j,N;
   N=A->N ;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      {
	A->M[i][j].real = 0.0 ;
	A->M[i][j].imag = 0.0 ;
      }
}
/* Sets Matrix A to Unity */
void UnitMatrix(Matrix *A)
{
  register int i,j,N;
  
   N=A->N ;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      {
	A->M[i][j].real = 1.0 ;
	A->M[i][j].imag = 0.0 ;
      }
}

/* Sets Matrix B to the Hermitian conjugate of A */
void HermitianConj(Matrix *A, Matrix *B)
{
  register int i,j,N;
  
   N=A->N ;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      {
	B->M[i][j].real =  A->M[j][i].real ;
	B->M[i][j].imag = -A->M[j][i].imag ;
      }
}

/* Computes the Schur decomposition of the hermitian matrix A so that  *
 * the (p,q) element is set to zero. p must be less that q.            *
 * (see Golub & Loan p. 446)                                           */
void HermitianSchur(Matrix *A, int p, int q, h_schur *schur)
{
  register double b, t, off, c, s, Apq_r, Apq_i ;
  
  schur->p = p ;
  schur->q = q ;
  Apq_r = A->M[p][q].real ;
  Apq_i = A->M[p][q].imag ;
  off = sqrt(Apq_r*Apq_r + Apq_i*Apq_i) ;
  if( off > DEPS_SCHUR )
    {
      off = 1.0/off ;
      b = .5*(A->M[q][q].real - A->M[p][p].real)*off ;
      t = sign(b)/(fabs(b) + sqrt(b*b+1)) ;
      schur->c = c = 1.0/sqrt(1+t*t) ;
      s = c*t*off ;
      schur->s.real = Apq_r * s;
      schur->s.imag = Apq_i * s;
    }
  else
    {
      schur->c = 1.0 ;
      schur->s.real = schur->s.imag = 0.0 ;
    }
}

/*   Multiplies from the Right the matrix A with the Shur rotation.   *
 *  A<-A*S                                                            */ 
void rightMultiplySchur(Matrix *A, h_schur *schur)
{
  register int i,N,p,q ;
  register double c, s_r,s_i ;
  register double Tp_r, Tp_i, Tq_r, Tq_i ;

  p  = schur->p ;
  q  = schur->q ;
  c  = schur->c ;
  s_r  = schur->s.real ;
  s_i  = schur->s.imag ;

  N = A->N ;

  for(i=0;i<N;i++)
    {
      Tp_r = A->M[i][p].real ; Tp_i = A->M[i][p].imag ; 
      Tq_r = A->M[i][q].real ; Tq_i = A->M[i][q].imag ; 
      
      A->M[i][p].real = Tp_r*c  - Tq_r*s_r - Tq_i*s_i ;
      A->M[i][p].imag = Tp_i*c  + Tq_r*s_i - Tq_i*s_r ;

      A->M[i][q].real = Tq_r*c  + Tp_r*s_r - Tp_i*s_i ;
      A->M[i][q].imag = Tq_i*c  + Tp_r*s_i + Tp_i*s_r ;
    }   
}

/*   Multiplies form the Left the matrix A with the Shur rotation.    *
 *  A<-S'*A  S' is the hermitian conjugate of S.                      */ 
void leftMultiplySchur(Matrix *A, h_schur *schur)
{
  register int i,N,p,q ;
  register double c, s_r,s_i ;
  register double Tp_r, Tp_i, Tq_r, Tq_i ;

  p  = schur->p ;
  q  = schur->q ;
  c  = schur->c ;
  s_r  = schur->s.real ;
  s_i  = schur->s.imag ;

  N = A->N ;
  
  for(i=0;i<N;i++)
    {
      Tp_r = A->M[p][i].real ; Tp_i = A->M[p][i].imag ; 
      Tq_r = A->M[q][i].real ; Tq_i = A->M[q][i].imag ; 
      
      A->M[p][i].real = Tp_r*c  - Tq_r*s_r + Tq_i*s_i ;
      A->M[p][i].imag = Tp_i*c  - Tq_r*s_i - Tq_i*s_r ;

      A->M[q][i].real = Tq_r*c  + Tp_r*s_r + Tp_i*s_i ;
      A->M[q][i].imag = Tq_i*c  - Tp_r*s_i + Tp_i*s_r ;
    }
}

/* Returns the Matrix A diagonalized. The Matrix V is such that        *
 * V'*A*V is diagonal                                                  */
void Jacobi(Matrix *A, Matrix *V, double tolerance)
{
  register int i,j,N ;
  h_schur S ;
  int itmax=1000;
  int iter = 0 ;
  register double eps,off,old_off ;
  N = A->N ;
  /* Initialize V to the unit matrix */
  for(i=0;i<N;i++)
    {
      V->M[i][i].real = 1.0 ;
      V->M[i][i].imag = 0.0 ;
      for(j=i+1;j<N;j++)
	V->M[i][j].real=V->M[i][j].imag=V->M[j][i].real=V->M[j][i].imag = 0.0 ;
    }
  
  eps = FrobeniusNorm(A)*tolerance ;
#ifdef DEBUG
  printf("In jacobi::Jacobi -- convergence eps: %g\n",eps) ;
  printf("%i off(A): %g\n", iter, OffDiag(A) ) ;
#endif
  old_off  = 1.0e+32 ;
  off = OffDiag(A) ;
  while((off>eps)&&(fabs(old_off-off)>eps)&&(iter<itmax))
    {
      iter++ ;
      for(i=0;i<N;i++)
	for(j=i+1;j<N;j++)
	  {
	    HermitianSchur(A, i, j, &S) ;
	    rightMultiplySchur(A, &S) ;
	    leftMultiplySchur(A, &S) ;
	    rightMultiplySchur(V, &S) ;
	  }
      old_off = off ;
      off = OffDiag(A) ;
#ifdef DEBUG
      printf("%i off(A): %g\n", iter, off ) ;
#endif
    }
if(iter==itmax)printf("Jacobi didn't converge in %d steps\n",iter);
}

 
/*
 Here is the code. The change is trivial. I just rearange the real part
of A->M[i][i] since the imaginary should be zero.
*/
/* Reorders eigenvectors so that they come with eigenvalues in      *     
 * ascending order. M->A[i][i].real is also rearanged.              */
void sort_eigenvectors(Matrix *A, Matrix *V)
{
  register int N , i, j, k ;
  register double fu ;
  register double_complex *bu ;
  double_complex *SaveTmp0 ;
  Matrix Tmp ;
  double *D;
  
  N = A->N ;
  Tmp = AllocateMatrix(N) ;
  /* Make V->M[i] the i th eigenvector.                                   *
   * Otherwise V->M[][i] is the eigenvector which makes reordering costly */
  HermitianConj(V, &Tmp) ;

  /* extract the eigenvalues */
  D = (double *) malloc(N*sizeof(double)) ;
  for(i=0;i<N;i++)
    D[i] = A->M[i][i].real ;
  
  /* Save Tmp.M[0] so that we can deallocate the memory Tmp occuppies  *
   * The Value of Tmp.M[0] may change with sorting                     */
  SaveTmp0 = Tmp.M[0] ;
  for(i=0;i<N;i++)
    for (j=1;j<N-i;j++)
      {
	k = j-1 ; 
	if(D[j]<D[k]) /* swicth them */
	  {
	    /* switch the eigenvalues */
	    fu = D[k] ;
	    D[k] = D[j] ;
	    D[j] = fu ;
	    /* switch the eigenvectors */
	    bu = Tmp.M[k] ;
	    Tmp.M[k] = Tmp.M[j] ;
	    Tmp.M[j] = bu ;
	  }
      }
  for(i=0;i<N;i++)
    A->M[i][i].real = D[i];
  free(D) ;
  HermitianConj(&Tmp,V) ;
  /*Restore the Tmp.M[0] for deallocation */
  Tmp.M[0] = SaveTmp0 ;
  deAllocate(&Tmp) ;
}

