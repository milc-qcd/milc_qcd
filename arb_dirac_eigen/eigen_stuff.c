/****** eigen_stuff.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* K.O. 8/99 Started. 
* MIMD version 6
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of ``any'' Hermetian matrix acting on data type wilson_vector 
*/
/**#define DEBUG**/
#define JACOBI_TOL 1.110223e-16
#define MINITER 5

#include "arb_dirac_eig_includes.h"    /* definitions files and prototypes */
#include <string.h>

void GramSchmidt(wilson_vector **vector, int Num, int parity) ;/* Not used */
void copy_Vector(wilson_vector *src, wilson_vector *res ) ; 
void norm2(wilson_vector *vec, double *norm, int parity); 
void dot_product(wilson_vector *vec1, wilson_vector *vec2, 
		 double_complex *dot, int parity) ;
void complex_vec_mult_sub(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2, int parity) ;
void complex_vec_mult_add(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2, int parity) ;
void double_vec_mult(double *a, wilson_vector *vec1, 
		     wilson_vector *vec2, int parity) ;
void double_vec_mult_sub(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2, int parity) ;
void double_vec_mult_add(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2, int parity) ;
void dax_p_by(double *a, wilson_vector *vec1, double *b, wilson_vector *vec2, 
	      int parity) ; 
void vec_plus_double_vec_mult(wilson_vector *vec1, double *a, wilson_vector *vec2, 
			      int parity) ;
void normalize(wilson_vector *vec,int parity) ;
void project_out(wilson_vector *vec, wilson_vector **vector, int Num, int parity);
void constructArray(wilson_vector **eigVec, Matrix *A, double *err, int parity) ;

void mult_spin_pseudoscalar(field_offset src, field_offset dest ) ;

/* Message tags to be used for the Matrix Vector multiplication */
static msg_tag *tags1[16],*tags2[16];
/* Temporary wilson_vector used for the squaring */ 
static wilson_vector *temp ;
/* flag indicating wether to start dslash. */
static int dslash_start = 1 ; /* 1 means start dslash */
/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. For this specific application it is the  *
 * -D_slash^2 of the KS fermions. We only compute on the "parity" sites.     *
 * Where parity can be EVEN, ODD, or ENENANDODD                              */
void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res, int parity)
{  
  register site *s;
  register  int i;
  int otherparity;
  /* store last source so that we know when to reinitialize the message tags */
  static wilson_vector *last_src=NULL ;
  if(dslash_start)
    {
#ifdef noDEBUG
      printf("dslash_start: old temp:%x\n", temp); 
#endif
      temp = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
#ifdef noDEBUG
      printf("dslash_start: new temp:%x\n", temp); 
      fflush(stdout) ;
#endif    
    }

  /*reinitialize the tags is the we have a new source */
  if(last_src != src)
    {
#ifdef noDEBUG
printf("dslash_start: cleaning tags. old src:%x new src:%x\n", last_src,src); 
#endif
/*      if(dslash_start)
	cleanup_gathers(tags1,tags2); */
      dslash_start = 1 ;
      last_src = src ;
    }
  switch(parity)
    {
    case EVEN:
      otherparity = ODD ;
      break ;
    case ODD:
      otherparity = EVEN ;
      break ;
    case EVENANDODD:
      otherparity = EVENANDODD ;
      break ;
    default:
      node0_printf("ERROR: wrong parity in eigen_stuff::Matrix_Vec_mult\n") ;
      terminate(1) ;
    }

  FORSOMEPARITY(i,s,parity){s->chi = src[i];}
  herm_delt(F_OFFSET(chi),F_OFFSET(psi));
  FORSOMEPARITY(i,s,parity){res[i]=s->psi;}

  dslash_start = 0 ;
}
/* Deallocates the tags and the temporaries the Matrix_Vec_mult needs */
void cleanup_Matrix()
{
  if(!dslash_start)
    {
      /* cleanup_gathers(tags1,tags2); 
	 cleanup_dslash_wtemps() ; */
      free(temp) ;
    }
  dslash_start = 1 ;
#ifdef DEBUG
  node0_printf("cleanup_Matrix(): done!\n") ; fflush(stdout) ;
#endif
}

/* Modified Gram-Schmidt orthonormalization pg. 219 Golub & Van Loan     *
 * Num is the number of vectors. vector are the vectors. They get        *
 * overwritten by the orthonormal vectors.                               *
 * parity is the parity on which we work on (EVEN,ODD,ENENANDODD).       */
void GramSchmidt(wilson_vector **vector, int Num, int parity)/* NEVER USED */
{
  register int i,j ;
  double norm ;
  double_complex cc ;

  for(i=0;i<Num;i++)
    {
      norm2(vector[i], &norm, parity);
      norm = 1.0/norm ;
      double_vec_mult(&norm, vector[i], vector[i], parity);
      for(j=i+1;j<Num;j++)
	{
	  dot_product(vector[i], vector[j], &cc, parity) ;
	  complex_vec_mult_sub(&cc, vector[i], vector[j], parity);
	}
    }
}

/*  Projects out the *vectors from the  vec. Num is the Number of vectors  *
 * and parity is the parity on which we work on.                           *
 * The vectors are assumed to be orthonormal.                              */
void project_out(wilson_vector *vec, wilson_vector **vector, int Num, int parity)
{
  register int i ;
  double_complex cc ;
  for(i=Num-1;i>-1;i--)
    {
      dot_product(vector[i], vec, &cc, parity) ;
      complex_vec_mult_sub(&cc, vector[i], vec, parity);
    }
}

/* Copies scr to res */
void copy_Vector(wilson_vector *src, wilson_vector *res)
{
  memcpy((void *)res, (void *)src, sites_on_node*sizeof(wilson_vector)) ;
}

/* Returns the 2-norm of a fermion vector */
void norm2(wilson_vector *vec, double *norm, int parity)
{
  register double n ;
  register site *s;
  register  int i;
  
  n=0 ; 
  FORSOMEPARITY(i,s,parity){
    n += magsq_wvec(&(vec[i]));
  }
  *norm = n ;
  g_doublesum(norm);
  *norm = sqrt(*norm) ;
}
 
/* Returns the dot product of two fermion vectors */
void dot_product(wilson_vector *vec1, wilson_vector *vec2, 
		   double_complex *dot, int parity) 
{
  register double re,im ;
  register site *s;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORSOMEPARITY(i,s,parity){
    cc = wvec_dot( &(vec1[i]), &(vec2[i]) );
    re += cc.real ;
    im += cc.imag ;
  }
  dot->real = re ; dot->imag = im ;
  g_dcomplexsum(dot);
}

/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
void complex_vec_mult_sub(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2, int parity)
{
  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= -(Real)(cc->real) ; sc.imag= -(Real)(cc->imag) ; 
  FORSOMEPARITY(i,s,parity){
    c_scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]),(&sc),&(vec2[i])) ;
  }
}

/* Returns vec2 = vec2 + cc*vec1   cc is a double complex   */
void complex_vec_mult_add(double_complex *cc, wilson_vector *vec1, 
			  wilson_vector *vec2, int parity)
{
  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; sc.imag= (Real)(cc->imag) ;
  FORSOMEPARITY(i,s,parity){
    c_scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]), (&sc),&(vec2[i])) ;
  }
}

/* Returns vec2 = vec2 - rr*vec1   rr is a double    */
void double_vec_mult_sub(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2, int parity)
{
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){
    scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]), -(Real)*rr, &(vec2[i]));
  }
}

/* Returns vec2 = vec2 + rr*vec1   rr is a double    */
void double_vec_mult_add(double *rr, wilson_vector *vec1,  
			 wilson_vector *vec2, int parity)
{
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){
    scalar_mult_add_wvec(&(vec2[i]), &(vec1[i]), (Real)*rr, &(vec2[i]));
  }
}

/* Returns vec2 = a*vec1   a is a double vec2 can be vec1*/
void double_vec_mult(double *a, wilson_vector *vec1, 
		     wilson_vector *vec2, int parity) 
{
  
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){ 
    scalar_mult_wvec( &(vec1[i]),(Real)*a, &(vec2[i])) ;
  }
}

/* Returns vec2 = vec1 + a*vec2 */
void vec_plus_double_vec_mult(wilson_vector *vec1, double *a, wilson_vector *vec2, 
			      int parity)
{
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){ 
    scalar_mult_wvec( &(vec2[i]), *a, &(vec2[i]) ) ;
    add_wilson_vector( &(vec1[i]), &(vec2[i]), &(vec2[i]) ) ;
  }
}

/* Returns vec1 = a*vec1 + b*vec2   a,b are double    */
void dax_p_by(double *a, wilson_vector *vec1, double *b, wilson_vector *vec2, 
	      int parity) 
{
  register site *s;
  register  int i;
  
  FORSOMEPARITY(i,s,parity){
    scalar_mult_wvec(&(vec1[i]), (Real)(*a), &(vec1[i])) ;
    scalar_mult_add_wvec(&(vec1[i]), &(vec2[i]), (Real)(*b), 
			       &(vec1[i])) ;
  }
}

/* normalizes the vecror vec. Work only on parity. */
void normalize(wilson_vector *vec,int parity)
{
  double norm ;
  norm2(vec,&norm,parity) ;
  norm = 1.0/norm ;
  double_vec_mult(&norm,vec,vec,parity) ;
}

int Rayleigh_min(wilson_vector *vec, wilson_vector **eigVec, Real Tolerance, 
		 Real RelTol, int Nvecs, int MaxIter, int Restart, int parity)
{
  int iter ;
  double beta, cos_theta, sin_theta ;
  double quot, P_norm, theta, real_vecMp, pMp, vec_norm ;
  double g_norm, old_g_norm, start_g_norm ;
  double_complex cc ;
  wilson_vector *Mvec, *grad, *P, *MP ;
  
  Mvec     = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;
  grad     = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;
  P        = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;
  MP       = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;

  project_out(vec, eigVec, Nvecs, parity);
  normalize(vec,parity) ; 
  Matrix_Vec_mult(vec,Mvec,parity) ;
  project_out(Mvec, eigVec, Nvecs, parity);
  
  /* Compute the quotient quot=vev*M*vec */
  dot_product(vec, Mvec, &cc, parity) ;
  /* quot is real since M is hermitian. quot = vec*M*vec */
  quot = cc.real ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- quot=%g,%g\n", quot,cc.imag) ;
#endif  
  /* Compute the grad=M*vec - quot*vec */
  copy_Vector(Mvec,grad) ;
  double_vec_mult_sub(&quot, vec, grad, parity) ;
  /* set P (the search direction) equal to grad */
  copy_Vector(grad,P) ;
  /* compute the norms of P and grad */
  norm2(P   , &P_norm, parity) ;
  norm2(grad, &g_norm, parity) ;
  start_g_norm = g_norm ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- g_norm=%g\n", g_norm) ;
#endif  

  iter = 0 ;
  while( (g_norm>Tolerance)&&
	 ( ((iter<MaxIter)&&(g_norm/start_g_norm>RelTol)) || (iter<MINITER) ) 
	 )
    {
      iter++ ;
      Matrix_Vec_mult(P,MP,parity) ;
      dot_product(vec, MP, &cc, parity) ; 
      real_vecMp = cc.real ;
      dot_product(P, MP, &cc, parity) ; 
      pMp = cc.real ; /*pMp is real */
      theta = 0.5*atan(2.0*real_vecMp/(quot*P_norm - pMp/P_norm)) ;
      sin_theta = sin(theta) ;
      cos_theta = cos(theta) ;
      if(sin_theta*cos_theta*real_vecMp>0)
	{ 
	  theta = theta - 0.5*PI ; /* chose the minimum not the maximum */
	  sin_theta = sin(theta) ; /* the sin,cos calls can be avoided */
	  cos_theta = cos(theta) ;
	}
      sin_theta = sin_theta/P_norm ;
      /* vec = cos(theta)*vec +sin(theta)*P/p_norm */
      dax_p_by(&cos_theta, vec,&sin_theta, P,parity) ; 
      /* Mvec = cos(theta)*Mvec +sin(theta)*MP/p_norm */
      dax_p_by(&cos_theta,Mvec,&sin_theta,MP,parity) ; 
      /* renormalize vec ... */
      if(iter%Restart == 0 ) {
#ifdef DEBUG
	node0_printf("Renormalizing...");
	norm2(vec,&vec_norm,parity) ;
	node0_printf("  norm: %g\n",1.0/vec_norm);
#endif
	/* Project vec on the orthogonal complement of eigVec */ 
	project_out(vec, eigVec, Nvecs, parity);
	normalize(vec,parity) ; 
	Matrix_Vec_mult(vec,Mvec,parity);
	/* Recompute the quotient */
	dot_product(vec, Mvec, &cc, parity) ;
	/* quot is real since M is hermitian. quot = vec*M*vec */
	quot = cc.real ;
	/* Recompute the grad */
	copy_Vector(Mvec,grad) ;
	double_vec_mult_sub(&quot, vec, grad, parity) ;
	norm2(grad,&g_norm,parity) ; /* recompute the g_norm */
	/* Project P on the orthogonal complement of eigVec */
	project_out(P, eigVec, Nvecs, parity);
	/* make P orthogonal to vec */
	dot_product(vec, P, &cc, parity) ;
	complex_vec_mult_sub(&cc, vec, P, parity);
	/* make P orthogonal to grad */
	dot_product(grad, P, &cc, parity) ;
	complex_vec_mult_sub(&cc, grad, P, parity);
	norm2(P,&P_norm,parity) ; /* recompute the P_norm */
      }
      dot_product(vec, Mvec, &cc, parity) ;
      /* quot is real since M is hermitian. quot = vec*M*vec */
      quot = cc.real ;
#ifdef DEBUG
      node0_printf("Rayleigh_min: %i, quot=%8g g=%8g b=%6g P:%6g\n",
		   iter,cc.real,g_norm,beta,P_norm) ;    
#endif      
      old_g_norm = g_norm ;
       
      copy_Vector(Mvec,grad) ;
      double_vec_mult_sub(&quot, vec, grad, parity) ;

      norm2(grad,&g_norm,parity) ;
      beta = cos_theta*g_norm*g_norm/(old_g_norm*old_g_norm) ;
      /* Cut off beta */
      if( beta>2.0 ) 
	beta = 2.0 ;
      dot_product(vec, P, &cc, parity) ;
      cc.real *= beta ; cc.imag *= beta ;
      /*      P_norm = 1.0/P_norm ;
	      double_vec_mult(&P_norm,  P,  P, parity) ; */
      vec_plus_double_vec_mult(grad, &beta, P, parity) ;/* P = grad + beta*P */
      complex_vec_mult_sub(&cc, vec, P, parity) ; /* P = P - cc*vec */
      norm2(P   , &P_norm, parity) ;
    } 
  project_out(vec, eigVec, Nvecs, parity);
  normalize(vec,parity) ;
  cleanup_Matrix() ;
  free(MP) ;
  free(P) ;
  free(grad) ;
  free(Mvec) ;
  iter++ ;
  return iter ;
}

/* Returns the projected matrix A and the error of each eigenvector */
void constructArray(wilson_vector **eigVec, Matrix *A, double *err, int parity)
{
  int i,j,Nvecs ;
  wilson_vector *res,*grad ;
  double_complex cc,Aij,Aji ;

  Nvecs = A->N ;
  res = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  grad = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  for(i=0;i<Nvecs;i++)
    {
      Matrix_Vec_mult(eigVec[i],res,parity) ;
      dot_product(res, eigVec[i], &cc, parity) ;
      A->M[i][i].real = cc.real ; A->M[i][i].imag = 0.0 ;
      copy_Vector(res, grad) ;
      double_vec_mult_sub(&cc.real, eigVec[i], grad, parity) ;
      norm2(grad, &err[i], parity) ;
      for(j=i+1;j<Nvecs;j++)
	{
	  dot_product(res, eigVec[j], &cc, parity) ;
	  Aij=cc ;
	  CONJG(cc, Aji) ;
	  dot_product(eigVec[j],res, &cc, parity) ;
	  CSUM(Aji,cc) ;
	  CONJG(cc, cc) ;
	  CSUM(Aij,cc) ;
	  CMULREAL(Aij,0.5,A->M[i][j]) ;
	  CMULREAL(Aji,0.5,A->M[j][i]) ;
	}
    }
  free(grad) ;
  free(res) ;
}


void RotateBasis(wilson_vector **eigVec, Matrix *V, int parity)
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
      FORSOMEPARITY(i,s,parity){
	clear_wvec(&Tmp[j][i]) ;
      }
      for(k=0;k<N;k++)
	complex_vec_mult_add(&V->M[k][j],eigVec[k],Tmp[j],parity) ;
    }
  
  /* Copy rotated basis to the eigVec and free temporaries */
  for(i=0;i<N;i++)
    {
      copy_Vector(Tmp[i], eigVec[i]) ; 
      free(Tmp[i]) ;
    }
  free(Tmp) ;
}

int Kalkreuter(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int Kiters, int parity)
{
  int total_iters=0 ;
  int j,k;
  Matrix Array,V ;
  register site *s;
  register  int i ;
  wilson_vector *vec, *res ;
  double max_error = 1.0e+10 ;
  double *grad, *err ;
  int iter = 0 ;

  /** Allocate the array **/
  Array = AllocateMatrix(Nvecs) ;  
  /** Allocate the Eigenvector matrix **/
  V = AllocateMatrix(Nvecs) ;  

  vec = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  grad = (double *) malloc(Nvecs*sizeof(double)) ;
  err = (double *) malloc(Nvecs*sizeof(double)) ;
		  
  /* Initiallize all the eigenvectors to a random vector 
  for(j=0;j<Nvecs;j++)
    {
      grad[j] = 1.0e+10 ;
      grsource_w(parity);  
      FORSOMEPARITY(i,s,parity){
	copy_wvec(&(s->g_rand),&(eigVec[j][i]));
      }
      eigVal[j] = 1.0e+16 ;
    }
*/
  for(j=0;j<Nvecs;j++){grad[j] = 1.0e+10 ;/* eigVal[j] = 1.0e+16 ;*/}

  while((max_error>Tolerance)&&(iter<Kiters))
    {
      iter++ ;
      for(j=0;j<Nvecs;j++)
	{
	  if(grad[j]>(Tolerance))
	    {
	      copy_Vector(eigVec[j],vec) ;
	      total_iters += Rayleigh_min(vec, eigVec, Tolerance, RelTol, 
					  j, MaxIter , Restart, parity) ;
	      copy_Vector(vec,eigVec[j]) ;
	    }
	}
      constructArray(eigVec, &Array, grad, parity) ;

#ifdef DEBUG
      node0_printf("Eigenvalues before diagonalization\n");
      for(i=0;i<Nvecs;i++)
	node0_printf("quot(%i) = %g |grad|=%g\n",i,Array.M[i][i].real,grad[i]);
#endif

      Jacobi(&Array, &V, JACOBI_TOL) ;
      sort_eigenvectors(&Array,&V) ;
      RotateBasis(eigVec,&V,parity) ;
      constructArray(eigVec, &Array, grad, parity) ;

      /* find the maximum error */
      max_error = 0.0 ;
      for(i=0;i<Nvecs;i++)
        {
          err[i] = eigVal[i] ;
          eigVal[i] = Array.M[i][i].real ;
          err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol) ;
          if(err[i]>max_error) 
            max_error = err[i] ;
        }
      /* old code--see kostas_29_feb.tex
	{
	  if(grad[i]>max_error) 
	    max_error = grad[i] ;
	  err[i] = eigVal[i] ;
	  eigVal[i] = Array.M[i][i].real ;
	  err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol) ;
	}
	*/


      node0_printf("\nEigenvalues after diagonalization at iteration %i\n",
		   iter);
      for(i=0;i<Nvecs;i++)
	node0_printf("quot(%i) = %g +/- %8e |grad|=%g\n",
		     i,eigVal[i],err[i],grad[i]);

    }
  
  for(i=0;i<Nvecs;i++)
    node0_printf("Eigenvalue(%i) = %g +/- %8e \n",i,eigVal[i],err[i]);

  /** Deallocate the arrays **/
  deAllocate(&V) ;
  deAllocate(&Array) ;
  free(err) ;
  free(grad) ;
  free(vec) ;
  return total_iters ;
}




